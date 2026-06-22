# batch_blast_verify.py

import os
import csv
import json
import glob
import math
import time
import datetime
import numpy as np
import requests
import pytz
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

from bioinf_packages.verify_funcs._compute_csv import (
    apply_cross_gene_consistency,
)
from bioinf_packages.verify_funcs._score_hit import (
    pident_to_llr,
    sigmoid,
    softmax_weights,
    score_hit,
    aggregate_hit_scores,
)
from bioinf_packages.verify_funcs._subsamples_to_blast import (
    _submit_blast,
    _parse_blast_xml,
    load_metadata,
    parse_subsample_header,
    _load_existing_csv,
    _append_to_csv,
    CSV_FIELDS,
)

BLAST_URL     = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
POLL_INTERVAL = 30    # seconds between polling rounds
EASTERN       = pytz.timezone("US/Eastern")


# ─────────────────────────────────────────────────────────────────────────────
# NCBI USAGE POLICY HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _is_ncbi_peak_hours() -> bool:
    """Return True if within NCBI peak hours (weekdays 5am-9pm Eastern)."""
    now     = datetime.datetime.now(EASTERN)
    weekday = now.weekday()   # 0=Monday, 6=Sunday
    hour    = now.hour
    if weekday >= 5:
        return False
    return 5 <= hour < 21


def _wait_for_off_peak(check_interval: int = 300) -> None:
    """Block until NCBI off-peak hours, printing status every check_interval seconds."""
    if not _is_ncbi_peak_hours():
        return

    now          = datetime.datetime.now(EASTERN)
    off_peak     = now.replace(hour=21, minute=0, second=0, microsecond=0)
    wait_seconds = (off_peak - now).total_seconds()

    print(f"\n  Current time (Eastern): {now.strftime('%A %H:%M')}")
    print(f"  NCBI peak hours are weekdays 5am-9pm Eastern.")
    print(f"  Off-peak starts at: {off_peak.strftime('%H:%M')} Eastern "
          f"({int(wait_seconds/3600)}h {int((wait_seconds%3600)/60)}m away)")
    print(f"  Waiting for off-peak hours...")

    while _is_ncbi_peak_hours():
        now       = datetime.datetime.now(EASTERN)
        remaining = (off_peak - now).total_seconds()
        print(f"    [{now.strftime('%H:%M')}] Still in peak hours — "
              f"{int(remaining/3600)}h {int((remaining%3600)/60)}m remaining...")
        time.sleep(check_interval)

    print(f"  Off-peak hours reached — proceeding with submission.")


# ─────────────────────────────────────────────────────────────────────────────
# JOB RESULT SCORING
# ─────────────────────────────────────────────────────────────────────────────

def _score_job_results(
        records: list,
        blast_results: dict,
        metadata: dict,
        taxonomy: dict,
        glossary: dict,
        label: str,
        rid: str,
        fasta_file: str,
) -> list:
    """
    Score BLAST results for one submitted job and return CSV rows.

    Parameters
    ----------
    records       : list   SeqRecord objects from the subsample FASTA
    blast_results : dict   query_id -> list of hit dicts from _parse_blast_xml
    metadata      : dict   (accession, gene) -> metadata dict
    taxonomy      : dict   from load_taxonomy()
    glossary      : dict   from load_glossary()
    label         : str    label to write into CSV rows
    rid           : str    BLAST RID for traceability
    fasta_file    : str    filename of the subsample FASTA

    Returns
    -------
    list of dicts, one per hit (or one NO_HITS row per subsample
    with no hits), ready to pass to _append_to_csv()
    """
    new_rows = []

    for record in records:
        sub_info  = parse_subsample_header(record.id)
        accession = sub_info["accession"]
        gene      = sub_info["gene"]
        taxon     = sub_info["taxon"]
        sample    = sub_info["sample"]

        meta_key = (accession, gene)
        parent   = metadata.get(meta_key, {
            "taxon":            taxon,
            "gene":             gene,
            "accession":        accession,
            "orientation":      sub_info.get("orientation", "+"),
            "seq_length":       sub_info.get("parent_length") or 0,
            "subsample_length": (
                (sub_info.get("subsample_end") or 0) -
                (sub_info.get("subsample_start") or 0)
            ),
            "n_samples": 0,
        })

        # Match query ID to BLAST result — BLAST sometimes truncates IDs
        hits_raw = blast_results.get(record.id, [])
        if not hits_raw:
            for blast_id, hits in blast_results.items():
                if (record.id.startswith(blast_id)
                        or blast_id.startswith(record.id[:30])):
                    hits_raw = hits
                    break

        if not hits_raw:
            new_rows.append({
                "label":                 label,
                "taxon":                 taxon,
                "gene":                  gene,
                "accession":             accession,
                "orientation":           parent.get("orientation", "+"),
                "parent_length":         parent.get("seq_length", 0),
                "subsample_length":      parent.get("subsample_length", 0),
                "n_samples_total":       parent.get("n_samples", 0),
                "sample":                sample,
                "subsample_start":       sub_info.get("subsample_start"),
                "subsample_end":         sub_info.get("subsample_end"),
                "hit_accession":         "NO_HITS",
                "hit_description":       "",
                "bitscore":              "",
                "pident":                "",
                "evalue":                "",
                "align_length":          "",
                "species_llr":           "",
                "gene_llr":              "",
                "total_hit_llr":         "",
                "subsample_posterior":   "",
                "subsample_total_llr":   "",
                "n_hits_this_subsample": 0,
                "rid":                   rid,
                "fasta_file":            fasta_file,
            })
            continue

        # Score hits — aggregate_hit_scores() uses no prior/k/beta
        # parameters; those are only needed during calibration scoring.
        # The CSV stores raw species_llr and gene_llr values which
        # calibrate_from_csv() and score_from_csv() use directly.
        scored_hits   = [
            score_hit(hit, taxon, gene, taxonomy, glossary)
            for hit in hits_raw
        ]
        subsample_agg = aggregate_hit_scores(scored_hits)

        for hit, scored in zip(hits_raw, scored_hits):
            new_rows.append({
                "label":                 label,
                "taxon":                 taxon,
                "gene":                  gene,
                "accession":             accession,
                "orientation":           parent.get("orientation", "+"),
                "parent_length":         parent.get("seq_length", 0),
                "subsample_length":      parent.get("subsample_length", 0),
                "n_samples_total":       parent.get("n_samples", 0),
                "sample":                sample,
                "subsample_start":       sub_info.get("subsample_start"),
                "subsample_end":         sub_info.get("subsample_end"),
                "hit_accession":         scored["accession"],
                "hit_description":       hit["description"][:200],
                "bitscore":              round(scored["bitscore"], 2),
                "pident":                round(scored["pident"], 4),
                "evalue":                scored["evalue"],
                "align_length":          scored["align_length"],
                "species_llr":           round(scored["species_llr"], 4),
                "gene_llr":              round(scored["gene_llr"], 4),
                "total_hit_llr":         round(scored["total_llr"], 4),
                "subsample_posterior":   round(
                    subsample_agg["posterior_prob"], 4),
                "subsample_total_llr":   round(
                    subsample_agg["total_llr"], 4),
                "n_hits_this_subsample": len(hits_raw),
                "rid":                   rid,
                "fasta_file":            fasta_file,
            })

    return new_rows


# ─────────────────────────────────────────────────────────────────────────────
# BATCH BLAST SUBMISSION AND POLLING
# ─────────────────────────────────────────────────────────────────────────────

def blast_subsample_directory(
        subsample_dir: str,
        metadata_csv: str,
        taxonomy: dict,
        glossary: dict,
        results_csv: str,
        label: str = "QUERY",
        file_prefix: str = "subsample_batch",
        filter_thresholds_json: str = None,
        evalue: float = 1e-10,
        max_hits: int = 20,
        database: str = "nt",
        entrez_query: str = "Vertebrata[Organism]",
        max_polls: int = 60,
        respect_peak_hours: bool = True,
        submission_interval: int = 30,
) -> dict:
    """
    Submit all subsample FASTA files in a directory to NCBI BLAST,
    accumulating results into a single CSV as jobs complete.

    All files are submitted first (one every submission_interval
    seconds), then results are polled concurrently and written to
    CSV as each job completes.

    Parameters
    ----------
    subsample_dir          : str   directory containing subsample FASTAs
    metadata_csv           : str   metadata CSV from prepare_subsample_fastas
    taxonomy               : dict  from load_taxonomy()
    glossary               : dict  from load_glossary()
    results_csv            : str   path to accumulate results CSV
    label                  : str   label for CSV rows. Default "QUERY".
    file_prefix            : str   FASTA filename prefix to glob
    filter_thresholds_json : str   optional per-gene filter thresholds
    evalue                 : float BLAST e-value threshold
    max_hits               : int   max hits per query sequence
    database               : str   BLAST database
    entrez_query           : str   Entrez filter
    max_polls              : int   max poll attempts per job
    respect_peak_hours     : bool  wait for NCBI off-peak if True
    submission_interval    : int   minimum seconds between submissions

    Returns
    -------
    dict with n_files_submitted, n_files_completed, n_files_failed,
              n_rows_added, results_csv, completed_rids, failed_files
    """
    pattern     = os.path.join(subsample_dir, f"{file_prefix}_*.fasta")
    fasta_files = sorted(glob.glob(pattern))

    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files matching '{pattern}'")

    print(f"\n{'='*60}")
    print(f"BLAST SUBSAMPLE DIRECTORY")
    print(f"{'='*60}")
    print(f"  Directory      : {subsample_dir}")
    print(f"  FASTA files    : {len(fasta_files)}")
    print(f"  Results CSV    : {results_csv}")
    print(f"  Label          : {label}")
    print(f"  Peak hrs check : {respect_peak_hours}")
    print(f"  Submit interval: {submission_interval}s")

    # Load per-gene filter thresholds if provided
    per_gene_thresholds = None
    if filter_thresholds_json and os.path.exists(filter_thresholds_json):
        with open(filter_thresholds_json) as f:
            thresh_data         = json.load(f)
            per_gene_thresholds = thresh_data.get("per_gene_thresholds")
        print(f"  Filter thresholds loaded for "
              f"{len(per_gene_thresholds)} genes")

    metadata                     = load_metadata(metadata_csv)
    existing_keys, existing_rows = _load_existing_csv(results_csv)
    print(f"  Existing CSV rows: {len(existing_rows)}")

    if respect_peak_hours:
        _wait_for_off_peak()

    # ── Phase 1: Submit all jobs ──────────────────────────────────────────────
    print(f"\n── Phase 1: Submitting {len(fasta_files)} jobs ──────────")

    submitted_jobs   = []
    last_submit_time = 0.0

    for i, fasta_file in enumerate(fasta_files, 1):

        if respect_peak_hours and _is_ncbi_peak_hours():
            print(f"\n  Peak hours detected before submission "
                  f"{i}/{len(fasta_files)} — waiting...")
            _wait_for_off_peak()

        elapsed = time.time() - last_submit_time
        if elapsed < submission_interval:
            wait = submission_interval - elapsed
            print(f"  Waiting {wait:.0f}s before next submission...")
            time.sleep(wait)

        records  = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            print(f"  [{i}/{len(fasta_files)}] SKIP "
                  f"{Path(fasta_file).name} — empty")
            continue

        total_bp = sum(len(str(r.seq)) for r in records)
        if total_bp > 900_000:
            print(f"  [{i}/{len(fasta_files)}] SKIP "
                  f"{Path(fasta_file).name} — "
                  f"{total_bp:,}bp exceeds 900k limit")
            continue
        if total_bp > 800_000:
            print(f"  [{i}/{len(fasta_files)}] WARNING "
                  f"{Path(fasta_file).name} — "
                  f"{total_bp:,}bp above 80% limit")

        with open(fasta_file) as f:
            fasta_content = f.read()

        try:
            rid = _submit_blast(
                fasta_content = fasta_content,
                evalue        = evalue,
                max_hits      = max_hits,
                database      = database,
                entrez_query  = entrez_query,
            )
            last_submit_time = time.time()
            submitted_jobs.append({
                "rid":        rid,
                "fasta_file": fasta_file,
                "records":    records,
                "total_bp":   total_bp,
                "submitted":  time.time(),
                "status":     "WAITING",
                "query_gene_map": {
                    r.id: parse_subsample_header(r.id).get("gene")
                    for r in records
                },
            })
            print(f"  [{i}/{len(fasta_files)}] SUBMITTED "
                  f"{Path(fasta_file).name} — "
                  f"RID={rid} ({total_bp:,}bp, {len(records)} seqs)")

        except Exception as e:
            print(f"  [{i}/{len(fasta_files)}] FAILED to submit "
                  f"{Path(fasta_file).name}: {e}")

    print(f"\n  Submitted {len(submitted_jobs)} jobs successfully")

    # ── Phase 2: Poll and save as results arrive ──────────────────────────────
    print(f"\n── Phase 2: Polling for results ─────────────────────────")
    print(f"  Polling every {POLL_INTERVAL}s per job, "
          f"max {max_polls} attempts each")

    completed_rids = []
    failed_files   = []
    n_rows_added   = 0
    poll_counts    = {job["rid"]: 0 for job in submitted_jobs}

    while any(j["status"] == "WAITING" for j in submitted_jobs):

        for job in submitted_jobs:
            if job["status"] != "WAITING":
                continue

            rid = job["rid"]
            poll_counts[rid] += 1

            if poll_counts[rid] > max_polls:
                print(f"  RID {rid}: exceeded {max_polls} polls "
                      f"— marking as failed")
                job["status"] = "FAILED"
                failed_files.append(job["fasta_file"])
                continue

            try:
                status_response = requests.get(BLAST_URL, params={
                    "CMD":           "Get",
                    "RID":           rid,
                    "FORMAT_TYPE":   "XML",
                    "FORMAT_OBJECT": "SearchInfo",
                }, timeout=60)
                status_text = status_response.text

                if "Status=WAITING" in status_text:
                    elapsed = int(time.time() - job["submitted"])
                    print(f"  RID {rid}: waiting "
                          f"({elapsed}s elapsed, "
                          f"poll {poll_counts[rid]}/{max_polls})")
                    continue

                elif "Status=FAILED" in status_text:
                    print(f"  RID {rid}: FAILED on NCBI servers")
                    job["status"] = "FAILED"
                    failed_files.append(job["fasta_file"])
                    continue

                elif "Status=UNKNOWN" in status_text:
                    print(f"  RID {rid}: UNKNOWN — may have expired")
                    job["status"] = "FAILED"
                    failed_files.append(job["fasta_file"])
                    continue

                elif "Status=READY" in status_text:
                    print(f"  RID {rid}: READY — fetching results...")

                    result_response = requests.get(BLAST_URL, params={
                        "CMD":           "Get",
                        "RID":           rid,
                        "FORMAT_TYPE":   "XML",
                        "FORMAT_OBJECT": "Alignment",
                    }, timeout=300)
                    xml_text = result_response.text

                    # Parse XML — may raise on malformed response.
                    # Catch XML errors specifically and retry next
                    # poll round rather than marking as failed.
                    try:
                        blast_results = _parse_blast_xml(
                            xml_text            = xml_text,
                            evalue_threshold    = evalue,
                            min_pident          = 70.0,
                            min_align_length    = 44,
                            per_gene_thresholds = per_gene_thresholds,
                            query_gene_map      = job["query_gene_map"],
                        )
                    except Exception as xml_err:
                        print(f"  RID {rid}: XML parse error "
                              f"({xml_err}) — will retry")
                        continue   # retry next polling round

                    new_rows = _score_job_results(
                        records       = job["records"],
                        blast_results = blast_results,
                        metadata      = metadata,
                        taxonomy      = taxonomy,
                        glossary      = glossary,
                        label         = label,
                        rid           = rid,
                        fasta_file    = Path(job["fasta_file"]).name,
                    )

                    added = _append_to_csv(
                        csv_path      = results_csv,
                        new_rows      = new_rows,
                        existing_keys = existing_keys,
                        existing_rows = existing_rows,
                    )
                    n_rows_added  += added
                    job["status"]  = "COMPLETE"
                    completed_rids.append(rid)

                    elapsed = int(time.time() - job["submitted"])
                    print(f"  RID {rid}: COMPLETE — "
                          f"{added} rows saved ({elapsed}s total)")

            except requests.exceptions.Timeout:
                print(f"  RID {rid}: request timed out — will retry")
                continue

            except Exception as e:
                print(f"  RID {rid}: poll error ({e}) — will retry")
                continue

        pending = sum(1 for j in submitted_jobs
                      if j["status"] == "WAITING")
        if pending > 0:
            print(f"\n  {pending} job(s) still pending — "
                  f"waiting {POLL_INTERVAL}s...")
            time.sleep(POLL_INTERVAL)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_complete = sum(1 for j in submitted_jobs
                     if j["status"] == "COMPLETE")
    n_failed   = sum(1 for j in submitted_jobs
                     if j["status"] == "FAILED")

    print(f"\n{'='*60}")
    print(f"BLAST DIRECTORY COMPLETE")
    print(f"{'='*60}")
    print(f"  Files submitted : {len(submitted_jobs)}")
    print(f"  Files completed : {n_complete}")
    print(f"  Files failed    : {n_failed}")
    print(f"  Total rows added: {n_rows_added}")
    if failed_files:
        print(f"\n  Failed files (can be resubmitted):")
        for f in failed_files:
            print(f"    {Path(f).name}")

    return {
        "n_files_submitted": len(submitted_jobs),
        "n_files_completed": n_complete,
        "n_files_failed":    n_failed,
        "n_rows_added":      n_rows_added,
        "results_csv":       results_csv,
        "completed_rids":    completed_rids,
        "failed_files":      failed_files,
    }


# ─────────────────────────────────────────────────────────────────────────────
# SCORE FROM CSV
# ─────────────────────────────────────────────────────────────────────────────


def score_from_csv(
        results_csv: str,
        calibration_json: str,
        output_path: str = None,
        label_filter: str = "QUERY",
        label_map: dict = None,
) -> dict:
    """
    Apply the calibrated Bayesian chimera detection model to BLAST
    hit results stored in a CSV file.

    Parameters
    ----------
    results_csv      : str   path to BLAST hit results CSV
    calibration_json : str   path to model_calibration.json
    output_path      : str   optional path to save chimera calls CSV
    label_filter     : str   only score rows with this label value.
                             Default "QUERY".
    label_map        : dict  optional mapping of raw label -> label_filter
    """
    with open(calibration_json, "r") as f:
        cal = json.load(f)

    best              = cal["best_params"]
    gene_thresholds   = cal.get("gene_prop_neg_thresholds", {})
    orphan_taxa       = cal.get("known_orphan_taxa", [])
    orphan_thresh     = cal.get("orphan_prop_neg_threshold", 0.8)
    max_poorly        = cal.get("max_poorly_scoring_for_rescue", 3)

    prior             = best["prior_prob"]
    pident_thresh     = best["pident_threshold"]
    k                 = best["k"]
    beta              = best["beta"]
    post_var_thresh   = best["posterior_var_threshold"]
    llr_var_thresh    = best["taxon_llr_var_threshold"]
    prop_neg_global   = best["prop_neg_threshold"]
    no_hits_thresh    = best.get("no_hits_threshold", 1.1)
    single_hit_thresh = cal.get("single_hit_taxon_llr_threshold", -2.0)

    print(f"\n{'='*60}")
    print(f"SCORE FROM CSV")
    print(f"{'='*60}")
    print(f"  Results CSV          : {results_csv}")
    print(f"  Calibration JSON     : {calibration_json}")
    print(f"  Label filter         : {label_filter}")
    print(f"  prior_prob           = {prior}")
    print(f"  pident_threshold     = {pident_thresh}")
    print(f"  k                    = {k}")
    print(f"  beta                 = {beta}")
    print(f"  posterior_var        > {post_var_thresh}")
    print(f"  taxon_llr_var        > {llr_var_thresh}")
    print(f"  prop_neg (global)    > {prop_neg_global}")
    print(f"  no_hits              > {no_hits_thresh}")
    print(f"  max_poorly_scoring   = {max_poorly}")
    print(f"  Orphan taxa          : {orphan_taxa}")

    # ── Load hit data ─────────────────────────────────────────────────────────
    subsample_hits = defaultdict(list)
    no_hit_samples = set()
    parent_taxon   = {}

    with open(results_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw_label = row.get("label", "").strip()
            if label_map is not None:
                matched = label_map.get(raw_label) == label_filter
            else:
                matched = raw_label == label_filter
            if not matched:
                continue

            acc   = row.get("accession", "").strip()
            gene  = row.get("gene", "").strip()
            taxon = row.get("taxon", "").strip()
            if not acc or not gene:
                continue

            parent_taxon[(acc, gene)] = taxon

            if row.get("hit_accession") == "NO_HITS":
                try:
                    no_hit_samples.add((acc, gene, row["sample"]))
                except KeyError:
                    pass
                continue

            try:
                sample      = row["sample"]
                bitscore    = float(row["bitscore"])
                pident      = float(row["pident"])
                species_llr = float(row["species_llr"])
                gene_llr    = float(row["gene_llr"])
            except (ValueError, KeyError, TypeError):
                continue

            subsample_hits[(acc, gene, sample)].append({
                "bitscore":    bitscore,
                "pident":      pident,
                "species_llr": species_llr,
                "gene_llr":    gene_llr,
            })

    if not parent_taxon:
        raise ValueError(
            f"No rows found with label='{label_filter}' in {results_csv}"
        )

    print(f"\n  Parent sequences found : {len(parent_taxon)}")
    print(f"  Subsample groups       : {len(subsample_hits)}")

    # ── Group subsamples by parent ────────────────────────────────────────────
    parent_samples = defaultdict(list)
    for (acc, gene, sample), hits in subsample_hits.items():
        parent_samples[(acc, gene)].append({
            "sample": sample,
            "hits":   hits,
        })

    seen = {(acc, gene, sd["sample"])
            for (acc, gene), sds in parent_samples.items()
            for sd in sds}
    for (acc, gene, sample) in no_hit_samples:
        if (acc, gene, sample) not in seen:
            parent_samples[(acc, gene)].append({
                "sample": sample,
                "hits":   [],
            })
    for (acc, gene, sample) in no_hit_samples:
        if (acc, gene) not in parent_samples:
            parent_samples[(acc, gene)].append({
                "sample": sample,
                "hits":   [],
            })

    # ── Compute per-parent metrics ────────────────────────────────────────────
    print(f"\nComputing per-parent metrics...")
    metrics  = {}
    prior_lo = math.log(prior / (1 - prior))
    orphan_set = set(orphan_taxa)

    for (acc, gene), samples in parent_samples.items():
        taxon         = parent_taxon.get((acc, gene), "")
        genus         = taxon.split("_")[0] if taxon else ""
        is_orphan     = genus in orphan_set
        prelim_thresh = 0.3 if is_orphan else 0.5

        all_taxon_llrs       = []
        sample_posts         = []
        any_single_flag      = False
        cumulative_log_odds  = prior_lo
        n_subsamples_total   = len(samples)
        n_subsamples_no_hits = 0

        for sample_data in samples:
            hits = sample_data["hits"]
            if not hits:
                n_subsamples_no_hits += 1
                continue

            bitscores    = [h["bitscore"] for h in hits]
            weights      = softmax_weights(bitscores, beta=beta)
            sample_llr   = 0.0
            s_log_odds   = prior_lo
            top_llr      = None
            top_bitscore = -1

            for h, w in zip(hits, weights):
                llr_p     = pident_to_llr(h["pident"],
                                          threshold=pident_thresh,
                                          k=k)
                llr_taxon = h["species_llr"]
                llr_gene  = h["gene_llr"]
                combined  = w * (llr_taxon + llr_gene + llr_p)
                sample_llr         += combined
                s_log_odds         += llr_taxon + llr_gene + llr_p
                all_taxon_llrs.append(llr_taxon)
                if h["bitscore"] > top_bitscore:
                    top_bitscore = h["bitscore"]
                    top_llr      = llr_taxon

            cumulative_log_odds += sample_llr
            sample_posts.append(sigmoid(s_log_odds))

            if top_llr is not None and top_llr < single_hit_thresh:
                any_single_flag = True

        post_var     = (float(np.var(sample_posts))
                        if len(sample_posts) > 1 else 0.0)
        llr_var      = (float(np.var(all_taxon_llrs))
                        if len(all_taxon_llrs) > 1 else 0.0)
        prop_neg     = (sum(1 for x in all_taxon_llrs if x < 0) /
                        len(all_taxon_llrs) if all_taxon_llrs else 0.0)
        prop_no_hits = (n_subsamples_no_hits / n_subsamples_total
                        if n_subsamples_total > 0 else 0.0)

        metrics[(acc, gene)] = {
            "posterior_variance":   post_var,
            "taxon_llr_variance":   llr_var,
            "prop_neg_taxon_hits":  prop_neg,
            "prop_no_hits":         prop_no_hits,
            "n_subsamples_no_hits": n_subsamples_no_hits,
            "n_subsamples_total":   n_subsamples_total,
            "mean_posterior":       (float(np.mean(sample_posts))
                                     if sample_posts else 0.0),
            "single_hit_flag":      any_single_flag,
            "cumulative_log_odds":  cumulative_log_odds,
            "label":                label_filter,
            "taxon":                taxon,
            "chimera_flag":         (any_single_flag
                                     or prop_neg > prelim_thresh
                                     or prop_no_hits > 0.5),
        }

    # ── Apply cross-gene rescue ───────────────────────────────────────────────
    print(f"Applying cross-gene consistency rescue...")
    metrics_rescued = apply_cross_gene_consistency(
        metrics,
        max_poorly_scoring_for_rescue = max_poorly,
        known_orphan_taxa             = orphan_taxa or None,
        orphan_prop_neg_threshold     = orphan_thresh,
    )

    # ── Make chimera calls ────────────────────────────────────────────────────
    print(f"Making chimera calls...")
    chimera_calls = {}
    flagged       = []

    for (acc, gene), m in metrics_rescued.items():
        prop_neg_thresh = gene_thresholds.get(gene, prop_neg_global)

        if m.get("cross_gene_rescue", False):
            chimera = False
            trigger = None
        else:
            if m["posterior_variance"] > post_var_thresh:
                chimera, trigger = True, "posterior_variance"
            elif m["taxon_llr_variance"] > llr_var_thresh:
                chimera, trigger = True, "taxon_llr_variance"
            elif m["prop_neg_taxon_hits"] > prop_neg_thresh:
                chimera, trigger = True, "prop_neg_taxon_hits"
            elif m["single_hit_flag"]:
                chimera, trigger = True, "single_hit"
            elif m["prop_no_hits"] > no_hits_thresh:
                chimera, trigger = True, "prop_no_hits"
            else:
                chimera, trigger = False, None

        chimera_calls[(acc, gene)] = {
            "accession":            acc,
            "gene":                 gene,
            "taxon":                m.get("taxon", ""),
            "chimera":              chimera,
            "trigger":              trigger,
            "prop_neg":             round(m["prop_neg_taxon_hits"], 4),
            "prop_no_hits":         round(m["prop_no_hits"], 4),
            "posterior_variance":   round(m["posterior_variance"], 6),
            "taxon_llr_variance":   round(m["taxon_llr_variance"], 4),
            "mean_posterior":       round(m["mean_posterior"], 4),
            "cross_gene_rescue":    m.get("cross_gene_rescue", False),
            "well_scoring_frac":    round(
                m.get("well_scoring_fraction", 0.0), 4),
            "n_subsamples_no_hits": m.get("n_subsamples_no_hits", 0),
            "n_subsamples_total":   m.get("n_subsamples_total", 0),
        }

        if chimera:
            flagged.append((acc, gene))

    flagged.sort()
    n_chimeras = len(flagged)
    n_genuine  = len(chimera_calls) - n_chimeras

    # ── Print summary ─────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"VERIFICATION COMPLETE")
    print(f"{'='*60}")
    print(f"  Sequences scored  : {len(chimera_calls)}")
    print(f"  Flagged chimeric  : {n_chimeras}")
    print(f"  Cleared genuine   : {n_genuine}")

    if flagged:
        print(f"\n  Chimeric sequences:")
        print(f"  {'Accession':<16} {'Gene':<10} {'Taxon':<28} "
              f"{'Trigger':<22} {'prop_neg':>9} {'prop_no_hits':>13}")
        print(f"  {'-'*100}")
        for (acc, gene) in flagged:
            c = chimera_calls[(acc, gene)]
            print(f"  {acc:<16} {gene:<10} {c['taxon']:<28} "
                  f"{str(c['trigger']):<22} "
                  f"{c['prop_neg']:>9.3f} "
                  f"{c['prop_no_hits']:>13.3f}")

    rescued = [(acc, gene) for (acc, gene), m
               in metrics_rescued.items()
               if m.get("cross_gene_rescue", False)]
    if rescued:
        print(f"\n  Sequences rescued by cross-gene logic: {len(rescued)}")
        for acc, gene in sorted(rescued):
            print(f"    {acc} | {gene}")

    # ── Save output CSV ───────────────────────────────────────────────────────
    summary_csv = None
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fields = [
            "accession", "gene", "taxon", "chimera", "trigger",
            "prop_neg", "prop_no_hits", "posterior_variance",
            "taxon_llr_variance", "mean_posterior",
            "cross_gene_rescue", "well_scoring_frac",
            "n_subsamples_no_hits", "n_subsamples_total",
        ]
        with open(output_path, "w", newline="",
                  encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields,
                                    extrasaction="ignore")
            writer.writeheader()
            for row in sorted(chimera_calls.values(),
                               key=lambda x: (x["accession"],
                                              x["gene"])):
                writer.writerow(row)
        summary_csv = output_path
        print(f"\n  Results saved: {output_path}")

    return {
        "chimera_calls": chimera_calls,
        "flagged":       flagged,
        "n_query_seqs":  len(chimera_calls),
        "n_chimeras":    n_chimeras,
        "n_genuine":     n_genuine,
        "summary_csv":   summary_csv,
        "orphan_taxa":   list(orphan_set)
    }

from bioinf_packages.verify_funcs._gene_parser import load_glossary
from bioinf_packages.verify_funcs._species_parser import load_taxonomy

#blast_subsample_directory(
#    subsample_dir          = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/query",
#    metadata_csv           = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/query/subsample_batch_metadata.csv",
#    taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
#    glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"),
#    results_csv            = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/blast_results.csv",
#    label                  = "QUERY",
#    filter_thresholds_json = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/filter_thresholds.json",
#)


results = score_from_csv(
    results_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/blast_results.csv",
    calibration_json = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/model_calibration.json",
    output_path = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/chimera_calls.csv",
    label_filter = "GENUINE",
)

print(f"\n{results['n_chimeras']} chimeric sequences detected "
      f"out of {results['n_query_seqs']} total")

for acc, gene in results["flagged"]:
    c = results["chimera_calls"][(acc, gene)]
    print(f"  CHIMERA: {acc} | {gene} | "
          f"trigger={c['trigger']} | "
          f"prop_neg={c['prop_neg']:.3f}")
