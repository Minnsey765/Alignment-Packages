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

from bioinf_packages.verify_funcs._compute_csv import apply_cross_gene_consistency
from bioinf_packages.verify_funcs._score_hit import (pident_to_llr, sigmoid, softmax_weights, score_hit, aggregate_hit_scores)

from bioinf_packages.verify_funcs._subsamples_to_blast import (
    _submit_blast,
    _poll_blast,
    _parse_blast_xml,
    load_metadata,
    parse_subsample_header,
    _load_existing_csv,
    _append_to_csv,
    CSV_FIELDS,
)

from Bio import SeqIO
from pathlib import Path
import numpy as np
# ─────────────────────────────────────────────────────────────────────────────
# NCBI USAGE POLICY HELPERS
# ─────────────────────────────────────────────────────────────────────────────

EASTERN = pytz.timezone("US/Eastern")
POLL_INTERVAL = 30   # seconds between submissions


def _is_ncbi_peak_hours() -> bool:
    """
    Return True if current time is within NCBI peak hours.

    NCBI requests that batch BLAST jobs are submitted outside of
    peak hours: weekdays 5am-9pm Eastern Standard Time.
    Weekend submissions are always permitted.

    Returns
    -------
    bool   True if currently in peak hours (should wait)
    """
    now_eastern = datetime.datetime.now(EASTERN)
    weekday     = now_eastern.weekday()   # 0=Monday, 6=Sunday
    hour        = now_eastern.hour        # 0-23

    if weekday >= 5:
        # Saturday or Sunday — always off-peak
        return False

    # Weekday — peak hours are 5am to 9pm (hour 5 to 20 inclusive)
    return 5 <= hour < 21


def _wait_for_off_peak(check_interval: int = 300) -> None:
    """
    Block until NCBI off-peak hours if currently in peak hours.

    Prints a message every check_interval seconds while waiting.

    Parameters
    ----------
    check_interval : int   seconds between peak hour checks while
                           waiting. Default 300 (5 minutes).
    """
    if not _is_ncbi_peak_hours():
        return

    now_eastern  = datetime.datetime.now(EASTERN)
    # Off-peak starts at 9pm (21:00) or next day at midnight
    # if already past 9pm (should not happen if in peak hours)
    if now_eastern.hour < 21:
        off_peak_today = now_eastern.replace(
            hour=21, minute=0, second=0, microsecond=0
        )
        wait_seconds = (off_peak_today - now_eastern).total_seconds()
    else:
        # Should not reach here given is_peak logic, but defensive
        return

    print(f"\n  Current time (Eastern): "
          f"{now_eastern.strftime('%A %H:%M')}")
    print(f"  NCBI peak hours are weekdays 5am-9pm Eastern.")
    print(f"  Off-peak starts at: "
          f"{off_peak_today.strftime('%H:%M')} Eastern "
          f"({int(wait_seconds/3600)}h "
          f"{int((wait_seconds%3600)/60)}m away)")
    print(f"  Waiting for off-peak hours...")

    while _is_ncbi_peak_hours():
        now_eastern = datetime.datetime.now(EASTERN)
        remaining   = (off_peak_today - now_eastern).total_seconds()
        print(f"    [{now_eastern.strftime('%H:%M')}] "
              f"Still in peak hours — "
              f"{int(remaining/3600)}h "
              f"{int((remaining%3600)/60)}m remaining...")
        time.sleep(check_interval)

    print(f"  Off-peak hours reached — proceeding with submission.")


# ─────────────────────────────────────────────────────────────────────────────
# BATCH BLAST WITH PARALLEL SUBMISSION
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

    Submission strategy
    -------------------
    All files are submitted first (one every submission_interval
    seconds), then results are polled concurrently and saved to CSV
    as each job completes. This is faster than sequential
    submit-wait-save for multiple files since NCBI processes jobs
    in parallel.

    NCBI usage policy
    -----------------
    If respect_peak_hours=True (default), submission is delayed until
    NCBI off-peak hours: weekday evenings after 9pm Eastern or
    weekends. This is in line with NCBI's guidelines for batch
    submissions. The function checks peak hours before each
    submission batch and will wait if necessary.

    A minimum of submission_interval seconds (default 30s) is
    enforced between successive submissions regardless of peak hours,
    in line with NCBI's rate limiting guidelines.

    Parameters
    ----------
    subsample_dir       : str   directory containing subsample FASTAs
    metadata_csv        : str   metadata CSV from prepare_subsample_fastas
    taxonomy            : dict  from load_taxonomy()
    glossary            : dict  from load_glossary()
    results_csv         : str   accumulating results CSV
    label               : str   label for CSV rows. Default "QUERY".
    file_prefix         : str   FASTA filename prefix to glob
    filter_thresholds_json: str optional per-gene filter thresholds
    evalue              : float BLAST e-value threshold
    max_hits            : int   max hits per query sequence
    database            : str   BLAST database
    entrez_query        : str   Entrez filter
    max_polls           : int   max poll attempts per job (2 min each)
    respect_peak_hours  : bool  wait for NCBI off-peak if True.
                                Default True.
    submission_interval : int   minimum seconds between submissions.
                                Default 30.

    Returns
    -------
    dict with:
        n_files_submitted  : int
        n_files_completed  : int
        n_files_failed     : int
        n_rows_added       : int
        results_csv        : str
        completed_rids     : list
        failed_files       : list
    """


    pattern     = os.path.join(
        subsample_dir, f"{file_prefix}_*.fasta"
    )
    fasta_files = sorted(glob.glob(pattern))

    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files matching '{pattern}'"
        )

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
    if filter_thresholds_json and os.path.exists(
            filter_thresholds_json):
        with open(filter_thresholds_json) as f:
            thresh_data         = json.load(f)
            per_gene_thresholds = thresh_data.get(
                "per_gene_thresholds"
            )
        print(f"  Filter thresholds loaded for "
              f"{len(per_gene_thresholds)} genes")

    # Load metadata and existing CSV once
    metadata                   = load_metadata(metadata_csv)
    existing_keys, existing_rows = _load_existing_csv(results_csv)
    print(f"  Existing CSV rows: {len(existing_rows)}")

    # Check / wait for off-peak before starting
    if respect_peak_hours:
        _wait_for_off_peak()

    # ── Phase 1: Submit all jobs ──────────────────────────────────────────────
    print(f"\n── Phase 1: Submitting {len(fasta_files)} jobs ──────────")

    submitted_jobs = []   # list of dicts with rid, fasta_file, records
    last_submit_time = 0.0

    for i, fasta_file in enumerate(fasta_files, 1):

        # Check peak hours before each submission
        if respect_peak_hours and _is_ncbi_peak_hours():
            print(f"\n  Peak hours detected before submission "
                  f"{i}/{len(fasta_files)} — waiting...")
            _wait_for_off_peak()

        # Enforce minimum interval between submissions
        elapsed = time.time() - last_submit_time
        if elapsed < submission_interval:
            wait = submission_interval - elapsed
            print(f"  Waiting {wait:.0f}s before next submission...")
            time.sleep(wait)

        # Load and validate FASTA
        records = list(SeqIO.parse(fasta_file, "fasta"))
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

        # Read FASTA content
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
            })
            print(f"  [{i}/{len(fasta_files)}] SUBMITTED "
                  f"{Path(fasta_file).name} — "
                  f"RID={rid} ({total_bp:,}bp, "
                  f"{len(records)} seqs)")

        except Exception as e:
            print(f"  [{i}/{len(fasta_files)}] FAILED to submit "
                  f"{Path(fasta_file).name}: {e}")

    print(f"\n  Submitted {len(submitted_jobs)} jobs successfully")

    # ── Phase 2: Poll and save as results arrive ──────────────────────────────
    print(f"\n── Phase 2: Polling for results ─────────────────────────")
    print(f"  Polling every {POLL_INTERVAL}s per job, "
          f"max {max_polls} attempts each")

    # Build query_gene_map for each job
    for job in submitted_jobs:
        job["query_gene_map"] = {
            r.id: parse_subsample_header(r.id).get("gene")
            for r in job["records"]
        }

    completed_rids = []
    failed_files   = []
    n_rows_added   = 0
    poll_counts    = {job["rid"]: 0 for job in submitted_jobs}

    # Keep polling until all jobs are done or max_polls exceeded
    while any(j["status"] == "WAITING" for j in submitted_jobs):

        for job in submitted_jobs:
            if job["status"] != "WAITING":
                continue

            rid        = job["rid"]
            poll_counts[rid] += 1

            if poll_counts[rid] > max_polls:
                print(f"  RID {rid}: exceeded {max_polls} polls "
                      f"— marking as failed")
                job["status"] = "FAILED"
                failed_files.append(job["fasta_file"])
                continue

            # Poll status
            try:
                status_response = requests.get(BLAST_URL, params={
                    "CMD":          "Get",
                    "RID":          rid,
                    "FORMAT_TYPE":  "XML",
                    "FORMAT_OBJECT":"SearchInfo",
                })
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

                    # Fetch XML
                    result_response = requests.get(BLAST_URL, params={
                        "CMD":           "Get",
                        "RID":           rid,
                        "FORMAT_TYPE":   "XML",
                        "FORMAT_OBJECT": "Alignment",
                    })
                    xml_text = result_response.text

                    # Parse with per-gene thresholds
                    blast_results = _parse_blast_xml(
                        xml_text            = xml_text,
                        evalue_threshold    = evalue,
                        min_pident          = 70.0,
                        min_align_length    = 44,
                        per_gene_thresholds = per_gene_thresholds,
                        query_gene_map      = job["query_gene_map"],
                    )

                    # Score and write to CSV
                    new_rows = _score_job_results(
                        records      = job["records"],
                        blast_results= blast_results,
                        metadata     = metadata,
                        taxonomy     = taxonomy,
                        glossary     = glossary,
                        label        = label,
                        rid          = rid,
                        fasta_file   = Path(job["fasta_file"]).name,
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
                          f"{added} rows saved "
                          f"({elapsed}s total)")

            except Exception as e:
                print(f"  RID {rid}: poll error — {e}")
                continue

        # Wait before next polling round if jobs still pending
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


def score_from_csv(
        results_csv: str,
        calibration_json: str,
        output_path: str = None,
        label_filter: str = "QUERY",
) -> dict:
    """
    Apply the calibrated Bayesian chimera detection model to BLAST
    hit results stored in a CSV file and return chimera calls for
    each (accession, gene) pair.

    Reads the calibration JSON to get all model parameters including
    gene-specific prop_neg thresholds, orphan taxa, and Bayesian
    scoring parameters. Applies cross-gene consistency rescue after
    computing per-gene metrics. Produces a final chimera call per
    (accession, gene) pair.

    This is the second step in live verification, called after
    blast_subsample_directory() has populated the results CSV.

    Parameters
    ----------
    results_csv      : str   path to BLAST hit results CSV from
                             blast_subsample_directory() or
                             subsamples_to_blast()
    calibration_json : str   path to model_calibration.json from
                             calibrate_from_csv()
    output_path      : str   optional path to save chimera calls as
                             CSV. If None, results are returned only.
    label_filter     : str   only score rows with this label value.
                             Use "QUERY" for live verification,
                             "GENUINE"/"ARTEFACT" for re-scoring
                             calibration data. Default "QUERY".

    Returns
    -------
    dict with:
        chimera_calls   : dict   (accession, gene) -> call dict with:
                                     chimera       : bool
                                     trigger       : str or None
                                     prop_neg      : float
                                     post_var      : float
                                     llr_var       : float
                                     mean_posterior: float
                                     cross_gene_rescue: bool
                                     taxon         : str
        flagged         : list   (accession, gene) pairs flagged as
                                 chimeric, sorted by accession then gene
        n_query_seqs    : int    total parent sequences scored
        n_chimeras      : int    number flagged as chimeric
        n_genuine       : int    number cleared as genuine
        summary_csv     : str    path to output CSV if written
    """
    # ── Load calibration parameters ───────────────────────────────────────────
    with open(calibration_json, "r") as f:
        cal = json.load(f)

    best             = cal["best_params"]
    gene_thresholds  = cal.get("gene_prop_neg_thresholds", {})
    orphan_taxa      = cal.get("known_orphan_taxa", [])
    orphan_thresh    = cal.get("orphan_prop_neg_threshold", 0.8)

    prior            = best["prior_prob"]
    pident_thresh    = best["pident_threshold"]
    k                = best["k"]
    beta             = best["beta"]
    post_var_thresh  = best["posterior_var_threshold"]
    llr_var_thresh   = best["taxon_llr_var_threshold"]
    prop_neg_global  = best["prop_neg_threshold"]

    print(f"\n{'='*60}")
    print(f"SCORE FROM CSV")
    print(f"{'='*60}")
    print(f"  Results CSV      : {results_csv}")
    print(f"  Calibration JSON : {calibration_json}")
    print(f"  Label filter     : {label_filter}")
    print(f"\n  Bayesian params:")
    print(f"    prior_prob             = {prior}")
    print(f"    pident_threshold       = {pident_thresh}")
    print(f"    k                      = {k}")
    print(f"    beta                   = {beta}")
    print(f"  Chimera thresholds:")
    print(f"    posterior_var          > {post_var_thresh}")
    print(f"    taxon_llr_var          > {llr_var_thresh}")
    print(f"    prop_neg (global)      > {prop_neg_global}")
    print(f"  Orphan taxa: {orphan_taxa}")

    # ── Load and group hit data ───────────────────────────────────────────────
    subsample_hits = defaultdict(list)
    parent_taxon   = {}

    with open(results_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = row.get("label", "").strip()
            if label_filter and label != label_filter:
                continue

            acc  = row.get("accession", "").strip()
            gene = row.get("gene", "").strip()
            taxon = row.get("taxon", "").strip()
            if not acc or not gene:
                continue

            parent_taxon[(acc, gene)] = taxon

            if row.get("hit_accession") == "NO_HITS":
                continue

            try:
                sample      = row["sample"]
                bitscore    = float(row["bitscore"])
                pident      = float(row["pident"])
                species_llr = float(row["species_llr"])
                gene_llr    = float(row["gene_llr"])
            except (ValueError, KeyError, TypeError):
                continue

            hit_key = (acc, gene, sample)
            subsample_hits[hit_key].append({
                "bitscore":    bitscore,
                "pident":      pident,
                "species_llr": species_llr,
                "gene_llr":    gene_llr,
            })

    if not parent_taxon:
        raise ValueError(
            f"No rows found with label='{label_filter}' in "
            f"{results_csv}"
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

    # ── Compute per-parent metrics ────────────────────────────────────────────
    print(f"\nComputing per-parent metrics...")
    metrics  = {}
    prior_lo = math.log(prior / (1 - prior))

    orphan_set      = set(orphan_taxa)
    prelim_orphan_t = 0.3   # lower preliminary threshold for orphan taxa

    for (acc, gene), samples in parent_samples.items():
        taxon          = parent_taxon.get((acc, gene), "")
        genus          = taxon.split("_")[0] if taxon else ""
        is_orphan      = genus in orphan_set
        prelim_thresh  = prelim_orphan_t if is_orphan else 0.5

        all_taxon_llrs      = []
        sample_posts        = []
        any_single_flag     = False
        cumulative_log_odds = prior_lo

        for sample_data in samples:
            hits = sample_data["hits"]
            if not hits:
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

            if (top_llr is not None
                    and top_llr < cal.get(
                        "single_hit_taxon_llr_threshold", -2.0)):
                any_single_flag = True

        post_var = (float(np.var(sample_posts))
                    if len(sample_posts) > 1 else 0.0)
        llr_var  = (float(np.var(all_taxon_llrs))
                    if len(all_taxon_llrs) > 1 else 0.0)
        prop_neg = (sum(1 for x in all_taxon_llrs if x < 0) /
                    len(all_taxon_llrs) if all_taxon_llrs else 0.0)

        metrics[(acc, gene)] = {
            "posterior_variance":  post_var,
            "taxon_llr_variance":  llr_var,
            "prop_neg_taxon_hits": prop_neg,
            "mean_posterior":      (float(np.mean(sample_posts))
                                    if sample_posts else 0.0),
            "single_hit_flag":     any_single_flag,
            "cumulative_log_odds": cumulative_log_odds,
            "label":               label_filter,
            "taxon":               taxon,
            "chimera_flag":        any_single_flag or
                                   prop_neg > prelim_thresh,
        }

    # ── Apply cross-gene rescue ───────────────────────────────────────────────
    print(f"Applying cross-gene consistency rescue...")
    metrics_rescued = apply_cross_gene_consistency(
        metrics,
        known_orphan_taxa         = orphan_taxa or None,
        orphan_prop_neg_threshold = orphan_thresh,
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
                chimera = True
                trigger = "posterior_variance"
            elif m["taxon_llr_variance"] > llr_var_thresh:
                chimera = True
                trigger = "taxon_llr_variance"
            elif m["prop_neg_taxon_hits"] > prop_neg_thresh:
                chimera = True
                trigger = "prop_neg_taxon_hits"
            elif m["single_hit_flag"]:
                chimera = True
                trigger = "single_hit"
            else:
                chimera = False
                trigger = None

        chimera_calls[(acc, gene)] = {
            "accession":         acc,
            "gene":              gene,
            "taxon":             m.get("taxon", ""),
            "chimera":           chimera,
            "trigger":           trigger,
            "prop_neg":          round(m["prop_neg_taxon_hits"], 4),
            "posterior_variance":round(m["posterior_variance"], 6),
            "taxon_llr_variance":round(m["taxon_llr_variance"], 4),
            "mean_posterior":    round(m["mean_posterior"], 4),
            "cross_gene_rescue": m.get("cross_gene_rescue", False),
            "well_scoring_frac": round(
                m.get("well_scoring_fraction", 0.0), 4),
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
              f"{'Trigger':<22} {'prop_neg':>9}")
        print(f"  {'-'*87}")
        for (acc, gene) in flagged:
            c = chimera_calls[(acc, gene)]
            print(f"  {acc:<16} {gene:<10} {c['taxon']:<28} "
                  f"{str(c['trigger']):<22} {c['prop_neg']:>9.3f}")

    rescued = [(acc, gene) for (acc, gene), m
               in metrics_rescued.items()
               if m.get("cross_gene_rescue", False)]
    if rescued:
        print(f"\n  Sequences rescued by cross-gene logic: "
              f"{len(rescued)}")
        for acc, gene in sorted(rescued):
            print(f"    {acc} | {gene}")

    # ── Save output CSV ───────────────────────────────────────────────────────
    summary_csv = None
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fields = [
            "accession", "gene", "taxon", "chimera", "trigger",
            "prop_neg", "posterior_variance", "taxon_llr_variance",
            "mean_posterior", "cross_gene_rescue", "well_scoring_frac",
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
    }

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
    Extracted from the main loop to keep blast_subsample_directory
    readable.
    """
    from bioinf_packages.verify_funcs._subsamples_to_blast import (
        parse_subsample_header,
    )
    from bioinf_packages.verify_funcs._score_hit import (
        score_hit,
        aggregate_hit_scores,
    )
    import numpy as np
    from collections import defaultdict

    new_rows              = []
    parent_subsample_data = defaultdict(list)

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
            "n_samples":        0,
        })

        hits_raw = blast_results.get(record.id, [])
        if not hits_raw:
            for blast_id, hits in blast_results.items():
                if (record.id.startswith(blast_id)
                        or blast_id.startswith(record.id[:30])):
                    hits_raw = hits
                    break

        if not hits_raw:
            new_rows.append({
                "label":               label,
                "taxon":               taxon,
                "gene":                gene,
                "accession":           accession,
                "orientation":         parent.get("orientation", "+"),
                "parent_length":       parent.get("seq_length", 0),
                "subsample_length":    parent.get("subsample_length",
                                                  0),
                "n_samples_total":     parent.get("n_samples", 0),
                "sample":              sample,
                "subsample_start":     sub_info.get("subsample_start"),
                "subsample_end":       sub_info.get("subsample_end"),
                "hit_accession":       "NO_HITS",
                "hit_description":     "",
                "bitscore":            "",
                "pident":              "",
                "evalue":              "",
                "align_length":        "",
                "species_llr":         "",
                "gene_llr":            "",
                "total_hit_llr":       "",
                "subsample_posterior": "",
                "subsample_total_llr": "",
                "n_hits_this_subsample": 0,
                "rid":                 rid,
                "fasta_file":          fasta_file,
            })
            continue

        scored_hits   = [
            score_hit(hit, taxon, gene, taxonomy, glossary)
            for hit in hits_raw
        ]
        subsample_agg = aggregate_hit_scores(scored_hits)

        parent_subsample_data[meta_key].append({
            "sample":      sample,
            "posterior":   subsample_agg["posterior_prob"],
            "total_llr":   subsample_agg["total_llr"],
            "scored_hits": scored_hits,
        })

        for hit, scored in zip(hits_raw, scored_hits):
            new_rows.append({
                "label":               label,
                "taxon":               taxon,
                "gene":                gene,
                "accession":           accession,
                "orientation":         parent.get("orientation", "+"),
                "parent_length":       parent.get("seq_length", 0),
                "subsample_length":    parent.get("subsample_length",
                                                  0),
                "n_samples_total":     parent.get("n_samples", 0),
                "sample":              sample,
                "subsample_start":     sub_info.get("subsample_start"),
                "subsample_end":       sub_info.get("subsample_end"),
                "hit_accession":       scored["accession"],
                "hit_description":     hit["description"][:200],
                "bitscore":            round(scored["bitscore"], 2),
                "pident":              round(scored["pident"], 4),
                "evalue":              scored["evalue"],
                "align_length":        scored["align_length"],
                "species_llr":         round(scored["species_llr"], 4),
                "gene_llr":            round(scored["gene_llr"], 4),
                "total_hit_llr":       round(scored["total_llr"], 4),
                "subsample_posterior": round(
                    subsample_agg["posterior_prob"], 4),
                "subsample_total_llr": round(
                    subsample_agg["total_llr"], 4),
                "n_hits_this_subsample": len(hits_raw),
                "rid":                 rid,
                "fasta_file":          fasta_file,
            })

    return new_rows

from bioinf_packages.verify_funcs._gene_parser import load_glossary
from bioinf_packages.verify_funcs._species_parser import load_taxonomy

blast_subsample_directory(
    subsample_dir          = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/query",
    metadata_csv           = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/query/subsample_batch_metadata.csv",
    taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
    glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"),
    results_csv            = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/blast_results.csv",
    label                  = "QUERY",
    filter_thresholds_json = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/filter_thresholds.json",
)