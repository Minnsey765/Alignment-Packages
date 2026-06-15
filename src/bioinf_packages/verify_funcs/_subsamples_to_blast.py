# subsamples_to_blast.py

import os
import csv
import time
import requests
import json

from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

try:
    from bioinf_packages.verify_funcs._score_hit import (
        score_hit,
        sigmoid,
        pident_to_llr,
        softmax_weights,
        aggregate_hit_scores,
    )
    from bioinf_packages.verify_funcs._species_parser import load_taxonomy
    from bioinf_packages.verify_funcs._gene_parser import load_glossary
except ImportError:
    from _score_hit import (
        score_hit,
        sigmoid,
        pident_to_llr,
        softmax_weights,
        aggregate_hit_scores,
    )
    from _species_parser import load_taxonomy
    from _gene_parser import load_glossary

BLAST_URL   = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
POLL_INTERVAL = 120   # seconds between status checks — NCBI ask for
                      # no more than one request per minute; 2 minutes
                      # is courteous and avoids rate limiting


# ─────────────────────────────────────────────────────────────────────────────
# METADATA LOADER
# ─────────────────────────────────────────────────────────────────────────────

def load_metadata(metadata_csv: str) -> dict:
    """
    Load the metadata CSV written by prepare_subsample_fastas() into a
    dict keyed by (accession, gene) tuples.

    The same accession number can appear more than once if the parent
    FASTA contained sequences for the same accession from different
    genes (e.g. a mitogenome accession contributing both CYTB and ND2
    sequences). Keying on (accession, gene) handles this correctly.

    Parameters
    ----------
    metadata_csv : str   path to {file_prefix}_metadata.csv written by
                         prepare_subsample_fastas()

    Returns
    -------
    dict keyed by (accession, gene) -> metadata dict with:
        taxon, gene, accession, orientation, seq_length,
        subsample_length, n_samples, file_indices
    """
    metadata = {}
    with open(metadata_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row["accession"].strip(), row["gene"].strip())
            metadata[key] = {
                "taxon":            row["taxon"].strip(),
                "gene":             row["gene"].strip(),
                "accession":        row["accession"].strip(),
                "orientation":      row.get("orientation", "+").strip(),
                "seq_length":       int(row["seq_length"]),
                "subsample_length": int(row["subsample_length"]),
                "n_samples":        int(row["n_samples"]),
                "file_indices":     row.get("file_indices", ""),
            }
    print(f"Loaded metadata for {len(metadata)} parent sequences "
          f"from {metadata_csv}")
    return metadata


# ─────────────────────────────────────────────────────────────────────────────
# HEADER PARSER
# ─────────────────────────────────────────────────────────────────────────────

def parse_subsample_header(header: str) -> dict:
    """
    Parse a subsample sequence header written by prepare_subsample_fastas().

    Format:
        Genus_species|gene|orien:+/-|accession:XXXX|Sample:N|
        parent_length:L|subsample_start:S|subsample_end:E

    Parameters
    ----------
    header : str   sequence ID from the subsample FASTA

    Returns
    -------
    dict with taxon, gene, orientation, accession, sample,
             parent_length, subsample_start, subsample_end
    """
    parts       = header.strip().split("|")
    taxon       = parts[0].strip()
    gene        = parts[1].strip() if len(parts) > 1 else None
    orientation = "+"
    accession   = None
    sample      = None
    parent_len  = None
    sub_start   = None
    sub_end     = None

    for p in parts:
        p = p.strip()
        if p.startswith("orien:"):
            orientation = p.split(":", 1)[1].strip()
        elif p.startswith("accession:"):
            accession = p.split(":", 1)[1].strip()
        elif p.startswith("Sample:"):
            try:
                sample = int(p.split(":", 1)[1].strip())
            except ValueError:
                pass
        elif p.startswith("parent_length:"):
            try:
                parent_len = int(p.split(":", 1)[1].strip())
            except ValueError:
                pass
        elif p.startswith("subsample_start:"):
            try:
                sub_start = int(p.split(":", 1)[1].strip())
            except ValueError:
                pass
        elif p.startswith("subsample_end:"):
            try:
                sub_end = int(p.split(":", 1)[1].strip())
            except ValueError:
                pass

    return {
        "taxon":            taxon,
        "gene":             gene,
        "orientation":      orientation,
        "accession":        accession,
        "sample":           sample,
        "parent_length":    parent_len,
        "subsample_start":  sub_start,
        "subsample_end":    sub_end,
    }


# ─────────────────────────────────────────────────────────────────────────────
# REMOTE BLAST SUBMISSION AND POLLING
# ─────────────────────────────────────────────────────────────────────────────

def _submit_blast(fasta_content: str,
                  evalue: float = 1e-10,
                  max_hits: int = 20,
                  database: str = "nt",
                  entrez_query: str = "Vertebrata[Organism]") -> str:
    """
    Submit a multi-sequence FASTA to NCBI remote BLAST and return
    the RID (Request ID) for polling.

    Parameters
    ----------
    fasta_content  : str   full content of the FASTA file as a string
    evalue         : float e-value threshold
    max_hits       : int   maximum hits per query sequence
    database       : str   BLAST database (default "nt")
    entrez_query   : str   Entrez filter to restrict results.
                           Default "Vertebrata[Organism]" — restricts
                           to vertebrates, which covers all your taxa
                           and common mammalian contaminants while
                           excluding bacteria, plants, fungi.

    Returns
    -------
    str   RID (Request ID) for polling
    """
    params = {
        "CMD":          "Put",
        "PROGRAM":      "blastn",
        "DATABASE":     database,
        "QUERY":        fasta_content,
        "HITLIST_SIZE": str(max_hits),
        "EXPECT":       str(evalue),
        "FORMAT_TYPE":  "XML",
        "ENTREZ_QUERY": entrez_query,
        "tool":         "bioinf_packages_chimera_check",
        "email":        "om380@cam.ac.uk",
    }

    response = requests.post(BLAST_URL, data=params)
    response.raise_for_status()

    rid = None
    for line in response.text.splitlines():
        if "RID =" in line:
            rid = line.split("=")[1].strip()
            break

    if not rid:
        raise RuntimeError(
            f"Failed to extract RID from BLAST submission response. "
            f"Response text: {response.text[:500]}"
        )

    print(f"  BLAST job submitted. RID: {rid}")
    return rid


def _poll_blast(rid: str,
                max_polls: int = 60) -> str:
    """
    Poll NCBI BLAST for results, waiting POLL_INTERVAL seconds between
    each check. Returns the XML result string when finished.

    Parameters
    ----------
    rid        : str   Request ID from _submit_blast()
    max_polls  : int   maximum number of status checks before giving up.
                       At 2 minutes per poll this is 2 hours by default.

    Returns
    -------
    str   raw XML result string
    """
    print(f"  Polling for results (every {POLL_INTERVAL}s, "
          f"max {max_polls} attempts = "
          f"{max_polls * POLL_INTERVAL // 60} minutes)...")

    for poll_num in range(1, max_polls + 1):
        time.sleep(POLL_INTERVAL)

        status_response = requests.get(BLAST_URL, params={
            "CMD":         "Get",
            "RID":         rid,
            "FORMAT_TYPE": "XML",
            "FORMAT_OBJECT":"SearchInfo",
        })
        status_response.raise_for_status()

        status_text = status_response.text
        if "Status=WAITING" in status_text:
            print(f"  [{poll_num}/{max_polls}] Still running...")
            continue
        elif "Status=FAILED" in status_text:
            raise RuntimeError(
                f"BLAST job {rid} failed on NCBI servers."
            )
        elif "Status=UNKNOWN" in status_text:
            raise RuntimeError(
                f"BLAST job {rid} expired or is unknown. "
                f"RIDs expire after 24 hours."
            )
        elif "Status=READY" in status_text:
            print(f"  [{poll_num}/{max_polls}] Results ready.")
            break
    else:
        raise RuntimeError(
            f"BLAST job {rid} did not complete within "
            f"{max_polls} polls ({max_polls * POLL_INTERVAL // 60} "
            f"minutes). Increase max_polls or check NCBI status."
        )

    # Fetch the actual XML results
    result_response = requests.get(BLAST_URL, params={
        "CMD":          "Get",
        "RID":          rid,
        "FORMAT_TYPE":  "XML",
        "FORMAT_OBJECT":"Alignment",
    })
    result_response.raise_for_status()
    return result_response.text


def _parse_blast_xml(xml_text: str,
                     evalue_threshold: float = 1e-10,
                     min_pident: float = 70.0,
                     min_align_length: int = 44,
                     per_gene_thresholds: dict = None,
                     query_gene_map: dict = None) -> dict:
    """
    Parse BLAST XML results, applying per-gene or global significance
    filters to remove non-significant hits before scoring.

    Parameters
    ----------
    xml_text              : str    raw XML from NCBI BLAST
    evalue_threshold      : float  global maximum e-value. Used when
                                   no per-gene threshold is available.
    min_pident            : float  global minimum percent identity.
    min_align_length      : int    global minimum alignment length.
    per_gene_thresholds   : dict   from filter_thresholds.json. Keys
                                   are gene names, values are dicts
                                   with recommended_evalue,
                                   recommended_pident,
                                   recommended_align_length.
                                   If provided, gene-specific thresholds
                                   override the global values.
    query_gene_map        : dict   maps query_id -> gene name so the
                                   parser knows which gene each query
                                   sequence belongs to. Built from the
                                   subsample FASTA headers.

    Returns
    -------
    dict keyed by query_id -> list of filtered hit dicts
    """
    from io import StringIO
    from Bio.Blast import NCBIXML

    results       = {}
    blast_records = list(NCBIXML.parse(StringIO(xml_text)))

    for record in blast_records:
        query_id = record.query.split()[0]
        hits     = []

        # Determine which gene this query belongs to
        gene = None
        if query_gene_map:
            gene = query_gene_map.get(query_id)
            # Also try partial match if BLAST truncated the ID
            if gene is None:
                for qid, g in query_gene_map.items():
                    if query_id in qid or qid[:30] in query_id:
                        gene = g
                        break

        # Get thresholds for this gene
        if (per_gene_thresholds and gene
                and gene in per_gene_thresholds):
            gene_thresh  = per_gene_thresholds[gene]
            ev_thresh    = gene_thresh.get("recommended_evalue",
                                           evalue_threshold)
            pid_thresh   = gene_thresh.get("recommended_pident",
                                           min_pident)
            aln_thresh   = gene_thresh.get("recommended_align_length",
                                           min_align_length)
        else:
            ev_thresh  = evalue_threshold
            pid_thresh = min_pident
            aln_thresh = min_align_length

        for alignment in record.alignments:
            for hsp in alignment.hsps:

                if hsp.expect > ev_thresh:
                    continue
                if hsp.align_length < aln_thresh:
                    continue

                pident = (hsp.identities / hsp.align_length * 100
                          if hsp.align_length > 0 else 0.0)

                if pident < pid_thresh:
                    continue

                hits.append({
                    "accession":    alignment.accession,
                    "description":  alignment.title,
                    "bitscore":     hsp.bits,
                    "pident":       pident,
                    "evalue":       hsp.expect,
                    "align_length": hsp.align_length,
                })

        results[query_id] = hits

    return results


# ─────────────────────────────────────────────────────────────────────────────
# CSV MANAGEMENT
# ─────────────────────────────────────────────────────────────────────────────

# All fields written to the output CSV — one row per BLAST hit
CSV_FIELDS = [
    # Parent sequence information
    "label",                # GENUINE or ARTEFACT
    "taxon",
    "gene",
    "accession",
    "orientation",
    "parent_length",
    "subsample_length",
    "n_samples_total",
    # Subsample information
    "sample",
    "subsample_start",
    "subsample_end",
    # BLAST hit information
    "hit_accession",
    "hit_description",
    "bitscore",
    "pident",
    "evalue",
    "align_length",
    # Scoring
    "species_llr",
    "gene_llr",
    "total_hit_llr",
    # Aggregate per-subsample
    "subsample_posterior",
    "subsample_total_llr",
    # Chimera indicators (computed per parent sequence once all
    # subsamples for that parent are available)
    "n_hits_this_subsample",
    "rid",                  # BLAST RID for traceability
    "fasta_file",           # which batch FASTA file this came from
]


def _load_existing_csv(csv_path: str) -> tuple:
    """
    Load an existing results CSV and return:
      - a set of (accession, gene, sample) tuples already recorded
        so we can avoid duplicating rows
      - the existing rows as a list of dicts

    Parameters
    ----------
    csv_path : str   path to existing CSV (may not exist yet)

    Returns
    -------
    tuple (existing_keys: set, existing_rows: list)
    """
    existing_keys = set()
    existing_rows = []

    if not Path(csv_path).exists():
        return existing_keys, existing_rows

    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            existing_rows.append(row)
            key = (
                row.get("accession", ""),
                row.get("gene", ""),
                row.get("sample", ""),
                row.get("hit_accession", ""),
            )
            existing_keys.add(key)

    print(f"  Existing CSV: {len(existing_rows)} rows already recorded")
    return existing_keys, existing_rows


def _append_to_csv(csv_path: str,
                   new_rows: list,
                   existing_keys: set,
                   existing_rows: list) -> int:
    """
    Append new rows to the results CSV, skipping any that are already
    present (identified by accession + gene + sample + hit_accession).

    Rewrites the entire file to ensure consistent column ordering.

    Parameters
    ----------
    csv_path      : str   path to CSV file
    new_rows      : list  list of dicts to add
    existing_keys : set   (accession, gene, sample, hit_accession)
                          tuples already in the CSV
    existing_rows : list  existing rows as list of dicts

    Returns
    -------
    int   number of new rows actually written
    """
    added = 0
    for row in new_rows:
        key = (
            row.get("accession", ""),
            row.get("gene", ""),
            str(row.get("sample", "")),
            row.get("hit_accession", ""),
        )
        if key not in existing_keys:
            existing_rows.append(row)
            existing_keys.add(key)
            added += 1

    # Rewrite entire file with all rows
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(existing_rows)

    return added


# ─────────────────────────────────────────────────────────────────────────────
# MAIN FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

def subsamples_to_blast(
        subsample_fasta: str,
        metadata_csv: str,
        taxonomy: dict,
        glossary: dict,
        results_csv: str,
        label: str,
        prior_prob: float = 0.9,
        pident_threshold: float = 90.0,
        k: float = 0.2,
        beta: float = 0.05,
        evalue: float = 1e-10,
        max_hits: int = 20,
        database: str = "nt",
        entrez_query: str = "Vertebrata[Organism]",
        max_polls: int = 60,
        filter_thresholds_json: str = None,    # ← new parameter
) -> dict:
    """
    Submit a subsample FASTA file (prepared by prepare_subsample_fastas)
    to NCBI remote BLAST, score each hit using score_hit(), and append
    results to a running CSV file.

    The CSV file accumulates results across multiple calls — one per
    batch FASTA file — so the full calibration dataset can be built up
    incrementally. Duplicate rows (same accession + gene + sample +
    hit_accession) are never written twice, which handles the case where
    a parent sequence's subsamples are split across multiple batch files.

    BP limits
    ---------
    NCBI remote BLAST has a practical limit of ~1,000,000bp per
    submission. This function warns at 800,000bp (80% capacity) and
    raises an error above 900,000bp to prevent failed submissions.

    Parameters
    ----------
    subsample_fasta  : str   path to a batch FASTA file written by
                             prepare_subsample_fastas(). Each sequence
                             in this file is one subsample.
    metadata_csv     : str   path to the metadata CSV written by
                             prepare_subsample_fastas() for the same
                             parent FASTA. Used to look up parent
                             sequence details (taxon, seq_length etc.)
                             by (accession, gene) key.
    taxonomy         : dict  from load_taxonomy()
    glossary         : dict  from load_glossary()
    results_csv      : str   path to the running results CSV. Created
                             if it does not exist; appended to if it
                             does. All calibration data accumulates here
                             across multiple calls.
    label            : str   "GENUINE" or "ARTEFACT" — recorded in the
                             CSV so the calibration function can
                             distinguish the two classes.
    prior_prob       : float prior P(genuine) for aggregate_hit_scores()
    pident_threshold : float pident contributing LLR of 0
    k                : float pident LLR steepness
    beta             : float softmax temperature over bitscores
    evalue           : float BLAST e-value threshold
    max_hits         : int   maximum BLAST hits per query sequence
    database         : str   BLAST database (default "nt")
    entrez_query     : str   Entrez filter. Default "Vertebrata[Organism]"
                             restricts to vertebrates which covers all
                             Eulipotyphla and common contaminants.
    max_polls        : int   maximum status poll attempts at 2min each.
                             Default 60 = 2 hours.
    filter_thresholds_json : str   optional path to filter_thresholds.json
                                   from calibrate_hit_filters(). If provided,
                                   per-gene evalue, pident and align_length
                                   thresholds are applied when filtering BLAST
                                   results. If None, global evalue and
                                   perc_identity parameters are used for all
                                   genes.

    Returns
    -------
    dict with:
        rid              : str   BLAST Request ID for traceability
        n_queries        : int   number of subsample sequences submitted
        total_bp         : int   total base pairs submitted
        n_hits_total     : int   total BLAST hits across all queries
        n_rows_added     : int   new rows written to results_csv
        n_rows_skipped   : int   rows skipped (already in CSV)
        results_csv      : str   path to the results CSV
        parent_summaries : dict  (accession, gene) -> per-parent summary
    """
    _start = time.time()

    # ── Load per-gene filter thresholds if provided ───────────────────────────
    per_gene_thresholds = None
    if filter_thresholds_json and os.path.exists(filter_thresholds_json):
        with open(filter_thresholds_json, "r") as f:
            thresh_data         = json.load(f)
            per_gene_thresholds = thresh_data.get("per_gene_thresholds")
        print(f"  Loaded per-gene thresholds for "
              f"{len(per_gene_thresholds)} genes from "
              f"{filter_thresholds_json}")
    elif filter_thresholds_json:
        print(f"  Warning: filter_thresholds_json not found at "
              f"{filter_thresholds_json} — using global thresholds")
        
    if isinstance(taxonomy, str):
        raise TypeError(
            "taxonomy must be a dict from load_taxonomy(), not a path."
        )
    if label not in ("GENUINE", "ARTEFACT"):
        raise ValueError(
            f"label must be 'GENUINE' or 'ARTEFACT', got '{label}'"
        )

    # ── Load subsample FASTA ──────────────────────────────────────────────────
    records = list(SeqIO.parse(subsample_fasta, "fasta"))
    if not records:
        raise ValueError(
            f"No sequences found in {subsample_fasta}"
        )

    print(f"\n{'='*60}")
    print(f"subsamples_to_blast()")
    print(f"{'='*60}")
    print(f"  Input FASTA  : {subsample_fasta}")
    print(f"  Label        : {label}")
    print(f"  Sequences    : {len(records)}")
    print(f"  Results CSV  : {results_csv}")

    # ── Check total base pairs ────────────────────────────────────────────────
    total_bp = sum(len(str(r.seq)) for r in records)
    print(f"  Total bp     : {total_bp:,}")

    if total_bp > 900_000:
        raise ValueError(
            f"Total base pairs ({total_bp:,}) exceeds 900,000bp safety "
            f"limit for NCBI remote BLAST submissions. Split your "
            f"subsample FASTA into smaller files using "
            f"prepare_subsample_fastas(max_seqs_per_file=...)."
        )
    elif total_bp > 800_000:
        print(f"  WARNING: {total_bp:,}bp is above 80% of the "
              f"recommended 1,000,000bp NCBI limit. Consider splitting "
              f"into smaller files to avoid timeouts.")

    # ── Load metadata ─────────────────────────────────────────────────────────
    metadata = load_metadata(metadata_csv)

    # ── Load existing CSV rows to avoid duplication ───────────────────────────
    existing_keys, existing_rows = _load_existing_csv(results_csv)
    print(f"  Existing CSV rows: {len(existing_rows)}")

    # ── Read FASTA content for submission ─────────────────────────────────────
    with open(subsample_fasta, "r") as f:
        fasta_content = f.read()

    # ── Submit to NCBI BLAST ──────────────────────────────────────────────────
    print(f"\nSubmitting to NCBI BLAST ({database}, "
          f"filter: '{entrez_query}')...")
    rid = _submit_blast(
        fasta_content = fasta_content,
        evalue        = evalue,
        max_hits      = max_hits,
        database      = database,
        entrez_query  = entrez_query,
    )

    # ── Poll for results ──────────────────────────────────────────────────────
    xml_text = _poll_blast(rid=rid, max_polls=max_polls)

    # ── Build query_gene_map from subsample headers ───────────────────────────
    # Maps each subsample sequence ID to its gene name so _parse_blast_xml
    # can apply the correct per-gene threshold to each query result
    query_gene_map = {}
    for record in records:
        sub_info = parse_subsample_header(record.id)
        if sub_info.get("gene"):
            query_gene_map[record.id] = sub_info["gene"]

    # ── Parse XML with per-gene thresholds ───────────────────────────────────
    # Replace the existing _parse_blast_xml call with this:
    print(f"\nParsing BLAST XML results...")
    blast_results = _parse_blast_xml(
        xml_text            = xml_text,
        evalue_threshold    = evalue,
        min_pident          = 70.0,     # permissive global fallback —
                                        # per-gene thresholds handle
                                        # gene-specific filtering
        min_align_length    = 44,
        per_gene_thresholds = per_gene_thresholds,
        query_gene_map      = query_gene_map,
    )
    print(f"  Query sequences with results: {len(blast_results)}")

    # ── Score hits and build CSV rows ─────────────────────────────────────────
    # Group subsamples by (accession, gene) so we can compute per-parent
    # aggregate statistics once all subsamples for a parent are available.

    print(f"\nScoring hits...")
    new_rows        = []
    n_hits_total    = 0

    # Track per-parent data for aggregate metrics
    parent_subsample_data = defaultdict(list)

    fasta_file_name = Path(subsample_fasta).name

    for record in records:
        sub_info = parse_subsample_header(record.id)

        accession = sub_info["accession"]
        gene      = sub_info["gene"]
        taxon     = sub_info["taxon"]
        sample    = sub_info["sample"]

        # Look up parent metadata by (accession, gene)
        meta_key = (accession, gene)
        parent   = metadata.get(meta_key)

        if parent is None:
            print(f"  Warning: no metadata for {accession} | {gene} — "
                  f"using header values only")
            parent = {
                "taxon":            taxon,
                "gene":             gene,
                "accession":        accession,
                "orientation":      sub_info["orientation"],
                "seq_length":       sub_info["parent_length"] or 0,
                "subsample_length": (
                    (sub_info["subsample_end"] or 0) -
                    (sub_info["subsample_start"] or 0)
                ),
                "n_samples":        0,
            }

        # Get BLAST hits for this subsample
        # The query ID in BLAST results may be truncated — match by
        # checking if the record ID starts with the blast result key
        hits_raw = blast_results.get(record.id, [])
        if not hits_raw:
            # Try partial match — BLAST sometimes truncates long IDs
            for blast_id, hits in blast_results.items():
                if record.id.startswith(blast_id) or \
                        blast_id.startswith(record.id[:30]):
                    hits_raw = hits
                    break

        n_hits_total += len(hits_raw)

        if not hits_raw:
            print(f"  No hits: {taxon} | {gene} | "
                  f"{accession} | Sample {sample}")
            # Write a row with no-hit information so the CSV is
            # complete and calibration can account for missing hits
            new_rows.append({
                "label":               label,
                "taxon":               taxon,
                "gene":                gene,
                "accession":           accession,
                "orientation":         parent["orientation"],
                "parent_length":       parent["seq_length"],
                "subsample_length":    parent["subsample_length"],
                "n_samples_total":     parent["n_samples"],
                "sample":              sample,
                "subsample_start":     sub_info["subsample_start"],
                "subsample_end":       sub_info["subsample_end"],
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
                "fasta_file":          fasta_file_name,
            })
            continue

        # Score hits using score_hit()
        scored_hits = [
            score_hit(hit, taxon, gene, taxonomy, glossary)
            for hit in hits_raw
        ]

        # Aggregate across hits for this subsample
        subsample_agg = aggregate_hit_scores(
            scored_hits,
            prior_prob       = prior_prob,
            pident_threshold = pident_threshold,
            k                = k,
            beta             = beta,
        )

        # Store for per-parent aggregate computation
        parent_subsample_data[meta_key].append({
            "sample":           sample,
            "posterior":        subsample_agg["posterior_prob"],
            "total_llr":        subsample_agg["total_llr"],
            "scored_hits":      scored_hits,
        })

        # Write one row per hit (not per subsample) so every piece
        # of information is available for calibration
        for hit, scored in zip(hits_raw, scored_hits):
            new_rows.append({
                "label":               label,
                "taxon":               taxon,
                "gene":                gene,
                "accession":           accession,
                "orientation":         parent["orientation"],
                "parent_length":       parent["seq_length"],
                "subsample_length":    parent["subsample_length"],
                "n_samples_total":     parent["n_samples"],
                "sample":              sample,
                "subsample_start":     sub_info["subsample_start"],
                "subsample_end":       sub_info["subsample_end"],
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
                "fasta_file":          fasta_file_name,
            })

    # ── Append to CSV ─────────────────────────────────────────────────────────
    print(f"\nUpdating results CSV...")
    n_rows_added = _append_to_csv(
        csv_path      = results_csv,
        new_rows      = new_rows,
        existing_keys = existing_keys,
        existing_rows = existing_rows,
    )
    n_rows_skipped = len(new_rows) - n_rows_added

    print(f"  New rows added   : {n_rows_added}")
    print(f"  Rows skipped     : {n_rows_skipped} (already in CSV)")
    print(f"  Total CSV rows   : {len(existing_rows)}")

    # ── Per-parent summaries ──────────────────────────────────────────────────
    import numpy as np
    import math

    parent_summaries = {}
    for meta_key, sub_data in parent_subsample_data.items():
        posteriors  = [s["posterior"] for s in sub_data]
        total_llrs  = [s["total_llr"] for s in sub_data]
        all_tax_llrs = [
            h["species_llr"]
            for s in sub_data
            for h in s["scored_hits"]
        ]

        post_var    = float(np.var(posteriors)) if len(posteriors) > 1 \
                      else 0.0
        llr_var     = float(np.var(all_tax_llrs)) if all_tax_llrs \
                      else 0.0
        prop_neg    = (sum(1 for x in all_tax_llrs if x < 0) /
                       len(all_tax_llrs) if all_tax_llrs else 0.0)

        parent_summaries[meta_key] = {
            "taxon":               metadata.get(meta_key, {}).get(
                "taxon", meta_key[0]),
            "accession":           meta_key[0],
            "gene":                meta_key[1],
            "n_subsamples_scored": len(sub_data),
            "mean_posterior":      round(float(np.mean(posteriors)), 4),
            "posterior_variance":  round(post_var, 6),
            "taxon_llr_variance":  round(llr_var, 4),
            "prop_neg_taxon_hits": round(prop_neg, 4),
            "chimera_flag":        post_var > 0.05 or llr_var > 2.0,
        }

    # ── Print summary ─────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"COMPLETE")
    print(f"{'='*60}")
    print(f"  RID                  : {rid}")
    print(f"  Queries submitted    : {len(records)}")
    print(f"  Total bp submitted   : {total_bp:,}")
    print(f"  Total BLAST hits     : {n_hits_total}")
    print(f"  New CSV rows added   : {n_rows_added}")
    print(f"  Parent seqs scored   : {len(parent_summaries)}")

    if parent_summaries:
        flagged = [v for v in parent_summaries.values()
                   if v["chimera_flag"]]
        print(f"  Chimera flagged      : {len(flagged)}/"
              f"{len(parent_summaries)}")
        
    elapsed = time.time() - _start
    print(f"\nTotal runtime: {int(elapsed // 60)}m {int(elapsed % 60)}s")


    return {
        "rid":              rid,
        "n_queries":        len(records),
        "total_bp":         total_bp,
        "n_hits_total":     n_hits_total,
        "n_rows_added":     n_rows_added,
        "n_rows_skipped":   n_rows_skipped,
        "results_csv":      results_csv,
        "parent_summaries": parent_summaries,
    }



#subsamples_to_blast(
#    subsample_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/genuine/subsample_batch_001.fasta",
#    metadata_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/genuine/subsample_batch_metadata.csv",
#    taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
#    glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"),
#    results_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
#    label = "GENUINE",
#)

subsamples_to_blast(
    subsample_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/artefact/subsample_batch_001.fasta",
    metadata_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/artefact/subsample_batch_metadata.csv",
    taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
    glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"),
    results_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
    label = "ARTEFACT",
    filter_thresholds_json= "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/filter_thresholds.json",
)

#from bioinf_packages.verify_funcs._compute_csv import calibrate_hit_filters
#calibrate_hit_filters(
#    results_csv="C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
#    output_path="C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/filter_thresholds.json"
#)
"""
Note the following post running the above code:
Gene precision values: ATP6 (0.67), COX2 (0.60), COX3 (0.71), ND3 (0.61)
These genes have very low evalues and pidents after blasting, even when scoring a correct hit.
This means that they are likely to be flagged as chimeras unless a good taxon_llr score is returned
"""


#subsamples_to_blast(
#    subsample_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/query/subsample_batch_001.fasta",
#    metadata_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/query/subsample_batch_metadata.csv",
#    taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
#    glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"),
#    results_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/blast_results.csv",
#    label = "QUERY",
#)