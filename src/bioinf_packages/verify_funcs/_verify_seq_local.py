# _verify_seq.py

import os
import re
import json
import random
import subprocess
import tempfile
import numpy as np

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def run_blast_local(query_fasta: str,
                    db: str,
                    out_file: str,
                    blast_bin: str = "blastn",
                    evalue: float = 1e-10,
                    max_target_seqs: int = 20,
                    num_threads: int = 4,
                    perc_identity: float = 70.0) -> str:
    """
    Run BLASTN locally against a local database (e.g. core_nt).

    Parameters
    ----------
    query_fasta     : str   path to query FASTA file
    db              : str   path to local BLAST database
                            e.g. "C:/blast/db/core_nt"
    out_file        : str   path for BLAST output (XML format)
    blast_bin       : str   path to blastn binary
                            e.g. "C:/blast/bin/blastn.exe"
    evalue          : float e-value threshold
    max_target_seqs : int   maximum hits to return per query
    num_threads     : int   CPU threads for BLAST
    perc_identity   : float minimum percent identity threshold

    Returns
    -------
    str   path to output XML file
    """
    cmd = [
        blast_bin,
        "-query",          query_fasta,
        "-db",             db,
        "-out",            out_file,
        "-outfmt",         "5",        # XML format for parsing
        "-evalue",         str(evalue),
        "-max_target_seqs",str(max_target_seqs),
        "-num_threads",    str(num_threads),
        "-perc_identity",  str(perc_identity),
        "-dust",           "no",       # disable dust filter for
                                       # short subsamples
    ]

    print(f"  Running BLAST: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"BLAST failed:\n{result.stderr}"
        )

    if not os.path.exists(out_file):
        raise FileNotFoundError(
            f"BLAST finished but output not found: {out_file}"
        )

    return out_file


def verify_seq_local(query_seq: str,
                     query_taxon: str,
                     query_gene: str,
                     query_accession: str,
                     query_orientation: str,
                     blast_db: str,
                     output_dir: str,
                     taxonomy: dict,
                     glossary: dict,
                     blast_bin: str = "blastn",
                     n_samples: int = None,
                     insert_fraction: float = 0.30,
                     target_detection_prob: float = 0.90,
                     evalue: float = 1e-10,
                     max_target_seqs: int = 20,
                     num_threads: int = 4,
                     perc_identity: float = 70.0,
                     flag_on_single_hit: bool = True,
                     single_hit_llr_threshold: float = -2.0) -> dict:
    """
    Verify a sequence by BLAST searching subsamples against a local
    BLAST database.

    With flag_on_single_hit=True (default), a single subsample
    returning a hit with taxon LLR below single_hit_llr_threshold
    flags the entire sequence as questionable. This matches the
    subsampling logic — 90% probability that at least one subsample
    overlaps a chimeric junction — so a single bad hit is meaningful
    signal, not noise.

    Parameters
    ----------
    query_seq               : str   full query sequence
    query_taxon             : str   e.g. "Solenodon_paradoxus"
    query_gene              : str   e.g. "12S_rRNA"
    query_accession         : str   e.g. "AY530070"
    query_orientation       : str   "+" or "-"
    blast_db                : str   path to local BLAST database
    output_dir              : str   directory to save output files
    taxonomy                : dict  from load_taxonomy()
    glossary                : dict  from load_glossary()
    blast_bin               : str   path to blastn executable
    n_samples               : int   number of subsamples. If None,
                                    calculated from recommended_n_samples()
    insert_fraction         : float expected chimeric insert fraction
    target_detection_prob   : float desired P(detect chimera)
    evalue                  : float BLAST e-value threshold
    max_target_seqs         : int   max BLAST hits per subsample
    num_threads             : int   BLAST CPU threads
    perc_identity           : float minimum percent identity for hits
    flag_on_single_hit      : bool  if True, a single subsample with
                                    taxon LLR below threshold flags the
                                    whole sequence. Default True.
    single_hit_llr_threshold: float taxon LLR below which a hit is
                                    considered taxonomically wrong.
                                    Default -2.0.

    Returns
    -------
    dict with per-sample BLAST results and sequence-level summary
    """
    from ._score_hit import (score_blast, recommended_subsample_length,
                              recommended_n_samples)

    os.makedirs(output_dir, exist_ok=True)

    seq_length = len(query_seq)
    sub_length = recommended_subsample_length(
        seq_length, insert_fraction=insert_fraction
    )

    if n_samples is None:
        n_samples = recommended_n_samples(
            seq_length            = seq_length,
            subsample_length      = sub_length,
            target_detection_prob = target_detection_prob,
            insert_fraction       = insert_fraction,
        )

    print(f"\nVerifying {query_taxon} | {query_gene} | "
          f"{query_accession}")
    print(f"Sequence: {seq_length}bp | "
          f"Subsample: {sub_length}bp | Samples: {n_samples}")
    print(f"BLAST db: {blast_db}")

    # ── Draw subsamples ───────────────────────────────────────────────────────
    subsamples  = []
    sample_meta = []

    for i in range(1, n_samples + 1):
        start  = random.randint(0, seq_length - sub_length)
        sample = query_seq[start:start + sub_length]
        header = (f"{query_taxon}|{query_gene}|"
                  f"orien:{query_orientation}|"
                  f"accession:{query_accession}|Sample{i}")
        subsamples.append(
            SeqRecord(Seq(sample), id=header, description="")
        )
        sample_meta.append({
            "sample": i,
            "start":  start,
            "end":    start + sub_length,
        })

    # ── Write subsamples to FASTA ─────────────────────────────────────────────
    query_fasta = os.path.join(
        output_dir,
        f"{query_taxon}_{query_accession}_subsamples.fasta"
    )
    SeqIO.write(subsamples, query_fasta, "fasta")

    # ── Run BLAST ─────────────────────────────────────────────────────────────
    blast_out = os.path.join(
        output_dir,
        f"{query_taxon}_{query_accession}_blast.xml"
    )
    run_blast_local(
        query_fasta     = query_fasta,
        db              = blast_db,
        out_file        = blast_out,
        blast_bin       = blast_bin,
        evalue          = evalue,
        max_target_seqs = max_target_seqs,
        num_threads     = num_threads,
        perc_identity   = perc_identity,
    )

    # ── Parse and score BLAST results ─────────────────────────────────────────
    blast_results = score_blast(
        blast_xml       = blast_out,
        query_taxon     = query_taxon,
        query_gene      = query_gene,
        taxonomy        = taxonomy,
        glossary        = glossary,
        sample_metadata = sample_meta,
    )

    # ── Apply single-hit flagging logic ───────────────────────────────────────
    # A single subsample with a clearly wrong taxonomic hit flags
    # the whole sequence. This is justified by the subsampling design
    # which gives 90% probability that at least one subsample
    # overlaps a chimeric region if one exists.
    #
    # Without this logic, a chimeric sequence could pass verification
    # if only 1-2 of 10 subsamples land in the contaminated region,
    # because the mean LLR across all subsamples would still look
    # acceptable.

    flagged_samples = []

    if flag_on_single_hit:
        for sample_result in blast_results.get("sample_scores", []):
            taxon_llr = sample_result.get("taxon_llr")
            if taxon_llr is not None and taxon_llr < single_hit_llr_threshold:
                flagged_samples.append(sample_result)
                print(f"  ⚠ Sample {sample_result['sample']} flagged: "
                      f"taxon_llr={taxon_llr:.3f} < "
                      f"{single_hit_llr_threshold} | "
                      f"top_hit={sample_result.get('top_taxon', 'N/A')} | "
                      f"pos={sample_result.get('start', '?')}-"
                      f"{sample_result.get('end', '?')}")

    # chimera_flag: True if ANY subsample has a bad taxonomic hit
    # (when flag_on_single_hit=True) OR if the aggregate score
    # indicates a chimera (existing logic)
    single_hit_flag = len(flagged_samples) > 0

    chimera_flag = (
        single_hit_flag
        or blast_results.get("chimera_flag", False)
    )

    # ── Build summary ─────────────────────────────────────────────────────────
    summary = {
        **blast_results,
        "query_taxon":              query_taxon,
        "query_gene":               query_gene,
        "query_accession":          query_accession,
        "query_orientation":        query_orientation,
        "n_samples":                n_samples,
        "subsample_length":         sub_length,
        "flag_on_single_hit":       flag_on_single_hit,
        "single_hit_llr_threshold": single_hit_llr_threshold,
        "n_flagged_samples":        len(flagged_samples),
        "flagged_sample_indices":   [s["sample"]
                                     for s in flagged_samples],
        "chimera_flag":             chimera_flag,
        "single_hit_flag":          single_hit_flag,
    }

    print(f"\n── BLAST summary ────────────────────────────────────────────")
    print(f"  Taxon             : {query_taxon}")
    print(f"  Gene              : {query_gene}")
    print(f"  Samples flagged   : "
          f"{len(flagged_samples)}/{n_samples}")
    print(f"  Single hit flag   : {single_hit_flag}")
    print(f"  Chimera flag      : {chimera_flag}")
    if flagged_samples:
        print(f"  Flagged positions :")
        for s in flagged_samples:
            print(f"    Sample {s['sample']} | "
                  f"pos {s.get('start','?')}-{s.get('end','?')} | "
                  f"taxon_llr={s.get('taxon_llr', 'N/A'):.3f} | "
                  f"top_hit={s.get('top_taxon', 'N/A')}")

    return summary


def verify_batch_local(query_fasta: str,
                        blast_db: str,
                        output_dir: str,
                        taxonomy: dict,
                        glossary: dict,
                        blast_bin: str = "blastn",
                        insert_fraction: float = 0.30,
                        target_detection_prob: float = 0.90,
                        evalue: float = 1e-10,
                        max_target_seqs: int = 20,
                        num_threads: int = 4,
                        perc_identity: float = 70.0,
                        flag_on_single_hit: bool = True,
                        single_hit_llr_threshold: float = -2.0) -> list:
    """
    Run verify_seq_local() on every sequence in a FASTA file.

    FASTA header format:
        >Species_name|gene|orien:+/-|accession:XXXX

    Parameters
    ----------
    query_fasta             : str   FASTA of sequences to verify
    blast_db                : str   path to local BLAST database
    output_dir              : str   parent directory for outputs
    taxonomy                : dict  from load_taxonomy()
    glossary                : dict  from load_glossary()
    blast_bin               : str   path to blastn executable
    insert_fraction         : float expected chimeric insert fraction
    target_detection_prob   : float desired P(detect chimera)
    evalue                  : float BLAST e-value threshold
    max_target_seqs         : int   max BLAST hits per subsample
    num_threads             : int   BLAST CPU threads
    perc_identity           : float minimum percent identity
    flag_on_single_hit      : bool  flag whole sequence on single bad
                                    subsample hit
    single_hit_llr_threshold: float taxon LLR threshold for bad hit

    Returns
    -------
    list of summary dicts, one per sequence
    """
    import csv

    os.makedirs(output_dir, exist_ok=True)

    records = list(SeqIO.parse(query_fasta, "fasta"))
    total   = len(records)
    results = []

    print(f"\nBatch local BLAST verification: {total} sequences")
    print(f"Database : {blast_db}")
    print(f"Output   : {output_dir}")
    print(f"Single-hit flagging: {flag_on_single_hit} "
          f"(threshold={single_hit_llr_threshold})\n")

    for i, record in enumerate(records, 1):
        parts       = record.description.split("|")
        taxon       = parts[0].strip()
        gene        = (parts[1].strip()
                       if len(parts) > 1 else "unknown")
        orientation = "+"
        accession   = None

        for p in parts:
            if p.startswith("orien:"):
                orientation = p.split(":", 1)[1].strip()
            if p.startswith("accession:"):
                accession = p.split(":", 1)[1].strip()

        if not accession:
            print(f"  [{i}/{total}] Skipping {taxon} — "
                  f"no accession in header")
            continue

        print(f"\n[{i}/{total}] {taxon} | {gene} | {accession}")

        # Create per-sequence output directory
        seq_dir = os.path.join(
            output_dir,
            f"{taxon}_{accession}"
        )
        os.makedirs(seq_dir, exist_ok=True)

        try:
            result = verify_seq_local(
                query_seq               = str(record.seq),
                query_taxon             = taxon,
                query_gene              = gene,
                query_accession         = accession,
                query_orientation       = orientation,
                blast_db                = blast_db,
                output_dir              = seq_dir,
                taxonomy                = taxonomy,
                glossary                = glossary,
                blast_bin               = blast_bin,
                insert_fraction         = insert_fraction,
                target_detection_prob   = target_detection_prob,
                evalue                  = evalue,
                max_target_seqs         = max_target_seqs,
                num_threads             = num_threads,
                perc_identity           = perc_identity,
                flag_on_single_hit      = flag_on_single_hit,
                single_hit_llr_threshold= single_hit_llr_threshold,
            )
            result["_header"] = record.description
            results.append(result)

            flag_str = "⚠ CHIMERA" if result["chimera_flag"] else "✓ OK"
            print(f"  {flag_str} | "
                  f"flagged_samples={result['n_flagged_samples']}/"
                  f"{result['n_samples']} | "
                  f"single_hit_flag={result['single_hit_flag']}")

        except Exception as e:
            print(f"  ✗ Failed: {e}")
            import traceback
            traceback.print_exc()

    # ── Write summary CSV ─────────────────────────────────────────────────────
    if results:
        summary_path = os.path.join(
            output_dir, "blast_verification_summary.csv"
        )
        exclude  = {"sample_scores", "hit_details"}
        csv_keys = [k for k in results[0].keys()
                    if k not in exclude]

        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=csv_keys,
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(results)

        print(f"\nSummary CSV: {summary_path}")

        flagged = [r for r in results if r["chimera_flag"]]
        print(f"\n{'⚠' if flagged else '✓'} "
              f"{len(flagged)}/{len(results)} sequences flagged")
        for r in flagged:
            print(f"  {r['query_taxon']} "
                  f"({r['query_accession']}) | "
                  f"{r.get('query_gene','?')} | "
                  f"flagged_samples="
                  f"{r['n_flagged_samples']}/{r['n_samples']}")

    return results


#from bioinf_packages.verify_funcs._verify_seq import verify_batch_local
from bioinf_packages.verify_funcs._species_parser import load_taxonomy
from bioinf_packages.verify_funcs._gene_parser import load_glossary

taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv")
import os
os.environ["BLASTDB"] = r"D:\blastdb\core_nt"

run_blast_local(query_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/raw_fastas/KF783046.fasta",
                    db = "core_nt",
                    out_file = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/test1.xml",
                    blast_bin = "C:/Program Files/NCBI/blast-2.17.0+/bin/blastn.exe",
                    evalue = 1e-10,
                    max_target_seqs = 20,
                    num_threads = 4,
                    perc_identity = 70.0) 

#results = verify_batch_local(
#    query_fasta              = "sequences_to_verify.fasta",
#    blast_db                 = "C:/blast/db/core_nt",
#    output_dir               = "blast_verification_output",
#    taxonomy                 = taxonomy,
#    glossary                 = glossary,
#    blast_bin                = "C:/blast/bin/blastn.exe",
#    num_threads              = 8,
#    flag_on_single_hit       = True,
#    single_hit_llr_threshold = -2.0,
#)

# Sequences that pass BLAST can be used as gold standard
# for distance verification calibration
#passed = [r for r in results if not r["chimera_flag"]]
#print(f"{len(passed)}/{len(results)} sequences passed BLAST verification")