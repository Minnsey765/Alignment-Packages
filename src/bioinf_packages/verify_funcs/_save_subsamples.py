# _prepare_subsamples.py

import os
import csv
import random

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# recommended_subsample_length — calculates how long each subsample
# should be based on sequence length and expected insert fraction
# recommended_n_samples — calculates how many subsamples to draw to
# achieve the target detection probability
from bioinf_packages.verify_funcs._score_hit import (
    recommended_subsample_length,
    recommended_n_samples,
)

def prepare_subsample_fastas(
        input_fasta: str,
        output_dir: str,
        max_seqs_per_file: int,
        max_files: int = None,
        insert_fraction: float = 0.30,
        target_detection_prob: float = 0.90,
        file_prefix: str = "subsample_batch",
) -> dict:
    """
    Draw subsamples from every sequence in a FASTA file and write them
    to batch FASTA files for remote BLAST querying.

    Each output file contains at most max_seqs_per_file subsample
    sequences, allowing you to control how large each BLAST query is.
    The number of subsamples drawn per parent sequence is calculated
    automatically from recommended_n_samples() based on sequence length.

    Output files are overwritten on each run to avoid duplication.

    Subsample headers contain all parent sequence information plus a
    sample index, in the format:
        >Genus_species|gene|orien:+/-|accession:XXXX|Sample:N|
         parent_length:L|subsample_start:S|subsample_end:E

    This means the output FASTA can be queried directly with BLAST and
    the results can be traced back to the parent sequence and subsample
    position without needing a separate metadata file.

    Parameters
    ----------
    input_fasta        : str   path to input FASTA file. Headers must
                               follow the pipe-delimited format:
                               >Genus_species|gene|orien:+/-|accession:XXXX
    output_dir         : str   directory to write subsample FASTA files.
                               Created if it does not exist. Existing
                               files matching file_prefix are overwritten.
    max_seqs_per_file  : int   maximum number of subsample sequences per
                               output file. Controls the size of each
                               BLAST query. Default 50.
    max_files          : int   maximum number of output files to write.
                               If None, processes all sequences in the
                               input FASTA. Use a small number (e.g. 2)
                               to test query speed before committing to
                               the full dataset.
    insert_fraction    : float expected chimeric insert as fraction of
                               sequence length. Passed to
                               recommended_subsample_length() and
                               recommended_n_samples(). Default 0.30.
    target_detection_prob: float desired P(detect chimera junction).
                               Passed to recommended_n_samples().
                               Default 0.90.
    file_prefix        : str   prefix for output filenames. Files are
                               named {file_prefix}_001.fasta,
                               {file_prefix}_002.fasta etc.
                               Default "subsample_batch".

    Returns
    -------
    dict with:
        output_files      : list  paths to all written FASTA files
        n_files           : int   number of files written
        n_parent_seqs     : int   number of parent sequences processed
        n_subsamples_total: int   total subsample sequences written
        skipped           : list  parent sequences skipped (too short
                                  or missing header fields)
        metadata          : list  one dict per parent sequence with:
                                    taxon, gene, accession,
                                    seq_length, subsample_length,
                                    n_samples, file_indices
                                  (which batch files contain its
                                  subsamples)
    """
    import random
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    os.makedirs(output_dir, exist_ok=True)

    # ── Delete any existing files matching the prefix to avoid duplication ────
    existing = list(Path(output_dir).glob(f"{file_prefix}_*.fasta"))
    for f in existing:
        f.unlink()
    if existing:
        print(f"Deleted {len(existing)} existing file(s) matching "
              f"'{file_prefix}_*.fasta'")

    # ── Parse input FASTA ─────────────────────────────────────────────────────
    records = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"\nInput FASTA        : {input_fasta}")
    print(f"Parent sequences   : {len(records)}")
    print(f"Max seqs per file  : {max_seqs_per_file}")
    print(f"Max files          : {max_files or 'all'}")
    print(f"Output directory   : {output_dir}")

    # ── Helper: parse pipe-delimited header ───────────────────────────────────
    def parse_header(description: str) -> dict:
        parts       = description.strip().split("|")
        taxon       = parts[0].strip()
        gene        = parts[1].strip() if len(parts) > 1 else None
        orientation = "+"
        accession   = None
        for p in parts:
            if p.startswith("orien:"):
                orientation = p.split(":", 1)[1].strip()
            if p.startswith("accession:"):
                accession = p.split(":", 1)[1].strip()
        return {
            "taxon":       taxon,
            "gene":        gene,
            "orientation": orientation,
            "accession":   accession,
        }

    # ── Generate all subsamples in order ──────────────────────────────────────
    # First build the full list of SeqRecord objects so we can split
    # them into files of max_seqs_per_file without partial parent seqs
    # being split awkwardly across files.

    all_subsamples = []   # list of SeqRecord
    metadata       = []   # one dict per parent sequence
    skipped        = []

    for record in records:
        info = parse_header(record.description)

        if not info["gene"] or not info["accession"]:
            print(f"  Skipping {record.description[:60]} "
                  f"— missing gene or accession in header")
            skipped.append(record.description)
            continue

        seq        = str(record.seq).replace("-", "")
        seq_length = len(seq)

        if seq_length < 50:
            print(f"  Skipping {info['taxon']} | {info['gene']} | "
                  f"{info['accession']} — sequence too short "
                  f"({seq_length}bp)")
            skipped.append(record.description)
            continue

        # Calculate recommended subsample parameters
        sub_length = recommended_subsample_length(
            seq_length, insert_fraction=insert_fraction
        )
        n_samples  = recommended_n_samples(
            seq_length            = seq_length,
            subsample_length      = sub_length,
            target_detection_prob = target_detection_prob,
            insert_fraction       = insert_fraction,
        )

        print(f"  {info['taxon']} | {info['gene']} | "
              f"{info['accession']} | "
              f"{seq_length}bp → {n_samples} subsamples "
              f"of {sub_length}bp")

        # Draw subsamples
        sample_records = []
        for i in range(1, n_samples + 1):
            start  = random.randint(0, seq_length - sub_length)
            sample = seq[start:start + sub_length]

            # Header encodes all information needed to trace this
            # subsample back to its parent and position
            header = (
                f"{info['taxon']}|{info['gene']}|"
                f"orien:{info['orientation']}|"
                f"accession:{info['accession']}|"
                f"Sample:{i}|"
                f"parent_length:{seq_length}|"
                f"subsample_start:{start}|"
                f"subsample_end:{start + sub_length}"
            )
            sample_records.append(
                SeqRecord(Seq(sample), id=header, description="")
            )

        all_subsamples.extend(sample_records)

        metadata.append({
            "taxon":            info["taxon"],
            "gene":             info["gene"],
            "accession":        info["accession"],
            "orientation":      info["orientation"],
            "seq_length":       seq_length,
            "subsample_length": sub_length,
            "n_samples":        n_samples,
            "file_indices":     [],  # filled in below
        })

    print(f"\nTotal subsamples generated : {len(all_subsamples)}")
    print(f"Splitting into files of    : {max_seqs_per_file} sequences")

    # ── Split into files ──────────────────────────────────────────────────────
    output_files     = []
    n_subsamples_out = 0
    file_number      = 1

    for batch_start in range(0, len(all_subsamples), max_seqs_per_file):

        # Stop if max_files reached
        if max_files is not None and file_number > max_files:
            print(f"\nReached max_files={max_files} — stopping.")
            break

        batch = all_subsamples[batch_start:batch_start + max_seqs_per_file]

        filename  = f"{file_prefix}_{file_number:03d}.fasta"
        filepath  = os.path.join(output_dir, filename)

        SeqIO.write(batch, filepath, "fasta")
        output_files.append(filepath)

        print(f"  Written: {filename} "
              f"({len(batch)} subsample sequences)")

        n_subsamples_out += len(batch)
        file_number      += 1

    # ── Fill in file_indices in metadata ─────────────────────────────────────
    # For each parent sequence, record which batch file(s) contain
    # its subsamples — useful for tracking results back to parents
    subsample_idx = 0
    for meta in metadata:
        n = meta["n_samples"]
        file_indices = set()
        for j in range(n):
            global_idx = subsample_idx + j
            # Which file does this subsample end up in?
            file_idx   = global_idx // max_seqs_per_file
            if max_files is None or file_idx < max_files:
                file_indices.add(file_idx + 1)  # 1-indexed
        meta["file_indices"] = sorted(file_indices)
        subsample_idx += n

    # ── Write metadata CSV ────────────────────────────────────────────────────
    metadata_path = os.path.join(output_dir, f"{file_prefix}_metadata.csv")
    import csv
    if metadata:
        with open(metadata_path, "w", newline="") as f:
            keys = list(metadata[0].keys())
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for row in metadata:
                # Convert list to string for CSV
                row_out = dict(row)
                row_out["file_indices"] = str(row["file_indices"])
                writer.writerow(row_out)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_parent_processed = len(metadata)

    print(f"\n── Summary ──────────────────────────────────────────────────")
    print(f"  Parent sequences processed : {n_parent_processed}")
    print(f"  Parent sequences skipped   : {len(skipped)}")
    print(f"  Total subsamples written   : {n_subsamples_out}")
    print(f"  Output files written       : {len(output_files)}")
    for f in output_files:
        n = len(list(SeqIO.parse(f, "fasta")))
        print(f"    {Path(f).name} — {n} sequences")
    print(f"  Metadata CSV               : {metadata_path}")

    if skipped:
        print(f"\n  Skipped sequences:")
        for s in skipped:
            print(f"    {s[:70]}")

    return {
        "output_files":       output_files,
        "n_files":            len(output_files),
        "n_parent_seqs":      n_parent_processed,
        "n_subsamples_total": n_subsamples_out,
        "skipped":            skipped,
        "metadata":           metadata,
        "metadata_csv":       metadata_path,
    }

#prepare_subsample_fastas(input_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta",
#                         output_dir = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/genuine",
#                         max_seqs_per_file = 2000,
#                         max_files = None,
#                         )

#prepare_subsample_fastas(input_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_seq/modifed_seqs.fasta",
#                         output_dir = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/artefact",
#                         max_seqs_per_file = 1000,
#                         max_files = None,
#                         )

prepare_subsample_fastas(input_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/all_queries.fasta",
                         output_dir = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/subsample_fastas/query",
                         max_seqs_per_file = 1000,
                         max_files = None,
                         )

