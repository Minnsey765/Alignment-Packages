# _big_fat_file_maker.py

import os
import difflib
import pandas as pd

from pathlib import Path
from Bio import SeqIO, AlignIO
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# ID NORMALISATION AND MATCHING
# ─────────────────────────────────────────────────────────────────────────────

def _normalise_id_for_matching(record_id: str) -> str:
    """
    Convert a sequence record ID from either short-ID or full
    pipe-delimited format into a plain Genus_species string
    suitable for fuzzy matching against species names.

    Short ID  : S01_Sorex_AB1234        → Sorex
    Full ID   : Sorex_araneus|CYTB|...  → Sorex_araneus
    Plain     : Sorex_araneus_AB1234    → Sorex_araneus
    """
    import re

    # Full pipe-delimited — take the taxon part before the first pipe
    if "|" in record_id:
        return record_id.split("|")[0].strip()

    # Short ID format: S01_Genus_AccSuffix or Q01_Genus_AccSuffix
    # Strip the leading S01_ / Q01_ prefix
    stripped = re.sub(r'^[SQ]\d+_', '', record_id)

    # Take first two underscore parts as Genus_species
    parts = stripped.split("_")
    if len(parts) >= 2:
        return "_".join(parts[:2])
    return parts[0]


def find_best_key(record_id: str, dictionary: dict):
    """
    Match a sequence record ID to the closest key in a dictionary.

    Tries, in order:
      1. Exact match (case-insensitive)
      2. Target contained in key, or key contained in target
      3. Genus-only match (first word of target matches start of key)
      4. Fuzzy match (difflib, cutoff 0.6)

    Returns the matched key string, or None if no match found.
    """
    target     = record_id.replace(" ", "_").lower()
    keys       = list(dictionary.keys())
    keys_lower = [k.lower() for k in keys]

    # 1. Exact match
    for i, key in enumerate(keys_lower):
        if target == key:
            return keys[i]

    # 2. Containment match
    for i, key in enumerate(keys_lower):
        if target in key or key in target:
            return keys[i]

    # 3. Genus-only match
    genus = target.split("_")[0]
    for i, key in enumerate(keys_lower):
        if key.startswith(genus + "_") or key == genus:
            return keys[i]

    # 4. Fuzzy match
    matches = difflib.get_close_matches(target, keys_lower,
                                         n=1, cutoff=0.6)
    if matches:
        return keys[keys_lower.index(matches[0])]

    return None


# ─────────────────────────────────────────────────────────────────────────────
# PER-GENE TRIMMING
# ─────────────────────────────────────────────────────────────────────────────

def trim_gene_alignment(alignment_path: str,
                         min_coverage: float = 0.5,
                         trim_ends_only: bool = True) -> dict:
    """
    Trim a single gene alignment (Nexus or FASTA) before concatenation.

    Should be called per-gene, not on the final supermatrix. Trimming
    the supermatrix would incorrectly penalise columns that are all
    gaps for a species that lacks a gene entirely, even if every
    species that has the gene is fully covered at that position.

    Parameters
    ----------
    alignment_path : str    path to aligned Nexus (.nex/.nexus) or
                            FASTA (.fasta/.fa) file
    min_coverage   : float  minimum fraction of non-gap characters
                            required to keep a column. Default 0.5
                            is appropriate for small (~12 taxa) diverse
                            datasets with distant outgroups. Use 0.8
                            for large datasets with closely related taxa.
    trim_ends_only : bool   if True (default), only remove leading and
                            lagging columns below min_coverage. Internal
                            gap columns (from genuine indels) are kept.
                            If False, all low-coverage columns are
                            removed including internal ones.

    Returns
    -------
    dict with:
        sequences       : dict  record_id -> trimmed sequence string
        original_length : int   alignment length before trimming
        trimmed_length  : int   alignment length after trimming
        columns_removed : int   number of columns removed
    """
    ext = Path(alignment_path).suffix.lower()
    fmt = "nexus" if ext in (".nex", ".nexus") else "fasta"

    alignment = AlignIO.read(alignment_path, fmt)
    n_seq     = len(alignment)
    aln_len   = alignment.get_alignment_length()

    # Per-column coverage: fraction of sequences with a real base
    coverage = np.array([
        1.0 - (alignment[:, i].count("-") / n_seq)
        for i in range(aln_len)
    ])

    if trim_ends_only:
        # Find the first and last column meeting min_coverage.
        # Keep everything between them — internal gaps are preserved.
        covered = coverage >= min_coverage

        if not covered.any():
            print(f"    WARNING: No columns meet min_coverage "
                  f"{min_coverage} — returning untrimmed alignment")
            keep_cols = list(range(aln_len))
            first, last = 0, aln_len
        else:
            first     = int(np.argmax(covered))
            last      = int(aln_len - np.argmax(covered[::-1]))
            keep_cols = list(range(first, last))

        print(f"    End-trim: {aln_len}bp → {len(keep_cols)}bp "
              f"(removed {first} leading + "
              f"{aln_len - last} trailing columns)")

    else:
        # Remove all columns below threshold including internal ones
        keep_cols = [i for i in range(aln_len)
                     if coverage[i] >= min_coverage]
        first, last = keep_cols[0] if keep_cols else 0, aln_len

        print(f"    Full-trim: {aln_len}bp → {len(keep_cols)}bp "
              f"({aln_len - len(keep_cols)} columns removed)")

    if len(keep_cols) < 50:
        print(f"    WARNING: Only {len(keep_cols)}bp retained. "
              f"Consider lowering min_coverage.")

    sequences = {}
    for record in alignment:
        trimmed_seq = "".join(str(record.seq)[i] for i in keep_cols)
        sequences[record.id] = trimmed_seq

    return {
        "sequences":       sequences,
        "original_length": aln_len,
        "trimmed_length":  len(keep_cols),
        "columns_removed": aln_len - len(keep_cols),
    }


# ─────────────────────────────────────────────────────────────────────────────
# MAIN FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

def big_fat_file_maker(raw_fol: str,
                        csv_path: str,
                        output: str,
                        datatype: str,
                        use_short_ids: bool = True,
                        min_coverage: float = 0.5,
                        trim_ends_only: bool = True) -> dict:
    """
    Concatenate per-gene aligned Nexus files into a single supermatrix
    Nexus file for multi-gene phylogenetic analysis.

    Missing genes for a species are gap-filled so every species has
    the same total sequence length. IQ-TREE treats gaps as missing
    data, so gap-filled regions contribute no information rather than
    misleading signal.

    Each gene alignment is end-trimmed before concatenation to remove
    columns caused by partial sequences that do not cover the full gene
    length. Internal gap columns (genuine indels) are preserved.

    Also writes an IQ-TREE partition file so each gene can be given
    its own substitution model in a partitioned analysis, which
    substantially improves tree accuracy for concatenated multi-gene
    data.

    Parameters
    ----------
    raw_fol        : str    directory containing aligned Nexus files.
                            If use_short_ids=True, looks for
                            *_short_aln.nex files. Otherwise looks for
                            *_aln.nex (excluding *_short_aln.nex).
    csv_path       : str    CSV with Species in first column and gene
                            names as remaining column headers. Species
                            names must match sequence headers closely
                            enough for fuzzy matching to work.
    output         : str    directory to write output files
    datatype       : str    Nexus datatype string, e.g. "dna"
    use_short_ids  : bool   if True (default), match against short-ID
                            Nexus files (_short_aln.nex). Set to False
                            to use full-header files (_aln.nex).
    min_coverage   : float  per-gene trimming threshold. Default 0.5
                            is appropriate for small datasets (~12 taxa)
                            with distant outgroups such as Xenopus.
                            Increase to 0.8 for large datasets.
    trim_ends_only : bool   if True (default), only trim leading and
                            lagging columns. Preserves internal indels.

    Returns
    -------
    dict with:
        nexus_file      : str    path to concatenated_aln.nex
        partition_file  : str    path to partitions.txt
        partition_info  : list   dicts with gene, start, end, length
        total_length    : int    total alignment length in bp
        n_genes         : int    number of genes successfully included
        n_taxa          : int    number of taxa in the supermatrix
    """
    csv           = pd.read_csv(csv_path)
    species_names = csv["Species"].tolist()
    gene_names    = csv.columns[1:].tolist()

    # Initialise supermatrix with empty strings for each species
    my_dict = {item: "" for item in species_names}

    folder         = Path(raw_fol)
    partition_info = []
    position       = 1

    print(f"\nBuilding supermatrix from {len(gene_names)} genes, "
          f"{len(species_names)} taxa")
    print(f"Source directory : {raw_fol}")
    print(f"Trimming         : min_coverage={min_coverage}, "
          f"trim_ends_only={trim_ends_only}\n")

    for name in gene_names:

        # ── Find the aligned Nexus file for this gene ─────────────────────────
        if use_short_ids:
            matching_files = list(folder.glob(f"*{name}*_short_aln.nex"))
        else:
            all_matches    = list(folder.glob(f"*{name}*_aln.nex"))
            matching_files = [f for f in all_matches
                              if "_short_" not in f.name]

        if not matching_files:
            print(f"  [{name}] No Nexus file found — skipping")
            continue

        if len(matching_files) > 1:
            print(f"  [{name}] Multiple files matched — "
                  f"using {matching_files[0].name}")

        file_path = matching_files[0]
        print(f"  [{name}] Reading {file_path.name}")

        # ── Trim before reading sequences ─────────────────────────────────────
        if min_coverage is not None:
            trim_result = trim_gene_alignment(
                alignment_path = str(file_path),
                min_coverage   = min_coverage,
                trim_ends_only = trim_ends_only,
            )
            raw_seqs    = trim_result["sequences"]
            gene_length = trim_result["trimmed_length"]
        else:
            # No trimming — read sequences directly
            raw_seqs = {
                r.id: str(r.seq)
                for r in SeqIO.parse(file_path, "nexus")
            }
            gene_length = (len(next(iter(raw_seqs.values())))
                           if raw_seqs else 0)

        if gene_length == 0:
            print(f"  [{name}] Empty alignment after trimming — skipping")
            continue

        # ── Match sequences to species list ───────────────────────────────────
        gene_seqs = {}
        for record_id, seq in raw_seqs.items():

            # Skip duplicate/copy entries
            if any(record_id.endswith(s)
                   for s in (".copy", "_copy", ".2", "_dup")):
                continue

            normalised_id = _normalise_id_for_matching(record_id)
            match         = find_best_key(normalised_id, my_dict)

            if match:
                # If a species already has a sequence for this gene
                # (e.g. from a duplicate alignment entry), keep the
                # first one
                if match not in gene_seqs:
                    gene_seqs[match] = seq
            else:
                print(f"    No match for: {record_id} "
                      f"(normalised: {normalised_id})")

        if not gene_seqs:
            print(f"  [{name}] No sequences matched — skipping")
            continue

        # ── Record partition info for IQ-TREE ─────────────────────────────────
        partition_info.append({
            "gene":   name,
            "start":  position,
            "end":    position + gene_length - 1,
            "length": gene_length,
        })
        position += gene_length

        # ── Append to supermatrix — gap-fill missing species ──────────────────
        for species in species_names:
            if species in gene_seqs:
                my_dict[species] += gene_seqs[species]
            else:
                my_dict[species] += "-" * gene_length

        # ── Coverage report ───────────────────────────────────────────────────
        n_present = sum(1 for sp in species_names if sp in gene_seqs)
        n_missing = len(species_names) - n_present
        print(f"  [{name}] {n_present}/{len(species_names)} species "
              f"found, {n_missing} gap-filled, "
              f"{gene_length}bp after trimming")

        if n_missing > 0:
            missing = [sp for sp in species_names
                       if sp not in gene_seqs]
            print(f"    Gap-filled: {', '.join(missing)}")

    # ── Validate supermatrix ──────────────────────────────────────────────────
    os.makedirs(output, exist_ok=True)

    lengths = [len(seq) for seq in my_dict.values()]
    if len(set(lengths)) != 1:
        raise ValueError(
            f"Sequences have unequal lengths after concatenation. "
            f"Unique lengths: {set(lengths)}"
        )

    nchar = lengths[0]
    ntax  = len(my_dict)

    if nchar == 0:
        raise ValueError(
            "Supermatrix is empty — no genes were successfully "
            "concatenated. Check that gene names in the CSV match "
            "file names in the alignment directory."
        )

    # ── Write Nexus supermatrix ───────────────────────────────────────────────
    output_file  = os.path.join(output, "concatenated_aln.nex")
    max_name_len = max(len(name) for name in my_dict.keys())

    print(f"\nWriting supermatrix Nexus: {ntax} taxa, {nchar}bp total, "
          f"{len(partition_info)} genes")

    with open(output_file, "w") as nex:
        nex.write("#NEXUS\n")
        nex.write("Begin DATA;\n")
        nex.write(f"    Dimensions ntax={ntax} nchar={nchar};\n")
        nex.write(f"    Format datatype={datatype} missing=? gap=-;\n")
        nex.write("    Matrix\n")
        for species, seq in my_dict.items():
            padded = species.ljust(max_name_len)
            nex.write(f"{padded} {seq}\n")
        nex.write("    ;\n")
        nex.write("End;\n")

    # ── Write IQ-TREE partition file ──────────────────────────────────────────
    # Allows IQ-TREE to fit a separate substitution model per gene
    # using: iqtree -s concatenated_aln.nex -p partitions.txt -m TEST
    partition_file = os.path.join(output, "partitions.txt")
    with open(partition_file, "w") as pf:
        for p in partition_info:
            pf.write(f"DNA, {p['gene']} = {p['start']}-{p['end']}\n")

    print(f"Partition file written  : {partition_file}")
    print(f"Supermatrix written     : {output_file}")
    print(f"Total alignment length  : {nchar}bp")
    print(f"Genes included          : {len(partition_info)}")
    print("Finished!")

    return {
        "nexus_file":     output_file,
        "partition_file": partition_file,
        "partition_info": partition_info,
        "total_length":   nchar,
        "n_genes":        len(partition_info),
        "n_taxa":         ntax,
    }

big_fat_file_maker(
    raw_fol      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_nexus/Framework_alignment/genes",
    csv_path     = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/MSA_framework.csv",
    output       = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_nexus/Framework_alignment/concatenated",
    datatype     = "dna",
    use_short_ids = True,   # use _short_aln.nex files
)