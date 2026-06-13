# _big_fat_file_maker.py

# _big_fat_file_maker.py

import os
import re
import json
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
    Convert a sequence record ID into a Genus_species string suitable
    for fuzzy matching against a species name list.

    Handles three formats:

    Full pipe-delimited header (from run_clustalo full output):
        Sorex_araneus|CYTB|orien:+|accession:AB1234  →  Sorex_araneus

    Short ID with resolved full header (preferred — via id_map):
        Same as above since id_map gives us the full header first.

    Short ID without id_map (fallback only):
        S01_Sorex_AB1234  →  Sorex
        Returns genus only because the second token is an accession
        suffix, not a species epithet. This is ambiguous when multiple
        species of the same genus are present — always provide an
        id_map to avoid this.
    """
    # Full pipe-delimited header — take taxon before first pipe
    if "|" in record_id:
        taxon = record_id.split("|")[0].strip()
        return taxon

    # Strip leading short-ID prefix S01_ or Q01_
    stripped = re.sub(r'^[SQ]\d+_', '', record_id)

    # Determine whether the second underscore-delimited token is a
    # species epithet or an accession suffix.
    # Accession suffixes start with a digit or are all-caps+digits.
    parts = stripped.split("_")
    if len(parts) >= 2:
        second       = parts[1]
        is_accession = (
            second[0].isdigit()
            or (second.isupper() and any(c.isdigit() for c in second))
            or len(second) <= 2   # very short — likely accession fragment
        )
        if is_accession:
            return parts[0]           # genus only
        else:
            return "_".join(parts[:2])  # Genus_species

    return parts[0]


def find_best_key(record_id: str, dictionary: dict):
    """
    Match a normalised sequence ID to the closest key in a dictionary.

    Matching is attempted in this order, stopping at the first hit:

    1. Exact match (case-insensitive)
    2. Containment — target in key, or key in target
    3. Genus-only — first word of target matches start of any key
    4. Fuzzy — difflib closest match (cutoff 0.6)

    Parameters
    ----------
    record_id  : str   normalised ID string from _normalise_id_for_matching
    dictionary : dict  species name dict (keys are species names)

    Returns
    -------
    str or None   matched key, or None if no match found
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

    # 3. Genus-only match — useful when normalisation returned only genus
    genus = target.split("_")[0]
    for i, key in enumerate(keys_lower):
        if key.startswith(genus + "_") or key == genus:
            return keys[i]

    # 4. Fuzzy match as last resort
    matches = difflib.get_close_matches(target, keys_lower,
                                         n=1, cutoff=0.6)
    if matches:
        return keys[keys_lower.index(matches[0])]

    return None


# ─────────────────────────────────────────────────────────────────────────────
# FILE NAME HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _get_name_variants(gene_name: str) -> list:
    """
    Generate lowercase filename variants for a gene name to handle
    differences between CSV column headers and actual file names.

    Examples
    --------
    ND4      →  [nd4]          exact only — avoids matching nd4l
    ND4L     →  [nd4l]
    rRNA_12S →  [rrna_12s, 12s_rrna]
    12S_rRNA →  [12s_rrna, rrna_12s]
    COX1     →  [cox1]
    CYTB     →  [cytb]
    """
    name_lower = gene_name.lower()
    variants   = [name_lower]

    # For two-part underscore names, also try reversing the parts
    # e.g. rrna_12s → 12s_rrna and vice versa
    parts = name_lower.split("_")
    if len(parts) == 2:
        reversed_name = "_".join(reversed(parts))
        if reversed_name != name_lower:
            variants.append(reversed_name)

    # Deduplicate preserving order
    return list(dict.fromkeys(variants))


def _find_nexus_file(folder: Path, gene_name: str,
                      use_short_ids: bool) -> Path:
    """
    Find the aligned Nexus file for a gene using exact name matching.

    Uses exact filename matching rather than glob substring matching
    to prevent ND4 matching ND4L, ND5 matching ND5L etc.

    Searches for:
        {variant}_short_aln.nex   (if use_short_ids=True)
        {variant}_aln.nex         (if use_short_ids=False)

    where variant is each entry from _get_name_variants(gene_name).

    Parameters
    ----------
    folder       : Path   directory to search
    gene_name    : str    gene name from CSV column header
    use_short_ids: bool   whether to look for short-ID Nexus files

    Returns
    -------
    Path or None   path to matched file, or None if not found
    """
    variants   = _get_name_variants(gene_name)
    all_nex    = list(folder.glob("*.nex"))
    all_names  = {f.name.lower(): f for f in all_nex}

    for variant in variants:
        if use_short_ids:
            target = f"{variant}_short_aln.nex"
        else:
            target = f"{variant}_aln.nex"

        if target in all_names:
            return all_names[target]

    return None


# ─────────────────────────────────────────────────────────────────────────────
# PER-GENE TRIMMING
# ─────────────────────────────────────────────────────────────────────────────

def trim_gene_alignment(alignment_path: str,
                         min_coverage: float = 0.5,
                         trim_ends_only: bool = True) -> dict:
    """
    Trim a single gene alignment (Nexus or FASTA) before concatenation.

    Call this per-gene before concatenation. Do not trim the final
    supermatrix — gap columns in the supermatrix may represent species
    missing an entire gene, not poorly aligned regions.

    Parameters
    ----------
    alignment_path : str    path to aligned Nexus (.nex/.nexus) or
                            FASTA (.fasta/.fa) file
    min_coverage   : float  minimum fraction of non-gap characters
                            required to define the shared region
                            endpoints. Default 0.5 is appropriate for
                            small (~12 taxa) datasets with distant
                            outgroups. Use 0.8 for larger datasets
                            with closely related taxa only.
    trim_ends_only : bool   if True (default), only remove leading and
                            lagging columns that fall below min_coverage.
                            Everything between the first and last
                            covered column is kept, preserving internal
                            indels that may be phylogenetically
                            informative. If False, all columns below
                            min_coverage are removed including internal
                            ones — use this only if you want to remove
                            highly gapped internal regions.

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
        covered = coverage >= min_coverage

        if not covered.any():
            print(f"    WARNING: No columns meet min_coverage "
                  f"{min_coverage} — returning untrimmed alignment")
            keep_cols   = list(range(aln_len))
            first, last = 0, aln_len
        else:
            # First and last column meeting the threshold define the
            # shared region — everything between is kept
            first     = int(np.argmax(covered))
            last      = int(aln_len - np.argmax(covered[::-1]))
            keep_cols = list(range(first, last))

        print(f"    End-trim: {aln_len}bp → {len(keep_cols)}bp "
              f"(removed {first} leading + "
              f"{aln_len - last} trailing columns)")

    else:
        # Remove all columns below threshold
        keep_cols = [i for i in range(aln_len)
                     if coverage[i] >= min_coverage]

        n_removed = aln_len - len(keep_cols)
        print(f"    Full-trim: {aln_len}bp → {len(keep_cols)}bp "
              f"({n_removed} columns removed)")

    if len(keep_cols) == 0:
        print(f"    WARNING: All columns removed — "
              f"returning untrimmed alignment")
        keep_cols = list(range(aln_len))
    elif len(keep_cols) < 30:
        print(f"    WARNING: Only {len(keep_cols)}bp retained. "
              f"Check alignment quality and consider lowering "
              f"min_coverage.")

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
                        trim_ends_only: bool = True,
                        trim: bool = True) -> dict:
    """
    Concatenate per-gene aligned Nexus files into a single supermatrix
    Nexus file for multi-gene phylogenetic analysis.

    For each gene, the function:
      1. Finds the aligned Nexus file using exact name matching
         (preventing ND4 from matching ND4L etc.)
      2. Loads the id_map JSON saved by run_clustalo() to resolve
         short IDs (S01_Sorex_AB1234) back to full species names
         (Sorex_araneus). This correctly handles multiple species
         from the same genus.
      3. End-trims the alignment to remove columns from partial
         sequences that don't cover the full gene length.
      4. Matches each sequence to the species list using fuzzy
         matching and appends to the supermatrix.
      5. Gap-fills species that are missing the gene.

    Also writes an IQ-TREE partition file (partitions.txt) so each
    gene can be assigned its own substitution model in a partitioned
    analysis:
        iqtree -s concatenated_aln.nex -p partitions.txt -m TEST -B 1000

    Parameters
    ----------
    raw_fol        : str    directory containing aligned Nexus files
                            and their corresponding id_map JSON files.
                            If use_short_ids=True, looks for
                            {gene}_short_aln.nex files. Otherwise
                            looks for {gene}_aln.nex files.
    csv_path       : str    CSV with Species in first column and gene
                            names as remaining column headers. Gene
                            names must match Nexus filenames exactly
                            (case-insensitive, word-order variants
                            tried automatically e.g. rRNA_12S →
                            12S_rRNA).
    output         : str    directory to write output files
    datatype       : str    Nexus datatype string e.g. "dna"
    use_short_ids  : bool   if True (default), read short-ID Nexus
                            files (_short_aln.nex) and resolve species
                            names via id_map JSON. Set False to use
                            full-header Nexus files (_aln.nex).
    min_coverage   : float  per-gene end-trimming coverage threshold.
                            Default 0.5 suits small (~12 taxa) datasets
                            with distant outgroups (e.g. Xenopus).
                            Use 0.8 for large datasets.
    trim_ends_only : bool   if True (default), only trim leading and
                            lagging columns — internal indels are kept.

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

    # Initialise supermatrix — empty string per species
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

        # ── Find the Nexus file for this gene ─────────────────────────────────
        file_path = _find_nexus_file(folder, name, use_short_ids)

        if file_path is None:
            print(f"  [{name}] No Nexus file found — skipping")
            variants = _get_name_variants(name)
            suffix   = "_short_aln.nex" if use_short_ids else "_aln.nex"
            print(f"    Looked for: "
                  f"{[v + suffix for v in variants]}")
            continue

        print(f"  [{name}] Reading {file_path.name}")

        # ── Load id_map to resolve short IDs to full species names ────────────
        # id_map maps S01_Sorex_AB1234 → Sorex_araneus|CYTB|...
        # Without this, multiple Sorex species would all normalise to
        # "Sorex" and arbitrarily match the first entry in the dict.
        id_map      = {}
        id_map_path = str(file_path).replace("_short_aln.nex",
                                              "_id_map.json")
        if use_short_ids:
            if os.path.exists(id_map_path):
                with open(id_map_path, "r") as f:
                    id_map = json.load(f)
                print(f"    id_map loaded: {len(id_map)} entries")
            else:
                print(f"    WARNING: no id_map found at {id_map_path}")
                print(f"    Species with multiple members of the same "
                      f"genus may be mismatched. Re-run run_clustalo() "
                      f"to regenerate id_map files.")

        # ── Trim alignment ────────────────────────────────────────────────────
        # ── Trim alignment ────────────────────────────────────────────────────
        if trim and min_coverage is not None:
            trim_result = trim_gene_alignment(
                alignment_path = str(file_path),
                min_coverage   = min_coverage,
                trim_ends_only = trim_ends_only,
            )
            raw_seqs    = trim_result["sequences"]
            gene_length = trim_result["trimmed_length"]
        else:
            if not trim:
                print(f"    Trimming disabled — using raw alignment")
            raw_seqs = {
                r.id: str(r.seq)
                for r in SeqIO.parse(file_path, "nexus")
            }
            gene_length = (len(next(iter(raw_seqs.values())))
                           if raw_seqs else 0)

        if gene_length == 0:
            print(f"  [{name}] Zero length after trimming — skipping")
            continue

        # ── Match each sequence to the species list ───────────────────────────
        gene_seqs = {}

        for record_id, seq in raw_seqs.items():

            # Skip duplicate/copy entries
            if any(record_id.endswith(s)
                   for s in (".copy", "_copy", ".2", "_dup")):
                continue

            # Resolve full header via id_map if available.
            # This recovers Sorex_araneus from S01_Sorex_AB1234
            # and correctly distinguishes between Sorex species.
            if use_short_ids and record_id in id_map:
                full_header   = id_map[record_id]
                normalised_id = _normalise_id_for_matching(full_header)
            else:
                # No id_map available — parse the ID directly.
                # Genus-only for short IDs, full name for pipe headers.
                normalised_id = _normalise_id_for_matching(record_id)

            match = find_best_key(normalised_id, my_dict)

            if match:
                if match not in gene_seqs:
                    gene_seqs[match] = seq
                else:
                    print(f"    Duplicate match for {match} "
                          f"from {record_id} — keeping first")
            else:
                print(f"    No match for: {record_id} "
                      f"(normalised: {normalised_id})")

        if not gene_seqs:
            print(f"  [{name}] No sequences matched — skipping")
            continue

        # ── Record partition for IQ-TREE ──────────────────────────────────────
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
            f"Unique lengths found: {set(lengths)}. "
            f"This indicates a bug in the gap-filling logic."
        )

    nchar = lengths[0]
    ntax  = len(my_dict)

    if nchar == 0:
        raise ValueError(
            "Supermatrix is empty — no genes were successfully "
            "concatenated. Check that gene names in the CSV match "
            "Nexus filenames in the alignment directory."
        )

    # ── Write Nexus supermatrix ───────────────────────────────────────────────
    output_file  = os.path.join(output, "concatenated_aln.nex")
    max_name_len = max(len(n) for n in my_dict.keys())

    print(f"\nWriting supermatrix Nexus: "
          f"{ntax} taxa, {nchar}bp total, "
          f"{len(partition_info)} genes")

    with open(output_file, "w") as nex:
        nex.write("#NEXUS\n")
        nex.write("Begin DATA;\n")
        nex.write(f"    Dimensions ntax={ntax} nchar={nchar};\n")
        nex.write(f"    Format datatype={datatype} "
                  f"missing=? gap=-;\n")
        nex.write("    Matrix\n")
        for species, seq in my_dict.items():
            padded = species.ljust(max_name_len)
            nex.write(f"{padded} {seq}\n")
        nex.write("    ;\n")
        nex.write("End;\n")

    # ── Write IQ-TREE partition file ──────────────────────────────────────────
    partition_file = os.path.join(output, "partitions.txt")
    with open(partition_file, "w") as pf:
        for p in partition_info:
            pf.write(
                f"DNA, {p['gene']} = {p['start']}-{p['end']}\n"
            )

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
    trim = False,
)