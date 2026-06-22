# big_fat_file_maker_v2.py
#
# Concatenate per-gene aligned NEXUS files into a single supermatrix.
#
# Accession extraction handles the header formats produced by restore_headers.py:
#   FASTA : Sorex_araneus||CYTB||NC_027963   (split on ||, last field)
#   NEXUS : Sorex_araneus__CYTB__NC_027963   (split on __, last field)
#   Legacy: Sorex_araneus_CYTB_NC_027963     (regex fallback)
#
# Functions:
#   big_fat_file_maker()         — main concatenation function
#   trim_gene_alignment()        — per-gene trimming (use if not manually trimmed)
#   write_partition_metadata()   — write gene positions + reading frames to JSON
#   generate_mrbayes_charsets()  — generate MrBayes charset block from metadata
#   _codon_positions()           — calculate codon position ranges

import os
import re
import json
import difflib
import pandas as pd

from pathlib import Path
from Bio import AlignIO
import numpy as np


# ===========================================================================
# Accession extraction
# ===========================================================================

def _extract_accession_from_header(record_id: str) -> str:
    """
    Extract the accession from a sequence record header.

    Handles the following formats (produced by restore_headers.py):

        FASTA medium : Sorex_araneus||CYTB||NC_027963
                       → split on "||", take last field → NC_027963

        NEXUS medium : Sorex_araneus__CYTB__NC_027963
                       → split on "__", take last field → NC_027963

        Full/legacy  : Sorex_araneus|CYTB|orien:+|accession:NC_027963
                       → split on "|", find "accession:" field → NC_027963

        Legacy underscore: Sorex_araneus_CYTB_NC_027963
                       → regex fallback → NC_027963

    Returns empty string if no real accession found (e.g. BraceEtAl).
    """
    # Double-pipe format (FASTA medium)
    if "||" in record_id:
        last = record_id.split("||")[-1].strip()
        if last and any(c.isdigit() for c in last):
            return last

    # Double-underscore format (NEXUS medium)
    if "__" in record_id:
        last = record_id.split("__")[-1].strip()
        if last and any(c.isdigit() for c in last):
            return last

    # Full format with "accession:" prefix (original run_clustalo headers)
    if "accession:" in record_id:
        for part in record_id.split("|"):
            if part.startswith("accession:"):
                val = part.split(":", 1)[1].strip()
                if any(c.isdigit() for c in val):
                    return val

    # Regex fallback: GenBank accession pattern anywhere in the string
    acc_match = re.search(
        r'(?<![A-Z0-9])([A-Z]{1,2}_?\d{5,9}(?:\.\d+)?)(?![A-Z0-9])',
        record_id
    )
    return acc_match.group(1) if acc_match else ""


def _normalise_accession(acc: str) -> str:
    """Strip version suffix (.1 etc.), whitespace, uppercase."""
    if not acc:
        return ""
    return re.sub(r'\.\d+$', '', acc.strip().upper())


def _build_accession_index(raw_seqs: dict) -> dict:
    """Build {normalised_accession: record_id} for fast lookup."""
    index = {}
    for record_id in raw_seqs:
        acc = _normalise_accession(_extract_accession_from_header(record_id))
        if acc:
            index[acc] = record_id
    return index


def _find_by_accession(target_acc: str, acc_index: dict, raw_seqs: dict) -> tuple:
    """Return (record_id, sequence) for the accession, or (None, None)."""
    norm = _normalise_accession(target_acc)
    if not norm:
        return None, None
    record_id = acc_index.get(norm)
    return (record_id, raw_seqs[record_id]) if record_id else (None, None)


# ===========================================================================
# Fuzzy name matching (fallback only)
# ===========================================================================

def _normalise_id_for_matching(record_id: str) -> str:
    """Extract Genus_species from a header for fuzzy matching."""
    # Double-pipe or double-underscore medium format — species is first field
    if "||" in record_id:
        return record_id.split("||")[0].strip()
    if "__" in record_id:
        return record_id.split("__")[0].strip()
    if "|" in record_id:
        return record_id.split("|")[0].strip()
    stripped = re.sub(r'^[SQ]\d+_', '', record_id)
    parts    = stripped.split("_")
    if len(parts) >= 2:
        second = parts[1]
        is_acc = (
            second[0].isdigit()
            or (second.isupper() and any(c.isdigit() for c in second))
            or len(second) <= 2
        )
        return parts[0] if is_acc else "_".join(parts[:2])
    return parts[0]


def find_best_key(record_id: str, dictionary: dict):
    """Match a normalised ID to the closest key in a dictionary."""
    target     = record_id.replace(" ", "_").lower()
    keys       = list(dictionary.keys())
    keys_lower = [k.lower() for k in keys]
    for i, key in enumerate(keys_lower):
        if target == key:
            return keys[i]
    for i, key in enumerate(keys_lower):
        if target in key or key in target:
            return keys[i]
    genus = target.split("_")[0]
    for i, key in enumerate(keys_lower):
        if key.startswith(genus + "_") or key == genus:
            return keys[i]
    matches = difflib.get_close_matches(target, keys_lower, n=1, cutoff=0.6)
    if matches:
        return keys[keys_lower.index(matches[0])]
    return None


# ===========================================================================
# File helpers
# ===========================================================================

def _get_name_variants(gene_name: str) -> list:
    """
    Lowercase filename variants for a gene name.
    Handles reversed two-part names (12S_rRNA ↔ rRNA_12S).
    Exact matching — avoids ND4 matching ND4L.
    """
    name_lower = gene_name.lower()
    variants   = [name_lower]
    parts      = name_lower.split("_")
    if len(parts) == 2:
        rev = "_".join(reversed(parts))
        if rev != name_lower:
            variants.append(rev)
    return list(dict.fromkeys(variants))


def _find_nexus_file(folder: Path, gene_name: str, suffix: str = "_aln.nex") -> Path:
    """Find the aligned NEXUS file for a gene using exact name matching."""
    variants  = _get_name_variants(gene_name)
    all_nex   = list(folder.glob("*.nex"))
    all_names = {f.name.lower(): f for f in all_nex}
    for variant in variants:
        if f"{variant}{suffix}" in all_names:
            return all_names[f"{variant}{suffix}"]
    return None


# ===========================================================================
# Special source detection
# ===========================================================================

_NON_ACCESSION_VALUES = {
    "brace et al.", "brace et al", "brace et al. (2016)",
    "yuan et al.", "yuan et al", "in-house", "unpublished",
}

def _is_special_source(cell_value: str) -> bool:
    return cell_value.strip().lower() in _NON_ACCESSION_VALUES


# ===========================================================================
# Per-gene trimming
# ===========================================================================

def trim_gene_alignment(alignment_path: str,
                         min_coverage: float = 0.5,
                         trim_ends_only: bool = True) -> dict:
    """
    Trim a single gene alignment before concatenation.
    Only trim individual gene files — never the final supermatrix.

    Parameters
    ----------
    alignment_path : str    path to .nex/.nexus or .fasta alignment
    min_coverage   : float  minimum non-gap fraction to define endpoints
    trim_ends_only : bool   if True, only remove leading/trailing columns

    Returns
    -------
    dict with sequences, original_length, trimmed_length, columns_removed
    """
    ext       = Path(alignment_path).suffix.lower()
    fmt       = "nexus" if ext in (".nex", ".nexus") else "fasta"
    alignment = AlignIO.read(alignment_path, fmt)
    n_seq     = len(alignment)
    aln_len   = alignment.get_alignment_length()

    coverage = np.array([
        1.0 - (alignment[:, i].count("-") / n_seq)
        for i in range(aln_len)
    ])

    if trim_ends_only:
        covered = coverage >= min_coverage
        if not covered.any():
            print(f"    WARNING: no columns meet min_coverage {min_coverage}")
            keep_cols   = list(range(aln_len))
            first, last = 0, aln_len
        else:
            first     = int(np.argmax(covered))
            last      = int(aln_len - np.argmax(covered[::-1]))
            keep_cols = list(range(first, last))
        print(f"    End-trim: {aln_len}bp → {len(keep_cols)}bp "
              f"(removed {first} leading + {aln_len - last} trailing)")
    else:
        keep_cols = [i for i in range(aln_len) if coverage[i] >= min_coverage]
        print(f"    Full-trim: {aln_len}bp → {len(keep_cols)}bp")

    if not keep_cols:
        keep_cols = list(range(aln_len))

    return {
        "sequences": {
            record.id: "".join(str(record.seq)[i] for i in keep_cols)
            for record in alignment
        },
        "original_length": aln_len,
        "trimmed_length":  len(keep_cols),
        "columns_removed": aln_len - len(keep_cols),
    }


# ===========================================================================
# Main concatenation function
# ===========================================================================

def big_fat_file_maker(raw_fol: str,
                        csv_path: str,
                        output: str,
                        datatype: str,
                        file_suffix: str = "_aln.nex",
                        min_coverage: float = 0.5,
                        trim_ends_only: bool = True,
                        trim: bool = False) -> dict:
    """
    Concatenate per-gene aligned NEXUS files into a single supermatrix
    NEXUS file for multi-gene phylogenetic analysis.

    Matching strategy (accession-first)
    ------------------------------------
    For each species/gene cell in the CSV:
    1. Accession match — the exact accession from the CSV cell is looked
       up in the alignment file headers. Handles medium-format headers
       (Sorex_araneus__CYTB__NC_027963) and legacy formats.
    2. Citation fallback — "Brace et al." etc. have no accession; falls
       back to genus+species name matching.
    3. Empty cell — gap-filled with '-' characters.

    Parameters
    ----------
    raw_fol        : str    directory with [gene]_aln.nex files
                            (output of loop_restore_headers)
    csv_path       : str    CSV with Species in column 0, gene names as
                            remaining headers, accessions as cell values
    output         : str    directory to write supermatrix and partition files
    datatype       : str    NEXUS datatype string e.g. "dna"
    file_suffix    : str    alignment file suffix (default "_aln.nex")
    min_coverage   : float  trimming threshold (only used if trim=True)
    trim_ends_only : bool   if True, only trim leading/trailing columns
    trim           : bool   if False (default), skip auto-trimming

    Returns
    -------
    dict with nexus_file, partition_file, partition_info,
             total_length, n_genes, n_taxa
    """
    csv            = pd.read_csv(csv_path)
    csv["Species"] = csv["Species"].str.strip()
    species_names  = csv["Species"].tolist()
    gene_names     = csv.columns[1:].tolist()

    my_dict        = {sp: "" for sp in species_names}
    folder         = Path(raw_fol)
    partition_info = []
    position       = 1

    print(f"\nBuilding supermatrix from {len(gene_names)} genes, "
          f"{len(species_names)} taxa")
    print(f"Source directory : {raw_fol}")
    print(f"Trimming         : {'disabled' if not trim else f'min_coverage={min_coverage}'}\n")

    for gene in gene_names:

        file_path = _find_nexus_file(folder, gene, suffix=file_suffix)
        if file_path is None:
            print(f"  [{gene}] No file found — skipping")
            print(f"    Looked for: "
                  f"{[v + file_suffix for v in _get_name_variants(gene)]}")
            continue

        print(f"  [{gene}] Reading {file_path.name}")

        # Load sequences
        if trim:
            result      = trim_gene_alignment(str(file_path), min_coverage, trim_ends_only)
            raw_seqs    = result["sequences"]
            gene_length = result["trimmed_length"]
        else:
            fmt      = "nexus" if file_path.suffix.lower() in (".nex", ".nexus") else "fasta"
            raw_seqs = {r.id: str(r.seq)
                        for r in AlignIO.read(str(file_path), fmt)}
            gene_length = len(next(iter(raw_seqs.values()))) if raw_seqs else 0

        if gene_length == 0:
            print(f"  [{gene}] Zero length — skipping")
            continue

        acc_index    = _build_accession_index(raw_seqs)
        gene_seqs    = {}
        need_fallback = []

        for species in species_names:
            row      = csv[csv["Species"] == species]
            cell_val = str(row[gene].values[0]).strip() if len(row) > 0 else ""

            if cell_val in ("", "nan", "NaN"):
                continue

            if _is_special_source(cell_val):
                need_fallback.append((species, cell_val))
                continue

            record_id, seq = _find_by_accession(cell_val, acc_index, raw_seqs)
            if seq is not None:
                gene_seqs[species] = seq
            else:
                print(f"    [warn] {species}: accession {cell_val} not found "
                      f"— trying name match")
                need_fallback.append((species, cell_val))

        # Name-based fallback
        for species, cell_val in need_fallback:
            if species in gene_seqs:
                continue

            matched_id = None
            sp_lower   = species.lower()

            for record_id in raw_seqs:
                rid_lower = record_id.lower()
                genus, *rest = sp_lower.split("_")
                epithet = rest[0] if rest else ""
                if genus in rid_lower and (not epithet or epithet in rid_lower):
                    matched_id = record_id
                    break

            if not matched_id:
                norm = _normalise_id_for_matching(species)
                best = find_best_key(norm, {rid: True for rid in raw_seqs})
                if best:
                    matched_id = best

            if matched_id:
                gene_seqs[species] = raw_seqs[matched_id]
                src = "citation" if _is_special_source(cell_val) else "name fallback"
                print(f"    [fallback/{src}] {species} → '{matched_id}'")
            else:
                print(f"    [no match] {species} (cell='{cell_val}') — gap-filled")

        if not gene_seqs:
            print(f"  [{gene}] No sequences matched — skipping")
            continue

        partition_info.append({
            "gene":   gene,
            "start":  position,
            "end":    position + gene_length - 1,
            "length": gene_length,
        })
        position += gene_length

        for species in species_names:
            my_dict[species] += gene_seqs.get(species, "-" * gene_length)

        n_present = len(gene_seqs)
        n_missing = len(species_names) - n_present
        print(f"  [{gene}] {n_present}/{len(species_names)} present, "
              f"{n_missing} gap-filled, {gene_length}bp")
        if n_missing > 0:
            missing = [sp for sp in species_names if sp not in gene_seqs]
            print(f"    Gap-filled: {', '.join(missing)}")

    # Validate
    os.makedirs(output, exist_ok=True)
    lengths = [len(seq) for seq in my_dict.values()]
    if len(set(lengths)) != 1:
        raise ValueError(f"Unequal sequence lengths: {set(lengths)}")

    nchar = lengths[0]
    ntax  = len(my_dict)

    if nchar == 0:
        raise ValueError("Supermatrix is empty — no genes concatenated.")

    # Write NEXUS supermatrix
    output_file  = os.path.join(output, "concatenated_aln.nex")
    max_name_len = max(len(n) for n in my_dict)

    print(f"\nWriting supermatrix: {ntax} taxa, {nchar}bp, {len(partition_info)} genes")

    with open(output_file, "w") as nex:
        nex.write("#NEXUS\n")
        nex.write("Begin DATA;\n")
        nex.write(f"    Dimensions ntax={ntax} nchar={nchar};\n")
        nex.write(f"    Format datatype={datatype} missing=? gap=-;\n")
        nex.write("    Matrix\n")
        for species, seq in my_dict.items():
            nex.write(f"{species.ljust(max_name_len)} {seq}\n")
        nex.write("    ;\n")
        nex.write("End;\n")

    # Write IQ-TREE / MrBayes partition file
    partition_file = os.path.join(output, "partitions.txt")
    with open(partition_file, "w") as pf:
        for p in partition_info:
            pf.write(f"DNA, {p['gene']} = {p['start']}-{p['end']}\n")

    print(f"Supermatrix written  : {output_file}")
    print(f"Partition file       : {partition_file}")
    print(f"Total length         : {nchar}bp")
    print(f"Genes included       : {len(partition_info)}")
    print("Finished!")

    return {
        "nexus_file":     output_file,
        "partition_file": partition_file,
        "partition_info": partition_info,
        "total_length":   nchar,
        "n_genes":        len(partition_info),
        "n_taxa":         ntax,
    }


# ===========================================================================
# Partition metadata and MrBayes charset generation
# ===========================================================================

def _codon_positions(start: int, end: int, frame: int) -> dict:
    """
    Calculate codon position column ranges within a gene block.

    Parameters
    ----------
    start : int   1-based start column in the supermatrix
    end   : int   1-based end column in the supermatrix
    frame : int   reading frame (1, 2, or 3), or None

    Returns
    -------
    dict {1: "start-end\\3", 2: ..., 3: ...} or None if frame unknown
    """
    if frame is None:
        return None
    offsets = {
        1: {1: 0, 2: 1, 3: 2},
        2: {1: 2, 2: 0, 3: 1},
        3: {1: 1, 2: 2, 3: 0},
    }
    if frame not in offsets:
        return None
    return {
        codon_pos: f"{start + offset}-{end}\\3"
        for codon_pos, offset in offsets[frame].items()
    }


def write_partition_metadata(partition_info: list,
                              output_dir: str,
                              reading_frames: dict = None) -> str:
    """
    Write a JSON file recording gene positions in the supermatrix,
    genetic codes, and codon position ranges for MrBayes.

    Parameters
    ----------
    partition_info  : list   from big_fat_file_maker() return dict
    output_dir      : str    where to write partition_metadata.json
    reading_frames  : dict   {gene_name: int} e.g. {"CYTB": 2, "COX1": 1}
                             Only for protein-coding genes; omit rRNA.

    Returns
    -------
    str   path to partition_metadata.json
    """
    MT_GENES   = {"ATP6","ATP8","COX1","COX2","COX3",
                  "CYTB","ND1","ND2","ND3","ND4","ND4L","ND5","ND6"}
    NU_GENES   = {"APOB","BDNF","BRCA1","GHR","RAG1","RAG2"}
    RRNA_GENES = {"12S_rRNA","16S_rRNA"}

    reading_frames = reading_frames or {}
    metadata = {
        "total_length": sum(p["length"] for p in partition_info),
        "n_genes":      len(partition_info),
        "genes":        [],
    }

    for p in partition_info:
        gene  = p["gene"]
        gtype = ("mitochondrial_protein_coding" if gene in MT_GENES else
                 "nuclear_protein_coding"        if gene in NU_GENES else
                 "rRNA"                          if gene in RRNA_GENES else
                 "unknown")
        code  = ("vertebrate_mitochondrial" if gene in MT_GENES else
                 "standard"                 if gene in NU_GENES else
                 None)
        frame = reading_frames.get(gene)
        metadata["genes"].append({
            "gene":            gene,
            "type":            gtype,
            "genetic_code":    code,
            "start":           p["start"],
            "end":             p["end"],
            "length":          p["length"],
            "reading_frame":   frame,
            "codon_positions": (
                _codon_positions(p["start"], p["end"], frame)
                if gtype != "rRNA" else None
            ),
        })

    out_path = Path(output_dir) / "partition_metadata.json"
    with open(out_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"Partition metadata written: {out_path}")
    return str(out_path)


def generate_mrbayes_charsets(metadata_path: str,
                               output_path: str = None) -> str:
    """
    Read partition_metadata.json and generate a MrBayes charset block
    defining which supermatrix columns belong to which gene and codon
    position.
 
    For protein-coding genes where a reading frame was recorded, each
    gene is split into three charsets using the \3 stepping syntax:
        charset CYTB_pos1 = 502-1640\3;   [every 3rd col from 502]
        charset CYTB_pos2 = 500-1640\3;   [every 3rd col from 500]
        charset CYTB_pos3 = 501-1640\3;   [every 3rd col from 501]
 
    This lets MrBayes apply different substitution models to each codon
    position — important because position 3 evolves much faster than
    positions 1 and 2 due to synonymous substitutions.
 
    For rRNA genes (12S, 16S) a single undivided charset is written:
        charset 12S_rRNA = 14500-15518;
 
    For protein-coding genes with no reading frame recorded, a single
    undivided charset is written with a warning comment.
 
    Parameters
    ----------
    metadata_path : str         path to partition_metadata.json
                                (written by write_partition_metadata)
    output_path   : str or None if provided, also writes the charset
                                block to this file for pasting into
                                your MrBayes NEXUS file
 
    Returns
    -------
    str   the complete charset block text (also printed to console)
    """
    with open(metadata_path) as f:
        metadata = json.load(f)
 
    lines = [
        "begin mrbayes;",
        "",
        "  [Gene positions in the concatenated supermatrix]",
        "  [Protein-coding genes split by codon position (1, 2, 3)]",
        "  [rRNA genes kept as single undivided blocks]",
        "",
        "  [define charsets]",
    ]
 
    for g in metadata["genes"]:
        gene  = g["gene"]
        start = g["start"]
        end   = g["end"]
        gtype = g["type"]
        frame = g.get("reading_frame")
        cp    = g.get("codon_positions")  # dict with keys "1","2","3" (JSON str keys)
 
        if gtype == "rRNA":
            # rRNA — single block, no codon positions
            lines.append(f"  charset {gene} = {start}-{end};")
 
        elif cp is not None:
            # Protein-coding with known reading frame — split by codon position
            # JSON stores integer keys as strings, so access as cp["1"] etc.
            lines.append(f"  charset {gene}_pos1 = {cp['1']};  "
                         f"[{gene} frame={frame} codon position 1]")
            lines.append(f"  charset {gene}_pos2 = {cp['2']};  "
                         f"[{gene} frame={frame} codon position 2]")
            lines.append(f"  charset {gene}_pos3 = {cp['3']};  "
                         f"[{gene} frame={frame} codon position 3]")
 
        else:
            # Protein-coding but reading frame not recorded — single block
            lines.append(f"  charset {gene} = {start}-{end};  "
                         f"[WARNING: reading frame not recorded — add to write_partition_metadata]")
 
    lines.append("")
    lines.append("end;")
 
    block = "\n".join(lines)
    print(block)
 
    if output_path:
        with open(output_path, "w") as f:
            f.write(block + "\n")
        print(f"\nCharset block written to: {output_path}")
 
    return block






"""
# ===========================================================================
# Example usage
# ===========================================================================
if __name__ == "__main__":

    result_relaxed = big_fat_file_maker(
    raw_fol      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/final_nex/relaxed",
    csv_path     = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/final_accession_relaxed.csv",
    output       = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/concatenated/relaxed",
    datatype     = "dna",
    file_suffix  = "_aln.nex",
    trim = False,
    )

    result_harsh = big_fat_file_maker(
        raw_fol      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/final_nex/harsh",
        csv_path     = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/final_accession_harsh.csv",
        output       = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/concatenated/harsh",
        datatype     = "dna",
        file_suffix  = "_aln.nex",
        trim = False,
    )


    write_partition_metadata(
        partition_info = result_relaxed["partition_info"],
        output_dir     = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/concatenated/relaxed",
        reading_frames = {
            "APOB": 2, "BDNF": 3, "BRCA1": 1, "GHR": 1, "RAG1": 1, "RAG2": 3,
            "ATP6": 3, "ATP8": 3, "COX1": 3,  "COX2": 3, "COX3": 3, "CYTB": 3,
            "ND1":  3, "ND2":  3, "ND3":  3,  "ND4":  3, "ND4L": 3, "ND5":  3,
            "ND6":  1,
        }
    )
    write_partition_metadata(
        partition_info = result_harsh["partition_info"],
        output_dir     = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/concatenated/harsh",
        reading_frames = {
            "APOB": 3, "BDNF": 3, "BRCA1": 3, "GHR": 1, "RAG1": 3, "RAG2": 3,
            "ATP6": 3, "ATP8": 3, "COX1": 3,  "COX2": 3, "COX3": 1, "CYTB": 1,
            "ND1":  1, "ND2":  1, "ND3":  3,  "ND4":  3, "ND4L": 3, "ND5":  3,
            "ND6":  1,
        }
    )

"""


generate_mrbayes_charsets(
    metadata_path= "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/concatenated/relaxed/partition_metadata.json",
    output_path= "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/concatenated/relaxed/mrbayes.txt",
    )

generate_mrbayes_charsets(
    metadata_path= "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/concatenated/harsh/partition_metadata.json",
    output_path= "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/concatenated/harsh/mrbayes.txt"
    )