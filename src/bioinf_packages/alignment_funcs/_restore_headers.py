# restore_headers.py
#
# Functions for restoring readable headers to alignment files after
# trimming in AliView/Mesquite, and looping over a directory of files.
#
# Header format produced (medium):
#   FASTA : >Sorex_araneus||CYTB||NC_027963
#   NEXUS : Sorex_araneus__CYTB__NC_027963
#
#   Double pipe (FASTA) or double underscore (NEXUS) separates fields.
#   This avoids colons (misread as sequence) and single underscores
#   (appear in species names). big_fat_file_maker_v2 extracts the
#   accession by splitting on || or __ and taking the last field.
#
# Functions:
#   restore_trimmed_headers()   — restore headers in a single file
#   loop_restore_headers()      — loop over a directory of trimmed files
#   _build_replacement()        — internal: build medium/full/species header
#   _replace_fasta_headers()    — internal: replace headers in FASTA
#   _replace_nexus_headers()    — internal: replace headers in NEXUS

import json
import re
import traceback
from pathlib import Path


# ===========================================================================
# Internal helpers
# ===========================================================================

def _build_replacement(full_id: str, fmt: str) -> str:
    """
    Convert a full pipe-delimited header (from run_clustalo id_map)
    to the chosen output format.

    Input (original id_map format):
        Sorex_araneus|CYTB|orien:+|accession:NC_027963

    Output:
        "medium"       → Sorex_araneus||CYTB||NC_027963
                         (NEXUS-safe: Sorex_araneus__CYTB__NC_027963)
        "species_only" → Sorex_araneus
        "full"         → Sorex_araneus|CYTB|orien:+|accession:NC_027963

    The double-pipe / double-underscore delimiter is unambiguous:
        - Double pipe never appears in FASTA sequence data or accessions
        - Double underscore never appears in species names or accessions
        - The accession is the last field with no prefix, so
          big_fat_file_maker_v2 extracts it by splitting on || or __
          and taking the final element.
    """
    if fmt == "full":
        return full_id

    parts = full_id.split("|")

    # Species: first pipe-segment e.g. "Sorex_araneus"
    species = parts[0].strip() if parts else full_id

    if fmt == "species_only":
        return species

    # Medium: extract gene (second segment, no colon) and accession
    gene = ""
    for p in parts[1:]:
        if ":" not in p and p.strip():
            gene = p.strip()
            break

    accession = ""
    for p in parts:
        if p.startswith("accession:"):
            accession = p.split(":", 1)[1].strip()
            break

    # Use double-pipe as delimiter — unambiguous in FASTA context.
    # _replace_nexus_headers converts || → __ for NEXUS compatibility.
    segments = [s for s in [species, gene, accession] if s]
    return "||".join(segments)


def _replace_fasta_headers(content: str, replacement_map: dict) -> str:
    """
    Replace short IDs on FASTA header lines (lines starting with >).
    Double-pipe delimiters are kept as-is in FASTA.
    """
    lines = content.splitlines(keepends=True)
    out   = []
    for line in lines:
        if line.startswith(">"):
            header_body = line[1:].rstrip("\n\r")
            tokens      = header_body.split()
            short_id    = tokens[0] if tokens else ""
            if short_id in replacement_map:
                new_header = replacement_map[short_id]
                rest       = " ".join(tokens[1:])
                line       = f">{new_header}" + (f" {rest}" if rest else "") + "\n"
        out.append(line)
    return "".join(out)


def _replace_nexus_headers(content: str, replacement_map: dict) -> str:
    """
    Replace short IDs in a NEXUS file (MATRIX block and TAXLABELS).

    Pipes are not allowed in NEXUS taxon names, so double-pipe (||)
    is converted to double-underscore (__). This gives headers like:
        Sorex_araneus__CYTB__NC_027963
    which big_fat_file_maker_v2 parses by splitting on __ and taking
    the last element as the accession.

    Replaces longest keys first to avoid partial matches.
    """
    for short_id, new_id in sorted(replacement_map.items(),
                                    key=lambda x: len(x[0]),
                                    reverse=True):
        pattern     = r'(?<![A-Za-z0-9_])' + re.escape(short_id) + r'(?![A-Za-z0-9_])'
        new_id_safe = new_id.replace("||", "__")  # pipes → double underscore for NEXUS
        content     = re.sub(pattern, new_id_safe, content)
    return content


# ===========================================================================
# Function 1 — restore headers in a single file
# ===========================================================================

def restore_trimmed_headers(trimmed_file: str,
                             id_map_path: str,
                             output_path: str,
                             header_format: str = "medium") -> str:
    """
    Restore readable headers in a trimmed alignment file, replacing short
    IDs (e.g. S01_Sorex_AB1234) with a cleaner format derived from the
    original full headers stored in the ID map JSON produced by run_clustalo.

    Header format options
    ---------------------
    "medium"  (default)
        FASTA : >Sorex_araneus||CYTB||NC_027963
        NEXUS : Sorex_araneus__CYTB__NC_027963
        Species, gene, and accession separated by double-pipe (FASTA) or
        double-underscore (NEXUS). No colons in the header so nothing is
        misread as part of the sequence. big_fat_file_maker_v2 extracts
        the accession by splitting on || or __ and taking the last field.

    "species_only"
        >Sorex_araneus
        Just the binomial — useful for tree viewers.

    "full"
        Restores the complete original header exactly as produced by
        run_clustalo. May be too long for some programs (>30 chars).

    Parameters
    ----------
    trimmed_file  : str   path to AliView/Mesquite-trimmed .nex or .fasta
                          file containing short IDs
    id_map_path   : str   path to the _id_map.json produced by run_clustalo
    output_path   : str   where to write the restored file
    header_format : str   "medium" (default), "species_only", or "full"

    Returns
    -------
    str   path to the restored output file
    """
    with open(id_map_path, "r") as f:
        id_map = json.load(f)   # short_id -> full_original_id

    replacement_map = {
        short_id: _build_replacement(full_id, header_format)
        for short_id, full_id in id_map.items()
    }

    with open(trimmed_file, "r") as f:
        content = f.read()

    is_nexus = "begin data" in content.lower() or "begin characters" in content.lower()

    if is_nexus:
        content = _replace_nexus_headers(content, replacement_map)
    else:
        content = _replace_fasta_headers(content, replacement_map)

    with open(output_path, "w") as f:
        f.write(content)

    print(f"  [restore] → {Path(output_path).name}  "
          f"({len(replacement_map)} IDs, format='{header_format}')")
    return output_path


# ===========================================================================
# Function 2 — loop over a directory of trimmed files
# ===========================================================================

def loop_restore_headers(trimmed_dir: str,
                          id_map_dir: str,
                          output_dir: str,
                          header_format: str = "medium") -> dict:
    """
    Loop through a directory of AliView/Mesquite-trimmed NEXUS files and
    restore readable headers using the ID map JSONs from run_clustalo.

    File naming conventions
    -----------------------
    Trimmed files  : [gene]_short_aln.nex   (in trimmed_dir)
    ID map files   : [gene]_id_map.json     (in id_map_dir)
    Output files   : [gene]_aln.nex         (in output_dir)

    Output files named [gene]_aln.nex so big_fat_file_maker finds them
    directly with file_suffix="_aln.nex".

    Parameters
    ----------
    trimmed_dir   : str   directory containing *_short_aln.nex files
    id_map_dir    : str   directory containing *_id_map.json files
                          (same directory as run_clustalo output)
    output_dir    : str   directory to write restored [gene]_aln.nex files
    header_format : str   "medium" (default), "species_only", or "full"

    Returns
    -------
    dict with restored, skipped, failed (lists of gene names), output_dir
    """
    trimmed_dir = Path(trimmed_dir)
    id_map_dir  = Path(id_map_dir)
    output_dir  = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    trimmed_files = sorted(trimmed_dir.glob("*_short_aln.nex"))

    if not trimmed_files:
        print(f"[loop_restore_headers] No *_short_aln.nex files found in {trimmed_dir}")
        return {"restored": [], "skipped": [], "failed": [], "output_dir": str(output_dir)}

    print(f"[loop_restore_headers]")
    print(f"  Trimmed dir  : {trimmed_dir}")
    print(f"  ID map dir   : {id_map_dir}")
    print(f"  Output dir   : {output_dir}")
    print(f"  Format       : {header_format}")
    print(f"  Found {len(trimmed_files)} trimmed files\n")

    restored, skipped, failed = [], [], []

    for trimmed_path in trimmed_files:
        gene        = trimmed_path.name.replace("_short_aln.nex", "")
        id_map_path = id_map_dir / f"{gene}_id_map.json"

        if not id_map_path.exists():
            print(f"  [skip]  {gene:<25} — id_map not found: {id_map_path.name}")
            skipped.append(gene)
            continue

        output_path = output_dir / f"{gene}_aln.nex"

        try:
            restore_trimmed_headers(
                trimmed_file  = str(trimmed_path),
                id_map_path   = str(id_map_path),
                output_path   = str(output_path),
                header_format = header_format,
            )
            restored.append(gene)

        except Exception as e:
            print(f"  [error] {gene:<25} — {e}")
            traceback.print_exc()
            failed.append(gene)

    print(f"\n  Done : {len(restored)} restored, "
          f"{len(skipped)} skipped, {len(failed)} failed")

    return {
        "restored":   restored,
        "skipped":    skipped,
        "failed":     failed,
        "output_dir": str(output_dir),
    }


# ===========================================================================
# Example usage
# ===========================================================================







# ===========================================================================
# Example usage
# ===========================================================================
if __name__ == "__main__":

    # ── Step 1: Restore headers for all trimmed files ──────────────────────
    loop_restore_headers(
        trimmed_dir   = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/edited_nexus/relaxed",       # <-- your Mesquite/AliView trimmed files
        id_map_dir    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/aligned_fastas/relaxed",            # <-- where run_clustalo saved id_map JSONs
        output_dir    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/final_nex/relaxed",       # <-- where to write restored files
        header_format = "medium",
    )
    loop_restore_headers(
        trimmed_dir   = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/edited_nexus/harsh",       # <-- your Mesquite/AliView trimmed files
        id_map_dir    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/aligned_fastas/harsh",            # <-- where run_clustalo saved id_map JSONs
        output_dir    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/final_nex/harsh",       # <-- where to write restored files
        header_format = "medium",
    )


    # ── Step 2 (optional): Write metadata with reading frames ──────────────
    # Fill in the reading frames you noted down from AliView
    # Only needed for protein-coding genes — leave out rRNA
    reading_frames = {
        "APOB":  1,    # <-- replace with your noted reading frames
        "BDNF":  1,
        "BRCA1": 1,
        "GHR":   1,
        "RAG1":  1,
        "RAG2":  1,
        "ATP6":  1,
        "ATP8":  1,
        "COX1":  1,
        "COX2":  1,
        "COX3":  1,
        "CYTB":  1,
        "ND1":   1,
        "ND2":   1,
        "ND3":   1,
        "ND4":   1,
        "ND4L":  1,
        "ND5":   1,
        "ND6":   1,
    }

    # After running big_fat_file_maker, call this with its partition_info output
    # write_partition_metadata(
    #     partition_info = result["partition_info"],
    #     output_dir     = "supermatrix/",
    #     reading_frames = reading_frames,
    # )

    # ── Step 3 (optional): Generate MrBayes charset block ─────────────────
    # generate_mrbayes_charsets("supermatrix/partition_metadata.json")