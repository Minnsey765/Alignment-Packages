#test if gene name or symbol is in a species definition
import re
import csv


# ─────────────────────────────────────────────────────────────────────────────
# GENE PARSER
# ─────────────────────────────────────────────────────────────────────────────
 
# LLR values for gene matching.
# ──────────────────────────────────────────────────────────────────────────────
# Again, *prior estimates* — replace with empirically fitted values once you
# have calibration data from your genuine/artefact hit sets.
# ──────────────────────────────────────────────────────────────────────────────
GENE_LLR = {
    "match":    2.5,   # hit description contains the expected gene → strong support
    "no_match": -1.5,  # gene absent or different → moderate negative evidence
    "unknown":  0.0,   # description unreadable / gene not annotated → neutral
}
 
 
def load_glossary(glossary_csv: str) -> dict:
    """
    Load the two-row gene glossary CSV into a symbol → canonical_name dict.
 
    The glossary CSV format (tab- or comma-separated, two rows):
        Row 0 (header): raw symbol variants used in GenBank descriptions
        Row 1          : corresponding canonical gene names
 
    Returns
    -------
    dict mapping raw_symbol (str, lowercase-stripped) -> canonical_name (str)
 
    Notes
    -----
    The original code did `[row for row in reader][0]` which read *only*
    the first data row (row index 0 after the header) as a dict whose keys
    are the header values.  That meant symbol_key was e.g.:
        {"12S ribosomal RNA": "12S_rRNA", "16S ribosomal RNA": "16S_rRNA", ...}
    which was actually correct for a two-row file used that way — but it
    relied on csv.DictReader treating row 0 as the header and row 1 as the
    single data row.  This reimplementation makes that intent explicit and
    adds case-insensitive lookup.
    """
    glossary = {}
    with open(glossary_csv, newline="", encoding="utf-8") as f:
        # Read all rows without a header assumption
        reader = csv.reader(f)
        rows = [row for row in reader if any(cell.strip() for cell in row)]
 
    if len(rows) < 2:
        raise ValueError(
            f"Glossary file '{glossary_csv}' must have at least two rows: "
            "symbols row and canonical names row."
        )
 
    symbols   = [cell.strip() for cell in rows[0]]
    canonical = [cell.strip() for cell in rows[1]]
 
    if len(symbols) != len(canonical):
        raise ValueError(
            "Glossary rows have different lengths "
            f"({len(symbols)} symbols vs {len(canonical)} canonical names)."
        )
 
    for sym, can in zip(symbols, canonical):
        if sym:
            glossary[sym.lower()] = can  # lowercase keys for case-insensitive lookup
 
    return glossary
 
 
# Mitochondrial genes that are routinely found in complete mitogenome
# records. If a hit description contains mitogenome keywords, any of
# these genes scores as a match regardless of whether the gene symbol
# appears explicitly in the description.
MITO_GENES = {
    "COX1", "COX2", "COX3",
    "ND1",  "ND2",  "ND3",  "ND4", "ND4L", "ND5", "ND6",
    "ATP6", "ATP8",
    "CYTB",
    "12S_rRNA", "16S_rRNA",
    # Common synonyms
    "12S",  "16S",
    "COI",  "COII", "COIII",
    "NAD1", "NAD2", "NAD3", "NAD4", "NAD4L", "NAD5", "NAD6",
}

# Keywords in BLAST hit descriptions that indicate a complete or
# partial mitochondrial genome record
def gene_parser(expected_gene: str,
                hit_description: str,
                glossary: dict) -> dict:
    """
    Score a BLAST hit description for gene identity.

    Priority order:
      1. Expected gene symbol or any synonym found → exact match (+2.0)
      2. Different gene symbol found → mismatch (-2.0), UNLESS both
         the expected and found genes are mitochondrial and the
         description indicates a mitochondrial record (they coexist
         in mitogenomes)
      3. No gene symbol found but mitochondrial origin confirmed
         AND expected gene is mitochondrial → weak match (+1.0)
      4. No gene symbol, not mitochondrial → unknown (0.0)

    Parameters
    ----------
    expected_gene   : str   canonical gene name e.g. "12S_rRNA", "CYTB"
    hit_description : str   full BLAST hit description string
    glossary        : dict  from load_glossary() — maps raw_symbol ->
                            canonical_name

    Returns
    -------
    dict with llr, match_type, matched_gene
    """
    desc_lower     = hit_description.lower()
    desc_upper     = hit_description.upper()
    expected_clean = expected_gene.strip().upper().replace("-", "_")

    # ── Build reverse lookup: canonical → set of all raw symbols ─────────────
    # The glossary maps raw_symbol -> canonical, so we invert it to
    # find all raw symbols that map to a given canonical gene name.

    def get_names_for_canonical(canonical: str) -> set:
        """
        Return all raw symbol strings that map to this canonical gene,
        plus short-form parts of the canonical name itself.
        """
        target = canonical.upper()
        names  = {target}

        for raw_symbol, canon_val in glossary.items():
            if isinstance(canon_val, str):
                if canon_val.upper() == target:
                    names.add(raw_symbol.upper())
            elif isinstance(canon_val, dict):
                c = canon_val.get("canonical", "")
                if c and c.upper() == target:
                    names.add(raw_symbol.upper())

        # Add parts of compound canonical names
        # 12S_rRNA → {12S, RRNA}
        # 16S_rRNA → {16S, RRNA}
        for part in target.replace("_", " ").split():
            if len(part) >= 2:
                names.add(part)

        return names

    # Build lookup for all canonical gene names present in glossary
    # so Priority 2 can check any gene efficiently
    all_canonicals = set()
    for canon_val in glossary.values():
        if isinstance(canon_val, str) and canon_val:
            all_canonicals.add(canon_val.upper())
        elif isinstance(canon_val, dict):
            c = canon_val.get("canonical", "")
            if c:
                all_canonicals.add(c.upper())

    expected_names = get_names_for_canonical(expected_clean)

    # ── Mitochondrial status (needed by Priority 2 and 3) ────────────────────
    MITO_TERMS   = ("mitochondri", "mitogenome")
    is_mito      = any(term in desc_lower for term in MITO_TERMS)
    is_mito_gene = (
        expected_clean in MITO_GENES
        or expected_clean.replace("_RRNA", "") in MITO_GENES
        or expected_clean.replace("RRNA_", "") in MITO_GENES
    )

    # ── Priority 1: expected gene symbol or synonym found ────────────────────
    for name in expected_names:
        if name and len(name) > 1 and name in desc_upper:
            return {
                "llr":          2.0,
                "match_type":   "exact",
                "matched_gene": name,
            }

    # ── Priority 2: different gene symbol found ───────────────────────────────
    for other_canonical in all_canonicals:
        if other_canonical == expected_clean:
            continue

        other_names = get_names_for_canonical(other_canonical)

        # Remove any names that overlap with the expected gene's names
        # to avoid cross-matching synonyms shared between genes
        other_names -= expected_names

        for name in other_names:
            if name and len(name) > 2 and name in desc_upper:

                other_is_mito = (
                    other_canonical in MITO_GENES
                    or other_canonical.replace("_RRNA", "")
                    in MITO_GENES
                )

                # Both mito genes in a mito record — coexist, not a
                # mismatch. Skip and keep checking.
                if is_mito_gene and other_is_mito and is_mito:
                    continue

                return {
                    "llr":          -2.0,
                    "match_type":   "mismatch",
                    "matched_gene": name,
                }

    # ── Priority 3: mitochondrial fallback ───────────────────────────────────
    if is_mito and is_mito_gene:
        return {
            "llr":          1.0,
            "match_type":   "mitochondrial",
            "matched_gene": None,
        }

    if is_mito and not is_mito_gene:
        return {
            "llr":          -2.0,
            "match_type":   "mismatch",
            "matched_gene": None,
        }

    # ── Priority 4: unknown ───────────────────────────────────────────────────
    return {
        "llr":          0.0,
        "match_type":   "unknown",
        "matched_gene": None,
    }
#print(gene_parser("APOB", "Anourosorex squamipes isolate 1 ApoB (ApoB) gene, partial cds", load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv")))