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
 
 
def gene_parser(expected_gene_symbol: str, hit_description: str, glossary: dict) -> dict:
    """
    Assess gene-label support for a BLAST hit and return a log-likelihood ratio.
 
    Strategy
    --------
    1.  Normalise the expected gene to its canonical form via the glossary.
    2.  Search the hit description for *all* glossary symbols (not just
        bracketed ones) using whole-word regex — GenBank descriptions use
        many formats, e.g.:
            "cytochrome b (cytb) gene"
            "NADH dehydrogenase subunit 2 (ND2)"
            "COX1 gene, partial cds"
            "cytochrome oxidase subunit I gene"
    3.  Collect every canonical gene name found in the description and check
        whether any matches the expected canonical name.
 
    Parameters
    ----------
    expected_gene_symbol : str   e.g. "CYTB", "cytb", "cytochrome b"
    hit_description      : str   the BLAST hit_def string
    glossary             : dict  output of load_glossary()
 
    Returns
    -------
    dict with keys:
        "expected_canonical"  : str   canonical form of the expected gene
        "found_canonicals"    : list  all canonical gene names found in description
        "match"               : bool  True if expected canonical is in found list
        "llr"                 : float log-likelihood ratio for this hit's gene label
    """
    # Resolve the expected symbol to its canonical name
    expected_canonical = glossary.get(expected_gene_symbol.lower())
    if expected_canonical is None:
        # Symbol not in glossary — can't make a judgement
        return {
            "expected_canonical": expected_gene_symbol,
            "found_canonicals":   [],
            "match":              False,
            "llr":                GENE_LLR["unknown"],
        }
 
    # Search description for every known symbol
    found_canonicals = []
    desc_lower = hit_description.lower()
 
    for raw_sym, canon in glossary.items():
        # Whole-word / whole-phrase match (handles multi-word symbols)
        pattern = rf'(?<!\w){re.escape(raw_sym)}(?!\w)'
        if re.search(pattern, desc_lower, re.IGNORECASE):
            if canon not in found_canonicals:
                found_canonicals.append(canon)
 
    match = expected_canonical in found_canonicals
 
    if not found_canonicals:
        llr = GENE_LLR["unknown"]   # description has no recognisable gene label
    elif match:
        llr = GENE_LLR["match"]
    else:
        llr = GENE_LLR["no_match"]
 
    return {
        "expected_canonical": expected_canonical,
        "found_canonicals":   found_canonicals,
        "match":              match,
        "llr":                llr,
    }

#print(gene_parser("APOB", "Anourosorex squamipes isolate 1 ApoB (ApoB) gene, partial cds", load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv")))