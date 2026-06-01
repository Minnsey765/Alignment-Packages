import re
import csv

def load_taxonomy(taxonomy_csv: str) -> dict:
    """
    Load the taxonomy CSV into a dict keyed by genus.
 
    Accepts comma- or tab-separated files, with or without a UTF-8 BOM.
    Column names are matched case-insensitively and with/without hyphens,
    so "Sub-family", "Subfamily", "sub_family" etc. all work.
 
    Expected columns (any order): Species, Genus, Sub-family, Family, Order
 
    Returns
    -------
    dict mapping genus (str) -> {"subfamily": str, "family": str, "order": str}
    """
    # ── 1. auto-detect delimiter ──────────────────────────────────────────────
    with open(taxonomy_csv, newline="", encoding="utf-8-sig") as f:  # utf-8-sig strips BOM
        sample = f.read(2048)
 
    delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
 
    # ── 2. normalised-key helper for fuzzy column matching ────────────────────
    def _norm(s: str) -> str:
        return re.sub(r"[\s\-_]", "", s).lower()
 
    COLUMN_MAP = {
        _norm("Genus"):      "genus",
        _norm("Sub-family"): "subfamily",
        _norm("Family"):     "family",
        _norm("Order"):      "order",
    }
 
    with open(taxonomy_csv, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
 
        if reader.fieldnames is None:
            raise ValueError(f"Could not read headers from '{taxonomy_csv}'")
 
        # Map actual header strings -> our internal keys
        # Skip empty strings — Excel often writes phantom trailing columns
        header_lookup = {}
        for raw_header in reader.fieldnames:
            if not raw_header or not raw_header.strip():
                continue
            norm = _norm(raw_header)
            if norm in COLUMN_MAP:
                header_lookup[COLUMN_MAP[norm]] = raw_header
 
        missing = [k for k in ("genus", "subfamily", "family", "order")
                   if k not in header_lookup]
        if missing:
            raise ValueError(
                f"Could not find columns {missing} in '{taxonomy_csv}'.\n"
                f"Headers found: {reader.fieldnames}\n"
                f"Check delimiter (detected: {repr(delimiter)}) and spelling."
            )
 
        # ── 3. load rows ──────────────────────────────────────────────────────
        taxonomy = {}
        for row in reader:
            genus = row[header_lookup["genus"]].strip()
            if genus and genus not in taxonomy:
                taxonomy[genus] = {
                    "subfamily": row[header_lookup["subfamily"]].strip(),
                    "family":    row[header_lookup["family"]].strip(),
                    "order":     row[header_lookup["order"]].strip(),
                }
 
    return taxonomy
 
 
def get_taxonomic_ranks(species_name: str, taxonomy: dict) -> dict | None:
    """
    Given a species_name string (e.g. "Sorex_araneus" or "Erinaceus_sp"),
    return the taxonomic ranks dict for its genus, or None if not found.
    """
    genus = species_name.split("_")[0]
    return taxonomy.get(genus, None)
 
 
# ─────────────────────────────────────────────────────────────────────────────
# SPECIES PARSER
# ─────────────────────────────────────────────────────────────────────────────
 
# Taxonomic rank hierarchy used for scoring.
# Higher index = more distantly related = weaker support for genuineness.
RANK_ORDER = ["genus", "subfamily", "family", "order"]
 
# Log-likelihood ratios (LLR) for each level of taxonomic match.
# These represent log( P(hit at rank X | genuine) / P(hit at rank X | artefact) ).
# Values are intentionally asymmetric: a genus-level match is strong positive
# evidence; matching only at order level is weak; no match is negative evidence.
# ──────────────────────────────────────────────────────────────────────────────
# These are *prior estimates* to be replaced once you have calibration data.
# Fit Beta distributions to genuine/artefact hit sets and recompute empirically.
# ──────────────────────────────────────────────────────────────────────────────
TAXONOMIC_LLR = {
    "species":   3.0,   # exact species match — strongest positive evidence
    "genus":     2.0,   # genus-level match — strong positive evidence
    "subfamily": 1.0,   # moderate positive evidence
    "family":    0.5,   # weak positive evidence
    "order":     0.1,   # very weak positive evidence
    "none":     -2.0,   # negative evidence (hit is from unrelated taxon)
}
 
 
def species_parser(species_name: str, hit_description: str, taxonomy: dict) -> dict:
    """
    Assess taxonomic support for a BLAST hit and return a log-likelihood ratio.
 
    Behaviour
    ---------
    1.  Extract the genus from species_name (handles both "Genus_species"
        and "Genus_sp" correctly — genus-only matching for "sp" entries).
    2.  Look up the full taxonomic chain (subfamily, family, order) from
        the taxonomy dict loaded via load_taxonomy().
    3.  Search the hit description for each rank in order from most to
        least specific, and return the LLR for the closest match found.
 
    Parameters
    ----------
    species_name    : str   e.g. "Sorex_araneus" or "Erinaceus_sp"
    hit_description : str   the BLAST hit_def string
    taxonomy        : dict  output of load_taxonomy()
 
    Returns
    -------
    dict with keys:
        "best_match_rank"  : str   one of "species", "genus", "subfamily",
                                   "family", "order", or "none"
        "matched_term"     : str   the actual string matched in the description
        "llr"              : float log-likelihood ratio for this hit's taxonomy
        "is_sp_level"      : bool  True if query is only identified to genus level
    """
    # Accept either a pre-loaded dict or a file path string
    if isinstance(taxonomy, str):
        taxonomy = load_taxonomy(taxonomy)
 
    parts = species_name.split("_")
    genus = parts[0]
    is_sp = len(parts) < 2 or parts[1].lower() == "sp"
 
    # Build ordered list of (rank_name, term_to_search_for), most specific first.
    # Species-level check is only meaningful when the query is identified to species.
    ranks_to_check = []
    if not is_sp:
        species_epithet = parts[1]
        genus_species = f"{genus} {species_epithet}"
        ranks_to_check.append(("species", genus_species))
    ranks_to_check.append(("genus", genus))
 
    tax_info = taxonomy.get(genus)
    if tax_info:
        ranks_to_check += [
            ("subfamily", tax_info["subfamily"]),
            ("family",    tax_info["family"]),
            ("order",     tax_info["order"]),
        ]
 
    # Search description for each rank (most specific first)
    best_rank = "none"
    matched_term = ""
    for rank_name, term in ranks_to_check:
        # Word-boundary search so "Sorex" doesn't match "Soriculus"
        if term and re.search(rf'\b{re.escape(term)}\b', hit_description, re.IGNORECASE):
            best_rank = rank_name
            matched_term = term
            break  # stop at the closest (most specific) match
 
    return {
        "best_match_rank": best_rank,
        "matched_term":    matched_term,
        "llr":             TAXONOMIC_LLR[best_rank],
        "is_sp_level":     is_sp,
    }


#taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
#print(list(taxonomy.items())[:3])

#print(load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"))
#print(species_parser("Anourosorex_squamipes", "Anourosorex squamipes isolate 1 ApoB (ApoB) gene, partial cds", load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")))