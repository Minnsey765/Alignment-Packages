import re
import csv
from pathlib import Path

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
 
 
def species_parser(species_name: str,
                   hit_description: str,
                   taxonomy: dict) -> dict:
    """
    Assess taxonomic support for a BLAST hit by comparing the query
    taxon's taxonomic chain against the hit taxon's taxonomic chain.

    Rather than searching the hit description string for taxonomic
    terms, this function extracts the genus from the hit description,
    looks it up in the taxonomy dict, and compares the two taxonomic
    chains. This correctly scores hits like:
        query=Erinaceus_europaeus, hit=Hylomys_suillus
    as a family-level match (both Erinaceidae) rather than "none".

    Atopogale and Solenodon are treated as the same genus since many
    GenBank records use Atopogale cubana for what is also known as
    Solenodon cubanus.

    Parameters
    ----------
    species_name    : str   e.g. "Sorex_araneus" or "Erinaceus_sp"
                            May include subspecies e.g.
                            "Solenodon_paradoxus_paradoxus" — the
                            third token is ignored.
    hit_description : str   BLAST hit description string e.g.
                            "Hylomys suillus mitochondrion, complete
                            genome" or "Sorex araneus BDNF mRNA"
    taxonomy        : dict  from load_taxonomy()

    Returns
    -------
    dict with:
        best_match_rank : str   "species", "genus", "subfamily",
                                "family", "order", or "none"
        matched_term    : str   the term that matched
        llr             : float log-likelihood ratio
        query_genus     : str   genus extracted from species_name
        hit_genus       : str   genus extracted from hit description
                                (None if not found in taxonomy)
        is_sp_level     : bool  True if query only identified to genus
    """
    if isinstance(taxonomy, str):
        taxonomy = load_taxonomy(taxonomy)

    # ── Parse query taxon ─────────────────────────────────────────────────────
    # Handle Genus_species, Genus_species_subspecies, Genus_sp
    # Subspecies token (third) is ignored — treat as species level
    parts   = species_name.replace(" ", "_").split("_")
    genus   = parts[0]
    is_sp   = len(parts) < 2 or parts[1].lower() in ("sp", "sp.")

    # Normalise Atopogale → Solenodon
    genus = _normalise_genus(genus)

    # Look up query taxonomy
    query_tax = taxonomy.get(genus, {})

    # ── Extract hit genus from description ────────────────────────────────────
    # BLAST descriptions typically start with "Genus species ..."
    # Extract the first word as a candidate genus and look it up
    hit_genus     = _extract_genus_from_description(hit_description,
                                                     taxonomy)
    hit_genus_raw = hit_genus   # before normalisation
    if hit_genus:
        hit_genus = _normalise_genus(hit_genus)

    hit_tax = taxonomy.get(hit_genus, {}) if hit_genus else {}

    # ── Compare taxonomic chains ──────────────────────────────────────────────
    # Work from most specific to least specific.
    # Species match: same genus + same species epithet in description
    best_rank    = "none"
    matched_term = ""

    if hit_genus and hit_genus == genus:
        # Same genus — check species epithet
        if not is_sp and len(parts) >= 2:
            species_epithet = parts[1]
            # Search description for genus + species
            genus_species = f"{genus} {species_epithet}"
            alt_gs        = f"{genus}_{species_epithet}"
            if (re.search(rf'\b{re.escape(genus_species)}\b',
                          hit_description, re.IGNORECASE)
                    or alt_gs.lower() in hit_description.lower()):
                best_rank    = "species"
                matched_term = genus_species
            else:
                best_rank    = "genus"
                matched_term = genus
        else:
            best_rank    = "genus"
            matched_term = genus

    elif hit_genus and query_tax and hit_tax:
        # Different genus — compare at subfamily, family, order
        q_subfamily = query_tax.get("subfamily", "")
        q_family    = query_tax.get("family", "")
        q_order     = query_tax.get("order", "")

        h_subfamily = hit_tax.get("subfamily", "")
        h_family    = hit_tax.get("family", "")
        h_order     = hit_tax.get("order", "")

        if (q_subfamily and h_subfamily
                and q_subfamily == h_subfamily):
            best_rank    = "subfamily"
            matched_term = q_subfamily
        elif q_family and h_family and q_family == h_family:
            best_rank    = "family"
            matched_term = q_family
        elif q_order and h_order and q_order == h_order:
            best_rank    = "order"
            matched_term = q_order
        else:
            best_rank    = "none"
            matched_term = ""

    elif hit_genus and not hit_tax:
        # Hit genus not in taxonomy — fall back to searching the
        # description string for the query's taxonomic terms, as
        # the original implementation did. This handles cases where
        # the hit is from a taxon not in our local taxonomy CSV.
        ranks_to_check = []
        if not is_sp and len(parts) >= 2:
            ranks_to_check.append(
                ("species", f"{genus} {parts[1]}")
            )
        ranks_to_check.append(("genus", genus))
        if query_tax:
            if query_tax.get("subfamily"):
                ranks_to_check.append(
                    ("subfamily", query_tax["subfamily"])
                )
            if query_tax.get("family"):
                ranks_to_check.append(
                    ("family", query_tax["family"])
                )
            if query_tax.get("order"):
                ranks_to_check.append(
                    ("order", query_tax["order"])
                )

        for rank_name, term in ranks_to_check:
            if term and re.search(rf'\b{re.escape(term)}\b',
                                   hit_description, re.IGNORECASE):
                best_rank    = rank_name
                matched_term = term
                break

    return {
        "best_match_rank": best_rank,
        "matched_term":    matched_term,
        "llr":             TAXONOMIC_LLR[best_rank],
        "query_genus":     genus,
        "hit_genus":       hit_genus,
        "is_sp_level":     is_sp,
    }


def _normalise_genus(genus: str) -> str:
    """
    Normalise genus synonyms to a canonical name.

    Atopogale is the newer name for the Cuban solenodon (previously
    Solenodon cubanus, now Atopogale cubana) but many GenBank records
    still use Solenodon. Treat both as Solenodon so they score
    correctly against each other and against the taxonomy CSV.
    """
    synonyms = {
        "Atopogale": "Solenodon",
    }
    return synonyms.get(genus, genus)

# Species-level synonyms — maps non-canonical name to canonical name.
# Used when two names refer to the same biological entity and should
# score identically against BLAST hits.
SPECIES_SYNONYMS = {
    "Crocidura_geldestani": "Crocidura_suaveolens",
    # add further synonyms here as needed
}

def _normalise_species(taxon: str) -> str:
    """
    Normalise species-level synonyms to a canonical binomial.

    Applied after _normalise_genus() so genus-level normalisation
    takes precedence. Input and output are Genus_species strings.
    """
    return SPECIES_SYNONYMS.get(taxon, taxon)


def _extract_genus_from_description(description: str,
                                     taxonomy: dict) -> str:
    """
    Extract a genus name from a BLAST hit description by scanning
    the first few words and checking against the taxonomy dict.

    BLAST descriptions typically start with "Genus species ...",
    but some start with accession numbers, strain names, or other
    prefixes. This function tries the first word, then the second,
    and returns the first one found in the taxonomy dict.

    If no match is found in the taxonomy dict, returns the first
    word of the description as a best guess (it may still be a
    valid genus not in our local taxonomy).

    Parameters
    ----------
    description : str   BLAST hit description
    taxonomy    : dict  from load_taxonomy()

    Returns
    -------
    str or None   best guess at the genus name
    """
    # Strip any leading accession or pipe characters
    clean = description.strip()
    for prefix in (">", "gb|", "ref|", "emb|", "dbj|"):
        if clean.startswith(prefix):
            clean = clean.split(None, 1)[-1].strip()

    words = clean.split()
    if not words:
        return None

    # Try first three words as candidate genus
    for word in words[:3]:
        # Strip trailing punctuation
        word_clean = re.sub(r'[^A-Za-z]', '', word)
        if len(word_clean) < 3:
            continue
        # Check in taxonomy dict (and normalise)
        candidate = _normalise_genus(word_clean.capitalize())
        if candidate in taxonomy:
            return candidate
        # Also check uncapitalised
        if word_clean.capitalize() in taxonomy:
            return word_clean.capitalize()

    # Fall back to first capitalised word as best guess
    for word in words[:5]:
        word_clean = re.sub(r'[^A-Za-z]', '', word)
        if len(word_clean) >= 3 and word_clean[0].isupper():
            return word_clean

    return words[0] if words else None


def flag_unknown_taxa(
        results_csv: str,
        taxonomy_csv: str,
        output_unknown_csv: str = None,
        auto_add_to_taxonomy: bool = True,
        min_hit_count: int = 1,
) -> dict:
    """
    Scan blast_results.csv for hit genera not present in taxonomy_data.csv
    and optionally add them as skeleton rows for manual curation.

    When auto_add_to_taxonomy=True, unknown genera are appended to
    taxonomy_data.csv with empty Sub-family, Family, and Order fields.
    Open taxonomy_data.csv afterwards and fill in the missing columns,
    then reload the taxonomy dict before rescoring.

    Parameters
    ----------
    results_csv          : str   path to blast_results.csv
    taxonomy_csv         : str   path to taxonomy_data.csv
    output_unknown_csv   : str   optional separate CSV of unknowns
                                 for review before adding to taxonomy
    auto_add_to_taxonomy : bool  if True, append skeleton rows for
                                 unknown genera to taxonomy_data.csv.
                                 Default True.
    min_hit_count        : int   only report/add genera that appear
                                 at least this many times as hits.
                                 Default 1 (report everything).

    Returns
    -------
    dict with:
        unknown_genera   : list   genus strings not in taxonomy
        hit_counts       : dict   genus -> total hit count
        rows_added       : int    rows appended to taxonomy_data.csv
    """
    from collections import Counter, defaultdict

    # ── Load known genera from taxonomy CSV ───────────────────────────────────
    known_genera    = set()
    existing_rows   = []
    taxonomy_fields = ["Species", "Genus", "Sub-family",
                       "Family", "Order"]

    with open(taxonomy_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        taxonomy_fields = reader.fieldnames or taxonomy_fields
        for row in reader:
            existing_rows.append(row)
            genus = row.get("Genus", "").strip()
            if genus:
                known_genera.add(genus)

    print(f"\n{'='*60}")
    print(f"FLAG UNKNOWN TAXA")
    print(f"{'='*60}")
    print(f"  Taxonomy CSV   : {taxonomy_csv}")
    print(f"  Known genera   : {len(known_genera)}")

    # ── Scan hit descriptions for unknown genera ──────────────────────────────
    genus_hit_counts    = Counter()
    genus_species_seen  = defaultdict(set)  # genus -> set of species epithets

    with open(results_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            hit_desc = row.get("hit_description", "").strip()
            hit_acc  = row.get("hit_accession", "").strip()

            if not hit_desc or hit_acc == "NO_HITS":
                continue

            # Extract genus from hit description.
            # BLAST descriptions typically start with
            # "Genus species ..." — take the first capitalised word.
            words = hit_desc.split()
            if not words:
                continue

            # First token is usually "Genus" if description starts
            # with binomial. Skip common non-genus prefixes.
            skip_prefixes = {
                "PREDICTED:", "UNVERIFIED:", "LOW", "QUALITY",
                "MULTISPECIES:", "TPA:", "TPA_inf:",
            }
            genus_candidate = None
            species_candidate = None

            for i, word in enumerate(words):
                if word in skip_prefixes or word.endswith(":"):
                    continue
                # A genus starts with a capital and is alphabetic
                if (word[0].isupper()
                        and word.replace("-", "").isalpha()
                        and len(word) > 2):
                    genus_candidate = word
                    # Species epithet is the next word if lowercase
                    if (i + 1 < len(words)
                            and words[i+1][0].islower()
                            and words[i+1].replace("-","").isalpha()):
                        species_candidate = words[i+1]
                    break

            if not genus_candidate:
                continue

            # Normalise using the same function as species_parser
            genus_norm = _normalise_genus(genus_candidate)

            if genus_norm not in known_genera:
                genus_hit_counts[genus_norm] += 1
                if species_candidate:
                    genus_species_seen[genus_norm].add(
                        species_candidate
                    )

    # Filter by minimum hit count
    unknown_genera = sorted(
        g for g, count in genus_hit_counts.items()
        if count >= min_hit_count
    )

    print(f"  Unknown genera : {len(unknown_genera)}")
    if min_hit_count > 1:
        print(f"  (min_hit_count = {min_hit_count})")

    if not unknown_genera:
        print(f"  All hit genera are in taxonomy_data.csv.")
        return {
            "unknown_genera": [],
            "hit_counts":     {},
            "rows_added":     0,
        }

    # ── Build skeleton rows for unknown genera ────────────────────────────────
    # One row per (genus, species) combination seen in hits.
    # If no species was extracted, add one row with empty species.
    new_rows = []
    print(f"\n  Unknown genera (genus: hit_count, species seen):")
    print(f"  {'-'*56}")

    for genus in unknown_genera:
        count   = genus_hit_counts[genus]
        species = sorted(genus_species_seen.get(genus, set()))
        print(f"  {genus:<30} hits={count:<6} "
              f"species=[{', '.join(species[:3])}"
              f"{'...' if len(species) > 3 else ''}]")

        if species:
            for sp in species:
                new_rows.append({
                    "Species":    sp,
                    "Genus":      genus,
                    "Sub-family": "",   # fill in manually
                    "Family":     "",   # fill in manually
                    "Order":      "",   # fill in manually
                })
        else:
            new_rows.append({
                "Species":    "",
                "Genus":      genus,
                "Sub-family": "",
                "Family":     "",
                "Order":      "",
            })

    # ── Write output unknown CSV for review ───────────────────────────────────
    if output_unknown_csv:
        Path(output_unknown_csv).parent.mkdir(
            parents=True, exist_ok=True
        )
        with open(output_unknown_csv, "w", newline="",
                  encoding="utf-8") as f:
            writer = csv.DictWriter(
                f, fieldnames=taxonomy_fields,
                extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(new_rows)
        print(f"\n  Unknown taxa written to: {output_unknown_csv}")
        print(f"  Fill in Sub-family, Family, Order columns,")
        print(f"  then set auto_add_to_taxonomy=True to append.")

    # ── Append skeleton rows to taxonomy CSV ──────────────────────────────────
    rows_added = 0
    if auto_add_to_taxonomy and new_rows:
        # Re-check which genera are already present to avoid
        # duplicates if this function is run more than once
        already_present = set()
        for row in existing_rows:
            g = row.get("Genus", "").strip()
            s = row.get("Species", "").strip()
            already_present.add((g, s))

        rows_to_add = [
            r for r in new_rows
            if (r["Genus"], r["Species"]) not in already_present
        ]

        if rows_to_add:
            with open(taxonomy_csv, "a", newline="",
                      encoding="utf-8") as f:
                writer = csv.DictWriter(
                    f, fieldnames=taxonomy_fields,
                    extrasaction="ignore"
                )
                writer.writerows(rows_to_add)
            rows_added = len(rows_to_add)
            print(f"\n  Added {rows_added} skeleton rows to "
                  f"{taxonomy_csv}")
            print(f"  Open the file and fill in Sub-family, "
                  f"Family, Order for each new genus.")
            print(f"  Then reload taxonomy with load_taxonomy() "
                  f"before rescoring.")
        else:
            print(f"\n  No new rows needed — all genera already "
                  f"in taxonomy CSV.")

    return {
        "unknown_genera": unknown_genera,
        "hit_counts":     dict(genus_hit_counts),
        "rows_added":     rows_added,
    }


"""
Run on the csv after retrieving the results from a blast run
"""
# Run after subsamples_to_blast() completes
#unknowns = flag_unknown_taxa(
#    results_csv          = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
#    taxonomy_csv         = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv",
#    output_unknown_csv   = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/unknown_taxa.csv",
#    auto_add_to_taxonomy = True,   # appends skeleton rows immediately
#    min_hit_count        = 3,      # ignore genera hit only once or
#                                   # twice — likely noise
#)

#taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
#print(list(taxonomy.items())[:3])

#print(load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"))
#print(species_parser("Anourosorex_squamipes", "Anourosorex squamipes isolate 1 ApoB (ApoB) gene, partial cds", load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")))