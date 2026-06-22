"""
fill_taxonomy.py

Two main functions:

1. fill_taxonomy_csv()
   Reads a CSV of unknown taxa, queries GBIF then Catalogue of Life (CoL)
   for each unique genus, and writes a new CSV with Sub-family, Family, and
   Order filled in. Empty Species cells are set to "sp" where a Genus exists.
   Returns a list of genera that could not be resolved by either source.

2. update_main_csv()
   Patches the main taxa CSV using a pre-filled source CSV, or by querying
   APIs directly for any rows that still have missing taxonomy columns.
   Returns a list of genera that could not be resolved.

Both functions only fill EMPTY cells by default (overwrite_existing=False),
so it is safe to rerun over a partially-filled file.
"""

import csv
import time
import requests
import shutil
from pathlib import Path

GBIF_MATCH_API = "https://api.gbif.org/v1/species/match"
GBIF_NAME_API  = "https://api.gbif.org/v1/species"
COL_SEARCH_API = "https://api.catalogueoflife.org/nameusage/search"

TARGET_RANKS = ("Sub-family", "Family", "Order")


# ─────────────────────────────────────────────────────────────────────────────
# Column-name resolver
# ─────────────────────────────────────────────────────────────────────────────

def _resolve_col(fieldnames: list, candidates: list):
    lower_map = {f.lower(): f for f in fieldnames}
    for c in candidates:
        hit = lower_map.get(c.lower())
        if hit:
            return hit
    return None


def _col(fieldnames, *candidates):
    col = _resolve_col(fieldnames, list(candidates))
    if col is None:
        raise ValueError(
            f"Could not find a column matching {candidates} in: {fieldnames}"
        )
    return col


# ─────────────────────────────────────────────────────────────────────────────
# GBIF
# ─────────────────────────────────────────────────────────────────────────────

def _query_gbif(genus: str, verbose: bool) -> dict:
    """Query GBIF species/match and parent-walk for Sub-family/Family/Order."""
    result = {rank: "" for rank in TARGET_RANKS}
    try:
        resp = requests.get(
            GBIF_MATCH_API,
            params={"name": genus, "rank": "GENUS", "strict": False},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        if verbose:
            print(f"    [GBIF WARN] {e}")
        return result

    if data.get("matchType", "NONE") == "NONE":
        return result

    gbif_key_map = {"Sub-family": "subfamily", "Family": "family", "Order": "order"}
    for col, key in gbif_key_map.items():
        val = data.get(key, "")
        if val:
            result[col] = val

    # Attempt parent-walk for subfamily if still missing
    if not result["Sub-family"] and "usageKey" in data:
        try:
            pr = requests.get(
                f"{GBIF_NAME_API}/{data['usageKey']}/parents", timeout=10
            )
            pr.raise_for_status()
            for parent in pr.json():
                if parent.get("rank", "").upper() == "SUBFAMILY":
                    result["Sub-family"] = parent.get("canonicalName", "")
                    break
        except requests.RequestException:
            pass

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Catalogue of Life
# ─────────────────────────────────────────────────────────────────────────────

def _query_col(genus: str, verbose: bool) -> dict:
    """
    Query the Catalogue of Life API for a genus and extract
    Sub-family, Family, and Order from the classification.
    """
    result = {rank: "" for rank in TARGET_RANKS}
    try:
        resp = requests.get(
            COL_SEARCH_API,
            params={"q": genus, "rank": "genus", "limit": 1},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        if verbose:
            print(f"    [CoL WARN] {e}")
        return result

    results = data.get("result", [])
    if not results:
        return result

    # Classification is a list of ancestor dicts with rank + name
    classification = results[0].get("classification", [])
    rank_map = {
        "subfamily": "Sub-family",
        "family":    "Family",
        "order":     "Order",
    }
    for entry in classification:
        rank_lower = entry.get("rank", "").lower()
        col_name   = rank_map.get(rank_lower)
        if col_name and not result[col_name]:
            result[col_name] = entry.get("name", "")

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Combined lookup: GBIF first, fill any gaps with CoL
# ─────────────────────────────────────────────────────────────────────────────

def query_genus(genus: str, verbose: bool = True, request_delay: float = 0.3) -> dict:
    """
    Try GBIF, then supplement missing ranks from Catalogue of Life.
    Returns a dict keyed by TARGET_RANKS.
    """
    genus = genus.strip()
    if not genus:
        return {rank: "" for rank in TARGET_RANKS}

    if verbose:
        print(f"  Querying GBIF...", end=" ", flush=True)
    result = _query_gbif(genus, verbose)

    missing = [r for r in TARGET_RANKS if not result[r]]
    if missing:
        if verbose:
            print(f"missing {missing}, trying CoL...", end=" ", flush=True)
        time.sleep(request_delay)
        col_result = _query_col(genus, verbose)
        for rank in missing:
            if col_result[rank]:
                result[rank] = col_result[rank]

    if verbose:
        filled = {k: v for k, v in result.items() if v}
        if filled:
            print(f"→ {filled}")
        else:
            print("→ no data found")

    return result


# ─────────────────────────────────────────────────────────────────────────────
# 1. fill_taxonomy_csv
# ─────────────────────────────────────────────────────────────────────────────

def fill_taxonomy_csv(
    input_path: str,
    output_path: str,
    overwrite_existing: bool = False,
    request_delay: float = 0.3,
    verbose: bool = True,
) -> list:
    """
    Read a CSV of unknown/new taxa, query GBIF + CoL for every unique genus,
    and write an enriched CSV with Sub-family, Family, Order filled in.
    Empty Species cells are set to 'sp' where a Genus exists.

    Safe to rerun over a partially-filled file — already-filled cells are
    skipped unless overwrite_existing=True.

    Returns
    -------
    List of genus names that could not be resolved by either API.
    """
    input_path  = Path(input_path)
    output_path = Path(output_path)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    with open(input_path, newline="", encoding="utf-8-sig") as f:
        reader     = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows       = list(reader)

    genus_col   = _col(fieldnames, "Genus",   "genus")
    species_col = _col(fieldnames, "Species", "species")

    for rank in TARGET_RANKS:
        if rank not in fieldnames:
            fieldnames.append(rank)

    # Only query genera that are actually missing at least one rank
    # (or all of them if overwrite_existing)
    def _needs_query(row):
        g = row.get(genus_col, "").strip()
        if not g:
            return False
        if overwrite_existing:
            return True
        return any(not row.get(rank, "").strip() for rank in TARGET_RANKS)

    unique_genera = {row[genus_col].strip() for row in rows if _needs_query(row)}
    cache = {}
    unresolved = []

    if verbose:
        print(f"Found {len(unique_genera)} genera needing lookup across {len(rows)} rows.")
        print("Querying GBIF + Catalogue of Life...\n")

    for i, genus in enumerate(sorted(unique_genera), 1):
        if verbose:
            print(f"[{i}/{len(unique_genera)}] {genus}")
        result = query_genus(genus, verbose=verbose, request_delay=request_delay)
        cache[genus] = result
        if not any(result.values()):
            unresolved.append(genus)
        time.sleep(request_delay)

    filled_count = 0
    sp_count     = 0

    for row in rows:
        genus = row.get(genus_col, "").strip()

        if genus and not row.get(species_col, "").strip():
            row[species_col] = "sp"
            sp_count += 1

        taxonomy = cache.get(genus, {})
        for rank in TARGET_RANKS:
            cell = row.get(rank, "").strip()
            if not cell or overwrite_existing:
                new_value = taxonomy.get(rank, "")
                if new_value:
                    row[rank] = new_value
                    filled_count += 1
                elif rank not in row:
                    row[rank] = ""

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if verbose:
        print(f"\nDone.")
        print(f"  Taxonomy cells filled : {filled_count}")
        print(f"  Species set to 'sp'   : {sp_count}")
        print(f"  Output                : {output_path}")
        if unresolved:
            print(f"\n  Could not resolve {len(unresolved)} genera (manual entry needed):")
            for g in unresolved:
                print(f"    - {g}")

    return unresolved


# ─────────────────────────────────────────────────────────────────────────────
# 2. update_main_csv
# ─────────────────────────────────────────────────────────────────────────────

def update_main_csv(
    main_path: str,
    source_path: str = None,
    backup: bool = True,
    overwrite_existing: bool = False,
    request_delay: float = 0.3,
    verbose: bool = True,
) -> list:
    """
    Update the main taxa CSV by filling rows with missing taxonomy columns.

    Modes
    -----
    source_path provided:
        Match rows by (Genus, Species) against the pre-filled source CSV.
        No API calls made.

    source_path omitted:
        Query GBIF + CoL directly for any genera still missing ranks.

    Safe to rerun — only empty cells are patched unless overwrite_existing=True.
    A .bak.csv backup is written before any changes.

    Returns
    -------
    List of genus names that could not be resolved (empty list in source mode).
    """
    main_path = Path(main_path)
    if not main_path.exists():
        raise FileNotFoundError(f"Main CSV not found: {main_path}")

    if backup:
        bak = main_path.with_suffix(".bak.csv")
        shutil.copy2(main_path, bak)
        if verbose:
            print(f"Backup written to: {bak}")

    with open(main_path, newline="", encoding="utf-8-sig") as f:
        reader     = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        main_rows  = list(reader)

    genus_col   = _col(fieldnames, "Genus",   "genus")
    species_col = _col(fieldnames, "Species", "species")

    for rank in TARGET_RANKS:
        if rank not in fieldnames:
            fieldnames.append(rank)

    patch_map  = {}   # (genus, species) → {rank: value}
    unresolved = []

    if source_path:
        src = Path(source_path)
        if not src.exists():
            raise FileNotFoundError(f"Source CSV not found: {src}")

        with open(src, newline="", encoding="utf-8-sig") as f:
            src_reader = csv.DictReader(f)
            src_fields = list(src_reader.fieldnames or [])
            src_rows   = list(src_reader)

        src_genus   = _col(src_fields, "Genus",   "genus")
        src_species = _col(src_fields, "Species", "species")

        for row in src_rows:
            g = row.get(src_genus,   "").strip()
            s = row.get(src_species, "").strip() or "sp"
            taxonomy = {
                rank: row.get(rank, "").strip()
                for rank in TARGET_RANKS
                if row.get(rank, "").strip()
            }
            if g and taxonomy:
                patch_map[(g, s)] = taxonomy
                patch_map.setdefault((g, "sp"), taxonomy)

        if verbose:
            print(f"Loaded {len(patch_map)} entries from source file.")

    else:
        # Find genera still missing at least one rank in the main CSV
        def _needs_query(row):
            g = row.get(genus_col, "").strip()
            if not g:
                return False
            if overwrite_existing:
                return True
            return any(not row.get(rank, "").strip() for rank in TARGET_RANKS)

        incomplete_genera = {
            row[genus_col].strip() for row in main_rows if _needs_query(row)
        }

        if verbose:
            print(
                f"Found {len(incomplete_genera)} genera with missing taxonomy. "
                f"Querying GBIF + CoL...\n"
            )

        gbif_cache = {}
        for i, genus in enumerate(sorted(incomplete_genera), 1):
            if verbose:
                print(f"[{i}/{len(incomplete_genera)}] {genus}")
            result = query_genus(genus, verbose=verbose, request_delay=request_delay)
            gbif_cache[genus] = result
            if not any(result.values()):
                unresolved.append(genus)
            time.sleep(request_delay)

        for row in main_rows:
            g = row.get(genus_col,   "").strip()
            s = row.get(species_col, "").strip() or "sp"
            if g and g in gbif_cache:
                patch_map[(g, s)] = gbif_cache[g]
                patch_map.setdefault((g, "sp"), gbif_cache[g])

    patched_cells = 0
    sp_count      = 0

    for row in main_rows:
        g = row.get(genus_col,   "").strip()
        s = row.get(species_col, "").strip()

        if g and not s:
            row[species_col] = "sp"
            s = "sp"
            sp_count += 1

        taxonomy = patch_map.get((g, s)) or patch_map.get((g, "sp"), {})

        for rank in TARGET_RANKS:
            cell = row.get(rank, "").strip()
            if not cell or overwrite_existing:
                new_value = taxonomy.get(rank, "")
                if new_value:
                    row[rank] = new_value
                    patched_cells += 1
                elif rank not in row:
                    row[rank] = ""

    with open(main_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(main_rows)

    if verbose:
        print(f"\nMain CSV updated in-place: {main_path}")
        print(f"  Taxonomy cells patched : {patched_cells}")
        print(f"  Species set to 'sp'    : {sp_count}")
        if unresolved:
            print(f"\n  Could not resolve {len(unresolved)} genera (manual entry needed):")
            for g in unresolved:
                print(f"    - {g}")

    return unresolved


# ─────────────────────────────────────────────────────────────────────────────
# 3. fallback_subfamily_from_family
# ─────────────────────────────────────────────────────────────────────────────

def fallback_subfamily_from_family(
    main_path: str,
    backup: bool = True,
    verbose: bool = True,
) -> int:
    """
    For any row in the main CSV where Sub-family is empty but Family is filled,
    derive a placeholder Sub-family by replacing the 'idae' suffix with 'inae'.

    If the Family name doesn't end in 'idae' it is used as-is, since the
    pattern doesn't apply cleanly.

    Safe to rerun — only touches rows where Sub-family is already empty.

    Parameters
    ----------
    main_path : Path to the main taxa CSV (edited in-place).
    backup    : Write a .bak.csv backup before making changes.
    verbose   : Print a summary to stdout.

    Returns
    -------
    Number of Sub-family cells filled.
    """
    main_path = Path(main_path)
    if not main_path.exists():
        raise FileNotFoundError(f"Main CSV not found: {main_path}")

    if backup:
        bak = main_path.with_suffix(".bak.csv")
        shutil.copy2(main_path, bak)
        if verbose:
            print(f"Backup written to: {bak}")

    with open(main_path, newline="", encoding="utf-8-sig") as f:
        reader     = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows       = list(reader)

    subfamily_col = _col(fieldnames, "Sub-family", "subfamily", "Subfamily")
    family_col    = _col(fieldnames, "Family",     "family")

    filled_count = 0

    for row in rows:
        if row.get(subfamily_col, "").strip():
            continue  # already filled — skip

        family = row.get(family_col, "").strip()
        if not family:
            continue  # nothing to derive from

        if family.lower().endswith("idae"):
            derived = family[:-4] + "inae"
        else:
            derived = family  # use as-is if suffix doesn't match

        row[subfamily_col] = derived
        filled_count += 1

    with open(main_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if verbose:
        print(f"Sub-family cells filled from Family name: {filled_count}")

    return filled_count



# ─────────────────────────────────────────────────────────────────────────────
# Edit paths below and run the file to execute
# ─────────────────────────────────────────────────────────────────────────────

#unresolved = fill_taxonomy_csv(
#    input_path="C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/unknown_taxa.csv",
#    output_path="C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/unknown_taxa_filled.csv",
#)

#unresolved += update_main_csv(
#    main_path="C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv",
#    source_path="C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/unknown_taxa_filled.csv",
#)

#if unresolved:
#    print(f"\nFull list of unresolved genera requiring manual input:")
#    for g in sorted(set(unresolved)):
#        print(f"  {g}")

#fallback_subfamily_from_family(main_path="C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")

    # Step 1: enrich the unknown taxa CSV via GBIF


    # Step 2: patch the main CSV using the filled file

#run this: python -m bioinf_packages.verify_funcs.fill_taxonomy

#do this if you have pre-filled unkown_taxa.csv
#python fill_taxonomy.py update --main "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv" --source "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/unknown_taxa.csv"
#do this if you don't have it (slower)
#python fill_taxonomy.py update --main "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"

#or both if you don't have unknown taxa filled.
#python -m bioinf_packages.verify_funcs.fill_taxonomy fill-and-update --input "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/unknown_taxa.csv" --main "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"