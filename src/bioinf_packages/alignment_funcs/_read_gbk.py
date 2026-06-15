import os
from Bio import SeqIO
import csv
import re

from ..dictionary_funcs.correction_finder import correction_finder


# ─────────────────────────────────────────────────────────────────────────────
# TAXON NAME NORMALISATION
# ─────────────────────────────────────────────────────────────────────────────

def normalise_taxon_name(taxon: str) -> str:
    """
    Normalise a taxon name to Genus_species, stripping subspecies
    and standardising separators.

    Examples
    --------
    Solenodon_paradoxus_paradoxus → Solenodon_paradoxus
    Solenodon paradoxus paradoxus → Solenodon_paradoxus
    Solenodon_paradoxus           → Solenodon_paradoxus  (unchanged)
    Atopogale_cubana              → Solenodon_cubana
    Atopogale cubana              → Solenodon_cubana

    Parameters
    ----------
    taxon : str   taxon name in any reasonable format

    Returns
    -------
    str   normalised Genus_species string
    """
    # Normalise separators — replace spaces and hyphens with underscores
    normalised = taxon.strip().replace(" ", "_").replace("-", "_")

    # Split and keep only genus and species (drop subspecies and beyond)
    parts = normalised.split("_")
    parts = [p for p in parts if p]  # remove empty tokens

    if len(parts) >= 2:
        genus   = parts[0]
        species = parts[1]
    elif len(parts) == 1:
        genus   = parts[0]
        species = None
    else:
        return taxon  # cannot parse — return unchanged

    # Normalise genus synonyms
    # Atopogale is the newer name for Cuban solenodon but many GenBank
    # records still use Solenodon cubanus — treat as the same genus
    genus_synonyms = {
        "Atopogale": "Solenodon",
    }
    genus = genus_synonyms.get(genus, genus)

    if species:
        return f"{genus}_{species}"
    else:
        return genus


# ─────────────────────────────────────────────────────────────────────────────
# GBK AUDIT
# ─────────────────────────────────────────────────────────────────────────────

def audit_gbk_features(data_folder: str,
                        min_length: int = 0,
                        save_path: str = None) -> dict:
    """
    Scan all .gbk files in a directory and report every feature type,
    qualifier key, symbol value, and sequence length found.

    Use this to determine which feature types and qualifier keys to
    include in read_gbk() for your specific dataset.

    Parameters
    ----------
    data_folder : str   directory containing .gbk files
    min_length  : int   only report features longer than this (default 0)
    save_path   : str   optional path to save report as CSV

    Returns
    -------
    dict with keys:
        feature_summary : dict   feature_type -> list of entry dicts
        type_counts     : dict   feature_type -> count across all files
        short_features  : list   features under 50bp
        long_features   : list   features 50bp and over
    """
    from collections import defaultdict

    feature_summary = defaultdict(list)
    type_counts     = defaultdict(int)
    short_features  = []
    long_features   = []

    gbk_files = [f for f in os.listdir(data_folder)
                 if f.endswith(".gbk") or f.endswith(".gb")]

    print(f"Found {len(gbk_files)} GBK files in {data_folder}\n")

    for filename in sorted(gbk_files):
        path      = os.path.join(data_folder, filename)
        accession = os.path.splitext(filename)[0]

        try:
            for record in SeqIO.parse(path, "genbank"):
                species = record.annotations.get("organism", "Unknown")

                for feature in record.features:

                    if feature.type in ("source",):
                        continue

                    qualifiers = feature.qualifiers
                    length     = (int(feature.location.end) -
                                  int(feature.location.start))

                    if length < min_length:
                        continue

                    symbol = (
                        qualifiers.get("standard_name", [None])[0] or
                        qualifiers.get("gene",          [None])[0] or
                        qualifiers.get("product",       [None])[0] or
                        qualifiers.get("note",          [None])[0] or
                        "NO_SYMBOL"
                    )

                    entry = {
                        "accession":    accession,
                        "species":      species,
                        "feature_type": feature.type,
                        "symbol":       symbol,
                        "length_bp":    length,
                        "start":        int(feature.location.start) + 1,
                        "end":          int(feature.location.end),
                        "orientation":  ("+"
                                         if feature.location.strand == 1
                                         else "-"),
                        "qualifiers":   list(qualifiers.keys()),
                    }

                    feature_summary[feature.type].append(entry)
                    type_counts[feature.type] += 1

                    if length < 50:
                        short_features.append(entry)
                    else:
                        long_features.append(entry)

        except Exception as e:
            print(f"  ✗ Failed to parse {filename}: {e}")

    # ── Print summary ─────────────────────────────────────────────────────────
    print("── Feature type counts ──────────────────────────────────────")
    for ftype, count in sorted(type_counts.items(),
                                key=lambda x: x[1], reverse=True):
        entries     = feature_summary[ftype]
        lengths     = [e["length_bp"] for e in entries]
        mean_len    = sum(lengths) / len(lengths) if lengths else 0
        unique_syms = set(e["symbol"] for e in entries)
        print(f"\n  {ftype:<20} count={count:<6} "
              f"mean_length={mean_len:.0f}bp")
        print(f"  {'':20} symbols: "
              f"{', '.join(sorted(unique_syms))}")

    print(f"\n── Length distribution ──────────────────────────────────────")
    print(f"  Features < 50bp  : {len(short_features)}")
    print(f"  Features >= 50bp : {len(long_features)}")

    if save_path:
        all_entries = [e for entries in feature_summary.values()
                       for e in entries]
        for e in all_entries:
            e["qualifiers"] = ";".join(e["qualifiers"])

        fieldnames = ["accession", "species", "feature_type",
                      "symbol", "length_bp", "start", "end",
                      "orientation", "qualifiers"]

        with open(save_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_entries)

        print(f"\nReport saved: {save_path}")

    return {
        "feature_summary": dict(feature_summary),
        "type_counts":     dict(type_counts),
        "short_features":  short_features,
        "long_features":   long_features,
    }


# ─────────────────────────────────────────────────────────────────────────────
# READ GBK
# ─────────────────────────────────────────────────────────────────────────────

def read_gbk(accession: str,
             data_folder: str,
             glossary: str,
             min_length: int = 100):
    """
    Read a GenBank (.gbk) file and extract coding genes.

    Feature type priority: CDS > rRNA > gene
    Excludes tRNA, mRNA, exon, misc_feature, D-loop, rep_origin,
    regulatory, intron, gap, repeat_region, and all other non-gene
    feature types.

    Species names are normalised via normalise_taxon_name() which:
      - strips subspecies tokens (Solenodon_paradoxus_paradoxus
        → Solenodon_paradoxus)
      - normalises Atopogale → Solenodon
      - replaces spaces and hyphens with underscores

    Parameters
    ----------
    accession   : str   GenBank accession number
    data_folder : str   directory containing .gbk files
    glossary    : str   path to gene symbol glossary CSV
    min_length  : int   minimum feature length in bp (default 100)

    Returns
    -------
    list of dicts with keys:
        symbol, begin, end, orientation, species, accession
    """
    with open(glossary, newline='', encoding='utf-8') as f:
        reader     = csv.DictReader(f)
        symbol_key = [row for row in reader][0]

    gbk_path = os.path.join(data_folder, f"{accession}.gbk")
    if not os.path.exists(gbk_path):
        raise FileNotFoundError(
            f"No .gbk file found at: {gbk_path}"
        )

    FEATURE_PRIORITY = {"CDS": 0, "rRNA": 1, "gene": 2}

    gene_list = []

    for record in SeqIO.parse(gbk_path, "genbank"):

        # ── Normalise species name ────────────────────────────────────────────
        # Strip subspecies tokens and normalise genus synonyms.
        # This ensures "Solenodon paradoxus paradoxus" becomes
        # "Solenodon_paradoxus" and "Atopogale cubana" becomes
        # "Solenodon_cubana" consistently across all downstream
        # functions.
        raw_species = record.annotations.get("organism", "Unknown")
        species     = normalise_taxon_name(raw_species)

        candidates = []

        for feature in record.features:

            if feature.type not in FEATURE_PRIORITY:
                continue

            qualifiers = feature.qualifiers

            symbol = (
                qualifiers.get("standard_name", [None])[0] or
                qualifiers.get("gene",          [None])[0] or
                qualifiers.get("product",       [None])[0]
            )

            if not symbol:
                continue
            if "trna" in symbol.lower() or symbol.startswith("trn"):
                continue
            if symbol == "NO_SYMBOL":
                continue

            feat_length = (int(feature.location.end) -
                           int(feature.location.start))
            if feat_length < min_length:
                continue

            try:
                symbol_correct = correction_finder(symbol, symbol_key,
                                                   glossary)
            except Exception:
                symbol_correct = symbol

            candidates.append({
                "priority":    FEATURE_PRIORITY[feature.type],
                "symbol":      symbol_correct,
                "begin":       int(feature.location.start) + 1,
                "end":         int(feature.location.end),
                "orientation": ("+"
                                if feature.location.strand == 1
                                else "-"),
                "species":     species,   # ← normalised name used here
                "accession":   accession,
            })

        # ── Deduplicate by symbol ─────────────────────────────────────────────
        # For each unique gene keep only the highest-priority feature
        # type (CDS > rRNA > gene).
        best = {}
        for c in candidates:
            key = c["symbol"].upper()
            if key not in best or c["priority"] < best[key]["priority"]:
                best[key] = c

        for c in best.values():
            gene_list.append({k: v for k, v in c.items()
                               if k != "priority"})

    return gene_list

#print(read_gbk("OZ208999.1", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info"))
#print(read_gbk("KC516842.1", "fasta_info"))
#print(read_gbk("NC_010298", "fasta_info"))
#print(read_gbk("KR710721", "fasta_info"))