import os
from Bio import SeqIO
import csv

#import symbol key and symbol formatting function

from ..dictionary_funcs.correction_finder import correction_finder

#function to find all possible labels used in .gbk files to determine which ones to use.
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
    min_length  : int   only report features longer than this (default 0 = all)
    save_path   : str   optional path to save report as CSV

    Returns
    -------
    dict with keys:
        "feature_summary" : dict   feature_type -> list of dicts with
                                   qualifier values and lengths
        "type_counts"     : dict   feature_type -> count across all files
        "short_features"  : list   features under 50bp (likely spurious)
        "long_features"   : list   features over 50bp (likely genuine genes)
    """
    import os
    import csv
    from Bio import SeqIO
    from collections import defaultdict

    feature_summary = defaultdict(list)
    type_counts     = defaultdict(int)
    short_features  = []
    long_features   = []

    gbk_files = [f for f in os.listdir(data_folder)
                 if f.endswith(".gbk") or f.endswith(".gb")]

    print(f"Found {len(gbk_files)} GBK files in {data_folder}\n")

    for filename in sorted(gbk_files):
        path    = os.path.join(data_folder, filename)
        accession = os.path.splitext(filename)[0]

        try:
            for record in SeqIO.parse(path, "genbank"):
                species = record.annotations.get("organism", "Unknown")

                for feature in record.features:

                    # Skip uninformative top-level features
                    if feature.type in ("source",):
                        continue

                    qualifiers = feature.qualifiers
                    length     = int(feature.location.end) - \
                                 int(feature.location.start)

                    if length < min_length:
                        continue

                    # Extract every possible name qualifier
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
                        "orientation":  "+" if feature.location.strand == 1
                                        else "-",
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
    print("── Feature type counts across all GBK files ──────────────────")
    for ftype, count in sorted(type_counts.items(),
                                key=lambda x: x[1], reverse=True):
        entries      = feature_summary[ftype]
        lengths      = [e["length_bp"] for e in entries]
        mean_len     = sum(lengths) / len(lengths) if lengths else 0
        unique_syms  = set(e["symbol"] for e in entries)
        print(f"\n  {ftype:<20} count={count:<6} "
              f"mean_length={mean_len:.0f}bp")
        print(f"  {'':20} symbols: {', '.join(sorted(unique_syms))}")

    print(f"\n── Length distribution ───────────────────────────────────────")
    print(f"  Features < 50bp  (likely spurious) : {len(short_features)}")
    print(f"  Features >= 50bp (likely genuine)  : {len(long_features)}")

    print(f"\n── Features >= 50bp by type ──────────────────────────────────")
    long_by_type = defaultdict(list)
    for e in long_features:
        long_by_type[e["feature_type"]].append(e)

    for ftype, entries in sorted(long_by_type.items()):
        unique_syms = sorted(set(e["symbol"] for e in entries))
        print(f"\n  {ftype} ({len(entries)} features):")
        for sym in unique_syms:
            sym_entries = [e for e in entries if e["symbol"] == sym]
            lengths     = [e["length_bp"] for e in sym_entries]
            print(f"    {sym:<30} n={len(sym_entries):<4} "
                  f"lengths: {min(lengths)}-{max(lengths)}bp")

    # ── Save to CSV if requested ──────────────────────────────────────────────
    if save_path:
        all_entries = [e for entries in feature_summary.values()
                       for e in entries]

        # Flatten qualifiers list to string for CSV
        for e in all_entries:
            e["qualifiers"] = ";".join(e["qualifiers"])

        fieldnames = ["accession", "species", "feature_type", "symbol",
                      "length_bp", "start", "end", "orientation",
                      "qualifiers"]

        with open(save_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_entries)

        print(f"\nFull report saved to: {save_path}")

    return {
        "feature_summary": dict(feature_summary),
        "type_counts":     dict(type_counts),
        "short_features":  short_features,
        "long_features":   long_features,
    }


#usage
#report = audit_gbk_features(
#    data_folder = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
#    min_length  = 0,          # report everything
#    save_path   = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/gbk_audit.csv"   # optional — saves full report to CSV
#)

# Then inspect specific feature types interactively:
#for entry in report["feature_summary"]["misc_feature"]:
#    if entry["length_bp"] >= 50:
#        print(f"{entry['accession']} | {entry['symbol']} | "
#              f"{entry['length_bp']}bp")


#read data file in gbk format
def read_gbk(accession: str, data_folder: str, glossary: str,
             min_length: int = 100):
    """
    Reads a GenBank (.gbk) file and extracts coding genes.

    Feature type priority: CDS > rRNA > gene
    Excludes: tRNA, mRNA, exon, misc_feature, D-loop, rep_origin,
              regulatory, intron, gap, repeat_region, and all other
              non-gene feature types.

    Parameters
    ----------
    accession   : str   GenBank accession number
    data_folder : str   directory containing .gbk files
    glossary    : str   path to gene symbol glossary CSV
    min_length  : int   minimum feature length in bp (default 100)
                        filters tRNA-sized fragments and sub-gene annotations

    Returns
    -------
    list of dicts with keys: symbol, begin, end, orientation,
                             species, accession
    """
    from ..dictionary_funcs.correction_finder import correction_finder

    with open(glossary, newline='', encoding='utf-8') as f:
        reader     = csv.DictReader(f)
        symbol_key = [row for row in reader][0]

    gbk_path = os.path.join(data_folder, f"{accession}.gbk")
    if not os.path.exists(gbk_path):
        raise FileNotFoundError(f"No .gbk file found at: {gbk_path}")

    # Priority order for feature types — CDS is most reliable,
    # rRNA is needed for 12S/16S, gene is fallback only
    FEATURE_PRIORITY = {"CDS": 0, "rRNA": 1, "gene": 2}

    gene_list = []

    for record in SeqIO.parse(gbk_path, "genbank"):
        species = record.annotations.get("organism", "Unknown")

        # Collect all valid features first, then deduplicate
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

            feat_length = int(feature.location.end) - \
                          int(feature.location.start)
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
                "orientation": "+" if feature.location.strand == 1
                               else "-",
                "species":     species.strip()
                                      .replace(" ", "_")
                                      .replace("-", "_"),
                "accession":   accession,
            })

        # ── Deduplicate by (symbol, accession) ───────────────────────
        # For each unique gene in this accession, keep only the
        # highest-priority feature type (CDS > rRNA > gene).
        # This handles cases where gene and CDS features exist for
        # the same locus with slightly different coordinates.
        best = {}   # symbol_upper -> candidate dict
        for c in candidates:
            key = c["symbol"].upper()
            if key not in best or c["priority"] < best[key]["priority"]:
                best[key] = c

        # Strip the priority field before adding to output
        for c in best.values():
            gene_list.append({k: v for k, v in c.items()
                               if k != "priority"})

    return gene_list

#print(read_gbk("OZ208999.1", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info"))
#print(read_gbk("KC516842.1", "fasta_info"))
#print(read_gbk("NC_010298", "fasta_info"))
#print(read_gbk("KR710721", "fasta_info"))