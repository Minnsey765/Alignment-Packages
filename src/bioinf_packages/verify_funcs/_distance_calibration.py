import os
import csv
import random
import numpy as np

from pathlib import Path
from Bio import SeqIO

def calibrate_distance_thresholds(
        genuine_fasta: str,
        artefact_fasta: str,
        scaffold_dir: str,
        taxonomy: dict,
        insert_fraction: float = 0.30,
        target_detection_prob: float = 0.90,
        k_nearest: int = 3,
        scaffold_extensions: tuple = (".fasta", ".fa", ".fas")
) -> dict:
    """
    Calibrate outlier_ratio_threshold and min_ratio_variance_flag for
    verify_by_distance() using known genuine and artefactual sequences.

    Runs verify_by_distance() on every sequence in both FASTA files,
    automatically matching each sequence to its corresponding scaffold
    FASTA by gene name extracted from the sequence header. Computes the
    distribution of mean_similarity_ratio, knn_similarity_ratio and
    ratio_variance for genuine and artefactual sequences separately,
    then recommends thresholds that maximise separation between the two
    groups.

    FASTA header formats accepted
    ------------------------------
    Genuine  : >Genus_species|gene|orien:+|accession:XXXX
    Artefact : >Genus_species|gene|mod_N|orien:+|accession:XXXX

    The gene field (second pipe-delimited token) is used to find the
    corresponding scaffold FASTA in scaffold_dir. For example a header
    containing |12S_rRNA| will look for 12S_rRNA.fasta in scaffold_dir.

    Parameters
    ----------
    genuine_fasta        : str    path to FASTA of known genuine sequences
    artefact_fasta       : str    path to FASTA of known artefactual
                                  sequences
    scaffold_dir         : str    directory containing per-gene scaffold
                                  FASTA files named {gene}.fasta
    taxonomy             : dict   from load_taxonomy()
    insert_fraction      : float  passed to verify_by_distance()
    target_detection_prob: float  passed to verify_by_distance()
    k_nearest            : int    passed to verify_by_distance()
    scaffold_extensions  : tuple  file extensions to search for scaffold
                                  files (default .fasta, .fa, .fas)

    Returns
    -------
    dict with:
        genuine_results     : list   raw verify_by_distance() outputs
        artefact_results    : list   raw verify_by_distance() outputs
        genuine_stats       : dict   mean/std/min/max of each metric
                                     across genuine sequences
        artefact_stats      : dict   mean/std/min/max of each metric
                                     across artefact sequences
        recommended_thresholds : dict
            outlier_ratio_threshold  : float
            min_ratio_variance_flag  : float
        per_gene_stats      : dict   stats broken down by gene
        summary_df_path     : str    path to saved CSV summary
    """
    import os
    import csv
    from pathlib import Path
    from collections import defaultdict

    if isinstance(taxonomy, str):
        raise TypeError(
            f"taxonomy must be a dict from load_taxonomy(), not a path."
        )

    scaffold_dir_path = Path(scaffold_dir)

    # ── Helper: find scaffold FASTA for a gene name ───────────────────────────

    def find_scaffold(gene_name: str):
        """
        Find scaffold FASTA for a gene name in scaffold_dir.
        Tries exact match and case-insensitive match.
        Returns path string or None.
        """
        for ext in scaffold_extensions:
            # Exact match
            candidate = scaffold_dir_path / f"{gene_name}{ext}"
            if candidate.exists():
                return str(candidate)

            # Case-insensitive match
            for f in scaffold_dir_path.iterdir():
                if (f.stem.lower() == gene_name.lower()
                        and f.suffix.lower() == ext.lower()):
                    return str(f)

        return None

    # ── Helper: parse header into components ──────────────────────────────────

    def parse_header(description: str) -> dict:
        """
        Parse pipe-delimited FASTA header.

        Genuine : Genus_species|gene|orien:+|accession:XXXX
        Artefact: Genus_species|gene|mod_N|orien:+|accession:XXXX

        Returns dict with taxon, gene, accession, is_modified.
        """
        parts       = description.strip().split("|")
        taxon       = parts[0].strip()
        gene        = parts[1].strip() if len(parts) > 1 else None
        accession   = None
        is_modified = False

        for p in parts:
            if p.startswith("accession:"):
                accession = p.split(":", 1)[1].strip()
            if p.startswith("mod_"):
                is_modified = True

        # Also check for mod_ in gene field position
        if gene and gene.startswith("mod_"):
            is_modified = True
            # In this case the real gene is in parts[1] before mod
            gene = parts[1].strip()

        return {
            "taxon":       taxon,
            "gene":        gene,
            "accession":   accession,
            "is_modified": is_modified,
        }

    # ── Helper: run verify_by_distance on one FASTA file ─────────────────────

    def run_on_fasta(fasta_path: str, label: str) -> list:
        """
        Run verify_by_distance() on every sequence in fasta_path.
        Skips sequences where no scaffold can be found.
        Returns list of result dicts with added fields:
            label, gene, scaffold_used
        """
        records = list(SeqIO.parse(fasta_path, "fasta"))
        results = []
        skipped = []

        print(f"\n{'='*60}")
        print(f"Processing {label} sequences: {len(records)} total")
        print(f"{'='*60}")

        for i, record in enumerate(records, 1):
            info = parse_header(record.description)

            if not info["gene"]:
                print(f"  [{i}/{len(records)}] Skipping "
                      f"{record.description[:50]} — no gene in header")
                skipped.append(record.description)
                continue

            if not info["accession"]:
                print(f"  [{i}/{len(records)}] Skipping "
                      f"{info['taxon']} — no accession in header")
                skipped.append(record.description)
                continue

            scaffold_path = find_scaffold(info["gene"])
            if scaffold_path is None:
                print(f"  [{i}/{len(records)}] Skipping "
                      f"{info['taxon']} | {info['gene']} — "
                      f"no scaffold found in {scaffold_dir}")
                skipped.append(record.description)
                continue

            print(f"\n[{label}] {i}/{len(records)} | "
                  f"{info['taxon']} | {info['gene']} | "
                  f"{info['accession']}")
            print(f"  Scaffold: {Path(scaffold_path).name}")

            try:
                result = verify_by_distance(
                    query_seq             = str(record.seq),
                    query_taxon           = info["taxon"],
                    query_accession       = info["accession"],
                    scaffold_fasta        = scaffold_path,
                    taxonomy              = taxonomy,
                    insert_fraction       = insert_fraction,
                    target_detection_prob = target_detection_prob,
                    # Use permissive thresholds during calibration —
                    # we want raw scores, not filtered results
                    outlier_ratio_threshold = 0.0,
                    min_ratio_variance_flag  = 0.0,
                    k_nearest               = k_nearest,
                )
                result["label"]         = label
                result["gene"]          = info["gene"]
                result["scaffold_used"] = scaffold_path
                result["header"]        = record.description
                results.append(result)

            except Exception as e:
                print(f"  ✗ Failed: {e}")
                import traceback
                traceback.print_exc()
                skipped.append(record.description)

        print(f"\n{label} complete: "
              f"{len(results)}/{len(records)} succeeded, "
              f"{len(skipped)} skipped")
        if skipped:
            print(f"  Skipped:")
            for s in skipped:
                print(f"    {s[:70]}")

        return results

    # ── Run on both files ─────────────────────────────────────────────────────

    genuine_results  = run_on_fasta(genuine_fasta,  "GENUINE")
    artefact_results = run_on_fasta(artefact_fasta, "ARTEFACT")

    if not genuine_results:
        raise ValueError("No genuine sequences were successfully scored.")
    if not artefact_results:
        raise ValueError("No artefact sequences were successfully scored.")

    # ── Compute statistics ────────────────────────────────────────────────────

    def compute_stats(results: list, label: str) -> dict:
        """
        Compute mean, std, min, max for key metrics across a set of
        verify_by_distance() results.
        """
        metrics = [
            "mean_similarity_ratio",
            "mean_knn_similarity_ratio",
            "min_mean_ratio",
            "min_knn_ratio",
            "mean_ratio_variance",
            "knn_ratio_variance",
            "mean_ingroup_identity",
            "mean_outgroup_identity",
            "prop_outliers",
        ]

        stats = {"label": label, "n": len(results)}
        for m in metrics:
            values = [r[m] for r in results
                      if r.get(m) is not None]
            if values:
                stats[f"{m}_mean"] = round(float(np.mean(values)), 4)
                stats[f"{m}_std"]  = round(float(np.std(values)),  4)
                stats[f"{m}_min"]  = round(float(np.min(values)),  4)
                stats[f"{m}_max"]  = round(float(np.max(values)),  4)
            else:
                stats[f"{m}_mean"] = None
                stats[f"{m}_std"]  = None
                stats[f"{m}_min"]  = None
                stats[f"{m}_max"]  = None

        return stats

    genuine_stats  = compute_stats(genuine_results,  "genuine")
    artefact_stats = compute_stats(artefact_results, "artefact")

    # ── Per-gene stats ────────────────────────────────────────────────────────

    def compute_per_gene_stats(genuine: list,
                                artefact: list) -> dict:
        """
        Break down statistics by gene so you can see which genes
        provide the best discrimination between genuine and artefact.
        """
        from collections import defaultdict

        gene_groups = defaultdict(lambda: {"genuine": [], "artefact": []})

        for r in genuine:
            gene_groups[r["gene"]]["genuine"].append(r)
        for r in artefact:
            gene_groups[r["gene"]]["artefact"].append(r)

        per_gene = {}
        for gene, groups in gene_groups.items():
            g_ratios = [r["mean_similarity_ratio"]
                        for r in groups["genuine"]
                        if r.get("mean_similarity_ratio") is not None]
            a_ratios = [r["mean_similarity_ratio"]
                        for r in groups["artefact"]
                        if r.get("mean_similarity_ratio") is not None]

            g_knn = [r["mean_knn_similarity_ratio"]
                     for r in groups["genuine"]
                     if r.get("mean_knn_similarity_ratio") is not None]
            a_knn = [r["mean_knn_similarity_ratio"]
                     for r in groups["artefact"]
                     if r.get("mean_knn_similarity_ratio") is not None]

            per_gene[gene] = {
                "n_genuine":            len(groups["genuine"]),
                "n_artefact":           len(groups["artefact"]),
                "genuine_mean_ratio":   round(np.mean(g_ratios), 4)
                                        if g_ratios else None,
                "artefact_mean_ratio":  round(np.mean(a_ratios), 4)
                                        if a_ratios else None,
                "genuine_knn_ratio":    round(np.mean(g_knn), 4)
                                        if g_knn else None,
                "artefact_knn_ratio":   round(np.mean(a_knn), 4)
                                        if a_knn else None,
                # Separation: how far apart are genuine and artefact
                # means? Larger = better discrimination for this gene
                "mean_ratio_separation": round(
                    np.mean(g_ratios) - np.mean(a_ratios), 4
                ) if g_ratios and a_ratios else None,
                "knn_ratio_separation": round(
                    np.mean(g_knn) - np.mean(a_knn), 4
                ) if g_knn and a_knn else None,
            }

        return per_gene

    per_gene_stats = compute_per_gene_stats(genuine_results,
                                             artefact_results)

    # ── Recommend thresholds ──────────────────────────────────────────────────
    # Threshold recommendation strategy:
    #
    # For outlier_ratio_threshold:
    #   Set at the midpoint between the minimum genuine mean_ratio and
    #   the maximum artefact mean_ratio. This maximises the gap between
    #   the two distributions.
    #   If there is no clean separation, set conservatively at the
    #   5th percentile of genuine scores (5% false positive rate).
    #
    # For min_ratio_variance_flag:
    #   Set at the midpoint between the maximum genuine variance and
    #   the minimum artefact variance.

    g_mean_ratios = [r["mean_similarity_ratio"]
                     for r in genuine_results
                     if r.get("mean_similarity_ratio") is not None]
    a_mean_ratios = [r["mean_similarity_ratio"]
                     for r in artefact_results
                     if r.get("mean_similarity_ratio") is not None]
    g_variances   = [r["mean_ratio_variance"]
                     for r in genuine_results
                     if r.get("mean_ratio_variance") is not None]
    a_variances   = [r["mean_ratio_variance"]
                     for r in artefact_results
                     if r.get("mean_ratio_variance") is not None]

    # Ratio threshold — midpoint between distributions if separated,
    # otherwise 5th percentile of genuine (conservative)
    g_min  = float(np.min(g_mean_ratios)) if g_mean_ratios else None
    a_max  = float(np.max(a_mean_ratios)) if a_mean_ratios else None
    g_5th  = float(np.percentile(g_mean_ratios, 5)) \
             if g_mean_ratios else None

    if g_min is not None and a_max is not None and g_min > a_max:
        # Clean separation — use midpoint
        ratio_threshold = round((g_min + a_max) / 2, 4)
        separation_note = (
            f"Clean separation: genuine min={g_min:.4f}, "
            f"artefact max={a_max:.4f}, "
            f"midpoint={ratio_threshold:.4f}"
        )
    else:
        # Overlapping distributions — use 5th percentile of genuine
        ratio_threshold = round(g_5th, 4) if g_5th else 0.95
        separation_note = (
            f"Overlapping distributions: using 5th percentile of "
            f"genuine scores ({ratio_threshold:.4f}). "
            f"Genuine min={g_min:.4f}, artefact max={a_max:.4f}. "
            f"Consider using a more variable gene for better "
            f"discrimination."
        )

    # Variance threshold
    g_var_max = float(np.max(g_variances)) if g_variances else None
    a_var_min = float(np.min(a_variances)) if a_variances else None

    if (g_var_max is not None and a_var_min is not None
            and a_var_min > g_var_max):
        variance_threshold = round((g_var_max + a_var_min) / 2, 6)
        var_note = (
            f"Genuine max var={g_var_max:.6f}, "
            f"artefact min var={a_var_min:.6f}, "
            f"midpoint={variance_threshold:.6f}"
        )
    else:
        variance_threshold = 0.001
        var_note = (
            f"Overlapping variance distributions — using default 0.001. "
            f"Variance alone may not discriminate well for this dataset."
        )

    recommended_thresholds = {
        "outlier_ratio_threshold":  ratio_threshold,
        "min_ratio_variance_flag":  variance_threshold,
        "separation_note":          separation_note,
        "variance_note":            var_note,
    }

    # ── Print calibration report ──────────────────────────────────────────────

    print(f"\n{'='*60}")
    print(f"CALIBRATION REPORT")
    print(f"{'='*60}")
    print(f"\nGenuine sequences  : {len(genuine_results)}")
    print(f"Artefact sequences : {len(artefact_results)}")

    print(f"\n── Mean similarity ratio ─────────────────────────────────")
    print(f"  Genuine  : mean={genuine_stats['mean_similarity_ratio_mean']:.4f} "
          f"± {genuine_stats['mean_similarity_ratio_std']:.4f} "
          f"(min={genuine_stats['mean_similarity_ratio_min']:.4f}, "
          f"max={genuine_stats['mean_similarity_ratio_max']:.4f})")
    print(f"  Artefact : mean={artefact_stats['mean_similarity_ratio_mean']:.4f} "
          f"± {artefact_stats['mean_similarity_ratio_std']:.4f} "
          f"(min={artefact_stats['mean_similarity_ratio_min']:.4f}, "
          f"max={artefact_stats['mean_similarity_ratio_max']:.4f})")

    print(f"\n── KNN similarity ratio ──────────────────────────────────")
    print(f"  Genuine  : mean={genuine_stats['mean_knn_similarity_ratio_mean']:.4f} "
          f"± {genuine_stats['mean_knn_similarity_ratio_std']:.4f}")
    print(f"  Artefact : mean={artefact_stats['mean_knn_similarity_ratio_mean']:.4f} "
          f"± {artefact_stats['mean_knn_similarity_ratio_std']:.4f}")

    print(f"\n── Ratio variance ────────────────────────────────────────")
    print(f"  Genuine  : mean={genuine_stats['mean_ratio_variance_mean']:.6f} "
          f"± {genuine_stats['mean_ratio_variance_std']:.6f}")
    print(f"  Artefact : mean={artefact_stats['mean_ratio_variance_mean']:.6f} "
          f"± {artefact_stats['mean_ratio_variance_std']:.6f}")

    print(f"\n── Per-gene separation (genuine - artefact mean ratio) ───")
    for gene, stats in sorted(per_gene_stats.items(),
                               key=lambda x: -(x[1]["mean_ratio_separation"]
                                               or 0)):
        sep = stats["mean_ratio_separation"]
        print(f"  {gene:<20} separation={sep:.4f}"
              f"  (genuine={stats['genuine_mean_ratio']}, "
              f"artefact={stats['artefact_mean_ratio']}, "
              f"n={stats['n_genuine']}g/{stats['n_artefact']}a)"
              if sep is not None else
              f"  {gene:<20} insufficient data")

    print(f"\n── Recommended thresholds ────────────────────────────────")
    print(f"  outlier_ratio_threshold : {ratio_threshold}")
    print(f"  min_ratio_variance_flag : {variance_threshold}")
    print(f"\n  Ratio threshold rationale:")
    print(f"    {separation_note}")
    print(f"\n  Variance threshold rationale:")
    print(f"    {var_note}")

    print(f"\n  Use these in verify_by_distance():")
    print(f"    verify_by_distance(")
    print(f"        ...,")
    print(f"        outlier_ratio_threshold = {ratio_threshold},")
    print(f"        min_ratio_variance_flag = {variance_threshold},")
    print(f"    )")

    # ── Save summary CSV ──────────────────────────────────────────────────────

    output_dir   = Path(scaffold_dir).parent / "calibration_output"
    output_dir.mkdir(exist_ok=True)
    summary_path = str(output_dir / "calibration_summary.csv")

    all_results = genuine_results + artefact_results
    csv_keys    = [
        "label", "gene", "query_taxon", "query_accession",
        "expected_order", "n_samples", "subsample_length",
        "mean_similarity_ratio", "mean_knn_similarity_ratio",
        "min_mean_ratio", "min_knn_ratio",
        "mean_ratio_variance", "knn_ratio_variance",
        "mean_ingroup_identity", "mean_outgroup_identity",
        "chimera_flag", "scaffold_used",
    ]

    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_keys,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_results)

    print(f"\nCalibration CSV saved: {summary_path}")

    return {
        "genuine_results":        genuine_results,
        "artefact_results":       artefact_results,
        "genuine_stats":          genuine_stats,
        "artefact_stats":         artefact_stats,
        "per_gene_stats":         per_gene_stats,
        "recommended_thresholds": recommended_thresholds,
        "summary_csv":            summary_path,
    }

def validate_chimera_detection(
        genuine_fasta: str,
        artefact_fasta: str,
        scaffold_dir: str,
        taxonomy: dict,
        output_dir: str,
        insert_fraction: float = 0.30,
        target_detection_prob: float = 0.90,
        k_nearest: int = 3,
        specific_ingroup_sd_threshold: float = 3.0,
        wrong_group_margin_threshold: float = 0.02,
        threshold_json: str = None,
        scaffold_extensions: tuple = (".fasta", ".fa", ".fas")
) -> dict:
    """
    Validate chimera detection by running verify_by_ingroup_distance()
    on known genuine and artefactual sequences and reporting how well
    the function discriminates between them.

    Uses the same input files and scaffold directory as
    calibrate_distance_thresholds(). Does not calibrate anything —
    just produces a CSV showing whether each sequence was correctly
    classified.

    FASTA header formats accepted
    ------------------------------
    Genuine  : >Genus_species|gene|orien:+|accession:XXXX
    Artefact : >Genus_species|gene|mod_N|orien:+|accession:XXXX

    Parameters
    ----------
    genuine_fasta                : str   known genuine sequences
    artefact_fasta               : str   known artefactual sequences
    scaffold_dir                 : str   directory of per-gene scaffold
                                         FASTA files named {gene}.fasta
    taxonomy                     : dict  from load_taxonomy()
    output_dir                   : str   directory to save CSV
    insert_fraction              : float passed to verify function
    target_detection_prob        : float passed to verify function
    k_nearest                    : int   passed to verify function
    specific_ingroup_sd_threshold: float passed to verify function
    scaffold_extensions          : tuple file extensions to search for

    Returns
    -------
    dict with:
        results        : list   all per-sequence result dicts
        summary_csv    : str    path to saved CSV
        true_positives : int    artefacts correctly flagged
        true_negatives : int    genuine sequences correctly not flagged
        false_positives: int    genuine sequences incorrectly flagged
        false_negatives: int    artefacts incorrectly not flagged
        sensitivity    : float  TP / (TP + FN)
        specificity    : float  TN / (TN + FP)
    """
    import csv
    from pathlib import Path

    if isinstance(taxonomy, str):
        raise TypeError(
            "taxonomy must be a dict from load_taxonomy(), not a path."
        )

    os.makedirs(output_dir, exist_ok=True)

    # ── Helper: find scaffold FASTA for a gene ────────────────────────────────

    def find_scaffold(gene_name: str):
        scaffold_path = Path(scaffold_dir)
        for ext in scaffold_extensions:
            candidate = scaffold_path / f"{gene_name}{ext}"
            if candidate.exists():
                return str(candidate)
            for f in scaffold_path.iterdir():
                if (f.stem.lower() == gene_name.lower()
                        and f.suffix.lower() == ext.lower()):
                    return str(f)
        return None

    # ── Helper: parse header ──────────────────────────────────────────────────

    def parse_header(description: str) -> dict:
        parts       = description.strip().split("|")
        taxon       = parts[0].strip()
        gene        = parts[1].strip() if len(parts) > 1 else None
        accession   = None
        is_modified = False

        for p in parts:
            if p.startswith("accession:"):
                accession = p.split(":", 1)[1].strip()
            if p.startswith("mod_"):
                is_modified = True

        return {
            "taxon":       taxon,
            "gene":        gene,
            "accession":   accession,
            "is_modified": is_modified,
        }

    # ── Run on one FASTA file ─────────────────────────────────────────────────

    def run_on_fasta(fasta_path: str, label: str) -> list:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        results = []

        print(f"\n{'='*60}")
        print(f"Processing {label}: {len(records)} sequences")
        print(f"{'='*60}")

        for i, record in enumerate(records, 1):
            info = parse_header(record.description)

            if not info["gene"] or not info["accession"]:
                print(f"  [{i}] Skipping {record.description[:50]} "
                      f"— missing gene or accession")
                continue

            scaffold_path = find_scaffold(info["gene"])
            if scaffold_path is None:
                print(f"  [{i}] Skipping {info['taxon']} | "
                      f"{info['gene']} — no scaffold found")
                continue

            print(f"\n  [{i}/{len(records)}] {info['taxon']} | "
                  f"{info['gene']} | {info['accession']} "
                  f"[{label}]")

            try:
                result = verify_by_ingroup_distance(
                    query_seq             = str(record.seq),
                    query_taxon           = info["taxon"],
                    query_accession       = info["accession"],
                    scaffold_fasta        = scaffold_path,
                    taxonomy              = taxonomy,
                    insert_fraction       = insert_fraction,
                    target_detection_prob = target_detection_prob,
                    k_nearest             = k_nearest,
                    specific_ingroup_sd_threshold =
                        specific_ingroup_sd_threshold,
                    wrong_group_margin_threshold=
                        wrong_group_margin_threshold,
                )

                # Add classification metadata
                result["label"]            = label
                result["gene"]             = info["gene"]
                result["header"]           = record.description
                result["scaffold_used"]    = scaffold_path
                result["expected_chimera"] = (label == "ARTEFACT")
                result["correctly_classified"] = (
                    (label == "ARTEFACT" and result["chimera_flag"])
                    or
                    (label == "GENUINE"  and not result["chimera_flag"])
                )
                result["classification"] = (
                    "TRUE_POSITIVE"  if label == "ARTEFACT"
                                        and result["chimera_flag"]
                    else "FALSE_NEGATIVE" if label == "ARTEFACT"
                                        and not result["chimera_flag"]
                    else "TRUE_NEGATIVE"  if label == "GENUINE"
                                        and not result["chimera_flag"]
                    else "FALSE_POSITIVE"
                )

                results.append(result)

                icon = "✓" if result["correctly_classified"] else "✗"
                print(f"    {icon} chimera_flag="
                      f"{result['chimera_flag']} | "
                      f"expected={'chimera' if label == 'ARTEFACT' else 'genuine'} | "
                      f"{result['classification']} | "
                      f"outliers={result['n_outliers']}/"
                      f"{result['n_samples']}")
                print(f"    group={result['query_specific_group']} | "
                      f"scoring_level={result['scoring_level']} | "
                      f"evidence={result['chimera_evidence'][:80]}")
                      
            except Exception as e:
                print(f"    ✗ Failed: {e}")
                import traceback
                traceback.print_exc()

        return results

    # ── Run on both files ─────────────────────────────────────────────────────

    genuine_results  = run_on_fasta(genuine_fasta,  "GENUINE")
    artefact_results = run_on_fasta(artefact_fasta, "ARTEFACT")
    all_results      = genuine_results + artefact_results

    # ── Compute performance metrics ───────────────────────────────────────────

    tp = sum(1 for r in all_results
             if r["classification"] == "TRUE_POSITIVE")
    tn = sum(1 for r in all_results
             if r["classification"] == "TRUE_NEGATIVE")
    fp = sum(1 for r in all_results
             if r["classification"] == "FALSE_POSITIVE")
    fn = sum(1 for r in all_results
             if r["classification"] == "FALSE_NEGATIVE")

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else None
    specificity = tn / (tn + fp) if (tn + fp) > 0 else None
    precision   = tp / (tp + fp) if (tp + fp) > 0 else None
    f1          = (2 * precision * sensitivity /
                   (precision + sensitivity)
                   if precision and sensitivity else None)

    # ── Per-gene breakdown ────────────────────────────────────────────────────

    from collections import defaultdict
    gene_breakdown = defaultdict(lambda: {
        "true_positive": 0,
        "true_negative": 0,
        "false_positive": 0,
        "false_negative": 0,
    })
    for r in all_results:
        gene_breakdown[r["gene"]][
            r["classification"].lower()
        ] += 1

    # ── Write CSV ─────────────────────────────────────────────────────────────

    csv_path = os.path.join(output_dir,
                             "chimera_detection_validation.csv")

    csv_keys = [
        "label",
        "classification",
        "correctly_classified",
        "gene",
        "query_taxon",
        "query_accession",
        "query_family",
        "expected_order",
        "chimera_flag",
        "n_outliers",
        "n_samples",
        "prop_outliers",
        "specific_ingroup_mean",
        "specific_ingroup_std",
        "specific_ingroup_variance",
        "outlier_threshold",
        "n_correct_family",
        "prop_correct_family",
        "chimera_evidence",
        "scaffold_used",
    ]

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_keys,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_results)

    # ── Print summary ─────────────────────────────────────────────────────────

    print(f"\n{'='*60}")
    print(f"VALIDATION SUMMARY")
    print(f"{'='*60}")
    print(f"  Genuine sequences  : {len(genuine_results)}")
    print(f"  Artefact sequences : {len(artefact_results)}")
    print(f"\n  TRUE POSITIVES  (artefacts correctly flagged)  : {tp}")
    print(f"  TRUE NEGATIVES  (genuine correctly not flagged) : {tn}")
    print(f"  FALSE POSITIVES (genuine incorrectly flagged)   : {fp}")
    print(f"  FALSE NEGATIVES (artefacts missed)              : {fn}")
    print(f"\n  Sensitivity (recall)    : "
          f"{sensitivity:.3f}" if sensitivity else "  Sensitivity: N/A")
    print(f"  Specificity             : "
          f"{specificity:.3f}" if specificity else "  Specificity: N/A")
    print(f"  Precision               : "
          f"{precision:.3f}" if precision else "  Precision: N/A")
    print(f"  F1 score                : "
          f"{f1:.3f}" if f1 else "  F1: N/A")

    if fp > 0:
        print(f"\n  False positives (genuine flagged as chimeric):")
        for r in all_results:
            if r["classification"] == "FALSE_POSITIVE":
                print(f"    {r['query_taxon']} ({r['query_accession']}) "
                      f"| {r['gene']} | "
                      f"{r['chimera_evidence'][:80]}")

    if fn > 0:
        print(f"\n  False negatives (artefacts not detected):")
        for r in all_results:
            if r["classification"] == "FALSE_NEGATIVE":
                print(f"    {r['query_taxon']} ({r['query_accession']}) "
                      f"| {r['gene']} | "
                      f"{r['chimera_evidence'][:80]}")

    print(f"\n  Per-gene breakdown:")
    print(f"  {'Gene':<15} {'TP':>4} {'TN':>4} "
          f"{'FP':>4} {'FN':>4} {'Sens':>6} {'Spec':>6}")
    print(f"  {'-'*50}")
    for gene, counts in sorted(gene_breakdown.items()):
        tp_g = counts["true_positive"]
        tn_g = counts["true_negative"]
        fp_g = counts["false_positive"]
        fn_g = counts["false_negative"]
        g_sens = (tp_g / (tp_g + fn_g)
                  if (tp_g + fn_g) > 0 else None)
        g_spec = (tn_g / (tn_g + fp_g)
                  if (tn_g + fp_g) > 0 else None)
        print(f"  {gene:<15} "
              f"{tp_g:>4} "
              f"{tn_g:>4} "
              f"{fp_g:>4} "
              f"{fn_g:>4} "
              f"{f'{g_sens:.3f}' if g_sens is not None else 'N/A':>6} "
              f"{f'{g_spec:.3f}' if g_spec is not None else 'N/A':>6}")

    print(f"\n  CSV saved: {csv_path}")

    return {
        "results":         all_results,
        "summary_csv":     csv_path,
        "true_positives":  tp,
        "true_negatives":  tn,
        "false_positives": fp,
        "false_negatives": fn,
        "sensitivity":     sensitivity,
        "specificity":     specificity,
        "precision":       precision,
        "f1":              f1,
        "gene_breakdown":  dict(gene_breakdown),
    }


def calibrate_ingroup_thresholds(
        genuine_fasta: str,
        artefact_fasta: str,
        scaffold_dir: str,
        taxonomy: dict,
        output_path: str,
        insert_fraction: float = 0.30,
        target_detection_prob: float = 0.90,
        k_nearest: int = 3,
        sd_threshold_range: tuple = (1.5, 2.0, 2.5, 3.0, 3.5, 4.0),
        margin_threshold_range: tuple = (0.0, 0.01, 0.02, 0.03, 0.05),
        scaffold_extensions: tuple = (".fasta", ".fa", ".fas"),
) -> dict:
    """
    Calibrate gene-specific thresholds for verify_by_ingroup_distance()
    by finding the combination of specific_ingroup_sd_threshold and
    wrong_group_margin_threshold that maximises F1 score for each gene
    separately.

    Runs verify_by_ingroup_distance() on all genuine and artefact
    sequences with permissive settings to collect raw scores, then
    sweeps threshold combinations per gene to find the best F1.

    Parameters
    ----------
    genuine_fasta         : str    known genuine sequences
    artefact_fasta        : str    known artefactual sequences
    scaffold_dir          : str    per-gene scaffold FASTA directory
    taxonomy              : dict   from load_taxonomy()
    output_path           : str    path to save calibration JSON.
                                   Pass this to verify_by_ingroup_distance
                                   as threshold_json=output_path to use
                                   gene-specific thresholds automatically.
    insert_fraction       : float  passed to verify function
    target_detection_prob : float  passed to verify function
    k_nearest             : int    passed to verify function
    sd_threshold_range    : tuple  SD multiplier values to sweep
    margin_threshold_range: tuple  wrong_group_margin values to sweep
    scaffold_extensions   : tuple  file extensions to search for

    Returns
    -------
    dict with:
        thresholds   : dict   gene -> {sd_threshold, margin_threshold,
                                       f1, sensitivity, specificity}
        summary_csv  : str    path to per-sequence results CSV
        threshold_json: str   path to saved threshold JSON
    """
    import json
    import csv
    from pathlib import Path
    from collections import defaultdict

    if isinstance(taxonomy, str):
        raise TypeError(
            "taxonomy must be a dict from load_taxonomy(), not a path."
        )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # ── Helper: find scaffold for a gene ──────────────────────────────────────
    def find_scaffold(gene_name: str):
        scaffold_path = Path(scaffold_dir)
        for ext in scaffold_extensions:
            candidate = scaffold_path / f"{gene_name}{ext}"
            if candidate.exists():
                return str(candidate)
            for f in scaffold_path.iterdir():
                if (f.stem.lower() == gene_name.lower()
                        and f.suffix.lower() == ext.lower()):
                    return str(f)
        return None

    # ── Helper: parse FASTA header ────────────────────────────────────────────
    def parse_header(description: str) -> dict:
        parts     = description.strip().split("|")
        taxon     = parts[0].strip()
        gene      = parts[1].strip() if len(parts) > 1 else None
        accession = None
        for p in parts:
            if p.startswith("accession:"):
                accession = p.split(":", 1)[1].strip()
        return {"taxon": taxon, "gene": gene, "accession": accession}

    # ── Step 1: collect raw per-subsample scores using permissive settings ────
    # We collect the raw scores here so we can sweep thresholds
    # without re-running the alignment for each combination.
    # Use sd_threshold=99 and margin=0 so nothing is flagged —
    # we just want the raw identity scores.

    print(f"\n{'='*60}")
    print(f"COLLECTING RAW SCORES FOR CALIBRATION")
    print(f"{'='*60}")

    raw_results = []  # list of dicts, one per sequence

    for fasta_path, label in [(genuine_fasta,  "GENUINE"),
                               (artefact_fasta, "ARTEFACT")]:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        print(f"\nProcessing {label}: {len(records)} sequences")

        for i, record in enumerate(records, 1):
            info = parse_header(record.description)
            if not info["gene"] or not info["accession"]:
                continue

            scaffold_path = find_scaffold(info["gene"])
            if scaffold_path is None:
                print(f"  [{i}] Skipping {info['taxon']} | "
                      f"{info['gene']} — no scaffold found")
                continue

            print(f"  [{i}/{len(records)}] {info['taxon']} | "
                  f"{info['gene']} [{label}]")

            try:
                # Run with permissive thresholds — collect raw scores
                result = verify_by_ingroup_distance(
                    query_seq             = str(record.seq),
                    query_taxon           = info["taxon"],
                    query_accession       = info["accession"],
                    scaffold_fasta        = scaffold_path,
                    taxonomy              = taxonomy,
                    insert_fraction       = insert_fraction,
                    target_detection_prob = target_detection_prob,
                    k_nearest             = k_nearest,
                    # Permissive — collect all raw scores
                    specific_ingroup_sd_threshold = 99.0,
                    wrong_group_margin_threshold  = 99.0,
                )

                raw_results.append({
                    "label":        label,
                    "gene":         info["gene"],
                    "taxon":        info["taxon"],
                    "accession":    info["accession"],
                    "is_artefact":  label == "ARTEFACT",
                    "sample_scores":result["sample_scores"],
                    "specific_ingroup_mean": result["specific_ingroup_mean"],
                    "specific_ingroup_std":  result["specific_ingroup_std"],
                    "query_specific_group":  result["query_specific_group"],
                    "scoring_level":         result["scoring_level"],
                })

            except Exception as e:
                print(f"  ✗ Failed: {e}")
                import traceback
                traceback.print_exc()

    if not raw_results:
        raise ValueError("No sequences were successfully scored.")

    # ── Step 2: sweep thresholds per gene ─────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"SWEEPING THRESHOLDS PER GENE")
    print(f"{'='*60}")

    # Group results by gene
    by_gene = defaultdict(list)
    for r in raw_results:
        by_gene[r["gene"]].append(r)

    gene_thresholds = {}

    for gene, seqs in sorted(by_gene.items()):
        genuine_seqs  = [s for s in seqs if not s["is_artefact"]]
        artefact_seqs = [s for s in seqs if s["is_artefact"]]

        if not genuine_seqs or not artefact_seqs:
            print(f"\n  {gene}: insufficient data "
                  f"(genuine={len(genuine_seqs)}, "
                  f"artefact={len(artefact_seqs)}) — skipping")
            continue

        print(f"\n  {gene}: {len(genuine_seqs)} genuine, "
              f"{len(artefact_seqs)} artefact")

        best_f1          = -1.0
        best_sensitivity = 0.0
        best_specificity = 0.0
        best_sd          = 2.0
        best_margin      = 0.02

        # Sweep all combinations of sd_threshold and margin_threshold
        for sd_thresh in sd_threshold_range:
            for margin_thresh in margin_threshold_range:

                tp = tn = fp = fn = 0

                for seq_result in seqs:
                    # Recompute outlier flags with these thresholds
                    scores       = seq_result["sample_scores"]
                    sp_mean      = seq_result["specific_ingroup_mean"]
                    sp_std       = seq_result["specific_ingroup_std"]
                    threshold    = sp_mean - sd_thresh * sp_std

                    any_outlier = False
                    for s in scores:
                        sp_id      = s["specific_ingroup_identity"]
                        sp_knn     = s["specific_ingroup_knn"]
                        nearest    = s["nearest_group_match"]
                        correct    = s["is_correct_group"]

                        # Recompute wrong group margin
                        if not correct and "per_group_identity" in s:
                            nearest_knn = s["per_group_identity"].get(
                                nearest, {}
                            ).get("knn_mean", sp_knn)
                            margin = nearest_knn - sp_knn
                        else:
                            margin = 0.0

                        wrong_flag = (
                            not correct
                            and margin > margin_thresh
                        )
                        is_outlier = (
                            sp_id < threshold or wrong_flag
                        )

                        if is_outlier:
                            any_outlier = True
                            break

                    predicted_chimera = any_outlier
                    actual_chimera    = seq_result["is_artefact"]

                    if predicted_chimera and actual_chimera:
                        tp += 1
                    elif not predicted_chimera and not actual_chimera:
                        tn += 1
                    elif predicted_chimera and not actual_chimera:
                        fp += 1
                    else:
                        fn += 1

                sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
                precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
                f1          = (2 * precision * sensitivity /
                               (precision + sensitivity)
                               if (precision + sensitivity) > 0
                               else 0.0)

                # Prefer higher F1; break ties by higher specificity
                # (we already have good sensitivity — want fewer FP)
                if (f1 > best_f1
                        or (f1 == best_f1
                            and specificity > best_specificity)):
                    best_f1          = f1
                    best_sensitivity = sensitivity
                    best_specificity = specificity
                    best_sd          = sd_thresh
                    best_margin      = margin_thresh

        gene_thresholds[gene] = {
            "sd_threshold":     best_sd,
            "margin_threshold": best_margin,
            "f1":               round(best_f1, 4),
            "sensitivity":      round(best_sensitivity, 4),
            "specificity":      round(best_specificity, 4),
            "n_genuine":        len(genuine_seqs),
            "n_artefact":       len(artefact_seqs),
        }

        print(f"    Best: sd={best_sd}, margin={best_margin} → "
              f"F1={best_f1:.3f}, "
              f"sensitivity={best_sensitivity:.3f}, "
              f"specificity={best_specificity:.3f}")

    # ── Step 3: save thresholds to JSON ───────────────────────────────────────
    with open(output_path, "w") as f:
        json.dump(gene_thresholds, f, indent=2)

    print(f"\n{'='*60}")
    print(f"CALIBRATION COMPLETE")
    print(f"{'='*60}")
    print(f"\n{'Gene':<15} {'SD':>5} {'Margin':>8} "
          f"{'F1':>6} {'Sens':>6} {'Spec':>6}")
    print(f"{'-'*50}")
    for gene, t in sorted(gene_thresholds.items()):
        print(f"{gene:<15} {t['sd_threshold']:>5} "
              f"{t['margin_threshold']:>8.3f} "
              f"{t['f1']:>6.3f} "
              f"{t['sensitivity']:>6.3f} "
              f"{t['specificity']:>6.3f}")

    print(f"\nThresholds saved to: {output_path}")

    return {
        "thresholds":    gene_thresholds,
        "threshold_json": output_path,
        "raw_results":   raw_results,
    }

from bioinf_packages.verify_funcs._verify_distance import (verify_by_distance, verify_by_ingroup_distance, verify_batch_by_ingroup_distance)
from bioinf_packages.verify_funcs._species_parser import load_taxonomy

#cal = calibrate_distance_thresholds(
#    genuine_fasta  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta",
#    artefact_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_seq/modifed_seqs.fasta",
#    scaffold_dir   = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas",
#    taxonomy       = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
#)

#print(cal)



# Step 1 — calibrate
cal = calibrate_ingroup_thresholds(
    genuine_fasta  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta",
    artefact_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_seq/modifed_seqs.fasta",
    scaffold_dir   = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas",
    taxonomy       = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
    output_path    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gene_thresholds.json",
)

# Step 2 — validate with calibrated thresholds
validation = validate_chimera_detection(
    genuine_fasta  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta",
    artefact_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_seq/modifed_seqs.fasta",
    scaffold_dir   = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas",
    taxonomy       = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
    output_dir     = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/validation_output",
    threshold_json = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gene_thresholds.json",
)

# Step 3 — use calibrated thresholds for real verification
results = verify_batch_by_ingroup_distance(
    query_fasta    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data/CYTB.fasta",
    scaffold_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_generation/fasta_data/CYTB.fasta",
    taxonomy       = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"),
    output_dir     = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/validation_output",
    threshold_json = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gene_thresholds.json",
)