# _compute_csv.py

import csv
import json
import math
import numpy as np
from pathlib import Path
from collections import defaultdict
import copy


from bioinf_packages.verify_funcs._score_hit import (
    softmax_weights,
    sigmoid,
    pident_to_llr,
)


# ─────────────────────────────────────────────────────────────────────────────
# HIT FILTER CALIBRATION
# ─────────────────────────────────────────────────────────────────────────────

def calibrate_hit_filters(results_csv: str,
                           output_path: str = None) -> dict:
    """
    Analyse the distribution of evalue and pident across BLAST hits
    to find appropriate filtering thresholds, with per-gene breakdown.

    Parameters
    ----------
    results_csv : str   path to blast_results.csv
    output_path : str   optional path to save recommendations as JSON

    Returns
    -------
    dict with recommended_evalue, recommended_pident,
              recommended_align_length, per_gene_thresholds
    """
    genuine_hits = []

    with open(results_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("hit_accession") == "NO_HITS":
                continue
            if row.get("label", "").strip() != "GENUINE":
                continue
            try:
                genuine_hits.append({
                    "evalue":       float(row["evalue"]),
                    "pident":       float(row["pident"]),
                    "align_length": int(row["align_length"]),
                    "bitscore":     float(row["bitscore"]),
                    "species_llr":  float(row["species_llr"]),
                    "gene":         row.get("gene", ""),
                    "taxon":        row.get("taxon", ""),
                    "accession":    row.get("accession", ""),
                })
            except (ValueError, KeyError):
                continue

    if not genuine_hits:
        raise ValueError("No GENUINE hits found in CSV.")

    correct_hits  = [h for h in genuine_hits if h["species_llr"] > 0]
    spurious_hits = [h for h in genuine_hits if h["species_llr"] <= 0]

    print(f"\nLoaded {len(genuine_hits)} genuine hits")
    print(f"  Taxonomically correct  : {len(correct_hits)}")
    print(f"  Taxonomically spurious : {len(spurious_hits)}")
    print(f"  Spurious rate          : "
          f"{len(spurious_hits)/len(genuine_hits)*100:.1f}%")

    by_gene = defaultdict(lambda: {"correct": [], "spurious": []})
    for h in correct_hits:
        by_gene[h["gene"]]["correct"].append(h)
    for h in spurious_hits:
        by_gene[h["gene"]]["spurious"].append(h)

    evalue_candidates = [1e-10, 1e-12, 1e-15, 1e-18, 1e-20,
                         1e-25, 1e-30, 1e-40, 1e-50]
    pident_candidates = [70.0, 75.0, 80.0, 82.0, 85.0,
                         88.0, 90.0, 92.0, 95.0]
    align_candidates  = [44, 50, 60, 70, 80]

    print(f"\n── Per-gene threshold optimisation ─────────────────────────")

    per_gene_thresholds = {}

    for gene, hits in sorted(by_gene.items()):
        correct  = hits["correct"]
        spurious = hits["spurious"]

        if not correct:
            continue

        n_correct  = len(correct)
        n_spurious = len(spurious)
        best_precision   = -1.0
        best_evalue      = 1e-10
        best_pident      = 70.0
        best_align       = 44
        best_c_retained  = n_correct
        best_s_retained  = n_spurious

        for ev in evalue_candidates:
            for pid in pident_candidates:
                for aln in align_candidates:
                    c_kept = sum(
                        1 for h in correct
                        if h["evalue"] <= ev
                        and h["pident"] >= pid
                        and h["align_length"] >= aln
                    )
                    s_kept = sum(
                        1 for h in spurious
                        if h["evalue"] <= ev
                        and h["pident"] >= pid
                        and h["align_length"] >= aln
                    )
                    if n_correct > 0 and c_kept < n_correct * 0.9:
                        continue
                    total_kept = c_kept + s_kept
                    if total_kept == 0:
                        continue
                    precision = c_kept / total_kept
                    if precision > best_precision:
                        best_precision  = precision
                        best_evalue     = ev
                        best_pident     = pid
                        best_align      = aln
                        best_c_retained = c_kept
                        best_s_retained = s_kept

        c_evalues = [h["evalue"] for h in correct]
        c_pidents = [h["pident"] for h in correct]

        per_gene_thresholds[gene] = {
            "recommended_evalue":       best_evalue,
            "recommended_pident":       best_pident,
            "recommended_align_length": best_align,
            "filter_precision":         round(best_precision, 4),
            "correct_hits_retained":    best_c_retained,
            "correct_hits_total":       n_correct,
            "spurious_hits_retained":   best_s_retained,
            "spurious_hits_total":      n_spurious,
            "pct_correct_retained":     round(
                best_c_retained / n_correct * 100, 1),
            "correct_evalue_max":       float(max(c_evalues)),
            "correct_pident_min":       float(min(c_pidents)),
        }

        print(f"\n  {gene} ({n_correct} correct, {n_spurious} spurious):")
        print(f"    Best evalue      : {best_evalue:.1e}")
        print(f"    Best pident      : {best_pident:.1f}%")
        print(f"    Best align_length: {best_align}bp")
        print(f"    Filter precision : {best_precision:.3f} "
              f"({best_c_retained}/{n_correct} correct retained, "
              f"{best_s_retained}/{n_spurious} spurious retained)")
        print(f"    Correct hit range: evalue up to "
              f"{max(c_evalues):.1e}, "
              f"pident down to {min(c_pidents):.1f}%")

    if per_gene_thresholds:
        global_evalue = max(
            v["recommended_evalue"]
            for v in per_gene_thresholds.values()
        )
        global_pident = min(
            v["recommended_pident"]
            for v in per_gene_thresholds.values()
        )
        global_align  = min(
            v["recommended_align_length"]
            for v in per_gene_thresholds.values()
        )
    else:
        global_evalue = 1e-10
        global_pident = 70.0
        global_align  = 50

    print(f"\n{'='*60}")
    print(f"RECOMMENDED GLOBAL THRESHOLDS")
    print(f"{'='*60}")
    print(f"  evalue        : {global_evalue:.1e}")
    print(f"  pident        : {global_pident:.1f}%")
    print(f"  align_length  : {global_align}bp")

    results = {
        "recommended_evalue":        global_evalue,
        "recommended_pident":        global_pident,
        "recommended_align_length":  global_align,
        "per_gene_thresholds":       per_gene_thresholds,
        "n_genuine_hits":            len(genuine_hits),
        "n_correct_hits":            len(correct_hits),
        "n_spurious_hits":           len(spurious_hits),
        "spurious_rate_pct":         round(
            len(spurious_hits) / len(genuine_hits) * 100, 2),
    }

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nSaved: {output_path}")

    return results


# ─────────────────────────────────────────────────────────────────────────────
# MITOCHONDRIAL GENE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

MITO_GENES_SET = {
    "COX1", "COX2", "COX3",
    "ND1",  "ND2",  "ND3",  "ND4", "ND4L", "ND5", "ND6",
    "ATP6", "ATP8", "CYTB",
    "12S_rRNA", "16S_rRNA",
    "12S_RRNA", "16S_RRNA",
}


def _is_mito_gene(gene: str) -> bool:
    """Check whether a gene name refers to a mitochondrial gene."""
    return (
        gene in MITO_GENES_SET
        or gene.upper() in MITO_GENES_SET
        or gene.upper().replace("_RRNA", "") in {"12S", "16S",
                                                   "RRNS", "RRNL"}
    )


# ─────────────────────────────────────────────────────────────────────────────
# CROSS-GENE CONSISTENCY
# ─────────────────────────────────────────────────────────────────────────────

def apply_cross_gene_consistency(
        metrics: dict,
        min_well_scoring_fraction: float = 0.6,
        well_scoring_threshold_prop_neg: float = 0.4,
        max_poorly_scoring_for_rescue: int = 3,
        # Do not rescue if more than this many genes are poorly
        # scoring — indicates a real chimera rather than noise.
        # Default 3: rescue fires only when <= 3 genes are flagged.
        known_orphan_taxa: list = None,
        orphan_prop_neg_threshold: float = 0.8,
) -> dict:

    result     = copy.deepcopy(metrics)
    orphan_set = set(known_orphan_taxa or [])

    def get_genus(parent_key: tuple) -> str:
        taxon = result.get(parent_key, {}).get("taxon", "")
        return taxon.split("_")[0] if taxon else ""

    # Group mito genes by (accession, label) — never compare across
    # labels. Artefact sequences share accession numbers with genuine
    # sequences they were derived from, so grouping by accession alone
    # would cause genuine well-scoring genes to rescue artefact flags.
    mito_by_acc = defaultdict(dict)
    for (acc, gene) in result:
        if _is_mito_gene(gene):
            label = result[(acc, gene)].get("label", "UNKNOWN")
            mito_by_acc[(acc, label)][gene] = (acc, gene)

    for (acc, label), gene_key_map in mito_by_acc.items():

        if len(gene_key_map) < 2:
            gene       = list(gene_key_map.keys())[0]
            parent_key = (acc, gene)
            result[parent_key].update({
                "cross_gene_rescue":      False,
                "cross_gene_evidence":    "Only one mito gene for "
                                          "this accession",
                "n_mito_genes_accession": 1,
                "n_well_scoring_mito":    0,
                "well_scoring_fraction":  0.0,
            })
            continue

        # Classify using chimera_flag — a gene not flagged is
        # well-scoring by definition
        well_scoring   = set()
        poorly_scoring = set()

        for gene, parent_key in gene_key_map.items():
            if not result[parent_key].get("chimera_flag", False):
                well_scoring.add(gene)
            else:
                poorly_scoring.add(gene)

        n_total   = len(gene_key_map)
        n_well    = len(well_scoring)
        frac_well = n_well / n_total

        for gene, parent_key in gene_key_map.items():
            m          = result[parent_key]
            is_poor    = gene in poorly_scoring
            is_flagged = m.get("chimera_flag", False)
            genus      = get_genus(parent_key)
            is_orphan  = genus in orphan_set

            result[parent_key]["n_mito_genes_accession"] = n_total
            result[parent_key]["n_well_scoring_mito"]    = n_well
            result[parent_key]["well_scoring_fraction"]  = round(
                frac_well, 4)

            if not is_poor or not is_flagged:
                result[parent_key]["cross_gene_rescue"]   = False
                result[parent_key]["cross_gene_evidence"] = (
                    f"Gene scores well — no rescue needed. "
                    f"{n_well}/{n_total} well-scoring for {acc}."
                )
                continue

            if is_orphan and label == "GENUINE":
                # Orphan rescue only applies to genuine sequences.
                # Artefact sequences from orphan taxa are still
                # chimeric — poor database coverage is not an excuse
                # for a genuinely chimeric insertion.
                result[parent_key]["chimera_flag"]        = False
                result[parent_key]["cross_gene_rescue"]   = True
                result[parent_key]["cross_gene_evidence"] = (
                    f"Flag REMOVED — orphan taxon ({genus}). "
                    f"Poor BLAST scores expected due to sparse "
                    f"database coverage."
                )
                print(f"  ✓ ORPHAN RESCUE: {acc} | {label} | "
                      f"{gene} | genus={genus}")

            elif is_orphan and label == "ARTEFACT":
                # Do not rescue artefact sequences even for orphan taxa
                result[parent_key]["cross_gene_rescue"]   = False
                result[parent_key]["cross_gene_evidence"] = (
                    f"Flag KEPT — orphan taxon ({genus}) but label "
                    f"is ARTEFACT. Chimeric insertion is genuine "
                    f"regardless of database coverage."
                )

            elif (frac_well >= min_well_scoring_fraction
                and len(poorly_scoring) <= max_poorly_scoring_for_rescue):
                result[parent_key]["chimera_flag"]        = False
                result[parent_key]["cross_gene_rescue"]   = True
                result[parent_key]["cross_gene_evidence"] = (
                    f"Flag REMOVED by cross-gene rescue. "
                    f"{n_well}/{n_total} well-scoring "
                    f"(frac={frac_well:.2f} >= "
                    f"{min_well_scoring_fraction}). "
                    f"Well-scoring: "
                    f"{', '.join(sorted(well_scoring))}."
                )
                print(f"  ✓ RESCUED: {acc} | {label} | {gene} | "
                      f"{n_well}/{n_total} well-scoring")

            else:
                result[parent_key]["cross_gene_rescue"]   = False
                result[parent_key]["cross_gene_evidence"] = (
                    f"Flag KEPT. Only {n_well}/{n_total} well-scoring "
                    f"(frac={frac_well:.2f} < "
                    f"{min_well_scoring_fraction})."
                )
                print(f"  ⚠ KEPT: {acc} | {label} | {gene} | "
                      f"{n_well}/{n_total} well-scoring")

    for parent_key in result:
        if "cross_gene_rescue" not in result[parent_key]:
            result[parent_key].update({
                "cross_gene_rescue":      False,
                "cross_gene_evidence":    "Non-mitochondrial gene",
                "n_mito_genes_accession": 0,
                "n_well_scoring_mito":    0,
                "well_scoring_fraction":  0.0,
            })

    return result


# ─────────────────────────────────────────────────────────────────────────────
# GENE-SPECIFIC PROP_NEG THRESHOLD CALIBRATION
# ─────────────────────────────────────────────────────────────────────────────

def calibrate_gene_prop_neg_thresholds(
        genuine_by_parent: dict,
        artefact_by_parent: dict,
        prior: float,
        pident_thresh: float,
        k: float,
        beta: float,
        compute_metrics_fn,
        candidate_thresholds: tuple = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
) -> dict:
    """
    Calibrate gene-specific prop_neg_taxon_hits thresholds.

    For each gene, finds the prop_neg threshold that maximises F1
    between genuine and artefact sequences at the given Bayesian
    scoring parameters.

    Parameters
    ----------
    genuine_by_parent   : dict      from calibrate_from_csv()
    artefact_by_parent  : dict      from calibrate_from_csv()
    prior               : float     best prior_prob from sweep
    pident_thresh       : float     best pident_threshold from sweep
    k                   : float     best k from sweep
    beta                : float     best beta from sweep
    compute_metrics_fn  : callable  the compute_parent_metrics function
    candidate_thresholds: tuple     prop_neg values to try per gene

    Returns
    -------
    dict   gene -> best prop_neg threshold
    """
    g_metrics = compute_metrics_fn(
        genuine_by_parent,  "GENUINE",
        prior, pident_thresh, k, beta
    )
    a_metrics = compute_metrics_fn(
        artefact_by_parent, "ARTEFACT",
        prior, pident_thresh, k, beta
    )

    g_by_gene = defaultdict(list)
    a_by_gene = defaultdict(list)

    for (acc, gene), m in g_metrics.items():
        g_by_gene[gene].append(m["prop_neg_taxon_hits"])
    for (acc, gene), m in a_metrics.items():
        a_by_gene[gene].append(m["prop_neg_taxon_hits"])

    all_genes       = set(g_by_gene.keys()) | set(a_by_gene.keys())
    gene_thresholds = {}

    print(f"\n── Gene-specific prop_neg thresholds ───────────────────────")
    print(f"  {'Gene':<15} {'G mean':>8} {'A mean':>8} "
          f"{'Threshold':>10} {'F1':>6}")
    print(f"  {'-'*53}")

    for gene in sorted(all_genes):
        g_vals = g_by_gene.get(gene, [])
        a_vals = a_by_gene.get(gene, [])

        if not g_vals or not a_vals:
            gene_thresholds[gene] = 0.5
            continue

        best_f1     = -1.0
        best_thresh = 0.5

        for thresh in candidate_thresholds:
            tn = sum(1 for v in g_vals if v <= thresh)
            fp = sum(1 for v in g_vals if v > thresh)
            tp = sum(1 for v in a_vals if v > thresh)
            fn = sum(1 for v in a_vals if v <= thresh)

            precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            f1          = (2 * precision * sensitivity /
                           (precision + sensitivity)
                           if (precision + sensitivity) > 0 else 0.0)

            if f1 > best_f1:
                best_f1     = f1
                best_thresh = thresh

        gene_thresholds[gene] = best_thresh
        print(f"  {gene:<15} "
              f"{np.mean(g_vals):>8.3f} "
              f"{np.mean(a_vals):>8.3f} "
              f"{best_thresh:>10.2f} "
              f"{best_f1:>6.3f}")

    return gene_thresholds



def calibrate_from_csv(
        results_csv: str,
        output_path: str,
        prior_prob_range: tuple = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
        pident_threshold_range: tuple = (85.0, 88.0, 90.0, 92.0, 95.0),
        k_range: tuple = (0.1, 0.2, 0.5, 1.0),
        beta_range: tuple = (0.01, 0.05, 0.1),
        posterior_var_threshold_range: tuple = (0.01, 0.02, 0.05,
                                                0.10, 0.15),
        taxon_llr_var_threshold_range: tuple = (0.5, 1.0, 2.0,
                                                3.0, 5.0),
        prop_neg_threshold_range: tuple = (0.3, 0.4, 0.5, 0.6, 0.7),
        no_hits_threshold_range: tuple = (0.3, 0.5, 0.7, 0.9),
        max_poorly_scoring_range: tuple = (1, 2, 3, 4, 5),
        flag_on_single_hit: bool = True,
        single_hit_taxon_llr_threshold: float = -2.0,
        min_well_scoring_fraction: float = 0.6,
        well_scoring_threshold_prop_neg: float = 0.4,
        known_orphan_taxa: list = None,
        orphan_prop_neg_threshold: float = 0.8,
        orphan_rescue_mode: str = "genuine_only",
        label_map: dict = None,
) -> dict:
    """
    Calibrate Bayesian model parameters from blast_results.csv.

    Genuine and artefact sequences are kept in completely separate
    data structures because artefact sequences share accession numbers
    with genuine sequences they were derived from.

    Orphan rescue is always applied to GENUINE sequences only.
    Artefact sequences from orphan taxa are never rescued — poor
    database coverage does not excuse a genuine chimeric insertion.

    orphan_rescue_mode controls the parameter sweep and mode
    comparison table:
        "genuine_only" — rescue applied to GENUINE during sweep
                         (recommended, used for the saved JSON)
        "none"         — no rescue during sweep (baseline comparison)

    label_map maps raw CSV label values to calibration classes.
    Default treats "GENUINE" as genuine and "ARTEFACT" as artefact.
    Use to incorporate additional label types, e.g.:
        label_map = {
            "GENUINE":       "GENUINE",
            "ARTEFACT":      "ARTEFACT",
            "REAL_ARTEFACT": "ARTEFACT",
            "CHIMERA":       "ARTEFACT",
        }
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # ── Build label map ───────────────────────────────────────────────────────
    if label_map is None:
        label_map = {
            "GENUINE":  "GENUINE",
            "ARTEFACT": "ARTEFACT",
        }
    genuine_labels  = {k for k, v in label_map.items()
                       if v == "GENUINE"}
    artefact_labels = {k for k, v in label_map.items()
                       if v == "ARTEFACT"}

    print(f"\n{'='*60}")
    print(f"CALIBRATE FROM CSV")
    print(f"{'='*60}")
    print(f"Loading: {results_csv}")
    print(f"Orphan rescue mode: {orphan_rescue_mode}")
    print(f"Genuine labels : {sorted(genuine_labels)}")
    print(f"Artefact labels: {sorted(artefact_labels)}")

    # ── Load data ─────────────────────────────────────────────────────────────
    genuine_hits     = defaultdict(list)
    artefact_hits    = defaultdict(list)
    genuine_parents  = set()
    artefact_parents = set()
    parent_taxon     = {}

    # Track which (acc, gene, sample) keys had no hits so
    # compute_parent_metrics can correctly compute prop_no_hits.
    # NO_HITS rows are recorded here as empty hit lists rather than
    # being skipped entirely.
    genuine_no_hit_samples  = set()   # (acc, gene, sample) with no hits
    artefact_no_hit_samples = set()

    with open(results_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw_label = row.get("label", "").strip()
            label     = label_map.get(raw_label)
            if label is None:
                continue

            acc   = row.get("accession", "").strip()
            gene  = row.get("gene", "").strip()
            taxon = row.get("taxon", "").strip()
            if not acc or not gene:
                continue

            parent_key = (acc, gene)
            if label == "GENUINE":
                genuine_parents.add(parent_key)
            else:
                artefact_parents.add(parent_key)

            parent_taxon[(acc, gene, label)] = taxon

            if row.get("hit_accession") == "NO_HITS":
                # Record the sample as a no-hit subsample so
                # prop_no_hits can be computed correctly
                try:
                    sample = row["sample"]
                except KeyError:
                    continue
                no_hit_key = (acc, gene, sample)
                if label == "GENUINE":
                    genuine_no_hit_samples.add(no_hit_key)
                else:
                    artefact_no_hit_samples.add(no_hit_key)
                continue

            try:
                sample      = row["sample"]
                bitscore    = float(row["bitscore"])
                pident      = float(row["pident"])
                species_llr = float(row["species_llr"])
                gene_llr    = float(row["gene_llr"])
            except (ValueError, KeyError, TypeError):
                continue

            hit = {
                "bitscore":    bitscore,
                "pident":      pident,
                "species_llr": species_llr,
                "gene_llr":    gene_llr,
            }
            hit_key = (acc, gene, sample)
            if label == "GENUINE":
                genuine_hits[hit_key].append(hit)
            else:
                artefact_hits[hit_key].append(hit)

    print(f"  Genuine parent sequences  : {len(genuine_parents)}")
    print(f"  Artefact parent sequences : {len(artefact_parents)}")
    print(f"  Genuine subsample groups  : {len(genuine_hits)}")
    print(f"  Artefact subsample groups : {len(artefact_hits)}")

    if not genuine_parents:
        raise ValueError("No GENUINE sequences found in CSV.")
    if not artefact_parents:
        raise ValueError("No ARTEFACT sequences found in CSV.")

    # ── Group subsamples by parent ────────────────────────────────────────────
    def group_by_parent(hit_dict: dict,
                        no_hit_samples: set) -> dict:
        """
        Group hits by (acc, gene) parent key. Also adds empty hit
        lists for no-hit subsamples so prop_no_hits is computed
        correctly — without this, subsamples with no hits are
        invisible to compute_parent_metrics.
        """
        grouped = defaultdict(list)

        # Add subsamples that had hits
        for (acc, gene, sample), hits in hit_dict.items():
            grouped[(acc, gene)].append({
                "sample": sample,
                "hits":   hits,
            })

        # Add subsamples that had no hits as empty lists
        seen_keys = {(acc, gene, sd["sample"])
                     for (acc, gene), samples in grouped.items()
                     for sd in samples}
        for (acc, gene, sample) in no_hit_samples:
            if (acc, gene, sample) not in seen_keys:
                grouped[(acc, gene)].append({
                    "sample": sample,
                    "hits":   [],
                })

        return dict(grouped)

    genuine_by_parent  = group_by_parent(genuine_hits,
                                         genuine_no_hit_samples)
    artefact_by_parent = group_by_parent(artefact_hits,
                                         artefact_no_hit_samples)

    # ── compute_parent_metrics ────────────────────────────────────────────────
    def compute_parent_metrics(parent_samples_dict: dict,
                                label: str,
                                prior: float,
                                pident_thresh: float,
                                k: float,
                                beta: float) -> dict:
        results = {}
        for parent_key, samples in parent_samples_dict.items():
            acc, gene           = parent_key
            all_taxon_llrs      = []
            sample_posts        = []
            any_single_flag     = False
            prior_log_odds      = math.log(prior / (1 - prior))
            cumulative_log_odds = prior_log_odds

            n_subsamples_total   = len(samples)
            n_subsamples_no_hits = 0

            for sample_data in samples:
                hits = sample_data["hits"]
                if not hits:
                    n_subsamples_no_hits += 1
                    continue

                bitscores    = [h["bitscore"] for h in hits]
                weights      = softmax_weights(bitscores, beta=beta)
                sample_llr   = 0.0
                s_log_odds   = math.log(prior / (1 - prior))
                top_llr      = None
                top_bitscore = -1

                for h, w in zip(hits, weights):
                    llr_p     = pident_to_llr(h["pident"],
                                              threshold=pident_thresh,
                                              k=k)
                    llr_taxon = h["species_llr"]
                    llr_gene  = h["gene_llr"]
                    combined  = w * (llr_taxon + llr_gene + llr_p)
                    sample_llr         += combined
                    s_log_odds         += llr_taxon + llr_gene + llr_p
                    all_taxon_llrs.append(llr_taxon)
                    if h["bitscore"] > top_bitscore:
                        top_bitscore = h["bitscore"]
                        top_llr      = llr_taxon

                cumulative_log_odds += sample_llr
                sample_posts.append(sigmoid(s_log_odds))

                if (flag_on_single_hit
                        and top_llr is not None
                        and top_llr < single_hit_taxon_llr_threshold):
                    any_single_flag = True

            post_var = (float(np.var(sample_posts))
                        if len(sample_posts) > 1 else 0.0)
            llr_var  = (float(np.var(all_taxon_llrs))
                        if len(all_taxon_llrs) > 1 else 0.0)
            prop_neg = (sum(1 for x in all_taxon_llrs if x < 0) /
                        len(all_taxon_llrs)
                        if all_taxon_llrs else 0.0)
            prop_no_hits = (n_subsamples_no_hits / n_subsamples_total
                            if n_subsamples_total > 0 else 0.0)

            taxon         = parent_taxon.get((acc, gene, label), "")
            taxon_genus   = taxon.split("_")[0] if taxon else ""
            is_orphan     = taxon_genus in set(known_orphan_taxa or [])
            prelim_thresh = 0.3 if is_orphan else 0.5

            # Preliminary chimera flag for cross-gene rescue.
            # prop_no_hits > 0.5 included so sequences where most
            # subsamples return no hits are eligible for rescue
            # assessment rather than being silently cleared.
            results[parent_key] = {
                "posterior_variance":   post_var,
                "taxon_llr_variance":   llr_var,
                "prop_neg_taxon_hits":  prop_neg,
                "prop_no_hits":         prop_no_hits,
                "n_subsamples_no_hits": n_subsamples_no_hits,
                "n_subsamples_total":   n_subsamples_total,
                "mean_posterior":       (float(np.mean(sample_posts))
                                         if sample_posts else 0.0),
                "single_hit_flag":      any_single_flag,
                "cumulative_log_odds":  cumulative_log_odds,
                "label":                label,
                "taxon":                taxon,
                "chimera_flag":         (any_single_flag
                                         or prop_neg > prelim_thresh
                                         or prop_no_hits > 0.5),
            }
        return results

    # ── Discrimination report ─────────────────────────────────────────────────
    print(f"\n── Discrimination report (neutral parameters) ───────────────")

    neutral_genuine  = compute_parent_metrics(
        genuine_by_parent,  "GENUINE",  0.5, 90.0, 0.2, 0.05
    )
    neutral_artefact = compute_parent_metrics(
        artefact_by_parent, "ARTEFACT", 0.5, 90.0, 0.2, 0.05
    )

    def report_separation(metric_name: str) -> dict:
        g_vals = [v[metric_name] for v in neutral_genuine.values()]
        a_vals = [v[metric_name] for v in neutral_artefact.values()]
        if not g_vals or not a_vals:
            return {}
        print(f"\n  {metric_name}:")
        print(f"    Genuine  : mean={np.mean(g_vals):.4f} "
              f"std={np.std(g_vals):.4f} "
              f"min={np.min(g_vals):.4f} "
              f"max={np.max(g_vals):.4f}")
        print(f"    Artefact : mean={np.mean(a_vals):.4f} "
              f"std={np.std(a_vals):.4f} "
              f"min={np.min(a_vals):.4f} "
              f"max={np.max(a_vals):.4f}")
        separation = abs(float(np.mean(a_vals)) -
                         float(np.mean(g_vals)))
        print(f"    Separation: {separation:.4f}")
        return {
            "genuine_mean":  float(np.mean(g_vals)),
            "genuine_std":   float(np.std(g_vals)),
            "artefact_mean": float(np.mean(a_vals)),
            "artefact_std":  float(np.std(a_vals)),
            "separation":    separation,
        }

    discrimination_report = {
        "posterior_variance":  report_separation("posterior_variance"),
        "taxon_llr_variance":  report_separation("taxon_llr_variance"),
        "prop_neg_taxon_hits": report_separation("prop_neg_taxon_hits"),
        "prop_no_hits":        report_separation("prop_no_hits"),
        "mean_posterior":      report_separation("mean_posterior"),
    }

    # ── Parameter sweep ───────────────────────────────────────────────────────
    total_combos = (len(prior_prob_range) *
                    len(pident_threshold_range) *
                    len(k_range) *
                    len(beta_range) *
                    len(prop_neg_threshold_range) *
                    len(posterior_var_threshold_range) *
                    len(taxon_llr_var_threshold_range) *
                    len(no_hits_threshold_range)*
                    len(max_poorly_scoring_range))

    print(f"\n{'='*60}")
    print(f"SWEEPING {total_combos} PARAMETER COMBINATIONS")
    print(f"{'='*60}")

    best_f1          = -1.0
    best_params      = {}
    best_sensitivity = 0.0
    best_specificity = 0.0
    sweep_results    = []

    orphan_for_genuine = (
        known_orphan_taxa
        if orphan_rescue_mode == "genuine_only"
        else None
    )

    for prior in prior_prob_range:
        for pident_thresh in pident_threshold_range:
            for k in k_range:
                for beta in beta_range:

                    g_metrics_base = compute_parent_metrics(
                        genuine_by_parent, "GENUINE",
                        prior, pident_thresh, k, beta
                    )
                    a_metrics_base = compute_parent_metrics(
                        artefact_by_parent, "ARTEFACT",
                        prior, pident_thresh, k, beta
                    )

                    for max_poorly in max_poorly_scoring_range:

                        g_metrics = apply_cross_gene_consistency(
                            g_metrics_base,
                            min_well_scoring_fraction       =
                                min_well_scoring_fraction,
                            well_scoring_threshold_prop_neg =
                                well_scoring_threshold_prop_neg,
                            max_poorly_scoring_for_rescue   =
                                max_poorly,
                            known_orphan_taxa               =
                                orphan_for_genuine,
                            orphan_prop_neg_threshold       =
                                orphan_prop_neg_threshold,
                        )
                        a_metrics = apply_cross_gene_consistency(
                            a_metrics_base,
                            min_well_scoring_fraction       =
                                min_well_scoring_fraction,
                            well_scoring_threshold_prop_neg =
                                well_scoring_threshold_prop_neg,
                            max_poorly_scoring_for_rescue   =
                                max_poorly,
                            known_orphan_taxa               = None,
                            orphan_prop_neg_threshold       =
                                orphan_prop_neg_threshold,
                        )

                        for prop_neg_thresh in prop_neg_threshold_range:
                            for post_var_thresh in \
                                    posterior_var_threshold_range:
                                for llr_var_thresh in \
                                        taxon_llr_var_threshold_range:
                                    for no_hits_thresh in \
                                            no_hits_threshold_range:

                                        tp = tn = fp = fn = 0

                                        for m in g_metrics.values():
                                            if m.get("cross_gene_rescue",
                                                    False):
                                                predicted = False
                                            else:
                                                predicted = (
                                                    m["posterior_variance"]
                                                    > post_var_thresh
                                                    or m["taxon_llr_variance"]
                                                    > llr_var_thresh
                                                    or m["prop_neg_taxon_hits"]
                                                    > prop_neg_thresh
                                                    or m["single_hit_flag"]
                                                    or m["prop_no_hits"]
                                                    > no_hits_thresh
                                                )
                                            if predicted: fp += 1
                                            else:         tn += 1

                                        for m in a_metrics.values():
                                            if m.get("cross_gene_rescue",
                                                    False):
                                                predicted = False
                                            else:
                                                predicted = (
                                                    m["posterior_variance"]
                                                    > post_var_thresh
                                                    or m["taxon_llr_variance"]
                                                    > llr_var_thresh
                                                    or m["prop_neg_taxon_hits"]
                                                    > prop_neg_thresh
                                                    or m["single_hit_flag"]
                                                    or m["prop_no_hits"]
                                                    > no_hits_thresh
                                                )
                                            if predicted: tp += 1
                                            else:         fn += 1

                                        sensitivity = (tp / (tp + fn)
                                                    if (tp + fn) > 0
                                                    else 0.0)
                                        specificity = (tn / (tn + fp)
                                                    if (tn + fp) > 0
                                                    else 0.0)
                                        precision   = (tp / (tp + fp)
                                                    if (tp + fp) > 0
                                                    else 0.0)
                                        f1          = (
                                            2 * precision * sensitivity /
                                            (precision + sensitivity)
                                            if (precision + sensitivity) > 0
                                            else 0.0
                                        )

                                        combo = {
                                            "prior_prob":
                                                prior,
                                            "pident_threshold":
                                                pident_thresh,
                                            "k":
                                                k,
                                            "beta":
                                                beta,
                                            "prop_neg_threshold":
                                                prop_neg_thresh,
                                            "posterior_var_threshold":
                                                post_var_thresh,
                                            "taxon_llr_var_threshold":
                                                llr_var_thresh,
                                            "no_hits_threshold":
                                                no_hits_thresh,
                                            "max_poorly_scoring_for_rescue":
                                                max_poorly,
                                            "tp":          tp,
                                            "tn":          tn,
                                            "fp":          fp,
                                            "fn":          fn,
                                            "sensitivity":
                                                round(sensitivity, 4),
                                            "specificity":
                                                round(specificity, 4),
                                            "precision":
                                                round(precision, 4),
                                            "f1":
                                                round(f1, 4),
                                        }
                                        sweep_results.append(combo)

                                        if (f1 > best_f1
                                                or (f1 == best_f1
                                                    and specificity >
                                                    best_specificity)):
                                            best_f1          = f1
                                            best_sensitivity = sensitivity
                                            best_specificity = specificity
                                            best_params      = dict(combo)

    # ── Gene-specific prop_neg thresholds ─────────────────────────────────────
    print(f"\nCalibrating gene-specific prop_neg thresholds...")
    gene_prop_neg_thresholds = calibrate_gene_prop_neg_thresholds(
        genuine_by_parent  = genuine_by_parent,
        artefact_by_parent = artefact_by_parent,
        prior              = best_params["prior_prob"],
        pident_thresh      = best_params["pident_threshold"],
        k                  = best_params["k"],
        beta               = best_params["beta"],
        compute_metrics_fn = compute_parent_metrics,
    )

    # ── Orphan taxa report ────────────────────────────────────────────────────
    orphan_set = set(known_orphan_taxa or [])

    def get_genus(taxon: str) -> str:
        return taxon.split("_")[0] if taxon else ""

    orphan_report = {
        "known_orphan_taxa":               list(orphan_set),
        "orphan_genuine_would_be_flagged": 0,
        "orphan_artefact_detected":        0,
        "orphan_genuine_rescued":          0,
        "note": (
            "Orphan rescue is applied to GENUINE sequences only. "
            "Artefact sequences from orphan taxa are never rescued — "
            "chimeric insertions are genuine errors regardless of "
            "database coverage for the source taxon."
        ),
    }

    if orphan_set:
        g_metrics_final = compute_parent_metrics(
            genuine_by_parent, "GENUINE",
            prior         = best_params["prior_prob"],
            pident_thresh = best_params["pident_threshold"],
            k             = best_params["k"],
            beta          = best_params["beta"],
        )
        a_metrics_final = compute_parent_metrics(
            artefact_by_parent, "ARTEFACT",
            prior         = best_params["prior_prob"],
            pident_thresh = best_params["pident_threshold"],
            k             = best_params["k"],
            beta          = best_params["beta"],
        )
        g_metrics_rescued = apply_cross_gene_consistency(
            g_metrics_final,
            min_well_scoring_fraction       = min_well_scoring_fraction,
            well_scoring_threshold_prop_neg = well_scoring_threshold_prop_neg,
            max_poorly_scoring_for_rescue   =
                best_params["max_poorly_scoring_for_rescue"],
            known_orphan_taxa               = known_orphan_taxa,
            orphan_prop_neg_threshold       = orphan_prop_neg_threshold,
        )

        pnt = best_params["prop_neg_threshold"]
        pvt = best_params["posterior_var_threshold"]
        lvt = best_params["taxon_llr_var_threshold"]
        nht = best_params["no_hits_threshold"]

        for (acc, gene), m in g_metrics_final.items():
            if get_genus(m.get("taxon", "")) not in orphan_set:
                continue
            flagged_without = (
                m["posterior_variance"]  > pvt
                or m["taxon_llr_variance"] > lvt
                or m["prop_neg_taxon_hits"] > pnt
                or m["single_hit_flag"]
                or m["prop_no_hits"] > nht
            )
            if flagged_without:
                orphan_report["orphan_genuine_would_be_flagged"] += 1
                if g_metrics_rescued.get(
                        (acc, gene), {}).get(
                        "cross_gene_rescue", False):
                    orphan_report["orphan_genuine_rescued"] += 1

        for (acc, gene), m in a_metrics_final.items():
            if get_genus(m.get("taxon", "")) not in orphan_set:
                continue
            if (m["posterior_variance"]  > pvt
                    or m["taxon_llr_variance"] > lvt
                    or m["prop_neg_taxon_hits"] > pnt
                    or m["single_hit_flag"]
                    or m["prop_no_hits"] > nht):
                orphan_report["orphan_artefact_detected"] += 1

        print(f"\n── Orphan taxa report ───────────────────────────────────────")
        print(f"  Orphan taxa                    : "
              f"{', '.join(sorted(orphan_set))}")
        print(f"  Genuine flagged without rescue : "
              f"{orphan_report['orphan_genuine_would_be_flagged']}")
        print(f"  Genuine rescued by orphan logic: "
              f"{orphan_report['orphan_genuine_rescued']}")
        print(f"  Artefact detected              : "
              f"{orphan_report['orphan_artefact_detected']}")

    # ── Mode comparison at fixed best params ──────────────────────────────────
    mode_comparison = {}
    print(f"\n── Orphan rescue mode comparison at best params ─────────────")
    print(f"  {'Mode':<16} {'TP':>4} {'TN':>4} {'FP':>4} {'FN':>4} "
          f"{'F1':>7} {'Sens':>7} {'Spec':>7}")
    print(f"  {'-'*60}")

    pnt = best_params["prop_neg_threshold"]
    pvt = best_params["posterior_var_threshold"]
    lvt = best_params["taxon_llr_var_threshold"]
    nht = best_params["no_hits_threshold"]

    for mode in ("genuine_only", "none"):
        of_g = known_orphan_taxa if mode == "genuine_only" else None

        gm = compute_parent_metrics(
            genuine_by_parent, "GENUINE",
            best_params["prior_prob"],
            best_params["pident_threshold"],
            best_params["k"],
            best_params["beta"],
        )
        am = compute_parent_metrics(
            artefact_by_parent, "ARTEFACT",
            best_params["prior_prob"],
            best_params["pident_threshold"],
            best_params["k"],
            best_params["beta"],
        )
        gm = apply_cross_gene_consistency(
            gm,
            min_well_scoring_fraction       = min_well_scoring_fraction,
            well_scoring_threshold_prop_neg = well_scoring_threshold_prop_neg,
            max_poorly_scoring_for_rescue   =
                best_params["max_poorly_scoring_for_rescue"],
            known_orphan_taxa               = of_g,
            orphan_prop_neg_threshold       = orphan_prop_neg_threshold,
        )
        am = apply_cross_gene_consistency(
            am,
            min_well_scoring_fraction       = min_well_scoring_fraction,
            well_scoring_threshold_prop_neg = well_scoring_threshold_prop_neg,
            max_poorly_scoring_for_rescue   =
                best_params["max_poorly_scoring_for_rescue"],
            known_orphan_taxa               = None,
            orphan_prop_neg_threshold       = orphan_prop_neg_threshold,
        )
        tp = tn = fp = fn = 0

        for m in gm.values():
            pred = (not m.get("cross_gene_rescue", False) and (
                m["posterior_variance"]  > pvt
                or m["taxon_llr_variance"] > lvt
                or m["prop_neg_taxon_hits"] > pnt
                or m["single_hit_flag"]
                or m["prop_no_hits"] > nht
            ))
            if pred: fp += 1
            else:    tn += 1

        for m in am.values():
            pred = (not m.get("cross_gene_rescue", False) and (
                m["posterior_variance"]  > pvt
                or m["taxon_llr_variance"] > lvt
                or m["prop_neg_taxon_hits"] > pnt
                or m["single_hit_flag"]
                or m["prop_no_hits"] > nht
            ))
            if pred: tp += 1
            else:    fn += 1

        sens = tp/(tp+fn) if (tp+fn) > 0 else 0.0
        spec = tn/(tn+fp) if (tn+fp) > 0 else 0.0
        prec = tp/(tp+fp) if (tp+fp) > 0 else 0.0
        f1   = (2*prec*sens/(prec+sens)
                if (prec+sens) > 0 else 0.0)

        mode_comparison[mode] = {
            "tp": tp, "tn": tn, "fp": fp, "fn": fn,
            "sensitivity": round(sens, 4),
            "specificity": round(spec, 4),
            "precision":   round(prec, 4),
            "f1":          round(f1, 4),
        }
        print(f"  {mode:<16} {tp:>4} {tn:>4} {fp:>4} {fn:>4} "
              f"{f1:>7.4f} {sens:>7.4f} {spec:>7.4f}")

    # ── Save JSON ─────────────────────────────────────────────────────────────
    calibration_out = {
        "best_params":                    best_params,
        "best_f1":                        round(best_f1, 4),
        "best_sensitivity":               round(best_sensitivity, 4),
        "best_specificity":               round(best_specificity, 4),
        "n_genuine":                      len(genuine_parents),
        "n_artefact":                     len(artefact_parents),
        "discrimination_report":          discrimination_report,
        "flag_on_single_hit":             flag_on_single_hit,
        "single_hit_taxon_llr_threshold": single_hit_taxon_llr_threshold,
        "gene_prop_neg_thresholds":       gene_prop_neg_thresholds,
        "known_orphan_taxa":              list(known_orphan_taxa or []),
        "orphan_prop_neg_threshold":      orphan_prop_neg_threshold,
        "orphan_taxa_report":             orphan_report,
        "orphan_rescue_mode":             orphan_rescue_mode,
        "mode_comparison":                mode_comparison,
        "label_map":                      label_map,
        "genuine_labels":                 sorted(genuine_labels),
        "artefact_labels":                sorted(artefact_labels),
    }

    with open(output_path, "w") as f:
        json.dump(calibration_out, f, indent=2)

    # ── Save sweep CSV ────────────────────────────────────────────────────────
    sweep_csv = output_path.replace(".json", "_sweep.csv")
    if sweep_results:
        import csv as csv_module
        with open(sweep_csv, "w", newline="") as f:
            writer = csv_module.DictWriter(
                f, fieldnames=list(sweep_results[0].keys())
            )
            writer.writeheader()
            writer.writerows(sweep_results)

    # ── Print summary ─────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"CALIBRATION COMPLETE")
    print(f"{'='*60}")
    print(f"  Genuine sequences  : {len(genuine_parents)}")
    print(f"  Artefact sequences : {len(artefact_parents)}")
    print(f"  Best F1            : {best_f1:.4f}")
    print(f"  Best sensitivity   : {best_sensitivity:.4f}")
    print(f"  Best specificity   : {best_specificity:.4f}")
    print(f"\n  Best parameters:")
    skip = {"tp", "tn", "fp", "fn", "sensitivity",
            "specificity", "precision", "f1"}
    for pname, pval in best_params.items():
        if pname not in skip:
            print(f"    {pname:<30} : {pval}")
    print(f"\n  JSON saved : {output_path}")
    print(f"  CSV saved  : {sweep_csv}")

    return {
        "best_params":              best_params,
        "best_f1":                  round(best_f1, 4),
        "best_sensitivity":         round(best_sensitivity, 4),
        "best_specificity":         round(best_specificity, 4),
        "sweep_results":            sweep_results,
        "discrimination_report":    discrimination_report,
        "calibration_json":         output_path,
        "sweep_csv":                sweep_csv,
        "genuine_by_parent":        genuine_by_parent,
        "artefact_by_parent":       artefact_by_parent,
        "genuine_parents":          genuine_parents,
        "artefact_parents":         artefact_parents,
        "gene_prop_neg_thresholds": gene_prop_neg_thresholds,
        "compute_parent_metrics":   compute_parent_metrics,
        "mode_comparison":          mode_comparison,
    }



def run_calibration_diagnostics(
        cal: dict,
        known_orphan_taxa: list = None,
        orphan_prop_neg_threshold: float = 0.8,
        flagged_accessions: list = None,
) -> dict:
    """
    Run diagnostics on calibration results from calibrate_from_csv().

    Prints false positive analysis and per-accession cross-gene
    rescue status for specified accessions. Call immediately after
    calibrate_from_csv().

    Parameters
    ----------
    cal                      : dict   return value of calibrate_from_csv()
    known_orphan_taxa        : list   orphan genera for rescue
    orphan_prop_neg_threshold: float  permissive threshold for orphans
    flagged_accessions       : list   accessions to inspect in detail.
                                      Default: all unique accessions
                                      in false positives.

    Returns
    -------
    dict with:
        false_positives  : list   FP sequences with trigger info
        fp_by_trigger    : dict   trigger -> count
        fp_by_gene       : dict   gene -> count
        accession_details: dict   accession -> gene metrics
    """
    best               = cal["best_params"]
    compute_metrics    = cal["compute_parent_metrics"]
    genuine_by_parent  = cal["genuine_by_parent"]

    prop_neg_global    = best["prop_neg_threshold"]
    post_var_thresh    = best["posterior_var_threshold"]
    llr_var_thresh     = best["taxon_llr_var_threshold"]
    gene_thresholds    = cal.get("gene_prop_neg_thresholds", {})

    # Recompute with best params
    g_metrics = compute_metrics(
        genuine_by_parent, "GENUINE",
        prior         = best["prior_prob"],
        pident_thresh = best["pident_threshold"],
        k             = best["k"],
        beta          = best["beta"],
    )
    g_rescued = apply_cross_gene_consistency(
        g_metrics,
        known_orphan_taxa         = known_orphan_taxa,
        orphan_prop_neg_threshold = orphan_prop_neg_threshold,
    )

    # Find false positives
    false_positives = []
    for (acc, gene), m in g_rescued.items():
        prop_neg_thresh = gene_thresholds.get(gene, prop_neg_global)

        if m.get("cross_gene_rescue", False):
            predicted = False
        else:
            predicted = (
                m["posterior_variance"]  > post_var_thresh
                or m["taxon_llr_variance"] > llr_var_thresh
                or m["prop_neg_taxon_hits"] > prop_neg_thresh
                or m["single_hit_flag"]
            )

        if predicted:
            trigger = (
                "post_var"   if m["posterior_variance"] > post_var_thresh
                else "llr_var"  if m["taxon_llr_variance"] > llr_var_thresh
                else "prop_neg" if m["prop_neg_taxon_hits"] > prop_neg_thresh
                else "single_hit"
            )
            false_positives.append({
                "accession":  acc,
                "gene":       gene,
                "taxon":      m.get("taxon", ""),
                "prop_neg":   round(m["prop_neg_taxon_hits"], 3),
                "post_var":   round(m["posterior_variance"], 4),
                "llr_var":    round(m["taxon_llr_variance"], 4),
                "rescued":    m.get("cross_gene_rescue", False),
                "trigger":    trigger,
            })

    # Summary counts
    from collections import Counter
    fp_by_trigger = dict(Counter(fp["trigger"] for fp in false_positives))
    fp_by_gene    = dict(Counter(fp["gene"]    for fp in false_positives))

    print(f"\nFalse positives: {len(false_positives)}")
    print(f"{'Accession':<16} {'Gene':<10} {'Taxon':<30} "
          f"{'prop_neg':>9} {'trigger'}")
    print(f"  {'-'*75}")
    for fp in sorted(false_positives, key=lambda x: x["gene"]):
        print(f"  {fp['accession']:<14} {fp['gene']:<10} "
              f"{fp['taxon']:<30} "
              f"{fp['prop_neg']:>9.3f} {fp['trigger']}")

    print(f"\nBy trigger: {fp_by_trigger}")
    print(f"By gene   : {fp_by_gene}")

    # Per-accession detail for specified accessions
    if flagged_accessions is None:
        flagged_accessions = list({fp["accession"]
                                   for fp in false_positives})

    accession_details = {}
    for acc in sorted(flagged_accessions):
        genes = {gene: m for (a, gene), m in g_rescued.items()
                 if a == acc and _is_mito_gene(gene)}
        if not genes:
            genes = {gene: m for (a, gene), m in g_rescued.items()
                     if a == acc}
        if not genes:
            print(f"\n  {acc}: NOT FOUND")
            continue

        print(f"\n{acc} ({len(genes)} genes):")
        for gene, m in sorted(genes.items()):
            prop_neg_thresh = gene_thresholds.get(gene, prop_neg_global)
            if m.get("cross_gene_rescue", False):
                flagged = False
            else:
                flagged = (
                    m["posterior_variance"]  > post_var_thresh
                    or m["taxon_llr_variance"] > llr_var_thresh
                    or m["prop_neg_taxon_hits"] > prop_neg_thresh
                    or m["single_hit_flag"]
                )
            print(f"  {gene:<10} prop_neg={m['prop_neg_taxon_hits']:.3f} "
                  f"taxon={m.get('taxon',''):<26} "
                  f"flagged={flagged} "
                  f"rescued={m.get('cross_gene_rescue', False)}")
        accession_details[acc] = genes

    return {
        "false_positives":   false_positives,
        "fp_by_trigger":     fp_by_trigger,
        "fp_by_gene":        fp_by_gene,
        "accession_details": accession_details,
        "g_metrics_rescued": g_rescued,
    }

"""
    prior_prob_range              = (0.5, 0.7, 0.9, 0.95),
    pident_threshold_range        = (90.0, 92.0, 95.0),
    k_range                       = (0.05, 0.1, 0.2),
    beta_range                    = (0.005, 0.01, 0.05),
    prop_neg_threshold_range      = (0.4, 0.5, 0.6),
    posterior_var_threshold_range = (0.03, 0.05, 0.07, 0.10),
    taxon_llr_var_threshold_range = (1.5, 2.0, 2.5, 3.0, 5.0),

"""
"""
cal = calibrate_from_csv(
    results_csv         = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
    output_path         = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/model_calibration.json",
    # Narrow ranges around best values
    prior_prob_range              = (0.5,),
    pident_threshold_range        = (95.0,),
    k_range                       = (0.05,),
    beta_range                    = (0.005,),
    prop_neg_threshold_range      = (0.4, 0.5, 0.6),
    posterior_var_threshold_range = (0.02, 0.03, 0.05),
    taxon_llr_var_threshold_range = (3.0, 5.0),
    no_hits_threshold_range       = (0.3, 0.5, 0.9),
    max_poorly_scoring_range      = (1, 2, 3, 4, 5),

    known_orphan_taxa             = ["Uropsilus", "Urotrichus",
                                    "Solenodon", "Condylura", 
                                    "Desmana", "Galemys", "Dymecodon", 
                                    "Nectogale", "Scapanulus", 
                                    "Parascaptor", "Podogymnura", 
                                    "Neohylomys", "Congosorex", 
                                    "Diplomesodon", "Scaptochirus", 
                                    "Chimarrogale", "Episoriculus", 
                                    "Neotetracus", "Soriculus", 
                                    "Scaptonyx", "Euroscaptor", 
                                    "Blarinella", "Anourosorex", 
                                    "Sylvisorex", "Myosorex"],
    orphan_prop_neg_threshold     = 0.8,
    orphan_rescue_mode= "genuine_only"
)
"""



"""
Potential other orphan taxa:
Chimarrogale, Episoriculus, Neotetracus, Soriculus, Scaptonyx, Euroscaptor, Blarinella, Anourosorex, Sylvisorex, Myosorex
"""
#make json files for all types of calibration method
#genuine only
#cal = calibrate_from_csv(
#    results_csv         = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
#    output_path         = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/model_calibration_genuineonly.json",
#    # Narrow ranges around best values
#    prior_prob_range              = (0.5, 0.7, 0.9, 0.95),
#    pident_threshold_range        = (90.0, 92.0, 95.0),
#    k_range                       = (0.05, 0.1, 0.2),
#    beta_range                    = (0.005, 0.01, 0.05),
#    prop_neg_threshold_range      = (0.4, 0.5, 0.6),
#    posterior_var_threshold_range = (0.03, 0.05, 0.07, 0.10),
#    taxon_llr_var_threshold_range = (1.5, 2.0, 2.5, 3.0, 5.0),
#    known_orphan_taxa             = ["Uropsilus", "Urotrichus",
#                                     "Solenodon", "Condylura"],
#    orphan_prop_neg_threshold     = 0.8,
#    orphan_rescue_mode="genuine_only",
#)

#neither
#cal = calibrate_from_csv(
#    results_csv         = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
#    output_path         = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/model_calibration_neither.json",
#    # Narrow ranges around best values
#    prior_prob_range              = (0.5, 0.7, 0.9, 0.95),
#    pident_threshold_range        = (90.0, 92.0, 95.0),
#    k_range                       = (0.05, 0.1, 0.2),
#    beta_range                    = (0.005, 0.01, 0.05),
#    prop_neg_threshold_range      = (0.4, 0.5, 0.6),
#    posterior_var_threshold_range = (0.03, 0.05, 0.07, 0.10),
#    taxon_llr_var_threshold_range = (1.5, 2.0, 2.5, 3.0, 5.0),
#    known_orphan_taxa             = ["Uropsilus", "Urotrichus",
#                                     "Solenodon", "Condylura"],
#    orphan_prop_neg_threshold     = 0.8,
#    orphan_rescue_mode="None",
#)
"""

#run calibration on standard cal
diag = run_calibration_diagnostics(
    cal                      = cal,
    known_orphan_taxa        = ["Uropsilus", "Urotrichus",
                                "Solenodon", "Condylura"],
    orphan_prop_neg_threshold = 0.8,
    flagged_accessions       = ["NC_002080", "NC_002391",
                                "NC_005034", "NC_023244"],
)

print(f"\nFP by trigger: {diag['fp_by_trigger']}")
print(f"FP by gene   : {diag['fp_by_gene']}")

"""