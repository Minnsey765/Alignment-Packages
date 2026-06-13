# _model_calibration.py

import os
import numpy as np
import json
import math
import csv

from pathlib import Path
from scipy.stats import beta as beta_dist
from Bio import SeqIO

from bioinf_packages.verify_funcs._score_hit import (
    score_blast,
    sigmoid,
    pident_to_llr,
    softmax_weights,
    recommended_subsample_length,
    recommended_n_samples,
)
from bioinf_packages.verify_funcs._species_parser import load_taxonomy
from bioinf_packages.verify_funcs._gene_parser import load_glossary


# ─────────────────────────────────────────────────────────────────────────────
# DISTRIBUTION FITTING HELPERS
# These are kept here for any downstream code that imports them from
# _model_calibration. They are not used by calibrate_priors() directly.
# ─────────────────────────────────────────────────────────────────────────────

def fit_beta_safe(data, feature_name):
    """
    Attempt Beta distribution fit with fallback strategies if the solver
    fails due to low variance or extreme values in the data.
    """
    data = np.clip(data, 1e-4, 1 - 1e-4)

    print(f"\n  {feature_name}: n={len(data)}, "
          f"mean={np.mean(data):.4f}, std={np.std(data):.4f}, "
          f"min={np.min(data):.4f}, max={np.max(data):.4f}")

    if len(data) < 5:
        print(f"  WARNING: Too few data points for {feature_name} "
              f"— using uniform prior (a=1, b=1)")
        return 1.0, 1.0

    if np.std(data) < 1e-4:
        print(f"  WARNING: Near-zero variance in {feature_name} "
              f"— using moment estimates")
        mu  = np.mean(data)
        var = max(np.var(data), 1e-6)
        a   = mu * (mu * (1 - mu) / var - 1)
        b   = (1 - mu) * (mu * (1 - mu) / var - 1)
        return max(a, 0.1), max(b, 0.1)

    try:
        a, b, _, _ = beta_dist.fit(data, floc=0, fscale=1)
        return a, b
    except Exception:
        print(f"  WARNING: MLE fit failed for {feature_name} "
              f"— falling back to method of moments")
        mu  = np.mean(data)
        var = max(np.var(data), 1e-6)
        a   = mu * (mu * (1 - mu) / var - 1)
        b   = (1 - mu) * (mu * (1 - mu) / var - 1)
        return max(a, 0.1), max(b, 0.1)


def fit_gaussian_safe(data, feature_name):
    """
    Fit a Gaussian with a fallback minimum std to prevent zero-variance
    issues when discrete LLR values are all identical in a class.
    """
    print(f"\n  {feature_name}: n={len(data)}, "
          f"mean={np.mean(data):.4f}, std={np.std(data):.4f}, "
          f"min={np.min(data):.4f}, max={np.max(data):.4f}")

    if len(data) < 5:
        print(f"  WARNING: Too few data points for {feature_name} "
              f"— using default mean=0, std=1")
        return 0.0, 1.0

    mean = float(np.mean(data))
    std  = float(np.std(data))

    if std < 1e-4:
        print(f"  WARNING: Near-zero variance in {feature_name} "
              f"— setting minimum std of 0.1")
        std = 0.1

    return mean, std


# ─────────────────────────────────────────────────────────────────────────────
# CALIBRATE PRIORS
# ─────────────────────────────────────────────────────────────────────────────

def calibrate_priors(genuine_fasta: str,
                     artefact_fasta: str,
                     taxonomy: dict,
                     glossary: dict,
                     output_path: str,
                     blast_db: str = None,
                     blast_bin: str = "blastn",
                     num_threads: int = 4,
                     evalue: float = 1e-10,
                     max_target_seqs: int = 20,
                     perc_identity: float = 70.0,
                     output_dir: str = None,
                     insert_fraction: float = 0.30,
                     target_detection_prob: float = 0.90,
                     prior_prob_range: tuple = (0.5, 0.6, 0.7,
                                                0.8, 0.9, 0.95),
                     pident_threshold_range: tuple = (85.0, 88.0,
                                                      90.0, 92.0,
                                                      95.0),
                     k_range: tuple = (0.1, 0.2, 0.5, 1.0),
                     beta_range: tuple = (0.01, 0.05, 0.1),
                     flag_on_single_hit: bool = True,
                     single_hit_taxon_llr_threshold: float = -2.0) -> dict:
    """
    Calibrate BLAST scoring priors using known genuine and artefactual
    sequences by sweeping combinations of prior_prob, pident_threshold,
    k, and beta to find the combination that maximises F1 score.

    Unlike the distance calibration which requires per-gene scaffold
    files and therefore limits sample size, this function uses all genes
    from both FASTA files simultaneously — BLAST finds its own hits from
    the database without needing a pre-built per-gene scaffold.

    The calibration runs BLAST once per sequence with neutral parameters
    to collect raw hit data (bitscore, pident, taxon LLR, gene LLR),
    then sweeps aggregation parameter combinations analytically without
    re-running BLAST. This makes the sweep fast regardless of how many
    parameter combinations are tested.

    FASTA header formats accepted
    ------------------------------
    Genuine  : >Genus_species|gene|orien:+|accession:XXXX
    Artefact : >Genus_species|gene|mod_N|orien:+|accession:XXXX

    Parameters
    ----------
    genuine_fasta                  : str   known genuine sequences
    artefact_fasta                 : str   known artefactual sequences
    taxonomy                       : dict  from load_taxonomy()
    glossary                       : dict  from load_glossary()
    output_path                    : str   path to save calibration JSON
    blast_db                       : str   path to local BLAST database.
                                           If None, uses remote BLAST via
                                           verify_seq() — not recommended
                                           for large calibration sets.
    blast_bin                      : str   path to blastn executable
    num_threads                    : int   BLAST CPU threads
    evalue                         : float BLAST e-value threshold
    max_target_seqs                : int   max hits per subsample
    perc_identity                  : float minimum percent identity
    output_dir                     : str   directory for BLAST XML output
                                           files. Required if blast_db set.
    insert_fraction                : float expected chimeric insert fraction
    target_detection_prob          : float desired P(detect chimera)
    prior_prob_range               : tuple prior_prob values to sweep
    pident_threshold_range         : tuple pident_threshold values to sweep
    k_range                        : tuple k (pident LLR steepness) values
    beta_range                     : tuple softmax temperature values
    flag_on_single_hit             : bool  apply single-hit flagging during
                                           sweep (recommended True)
    single_hit_taxon_llr_threshold : float taxon LLR below which a single
                                           subsample triggers chimera flag

    Returns
    -------
    dict with:
        best_params      : dict   best prior_prob, pident_threshold, k, beta
        best_f1          : float
        best_sensitivity : float
        best_specificity : float
        sweep_results    : list   all parameter combinations with metrics
        raw_results      : list   per-sequence score_blast() raw outputs
        calibration_json : str    path to saved JSON
        sweep_csv        : str    path to saved sweep CSV
    """
    if isinstance(taxonomy, str):
        raise TypeError(
            "taxonomy must be a dict from load_taxonomy(), not a path."
        )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

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

    # ── Step 1: collect raw BLAST scores with neutral parameters ──────────────
    # Run BLAST once per sequence using neutral priors so we store all
    # hit-level data. The parameter sweep then operates analytically on
    # these stored scores without re-running BLAST.
    print(f"\n{'='*60}")
    print(f"COLLECTING RAW BLAST SCORES FOR CALIBRATION")
    print(f"{'='*60}")

    raw_results = []

    for fasta_path, label in [(genuine_fasta,  "GENUINE"),
                               (artefact_fasta, "ARTEFACT")]:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        print(f"\nProcessing {label}: {len(records)} sequences")

        for i, record in enumerate(records, 1):
            info = parse_header(record.description)

            if not info["gene"] or not info["accession"]:
                print(f"  [{i}] Skipping — missing gene or accession "
                      f"in header: {record.description[:60]}")
                continue

            seq_dir = None
            if output_dir:
                seq_dir = os.path.join(
                    output_dir,
                    f"{label}_{info['taxon']}_{info['accession']}"
                )

            print(f"\n  [{i}/{len(records)}] "
                  f"{info['taxon']} | {info['gene']} | "
                  f"{info['accession']} [{label}]")

            try:
                result = score_blast(
                    query_accession       = info["accession"],
                    seq                   = str(record.seq),
                    taxon                 = info["taxon"],
                    gene                  = info["gene"],
                    taxonomy              = taxonomy,
                    glossary              = glossary,
                    prior_prob            = 0.5,    # neutral
                    pident_threshold      = 90.0,   # neutral
                    k                     = 0.2,    # neutral
                    beta                  = 0.05,   # neutral
                    insert_fraction       = insert_fraction,
                    target_detection_prob = target_detection_prob,
                    blast_db              = blast_db,
                    blast_bin             = blast_bin,
                    num_threads           = num_threads,
                    evalue                = evalue,
                    max_target_seqs       = max_target_seqs,
                    perc_identity         = perc_identity,
                    output_dir            = seq_dir,
                    flag_on_single_hit    = False,  # collect raw only
                )

                result["label"]       = label
                result["is_artefact"] = (label == "ARTEFACT")
                raw_results.append(result)

                print(f"    ✓ {len(result['sample_results'])} samples | "
                      f"{len(result['hit_details'])} total hits")

            except Exception as e:
                print(f"    ✗ Failed: {e}")
                import traceback
                traceback.print_exc()

    if not raw_results:
        raise ValueError("No sequences were successfully scored. "
                         "Check BLAST database path and input FASTA.")

    genuine_n  = sum(1 for r in raw_results if not r["is_artefact"])
    artefact_n = sum(1 for r in raw_results if r["is_artefact"])
    print(f"\nCollected: {genuine_n} genuine, {artefact_n} artefact")

    if genuine_n == 0 or artefact_n == 0:
        raise ValueError(
            f"Need both genuine and artefact sequences. "
            f"Got genuine={genuine_n}, artefact={artefact_n}."
        )

    # ── Step 2: print feature distribution diagnostics ────────────────────────
    print(f"\n{'='*60}")
    print(f"FEATURE DISTRIBUTIONS")
    print(f"{'='*60}")

    def extract_features(results: list) -> dict:
        pidents, bitscores, llr_taxons, llr_genes = [], [], [], []
        prop_negs, taxon_vars = [], []
        for r in results:
            for h in r.get("hit_details", []):
                if h.get("gated_out", False):
                    continue
                pidents.append(h.get("pident", 0))
                bitscores.append(h.get("bitscore", 0))
                llr_taxons.append(h.get("llr_taxon", 0))
                llr_genes.append(h.get("llr_gene", 0))
            prop_negs.append(r.get("prop_neg_taxon_hits", 0.0))
            taxon_vars.append(r.get("taxon_llr_variance", 0.0))
        return {
            "pident":              np.array(pidents),
            "bitscore":            np.array(bitscores),
            "llr_taxon":           np.array(llr_taxons),
            "llr_gene":            np.array(llr_genes),
            "prop_neg_taxon_hits": np.array(prop_negs),
            "taxon_llr_variance":  np.array(taxon_vars),
        }

    g_feat = extract_features(
        [r for r in raw_results if not r["is_artefact"]]
    )
    a_feat = extract_features(
        [r for r in raw_results if r["is_artefact"]]
    )

    print(f"\n{'Feature':<24} {'Class':<10} {'n':<6} {'mean':<8} "
          f"{'std':<8} {'min':<8} {'max':<8}")
    print("─" * 76)

    for feature in ["pident", "bitscore", "llr_taxon", "llr_gene",
                    "prop_neg_taxon_hits", "taxon_llr_variance"]:
        for lbl, feat in [("genuine", g_feat), ("artefact", a_feat)]:
            d = feat[feature]
            if len(d) == 0:
                print(f"{feature:<24} {lbl:<10} NO DATA")
            else:
                print(f"{feature:<24} {lbl:<10} {len(d):<6} "
                      f"{np.mean(d):<8.3f} {np.std(d):<8.3f} "
                      f"{np.min(d):<8.3f} {np.max(d):<8.3f}")
        print()

    # ── Step 3: sweep parameter combinations analytically ─────────────────────
    print(f"\n{'='*60}")
    print(f"SWEEPING PARAMETER COMBINATIONS")
    print(f"{'='*60}")

    total_combos = (len(prior_prob_range) *
                    len(pident_threshold_range) *
                    len(k_range) *
                    len(beta_range))
    print(f"Total combinations to sweep: {total_combos}")

    best_f1          = -1.0
    best_params      = {}
    best_sensitivity = 0.0
    best_specificity = 0.0
    sweep_results    = []

    for prior in prior_prob_range:
        for pident_thresh in pident_threshold_range:
            for k in k_range:
                for beta in beta_range:

                    tp = tn = fp = fn = 0

                    for seq_result in raw_results:

                        prior_log_odds      = math.log(
                            prior / (1 - prior)
                        )
                        cumulative_log_odds = prior_log_odds
                        any_single_hit_flag = False
                        sample_posts        = []
                        all_taxon_llrs      = []

                        for sample in seq_result.get("sample_results",
                                                     []):
                            hit_details = sample.get("hit_details", [])
                            if not hit_details:
                                continue

                            bitscores = [h["bitscore"]
                                         for h in hit_details]
                            weights   = softmax_weights(bitscores,
                                                        beta=beta)

                            sample_llr   = 0.0
                            top_llr      = None
                            top_bitscore = -1
                            s_log_odds   = math.log(prior / (1 - prior))

                            for h, w in zip(hit_details, weights):
                                llr_pident = pident_to_llr(
                                    h["pident"],
                                    threshold=pident_thresh,
                                    k=k
                                )
                                llr_taxon = h["llr_taxon"]
                                llr_gene  = h["llr_gene"]
                                combined  = w * (llr_taxon
                                                 + llr_gene
                                                 + llr_pident)
                                sample_llr += combined
                                s_log_odds += (llr_taxon
                                               + llr_gene
                                               + llr_pident)

                                all_taxon_llrs.append(llr_taxon)

                                if h["bitscore"] > top_bitscore:
                                    top_bitscore = h["bitscore"]
                                    top_llr      = llr_taxon

                            cumulative_log_odds += sample_llr
                            sample_posts.append(sigmoid(s_log_odds))

                            # Single-hit check
                            if (flag_on_single_hit
                                    and top_llr is not None
                                    and top_llr 
                                    single_hit_taxon_llr_threshold):
                                any_single_hit_flag = True

                        # Aggregate chimera flag
                        post_var = (float(np.var(sample_posts))
                                    if len(sample_posts) > 1 else 0.0)
                        llr_var  = (float(np.var(all_taxon_llrs))
                                    if len(all_taxon_llrs) > 1
                                    else 0.0)

                        agg_flag          = (post_var > 0.05
                                             or llr_var > 2.0)
                        predicted_chimera = (agg_flag
                                             or any_single_hit_flag)
                        actual_chimera    = seq_result["is_artefact"]

                        if predicted_chimera and actual_chimera:
                            tp += 1
                        elif not predicted_chimera and not actual_chimera:
                            tn += 1
                        elif predicted_chimera and not actual_chimera:
                            fp += 1
                        else:
                            fn += 1

                    sensitivity = (tp / (tp + fn)
                                   if (tp + fn) > 0 else 0.0)
                    specificity = (tn / (tn + fp)
                                   if (tn + fp) > 0 else 0.0)
                    precision   = (tp / (tp + fp)
                                   if (tp + fp) > 0 else 0.0)
                    f1          = (2 * precision * sensitivity
                                   / (precision + sensitivity)
                                   if (precision + sensitivity) > 0
                                   else 0.0)

                    combo = {
                        "prior_prob":       prior,
                        "pident_threshold": pident_thresh,
                        "k":                k,
                        "beta":             beta,
                        "tp":               tp,
                        "tn":               tn,
                        "fp":               fp,
                        "fn":               fn,
                        "sensitivity":      round(sensitivity, 4),
                        "specificity":      round(specificity, 4),
                        "precision":        round(precision, 4),
                        "f1":               round(f1, 4),
                    }
                    sweep_results.append(combo)

                    # Prefer higher F1; break ties by higher specificity
                    if (f1 > best_f1
                            or (f1 == best_f1
                                and specificity > best_specificity)):
                        best_f1          = f1
                        best_sensitivity = sensitivity
                        best_specificity = specificity
                        best_params      = {
                            "prior_prob":       prior,
                            "pident_threshold": pident_thresh,
                            "k":                k,
                            "beta":             beta,
                        }

    # ── Step 4: save results ──────────────────────────────────────────────────
    calibration_out = {
        "best_params":                  best_params,
        "best_f1":                      round(best_f1, 4),
        "best_sensitivity":             round(best_sensitivity, 4),
        "best_specificity":             round(best_specificity, 4),
        "n_genuine":                    genuine_n,
        "n_artefact":                   artefact_n,
        "flag_on_single_hit":           flag_on_single_hit,
        "single_hit_taxon_llr_threshold": single_hit_taxon_llr_threshold,
    }

    with open(output_path, "w") as f:
        json.dump(calibration_out, f, indent=2)

    sweep_csv = output_path.replace(".json", "_sweep.csv")
    if sweep_results:
        with open(sweep_csv, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=list(sweep_results[0].keys())
            )
            writer.writeheader()
            writer.writerows(sweep_results)

    # ── Step 5: print summary ─────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"CALIBRATION COMPLETE")
    print(f"{'='*60}")
    print(f"  Genuine sequences   : {genuine_n}")
    print(f"  Artefact sequences  : {artefact_n}")
    print(f"  Best parameters     :")
    for pname, pval in best_params.items():
        print(f"    {pname:<22} : {pval}")
    print(f"  Best F1             : {best_f1:.4f}")
    print(f"  Best sensitivity    : {best_sensitivity:.4f}")
    print(f"  Best specificity    : {best_specificity:.4f}")
    print(f"\n  JSON saved          : {output_path}")
    print(f"  Sweep CSV saved     : {sweep_csv}")
    print(f"\n  Use these in score_blast():")
    print(f"    score_blast(")
    for pname, pval in best_params.items():
        print(f"        {pname:<22} = {pval},")
    print(f"    )")

    return {
        "best_params":      best_params,
        "best_f1":          round(best_f1, 4),
        "best_sensitivity": round(best_sensitivity, 4),
        "best_specificity": round(best_specificity, 4),
        "sweep_results":    sweep_results,
        "raw_results":      raw_results,
        "calibration_json": output_path,
        "sweep_csv":        sweep_csv,
    }


# ─────────────────────────────────────────────────────────────────────────────
# RUN CALIBRATION  (convenience wrapper)
# ─────────────────────────────────────────────────────────────────────────────

def run_calibration(genuine_fasta: str,
                    artefact_fasta: str,
                    taxonomy: dict,
                    glossary: dict,
                    output_path: str,
                    blast_db: str,
                    blast_bin: str = "blastn",
                    num_threads: int = 4,
                    evalue: float = 1e-10,
                    max_target_seqs: int = 20,
                    perc_identity: float = 70.0,
                    output_dir: str = None,
                    insert_fraction: float = 0.30,
                    target_detection_prob: float = 0.90,
                    flag_on_single_hit: bool = True,
                    single_hit_taxon_llr_threshold: float = -2.0) -> dict:
    """
    Convenience wrapper around calibrate_priors() with sensible defaults.

    Runs the full calibration pipeline:
      1. BLAST all genuine and artefact sequences (once, with neutral
         parameters)
      2. Print feature distribution diagnostics
      3. Sweep parameter combinations analytically
      4. Save best parameters to JSON

    Parameters
    ----------
    genuine_fasta                  : str   known genuine sequences FASTA
    artefact_fasta                 : str   known artefactual sequences FASTA
    taxonomy                       : dict  from load_taxonomy()
    glossary                       : dict  from load_glossary()
    output_path                    : str   path to save calibration JSON
    blast_db                       : str   path to local BLAST database
    blast_bin                      : str   path to blastn executable
    num_threads                    : int   BLAST CPU threads
    evalue                         : float BLAST e-value threshold
    max_target_seqs                : int   max hits per subsample
    perc_identity                  : float minimum percent identity
    output_dir                     : str   directory for BLAST XML outputs
    insert_fraction                : float expected chimeric insert fraction
    target_detection_prob          : float desired P(detect chimera)
    flag_on_single_hit             : bool  apply single-hit flagging
    single_hit_taxon_llr_threshold : float single-hit taxon LLR threshold

    Returns
    -------
    dict   output of calibrate_priors()
    """
    return calibrate_priors(
        genuine_fasta                  = genuine_fasta,
        artefact_fasta                 = artefact_fasta,
        taxonomy                       = taxonomy,
        glossary                       = glossary,
        output_path                    = output_path,
        blast_db                       = blast_db,
        blast_bin                      = blast_bin,
        num_threads                    = num_threads,
        evalue                         = evalue,
        max_target_seqs                = max_target_seqs,
        perc_identity                  = perc_identity,
        output_dir                     = output_dir,
        insert_fraction                = insert_fraction,
        target_detection_prob          = target_detection_prob,
        flag_on_single_hit             = flag_on_single_hit,
        single_hit_taxon_llr_threshold = single_hit_taxon_llr_threshold,
    )




#usage
taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv")


print(run_calibration(
    genuine_fasta  = "gold_set/genuine_sequences.fasta",
    artefact_fasta = "gold_set/artefact_sequences.fasta",
    taxonomy       = taxonomy,
    glossary       = glossary,
    output_path    = "calibration_output/blast_priors.json",
    blast_db       = "C:/blast/db/core_nt",
    blast_bin      = "C:/blast/bin/blastn.exe",
    num_threads    = 8,
    output_dir     = "calibration_output/blast_runs",
))