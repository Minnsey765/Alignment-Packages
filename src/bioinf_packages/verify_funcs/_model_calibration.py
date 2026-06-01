# _model_calibration.py
import numpy as np
import json
import re
from scipy.stats import beta as beta_dist, norm as norm_dist
from Bio import SeqIO

from _score_hit import score_blast
from _species_parser import load_taxonomy
from _gene_parser import load_glossary


# ─────────────────────────────────────────────────────────────────────────────
# DISTRIBUTION FITTING HELPERS
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

def calibrate_priors(genuine_blast_results: list,
                     artefact_blast_results: list,
                     save_path: str = None) -> dict:
    """
    Fit empirical likelihood distributions from labelled score_blast()
    outputs and return calibrated LLR functions for each signal.

    Parameters
    ----------
    genuine_blast_results  : list   score_blast() outputs for genuine seqs
    artefact_blast_results : list   score_blast() outputs for artefact seqs
    save_path              : str    optional path to save JSON params

    Returns
    -------
    dict with keys:
        "params"      : fitted distribution parameters per feature
        "llr_funcs"   : callable LLR functions per feature
        "prior_prob"  : empirical prior P(genuine)
        "summary"     : mean/std per feature per class
    """

    def extract_features(blast_results):
        features = {
            "pident":              [],
            "bitscore":            [],
            "llr_taxon":           [],
            "llr_gene":            [],
            "prop_neg_taxon_hits": [],
            "taxon_llr_variance":  [],
        }
        for result in blast_results:
            # Hit-level features — one value per hit
            for hit in result["hit_details"]:
                if hit.get("gated_out", False):
                    continue
                features["pident"].append(hit["pident"] / 100)
                features["bitscore"].append(hit["bitscore"])
                features["llr_taxon"].append(hit["llr_taxon"])
                features["llr_gene"].append(hit["llr_gene"])

            # Sequence-level features — one value per sequence
            features["prop_neg_taxon_hits"].append(
                result.get("prop_neg_taxon_hits", 0.0)
            )
            features["taxon_llr_variance"].append(
                result.get("taxon_llr_variance", 0.0)
            )

        return {k: np.array(v) for k, v in features.items()}

    genuine_features  = extract_features(genuine_blast_results)
    artefact_features = extract_features(artefact_blast_results)

    # ── Empirical prior ───────────────────────────────────────────────────────
    n_genuine  = len(genuine_blast_results)
    n_artefact = len(artefact_blast_results)

    if n_genuine + n_artefact == 0:
        raise ValueError("Both result lists are empty — cannot calibrate.")

    prior_prob = n_genuine / (n_genuine + n_artefact)
    print(f"Genuine sequences:  {n_genuine}")
    print(f"Artefact sequences: {n_artefact}")
    print(f"Empirical prior P(genuine): {prior_prob:.3f}")

    params    = {}
    llr_funcs = {}
    summary   = {}

    def clip_beta(arr):
        return np.clip(arr, 1e-6, 1 - 1e-6)

    # ── pident (Beta) ─────────────────────────────────────────────────────────
    g_pident = clip_beta(genuine_features["pident"])
    a_pident = clip_beta(artefact_features["pident"])

    ag, bg = fit_beta_safe(g_pident, "pident (genuine)")
    aa, ba = fit_beta_safe(a_pident, "pident (artefact)")

    params["pident"] = {"genuine":  {"a": ag, "b": bg},
                        "artefact": {"a": aa, "b": ba}}

    def llr_pident(pident_raw):
        x = np.clip(pident_raw / 100, 1e-6, 1 - 1e-6)
        return float(np.log(beta_dist.pdf(x, ag, bg) /
                            beta_dist.pdf(x, aa, ba)))

    llr_funcs["pident"] = llr_pident

    summary["pident"] = {
        "genuine_mean":  float(np.mean(genuine_features["pident"]) * 100),
        "genuine_std":   float(np.std(genuine_features["pident"])  * 100),
        "artefact_mean": float(np.mean(artefact_features["pident"]) * 100),
        "artefact_std":  float(np.std(artefact_features["pident"])  * 100),
    }

    # ── bitscore (Beta) ───────────────────────────────────────────────────────
    all_bitscores = np.concatenate([genuine_features["bitscore"],
                                    artefact_features["bitscore"]])
    bitscore_max  = float(np.max(all_bitscores))

    g_bits = clip_beta(genuine_features["bitscore"]  / bitscore_max)
    a_bits = clip_beta(artefact_features["bitscore"] / bitscore_max)

    agb, bgb = fit_beta_safe(g_bits, "bitscore (genuine)")
    aab, bab = fit_beta_safe(a_bits, "bitscore (artefact)")

    params["bitscore"] = {"genuine":      {"a": agb, "b": bgb},
                          "artefact":     {"a": aab, "b": bab},
                          "max_observed": bitscore_max}

    def llr_bitscore(bitscore_raw):
        x = np.clip(bitscore_raw / bitscore_max, 1e-6, 1 - 1e-6)
        return float(np.log(beta_dist.pdf(x, agb, bgb) /
                            beta_dist.pdf(x, aab, bab)))

    llr_funcs["bitscore"] = llr_bitscore

    summary["bitscore"] = {
        "genuine_mean":  float(np.mean(genuine_features["bitscore"])),
        "genuine_std":   float(np.std(genuine_features["bitscore"])),
        "artefact_mean": float(np.mean(artefact_features["bitscore"])),
        "artefact_std":  float(np.std(artefact_features["bitscore"])),
    }

    # ── llr_taxon (Gaussian) ──────────────────────────────────────────────────
    mg,  sg  = fit_gaussian_safe(genuine_features["llr_taxon"],
                                  "llr_taxon (genuine)")
    ma,  sa  = fit_gaussian_safe(artefact_features["llr_taxon"],
                                  "llr_taxon (artefact)")

    params["llr_taxon"] = {"genuine":  {"mean": mg, "std": sg},
                           "artefact": {"mean": ma, "std": sa}}

    def llr_taxon_func(llr_taxon_val):
        g_pdf = norm_dist.pdf(llr_taxon_val, mg, sg)
        a_pdf = norm_dist.pdf(llr_taxon_val, ma, sa)
        if a_pdf == 0:
            return 10.0
        return float(np.log(g_pdf / a_pdf))

    llr_funcs["llr_taxon"] = llr_taxon_func

    summary["llr_taxon"] = {
        "genuine_mean":  mg,  "genuine_std":  sg,
        "artefact_mean": ma,  "artefact_std": sa,
    }

    # ── llr_gene (Gaussian) ───────────────────────────────────────────────────
    mg2, sg2 = fit_gaussian_safe(genuine_features["llr_gene"],
                                  "llr_gene (genuine)")
    ma2, sa2 = fit_gaussian_safe(artefact_features["llr_gene"],
                                  "llr_gene (artefact)")

    params["llr_gene"] = {"genuine":  {"mean": mg2, "std": sg2},
                          "artefact": {"mean": ma2, "std": sa2}}

    def llr_gene_func(llr_gene_val):
        g_pdf = norm_dist.pdf(llr_gene_val, mg2, sg2)
        a_pdf = norm_dist.pdf(llr_gene_val, ma2, sa2)
        if a_pdf == 0:
            return 10.0
        return float(np.log(g_pdf / a_pdf))

    llr_funcs["llr_gene"] = llr_gene_func

    summary["llr_gene"] = {
        "genuine_mean":  mg2,  "genuine_std":  sg2,
        "artefact_mean": ma2,  "artefact_std": sa2,
    }

    # ── prop_neg_taxon_hits (Beta — sequence level) ───────────────────────────
    g_prop = clip_beta(genuine_features["prop_neg_taxon_hits"])
    a_prop = clip_beta(artefact_features["prop_neg_taxon_hits"])

    agp, bgp = fit_beta_safe(g_prop, "prop_neg_taxon (genuine)")
    aap, bap = fit_beta_safe(a_prop, "prop_neg_taxon (artefact)")

    params["prop_neg_taxon_hits"] = {
        "genuine":  {"a": agp, "b": bgp},
        "artefact": {"a": aap, "b": bap},
    }

    def llr_prop_neg_taxon(prop):
        x = np.clip(prop, 1e-6, 1 - 1e-6)
        return float(np.log(beta_dist.pdf(x, agp, bgp) /
                            beta_dist.pdf(x, aap, bap)))

    llr_funcs["prop_neg_taxon_hits"] = llr_prop_neg_taxon

    summary["prop_neg_taxon_hits"] = {
        "genuine_mean":  float(np.mean(genuine_features["prop_neg_taxon_hits"])),
        "genuine_std":   float(np.std(genuine_features["prop_neg_taxon_hits"])),
        "artefact_mean": float(np.mean(artefact_features["prop_neg_taxon_hits"])),
        "artefact_std":  float(np.std(artefact_features["prop_neg_taxon_hits"])),
    }

    # ── taxon_llr_variance (Gaussian — sequence level) ────────────────────────
    mg3, sg3 = fit_gaussian_safe(genuine_features["taxon_llr_variance"],
                                  "taxon_llr_variance (genuine)")
    ma3, sa3 = fit_gaussian_safe(artefact_features["taxon_llr_variance"],
                                  "taxon_llr_variance (artefact)")

    params["taxon_llr_variance"] = {"genuine":  {"mean": mg3, "std": sg3},
                                    "artefact": {"mean": ma3, "std": sa3}}

    def llr_taxon_variance_func(var_val):
        g_pdf = norm_dist.pdf(var_val, mg3, sg3)
        a_pdf = norm_dist.pdf(var_val, ma3, sa3)
        if a_pdf == 0:
            return 10.0
        return float(np.log(g_pdf / a_pdf))

    llr_funcs["taxon_llr_variance"] = llr_taxon_variance_func

    summary["taxon_llr_variance"] = {
        "genuine_mean":  mg3,  "genuine_std":  sg3,
        "artefact_mean": ma3,  "artefact_std": sa3,
    }

    # ── Save to JSON ──────────────────────────────────────────────────────────
    if save_path:
        with open(save_path, "w") as f:
            json.dump({"params": params, "prior_prob": prior_prob,
                       "summary": summary}, f, indent=2)
        print(f"\nCalibration parameters saved to {save_path}")

    # ── Print summary ─────────────────────────────────────────────────────────
    print("\n── Calibration summary ──────────────────────────────────────")
    for feature, stats in summary.items():
        print(f"\n{feature}:")
        print(f"  Genuine:  mean={stats['genuine_mean']:.3f}  "
              f"std={stats['genuine_std']:.3f}")
        print(f"  Artefact: mean={stats['artefact_mean']:.3f}  "
              f"std={stats['artefact_std']:.3f}")

    return {
        "params":     params,
        "llr_funcs":  llr_funcs,
        "prior_prob": prior_prob,
        "summary":    summary,
    }


# ─────────────────────────────────────────────────────────────────────────────
# RUN CALIBRATION
# ─────────────────────────────────────────────────────────────────────────────

def run_calibration(genuine_fasta: str,
                    artefact_fasta: str,
                    taxonomy: dict,
                    glossary: dict,
                    save_path: str = None,
                    prior_prob: float = 0.95,
                    pident_threshold: float = 97.0,
                    k: float = 0.3,
                    beta: float = 0.01,
                    insert_fraction: float = 0.30,
                    target_detection_prob: float = 0.90) -> dict:
    """
    Parse two labelled FASTA files, run score_blast() on every sequence,
    diagnose feature distributions, and pass results to calibrate_priors().

    Note: subsample length (n) and number of samples are now calculated
    automatically inside score_blast() based on sequence length, so they
    are no longer parameters here.

    FASTA header format (pipe-delimited):
        >Species_name|gene_symbol|orien:+/-|accession:XXXXXXX

    Parameters
    ----------
    genuine_fasta         : str    path to FASTA of confirmed genuine seqs
    artefact_fasta        : str    path to FASTA of confirmed artefact seqs
    taxonomy              : dict   from load_taxonomy()
    glossary              : dict   from load_glossary()
    save_path             : str    optional path to save calibration JSON
    prior_prob            : float  prior P(genuine) passed to score_blast()
    pident_threshold      : float  pident threshold passed to score_blast()
    k                     : float  pident LLR steepness
    beta                  : float  softmax temperature
    insert_fraction       : float  expected insert fraction for n_samples
                                   calculation (default 0.30)
    target_detection_prob : float  desired P(detect chimera) (default 0.90)

    Returns
    -------
    dict   output of calibrate_priors()
    """

    def parse_header(header: str) -> dict:
        """
        Parse a pipe-delimited FASTA header into its components.
        '>Talpa_europaea|12S_rRNA|orien:+|accession:NC_002391'
        """
        parts     = header.strip().split("|")
        species   = parts[0].strip()
        gene      = parts[1].strip() if len(parts) > 1 else None
        accession = None
        for part in parts:
            if part.startswith("accession:"):
                accession = part.split(":", 1)[1].strip()
                break
        return {"species": species, "gene": gene, "accession": accession}

    def process_fasta(fasta_path: str, label: str) -> list:
        """
        Run score_blast() on every record in a FASTA file.
        Returns a list of score_blast() output dicts.
        """
        results = []
        records = list(SeqIO.parse(fasta_path, "fasta"))
        total   = len(records)

        print(f"\nProcessing {label} set: {total} sequences from {fasta_path}")

        for i, record in enumerate(records, 1):
            meta      = parse_header(record.description)
            species   = meta["species"]
            gene      = meta["gene"]
            accession = meta["accession"]
            seq       = str(record.seq)

            print(f"  [{i}/{total}] {species} | {gene} | {accession}")

            try:
                result = score_blast(
                    query_accession      = accession,
                    seq                  = seq,
                    taxon                = species,
                    gene                 = gene,
                    taxonomy             = taxonomy,
                    glossary             = glossary,
                    prior_prob           = prior_prob,
                    pident_threshold     = pident_threshold,
                    k                    = k,
                    beta                 = beta,
                    insert_fraction      = insert_fraction,
                    target_detection_prob= target_detection_prob,
                )
                result["_meta"] = meta
                results.append(result)

            except Exception as e:
                print(f"    ✗ Failed: {e} — skipping")
                import traceback
                traceback.print_exc()

        print(f"  Done: {len(results)}/{total} succeeded")
        return results

    def diagnose_calibration_data(genuine_results: list,
                                  artefact_results: list):
        """Print feature distribution summary before fitting."""

        def extract(results):
            pidents, bitscores, llr_taxons, llr_genes = [], [], [], []
            prop_negs, taxon_vars = [], []
            for r in results:
                for h in r["hit_details"]:
                    if h.get("gated_out", False):
                        continue
                    pidents.append(h["pident"])
                    bitscores.append(h["bitscore"])
                    llr_taxons.append(h["llr_taxon"])
                    llr_genes.append(h["llr_gene"])
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

        g = extract(genuine_results)
        a = extract(artefact_results)

        print(f"\n{'Feature':<24} {'Class':<10} {'n':<6} {'mean':<8} "
              f"{'std':<8} {'min':<8} {'max':<8}")
        print("─" * 76)

        features = ["pident", "bitscore", "llr_taxon", "llr_gene",
                    "prop_neg_taxon_hits", "taxon_llr_variance"]

        for feature in features:
            for label, data in [("genuine", g[feature]),
                                 ("artefact", a[feature])]:
                if len(data) == 0:
                    print(f"{feature:<24} {label:<10} NO DATA")
                    continue
                print(f"{feature:<24} {label:<10} {len(data):<6} "
                      f"{np.mean(data):<8.3f} {np.std(data):<8.3f} "
                      f"{np.min(data):<8.3f} {np.max(data):<8.3f}")
            print()

    # ── Run both sets ─────────────────────────────────────────────────────────
    genuine_results  = process_fasta(genuine_fasta,  label="genuine")
    artefact_results = process_fasta(artefact_fasta, label="artefact")

    # ── Diagnose before calibrating ───────────────────────────────────────────
    print("\nDiagnosing feature distributions...")
    diagnose_calibration_data(genuine_results, artefact_results)

    # ── Calibrate ─────────────────────────────────────────────────────────────
    calibration = calibrate_priors(genuine_results, artefact_results,
                                   save_path=save_path)

    return calibration





#usage
taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv")


#run callibration
calibration = run_calibration(
    genuine_fasta  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta",
    artefact_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_seq/modifed_seqs.fasta",
    taxonomy       = taxonomy,
    glossary       = glossary,
    save_path      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/callibration.json",
)