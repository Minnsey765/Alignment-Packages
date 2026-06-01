# _score_hit.py
import csv
import re
import math
import numpy as np

from _gene_parser import gene_parser, load_glossary
from _species_parser import load_taxonomy, species_parser
from _verify_seq import verify_seq


# ─────────────────────────────────────────────────────────────────────────────
# INDIVIDUAL HIT SCORER
# ─────────────────────────────────────────────────────────────────────────────

def score_hit(hit: dict,
              species_name: str,
              expected_gene: str,
              taxonomy: dict,
              glossary: dict) -> dict:
    """
    Score a single BLAST hit dict (as returned by verify_seq) and return
    all component LLRs ready for Bayesian aggregation.

    Parameters
    ----------
    hit           : dict   one entry from verify_seq()["hits"]
    species_name  : str    e.g. "Sorex_araneus"
    expected_gene : str    e.g. "CYTB"
    taxonomy      : dict   from load_taxonomy()
    glossary      : dict   from load_glossary()

    Returns
    -------
    dict with keys:
        "accession"       : str
        "bitscore"        : float
        "pident"          : float
        "evalue"          : float
        "align_length"    : int
        "species_result"  : dict  (full output of species_parser)
        "gene_result"     : dict  (full output of gene_parser)
        "species_llr"     : float
        "gene_llr"        : float
        "total_llr"       : float
    """
    # Truncate at first ">" to strip concatenated secondary hits
    desc = hit["description"].split(">")[0].strip()

    sp_result   = species_parser(species_name, desc, taxonomy)
    gene_result = gene_parser(expected_gene, desc, glossary)

    total_llr = sp_result["llr"] + gene_result["llr"]

    return {
        "accession":      hit["accession"],
        "bitscore":       hit["bitscore"],
        "pident":         hit["pident"],
        "evalue":         hit["evalue"],
        "align_length":   hit["align_length"],
        "species_result": sp_result,
        "gene_result":    gene_result,
        "species_llr":    sp_result["llr"],
        "gene_llr":       gene_result["llr"],
        "total_llr":      total_llr,
    }


# ─────────────────────────────────────────────────────────────────────────────
# AGGREGATION HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def sigmoid(x):
    return 1 / (1 + math.exp(-x))


def pident_to_llr(pident, threshold, k):
    """Linear LLR centred on threshold. Positive above, negative below."""
    return k * (pident - threshold)


def softmax_weights(bitscores, beta):
    """Softmax over bitscores. Beta controls how much top hits dominate."""
    exp_scores = [math.exp(beta * b) for b in bitscores]
    total      = sum(exp_scores)
    return [e / total for e in exp_scores]


def aggregate_hit_scores(scored_hits, prior_prob, pident_threshold, k, beta):
    """
    Aggregate scored hits into a posterior probability of genuineness.

    Parameters
    ----------
    scored_hits       : list of dicts from score_hit()
    prior_prob        : float  prior P(genuine) before seeing BLAST results
    pident_threshold  : float  pident value that contributes LLR of 0
    k                 : float  steepness of pident->LLR mapping
    beta              : float  softmax temperature over bitscores

    Returns
    -------
    dict with per-hit breakdown and final posterior probability
    """
    prior_log_odds = math.log(prior_prob / (1 - prior_prob))

    bitscores = [h["bitscore"] for h in scored_hits]
    weights   = softmax_weights(bitscores, beta=beta)

    hit_details = []
    total_llr   = 0.0

    for hit, w in zip(scored_hits, weights):
        llr_pident   = pident_to_llr(hit["pident"],
                                     threshold=pident_threshold, k=k)
        llr_taxon    = hit["species_llr"]
        llr_gene     = hit["gene_llr"]
        combined_llr = w * (llr_taxon + llr_gene + llr_pident)
        total_llr   += combined_llr

        hit_details.append({
            "accession":    hit["accession"],
            "weight":       round(w, 4),
            "bitscore":     hit["bitscore"],
            "pident":       hit["pident"],
            "align_length": hit["align_length"],
            "llr_taxon":    round(llr_taxon, 4),
            "llr_gene":     round(llr_gene, 4),
            "llr_pident":   round(llr_pident, 4),
            "combined_llr": round(combined_llr, 4),
        })

    posterior_log_odds = prior_log_odds + total_llr
    posterior_prob     = sigmoid(posterior_log_odds)

    return {
        "prior_prob":         prior_prob,
        "prior_log_odds":     round(prior_log_odds, 4),
        "total_llr":          round(total_llr, 4),
        "posterior_log_odds": round(posterior_log_odds, 4),
        "posterior_prob":     round(posterior_prob, 4),
        "hit_details":        hit_details,
    }


# ─────────────────────────────────────────────────────────────────────────────
# SAMPLING HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def recommended_subsample_length(seq_length: int,
                                  insert_fraction: float = 0.30,
                                  min_length: int = 50,
                                  max_length: int = 200) -> int:
    """
    Recommend a subsample length that is small enough to sit mostly within
    the chimeric insert when it lands there, but long enough for reliable
    BLAST results.

    Target is half the expected insert size, clamped between min_length
    and max_length. For very short sequences the floor is reduced
    proportionally to avoid the subsample being larger than the insert.

    Parameters
    ----------
    seq_length      : int   full sequence length in bp
    insert_fraction : float expected insert size as fraction of seq_length
    min_length      : int   absolute minimum subsample length (default 50)
    max_length      : int   absolute maximum subsample length (default 200)
    """
    insert_size = seq_length * insert_fraction

    # Target: half the insert size so subsample can sit mostly within it
    recommended = int(insert_size * 0.5)

    # Floor scales with insert size — never let subsample exceed insert
    # For very short genes, accept shorter subsamples (worse BLAST but
    # better detection geometry)
    dynamic_floor = max(min_length, int(insert_size * 0.3))

    return int(np.clip(recommended, dynamic_floor, max_length))


def recommended_n_samples(seq_length: int,
                           subsample_length: int,
                           target_detection_prob: float = 0.90,
                           insert_fraction: float = 0.30,
                           min_overlap_fraction: float = 0.50,
                           max_samples: int = 15) -> int:
    """
    Calculate how many subsamples are needed to achieve a target probability
    of at least one subsample landing >50% within the chimeric insert.

    Capped at max_samples to prevent impractical numbers of BLAST calls
    for short sequences where the geometry is unfavourable.

    Parameters
    ----------
    seq_length            : int    full sequence length in bp
    subsample_length      : int    length of each subsample
    target_detection_prob : float  desired P(detect) e.g. 0.90
    insert_fraction       : float  expected insert size as fraction of L
    min_overlap_fraction  : float  minimum overlap fraction for detection
    max_samples           : int    hard cap on number of samples (default 15)
    """
    import math

    insert_size       = seq_length * insert_fraction
    min_overlap_bp    = subsample_length * min_overlap_fraction
    detectable_window = max(0.0, insert_size - min_overlap_bp)
    p_single_detect   = min(detectable_window / seq_length, 0.999)

    if p_single_detect <= 0:
        print(f"  WARNING: Sequence too short for reliable chimera detection "
              f"(insert={insert_size:.0f}bp, subsample={subsample_length}bp) "
              f"— using minimum {max_samples} samples")
        return max_samples

    n_samples = math.ceil(
        math.log(1 - target_detection_prob) /
        math.log(1 - p_single_detect)
    )

    if n_samples > max_samples:
        # Calculate what detection probability we actually achieve at the cap
        actual_prob = 1 - (1 - p_single_detect) ** max_samples
        print(f"  NOTE: {n_samples} samples needed for "
              f"{target_detection_prob*100:.0f}% detection but capped at "
              f"{max_samples} — actual P(detect) = {actual_prob*100:.1f}%")
        return max_samples

    return max(n_samples, 3)   # minimum 3

# ─────────────────────────────────────────────────────────────────────────────
# TOP-LEVEL SCORER
# ─────────────────────────────────────────────────────────────────────────────

def score_blast(query_accession, seq: str, taxon: str, gene: str,
                taxonomy, glossary, prior_prob, pident_threshold, k, beta,
                n: int = None,
                insert_fraction: float = 0.30,
                target_detection_prob: float = 0.90) -> dict:
    """
    Run verify_seq() with automatically scaled subsample length and number
    of samples based on sequence length and expected insert size, then
    score all hits and return an aggregated posterior.

    Parameters
    ----------
    query_accession      : str    accession of the query (for self-hit removal)
    seq                  : str    full query sequence
    taxon                : str    e.g. "Sorex_araneus"
    gene                 : str    e.g. "CYTB" or "CYTB_modified_1"
    taxonomy             : dict   from load_taxonomy()
    glossary             : dict   from load_glossary()
    prior_prob           : float  prior P(genuine)
    pident_threshold     : float  pident value contributing LLR of 0
    k                    : float  steepness of pident LLR mapping
    beta                 : float  softmax temperature over bitscores
    n                    : int    subsample length — if None, calculated
                                  automatically from sequence length
    insert_fraction      : float  expected insert size as fraction of
                                  sequence length (default 0.30)
    target_detection_prob: float  desired P(detect chimera) (default 0.90)

    Returns
    -------
    dict with posterior probability, chimera indicators, and full
    per-sample and per-hit breakdowns
    """
    # ── Strip _modified_N suffix so gene parser can resolve the gene name ─────
    clean_gene = re.sub(r'_modified(_\d+)?$', '', gene)

    # ── Auto-scale subsample length and n_samples ─────────────────────────────
    if n is None:
        n = recommended_subsample_length(len(seq),
                                         insert_fraction=insert_fraction)

    n_samples = recommended_n_samples(
        seq_length            = len(seq),
        subsample_length      = n,
        target_detection_prob = target_detection_prob,
        insert_fraction       = insert_fraction,
    )

    print(f"  Sequence: {len(seq)}bp | "
          f"Subsample: {n}bp | "
          f"Samples: {n_samples}")

    # ── Run BLAST across all subsamples ───────────────────────────────────────
    blast_result = verify_seq(seq, n, query_accession,
                              n_samples=n_samples)

    # ── Score each sample independently ──────────────────────────────────────
    sample_results    = []
    sample_posteriors = []

    for sample in blast_result["samples"]:
        if not sample["hits"]:
            print(f"  Sample {sample['sample_index'] + 1}: "
                  f"no hits — skipping")
            continue

        scored_hits = [
            score_hit(hit, taxon, clean_gene, taxonomy, glossary)
            for hit in sample["hits"]
        ]

        sample_agg = aggregate_hit_scores(
            scored_hits, prior_prob, pident_threshold, k, beta
        )
        sample_agg["sample_index"] = sample["sample_index"]
        sample_agg["query_start"]  = sample["query_start"]

        sample_posteriors.append(sample_agg["posterior_prob"])
        sample_results.append(sample_agg)

    if not sample_results:
        raise ValueError("No samples returned any BLAST hits.")

    # ── Aggregate log-odds across samples ─────────────────────────────────────
    # Prior applied once. Each sample contributes an independent LLR update.
    prior_log_odds      = math.log(prior_prob / (1 - prior_prob))
    cumulative_log_odds = prior_log_odds

    for s in sample_results:
        cumulative_log_odds += s["total_llr"]

    final_posterior = sigmoid(cumulative_log_odds)

    # ── Chimera indicators ────────────────────────────────────────────────────
    all_taxon_llrs = [
        h["llr_taxon"]
        for s in sample_results
        for h in s["hit_details"]
        if not h.get("gated_out", False)
    ]

    posterior_variance = float(np.var(sample_posteriors)) \
                         if len(sample_posteriors) > 1 else 0.0
    taxon_llr_variance = float(np.var(all_taxon_llrs)) \
                         if len(all_taxon_llrs) > 1 else 0.0
    prop_neg_taxon     = (sum(1 for x in all_taxon_llrs if x < 0) /
                          len(all_taxon_llrs)) if all_taxon_llrs else 0.0

    chimera_flag = posterior_variance > 0.05 or taxon_llr_variance > 2.0

    return {
        "query_accession":     query_accession,
        "taxon":               taxon,
        "gene":                gene,
        "seq_length":          len(seq),
        "subsample_length":    n,
        "n_samples":           n_samples,
        "sample_results":      sample_results,
        "sample_posteriors":   sample_posteriors,
        "posterior_variance":  round(posterior_variance, 6),
        "taxon_llr_variance":  round(taxon_llr_variance, 4),
        "prop_neg_taxon_hits": round(prop_neg_taxon, 4),
        "chimera_flag":        chimera_flag,
        "cumulative_log_odds": round(cumulative_log_odds, 4),
        "final_posterior":     round(final_posterior, 4),
        # Kept for backward compatibility with calibrate_priors()
        "prior_prob":          prior_prob,
        "prior_log_odds":      round(prior_log_odds, 4),
        "total_llr":           round(cumulative_log_odds - prior_log_odds, 4),
        "posterior_log_odds":  round(cumulative_log_odds, 4),
        "posterior_prob":      round(final_posterior, 4),
        "hit_details":         [h for s in sample_results
                                for h in s["hit_details"]],
    }

#print(score_blast("GU981106", "CACTTCCTTTGGATATGYTTGATGTGTTTTTGAATCATAATATCAATTCCTTTCTGAGGCAAGTTGAGAAGGTCAGAGATGAGGCATTGGTTCTTGTTATTCAATCCTATAATGAAGCAAAAATGAAATTTGATGAGCATAAGGTTGAAAAATCTATCACCCAACAACGAAAGACCTTTCAAATTCCAGGGTACACCATTCCTGTTGTTAATGTCGAAGTGTCTCCATTCACAGTAGAGATGTTTCCATTTGGTTATGTGATCCCAAAGGAGGTCAGCACCCCAAAGTTCACCATCCTGGGTTCTGGTTTCTCTGTGCCTTCCTATACTTTAGTCCTGCCCTTTCTAGAACTACCAGCTCTTCATATCCCTAAGTTTCTTGAGCTTTCTTTTCCAGACTTCAAAGTATCGAGTATCCCAAGGAATATTTTCATTCCAGCCCTGGGAAATGTTACATATGATTTTTCCTTTAAGTCAAGTGTCATTACACTGAATGCCAATGCTGGACTTTAT", 50,"Anourosorex_squamipes", "APOB", load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"), load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"), 0.95, 97.0, 0.3, 0.1))