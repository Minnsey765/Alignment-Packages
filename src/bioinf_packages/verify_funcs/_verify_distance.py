# _distance_verify.py

import os
import csv
import random
import numpy as np

from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner

try:
    from ..verify_funcs._score_hit import (recommended_subsample_length,
                                            recommended_n_samples)
    from ..verify_funcs._species_parser import load_taxonomy
except ImportError:
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from verify_funcs._score_hit import (recommended_subsample_length,
                                          recommended_n_samples)
    from verify_funcs._species_parser import load_taxonomy


# ─────────────────────────────────────────────────────────────────────────────
# TAXONOMY HELPERS
# copied here so this module is self-contained and does not depend
# on _phylo_verify.py being importable
# ─────────────────────────────────────────────────────────────────────────────

def get_expected_order(taxon: str, taxonomy: dict) -> str:
    """
    Return the taxonomic order for a Genus_species taxon string.

    Parameters
    ----------
    taxon    : str   e.g. "Solenodon_paradoxus"
    taxonomy : dict  loaded dict from load_taxonomy(), not a file path

    Returns
    -------
    str or None
    """
    if isinstance(taxonomy, str):
        raise TypeError(
            f"taxonomy must be a dict loaded by load_taxonomy(), "
            f"not a file path string. "
            f"Call taxonomy = load_taxonomy('{taxonomy}') first."
        )

    genus    = taxon.split("_")[0]
    tax_info = taxonomy.get(genus)
    if tax_info:
        return tax_info.get("order")
    return None


def get_ingroup_genera(order: str, taxonomy: dict) -> set:
    """
    Return all genera belonging to the given order in the taxonomy.
    """
    return {
        genus for genus, info in taxonomy.items()
        if info.get("order") == order
    }


# ─────────────────────────────────────────────────────────────────────────────
# PAIRWISE ALIGNER SETUP
# ─────────────────────────────────────────────────────────────────────────────

def _build_semiglobal_aligner() -> PairwiseAligner:
    """
    Build a semi-global PairwiseAligner suitable for comparing a short
    query subsample against full-length scaffold sequences.

    Semi-global means:
      - The entire query must align (no free end gaps on query side)
      - The target (scaffold) may extend beyond the query at either
        end without penalty (free end gaps on target side)

    This is appropriate for comparing a ~100bp subsample against a
    ~1000bp scaffold sequence. The subsample can "float" within the
    scaffold and align to whichever region it best matches, but every
    position of the subsample must participate in the alignment.

    This avoids the problem with local alignment (BLAST-style) where
    a short high-scoring patch can make two sequences appear more
    similar than they really are.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"

    aligner.match_score    =  1.0
    aligner.mismatch_score = -1.0

    # Updated attribute names for BioPython >= 1.82
    aligner.open_deletion_score    = -2.0
    aligner.extend_deletion_score  = -0.5

    aligner.open_left_insertion_score    = 0.0
    aligner.extend_left_insertion_score  = 0.0
    aligner.open_right_insertion_score   = 0.0
    aligner.extend_right_insertion_score = 0.0

    aligner.open_internal_insertion_score   = -2.0
    aligner.extend_internal_insertion_score = -0.5

    return aligner


# ─────────────────────────────────────────────────────────────────────────────
# PAIRWISE IDENTITY
# ─────────────────────────────────────────────────────────────────────────────

def semiglobal_identity(query: str,
                         target: str,
                         aligner: PairwiseAligner) -> float:
    """
    Compute semi-global pairwise identity between a query subsample
    and a target scaffold sequence.

    Identity is the fraction of query positions that are identical
    to their aligned target position, normalised by query length.
    This means a 100bp query that perfectly matches 100bp of a
    1000bp target scores 1.0, not 0.1.

    Parameters
    ----------
    query   : str              query subsample (no gaps)
    target  : str              scaffold sequence (gaps stripped)
    aligner : PairwiseAligner  from _build_semiglobal_aligner()

    Returns
    -------
    float   identity in range [0, 1]
    """
    if not query or not target:
        return 0.0

    try:
        alignments = aligner.align(target, query)
        if not alignments:
            return 0.0

        best = alignments[0]

        # Count identical positions
        # aligned_sequences gives the two rows of the alignment
        aligned = best.aligned

        # Use the alignment coordinates to count matches
        n_identical = 0
        for (t_start, t_end), (q_start, q_end) in zip(
            aligned[0], aligned[1]
        ):
            t_block = target[t_start:t_end]
            q_block = query[q_start:q_end]
            n_identical += sum(t == q for t, q in
                               zip(t_block, q_block))

        return n_identical / len(query)

    except Exception as e:
        return 0.0


def verify_by_distance(query_seq: str,
                        query_taxon: str,
                        query_accession: str,
                        scaffold_fasta: str,
                        taxonomy: dict,
                        insert_fraction: float = 0.30,
                        target_detection_prob: float = 0.90,
                        outlier_ratio_threshold: float = 0.95,
                        min_ratio_variance_flag: float = 0.001,
                        k_nearest: int = 3) -> dict:
    """
    Verify a single sequence by comparing random subsamples against
    scaffold sequences using semi-global pairwise alignment.

    For each subsample two similarity metrics are computed:

    Mean similarity ratio
        Ingroup mean identity / outgroup mean identity across all
        scaffold sequences. Captures overall ingroup affinity.

    K-nearest similarity ratio
        Mean identity to the k most similar ingroup sequences /
        mean identity to the k most similar outgroup sequences.
        More robust when the scaffold is phylogenetically
        heterogeneous — a genuine sequence from an underrepresented
        family may have low mean ingroup identity but should still
        have close relatives in the ingroup.

    A subsample is flagged as an outlier only when BOTH metrics fall
    below outlier_ratio_threshold. Requiring both to agree reduces
    false positives from scaffold heterogeneity while maintaining
    sensitivity to genuine contamination.

    High variance in similarity ratio across subsamples is a chimera
    indicator even when the mean ratio looks acceptable — some
    subsamples land in the contaminated region and some do not.

    Parameters
    ----------
    query_seq               : str   full query gene sequence (no gaps)
    query_taxon             : str   e.g. "Solenodon_paradoxus"
    query_accession         : str   e.g. "AY530070"
    scaffold_fasta          : str   path to scaffold FASTA file.
                                    Aligned or unaligned — gaps are
                                    stripped before alignment so
                                    either works.
    taxonomy                : dict  loaded dict from load_taxonomy().
                                    Must be a dict, not a file path.
    insert_fraction         : float expected chimeric insert as
                                    fraction of total sequence length
                                    (default 0.3)
    target_detection_prob   : float desired P(detect chimera junction)
                                    (default 0.9)
    outlier_ratio_threshold : float a subsample is only flagged as an
                                    outlier when BOTH mean_ratio and
                                    knn_ratio fall below this value.
                                    Default 0.95 — the subsample must
                                    be meaningfully more similar to
                                    outgroup than ingroup on both
                                    metrics to be considered suspicious.
                                    Prevents ratio=0.9994 noise from
                                    triggering false positives.
    min_ratio_variance_flag : float chimera_flag is set if variance
                                    exceeds this even without genuine
                                    outliers. Default 0.001.
    k_nearest               : int   number of nearest neighbours for
                                    knn_ratio metric (default 3)

    Returns
    -------
    dict with:
        query_accession          : str
        query_taxon              : str
        expected_order           : str
        n_ingroup_scaffold       : int
        n_outgroup_scaffold      : int
        n_samples                : int
        subsample_length         : int
        n_outliers               : int   samples where both metrics
                                         are below threshold
        prop_outliers            : float
        mean_similarity_ratio    : float mean of per-sample mean ratios
        mean_knn_similarity_ratio: float mean of per-sample knn ratios
        min_mean_ratio           : float
        min_knn_ratio            : float
        mean_ratio_variance      : float key chimera indicator
        knn_ratio_variance       : float
        mean_ingroup_identity    : float
        mean_outgroup_identity   : float
        chimera_flag             : bool
        sample_scores            : list of per-sample dicts
    """
    # Guard against passing file path instead of loaded dict
    if isinstance(taxonomy, str):
        raise TypeError(
            f"taxonomy must be a dict from load_taxonomy(), "
            f"not a file path. "
            f"Call: taxonomy = load_taxonomy('{taxonomy}')"
        )

    seq_length = len(query_seq)

    # ── Calculate recommended subsample parameters ────────────────────────────
    n = recommended_subsample_length(
        seq_length, insert_fraction=insert_fraction
    )
    n_samples = recommended_n_samples(
        seq_length            = seq_length,
        subsample_length      = n,
        target_detection_prob = target_detection_prob,
        insert_fraction       = insert_fraction,
    )

    print(f"\nVerifying {query_taxon} ({query_accession}) "
          f"by pairwise distance")
    print(f"Sequence: {seq_length}bp | "
          f"Subsample: {n}bp | Samples: {n_samples}")

    # ── Load and classify scaffold sequences ──────────────────────────────────
    scaffold_records = list(SeqIO.parse(scaffold_fasta, "fasta"))

    if not scaffold_records:
        raise FileNotFoundError(
            f"No sequences found in scaffold FASTA: {scaffold_fasta}"
        )

    expected_order = get_expected_order(query_taxon, taxonomy)
    ingroup_genera = (get_ingroup_genera(expected_order, taxonomy)
                      if expected_order else set())

    def get_genus(record):
        return record.id.split("|")[0].split("_")[0]

    # Strip alignment gaps — works whether FASTA is aligned or not
    ingroup_seqs  = [str(r.seq).replace("-", "")
                     for r in scaffold_records
                     if get_genus(r) in ingroup_genera]
    outgroup_seqs = [str(r.seq).replace("-", "")
                     for r in scaffold_records
                     if get_genus(r) not in ingroup_genera]

    print(f"Expected order   : {expected_order}")
    print(f"Ingroup genera   : {', '.join(sorted(ingroup_genera))}")
    print(f"Ingroup scaffold : {len(ingroup_seqs)} sequences")
    print(f"Outgroup scaffold: {len(outgroup_seqs)} sequences")

    if not ingroup_seqs:
        raise ValueError(
            f"No ingroup sequences found for order '{expected_order}' "
            f"in {scaffold_fasta}. Check that your scaffold contains "
            f"taxa from your taxonomy CSV."
        )

    # ── Build aligner once, reuse for all comparisons ─────────────────────────
    aligner = _build_semiglobal_aligner()

    # ── Define helper functions ───────────────────────────────────────────────

    def mean_identity_to_group(subsample: str,
                                seqs: list) -> tuple:
        """
        Mean and std of semi-global identity between subsample and
        all sequences in a group.
        """
        if not seqs:
            return 0.0, 0.0
        identities = [semiglobal_identity(subsample, t, aligner)
                      for t in seqs]
        return float(np.mean(identities)), float(np.std(identities))

    def top_k_identity(subsample: str,
                       seqs: list,
                       k: int) -> tuple:
        """
        Mean and std of identity to the k most similar sequences in
        a group. More robust than overall mean when the group contains
        phylogenetically heterogeneous sequences, because it focuses
        on the closest relatives rather than being dragged down by
        distant members of the same group.
        """
        if not seqs:
            return 0.0, 0.0
        identities = sorted(
            [semiglobal_identity(subsample, t, aligner)
             for t in seqs],
            reverse=True
        )
        top_k = identities[:min(k, len(identities))]
        return float(np.mean(top_k)), float(np.std(top_k))

    # ── Score each subsample ──────────────────────────────────────────────────
    sample_scores = []

    for i in range(1, n_samples + 1):
        start     = random.randint(0, seq_length - n)
        subsample = query_seq[start:start + n]

        # Mean identity across all scaffold sequences in each group
        ig_mean, ig_std = mean_identity_to_group(subsample,
                                                  ingroup_seqs)
        og_mean, og_std = mean_identity_to_group(subsample,
                                                  outgroup_seqs)

        # K-nearest identity — focuses on closest relatives
        ig_knn_mean, ig_knn_std = top_k_identity(subsample,
                                                   ingroup_seqs,
                                                   k_nearest)
        og_knn_mean, og_knn_std = top_k_identity(subsample,
                                                   outgroup_seqs,
                                                   k_nearest)

        # Similarity ratios — values > 1 mean more similar to ingroup
        mean_ratio = (ig_mean / og_mean
                      if og_mean > 0 else float("inf"))
        knn_ratio  = (ig_knn_mean / og_knn_mean
                      if og_knn_mean > 0 else float("inf"))

        # Outlier only when BOTH metrics agree the subsample is more
        # similar to outgroup than ingroup by a meaningful margin.
        # Requiring both to agree reduces false positives from scaffold
        # heterogeneity while maintaining sensitivity to contamination.
        is_outlier = (
            mean_ratio < outlier_ratio_threshold
            and knn_ratio < outlier_ratio_threshold
        )

        score = {
            "sample":                i,
            "start":                 start,
            "end":                   start + n,
            "ingroup_identity":      round(ig_mean, 4),
            "ingroup_identity_std":  round(ig_std, 4),
            "outgroup_identity":     round(og_mean, 4),
            "outgroup_identity_std": round(og_std, 4),
            "mean_similarity_ratio": round(mean_ratio, 4),
            "ingroup_knn_identity":  round(ig_knn_mean, 4),
            "ingroup_knn_std":       round(ig_knn_std, 4),
            "outgroup_knn_identity": round(og_knn_mean, 4),
            "outgroup_knn_std":      round(og_knn_std, 4),
            "knn_similarity_ratio":  round(knn_ratio, 4),
            "is_outlier":            is_outlier,
        }
        sample_scores.append(score)

        flag = "⚠ OUTLIER" if is_outlier else "✓"
        print(f"  {flag} Sample {i}/{n_samples} | "
              f"pos {start}-{start+n} | "
              f"mean: ingroup={ig_mean:.3f}±{ig_std:.3f} "
              f"outgroup={og_mean:.3f}±{og_std:.3f} "
              f"ratio={mean_ratio:.3f} | "
              f"knn: ingroup={ig_knn_mean:.3f} "
              f"outgroup={og_knn_mean:.3f} "
              f"ratio={knn_ratio:.3f}")

    # ── Sequence-level summary ────────────────────────────────────────────────
    mean_ratios = [s["mean_similarity_ratio"] for s in sample_scores
                   if s["mean_similarity_ratio"] != float("inf")]
    knn_ratios  = [s["knn_similarity_ratio"] for s in sample_scores
                   if s["knn_similarity_ratio"] != float("inf")]
    outliers    = [s for s in sample_scores if s["is_outlier"]]

    mean_ratio_variance = (round(float(np.var(mean_ratios)), 4)
                           if mean_ratios else None)
    knn_ratio_variance  = (round(float(np.var(knn_ratios)), 4)
                           if knn_ratios else None)

    # chimera_flag is set when:
    #   - at least one subsample is a genuine outlier on both metrics, OR
    #   - ratio variance is high even without clear outliers
    #     (mixed signal — some subsamples ingroup, some outgroup)
    chimera_flag = (
        len(outliers) > 0
        or (mean_ratio_variance is not None
            and mean_ratio_variance > min_ratio_variance_flag)
        or (knn_ratio_variance is not None
            and knn_ratio_variance > min_ratio_variance_flag)
    )

    summary = {
        "query_accession":           query_accession,
        "query_taxon":               query_taxon,
        "expected_order":            expected_order,
        "n_ingroup_scaffold":        len(ingroup_seqs),
        "n_outgroup_scaffold":       len(outgroup_seqs),
        "n_samples":                 n_samples,
        "subsample_length":          n,
        "k_nearest":                 k_nearest,
        "outlier_ratio_threshold":   outlier_ratio_threshold,
        "n_outliers":                len(outliers),
        "prop_outliers":             (round(len(outliers) / n_samples,
                                           4) if n_samples else None),
        "mean_similarity_ratio":     (round(float(np.mean(mean_ratios)),
                                           4) if mean_ratios else None),
        "mean_knn_similarity_ratio": (round(float(np.mean(knn_ratios)),
                                           4) if knn_ratios else None),
        "min_mean_ratio":            (round(float(np.min(mean_ratios)),
                                           4) if mean_ratios else None),
        "max_mean_ratio":            (round(float(np.max(mean_ratios)),
                                           4) if mean_ratios else None),
        "min_knn_ratio":             (round(float(np.min(knn_ratios)),
                                           4) if knn_ratios else None),
        "max_knn_ratio":             (round(float(np.max(knn_ratios)),
                                           4) if knn_ratios else None),
        "mean_ratio_variance":       mean_ratio_variance,
        "knn_ratio_variance":        knn_ratio_variance,
        "mean_ingroup_identity":     round(float(np.mean(
                                       [s["ingroup_identity"]
                                        for s in sample_scores])), 4),
        "mean_outgroup_identity":    round(float(np.mean(
                                       [s["outgroup_identity"]
                                        for s in sample_scores])), 4),
        "mean_ingroup_knn_identity": round(float(np.mean(
                                       [s["ingroup_knn_identity"]
                                        for s in sample_scores])), 4),
        "mean_outgroup_knn_identity":round(float(np.mean(
                                       [s["outgroup_knn_identity"]
                                        for s in sample_scores])), 4),
        "chimera_flag":              chimera_flag,
        "sample_scores":             sample_scores,
    }

    print(f"\n── Distance summary ─────────────────────────────────────────")
    print(f"  Expected order               : {expected_order}")
    print(f"  Mean ingroup identity        : "
          f"{summary['mean_ingroup_identity']}")
    print(f"  Mean outgroup identity       : "
          f"{summary['mean_outgroup_identity']}")
    print(f"  Mean similarity ratio        : "
          f"{summary['mean_similarity_ratio']}")
    print(f"  Mean knn similarity ratio    : "
          f"{summary['mean_knn_similarity_ratio']}")
    print(f"  Min mean ratio               : "
          f"{summary['min_mean_ratio']}")
    print(f"  Min knn ratio                : "
          f"{summary['min_knn_ratio']}")
    print(f"  Mean ratio variance          : {mean_ratio_variance}")
    print(f"  Knn ratio variance           : {knn_ratio_variance}")
    print(f"  Outlier samples (both < "
          f"{outlier_ratio_threshold})   : "
          f"{len(outliers)}/{n_samples}")
    print(f"  Chimera flag                 : {chimera_flag}")

    return summary


# ─────────────────────────────────────────────────────────────────────────────
# BATCH VERIFICATION
# ─────────────────────────────────────────────────────────────────────────────

def verify_batch_by_distance(query_fasta: str,
                              scaffold_fasta: str,
                              taxonomy: dict,
                              output_dir: str,
                              insert_fraction: float = 0.30,
                              target_detection_prob: float = 0.90
                              ) -> list:
    """
    Run verify_by_distance() on every sequence in a FASTA file.

    FASTA header format expected:
        >Species_name|gene|orien:+/-|accession:XXXX

    Parameters
    ----------
    query_fasta           : str   path to FASTA of sequences to verify
    scaffold_fasta        : str   path to scaffold FASTA (aligned or
                                  unaligned)
    taxonomy              : dict  loaded dict from load_taxonomy()
    output_dir            : str   directory to save summary CSV and
                                  per-sequence sample score CSVs
    insert_fraction       : float expected chimeric insert fraction
    target_detection_prob : float desired P(detect chimera)

    Returns
    -------
    list of summary dicts, one per sequence. Writes:
        {output_dir}/distance_verification_summary.csv
        {output_dir}/{taxon}_{accession}_sample_scores.csv
    """
    if isinstance(taxonomy, str):
        raise TypeError(
            f"taxonomy must be a dict from load_taxonomy(), "
            f"not a file path. "
            f"Call: taxonomy = load_taxonomy('{taxonomy}')"
        )

    os.makedirs(output_dir, exist_ok=True)

    records = list(SeqIO.parse(query_fasta, "fasta"))
    total   = len(records)
    results = []

    print(f"\nBatch distance verification: {total} sequences")
    print(f"Scaffold : {scaffold_fasta}")
    print(f"Output   : {output_dir}\n")

    for i, record in enumerate(records, 1):

        # ── Parse FASTA header ────────────────────────────────────────────────
        parts       = record.description.split("|")
        taxon       = parts[0].strip()
        gene        = parts[1].strip() if len(parts) > 1 else "unknown"
        orientation = "+"
        accession   = None

        for p in parts:
            if p.startswith("orien:"):
                orientation = p.split(":", 1)[1].strip()
            if p.startswith("accession:"):
                accession = p.split(":", 1)[1].strip()

        if not accession:
            print(f"  [{i}/{total}] Skipping {taxon} — "
                  f"no accession in header")
            continue

        print(f"\n[{i}/{total}] {taxon} | {gene} | {accession}")

        try:
            result = verify_by_distance(
                query_seq             = str(record.seq),
                query_taxon           = taxon,
                query_accession       = accession,
                scaffold_fasta        = scaffold_fasta,
                taxonomy              = taxonomy,
                insert_fraction       = insert_fraction,
                target_detection_prob = target_detection_prob,
            )
            result["gene"]    = gene
            result["_header"] = record.description
            results.append(result)

            # ── Save per-sequence sample scores to CSV ────────────────────────
            scores_path = os.path.join(
                output_dir,
                f"{taxon}_{accession}_sample_scores.csv"
            )
            if result["sample_scores"]:
                score_keys = list(result["sample_scores"][0].keys())
                with open(scores_path, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=score_keys)
                    writer.writeheader()
                    writer.writerows(result["sample_scores"])
                print(f"  Sample scores: {scores_path}")

            print(f"  ✓ chimera_flag={result['chimera_flag']} | "
                  f"mean_ratio={result['mean_similarity_ratio']} | "
                  f"variance={result['ratio_variance']} | "
                  f"outliers={result['n_outliers']}/"
                  f"{result['n_samples']}")

        except Exception as e:
            print(f"  ✗ Failed: {e}")
            import traceback
            traceback.print_exc()

    # ── Write summary CSV ─────────────────────────────────────────────────────
    if results:
        summary_path = os.path.join(
            output_dir, "distance_verification_summary.csv"
        )
        summary_keys = [k for k in results[0].keys()
                        if k != "sample_scores"]

        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summary_keys)
            writer.writeheader()
            for r in results:
                writer.writerow({k: r[k] for k in summary_keys})

        print(f"\nSummary CSV: {summary_path}")

        flagged = [r for r in results if r["chimera_flag"]]
        if flagged:
            print(f"\n⚠  {len(flagged)} sequence(s) flagged:")
            for r in flagged:
                print(f"  {r['query_taxon']} "
                      f"({r['query_accession']}) | "
                      f"gene={r.get('gene','?')} | "
                      f"outliers={r['n_outliers']}/"
                      f"{r['n_samples']} | "
                      f"ratio={r['mean_similarity_ratio']} | "
                      f"var={r['ratio_variance']}")
        else:
            print(f"\n✓ No sequences flagged")

    print(f"\nBatch complete: {len(results)}/{total} succeeded")
    return results


def verify_by_ingroup_distance(
        query_seq: str,
        query_taxon: str,
        query_accession: str,
        scaffold_fasta: str,
        taxonomy: dict,
        insert_fraction: float = 0.30,
        target_detection_prob: float = 0.90,
        k_nearest: int = 3,
        min_specific_ingroup_identity: float = None,
        specific_ingroup_sd_threshold: float = 3.0,
        wrong_group_margin_threshold: float = 0.02,
        known_orphan_taxa: list = None,
        orphan_sd_threshold: float = 3.0,
        threshold_json: str = None) -> dict:
    """
    Verify a sequence by testing how well each subsample scores
    against the ingroup scaffold only, broken down by taxonomic
    grouping.

    A genuine sequence should score highest against its own
    taxonomic group (subfamily > family > order, whichever is
    available in the scaffold) and progressively lower against
    more distant ingroup groups. A chimeric sequence will score
    inconsistently — some subsamples will score well against the
    expected group and others will score poorly because they fall
    in the contaminated region and match a different group better.

    Taxonomic level used for scoring is determined hierarchically:
        1. Subfamily — if scaffold contains sequences from the
                       query's subfamily
        2. Family    — if no subfamily match in scaffold
        3. Order     — if no family match in scaffold (with warning)

    This means Podogymnura (Galericinae, Erinaceidae) will score
    against Erinaceidae if Galericinae is absent from the scaffold,
    rather than falling back to all Eulipotyphla.

    Because the subsampling method gives a 90% probability that at
    least one subsample overlaps a chimeric junction, the chimera
    flag is set if ANY single subsample is an outlier.

    Parameters
    ----------
    query_seq                    : str   full query sequence
    query_taxon                  : str   e.g. "Podogymnura_truei"
    query_accession              : str   e.g. "JN414025"
    scaffold_fasta               : str   path to scaffold FASTA
                                         (aligned or unaligned —
                                         gaps are stripped)
    taxonomy                     : dict  from load_taxonomy(). Must
                                         contain family and ideally
                                         subfamily columns.
    insert_fraction              : float expected chimeric insert
                                         fraction (default 0.3)
    target_detection_prob        : float desired P(detect chimera)
                                         (default 0.9)
    k_nearest                    : int   nearest neighbours within
                                         each group for scoring
                                         (default 3)
    min_specific_ingroup_identity: float fixed minimum identity
                                         threshold for the specific
                                         ingroup. If None, threshold
                                         is auto-calibrated as
                                         mean - specific_ingroup_sd_threshold
                                         * std across all subsamples.
    specific_ingroup_sd_threshold: float SDs below mean to set the
                                         auto-calibrated threshold.
                                         Default 2.0. Lower = more
                                         sensitive, higher = more
                                         conservative.
    known_orphan_taxa            : list  genus names with no close
                                         relatives in scaffold.
                                         These use a more lenient
                                         threshold and do not require
                                         correct family matching.
                                         e.g. ["Solenodon", "Atopogale"]
    orphan_sd_threshold          : float SD threshold for orphan taxa.
                                         Default 3.0 (more lenient).

    Returns
    -------
    dict with:
        query_accession           : str
        query_taxon               : str
        query_specific_group      : str   the taxonomic group used as
                                          specific ingroup
        scoring_level             : str   "subfamily", "family", or
                                          "order"
        fallback_used             : bool  True if preferred level was
                                          not in scaffold
        is_orphan                 : bool
        expected_order            : str
        n_samples                 : int
        subsample_length          : int
        outlier_threshold         : float
        n_outliers                : int
        prop_outliers             : float
        n_correct_group           : int   subsamples whose nearest
                                          group match is correct
        prop_correct_group        : float
        specific_ingroup_mean     : float mean identity to specific
                                          ingroup across all subsamples
        specific_ingroup_std      : float
        specific_ingroup_variance : float key chimera indicator
        per_group_means           : dict  group -> mean knn identity
        chimera_flag              : bool  True if any outlier
        chimera_evidence          : str   human-readable explanation
        scoring_groups            : dict  group -> list of scaffold IDs
        sample_scores             : list  per-subsample dicts
    """
    if isinstance(taxonomy, str):
        raise TypeError(
            f"taxonomy must be a dict from load_taxonomy(), "
            f"not a file path. "
            f"Call: taxonomy = load_taxonomy('path')"
        )

    seq_length = len(query_seq)

    # ── Recommended subsample parameters ─────────────────────────────────────
    n = recommended_subsample_length(
        seq_length, insert_fraction=insert_fraction
    )
    n_samples = recommended_n_samples(
        seq_length            = seq_length,
        subsample_length      = n,
        target_detection_prob = target_detection_prob,
        insert_fraction       = insert_fraction,
    )

    print(f"\nVerifying {query_taxon} ({query_accession}) "
          f"by ingroup distance")
    print(f"Sequence: {seq_length}bp | "
          f"Subsample: {n}bp | Samples: {n_samples}")

    # ── Determine query taxonomic levels ──────────────────────────────────────
    query_genus = query_taxon.split("_")[0]
    query_tax   = taxonomy.get(query_genus, {})

    expected_order  = query_tax.get("order")
    query_order = expected_order
    query_subfamily = (query_tax.get("sub-family")
                       or query_tax.get("Sub-family")
                       or query_tax.get("subfamily"))
    query_family    = (query_tax.get("family")
                       or query_tax.get("Family"))

    is_orphan = (known_orphan_taxa is not None
                 and query_genus in known_orphan_taxa)

    effective_sd_threshold = (orphan_sd_threshold
                              if is_orphan
                              else specific_ingroup_sd_threshold)

    print(f"Query genus      : {query_genus}")
    print(f"Query subfamily  : {query_subfamily or 'not in taxonomy'}")
    print(f"Query family     : {query_family or 'not in taxonomy'}")
    print(f"Expected order   : {expected_order}")
    print(f"Is orphan taxon  : {is_orphan}")

    if not query_family and not query_subfamily and not expected_order:
        raise ValueError(
            f"Cannot determine any taxonomic level for genus "
            f"'{query_genus}'. Check taxonomy CSV."
        )

    # ── Load scaffold sequences ───────────────────────────────────────────────
    scaffold_records = list(SeqIO.parse(scaffold_fasta, "fasta"))
    if not scaffold_records:
        raise FileNotFoundError(
            f"No sequences found in scaffold: {scaffold_fasta}"
        )

    # ── Helper: extract genus from record ID ──────────────────────────────────
    def get_genus_from_id(record_id: str) -> str:
        return record_id.split("|")[0].split("_")[0]

    # ── Helper: get taxonomic levels for a genus ──────────────────────────────
    def get_tax_levels(genus: str) -> dict:
        tax = taxonomy.get(genus, {})
        return {
            "subfamily": (tax.get("sub-family")
                          or tax.get("Sub-family")
                          or tax.get("subfamily")),
            "family":    (tax.get("family")
                          or tax.get("Family")),
            "order":     tax.get("order"),
        }

    # ── Group scaffold sequences by all taxonomic levels ──────────────────────
    from collections import defaultdict

    subfamily_groups = defaultdict(list)
    family_groups    = defaultdict(list)
    order_groups     = defaultdict(list)

    for r in scaffold_records:
        genus  = get_genus_from_id(r.id)
        levels = get_tax_levels(genus)
        pair   = (r.id, str(r.seq).replace("-", ""))

        if levels["subfamily"]:
            subfamily_groups[levels["subfamily"]].append(pair)
        if levels["family"]:
            family_groups[levels["family"]].append(pair)
        if levels["order"]:
            order_groups[levels["order"]].append(pair)

    # ── Determine specific ingroup using hierarchical fallback ────────────────
    specific_ingroup_seqs  = []
    specific_ingroup_label = None
    scoring_level          = None
    fallback_used          = False
    scoring_groups         = {}

    if query_subfamily and query_subfamily in subfamily_groups:
        # Best case — subfamily present in scaffold
        specific_ingroup_seqs  = subfamily_groups[query_subfamily]
        specific_ingroup_label = query_subfamily
        scoring_level          = "subfamily"
        scoring_groups         = dict(subfamily_groups)

    elif query_subfamily and query_family and \
            query_family in family_groups:
        # Subfamily not in scaffold, fall back to family
        specific_ingroup_seqs  = family_groups[query_family]
        specific_ingroup_label = query_family
        scoring_level          = "family"
        scoring_groups         = dict(family_groups)
        fallback_used          = True
        print(f"\n  Note: subfamily '{query_subfamily}' not found in "
              f"scaffold — falling back to family '{query_family}' "
              f"({len(specific_ingroup_seqs)} sequences). "
              f"This is expected e.g. for Galericinae → Erinaceidae.")

    elif query_family and query_family in family_groups:
        # No subfamily in taxonomy, use family directly
        specific_ingroup_seqs  = family_groups[query_family]
        specific_ingroup_label = query_family
        scoring_level          = "family"
        scoring_groups         = dict(family_groups)

    elif query_order and query_order in order_groups:
        # No family match — fall back to order (least desirable)
        specific_ingroup_seqs  = order_groups[query_order]
        specific_ingroup_label = query_order
        scoring_level          = "order"
        scoring_groups         = dict(order_groups)
        fallback_used          = True
        print(f"\n  Warning: no family or subfamily match found for "
              f"'{query_genus}' in scaffold — falling back to order "
              f"'{query_order}'. Chimera detection will be less "
              f"precise. Consider adding sequences from the same "
              f"family to the scaffold.")

    else:
        raise ValueError(
            f"No scaffold sequences found for '{query_genus}' at "
            f"any taxonomic level (subfamily={query_subfamily}, "
            f"family={query_family}, order={expected_order}). "
            f"Check scaffold contains Eulipotyphla sequences and "
            f"taxonomy CSV is complete."
        )

    print(f"\nSpecific ingroup : {specific_ingroup_label} "
          f"[{scoring_level}] "
          f"({len(specific_ingroup_seqs)} sequences)")
    print(f"Fallback used    : {fallback_used}")

    print(f"\nScaffold groups at {scoring_level} level:")
    for grp, members in sorted(scoring_groups.items()):
        marker = " ◄ QUERY" if grp == specific_ingroup_label else ""
        print(f"  {grp}: {len(members)} sequences{marker}")
        for rid, _ in members:
            print(f"    {rid}")

    # ── Build aligner ─────────────────────────────────────────────────────────
    aligner = _build_semiglobal_aligner()

    # ── Score helper ──────────────────────────────────────────────────────────
    def score_against_pairs(subsample: str, pairs: list) -> dict:
        """
        Score subsample against a list of (id, ungapped_seq) pairs.
        Returns mean, std, knn mean, nearest id and identity.
        """
        if not pairs:
            return {
                "mean":             0.0,
                "std":              0.0,
                "knn_mean":         0.0,
                "nearest_id":       None,
                "nearest_identity": 0.0,
            }

        identities = [
            (semiglobal_identity(subsample, seq, aligner), rid)
            for rid, seq in pairs
        ]
        identities.sort(reverse=True)

        vals     = [v for v, _ in identities]
        top_k    = vals[:min(k_nearest, len(vals))]
        best_val, best_id = identities[0]

        return {
            "mean":             float(np.mean(vals)),
            "std":              float(np.std(vals)),
            "knn_mean":         float(np.mean(top_k)),
            "nearest_id":       best_id,
            "nearest_identity": best_val,
        }

    # ── Score each subsample ──────────────────────────────────────────────────
    sample_scores = []

    for i in range(1, n_samples + 1):
        start     = random.randint(0, seq_length - n)
        subsample = query_seq[start:start + n]

        # Score against specific ingroup
        sp = score_against_pairs(subsample, specific_ingroup_seqs)

        # Score against every group at the chosen scoring level
        per_group = {}
        for grp, pairs in scoring_groups.items():
            g = score_against_pairs(subsample, pairs)
            per_group[grp] = {
                "mean":             round(g["mean"], 4),
                "knn_mean":         round(g["knn_mean"], 4),
                "nearest_id":       g["nearest_id"],
                "nearest_identity": round(g["nearest_identity"], 4),
            }

        # Which group does this subsample match best by knn?
        nearest_group  = max(per_group.keys(),
                             key=lambda g: per_group[g]["knn_mean"])
        is_correct     = nearest_group == specific_ingroup_label

        sample_scores.append({
            "sample":                    i,
            "start":                     start,
            "end":                       start + n,
            "specific_ingroup_identity": round(sp["mean"], 4),
            "specific_ingroup_std":      round(sp["std"], 4),
            "specific_ingroup_knn":      round(sp["knn_mean"], 4),
            "nearest_specific_id":       sp["nearest_id"],
            "nearest_specific_identity": round(sp["nearest_identity"],
                                               4),
            "per_group_identity":        per_group,
            "nearest_group_match":       nearest_group,
            "is_correct_group":          is_correct,
            "is_outlier":                False,  # set below
        })

    # ── Auto-calibrate or apply fixed outlier threshold ───────────────────────
    specific_ids = [s["specific_ingroup_identity"]
                    for s in sample_scores]
    sp_mean = float(np.mean(specific_ids))
    sp_std  = float(np.std(specific_ids))

    if min_specific_ingroup_identity is not None:
        threshold = min_specific_ingroup_identity
        print(f"\nUsing fixed threshold  : {threshold:.4f}")
    else:
        threshold = sp_mean - effective_sd_threshold * sp_std
        print(f"\nAuto-calibrated threshold: "
              f"{sp_mean:.4f} - {effective_sd_threshold} × "
              f"{sp_std:.4f} = {threshold:.4f}")

    # ── Apply outlier flag ─────────────────────────────────────────────
    for s in sample_scores:

        # Compute how much better the nearest wrong group scores
        # compared to the specific ingroup
        specific_knn = s["specific_ingroup_knn"]
        nearest_grp  = s["nearest_group_match"]

        if nearest_grp != specific_ingroup_label:
            nearest_knn = s["per_group_identity"][nearest_grp]["knn_mean"]
            wrong_group_margin = nearest_knn - specific_knn
        else:
            wrong_group_margin = 0.0

        s["wrong_group_margin"] = round(wrong_group_margin, 4)

        if is_orphan:
            s["is_outlier"] = (
                s["specific_ingroup_identity"] < threshold
            )
        else:
            # Wrong group match only flagged if the margin is
            # meaningful — prevents noise from triggering false
            # positives when two groups score nearly identically.
            # A genuine chimera will show a clear preference for
            # the wrong group, not a marginal 0.001 difference.
            wrong_group_flag = (
                not s["is_correct_group"]
                and wrong_group_margin > wrong_group_margin_threshold
            )
            s["is_outlier"] = (
                s["specific_ingroup_identity"] < threshold
                or wrong_group_flag
            )
    # ── Print per-sample results ──────────────────────────────────────────────
    for s in sample_scores:
        flag = "⚠ OUTLIER" if s["is_outlier"] else "✓"
        group_str = " | ".join(
            f"{g}={v['knn_mean']:.3f}"
            for g, v in sorted(s["per_group_identity"].items())
        )
        print(f"  {flag} Sample {s['sample']}/{n_samples} | "
              f"pos {s['start']}-{s['end']} | "
              f"specific={s['specific_ingroup_knn']:.3f} | "
              f"nearest_group={s['nearest_group_match']} | "
              f"correct={s['is_correct_group']}")
        print(f"    per-group knn: {group_str}")
        print(f"    nearest specific match: "
              f"{s['nearest_specific_id']} "
              f"({s['nearest_specific_identity']:.3f})")

    # ── Sequence-level summary ────────────────────────────────────────────────
    outliers     = [s for s in sample_scores if s["is_outlier"]]
    correct      = [s for s in sample_scores if s["is_correct_group"]]
    sp_variance  = round(float(np.var(specific_ids)), 6)

    # Per-group mean knn identity across all subsamples
    per_group_means = {}
    for grp in scoring_groups:
        vals = [s["per_group_identity"][grp]["knn_mean"]
                for s in sample_scores
                if grp in s["per_group_identity"]]
        per_group_means[grp] = (round(float(np.mean(vals)), 4)
                                 if vals else None)

    # chimera_flag: True if ANY subsample is an outlier
    chimera_flag = len(outliers) > 0

    # Build human-readable evidence
    if not chimera_flag:
        evidence = (
            f"All {n_samples} subsamples score consistently with "
            f"{specific_ingroup_label} [{scoring_level}] "
            f"(mean={sp_mean:.4f}, variance={sp_variance:.6f})"
        )
    else:
        outlier_positions = [(s["start"], s["end"])
                             for s in outliers]
        wrong_groups      = list({
            s["nearest_group_match"]
            for s in outliers
            if not s["is_correct_group"]
        })
        evidence = (
            f"{len(outliers)}/{n_samples} subsamples flagged. "
            f"Positions: {outlier_positions}. "
        )
        if wrong_groups:
            evidence += (
                f"These subsamples scored highest against: "
                f"{', '.join(wrong_groups)} instead of "
                f"{specific_ingroup_label}."
            )
        else:
            evidence += (
                f"Subsamples scored unusually low against "
                f"{specific_ingroup_label} "
                f"(threshold={threshold:.4f})."
            )

    summary = {
        "query_accession":          query_accession,
        "query_taxon":              query_taxon,
        "query_specific_group":     specific_ingroup_label,
        "scoring_level":            scoring_level,
        "fallback_used":            fallback_used,
        "is_orphan":                is_orphan,
        "expected_order":           expected_order,
        "n_samples":                n_samples,
        "subsample_length":         n,
        "k_nearest":                k_nearest,
        "outlier_threshold":        round(threshold, 4),
        "n_outliers":               len(outliers),
        "prop_outliers":            round(len(outliers) / n_samples,
                                         4),
        "n_correct_group":          len(correct),
        "prop_correct_group":       round(len(correct) / n_samples,
                                         4),
        "specific_ingroup_mean":    round(sp_mean, 4),
        "specific_ingroup_std":     round(sp_std, 4),
        "specific_ingroup_variance":sp_variance,
        "per_group_means":          per_group_means,
        "chimera_flag":             chimera_flag,
        "chimera_evidence":         evidence,
        "scoring_groups":           {
            grp: [rid for rid, _ in pairs]
            for grp, pairs in scoring_groups.items()
        },
        "sample_scores":            sample_scores,
    }

    print(f"\n── Ingroup distance summary ─────────────────────────────────")
    print(f"  Specific ingroup         : "
          f"{specific_ingroup_label} [{scoring_level}]")
    print(f"  Fallback used            : {fallback_used}")
    print(f"  Is orphan                : {is_orphan}")
    print(f"  Specific ingroup mean    : {sp_mean:.4f}")
    print(f"  Specific ingroup std     : {sp_std:.4f}")
    print(f"  Specific ingroup variance: {sp_variance:.6f}")
    print(f"  Outlier threshold        : {threshold:.4f}")
    print(f"  Correct group matches    : "
          f"{len(correct)}/{n_samples}")
    print(f"  Outlier subsamples       : "
          f"{len(outliers)}/{n_samples}")
    print(f"  Chimera flag             : {chimera_flag}")
    print(f"  Evidence                 : {evidence}")
    print(f"\n  Per-group mean knn identity:")
    for grp, mean in sorted(per_group_means.items(),
                             key=lambda x: -(x[1] or 0)):
        marker = " ◄ EXPECTED" if grp == specific_ingroup_label else ""
        print(f"    {grp:<30} {mean:.4f}{marker}")

    return summary


def verify_batch_by_ingroup_distance(
        query_fasta: str,
        scaffold_fasta: str,
        taxonomy: dict,
        output_dir: str,
        insert_fraction: float = 0.30,
        target_detection_prob: float = 0.90,
        k_nearest: int = 3,
        min_specific_ingroup_identity: float = None,
        specific_ingroup_sd_threshold: float = 3.0,
        wrong_group_margin_threshold: float = 0.02,
        known_orphan_taxa: list = None,
        orphan_sd_threshold: float = 3.0,
        threshold_json: str = None) -> dict:
    """
    Run verify_by_ingroup_distance() on every sequence in a FASTA.

    FASTA header format:
        >Genus_species|gene|orien:+/-|accession:XXXX

    Parameters
    ----------
    query_fasta                  : str   FASTA of sequences to verify
    scaffold_fasta               : str   scaffold FASTA for the gene
    taxonomy                     : dict  from load_taxonomy()
    output_dir                   : str   directory for output CSVs
    insert_fraction              : float passed to verify function
    target_detection_prob        : float passed to verify function
    k_nearest                    : int   passed to verify function
    specific_ingroup_sd_threshold: float passed to verify function

    Returns
    -------
    list of summary dicts, one per sequence
    """
    import json

    # ── Load gene-specific thresholds if provided ─────────────────────────────
    gene_thresholds = {}
    if threshold_json and os.path.exists(threshold_json):
        with open(threshold_json, "r") as f:
            gene_thresholds = json.load(f)
        print(f"Loaded gene-specific thresholds: "
              f"{len(gene_thresholds)} genes from {threshold_json}")
    elif threshold_json:
        print(f"Warning: threshold_json path not found: "
              f"{threshold_json} — using default thresholds")
                
    os.makedirs(output_dir, exist_ok=True)

    records = list(SeqIO.parse(query_fasta, "fasta"))
    total   = len(records)
    results = []

    print(f"\nBatch ingroup distance verification: {total} sequences")
    print(f"Scaffold : {scaffold_fasta}")
    print(f"Output   : {output_dir}\n")

    for i, record in enumerate(records, 1):
            parts     = record.description.split("|")
            taxon     = parts[0].strip()
            gene      = parts[1].strip() if len(parts) > 1 else "unknown"
            accession = None

            for p in parts:
                if p.startswith("accession:"):
                    accession = p.split(":", 1)[1].strip()

            if not accession:
                print(f"  [{i}/{total}] Skipping {taxon} — "
                    f"no accession in header")
                continue

            print(f"\n[{i}/{total}] {taxon} | {gene} | {accession}")

            # ── Look up gene-specific thresholds if available ─────────────────────
            if gene in gene_thresholds:
                sd_thresh     = gene_thresholds[gene]["sd_threshold"]
                margin_thresh = gene_thresholds[gene]["margin_threshold"]
                print(f"  Using calibrated thresholds for {gene}: "
                    f"sd={sd_thresh}, margin={margin_thresh}")
            else:
                sd_thresh     = specific_ingroup_sd_threshold
                margin_thresh = wrong_group_margin_threshold
                if gene_thresholds:
                    print(f"  No calibrated threshold for {gene} — "
                        f"using defaults: "
                        f"sd={sd_thresh}, margin={margin_thresh}")

            try:
                result = verify_by_ingroup_distance(
                    query_seq                     = str(record.seq),
                    query_taxon                   = taxon,
                    query_accession               = accession,
                    scaffold_fasta                = scaffold_fasta,
                    taxonomy                      = taxonomy,
                    insert_fraction               = insert_fraction,
                    target_detection_prob         = target_detection_prob,
                    k_nearest                     = k_nearest,
                    specific_ingroup_sd_threshold = sd_thresh,
                    wrong_group_margin_threshold  = margin_thresh,
                    known_orphan_taxa             = known_orphan_taxa,
                    orphan_sd_threshold           = orphan_sd_threshold,
                )
                result["gene"]    = gene
                result["_header"] = record.description
                results.append(result)

                print(f"  ✓ chimera_flag={result['chimera_flag']} | "
                    f"group={result['query_specific_group']} | "
                    f"level={result['scoring_level']} | "
                    f"outliers={result['n_outliers']}/"
                    f"{result['n_samples']} | "
                    f"evidence={result['chimera_evidence'][:80]}")

            except Exception as e:
                print(f"  ✗ Failed: {e}")
                import traceback
                traceback.print_exc()

    # ── Write summary CSV ─────────────────────────────────────────────────────
    if results:
        summary_path = os.path.join(
            output_dir,
            "ingroup_verification_summary.csv"
        )
        exclude = {"sample_scores", "family_groups",
                   "per_family_means"}
        keys    = [k for k in results[0].keys()
                   if k not in exclude]

        import csv
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for r in results:
                writer.writerow({k: r[k] for k in keys})

        print(f"\nSummary: {summary_path}")

        flagged = [r for r in results if r["chimera_flag"]]
        print(f"\n{'⚠' if flagged else '✓'} "
              f"{len(flagged)}/{len(results)} sequences flagged")
        for r in flagged:
            print(f"  {r['query_taxon']} ({r['query_accession']}) | "
                  f"{r.get('gene','?')} | "
                  f"{r['chimera_evidence'][:100]}")

    return results


#print(verify_by_ingroup_distance("ACTTTGGTTGCATGAAGGCTGCCCCCATGAAAGAAGCACACTTCCGAGGACAAGGCAGCTTGGCCTACCCAGGTCTGCGGACCCATGGGACTCTGGAGAGCGCAAATGGGCCCAAGGCAAGTTCAAGAGACCTGGCGTTGGCTAGCACTTTTGAACATGTGCTGGAAGAGCTGTTGGACGAGGACCAGAAGATTCGTCCCCATGAAGAAACCCCTAAGGACGCGGACTTGTATACTTCCCGAGTGATGCTCAGCAGTCAAGTGCCTTTGGAGCCACCACTTCTCTTTCTGCTTGAGGAATACAAAAATTACCTGGATGCTGCAAACATGTCGATGAGGGTCCGACGCCACTCCGACCCTGCCCGCCGTGGGGAGCTGAGCGTGTGCGACAGCGTTAGCCAGTGGGTGACAGCAGCAGATAAAAAGACTGCAGTGGACATGTCGGGCGGGACGGTCACGGTCCTGGAAAAGGTCCCTGTGTCCAAAGGCCAACTGAAGCAGTACTTCTACGAGACCAGGTGCAATCCCCTGGGTTTCACGAAGGAAGGCTGCAGG",
#                         "Podogymnura_truei",
#                         "JN633375",
#                         "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas/BDNF.fasta",
#                         load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
#))


#print(verify_by_ingroup_distance("ACTTCGGTTGCATGAAGGCTGCCCCCATGAAAGAAGCCAGTGTCCGAGGACCAGGCAGCTTGGCCTACCCAGGTGTGCGGACCCATGGGACTCTGGAGAGCGTGAATGGGCCCAAGGCAGGTTCGAGAGGCCTGACTTTGGCTGACACTTTTGAACACGTGATAGAAGAGCTCCTGGATGAGGACCAGAAAGTCCGGCCCCACGAAGAGAACAATAAGGACGCGGACTTGTACACCTCCAGGGTGATGCTCAGTAGTCAAGTGCCTTTGGAGCCGCCTCTCCTCTTTCTGCTCGAGGAATACAAAAATTACCTGGATGCTGCAAACATGTCTATGCGGGTCCGGCGCCACTCCGACCCTGCCCGCCGCGGGGAGCTGAGCGTGTGCGACAGCATTAGCGAGTGGGTGACGGCGGCGGATAAAAAGACTGCAGTGGACATGTCGGGCGGGACGGTGACGGTCCTGGAGAAAGTCCCTGTATCGAAAGGCCAACTGAAGCAGTACTTCTACGAGACCAAGTGCAATCCCATGGGTTACACAAAGGAGGGCTGCAGG",
#                         "Uropsilus_soricipes",
#                         "KF778036",
#                         "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas/BDNF.fasta",
#                         load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
#))


#print(verify_by_distance("ACTTCGGTTGCATGAAGGCTGCCCCCATGAAAGAAGCCAGCGTCCGAGGACAAGGCAGCTTGGCCTACCCAGGTGTGCGGACCCATGGGACTCTGGAGAGTGTGAATGGGCCCAAGGCAGGTGCCAGAGGCCTGACGTCCTTGGCTGACACTTTTGAACACGTGATCGAAGAGCTGTTGGAAGAGGACCAGAAAGTTCGTCCCCATGAAGAAACCAATAAGGACGCGGACTTGTACACTTCCCGGGTGATGCTGAGTAGTCAAGTGCCTTTGGAGCCTCCTCTTCTCTTTCTGCTGGAGGAATACAAAAATTACCTGGATGCTGCAAACATGTCCATGAGGGTCCGGCGCCACTCCGACCCCGCCCGCCGCGGGGAGCTGAGCGTGTGTGACAGCATCAGCGAGTGGGTGACAGCAGCGGATAAAAAGACTGCAGTGGACATGTCGGGCGGGACGGTCACTGTCCTGGAAAAAGTCCCTGTATCCAAAGGCCAACTGAAGCAGTACTTCTACGAGACCAAGTGCAATCCCATGGGTTACACGAAGGAGGGCTGCAGG",
#                         "Solenodon_paradoxus",
#                         "AY530070",
#                         "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas/BDNF.fasta",
#                         load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
#))

#print(verify_by_distance("ATAATTAAAGGTCTGGTCCCAGCCTTCCTATTTTCTATTAGTAGAATTACACATGCAAGTATCAGCTACCCAGTGCGAATGCCCTCTAACCCTACCATTAATAGGTGTAAAGGAGCGGATATCAAGTACACACATATGTTGCTAATGACATCTTGCTTAACCACACCCCCACGGGAAACAGCAGTGATAAATATTGAGCTATAAACGAAAGTTTGACTAAGCCATATTAATTTAGGGTTGGTAAATCTCGTCCGAGCCACCGCGGTCATACGATTAACCCATGAGAATAGGAAATCGGCGTAAAGAGTGTTTAGGATATTAATGTAATGAAATTAAAAAATGACTTAGCTGTAAAAAGCTCATTTCATAAATAAAAACATCTACAAAAGTGATTTCATAGGATCTTATTACACGTGAGCTAAGACCCAAACTAGGATTAGATACCCTATTATGCTTAGCCCTAAACTTAGACAGTTACTATTTGCCAGAGAACTACTAGCCATAGCTTAAAACTCAAAGGACTTGGCGGTACTTTATATCCATCTAGAGGAGCCTGTTCTATAATCGATAAACCCCGCTCTACCTCACCATCTCTTGCTAATTCAGCCTATATACCGCCATCTTCAGCAAACCCTAAAAAGGTATTAAAGTAAGCAAAAGAATCAAACATAAAAACGTTAGGTCAAGGTGTAGCCAATGAAATGGGAAGAAATGGGCTACATTTTCTTATAAAAGAACATTACTATACCCTTTATGAAACTAAAGGATTAAGGAGGATTTAGTAGTAAATTAAGAATAGAGAGCTTAATTGAATTGAGCAATTTGGCAATGAAGCATGCACACACCGCCCGTCACCCTCTTCAAGCATATTAAGTCACCAACCTATATAATTAATGTTATAATGATAATCACATGCAAGAAGAGATAAGTCGTAACAAGGTAAGTATACTGGAAAGTGTACTTGGATTAT",
#                         "Erinaceus_europaeus",
#                         "NC_002080",
#                         "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas/BDNF.fasta",
#                         load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
#))