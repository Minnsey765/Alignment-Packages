# _distance_verify.py

import os
import csv
import random
import numpy as np

from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner

from bioinf_packages.verify_funcs._build_reference_fastas import (
    build_reference_fastas,
    build_combined_fasta,
    add_calibration_to_qualified,
    _parse_pipe_header,
    resolve_accession_keys,
)
from bioinf_packages.verify_funcs._qualify_cleared import (
    build_gbk_index,
    qualify_cleared_sequences,
)

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
# ─────────────────────────────────────────────────────────────────────────────

def get_expected_order(taxon: str, taxonomy: dict) -> str:
    if isinstance(taxonomy, str):
        raise TypeError(
            f"taxonomy must be a dict from load_taxonomy(), "
            f"not a file path."
        )
    genus    = taxon.split("_")[0]
    tax_info = taxonomy.get(genus)
    if tax_info:
        return tax_info.get("order")
    return None


def get_ingroup_genera(order: str, taxonomy: dict) -> set:
    return {
        genus for genus, info in taxonomy.items()
        if info.get("order") == order
    }


# ─────────────────────────────────────────────────────────────────────────────
# PAIRWISE ALIGNER
# ─────────────────────────────────────────────────────────────────────────────

def _build_semiglobal_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode                         = "global"
    aligner.match_score                  =  1.0
    aligner.mismatch_score               = -1.0
    aligner.open_deletion_score          = -2.0
    aligner.extend_deletion_score        = -0.5
    aligner.open_left_insertion_score    =  0.0
    aligner.extend_left_insertion_score  =  0.0
    aligner.open_right_insertion_score   =  0.0
    aligner.extend_right_insertion_score =  0.0
    aligner.open_internal_insertion_score   = -2.0
    aligner.extend_internal_insertion_score = -0.5
    return aligner


def semiglobal_identity(query: str,
                         target: str,
                         aligner: PairwiseAligner) -> float:
    if not query or not target:
        return 0.0
    try:
        alignments = aligner.align(target, query)
        if not alignments:
            return 0.0
        best    = alignments[0]
        aligned = best.aligned
        n_identical = 0
        for (t_start, t_end), (q_start, q_end) in zip(
                aligned[0], aligned[1]):
            t_block = target[t_start:t_end]
            q_block = query[q_start:q_end]
            n_identical += sum(t == q for t, q in
                               zip(t_block, q_block))
        return n_identical / len(query)
    except Exception:
        return 0.0


# ─────────────────────────────────────────────────────────────────────────────
# ORDER-LEVEL DISTANCE VERIFICATION (legacy — kept for reference)
# ─────────────────────────────────────────────────────────────────────────────

def verify_by_distance(query_seq, query_taxon, query_accession,
                        scaffold_fasta, taxonomy,
                        insert_fraction=0.30,
                        target_detection_prob=0.90,
                        outlier_ratio_threshold=0.95,
                        min_ratio_variance_flag=0.001,
                        k_nearest=3):
    # [function body unchanged from your current file]
    pass


def verify_batch_by_distance(query_fasta, scaffold_fasta, taxonomy,
                              output_dir, insert_fraction=0.30,
                              target_detection_prob=0.90):
    # [function body unchanged from your current file]
    pass


def verify_by_ingroup_distance(query_seq, query_taxon, query_accession,
                                scaffold_fasta, taxonomy,
                                insert_fraction=0.30,
                                target_detection_prob=0.90,
                                k_nearest=3,
                                min_specific_ingroup_identity=None,
                                specific_ingroup_sd_threshold=3.0,
                                wrong_group_margin_threshold=0.02,
                                known_orphan_taxa=None,
                                orphan_sd_threshold=3.0,
                                threshold_json=None):
    # [function body unchanged from your current file]
    pass


def verify_batch_by_ingroup_distance(query_fasta, scaffold_fasta,
                                      taxonomy, output_dir,
                                      insert_fraction=0.30,
                                      target_detection_prob=0.90,
                                      k_nearest=3,
                                      min_specific_ingroup_identity=None,
                                      specific_ingroup_sd_threshold=3.0,
                                      wrong_group_margin_threshold=0.02,
                                      known_orphan_taxa=None,
                                      orphan_sd_threshold=3.0,
                                      threshold_json=None):
    # [function body unchanged from your current file]
    pass


# ─────────────────────────────────────────────────────────────────────────────
# FAMILY DISTANCE VERIFICATION
# ─────────────────────────────────────────────────────────────────────────────

def verify_subsamples_by_family_distance(
        subsample_seqs: list,
        query_taxon: str,
        gene: str,
        reference_sequences: dict,
        taxonomy: dict,
        family_distance_threshold: float = 0.15,
        subfamily_z_threshold: float = 2.5,
        min_references_family: int = 3,
) -> dict:
    """
    Check each subsample independently for family grouping.

    Confidence tiers
    ----------------
    CONFIDENT_CHIMERA
        One or more subsamples do not group within the correct family.
    PROBABLE_CHIMERA
        All subsamples within family but high distance variance.
    NOT_DETECTED
        All subsamples group consistently within the family.
    INSUFFICIENT_DATA
        Too few family references or family not in taxonomy.
    """
    import numpy as np

    def p_distance(seq_a, seq_b):
        comparable = [(a, b) for a, b in zip(seq_a, seq_b)
                      if a != '-' and b != '-']
        if not comparable:
            return None
        return sum(1 for a, b in comparable
                   if a != b) / len(comparable)

    def get_family(taxon):
        genus = taxon.split("_")[0] if "_" in taxon else taxon
        tax   = taxonomy.get(genus, {})
        return (tax.get("family") or tax.get("Family")
                or tax.get("FAMILY") or "").strip()

    def get_subfamily(taxon):
        genus = taxon.split("_")[0] if "_" in taxon else taxon
        tax   = taxonomy.get(genus, {})
        return (tax.get("subfamily") or tax.get("Subfamily")
                or tax.get("sub-family") or "").strip()

    query_family    = get_family(query_taxon)
    query_subfamily = get_subfamily(query_taxon)

    if not query_family:
        return {
            "confidence":      "INSUFFICIENT_DATA",
            "note":            (f"Family not found in taxonomy for "
                                f"'{query_taxon.split('_')[0]}'. "
                                f"Check taxonomy CSV column name."),
            "query_family":    query_family,
            "query_subfamily": query_subfamily,
            "query_taxon":     query_taxon,
            "gene":            gene,
        }

    family_refs = {
        t: s for t, s in reference_sequences.items()
        if t != query_taxon and get_family(t) == query_family
    }

    if len(family_refs) < min_references_family:
        return {
            "confidence":      "INSUFFICIENT_DATA",
            "note":            (f"Only {len(family_refs)} family "
                                f"references for {query_family} "
                                f"(need {min_references_family})."),
            "query_family":    query_family,
            "query_subfamily": query_subfamily,
            "query_taxon":     query_taxon,
            "gene":            gene,
            "n_family_refs":   len(family_refs),
        }

    subsample_results = []
    n_out_of_family   = 0

    for i, sub_seq in enumerate(subsample_seqs):
        dists = {}
        for ref_taxon, ref_seq in family_refs.items():
            d = p_distance(sub_seq, ref_seq)
            if d is not None:
                dists[ref_taxon] = d

        if not dists:
            continue

        min_dist      = min(dists.values())
        min_dist_ref  = min(dists, key=dists.get)
        mean_dist     = float(np.mean(list(dists.values())))
        out_of_family = min_dist > family_distance_threshold

        if out_of_family:
            n_out_of_family += 1

        subsample_results.append({
            "subsample_index": i,
            "min_dist":        round(min_dist, 4),
            "min_dist_ref":    min_dist_ref,
            "mean_dist":       round(mean_dist, 4),
            "out_of_family":   out_of_family,
            "distances":       {t: round(d, 4)
                                for t, d in sorted(dists.items())},
        })

    if not subsample_results:
        return {
            "confidence":      "INSUFFICIENT_DATA",
            "note":            "No comparable sites found.",
            "query_family":    query_family,
            "query_subfamily": query_subfamily,
            "query_taxon":     query_taxon,
            "gene":            gene,
        }

    prop_out = n_out_of_family / len(subsample_results)
    min_dists = [r["min_dist"] for r in subsample_results]
    dist_var  = float(np.var(min_dists))  if len(min_dists) > 1 else 0.0
    dist_std  = float(np.std(min_dists))  if len(min_dists) > 1 else 0.0

    if n_out_of_family > 0:
        confidence = "CONFIDENT_CHIMERA"
        note = (f"{n_out_of_family}/{len(subsample_results)} subsamples "
                f"do not group within {query_family} "
                f"(min_dist > {family_distance_threshold:.4f}).")
    elif dist_std > subfamily_z_threshold * 0.05:
        confidence = "PROBABLE_CHIMERA"
        note = (f"All subsamples within {query_family} but high "
                f"variance (std={dist_std:.4f}).")
    else:
        confidence = "NOT_DETECTED"
        note = (f"All {len(subsample_results)} subsamples group "
                f"within {query_family} (std={dist_std:.4f}).")

    return {
        "confidence":         confidence,
        "prop_out_of_family": round(prop_out, 4),
        "n_out_of_family":    n_out_of_family,
        "n_subsamples":       len(subsample_results),
        "n_family_refs":      len(family_refs),
        "dist_variance":      round(dist_var, 4),
        "dist_std":           round(dist_std, 4),
        "subsample_results":  subsample_results,
        "query_family":       query_family,
        "query_subfamily":    query_subfamily,
        "query_taxon":        query_taxon,
        "gene":               gene,
        "note":               note,
    }


def calibrate_family_distance_thresholds(
        qualified: dict,
        alignment_folder: str,
        taxonomy: dict,
        min_taxa_per_family: int = 3,
        percentile_threshold: float = 95.0,
        output_path: str = None,
        file_pattern: str = "*.fasta",
) -> dict:
    """
    Calibrate family distance thresholds from cleared sequences.

    Family membership is derived from alignment file headers, not
    from the qualified dict keys (which may be accessions).

    Parameters
    ----------
    qualified           : dict   taxon or accession -> list of genes.
                                 Not used for family grouping —
                                 alignment headers are used instead.
    alignment_folder    : str    directory with aligned FASTA files
    taxonomy            : dict   from load_taxonomy()
    min_taxa_per_family : int    minimum taxa per family
    percentile_threshold: float  percentile for threshold (default 95)
    output_path         : str    optional JSON output path
    file_pattern        : str    glob pattern for alignment files.
                                 Default "*.fasta". Use "*_aln.fasta"
                                 to exclude short_aln files when both
                                 are in the same folder.

    Returns
    -------
    dict with family_distance_threshold, dist_variance_threshold,
         per_gene_thresholds, per_family_thresholds,
         calibration_distances, n_calibration_sequences,
         n_families, family_stats, percentile_used
    """
    import numpy as np
    import json
    from collections import defaultdict

    rrna_map = {"12S_RRNA": "12S_rRNA", "16S_RRNA": "16S_rRNA"}

    def get_family(taxon):
        genus = taxon.split("_")[0] if "_" in taxon else taxon
        tax   = taxonomy.get(genus, {})
        return (tax.get("family") or tax.get("Family")
                or tax.get("FAMILY") or "").strip()

    def p_distance(seq_a, seq_b):
        comparable = [(a, b) for a, b in zip(seq_a, seq_b)
                      if a != '-' and b != '-']
        if not comparable:
            return None
        return sum(1 for a, b in comparable
                   if a != b) / len(comparable)

    sample_genus = next(iter(taxonomy), None)
    if sample_genus:
        print(f"\n  Taxonomy key check for '{sample_genus}': "
              f"{taxonomy[sample_genus]}")

    print(f"\n{'='*60}")
    print(f"CALIBRATE FAMILY DISTANCE THRESHOLDS")
    print(f"{'='*60}")
    print(f"  Alignment folder    : {alignment_folder}")
    print(f"  Qualified entries   : {len(qualified)}")
    print(f"  Percentile          : {percentile_threshold}")
    print(f"  File pattern        : {file_pattern}")

    aln_folder = Path(alignment_folder)
    alignments = {}

    for aln_file in sorted(aln_folder.glob(file_pattern)):
        gene = (aln_file.stem
                .replace("_aligned", "")
                .replace("_aln", "")
                .replace("_short", "")
                .upper())
        gene = rrna_map.get(gene, gene)

        seqs = {}
        try:
            for record in SeqIO.parse(str(aln_file), "fasta"):
                taxon_key = None
                parts     = record.id.split("_")
                for i in range(len(parts) - 1):
                    if (parts[i] and parts[i][0].isupper()
                            and parts[i].isalpha()
                            and parts[i+1]
                            and parts[i+1][0].islower()
                            and parts[i+1].isalpha()):
                        taxon_key = f"{parts[i]}_{parts[i+1]}"
                        break
                if taxon_key is None:
                    taxon_key = (f"{parts[0]}_{parts[1]}"
                                 if len(parts) >= 2 else record.id)
                seqs[taxon_key] = str(record.seq)
        except Exception as e:
            print(f"  Warning: could not load {aln_file.name}: {e}")
            continue

        if seqs:
            alignments[gene] = seqs
            print(f"  Loaded {gene:<15} "
                  f"{len(seqs)} taxa from {aln_file.name}")

    print(f"\n  Total alignment files loaded: {len(alignments)}")
    if not alignments:
        raise FileNotFoundError(
            f"No FASTA files matching '{file_pattern}' in "
            f"{alignment_folder}"
        )

    # Group taxa by family from alignment headers
    alignment_taxa = set()
    for seqs in alignments.values():
        alignment_taxa.update(seqs.keys())

    print(f"\n  Unique taxa in alignments: {len(alignment_taxa)}")

    family_taxa    = defaultdict(set)
    unknown_genera = set()

    for taxon in alignment_taxa:
        family = get_family(taxon)
        if family:
            family_taxa[family].add(taxon)
        else:
            unknown_genera.add(
                taxon.split("_")[0] if "_" in taxon else taxon
            )

    if unknown_genera:
        print(f"  Warning: {len(unknown_genera)} genera not in "
              f"taxonomy: {', '.join(sorted(unknown_genera)[:10])}")

    eligible_families = {
        fam: taxa for fam, taxa in family_taxa.items()
        if len(taxa) >= min_taxa_per_family
    }

    print(f"\n  Families with >= {min_taxa_per_family} taxa: "
          f"{len(eligible_families)}")
    for fam, taxa in sorted(eligible_families.items()):
        print(f"    {fam:<30} {len(taxa)} taxa: "
              f"{', '.join(sorted(taxa))}")

    if not eligible_families:
        raise ValueError(
            f"No families have >= {min_taxa_per_family} taxa in the "
            f"alignments. Taxonomy column keys found: "
            f"{list(next(iter(taxonomy.values()), {}).keys())}. "
            f"Add more sequences to the reference alignment or lower "
            f"min_taxa_per_family."
        )

    all_min_distances    = []
    per_gene_distances   = defaultdict(list)
    per_family_distances = defaultdict(list)
    all_variances        = []
    n_computed           = 0
    family_stats         = {}

    for family, taxa in sorted(eligible_families.items()):
        family_min_dists = []

        for taxon in sorted(taxa):
            genes_in_aln = [
                g for g, seqs in alignments.items()
                if taxon in seqs
            ]

            for gene_upper in genes_in_aln:
                aln       = alignments[gene_upper]
                query_seq = aln.get(taxon)
                if query_seq is None:
                    continue

                family_dists = []
                for other_taxon in taxa:
                    if other_taxon == taxon:
                        continue
                    other_seq = aln.get(other_taxon)
                    if other_seq is None:
                        continue
                    d = p_distance(query_seq, other_seq)
                    if d is not None:
                        family_dists.append(d)

                if not family_dists:
                    continue

                min_d = min(family_dists)
                std_d = float(np.std(family_dists))

                all_min_distances.append(min_d)
                per_gene_distances[gene_upper].append(min_d)
                per_family_distances[family].append(min_d)
                family_min_dists.append(min_d)
                if std_d > 0:
                    all_variances.append(std_d)
                n_computed += 1

        if family_min_dists:
            family_stats[family] = {
                "n_sequences":   len(family_min_dists),
                "mean_min_dist": round(float(np.mean(
                    family_min_dists)), 4),
                "std_min_dist":  round(float(np.std(
                    family_min_dists)), 4),
                "max_min_dist":  round(float(np.max(
                    family_min_dists)), 4),
                f"p{int(percentile_threshold)}": round(
                    float(np.percentile(
                        family_min_dists, percentile_threshold)), 4),
            }

    print(f"\n  Computed {n_computed} within-family distances "
          f"across {len(eligible_families)} families")

    if not all_min_distances:
        raise ValueError(
            "No distances computed — check alignment taxon names "
            "match taxonomy CSV."
        )

    global_threshold = float(np.percentile(
        all_min_distances, percentile_threshold
    ))

    per_gene_thresholds = {}
    for gene, dists in sorted(per_gene_distances.items()):
        if len(dists) >= 5:
            per_gene_thresholds[gene] = round(
                float(np.percentile(dists, percentile_threshold)), 4
            )

    per_family_thresholds = {}
    for family, dists in sorted(per_family_distances.items()):
        if len(dists) >= 5:
            per_family_thresholds[family] = round(
                float(np.percentile(dists, percentile_threshold)), 4
            )

    dist_variance_threshold = (
        float(np.percentile(all_variances, percentile_threshold))
        if all_variances else 0.05
    )

    print(f"\n── Calibrated thresholds ────────────────────────────────")
    print(f"  Global threshold : {global_threshold:.4f}")
    print(f"  Variance thresh  : {dist_variance_threshold:.4f}")
    print(f"\n  Per-gene thresholds:")
    for gene, thresh in sorted(per_gene_thresholds.items()):
        print(f"    {gene:<15} {thresh:.4f}")
    print(f"\n  Per-family thresholds:")
    for fam, thresh in sorted(per_family_thresholds.items()):
        print(f"    {fam:<30} {thresh:.4f}")
    print(f"\n  n={len(all_min_distances)} "
          f"mean={np.mean(all_min_distances):.4f} "
          f"max={np.max(all_min_distances):.4f} "
          f"p{int(percentile_threshold)}={global_threshold:.4f}")

    result = {
        "family_distance_threshold":  round(global_threshold, 4),
        "dist_variance_threshold":    round(dist_variance_threshold, 4),
        "per_gene_thresholds":        per_gene_thresholds,
        "per_family_thresholds":      per_family_thresholds,
        "calibration_distances":      [round(d, 4)
                                       for d in all_min_distances],
        "n_calibration_sequences":    n_computed,
        "n_families":                 len(eligible_families),
        "family_stats":               family_stats,
        "percentile_used":            percentile_threshold,
        "alignment_folder":           str(alignment_folder),
        "min_taxa_per_family":        min_taxa_per_family,
    }

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            import json as _json
            _json.dump(result, f, indent=2)
        print(f"\n  Thresholds saved: {output_path}")

    return result


def get_subsamples_from_csv(
        results_csv: str,
        taxon: str,
        gene: str,
        label_filter: str = "GENUINE",
) -> dict:
    """
    Extract subsample coordinate info for a (taxon, gene) pair.

    Tries matching on taxon column first, then accession column,
    then accession without version suffix (e.g. LC124901.1 ->
    LC124901).
    """
    def _scan(f, match_field, match_value):
        f.seek(0)
        reader     = csv.DictReader(f)
        subsamples = []
        seen       = set()
        for row in reader:
            if row.get("label", "").strip() != label_filter:
                continue
            if row.get(match_field, "").strip() != match_value:
                continue
            if row.get("gene", "").strip() != gene:
                continue
            if row.get("hit_accession", "").strip() == "NO_HITS":
                continue
            sample = row.get("sample", "").strip()
            if sample in seen:
                continue
            seen.add(sample)
            try:
                start = int(row["subsample_start"])
                end   = int(row["subsample_end"])
            except (KeyError, ValueError, TypeError):
                continue
            subsamples.append({
                "sample":          sample,
                "subsample_start": start,
                "subsample_end":   end,
                "taxon_in_csv":    row.get("taxon", "").strip(),
                "accession":       row.get("accession", "").strip(),
            })
        return subsamples

    with open(results_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        subsamples = _scan(f, "taxon", taxon)
        if not subsamples:
            subsamples = _scan(f, "accession", taxon)
        if not subsamples:
            taxon_base = taxon.split(".")[0]
            if taxon_base != taxon:
                subsamples = _scan(f, "accession", taxon_base)

    return {
        "subsamples":   subsamples,
        "n_subsamples": len(subsamples),
        "taxon":        taxon,
        "gene":         gene,
    }


def extract_subsample_sequences(
        subsample_info: dict,
        fasta_folder: str,
        taxon: str,
        gene: str,
) -> list:
    """
    Extract nucleotide sequences for each subsample window from the
    original per-gene FASTA files using coordinates from the CSV.

    fasta_folder should point to the original per-gene FASTA files
    (e.g. fasta_data/CYTB.fasta) not the subsample batch files.
    """
    import glob

    taxon_in_csv = None
    accession    = None
    if subsample_info.get("subsamples"):
        first        = subsample_info["subsamples"][0]
        taxon_in_csv = first.get("taxon_in_csv", "").strip() or None
        accession    = first.get("accession", "").strip() or None

    search_terms = []
    if taxon_in_csv:
        search_terms.append(taxon_in_csv)
    if accession:
        search_terms.append(accession)
    if taxon not in search_terms:
        search_terms.append(taxon)

    parent_seq = None
    search_patterns = [
        f"{fasta_folder}/**/*{gene}*.fasta",
        f"{fasta_folder}/*{gene}*.fasta",
        f"{fasta_folder}/**/*.fasta",
        f"{fasta_folder}/*.fasta",
    ]

    seen_files = set()
    for pattern in search_patterns:
        if parent_seq:
            break
        for match in glob.glob(pattern, recursive=True):
            if match in seen_files:
                continue
            seen_files.add(match)
            try:
                for record in SeqIO.parse(match, "fasta"):
                    record_str = record.id + " " + record.description
                    if any(term in record_str
                           for term in search_terms if term):
                        parent_seq = str(record.seq).replace("-", "")
                        break
            except Exception:
                continue
            if parent_seq:
                break

    if parent_seq is None:
        return []

    sequences = []
    for sub in subsample_info["subsamples"]:
        start = sub["subsample_start"] - 1
        end   = sub["subsample_end"]
        if start >= 0 and end <= len(parent_seq) and start < end:
            sequences.append(parent_seq[start:end])

    return sequences


def profile_align_subsample(
        subsample_seq: str,
        reference_alignment: dict,
) -> str:
    """
    Insert a short subsample into an existing MSA by pairwise
    alignment to the reference consensus.
    """
    from collections import Counter

    try:
        from Bio import pairwise2
    except ImportError:
        return subsample_seq

    if not reference_alignment:
        return subsample_seq

    ref_seqs  = list(reference_alignment.values())
    aln_len   = len(ref_seqs[0])
    consensus = []

    for pos in range(aln_len):
        chars = [s[pos] for s in ref_seqs
                 if pos < len(s) and s[pos] != '-']
        if chars:
            consensus.append(Counter(chars).most_common(1)[0][0])
        else:
            consensus.append('N')

    consensus_str = "".join(consensus)

    alignments = pairwise2.align.localms(
        subsample_seq.upper(), consensus_str.upper(),
        2, -1, -2, -0.5, one_alignment_only=True,
    )

    if not alignments:
        return subsample_seq

    aligned_sub, aligned_ref, score, begin, end = alignments[0]

    ref_pos     = 0
    sub_aligned = []

    for sub_char, ref_char in zip(aligned_sub, aligned_ref):
        if ref_char == '-':
            continue
        elif sub_char == '-':
            sub_aligned.append((ref_pos, '-'))
            ref_pos += 1
        else:
            sub_aligned.append((ref_pos, sub_char))
            ref_pos += 1

    result = ['-'] * aln_len
    for col, char in sub_aligned:
        if col < aln_len:
            result[col] = char

    return "".join(result)


def run_family_distance_verification(
        flagged_list: list,
        results_csv: str,
        fasta_folder: str,
        alignment_folder: str,
        taxonomy: dict,
        thresholds: dict,
        label_filter: str = "GENUINE",
        output_path: str = None,
) -> dict:
    """
    Run family distance verification on flagged sequences.

    Accepts ALL flagged sequences (LIKELY_TP, REVIEW, LIKELY_FP,
    ORPHAN_FP) from the BLAST model — not just REVIEW. This allows
    the family distance test to rescue false positives from any
    category and confirm true positives.

    For sequences where INSUFFICIENT_DATA is returned (typically
    orphan taxa like Solenodon with no family representatives in the
    reference alignment), the result is recorded but no chimera call
    is made by this function — those sequences should be handled
    separately using PubMed qualification.

    Parameters
    ----------
    flagged_list     : list   (taxon, gene) tuples — all flagged
                              sequences from result["flagged"]
    results_csv      : str    path to blast_results.csv
    fasta_folder     : str    original per-gene FASTA folder
                              (NOT subsample batch files)
    alignment_folder : str    per-gene reference MSA folder
    taxonomy         : dict   from load_taxonomy()
    thresholds       : dict   from calibrate_family_distance_thresholds()
    label_filter     : str    label in blast_results.csv
    output_path      : str    optional JSON output path

    Returns
    -------
    dict with results, confident_chimera, probable_chimera,
              not_detected, insufficient
    """
    import json

    rrna_map = {"12S_RRNA": "12S_rRNA", "16S_RRNA": "16S_rRNA"}

    print(f"\n{'='*60}")
    print(f"FAMILY DISTANCE VERIFICATION")
    print(f"{'='*60}")
    print(f"  Flagged sequences: {len(flagged_list)}")
    print(f"  Alignment folder : {alignment_folder}")
    print(f"  FASTA folder     : {fasta_folder}")
    print(f"  Label filter     : {label_filter}")

    # Load reference alignments once per gene
    aln_folder     = Path(alignment_folder)
    ref_alignments = {}

    genes_needed = set()
    for _, gene in flagged_list:
        g = gene.upper()
        genes_needed.add(rrna_map.get(g, g))

    for gene_upper in sorted(genes_needed):
        found_path = None
        for candidate in [
            aln_folder / f"{gene_upper}_aln.fasta",
            aln_folder / f"{gene_upper}_aligned.fasta",
            aln_folder / f"{gene_upper}.fasta",
            aln_folder / f"{gene_upper}_aln.fa",
            aln_folder / f"{gene_upper}.fa",
            aln_folder / f"{gene_upper.replace('_RRNA', '_rRNA')}_aln.fasta",
            aln_folder / f"{gene_upper.replace('_RRNA', '_rRNA')}.fasta",
        ]:
            if candidate.exists():
                found_path = candidate
                break

        if found_path is None:
            print(f"  Warning: no alignment for {gene_upper}")
            continue

        seqs = {}
        for record in SeqIO.parse(str(found_path), "fasta"):
            taxon_key = None
            parts     = record.id.split("_")
            for i in range(len(parts) - 1):
                if (parts[i] and parts[i][0].isupper()
                        and parts[i].isalpha()
                        and parts[i+1]
                        and parts[i+1][0].islower()
                        and parts[i+1].isalpha()):
                    taxon_key = f"{parts[i]}_{parts[i+1]}"
                    break
            if taxon_key is None:
                taxon_key = (f"{parts[0]}_{parts[1]}"
                             if len(parts) >= 2 else record.id)
            seqs[taxon_key] = str(record.seq)

        ref_alignments[gene_upper] = seqs
        print(f"  Loaded {gene_upper:<15} "
              f"{len(seqs)} taxa from {found_path.name}")

    # Run verification
    results           = {}
    confident_chimera = []
    probable_chimera  = []
    not_detected      = []
    insufficient      = []

    for taxon, gene in sorted(flagged_list):
        gene_upper = rrna_map.get(gene.upper(), gene.upper())
        ref_aln    = ref_alignments.get(gene_upper, {})

        sub_info = get_subsamples_from_csv(
            results_csv  = results_csv,
            taxon        = taxon,
            gene         = gene,
            label_filter = label_filter,
        )

        if sub_info["n_subsamples"] == 0:
            print(f"  {taxon} | {gene}: no subsamples in CSV")
            results[(taxon, gene)] = {
                "confidence": "INSUFFICIENT_DATA",
                "note":       "No subsamples found in CSV",
            }
            insufficient.append((taxon, gene))
            continue

        # Resolve actual taxon name (may differ if keyed by accession)
        actual_taxon = taxon
        if sub_info["subsamples"]:
            csv_taxon = sub_info["subsamples"][0].get(
                "taxon_in_csv", ""
            ).strip()
            if (csv_taxon and csv_taxon != taxon
                    and "_" in csv_taxon
                    and csv_taxon[0].isupper()
                    and not csv_taxon[0].isdigit()):
                actual_taxon = csv_taxon
                print(f"  Note: '{taxon}' → '{actual_taxon}'")

        sub_seqs = extract_subsample_sequences(
            subsample_info = sub_info,
            fasta_folder   = fasta_folder,
            taxon          = actual_taxon,
            gene           = gene,
        )

        if not sub_seqs:
            print(f"  {taxon} | {gene}: sequences not found "
                  f"in FASTA folder")
            results[(taxon, gene)] = {
                "confidence": "INSUFFICIENT_DATA",
                "note":       "Sequences not found in FASTA folder.",
            }
            insufficient.append((taxon, gene))
            continue

        aligned_subs = (
            [profile_align_subsample(seq, ref_aln) for seq in sub_seqs]
            if ref_aln else sub_seqs
        )

        if not ref_aln:
            print(f"  Warning: no alignment for {gene_upper} — "
                  f"using unaligned distances")

        family_threshold = (
            thresholds.get("per_gene_thresholds", {})
                      .get(gene_upper,
                           thresholds["family_distance_threshold"])
        )
        variance_threshold = thresholds.get(
            "dist_variance_threshold", 0.05
        )

        result = verify_subsamples_by_family_distance(
            subsample_seqs            = aligned_subs,
            query_taxon               = actual_taxon,
            gene                      = gene,
            reference_sequences       = ref_aln,
            taxonomy                  = taxonomy,
            family_distance_threshold = family_threshold,
            subfamily_z_threshold     = variance_threshold,
        )

        results[(taxon, gene)] = result
        confidence = result["confidence"]

        if confidence == "CONFIDENT_CHIMERA":
            confident_chimera.append((taxon, gene))
        elif confidence == "PROBABLE_CHIMERA":
            probable_chimera.append((taxon, gene))
        elif confidence == "NOT_DETECTED":
            not_detected.append((taxon, gene))
        else:
            insufficient.append((taxon, gene))

        print(f"  {taxon:<40} {gene:<12} → {confidence}")
        if result.get("note"):
            print(f"    {result['note']}")

    print(f"\n── Summary ──────────────────────────────────────────────")
    print(f"  CONFIDENT_CHIMERA : {len(confident_chimera)}")
    print(f"  PROBABLE_CHIMERA  : {len(probable_chimera)}")
    print(f"  NOT_DETECTED      : {len(not_detected)}")
    print(f"  INSUFFICIENT_DATA : {len(insufficient)}")

    if confident_chimera:
        print(f"\n  Confirmed chimeras (exclude):")
        for t, g in sorted(confident_chimera):
            print(f"    {t} | {g}")

    if not_detected:
        print(f"\n  Cleared by family distance (retain):")
        for t, g in sorted(not_detected):
            print(f"    {t} | {g}")

    output = {
        "results":           {f"{t}|{g}": v
                              for (t, g), v in results.items()},
        "confident_chimera": [list(x) for x in confident_chimera],
        "probable_chimera":  [list(x) for x in probable_chimera],
        "not_detected":      [list(x) for x in not_detected],
        "insufficient":      [list(x) for x in insufficient],
    }

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(output, f, indent=2, default=str)
        print(f"\n  Results saved: {output_path}")

    return output


def build_final_dataset(
        blast_result: dict,
        verification: dict,
        qual: dict,
        gbk_index: dict,
        gbk_folder: str,
        output_path: str = None,
) -> dict:
    """
    Assemble the final dataset from BLAST and family distance
    verification results, with PubMed qualification for all sequences.

    Logic
    -----
    Cleared by BLAST + PubMed qualified  → include (high confidence)
    Cleared by BLAST + no PubMed         → include with flag
    Flagged by BLAST + NOT_DETECTED      → include (rescued by distance)
    Flagged by BLAST + PROBABLE_CHIMERA  → exclude conservatively
    Flagged by BLAST + CONFIDENT_CHIMERA → exclude
    Flagged by BLAST + INSUFFICIENT_DATA → include if PubMed, else exclude
      (handles orphan taxa like Solenodon that cannot be distance-tested)

    PubMed information is recorded for ALL sequences in the dataset
    regardless of how they were qualified, providing an additional
    quality indicator for downstream use.

    Parameters
    ----------
    blast_result   : dict   from score_from_csv()
    verification   : dict   from run_family_distance_verification()
    qual           : dict   from qualify_cleared_sequences()
    gbk_index      : dict   from build_gbk_index()
    gbk_folder     : str    directory containing .gbk files
    output_path    : str    optional path to save final dataset CSV

    Returns
    -------
    dict with:
        dataset        : list   dicts with taxon, gene, source,
                                pubmed, confidence fields
        included       : list   (taxon, gene) pairs in final dataset
        excluded       : list   (taxon, gene) pairs excluded
        pubmed_status  : dict   (taxon, gene) -> bool (has PubMed)
        summary        : dict   counts by source
    """
    from collections import defaultdict

    def check_pubmed(taxon: str) -> bool:
        """Check if any GBK file for this taxon has a PubMed ref."""
        # taxon may actually be an accession string — try it directly first
        candidates = []
        
        # Direct accession lookup (taxon IS the accession)
        direct = Path(gbk_folder) / f"{taxon}.gbk"
        if not direct.exists():
            direct = Path(gbk_folder) / f"{taxon}.gb"
        if direct.exists():
            candidates.append(direct)
        
        # Also try via gbk_index (for taxon-name keyed entries)
        for acc in gbk_index.get(taxon, []):
            for ext in ("gbk", "gb"):
                p = Path(gbk_folder) / f"{acc}.{ext}"
                if p.exists():
                    candidates.append(p)
        
        for gbk_path in candidates:
            try:
                for record in SeqIO.parse(str(gbk_path), "genbank"):
                    for ref in record.annotations.get("references", []):
                        if getattr(ref, "pubmed_id", "").strip():
                            return True
            except Exception:
                continue
        return False

    print(f"\n{'='*60}")
    print(f"BUILD FINAL DATASET")
    print(f"{'='*60}")

    # Index verification results for fast lookup
    not_detected      = set(
        tuple(x) for x in verification.get("not_detected", [])
    )
    confident_chimera = set(
        tuple(x) for x in verification.get("confident_chimera", [])
    )
    probable_chimera  = set(
        tuple(x) for x in verification.get("probable_chimera", [])
    )
    insufficient      = set(
        tuple(x) for x in verification.get("insufficient", [])
    )

    all_flagged = set(tuple(x) for x in blast_result.get("flagged", []))
    all_cleared = set(
        (taxon, gene)
        for (taxon, gene), call in blast_result["chimera_calls"].items()
        if not call["chimera"]
    )

    # PubMed cache — check each taxon once
    pubmed_cache  = {}
    pubmed_status = {}

    def get_pubmed(taxon):
        if taxon not in pubmed_cache:
            pubmed_cache[taxon] = check_pubmed(taxon)
        return pubmed_cache[taxon]

    dataset  = []
    included = []
    excluded = []
    summary  = defaultdict(int)

    # Process cleared sequences
    for taxon, gene in sorted(all_cleared):
        has_pubmed = get_pubmed(taxon)
        pubmed_status[(taxon, gene)] = has_pubmed
        is_qualified = gene in qual.get("qualified", {}).get(taxon, [])

        source     = "blast_cleared_pubmed" if has_pubmed \
                     else "blast_cleared"
        confidence = "HIGH" if has_pubmed else "MEDIUM"

        dataset.append({
            "taxon":      taxon,
            "gene":       gene,
            "included":   True,
            "source":     source,
            "confidence": confidence,
            "pubmed":     has_pubmed,
            "qualified":  is_qualified,
        })
        included.append((taxon, gene))
        summary[source] += 1

    # Process flagged sequences
    for taxon, gene in sorted(all_flagged):
        has_pubmed = get_pubmed(taxon)
        pubmed_status[(taxon, gene)] = has_pubmed
        pair = (taxon, gene)

        if pair in not_detected:
            # Family distance says genuine — retain
            source     = "flagged_blast_cleared_distance"
            confidence = "MEDIUM"
            include    = True

        elif pair in confident_chimera:
            # Both BLAST and distance agree — exclude
            source     = "confirmed_chimera"
            confidence = "HIGH_CHIMERA"
            include    = False

        elif pair in probable_chimera:
            # Probable chimera — exclude conservatively
            source     = "probable_chimera"
            confidence = "PROBABLE_CHIMERA"
            include    = False

        elif pair in insufficient:
            # Could not be distance-tested (orphan taxon etc.)
            # Retain if PubMed-backed, exclude otherwise
            if has_pubmed:
                source     = "flagged_blast_pubmed_retained"
                confidence = "LOW"
                include    = True
            else:
                source     = "flagged_blast_no_evidence"
                confidence = "EXCLUDE"
                include    = False

        else:
            # Not in verification results at all — exclude conservatively
            source     = "flagged_blast_unverified"
            confidence = "EXCLUDE"
            include    = False

        dataset.append({
            "taxon":      taxon,
            "gene":       gene,
            "included":   include,
            "source":     source,
            "confidence": confidence,
            "pubmed":     has_pubmed,
            "qualified":  False,
        })

        if include:
            included.append(pair)
            summary[source] += 1
        else:
            excluded.append(pair)
            summary[source] += 1

    # Print summary
    print(f"\n  Total sequences evaluated : "
          f"{len(all_cleared) + len(all_flagged)}")
    print(f"  Included in final dataset : {len(included)}")
    print(f"  Excluded                  : {len(excluded)}")
    print(f"\n  By source:")
    for source, count in sorted(summary.items()):
        print(f"    {source:<45} {count}")

    pubmed_included = sum(1 for t, g in included
                          if pubmed_status.get((t, g), False))
    print(f"\n  PubMed-backed in final dataset: "
          f"{pubmed_included}/{len(included)} "
          f"({100*pubmed_included/len(included):.1f}%)"
          if included else "")

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fields = ["taxon", "gene", "included", "source",
                  "confidence", "pubmed", "qualified"]
        with open(output_path, "w", newline="",
                  encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(dataset)
        print(f"\n  Dataset saved: {output_path}")

    sample = list(blast_result["chimera_calls"].keys())[:5]
    for key in sample:
        taxon = key[0] if isinstance(key, tuple) else key
        p = Path(gbk_folder) / f"{taxon}.gbk"
        print(f"{taxon} → exists: {p.exists()}")
    return {
        "dataset":       dataset,
        "included":      included,
        "excluded":      excluded,
        "pubmed_status": pubmed_status,
        "summary":       dict(summary),
    }




from bioinf_packages.verify_funcs.batch_blast_verify import score_from_csv
from bioinf_packages.alignment_funcs._extract_accessions import extract_accessions
from bioinf_packages.verify_funcs._qualify_cleared import (extract_verification_lists, qualify_cleared_sequences, build_gbk_index, assess_flagged_sequences)
from bioinf_packages.verify_funcs._species_parser import load_taxonomy

results = score_from_csv(
    results_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/blast_results.csv",
    calibration_json = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/model_calibration.json",
    output_path = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/chimera_calls.csv",
    label_filter = "GENUINE",
)
taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")
lists = extract_verification_lists(result=results)

calibration_set = extract_accessions("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/callibration_accessions.csv")
print(calibration_set)
index = build_gbk_index(
    gbk_folder="C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    glossary ="C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv",
)

qual = qualify_cleared_sequences(
    cleared_by_taxon       = lists["cleared_by_taxon"],
    gbk_folder             = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    calibration_accessions = calibration_set,
    gbk_index = index,
    glossary = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"
)

assessment = assess_flagged_sequences(
    result     = results,
    qual       = qual,
    gbk_folder = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    glossary   = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv",
    gbk_index  = index,
)

# Step 1: combine all FASTAs into one master file  
#combined = build_combined_fasta(
#    fasta_paths = [
#        "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/all_queries.fasta",
#        "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta",
#    ],
#    output_path = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/all_sequences_combined.fasta",
#)

from collections import defaultdict

# Build a broader qualified set for calibration purposes only
# — includes all cleared sequences plus calibration set,
# without requiring PubMed validation
calibration_qualified = {}

# Add all cleared sequences
for taxon, genes in lists["cleared_by_taxon"].items():
    calibration_qualified[taxon] = list(genes)

# Step 2: add calibration sequences to qualified dict
# (they were never scored by score_from_csv so are not in qual)
qualified_with_cal = add_calibration_to_qualified(
    qualified          = calibration_qualified,
    calibration_set    = calibration_set,
    gbk_index          = index,
    calibration_fasta  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta",
)
print(f"Calibration qualified taxa: {len(calibration_qualified)}")

# Check family distribution
family_counts = defaultdict(list)
for taxon in calibration_qualified:
    genus  = taxon.split("_")[0]
    family = taxonomy.get(genus, {}).get("Family", "Unknown")
    family_counts[family].append(taxon)

print("Taxa per family:")
for fam, taxa in sorted(family_counts.items()):
    print(f"  {fam:<30} {len(taxa)} taxa")

from collections import defaultdict

# Diagnostic — print how many taxa per family
family_counts = defaultdict(list)
for taxon in qualified_with_cal:
    genus  = taxon.split("_")[0]
    family = taxonomy.get(genus, {}).get("Family", "Unknown")
    family_counts[family].append(taxon)

print("Taxa per family in qualified_with_cal:")
for fam, taxa in sorted(family_counts.items()):
    print(f"  {fam:<30} {len(taxa)} taxa: "
          f"{', '.join(sorted(taxa))}")

print(taxonomy.get("Erinaceus", {}))
print(taxonomy.get("Sorex", {}))
print(taxonomy.get("Talpa", {}))

sample = next(iter(taxonomy.values()))
print("Taxonomy columns:", list(sample.keys()))
    
# Step 3: build per-gene reference FASTAs

# Resolve accession-keyed cleared_by_taxon to taxon names
#cleared_by_taxon_resolved = resolve_accession_keys(
#    cleared_by_taxon = lists["cleared_by_taxon"],
#    results_csv      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/blast_results.csv",
#)
#fasta_info = build_reference_fastas(
#    cleared_by_taxon  = cleared_by_taxon_resolved,  # not qualified
#    combined_fasta    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/all_sequences_combined.fasta",
#    output_dir        = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/fasta_data",
#    calibration_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta",
#    min_taxa_per_gene = 3,
#)

# Step 4: run ClustalOmega on each file in fasta_data/
# Step 4.5: split short_aln and aln into separate folders.

# ── Step 5: Calibrate thresholds from qualified sequences ─────────────────────
# Uses the aligned reference FASTAs you just built.
# alignment_folder should contain [gene]_aln.fasta files.
# The function looks for *.fasta files so it will find them.
# However the current code strips "_aligned" from the stem to get
# the gene name — "_aln" needs to be handled the same way.


print("Calling calibrate_family_distance_thresholds...")
thresholds_long = calibrate_family_distance_thresholds(
    qualified        = calibration_qualified,   # still needed for the print
    alignment_folder = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/nexus_data/long",
    taxonomy         = taxonomy,
    min_taxa_per_family = 3,
    percentile_threshold = 95.0,
    output_path      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/fam_dist_thresh_long.json",
)


#thresholds_short = calibrate_family_distance_thresholds(
#    qualified         = calibration_qualified,
#    alignment_folder  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/nexus_data/short",
#    taxonomy          = taxonomy,
#    min_taxa_per_family = 3,
#    percentile_threshold = 95.0,
#    output_path       = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/fam_dist_thresh_short.json",
#)

# ── Step 6: Run verification on REVIEW sequences ──────────────────────────────
verification = run_family_distance_verification(
    flagged_list     = list(results["flagged"]),
    results_csv      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/blast_results.csv",
    fasta_folder     = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data",   # original per-gene FASTAs
    alignment_folder = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/nexus_data/long",
    taxonomy         = taxonomy,
    thresholds       = thresholds_long,
    label_filter     = "GENUINE",
    output_path      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/verify_family/verification_results.json",
)


# 3 — Build final dataset
final = build_final_dataset(
    blast_result = results,
    verification = verification,
    qual         = qual,
    gbk_index    = index,
    gbk_folder   = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    output_path  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/final_dataset.csv",
)


"""
Summary looks reasonable:

51 CONFIDENT_CHIMERA, 56 PROBABLE_CHIMERA, 16 INSUFFICIENT_DATA, 0 NOT_DETECTED
Solenodon sequences (AY451972, AY530070, AY530075, AY530080, JN414026, JN414741, KU697358, LC124832, LC124888, LC124921, LC124955) all correctly returning INSUFFICIENT_DATA as expected — good to manually retain those
Erinaceidae BDNF references were also 0 (AY986746/Echinosorex_gymnura, AY986748/Paraechinus_aethiopicus, JN633319/Podogymnura_truei, JN633375/Podogymnura_truei, OR554407/Hylomys_suillus) — this makes sense since Erinaceidae BDNF wasn't represented in your gold set or cleared sequences enough to build references for that gene

One thing to verify — NC_042734 (Blarina_brevicauda) is being verified as a mitogenome with many flagged genes. Looking at the results, all its CONFIDENT hits show min_dist > 0.1649 for COX1, which is the threshold. The subsamples that fail are borderline (min_dist 0.175-0.190), which is genuinely ambiguous for Blarina — it's a highly divergent Soricinae. Worth keeping in mind when interpreting.
One biological flag — NC_002808 (Echinosorex_gymnura) is flagged CONFIDENT across 12 genes with very high prop_out_of_family values (up to 0.8 for ND2). That's an extremely strong chimera signal across essentially the entire mitogenome, suggesting this is a genuine assembly-level problem or misidentification in GenBank.
Everything looks mechanically correct. The INSUFFICIENT_DATA cases are all genuinely data-limited taxa (Solenodontidae, sparse Erinaceidae genes), not errors. You're good to run build_final_dataset().
"""

"""
Final list removed after manual review (relaxed):
KX754645 - Dymecodon BDNF
MK410422 - Dymecodon COX1
AY986746 - Echinosorex BDNF
MG973451 - Episoriculus RAG1
EU122210 - Euroscaptor 16S rRNA
OR554407 - Hylomys BDNF
AY121754 - Hylomys BRCA1
KF783055 - Neotetracus RAG1
KX754649 - Neurotrichus BDNF
MZ150484 - Neurotrichus COX1
LC124967 - Notiosorex RAG1
MN061469 - Notiosorex COX1
AY986748 - Paraechinus BDNF
JN414025 - Podogymnura APOB
AF434829 - Podogymnura ND2
AY170059 - Scalopus ND2
HM902768 - Scalopus COX1
AF284007 - Scalopus BRCA1
DQ630353 - Notiosorex 16S rRNA
AF069539 - Scalopus 12S rRNA
EF027282 - Sorex alpinus 12S rRNA

Problem Taxa Still Included:
Echinosorex gymnura (mitogenome)
Mogera wogura (mitogenome)
Myosorex kihaulei (mitogenome)

Final list removed after manual review (harsh):
GU981136 - Anourosorex yamashinai ATP6
OM220058 - Anourosorex yamashinai COX1
GU981349 - Anourosorex yamashinai ND4
GU981395 - Anourosorex yamashinai ND5
MG973425 - Congosorex CYTB
GU473580 - Congosorex ND2
AY121755 - Crocidura russula BRCA1
LC124966 - Cryptotis RAG1
OP855695 - Myosorex COX1
OP853578 - Myosorex ND4
OP853565 - Mysorex ND5
AY691822 - Myosorex 12S rRNA
MZ217184 - Neurotrichus COX3
AF434834 - Parascalops ND2
MZ217186 - Scalopus COX3
AF250464 - Sorex monticolus ND4
KC113264 - Sorex araneus RAG1
NC_002808 - Echinosorex mitogenome
KX754645 - Dymecodon BDNF
MK410422 - Dymecodon COX1
AY986746 - Echinosorex BDNF
MG973451 - Episoriculus RAG1
EU122210 - Euroscaptor 16S rRNA
OR554407 - Hylomys BDNF
AY121754 - Hylomys BRCA1
KF783055 - Neotetracus RAG1
KX754649 - Neurotrichus BDNF
MZ150484 - Neurotrichus COX1
LC124967 - Notiosorex RAG1
MN061469 - Notiosorex COX1
AY986748 - Paraechinus BDNF
JN414025 - Podogymnura APOB
AF434829 - Podogymnura ND2
AY170059 - Scalopus ND2
HM902768 - Scalopus COX1
AF284007 - Scalopus BRCA1
DQ630353 - Notiosorex 16S rRNA
AF069539 - Scalopus 12S rRNA
EF027282 - Sorex alpinus 12S rRNA

Problem taxa still included:
Echinosorex gymnura (mitogenome)
Mogera wogura (mitogenome)

"""