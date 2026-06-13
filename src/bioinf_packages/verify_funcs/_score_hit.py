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

def score_blast(query_accession: str,
                seq: str,
                taxon: str,
                gene: str,
                taxonomy: dict,
                glossary: dict,
                prior_prob: float,
                pident_threshold: float,
                k: float,
                beta: float,
                n: int = None,
                insert_fraction: float = 0.30,
                target_detection_prob: float = 0.90,
                blast_db: str = None,
                blast_bin: str = "blastn",
                num_threads: int = 4,
                evalue: float = 1e-10,
                max_target_seqs: int = 20,
                perc_identity: float = 70.0,
                output_dir: str = None,
                flag_on_single_hit: bool = True,
                single_hit_taxon_llr_threshold: float = -2.0) -> dict:
    """
    Run BLAST across subsamples of a query sequence, score all hits,
    and return an aggregated posterior probability of genuineness.

    With flag_on_single_hit=True (default), a single subsample whose
    top hit has a taxon LLR below single_hit_taxon_llr_threshold flags
    the entire sequence as a chimera candidate regardless of the
    aggregate posterior. This matches the subsampling design which gives
    90% probability that at least one subsample overlaps a chimeric
    junction — so a single clearly wrong taxonomic hit is meaningful
    signal.

    Supports both local BLAST (blast_db + blast_bin) and the legacy
    verify_seq() remote BLAST path. If blast_db is provided, local
    BLAST is used. Otherwise falls back to verify_seq().

    Parameters
    ----------
    query_accession                : str   accession for self-hit removal
    seq                            : str   full query sequence
    taxon                          : str   e.g. "Sorex_araneus"
    gene                           : str   e.g. "CYTB" or "CYTB_modified_1"
    taxonomy                       : dict  from load_taxonomy()
    glossary                       : dict  from load_glossary()
    prior_prob                     : float prior P(genuine)
    pident_threshold               : float pident contributing LLR of 0
    k                              : float steepness of pident LLR mapping
    beta                           : float softmax temperature over bitscores
    n                              : int   subsample length. If None,
                                           calculated automatically.
    insert_fraction                : float expected insert size as fraction
                                           of sequence length (default 0.30)
    target_detection_prob          : float desired P(detect chimera)
                                           (default 0.90)
    blast_db                       : str   path to local BLAST database.
                                           If None, uses remote verify_seq().
    blast_bin                      : str   path to blastn executable
    num_threads                    : int   BLAST CPU threads
    evalue                         : float BLAST e-value threshold
    max_target_seqs                : int   max hits per subsample
    perc_identity                  : float minimum percent identity
    output_dir                     : str   directory for BLAST output files.
                                           Required when blast_db is set.
    flag_on_single_hit             : bool  if True, any subsample with
                                           taxon_llr below threshold flags
                                           whole sequence (default True)
    single_hit_taxon_llr_threshold : float taxon LLR below which a single
                                           subsample triggers chimera flag.
                                           Default -2.0.

    Returns
    -------
    dict with:
        query_accession         : str
        taxon                   : str
        gene                    : str
        seq_length              : int
        subsample_length        : int
        n_samples               : int
        sample_results          : list  per-sample aggregated scores
        sample_posteriors       : list  per-sample posterior probabilities
        posterior_variance      : float variance across sample posteriors
        taxon_llr_variance      : float variance of taxon LLRs across hits
        prop_neg_taxon_hits     : float fraction of hits with negative LLR
        chimera_flag            : bool  True if aggregate OR single-hit
                                        flag triggered
        single_hit_flag         : bool  True if any single subsample
                                        triggered the flag
        n_flagged_samples       : int   number of subsamples flagged
        flagged_sample_indices  : list  which sample numbers were flagged
        flagged_sample_details  : list  dicts with position and top hit
                                        for each flagged subsample
        cumulative_log_odds     : float
        final_posterior         : float
        prior_prob              : float
        prior_log_odds          : float
        total_llr               : float
        posterior_log_odds      : float
        posterior_prob          : float
        hit_details             : list  all hit dicts across all samples
    """
    import tempfile

    # ── Strip _modified_N suffix for gene parser ──────────────────────────────
    clean_gene = re.sub(r'_modified(_\d+)?$', '', gene)

    # ── Auto-scale subsample parameters ──────────────────────────────────────
    if n is None:
        n = recommended_subsample_length(
            len(seq), insert_fraction=insert_fraction
        )

    n_samples = recommended_n_samples(
        seq_length            = len(seq),
        subsample_length      = n,
        target_detection_prob = target_detection_prob,
        insert_fraction       = insert_fraction,
    )

    print(f"  Sequence: {len(seq)}bp | "
          f"Subsample: {n}bp | Samples: {n_samples}")

    # ── Run BLAST — local or remote ───────────────────────────────────────────
    if blast_db is not None:
        # Local BLAST via blastn subprocess
        if output_dir is None:
            raise ValueError(
                "output_dir must be provided when using local BLAST "
                "(blast_db is set)."
            )
        os.makedirs(output_dir, exist_ok=True)

        # Draw subsamples and write to FASTA
        import random
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq as BioSeq

        subsamples  = []
        sample_meta = []
        for i in range(1, n_samples + 1):
            start  = random.randint(0, len(seq) - n)
            sample = seq[start:start + n]
            header = (f"{taxon}|{clean_gene}|accession:{query_accession}"
                      f"|Sample{i}")
            subsamples.append(
                SeqRecord(BioSeq(sample), id=header, description="")
            )
            sample_meta.append({
                "sample_index": i - 1,
                "sample":       i,
                "query_start":  start,
                "query_end":    start + n,
            })

        query_fasta = os.path.join(
            output_dir,
            f"{taxon}_{query_accession}_subsamples.fasta"
        )
        SeqIO.write(subsamples, query_fasta, "fasta")

        blast_out = os.path.join(
            output_dir,
            f"{taxon}_{query_accession}_blast.xml"
        )
        run_blast_local(
            query_fasta     = query_fasta,
            db              = blast_db,
            out_file        = blast_out,
            blast_bin       = blast_bin,
            evalue          = evalue,
            max_target_seqs = max_target_seqs,
            num_threads     = num_threads,
            perc_identity   = perc_identity,
        )

        # Parse XML results into the same format as verify_seq()
        blast_result = _parse_blast_xml_to_samples(
            blast_xml    = blast_out,
            sample_meta  = sample_meta,
            accession    = query_accession,
        )

    else:
        # Legacy remote BLAST via verify_seq()
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
        sample_agg["query_start"]  = sample.get("query_start", None)
        sample_agg["query_end"]    = sample.get("query_end", None)
        sample_agg["sample"]       = sample["sample_index"] + 1

        # Store top hit taxon LLR for single-hit flagging
        if scored_hits:
            top_hit = max(scored_hits, key=lambda h: h["bitscore"])
            sample_agg["top_taxon_llr"]   = top_hit["species_llr"]
            sample_agg["top_taxon"]        = (
                top_hit["species_result"].get("matched_taxon", "unknown")
            )
            sample_agg["top_gene_llr"]    = top_hit["gene_llr"]
            sample_agg["top_accession"]   = top_hit["accession"]
        else:
            sample_agg["top_taxon_llr"]  = None
            sample_agg["top_taxon"]       = None
            sample_agg["top_gene_llr"]   = None
            sample_agg["top_accession"]  = None

        sample_posteriors.append(sample_agg["posterior_prob"])
        sample_results.append(sample_agg)

    if not sample_results:
        raise ValueError(
            f"No samples returned any BLAST hits for "
            f"{taxon} | {gene} | {query_accession}."
        )

    # ── Aggregate log-odds across all samples ─────────────────────────────────
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

    posterior_variance = (float(np.var(sample_posteriors))
                          if len(sample_posteriors) > 1 else 0.0)
    taxon_llr_variance = (float(np.var(all_taxon_llrs))
                          if len(all_taxon_llrs) > 1 else 0.0)
    prop_neg_taxon     = (sum(1 for x in all_taxon_llrs if x < 0) /
                          len(all_taxon_llrs)
                          if all_taxon_llrs else 0.0)

    # Aggregate chimera flag — existing variance-based logic
    aggregate_chimera_flag = (
        posterior_variance > 0.05
        or taxon_llr_variance > 2.0
    )

    # ── Single-hit flagging ───────────────────────────────────────────────────
    # Flag any subsample whose top hit has a taxon LLR below the threshold.
    # This catches chimeras where only 1-2 subsamples land in the
    # contaminated region — the aggregate posterior would look fine but
    # the individual subsample clearly hits the wrong taxon.
    flagged_samples        = []
    flagged_sample_details = []

    if flag_on_single_hit:
        for s in sample_results:
            top_llr = s.get("top_taxon_llr")
            if top_llr is not None and top_llr < single_hit_taxon_llr_threshold:
                flagged_samples.append(s["sample"])
                flagged_sample_details.append({
                    "sample":        s["sample"],
                    "query_start":   s.get("query_start"),
                    "query_end":     s.get("query_end"),
                    "top_taxon_llr": round(top_llr, 4),
                    "top_taxon":     s.get("top_taxon", "unknown"),
                    "top_gene_llr":  s.get("top_gene_llr"),
                    "top_accession": s.get("top_accession"),
                })
                print(f"  ⚠ Sample {s['sample']} flagged: "
                      f"taxon_llr={top_llr:.3f} < "
                      f"{single_hit_taxon_llr_threshold} | "
                      f"top_taxon={s.get('top_taxon', 'N/A')} | "
                      f"pos={s.get('query_start','?')}-"
                      f"{s.get('query_end','?')}")

    single_hit_flag = len(flagged_samples) > 0
    chimera_flag    = aggregate_chimera_flag or single_hit_flag

    # ── Print summary ─────────────────────────────────────────────────────────
    flag_str = "⚠ CHIMERA" if chimera_flag else "✓ OK"
    print(f"  {flag_str} | "
          f"posterior={final_posterior:.4f} | "
          f"post_var={posterior_variance:.4f} | "
          f"taxon_llr_var={taxon_llr_variance:.4f} | "
          f"prop_neg={prop_neg_taxon:.3f} | "
          f"flagged_samples={len(flagged_samples)}/{len(sample_results)}")

    return {
        "query_accession":              query_accession,
        "taxon":                        taxon,
        "gene":                         gene,
        "seq_length":                   len(seq),
        "subsample_length":             n,
        "n_samples":                    n_samples,
        "sample_results":               sample_results,
        "sample_posteriors":            sample_posteriors,
        "posterior_variance":           round(posterior_variance, 6),
        "taxon_llr_variance":           round(taxon_llr_variance, 4),
        "prop_neg_taxon_hits":          round(prop_neg_taxon, 4),
        "aggregate_chimera_flag":       aggregate_chimera_flag,
        "single_hit_flag":              single_hit_flag,
        "chimera_flag":                 chimera_flag,
        "n_flagged_samples":            len(flagged_samples),
        "flagged_sample_indices":       flagged_samples,
        "flagged_sample_details":       flagged_sample_details,
        "flag_on_single_hit":           flag_on_single_hit,
        "single_hit_taxon_llr_threshold": single_hit_taxon_llr_threshold,
        "cumulative_log_odds":          round(cumulative_log_odds, 4),
        "final_posterior":              round(final_posterior, 4),
        "prior_prob":                   prior_prob,
        "prior_log_odds":               round(prior_log_odds, 4),
        "total_llr":                    round(
                                          cumulative_log_odds -
                                          prior_log_odds, 4),
        "posterior_log_odds":           round(cumulative_log_odds, 4),
        "posterior_prob":               round(final_posterior, 4),
        "hit_details":                  [h for s in sample_results
                                         for h in s["hit_details"]],
    }


def _parse_blast_xml_to_samples(blast_xml: str,
                                 sample_meta: list,
                                 accession: str) -> dict:
    """
    Parse a BLAST XML output file into the same dict format that
    verify_seq() returns, so score_blast() can handle both local
    and remote BLAST with identical downstream code.

    Parameters
    ----------
    blast_xml   : str   path to BLAST XML output file
    sample_meta : list  dicts with sample_index, query_start, query_end
                        for each subsample — in the same order they were
                        written to the query FASTA
    accession   : str   query accession for self-hit removal

    Returns
    -------
    dict with "samples" key matching verify_seq() output format
    """
    from Bio.Blast import NCBIXML

    samples = []

    with open(blast_xml, "r") as f:
        blast_records = list(NCBIXML.parse(f))

    for idx, (record, meta) in enumerate(
        zip(blast_records, sample_meta)
    ):
        hits = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:

                # Remove self-hits by accession
                hit_acc = alignment.accession
                if accession and accession in hit_acc:
                    continue

                hits.append({
                    "accession":    hit_acc,
                    "description":  alignment.title,
                    "bitscore":     hsp.bits,
                    "pident":       (hsp.identities /
                                     hsp.align_length * 100
                                     if hsp.align_length > 0 else 0.0),
                    "evalue":       hsp.expect,
                    "align_length": hsp.align_length,
                })

        samples.append({
            "sample_index": meta["sample_index"],
            "query_start":  meta["query_start"],
            "query_end":    meta["query_end"],
            "hits":         hits,
        })

    return {"samples": samples}
#print(score_blast("GU981106", "CACTTCCTTTGGATATGYTTGATGTGTTTTTGAATCATAATATCAATTCCTTTCTGAGGCAAGTTGAGAAGGTCAGAGATGAGGCATTGGTTCTTGTTATTCAATCCTATAATGAAGCAAAAATGAAATTTGATGAGCATAAGGTTGAAAAATCTATCACCCAACAACGAAAGACCTTTCAAATTCCAGGGTACACCATTCCTGTTGTTAATGTCGAAGTGTCTCCATTCACAGTAGAGATGTTTCCATTTGGTTATGTGATCCCAAAGGAGGTCAGCACCCCAAAGTTCACCATCCTGGGTTCTGGTTTCTCTGTGCCTTCCTATACTTTAGTCCTGCCCTTTCTAGAACTACCAGCTCTTCATATCCCTAAGTTTCTTGAGCTTTCTTTTCCAGACTTCAAAGTATCGAGTATCCCAAGGAATATTTTCATTCCAGCCCTGGGAAATGTTACATATGATTTTTCCTTTAAGTCAAGTGTCATTACACTGAATGCCAATGCTGGACTTTAT", 50,"Anourosorex_squamipes", "APOB", load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"), load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"), 0.95, 97.0, 0.3, 0.1))