#"c:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_verify"


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

def plot_blast_result(score_blast_result: dict, sequence_label: str, save_dir: str):
    """
    Visualise the output of score_blast() for a single sequence.

    Parameters
    ----------
    score_blast_result : dict   output of score_blast()
    sequence_label     : str    e.g. "Anourosorex_squamipes_seq1"
    """
    hit_details  = score_blast_result["hit_details"]
    valid_hits   = [h for h in hit_details if not h.get("gated_out", False)]
    posterior    = score_blast_result["posterior_prob"]
    total_llr    = score_blast_result["total_llr"]

    ranks        = list(range(1, len(valid_hits) + 1))
    pidents      = [h["pident"]       for h in valid_hits]
    bitscores    = [h["bitscore"]     for h in valid_hits]
    llr_taxon    = [h["llr_taxon"]    for h in valid_hits]
    llr_gene     = [h["llr_gene"]     for h in valid_hits]
    llr_pident   = [h["llr_pident"]   for h in valid_hits]
    combined_llr = [h["combined_llr"] for h in valid_hits]
    weights      = [h["weight"]       for h in valid_hits]

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(f"{sequence_label}  |  P(genuine) = {posterior:.3f}  |  total LLR = {total_llr:.3f}",
                 fontsize=13, fontweight="bold")

    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    # ── Panel 1: pident by hit rank ───────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.bar(ranks, pidents, color="steelblue")
    ax1.axhline(97, color="red", linestyle="--", linewidth=1, label="97% threshold")
    ax1.set_xlabel("Hit rank")
    ax1.set_ylabel("% Identity")
    ax1.set_title("Sequence identity per hit")
    ax1.set_xticks(ranks)
    ax1.legend(fontsize=8)

    # ── Panel 2: bitscore by hit rank ─────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.bar(ranks, bitscores, color="darkorange")
    ax2.set_xlabel("Hit rank")
    ax2.set_ylabel("Bitscore")
    ax2.set_title("Bitscore per hit")
    ax2.set_xticks(ranks)

    # ── Panel 3: softmax weights ───────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.bar(ranks, weights, color="mediumpurple")
    ax3.set_xlabel("Hit rank")
    ax3.set_ylabel("Softmax weight")
    ax3.set_title("Hit weights")
    ax3.set_xticks(ranks)

    # ── Panel 4: LLR decomposition per hit (stacked bar) ─────────────────────
    ax4 = fig.add_subplot(gs[1, :2])
    bar_width = 0.25
    x = np.array(ranks)
    ax4.bar(x - bar_width, llr_taxon,  width=bar_width, label="Taxon LLR",  color="steelblue")
    ax4.bar(x,             llr_gene,   width=bar_width, label="Gene LLR",   color="darkorange")
    ax4.bar(x + bar_width, llr_pident, width=bar_width, label="Pident LLR", color="seagreen")
    ax4.axhline(0, color="black", linewidth=0.8)
    ax4.set_xlabel("Hit rank")
    ax4.set_ylabel("LLR contribution")
    ax4.set_title("Evidence decomposition per hit")
    ax4.set_xticks(ranks)
    ax4.legend(fontsize=8)

    # ── Panel 5: weighted combined LLR per hit ────────────────────────────────
    ax5 = fig.add_subplot(gs[1, 2])
    colours = ["seagreen" if v >= 0 else "firebrick" for v in combined_llr]
    ax5.bar(ranks, combined_llr, color=colours)
    ax5.axhline(0, color="black", linewidth=0.8)
    ax5.set_xlabel("Hit rank")
    ax5.set_ylabel("Weighted LLR")
    ax5.set_title("Combined LLR per hit\n(weight × sum of signals)")
    ax5.set_xticks(ranks)

    plt.savefig(os.path.join(save_dir, f"{sequence_label}_blast_score.png"), 
                dpi=150, bbox_inches="tight")
    plt.show()

plot_blast_result({'prior_prob': 0.95, 'prior_log_odds': 2.9444, 'total_llr': 5.095, 'posterior_log_odds': 8.0394, 'posterior_prob': 0.9997, 'hit_details': [{'accession': 'GU981107', 'weight': 0.1573, 'bitscore': 91.4547, 'pident': 100.0, 'align_length': 50, 'llr_taxon': 2.0, 'llr_gene': 2.5, 'llr_pident': 0.9, 'combined_llr': 0.8496}, {'accession': 'DQ630183', 'weight': 0.1573, 'bitscore': 91.4547, 'pident': 100.0, 'align_length': 50, 'llr_taxon': 2.0, 'llr_gene': 2.5, 'llr_pident': 0.9, 'combined_llr': 0.8496}, {'accession': 'KP995386', 'weight': 0.1573, 'bitscore': 91.4547, 'pident': 100.0, 'align_length': 50, 'llr_taxon': 3.0, 'llr_gene': 2.5, 'llr_pident': 0.9, 'combined_llr': 1.0069}, {'accession': 'GU981106', 'weight': 0.1573, 'bitscore': 91.4547, 'pident': 100.0, 'align_length': 50, 'llr_taxon': 3.0, 'llr_gene': 2.5, 'llr_pident': 0.9, 'combined_llr': 1.0069}, {'accession': 'LC124839', 'weight': 0.1573, 'bitscore': 91.4547, 'pident': 100.0, 'align_length': 50, 'llr_taxon': 2.0, 'llr_gene': 2.5, 'llr_pident': 0.9, 'combined_llr': 0.8496}, {'accession': 'DQ630123', 'weight': 0.1097, 'bitscore': 87.848, 'pident': 98.0, 'align_length': 50, 'llr_taxon': 3.0, 'llr_gene': 2.5, 'llr_pident': 0.3, 'combined_llr': 0.6362}, {'accession': 'MW670247', 'weight': 0.0259, 'bitscore': 73.4211, 'pident': 92.0, 'align_length': 50, 'llr_taxon': -2.0, 'llr_gene': 2.5, 'llr_pident': -1.5, 'combined_llr': -0.0259}, {'accession': 'MW670278', 'weight': 0.0259, 'bitscore': 73.4211, 'pident': 92.0, 'align_length': 50, 'llr_taxon': -2.0, 'llr_gene': 2.5, 'llr_pident': -1.5, 'combined_llr': -0.0259}, {'accession': 'MW670250', 'weight': 0.0259, 'bitscore': 73.4211, 'pident': 92.0, 'align_length': 50, 'llr_taxon': -2.0, 'llr_gene': 2.5, 'llr_pident': -1.5, 'combined_llr': -0.0259}, {'accession': 'MW670299', 'weight': 0.0259, 'bitscore': 73.4211, 'pident': 92.0, 'align_length': 50, 'llr_taxon': -2.0, 'llr_gene': 2.5, 'llr_pident': -1.5, 'combined_llr': -0.0259}]}, "anoSqu_APOB", "c:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_verify")