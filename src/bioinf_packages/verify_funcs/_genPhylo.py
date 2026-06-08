from _phylo_verify import build_scaffold_tree
from pathlib import Path
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json


def plot_scaffold_tree(treefile: str,
                       id_map_path: str = None,
                       save_path: str = None):
    """
    Plot a scaffold tree, restoring full headers if an id_map is provided.

    Parameters
    ----------
    treefile    : str   path to IQ-TREE .treefile (Newick)
    id_map_path : str   path to _id_map.json — if provided, short IDs
                        are replaced with full headers for readable labels
    save_path   : str   path to save figure (default: same dir as treefile)
    """
    tree = Phylo.read(treefile, "newick")

    # Restore full headers if id_map provided
    if id_map_path and Path(id_map_path).exists():
        with open(id_map_path, "r") as f:
            id_map = json.load(f)
        for clade in tree.get_terminals():
            if clade.name and clade.name in id_map:
                # Extract just genus_species and accession for readability
                full_id = id_map[clade.name]
                parts   = full_id.split("|")
                taxon   = parts[0]
                acc     = next((p.split(":",1)[1] for p in parts
                                if p.startswith("accession:")), "")
                clade.name = f"{taxon}_{acc}" if acc else taxon

    n_terminals = len(tree.get_terminals())
    fig, ax     = plt.subplots(figsize=(14, max(8, n_terminals * 0.4)))

    Phylo.draw(tree, axes=ax, do_show=False)

    gene_name = Path(treefile).stem
    ax.set_title(f"Scaffold tree — {gene_name}", fontsize=12)
    ax.set_xlabel("Branch length (substitutions per site)")
    plt.tight_layout()

    if save_path is None:
        save_path = str(treefile).replace(".treefile", "_tree.png")

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {save_path}")
    return save_path

scaffold_trees = {}

for nex_file in Path("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_nexus/Framework_alignment/genes").glob("*_short_aln.nex"):
    gene_name = nex_file.stem.replace("_short_aln", "")
    print(f"\nBuilding scaffold tree for {gene_name}...")

    result = build_scaffold_tree(
        scaffold_nexus    = str(nex_file),
        scaffold_tree_dir = f"C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_trees/scaffold_trees/{gene_name}",
        iqtree_bin        = "C:/Program Files/iqtree-3.0.1-Windows/bin/iqtree3.exe",
        seed              = 12345,
    )
    scaffold_trees[gene_name] = result
    print(f"  Tree: {result['treefile']}")


for gene_name, result in scaffold_trees.items():

    # Find the id_map — it's in the Framework_alignment folder
    id_map_path = (f"C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_nexus/Framework_alignment/genes"
                   f"{gene_name}_id_map.json")

    plot_scaffold_tree(
        treefile    = result["treefile"],
        id_map_path = id_map_path,
        save_path   = f"C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_trees/scaffold_trees/{gene_name}/{gene_name}_tree.png"
    )