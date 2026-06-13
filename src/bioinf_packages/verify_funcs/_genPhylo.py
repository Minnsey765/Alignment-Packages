from _phylo_verify import (build_scaffold_tree, plot_scaffold_tree)
from pathlib import Path
from bioinf_packages.verify_funcs._species_parser import load_taxonomy


concatenated_nex = (
    "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/"
    "verify_data/msa_verify/msaVerify_nexus/Framework_alignment/"
    "concatenated/concatenated_aln.nex"
)

scaffold_tree_dir = (
    "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/"
    "verify_data/msa_verify/msaVerify_nexus/Framework_alignment/"
    "concatenated/scaffold_tree"
)

taxonomy_path = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv"
iqtree_bin    = "C:/Program Files/iqtree-3.0.1-Windows/bin/iqtree3.exe"

# ── Load taxonomy (needed for constraint tree) ────────────────────────────────

taxonomy = load_taxonomy(taxonomy_path)

# ── Build scaffold tree ───────────────────────────────────────────────────────
# The concatenated Nexus already uses plain species names (not short IDs)
# so IQ-TREE can read it directly. No id_map needed here.
# Pass taxonomy and the Nexus path to generate a constraint tree that
# enforces order-level monophyly and prevents long-branch attraction.

result = build_scaffold_tree(
    scaffold_nexus         = concatenated_nex,
    scaffold_tree_dir      = scaffold_tree_dir,
    taxonomy               = taxonomy,
    scaffold_fasta_aligned = None,   # concatenated Nexus has plain names
                                     # so no FASTA needed for constraint —
                                     # see note below
    iqtree_bin             = iqtree_bin,
    model                  = "TEST",   # or "TEST" for auto model selection
    bootstrap              = 1000,
    n_threads              = 4,
    seed                   = 12345,
)

print(f"\nScaffold tree: {result['treefile']}")

# ── Plot scaffold tree ────────────────────────────────────────────────────────
# The concatenated Nexus uses full species names so no id_map is needed —
# the labels are already readable.

plot_scaffold_tree(
    treefile    = result["treefile"],
    id_map_path = None,          # not needed — names are already full
    save_path   = str(
        Path(scaffold_tree_dir) / "concatenated_scaffold_tree.png"
    ),
)





#scaffold_trees = {}

#for nex_file in Path("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_nexus/Framework_alignment/genes").glob("*_short_aln.nex"):
#    gene_name = nex_file.stem.replace("_short_aln", "")
#    print(f"\nBuilding scaffold tree for {gene_name}...")

#    result = build_scaffold_tree(
#        scaffold_nexus    = str(nex_file),
#        scaffold_tree_dir = f"C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_trees/scaffold_trees/{gene_name}",
#        iqtree_bin        = "C:/Program Files/iqtree-3.0.1-Windows/bin/iqtree3.exe",
#        seed              = 12345,
#    )
#    scaffold_trees[gene_name] = result
#    print(f"  Tree: {result['treefile']}")


#for gene_name, result in scaffold_trees.items():

    # Find the id_map — it's in the Framework_alignment folder
#    id_map_path = (f"C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_nexus/Framework_alignment/genes"
#                   f"{gene_name}_id_map.json")

#    plot_scaffold_tree(
#        treefile    = result["treefile"],
#        id_map_path = id_map_path,
#        save_path   = f"C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_trees/scaffold_trees/{gene_name}/{gene_name}_tree.png"
#    )