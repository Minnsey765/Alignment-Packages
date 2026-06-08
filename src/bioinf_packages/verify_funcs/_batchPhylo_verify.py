from _phylo_verify import verify_by_phylogeny, verify_batch_by_phylogeny
from _species_parser import load_taxonomy

taxonomy = load_taxonomy("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/taxonomy_data.csv")

# ── Single sequence ───────────────────────────────────────────────────────────
result = verify_by_phylogeny(
    query_seq         = "ACTTCGGTTGCATGAAGGCTGCCCCCATGAAAGAAGCCAGCGTCCGAGGACAAGGCAGCTTGGCCTACCCAGGTGTGCGGACCCATGGGACTCTGGAGAGTGTGAATGGGCCCAAGGCAGGTGCCAGAGGCCTGACGTCCTTGGCTGACACTTTTGAACACGTGATCGAAGAGCTGTTGGAAGAGGACCAGAAAGTTCGTCCCCATGAAGAAACCAATAAGGACGCGGACTTGTACACTTCCCGGGTGATGCTGAGTAGTCAAGTGCCTTTGGAGCCTCCTCTTCTCTTTCTGCTGGAGGAATACAAAAATTACCTGGATGCTGCAAACATGTCCATGAGGGTCCGGCGCCACTCCGACCCCGCCCGCCGCGGGGAGCTGAGCGTGTGTGACAGCATCAGCGAGTGGGTGACAGCAGCGGATAAAAAGACTGCAGTGGACATGTCGGGCGGGACGGTCACTGTCCTGGAAAAAGTCCCTGTATCCAAAGGCCAACTGAAGCAGTACTTCTACGAGACCAAGTGCAATCCCATGGGTTACACGAAGGAGGGCTGCAGG",
    query_taxon       = "Solenodon_paradoxus",
    query_gene        = "BDNF",
    query_accession   = "AY530070",
    query_orientation = "+",
    scaffold_fasta    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas/BDNF.fasta",
    output_dir        = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_nexus",   # subdir created automatically
    taxonomy          = taxonomy,
    iqtree_bin        = "C:/Program Files/iqtree-3.0.1-Windows/bin/iqtree3.exe",        # or "iqtree" if that's your install
)
# Output saved to: phylo_output/Solenodon_paradoxus_LC644162_phylo/

# ── Batch ─────────────────────────────────────────────────────────────────────
#results = verify_batch_by_phylogeny(
#    query_fasta    = "sequences_to_verify.fasta",
#    scaffold_fasta = "scaffold_12S.fasta",
#    output_dir     = "phylo_output",
#    taxonomy       = taxonomy,
#    iqtree_bin     = "iqtree2",
#    interval       = 30,
#)

# Summary across all sequences
#for r in results:
#    print(f"{r['query_taxon']} {r['query_accession']}: "
#          f"chimera={r['chimera_flag']} "
#          f"ratio={r['mean_placement_ratio']}")