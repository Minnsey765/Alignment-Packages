import csv

my_dict = {
    "12S ribosomal RNA": ["12S_rRNA"],
    "16S ribosomal RNA": ["16S_rRNA"],
    "APOB" : ["APOB"],
    "ATP6" : ["ATP6"],
    "ATP8" : ["ATP8"],
    "ApoB" : ["APOB"],
    "BDNF" : ["BDNF"],
    "BRCA1" : ["BRCA1"],
    "COI" : ["COX1"],
    "COX1" : ["COX1"],
    "COX2" : ["COX2"],
    "COX3" : ["COX3"],
    "CYTB" : ["CYTB"],
    "Cytb" : ["CYTB"],
    "GHR" : ["GHR"],
    "Ghr" : ["GHR"],
    "NADH2" : ["ND2"],
    "ND1" : ["ND1"],
    "ND2" : ["ND2"],
    "ND3" : ["ND3"],
    "ND4" : ["ND4"],
    "ND4L" : ["ND4L"],
    "ND5" : ["ND5"],
    "ND6" : ["ND6"],
    "RAG1" : ["RAG1"],
    "RAG2" : ["RAG2"],
    "Rag1" : ["RAG1"],
    "Rag2" : ["RAG2"],
    "apoB" : ["APOB"],
    "apoB-100" : ["APOB"],
    "apob" : ["APOB"],
    "apolipoprotein B" : ["APOB"],
    "bdnf" : ["BDNF"],
    "brca1" : ["BRCA1"],
    "breast cancer susceptibility 1" : ["BRCA1"],
    "cytb" : ["CYTB"],
    "cytochrome b" : ["CYTB"],
    "cytochrome c oxidase subunit I" : ["COX1"],
    "growth hormone receptor" : ["GHR"],
    "l-rRNA" : ["16S_rRNA"],
    "large subunit ribosomal RNA" : ["16S_rRNA"],
    "nad2" : ["ND2"],
    "rag1" : ["RAG1"],
    "rag2" : ["RAG2"],
    "recombination activating protein 2" : ["RAG2"],
    "rrn12" : ["12S_rRNA"],
    "rrn16" : ["16S_rRNA"],
    "rrnl" : ["16S_rRNA"],
    "rrns" : ["12S_rRNA"],
    "s-rRNA" : ["12S_rRNA"],
    "small subunit ribosomal RNA" : ["12S_rRNA"]
}

with open("src/bioinf_packages/dictionary_funcs/glossary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, my_dict.keys())
    w.writeheader()

    #stop it from writing "['WXYZ']" as a string into the dictionary
    rm_brkt = {
        k : v[0] if isinstance(v,list) and v else ""
        for k, v in my_dict.items()
    }
    w.writerow(rm_brkt)