#lowk a useless function
#generate a match score for hits found in verify seq

#try: #Import from package (elsewhere)
from ._verify_seq import verify_seq
from ._species_parser import species_parser
from ._gene_parser import gene_parser
#except ImportError: #Import locally
#    from bioinf_packages.verify_funcs._verify_seq import verify_seq
#    from bioinf_packages.verify_funcs._species_parser import species_parser
#    from bioinf_packages.verify_funcs._gene_parser import gene_parser
from Bio import SeqIO
import os

def loop_verify(fasta_file: str, n: int, output: str):
    fasta = SeqIO.parse(fasta_file, "fasta")
    gene_name = os.path.splitext(os.path.basename(fasta_file))[0]
    full_path = os.path.join(output,f"{gene_name}.txt")
    #set match score dictionary
    match_scores = {}
    for record in fasta:
        seq = record.seq
        #print(type(seq))
        species,gene,orien,identity = record.description.split("|")
        useless,accession = identity.split(":")
        #print(species)
        #print(gene)
        #print(accession)
        matches = 0
        print(f"Collecting Blast Data For: {accession}")
        verify_dic = verify_seq(seq,n)
        print(f"Collected Blast Data For {accession}")
        hit_ids = verify_dic["Accessions"]
        hit_des = verify_dic["Hit_Des"]
        #loop through hit ids 
        for id in hit_ids:
            if accession in id:
                matches = matches +1
        #loop through descriptions (species and gene)
        for des in hit_des:
            matches = matches + species_parser(species, des) + gene_parser(gene, des)
        match_score = matches #max points for given record = 4 (1 for match genus, species, gene, & accession)
                                                 #max points for all is therefore 3*number of records/hits + 1 (only accession match)
        print(f"Calculated Match Score for {accession} ({match_score})")
        #dynamically make key names for each accession and add match score
        match_scores[accession] = match_score
        print(f"Added to Dictionary")
    print("Finished! Writing to .txt file...")
    #write to file
    with open(full_path, "w") as f:
        for key, value in match_scores.items():
            f.write(f"{key}\t{value}\n")

    return(match_scores)

print(loop_verify("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data/APOB.fasta",50,"C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_verify"))