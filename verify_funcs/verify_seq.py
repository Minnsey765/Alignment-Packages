#verify whether sequences are reliable

#from overlord import overlord_function
from Bio import SeqIO
import random
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "om380@cam.ac.uk"

def verify_seq(seq: str, n: int):
    #random start point
    start = random.randint(0, len(seq) - n)
    sample = seq[start:start+n]
    #print(sample)
    # Run blast (here: nucleotide blast against nt database)
    result_handle = NCBIWWW.qblast(
        program="blastn",          # or "blastp" for proteins
        database="nt",             # genbank nucleotide database
        sequence=sample
    )
    # Save blast results
    with open("my_blast.xml", "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()

    results = {
        "Accessions" : [],
        "Hit_Des" : [],
    }
    # Parse blast results
    with open("my_blast.xml") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for record in blast_records:
            for alignment in record.alignments[:5]:  # show top 5 
                #print(f"Hit ID: {alignment.hit_id}")
                results["Accessions"].append(alignment.accession) #add hit IDs to dictionary
                #print(f"Description: {alignment.hit_def}")
                results["Hit_Des"].append(alignment.hit_def) #add hit descriptions to dictionary
                #print(f"Length: {alignment.length}")
                #print("------")
    return(results,start)

#print(verify_seq("CACTTCCTTTGGATATGYTTGATGTGTTTTTGAATCATAATATCAATTCCTTTCTGAGGCAAGTTGAGAAGGTCAGAGATGAGGCATTGGTTCTTGTTATTCAATCCTATAATGAAGCAAAAATGAAATTTGATGAGCATAAGGTTGAAAAATCTATCACCCAACAACGAAAGACCTTTCAAATTCCAGGGTACACCATTCCTGTTGTTAATGTCGAAGTGTCTCCATTCACAGTAGAGATGTTTCCATTTGGTTATGTGATCCCAAAGGAGGTCAGCACCCCAAAGTTCACCATCCTGGGTTCTGGTTTCTCTGTGCCTTCCTATACTTTAGTCCTGCCCTTTCTAGAACTACCAGCTCTTCATATCCCTAAGTTTCTTGAGCTTTCTTTTCCAGACTTCAAAGTATCGAGTATCCCAAGGAATATTTTCATTCCAGCCCTGGGAAATGTTACATATGATTTTTCCTTTAAGTCAAGTGTCATTACACTGAATGCCAATGCTGGACTTTAT", 50))