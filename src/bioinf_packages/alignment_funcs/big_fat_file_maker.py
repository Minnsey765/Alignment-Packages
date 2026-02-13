#function to concatenate all sequences into a single sequence to be aligned at once

import os
import sys
try: #Import from package (elsewhere)
    from ._matcher import find_best_key
except ImportError: #import locally
    from _matcher import find_best_key
    
import pandas as pd
from Bio import SeqIO
from pathlib import Path

def big_fat_file_maker(raw_fol: str, csv_path: str, output: str, datatype: str):

    #make list of all species names
    csv = pd.read_csv(csv_path)
    species_names = csv["Species"].tolist()
    #make list of all gene names (ignore species column)
    gene_names = csv.columns[1:].tolist()

    #convert list to dictionary where entries in list are dictionary keys
    my_dict = {item: "" for item in species_names}
    folder = Path(raw_fol)
    for name in gene_names:
        matching_files = list(folder.glob(f"*{name}_aln*.nex"))
        
        if not matching_files:
            print(f"No NEXUS file found for {name}")
            continue

        for file_path in matching_files:
            print(f"Reading {file_path.name}")
            gene_seqs = {}
            for record in SeqIO.parse(file_path, "nexus"):
                #ignore extra genes for a single species
                if record.id.endswith(".copy"):
                    continue
                    #print(f">{record.id}")

                #use matching function to locate correct key
                match = find_best_key(record.id, my_dict)
                if match:
                    gene_seqs[match] = str(record.seq)
                else:
                    #skip if no match
                    continue
            # Determine alignment length for this gene
            if gene_seqs:
                gene_length = len(next(iter(gene_seqs.values())))
            else:
                gene_length = 0

            # Append sequences or gaps to my_dict to maintain alignment
            for species in species_names:
                if species in gene_seqs:
                    my_dict[species] += gene_seqs[species]
                else:
                    my_dict[species] += "-" * gene_length
                #print(len(my_dict[species]))
    
    
    os.makedirs(output, exist_ok=True)

    #write actual nexus file
    # Check that all sequences have the same length
    lengths = [len(seq) for seq in my_dict.values()]
    if len(set(lengths)) != 1:
        raise ValueError("All sequences must have the same length!")
    
    nchar = lengths[0]
    ntax = len(my_dict)
    output_file = f"{output}/concatenated_aln.nex"

    #normalise name lengths for alignment
    #find longest name length
    max_name_len = max(len(name) for name in my_dict.keys())

    print("Writing to Nexus")
    with open(output_file, "w") as nex:
        nex.write("#NEXUS\n")
        nex.write("Begin DATA;\n")
        nex.write(f"    Dimensions ntax={ntax} nchar={nchar};\n")
        nex.write(f"    Format datatype={datatype} missing=? gap=-;\n")
        nex.write("    Matrix\n")
        
        # Write each species and its sequence
        for species, seq in my_dict.items():
            padded_name = species.ljust(max_name_len) # pad with spaces
            nex.write(f"{padded_name} {seq}\n")
        
        nex.write("    ;\n")
        nex.write("End;\n")
    print("Finished!")
    #return(my_dict)

    


#(big_fat_file_maker("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/edited_nexus", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/newGenBank.csv", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/edited_nexus", "dna"))