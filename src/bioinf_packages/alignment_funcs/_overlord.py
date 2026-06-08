#pip install -r requirements.txt
import os
from Bio import Entrez
import sys
from Bio import SeqIO
import pandas as pd
import csv

# IMPORTANT: You need to set your email for NCBI Entrez API
Entrez.email = "om380@cam.ac.uk"



#find all new gene symbols at once - call in overlord
#write an exception
class GlossaryIncompleteError(Exception):
    pass


def update_glossary(accessions, meta_folder, glossary):

    # load glossary
    with open(glossary, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        symbol_key = [row for row in reader][0]

    missing = set()

    for accession in accessions:

        gbk_path = os.path.join(meta_folder, f"{accession}.gbk")

        if not os.path.exists(gbk_path):
            continue

        for record in SeqIO.parse(gbk_path, "genbank"):

            for feature in record.features:

                if feature.type not in ("CDS", "rRNA", "misc_feature"):
                    continue

                qualifiers = feature.qualifiers

                symbol = qualifiers.get("standard_name", [None])[0]

                if not symbol:
                    symbol = qualifiers.get("gene", [None])[0]

                if not symbol:
                    symbol = qualifiers.get("product", [None])[0]

                if not symbol:
                    continue

                if "trna" in symbol.lower():
                    continue

                if symbol not in symbol_key:
                    missing.add(symbol)

    if missing:

        # add all missing terms
        for symbol in sorted(missing):
            symbol_key[symbol] = ""

        # rewrite glossary
        with open(glossary, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=symbol_key.keys())
            writer.writeheader()
            writer.writerow(symbol_key)

        raise GlossaryIncompleteError(
            f"{len(missing)} new symbols added to glossary. "
            "Please fill in their standard names and rerun."
        )

    return True

#take a list of accession numbers and cycle through repeating steps 1-3 appending all outputs to a single dictionary
"""
- collect metadata for each accession (efetch_gene)
- read metadata into correct format (read_gbk)
- read fasta file in conjunction with metadata and generate dictionary with sequences for specific accession
- append this dictionary to existing dictionary for all accession numbers
"""

from ._extract_accessions import extract_accessions
from ._fetch_fasta import fetch_fasta
from ._efetch_gene import efetch_gene
from ._read_gbk import read_gbk
from ._seq_finder import seq_finder


def overlord_function(raw_file: str, meta_file: str, csv_path: str, glossary: str):
    #get list of accession numbers
    accessions = extract_accessions(csv_path)

    #set empty dictionary
    dict_list = []
    #n=1
    #loop through accessions to get gene symbols and .gbks
    for value in accessions:
        #n=n+1

        #install meta data for accessions
        efetch_gene(value, meta_file)
    
    #scan for unknown gene symbols
    update_glossary(accessions, meta_file, glossary)

    #loop through accesions to do everything else
    for value in accessions:
        #install fasta files for accessions
        fetch_fasta(value, raw_file)
    
        #read meta data
        meta_data = read_gbk(value, meta_file, glossary)
        
        #add sequences to meta data
        mini_dict = seq_finder(meta_data, raw_file)

        #loop through mini dict incase it's a list with multiple values and append them
        for i in range(len(mini_dict)):
            #print(n)
            #print((mini_dict))
            #append to main dictionary
            dict_list.append(mini_dict[i])
    
    return(dict_list)
            

#test1 = (overlord_function("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/raw_fastas", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/newGenBank.csv"))
#print(test1)
#print([d for d in test1 if d.get("symbol") == "small subunit ribosomal RNA"])