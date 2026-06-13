#generate a database of sequences that are correct (gold set)

import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#from _assess_peer_review import loop_peer_review
from ..alignment_funcs._extract_accessions import extract_accessions


def _strip_version(acc):
    """
    Removes trailing version like .1, .2, etc.
    Example: AF069539.1 -> AF069539
    """
    return acc.split(".")[0]


def make_gold_set(accession_list, input_dir, output_dir, output_filename):

    """
    accession_list can be either:

    [
        {"accession": "AF069539.1"},
        {"accession": "AY123456.2"}
    ]

    or

    [
        "AF069539.1",
        "AY123456.2"
    ]
    """

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)

    # Normalize accessions
    accessions = set()

    for item in accession_list:
        if isinstance(item, dict):
            if "accession" in item:
                accessions.add(_strip_version(item["accession"]))
        else:
            accessions.add(_strip_version(item))

    matched_records = []

    for filename in os.listdir(input_dir):
        if not filename.endswith(".fasta"):
            continue

        path = os.path.join(input_dir, filename)

        for record in SeqIO.parse(path, "fasta"):
            header = record.description

            # Try to extract accession field explicitly if present
            # e.g. accession:AF069539.1
            match = re.search(r"accession:([A-Za-z0-9_.]+)", header)

            if match:
                acc = _strip_version(match.group(1))
            else:
                # fallback: search anywhere in header
                # still normalize by stripping version if found
                tokens = header.split("|")
                acc = None
                for t in tokens:
                    if "accession:" in t:
                        acc = _strip_version(t.split(":")[1])
                        break

                if acc is None:
                    continue

            if acc in accessions:
                matched_records.append(record)

    SeqIO.write(matched_records, output_path, "fasta")

    print(f"Found {len(matched_records)} matching sequences")
    print(f"Wrote output to {output_path}")



#make gold set using pubmed 
#make_gold_set(loop_peer_review("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info"), "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq", "goldset_seqs.fasta")

#make gold set using pre-picked accession numbers
#make_gold_set(extract_accessions("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_generation/artefact_sources.csv"), "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq", "goldset_seqs.fasta")