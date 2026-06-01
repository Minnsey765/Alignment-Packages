import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


def parse_genus_species(description: str):
    match = re.search(r"([A-Z][a-z]+)\s([a-z]+)", description)
    if match:
        return match.group(1), match.group(2)
    return None, None


def parse_accession(description: str):
    match = re.search(r"accession[:\| ]?([A-Za-z0-9_.]+)", description)
    if match:
        return match.group(1)
    return description.split()[0]


def parse_query_header(header: str):
    acc_match = re.search(r"accession:([A-Za-z0-9_.]+)", header)
    accession = acc_match.group(1) if acc_match else header.split()[0]
    return accession


def blast_fasta_compact(
    fasta_file,
    output_dir,
    top_hits=10
):

    os.makedirs(output_dir, exist_ok=True)

    rows = []

    for record in SeqIO.parse(fasta_file, "fasta"):

        query_seq = str(record.seq)
        query_acc = parse_query_header(record.description)

        print(f"\nBLASTing {query_acc}")

        row = {"query_accession": query_acc}

        try:
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=query_seq
            )

            blast_records = NCBIXML.parse(result_handle)

            hits = []

            for blast_record in blast_records:
                for alignment in blast_record.alignments[:top_hits]:

                    if not alignment.hsps:
                        continue

                    accession = parse_accession(alignment.hit_def)
                    genus, species = parse_genus_species(alignment.hit_def)

                    hit_str = f"{accession}|{genus} {species}"
                    hits.append(hit_str)

            result_handle.close()

            # ensure fixed number of columns
            for i in range(top_hits):
                row[f"hit{i+1}"] = hits[i] if i < len(hits) else None

        except Exception as e:
            print(f"Failed {query_acc}: {e}")
            row["error"] = str(e)

            for i in range(top_hits):
                row[f"hit{i+1}"] = None

        rows.append(row)

    df = pd.DataFrame(rows)

    out_path = os.path.join(output_dir, "blast_compact_results.csv")
    df.to_csv(out_path, index=False)

    print(f"\nSaved results to {out_path}")

    return df

#blast_fasta_compact("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq/goldset_seqs.fasta", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/gold_seq",10)