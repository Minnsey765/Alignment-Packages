from Bio import SeqIO
import os

#try: #import from package (elsewhere)
#    from ._read_gbk import read_gbk
#except ImportError: #import locally
#    from _read_gbk import read_gbk


#takes output from read_gbk and replaces start and stop with actual sequence
def seq_finder(meta_data: list[dict], raw_file: str):
    # identify correct fasta sequence
    if not meta_data:
        return []
    accession = meta_data[0].get("accession")
    # Look for any fasta file containing accession in its filename
    fasta_path = None
    for fname in os.listdir(raw_file):
        if accession in fname and fname.lower().endswith((".fasta", ".fa", ".fna")):
            fasta_path = os.path.join(raw_file, fname)
            break

    if fasta_path is None:
        raise FileNotFoundError(f"⚠️ No FASTA file containing accession {accession} found in {raw_file}")
    
    # Load the FASTA sequence
    record = SeqIO.read(fasta_path, "fasta")
    seq = record.seq

    updated_genes = []
    for gene in meta_data:
        start = gene.get("begin")
        end = gene.get("end")
        #strand = gene.get("orientation", 1)
        if start is None or end is None:
            continue  # skip if missing coords
        
        # Biopython uses 0-based indexing, GenBank is 1-based
        subseq = seq[start:end]  
        # Reverse complement if on the minus strand
        #if strand == -1:
        #    subseq = subseq.reverse_complement()

        # Build new dictionary: keep symbol, species, etc., but replace coords
        new_gene = gene.copy()
        new_gene.pop("begin", None)
        new_gene.pop("end", None)
        new_gene["sequence"] = str(subseq)

        updated_genes.append(new_gene)

    return updated_genes

#print(seq_finder(read_gbk("KR711200", "fasta_info"), "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/raw_fastas"))