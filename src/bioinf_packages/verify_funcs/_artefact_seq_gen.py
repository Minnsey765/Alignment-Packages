# take two sequences, A and B, and replace A with subsequnces from B

import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def replace_subsequences(A: str, B: str,
                         min_len_frac: float = 0.20,
                         max_len_frac: float = 0.40):
    """
    Replace a single random subsequence in A with the corresponding
    region from B. Insert length is 20-40% of A's total length.

    Parameters
    ----------
    A             : str    target sequence to modify
    B             : str    donor sequence to take insert from
    min_len_frac  : float  minimum insert as fraction of len(A) (default 0.20)
    max_len_frac  : float  maximum insert as fraction of len(A) (default 0.40)

    Returns
    -------
    str   modified version of A
    list  replacement events as [(start, end, length)]
    """
    if len(A) == 0 or len(B) == 0:
        raise ValueError("A and B must be non-empty.")

    # Insert length as fraction of A's length
    min_len = max(1, int(len(A) * min_len_frac))
    max_len = int(len(A) * max_len_frac)

    # Can only replace up to the length of the shorter sequence
    max_common_pos = min(len(A), len(B))

    if max_common_pos < min_len:
        raise ValueError(
            f"Sequences too short for a {min_len_frac*100:.0f}% insert: "
            f"A={len(A)}bp, B={len(B)}bp, min_insert={min_len}bp"
        )

    # Clamp max_len to what's available
    max_len = min(max_len, max_common_pos)

    if max_len < min_len:
        max_len = min_len

    length = random.randint(min_len, max_len)

    max_start = max_common_pos - length
    if max_start < 0:
        raise ValueError("Cannot fit insert — sequences too short.")

    start = random.randint(0, max_start)
    end   = start + length

    # Build modified sequence
    result = A[:start] + B[start:end] + A[end:]
    events = [(start, end, length)]

    return result, events


def generate_artefact_seqs(input_dir: str,
                           output_dir: str,
                           output_filename: str,
                           min_len_frac: float = 0.20,
                           max_len_frac: float = 0.40):
    """
    For each FASTA file in input_dir:
        - First sequence is the donor (B)
        - All subsequent sequences are targets (A) to be modified

    Each A sequence gets one chimeric replacement event.
    Modified sequences are written to a single output FASTA.

    Header format (modified from original):
        species|gene_modified_N|orien:+/-|accession:XXXX
    Where N is the index of the A sequence within that file (1-based).

    Parameters
    ----------
    input_dir        : str    directory containing input FASTA files
    output_dir       : str    directory to write output FASTA
    output_filename  : str    filename for output FASTA
    min_len_frac     : float  minimum insert length as fraction of sequence
    max_len_frac     : float  maximum insert length as fraction of sequence
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)

    all_modified_records = []

    for filename in sorted(os.listdir(input_dir)):
        if not filename.endswith(".fasta"):
            continue

        path    = os.path.join(input_dir, filename)
        records = list(SeqIO.parse(path, "fasta"))

        if len(records) < 2:
            print(f"  Skipping {filename} — needs at least 2 sequences "
                  f"(found {len(records)})")
            continue

        # First record is always the donor
        rec_B = records[0]
        seq_B = str(rec_B.seq)

        # All subsequent records are targets
        target_records = records[1:]

        gene_name = os.path.splitext(filename)[0]

        for idx, rec_A in enumerate(target_records, start=1):
            seq_A = str(rec_A.seq)

            try:
                modified_A, events = replace_subsequences(
                    seq_A, seq_B,
                    min_len_frac=min_len_frac,
                    max_len_frac=max_len_frac,
                )
            except ValueError as e:
                print(f"  Skipping sequence {idx} in {filename}: {e}")
                continue

            # ── Parse and rebuild header ──────────────────────────────────────
            # Expected: species|gene|orien:+/-|accession:XXXX
            header_parts = rec_A.description.split("|")

            if len(header_parts) < 4:
                print(f"  Warning: unexpected header in {filename}: "
                      f"{rec_A.description} — skipping")
                continue

            species  = header_parts[0]
            rest     = header_parts[2:]   # orien and accession fields

            # Number the gene field so multiple artefacts from the same
            # source sequence are distinguishable
            new_gene_field = f"{gene_name}|mod_{idx}"
            new_header     = "|".join([species, new_gene_field] + rest)

            start, end, length = events[0]
            print(f"  {filename} | target {idx}/{len(target_records)} | "
                  f"insert: {start}-{end} ({length}bp, "
                  f"{length/len(seq_A)*100:.1f}% of sequence)")

            all_modified_records.append(
                SeqRecord(
                    Seq(modified_A),
                    id=new_header,
                    description=""
                )
            )

    SeqIO.write(all_modified_records, output_path, "fasta")
    print(f"\nWrote {len(all_modified_records)} sequences to {output_path}")

#generate_artefact_seqs("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_generation/fasta_data", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_seq", "modifed_seqs.fasta", 0.20, 0.40)