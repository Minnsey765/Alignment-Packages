from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import os

def verify_seq(seq: str, n: int, query_accession: str,
               n_samples: int = 3,
               db: str = "nt",
               blast_bin: str = "blastn",
               num_threads: int = 4) -> dict:
    """
    Local BLAST+ version of verify_seq(). Requires BLAST+ installed
    and the nt database downloaded and indexed.

    Parameters
    ----------
    db          : str   path to BLAST database (e.g. "/data/blastdb/nt"
                        or just "nt" if it's on your BLASTDB path)
    blast_bin   : str   path to blastn binary if not on PATH
    num_threads : int   threads for BLAST+ (default 4)
    """
    import random, math

    if n > len(seq):
        raise ValueError(
            f"Subsample length {n} exceeds sequence length {len(seq)}."
        )

    samples  = []
    all_hits = []

    for sample_idx in range(n_samples):

        start  = random.randint(0, len(seq) - n)
        sample = seq[start:start + n]

        print(f"  Sample {sample_idx + 1}/{n_samples} | "
              f"position {start}-{start+n} | running BLAST+...")

        # Write subsample to temporary FASTA
        query_file  = f"query_{sample_idx}.fasta"
        output_file = f"blast_out_{sample_idx}.xml"

        with open(query_file, "w") as f:
            f.write(f">query_{sample_idx}\n{sample}\n")

        # Run BLAST+
        blastn = NcbiblastnCommandline(
            cmd=blast_bin,
            query=query_file,
            db=db,
            outfmt=5,           # XML output
            out=output_file,
            num_alignments=10,
            num_descriptions=10,
            num_threads=num_threads,
        )
        stdout, stderr = blastn()

        if stderr:
            print(f"  BLAST+ stderr: {stderr}")

        # Parse results — identical to before
        sample_hits = []

        with open(output_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)

            for record in blast_records:
                for alignment in record.alignments[:10]:

                    if query_accession and \
                       alignment.accession == query_accession:
                        continue

                    if not alignment.hsps:
                        continue

                    best_hsp = max(alignment.hsps, key=lambda h: h.bits)
                    pident   = (best_hsp.identities /
                                best_hsp.align_length) * 100

                    if pident == 100.0 and best_hsp.align_length == n:
                        continue

                    hit = {
                        "sample_index": sample_idx,
                        "accession":    alignment.accession,
                        "description":  alignment.hit_def,
                        "bitscore":     best_hsp.bits,
                        "evalue":       best_hsp.expect,
                        "pident":       pident,
                        "align_length": best_hsp.align_length,
                    }

                    sample_hits.append(hit)
                    all_hits.append(hit)

        # Clean up temp files
        os.remove(query_file)
        os.remove(output_file)

        print(f"  ✓ Sample {sample_idx + 1} | {len(sample_hits)} hits")

        samples.append({
            "sample_index": sample_idx,
            "query_start":  start,
            "query_seq":    sample,
            "hits":         sample_hits,
        })

    return {
        "samples":        samples,
        "all_hits":       all_hits,
        "n_samples":      n_samples,
        "sample_length":  n,
    }


# Local nt database on HPC
result = verify_seq("GCTA...", n=200, query_accession="XX_00000",
                    n_samples=3,
                    db="/path/to/blastdb/nt",
                    num_threads=8)

# Or if nt is on your BLASTDB environment variable
result = verify_seq("GCTA...", n=200, query_accession="XX_00000",
                    n_samples=3,
                    db="nt")