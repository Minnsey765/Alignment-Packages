#verify whether sequences are reliable

from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import random
import time

Entrez.email = "om380@cam.ac.uk"
#API_KEY = "cd463c1a96153058c1010d089da7d2359907"


def blast_with_retry(sequence, max_retries=3, timeout_secs=600):
    """
    Submit a BLAST query with automatic retry on timeout.
    
    Parameters
    ----------
    sequence    : str   nucleotide sequence to BLAST
    max_retries : int   number of times to retry before giving up
    timeout_secs: int   seconds to wait before considering it timed out
    """
    for attempt in range(1, max_retries + 1):
        try:
            print(f"    BLAST attempt {attempt}/{max_retries}...")
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=sequence,
                hitlist_size=10,
            )
            return result_handle

        except Exception as e:
            print(f"    Attempt {attempt} failed: {e}")
            if attempt < max_retries:
                wait = 30 * attempt   # back off: 30s, 60s, 90s
                print(f"    Waiting {wait}s before retry...")
                time.sleep(wait)
            else:
                raise RuntimeError(
                    f"BLAST failed after {max_retries} attempts"
                ) from e
            
        
def verify_seq(seq: str, n: int, query_accession: str,
               n_samples: int = 3) -> dict:
    """
    Draw n_samples independent random subsequences from seq, BLAST each,
    and return all hits together with per-sample metadata.

    Multiple subsamples serve two purposes:
        1. Increases the probability that at least one subsample spans a
           chimeric junction
        2. Variance in results across samples is itself a chimera signal —
           genuine sequences produce consistent hits; chimeras produce
           bimodal hit distributions

    Parameters
    ----------
    seq             : str   full query sequence
    n               : int   subsample length in bp
    query_accession : str   accession of the query sequence (for self-hit
                            removal)
    n_samples       : int   number of independent subsamples (default 3)

    Returns
    -------
    dict with keys:
        "samples"      : list of per-sample dicts, each containing:
                            "sample_index"  : int
                            "query_start"   : int
                            "query_seq"     : str
                            "hits"          : list of hit dicts
        "all_hits"     : flat list of all hits across all samples, with
                         "sample_index" added to each hit for traceability
        "n_samples"    : int
        "sample_length": int
    """
    if n > len(seq):
        raise ValueError(
            f"Subsample length {n} exceeds sequence length {len(seq)}."
        )

    samples  = []
    all_hits = []

    for sample_idx in range(n_samples):

        # ── 1. Draw random subsequence ────────────────────────────────────────
        start  = random.randint(0, len(seq) - n)
        sample = seq[start:start + n]

        print(f"  Sample {sample_idx + 1}/{n_samples} | "
              f"position {start}-{start + n} | fetching BLAST results...")

        # ── 2. Run BLAST ──────────────────────────────────────────────────────
        result_handle = blast_with_retry(sample, max_retries=3)

        # Be polite to NCBI — wait between requests
        if sample_idx < n_samples - 1:
            time.sleep(10)

        # ── 3. Save XML ───────────────────────────────────────────────────────
        xml_file = f"blast_sample_{sample_idx}.xml"
        with open(xml_file, "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

        # ── 4. Parse results ──────────────────────────────────────────────────
        sample_hits = []

        with open(xml_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)

            for record in blast_records:
                for alignment in record.alignments[:10]:

                    # Skip self-hit by accession
                    if query_accession and \
                       alignment.accession == query_accession:
                        continue

                    if not alignment.hsps:
                        continue

                    best_hsp = max(alignment.hsps, key=lambda h: h.bits)
                    pident   = (best_hsp.identities /
                                best_hsp.align_length) * 100

                    # Skip perfect full-length matches (duplicate accessions)
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

        print(f"  ✓ Sample {sample_idx + 1} | {len(sample_hits)} hits returned")

        samples.append({
            "sample_index": sample_idx,
            "query_start":  start,
            "query_seq":    sample,
            "hits":         sample_hits,
        })

    return {
        "samples":       samples,
        "all_hits":      all_hits,
        "n_samples":     n_samples,
        "sample_length": n,
    }

#print(verify_seq("CACTTCCTTTGGATATGYTTGATGTGTTTTTGAATCATAATATCAATTCCTTTCTGAGGCAAGTTGAGAAGGTCAGAGATGAGGCATTGGTTCTTGTTATTCAATCCTATAATGAAGCAAAAATGAAATTTGATGAGCATAAGGTTGAAAAATCTATCACCCAACAACGAAAGACCTTTCAAATTCCAGGGTACACCATTCCTGTTGTTAATGTCGAAGTGTCTCCATTCACAGTAGAGATGTTTCCATTTGGTTATGTGATCCCAAAGGAGGTCAGCACCCCAAAGTTCACCATCCTGGGTTCTGGTTTCTCTGTGCCTTCCTATACTTTAGTCCTGCCCTTTCTAGAACTACCAGCTCTTCATATCCCTAAGTTTCTTGAGCTTTCTTTTCCAGACTTCAAAGTATCGAGTATCCCAAGGAATATTTTCATTCCAGCCCTGGGAAATGTTACATATGATTTTTCCTTTAAGTCAAGTGTCATTACACTGAATGCCAATGCTGGACTTTAT", 50, "GU981106"))