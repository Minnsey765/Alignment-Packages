#align sequences using clustal omega
import os
import re
import json
import time
import requests
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"


def build_short_ids(records: list) -> tuple:
    """
    Replace full pipe-delimited headers with short safe IDs for
    ClustalOmega (which truncates IDs longer than 30 characters),
    and return a lookup dict to restore them afterwards.

    Short ID format:
        Scaffold : "S{N}_{Genus}_{AccSuffix}"   e.g. "S01_Sorex_AB1234"
        Query    : "Q{N}_Smp{M}_{AccSuffix}"    e.g. "Q01_Smp3_LC6441"

    Parameters
    ----------
    records : list   SeqRecord objects with full headers as .id

    Returns
    -------
    tuple of:
        short_records : list   SeqRecord objects with short IDs
        id_map        : dict   short_id -> original_full_id
        reverse_map   : dict   original_full_id -> short_id
    """
    short_records = []
    id_map        = {}
    reverse_map   = {}
    scaffold_n    = 0
    query_n       = 0

    for record in records:
        full_id = record.id
        parts   = full_id.split("|")
        is_query = "Sample" in full_id

        # Extract accession from header
        acc = ""
        for p in parts:
            if p.startswith("accession:"):
                acc = p.split(":", 1)[1].strip()
                break

        if is_query:
            sample_match = re.search(r'Sample(\d+)', full_id)
            sample_num   = sample_match.group(1) if sample_match else "X"
            query_n     += 1
            short_id     = f"Q{query_n:02d}_Smp{sample_num}_{acc[:6]}"
        else:
            genus      = parts[0].split("_")[0][:8] if parts else full_id[:8]
            scaffold_n += 1
            short_id    = f"S{scaffold_n:02d}_{genus}_{acc[-6:]}"

        # Ensure uniqueness
        base_short = short_id
        counter    = 1
        while short_id in id_map:
            short_id = f"{base_short}_{counter}"
            counter += 1

        id_map[short_id]     = full_id
        reverse_map[full_id] = short_id

        short_records.append(SeqRecord(
            record.seq, id=short_id, description=""
        ))

    return short_records, id_map, reverse_map


def restore_headers_in_file(file_path: str,
                             id_map: dict,
                             output_path: str = None) -> str:
    """
    Replace short IDs in any text file (Nexus, Newick, FASTA) with
    original full headers via string replacement.

    Parameters
    ----------
    file_path   : str   path to file with short IDs
    id_map      : dict  short_id -> full_id
    output_path : str   output path (default: overwrites input)

    Returns
    -------
    str   path to restored file
    """
    with open(file_path, "r") as f:
        content = f.read()

    # Replace longest keys first to avoid partial replacements
    for short_id, full_id in sorted(id_map.items(),
                                     key=lambda x: len(x[0]),
                                     reverse=True):
        content = content.replace(short_id, full_id)

    output_path = output_path or file_path
    with open(output_path, "w") as f:
        f.write(content)

    return output_path




def run_clustalo(fasta_file: str,
                 out_fold: str,
                 email: str = "om380@cam.ac.uk",
                 max_polls: int = 120,
                 poll_interval: int = 10) -> dict:

    os.makedirs(out_fold, exist_ok=True)
    base = os.path.splitext(os.path.basename(fasta_file))[0]

    # ── Remap headers ─────────────────────────────────────────────────────────
    records                  = list(SeqIO.parse(fasta_file, "fasta"))
    short_records, id_map, _ = build_short_ids(records)

    map_path = os.path.join(out_fold, f"{base}_id_map.json")
    with open(map_path, "w") as f:
        json.dump(id_map, f, indent=2)
    print(f"ID map saved: {map_path}")

    short_fasta = os.path.join(out_fold, f"{base}_short.fasta")
    SeqIO.write(short_records, short_fasta, "fasta")

    with open(short_fasta, "r") as f:
        fasta_data = f.read()

    # ── Submit ONCE ───────────────────────────────────────────────────────────
    print(f"\nSubmitting to EBI ClustalOmega...")

    params = {
        "sequence": fasta_data,
        "stype":    "dna",
        "email":    email,
        "outfmt":   "nexus",
        "order":    "aligned",
    }
    response = requests.post(BASE_URL + "/run", data=params)
    response.raise_for_status()
    job_id = response.text.strip()
    print(f"  Job ID: {job_id}")

    # ── Poll ──────────────────────────────────────────────────────────────────
    status = "RUNNING"
    polls  = 0
    while status in ("RUNNING", "PENDING") and polls < max_polls:
        time.sleep(poll_interval)
        status = requests.get(
            BASE_URL + f"/status/{job_id}"
        ).text.strip()
        polls += 1
        print(f"  [{polls}] {status}")

    if status != "FINISHED":
        raise RuntimeError(
            f"ClustalOmega job {job_id} ended with status: {status}"
        )

    # ── Fetch both formats ────────────────────────────────────────────────────
    # ↓↓↓ THIS IS THE SECTION THAT REPLACES YOUR EXISTING FETCH LOOP ↓↓↓
    results = {}

    for fmt, result_key, suffix in [
        ("nexus", "aln-nexus", "_aln.nex"),
        ("fasta", "aln-fasta", "_aln.fasta"),
    ]:
        print(f"  Fetching {fmt} result...")

        short_out = os.path.join(out_fold, f"{base}_short{suffix}")
        full_out  = os.path.join(out_fold, f"{base}{suffix}")

        try:
            result = requests.get(BASE_URL + f"/result/{job_id}/{result_key}")
            result.raise_for_status()

            with open(short_out, "w") as f:
                f.write(result.text)

            restore_headers_in_file(short_out, id_map, full_out)
            print(f"  Short-ID ({fmt}): {short_out}")
            print(f"  Full headers ({fmt}): {full_out}")

        except requests.exceptions.HTTPError as e:
            if fmt == "fasta":
                # Fall back: convert nexus to fasta locally
                print(f"  EBI fasta fetch failed ({e}) — "
                      f"converting nexus to fasta locally")
                nexus_short = results["nexus"]["short"]
                alignment   = AlignIO.read(nexus_short, "nexus")
                AlignIO.write(alignment, short_out, "fasta")
                restore_headers_in_file(short_out, id_map, full_out)
                print(f"  Converted locally — Short-ID (fasta): {short_out}")
                print(f"  Full headers (fasta): {full_out}")
            else:
                raise
        # ↑↑↑ END OF REPLACEMENT SECTION ↑↑↑

        results[fmt] = {"short": short_out, "full": full_out}

    alignment        = AlignIO.read(results["fasta"]["full"], "fasta")
    alignment_length = alignment.get_alignment_length()

    return {
        "nexus_alignment":  results["nexus"]["full"],
        "fasta_alignment":  results["fasta"]["full"],
        "short_nexus":      results["nexus"]["short"],
        "short_fasta":      results["fasta"]["short"],
        "id_map":           id_map,
        "id_map_path":      map_path,
        "n_sequences":      len(alignment),
        "alignment_length": alignment_length,
    }

#run_clustalo("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data/BDNF.fasta", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/aligned_fastas")

def loop_clustal(fasta_dir: str,
                 output_dir: str,
                 interval: int = 30,
                 email: str = "om380@cam.ac.uk"):
    """
    Run ClustalOmega on all FASTA files in a directory.

    Parameters
    ----------
    fasta_dir  : str   directory containing .fasta files
    output_dir : str   directory to write alignments
    interval   : int   seconds between submissions (default 30)
    email      : str   EBI API email
    """
    folder = Path(fasta_dir)
    fastas = list(folder.glob("*.fasta"))
    print(f"Found {len(fastas)} FASTA files in {fasta_dir}")

    for i, fasta_file in enumerate(fastas, 1):
        print(f"\n[{i}/{len(fastas)}] {fasta_file.name}")
        try:
            run_clustalo(str(fasta_file), output_dir, email=email)
        except Exception as e:
            print(f"  ✗ Failed: {e}")
            import traceback
            traceback.print_exc()
        if i < len(fastas):
            print(f"Waiting {interval}s...")
            time.sleep(interval)

# Usage
loop_clustal("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_nexus/Framework_alignment/genes")