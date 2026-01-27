import os
from Bio import Entrez
import sys

# IMPORTANT: You need to set your email for NCBI Entrez API
Entrez.email = "om380@cam.ac.uk"

#find fasta based on accession and directory to save it to
def fetch_fasta(accession: str, target_dir: str) -> str:
    """
    Download FASTA file for a given accession number if it does not already exist.

    Args:
        accession (str): NCBI accession number.
        target_dir (str): Directory to store the fasta file.

    Returns:
        str: Path to the fasta file.
    """
    os.makedirs(target_dir, exist_ok=True)

    fasta_path = os.path.join(target_dir, f"{accession}.fasta")

    # Check if file already exists
    if os.path.exists(fasta_path):
        print(f"✅ FASTA already exists: {fasta_path}")
        return fasta_path

    # Download fasta
    try:
        print(f"⬇️ Downloading FASTA for {accession}...")
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
        with open(fasta_path, "w") as f:
            f.write(handle.read())
        handle.close()
        print(f"✅ Download complete: {fasta_path}")
    except Exception as e:
        raise RuntimeError(f"❌ Failed to fetch FASTA for {accession}: {e}")

    return fasta_path

#print(fetch_fasta("NC_024563.1", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/raw_fastas"))
# Make callable from command line
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: python {sys.argv[0]} <accession> <output_directory>")
        sys.exit(1)

    accession = sys.argv[1]
    target_dir = sys.argv[2]

    fetch_fasta(accession, target_dir)