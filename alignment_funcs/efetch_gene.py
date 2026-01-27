#fetch gene with non-ncbi supported accession number (KC, AC, etc)
import os
import sys
#pip install biopython in venv
from Bio import Entrez

# Required by NCBI Entrez
Entrez.email = "om380@cam.ac.uk"

#install meta data for given accession number
def efetch_gene(accession: str, target_dir: str):
    """
    Fetches a gene sequence by accession using Entrez efetch.
    Saves it as a GenBank (.gbk) file inside target_dir.
    Skips download if file already exists.
    """
    os.makedirs(target_dir, exist_ok=True)

    # Save as <accession>.gbk inside target_dir
    gbk_path = os.path.join(target_dir, f"{accession}.gbk")

    # ✅ Skip download if file already exists
    if os.path.exists(gbk_path):
        print(f"⚠️ GenBank file already exists: {gbk_path}, skipping download.")
        return gbk_path

    try:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        with open(gbk_path, "w") as f:
            f.write(handle.read())
        handle.close()

        print(f"✅ GenBank accession {accession} downloaded to {gbk_path}")
        return gbk_path

    except Exception as e:
        print(f"❌ efetch failed for {accession}: {e}")
        return None

# Allow function to be run directly from command line
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: python {sys.argv[0]} <gene_accession> <output_directory>")
        sys.exit(1)

    accession = sys.argv[1]
    output_dir = sys.argv[2]

    efetch_gene(accession, output_dir)