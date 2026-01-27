#align sequences using clustal omega
import requests
import time
import os

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"

def run_clustalo(fasta_file, out_fold):

    #generate name for output
    #get rid of path and .fasta bit
    base = os.path.splitext(os.path.basename(fasta_file))[0]
    new_name = f"{base}_aln.nex"

    out_file = os.path.join(out_fold, new_name)
    # Read FASTA
    with open(fasta_file, "r") as f:
        fasta_data = f.read()
    
    # Submit job (do NOT include outfmt="nexus")
    params = {
        "sequence": fasta_data,
        "stype": "dna",   # or "protein"
        "email": "om380@cam.ac.uk",
        "outfmt": "nexus",
        "order": "aligned"
    }
    response = requests.post(BASE_URL + "/run", data=params)
    response.raise_for_status()
    job_id = response.text.strip()
    print(f"Job submitted. ID: {job_id}")
    
    # Poll for status
    status = "RUNNING"
    i = 0
    while status in ["RUNNING", "PENDING"]:
        i = i+1
        time.sleep(10)
        status = requests.get(BASE_URL + f"/status/{job_id}").text.strip()
        print(f"Job status: {status}")
        #dont get stuck in never ending loop
        if i == 100:
            break
    
    if status != "FINISHED":
        raise RuntimeError(f"Job failed with status: {status}")
    
    # Fetch results as NEXUS
    result = requests.get(BASE_URL + f"/result/{job_id}/aln-nexus")
    result.raise_for_status()
    
    # Save to file
    with open(out_file, "w+") as f:
        f.write(result.text)
    print(f"Alignment saved to {out_file}")

run_clustalo("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data/BDNF.fasta", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/aligned_fastas")