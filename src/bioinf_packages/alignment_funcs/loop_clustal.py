#run through all fasta files and align them

try: #import from package (elsewhere)
    from ._clustalOmega import run_clustalo
except ImportError: #Import locally
    from _clustalOmega import run_clustalo

import time
from pathlib import Path

def loop_clustal(fasta_dir: str, output_dir: str, interval):
    #loop through fasta files in directory
    folder = Path(fasta_dir)
    
    # Loop through all fasta files in the folder
    for fasta_file in folder.glob("*.fasta"):
        run_clustalo(fasta_file, output_dir)
        #print(f"Waiting {interval} seconds before next file...")
        #time.sleep(interval)  # pause for 30 seconds

# Usage
#loop_clustal("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/aligned_fastas", interval=30)