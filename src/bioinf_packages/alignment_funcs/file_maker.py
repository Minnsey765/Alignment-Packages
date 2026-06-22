#take large dictionary and write into fasta files
import os
import sys


from bioinf_packages.alignment_funcs._overlord import overlord_function

def file_maker(sort_crit: str, raw_file: str, meta_file: str, csv_path: str, output: str, glossary: str):
    #generate dictionary
    dict = overlord_function(raw_file, meta_file, csv_path, glossary)
    #dictionary list is apparently list of lists with single dictionary entry hence [0][0]
    #print((dict[0]))

    #make sure outfile exists
    os.makedirs(output, exist_ok=True)

    #extract all unique values in sort criterior dictionary key
    values = set()
    for entry in dict:
        if sort_crit in entry and entry[sort_crit]: #only add if key exists and is not empty
            values.add(entry[sort_crit])
    keys = sorted(values)

    #paste data into fasta files with sort_crit in common
    for key in keys:
        #find list containing all dictionaries with key value in common
        fasta_entries = [d for d in dict if d.get(sort_crit) == key]
        #make fasta file for each key value
        fasta_path = os.path.join(output, f"{key}.fasta")

        #make sure it doesn't already exist etc etc
        if os.path.exists(fasta_path):
            print(f"⚠️ File already exists: {fasta_path}")
        else:
            with open(fasta_path, "w") as f: # w mode clears file
                pass

        #loop through list of entries (dicts)
        for n in fasta_entries:
            #generate fasta sequence name and sequence
            # possible should change orientation to be further down so any '-' doesn't interfere with reading nexus file
            header = f"{n["species"]}|{n["symbol"]}|orien:{n["orientation"]}|accession:{n["accession"]}"
            seq = n["sequence"]

            #append fasta entries to new file
            with open(fasta_path, "a") as f:  # "a" = append mode
                f.write(f">{header}\n")
                f.write(seq + "\n")  # write sequence in one line
            print(f"➕ Added entry '{header}' to {fasta_path}")

# Allow function to be run directly from command line
#if __name__ == "__main__":
#    if len(sys.argv) < 3:
#        print(f"Usage: python {sys.argv[0]} <sort_criterion> <raw_data_directory> <meta_data_directory> <csv_path> <output_directory>")
#        sys.exit(1)

#    sort_crit = sys.argv[1]
#    raw_file = sys.argv[2]
#    meta_file = sys.argv[3]
#    csv_path = sys.argv[4]
#    output = sys.argv[5]


#    file_maker(sort_crit, raw_file, meta_file, csv_path, output)


#file_maker("symbol", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/raw_fastas", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/newGenBank.csv", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data", "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv")


#making MSA framework for alignment testing
#file_maker("symbol", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/raw_fastas", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/MSA_framework.csv", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/msa_verify/msaVerify_fastas", "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv")

#making mus musculus fasta files for making artefactual sequences

# run this in terminal: python -m bioinf_packages.alignment_funcs.file_maker
#file_maker("symbol", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_generation/raw_fastas", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_generation/fasta_info", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_generation/artefact_sources.csv", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/artefact_generation/fasta_data", "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv")

#relaxed
file_maker(
    sort_crit = "symbol",
    raw_file = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/raw_fastas",
    meta_file = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    csv_path = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/final_accession_relaxed.csv",
    output = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data_clean/relaxed",
    glossary = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv",
)
#harsh
file_maker(
    sort_crit = "symbol",
    raw_file = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/raw_fastas",
    meta_file = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    csv_path = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/final_accession_harsh.csv",
    output = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data_clean/harsh",
    glossary = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv",
)