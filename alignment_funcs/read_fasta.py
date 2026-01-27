#obselete function
import os
def read_fasta(directory):
    sequences = {}

    for filename in os.listdir(directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            file_path = os.path.join(directory, filename)
            seq_id = None
            seq_list = []
            with open(file_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):  # Header line
                        if seq_id:  # Save previous sequence
                            sequences[seq_id] = "".join(seq_list)
                        end = line.index(' ')
                        seq_id = line[1:end]  # Remove ">"
                        seq_list = []
                    else:
                        seq_list.append(line)  # Add sequence line

                # Add last sequence
                if seq_id:
                    sequences[seq_id] = "".join(seq_list)

    return sequences

print(read_fasta('C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/raw_fastas')["NC_024563.1"])
