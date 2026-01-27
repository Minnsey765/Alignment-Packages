# Overview

Python packages for alignment of sequences

# Development Environment

## Requirements

Python3

## Configuration

1. Create a virtual environment.

    python(3) -m venv venv

2. Activate the virtual environment.

e.g. in powershell: 

    .\venv\Scripts\activate.ps1

e.g in zsh/bash/sh

    .\venv\Scripts\activate


https://docs.python.org/3/library/venv.html 

3. Install dependencies

    pip(3) install -r requirements.txt


# Contributing

## Tests

Tests are written using [PyTest](https://docs.pytest.org/en/stable/index.html)

### Run Tests

In the root foler

    pytest -v


# Usage

Step by Step:

0. into terminal: pip(3) install -r path/to/requirements.txt

1. generates fasta files of sequences corresponding to a gene, species, or orientation

Run file_maker(sort_crit, raw_file, meta_file, csv_path, output)
- sort_crit = dictionary key to generate fasta files by: symbol, orientation, species (sequence & accession)
- raw_file = raw fasta data folder (aligment/raw_fastas)
- meta_file = meta data folder (aligment/code/fasta_info)
- csv_path = path to csv containing accession numbers as str (project/newGenBank.csv)
- output = output directory path as str (aligment/code/fasta_data)

2. sends fasta files generated in step 1 to clustal omega to be aligned and saves nexus files to a folder (aligned_fastas)

Run loop_clustal(fasta_dir, output_dir, interval)
- fasta_dir = folder containing sorted fasta files (fasta_data)
- output_dir = folder containing aligned nexus files (aligned_fastas)
- interval = wait period before sending another request to clustalOmega (not used)

2.5 Manually trim aligned nexus files for each gene (species/orientation) using mesquite
- save trimmed edit as new nexus file in new folder (edited_nexus)
- edit each nexus file by deleting the top and bottom sections mesquite adds to match nexus file format produced by clustalOmega
- edit each nexus file by replacing {top} with "begin data;", the {bottom} with "end;" and adding "NTAX=XY" and "NCHAR=WXYZ" to the dimensions line

3. Concatenate all aligned sequences into a single sequence for each species (gene, or orientation)

Run big_fat_file_maker(raw_fol, csv_path, output, datatype)
- raw_fol = folder containing nexus alignments (aligned_fastas or edited_nexus)
- csv_path = path to csv file containing accession numbers (Project/newGenBank.csv)
- output = folder to output concatenated alignment (usually same as raw_fol - aligned_fastas or edited_nexus)
- datatype = datatype to align by (dna)

4. Method one of verifying sequences from genbank
Run loop_verify(fasta_file, n, output)
- fasta_file = fasta file path of a given gene, e.g. APOB.fasta 
- n = number of entries returned
- output = folder to dump scores for each fasta file (verify_data)
