from bioinf_packages.alignment_funcs._overlord import overlord_function
import random
import re
import csv


def test_overlord_returnsExpected():
    acc_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/tests/bioinf_packages/alignment_funcs/test.csv"
    meta = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/tests/bioinf_packages/alignment_funcs/test_meta_data"
    raw = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/tests/bioinf_packages/alignment_funcs/test_raw_data"
    glossary = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/tests/bioinf_packages/alignment_funcs/test_glossary.csv"
    dict = overlord_function(raw, meta, acc_csv, glossary)
    #note that dict is a list of dictionaries not a dictionary itself

    #assert correct no. entries in dictionary
    assert len(dict) == 23

    #randomly select an entry
    i = random.randint(0,22)
    print(dict[i])

    #assert that the gene symbol used exists in the glossary
    with open(glossary, newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)
        second_row = rows[1]
        valid_genes = set(second_row) #rm duplicates
    assert dict[i]['symbol'] in valid_genes

    #assert orientation is either + or -
    assert dict[i]['orientation'] in ("+","-")

    #assert species contains a genus and species separated by "_"
    spec_pattern = r"^[A-Z][a-zA-Z]*_[a-zA-Z]+$"
    assert re.fullmatch(spec_pattern, dict[i]['species']), \
        f"Invalid format at index {i}: {dict[i]['species']!r}"
    
    #assert that an accession number exists
    acc_pattern = r"^[A-Z]{2}_?\d+(?:_\d+)?$"
    assert re.fullmatch(acc_pattern, dict[i]['accession']), \
        f"Invalid format at index {i}: {dict[i]['accession']!r}"
    
    #assert sequence exists in dictionary
    seq_pattern = r"^[GATC]+$"
    assert re.fullmatch(seq_pattern, dict[i]['sequence'][i]), \
        f"Invalid format at index {i}: {dict[i]['sequence']!r}"

    #assert False
