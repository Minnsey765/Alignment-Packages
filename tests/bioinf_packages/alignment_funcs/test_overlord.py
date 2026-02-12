from bioinf_packages.alignment_funcs.overlord import overlord_function
import os



def test_overlord_returnsExpected():
    csv = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/tests/bioinf_packages/alignment_funcs/test.csv"
    meta = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/tests/bioinf_packages/alignment_funcs/test_meta_data"
    raw = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/tests/bioinf_packages/alignment_funcs/test_raw_data"
    glossary = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/tests/bioinf_packages/alignment_funcs/test_glossary.csv"
    dict = overlord_function(raw, meta, csv, glossary)
    print(dict[0])
    #assert correct no. entries in dictionary
    assert len(dict) == 23
    assert dict[0]['symbol'] == "APOB"
    assert dict[0]['orientation'] == "+"
    assert dict[0]['species'] == "Anourosorex_squamipes"
    assert dict[0]['accession'] == "GU981106"
    assert False
