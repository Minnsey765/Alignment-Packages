#find correct symbol from incorrect symbol input and overwrite glossary to input missing symbol

import csv 

#write an exception
class GlossaryIncompleteError(Exception):
    pass


def correction_finder(incorrect: str, key: dict, glossary: str):
    #does the symbol exist in the dictionary?
    value = key.get(incorrect)
    #if it doesn't
    if value is None:
        #add new symbol to dictionary
        key[incorrect] = []
        #overwrite old glossary csv
        with open(glossary, "w", newline="") as f:
            w = csv.DictWriter(f, key.keys())
            w.writeheader()
            w.writerow(key)
        print(f"No match found for '{incorrect}'")
        print(f"Writing new entry ('{incorrect}') into '{glossary}'")
        raise GlossaryIncompleteError(
            "Before continuing, glossary must be filled in with new term identified."
        )
    #if it does
    else:
    #return as string
        correct = key[incorrect]
    return correct

#from write_dict_to_csv import my_dict
#print(correction_finder("SNL", my_dict, "dictionary_funcs/mycsv.csv"))