#matching a string to a key in a dictionary that aren't complete matches

import difflib

def find_best_key(record_id, dict):
    #normalise id (spaces to underscores)
    target = record_id.replace(" ", "_")
    #try exact matches
    for key in dict.keys():
        if target in key.lower():
            return key
        
    #try fuzzy matching
    keys = list(dict.keys())
    matches = difflib.get_close_matches(target, keys, n=1, cutoff=0.6)
    if matches:
        return matches[0]
    return None