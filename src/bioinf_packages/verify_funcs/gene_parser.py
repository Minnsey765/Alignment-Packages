#test if gene name or symbol is in a species definition
import re
from symbol_dictionary import symbol_key
def gene_parser(gene_symbol, des):
    #find gene symbol in description
    match = re.search(r'\(([^)]+)\)', des) #look for bracketed strings in description
    gene = match.group(1) if match else None
    #print(gene)
    matches = 0
    if gene in symbol_key:
        matches = matches + 1
    return(matches)

#print(gene_parser("APOB", "Anourosorex squamipes isolate 1 ApoB (ApoB) gene, partial cds"))