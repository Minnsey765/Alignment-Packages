#test if species name exists in a hit definition

def species_parser(species_name, des):
    genus,species = species_name.split("_") #split name at underscore
    genus_species = f"{genus} {species}" #concatenate genus species name with a space instead of underscore
    matches = 0
    #test and assign points
    if genus in des:
        matches = matches + 1
    if genus_species in des:
        matches = matches + 1
    return(matches)

#print(species_parser("Anourosorex_yamashinai", "Anourosorex squamipes isolate 1 ApoB (ApoB) gene, partial cds"))