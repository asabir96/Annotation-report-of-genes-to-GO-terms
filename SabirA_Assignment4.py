import sys
import re

def split_terms(filename):
    ''' 
    Split the entire split into individual terms
    Args:
        filename: a file path to a GO terms file
     '''
    terms = []
    try:
        with open(filename) as f:
            contents = f.read()
            split_terms = contents.split("[Term]")
            return split_terms[1:-1]
    except FileNotFoundError:
        print("File Not Found!")
    return terms


def map_protein_to_go(GAF):
    '''
    Builds a dictionary in order to map each individual protein, a key, to a list of its associated GO terms, the values.
    Args: 
        GAF: a file path to a GO annotations (GAF) file
    '''
    dict1 = {}
    try:
        with open(GAF) as f:
            for line in f:
                if not line.startswith("!"):
                    columns = line.split("\t")
                    protein = columns[1]
                    GOID = columns[3] if "GO:" in columns[3] else columns[4]
                    if protein in dict1:
                        dict1[protein].append(GOID)
                    else:
                        dict1[protein] = [GOID]
    except FileNotFoundError:
        print("File Not Found!")
    return(dict1)

def parse_go_term(term):
    '''
    Grab only the GO IDs from the GO Term and return the associated ID and is_a values in a collection
    Args:
        term: a single GO term from the larger GO file
    '''
    is_a = []
    id = None
    if term:
        for line in term.split("\n"):
            linesplit = line.split()
            if "id: GO:" in line:
                id = linesplit[1]
            if "is_a" in line:
                is_a.append(linesplit[1])
    return {id: is_a}

def find_parent_terms(go_id, go_dict):
    '''
    Takes in two parameters and looks for all parents terms of a GO term and returns them as a collection.
    Args:
        go_id: a singular GO ID
        go_dict: a dictionary of GO terms
    '''
    parent_terms = []
    for goid in go_dict.get(go_id, []):
        parent_terms.append(goid)
        parsed_parents = find_parent_terms(goid, parse_go_term(go_dict.get(goid, None)))
        if parsed_parents:
            parent_terms.extend(parsed_parents)
    parent_terms = list(set(parent_terms))
    return parent_terms

def mapping_to_term(terms):
    '''
    Obtain a new dictionary for output
    '''
    mapping_the_terms = {}
    for items in terms:
        protein_parent_term = parse_go_term(items)
        mapping_the_terms[next(iter(protein_parent_term))] = items # next(iter(protein_parent_term))
    return mapping_the_terms

def main():
    '''
    main function for output
    Run of command line: Python3 *NameOfProgram.py* <GO file> <GAF file>
    '''
    if len(sys.argv) > 2:
        input_terms = sys.argv[1]
        input_annotations = sys.argv[2]
        terms = split_terms(input_terms)
        mapped_protein = map_protein_to_go(input_annotations)

        mapped_term = mapping_to_term(terms)
        testing = {'A0A075B6K4': mapped_protein['A0A075B6K4']}
        if len(sys.argv) == 3:
            outfile = "results.tsv"
        else:
            outfile = sys.argv[3]

        with open(outfile, 'w+') as annotated_result:
            for prot, updated_goids in testing.items():
                annotated_result.write("\n"+prot+"\t")
                for new_goid in updated_goids:
                    combine = find_parent_terms(new_goid, parse_go_term(mapped_term.get(new_goid)))
                    proteinstring = " ".join(combine)
                    print("{0} {1}".format(new_goid, proteinstring))
                    annotated_result.write(proteinstring+"\n")
        print(
            find_parent_terms(
                "GO:0002250", parse_go_term(mapped_term['GO:0002250'])))
    else:
        print("Error! Not enough command line arguments")


if __name__=="__main__":
    main()