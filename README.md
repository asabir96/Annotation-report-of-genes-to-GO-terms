# Annotation-report-of-genes-to-GO-terms
Parse and grab objects from GO (gene ontology) and GAF (gene associated file) to create an create an annotation report of human genes and their list of associating GO terms


split_terms() function accepts a parameter: filename(a file path to a GO terms file)
This function will open the file, split it into individual GO terms, and return the separated terms as a list

map_protein_to_go() function accepts a parameter: filename(a file path to a GO annotations file (the GAF))
This function will open the file, build the mapping relationship between the protein ID (DB Object ID) and its list of associating GO terms

parse_go_term() function will accepts a parameter: term(a single GO term)
This function will parse for:
-ID
-is_a
and grab only the GO IDs from each element and once parsed, return the ID and is_a values

find_parent_terms() function will accept two parameters: go_id(a single GO ID), go_dict(a dictionary of GO terms)
This function will recursively look for the parent terms of the GO term and return them as a collection.

map GO terms to their direct parent GO terms and return results in a TSV output file
