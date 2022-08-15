obo2csv go.obo is_a.csv name.csv alt_id.csv
    convert obo format ontology definition to tab-delimited text

Input:
    go.obo     - obo format ontology definition

Output:
    is_a.csv   - GO_ID Aspect is_a_direct is_a_indirect
    name.csv   - GO_ID Aspect name
    alt_id.csv - alt_id GO_ID

*************************************************

backpropagate input.txt is_a.csv alt_id.csv output.txt
    back-propagate parent GO terms for the list of input GO terms

Input:
    input.txt  - list of input GO terms. GO terms separated by comma in
                 the same line will be back-propagate together. GO terms
                 at differnet lines will be be back-propagated separately.
    is_a.csv   - GO_ID Aspect is_a_direct is_a_indirect
    alt_id.csv - alt_id GO_ID

Output:
    output.txt - list of input GO terms and their parent terms
