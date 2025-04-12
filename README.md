# chem_analyzer
create a graph representation of a molecule by providing an sdf file or through the command line, generate molecular fingerprints and check for substructures

## Python Usage
1. clone repository: "git clone git@github.com:msse-chem274A-2024/final-assignment-JXDolan.git" 
2. create environments: "make environment"
3. activate environment: "conda activate chem_analyzer_env"
4. run pytests: "make tests"

5. perform substructure search from terminal constructing by lists: 
        python molecule.py \
        --molecule_atoms "[('C1', 'C'), ('C2', 'C'), ('O1', 'O')]" \
        --molecule_bonds "[('C1', 'C2', 1), ('O1', 'C2', 1)]" \
        --substructure_atoms "[('B1', 'B'), ('O1', 'O')]" \
        --substructure_bonds "[('B1', 'O1', 1)]"

6. perfom substructure search from terminal constructing by sdf:
    python molecule.py\
    --substructure_sdf "sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf" \
    --substructure_sdf_include_h "False"\
    --molecule_sdf "sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf" \
    --molecule_sdf_include_h "True" 

