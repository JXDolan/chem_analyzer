import pytest
from molecule import Molecule
from provided import parse_sdf
import networkx as nx
import os

#test that default constructor is working
@pytest.mark.parametrize("input_atoms, input_bonds, expected_atoms, expected_bonds",[
    # h2o test
    ([('H1', 'H'), ('O1', 'O'), ('H2', 'H')],
     [('H1', 'O1', 1), ('O1', 'H2', 1)],
     [('H1', {"label": "H"}), ('O1', {"label": "O"}), ('H2', {"label": "H"})],
     [('H1', 'O1', {"bond_order": 1, "aromatic": "False"}), ('O1', 'H2', {"bond_order": 1, "aromatic": "False"})]),
    
    # Ethanol skeletal representation
    ([('C1', 'C'), ('C2', 'C'), ('O1', 'O')],
     [('C1', 'C2', 1), ('O1', 'C2', 1)],
     [('C1', {"label": "C"}), ('C2', {"label": "C"}), ('O1', {"label": "O"})],
     [('C1', 'C2', {"bond_order": 1, "aromatic": "False"}), ('O1', 'C2', {"bond_order": 1, "aromatic": "False"})])
])
def test_constructor(input_atoms, input_bonds, expected_atoms, expected_bonds):
    mol = Molecule(input_atoms, input_bonds)
    
    assert list(mol.nodes(data=True)) == expected_atoms

    # Assert edges (bonds) are added correctly with bond orders
    for u,v,edge_data in expected_bonds:
        assert mol.has_edge(u,v)
        assert mol.get_edge_data(u,v) == edge_data

#test that sdf constructor from_sdf function is working
def test_sdf_constructor():
    test_sdf = Molecule.from_sdf("sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf", include_hydrogen = True)

    #check for expected number of nodes
    assert (len(list(test_sdf.nodes)) == 13)

    #check for the expected number of edges
    assert (len(list(test_sdf.edges)) == 13)

    #check that all expected nodes are present
    assert list(test_sdf.nodes(data=True)) == [
    ('C1', {"label" : "C"}),
    ('C2', {"label" : "C"}),
    ('N1', {"label" : "N"}),
    ('C3', {"label" : "C"}),
    ('C4', {"label" : "C"}),
    ('C5', {"label" : "C"}),
    ('C6', {"label" : "C"}),
    ('C7', {"label" : "C"}),
    ('C8', {"label" : "C"}),
    ('C9', {"label" : "C"}),
    ('C10', {"label" : "C"}),
    ('C11', {"label" : "C"}),
    ('Br1', {"label" : "Br"}),
    ]

    expected_bonds = [
    ('C1', 'C2'),
    ('C2', 'N1'),
    ('N1', 'C3'),
    ('N1', 'C4'),
    ('N1', 'C5'),
    ('C5', 'C6'),
    ('C6', 'C7'),
    ('C7', 'C8'),
    ('C8', 'C9'),
    ('C9', 'C10'),
    ('C10', 'C11'),
    ('C11', 'Br1'),
    ('C11', 'C6'),]

    #check that all expected bonds are in the graph
    for u,v in expected_bonds:
        assert test_sdf.has_edge(u,v)

#test visualization function sdf files
@pytest.mark.parametrize("input_value", ['sdf/AAKJLRGGTJKAMG-UHFFFAOYSA-N.sdf', 'sdf/AAOVKJBEBIDNHE-UHFFFAOYSA-N.sdf', 'sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf'])
def test_visualize(input_value):
    #delete png to make sure new one is being created each time
    if os.path.exists("molecule_visualization.png"):
        os.remove("molecule_visualization.png")
    mol = Molecule.from_sdf(input_value, include_hydrogen= False)
    mol.visualize()
    assert os.path.exists("molecule_visualization.png")  

#test visualize functions from lists 
@pytest.mark.parametrize(
    "atoms, bonds",
    [
        ([('H1', 'H'), ('O1', 'O'), ('H2', 'H')],
         [('H1', 'O1', 1), ('O1', 'H2', 1)]),
        ([('C1', 'C'), ('C2', 'C'), ('C3', 'C'), ('C4', 'C'), ('C5', 'C'), ('C6', 'C')],
         [('C1', 'C2', 2), ('C2', 'C3', 1), ('C3', 'C4', 2), 
          ('C4', 'C5', 1), ('C5', 'C6', 2), ('C6', 'C1', 1)])
])
def test_visualize_with_lists(atoms, bonds):
    """
    Tests the visualize function for molecules built using lists of atoms and bonds.
    """
    #delete png to make sure new one is being created each time
    if os.path.exists("molecule_visualization.png"):
        os.remove("molecule_visualization.png")
    mol = Molecule(atoms, bonds)
    mol.visualize()  
    assert os.path.exists("molecule_visualization.png")            


#testing aromatic function
@pytest.mark.parametrize("input_molecule, expected_bool", [
    ("sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf", True), 
    ("sdf/ZTVQQQVZCWLTDF-UHFFFAOYSA-N.sdf", True),  
    ("sdf/ZSBOMTDTBDDKMP-OAHLLOKOSA-N.sdf", False),
    (Molecule([('H1', 'H'), ('O1', 'O'), ('H2', 'H')], [('H1', 'O1', 1), ('H2', 'O1', 1)]), False),
    (Molecule([('B1', 'B'), ('B2', 'B'), ('B3', 'B'), ('B4', 'B'), ('B5', 'B'), ('B6', 'B')], [('B1', 'B2', 2), ('B2', 'B3', 1), ('B3', 'B4', 2), ('B4', 'B5', 1), ('B5', 'B6', 2), ('B6', 'B1', 1)]), False) 
])
def test_aromatic(input_molecule, expected_bool):
    if isinstance(input_molecule, str):
        mol = Molecule.from_sdf(input_molecule, include_hydrogen= False)
    else: 
        mol = input_molecule
    bool_result = mol.find_aromatic_rings()
    assert bool_result == expected_bool


#test fingerprint function, ensure that hydroxyl and carbonyl have different fingerprints
def test_fingerprint():
    #carbonyl and hydroxyl should be different
    #make hydroxyl
    hydroxyl = Molecule([('C1','C'),('O1','O'),('H1', 'H')],[('C1','O1', 1),('O1','H1',1)])
    hydroxyl_fp = hydroxyl.make_fingerprint(2048, 7)
    #make carbonyl
    carbonyl = Molecule([('C1','C'),('O1','O')],[('C1','O1', 2)])
    carbonyl_fp = carbonyl.make_fingerprint(2048,7)
    #check that they have different fingerprints
    hydroxl_v_carbonyl_same = (carbonyl_fp == hydroxyl_fp)
    assert hydroxl_v_carbonyl_same == False


#test == function
@pytest.mark.parametrize("graph1, graph2, expected_bool", 
    [
        #same molecule but constructed in different order
        (  
        Molecule([('H1', 'H'), ('O1', 'O'), ('H2', 'H')], [('H1', 'O1', 1), ('O1', 'H2', 1)]),
        Molecule([('H2', 'H'),('H1', 'H'), ('O1', 'O')],[('O1', 'H2', 1),('H1', 'O1', 1)]),
        True),
         #same sdf file
        (Molecule.from_sdf("sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf", include_hydrogen = False),
        Molecule.from_sdf("sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf", include_hydrogen = False),
        True),
        #different sdf files
        (
            Molecule.from_sdf("sdf/ZTVQQQVZCWLTDF-UHFFFAOYSA-N.sdf", include_hydrogen = False),
            Molecule.from_sdf("sdf/ZSBOMTDTBDDKMP-OAHLLOKOSA-N.sdf", include_hydrogen = False),
            False
        ),
        #bond_order different
        (  
        Molecule([('H1', 'H'), ('O1', 'O'), ('H2', 'H')], [('H1', 'O1', 1), ('O1', 'H2', 1)]),
        Molecule([('H2', 'H'),('H1', 'H'), ('O1', 'O')],[('O1', 'H2', 2),('H1', 'O1', 3)]),
        False
        ),
        #atoms involved in bonds different
        (  
        Molecule([('H1', 'H'), ('O1', 'O'), ('H2', 'H')], [('H1', 'O1', 1), ('O1', 'H2', 1)]),
        Molecule([('H2', 'H'),('H1', 'H'), ('O1', 'O')],[('O1', 'H2', 2),('H1', 'H2', 3)]),
        False
        ),


    ])
def test_eq(graph1, graph2, expected_bool):
    result = (graph1==graph2)
    assert result == expected_bool


# Test the substructure function
@pytest.mark.parametrize("molecule, substructure, expected_bool", 
    [
        # Water molecule contains a hydroxyl substructure
        (
            Molecule([('H1', 'H'), ('O1', 'O'), ('H2', 'H')], [('H1', 'O1', 1), ('O1', 'H2', 1)]), 
            Molecule([('O2', 'O'), ('H1', 'H')], [('H1', 'O2', 1)]), 
            True
        ),
        # Another test case where substructure is not present
        (
            Molecule([('C1', 'C'), ('O1', 'O')], [('C1', 'O1', 2)]), 
            Molecule([('O2', 'O'), ('H1', 'H')], [('H1', 'O2', 1)]), 
            False
        ),

        #check for benzene
        (
            Molecule.from_sdf("sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf", include_hydrogen = True),
            Molecule([('C1', 'C'), ('C2', 'C'), ('C3', 'C'),('C4', 'C'), ('C5', 'C'), ('C6', 'C')], [('C1', 'C2', 2), ('C2', 'C3', 1), ('C3', 'C4', 2),('C4', 'C5', 1), ('C5', 'C6', 2), ('C6', 'C1', 1)]),
            True
        )
    ])
def test_substructure(molecule, substructure, expected_bool):
    result_bool = molecule.substructure_screen(substructure)
    assert result_bool == expected_bool
    
def test_fingerprint_equivalance():
    #resonance 
    pyridine1 = Molecule([('C1', 'C'), ('C2', 'C'), ('C3', 'C'), ('N1', 'N'), ('C4', 'C'), ('C5', 'C')], [('C1', 'C2', 2), ('C2', 'C3', 1), ('C3', 'N1', 2),('N1', 'C4', 1), ('C4', 'C5', 2), ('C5', 'C1', 1)])
    pyridine2 = Molecule([('C1', 'C'), ('C2', 'C'), ('C3', 'C'), ('N1', 'N'), ('C4', 'C'), ('C5', 'C')], [('C1', 'C2', 1), ('C2', 'C3', 2), ('C3', 'N1', 1),('N1', 'C4', 2), ('C4', 'C5', 1), ('C5', 'C1', 2)])

    p1_fp = pyridine1.make_fingerprint()
    p2_fp = pyridine2.make_fingerprint()

    is_same = (p1_fp == p2_fp)

    assert (is_same == True)
