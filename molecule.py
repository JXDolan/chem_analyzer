import networkx as nx
from typing import Optional, Tuple, List
from provided import parse_sdf
import matplotlib.pyplot as plt
import random
import networkx.algorithms.isomorphism as iso
import argparse
import ast
import operator

#inherits from networkx Graph class
class Molecule(nx.Graph):
    
    #constructor 
    #atoms: list[tuples(str,str)], list[(node identifier (str), atomic symbol (str))], water example: [('H1','H'),('O1','O'),('H2','H')]
    #bonds: list[tuples(str,str,int)], list[(node identifier of first atom in bond (str), node identifier of second atom in bond (str), bond order (int) )]bonds example: [('H1','O1',1),('H2','O1',1)] 
    def __init__(self, atoms: list, bonds: list ):
        '''
        constructor for Molecule objects

        parameters:
        -----------
        atoms (list[tuples(str,str)]): list of the atoms in the molecule. 
        The first str in the tuple is the node ID (atomic symbol and node number), the second str in the tuple is the label (atomic symbol)
        example atoms list for water [('H1','H'),('O1','O'),('H2','H')]

        bonds (list[tuple(str,str,int)]): list of bonds in the molecule. 
        The first str in tuple is the node ID for the first atom in the bond, the second str in the tuple is the node ID for the second atom in the bond,
        the int in the tuple represents the bond order. 
        '''
        super().__init__()

        #defensive check that all node_ids are unique 
        node_id_list = []
        symbols_list = []
        for node_id, symbols in atoms:
            node_id_list.append(node_id)
            symbols_list.append(symbols)

        #make node_id_list into a set
        unique_node_ids = set(node_id_list)

        #check that node_id symbol and labels match
        for i in range(len(symbols_list)):
            if symbols_list[i] in node_id_list[i]:
                continue
            else:
                raise ValueError(f"{node_id_list[i]} and {symbols_list[i]} do not match. Expecting matching automic symbols in node id and label ex: ('C1','C')")
            
        #if there are only unique node ids in the atom, and symbols match, create nodes in molecule graph
        if len(unique_node_ids) == len(node_id_list):
            #add nodes to molecule graph
            for node_id, symbol in atoms:
                if ((type(node_id) == str) and (type(symbol) == str)): #defensive checks that node_id and symbol are strings
                    self.add_node(node_id, label = symbol)
                else:
                    raise ValueError('invalid node ID or atomic symbol. atom must be [(str,str)]')
        else: 
            raise ValueError('all node ids must be unique') 
            
        #add bonds
        atom_ids = {id for id, _ in atoms}
        for bond in bonds:
            #defensive check that bonds are being created from existing nodes and are not connecting the same node and tuple is (str,str,int)
            if (bond[0] in atom_ids and bond[1] in atom_ids) and (bond[0] != bond[1]) and (type(bond[0]) == str) and (type(bond[1]) == str) and (type(bond[2] == int)):
             self.add_edge(bond[0], bond[1], bond_order=bond[2], aromatic = "False")
            else:
                #print(f"Invalid bond: {bond}. Must be a tuple(str,str,int)")
                raise ValueError(f"Invalid bond: {bond}. Must be a tuple(str,str,int)")

    
    #constructor using sdf file 
    @classmethod
    def from_sdf(
        cls, filename: str, include_hydrogen: Optional[bool] = False
    ) ->"Molecule":
        '''
        creates a Molecule from a sdf file

        Parameters:
        ----------
        filename (str): name of sdf file
        include_hydrogen (bool): include hydrogens in the graph (True), do not include hydrogens in graph (False)

        Returns:
        ---------
        Molecule constructed from atoms and bonds
        atoms (list[tuple(str,str)]): list of tuples, first tuple is the node ID (atomic symbol+node number), the second tuple is the atomic symbol
        bonds (list[tuple(str,str,int)]): list of tuples, first tuple is first atom node ID, second tuple is the second atom in the bond's node ID, third tuple is the bond order
        '''
        names_and_elements_dict, bonds = parse_sdf(filename, include_hydrogen)
        #get list of elements from names_and_elements_dict
        atoms = list(names_and_elements_dict.items())
        return cls(atoms, bonds)
    

    #create a visual representation of the molecule
    def visualize(self):
        '''
        creates a visual representation coloring the nodes using CPK coloring

        Returns:
        none
            draws image using matplotlib and saves it
        '''
        #cpk color dictionary
        cpk = {
            'H' : 'white',
            'C' : 'black',
            'N' : 'blue',
            'O' : 'red',
            'F' : 'green',
            'Cl': 'green',
            'Br': '#8b0000', #hexidecimal for dark red
            'I' : '#9400D3', #dark violet
            'He': '#00FFFF', #cyan
            'Ne': '#00FFFF',
            'Ar': '#00FFFF',
            'Kr': '#00FFFF',
            'Xe': '#00FFFF',
            'Rn': '#00FFFF',
            'P' : 'orange',
            'S' : 'yellow',
            'B' : '#F5F5DC', #beige
            'Li': '#7F00FF', #violet
            'Na': '#7F00FF',
            'K' : '#7F00FF',
            'Rb': '#7F00FF',
            'Cs': '#7F00FF',
            'Fr': '#7F00FF',
            'Ti': '#808080', #grey
            'Fe': '#FF8C00'
                    }

        #position = nx.spring_layout(self, k = 0.05)
        position = nx.kamada_kawai_layout(self)
        
        labels = nx.get_node_attributes(self, "label")

        edge_labels = nx.get_edge_attributes(self, "bond_order")

        #print(edge_labels)
        
        colors = []
        for nodes in self.nodes():
            if labels[nodes] in cpk:
                colors.append(cpk[labels[nodes]])
            else:
                #if not in cpk set to beige
                colors.append('pink')
            pass

        #print(colors) 
        plt.figure(figsize=(8, 8))
        nx.draw_networkx_nodes(self, position, node_size = 1000, node_color = colors)
        nx.draw_networkx_edges(self, position, width = 1)
        nx.draw_networkx_edge_labels(self, position, edge_labels)
        nx.draw_networkx_labels(self, position, labels, font_color= 'orange') 
        plt.savefig("molecule_visualization.png")
 
        # clearing the current plot
        plt.clf()

    #check if the molecule has any aromatic rings
    def find_aromatic_rings(self):
        '''
        detects aromatic rings, returns boolean, True if the molecule has aromatic rings, False if the molecule does not have aromatic rings

        Returns
        ------
        (bool): True if aromatic, False if not aromatic
        cycles_list: bonds involed in aromatic cycle
        '''
        try:
            #is a list of bonds that form the cycle
            cycles_list = list(nx.find_cycle(self))
            #print(f'Ring detected!\nedges in the ring: {cycles_list}')
        except nx.NetworkXNoCycle:
            print("No rings detected in molecule")
            return False

        #check if the ring contains C, N, O, or S
        must_contain = ["C","N","O","S"]
        counter = 0

        for bond in cycles_list:
            atom1_label = self.nodes[bond[0]]['label']
            atom2_label = self.nodes[bond[1]]['label']
            if atom1_label in must_contain or atom2_label in must_contain:
                counter +=1
                
        if counter < 1:
            print("Consider only rings containing C, N, O, and S atoms. C, N, O, and S not present in this ring")
            return False

        #get the bond order for each edge in the cycle
        bond_order_list = []
        for u,v in cycles_list:
            #print(self.get_edge_data(u,v))
            bond_order_list.append(self.get_edge_data(u,v)['bond_order'])

        #print(f'bond orders in ring{bond_order_list}')

        #check if fully conjugated
        for i in range (len(bond_order_list)-1):
            current_bond_order = bond_order_list[i]
            next_bond_order = bond_order_list[i+1]
            if current_bond_order != next_bond_order:
                continue
            else:
                print("molecule is not fully conjugated")
                return False
        
        #print("molecule is fully conguated")

        #calculate number of pi electrons so we can calculate n 
        num_pi_electrons = bond_order_list.count(2) ** 2 #pi = number of double bonds ** 2

        #solve for n
        n = (num_pi_electrons -2) / 4

        #n must be equal to or greater than zero to be considered aromatic as well as pass all above checks
        if n >= 0:
            print('molecule is aromatic!')
            #change aromtic to True in edges/bonds in aromatic cycle
            for u,v in cycles_list:
                self[u][v]['aromatic'] = True
            return True

    
    def make_fingerprint(self, fingerprint_size = 2048, path_length =7):
        '''
        makes a fingerprint for the molecule

        parameters:
        -----------
        fingerprint_size (int): number of bits in the fingerprint. Default = 2048 bits, this is what rdkit uses
        path_length (int): the maximum path length. Default value = 7, standard for rdkit

        Returns:
        ---------
        fingerprint (list[bool]): list of booleans
        '''
        #defensive checks
        if path_length > 7:
            raise ValueError('Path length values must be between 1 and 7')
        
        if fingerprint_size < 512:
            raise ValueError('fingerprint_size must be greater than 512')
        
        num_bits_per_hash = 2
        fingerprint = [0] * fingerprint_size 
        self.find_aromatic_rings() #check for aromatic rings so edge attribute can be set to aromatic : True 

        #for all atoms in graph/molecule
        for atom in self.nodes():
            #print(f'atom: {atom}')
            for target in self.nodes():
                for path in nx.all_simple_paths(self, source = atom, target = target, cutoff = path_length): #find paths between atom and all other atoms
                    bond_info = []
                    if len(path) > 1: #if there is more than one atom in the path, check if the bonds between atoms are aromatic
                        for i in range(len(path)-1):
                            aromaticity = self[path[i]][path[i+1]]['aromatic']
                            bond_order = self[path[i]][path[i+1]]['bond_order']
                            bond_info.append(f'{aromaticity}{bond_order}')
                        path_labels = [self.nodes[atom]['label'] for atom in path]
                        path_str = ' '.join(path_labels)
                        arom_str = ' '.join(bond_info)
                        path_str = path_str + " " + arom_str #make a path string with aromatic information
                        #print(path_str)
                    else:
                        path_labels = [self.nodes[atom]['label'] for atom in path]
                        path_str = ' '.join(path_labels) #list is not hashable, must make path into something that is hashable like a string
                        #print(path_str)
                    seed = hash(path_str)
                    #print(f'seed: {seed}')
                    random.seed(seed)
                    random_ints = [random.randint(0,fingerprint_size -1) for i in range(num_bits_per_hash)]
                    #print(f'random_ints: {random_ints}')
                    for rInt in random_ints:
                        index = rInt % fingerprint_size
                        fingerprint[index] = True
        return fingerprint
    
    def __eq__(self,other):
        '''
        checks if two molecules are the same

        parameters:
        ----------
        self (Molecule): molecule 1
        other (Molecule): molecule 2

        Returns:
        ---------
        bool: True when molecule 1 and molecule 2 are the same, False when molecule 1 and molecule 2 are not the same
        '''
        em = iso.generic_edge_match(['bond_order', 'aromatic'], [1, False], [operator.eq,operator.eq])
        nm = iso.categorical_node_match('label', "MISSING SYMBOL" )
        return nx.is_isomorphic(self,other, node_match = nm, edge_match = em)

    def substructure_screen(self, substructure):
        '''
        checks if self Molecule contains substructure
        parameters:
        -----------
        self (Molecule): Molecule that we are checking if it contains substructure
        substructure (Moleucle): Molecule we are checking if it exists within the self 

        Retruns:
        --------
        (bool): False if substructure is not in self, True is substructure is in self
        '''
        em = iso.generic_edge_match(['bond_order', 'aromatic'], [1, False], [operator.eq,operator.eq])
        nm = iso.categorical_node_match('label', "MISSING SYMBOL" )
        GM = iso.GraphMatcher(self,substructure, node_match = nm, edge_match = em)
        #print(GM.subgraph_is_isomorphic())
        return GM.subgraph_is_isomorphic()

        
        
def main():
    parser = argparse.ArgumentParser(description = "check if a substructure is present in a molecule")
    #creating molecule and substructure using default constructor
    #argument to get atoms if making Molecule using default constructor
    parser.add_argument("--molecule_atoms",
                        type = str,
                        help = "atoms in a molecule is a list of tuples (str,str): H2O atoms example [('H1', 'H'), ('O1', 'O'), ('H2', 'H')]"
                        )
    #argument to get the bonds if making the Molecule using default constructor
    parser.add_argument("--molecule_bonds",
                        type = str,
                        help = "bonds in a molecule is a list of tuples (str,str,int): H2O bonds example [('H1', 'O1', 1), ('O1', 'H2', 1)")
    #add argument to get atoms if making the substructure moleucle with default constructor
    parser.add_argument("--substructure_atoms",
                        type = str,
                        help = "atoms in the substructure molecule is a list of tuples (str,str): H2O atoms example [('H1', 'H'), ('O1', 'O'), ('H2', 'H')]")
    #add argument to get the bonds if making the substructure molecule with default constructor
    parser.add_argument("--substructure_bonds",
                        type = str,
                        help = "bonds in the substructure molecule is a list of tuples (str,str): H2O atoms example [('H1', 'O1', 1), ('O1', 'H2', 1)")
    
    #creating the molecule and substructure by passing in sdf file
    parser.add_argument("--molecule_sdf",
                        type = str,
                        help = "molecule_sdf is a filepath to the sdf you want to create the molecule from. Example sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf ")
    #include hydrogens for molecule
    parser.add_argument("--molecule_sdf_include_h",
                        type =str,
                        help = "molecule_sdf_include_h is a boolean, if you want to include hydrogens in the construction of the Molecule (True), if you do not want to include hydrogen (False)")
    #creating the substructure by passing in sdf file
    parser.add_argument("--substructure_sdf",
                        type = str,
                        help = "substructure_sdf is a file path to the sdf file you want to use to create the substructure from. Example df/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf ")
    parser.add_argument("--substructure_sdf_include_h",
                        type = str,
                        help = "substructure_sdf_include_h is a boolean, if you want to include hydrogens in the construction of the Molecule (True), if you do not want to include hydrogen (False)")
    
    args =parser.parse_args()

    #creating molecule through sdf or lists
    if args.molecule_bonds and args.molecule_atoms:
        mol = Molecule(ast.literal_eval(args.molecule_atoms),ast.literal_eval(args.molecule_bonds))
    else:
        include_h = args.molecule_sdf_include_h == "True"
        mol = Molecule.from_sdf(args.molecule_sdf, include_h)

    #creating substructure through sdf or lists
    if args.substructure_bonds and args.substructure_atoms:
        substructure = Molecule(ast.literal_eval(args.substructure_atoms), ast.literal_eval(args.substructure_bonds))
    else:
        include_h = args.substructure_sdf_include_h == "True"
        substructure = Molecule.from_sdf(args.substructure_sdf, include_h)

    result = mol.substructure_screen(substructure)

    if result == True:
        print('substructure found in molecule')
    else:
        print('substructure is not present in molecule')

    
if __name__ == "__main__":
    main()
    '''
    #testing that argparse is working when creating molecule and substructure from lists (should be False)
    python molecule.py \
    --molecule_atoms "[('C1', 'C'), ('C2', 'C'), ('O1', 'O')]" \
    --molecule_bonds "[('C1', 'C2', 1), ('O1', 'C2', 1)]" \
    --substructure_atoms "[('B1', 'B'), ('O1', 'O')]" \
    --substructure_bonds "[('B1', 'O1', 1)]"
    '''       

    '''
    #testing argparse with sdf molecule and string substructure, substructure = benzene (should be true)
    python molecule.py \
    --molecule_sdf "sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf" \
    --molecule_sdf_include_h "False" \
    --substructure_atoms "[('C1', 'C'), ('C2', 'C'), ('C3', 'C'),('C4', 'C'), ('C5', 'C'), ('C6', 'C')]" \
    --substructure_bonds "[('C1', 'C2', 2), ('C2', 'C3', 1), ('C3', 'C4', 2),('C4', 'C5', 1), ('C5', 'C6', 2), ('C6', 'C1', 1)]"

    '''    

    '''
    #testing argparse with list molecule and sdf substructure (should be false)
    python molecule.py \
    --molecule_atoms "[('C1', 'C'), ('C2', 'C'), ('C3', 'C'),('C4', 'C'), ('C5', 'C'), ('C6', 'C')]" \
    --molecule_bonds "[('C1', 'C2', 2), ('C2', 'C3', 1), ('C3', 'C4', 2),('C4', 'C5', 1), ('C5', 'C6', 2), ('C6', 'C1', 1)]" \
    --substructure_sdf "sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf" \
    --substructure_sdf_include_h "False"
    ''' 

    '''
    #testing argpase with sdf molecule and sdf substructure (should be true)
    python molecule.py\
    --substructure_sdf "sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf" \
    --substructure_sdf_include_h "False"\
    --molecule_sdf "sdf/AAQOQKQBGPPFNS-UHFFFAOYSA-N.sdf" \
    --molecule_sdf_include_h "True"
    '''           

            
