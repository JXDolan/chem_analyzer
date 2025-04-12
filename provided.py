"""
provided functions.
"""

from typing import Optional

def generate_names(elements: list[str]) -> list[str]:
    """
    Generate a list of atom names from a list of elements.

    Parameters
    ----------
    elements
        An unnumbered list of elements.
    
    Returns
    -------
    names 
        A list of atom names with the format "element + number".
    """

    element_count = {x: 0 for x in set(elements)}
    names = []

    for element in elements:
        element_count[element] += 1
        names.append(f"{element}{element_count[element]}")

    return names

def parse_sdf(
    filename: str, include_hydrogen: Optional[bool] = False
) -> tuple[dict, list]:
    """
    Read an sdf file and return the atom and bond info

    Parameters
    ----------
    filename 
        The name of the file to analyze
    include_hydrogen
        Controls whether information about hydrogens is returned.

    Returns
    -------
    names_and_elements
        A dictionary where the keys are the atom names and the values are the
        elements.
    bonds
        A list of tuples where the first two numbers represent the atom indices
        and the third represents the bond order.
    """

    try:
        with open(filename) as f:
            data = [x.strip() for x in f.readlines()]
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Could not find {filename}. Please check the path and try again."
        )

    # read the number of atoms and bonds from the SDF file
    num_atoms = int(data[3].split()[0])
    num_bonds = int(data[3].split()[1])

    # define the end of the atoms and bonds sections
    atom_end = 4 + num_atoms
    bond_end = 4 + num_atoms + num_bonds

    # grab atom and bonds sections from file.
    atoms = data[4:atom_end]
    bonds = data[atom_end:bond_end]

    # check for consistency - does what we retrieved match what is in the top of the file?
    if num_atoms != len(atoms):
        raise ValueError(
            f"Error reading {filename}. The number of atoms does not match the number of atoms in the file."
        )
    if num_bonds != len(bonds):
        raise ValueError(
            f"Error reading {filename}. The number of bonds does not match the number of bonds in the file."
        )

    # extract the elements from the atoms section
    elements = [x.split()[3] for x in atoms]

    # generate unique name for the atoms
    names = generate_names(elements)

    # create a dictionary of atom names and elements  
    names_and_elements = {k: v for k, v in zip(names, elements)}

    # Get list of bonds - format [ (atom1, atom2, bond_order) ...]
    bond_names = [
        (names[int(x.split()[0]) - 1], names[int(x.split()[1]) - 1], int(x.split()[2]))
        for x in bonds
    ]

    # remove hydrogens.
    if not include_hydrogen:
        # remove bonds between hydrogen atoms and other elements
        # since we are using atom names instead of indices, this is
        # a lot cleaner. 
        bond_names = [
            x
            for x in bond_names
            if names_and_elements[x[0]] != "H" and names_and_elements[x[1]] != "H"
        ]

        # remove hydrogen atoms
        names_and_elements = {k: v for k, v in names_and_elements.items() if v != "H"}

    return names_and_elements, bond_names

