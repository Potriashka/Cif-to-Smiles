from rdkit import Chem
from CifFile import ReadCif

def cif_to_smiles_rdkit(cif_path, include_hydrogens=False):
    """
    Convert a CIF file with explicit bonds into a valid SMILES string using RDKit.
    
    Parameters:
        cif_path (str): Path to the CIF file.
        include_hydrogens (bool): Whether to include H atoms in the output SMILES.
        
    Returns:
        str: Valid SMILES string.
    """
    cf = ReadCif(cif_path)
    block_name = list(cf.keys())[0]
    block = cf[block_name]

    # --- Read atoms ---
    atom_labels = list(block["_atom_site_label"])
    elements = list(block["_atom_site_type_symbol"])
    label_to_index = {label: i for i, label in enumerate(atom_labels)}

    # --- Create RDKit molecule ---
    mol = Chem.RWMol()
    atom_indices = []

    for el in elements:
        if not include_hydrogens and el == "H":
            atom_indices.append(None)  # skip H in indexing
            continue
        atom_idx = mol.AddAtom(Chem.Atom(el))
        atom_indices.append(atom_idx)

    # --- Add bonds from CIF ---
    if "_geom_bond_atom_site_label_1" in block and "_geom_bond_atom_site_label_2" in block:
        bonds_1 = block["_geom_bond_atom_site_label_1"]
        bonds_2 = block["_geom_bond_atom_site_label_2"]

        for a1, a2 in zip(bonds_1, bonds_2):
            i = label_to_index[a1]
            j = label_to_index[a2]
            # Skip bonds involving H if desired
            if atom_indices[i] is None or atom_indices[j] is None:
                continue
            mol.AddBond(atom_indices[i], atom_indices[j], Chem.BondType.SINGLE)

    # --- Sanitize molecule ---
    Chem.SanitizeMol(mol)

    # --- Convert to SMILES ---
    smiles = Chem.MolToSmiles(mol)
    return smiles

# ----------------------------
# Run as script
# ----------------------------

if __name__ == "__main__":
    cif_file = 'is_m5.cif' # replace with your CIF path
    smiles = cif_to_smiles_rdkit(cif_file, include_hydrogens=False)

    with open("output.smiles", "w") as f:
        f.write(smiles + "\n")

    print("Valid SMILES output:", smiles)
