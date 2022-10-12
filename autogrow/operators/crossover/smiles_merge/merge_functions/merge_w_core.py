"""
This script is handles the merging of two molecules and cleans up
the molecule.
"""
import __future__

import copy
import rdkit
from rdkit import Chem

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH



# HANDLE THE MERGING OF THE R-groups and the MCS And cleanup
# Final steps
def merge_smiles_with_core(rs_chosen_smiles, mcs_mol):
    """
    This function runs most of the ligand merger portion of the script. It
    merges the chosen R-groups (rs_chosen_smiles) with the common core
    (mcs_mol) by bonding anchors in the core to the atoms bound to the
    respective anchor in the R-group fragment.

    Variables used in this function
        anchor_to_connection_dict: dict of anchor to connected atoms and bond
        types.

        example anchor_to_connection_dict = {10008: [1016, rdkit.Chem.rdchem.BondType.AROMATIC],
                                            10007: [1013, rdkit.Chem.rdchem.BondType.AROMATIC]})

        mol_frag, keys are isotope label, value is Idx.

        example mol_frag_iso_to_idx_dict = {10007: 0, 10008: 8, 1013: 1, 1014: 3, 1015: 5, 1016: 7,
                                           1017: 9, 1018: 6, 1019: 4, 1020: 2})

    Inputs:
    :param list rs_chosen_smiles: A list containing the  SMILES strings for
        the chosen R-groups to add
    :param rdkit.Chem.rdchem.Mol mcs_mol: an rdkit molecule representing the
        Most common substructure (MCS) which will be expanded by adding R-groups
        to make the child molecule

    Returns:
    :returns: rdkit.Chem.rdchem.Mol rw_core_merg: The child molecule with the
        added R-groups built onto the mcs_mol. returns None if the process fails
        or if a None-type makes it through.
    """

    # convert to RWMOL class of molecule which are able to add and remove
    # bonds. RWMOL class is the Read and Write-Mol Class in rdkit.
    rw_core_merg = Chem.RWMol(mcs_mol)

    # sanitize the mol_frag
    rw_core_merg = MOH.check_sanitization(rw_core_merg)
    if rw_core_merg is None:
        # ("rw_core_merg failed to be sanitizable (merge_smiles_with_core)")
        return None

    for r_groups in rs_chosen_smiles:

        for frag in r_groups:

            # empty dicts which will only correspond to individual fragments
            # for each R-group. Dicts will be passed into functions within the
            # for loop which should handle the merging of each R-group to the
            # MCS. Note that making mols with the sanitize=True will kill most
            # of these fragmented mols So always Chem.MolFromSmiles(frag,
            # sanitize = False)

            # make a rdkit mol out of the smiles string of the R-group frag
            mol_frag = Chem.MolFromSmiles(frag, sanitize=False)

            # Try to sanitize the mol_frag
            # It often fails but lets try on some
            mol_frag_copy = copy.deepcopy(mol_frag)
            mol_frag = MOH.check_sanitization(mol_frag)
            if mol_frag is None:
                # ("mol_frag failed to be sanitizable
                # (merge_smiles_with_core)"). It often fails so if it didn't
                # work lets just go back to it. Most often it fails when there
                # is a broken ring group...
                mol_frag = mol_frag_copy

            # make dict of anchor to connected atoms and bond types.
            # example anchor_to_connection_dict = {
            #           10008: [1016, rdkit.Chem.rdchem.BondType.AROMATIC],
            #           10007: [1013, rdkit.Chem.rdchem.BondType.AROMATIC]})
            anchor_to_connection_dict = make_anchor_to_bonds_and_type_for_frag(mol_frag)

            # Make Dict of all atoms in mol_frag, keys are isotope label,
            # value is Idx.
            # example mol_frag_iso_to_idx_dict = {10007: 0, 10008: 8, 1013: 1,
            #   1014: 3, 1015: 5, 1016: 7, 1017: 9, 1018: 6, 1019: 4, 1020: 2})
            mol_frag_iso_to_idx_dict = make_dict_all_atoms_iso_to_idx_dict(mol_frag)

            # Make list of atoms idx to remove the anchor atoms from the frag.
            # this is necessary to prevent the redundancy of multiple of the
            # same anchors once the mol_frag and the core are merged.
            anchors_idxs_to_remove = []
            for anchors in list(anchor_to_connection_dict.keys()):
                # anchors are iso numbers of 10,000 or higher
                idx_val = mol_frag_iso_to_idx_dict[anchors]
                anchors_idxs_to_remove.append(idx_val)

            # remove the anchor atoms from mol_frag
            mol_frag = MOH.remove_atoms(mol_frag, anchors_idxs_to_remove)

            # Merge the frag with the core. ie) mol1="CCC" and mol2="CCCCCC"
            # mol3 = Chem.CombineMols(mol1,mol2); mol3 == "CCC.CCCCCC"
            rw_core_merg = Chem.CombineMols(rw_core_merg, mol_frag)

            # convert to RWMOL class of molecule which are able to add and
            # remove bonds
            rw_core_merg = Chem.RWMol(rw_core_merg)

            # make a dictionary of every atom in rw_core_merg with Iso as the
            # key and the Idx as its value
            core_merg_iso_to_idx_dict = make_dict_all_atoms_iso_to_idx_dict(rw_core_merg)

            for anchor_atom_iso in list(anchor_to_connection_dict.keys()):
                # Idx of the anchor in merged core
                idx_for_anchor = core_merg_iso_to_idx_dict[anchor_atom_iso]
                # unpack list of atom to connect idx and bond types
                list_of_atom_idx, list_of_bond_types = unpack_lists_of_atoms_and_bond_type(
                    anchor_to_connection_dict, anchor_atom_iso, core_merg_iso_to_idx_dict
                )

                # USING THE BOND INFO DRAW BONDS TO MAKE NEW CHILD MOL. THIS
                # IS AN ITERATIVE PROCESS FOR EACH R-GROUP.
                for ai, bt in zip(list_of_atom_idx, list_of_bond_types):
                    # ai is the atom Idx
                    # bt is the bondtype
                    try:
                        rw_core_merg.AddBond(idx_for_anchor, ai, bt)
                    except:
                        return None

    return rw_core_merg


def make_anchor_to_bonds_and_type_for_frag(mol_frag):
    """
    Create a dictionary w anchor atoms as keys.

    for each key, the items are broken into lists of lists with the 1st number
    of each as the isotope of the atom bound and the second value as the bond
    type to recreat bonds later to merge..

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_frag: an R-group which was converted into
        a mol

    Returns:
    :returns: dict anchor_to_connection_dict: a dictionary of anchor atom
        isolabels as keys, item is a list of the isolabel of the atom the key is
        bound to and the bond type. ie. anchor_to_connection_dict[10007] =
        [1004,Chem.BondType.AROMATIC]
    """

    anchor_to_connection_dict = {}
    isos_anchors_idxs_to_remove = []

    for atom in mol_frag.GetAtoms():
        if atom.GetIsotope() > 9999:  # if atom is an anchor
            iso_anchor = atom.GetIsotope()  # isotopes of the anchor atom
            anchor_idx = atom.GetIdx()  # get that atoms idx
            isos_anchors_idxs_to_remove.append(
                anchor_idx
            )  # append to list to remove later

            # empty lists for subloop
            connection_iso_idx_list = (
                []
            )  # list of isotope number for atoms connected to an anchor
            bond_type_list = (
                []
            )  # list of bond types in the same order as connection_iso_idx_list

            neighbor = (
                atom.GetNeighbors()
            )  # all neighbor atoms to the Atom from above loop
            for x in neighbor:  # Atoms which are neighbors of anchor
                iso_neighbor_atom = x.GetIsotope()
                neighbor_bond_idx = x.GetIdx()
                connection_iso_idx_list.append(iso_neighbor_atom)

                # get bond type between anchor and connected atoms
                bond_object = mol_frag.GetBondBetweenAtoms(
                    anchor_idx, neighbor_bond_idx
                )
                bond_type = bond_object.GetBondType()
                bond_type_list.append(bond_type)

            for i, j in zip(connection_iso_idx_list, bond_type_list):
                list_of_atom_and_bond = [i, j]
                if iso_anchor in list(anchor_to_connection_dict.keys()):
                    tmp = anchor_to_connection_dict[iso_anchor]
                    tmp.append(list_of_atom_and_bond)
                    anchor_to_connection_dict[iso_anchor] = tmp
                else:
                    anchor_to_connection_dict[iso_anchor] = [list_of_atom_and_bond]

    return anchor_to_connection_dict


def make_dict_all_atoms_iso_to_idx_dict(mol):
    """
    Make a dictionary of every atom in a molecule with Iso as the key and the
    Idx as its value.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit molecule

    Return
    :returns: dict mol_iso_to_idx_dict: a dictionary of the iso-label of every
        atom in the mol as the keys and the idx of that atom in the mol object.
        ie) {1008: 7, 1009: 8, 1003: 4, 1004: 3, 1010: 9, 1006: 5, 1007: 6, 10000:
        0, 10001: 1, 10002: 2, 1005: 10}
    """

    mol_iso_to_idx_dict = {}
    for atom in mol.GetAtoms():
        iso = atom.GetIsotope()
        idx = atom.GetIdx()
        mol_iso_to_idx_dict[iso] = idx
    return mol_iso_to_idx_dict


def unpack_lists_of_atoms_and_bond_type(anchor_to_connection_dict, anchor_atom_iso,
                                        core_merg_iso_to_idx_dict):
    """
    Iterate through all atoms which will be bound to the anchor and unpackage
    all the bond types in a list.

    Inputs:
    :param dict anchor_to_connection_dict: a dictionary of anchor isotope
        labels as keys and a lists as the items. these lists have 2 variables, the
        1st is the atom iso-label of the atom connected to an anchor, and the
        second variable is the rdkit bond type. ie) {10004: [1007,
        rdkit.Chem.rdchem.BondType.SINGLE]}
    :param int anchor_atom_iso: the interger of an anchor atom's isotope
        label. ie) 10004
    :param dict core_merg_iso_to_idx_dict: a dictionary of atom's isotope
        labels as keys and their corresponding Idx as the items. ie) {1008: 14,
        1014: 11, 1009: 15, 1010: 7, 1007: 13, 10000: 0, 10001: 1, 10002: 2,
        10003: 3, 10004: 4}

    Returns:
    :returns: list list_of_atom_idx: a list containing the atom idx. ie) [13]
    :returns: list list_of_bond_types: a list containing the bond types of
        bonds connected to an anchor atom. ie) [rdkit.Chem.rdchem.BondType.SINGLE]
    """

    list_of_atom_idx = []
    list_of_bond_types = []
    # this if statement determines if there are multiple connections to the
    # anchor and how to unpack the dictionary. if there are 2 connections it
    # will look like the len is 2 (with 2 sub lists) which is the same as if
    # there is there is only 1 connection with two parts in the list. so if
    # type(anchor_to_connection_dict[atom_iso][0]) = int then theres only 1
    # connection. if type(anchor_to_connection_dict[atom_iso][0]) = list then
    # there are multiple lists of lists.
    connection_list = anchor_to_connection_dict[anchor_atom_iso]

    if type(connection_list[0]) == int:
        # get the atom iso label
        atom_iso = anchor_to_connection_dict[anchor_atom_iso][0]

        # get the atoms idx in the merged core
        atom_idx = core_merg_iso_to_idx_dict[atom_iso]

        # get bond type
        bond_type = anchor_to_connection_dict[anchor_atom_iso][1]

        # append to the lists
        list_of_atom_idx.append(atom_idx)
        list_of_bond_types.append(bond_type)

    elif type(connection_list[0]) == list:
        if type(connection_list[0][0]) == int:

            for i in range(0, len(connection_list)):
                # get the atom iso label
                atom_iso = anchor_to_connection_dict[anchor_atom_iso][i][0]

                # get the atoms idx in the merged core
                atom_idx = core_merg_iso_to_idx_dict[atom_iso]

                # get bond type
                bond_type = anchor_to_connection_dict[anchor_atom_iso][i][1]

                # append to the lists
                list_of_atom_idx.append(atom_idx)
                list_of_bond_types.append(bond_type)

    return list_of_atom_idx, list_of_bond_types


def remove_all_isolabels(rw_core_merg):
    """
    Remove all the isotope labels from a molecule. One of the finalizing steps
    used when LigSmiles is nearly complete.

    Inputs:
    :param rdkit.Chem.rdchem.RWMol rw_core_merg: a read-write rdkit molecule
        of the child molecule (after R-groups have been added) to have the
        isotopes to be removed

    Returns:
    :returns: rdkit.Chem.rdchem.RWMol rw_core_merg: the read-write rdkit
        molecule of the child molecule with all the isotope labels removed
    """

    # None's often end up in a pipeline use of RDKit so we handle this data
    # type as return None instead of raise TypeError
    if rw_core_merg is None:
        return None

    # If mol is wrong data type (excluding None) raise TypeError
    if (
            type(rw_core_merg) != rdkit.Chem.rdchem.Mol
            and type(rw_core_merg) != rdkit.Chem.rdchem.RWMol
    ):
        printout = "rw_core_merg is the wrong data type. \n"
        printout = (
            printout
            + "Input should be a rdkit.Chem.rdchem.Mol or rdkit.Chem.rdchem.RWMol\n"
        )
        printout = printout + "Input mol was {} type.".format(type(rw_core_merg))
        raise TypeError(printout)

    for atom in rw_core_merg.GetAtoms():
        atom.SetIsotope(0)

    return rw_core_merg
