"""This handles alignment and ring breaks for Crossover"""
import __future__

import random
import copy

import rdkit
from rdkit import Chem

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH




# handle_mcs_alignments_labeling_and_cyclicbreaks
def handle_mcs_align_labeling_and_cyclicbreaks(mol_1, mol_2, mcs_mol):
    """
    This will take 2 aligned molecules and pick a specific alignment, check
    for any ring and cyclic breaks and fragmentation removing any atoms from
    the common core which caused those issues, renumber and isotope atoms in
    mol_1, mol_2, mcs_mol to be tractable, and consistent amongst the three.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: an rdkit molecule
    :param rdkit.Chem.rdchem.Mol mol_2: an rdkit molecule
    :param Chem.MolFromSmarts(mcs_res_SMART) mcs_mol: a result object from
        rdkits MCS function

    Returns:
    :returns: rdkit.Chem.rdchem.Mol mol_1: an rdkit molecule isolabeled and
        renumbered or None if it fails
    :returns: rdkit.Chem.rdchem.Mol mol_2: an rdkit molecule isolabeled  and
        renumbered or None if it fails
    :returns: mcs_res_Mol mcs_mol: an MCS result isolabeled and no breaks or
        None if it fails
    """

    if type(mol_1) is not rdkit.Chem.rdchem.Mol:
        return None, None, None
    if type(mol_2) is not rdkit.Chem.rdchem.Mol:
        return None, None, None
    if type(mcs_mol) is not rdkit.Chem.rdchem.Mol:
        return None, None, None

    # Set Alignment, Isotope label, and Handle breaks
    picked_alignment = pick_mcs_alignment(mol_1, mol_2, mcs_mol)
    if picked_alignment is None:
        return None, None, None

    # Isotope label the MCS core
    index_tuple = add_mcs_isolabels(mol_1, mol_2, mcs_mol, picked_alignment)

    ## Reorder all atoms within ligands to be consistent
    mol_1 = renumber_to_mcs(mol_1, picked_alignment[0])
    mol_2 = renumber_to_mcs(mol_2, picked_alignment[1])

    #### CHECK FOR RING BREAKS AND FRAGMENTATION BEFORE CONTINUING FORWARD
    # SYSTEM ERRORS MAY OCCUR HERE WHEN MCS IS BAD
    # Check cyclic ring breaks, if break remove atoms from common core
    mcs_mol, new_index, are_there_breaks = check_cyclic_breaks(
        index_tuple, mol_1, mol_2, mcs_mol
    )

    # prevents None type errors
    if mcs_mol is None:
        return None, None, None

    # Check cyclic ring breaks, if break remove atoms from common core
    if are_there_breaks is True:
        # Run cyclic ring break again, if there are still breaks create sys
        # error and abort merge
        mcs_m3, new_index_2, are_there_breaks_2 = check_cyclic_breaks(
            new_index, mol_1, mol_2, mcs_mol
        )
        if are_there_breaks_2 is True:
            return None, None, None

    # prevents None type errors
    if mcs_mol is None:
        return None, None, None

    # Check for any additional fragmentation issues
    confirm_no_frag = Chem.GetMolFrags(mcs_mol, asMols=True, sanitizeFrags=False)
    if len(confirm_no_frag) != 1:
        return None, None, None

    # Renumber all atoms within both ligands to be consistent with the core
    mol_1 = renumber_to_mcs(mol_1, new_index[0])
    mol_2 = renumber_to_mcs(mol_2, new_index[1])

    # Add correct isolabels to the core atoms. requires 1 final fake alignment
    # tuple of range(0, len(mcs_mol.GetNumAtoms()))
    final_alignment = (new_index[2], new_index[2], new_index[2])
    # Isotope label the core
    index_tuple = add_mcs_isolabels(mol_1, mol_2, mcs_mol, final_alignment)

    # label the 1st atom in R-group with its idx + (1000 for lig1 and 2000 for
    # lig 2) the 1000 vs 2000 label will be used to deconvelute later. these
    # will be the iso values assigned to the 1st atom in an R-group.
    add_r_atom_isolabels(mol_1, mol_2)

    return mol_1, mol_2, mcs_mol


##########################
# CYCLE BREAKS AND REMOVE BAD ATOMS FUNCTIONS
##########################
# USE THESE FUNCTION TO TEST FOR CYCLIC BREYCLE BREAK AND REMOVE BAD ATOMS
# FUNCTIONS
##########################

def check_cyclic_breaks(alignment_tuple, mol_1, mol_2, core):
    """
    Check for cyclic breaks and fixes them.

    Fixing cyclic breaks is handled by removing the atoms from the core (MCS
    substructure). Removing atoms will often cause fragmentation (ie. removing
    a carbon from a methyl will leave 3 fragmented hydrogens). Fragmented
    atoms need to be removed. The largest fragment is assumed to be the core
    of interest which we don't want to remove.

    Inputs:
    :param tuple alignment_tuple: a tuple with the atoms with IDX which match
        in the same order ie. alignment_tuple[0][0] is the IDX for the atom in
        mol_1 which is matched to alignment_tuple[1][0] and alignment_tuple[2][0]
        (for mol_2) alignment_tuple[3] (for core) is the reference index which
        ranges from 0 to the number of atoms in the order (0,1,2,3,4...n)
    :param rdkit.Chem.rdchem.Mol mol_1: rdkit mol for ligand 1
    :param rdkit.Chem.rdchem.Mol mol_2: rdkit mol for ligand 2
    :param rdkit.Chem.rdchem.Mol core: rdkit mol for shared common core
        between mol_1 and mol_2

    Returns:
    :returns: rdkit.Chem.rdchem.Mol core: the original core rdkit mol returned
        if did_a_ring_break is False
    :returns: tuple alignment_tuple: the unaltered input param alignment_tuple
        returned if did_a_ring_break is False
    :returns: bool did_a_ring_break: True if the ring broke and was fixed;
        False if there were no breaks and required no modifications to be made to
        the alignment_tuple or core
    :returns: rdkit.Chem.rdchem.Mol new_core: the modified core rdkit mol
        returned if did_a_ring_break is True
    :returns: tuple new_align_tuple: the modified alignment_tuple returned if
        did_a_ring_break is True
    :returns: bool None: returns 3 Nones if it failed to fix the cyclic breaks
    """

    if type(mol_1) is not rdkit.Chem.rdchem.Mol:
        return None, None, None
    if type(mol_2) is not rdkit.Chem.rdchem.Mol:
        return None, None, None
    if type(core) is not rdkit.Chem.rdchem.Mol:
        return None, None, None

    mcs_ringbreak_idx = []
    mol_1_ringbreak_idx = []
    mol_2_ringbreak_idx = []
    for l1, l2, c1 in zip(alignment_tuple[0], alignment_tuple[1], alignment_tuple[2]):
        atom1 = mol_1.GetAtomWithIdx(l1)
        atom2 = mol_2.GetAtomWithIdx(l2)
        atom_c = core.GetAtomWithIdx(c1)

        # ring breaks can occur when an atom in either lig is a ring atom
        # but the common substructure has that as a non-ring atom
        if atom_c.IsInRing() is False and (
                atom1.IsInRing() is True or atom2.IsInRing() is True
        ):
            mcs_ringbreak_idx.append(l1)
            mol_1_ringbreak_idx.append(l2)
            mol_2_ringbreak_idx.append(c1)

    if len(mcs_ringbreak_idx) > 0:
        new_align_list_l1 = []
        new_align_list_l2 = []
        new_align_list_c1 = []

        # THIS IS A BIT COMPLEX HERE SO THIS IS THE IDEA: Using the list of
        # ringbreak idx's we will delete those ringbreak atoms from the core.
        # This will change the idx's within the core which is why we've
        # iso-labeled the core in previous steps. Once we delete from the core
        # we need to detemine the idx's of fragmented atoms but these will be
        # the idx's in the original core, as determined by the isotope labels.
        # ie. idx_of_atom_in_original_core = atom.GetIsotope()-10000 (We keep
        # the largest fragment as the future core).

        # After we have a list of all the atom idx's (in the original core)
        # for all the atoms which are ringbreaks and all the atoms that will
        # cause cyclicbreaks/fragmentation we will then delete all those atoms
        # at once, thus preventing issues of the idx's changing as we delete
        # and issues of fragmentation.

        # Make a copy of the core to test fragments
        temp_core_removed_breaks = copy.deepcopy(core)

        # delete atoms causing ringbreaks
        temp_core_removed_breaks = MOH.remove_atoms(
            temp_core_removed_breaks, mcs_ringbreak_idx
        )
        if temp_core_removed_breaks is None:
            return None, None, None
        # check for fragmentation and add the smallest fragments indexes to
        # the to delete list
        all_atoms_to_delete = ringbreak_frag_handling(
            temp_core_removed_breaks, mcs_ringbreak_idx
        )
        if all_atoms_to_delete is None:
            return None, None, None

        # Now work on the original core. THIS WILL BE THE OFFICIAL NEW CORE.
        # delete any cyclic breaks or anything connected to a cyclic break
        # which would fragment only delete from core mol
        new_core = MOH.remove_atoms(core, all_atoms_to_delete)
        if new_core is None:
            return None, None, None
        new_core = MOH.check_sanitization(new_core)
        if new_core is None:
            return None, None, None

        # now that we've made a new core, the idx's are different so we need
        # to relabel mol_1 and mol_2

        # remove the Iso-labels from lig 1 and 2 for anything deleted
        remove_iso_labels(mol_1, all_atoms_to_delete)
        remove_iso_labels(mol_2, all_atoms_to_delete)

        # make a new index series for comparing mol_1 and mol_2 to the core.
        # this is done using the original indexing and the atoms which were
        # removed from mcs.
        count = 0
        for l1, l2, c1 in zip(
                alignment_tuple[0], alignment_tuple[1], alignment_tuple[2]
        ):
            if c1 not in all_atoms_to_delete:
                new_align_list_l1.append(l1)
                new_align_list_l2.append(l2)
                new_align_list_c1.append(count)
                count = count + 1
        new_align_tuple = (new_align_list_l1, new_align_list_l2, new_align_list_c1)
        did_a_ring_break = True
        return new_core, new_align_tuple, did_a_ring_break

    # len(mcs_ringbreak_idx) less than or equal to 0
    did_a_ring_break = False
    return core, alignment_tuple, did_a_ring_break


def ringbreak_frag_handling(new_core, mcs_ringbreak_idx):
    """
    This takes a rdkit mol of the core (after atoms were removed to resolve
    cyclic breaks). It then tests and handles ringbreaks and fragmentation.

    if an atom in the core was deleted, if it had additional atoms attached to
    it it can cause fragmenation. so this step handles that and takes the
    largest fragment as the new common core.

    Inputs:
    :param rdkit.Chem.rdchem.Mol new_core: common core mol object
    :param list mcs_ringbreak_idx: list of the idx's of the common core for
        iso labels and later adjustment iso labels and idx numbers

    Returns:
    :returns: list iso_core_frag_list: list of the idx's of common core; same
        as mcs_ringbreak_idx unless there was fragmentation that needed to be
        handled.
    """

    ##########################
    # Check for fragmentation, if core is fragmented than additional
    # processing is required to find the largest frag and to then reassign the
    # indexes in the ligs and core

    check_fragmentation = Chem.GetMolFrags(new_core, asMols=True, sanitizeFrags=False)
    num_frag_len = len(check_fragmentation)
    iso_core_frag_list = copy.deepcopy(mcs_ringbreak_idx)

    if num_frag_len > 1:
        # determine the largest fragment in the list of frags
        largest_frag, largest_frag_index_num = find_biggest_frag(check_fragmentation)

        # the core is now the largest fragment
        core = check_fragmentation[largest_frag_index_num]

        # make a list without the largest fragment
        list_frag_mols = []
        list_of_frag_idxs = range(0, len(check_fragmentation))

        for i in list_of_frag_idxs:
            if i == largest_frag_index_num:
                continue

            frag = check_fragmentation[int(i)]
            list_frag_mols.append(frag)

        # get the idx for all atoms in all frags EXCEPT THE LARGEST FRAG.
        # these will be the idx's of the original common core, before deleting
        # things which will be identified using the Isolabels we added before.
        # We will be deleting these atoms shortly
        for frag in list_frag_mols:

            # get all atom idx's (for the original unaltered common_core) in
            # the frag based on the Iso-labels we added before
            for atoms in frag.GetAtoms():
                index_val = atoms.GetIsotope() - 10000
                iso_core_frag_list.append(index_val)

        # Remove redundancy
        iso_core_frag_list = list(set(iso_core_frag_list))

        return iso_core_frag_list

    # if no fragmentation occured
    return iso_core_frag_list


#
def find_biggest_frag(frag_mols_obj):
    """
    This will take a frag mol object and return the largest fragment and the
    index in frag_mols_obj.

    Inputs:
    :param tuple frag_mols_obj: A tuple containing all the fragments of an
        rdkit mol.

    Returns:
    :returns: rdkit.Chem.rdchem.Mol frag_mols_obj: The largest rdkit mol obj
        in the provided tuple
    :returns: int idx_of_max: the idx number of the largest rdkit mol obj in
        the provided tuple.
    """

    if len(frag_mols_obj) > 1:

        idx_of_max = None
        num_atoms_max = None

        for i in range(0, len(frag_mols_obj)):
            frag = frag_mols_obj[i]
            atom_count = frag.GetNumAtoms()
            if num_atoms_max is None:
                idx_of_max = i
                num_atoms_max = atom_count
            elif num_atoms_max < atom_count:
                idx_of_max = i
                num_atoms_max = atom_count
            else:
                continue

        return frag_mols_obj, idx_of_max

    else:
        return frag_mols_obj, 0


def remove_iso_labels(mol, list_of_idx_to_remove):
    """
    Given a molecule and a set of idx numbers this will remove the isotope
    labels for atoms with the idx numbers in the list.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: rdkit mol which needs the iso labeles
        removed
    :param list list_of_idx_to_remove: a list of the idx's which need the
        isolabels removed
    """

    for i in list_of_idx_to_remove:
        atom = mol.GetAtomWithIdx(i)
        atom.SetIsotope(0)


def add_r_atom_isolabels(mol_1, mol_2):
    """
    Label the 1st atom in R-group with its idx + (1000 for lig1 and 2000 for
    lig 2) the 1000 vs 2000 label will be used to deconvelute later. These
    will be the iso values assigned to the 1st atom in an R-group

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: rdkit mol for ligand 1
    :param rdkit.Chem.rdchem.Mol mol_2: rdkit mol for ligand 2
    """

    # isotope label mol_1
    for atom in mol_1.GetAtoms():
        if atom.GetIsotope() < 9999:
            atom_iso = atom.GetIdx() + 1000
            atom.SetIsotope(atom_iso)

    # isotope label mol_2
    for atom in mol_2.GetAtoms():
        if atom.GetIsotope() < 9999:
            atom_iso = atom.GetIdx() + 2000
            atom.SetIsotope(atom_iso)


# alignment and labeling
def pick_mcs_alignment(mol_1, mol_2, common_core):
    """
    This will take the common substructure (aka the common_core) and find
    every match to it within each of the two ligands and produce to list of
    tuples for the atom index (Idx) relative to the indexing of the
    common_core substrure.

    It will then randomly pick combination of an alignment for mol_1 and
    mol_2.

    This should pick the alignments for future numbering.

    This should return picked_alignment (tuple)
        picked_alignment[0] is alignment for mol_1
        picked_alignment[1] is alignment for mol_2
        picked_alignment[2] is the atom idx of the common_core

        ie. picked_alignment[0][0] is the atom IDx for the atom in mol_1 which
        corresponds to the atom in the common substructure index as 0

        picked_alignment[0][1] is the atom IDx for the atom in mol_1 which
        corresponds to the atom in the common substructure index as 1

        picked_alignment[1][0] is the atom IDx for the atom in mol_2 which
        corresponds to the atom in the common substructure index as 0

        picked_alignment[1][1] is the atom IDx for the atom in mol_2 which
        corresponds to the atom in the common substructure index as 1

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: rdkit mol for ligand 1
    :param rdkit.Chem.rdchem.Mol mol_2: rdkit mol for ligand 2
    :param rdkit.Chem.rdchem.Mol common_core: rdkit mol for the shared core
        between mol_1 and mol_2

    Returns:
    :returns: tuple picked_alignment: tuple with 3 sub lists of the atoms IDx
        which correspond to their respective mol object in the order that they
        match the atom in the same index in the sublist for all 3 mols
        (mol_1,mol_2,common_core)
    """

    # Get the substructure match for the MCS within each ligand
    mol_1_match_idx = mol_1.GetSubstructMatches(
        common_core, uniquify=False, maxMatches=10
    )
    mol_2_match_idx = mol_2.GetSubstructMatches(
        common_core, uniquify=False, maxMatches=10
    )

    all_drug_pairings = []
    for mol_1_match in mol_1_match_idx:
        for mol_2_match in mol_2_match_idx:
            all_drug_pairings.append((mol_1_match, mol_2_match))

    # Check that it is a list of tuples. otherwise the random.choice function
    # breaks in python 2.7
    if type(all_drug_pairings) != list:
        return None
    else:
        if len(all_drug_pairings) == 0 or type(all_drug_pairings[0]) != tuple:
            return None
        if len(all_drug_pairings[0]) == 0:
            return None

    # chose an alignment
    alignment_choice = random.choice(all_drug_pairings)

    # add the common_core atom idx to the tuple
    substruc_idx = tuple(list(range(0, len(common_core.GetAtoms()))))
    picked_alignment = (alignment_choice[0], alignment_choice[1], substruc_idx)

    return picked_alignment


###### Index and convert Ligs w isotope labels
def add_mcs_isolabels(mol_1, mol_2, common_core, picked_alignment):
    """
    This will modify every atom in mol_1, mol_2, and the common_core to have
    the same isotope labels based on the index of the common_core atoms.

    Isotope number is set as the index number for the common_core atoms + 10,000

    Input
    :param rdkit.Chem.rdchem.Mol mol_1: an rdkit molecule
    :param rdkit.Chem.rdchem.Mol mol_2: an rdkit molecule
    :param Chem.MolFromSmarts(mcs_res_SMART) common_core: mcs_res_Mol
    :param tupple picked_alignment: a tuple with 3 subtuples for
        mol_1,mol_2,common_core. The numbers within the sublist is the atom IDx
        for a given atom in a ligand. The sublist index for each atom in for
        picked_alignment corresponds to the Idx of that atoms match in the
        commmon_core. ie picked_alignment[1][3] = 10; thus the IDx atom in mol_2
        which corresponds to the 3rd atom in the common_core is 10.

    Returns:
    :returns: tuple final_index: tuple with three sublists which are the same
        as picked_alignment[2]).
    """
    i = 0
    index_list = []
    for lig1, lig2, c1 in zip(
            picked_alignment[0], picked_alignment[1], picked_alignment[2]
    ):
        atom1 = mol_1.GetAtomWithIdx(lig1)
        atom2 = mol_2.GetAtomWithIdx(lig2)
        atom_c = common_core.GetAtomWithIdx(c1)

        atom1.SetIsotope(10000 + i)
        atom2.SetIsotope(10000 + i)
        atom_c.SetIsotope(10000 + i)
        index_list.append(i)
        i = i + 1
    final_index = index_list, index_list, index_list
    return final_index


def renumber_to_mcs(mol, tuple_order_list):
    """
    This renumbers the indexes of the atoms in a lig to that of the MCS and
    returns a renumbered atom.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit molecule
    :param tuple tuple_order_list: a tuple with 3 subtuples for
        mol_1,mol_2,common_core. The numbers within the sublist is the atom IDx
        for a given atom in a ligand. The sublist index for each atom in for
        tuple_order_list corresponds to the Idx of that atoms match in the
        commmon_core. ie tuple_order_list[1][3] = 10; thus the IDx atom in mol_2
        which corresponds to the 3rd atom in the common_core is 10.

    :returns:
    :returns: rdkit.Chem.rdchem.Mol mol: the same rdkit molecule but with the
        atom idx's renumbered to be consistent with the common core, as provided
        by the tuple_order_list
    """

    # make a tuple the same length as the mol the beggining needs to be the
    # same as the aligment numbering to the MCS and there can be no redundancy
    # in the tuple.
    full_tuple_order = [x for x in tuple_order_list]
    num = 0
    while num < len(mol.GetAtoms()):
        if num in full_tuple_order:
            num = num + 1
        else:
            full_tuple_order.append(num)
            num = num + 1

    # reorder the Idx's of the mol to the full_tuple_order
    mol = Chem.RenumberAtoms(mol, full_tuple_order)
    return mol
