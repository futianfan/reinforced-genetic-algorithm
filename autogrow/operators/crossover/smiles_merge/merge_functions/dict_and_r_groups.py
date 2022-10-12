"""
Dictionary and Dictionary handling functions
"""
import __future__

import copy

import rdkit
from rdkit import Chem

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.crossover.smiles_merge.merge_functions.mapping_class as mapping_class


def handle_dicts_and_select_b_groups(mol_1, mol_2, mcs_mol):
    """
    this takes 3 rdkit.Chem.rdchem.Mol objects 1 for lig_1,lig_2, and the
    common core(mcs_mol). It creates all the necessary dictionaries, mapping,
    and selects the ligands that will be added to make the final molecule.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: rdkit mol for ligand 1
    :param rdkit.Chem.rdchem.Mol mol_2: rdkit mol for ligand 2
    :param rdkit.Chem.rdchem.Mol mcs_mol: rdkit mol for shared common core
        between mol_1 and mol_2

    Returns:
    :returns: list rs_chosen_smiles: smiles strings for the the R groups which
        correspond to the chosen B's returns None if it fails
    """

    # Confirm that mcs_mol can be replaced in mol_1 and mol_2 around 0.8% of
    # the time this function fails so we will filter this 1st
    they_pass = check_replace_mol(mol_1, mol_2, mcs_mol)
    if they_pass is False:
        return None

    r_smiles_dict_1, b_to_r_master_dict_1, b_to_anchor_master_dict_1 = mol_handling_of_fragmenting_labeling_and_indexing(
        mol_1, mcs_mol, 1
    )
    # check that this worked (ie if it failed they will return None)
    if r_smiles_dict_1 is None:
        return None

    r_smiles_dict_2, b_to_r_master_dict_2, b_to_anchor_master_dict_2 = mol_handling_of_fragmenting_labeling_and_indexing(
        mol_2, mcs_mol, 2
    )
    # check that this worked (ie if it failed they will return None)
    if r_smiles_dict_2 is None:
        return None

    # Merge b_to_anchor_master_dict into 1 master dictionary of B_to_anchors.
    # the keys will be all be unique so we can add these dictionaries together
    # without worry of overrighting an entry. We will invert the dict after to
    # get anchors as the keys and the B's as the items. example
    # b_to_anchor_master {'1B1':[10008,10007],'1B2':[10000],'1B3':[10006],
    # '2B3':[10006,10007],'2B2':[10000],'2B1':[10008]}
    b_to_anchor_master = b_to_anchor_master_dict_1
    for i in list(b_to_anchor_master_dict_2.keys()):
        b_to_anchor_master[i] = b_to_anchor_master_dict_2[i]

    # Invert b_dictionary to produce a master I dictionary. example
    # anchor_to_b_master = {10008:['1B1','2B1'],10000:['1B2','2B2']}
    anchor_to_b_master = invert_dictionary(b_to_anchor_master)

    bs_chosen = mapping_class.run_mapping(b_to_anchor_master, anchor_to_b_master)

    # Get the R groups which correspond to the chosen B's
    # ['1R1', '1R5', '2R2']
    rs_chosen = get_rs_chosen_from_bs(
        bs_chosen, b_to_r_master_dict_1, b_to_r_master_dict_2
    )

    # Get the smiles strings for the the R groups which correspond to the
    # chosen B's
    rs_chosen_smiles = get_rs_chosen_smiles(rs_chosen, r_smiles_dict_1, r_smiles_dict_2)

    return rs_chosen_smiles


def mol_handling_of_fragmenting_labeling_and_indexing(mol, mcs_mol, lig_number):
    """
    This takes an rdkit mol for a ligand and 1 for the mcs_mol. It fragments
    the ligand by replacing the MCS. and it determines which anchors are in
    each fragment. These fragments are our R-groups and the assignment of
    anchors. is how we determine which R-group goes where relative to the MCS.

    lig_number  int    is the number of the ligand that is mol
                       ie if mol is mol_1 lig_number = 1
                       ie if mol is mol_2 lig_number = 2


    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit mol (either mol_1 or mol_2)
    :param rdkit.Chem.rdchem.Mol mcs_mol: rdkit mol for shared common core
        between mol_1 and mol_2
    :param int lig_number: an int either 1 or 2 for (mol_1 or mol_2
        respectively)

    Returns:
    :returns: dict r_smiles_dictionary: a dictionary of the R-groups which
        branch off the common core keys are the R-groups; items are the SMILES
        strings of that R-groups returns None if it fails. Example: {'1R1':
        '[10003*][1007N]=[1013O]', '1R2': '[10000*][1011CH2]=[1008O]'}
    :returns: dict b_to_r_master_dict: A dictionary which tracks the R groups
        which belong to a B-group keys are the B-groups; items are the R-groups
        which belong to the B-group. returns None if it fails. Example: {'1B1':
        ['1R2'], '1B2': ['1R1']}
    :returns: dict b_to_anchor_master_dict: A dictionary which tracks the iso
        label of the anchor atoms for B-group. keys are the B-groups; items are
        the iso label of the anchor atoms for B-group returns None if it fails.
        Example:{'1B1': [10000], '1B2': [10003]}
    """

    # Find which MCS Atoms Rs branch from Function to find all neighbors for a
    # set of molecules touching an Isolabeled core
    mcs_touches = get_atoms_touch_mcs(mol)

    # invert dictionary
    lig_r_atoms_touch_mcs = invert_dictionary(mcs_touches)

    # remove the Core atoms from each ligand this gives us the R-groups
    replace_core = r_group_list(mol, mcs_mol)
    if replace_core is None:
        # replace_core failed to handle fragments"
        return None, None, None

    replace_core = replace_core_mol_dummy_atoms(mol, mcs_mol, replace_core)
    if replace_core is None:
        # replace_core failed to handle fragments"
        return None, None, None

    # A single anchor (isotope label) may now be present multiple times in the
    # replace_core_mols as they are fragmented replace_frag_w_anchor_isolabels
    # can return a None if failed so lets check for None before we move on
    if replace_core is None:
        # replace_core failed to handle fragments"
        return None, None, None

    # MAKE NEW MOL FRAGS FROM LABELED replace_core
    mol_frags = Chem.GetMolFrags(replace_core, asMols=True, sanitizeFrags=False)
    list_r_groups = []
    i = 0
    while i < len(mol_frags):
        val = Chem.MolToSmiles(mol_frags[i], isomericSmiles=True)
        list_r_groups.append(val)
        i = i + 1

    # Generate all the R-libraries with full R-groups using the index of its
    # respective Lig r_chain_dictionary is the master dictionary for R-groups
    r_chain_dictionary, r_smiles_dictionary = r_groups_dict(mol_frags, lig_number)

    # r_dict is a secondary dictionary for searching I's in R's. this
    # dictionary is limited to only the R-group and anchor(I).
    r_dict = get_r_dict(r_chain_dictionary, lig_r_atoms_touch_mcs)

    # make inversion of r_dict. keys are the Anchor atom iso_labels while the
    # items are the R-group numbers which are attached to that anchor atom.
    # Example: {10008: ['2R3'], 10000: ['2R2'], 10006: ['2R1'], 10007:
    # ['2R1']}
    i_dict = invert_dictionary(r_dict)

    """
    B-dictionaries:
    Ligmerge will randomly select R-groups to append to a shared common core
    from 2 separate ligands but so anchor atoms in the common core may have
    more than 1 R-group attached to it.

    ie. if an anchor carbon has di-methyls attached to it (which are not part
        of the shared core) are these di-methyls 2 separate R-groups or is the
        contextual chemical environment created by having both methyls
        different from having one alone and thus should be treated as a single
        R-group which just happen to branch. This author would argue that
        context is important here and so this version of Ligmerge treats
        anything attached to an anchor atom in the common core as a singular
        contextual functional group which shall be referred to as a B-groups.
    ie. a B-group consists of 1 or more R-groups which are attached to an
        anchor atom in the shared common core. This makes a significant
        difference in how we select for which pieces are added to build our
        child molecule. Additionally this has significance in the decision
        tree use to build a child molecule. An R/B group can be connected to
        multiple anchor atoms so once we chose a B group we will need to know
        which anchor atoms are affected by that decision. This is something
        handled more in the Mapping class, but this is why the nomenclature
        change from R-groups to B-groups and why the next several steps are
        important.
    make_b_dictionaries (B is the name we gave to R-groups sets)
    """

    b_to_r_master_dict, b_to_anchor_master_dict = make_b_dic(i_dict, r_dict, lig_number)

    return r_smiles_dictionary, b_to_r_master_dict, b_to_anchor_master_dict


def check_replace_mol(mol_1, mol_2, mcs_mol):
    """
    Confirm that mcs_mol can be replaced in mol_1 and mol_2 around 0.8% of the
    time this function fails so we will filter this 1st

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: an rdkit mol
    :param rdkit.Chem.rdchem.Mol mol_2: an rdkit mol
    :param rdkit.Chem.rdchem.Mol mcs_mol: rdkit mol for shared common core
        between mol_1 and mol_2

    Returns:
    :returns: bool True/False: Returns True if it passes for both mol_1 and
        mol_2 returns False if either fails.
    """

    temp = r_group_list(mol_1, mcs_mol)
    if temp is None:
        return False
    temp = r_group_list(mol_2, mcs_mol)
    if temp is None:
        return False
    return True


###########################################################
###########################################################
# HANDLE THE OBTAINING THE R-Groups for a given mol


def r_group_list(mol, core_mol):
    """
    This takes a mol and the common core and finds all the R-groups by
    replacing the atoms in the ligand (which make up the common core) with
    nothing.

    This fragments the ligand and from those fragments we are able to
    determine what our R-groups are. for any common core atom which touched
    the fragment a * will replace that atom in the fragments.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit molecule
    :param rdkit.Chem.rdchem.Mol core_mol: an rdkit molecule for the shared
        common core

    Returns:
    :returns: rdkit.Chem.rdchem.Mol replace_core_mol: an rdkit molecule with
        the common core removed from a ligand fragments the mol which can be used
        to make lists of R-groups
    """

    # This returns all the mol frags for a particular compound against the
    # core molecule
    replace_core_mol = Chem.ReplaceCore(
        mol, core_mol, labelByIndex=True, replaceDummies=True, requireDummyMatch=False
    )

    if len(replace_core_mol.GetAtoms()) == 0:
        # This means that the mol either did not contain the core_mol or the
        # core_mol is the same mol as the mol. ie) if mol_string
        # ="[10000N-]=[10001N+]=[10002N][10003CH]1[10004O][10005CH]([10006CH2][10007OH])[10008CH]([10013OH])[10009CH]([10012OH])[10010CH]1[10011OH]"
        # and core_string
        # ="[10000NH]=[10001N+]=[10002N][10003CH]1[10004O][10005CH]([10006CH2][10007OH])[10008CH]([10013OH])[10009CH]([10012OH])[10010CH]1[10011OH]"
        # the only difference is the H's which means it can be replaced within
        # because its the same mol This is rare but does occur.
        return None

    return replace_core_mol

def replace_core_mol_dummy_atoms(mol, mcs, replace_core_mol):
    """
    This function will replace the dummy atoms (*) with the isotope label from
    the core atoms. example:
        mol = Chem.MolFromSmiles("[10000N-]=[10001N+]=[10002N][10003CH2][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
        mcs = Chem.MolFromSmiles("[10003CH3][10002N]=[10001N+]=[10000NH]")
        replace_core = Chem.MolFromSmiles("[3*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")

        resulting replace_core = '[10003*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]'

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit molecule
    :param rdkit.Chem.rdchem.Mol mcs: an rdkit molecule for the shared common
        core
    :param rdkit.Chem.rdchem.Mol replace_core_mol: the mol with the MCS
        anchors labeled with * and an isotope label of the idx of the core anchor
        atom

    Returns:
    :returns: rdkit.Chem.rdchem.Mol replace_core_mol: an rdkit molecule with
        the common core removed from a ligand fragments the mol which can be used
        to make lists of R-groups. The * atoms will be isotope labeled with the
        isotope label from the core.
    """


    replace_core_mol_original = copy.deepcopy(replace_core_mol)
    anchor_dict = {}
    anchor_to_set_dict = {}
    for atom in replace_core_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            anchor_iso = atom.GetIsotope() + 10000
            neighbors = atom.GetNeighbors()
            tmp = []
            for n_atom in neighbors:
                tmp.append(n_atom.GetIsotope())
            anchor_dict[anchor_iso] = tmp

            anchor_to_set_dict[atom.GetIdx()] = anchor_iso

    for idx in list(anchor_to_set_dict.keys()):

        atom = replace_core_mol.GetAtomWithIdx(idx)
        anchor_iso = anchor_to_set_dict[idx]
        atom.SetIsotope(anchor_iso)

    return replace_core_mol


def r_groups_dict(mol_frags, lig_number_for_multiplier):
    """
    given a set of mol_frags and the ligand_number (ie. 1 for mol_1 and 2 for
    mol_2) this will make dictionaries of all the Rgroup and all the smiles
    for each Rgroup

    Input
    :param rdkit.Chem.rdchem.Mol mol_frags: a rdkit mol containing fragments
    :param int lig_number_for_multiplier: an int either 1 for mol_1 or 2 for
        mol_2, used to make labels which are traceable to the ligand being used

    Returns:
    :returns: dict r_chain_dictionary: a dictionary with the R-groups and the
        anchor atoms they connect to ie) {'1R1':[13,14],'1R2':[21,22],'1R3':[25]}
    :returns: dict r_smiles_dictionary: a dictionary with the R-groups and the
        SMILES strings of those groups ie
        {'1R1':'[1*]:[1013c]([1020H])[1014c]([1019H])[1015c]([1018H])[1016c](:[2*])[1017H]',
        '1R2':'[3*][1024C]([1026H])([1027H])[1023N] = [1022N+] = [1021N-]',
        '1R3':'[4*][1025O][1029H]'}
    """

    num_frags = len(mol_frags)
    i = 0
    r_chain_dictionary = {}
    r_smiles_dictionary = {}
    k = int(lig_number_for_multiplier)
    while i < num_frags:
        frag = mol_frags[i]
        r_list_temp = []
        r_list_smiles = Chem.MolToSmiles(frag, isomericSmiles=True)
        for atoms in frag.GetAtoms():
            iso = atoms.GetIsotope()
            if 3000 > iso > 100:
                r_list_temp.append(iso - (1000 * k))
                atoms.SetIsotope(0)
            if iso > 3000:
                name = "I{}".format(iso - 10000)
                r_list_temp.append(iso)
            lig_num_r_r_num = "{}R{}".format(k, i + 1)
            r_chain_dictionary[lig_num_r_r_num] = r_list_temp
            r_smiles_dictionary[lig_num_r_r_num] = r_list_smiles
        i = i + 1

    return r_chain_dictionary, r_smiles_dictionary


def get_r_dict(r_chain_dict, lig_r_atom_touch_mcs):
    """
    This will take the r_chain_dict and the dict of all the atoms which touch
    the core and return a dict of Rs groups as keys and their nodes as values

    Inputs:
    :param dict r_chain_dict: dict of all the atom isolabels for in an
        R-group. keys are R-groups;  items are iso-labels of atoms in the R-group.
        ie) {'1R1': [3, 4, 5, 6, 7, 8, 9, 10, 11, 10000]}
    :param dict lig_r_atom_touch_mcs: dict of all the atoms which directly
        touch the core and what anchor they touch. keys are atom isolabels of
        atoms touching an anchor; items are iso-labels of anchor atoms. ie) {3:
        [10000]}

    Returns:
    :returns: dict r_s_dict:  dictionary of R-groups and anchor atoms they are
        connected to. keys are R-groups. items are isolabel of anchor atoms. ie)
        {'1R1': [10000]}
    """

    r_s_dict = {}
    for key in list(r_chain_dict.keys()):
        temp_r_list = r_chain_dict[key]
        node_list = []
        for atom in r_chain_dict[key]:
            for key_id in list(lig_r_atom_touch_mcs.keys()):
                if atom == key_id:
                    for x in lig_r_atom_touch_mcs[key_id]:
                        node_list.append(x)
                    r_s_dict[key] = node_list

    return r_s_dict


##########
# Mapping functions and finding neighbors
#########
def get_idx_using_unique_iso(mol, iso_val):
    """
    This function takes a value for an isotope label and finds the atom in a
    mol which has that isotope label. This assumes there is only 1 atom in a
    mol with the same isotope value

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: a molecule whose atom's have unique
        isotope labels
    :param int iso_val:  the isotope value to search by

    Returns:
    :returns: int idx:  the Idx index number of the atom whose isotope label
        is the same as iso_val. Returns None if iso_val not in mol.
    """

    for atom in mol.GetAtoms():
        if atom.GetIsotope() == iso_val:
            idx = atom.GetIdx()
            return idx
    return None


def make_b_dic(i_dictionary, r_dict_num, lig_number):
    """
    This generates the dictionaries for the B-groups. one is to track the
    R-groups which a B-group represents (this is the b_to_r_master_dict). one
    is to track the anchor atoms a B-group branches from (this is the
    b_to_anchor_master_dict).

    Inputs:
    :param dict i_dictionary:dictionary for R groups bound to nodes (aka I's).
        ie) {'10008':[1R1,1R2],'10009':[1R2,1R3]}
    :param dict r_dict_num: dictionary for anchors which are attached to an R
        group. ie) {'1R1':[10008],'1R2':[10008,10009],'1R3':[10009]}
    :param int lig_number: an int either 1 or 2 for (mol_1 or mol_2
        respectively)

    Returns:
    :returns: dict b_to_r_master_dict: key is unique B-name and the R-groups
        it represents. example {'1B1':['1R1'],'1B2':['1R2','1R3','1R4'],'1B3':
        ['1R5']}
    :returns: dict b_to_anchor_master_dict: key is unique B-name and items are
        anchors that B connects to. example
        {'1B1':[10008,10007],'1B2':[10000],'1B3':[10006]}
    """

    k = lig_number
    b_to_r_master_dict = {}
    b_to_anchor_master_dict = {}
    counter = 1
    anchor_list = list(i_dictionary.keys())
    # anchor_list = [10008, 10000, 10006, 10007]

    while len(anchor_list) > 0:
        anchor = anchor_list[0]
        B_key = "{}B{}".format(k, counter)
        temp_r_list = []
        temp_anchor_list = []

        for Rs in i_dictionary[anchor]:
            # example Rs in i_dictionary[anchor]: '1R1')
            temp_r_list.append(Rs)
            r_dict_i = r_dict_num[Rs]
            for I in r_dict_i:
                # example Rs in i_dictionary[anchor]: '1R1')
                temp_anchor_list.append(I)
        # remove any redundancies in the list by list(set(list_of_things))
        temp_anchor_list = list(set(temp_anchor_list))
        temp_r_list = list(set(temp_r_list))

        # make new B-group entry in the dictionaries
        b_to_r_master_dict[B_key] = temp_r_list  # This B-represents these R-groups
        b_to_anchor_master_dict[
            B_key
        ] = temp_anchor_list  # This B connects to these anchor atoms

        counter = counter + 1

        # make a list of atoms to remove in the next iteration if they are in
        # both the temp_anchor_list and anchor_list
        for i in temp_anchor_list:
            if i in anchor_list:
                anchor_list.remove(i)

    # example
    # b_to_r_master_dict:{'1B1':['1R1'],'1B2':['1R2','1R3','1R4'],'1B3':
    # ['1R5']}. example
    # b_to_anchor_master_dict:{'1B1':[10008,10007],'1B2':[10000],'1B3':[10006]}.
    return b_to_r_master_dict, b_to_anchor_master_dict


def invert_dictionary(old_dic):
    """
    This will invert any dictionary so that the keys are the values and the
    values are the keys.

    Inputs:
    :param dict old_dic: a dictionary to invert

    Returns:
    :returns: dict inverted_dic: old_dict dict inverted so the keys are the
        items and the items are the keys
    """

    # inverted_dic = {}
    # for k, v in old_dic.iteritems():
    # keys = inverted_dic.setdefault(v, [])
    # keys.append(k)
    values = set([a for b in list(old_dic.values()) for a in b])
    values = list(values)
    inverted_dic = dict(
        (new_key, [key for key, value in list(old_dic.items()) if new_key in value])
        for new_key in values
    )

    return inverted_dic


def get_atoms_touch_mcs(mol):
    """
    Function to find all neighbors for a set of molecules touching. Isolabeled
    core atoms.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: isolabeled with atoms in the core having
        isotope. labels set as their idx number + 10000 and atoms not shared in
        the common core isotope labels set as:
            for lig_1: atom idx number + 1000
            for lig_1: atom idx number + 2000


    Returns:
    :returns: dict mcs_touches dict:  a dictionary with keys being the isotope
        label of core atoms and the items being the idx's of all non-core atoms
        which touch it. If a core atom touch no non-core atoms it will not be
        added to the dictionary.
    """

    mcs_touches = {}
    all_atoms = mol.GetAtoms()

    for atom in all_atoms:
        # find all atoms in the mol which are also in the Core using iso
        # labels

        iso = atom.GetIsotope()
        if iso > 9999:
            # then its a core atom
            neighbors = atom.GetNeighbors()
            values = []

            for neighbor_atom in neighbors:
                # compile list of the Indexes of all neighbor of the atom
                Idx_neighbor = neighbor_atom.GetIdx()

                # Select for only neighbors which are not in the core using
                # iso
                iso_neighbor_x = neighbor_atom.GetIsotope()
                if iso_neighbor_x < 9999:
                    # Then this is an atom which is not in the core but
                    # touches the core
                    idx_of_neighbor = neighbor_atom.GetIdx()
                    values.append(idx_of_neighbor)
                    mcs_touches[iso] = values

    return mcs_touches


##########
# Handling after B-groups are chosen
##########
def get_rs_chosen_from_bs(bs_chosen, b_to_r_master_dict_1, b_to_r_master_dict_2):
    """
    this function returns a list of R-groups chosen based on the list of
    chosen B's. It requires the b_to_r_master_dict_1 for both ligands to
    function.

    Inputs:
    :param list bs_chosen: A list of the chosen B-groups. ie) ['1B1', 1B2',
        '2B3']
    :param dict b_to_r_master_dict_1: a Dictionary to reference B and R-groups
        from mol_1. keys are names of B-groups; items are R-groups that a B-group
        represents. ie) {'1B1':['1R1'],'1B2':['1R2','1R3','1R4'],'1B3': ['1R5']}
    :param dict b_to_r_master_dict_2: a Dictionary to reference B and R-groups
        from mol_2. keys are names of B-groups; items are R-groups that a B-group
        represents. ie) {'2B1':['2R1'],'2B2':['2R2','2R3','2R4'],'2B3':
        ['2R5','2R6]}

    Returns:
    :returns: list rs_chosen: a list containing all the R-groups represented
        by the chosen B-groups. ie) ['1R1', '1R2', '1R3','1R4', '2R5', '2R6']
    """

    rs_chosen = []
    for B in bs_chosen:
        Rs_for_the_B = []
        lig_number = B[0]
        B_number = B[2]
        if lig_number == str(1):
            for i in b_to_r_master_dict_1[B]:
                Rs_for_the_B.append(i)

        elif lig_number == str(2):
            for i in b_to_r_master_dict_2[B]:
                Rs_for_the_B.append(i)
        for i in Rs_for_the_B:
            rs_chosen.append(i)

    # rs_chosen looks like ['1R1', '1R5', '2R2']
    return rs_chosen


def get_rs_chosen_smiles(rs_chosen, r_smiles_dict_1, r_smiles_dict_2):
    """
    This function returns a list of SMILES strings for every R-group chosen.
    It requires the R_smile_dictionary for both ligands to function.

    Inputs:
    :param list rs_chosen: A list of the chosen R-groups which will be used to
        generate a new mol. ie) ['2R2', '1R1']
    :param dict r_smiles_dict_1: A dictionary which has can find the SMILES
        string for each R-group of Ligand 1. ie) {'1R1':
        '[10006*][1009N]=[1008N+]=[1007N-]'}
    :param dict r_smiles_dict_2: A dictionary which has can find the SMILES
        string for each R-group of Ligand 2. ie) {'2R2': '[10006*][2009OH]',
        '2R1': '[10003*][2007CH2][2008OH]'}


    Returns:
    :returns: list rs_chosen_smiles: A list of all the SMILES string which are
        to be added to make the child ligand. Each SMILES is a sublist.
        ie)[['[10006*][1009N]=[1008N+]=[1007N-]'],['[10006*][2009OH]']]
    """

    rs_chosen_smiles = []
    for R in rs_chosen:
        Rs_for_the_R = []
        lig_number = R[0]
        R_number = R[2]
        if lig_number == str(1):
            Rs_for_the_R.append(r_smiles_dict_1[R])
        elif lig_number == str(2):
            Rs_for_the_R.append(r_smiles_dict_2[R])

        rs_chosen_smiles.append(Rs_for_the_R)

    return rs_chosen_smiles
