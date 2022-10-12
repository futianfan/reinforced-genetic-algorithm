"""
This script will fragment a .smi
Example Run:
python fragmenter_of_smi_mol.py \
    --smi_file autogrow4/source_compounds/PARPi.smi
"""
import itertools
import copy
import random
import os
import argparse

import rdkit
import rdkit.Chem as Chem
from rdkit.Chem.BRICS import BRICSDecompose
from rdkit import RDLogger

# Turn off warnings
RDLogger.DisableLog("rdApp.*")

import support_scripts.Multiprocess as mp
import support_scripts.mol_object_handling as MOH


def get_atom_w_iso_num(mol, iso_num):
    """
    Find all permutations of bonds to cut on a molecule.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: any rdkit mol
    :param int iso_num: the isotope number to get index

    Returns:
    :returns: int atom_idx: the atom.GetIdx() of the atom with the
        isotope label. If not in atom return None
    """
    for atom in mol.GetAtoms():
        if atom.GetIsotope() == iso_num:
            return atom.GetIdx()
        continue
    return None


#
def label_iso_num_w_idx(mol):
    """
    Find all permutations of bonds to cut on a molecule.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: any rdkit mol

    Returns:
    :returns: rdkit.Chem.rdchem.Mol mol: a rdkit mol
    """
    for atom in mol.GetAtoms():
        atom.SetIsotope(atom.GetIdx())
    return mol


#
def get_rot_bond_permutations_to_cut(mol, c_c_bonds_off=False):
    """
    Find all permutations of bonds to cut on a molecule.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: any rdkit mol

    :param bool c_c_bonds_off: whether to fragment C-C bonds
    Returns:
    :returns: list permutations_of_bonds_to_remove: list of bonds to cut for a mol
    """
    rotatable_bond = Chem.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]")
    rotatable_bonds_set = mol.GetSubstructMatches(rotatable_bond)

    rotatable_bonds_to_frag = []
    for rot_bond in rotatable_bonds_set:
        atom1 = mol.GetAtomWithIdx(rot_bond[0])
        atom2 = mol.GetAtomWithIdx(rot_bond[1])
        atom_isos = [atom1.GetIsotope(), atom2.GetIsotope()]
        bond = mol.GetBondBetweenAtoms(rot_bond[0], rot_bond[1])
        if bond.GetIsAromatic() is True:
            continue
        # Remove any bonds including Hydrogen
        if atom1.GetAtomicNum() == 1 or atom2.GetAtomicNum() == 1:
            continue
        # Remove any C-C single bonds
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
            if c_c_bonds_off is True:
                continue

            rotatable_bonds_to_frag.append(atom_isos)

        rotatable_bonds_to_frag.append(atom_isos)
    permutations_of_bonds_to_remove = []

    for i in range(1, len(rotatable_bonds_to_frag) + 1):
        #  xrange will return the values 1, 2, 3, 4 in this loop
        temp_perm_list = list(itertools.combinations(rotatable_bonds_to_frag, i))
        permutations_of_bonds_to_remove.extend(temp_perm_list)
    return permutations_of_bonds_to_remove


#
def remove_atoms(mol, list_of_idx_to_remove):
    """
    This function removes atoms from an rdkit mol based on
    a provided list. The RemoveAtom function in Rdkit requires
    converting the mol to an more editable version of the rdkit mol
    object (Chem.EditableMol).

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: any rdkit mol
    :param list list_of_idx_to_remove: a list of idx values to remove
                                        from mol
    Returns:
    :returns: rdkit.Chem.rdchem.Mol new_mol: the rdkit mol as input but with
                                            the atoms from the list removed
    """

    if mol is None:
        return None

    try:
        atoms_to_remove = list_of_idx_to_remove
        atoms_to_remove.sort(reverse=True)
    except:
        return None

    try:
        em1 = Chem.EditableMol(mol)
        for atom in atoms_to_remove:
            em1.RemoveAtom(atom)

        new_mol = em1.GetMol()

        return new_mol
    except:
        return None


#
def get_brics_permutations(mol, min_frag_size=3):
    """
    Fragment a mol using BRICS methods.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: any rdkit mol
    :param int min_frag_size: minimum size of fragments to keep
    Returns:
    :returns: list clean_frag_list: list of fragmented SMILES
    """
    res = list(BRICSDecompose(mol, returnMols=True, minFragmentSize=min_frag_size))
    smis = [Chem.MolToSmiles(x, True) for x in res]

    # Get larger pieces
    res = list(
        BRICSDecompose(
            mol, returnMols=True, keepNonLeafNodes=True, minFragmentSize=min_frag_size
        )
    )
    smis.extend([Chem.MolToSmiles(x, True) for x in res])
    clean_frag_list = []
    for x in res:
        list_to_remove = []
        for i in x.GetAtoms():
            if i.GetAtomicNum() == 0:
                list_to_remove.append(i.GetIdx())

        x = remove_atoms(x, list_to_remove)

        for atom in x.GetAtoms():
            atom.SetIsotope(0)

        clean_frag_list.append(Chem.MolToSmiles(x))
    list(set(list(clean_frag_list)))

    return clean_frag_list


#
def remove_bonds(mol, list_of_atomiso_bondsets_to_remove):
    """
    This function removes bond from an rdkit mol based on
    a provided list. This list is a list of sets, with each set containing
    two atoms with the isotope label of that atom. Using Isotopes is to ensure
    that atom Idx dont change.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: any rdkit mol
    :param list list_of_atomiso_bondsets_to_remove: a list of idx values to remove
                                        from mol
    Returns:
    :returns: rdkit.Chem.rdchem.Mol new_mol: the rdkit mol as input but with
                                            the atoms from the list removed
    """
    # None's often end up in a pipeline use of RDKit so we handle this data type as return None
    # instead of raise TypeError
    if mol is None:
        return None

    # If mol is wrong data type (excluding None) raise TypeError
    if type(mol) != rdkit.Chem.rdchem.Mol and type(mol) != rdkit.Chem.rdchem.RWMol:
        printout = "mol is the wrong data type. \n"
        printout = printout + "Input should be a rdkit.Chem.rdchem.Mol\n"
        printout = printout + "Input mol was {} type.".format(type(mol))
        raise TypeError(printout)
    new_mol = copy.deepcopy(mol)
    if len(list_of_atomiso_bondsets_to_remove) == 0:
        return None
    for atomiso_bondsets in list_of_atomiso_bondsets_to_remove:
        if len(atomiso_bondsets) == 0:
            continue
        if len(atomiso_bondsets) != 2:
            printout = "list_of_atomiso_bondsets_to_remove needs to be 2 isolabels for the atoms"
            raise TypeError(printout)

        atom_1_idx = int(get_atom_w_iso_num(new_mol, atomiso_bondsets[0]))
        atom_2_idx = int(get_atom_w_iso_num(new_mol, atomiso_bondsets[1]))

        try:
            new_mol = Chem.FragmentOnBonds(
                new_mol, [atom_1_idx, atom_2_idx], addDummies=False
            )
        except:
            return None

        new_mol = MOH.check_sanitization(new_mol)
        if new_mol is None:
            return None
    new_mol = MOH.check_sanitization(new_mol)
    if new_mol is None:
        return None
    return new_mol


#
def make_list_of_all_unique_frags(fragment_list):
    """
    This function takes a list of all molecules after fragmentation and separates the
    the fragments into individual rdkit mol objects, sanitizes each, removes isotopes
    and converts them into a SMILES string. The SMILES are compiled into a list,
    and then redundant strings are reduced to a single entry.

    It returns a list of all unique sanitized canonical SMILES for every fragment made
    from all permutations of bond breaking.

    Inputs:
    :param list fragment_list: list of fragmented rdkit mols which haven't been separated
        yet

    Returns:
    :returns: list clean_frag_list: List of unique sanitized SMILES strings from all objects
                in fragment_list. Isotope labels are also removed here.
    """
    clean_frag_list = []
    for fragments in fragment_list:
        frags = Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=False)
        for frag in frags:
            frag = MOH.check_sanitization(frag)
            if frag is None:
                continue

            # Remove those under 2 atoms minimum
            list_mol_atoms = frag.GetAtoms()
            if len(list_mol_atoms) < 3:
                continue

            for atom in frag.GetAtoms():
                atom.SetIsotope(0)
            clean_frag_list.append(
                Chem.MolToSmiles(frag, isomericSmiles=True, canonical=True)
            )
        list(set(list(clean_frag_list)))

    return clean_frag_list


#
def make_unique_lig_id(parent_lig_name, current_lig_list):
    """
    This will make a ligand name from the parent name. Keep start names simple.

    Format of names:
        - str(parent_lig_name) + "_Frag_" + str(random_int)

    Inputs:
    :param str parent_lig_name: str of the ligand Id for the parent mol
    :param list current_lig_list: the list of names already taken

    Returns:
    :returns: str unique_lig_id: A unique ID/name for the child ligand.
    """
    if type(parent_lig_name) != str:
        raise Exception("Ligand ID's to seed this must have Unique string IDs")
    parent_lig_name = parent_lig_name.replace(" ", "")
    picked_name = False
    while picked_name is False:
        random_int = random.choice(range(100000, 999999))
        unique_lig_id = str(parent_lig_name) + "_Frag_" + str(random_int)
        if unique_lig_id in current_lig_list:
            continue
        picked_name = True
        break
    return unique_lig_id


#
def make_frag_list_for_one_mol(
    mol_info, frags_per_seed_lig, run_brics, run_frag, c_c_bonds_off=False
):
    """
    This will take a ligand string and ID encased in the list mol_info.
    This will then be fragmented along all non Carbon-carbon rotatable bonds which
    are not aromatic.

    It will make all permutations of all potential bond breaks, reduce to only unique
    fragments and than pick the number of chosen fragments. Then it will create unique ID's
    for each and return a list of lists containing the chosen unique fragments.

    Inputs:
    :param list mol_info: list containing [mol_string, mol_id]
                mol_info[0] = the SMILE string of the parent mol
                mol_info[1] = the Unique ID of the parent mol
    :param bool run_brics: whether to fragment using BRICS method
    :param bool run_frag: whether to fragment all bonds
    :param bool c_c_bonds_off: whether to fragment C-C bonds

    Returns:
    :returns: list final_frag_list: A list of lists containing the chosen unique fragments.
            final_frag_list[0] = [SMILE, mol_id]
    """
    mol_smile = mol_info[0]
    lig_id = mol_info[1]

    mol = Chem.MolFromSmiles(mol_smile, sanitize=False)
    mol = MOH.check_sanitization(mol)
    if mol is None:
        printout = "\nMolecule {} failed to sanitize. \
                    Could not make any fragments from it".format(
            lig_id
        )
        raise Exception(printout)
    mol_smile = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

    mol = label_iso_num_w_idx(mol)
    mol_copy = copy.deepcopy(mol)
    bonds_to_remove_permutations = get_rot_bond_permutations_to_cut(
        mol_copy, c_c_bonds_off
    )

    fragment_list = []
    for bond_set_to_del in bonds_to_remove_permutations:
        mol_copy = copy.deepcopy(mol)
        x = remove_bonds(mol_copy, bond_set_to_del)
        if x is None:
            continue
        fragment_list.append(x)

    clean_frag_list = []
    if run_frag is True:
        clean_frag_list = make_list_of_all_unique_frags(fragment_list)
        clean_frag_list = list(set(clean_frag_list))

    if run_brics is True:
        mol_copy = copy.deepcopy(mol)
        bric_mols = get_brics_permutations(mol_copy, min_frag_size=3)

        clean_frag_list.extend(bric_mols)
        clean_frag_list = list(set(clean_frag_list))

    if len(clean_frag_list) == 0:
        printout = "\nNo fragments were made for {}.\n".format(lig_id)
        print(printout)
        return [[mol_smile, lig_id]]

    # Pick the number of ligands to make
    final_frag_list = [[mol_smile, lig_id]]

    if frags_per_seed_lig == -1:
        printout = "\nFor {}: {} fragmented were made.".format(
            lig_id, len(clean_frag_list)
        )
        print(printout)
        for frag in clean_frag_list:
            unique_lig_id = make_unique_lig_id(lig_id, final_frag_list)
            temp_frag_info = [frag, unique_lig_id]
            final_frag_list.append(temp_frag_info)
    return final_frag_list


#
def get_ligands_from_smi(smi_file):
    """
    Get the ligands from the smi_file

    Inputs:
    :param str smi_file: Path to smiles file

    Returns:
    :returns: list list_of_ligands: A list of lists containing the chosen unique fragments.
            final_frag_list[0] = [SMILE, mol_id]
    """
    list_of_ligands = []
    with open(smi_file, "r") as smiles_file:
        line_counter = 0
        for line in smiles_file:
            line_counter = line_counter + 1
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts separated by 4-spaces
            if len(parts) == 1:
                parts = line.split(
                    "    "
                )  # split line into parts separated by 4-spaces

            if len(parts) == 2 or len(parts) > 2:
                mol_string = parts[0]
                mol_id = parts[1]
                if type(mol_id) != str:
                    print(
                        "Miss Formatted within .SMI. Line number {}".format(
                            str(line_counter)
                        )
                    )
                    continue

                try:
                    mol = Chem.MolFromSmiles(mol_string, sanitize=False)
                except:
                    print(
                        "Miss Formatted within .SMI. Line number {}".format(
                            str(line_counter)
                        )
                    )
                    continue
                mol = MOH.check_sanitization(mol)
                if mol is None:
                    continue

                mol_smile = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                mol_info = [mol_smile, mol_id]
                list_of_ligands.append(mol_info)

            else:
                continue
    print(
        "Was able to import and sanitize {} \
          ligands from the .smi.".format(
            len(list_of_ligands)
        )
    )
    if line_counter != len(list_of_ligands):
        print(
            "\t Failed to sanitize/import \
              {} ligands from the .smi".format(
                line_counter - len(list_of_ligands)
            )
        )
    print("########")

    return list_of_ligands


#
def run_fragmentation_main(vars):
    """
    This runs the fragmenter.

    Inputs:
    :param dict vars: variable with all of the user variables
    """
    print(vars)
    smi_file = vars["smi_file"]
    output_smi_file = vars["output_smi_file"]
    run_brics = vars["run_brics"]
    frags_per_seed_lig = vars["frags_per_seed_lig"]
    run_frag = vars["run_frag"]
    c_c_bonds_off = vars["c_c_bonds_off"]
    number_of_processors = vars["number_of_processors"]

    print("")
    print("STARTING FRAGMENTER")
    print("frags_per_seed_lig: ", frags_per_seed_lig)
    print("smi_file: ", smi_file)
    print("########")
    print("Importing .smi file")
    list_of_ligands = get_ligands_from_smi(smi_file)

    # create a set of jobs to multithread the fragmentation
    job_input = [
        [mol_info, frags_per_seed_lig, run_brics, run_frag, c_c_bonds_off]
        for mol_info in list_of_ligands
    ]
    job_input = [tuple(x) for x in job_input]
    list_of_ligands = None
    output = mp.multi_threading(
        job_input, number_of_processors, make_frag_list_for_one_mol
    )

    print("Finish multithread\n")
    #

    output = [x for x in output if x is not None]
    output = [x for x in output if x is not ""]

    initial_output_reduce = []

    for x in output:
        initial_output_reduce.extend(x)
    output = None

    initial_output_reduce = [x for x in initial_output_reduce if x[0] != ""]
    initial_output_reduce = [x for x in initial_output_reduce if x[1] != ""]

    # Reduce smile redundancies:
    smiles_list = []
    output_reduce = []
    for x in initial_output_reduce:
        if x[0] in smiles_list:
            continue
        output_reduce.append(x)
        smiles_list.append(x[0])

    final_mol_list = []
    master_smile_list = []
    master_id_list = []
    for x in output_reduce:
        temp_smile = x[0]
        temp_id = x[1]
        if temp_smile in master_smile_list:
            continue
        if temp_id in master_id_list:
            continue

        # Append to master lists and final_mol_list
        final_mol_list.append(x)
        master_smile_list.append(temp_smile)
        master_id_list.append(temp_id)

    # convert list of mols to a print statement
    printout = ""
    for x in final_mol_list:
        printout = printout + x[0] + "\t" + x[1] + "\n"

    print("####")
    print("\nSaving list to file")
    with open(output_smi_file, "w") as f:
        f.write(printout)

    print("Number of parent ligands:         {}".format(len(job_input)))
    print(
        "Number of new fragmented ligands: {}".format(
            len(final_mol_list) - len(job_input)
        )
    )
    print("Total number ligs in output file: {}".format(len(final_mol_list)))


def convert_to_bool(val):
    """
    Converts an integer, string, or boolean to the appropriate boolean.

    Inputs:
    :param int|bool|str val: The value to be converted.
    Returns:
    :returns: bool: The equivalent boolean.
    """

    return (
        True
        if val in [True, 1] or (type(val) is str and val.upper() == "TRUE")
        else False
    )


def process_inputs(inputs):
    """
    This will handle processing all parameters.

    Inputs:
    :params dict inputs: dictionary of argparse parameters
    Returns:
    :returns: dict inputs: dictionary of argparse parameters
    """
    # check input smi
    smi_file = os.path.abspath(inputs["smi_file"])
    if os.path.exists(smi_file) is False:
        raise Exception("\n.SMI file not found.\n")
    if os.path.isfile(smi_file) is False:
        raise Exception("\n.SMI file not found.\n")
    inputs["smi_file"] = smi_file

    # check output_smi_file
    if "output_smi_file" in inputs.keys():
        output_smi_file = inputs["output_smi_file"]
    else:
        output_smi_file = None

    if output_smi_file is None:
        output_smi_file = smi_file + "_Fragmented.smi"
    output_smi_file = os.path.abspath(output_smi_file)

    if os.path.exists(os.path.dirname(output_smi_file)) is False:
        raise Exception("\n.directory for output_smi_file not found.\n")
    if os.path.isfile(smi_file) is False:
        raise Exception("\n.SMI file not found.\n")
    inputs["output_smi_file"] = output_smi_file

    if "frags_per_seed_lig" in inputs.keys():
        inputs["frags_per_seed_lig"] = int(inputs["frags_per_seed_lig"])
    else:
        inputs["frags_per_seed_lig"] = -1

    if "run_brics" in inputs.keys():
        inputs["run_brics"] = convert_to_bool(inputs["run_brics"])
    else:
        inputs["run_brics"] = True

    if "run_frag" in inputs.keys():
        inputs["run_frag"] = convert_to_bool(inputs["run_brics"])
    else:
        inputs["run_frag"] = True

    if "c_c_bonds_off" in inputs.keys():
        inputs["c_c_bonds_off"] = convert_to_bool(inputs["run_brics"])
    else:
        inputs["c_c_bonds_off"] = True

    if "number_of_processors" in inputs.keys():
        inputs["number_of_processors"] = int(inputs["run_brics"])
    else:
        inputs["number_of_processors"] = True

    return inputs


# Argument parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument(
    "--smi_file",
    required=True,
    default=None,
    help="Path to tab-delineated .smi file to fragment",
)
PARSER.add_argument(
    "--output_smi_file",
    "-o",
    type=str,
    default=None,
    help="Path to output tab-delineated .smi file of fragments. \
    If not provided it will play a file in the same directory as smi_file \
    titled smi_file + _Fragmented.smi",
)
PARSER.add_argument(
    "--frags_per_seed_lig",
    type=int,
    required=False,
    default=-1,
    help="Number of fragments to create per input SMILES. \
    default is -1 which mean all possible fragments.",
)
PARSER.add_argument(
    "--run_brics",
    type=bool,
    required=False,
    default=True,
    help="Whether to fragment ligands using BRICS fragmentation. This fragments \
    along synthesizable bonds. Default is True.",
)
PARSER.add_argument(
    "--run_frag",
    type=bool,
    required=False,
    default=True,
    help="Whether to fragment ligands over all rotatable bonds. Default is True.",
)
PARSER.add_argument(
    "--c_c_bonds_off",
    type=bool,
    required=False,
    default=True,
    help="Whether to exclude fragmenting carbon-carbon single bonds. Default is True. \
    If True it will ignore fragments on C-C bonds; if False it will fragment.",
)
PARSER.add_argument(
    "--number_of_processors",
    "-p",
    type=int,
    metavar="N",
    default=-1,
    help="Number of processors to use for parallel calculations. \
    Set to -1 for all available CPUs.",
)


ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = process_inputs(ARGS_DICT)

run_fragmentation_main(ARGS_DICT)
print("Fragments located at: {}".format(ARGS_DICT["output_smi_file"]))
print("Finished")
