"""
Run a single crossover event.
"""
import __future__

import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.crossover.smiles_merge.merge_functions.merge_w_core as MWC
import autogrow.operators.crossover.smiles_merge.merge_functions.dict_and_r_groups as DnR
import autogrow.operators.crossover.smiles_merge.merge_functions.alignment_and_breaks as AnB
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH

def process_ligand_new_mol(ligand_new_mol):
    """
    This function processes the ligand_new_mol.
    It either returns the SMILES string of ligand_new_mol (ligand_new_smiles)
    or None if it failed at any step.

    Inputs:
    :param str lig_string_1: smile string for lig 1

    Returns:
    :returns: str ligand_new_smiles: either returns the SMILES
        string of ligand_new_mol or None if it failed at any step.
    """

    ligand_new_mol = MOH.check_sanitization(ligand_new_mol)
    if ligand_new_mol is None:
        return None

    # REMOVE ALL THE ISOTOPES IN THE NEW MOLECULE
    ligand_new_mol_final = MWC.remove_all_isolabels(ligand_new_mol)

    # Remove any fragments incase 1 made it through
    ligand_new_mol_final = MOH.handle_frag_check(ligand_new_mol_final)
    if ligand_new_mol_final is None:
        return None

    # Make sure there are no unassigned atoms which made it through. These are
    # very unlikely but possible
    ligand_new_mol_final = MOH.check_for_unassigned_atom(ligand_new_mol_final)
    if ligand_new_mol_final is None:
        return None

    ligand_new_mol = MOH.check_sanitization(ligand_new_mol_final)
    if ligand_new_mol is None:
        return None

    ligand_new_smiles = Chem.MolToSmiles(ligand_new_mol, isomericSmiles=True)

    return ligand_new_smiles

def run_main_smiles_merge(vars, lig_string_1, lig_string_2):
    """
    This runs the main script for SmileMerge.

    Inputs:
    :param dict vars: User variables which will govern how the programs runs

    :param str lig_string_1: smile string for lig 1
    :param str lig_string_2: smile string for lig 2. example: lig_string_1 =
        "[N-] = [N+] = NCC(O)COc1cccc2ccccc12"; example: lig_string_2 = "C#
        CCOc1ccc2ccccc2c1CO"

    Returns:
    :returns: str ligand_new_smiles: smile string for the child ligand derived
        from lig_1 and lig_2. Returns None if it failed at any point.
    """

    # lig_string_1 = "[N-] = [N+] = NCC(O)COc1cccc2ccccc12"
    # lig_string_2 = "C# CCOc1ccc2ccccc2c1CO"
    # lig_string_1 = "C1 = CC = CC = C1"
    # Sanitize
    lig_smile_1 = Chem.MolFromSmiles(lig_string_1, sanitize=False)
    lig_smile_2 = Chem.MolFromSmiles(lig_string_2, sanitize=False)

    # Sanitize, deprotanate, and reprotanate both molecules
    mol_1 = MOH.check_sanitization(lig_smile_1)
    mol_2 = MOH.check_sanitization(lig_smile_2)
    if mol_1 is None or mol_2 is None:
        return False

    protanate_step = vars["protanate_step"]
    mol_1 = MOH.handleHs(lig_smile_1, protanate_step)
    mol_2 = MOH.handleHs(lig_smile_2, protanate_step)

    # check that handleHs() worked for both molecules if fail move on to next
    # pair of molecules
    if mol_1 is None or mol_2 is None:
        return None

    # make a list of the two rdkit.Chem.rdchem.Mol objects
    mols = [mol_1, mol_2]

    # Use the below mcs_H function for Most Common Substructure searching.
    # This will prevent broken rings.
    mcs_results = rdFMCS.FindMCS(
        mols,
        matchValences=False,
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
        timeout=vars["max_time_mcs_thorough"],
    )

    if mcs_results.canceled is True:
        return None

    # confirm that this meets the minimum number of matching atoms
    if mcs_results.numAtoms < vars["min_atom_match_mcs"]:
        return None

    ### Convert mcs_res from into usable and referable forms
    mcs_mol = Chem.MolFromSmarts(mcs_results.smartsString)

    # handle_mcs_align_labeling_and_cyclicbreaks
    mol_1, mol_2, mcs_mol = AnB.handle_mcs_align_labeling_and_cyclicbreaks(
        mol_1, mol_2, mcs_mol
    )

    # confirm that this meets the minimum number of matching atoms
    if mol_1 is None or mol_2 is None or mcs_mol is None:
        return None
    # This is a very big step. This is where all of the atoms are fragmented
    # to determine what are the R-groups options. R-groups are consolidated
    # into B-groups for each molecule individually, the 2 dictionaries of
    # B-groups are consolidated into a master B-group which is fed into the
    # Mapping class which will randomly select a set of non-clashing B-groups,
    # such that no 2 chosen B's will have connections to the same anchors
    # (anchors are atoms in the common core). It returns a list of
    # smilestrings for the chosen R-groups. It returns None if a step failed.

    rs_chosen_smiles = DnR.handle_dicts_and_select_b_groups(mol_1, mol_2, mcs_mol)
    if rs_chosen_smiles is None:
        return None

    # MERGE ALL CHOSEN R-GROUPS WITH THE CORE
    ligand_new_mol = MWC.merge_smiles_with_core(rs_chosen_smiles, mcs_mol)
    if ligand_new_mol is None:
        return None

    # ligand_new_smiles is either a SMILES string if processing works
    # or None if processing fails
    ligand_new_smiles = process_ligand_new_mol(ligand_new_mol)

    return ligand_new_smiles






