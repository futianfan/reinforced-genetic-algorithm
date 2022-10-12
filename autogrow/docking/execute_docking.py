"""
This script handles the docking and file conversion for docking.
"""
import __future__
import os

from autogrow.docking.docking_class.get_child_class import get_all_subclasses

from autogrow.docking.docking_class.docking_class_children import *
from autogrow.docking.docking_class.parent_dock_class import ParentDocking
# from autogrow.docking.docking_class.docking_class_children \
#                           import VinaDocking, QuickVina2Docking

from autogrow.docking.docking_class.docking_file_conversion import *
from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter
# from autogrow.docking.docking_class.docking_file_conversion \
#                           import convert_with_obabel, convert_with_mgltools


def pick_docking_class_dict(dock_choice):
    """
    This will retrieve all the names of every child class of the parent class
    ParentDocking

    Inputs:
    :param list dock_choice: List with the User specified docking choices

    Returns:
    :returns: object child_dict[dock_choice]: the class for running the chosen
        docking method
    """

    children = get_all_subclasses(ParentDocking)

    child_dict = {}
    for child in children:
        child_name = child.__name__
        child_dict[child_name] = child

    return child_dict[dock_choice]


def pick_run_conversion_class_dict(conversion_choice):
    """
    This will retrieve all the names of every child class of the parent class
    ParentDocking

    Inputs:
    :param list conversion_choice: List with the User specified docking
        choices

    Returns:
    :returns: object child_dict[conversion_choice]: the class for running the
        chosen docking method
    """

    children = get_all_subclasses(ParentPDBQTConverter)

    child_dict = {}
    for child in children:
        child_name = child.__name__
        child_dict[child_name] = child

    return child_dict[conversion_choice]


def run_docking_common(vars, current_gen_int, current_generation_dir,
                       smile_file_new_gen):
    """
    where is SMILES -> pdb 
    
    A. pdb -> pdbqt

    B. run docking 

    ---------------------

    This section runs the functions common to all Docking programs.

    IF ONE INCORPORATES A NEW DOCKING SOFTWARE, CONFIRM THAT ITS INPUT/OUTPUTS
    CONFORM TO THIS SECTION.

    ############## VERY IMPORTANT SECTION ########################

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param int current_gen_int: the interger of the current generation indexed
        to zero
    :param str current_generation_dir: the current generation directory to
        find the subfolder with pdb files
    :param str smile_file_new_gen: the name of the file containing the
        molecules in the new population

    Returns:
    :returns: str unweighted_ranked_smile_file: the name of the
        unweighted-ranked SMILES with their docking score
    """

    # Get directory string of PDB files for Ligands
    current_generation_pdb_dir = current_generation_dir + "PDBs" + os.sep

    dock_choice = vars["dock_choice"]
    conversion_choice = vars["conversion_choice"]
    receptor = vars["filename_of_receptor"]

    # Use a temp vars dict so you don't put mpi multiprocess info through
    # itself...
    temp_vars = {}
    for key in list(vars.keys()):
        if key == "parallelizer":
            continue
        temp_vars[key] = vars[key]

    file_conversion_class_object = pick_run_conversion_class_dict(conversion_choice)
    file_conversion_class_object = file_conversion_class_object(temp_vars, receptor, test_boot=False)

    dock_class = pick_docking_class_dict(dock_choice)
    docking_object = dock_class(temp_vars, receptor, file_conversion_class_object, test_boot=False)

    if vars["docking_executable"] is None:
        docking_executable = docking_object.get_docking_executable_file(temp_vars)
        vars["docking_executable"] = docking_executable
    ##### vina or Qvina 

    # Find PDB's
    pdbs_in_folder = docking_object.find_pdb_ligands(current_generation_pdb_dir)
    job_input_convert_lig = tuple([tuple([docking_object, pdb]) for pdb in pdbs_in_folder])

    ##############################################################
    ##### part A. ########
    # print("####################")
    # print("Convert Ligand from PDB to PDBQT format")
    smiles_names_failed_to_convert = vars["parallelizer"].run(job_input_convert_lig, lig_convert_multithread)

    ########### print #############
    # deleted_smiles_names_list_convert = [x for x in smiles_names_failed_to_convert if x is not None]
    # deleted_smiles_names_list_convert = list(set(deleted_smiles_names_list_convert))
    # if len(deleted_smiles_names_list_convert) != 0:
    #     print("THE FOLLOWING LIGANDS WHICH FAILED TO CONVERT:")
    #     print(deleted_smiles_names_list_convert)
    # print("####################")


    # Docking the ligands which converted to PDBQT Find PDBQT's
    pdbqts_in_folder = docking_object.find_converted_ligands(current_generation_pdb_dir)

    job_input_dock_lig = tuple([tuple([docking_object, pdbqt]) for pdbqt in pdbqts_in_folder])
    #############################################################
    ###### part B. #######
    # print("####################")
    # print("Docking Begun")
    #############################
    ########### main ############
    #############################
    smiles_names_failed_to_dock = vars["parallelizer"].run(job_input_dock_lig, run_dock_multithread)
    # print("Docking Completed")
    # print("####################")


    ######################################
    ############### print ##############
    deleted_smiles_names_list_dock = [x for x in smiles_names_failed_to_dock if x is not None]
    deleted_smiles_names_list_dock = list(set(deleted_smiles_names_list_dock))
    # print("THE FOLLOWING LIGANDS WHICH FAILED TO DOCK:", deleted_smiles_names_list_dock)
    # print("####################")
    deleted_smiles_names_list = deleted_smiles_names_list_convert + deleted_smiles_names_list_dock
    if len(deleted_smiles_names_list) != 0:
        pass 
        # print("\nTHE FOLLOWING LIGANDS WHERE DELETED FOR FAILURE TO CONVERT OR DOCK:")
        # print(deleted_smiles_names_list)

    ###################################################
    ############## part B2. retrieve the results 
    # print("#################### save results #####################")
    # print("\nBegin Ranking and Saving results")
    unweighted_ranked_smile_file = docking_object.rank_and_save_output_smi(vars, current_generation_dir, current_gen_int, smile_file_new_gen, deleted_smiles_names_list)
    # print("\nCompleted Ranking and Saving results")

    return unweighted_ranked_smile_file


def lig_convert_multithread(docking_object, pdb):
    """
    Run the ligand conversion of a single molecule. If it failed
    failed_smiles_name will be a string of the SMILE which failed to convert
    If it converts failed_smiles_name will be a None.

    Inputs:
    :param object docking_object: the class for running the chosen docking
        method
    :param str pdb: the path to the pdb of a molecule

    Returns:
    :returns: list failed_smiles_name: if the molecule failed to convert to
        final format. (ie. pdbqt conversion fail)
    """

    failed_smiles_name = docking_object.run_ligand_handling_for_docking(pdb)
    return failed_smiles_name


def run_dock_multithread(docking_object, pdb):
    """
    Run the docking of a single molecule.

    Inputs:
    :param object docking_object: the class for running the chosen docking
        method
    :param str pdb: the path to the pdb of a molecule

    Returns:
    :returns: list failed_smiles_names: any smiles which were deleted (ie.
        docking failed)
    """

    # print("Attempt to Dock: ", pdb)
    failed_smiles_names = docking_object.run_dock(pdb)
    # print('------------- run_dock_multithread in execute_docking.py -----------')
    return failed_smiles_names



