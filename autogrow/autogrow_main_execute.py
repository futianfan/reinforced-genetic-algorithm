"""
Top level for running AutoGrow.
Runs all population generation (operations) and docking.
Runs plotting at end.
"""
import __future__

import os
import glob
import sys
import shutil

import autogrow.docking.execute_docking as DockingClass
import autogrow.operators.operations as operations
import autogrow.docking.concatenate_files as concatenate_files
from autogrow.model import Ligand2D

def main_execute(vars):
    """
    This function takes the user variables and runs Autogrow

    Inputs:
    :param dict vars: dict of user variables which will govern how the
        programs runs
    """

    # Unpack necessary variables
    # output_directory is the root output folder for the run
    output_directory = vars["output_directory"]
    num_gens_to_make = vars["num_generations"]

    # Determine what was the last completed generation in the Run directory
    last_generation = determine_current_gen(output_directory)
    if last_generation is None:
        # Check to see if there's a Run 0 based on the seed.
        if vars["use_docked_source_compounds"] is True:
            # This will assess and rank the source compounds prior to
            # generation 1. Thus using the source compounds as a generation 0
            starting_generation_num = 0
        else:
            starting_generation_num = 1

    else:
        starting_generation_num = last_generation + 1

    if starting_generation_num > num_gens_to_make:
        print("This simulation has already been completed to the user defined number \
                of generations. Please check your user variables.")
        raise Exception("This simulation has already been completed to the user defined number \
                of generations. Please check your user variables.")

    ##########################################################################################
    ##########################################################################################
    ### policy network 
    mutate_ligand_select_policy_net = Ligand2D() 
    mutate_reaction_select_policy_net = Ligand2D()
    crossover_ligand1_policy_net = Ligand2D() 
    crossover_ligand2_policy_net = Ligand2D()

    opt1 = torch.optim.Adam(mutate_ligand_select_policy_net.parameters(), lr=1e-3)
    opt2 = torch.optim.Adam(mutate_reaction_select_policy_net.parameters(), lr=1e-3)
    opt3 = torch.optim.Adam(crossover_ligand1_policy_net.parameters(), lr=1e-3)
    opt4 = torch.optim.Adam(crossover_ligand2_policy_net.parameters(), lr=1e-3)
    ##########################################################################################
    ##########################################################################################
    # This is the main loop which will control and execute all commands This
    # is broken into 3 main sections:
    # 1)  operations which populating the new generation with ligands which
    #     both pass the userdefined filter and convert from 1D smiles to 3D
    #     PDB
    # 2)  Docking which handles converting from PDBs to Docking specific
    #     formats and running the actual Docking simulations
    # 3)  Ranking the generation based on the Docking scores
    for current_generation_number in range(starting_generation_num, num_gens_to_make+1):
        sys.stdout.flush()

        # Get directory for smi to go
        current_generation_dir = vars["output_directory"] + "generation_{}{}".format(current_generation_number, os.sep)
        print(current_generation_dir)
        sys.stdout.flush()


        ##### 0-th generation 
        if current_generation_number == 0 and vars["use_docked_source_compounds"] is True:
            if os.path.exists(current_generation_dir + os.sep + "generation_0_ranked.smi") is True:
                continue

            ##################################################
            #################### main  #######################
            ##################################################
            ##### A. populate_generation 
            already_docked, smile_file_new_gen, new_gen_ligands_list = operations.populate_generation_zero(vars, generation_num=0)
            sys.stdout.flush()

            if already_docked is False:
                # Run file conversions of PDB to docking specific file type
                # and Begin Docking unweighted_ranked_smile_file is the file
                # name where the unweighted ranked but score .smi file resides
                ##########################################
                ##### B. docking 
                ##########################################
                unweighted_ranked_smile_file = DockingClass.run_docking_common(
                    vars, current_generation_number,
                    current_generation_dir, smile_file_new_gen)

        else:
            ##################################################
            #################### main  #######################
            ##################################################
            ##### A. populate_generation 
            smile_file_new_gen, new_gen_ligands_list = operations.populate_generation(vars, current_generation_number, 
                                                                                    mutate_ligand_select_policy_net, mutate_reaction_select_policy_net, 
                                                                                    crossover_ligand1_policy_net, crossover_ligand2_policy_net, )
            ## smiles -> sdf -> pdb 
            sys.stdout.flush()
            if new_gen_ligands_list is None:
                raise ValueError("Population failed to make enough mutants or crossovers... \
                                    Errors could include not enough diversity, too few seeds to the generation, \
                                    the seed mols are unable to cross-over due to lack of similarity,\
                                    or all of the seed lack functional groups for performing reactions.")

            # Run file conversions of PDB to docking specific file type and
            # Begin Docking unweighted_ranked_smile_file is the file name
            # where the unweighted ranked but score .smi file resides
            ##########################################
            ##### B. docking 
            ##########################################
            ## pdb -> pdbqt 
            unweighted_ranked_smile_file = DockingClass.run_docking_common(vars, current_generation_number, current_generation_dir, smile_file_new_gen)

        # Delete all temporary files; Skip if in Debugging Mode
        if vars["debug_mode"] is False:
            print("Deleting temporary files and directories")
            files_to_del = []
            folders_to_del = ["{}{}3D_SDFs{}".format(current_generation_dir, os.sep, os.sep), "{}{}3D_SDFs{}log{}".format(current_generation_dir, os.sep, os.sep, os.sep), "{}{}gypsum_submission_files{}".format(current_generation_dir, os.sep, os.sep)]
            for folder in folders_to_del:
                if os.path.exists(folder) is False:
                    continue
                files_to_del.extend(glob.glob(folder+"*"))

            job_input = tuple([tuple([x]) for x in files_to_del if os.path.isfile(x) is True])
            vars["parallelizer"].run(job_input, delete_temporary_files_and_folders)
            # Delete Folders in an ordered manor incase folders are nested
            for i in range(0, len(folders_to_del)):
                delete_temporary_files_and_folders(folders_to_del[i])

        sys.stdout.flush()
        if vars["reduce_files_sizes"] is True:
            # Reduce the files in the PDBs folder to a single compiled file.
            # This reduces the data size And makes it easier to transfer the
            # data
            pdbs_folder = "{}{}PDBs{}".format(current_generation_dir, os.sep, os.sep)
            if os.path.exists(pdbs_folder) is True:
                concatenate_files.run_concatenation(vars["parallelizer"], pdbs_folder)
            else:
                print("\nNo PDB folder to concatenate and compress. This is likely generation 0 seeded with a Ranked .smi file.\n")
        print("")
        print("Finished generation ", current_generation_number)

        sys.stdout.flush()

    if vars["generate_plot"] is True:
        matplotlib_is_callable = False
        try:
            import matplotlib
            matplotlib_is_callable = True
        except:
            matplotlib_is_callable = False
        if matplotlib_is_callable is False:
            print("Can not make figure as matplotlib is not installed")
        else:
            print("Plotting")
            import autogrow.plotting.generate_line_plot as plot
            plot.generate_figures(vars)

    sys.stdout.flush()
#

def determine_current_gen(output_directory):
    """
    Check if there has been any previous runs in the output directory. Returns
    an integer of the last completed generation folder. The last completed
    generation folder will be what seeds the next generation. If no previous
    runs exist which completed (have a ranked.smi file) then it returns a None
    which causes the program to start off at generation 0 using the
    source_compound_file to seed generation 1.

    If the last two generation folders were incomplete (ie both lack a
    ranked.smi file) then we will raise Exception.

    Additionally, if a generation failed to complete in a previous attempt,
    than that generation directory will be renamed so that we can make a new
    generation in its place without losing that data

    -ie if a failed generation directory was named PATH/generation_3 it will
    be rename Path/generation_3_Failed_0

    -if Path/generation_3_Failed_0 already exists it will be name
    Path/generation_3_Failed_1 or so on until unique

    Inputs:
    :param str output_directory: is the path of the Run folder within root
        output folder.

    Returns:
    :returns: int last_gen_number: the int of the last generation number or
        None if no previous generations were completed.
    """

    folder_path_gen = output_directory + "generation_"

    for tries in range(2):
        if tries == 2:
            print("We are in the following directory:", output_directory)

            raise Exception("The last 2 generations in this Run have failed to complete. \
                            Please check that the Run folder that there is something to continue off of.")

        last_gen_number = find_last_generation(folder_path_gen)
        if last_gen_number is None:
            # There are no previous runs in this directory
            return None

        # A previous run exists. The number of the last run.
        folder_path = "{}{}".format(folder_path_gen, last_gen_number)

        is_completed = determine_if_gen_completed(folder_path, last_gen_number)

        if is_completed is True:
            # The last generation (last_gen_number) completed and we will
            # continue our run from that
            return last_gen_number

        # The last generation in the folder crashed before completing.
        # So we will rename the directory by appending _FAILED to the
        # folder name

        printout = "Generation {} in {} failed in the previous simulation.".format(last_gen_number, folder_path)
        print(printout)

        counter = 0
        dir_exists = True
        while dir_exists is True:
            failed_folder_rename = "{}_FAILED".format(folder_path)
            failed_folder_rename_count = "{}_{}".format(failed_folder_rename, counter)

            if os.path.isdir(failed_folder_rename_count) is True:
                counter = counter + 1
            else:
                dir_exists = False

        os.rename(folder_path, failed_folder_rename_count)
        printout = "Renaming folder: {} \
                    to: {}".format(folder_path, failed_folder_rename)
        print(printout)

###################################
### main 
###################################






def find_last_generation(folder_path_string_no_gen):
    """
    This will take a folder path which is missing an interger at the end, and
    find if there are any folders which exist with that path with an interger.
    If there are it will return the highest interger that when added to the
    path exists as a directory.

    If no directories exist with that path+0 then we return None. This causes
    the starting generation of this attempt to run to be generation_0.
    Starting fresh.

    folder_path_string_no_gen = output_directory + "generation_"

    Inputs:
    :param str folder_path_string_no_gen: the folder to check.

    Returns:
    :returns: int last_gen_number: the int of the last generation number or
        None if no previous runs.
    """

    path_exists = True
    i = 1
    while path_exists is True:
        folder_path = "{}{}{}".format(folder_path_string_no_gen, i, os.sep)
        if os.path.exists(folder_path):
            i = i + 1

        else:
            path_exists = False

    if i == 1:
        # Check to see if there's a Run 0 based on the seed.
        i = 0
        folder_path = "{}{}{}".format(folder_path_string_no_gen, i, os.sep)
        if os.path.exists(folder_path) is False:
            return None

        # There are no previous runs in this directory
        last_gen_number = None

    else:
        last_gen_number = i - 1

    return last_gen_number
#

def determine_if_gen_completed(gen_dir_path, gen_number):
    """
    Check if this generation has completed or if it failed. Every generation
    which completes has a .smi file title generation_0_ranked.smi (with the
    number of the generation between the word generation and ranked).
    -If a Run failed due to either a hard crash or a soft crash, there should
        not be a ranked .smi file.

    Inputs:
    :param str gen_dir_path: is the path of the generation folder within a Run
        folder.
    :param int gen_number: The generation number of the folder.

    Returns:
    :returns: bool os.path.isfile(file_path): Returns True if the gen_dir_path
        has a ranked.smi file. Returns False if the gen_dir_path does not have a
        ranked.smi file
    """

    ranked_file_name = "generation_{}_ranked.smi".format(gen_number)
    file_path = "{}{}{}".format(gen_dir_path, os.sep, ranked_file_name)

    return os.path.isfile(file_path)
#

def delete_temporary_files_and_folders(file_or_folder):
    """
    This deletes all temporary files.

    Inputs:
    :param str file_or_folder: the file or folder to delete

    """
    if os.path.exists(file_or_folder) is True:
        if os.path.isdir(file_or_folder) is True:
            try:
                shutil.rmtree(file_or_folder)
            except:
                pass
        else:
            try:
                os.remove(file_or_folder)
            except:
                pass

        # If it failed to delete try via bash command
        if os.path.exists(file_or_folder) is True:
            command = "rm -rf {}".format(file_or_folder)
            try:
                os.system(command)
            except:
                pass
    else:
        pass
#
