"""user_vars
This should contain the functions for defining input variables.
Both the default variables and the user input variables.
This should also validate them.
"""

import __future__

import os
import copy
import datetime
import json
import sys
import platform
from shutil import copyfile


def program_info():
    """
    Get the program version number, etc.

    Returns:
    :returns: str program_output: a string for the print of the program information
    """
    program_output = "\nAutoGrow Version 4.0.3\n"
    program_output = program_output + " ================== \n"
    program_output = (
        program_output
        + "If you use AutoGrow 4.0.3 in your research, please cite the following reference:\n"
    )
    program_output = program_output + "Spiegel, J.O., Durrant, J.D. \n"
    program_output = program_output + "AutoGrow4: an open-source genetic algorithm "
    program_output = program_output + "for de novo drug design and lead optimization. \n"
    program_output = program_output + "J Cheminform 12, 25 (2020). \n"
    program_output = program_output + "[doi: 10.1186/s13321-020-00429-4]\n"
    program_output = program_output + " ================== \n\n"

    return program_output


#
def save_vars_as_json(vars):
    """
    This function saves the vars dictionary as a json file. This can be used
    later to track experiments and is necessary for several of the utility
    scripts.
    It saves all variables except the parallelizer class object.

    It saves the file to the output_directory + "vars.json"
        -If AutoGrow has been run multiple times for the same directory it
        will save the new vars file as append a number to the file name
        starting with 2. The util scripts will only look at the original "vars.json"
            ie) output_directory + "vars_2.json"

    Inputs:
    :param dict vars: dict of user variables which will govern how the programs runs
    """
    output_directory = vars["output_directory"]

    vars_file = output_directory + os.sep + "vars.json"
    if os.path.exists(vars_file):
        # vars.json already exist
        # lets make the next file
        path_exists = True
        i = 2
        while path_exists is True:
            vars_file = "{}{}vars_{}.json".format(output_directory, os.sep, i)
            if os.path.exists(vars_file):
                i = i + 1
            else:
                path_exists = False

    temp_vars = {}
    for k in vars.keys():
        if "parallelizer" in k or k == "filter_object_dict":
            continue

        temp_vars[k] = copy.deepcopy(vars[k])

    with open(vars_file, "w") as fp:
        json.dump(temp_vars, fp, indent=4)


def multiprocess_handling(vars):
    """
    This function handles the multiprocessing functions. It establishes a Paralellizer object
    and adds it to the vars dictionary.
    Inputs:
    :param dict vars: dict of user variables which will govern how the programs runs
    Returns:
    :returns: dict vars: dict of user variables which will govern how the programs runs
    """

    # Handle Serial overriding number_of_processors
    # serial fixes it to 1 processor
    if vars["multithread_mode"].lower() == "serial":
        vars["multithread_mode"] = "serial"
        if vars["number_of_processors"] != 1:
            print(
                "Because --multithread_mode was set to serial, "
                + "this will be run on a single processor."
            )
        vars["number_of_processors"] = 1

    # Handle mpi errors if mpi4py isn't installed
    if vars["multithread_mode"].lower() == "mpi":
        vars["multithread_mode"] = "mpi"
        try:
            import mpi4py
        except:
            printout = "mpi4py not installed but --multithread_mode is set to"
            printout = printout + " mpi. \n Either install mpi4py or switch "
            printout = printout + "multithread_mode to multithreading or serial"
            raise ImportError(printout)

        try:
            import func_timeout
            from func_timeout import func_timeout, FunctionTimedOut
        except:
            printout = "func_timeout not installed but --multithread_mode is "
            printout = printout + "set to mpi. \n Either install func_timeout "
            printout = printout + "or switch multithread_mode to"
            printout = printout + " multithreading or serial"
            raise ImportError(printout)

    # # # launch mpi workers
    if vars["multithread_mode"] == "mpi":
        # Avoid EOF error
        from autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import (
            Parallelizer,
        )

        vars["parallelizer"] = Parallelizer(
            vars["multithread_mode"], vars["number_of_processors"]
        )

        if vars["parallelizer"] is None:
            printout = "EOF ERRORS FAILED TO CREATE A PARALLIZER OBJECT"
            print(printout)
            raise Exception(printout)

    else:
        # Lower level mpi (ie making a new Parallelizer within an mpi)
        #   has problems with importing the MPI environment and mpi4py
        #   So we will flag it to skip the MPI mode and just go to multithread/serial
        # This is a saftey precaution
        from autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import Parallelizer

        vars["parallelizer"] = Parallelizer(
            vars["multithread_mode"], vars["number_of_processors"], True
        )

    return vars

def test_docking_executables(vars, vina_exe, qvina2_exe):
    """
    This will test if docking executables are compatible with OS.
    This is only required for MacOS.

    Test will output the version of Vina and QVina2.1 executables to txt file
    in the root_output_folder (docking_exe_MACOS_test.txt)
    If both executables are compatible with this MacOS there should be the following
    2 lines in the txt file:
        AutoDock Vina 1.1.2 (May 11, 2011)
        QuickVina 2.1 (24 Dec, 2017)

    Returns True if both work and returns False.

    Inputs:
    :param dict vars: dict of user variables which will govern how the programs runs
    :param str vina_exe: path to vina executable
    :param str qvina2_exe: path to quick vina 2 executable
    Returns:
    :returns: bool bool: returns True if both docking executables work; False if either fails
    """
    test_vina_outfile = vars["root_output_folder"] + os.sep + "docking_exe_MACOS_test.txt"
    try:
        command = "{} --version > {arg_2} 2>> {arg_2}".format(vina_exe, arg_2=test_vina_outfile)
        os.system(command)
        command = "{} --version >> {arg_2} 2>> {arg_2}".format(qvina2_exe, arg_2=test_vina_outfile)
        os.system(command)
    except:
        printout = "Docking executables could not be found."
        # is not compatible on this OS. \nPlease use docker "
        return False

    with open(test_vina_outfile, "r") as test_file:
        lines = test_file.readlines()
    if "AutoDock Vina 1.1.2" not in lines[0]:
        printout = "Vina docking is not compatible on this OS. \nPlease use docker or "
        printout = printout + "try provide a Vina executable compatible with the OS.\n"
        print(printout)
        if vars["dock_choice"] == "VinaDocking":
            return False

    if "QuickVina 2.1" not in lines[1]:
        printout = "QuickVina 2.1 docking is not compatible on this OS. \nPlease use docker"
        printout = printout + " or try provide a Vina executable compatible with the OS.\n"
        print(printout)
        if vars["dock_choice"] == "QuickVina2Docking":
            return False
    return True

def run_macos_notarization(vars):
    """
    This function runs notarization on vina and qvina2 docking.
    This is important for MacOS newer than 10.15 and newer than

    For MacOS newer than 10.15, this will require an internet connection.

    Inputs:
    :param dict vars: dict of user variables which will govern how the programs runs
    """
    current_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
    vina_exe = current_dir + os.sep.join(["docking", "docking_executables", "vina", \
                                          "autodock_vina_1_1_2_mac", "bin", "vina"])
    qvina2_exe = current_dir + os.sep.join(["docking", "docking_executables", \
                                            "q_vina_2", "q_vina_2_1_mac", "qvina2.1"])

    # Check executables exist. raise exception if not
    if os.path.exists(vina_exe) is False or os.path.exists(qvina2_exe) is False:
        printout = "Docking executables could not be found."
        raise Exception(printout)

    both_docking_exe_work = test_docking_executables(vars, vina_exe, qvina2_exe)

    if both_docking_exe_work is False:
        # Ensure permissions are unrestricted
        try:
            command = "chmod -R a+rwx {}".format(vina_exe)
            os.system(command)
            command = "chmod -R a+rwx {}".format(qvina2_exe)
            os.system(command)
        except:
            printout = "Permissions could not be adjusted on docking files."
            print(printout)
            raise Exception(printout)

        # Check Platform information
        mac_version = platform.mac_ver()[0].split(".")
        if int(mac_version[0]) < 10:
            printout = "We do not provide support for MacOS less than 10.7.\n" + \
                "Please run using docker version of AutoGrow."
            print(printout)
            raise Exception(printout)

        if int(mac_version[0]) == 10:
            if int(mac_version[1]) < 7:
                printout = "We do not support for MacOS less than 10.7.\n" + \
                    "Please run using docker version of AutoGrow."
                print(printout)
                raise Exception(printout)

            if int(mac_version[1]) > 15:
                # 10.15 is Catalina which requires notarizing docking software

                printout = "We have not tested MacOS higher than 10.15.\n" + \
                    "Please run using docker version of AutoGrow."
                print(printout)
                raise Exception(printout)

            try:
                command = "xattr -w com.apple.quarantine {}".format(vina_exe)
                os.system(command)
                command = "xattr -w com.apple.quarantine {}".format(qvina2_exe)
                os.system(command)
            except:
                printout = "Please install xattr. Can be installed using the command:"
                printout = printout  + "\n\tpip install xattr"
                print(printout)
                raise Exception(printout)

############################################
###### Variables Handlining Settings #######
############################################
def check_for_required_inputs(input_params):
    """
    Confirm all the required inputs were provided.

    Required Variables go here.

    Inputs:
    :param dict input_params: The parameters. A dictionary of {parameter name: value}.
    """
    keys_from_input = list(input_params.keys())

    list_of_required_inputs = [
        "filename_of_receptor",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "root_output_folder",
        "source_compound_file"
    ]

    missing_variables = []
    for variable in list_of_required_inputs:
        if variable in keys_from_input:
            continue
        missing_variables.append(variable)

    if len(missing_variables) != 0:
        printout = "\nRequired variables are missing from the input. A description \
            of each of these can be found by running python ./RunAutogrow -h"
        printout = printout + "\nThe following required variables are missing: "
        for variable in missing_variables:
            printout = printout + "\n\t" + variable
        print("")
        print(printout)
        print("")
        raise NotImplementedError("\n" + printout + "\n")

    # Make sure the dimmensions are in floats. If in int convert to float.
    for x in ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]:
        if type(input_params[x]) == float:
            continue
        if type(input_params[x]) == int:
            input_params[x] = float(input_params[x])
        else:
            printout = "\n{} must be a float value.\n".format(x)
            print(printout)
            raise Exception(printout)

    # Check Docking Exhaustiveness and modes...
    if "docking_exhaustiveness" in list(input_params.keys()):
        if input_params["docking_exhaustiveness"] == "None":
            input_params["docking_exhaustiveness"] = None
        if input_params["docking_exhaustiveness"] is not None:

            try:
                input_params["docking_exhaustiveness"] = int(
                    input_params["docking_exhaustiveness"]
                )
            except:
                pass
            if (
                    type(input_params["docking_exhaustiveness"]) != int
                    and type(input_params["docking_exhaustiveness"]) != float
            ):
                raise Exception(
                    "docking_exhaustiveness needs to be an interger. \
                    If you do not know what to use, leave this blank and the \
                    default for the docking software will be used."
                )
    if "docking_num_modes" in list(input_params.keys()):
        if input_params["docking_num_modes"] == "None":
            input_params["docking_num_modes"] = None
        if input_params["docking_num_modes"] is not None:
            try:
                input_params["docking_num_modes"] = int(
                    input_params["docking_num_modes"]
                )
            except:
                pass

            if (
                    type(input_params["docking_num_modes"]) != int
                    and type(input_params["docking_num_modes"]) != float
            ):
                raise Exception(
                    "docking_num_modes needs to be an interger. \
                    If you do not know what to use, leave this blank and the \
                    default for the docking software will be used."
                )

    # Check numbers which may be defined by first generation
    if "top_mols_to_seed_next_generation_first_generation" not in list(
            input_params.keys()
    ):
        if "top_mols_to_seed_next_generation" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["top_mols_to_seed_next_generation"] = 10
            input_params["top_mols_to_seed_next_generation_first_generation"] = 10
        else:
            input_params[
                "top_mols_to_seed_next_generation_first_generation"
            ] = input_params["top_mols_to_seed_next_generation"]

    if "number_of_crossovers_first_generation" not in list(input_params.keys()):
        if "number_of_crossovers" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_of_crossovers"] = 10
            input_params["number_of_crossovers_first_generation"] = 10
        else:
            input_params["number_of_crossovers_first_generation"] = input_params[
                "number_of_crossovers"
            ]

    if "number_of_mutants_first_generation" not in list(input_params.keys()):
        if "number_of_mutants" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_of_mutants"] = 10
            input_params["number_of_mutants_first_generation"] = 10
        else:
            input_params["number_of_mutants_first_generation"] = input_params[
                "number_of_mutants"
            ]

    if "number_elitism_advance_from_previous_gen_first_generation" not in list(
            input_params.keys()
    ):
        if "number_elitism_advance_from_previous_gen" not in list(input_params.keys()):
            # Use defined default of 10
            input_params["number_elitism_advance_from_previous_gen"] = 10
            input_params[
                "number_elitism_advance_from_previous_gen_first_generation"
            ] = 10
        else:
            input_params[
                "number_elitism_advance_from_previous_gen_first_generation"
            ] = input_params["number_elitism_advance_from_previous_gen"]

    #######################################
    # Check that all required files exist #
    #######################################

    # convert paths to abspath, in case necessary
    input_params["filename_of_receptor"] = os.path.abspath(
        input_params["filename_of_receptor"]
    )
    input_params["root_output_folder"] = os.path.abspath(
        input_params["root_output_folder"]
    )
    input_params["source_compound_file"] = os.path.abspath(
        input_params["source_compound_file"]
    )

    # Check filename_of_receptor exists
    if os.path.isfile(input_params["filename_of_receptor"]) is False:
        raise NotImplementedError(
            "Receptor file can not be found. File must be a .PDB file."
        )
    if ".pdb" not in input_params["filename_of_receptor"]:
        raise NotImplementedError("filename_of_receptor must be a .PDB file.")

    # Check root_output_folder exists
    if os.path.exists(input_params["root_output_folder"]) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(input_params["root_output_folder"])
        except:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

        if os.path.exists(input_params["root_output_folder"]) is False:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

    if os.path.isdir(input_params["root_output_folder"]) is False:
        raise NotImplementedError(
            "root_output_folder is not a directory. \
            Check your input parameters."
        )

    # Check source_compound_file exists
    if os.path.isfile(input_params["source_compound_file"]) is False:
        raise NotImplementedError(
            "source_compound_file can not be found. \
            File must be a tab delineated .smi file."
        )
    if ".smi" not in input_params["source_compound_file"]:
        raise NotImplementedError(
            "source_compound_file must be a \
            tab delineated .smi file."
        )



def determine_bash_timeout_vs_gtimeout():
    """
    This function tests whether we should use the BASH command "timeout" (for linux)
     or the coreutils function "gtimeout" for MacOS which can be obtained
     through homebrew

    Returns:
    :returns: str timeout_option: A string either "timeout" or "gtimeout" describing
     whether the bash terminal is able to use the bash function timeout or gtimeout
    """

    if sys.platform.lower() in ["linux", "linux2"]:
        # Should be true and default installed in all Linux machines
        return "timeout"

    command = 'timeout 1 echo " "'
    # Running the os.system command for command will return 0,1, or 32512
    # 0 means that the timeout function works (most likely this is a linux os)
    # 32512 means that the timeout function DOES NOT Work (most likely this is MacOS)

    try:  # timeout or gtimeout
        timeout_result = os.system("g" + command)
    except:
        raise Exception(
            "Something is very wrong. This OS may not be supported \
            by Autogrow or you may need to execute through Bash."
        )
    if timeout_result == 0:
        timeout_option = "gtimeout"
        return timeout_option
    print("gtimeout failed to run, we will check timeout")

    try:  # timeout or gtimeout
        timeout_result = os.system(command)
    except:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
            Autogrow or you may need to execute through Bash."
        )

    if timeout_result == 0:
        timeout_option = "timeout"
        return timeout_option

    printout = "Need to install GNU tools for Bash to work. \n"
    printout = (
        printout
        + "This is essential to use Bash Timeout function in Autogrow. \n"
    )
    printout = printout + "\t This will require 1st installing homebrew. \n"
    printout = printout + "\t\t Instructions found at: https://brew.sh/ \n"
    printout = printout + "\t Once brew is installed, please run:"
    printout = printout + " sudo brew install coreutils \n\n"
    print(printout)
    raise Exception(printout)


def check_dependencies():
    """
    This function will try to import all the installed dependencies that will be
    used in Autogrow. If it fails to import it will raise an ImportError
    """

    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option not in ["timeout", "gtimeout"]:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
        Autogrow or you may need to execute through Bash."
        )

    try:
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem.Draw import PrepareMolForDrawing
        from rdkit.Chem import rdFMCS
        from rdkit.Chem import FilterCatalog
        from rdkit.Chem.FilterCatalog import FilterCatalogParams
        import rdkit.Chem.Lipinski as Lipinski
        import rdkit.Chem.Crippen as Crippen
        import rdkit.Chem.Descriptors as Descriptors
        import rdkit.Chem.MolSurf as MolSurf

    except:
        print("You need to install rdkit and its dependencies.")
        raise ImportError("You need to install rdkit and its dependencies.")

    # molvs is prepackaged within gypsum_dl
    # try:
    #     from molvs import standardize_smiles as ssmiles
    # except:
    #     print("You need to install molvs and its dependencies.")
    #     raise ImportError("You need to install molvs and its dependencies.")

    try:
        import numpy
    except:
        print("You need to install numpy and its dependencies.")
        raise ImportError("You need to install numpy and its dependencies.")

    try:
        from scipy.cluster.vq import kmeans2
    except:
        print("You need to install scipy and its dependencies.")
        raise ImportError("You need to install scipy and its dependencies.")

    try:
        import os
        import sys
        import glob
        import subprocess
        import multiprocessing
        import time
    except:
        print(
            "Missing a Python Dependency. Could be import: os,sys,glob,\
                subprocess,multiprocess, time."
        )
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
                os,sys,glob,subprocess,multiprocess, time."
        )

    try:
        import copy
        import random
        import string
        import math
    except:
        print("Missing a Python Dependency. Could be import: copy,random, string,math")
        raise ImportError(
            "Missing a Python Dependency. Could be import: copy,random, string,math"
        )

    try:
        from collections import OrderedDict
        import webbrowser
        import argparse
        import itertools
        import unittest
    except:
        print(
            "Missing a Python Dependency. Could be import: collections,\
            webbrowser,argparse,itertools,unittest"
        )
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
            collections,webbrowser,argparse,itertools,unittest"
        )

    try:
        import textwrap
        import pickle
        import json
    except:
        print("Missing a Python Dependency. Could be import: textwrap, pickle,json")
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
            textwrap, pickle,json"
        )


def define_defaults():
    """
    Sets the command-line parameters to their default values.

    Returns:
    :returns: dict vars: a dictionary of all default variables
    """

    vars = {}

    # where we are currently (absolute filepath from route)
    # used for relative pathings
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Some variables which can be manually replaced but defaults
    # point to prepackaged locations.
    ## Neural Network executable for scoring binding
    vars["nn1_script"] = os.path.join(
        script_dir, "docking", "scoring", "nn_score_exe", "nnscore1", "NNScore.py"
    )
    # Example: vars['nn1_script'] =
    #    "/PATH/autogrow4/autogrow/docking/scoring/nn_score_exe/nnscore1/NNScore.py"

    vars["nn2_script"] = os.path.join(
        script_dir, "docking", "scoring", "nn_score_exe", "nnscore2", "NNScore2.py"
    )
    # Example: vars['nn2_script'] =
    #    "/PATH/autogrow4/autogrow/docking/scoring/nnscore2/NNScore2.py"

    #### OPTIONAL FILE-LOCATION VARIABLES ####
    # (RECOMMEND SETTING TO "" SO AUTOGROW CAN AUTOLOCATE THESE FILES)#

    # PARSER.add_argument('--conversion_choice', choices
    #    = ["MGLTools","obabel"], default="MGLTools",
    vars["conversion_choice"] = "MGLToolsConversion"
    vars["obabel_path"] = "obabel"
    vars["custom_conversion_script"] = ""
    # vars['prepare_ligand4.py'] =
    #   "/PATH/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    vars["prepare_ligand4.py"] = ""
    # vars['prepare_receptor4.py'] =
    #   "/PATH/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
    vars["prepare_receptor4.py"] = ""
    # vars['mgl_python'] = "/PATH/MGLTools-1.5.4/bin/pythonsh"
    vars["mgl_python"] = ""

    # Crossover function
    vars["start_a_new_run"] = False
    vars["max_time_mcs_prescreen"] = 1
    vars["max_time_mcs_thorough"] = 1
    vars["min_atom_match_mcs"] = 4
    vars["protanate_step"] = False

    # Mutation Settings
    vars["rxn_library"] = "click_chem_rxns"
    vars["rxn_library_file"] = ""
    vars["function_group_library"] = ""
    vars["complementary_mol_directory"] = ""

    # processors
    vars["number_of_processors"] = 1
    vars["multithread_mode"] = "multithreading"

    # Genetic Algorithm Components
    vars["selector_choice"] = "Roulette_Selector"
    vars["tourn_size"] = 0.1

    # Seeding next gen and diversity
    vars["top_mols_to_seed_next_generation_first_generation"] = 10
    vars["top_mols_to_seed_next_generation"] = 10
    vars["diversity_mols_to_seed_first_generation"] = 10
    vars["diversity_seed_depreciation_per_gen"] = 2

    # Populations settings
    vars["filter_source_compounds"] = True
    vars["use_docked_source_compounds"] = True
    vars["num_generations"] = 10
    vars["number_of_crossovers_first_generation"] = 10
    vars["number_of_mutants_first_generation"] = 10
    vars["number_of_crossovers"] = 10
    vars["number_of_mutants"] = 10
    vars["number_elitism_advance_from_previous_gen"] = 10
    vars["number_elitism_advance_from_previous_gen_first_generation"] = 10
    vars["redock_elite_from_previous_gen"] = False

    # Filters
    vars["LipinskiStrictFilter"] = False
    vars["LipinskiLenientFilter"] = False
    vars["GhoseFilter"] = False
    vars["GhoseModifiedFilter"] = False
    vars["MozziconacciFilter"] = False
    vars["VandeWaterbeemdFilter"] = False
    vars["PAINSFilter"] = False
    vars["NIHFilter"] = False
    vars["BRENKFilter"] = False
    vars["No_Filters"] = False
    vars["alternative_filter"] = None

    # docking
    vars["dock_choice"] = "QuickVina2Docking"
    vars["docking_executable"] = None
    vars["docking_exhaustiveness"] = None
    vars["docking_num_modes"] = None
    vars["docking_timeout_limit"] = 120
    vars["custom_docking_script"] = ""

    # scoring
    vars["scoring_choice"] = "VINA"
    vars["rescore_lig_efficiency"] = False
    vars["custom_scoring_script"] = ""

    # gypsum # max variance is the number of conformers made per ligand
    vars["max_variants_per_compound"] = 3
    vars["gypsum_thoroughness"] = 3
    vars["min_ph"] = 6.4
    vars["max_ph"] = 8.4
    vars["pka_precision"] = 1.0
    vars["gypsum_timeout_limit"] = 10

    # Other vars
    vars["debug_mode"] = False
    vars["reduce_files_sizes"] = False
    vars["generate_plot"] = True
    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option in  ["timeout", "gtimeout"]:
        vars["timeout_vs_gtimeout"] = timeout_option
    else:
        raise Exception("Something is very wrong. This OS may not be supported by \
                        Autogrow or you may need to execute through Bash.")

    return vars


############################################
######## Input Handlining Settings #########
############################################
def convert_json_params_from_unicode(params_unicode):
    """
    Set the parameters that will control this ConfGenerator object.

    :param dict params_unicode: The parameters. A dictionary of {parameter name:
                value}.
    Returns:
    :returns: dict params: Dictionary of User variables
    """
    # Also, rdkit doesn't play nice with unicode, so convert to ascii

    # Because Python2 & Python3 use different string objects, we separate their
    # usecases here.
    params = {}
    if sys.version_info < (3,):
        for param in params_unicode:
            val = params_unicode[param]
            if isinstance(val, unicode):
                val = str(val).encode("utf8")
            key = param.encode("utf8")
            params[key] = val
    else:
        for param in params_unicode:
            val = params_unicode[param]
            key = param
            params[key] = val
    return params


def check_value_types(vars, argv):
    """
    This checks that all the user variables loaded in use that same or comparable
    datatypes as the defaults in vars. This prevents type issues later in the
    simulation.

    Given the many uservars and the possibility for intentional differences,
    especially as the program is developed, this function tries to be
    NOT OPINIONATED, only correcting for several obvious and easy to correct issues
    of type discrepancies occur between argv[key] and vars[key]
        ie
            1) argv[key] = "true" and vars[key] = False
                this script will not change argv[key] to False... it will
                convert "true" to True
                ---> argv[key]=True
            2) argv[key] = "1.01" and vars[key] = 2.1
                this script will change argv[key] from "1.01" to float(1.01)

    Inputs:
    :param dict vars: Dictionary of program defaults, which will later be
        overwritten by argv values
    :param dict argv: Dictionary of User specified variables
    Returns:
    :returns: dict vars: Dictionary of program defaults, which will later
        be overwritten by argv values
    :returns: dict argv: Dictionary of User specified variables
    """
    for key in list(argv.keys()):
        if key not in list(vars.keys()):
            # Examples may be things like filename_of_receptor or
            # dimensions of the docking box
            #   Just skip these
            continue

        if type(argv[key]) != type(vars[key]):
            # Several variable default is None which means checks are
            # processed elsewhere...
            if vars[key] is None:
                # check argv[key] is "none" or "None"
                if type(argv[key]) == str:
                    if argv[key].lower() == "none":
                        argv[key] = None
                else:
                    continue

            # Handle number types
            elif type(vars[key]) == int or type(vars[key]) == float:
                if type(argv[key]) == int or type(argv[key]) == float:
                    # this is fine
                    continue
                elif type(argv[key]) == str:
                    try:
                        temp_item = float(argv[key])
                        if type(temp_item) == float:
                            argv[key] = temp_item
                        else:
                            printout = "This parameter is the wrong type.\n \t Check : "
                            printout = printout + "{} type={}\n".format(
                                key, type(argv[key])
                            )
                            printout = printout + "\t Should be type={}\n\t".format(
                                type(vars[key])
                            )
                            printout = (
                                printout
                                + "Please check Autogrow documentation using -h"
                            )
                            raise IOError(printout)
                    except:
                        printout = "This parameter is the wrong type. \n \t Check :"
                        printout = printout + " {} type={}\n".format(
                            key, type(argv[key])
                        )
                        printout = printout + "\t Should be type={}\n\t".format(
                            type(vars[key])
                        )
                        printout = (
                            printout + "Please check Autogrow documentation using -h"
                        )
                        raise IOError(printout)
                else:
                    printout = "This parameter is the wrong type. \n \t Check :"
                    printout = printout + " {} type={}\n".format(key, type(argv[key]))
                    printout = printout + "\t Should be type={}\n\t".format(
                        type(vars[key])
                    )
                    printout = printout + "Please check Autogrow documentation using -h"
                    raise IOError(printout)
            elif type(vars[key]) == bool:
                if argv[key] is None:
                    # Do not try to handle this. May make sense.
                    continue
                if type(argv[key]) == str:
                    if argv[key].lower() in ["true", "1"]:
                        argv[key] = True
                    elif argv[key].lower() in ["false", "0"]:
                        argv[key] = False
                    elif argv[key].lower() in ["none"]:
                        argv[key] = None
                    else:
                        printout = "This parameter is the wrong type. \n \t Check :"
                        printout = printout + " {} type={}\n".format(
                            key, type(argv[key])
                        )
                        printout = printout + "\t Should be type={}\n\t".format(
                            type(vars[key])
                        )
                        printout = (
                            printout + "Please check Autogrow documentation using -h"
                        )
                        raise IOError(printout)
                else:
                    printout = "This parameter is the wrong type. \n \t Check :"
                    printout = printout + " {} type={}\n".format(key, type(argv[key]))
                    printout = printout + "\t Should be type={}\n\t".format(
                        type(vars[key])
                    )
                    printout = printout + "Please check Autogrow documentation using -h"
                    raise IOError(printout)
    return vars, argv


def load_in_commandline_parameters(argv):
    """
    Load in the command-line parameters

    Inputs:
    :param dict argv: Dictionary of User specified variables

    Returns:
    :returns: dict vars: Dictionary of User variables
    :returns: str printout: a string to be printed to screen and saved to output file
    """

    vars = define_defaults()

    # Load the parameters from the json
    if "json" in argv:
        json_vars = json.load(open(argv["json"]))
        json_vars = convert_json_params_from_unicode(json_vars)
        check_for_required_inputs(json_vars)
        vars, json_vars = check_value_types(vars, json_vars)
        for key in list(json_vars.keys()):
            vars[key] = json_vars[key]

    else:
        check_for_required_inputs(argv)
        argv = handle_custom_inputs_if_argparsed(argv)
        vars, argv = check_value_types(vars, argv)
        for key in list(argv.keys()):
            vars[key] = argv[key]

    vars = multiprocess_handling(vars)

    printout = "(RE)STARTING AUTOGROW 4.0: " + str(datetime.datetime.now())
    printout = printout + program_info()
    printout = (
        printout + "\nUse the -h tag to get detailed help regarding program usage.\n"
    )
    print(printout)
    sys.stdout.flush()
    # Check all Dependencies are installed
    check_dependencies()

    vars = filter_choice_handling(vars)

    ###########################################
    ########## Check variables Exist ##########
    ###########################################

    # Check if Custom docking option if so there's a few things which
    # need to also be specified
    # if not lets flag the error
    if vars["dock_choice"] == "Custom":
        if vars["docking_executable"] is None:
            raise ValueError(
                "TO USE Custom DOCKING OPTION, MUST SPECIFY THE \
                PATH TO THE docking_executable AND THE DOCKING_CLASS"
            )
        if os.path.exists(vars["docking_executable"]) is False:
            raise ValueError(
                "Custom docking_executable could not be found at:\
                {}".format(
                    vars["docking_executable"]
                )
            )
        if (
                type(vars["custom_docking_script"]) != list
                or os.path.exists(vars["custom_docking_script"][1]) is not True
        ):
            raise ValueError(
                "TO USE Custom DOCKING OPTION, MUST SPECIFY THE \
                PATH TO THE Custom DOCKING SCRIPT"
            )
    if vars["dock_choice"] in ["VinaDocking", "QuickVina2Docking"]:
        if sys.platform.lower() == "darwin":
            # Some MacOS require docking software to be notarized.
            # This will require an internet signal
            run_macos_notarization(vars)

    if vars["conversion_choice"] == "Custom":
        if (
                type(vars["custom_conversion_script"]) != list
                or os.path.exists(vars["custom_conversion_script"][1]) is not True
        ):

            raise ValueError(
                "TO USE Custom conversion_choice OPTION, \
                MUST SPECIFY THE PATH TO THE custom Conversion SCRIPT"
            )

    if vars["scoring_choice"] == "Custom":
        if (
                type(vars["custom_scoring_script"]) != list
                or os.path.exists(vars["custom_scoring_script"][1]) is not True
        ):

            raise ValueError(
                "TO USE custom scoring_choice OPTION, \
                MUST SPECIFY THE PATH TO THE Custom SCORING SCRIPT"
            )

    if (
            vars["conversion_choice"] == "Custom"
            or vars["dock_choice"] == "Custom"
            or vars["scoring_choice"] == "Custom"
    ):
        vars = handle_custom_dock_and_conversion_scoring_options(vars)

    # Mutation Settings
    if vars["rxn_library"] == "Custom":
        if vars["rxn_library_file"] == "" or vars["function_group_library"] == "":
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY \
                 THE PATH TO THE REACTION LIBRARY USING INPUT PARAMETER rxn_library"
            )
        if os.path.exists(vars["rxn_library_file"]) is False:
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY \
                THE PATH TO THE REACTION LIBRARY USING INPUT PARAMETER rxn_library"
            )

        if vars["complementary_mol_directory"] == "":
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY THE PATH \
                TO THE REACTION LIBRARY USING INPUT PARAMETER function_group_library"
            )
        if os.path.isdir(vars["complementary_mol_directory"]) is False:
            raise ValueError(
                "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY THE PATH \
                TO THE REACTION LIBRARY USING INPUT PARAMETER complementary_mol_directory"
            )
    else:  # Using default settings
        if vars["rxn_library_file"] != "":
            raise ValueError(
                "You have selected a Custom rxn_library_file group \
            library but not chosen to use the Custom option for rxn_library. \
            Please use either the provided rxn_library options or chose the Custom \
            option for rxn_library"
            )
        if vars["function_group_library"] != "":
            raise ValueError(
                "You have selected a Custom function_group_library but \
            not chosen to use the Custom option for rxn_library. Please use \
            either the provided rxn_library options or chose the Custom option \
            for rxn_library"
            )
        if vars["complementary_mol_directory"] != "":
            raise ValueError(
                "You have selected a Custom complementary_mol_directory\
            but not chosen to use the Custom option for rxn_library. \
            Please use either the provided rxn_library options or chose the Custom\
            option for rxn_library"
            )

    # Check if the Operating System is Windows, if so turn off Multiprocessing.
    if os.name == "nt" or os.name == "ce":
        # so it's running under windows. multiprocessing disabled
        vars["number_of_processors"] = 1
        printout = (
            printout
            + "\nWARNING: Multiprocessing is disabled on\
            windows machines.\n"
        )

    # convert paths to abspath, in case necessary
    vars["nn1_script"] = os.path.abspath(vars["nn1_script"])
    vars["nn2_script"] = os.path.abspath(vars["nn2_script"])

    # make sure directories end in os.sep
    if vars["root_output_folder"][-1] != os.sep:
        vars["root_output_folder"] = vars["root_output_folder"] + os.sep

    # If MGLTools is being used handle its paths
    if vars["conversion_choice"] == "MGLToolsConversion":
        if "mgltools_directory" not in vars.keys():
            printout = "\nmgltools_directory was not provide but conversion_choice"
            printout = printout + " is set to MGLToolsConversion. Please "
            printout = printout + " provide the path to the mgltools_directory\n"
            print(printout)
            raise NotImplementedError(printout)

        vars["mgltools_directory"] = os.path.abspath(
            vars["mgltools_directory"]
        )
        if os.path.exists(vars["mgltools_directory"]) is False:
            raise NotImplementedError("mgltools_directory does not exist")
        if os.path.isdir(vars["mgltools_directory"]) is False:
            raise NotImplementedError(
                "mgltools_directory is not a directory. \
            Check your input parameters."
            )
        if vars["mgltools_directory"][-1] != os.sep:
            vars["mgltools_directory"] = vars["mgltools_directory"] + os.sep

        # find other mgltools-related scripts
        if vars["prepare_ligand4.py"] == "":
            vars["prepare_ligand4.py"] = (
                vars["mgltools_directory"]
                + "MGLToolsPckgs"
                + os.sep
                + "AutoDockTools"
                + os.sep
                + "Utilities24"
                + os.sep
                + "prepare_ligand4.py"
            )
        if vars["prepare_receptor4.py"] == "":
            vars["prepare_receptor4.py"] = (
                vars["mgltools_directory"]
                + "MGLToolsPckgs"
                + os.sep
                + "AutoDockTools"
                + os.sep
                + "Utilities24"
                + os.sep
                + "prepare_receptor4.py"
            )
        if vars["mgl_python"] == "":
            vars["mgl_python"] = vars["mgltools_directory"] + "bin" + os.sep + "pythonsh"

    # More Handling for Windows OS
    # convert path names with spaces if this is windows
    if os.name == "nt" or os.name == "ce":
        # so it's running under windows. multiprocessing disabled

        if " " in vars["filename_of_receptor"]:
            vars["filename_of_receptor"] = '"' + vars["filename_of_receptor"] + '"'
        if " " in vars["root_output_folder"]:
            vars["root_output_folder"] = '"' + vars["root_output_folder"] + '"'
        if " " in vars["nn1_script"]:
            vars["nn1_script"] = '"' + vars["nn1_script"] + '"'
        if " " in vars["nn2_script"]:
            vars["nn2_script"] = '"' + vars["nn2_script"] + '"'
        # If MGLTools is being used handle its paths
        if vars["conversion_choice"] == "MGLToolsConversion":
            if " " in vars["mgltools_directory"]:
                vars["mgltools_directory"] = '"' + vars["mgltools_directory"] + '"'
            if " " in vars["prepare_ligand4.py"]:
                vars["prepare_ligand4.py"] = '"' + vars["prepare_ligand4.py"] + '"'
            if " " in vars["prepare_receptor4.py"]:
                vars["prepare_receptor4.py"] = '"' + vars["prepare_receptor4.py"] + '"'
            if " " in vars["mgl_python"]:
                vars["mgl_python"] = '"' + vars["mgl_python"] + '"'

    # output the paramters used
    printout = printout + "\nPARAMETERS" + "\n"
    printout = printout + " ========== " + "\n"

    # Make sure scripts and executables exist
    # If MGLTools is being used handle its paths
    if vars["conversion_choice"] == "MGLToolsConversion":
        if not os.path.exists(vars["prepare_ligand4.py"]) and not os.path.exists(
                vars["prepare_ligand4.py"].replace('"', "")
        ):
            printout = (
                printout
                + "\nERROR: Could not find prepare_ligand4.py at "
                + vars["prepare_ligand4.py"]
                + "\n"
            )
            print(printout)
            raise NotImplementedError(printout)
        if not os.path.exists(vars["prepare_receptor4.py"]) and not os.path.exists(
                vars["prepare_receptor4.py"].replace('"', "")
        ):
            printout = (
                printout
                + "\nERROR: Could not find prepare_receptor4.py at "
                + vars["prepare_receptor4.py"]
                + "\n"
            )
            print(printout)
            raise NotImplementedError(printout)
        if not os.path.exists(vars["mgl_python"]) and not os.path.exists(
                vars["mgl_python"].replace('"', "")
        ):
            printout = (
                printout
                + "\nERROR: Could not find pythonsh at "
                + vars["mgl_python"]
                + "\n"
            )
            print(printout)
            raise NotImplementedError(printout)

    if not os.path.exists(vars["nn1_script"]) and not os.path.exists(
            vars["nn1_script"].replace('"', "")
    ):
        printout = (
            printout
            + "\nERROR: Could not find "
            + os.path.basename(vars["nn1_script"])
            + " at "
            + vars["nn1_script"]
            + "\n"
        )
        print(printout)
        raise NotImplementedError(printout)
    if not os.path.exists(vars["nn2_script"]) and not os.path.exists(
            vars["nn2_script"].replace('"', "")
    ):
        printout = (
            printout
            + "\nERROR: Could not find "
            + os.path.basename(vars["nn2_script"])
            + " at "
            + vars["nn2_script"]
            + "\n"
        )
        print(printout)
        raise NotImplementedError(printout)
    if not os.path.exists(vars["filename_of_receptor"]):
        printout = (
            printout
            + '\nERROR: There receptor file does not exist: "'
            + vars["filename_of_receptor"]
            + '".'
            + "\n"
        )
        print(printout)
        raise NotImplementedError(printout)

    # CHECK THAT NN1/NN2 are using only traditional Vina Docking
    if vars["scoring_choice"] == "NN1" or vars["scoring_choice"] == "NN2":
        if vars["dock_choice"] != "VinaDocking":
            printout = "\n\nNeural Networks 1 and 2 (NN1/NN2) are trained on data "
            printout = printout + "using PDBQT files converted by MGLTools \n"
            printout = printout + "and docked using Autodock Vina 1.1.2.\n"
            printout = (
                printout
                + "\nUsing conversion or docking software besides" +
                " these will not work. \n"
            )
            printout = (
                printout
                + "\nPlease switch dock_choice option to VinaDocking" +
                " or deselect NN1/NN2 as the scoring_choice.\n"
            )
            print(printout)
            raise Exception(printout)

        # IF ALTERNATIVE CONVERSION OF PDB2PDBQT CHECK THAT NN1/NN2 are using only MGLTOOLS
        if vars["conversion_choice"] != "MGLToolsConversion":
            printout = "\n\nNeural Networks 1 and 2 (NN1/NN2) are trained on data "
            printout = printout + "using PDBQT files converted by MGLTools \n"
            printout = printout + "and docked using Autodock Vina 1.1.2.\n"
            printout = (
                printout
                + "\nUsing conversion or docking software besides" +
                " these will not work. \n"
            )
            printout = (
                printout
                + "Please switch conversion_choice option to MGLToolsConversion" +
                " or deselect NN1/NN2 as the scoring_choice.\n"
            )
            print(printout)
            raise Exception(printout)

    # Check if the user wants to continue a run or start a new run.
    # Make new run directory if necessary. return the Run folder path
    # The run folder path will be where we place our generations and output files
    vars["output_directory"] = set_run_directory(
        vars["root_output_folder"], vars["start_a_new_run"]
    )

    # Save variables in vars dict to a .json file for later usage and reference
    # It saves the file to the output_directory + "vars.json"
    # -If AutoGrow has been run multiple times for the same directory it
    # will save the new vars file as append a number to the file name
    # starting with 2. The util scripts will only look at the original "vars.json"
    #     ie) output_directory + "vars_2.json"
    save_vars_as_json(vars)

    return vars, printout


############################################
######### File Handlining Settings #########
############################################
def find_previous_runs(folder_name_path):
    """
    This will check if there are any previous runs in the output directory.
        - If there are it will return the interger of the number label of the last Run folder path.
            - ie if there are folders Run_0, Run_1, Run_2 the function will return int(2)
        - If there are no previous Run folders it returns None.

    Inputs:
    :param str folder_name_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files

    Returns:
    :returns: int last_run_number: the int of the last run number or None if no previous runs.
    """

    path_exists = True
    i = 0
    while path_exists is True:
        folder_path = "{}{}{}".format(folder_name_path, i, os.sep)
        if os.path.exists(folder_path):
            i = i + 1
        else:
            path_exists = False

    if i == 0:
        # There are no previous runs in this directory
        last_run_number = None
        return None

    # A previous run exists. The number of the last run.
    last_run_number = i - 1
    return last_run_number


def set_run_directory(root_folder_path, start_a_new_run):
    """
    Determine and make the folder for the run directory.
        If start_a_new_run is True    Start a frest new run.
            -If no previous runs exist in the root_folder_path then make a new
                folder named root_folder_path + "Run_0"
            -If there are previous runs in the root_folder_path then make a
                new folder incremental increasing the name by 1 from the last
                run in the same output directory.
        If start_a_new_run is False    Find the last run folder and return that path
            -If no previous runs exist in the root_folder_path then make a new
            folder named root_folder_path + "Run_0"

    Inputs:
    :param str root_folder_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files
    :param bol start_a_new_run: True or False to determine if we continue from
        the last run or start a new run
        - This is set as a vars["start_a_new_run"]
        - The default is vars["start_a_new_run"] = True
    Returns:
    :returns: str folder_path: the string of the newly created directory for
        puting output folders
    """

    folder_name_path = root_folder_path + "Run_"
    print(folder_name_path)

    last_run_number = find_previous_runs(folder_name_path)

    if last_run_number is None:
        # There are no previous simulation runs in this directory
        print("There are no previous runs in this directory.")
        print("Starting a new run named Run_0.")

        # make a folder for the new generation
        run_number = 0
        folder_path = "{}{}{}".format(folder_name_path, run_number, os.sep)
        os.makedirs(folder_path)

    else:
        if start_a_new_run is False:
            # Continue from the last simulation run
            run_number = last_run_number
            folder_path = "{}{}{}".format(folder_name_path, last_run_number, os.sep)
        else:  # start_a_new_run is True
            # Start a new fresh simulation
            # Make a directory for the new run by increasing run number by +1
            # from last_run_number
            run_number = last_run_number + 1
            folder_path = "{}{}{}".format(folder_name_path, run_number, os.sep)
            os.makedirs(folder_path)

    print("The Run number is: ", run_number)
    print("The Run folder path is: ", folder_path)
    print("")
    return folder_path


############################################
########   Custom Option Settings   ########
############################################
def handle_custom_inputs_if_argparsed(input_params):
    """
    There are several Custom options such as filters, docking software
    which take a list of information. Because Filters can use multiple options
    at once it takes a list of list information.
    This function is used to properly import and parse those user variables if
     using the commandline argparse

    This function will handle those if there are used and return
    the modified input_params dict

    Inputs:
    :param dict input_params: The parameters. A dictionary of
        {parameter name: value}.
    Returns:
    :returns: dict input_params: The parameters. A dictionary of
        {parameter name: value}.
    """

    # Custom Filters
    if "alternative_filter" not in input_params.keys():
        input_params["alternative_filter"] = None
    if (
            input_params["alternative_filter"] is not None
            and input_params["alternative_filter"] != []
    ):
        orginal = input_params["alternative_filter"][0]
        orginal = orginal.replace("[[", "[").replace("]]", "]")
        new_alternative_filter = []
        for custom_filter in orginal.split("]"):
            custom_filter = custom_filter.replace("[", "").replace("]", "")
            custom_filter = [x for x in custom_filter.split(",") if x != ""]
            if len(custom_filter) == 2:
                new_alternative_filter.append(custom_filter)
        input_params["alternative_filter"] = new_alternative_filter

    # custom_conversion_script
    if "custom_conversion_script" not in input_params.keys():
        input_params["custom_conversion_script"] = None
    if input_params["custom_conversion_script"] not in [None, [], "", "[]"]:
        orginal = input_params["custom_conversion_script"][0]
        orginal = orginal.replace("[[", "[").replace("]]", "]")
        new_alternative_conversion = []
        for custom_converter in orginal.split("]"):
            custom_converter = custom_converter.replace("[", "").replace("]", "")
            custom_converter = [x for x in custom_converter.split(",") if x != ""]
            if len(custom_converter) == 2:
                new_alternative_conversion.append(custom_converter)
        input_params["custom_conversion_script"] = new_alternative_conversion

    # custom_docking_script
    if "custom_docking_script" not in input_params.keys():
        input_params["custom_docking_script"] = None
    if input_params["custom_docking_script"] not in [None, [], "", "[]"]:
        orginal = input_params["custom_docking_script"][0]
        orginal = orginal.replace("[[", "[").replace("]]", "]")
        new_alternative_docking = []
        for custom_docking in orginal.split("]"):
            custom_docking = custom_docking.replace("[", "").replace("]", "")
            custom_docking = [x for x in custom_docking.split(",") if x != ""]
            if len(custom_docking) == 2:
                new_alternative_docking.append(custom_docking)
        input_params["custom_docking_script"] = new_alternative_docking

    # Custom_Scoring script
    if "custom_scoring_script" not in input_params.keys():
        input_params["custom_scoring_script"] = None
    if input_params["custom_scoring_script"] not in [None, [], "", "[]"]:
        orginal = input_params["custom_scoring_script"][0]
        orginal = orginal.replace("[[", "[").replace("]]", "]")
        new_alternative_scoring = []
        for custom_scoring in orginal.split("]"):
            custom_scoring = custom_scoring.replace("[", "").replace("]", "")
            custom_scoring = [x for x in custom_scoring.split(",") if x != ""]
            if len(custom_scoring) == 2:
                new_alternative_scoring.append(custom_scoring)
        input_params["custom_scoring_script"] = new_alternative_scoring

    return input_params


#
def handle_alternative_filters(vars, filter_list):
    """
    This will handle Custom Filters

    Inputs:
    :param dict vars: Dictionary of User variables
    :param list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow the same
        nomenclature and enter that filter in the user vars["alternative_filters"]
        as the name of that class and place that file in the same folder as the
        other filter classes.

    Returns:
    :returns: list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow the same
        nomenclature and enter that filter in the user vars["alternative_filters"]
        as the name of that class and place that file in the same folder as the
        other filter classes.
    """
    if vars["alternative_filter"] is not None:
        if type(vars["alternative_filter"]) != list:
            raise Exception(
                "If you want to add Custom filters to the filter \
                child classes Must be a list of lists \
                [[name_filter1, Path/to/name_filter1.py],[name_filter2, Path/to/name_filter2.py]]"
            )
        if type(vars["alternative_filter"][0]) != list:
            print(vars["alternative_filter"])
            raise Exception(
                "If you want to add Custom filters to the filter \
                child classes Must be a list of lists \
                [[name_filter1, Path/to/name_filter1.py],[name_filter2, Path/to/name_filter2.py]]"
            )

        full_children_dict = make_complete_children_dict("filter")
        scripts_to_copy = []
        for custom_class in vars["alternative_filter"]:
            if custom_class[0] not in full_children_dict.keys():
                if os.path.exists(custom_class[1]) is False:
                    # Check that the path to the original script exists.
                    raise Exception(
                        "File can not be found for alternative_filter \
                        {}\n If you want to add Custom filters to the filter child \
                        classes Must be a list of lists \
                        [[name_filter1, Path/to/name_filter1.py],\
                        [name_filter2, Path/to/name_filter2.py]]".format(custom_class[1])
                    )

                new_file = os.sep.join(
                    [
                        os.path.abspath(os.path.dirname(__file__)),
                        "operators",
                        "filter",
                        "filter_classes",
                        "filter_children_classes",
                        os.path.basename(custom_class[0]) + ".py",
                    ]
                )

                if os.path.exists(new_file) is True:
                    # File has been copied to proper dir but is not being found by the code
                    printout = "A copy of the custom script {} has been moved \
                        to {}\n".format(custom_class[1], new_file)
                    printout = (
                        printout
                        + "Unfortunately this could not be  \
                        imported by the filter module."
                    )
                    printout = (
                        printout
                        + "Please check the file naming \
                        corresponding to: {}\n\n".format(
                            custom_class
                        )
                    )
                    print(printout)
                    raise Exception(printout)

                # Add to list of scripts to copy into the filter folder
                scripts_to_copy.append([custom_class[1], new_file])

            filter_list.append(custom_class[0])
        if len(scripts_to_copy) != 0:
            for filter_info in scripts_to_copy:
                print("copying Custom class file into the FilterClasses folder:")
                print(
                    "\t Copying : {}\n\t New file: {}\n".format(
                        custom_class[1], new_file
                    )
                )
                print(
                    "AutoGrow will need to be restarted once all custom scripts \
                    have been copied to their required location."
                )
                print(
                    "This is done once so if the script needs to be changed \
                    please either remove or replace the script within the \
                    FilterClasses folder."
                )
                print(
                    "Please ensure you unit test this code properly before \
                    incorporating.\n"
                )
                copyfile(filter_info[0], filter_info[1])

            print(
                "\n########################################"
                + "#####################################"
            )
            print("AutoGrow has incorporated the custom files into"
                  + " the filter Module.")
            print(
                " AutoGrow needs to be restarted and should now "
                + "be able to run custom scripts."
            )
            print("Please ensure you unit test this code properly before incorporating.")
            print(
                "#####################################"
                + "########################################\n"
            )
            # Technically Exit intentionally but maybe should be a raise Exception
            sys.exit(0)
    return filter_list

#
def make_complete_children_dict(purpose_of_object):
    """
    This will retrieve all the names of every child class of the parent class
    This can be either filter, parent_pdbqt_converter, ParentDocking,
    or ParentScoring

    Inputs:
    :param str purpose_of_object: either filter, parent_pdbqt_converter,
        ParentDocking, or ParentScoring
    Returns:
    :returns: dict child_dict: Dictionary of all the class objects for either
        Filtering, docking, Dockingfile conversion or scoring
    """
    if purpose_of_object == "filter":
        import autogrow.operators.filter.filter_classes.filter_children_classes
        from autogrow.operators.filter.filter_classes.parent_filter_class import ParentFilter as parent_object
        from autogrow.operators.filter.filter_classes.get_child_filter_class import get_all_subclasses

    elif purpose_of_object == "parent_pdbqt_converter":
        import autogrow.docking.docking_class.docking_file_conversion
        from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter as parent_object
        from autogrow.docking.docking_class.get_child_class import get_all_subclasses
    elif purpose_of_object == "ParentDocking":
        import autogrow.docking.docking_class.docking_class_children
        from autogrow.docking.docking_class.parent_dock_class import ParentDocking as parent_object
        from autogrow.docking.docking_class.get_child_class import get_all_subclasses
    elif purpose_of_object == "ParentScoring":
        import autogrow.docking.scoring.scoring_classes.scoring_functions
        from autogrow.docking.scoring.scoring_classes.parent_scoring_class import ParentScoring as parent_object
        from autogrow.docking.docking_class.get_child_class import get_all_subclasses

    children = get_all_subclasses(parent_object)
    child_dict = {}
    for child in children:
        child_object = child()
        child_name = child_object.get_name()
        child_dict[child_name] = child_object

    return child_dict


#
def handle_custom_conversion_script(vars):
    """
    This will handle Custom Conversion_scripts

    Inputs:
    :param dict vars: Dictionary of User variables
    Returns:
    :returns: dict vars: Dictionary of User variables modified with
        the vars["conversion_choice"] set to the new custom conversion_choice
    :returns: bool need_restart: If True AutoGrow will need to be restarted
        after all other files are incorporated
    :returns: str printout: "" or a message to be print prior to being
        restarted if needed
    """
    need_restart = False
    printout = ""
    if vars["custom_conversion_script"] is not None:
        if type(vars["custom_conversion_script"]) != list:
            print(vars["custom_conversion_script"])
            raise Exception(
                "If you want to add Custom Conversion_script \
                to the Conversion_script child classes Must be a list of \
                [name_Conversion_script1, Path/to/name_Conversion_script1.py]"
            )
        if type(vars["custom_conversion_script"][0]) != str:
            print("")
            print(vars["custom_conversion_script"])
            print("")
            raise Exception(
                "If you want to add Custom Conversion_script \
                to the Conversion_script child classes Must be a list of \
                [name_Conversion_script1, Path/to/name_Conversion_script1.py]"
            )

        full_children_dict = make_complete_children_dict("parent_pdbqt_converter")
        custom_class = vars["custom_conversion_script"]
        if custom_class[0] not in full_children_dict.keys():
            if os.path.exists(custom_class[1]) is False:
                print(custom_class)
                raise Exception(
                    "File can not be found for custom_conversion_script \
                    {}\n If you want to add Custom Conversion_scripts to the \
                    Conversion_script child classes Must be a list of \
                    [name_Conversion_script1, Path/to/name_Conversion_script1.py]".format(
                        custom_class[1]
                    )
                )

            new_file = os.sep.join(
                [
                    os.path.abspath(os.path.dirname(__file__)),
                    "docking",
                    "docking_class",
                    "docking_file_conversion",
                    os.path.basename(custom_class[0]) + ".py",
                ]
            )

            if os.path.exists(new_file) is True:
                # File has been copied to proper dir but is not being found by the code
                printout = "A copy of the custom script {} has been moved \
                    to {}\n".format(custom_class[1], new_file)
                printout = (
                    printout
                    + "Unfortunately this could not be \
                    imported by the Conversion_script module."
                )
                printout = (
                    printout
                    + "Please check the file naming corresponding \
                    to: {}\n\n".format(
                        custom_class
                    )
                )
                print(printout)
                raise Exception(printout)

            # Add copy the script to the docking_file_conversion folder
            print("copying Custom class file into the Conversion_script folder:")
            print(
                "\t Copying : {}\n\t New file: {}\n".format(
                    custom_class[1], new_file
                )
            )
            print(
                "AutoGrow will need to be restarted once the custom script \
                has been copied to their required location."
            )
            print(
                "This is done once so if the script needs to be changed \
                please either remove or replace the script within the \
                docking_file_conversion folder."
            )
            print(
                "Please ensure you unit test this code properly before \
                incorporating."
            )
            copyfile(custom_class[1], new_file)

            printout = (
                printout
                + "\n#########################################"
                + "####################################"
            )
            printout = (
                printout
                + "AutoGrow has incorporated the custom files into "
                + "the docking_file_conversion Module."
            )
            printout = (
                printout
                + "AutoGrow needs to be restarted and should now be "
                + "able to run custom scripts."
            )
            printout = (
                printout
                + "Please ensure you unit test this code properly "
                + "before incorporating."
            )
            printout = (
                printout
                + "#########################################"
                + "####################################\n"
            )

            need_restart = True

        vars["conversion_choice"] = custom_class[0]
    return vars, need_restart, printout


#
def handle_custom_docking_script(vars):
    """
    This will handle Custom Docking_scripts

    Inputs:
    :param dict vars: Dictionary of User variables
    Returns:
    :returns: dict vars: Dictionary of User variables modified with the
        vars["dock_choice"] set to the new custom dock_choice
    :returns: bool need_restart: If True AutoGrow will need to be estarted
         after all other files are incorporated
    :returns: str printout: "" or a message to be print prior to being
        restarted if needed
    """
    need_restart = False
    printout = ""
    if vars["custom_docking_script"] is not None:
        if type(vars["custom_docking_script"]) != list:
            print(vars["custom_docking_script"])
            raise Exception(
                "If you want to add Custom Docking_script to the \
                Docking_script child classes Must be a list of \
                [name_Docking_script1, Path/to/name_Docking_script1.py]"
            )
        if type(vars["custom_docking_script"][0]) != str:
            print("")
            print(vars["custom_docking_script"])
            print("")
            raise Exception(
                "If you want to add Custom Docking_script to the \
                Docking_script child classes Must be a list of \
                [name_Docking_script1, Path/to/name_Docking_script1.py]"
            )

        full_children_dict = make_complete_children_dict("ParentDocking")
        custom_class = vars["custom_docking_script"]
        if custom_class[0] not in full_children_dict.keys():
            if os.path.exists(custom_class[1]) is False:
                print(custom_class)
                raise Exception(
                    "File can not be found for custom_docking_script \
                    {}\n If you want to add Custom Docking_scripts to the \
                    Docking_script child classes Must be a list of \
                    [name_Docking_script1, Path/to/name_Docking_script1.py]".format(
                        custom_class[1]
                    )
                )

            new_file = os.sep.join(
                [
                    os.path.abspath(os.path.dirname(__file__)),
                    "docking",
                    "docking_class",
                    "docking_class_children",
                    os.path.basename(custom_class[0]) + ".py",
                ]
            )

            if os.path.exists(new_file) is True:
                # File has been copied to proper dir but is not being found by the code
                printout = "A copy of the custom script {} has been moved \
                    to {}\n".format(
                        custom_class[1], new_file
                    )
                printout = (
                    printout
                    + "Unfortunately this could not be imported \
                    by the docking module."
                )
                printout = (
                    printout
                    + "Please check the file naming corresponding \
                    to: {}\n\n".format(
                        custom_class
                    )
                )
                print(printout)
                raise Exception(printout)

            # Add copy the script to the children folder
            print("copying Custom class file into the children folder:")
            print(
                "\t Copying : {}\n\t New file: {}\n".format(
                    custom_class[1], new_file
                )
            )
            print(
                "AutoGrow will need to be restarted once the custom \
                script has been copied to their required location."
            )
            print(
                "This is done once so if the script needs to be changed \
                please either remove or replace the script within the \
                children folder."
            )
            print(
                "Please ensure you unit test this code properly before incorporating."
            )
            copyfile(custom_class[1], new_file)

            printout = (
                printout
                + "\n############################################"
                + "#################################"
            )
            printout = (
                printout
                + "AutoGrow has incorporated the custom files into the children Module."
            )
            printout = (
                printout
                + "AutoGrow needs to be restarted and should now be able to run custom scripts."
            )
            printout = (
                printout
                + "Please ensure you unit test this code properly before incorporating."
            )
            printout = (
                printout
                + "##############################################"
                + "###############################\n"
            )

            need_restart = True

        vars["dock_choice"] = custom_class[0]
    return vars, need_restart, printout


#
def handle_custom_scoring_script(vars):
    """
    This will handle Custom scoring_scripts

    Inputs:
    :param dict vars: Dictionary of User variables
    Returns:
    :returns: dict vars: Dictionary of User variables modified with the
        vars["dock_choice"] set to the new custom dock_choice
    :returns: bool need_restart: If True AutoGrow will need to be restarted
        after all other files are incorporated
    :returns: str printout: "" or a message to be print prior to
        being restarted if needed
    """
    need_restart = False
    printout = ""
    if vars["custom_scoring_script"] is not None:
        if type(vars["custom_scoring_script"]) != list:
            print(vars["custom_scoring_script"])
            raise Exception(
                "If you want to add Custom scoring_script \
                to the scoring_script child classes Must be a list of \
                [name_scoring_script1, Path/to/name_scoring_script1.py]"
            )
        if type(vars["custom_scoring_script"][0]) != str:
            print("")
            print(vars["custom_scoring_script"])
            print("")
            raise Exception(
                "If you want to add Custom scoring_script \
                to the scoring_script child classes Must be a list of \
                [name_scoring_script1, Path/to/name_scoring_script1.py]"
            )

        full_children_dict = make_complete_children_dict("ParentScoring")
        custom_class = vars["custom_scoring_script"]
        if custom_class[0] not in full_children_dict.keys():
            if os.path.exists(custom_class[1]) is False:
                print(custom_class)
                raise Exception(
                    "File can not be found for custom_scoring_script \
                    {}\n If you want to add Custom scoring_scripts to the \
                    scoring_script child classes Must be a list of \
                    [name_scoring_script1, Path/to/name_scoring_script1.py]".format(
                        custom_class[1]
                    )
                )

            new_file = os.sep.join(
                [
                    os.path.abspath(os.path.dirname(__file__)),
                    "docking",
                    "scoring",
                    "scoring_classes",
                    "scoring_functions",
                    os.path.basename(custom_class[0]) + ".py",
                ]
            )

            if os.path.exists(new_file) is True:
                # File has been copied to proper dir but is not being found by the code
                printout = "A copy of the custom script {} has been moved to {}\n".format(
                    custom_class[1], new_file
                )
                printout = (
                    printout
                    + "Unfortunately this could not be imported by the scoring module."
                )
                printout = (
                    printout
                    + "Please check the file naming corresponding to: {}\n\n".format(
                        custom_class
                    )
                )
                print(printout)
                raise Exception(printout)

            # Add copy the script to the scoring_choices folder
            print("copying Custom class file into the scoring_choices folder:")
            print(
                "\t Copying : {}\n\t New file: {}\n".format(
                    custom_class[1], new_file
                )
            )
            print(
                "AutoGrow will need to be restarted once the custom script \
                has been copied to their required location."
            )
            print(
                "This is done once so if the script needs to be changed \
                please either remove or replace the script within \
                the scoring_choices folder."
            )
            print(
                "Please ensure you unit test this code properly before incorporating."
            )
            copyfile(custom_class[1], new_file)

            printout = "\n#######################################"
            printout = printout + "######################################"
            printout = (
                printout
                + "AutoGrow has incorporated the custom files into the scoring Module."
            )
            printout = (
                printout
                + "AutoGrow needs to be restarted and should now be able to run custom scripts."
            )
            printout = (
                printout
                + "Please ensure you unit test this code properly before incorporating."
            )
            printout = (
                printout
                + "##############################################"
                + "###############################\n"
            )

            need_restart = True

        vars["scoring_choice"] = custom_class[0]
    return vars, need_restart, printout


#
def handle_custom_dock_and_conversion_scoring_options(vars):
    """
    This function handles selecting the user defined Custom options
    for Custom docking Conversion, and scoring scripts.

    Inputs:
    :param dict vars: Dictionary of User variables
    Returns:
    :returns: dict vars: Dictionary of User variables with the added options
    """
    master_need_restart = False
    master_printout = ""
    if vars["conversion_choice"] == "Custom":
        vars, need_restart, printout = handle_custom_conversion_script(vars)
        if need_restart is True:
            master_need_restart = True
            master_printout = master_printout + printout
    if vars["dock_choice"] == "Custom":
        vars, need_restart, printout = handle_custom_docking_script(vars)
        if need_restart is True:
            master_need_restart = True
            master_printout = master_printout + printout
    if vars["scoring_choice"] == "Custom":
        vars, need_restart, printout = handle_custom_scoring_script(vars)
        if need_restart is True:
            master_need_restart = True
            master_printout = master_printout + printout

    if master_need_restart is True:
        print(master_printout)
        sys.exit(
            0
        )  # Technically Exit intentionally but maybe should be a raise Exception

    return vars


############################################
######## Filter Handlining Settings ########
############################################
def filter_choice_handling(vars):
    """
    This function handles selecting the user defined Ligand filters.

    Inputs:
    :param dict vars: Dictionary of User variables
    Returns:
    :returns: dict vars: Dictionary of User variables with the
        chosen_ligand_filters added
    """
    if "No_Filters" in list(vars.keys()):
        if vars["No_Filters"] is True:
            chosen_ligand_filters = None
        else:
            chosen_ligand_filters, vars = picked_filters(vars)
    else:
        chosen_ligand_filters, vars = picked_filters(vars)
    vars["chosen_ligand_filters"] = chosen_ligand_filters

    import autogrow.operators.filter.execute_filters as Filter


    # get child filter class object function dictionary
    vars["filter_object_dict"] = Filter.make_run_class_dict(chosen_ligand_filters)

    return vars


#
def picked_filters(vars):
    """
    This will take the user vars and return a list of the filters
    which a molecule must pass to move into the next generation.

    Inputs:
    :param dict vars: Dictionary of User variables
    Returns:
    :returns: list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow
        the same nomenclature and enter that filter in the user
        vars["alternative_filters"] as the name of that class and place
        that file in the same folder as the other filter classes.
    """
    filter_list = []
    vars_keys = list(vars.keys())

    if "LipinskiStrictFilter" in vars_keys:
        if vars["LipinskiStrictFilter"] is True:
            filter_list.append("LipinskiStrictFilter")
    else:
        vars["LipinskiStrictFilter"] = False

    if "LipinskiLenientFilter" in vars_keys:
        if vars["LipinskiLenientFilter"] is True:
            filter_list.append("LipinskiLenientFilter")
    else:
        vars["LipinskiLenientFilter"] = False

    if "GhoseFilter" in vars_keys:
        if vars["GhoseFilter"] is True:
            filter_list.append("GhoseFilter")
    else:
        vars["GhoseFilter"] = False

    if "GhoseModifiedFilter" in vars_keys:
        if vars["GhoseModifiedFilter"] is True:
            filter_list.append("GhoseModifiedFilter")
    else:
        vars["GhoseModifiedFilter"] = False

    if "MozziconacciFilter" in vars_keys:
        if vars["MozziconacciFilter"] is True:
            filter_list.append("MozziconacciFilter")
    else:
        vars["MozziconacciFilter"] = False

    if "VandeWaterbeemdFilter" in vars_keys:
        if vars["VandeWaterbeemdFilter"] is True:
            filter_list.append("VandeWaterbeemdFilter")
    else:
        vars["VandeWaterbeemdFilter"] = False

    if "PAINSFilter" in vars_keys:
        if vars["PAINSFilter"] is True:
            filter_list.append("PAINSFilter")
    else:
        vars["PAINSFilter"] = False

    if "NIHFilter" in vars_keys:
        if vars["NIHFilter"] is True:
            filter_list.append("NIHFilter")
    else:
        vars["NIHFilter"] = False

    if "BRENKFilter" in vars_keys:
        if vars["BRENKFilter"] is True:
            filter_list.append("BRENKFilter")
    else:
        vars["BRENKFilter"] = False

    if "alternative_filter" in vars_keys:
        filter_list = handle_alternative_filters(vars, filter_list)
    else:
        vars["alternative_filter"] = None

    # if there is no user specified ligand filters but they haven't set
    # filters to None ---> set filter to default of LipinskiLenientFilter.
    if len(filter_list) == 0:
        vars["LipinskiLenientFilter"] = True
        filter_list.append("LipinskiLenientFilter")

    return filter_list, vars
