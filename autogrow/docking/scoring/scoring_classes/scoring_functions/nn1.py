"""
This script contains the class NN1 that rescores Vina type docking
using the program NNScore1.
"""
import __future__

import glob
import os
import sys


from autogrow.docking.scoring.scoring_classes.parent_scoring_class import ParentScoring
from autogrow.docking.scoring.scoring_classes.scoring_functions.vina import VINA


class NN1(VINA):
    """
    This will Score a given ligand for its binding affinity based on VINA or
    QuickVina02 type docking.

    Inputs:
    :param class ParentFilter: a parent class to initialize off of.
    """

    def __init__(self, vars=None, smiles_dict=None, test_boot=True):
        """
        This will take vars and a list of smiles.

        Inputs:
        :param dict vars: Dictionary of User variables
        :param dict smiles_dict: a dict of ligand info of SMILES, IDS, and
            short ID
        :param bool test_boot: used to initialize class without objects for
            testing purpose
        """

        if test_boot is False:
            self.vars = vars

            self.smiles_dict = smiles_dict
            print("")
            print("######################")
            print("Running NN1 rescoring on vina files")


    #######################
    # Executed by the Execute_Scoring.py script
    #######################
    def find_files_to_score(self, file_path):
        """
        Find all files of the appropriate file format within the dir. For this
        class its .pdbqt.vina files.

        ALL SCORING FUNCTIONS MUST HAVE THIS FUNCTION.

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list list_of_NN1_files: list of all files to be scored
            within the dir
        """

        self.file_path = file_path

        list_of_docked_files = []

        list_of_docked_files = glob.glob(file_path + "*.pdbqt.vina")

        return list_of_docked_files

    def run_rescoring(self, vina_output_file):
        """
        Run the NN1 scoring on all of these files. Return a list of a rescored
        file with NN1 ie *.pdbqt.vina.nn1

        Inputs:
        :param str vina_output_file: Path to a vina output file to be rescored

        Returns:
        :returns: list results of the rescoring function: [file_path,
            it_rescored]. [PATH, True] means it passed. [PATH, False] means it
            failed a results of all NN1 files.
        """

        result_of_rescore = run_nn_rescoring(self.vars, vina_output_file)

        return result_of_rescore

    def run_scoring(self, file_path):
        """
        Get all relevant scoring info and return as a list

        This is required for all Scoring Functions. Additional manipulations
        may go here but there are none for this script..

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list list_of_lig_data: information about the scored ligand.
            Score is last index (ie. [lig_id_shortname, any_details,
            fitness_score_to_use] )
        """

        if os.path.exists(file_path):
            lig_info = self.get_score_from_a_file(file_path)
            return lig_info

        # file_path does not exist
        return None

    def get_score_from_a_file(self, file_path):
        """
        Make a list of a ligands information including its docking score.

        Because a higher score is better for both NNScore functions, but
        AutoGrow4 selects based on most negative score, we multiple each NNScore
        value by -1.0. This ensures that the best score is the most negative
        score.

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list lig_info: a list of the ligands short_id_name and the
            docking score from the best pose.
        """

        if ".nn1" not in file_path:
            if ".vina" in file_path:
                file_path = file_path + ".nn1"
                if os.path.exists(file_path) is False:
                    return None
            else:
                return None
        if os.path.exists(file_path) is False:
            return None
        # grab the index of the ligand for the score
        basefile = os.path.basename(file_path)
        ligand_pose = basefile.split(".pdbqt.vina.nn1")[0]
        basefile_split = basefile.split("__")
        ligand_short_name = basefile_split[0]

        score = None

        with open(file_path, "r") as f:
            for line in f.readlines():

                if "Best score: " in line:
                    try:
                        tmp = line.split(" ")

                        # because for both NNScore functions, a higher score is better
                        # multiply score to be negative
                        temp_score = tmp[2].replace(",", "")
                        temp_score = float(temp_score) * float(-1.0)
                    except:
                        continue
                    if score is None:
                        score = temp_score
                    else:
                        if score > temp_score:
                            score = temp_score
                        else:
                            continue

                    # Pose info not used here but could be useful in some
                    # context
                    # tmp = line.replace(")",'')
                    # tmp = tmp.split(', MODEL ')
                    # vina_pose_to_use = int(tmp[1])

        if score is None:
            # This file lacks a pose to use
            return None

        lig_info = [ligand_short_name, ligand_pose, score]

        # Obtain additional file info
        lig_info = self.merge_smile_info_w_affinity_info(lig_info)

        if lig_info is None:
            return None
        lig_info = [str(x) for x in lig_info]

        return lig_info


###Outside class for multithreading
# Run NN1 rescoring
def run_nn_rescoring(vars, vina_output_file):
    """
    This will run NN1 on all of the vina files in the list. This is outside
    the class for multifunction purposes

    Returns A list containing the file name as item 1 and whether it passed as
    item 2. [PATH, True] means it passed. [PATH, False] means it failed a
    results of all NN1 files.

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param str vina_output_file: Path to a vina output file to be rescored

    Returns:
    :returns: list results of the rescoring function: [file_path,
        it_rescored]. [PATH, True] means it passed. [PATH, False] means it failed
        a results of all NN1 files.
    """

    if vina_output_file is None:
        return None
    # Unpackage vars
    receptor = vars["filename_of_receptor"] + "qt"
    nn1_executable = vars["nn1_script"]
    networks_dir = (
        os.path.dirname(vars["nn1_script"])
        + os.sep
        + "networks"
        + os.sep
        + "top_3_networks"
        + os.sep
    )

    nn1_output = vina_output_file + ".nn1"
    # sys.executable is the path to python executable
    torun = (
        sys.executable + " "
        + nn1_executable
        + " -receptor "
        + receptor
        + " -vina_output "
        + vina_output_file
        + " -networks_dir "
        + networks_dir
        + " > "
        + nn1_output
    )

    # A list containing the file name as item 1 and whether it passed as item
    # 2
    results = execute_nn_scoring(torun, nn1_output)

    # Will be None if it passed. A list containing the file name as item 1 and
    # whether it passed as item 2. [PATH, True] means it passed. [PATH, False]
    # means it failed.
    return results


def execute_nn_scoring(command, file_path):
    """
    Run an individual NN scoring function.

    returns A list containing the file name as item 1 and whether it passed as
    item 2. [PATH, True] means it passed. [PATH, False] means it failed.

    Inputs:
    :param str command: the rescoring bash style command to execute
    :param str file_path: Path to a vina output file to be rescored

    Returns:
    :returns: list results of the rescoring function: [file_path,
        it_rescored]. [PATH, True] means it passed. [PATH, False] means it failed
        a results of all NN1 files.
    """

    try:
        os.system(command)
        it_rescored = confirm_file_has_scoring(file_path)
    except:
        return [file_path, False]
    return [file_path, it_rescored]


def confirm_file_has_scoring(file_path):
    """
    Check the file has a rescore value

    Inputs:
    :param str file_path: Path to a vina output file to be rescored
    Returns:
    :returns: bol has_scoring: True if has score;
        False if no score found
    """

    if os.path.exists(file_path) is False:
        return False
    with open(file_path, "r") as f:
        has_scoring = False
        for line in f.readlines():

            if "Best Score:" in line or "Best score:":
                has_scoring = True
                return has_scoring
    return has_scoring
