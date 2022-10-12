"""
This script will convert a docked .pdbqt.vina file into separate .pdb file.

This is done by splitting up a single .pdbqt.vina into separate .pdbqt
files for each docked pose.
Then it removes a column of the .pdbqt and saves as a .pdb file.

If variable --max_num_of_poses is not set it will convert all poses.
    If --max_num_of_poses == 1 it will only convert the top docked pose to .pdb
    If --max_num_of_poses == 2 it will only convert the top 2 docked poses to .pdb
    If --max_num_of_poses == 10 but there only 8 poses it will convert the 8 poses and stop


If --max_docking_score is not set it will convert all poses to .pdb;
    If --max_docking_score == -10.0 it will only convert poses with docking scores less than
        or equal to -10.0 (Remember docking scores are better when more negative)

--max_docking_score and --max_num_of_poses work as AND type operators.
    Example:
        --max_docking_score == -11.4 and --max_num_of_poses==5
        It will take the top 5 poses as long as they also have docking scores <=-11.4

Example submit:
python PATH/autogrow4/accessory_scripts/convert_vina_docked_pdbqt_to_pdbs.py \
    --vina_docked_pdbqt_file \
        PATH/Run_1/Run_0/generation_30/PDBs/Gen_30_Cross_313228__1.pdbqt.vina \
    --output_folder PATH/outfolder/ \
    --max_num_of_poses 1 --number_of_processors -1
"""

import os
import copy
import argparse
import glob

import support_scripts.Multiprocess as mp

def run_conversion_for_a_vina_file(vina_file, output_folder, max_num_of_poses,
                                   max_docking_score, min_docking_score):
    """
    This script will convert .pdbqt.vina file to multiple .pdb files.

    If criteria such as  max_num_of_poses, max_docking_score,
    min_docking_score a pose must meet all criteria to be converted.

    Inputs:
    :param str vina_file: Path to vina file to convert
    :param str output_folder: Path to vina file to output folder
    :param int max_num_of_poses: Max number of poses to convert to pdb
    :param float max_docking_score: Most positive docking score to be converted; float or None
    :param float min_docking_score: Most negative docking score to be converted; float or None
    """

    if os.path.exists(vina_file) is False:
        raise Exception("CANT FIND FILE:", vina_file)

    if os.path.exists(output_folder) is False:
        raise Exception("CANT FIND outfolder:", output_folder)

    short_name = os.path.basename(vina_file).replace(".pdbqt.vina", "")

    with open(vina_file, "r") as f:
        pose_number = 1
        printout_list = []
        score = 0.0
        terminate_run = False
        write_pose = True
        for line in f.readlines():
            if pose_number > max_num_of_poses and max_num_of_poses != -1:
                # break if hit max number of poses
                break
            if terminate_run is True:
                break

            if "REMARK VINA RESULT" in line:
                write_pose = True
                if max_docking_score is None and min_docking_score is None:
                    printout_list.append(line)
                else:

                    temp_line = copy.deepcopy(line)
                    for i in range(10):
                        temp_line = temp_line.replace("  ", " ")
                    temp_line = temp_line.split("RESULT:")[1]
                    temp_line = [x for x in temp_line.split(" ") if x != "" and x != " "]
                    try:
                        score = float(temp_line[0])
                    except:
                        raise Exception("Score not in remark line for {}".format(vina_file))

                    if max_docking_score is not None:
                        if score > max_docking_score:
                            terminate_run = True
                            write_pose = False
                            break


                    if min_docking_score is not None:
                        if score < min_docking_score:
                            # This score is bellow the minimum but the
                            # poses after may not be docked as well.
                            # Normally this should be a stop but may
                            # be useful for studying poor poses...
                            write_pose = False

                    printout_list.append(line)

            elif "ENDMDL" in line:
                if write_pose is True:
                    printout_list.append(line)

                    # convert list of pdbqt info to
                    # .pdb format by removing the partial charge info in ATOM line
                    printout_pdb = convert_pdbqt_to_pdb(printout_list)

                    # write to a file
                    outfile = output_folder + os.sep + short_name +\
                        "_pose_{}.pdb".format(pose_number)

                    with open(outfile, "w") as f:
                        f.write(printout_pdb)

                # Reset variables for the next iteration

                printout_list = []
                score = 0.0
                terminate_run = False
                write_pose = True
                printout_pdb = ""

                # update the counter of the pose number
                pose_number += 1
            else:
                printout_list.append(line)
#
def convert_pdbqt_to_pdb(list_of_lines):
    """
    Converts lines from a pdbqt.vina file pose to pdb format.
    Inputs:
    :param list list_of_lines: list of lines of a docked pdbqt pose
    Returns:
    :returns: str printout: A string for a .pdb to write to a file
    """
    printout = ""
    line_index_range = [x for x in range(0, 61)] + [x for x in range(70, 80)]

    for line in list_of_lines:
        if "ATOM" in line or "HETATM" in line:
            short_line = ""
            for i in line_index_range:
                # print(i)
                if i >= len(line):
                    continue

                short_line = short_line + line[i]

            printout = printout + short_line
        elif "REMARK                            x       y       z     vdW  Elec" + \
            "       q    Type" in line \
            or "REMARK                         _______ _______ _______ _____ _____" + \
            "    ______ ____" in line:
            short_line = ""
            for i in line_index_range:
                # print(i)
                if i >= len(line):
                    continue

                short_line = short_line + line[i]

            printout = printout + short_line + "\n"
        else:
            printout = printout + line
    return printout
#
def start_run_main(vars):
    """
    This will run the main arguments for the script.

    Inputs:
    :param dict vars: dictionary of user variables.
    """
    output_folder = vars["output_folder"] + os.sep
    max_num_of_poses = vars["max_num_of_poses"]
    max_docking_score = vars["max_docking_score"]
    min_docking_score = vars["min_docking_score"]



    vina_docked_pdbqt_file = vars["vina_docked_pdbqt_file"]
    if os.path.isfile(vina_docked_pdbqt_file) is True:

        run_conversion_for_a_vina_file(vina_docked_pdbqt_file, output_folder,
                                       max_num_of_poses, max_docking_score,
                                       min_docking_score)

    else:

        # vina_docked_pdbqt_file is a folder run for all .pdbqt.vina files
        pdbqt_files = glob.glob(vina_docked_pdbqt_file + "*.pdbqt.vina")
        pdbqt_files.extend(glob.glob(vina_docked_pdbqt_file + "*.PDBQT.vina"))
        pdbqt_files.extend(glob.glob(vina_docked_pdbqt_file + "*.pdbqt.VINA"))
        pdbqt_files.extend(glob.glob(vina_docked_pdbqt_file + "*.PDBQT.VINA"))
        pdbqt_files = list(set(pdbqt_files))
        if len(pdbqt_files) == 0:
            printout = "No .pdbqt.vina were found at: {}".format(vina_docked_pdbqt_file)
            raise Exception(printout)
        job_input = tuple([tuple([vina_docked_pdbqt_file, output_folder, max_num_of_poses,
                                  max_docking_score, min_docking_score]) \
                                 for vina_docked_pdbqt_file in pdbqt_files])
        # run convert in multithread
        mol_usable_list = mp.multi_threading(job_input, -1, run_conversion_for_a_vina_file)
#
def get_arguments_from_argparse(args_dict):
    """
    This function handles the arg parser arguments for the script.

    Inputs:
    :param dict args_dict: dictionary of parameters
    Returns:
    :returns: dict args_dict: dictionary of parameters
    """


    # Argument handling
    if type(args_dict["vina_docked_pdbqt_file"]) != str:
        raise Exception("provided vina_docked_pdbqt_file must be either a docked vina file or \
            a directory of docked vina files.")

    if type(args_dict["output_folder"]) != str:
        raise Exception("provided output_folder must be a directory.")

    #  argument_handling
    if os.path.exists(args_dict["vina_docked_pdbqt_file"]) is False:
        raise Exception("provided vina_docked_pdbqt_file can not be found.")
    if ".pdbqt.vina" not in args_dict["vina_docked_pdbqt_file"].lower():
        if os.path.isdir(args_dict["vina_docked_pdbqt_file"]) is False:
            raise Exception("provided vina_docked_pdbqt_file must be either a docked vina file \
            containing .pdbqt.vina in file name or a directory of docked vina files.")

        args_dict["vina_docked_pdbqt_file"] = args_dict["vina_docked_pdbqt_file"] + os.sep

    if os.path.exists(args_dict["output_folder"]) is False:
        try:
            os.mkdir(args_dict["output_folder"])
        except:
            pass
        if os.path.exists(args_dict["output_folder"]) is False:
            raise Exception("output_folder could not be made or found.")
    else:
        if os.path.isdir(args_dict["output_folder"]) is False:
            raise Exception("output_folder needs to be a directory.")

        args_dict["output_folder"] = os.path.abspath(args_dict["output_folder"]) + os.sep

    # handle max_num_of_poses
    if args_dict["max_num_of_poses"] is not None:
        if type(args_dict["max_num_of_poses"]) != float and \
            type(args_dict["max_num_of_poses"]) != int:
            raise Exception("max_num_of_poses must be a int or None")
        if type(args_dict["max_num_of_poses"]) == float:
            args_dict["max_num_of_poses"] = int(args_dict["max_num_of_poses"])
    # handle max_docking_score
    if args_dict["max_docking_score"] is not None:
        if type(args_dict["max_docking_score"]) != float or \
            type(args_dict["max_docking_score"]) != int:
            raise Exception("max_docking_score must be a float or None")
        if type(args_dict["max_docking_score"]) == int:
            args_dict["max_docking_score"] = float(args_dict["max_num_of_poses"])
    # handle min_docking_score
    if args_dict["min_docking_score"] is not None:
        if type(args_dict["min_docking_score"]) != float or \
            type(args_dict["min_docking_score"]) != int:
            raise Exception("min_docking_score must be a float or None")
        if type(args_dict["min_docking_score"]) == int:
            args_dict["min_docking_score"] = float(args_dict["max_num_of_poses"])

    return args_dict
#


# Argument parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument(
    '--vina_docked_pdbqt_file', '-f',
    required=True, default=None,
    help='Path to .pdbqt.vina file to split into 1 .pdb file per pose that matches all criteria. \
    if this is a directory it will convert all of the files with the extension .pdbqt.vina'
)

PARSER.add_argument(
    '--output_folder', '-o', type=str, default=None,
    help='Path to folder where the .pdb files will be placed. \
    Files will be the basename of the docked file with _pose_{pose_number}.pdb \
    replacing the extension .pdbqt.vina.'
)
PARSER.add_argument(
    '--max_num_of_poses', type=int, required=False, default=-1,
    help='Each docked file will have 1 or more poses of the ligand. This setting \
    controls how many are converted. default is -1 which means all poses possible. \
    max_num_of_poses=1 means only the best docked pose will be converted. \
    If additional criteria like max_docking_score is applied a pose must meet both criteria \
    to be converted. ie) if max_num_of_poses= 5 and max_docking_score=-13.0 \
    for a pose to be converted it must be between the 1st and 5th pose in the file and \
    must have docked with a score less than or equal to -13.0.'
)

PARSER.add_argument(
    '--max_docking_score', type=float, required=False, default=None,
    help='The most positive docking score to be converted. (More negative scores \
    are predicted to bind better). If additional criteria such as \
    max_num_of_poses is applied a pose must meet both criteria \
    to be converted. ie) if max_num_of_poses= 5 and max_docking_score=-13.0 \
    for a pose to be converted it must be between the 1st and 5th pose in the file and \
    must have docked with a score less than or equal to -13.0.'
)

PARSER.add_argument(
    '--min_docking_score', type=float, required=False, default=None,
    help='The most negative docking score to be converted. (More negative scores \
    are predicted to bind better). If additional criteria such as \
    max_num_of_poses is applied a pose must meet both criteria \
    to be converted. ie) if min_docking_score= -15.0 and max_docking_score=-13.0 \
    for a pose to be converted it must: \
    -13.0. <= docking score <= -15.0'
)

PARSER.add_argument(
    "--number_of_processors",
    "-p",
    type=int,
    metavar="N",
    default=-1,
    help="Number of processors to use for parallel calculations.\
         Set to -1 for all available CPUs."
)


ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)

start_run_main(ARGS_DICT)
print("finished")
