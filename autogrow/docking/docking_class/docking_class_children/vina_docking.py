"""
The child classes from ParentExample
"""
import __future__

import os
import sys
import glob

import autogrow.docking.delete_failed_mol as Delete
import autogrow.docking.ranking.ranking_mol as Ranking
from autogrow.docking.docking_class.parent_dock_class import ParentDocking
import autogrow.docking.scoring.execute_scoring_mol as Scoring


class VinaDocking(ParentDocking):
    """
    RUN VINA DOCKING

    Inputs:
    :param class ParentDocking: Parent docking class to inherit from
    """

    def __init__(self, vars=None, receptor_file=None,
                 file_conversion_class_object=None, test_boot=True):
        """
        get the specifications for Vina/QuickVina2 from vars load them into
        the self variables we will need and convert the receptor to the proper
        file format (ie pdb-> pdbqt)

        Inputs:
        :param dict vars: Dictionary of User variables
        :param str receptor_file: the path for the receptor pdb
        :param obj file_conversion_class_object: object which is used to
            convert files from pdb to pdbqt
        :param bool test_boot: used to initialize class without objects for
            testing purpose
        """

        if test_boot is False:

            self.vars = vars
            self.debug_mode = vars["debug_mode"]
            self.file_conversion_class_object = file_conversion_class_object

            # VINA SPECIFIC VARS
            receptor_file = vars["filename_of_receptor"]
            # mgl_python = vars["mgl_python"]
            # receptor_template = vars["prepare_receptor4.py"]
            # number_of_processors = vars["number_of_processors"]
            # docking_executable = vars["docking_executable"]

            ###########################

            self.receptor_pdbqt_file = receptor_file + "qt"

            self.vars["docking_executable"] = self.get_docking_executable_file(
                self.vars
            )

    def run_ligand_handling_for_docking(self, pdb_file):
        """
        this function converts the ligands from PDB to PDBQT format. Returns
        NONE if it worked and the name if it failed to convert.

        Inputs:
        :param str pdb_file: the pdb file of a ligand to format, dock and
            score

        Returns:
        :returns: str smile_name: name of smiles if it failed to dock returns
            None if it docked properly
        """

        # convert ligands to pdbqt format
        # log("\nConverting ligand PDB files to PDBQT format...")
        did_it_convert, smile_name = self.file_conversion_class_object.convert_ligand_pdb_file_to_pdbqt(pdb_file)

        if did_it_convert is False:
            # conversion failed
            return smile_name
        # Conversion pass. Return None
        # only return failed smile_names which will be handled later
        return None

    def run_dock(self, pdbqt_filename):
        """
        this function runs the docking. Returns None if it worked and the name
        if it failed to dock.

        Inputs:
        :param str pdbqt_filename: the pdbqt file of a ligand to dock and
            score

        Returns:
        :returns: str smile_name: name of smiles if it failed to dock returns
            None if it docked properly
        """

        # print(" -------- Docking compounds using AutoDock Vina...")
        self.dock_ligand(pdbqt_filename)

        # check that it docked
        pdb_filename = pdbqt_filename.replace("qt", "")

        did_it_dock, smile_name = self.check_docked(pdb_filename)

        if did_it_dock is False:
            # Docking failed

            if smile_name is None:
                # print("Missing pdb and pdbqt files for : ", pdbqt_filename)
                pass 

            return smile_name

        return None

    #######################################
    # STUFF DONE BY THE INIT
    ##########################################
    def get_docking_executable_file(self, vars):
        """
        This retrieves the docking executable files Path.

        Inputs:
        :param dict vars: Dictionary of User variables

        Returns:
        :returns: str docking_executable: String for the docking executable
            file path
        """

        if vars["docking_executable"] is None:
            # get default docking_executable for vina
            script_dir = str(os.path.dirname(os.path.realpath(__file__)))
            docking_executable_directory = (
                script_dir.split(os.sep + "docking_class")[0]
                + os.sep
                + "docking_executables"
                + os.sep
            )

            if sys.platform == "linux" or sys.platform == "linux2":
                # Use linux version of Autodock Vina
                docking_executable = (
                    docking_executable_directory
                    + "vina"
                    + os.sep
                    + "autodock_vina_1_1_2_linux_x86"
                    + os.sep
                    + "bin"
                    + os.sep
                    + "vina"
                )

            elif sys.platform == "darwin":
                # Use OS X version of Autodock Vina
                docking_executable = (
                    docking_executable_directory
                    + "vina"
                    + os.sep
                    + "autodock_vina_1_1_2_mac"
                    + os.sep
                    + "bin"
                    + os.sep
                    + "vina"
                )

            elif sys.platform == "win32":
                # Windows...
                raise Exception("Windows is currently not supported")
            else:
                raise Exception("This OS is currently not supported")

        else:
            # if user specifies a different vina executable
            docking_executable = vars["docking_executable"]

        if os.path.exists(docking_executable) is False:
            printout = "Docking executable could not be found at: "
            printout = printout + "{}".format(docking_executable)
            print(printout)
            raise Exception(printout)

        return docking_executable

    # Finding PDBs for ligands in a folder
    def find_pdb_ligands(self, current_generation_pdb_dir):
        """
        This finds all the pdb files of ligands in a directory

        Inputs:
        :param str current_generation_pdb_dir: the dir path which contains the
            pdb files of ligands to be converted

        Returns:
        :returns: list pdbs_in_folder: a list of all PDB's in the dir
        """

        # make list of every pdb in the current generations pdb folder
        pdbs_in_folder = []
        for filename in glob.glob(current_generation_pdb_dir + "*.pdb"):
            pdbs_in_folder.append(filename)

        return pdbs_in_folder

    # Find ligands which converted to PDBQT
    def find_converted_ligands(self, current_generation_pdb_dir):
        """
        This finds all the pdbqt files of ligands in a directory

        Inputs:
        :param str current_generation_pdb_dir: the dir path which contains the
            pdbqt files of ligands to be docked

        Returns:
        :returns: list pdbqts_in_folder: a list of all PDBqt's in the dir
        """

        # make list of every pdbqt in the current generations pdb folder
        pdbqts_in_folder = []
        for filename in glob.glob(current_generation_pdb_dir + "*.pdbqt"):
            pdbqts_in_folder.append(filename)

        return pdbqts_in_folder

    #######################################
    # DOCK USING VINA                     #
    #######################################
    def dock_ligand(self, lig_pdbqt_filename):
        """
        Dock the ligand pdbqt files in a given directory using AutoDock Vina

        Inputs:
        :param str lig_pdbqt_filename: the ligand pdbqt filename
        """
        vars = self.vars
        timeout_option = vars["timeout_vs_gtimeout"]
        docking_timeout_limit = vars["docking_timeout_limit"]
        # do the docking of the ligand Run with a timeout_option limit.
        # Default setting is 5 minutes. This is excessive as most things run
        # within 30seconds This will prevent stalling out. timeout or gtimeout
        torun = (
            "{} {} {}".format(timeout_option, docking_timeout_limit, vars["docking_executable"])
            + " --center_x "
            + str(vars["center_x"])
            + " --center_y "
            + str(vars["center_y"])
            + " --center_z "
            + str(vars["center_z"])
            + " --size_x "
            + str(vars["size_x"])
            + " --size_y "
            + str(vars["size_y"])
            + " --size_z "
            + str(vars["size_z"])
            + " --receptor "
            + self.receptor_pdbqt_file
            + " --ligand "
            + lig_pdbqt_filename
            + " --out "
            + lig_pdbqt_filename
            + ".vina --cpu 1"
        )

        # Add optional user variables additional variable
        if (
                vars["docking_exhaustiveness"] is not None
                and vars["docking_exhaustiveness"] != "None"
        ):
            if (
                    type(vars["docking_exhaustiveness"]) == int
                    or type(vars["docking_exhaustiveness"]) == float
            ):
                torun = (
                    torun
                    + " --exhaustiveness "
                    + str(int(vars["docking_exhaustiveness"]))
                )
        if vars["docking_num_modes"] is not None and vars["docking_num_modes"] != "None":
            if (
                    type(vars["docking_num_modes"]) == int
                    or type(vars["docking_num_modes"]) == float
            ):
                torun = torun + " --num_modes " + str(int(vars["docking_num_modes"]))

        # Add output line MUST ALWAYS INCLUDE THIS LINE
        torun = (
            torun
            + " >>"
            + lig_pdbqt_filename
            + "_docking_output.txt "
            + " 2>>"
            + lig_pdbqt_filename
            + "_docking_output.txt"
        )

        # print("\tDocking: {}".format(lig_pdbqt_filename))
        results = self.execute_docking_vina(torun)

        if results is None or results is None or results == 256:
            made_changes = self.replace_atoms_not_handled_by_forcefield(
                lig_pdbqt_filename
            )
            if made_changes is True:
                results = self.execute_docking_vina(torun)
                if results == 256 or results is None:
                    # print(
                    #     "\nLigand failed to dock after corrections: {}\n".format(
                    #         lig_pdbqt_filename
                    #     )
                    # )
                    pass 
            else:
                # print("\tFinished Docking: {}".format(lig_pdbqt_filename))
                pass 
        else:
            pass 
            # print("\tFinished Docking: {}".format(lig_pdbqt_filename))

    def replace_atoms_not_handled_by_forcefield(self, lig_pdbqt_filename):
        """
        Replaces atoms not handled by the forcefield to prevent errors. Atoms
        include B and Si.

        Inputs:
        :param str lig_pdbqt_filename: the ligand pdbqt filename

        Returns:
        :returns: bool retry: If True it will be ligand will be redocked, if
            False its dones and wont be docked again.
        """

        # VINA/QuickVINA and MGL have problems with the forcefields for
        # certain atom types To correct this, Autodock Vina suggests replacing
        # the
        atoms_to_replace = [
            "B \n",
            "B\n",
            "Si \n",
            "Si\n",
        ]  # add the \n at the end so we replace the end portion of the line
        printout_of_file = ""
        printout_info = ""
        retry = False
        line_count = 0
        with open(lig_pdbqt_filename, "r") as f:
            for line in f.readlines():
                line_count = line_count + 1
                if "HETATM" in line:
                    for x in atoms_to_replace:
                        if x in line:
                            line = line.replace(x, "A \n")
                            retry = True

                            printout_info = (
                                printout_info
                                + "Changing '{}' to 'A ' in line: {} of {}".format(
                                    str(x.strip()), line_count, lig_pdbqt_filename
                                )
                            )  # x Need to remove whitespaces on both ends
                printout_of_file = printout_of_file + line

        if retry is True:
            print(printout_info)
            with open(lig_pdbqt_filename, "w") as f:
                f.write(printout_of_file)
        else:
            printout_info = "\nCheck the docking message for 'Parse error on'"
            printout_info = (
                printout_info
                + "\n\t This ligand failed to dock. Please check that all "
                + "atoms are covered by the docking forcefield"
            )
            printout_info = (
                printout_info
                + "\n\t Any atoms not covered by the forcefield should be "
                + "added to atoms_to_replace in the function "
                + "replace_atoms_not_handled_by_forcefield"
            )
            printout_info = printout_info + "\n\t Verify for this ligand: {}\n".format(
                lig_pdbqt_filename
            )
            print(printout_info)
        return retry

    def execute_docking_vina(self, command):
        """
        Run a single docking execution command

        Inputs:
        :param str command: string of command to run.

        Returns:
        :returns: int result: the exit output for the command. If its None of
            256 it failed.
        """

        try:
            result = os.system(command)
        except:
            result = None
            print("Failed to execute: " + command)
        return result

    def check_docked(self, pdb_file):
        """
        given a pdb_file name, test if a pdbqt.vina was created. If it failed
        to dock delete the file pdb and pdbqt file for that ligand -then
        return false

        if it docked properly return True

        Inputs:
        :param str pdb_file: pdb file path

        Returns:
        :returns: bool bool: false if not vina was unsuccessful
        :returns: str smile_name: name of the pdb file
        """

        if not os.path.exists(pdb_file):
            # PDB file doesn't exist
            return False, None
        smile_name = self.file_conversion_class_object.get_smile_name_from_pdb(
            pdb_file
        )
        if not os.path.exists(pdb_file + "qt.vina"):
        # so this pdbqt.vina file didn't exist
            if self.debug_mode is False:
                # print("Docking unsuccessful: Deleting "
                #         + os.path.basename(pdb_file) + "...")

                # REMOVE Failed molecules. Delete ones that were not docked
                # successfully
                Delete.delete_all_associated_files(pdb_file)
                # # delete pdbqt_file
                pdbqt_file = pdb_file + "qt"
                Delete.delete_all_associated_files(pdbqt_file)

                return False, smile_name

            # Failed to dock but in debug mode
            # print("Docking unsuccessful: " + os.path.basename(pdb_file) + "...")
            return False, smile_name

        # Successfully docked
        return True, smile_name

    ##########################################
    # Convert the dock outputs to a usable formatted .smi file
    # This is mandatory for all Docking classes but
    # implementation and approach varies by docking and scoring choice
    ##########################################

    def rank_and_save_output_smi(self, vars, current_generation_dir,
                                 current_gen_int, smile_file,
                                 deleted_smiles_names_list):
        """
        Given a folder with PDBQT's, rank all the SMILES based on docking
        score (High to low). Then format it into a .smi file. Then save the
        file.

        Inputs:
        :param dict vars: vars needs to be threaded here because it has the
            paralizer object which is needed within Scoring.run_scoring_common
        :param str current_generation_dir: path of directory of current
            generation
        :param int current_gen_int: the interger of the current generation
            indexed to zero
        :param str smile_file:  File path for the file with the ligands for
            the generation which will be a .smi file
        :param list deleted_smiles_names_list: list of SMILES which may have
            failed the conversion process

        Returns:
        :returns: str output_ranked_smile_file: the path of the output ranked
            .smi file
        """

        # Get directory string of PDB files for Ligands
        folder_with_pdbqts = current_generation_dir + "PDBs" + os.sep

        # Run any compatible Scoring Function
        smiles_list = Scoring.run_scoring_common(vars, smile_file, folder_with_pdbqts)

        # Before ranking these we need to handle Pass-Through ligands from the
        # last generation If it's current_gen_int==1 or if
        # vars['redock_elite_from_previous_gen'] is True -Both of these states
        # dock all ligands from the last generation so all of the pass-through
        # lig are already in the PDB's folder thus they should be accounted
        # for in smiles_list If vars['redock_elite_from_previous_gen'] is False
        # and current_gen_int != 1 - We need to append the scores form the
        # last gen to smiles_list

        # Only add these when we haven't already redocked the ligand
        if (
                self.vars["redock_elite_from_previous_gen"] is False
                and current_gen_int != 0
        ):
            # Go to previous generation folder
            prev_gen_num = str(current_gen_int - 1)
            run_folder = self.vars["output_directory"]
            previous_gen_folder = run_folder + "generation_{}{}".format(
                str(prev_gen_num), os.sep
            )
            ranked_smi_file_prev_gen = (
                previous_gen_folder
                + "generation_{}_ranked.smi".format(str(prev_gen_num))
            )

            # Also check sometimes Generation 1 won't have a previous
            # generation to do this with and sometimes it will
            if (
                    current_gen_int == 1
                    and os.path.exists(ranked_smi_file_prev_gen) is False
            ):
                pass
            else:
                print("Getting ligand scores from the previous generation")

                # Shouldn't happen but to be safe.
                if os.path.exists(ranked_smi_file_prev_gen) is False:
                    raise Exception(
                        "Previous generation ranked .smi file does not exist. "
                        + "Check if output folder has been moved"
                    )

                # Get the data for all ligands from previous generation ranked
                # file
                prev_gen_data_list = Ranking.get_usable_format(ranked_smi_file_prev_gen)

                # Get the list of pass through ligands
                current_gen_pass_through_smi = (
                    current_generation_dir
                    + "SeedFolder{}Chosen_Elite_To_advance_Gen_{}.smi".format(
                        os.sep, str(current_gen_int)
                    )
                )
                pass_through_list = Ranking.get_usable_format(
                    current_gen_pass_through_smi
                )

                # Convert lists to searchable Dictionaries.
                prev_gen_data_dict = Ranking.convert_usable_list_to_lig_dict(
                    prev_gen_data_list
                )

                pass_through_data = []
                for lig in pass_through_list:
                    smile_plus_id = str(lig[0] + lig[1])
                    lig_data = prev_gen_data_dict[smile_plus_id]
                    lig_info_remove_diversity = [
                        lig_data[x] for x in range(0, len(lig_data) - 1)
                    ]
                    pass_through_data.append(lig_info_remove_diversity)

                smiles_list.extend(pass_through_data)

        # Output format of the .smi file will be: SMILES    Full_lig_name
        # shorthandname   ...AnyCustominfo... Fitness_metric  diversity
        # Normally the docking score is the fitness metric but if we use a
        # Custom metric than dock score gets moved to index -3 and the new
        # fitness metric gets -2

        # sort list by the affinity of each sublist (which is the last index
        # of sublist)
        smiles_list.sort(key=lambda x: float(x[-1]), reverse=False)

        # score the diversity of each ligand compared to the rest of the
        # ligands in the group this adds on a float in the last column for the
        # sum of pairwise comparisons the lower the diversity score the more
        # unique a molecule is from the other mols in the same generation
        smiles_list = Ranking.score_and_append_diversity_scores(smiles_list)

        # name for the output file
        output_ranked_smile_file = smile_file.replace(".smi", "") + "_ranked.smi"

        # save to a new output smiles file. ie. save to ranked_smiles_file

        with open(output_ranked_smile_file, "w") as output:
            for ligand_info_list in smiles_list:
                str_ligand_info_list = [str(x) for x in ligand_info_list]
                output_line = "\t".join(str_ligand_info_list) + "\n"
                output.write(output_line)

        return output_ranked_smile_file
