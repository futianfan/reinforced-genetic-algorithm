"""
This function handles the scoring/rescoring of docked molecules.
"""

import __future__

from autogrow.docking.docking_class.get_child_class import get_all_subclasses

# importing scoring_functions is necessary to find rescoring modules
import autogrow.docking.scoring.scoring_classes.scoring_functions
from autogrow.docking.scoring.scoring_classes.parent_scoring_class import ParentScoring


def pick_run_class_dict(scoring_choice):
    """
    This will retrieve all the names of every child class of the parent class
    ParentScoring

    Inputs:
    :param list scoring_choice: List with the User specified scoring choices

    Returns:
    :returns: object child_dict[scoring_choice]: the class for running the
        chosen scoring method
    """

    children = get_all_subclasses(ParentScoring)

    child_dict = {}
    for child in children:
        child_name = child.__name__
        child_dict[child_name] = child
    return child_dict[scoring_choice]


############
############
def run_scoring_common(vars, smile_file, folder_to_search):
    """
    This section runs the functions common to all scoring functions.

    IF ONE INCORPORATES A NEW DOCKING OR SCORING SOFTWARE, CONFIRM THAT ITS
    INPUT/OUTPUTS CONFORM TO THIS SECTION.

    ############## VERY IMPORTANT SECTION########################

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param str smile_file: File path for the file with the ligands for the
        generation which will be a .smi file
    :param str folder_to_search: a directory to search containing docked
        ligands

    Returns:
    :returns: dict lig_dict: a dictionary where the keys are the ligand
        shorthand names and the items are a list of any information about the
        ligand with the fitness measure as the -1 idx in each list
    """

    # Retrieve a list of all files with the proper information within
    # folder_to_search
    scoring_choice = vars["scoring_choice"]
    scoring_class = pick_run_class_dict(scoring_choice)

    # Make a dictionary of all ligands which may have been docked. This will
    # be taken from the .smi file gen_n_to_convert.smi file Keys are the
    # shortened name of the ligand id and the item is the SMILES string and
    # the ligand id full name
    smiles_dict = make_dict_of_smiles(smile_file)
    # print('------------ smiles_dict', smiles_dict)
    # Use a temp vars dict so you don't put mpi multiprocess info through
    # itself...
    temp_vars = {}
    for key in list(vars.keys()):
        if key == "parallelizer":
            continue
        temp_vars[key] = vars[key]

    # Initialize the scoring class
    scoring_object = scoring_class(temp_vars, smiles_dict, test_boot=False)

    # Find all the files that need to be scored. Because this depends on the
    # scoring function these may be different file types
    files_to_score = scoring_object.find_files_to_score(folder_to_search)
    # print('---------- files_to_score', files_to_score)
    # Run Rescoring If applicable (All classes should have this even if its
    # just returning None)
    files_to_score = run_rescoring(vars, scoring_object, files_to_score)
    files_to_score = [x for x in files_to_score if x is not None]

    # Determine if the values from the Scoring function should be adjusted by
    # the number of non-H atoms in the ligand ie) rescoring_lig_efficiency
    rescore_lig_efficiency = vars["rescore_lig_efficiency"]

    # create rescore_lig_efficiency object if needed
    if rescore_lig_efficiency is True:
        rescore_lig_efficiency_class = pick_run_class_dict("LigEfficiency")
        # Initialize the scoring class
        rescore_lig_efficiency_scoring_object = rescore_lig_efficiency_class(
            temp_vars, smiles_dict, test_boot=False
        )
    else:
        rescore_lig_efficiency_scoring_object = None

    job_input_files_to_score = tuple(
        [
            tuple(
                [
                    scoring_object,
                    file_path,
                    rescore_lig_efficiency,
                    rescore_lig_efficiency_scoring_object,
                ]
            )
            for file_path in files_to_score
        ]
    )

    # Format for list_of_raw_data must be [lig_id_shortname, any_details,
    # fitness_score_to_use]
    list_of_list_of_lig_data = vars["parallelizer"].run(
        job_input_files_to_score, score_files_multithread
    )

    # Convert all list_of_list_of_lig_data to a searchable dictionary This
    # removes redundancies from multiple conformations and any Nones. This is
    # left as a dictionary because it could be an area for expansion.
    lig_dict = make_lig_score_dictionary(list_of_list_of_lig_data)

    # Convert Dictionary to the final list form
    list_of_smiles_with_scores = []
    for key in list(lig_dict.keys()):

        lig_info = lig_dict[key]
        if lig_info is None:
            continue
        lig_info = [str(x) for x in lig_info]

        # Make sure every value is a string

        list_of_smiles_with_scores.append(lig_info)

    # Make sure every value is

    return list_of_smiles_with_scores


#
def run_rescoring(vars, scoring_object, files_to_score):
    """
    Run a rescoring function on docked ligands.

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param object scoring_object: the class for running the chosen scoring
        method
    :param list files_to_score: list of files to rescores

    Returns:
    :returns: list completed_rescore: a list of all ligands which passed a
        scoring function.
    """

    files_to_score = [x for x in files_to_score if x is not None]

    # Run Rescoring If applicable (All classes should have this even if its
    # just returning None)
    job_input_files_to_score = tuple(
        [tuple([file_path, scoring_object]) for file_path in files_to_score]
    )

    # Format for list_of_raw_data must be [lig_id_shortname, any_details,
    # fitness_score_to_use]
    results_rescore = vars["parallelizer"].run(
        job_input_files_to_score, rescore_single_file
    )

    if len(results_rescore) == 0:
        return files_to_score
    if results_rescore[0] == "Not Applicable":
        return files_to_score

    results_rescore = [x for x in results_rescore if x is not None]
    completed_rescore = [x[0] for x in results_rescore if x[1] is True]
    failed_to_rescore = [x[0] for x in results_rescore if x[1] is False]

    # print fails which made it through the try statement but failed to
    # produce an output file.
    if len(failed_to_rescore) != 0:
        print("The following files failed to be rescored: ")
        print(failed_to_rescore)
    else:
        print("All rescoring attempts were successful")
        pass 

    if len(completed_rescore) == 0:
        printout = (
            "All Rescoring attempts failed to create output files. No data to report."
        )
        raise Exception(printout)

    print("Finished rescoring")
    print("######################\n")

    return completed_rescore


#
def rescore_single_file(file_path, scoring_object):
    """
    Run scoring_object.run_rescoring through this function so multithread
    doesn't break.

    Inputs:
    :param str file_path: Path to a vina output file to be rescored
    :param object scoring_object: object that rescores such as an NN1 or NN2
        class object

    Returns:
    :returns: list results of a rescoring function: [file_path, it_rescored]
        [PATH, True] means it passed [PATH, False] means it failed a results of
        all NN1 files
    """

    return scoring_object.run_rescoring(file_path)


def make_lig_score_dictionary(list_of_list_of_lig_data):
    """
    Given a list of ligands with the scoring data make a dictionary.

    This will also reduce multiple Conformers and Poses down to the top score.

    # REQUIRES THE BEST SCORE TO BE MOST NEGATIVE

    Inputs:
    :param list list_of_list_of_lig_data: a list of lists containing all info
        on ligands after scoring [[SMILES, id, Shortid, additional_info...,
        fitness_score], [SMILES, id, Shortid, additional_info..., fitness_score]]

    Returns:
    :returns: dict lig_dict: a dictionary containing the information of all
        ligands from list_of_list_of_lig_data this dictionary has the short_id as
        the key and the item is a list from list_of_list_of_lig_data for that
        ligand. Because there can be multiple files for multiple conformations of
        a given ligand, this will reduce down multiple confirmations to a single
        ligand with the most negative fitness score.
    """

    lig_dict = {}
    for lig in list_of_list_of_lig_data:
        if lig is None:
            continue
        lig_short_id = str(lig[2])
        fitness_score = float(lig[-1])

        # Handle if Multiple Conformers and Poses
        if lig_short_id in lig_dict.keys():
            if float(lig_dict[lig_short_id][-1]) > float(fitness_score):
                lig_dict[lig_short_id] = lig
            else:
                continue
        else:
            lig_dict[lig_short_id] = lig

    return lig_dict


def score_files_multithread(scoring_object, file_path, rescore_lig_efficiency,
                            lig_efficiency_scoring_object):
    """
    Run the scoring of a single molecule.

    Inputs:
    :param object scoring_object: the class for running the chosen scoring
        method
    :param str file_path: the path to the file to be scored
    :param bol rescore_lig_efficiency: if True than the final score will be
        adjusted to the ligand efficieny score, otherwise it will remain the
        output of the scoring_object
    :param object lig_efficiency_scoring_object: the class for running the
        rescoring by ligand efficieny

    Returns:
    :returns: list list_of_lig_data: information about the scored ligand.
        Score is last index (ie. [SMILES, lig_id, lig_id_shortname, any_details,
        fitness_score_to_use] )
    """

    list_of_lig_data = scoring_object.run_scoring(file_path)
    if rescore_lig_efficiency is True:

        list_of_lig_data = lig_efficiency_scoring_object.get_lig_efficiency_rescore_from_a_file(
            file_path, list_of_lig_data
        )

    return list_of_lig_data


############
############


def make_dict_of_smiles(smile_file):
    """
    This will take a .smi file and make a dictionary with all of the info
    about the smiles. This list won't have scores yet but will have all of the
    string info. This will be used later to search through.

    The keys will be the shorthand id for each ligand.

    Inputs:
    :param str smile_file: the path for the receptor pdb

    Returns:
    :return dict smiles_dict: a list of ligand info before docking
    """

    smiles_dict = {}
    # load smile file and convert to list with index
    with open(str(smile_file), "r") as smi:

        # index for this for-loop is zero based indexing while ligand naming
        # index is one-based-indexing so we will use index = index+1
        for index, line in enumerate(smi.readlines()):
            # split_line = line.replace("\n", "").split("\t")
            # ligand_name = line.replace("\n", "").split("\t")[1]  
            split_line = line.replace("\n", "").split()
            ligand_name = line.replace("\n", "").split()[1] 
            # ligand_name should be something like: ['CCC', '(ZINC123+ZINC345)Gen_0_Cross_99571'] or ['CCC', 'ZINC123']

            if len(ligand_name.split(")")) == 2:
                lig_name_short = ligand_name.split(")")[1]
            elif len(ligand_name.split(")")) == 1:
                lig_name_short = ligand_name

            smiles_dict[lig_name_short] = split_line

    return smiles_dict
