"""
Plots a line plot of the average score for each generation of AutoGrow run.

Example submit:
    python autogrow4/accessory_scripts/plot_autogrow_run.py\
        -i $PATH/Run_1/Run_0/ \
        --plot_reference_lines [['Olaparib Score',-12.8,'y'],\
            ['Niraparib',-10.7,'k'],['NAD/NADH',-10.3,'purple'],\
                ['ADP-ribose',-9.3,'maroon']]
"""
import __future__

import os
import glob
import json
import copy
import argparse

import matplotlib
import matplotlib.pyplot as plt

def get_usable_format(infile):
    """
    This code takes a string for an file which is formatted as an .smi file. It
    opens the file and reads in the components into a usable list.

    The .smi must follow the following format for each line:
        MANDATORY INFO
            part 1 is the SMILES string
            part 2 is the SMILES name/ID

        Optional info
            part -1 (the last piece of info) is the SMILES diversity score
                relative to its population
            part -2 (the second to last piece of info) is the fitness metric
                for evaluating
                - For default setting this is the Docking score
                - If you add a unique scoring function Docking score should be
                    -3 and that score function should be -2

            Any other information MUST be between part 2 and part -2 (this
            allows for the expansion of features without disrupting the rest of the code)

    Inputs:
    :param str infile: the string of the PATHname of a formatted .smi file to
        be read into the program

    Returns:
    :returns: list usable_list_of_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow
    """

    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_list_of_smiles = []

    if os.path.exists(infile) is False:
        print("\nFile of Source compounds does not exist: {}\n".format(infile))
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts separated by 4-spaces
            if len(parts) == 1:
                parts = line.split(
                    "    "
                )  # split line into parts separated by 4-spaces

            choice_list = []
            for i in range(0, len(parts)):
                choice_list.append(parts[i])
            usable_list_of_smiles.append(choice_list)

    return usable_list_of_smiles


def get_average_score_per_gen(infolder, folder_list):
    """
    This script will get the average docking score from the ranked .smi file
    from each generation.

    Inputs:
    :param str infolder: the path of the folder which has all of the
        generation folders
    :param list folder_list: a list of generation folders for each generation
        within infolder

    Returns:
    :returns: list usable_list_of_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow
    """

    average_affinity_dict = {}
    for gen_folder in folder_list:
        gen_folder_name = infolder + gen_folder + os.sep
        ranked_file = glob.glob(gen_folder_name + "*_ranked.smi")

        for rank_file in ranked_file:
            # write as a tab delineated .smi file
            with open(rank_file, "r") as f:
                gen_affinity_sum = float(0.0)
                num_lines_counter = float(0.0)
                for line in f:
                    line = line.replace("\n", "")
                    # split line into parts separated by 4-spaces
                    parts = line.split("\t")

                    choice_list = []
                    for i in range(0, len(parts)):
                        choice_list.append(parts[i])

                    gen_affinity_sum = gen_affinity_sum + float(choice_list[-2])
                    num_lines_counter = num_lines_counter + float(1.0)

            gen_affinity_average = gen_affinity_sum / num_lines_counter


            gen_num = os.path.basename(rank_file).split("_")[1]
            gen_name = "generation_{}".format(gen_num)
            average_affinity_dict[gen_name] = gen_affinity_average

    print_gens(average_affinity_dict)
    return average_affinity_dict

def get_average_top_score_per_gen(infolder, folder_list, top_score_per_gen):
    """
    This script will get the average docking score of the top N number of
    ligands ranked .smi file from each generation.

    Inputs:
    :param str infolder: the path of the folder which has all of the
        generation folders
    :param list folder_list: a list of generation folders for each generation
        within infolder
    :param int top_score_per_gen: the number of ligands to determine the
        average score. ie) if top_score_per_gen=50 it will return the average of
        the top 50 scores.

    Returns:
    :returns: dict average_affinity_dict: dictionary of average affinity
        scores for top_score_per_gen number of ligands
    """

    average_affinity_dict = {}

    for gen_folder in folder_list:
        gen_folder_name = infolder + gen_folder + os.sep
        ranked_file = glob.glob(gen_folder_name + "*_ranked.smi")

        for rank_file in ranked_file:
            # Check number of lines
            num_lines = 0
            with open(rank_file, "r") as rf:
                for line in rf:
                    num_lines = num_lines + 1

            if num_lines >= top_score_per_gen:
                # read as a tab delineated .smi file
                with open(rank_file, "r") as f:
                    gen_affinity_sum = float(0.0)

                    for i, line in enumerate(f.readlines()):
                        if i >= top_score_per_gen:
                            break
                        line = line.replace("\n", "")
                        parts = line.split(
                            "\t"
                        )  # split line into parts separated by 4-spaces

                        choice_list = []
                        for j in range(0, len(parts)):
                            choice_list.append(parts[j])

                        gen_affinity_sum = gen_affinity_sum + float(choice_list[-2])

                    gen_affinity_average = gen_affinity_sum / top_score_per_gen

                    gen_num = os.path.basename(rank_file).split("_")[1]
                    gen_name = "generation_{}".format(gen_num)
                    average_affinity_dict[gen_name] = gen_affinity_average

            else:
                gen_num = os.path.basename(rank_file).split("_")[1]
                gen_name = "generation_{}".format(gen_num)
                average_affinity_dict[gen_name] = "N/A"

    print_gens(average_affinity_dict)
    return average_affinity_dict

def print_gens(average_affinity_dict):
    """
    This prints out the average scores for each generation

    Inputs:
    :param dict average_affinity_dict: dictionary of average affinity scores
        for top_score_per_gen number of ligands
    """

    print("generation_number              average affinity score")
    affinity_keys = list(average_affinity_dict.keys())
    affinity_keys.sort(key=lambda x: int(x.split("_")[1]))
    for gen in affinity_keys:
        print(gen, "                  ", average_affinity_dict[gen])

def make_graph(dictionary):
    """
    Because some generations may not have 50 ligands this basically checks to see if
    theres enough ligands and prepares lists to be plotted

    Inputs:
    :param dict dictionary: dictionary of average affinity scores for
        top_score_per_gen number of ligands
    Returns:
    :returns: list list_generations: list of ints for each generation to be plotted.
        if a generation lacks ligands to generate the average it will return "N/A"
    :returns: list list_of_scores: list of averages for each generation;
        if a generation lacks ligands to generate the average it will return "N/A"
    """
    list_generations = []
    list_of_gen_names = []
    list_of_scores = []
    #print(dictionary)

    for key in dictionary.keys():
        #print(key)
        list_of_gen_names.append(key)

        score = dictionary[key]
        list_of_scores.append(score)

        gen = key.replace("generation_", "")

        gen = int(gen)
        list_generations.append(gen)
        list_of_gen_names.append(key)

    for i in list_of_scores:
        if i == "N/A":
            return None, None

    return list_generations, list_of_scores

def run_plotter(vars, dict_of_averages, outfile):
    """
    This plots the averages into a matplotlib figure. It will require you to
    answer questions about titles and labels

    Inputs:
    :param dict vars: dict of user variables which will govern how the
        programs runs
    :param dict dict_of_averages: a dictionary of dictionaries containing the
        average of each generation for the top 50,20, 10, and 1 ligand(s) and the
        overall average for each generation.
    :param str outfile: Path for the output file for the plot
    """

    average_affinity_dict = dict_of_averages["average_affinity_dict"]
    top_fifty_dict = dict_of_averages["top_fifty_dict"]
    top_twenty_dict = dict_of_averages["top_twenty_dict"]
    top_ten_dict = dict_of_averages["top_ten_dict"]
    top_one_dict = dict_of_averages["top_one_dict"]

    # print("Graphing Overall Average")
    list_generations_average, list_of_scores_average = make_graph(average_affinity_dict)
    # print("Graphing top_fifty_dict")
    print_fifty = True
    for key in top_fifty_dict.keys():
        if top_fifty_dict[key] == "N/A":
            print_fifty = False
    if print_fifty is True:
        list_generations_fifty, list_of_scores_fifty = make_graph(top_fifty_dict)
    # print("Graphing top_fifty_dict")
    print_twenty = True
    for key in top_twenty_dict.keys():
        if top_twenty_dict[key] == "N/A":
            print_twenty = False
    if print_twenty is True:
        list_generations_twenty, list_of_scores_twenty = make_graph(top_twenty_dict)

    # print("Graphing top_ten_dict")
    list_generations_ten, list_of_scores_ten = make_graph(top_ten_dict)
    # print("Graphing top_one_dict")
    list_generations_one, list_of_scores_one = make_graph(top_one_dict)
    # print("")

    ax = plt.subplot(111)

    ax.plot(
        list_generations_average, list_of_scores_average, color="b", label="Average"
    )
    if print_fifty is True:
        ax.plot(list_generations_fifty, list_of_scores_fifty, color="c", label="Top 50")

    if print_twenty is True:
        ax.plot(
            list_generations_twenty, list_of_scores_twenty, color="m", label="Top 20"
        )
    ax.plot(list_generations_ten, list_of_scores_ten, color="g", label="Top 10")
    ax.plot(list_generations_one, list_of_scores_one, color="r", label="Top 1")

    if vars["plot_reference_lines"] is not None:
        for ref_info in vars["plot_reference_lines"]:
            ax.axhline(y=ref_info[1], color=ref_info[2], linestyle=':', label=ref_info[0])

    ax.set_ylim()

    receptor_name = os.path.basename(vars["filename_of_receptor"])
    scoring_type = vars["scoring_choice"]
    docking_type = vars["scoring_choice"]
    num_lig = (
        int(vars["number_of_mutants"])
        + int(vars["number_of_crossovers"])
        + int(vars["number_elitism_advance_from_previous_gen"])
    )
    number_of_conf_per_lig = str(vars["max_variants_per_compound"])

    # Get Customizations
    title_of_figure = "{} Scores for {} using {}".format(
        scoring_type, receptor_name, docking_type
    )
    plt.title(title_of_figure, fontweight="semibold")

    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.274), fontsize="small")
    number_of_lig_per_gen = str(num_lig)

    output = (
        str(number_of_lig_per_gen)
        + " lig/gen"
        + "\n"
        + str(number_of_conf_per_lig)
        + " variants/lig"
    )

    plt.text(
        5.4, -8.5, output, bbox=dict(facecolor="white", alpha=0.5), fontsize="small"
    )

    # legend1 = plt.legend([lines[i].get_label() for i in range(0, lines_leg)],
    #           loc='center left', bbox_to_anchor=(1, 0.274),fontsize='small')
    # legend2 = plt.legend([output],loc='center left',
    #           bbox_to_anchor=(1, 0.774),fontsize='small')
    # # help(plt.legend)
    # ax.add_artist(legend1)
    # ax.add_artist(legend2)

    ax.set_ylim()

    if "VINA" in str(scoring_type):
        y_label = "Docking Affinity (kcal/mol"
    else:
        y_label = "Fitness Score"

    plt.ylabel(y_label, fontweight="semibold")

    plt.xlabel("Generation Number", fontweight="semibold")

    plt.savefig(outfile, bbox_inches="tight", \
        foramt=vars["outfile_format"], dpi=1000)

def print_data_table(infolder, folder_list):
    """
    This function takes a folder of an Autogrow Run and a list of all folders
    within the infolder, and finds the average of each generation, the average
    of the top 50,20, 10, and 1 ligand(s) in each generation.

    It prints the average docking score values in a table and returns that
    information as a dictionary of dictionaries.

    Inputs:
    :param str infolder: a string for the file path to a directory containing
        an Autogrow run. ie) "PATH/Run_0/"
    :param list folder_list: a list of every generation folders within the
        infolder

    Returns
    :returns: dict dict_of_averages: a dictionary of dictionaries containing
        the average of each generation for the top 50,20, 10, and 1 ligand(s) and
        the overall average for each generation.
    """

    print("Overall Scoring Average for all Compounds")
    average_affinity_dict = get_average_score_per_gen(infolder, folder_list)
    print("")
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", 50)
    top_fifty_dict = get_average_top_score_per_gen(infolder, folder_list, 50)
    print("")
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", 20)
    top_twenty_dict = get_average_top_score_per_gen(infolder, folder_list, 20)
    print("")
    print("Average for Top Scoring Compounds")
    print("Number of top scoring compounds: ", 10)
    top_ten_dict = get_average_top_score_per_gen(infolder, folder_list, 10)
    print("")
    print("Best Score per generation")
    print("Number of top scoring compounds: ", 1)
    top_one_dict = get_average_top_score_per_gen(infolder, folder_list, 1)
    print("")
    print("")
    dict_of_averages = {}
    dict_of_averages["average_affinity_dict"] = average_affinity_dict
    dict_of_averages["top_fifty_dict"] = top_fifty_dict
    dict_of_averages["top_twenty_dict"] = top_twenty_dict
    dict_of_averages["top_ten_dict"] = top_ten_dict
    dict_of_averages["top_one_dict"] = top_one_dict

    return dict_of_averages

def generate_figures(vars):
    """
    This runs everything to make a line plot of the results of an Autogrow
    simulation.

    Inputs:
    :param dict vars: dict of user variables which will govern how the
        programs runs
    """

    infolder = vars["infolder"]
    outfile = vars["outfile"]

    all_folders_list = [
        f for f in sorted(os.listdir(infolder)) if os.path.isdir(infolder + f)
    ]
    folder_list = []
    for folder in all_folders_list:
        if folder != "Data" and len(folder.split("_")) == 2:
            folder_list.append(folder)

    folder_list.sort(key=lambda x: int(x.split("_")[1]))

    dict_of_averages = print_data_table(infolder, folder_list)
    run_plotter(vars, dict_of_averages, outfile)

######## Handle Variables #####
def retrieve_vars_dict(autogrow_vars_json):
    """
    This will retrieve a variable dictionary from a AutoGrow vars json file.

    Inputs:
    :param str autogrow_vars_json: path to AutoGrow json variable file
    Returns:
    :returns: dict vars: a dictionary of variable to use
    """
    if os.path.exists(autogrow_vars_json) is False:
        raise Exception("variable file could not be found. It should be the \
            vars.json file written by AutoGrow in the output folder of the run.")
    try:
        with open(autogrow_vars_json, "r") as f:
            vars = json.load(f)
    except:
        raise Exception("variable file would not import. It should be the \
            vars.json file written by AutoGrow in the output folder of the run.")
    return vars
#
def process_inputs(inputs):
    """
    This will handle processing all parameters.

    inputs:
    :params dict inputs: dictionary of argparse parameters
    Returns:
    :returns: dict vars_dict: dictionary of argparse parameters
    """

    # handle input information
    inputs["infolder"] = os.path.abspath(inputs["infolder"]) + os.sep
    if os.path.exists(inputs["infolder"]) is False:
        raise Exception("Input folder {} does not\
            exist.".format(inputs["infolder"]))

    # get vars dict from last run
    inputs["vars_json"] = inputs["infolder"] + "vars.json"
    if os.path.exists(inputs["vars_json"]) is False:
        raise Exception("Input folder {} does not contain the vars.json file \
            necessary to run script. Please make sure the vars.json is in the \
            folder.".format(inputs["infolder"]))

    try:
        with open(inputs["vars_json"], "r") as f:
            vars_dict = json.load(f)
    except:
        raise Exception("variable file would not import. It should be the \
            vars.json file written by AutoGrow in the output folder of the run.")


    if "outfile_format" in inputs.keys():
        if inputs["outfile_format"] is None:
            inputs["outfile_format"] = "svg"
        if inputs["outfile_format"].lower() not in ["svg", "png", "jpg", "pdf"]:
            raise Exception("outfile_format not a valid format")

    if "outfile" in inputs.keys():
        if inputs["outfile"] is not None:
            if os.path.dirname(inputs["outfile"]) is False:
                try:
                    os.mkdir(os.path.dirname(inputs["outfile"]))
                except:
                    raise Exception("outfile directory does not exist")
            if os.path.dirname(inputs["outfile"]) is False:
                raise Exception("outfile directory does not exist")
        else:
            inputs["outfile"] = inputs["infolder"] + os.sep + \
                "data_line_plot." + inputs["outfile_format"]

    else:
        inputs["outfile"] = inputs["infolder"] + os.sep + \
            "data_line_plot." + inputs["outfile_format"]

    # update --plot_reference_lines
    if "plot_reference_lines" not in inputs.keys():
        inputs["plot_reference_lines"] = None

    if inputs["plot_reference_lines"] is not None:
        # names of all matplotlib color options
        matplot_colors = matplotlib.colors.get_named_colors_mapping().keys()

        ref_lines = inputs["plot_reference_lines"].replace("[[", "[").replace("]]", "]")
        ref_lines = ref_lines.split("],")
        ref_lines = [ref.replace("]", "").replace("[", "").split(",") for ref in ref_lines]

        new_ref_lines = []
        failed_io = False
        for ref_info in ref_lines:
            if len(ref_info) != 3:
                failed_io = True
                break

            # make new list with 1st item the str name
            temp_ref_lines = [str(ref_info[0])]

            try:
                temp_ref_lines.append(float(ref_info[1]))
            except:
                failed_io = True
                break
            if str(ref_info[2]) not in matplot_colors:
                print("COULD NOT FIND COLOR: " + str(ref_info[2]))
                failed_io = True
                break
            temp_ref_lines.append(str(ref_info[2]))
            new_ref_lines.append(temp_ref_lines)

        if failed_io is True:
            printout = "\n --plot_reference_lines must be list of lists where each "
            printout = printout + "sublist has three pieces of information in this "
            printout = printout + "order:\n\t [name, value, matplotlib_color]\n"
            printout = printout + "more details can be found using the -h option\n"
            print(printout)
            raise Exception(printout)
        inputs["plot_reference_lines"] = new_ref_lines
    # overwrite and return vars_dict with input commands
    for key in inputs.keys():
        vars_dict[key] = inputs[key]

    return vars_dict
#

######################################
######################################
######################################
PARSER = argparse.ArgumentParser()

# Get needed info
PARSER.add_argument(
    "--outfile",
    "-o",
    metavar="param.outfile",
    required=False,
    default=None,
    help="Path to folder to output files. It will be created if does not exist. \
    If not provide it will be placed in the infolder/data_line_plot.svg",
)
PARSER.add_argument(
    "--outfile_format",
    metavar="param.outfile_format",
    type=str, default="svg",
    choices=["svg", "png", "jpg", "pdf"],
    help="The type of file for figure to be exported as default is .svg file.",
)
PARSER.add_argument(
    "--infolder",
    "-i",
    metavar="param.infolder",
    required=True,
    help="Path to input folder containing the AutoGrow run. This should be the \
        top folder which contains the vars.json file.",
)
PARSER.add_argument(
    "--plot_reference_lines",
    default=None,
    help="This will be a list of lists, with each sublist being a different \
        dotted-line reference to plot. For each sublist the order of \
        information should be: [name, value, matplotlib_color] \
        For example a [['Olaparib score',-12.8,'y'],['Niraparib score',-10.7,'k']] \
        will add horizontal dotted lines at -12.8 (yellow) and -10.7 (black) \
        with Olaparib and Niraparib added to the legend. \
        Spaces must be within quotes and not be between variables. \
        matplotlib colors can be found with mcolors.get_named_colors_mapping().keys()"
)


ARGSDICT = vars(PARSER.parse_args())

# copying ARGSDICT so we can delete out of while iterating through the
# original ARGSDICT
INPUTS = copy.deepcopy(ARGSDICT)

for k, v in ARGSDICT.items():
    if v is None:
        del INPUTS[k]

USER_VARS = process_inputs(INPUTS)


generate_figures(USER_VARS)

print("FINISHED {}".format(USER_VARS["outfile"]))

print("finished")
