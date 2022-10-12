"""
This script creates figures for all ligands which parented a given ligand.

All compounds for the entire AutoGrow run will be compiled into a dictionary \
which is used to search when tracing lineages. We pickle these dictionaries so \
that if this script is run multiple times these dictionaries do not need to be \
recreated. For this reason the 1st time running this script on a data set will \
take longer than future runs.
"""
import os
import sys
import glob
import argparse
import json
import copy
import pickle

import matplotlib.pyplot as plt

import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import Draw, AllChem
from PIL import Image

#Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog('rdApp.*')
##################################################################
##################################################################
########### BASIC OPERATIONS #####################################
##################################################################
##################################################################
def get_obj_from_pickle_file(file_path):
    """
    This functions retrieves objects from a pickle_file
    Inputs:
    :param str file_path: path to pickle File
    Returns:
    :returns: unknown objects: object(s) from a pickle file
    """
    with open(file_path, 'rb') as handle:
        objects = pickle.load(handle)
    return objects

def write_pickle_to_file(file_path, obj):
    """
    This functions pickles an object into a pickle_file
    Inputs:
    :param str file_path: path to output pickle File
    :param unknown obj: object(s) to pickle
    """
    with open(file_path, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

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

#####################################################################
# Make images
#####################################################################

def get_image_dimensions(imagefile):
    """
    Helper function that returns the image dimensions.

    :param: imagefile str (path to image)
    :return dict (of the form: {width:<int>, height=<int>, size_bytes=<size_bytes>)
    """
    # Inline import for PIL because it is not a common library
    with Image.open(imagefile) as img:
        # Calculate the width and hight of an image
        width, height = img.size

    # calculate the size in bytes
    size_bytes = os.path.getsize(imagefile)

    return dict(width=width, height=height, size_bytes=size_bytes)

def get_grid_img(img_files_list, list_printout_info, result_grid_filename):
    """
    This will plot a row of imgs and save them to a file.

    Inputs:
    :param list img_files_list: list of paths to img file for each subplot
    :param list list_printout_info: list of info to add as
        caption to each subplot in order
    :param str result_grid_filename: path to outfile
    """

    images_count = len(img_files_list)
    dimmension_dict = get_image_dimensions(img_files_list[0])
    width = dimmension_dict["width"] / 30
    height = dimmension_dict["height"]/ 30
    # size_bytes = dimmension_dict["size_bytes"]

    fig, axs_list = plt.subplots(1, images_count, figsize=(width*images_count, height))

    if len(img_files_list) == 1:
        sub_ax = axs_list
        image_filename = img_files_list[0]
        printout = list_printout_info[0]
        plt_image = plt.imread(os.path.abspath(image_filename), printout)
        sub_ax.imshow(plt_image)

        sub_ax.set_title(printout, fontsize=40, fontweight="bold")
        sub_ax.grid(False)
        sub_ax.axis(False)
        sub_ax.autoscale_view('tight')

    else:
        for sub_ax, image_filename, printout in zip(axs_list, img_files_list, list_printout_info):
            plt_image = plt.imread(os.path.abspath(image_filename), printout)
            sub_ax.imshow(plt_image)

            sub_ax.set_title(printout, fontsize=40)
            sub_ax.grid(False)
            sub_ax.axis(False)
            sub_ax.autoscale_view('tight')

            del plt_image

    plt.savefig(result_grid_filename)
    del fig

def make_single_image_files(vars, lineage_dict, mol_dict):
    """
    This function will create individual image files for each ligand in
    an ancestry

    Inputs:
    :param dict vars: dictionary of variable to use
    :param dict lineage: a dict of lists of ancestors where the keys are the
        generation number relative to the creation of mol_name and the lists
        are each parent ligands full-length name, or None if there aren't two
        parents
    :param str mol_name: full-length name of child ligand to find parents.

    Returns:
    :returns: dict lineage: a dict of lists of ancestors where the keys are the
        generation number relative to the creation of mol_name and the lists
        are each parent ligands full-length name, or None if there aren't two
        parents
    """
    if len(list(lineage_dict.keys())) <= 6:
        img_size = 500
    else:
        img_size = 250



    # make single img files for each ligand
    # make a blank None image used later for spacers
    mol_none = Chem.MolFromSmiles("")
    img = Draw.MolsToGridImage([mol_none], molsPerRow=1, subImgSize=(img_size, img_size))
    img_file_name = vars["single_image_folder"]+ "None.png"
    img.save(img_file_name)

    for mol_name in mol_dict.keys():
        mol = copy.deepcopy(mol_dict[mol_name][-1])
        tmp = AllChem.Compute2DCoords(mol)

        img = Draw.MolsToGridImage([mol], molsPerRow=1,
                                   subImgSize=(img_size, img_size))
        img_file_name = vars["single_image_folder"]+ mol_name + ".png"
        img.save(img_file_name)
        del tmp

#
def make_image_files(vars, lineage_dict, mol_dict):
    """
    This function will create individual image files for each ligand in
    an ancestry and a full ancestry

    Inputs:
    :param dict vars: dictionary of variable to use
    :param dict lineage: a dict of lists of ancestors where the keys are the
        generation number relative to the creation of mol_name and the lists
        are each parent ligands full-length name, or None if there aren't two
        parents
    :param str mol_name: full-length name of child ligand to find parents.

    Returns:
    :returns: dict lineage: a dict of lists of ancestors where the keys are the
        generation number relative to the creation of mol_name and the lists
        are each parent ligands full-length name, or None if there aren't two
        parents
    """

    # make individual ligand images
    make_single_image_files(vars, lineage_dict, mol_dict)

    if os.path.exists(vars["ancestry_image_folder"]) is False:
        os.mkdir(vars["ancestry_image_folder"])

    for gen_num in lineage_dict.keys():
        result_grid_filename = vars["ancestry_image_folder"] + \
            str(gen_num) + ".png"
        lineage_name_list = lineage_dict[gen_num]
        img_files_list = []
        list_printout_info = []
        for mol_name in lineage_name_list:
            if mol_name is None:
                img_files_list.append(vars["single_image_folder"] + \
                "None.png")
                list_printout_info.append("")
            else:
                img_files_list.append(vars["single_image_folder"] + \
                    mol_name + ".png")


                # set properties
                if mol_dict[mol_name][4] is None:
                    printout = str(mol_dict[mol_name][2]) + "\nVina: " \
                        +  "COMP kcal/mol"
                else:
                    printout = str(mol_dict[mol_name][2]) + "\nVina: " \
                        + str(mol_dict[mol_name][4])+ " kcal/mol"

                list_printout_info.append(printout)


        get_grid_img(img_files_list, list_printout_info, result_grid_filename)


#####################################################################
# get parents for a ligand
#####################################################################

def get_parents_full_names(child_name, master_shortname_mol_dict):
    """
    Get full-length names for each parent for a given child ligand.
    Will return as list of names for parents. This will always be a list of 2.
    There are three options of what is returned:
        1) child ligand has no parents: ie) source/complementary ligand)
            will return [None, None]
        2) child ligand has 1 parent: ie) single reactant mutant
            will return ["parent_1_name", None]
        3) child ligand has 2 parent: ie) crossover or two reactant mutation
            will return ["parent_1_name", "parent_2_name"]

    Inputs:
    :param str child_name: full-length name of child ligand to find parents.
        ligand from the AutoGrow run. keys are full-length name of the ligands.
    :param dict master_shortname_mol_dict: dictionary where keys are
        shorthand names and the items the full-length name.
    Returns:
    :returns: list parent_list: a list of string or Nones for each parent.
        1) child ligand has no parents: ie) source/complementary ligand)
            will return [None, None]
        2) child ligand has 1 parent: ie) single reactant mutant
            will return ["parent_1_name", None]
        3) child ligand has 2 parent: ie) crossover or two reactant mutation
            will return ["parent_1_name", "parent_2_name"]
    """
    # Handle if no parents
    if "(" not in child_name and ")" not in child_name:
        return [None, None]

    parents_info = child_name.split(")")[0].replace("(", "")
    # Handle single parent cases
    if "+" not in parents_info:
        parent_1_short = parents_info
        if parent_1_short not in master_shortname_mol_dict.keys():
            raise Exception("a parent is not in master_shortname_mol_dict " \
                + "this means that the dictionary is missing information on" \
                + " a ligand. missing parrent is: {}".format(parent_1_short))

        parent_1_name = master_shortname_mol_dict[parent_1_short]
        return [parent_1_name, None]


    parent_1_short = parents_info.split("+")[0]
    parent_2_short = parents_info.split("+")[1]

    if parent_1_short not in master_shortname_mol_dict.keys():
        raise Exception("a parent is not in master_shortname_mol_dict " \
            + "this means that the dictionary is missing information on" \
            + " a ligand. missing parrent is: {}".format(parent_1_short))

    if parent_2_short not in master_shortname_mol_dict.keys():
        raise Exception("a parent is not in master_shortname_mol_dict " \
            + "this means that the dictionary is missing information on" \
            + " a ligand. missing parrent is: {}".format(parent_2_short))

    parent_1_name = master_shortname_mol_dict[parent_1_short]
    parent_2_name = master_shortname_mol_dict[parent_2_short]

    return [parent_1_name, parent_2_name]
#
def get_all_ancestors(mol_name, master_shortname_mol_dict):
    """
    This function will obtain all ancestors and store them in a
    diction where the key is the generation number (relative to the creation of the
    mol the user requested) and the items are lists of full-length molecule names.

    These lists must be ordered and will contain None as a place holder.
    Each previous generation will be double the length of its successors, even if
    multiple entries are None. This way we can create a lineage tree.

    Gen 0: [ A,    None,      B,     C,      D,    None,     None,     None]
             |      |         |      |       |      |          |      |
              ______           ______         ______            ______
                |                 |             |                  |
    Gen 1: [ (A)E,              (B+C)F        (D)G,             None       ]
              |                    |            |                    |
               ____________________              ____________________
                        |                                   |
    Gen 2: [         (E+F)H,                              (G)I             ]
                        |                                   |
                        _____________________________________
                                        |
    Gen 3: [                        (H+I)J                                 ]

    This tree would return:
    {
        0:[A, None, B, C, D, None, None, None],
        1:[(A)E, (B+C)F (D)G, None],
        2:[(E+F)H, (G)I],
        3:[(H+I)J],
    }


    Inputs:
    :param str mol_name: full-length name of child ligand to find parents.
    :param dict master_mol_dict: dictionary containing the information from every
        ligand from the AutoGrow run. keys are full-length name of the ligands.

    Returns:
    :returns: dict lineage: a dict of lists of ancestors where the keys are the
        generation number relative to the creation of mol_name and the lists
        are each parent ligands full-length name, or None if there aren't two
        parents
    """
    if ")Gen_" not in mol_name:
        raise Exception("mol_name provided either does not have parents " \
            +"and is likely from source compound list.")
    start_generation = int(mol_name.split(")Gen_")[-1].split("_")[0])
    lineage_dictionary = {}

    lineage_dictionary[start_generation] = [mol_name]
    # check that parents exist for main mol
    parents_to_check = get_parents_full_names(mol_name, master_shortname_mol_dict)
    if parents_to_check == [None, None]:
        raise Exception("mol_name provided either does not have parents " \
            +"and is likely from source compound list.")

    lineage_dictionary[start_generation - 1] = parents_to_check

    current_gen = start_generation - 2
    while current_gen >= 0:
        grand_parent_list = []


        for parent in parents_to_check:
            if parent is None:
                grand_parent_list.extend([None, None])
            else:
                parent_list = get_parents_full_names(parent,
                                                     master_shortname_mol_dict)

                grand_parent_list.extend(parent_list)

        parents_to_check = grand_parent_list
        lineage_dictionary[current_gen] = grand_parent_list

        current_gen = current_gen -1
        if list(set(parents_to_check)) == [None]:
            # All ancestors are None
            break

    return lineage_dictionary
#
#####################################################################
# make/retrive pickle dictionaries
#####################################################################

def make_master_shortname_mol_dict(vars, master_mol_dict):
    """
    Create and save a dictionary that can look up a full-length name using
    a shorthand name; keys are shorthand name w full-length names as items

    Save to a pickle file named master_shortname_mol_dict_pickle.

    Inputs:
    :param dict vars: dictionary of variable to use
    :param dict master_mol_dict: master dictionary with
        the long names of ligands as keys
    Returns:
    :returns: dict master_shortname_mol_dict: dictionary where keys are
        shorthand names and the items the full-length name.
    """
    master_shortname_mol_dict = {}
    for mol_entry in master_mol_dict.keys():
        short_name = master_mol_dict[mol_entry][2]
        master_shortname_mol_dict[short_name] = mol_entry

    # Write to pickle file
    master_shortname_mol_dict_pickle = vars["master_shortname_mol_dict_pickle"]
    write_pickle_to_file(master_shortname_mol_dict_pickle, master_shortname_mol_dict)

    return master_shortname_mol_dict
#
def merge_comp_and_ranked_dicts(vars):
    """
    Merge two dictions into one.

    Save to a pickle file named master_mol_dict_pickle.

    Inputs:
    :param dict vars: dictionary of variable to use
    Returns:
    :returns: dict master_mol_dict: master dictionary with all ligands entered
    """
    ranked_mol_dict_pickle = vars["ranked_mol_dict_pickle"]
    comp_dict_pickle = vars["comp_dict_pickle"]

    comp_mol_dict = get_obj_from_pickle_file(comp_dict_pickle)

    master_mol_dict = copy.deepcopy(comp_mol_dict)
    del comp_mol_dict

    ranked_mol_dict = get_obj_from_pickle_file(ranked_mol_dict_pickle)
    # Since there shouldn't be docking information in the comp_mol_dict
    # we can feel free to overwrite any duplicate entries. Any duplicate
    # entries would be caused by a ligand being in both the source ligands
    # and the complementary molecule library.
    for mol_entry in ranked_mol_dict.keys():
        master_mol_dict[mol_entry] = ranked_mol_dict[mol_entry]

    del ranked_mol_dict

    # Write to pickle file
    master_mol_dict_pickle = vars["master_mol_dict_pickle"]
    write_pickle_to_file(master_mol_dict_pickle, master_mol_dict)

    return master_mol_dict
#
def make_comp_mol_dict(vars):
    """
    Create and pickled dictionary of all complementary molecules.

    These mol_dicts can be quite large and memory intensive so we will create
    and save as pickle. And reopen later to minimize memory overhead.

    Inputs:
    :param dict vars: dictionary of variable to use
    """
    # Add complementary mols from reactions
    # Only valid for autoclickchem rxns

    comp_smi_list = glob.glob(vars["complementary_mol_directory"] + os.sep + "*.smi")
    if len(comp_smi_list) == 0:
        raise Exception("No .smi files found for complementary_mol_directory.\n" + \
            "please check: {}".format(vars["complementary_mol_directory"]))
    comp_dict = {}
    for smi in comp_smi_list:
        comp_mol_list = get_usable_format(smi)
        for mol_entry in comp_mol_list:
            comp_dict[mol_entry[1]] = mol_entry[0]
        del comp_mol_list
    del comp_smi_list

    #Make this match those with scores
    for mol_name in comp_dict.keys():
        mol = Chem.MolFromSmiles(comp_dict[mol_name])
        temp_info = [comp_dict[mol_name], mol_name, mol_name, mol_name, None, None, mol]

        comp_dict[mol_name] = temp_info

    comp_dict_pickle = vars["comp_dict_pickle"]
    write_pickle_to_file(comp_dict_pickle, comp_dict)
    del comp_dict
#
def make_ranked_files_mol_dict(vars):
    """
    Create and pickle dictionary of all ranked ligands.

    These mol_dicts can be quite large and memory intensive so we will create
    and save as pickle. And reopen later to minimize memory overhead.

    Inputs:
    :param dict vars: dictionary of variable to use
    """
    dir_w_all_gens = vars["input_dir"]
    ranked_file_list = glob.glob(str(os.sep).join([dir_w_all_gens,
                                                   "generation_*",
                                                   "generation_*_ranked.smi"]))
    if len(ranked_file_list) == 0:
        raise Exception("Could not find any ranked.smi files within the input_dir.\
            Please make sure input_dir has folders named 'generation_' + int that \
            contain .smi files named 'generation_{}_ranked.smi'.format(int) ")

    source_compound_file = vars["source_compound_file"]
    ranked_file_list = [x for x in list(set(ranked_file_list)) if x is not source_compound_file]

    # make a list of molecule information
    mol_list = []
    for i in ranked_file_list:
        mol_list.extend(get_usable_format(i))

    # Add Source compounds
    source_compound_list = get_usable_format(source_compound_file)
    len_of_each_mol_info = [len(x) for x in source_compound_list]
    if len(list(set(len_of_each_mol_info))) != 1:
        print(list(set(len_of_each_mol_info)))
        raise Exception("The source compound file is inconsistently with the number\
            of columns per line. Please correct this so that each line has the \
            same number of columns.")

    # If the source compounds were previously docked keep all info because there is docking info
    if vars["use_docked_source_compounds"] is True:
        mol_list.extend(mol_list)
    else:
        # If there are only two columns we assume it is the SMILES and the source name
        # Otherwise we print a message saying we are ignoring any additional information
        # because the --use_docked_source_compounds==False
        if list(set(len_of_each_mol_info))[0] != 2:
            print("\nWARNING: There are multiple columns within the source \
                compound file ({}), but --use_docked_source_compounds is set \
                to False. You will also need to delete the pickled dictionaries \
                produced by this script before re-running the script. \n\
                We will ignore any information other than the first \
                two columns in the source compound file. This may mean that we \
                ignore docking scores or use the full-length names of compounds \
                in generation zero.\n".format(source_compound_file))
        new_source_compound_list = []
        for mol_info in source_compound_list:
            temp_info = [mol_info[0], mol_info[1], mol_info[1], mol_info[1], None, None]
            new_source_compound_list.append(temp_info)
        mol_list.extend(new_source_compound_list)

    new_list = []
    for x in mol_list:
        temp = []
        for y in x:
            try:
                y = float(y)
                temp.append(y)
            except:
                temp.append(y)
        temp.append(Chem.MolFromSmiles(temp[0]))
        new_list.append(temp)

    del mol_list
    mol_dict = {}
    for mol_entry in new_list:
        if mol_entry[1] in mol_dict.keys():

            # source compounds may be unranked or may
            # be ranked from advancement but not in the source compound .smi
            # so we check to try to keep the ranked version if it exists.
            if type(mol_entry[-3]) not in [float, int]:
                # the current line is not ranked
                continue
            #
            if type(mol_dict[mol_entry[1]][-3]) not in [float, int]:
                # the current line is ranked but the previously entered version
                # was not. Lets overwrite it
                mol_dict[mol_entry[1]] = mol_entry

            if float(mol_entry[-3]) < float(mol_dict[mol_entry[1]][-3]):
                # Both entries of the ligand are ranked. Let us take
                # the better ranked version of the two
                mol_dict[mol_entry[1]] = mol_entry
            else:
                # the previous version was ranked better. Lets keep that.
                continue
        else:
            # new entry lets add it to the dictionary.
            mol_dict[mol_entry[1]] = mol_entry

    del new_list
    # Write to pickle file
    ranked_mol_dict_pickle = vars["ranked_mol_dict_pickle"]
    write_pickle_to_file(ranked_mol_dict_pickle, mol_dict)
    del mol_dict
#
def get_mol_dict(vars):
    """
    Retrieve pickled dictionary of all ligand information if it already has
    been created, or it creates the dictionary of ligand information and saves
    it as a pickle file.

    Inputs:
    :param dict vars: dictionary of variable to use
    Returns:
    :returns: dict master_mol_dict: dictionary containing the information from every
        ligand from the AutoGrow run. keys are full-length name of the ligands.
    :returns: dict master_shortname_mol_dict: dictionary where keys are
        shorthand names and the items the full-length name.
    """
    ranked_mol_dict_pickle = vars["ranked_mol_dict_pickle"]
    comp_dict_pickle = vars["comp_dict_pickle"]
    master_mol_dict_pickle = vars["master_mol_dict_pickle"]
    master_shortname_mol_dict_pickle = vars["master_shortname_mol_dict_pickle"]

    if os.path.exists(master_mol_dict_pickle) is False:
        # Handle ranked files
        if os.path.exists(ranked_mol_dict_pickle) is False:
            # ranked_mol_dict_pickle does not exist. Need to make and pickle
            # Will reopen later to minimize memory overhead
            print("Creating ranked_files_mol_dict")
            make_ranked_files_mol_dict(vars)

        # Handle complementary molecules
        if os.path.exists(comp_dict_pickle) is False:
            # comp_dict_pickle does not exist. Need to make and pickle
            # Will reopen later to minimize memory overhead
            print("Creating comp_mol_dict")
            make_comp_mol_dict(vars)

        # merge molecule dictionaries together to create a master dictionary
        print("Creating master_mol_dict")
        master_mol_dict = merge_comp_and_ranked_dicts(vars)

    else:

        # The master dictionary already exists so we don't need to make the
        # smaller dictionaries.
        # Just get it from the pickle file
        print("Getting master_mol_dict from pickle file")
        master_mol_dict = get_obj_from_pickle_file(master_mol_dict_pickle)

    if os.path.exists(master_shortname_mol_dict_pickle) is False:
        # Need to create the master_shortname_mol_dict
        master_shortname_mol_dict = make_master_shortname_mol_dict(vars, master_mol_dict)
    else:
        # The master_shortname_mol_dict
        # Get it from the pickle file
        print("Getting master_shortname_mol_dict from pickle file")
        master_shortname_mol_dict = get_obj_from_pickle_file(master_shortname_mol_dict_pickle)

    return master_mol_dict, master_shortname_mol_dict

##################################################################

#####################################################################
# I/O
#####################################################################
def get_full_length_mol_name(vars, master_mol_dict, master_shortname_mol_dict):
    """
    Get full-length mol_name and make sure that it is in the master_mol_dict

    Inputs:
    :param dict vars: dictionary of variable to use

    :param dict master_mol_dict: dictionary containing the information from every
        ligand from the AutoGrow run. keys are full-length name of the ligands.
    :param dict master_shortname_mol_dict: dictionary where keys are
        shorthand names and the items the full-length name.
    Returns:
    :returns: str mol_name: full-length name of ligand.
    """
    if vars["mol_name"] not in master_mol_dict.keys():
        if vars["mol_name"] not in master_shortname_mol_dict.keys():

            # may be a gypsum variant with '__{}'.format(num) at the end
            # ie Gen_5_Mutant_46_684401 could be represented as Gen_5_Mutant_46_684401__1
            if "__" in vars["mol_name"]:
                test_name = vars["mol_name"].split("__")[0]
                if (test_name not in master_shortname_mol_dict.keys() and \
                        test_name not in master_mol_dict.keys()):
                    printout = "mol_name provided not found in shorthand or" \
                        + "full-length dictionaries. Please check that mol_name is in" \
                        + "the AutoGrow run tested. \n" \
                        + "Name provided is :\n\t{}".format(vars["mol_name"] \
                        + "\nName should look like is :" \
                        + "\n\t  (Gen_2_Mutant_7_97143)Gen_4_Mutant_7_802531" \
                        + "\n\t\t or \n\t Gen_4_Mutant_7_802531")
                    print(printout)

                    raise Exception(printout)

                if test_name in master_shortname_mol_dict.keys():
                    mol_name = master_shortname_mol_dict[test_name]
                else:
                    mol_name = test_name
                del test_name
            else:
                printout = "mol_name provided not found in shorthand or" \
                    + "full-length dictionaries. Please check that mol_name is in" \
                    + "the AutoGrow run tested. \n" \
                    + "Name provided is :\n\t{}".format(vars["mol_name"] \
                    + "\nName should look like is :" \
                    + "\n\t  (Gen_2_Mutant_7_97143)Gen_4_Mutant_7_802531" \
                    + "\n\t\t or \n\t Gen_4_Mutant_7_802531")
                print(printout)

                raise Exception(printout)

        mol_name = master_shortname_mol_dict[vars["mol_name"]]
    else:
        # original name is already full-length name
        mol_name = vars["mol_name"]

    return mol_name
#
def run_purge_previous_pickled_files(vars):
    """
    This will delete previously created pickled files within the input_dir.
    The four files it will delete are:
    `$input_dir/comp_dict_pickle`, `$input_dir/master_mol_dict_pickle`,
    `$input_dir/master_shortname_mol_dict_pickle`,
    and `$input_dir/ranked_mol_dict_pickle`.

    These files save time when you are tracing the lineage of multiple
    compounds, however purging these files may be helpful for space saving
    or if it had been previously run with an invalid input variable.=

    Following file deletion the program will terminate.

    inputs:
    :params vars inputs: dictionary of argparse parameters
    """
    print("\nDELETING PREVIOUSLY GENERATED PICKLED FILES.\n")
    input_dir = vars["input_dir"] + os.sep
    if os.path.exists(input_dir) is False:
        raise Exception("Input folder {} does not\
            exist.".format(input_dir))
    for file_name in ["comp_dict_pickle", "master_mol_dict_pickle",
                      "master_shortname_mol_dict_pickle", "ranked_mol_dict_pickle"]:
        file_path = input_dir + file_name
        if os.path.exists(file_path) is False:
            printout = "Could not delete {} file".format(file_name)
            printout = printout + " as it was not located at:\n\t {}\n".format(file_path)
            print(printout)
        else:
            try:
                os.remove(file_path)
                print("Deleted: {}".format(file_path))
            except:
                printout = "WARNING: Could not delete {} file.\n".format(file_name)
                printout = printout + "\tPlease check file permissions of:"
                printout = printout + "\n\t\t {}\n".format(file_path)
                print(printout)
            # Check that is was successfully deleted
            if os.path.exists(file_path) is False:
                print("Deleted: {}".format(file_path))

    print("Attempt to delete files completed.")
    sys.exit(0)

def process_inputs(inputs):
    """
    This will handle processing all parameters.

    inputs:
    :params dict inputs: dictionary of argparse parameters
    Returns:
    :returns: dict inputs: dictionary of argparse parameters
    """

    # handle input information
    inputs["input_dir"] = os.path.abspath(inputs["input_dir"]) + os.sep
    if os.path.exists(inputs["input_dir"]) is False:
        raise Exception("Input folder {} does not\
            exist.".format(inputs["input_dir"]))

    # get vars dict from last run
    inputs["vars_json"] = inputs["input_dir"] + "vars.json"
    if os.path.exists(inputs["vars_json"]) is False:
        raise Exception("Input folder {} does not contain the vars.json file \
            necessary to run script. Please make sure the vars.json is in the \
            folder.".format(inputs["input_dir"]))

    try:
        with open(inputs["vars_json"], "r") as f:
            vars_dict = json.load(f)
    except:
        raise Exception("variable file would not import. It should be the \
            vars.json file written by AutoGrow in the output folder of the run.")
    if inputs["complementary_mol_directory"] in ["", None]:

        # Get complementary_mol_directory from vars.json
        if vars_dict["complementary_mol_directory"] not in ["", None]:
            if os.path.exists(vars_dict["complementary_mol_directory"]):
                inputs["complementary_mol_directory"] = vars_dict["complementary_mol_directory"]
            else:
                # Can not find the one used in vars.json list
                raise Exception("Please provide path to complementary_mol_directory. \
                    vars.json file lists custom path to complementary_mol_director={} \
                    but this directory can not be found. Please provide path using \
                    --complementary_mol_directory $PATH/TO/complementary_mol_directory/")

        # Get complementary_mol_directory from vars.json
        elif vars_dict["rxn_library"].lower() in ["click_chem_rxns", "robust_rxns", "all_rxns"]:
            dir_above_script_dir = str(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

            complementary_mol_directory = str(os.sep).join([
                dir_above_script_dir, "autogrow", "operators", "mutation",
                "smiles_click_chem", "reaction_libraries",
                vars_dict["rxn_library"].lower(), "complementary_mol_dir"])

            complementary_mol_directory = os.path.abspath(complementary_mol_directory)
            if os.path.exists(complementary_mol_directory) is False:
                raise Exception("Please provide path to complementary_mol_directory. \
                    Could not find the location of the directory")

            inputs["complementary_mol_directory"] = complementary_mol_directory + os.sep
        else:
            raise Exception("Please provide path to complementary_mol_directory. \
                Could not find the location of the directory")
    else: # complementary_mol_directory was provided.
        inputs["complementary_mol_directory"] = \
            os.path.abspath(inputs["complementary_mol_directory"]) + os.sep

        if os.path.exists(inputs["complementary_mol_directory"]) is False:
            # Can not find the one used in vars.json list
            raise Exception("Please provide path to complementary_mol_directory. \
                provided path could not be found " \
                    + ":\n\t{}".format(inputs["complementary_mol_directory"]))


        if len(glob.glob(inputs["complementary_mol_directory"] + "*.smi")) == 0:
            sub_dir = inputs["complementary_mol_directory"] + os.sep \
                + "complementary_mol_dir" + os.sep
            if len(glob.glob(sub_dir + "*.smi")) == 0:
                raise Exception("Please provide path to complementary_mol_directory. " \
                    + "provided path had no .smi files: " \
                    + "\n\t{}".format(inputs["complementary_mol_directory"]))
            # They provided 1 directory up...
            inputs["complementary_mol_directory"] = sub_dir

    if "use_docked_source_compounds" not in inputs.keys() or \
                inputs["use_docked_source_compounds"] in ["", None]:
        # Get whether they used use_docked_source_compounds from vars.json
        if "use_docked_source_compounds" in vars_dict:
            if vars_dict["use_docked_source_compounds"] in [True, False]:
                inputs["use_docked_source_compounds"] = vars_dict["use_docked_source_compounds"]
            else:
                raise Exception("Please provide the --use_docked_source_compounds setting \
                    used during the run. We could not auto-detect from the vars file.")
        else:
            raise Exception("Please provide the --use_docked_source_compounds setting \
                used during the run. We could not auto-detect from the vars file.")
    else:
        if inputs["use_docked_source_compounds"] in [True, "true", "True"]:
            inputs["use_docked_source_compounds"] = True
        elif inputs["use_docked_source_compounds"] in [False, "false", "False"]:
            inputs["use_docked_source_compounds"] = False
        else:
            raise Exception("Please check the --use_docked_source_compounds setting provided." \
                " --use_docked_source_compounds should be True or False.")

    # Handle output directory
    inputs["output_dir"] = os.path.abspath(inputs["output_dir"]) + os.sep
    if os.path.exists(inputs["output_dir"]) is False:

        try:
            os.mkdir(inputs["output_dir"])
            print("Made the output dir at: {}".format((inputs["output_dir"])))
        except:
            pass
        if os.path.exists(inputs["output_dir"]) is False:
            raise Exception("Output folder {} does not\
                exist.".format(inputs["output_dir"]))

    # handle source_compound_file .smi file
    if type(inputs["source_compound_file"]) is not str or inputs["source_compound_file"] == "":
        raise Exception("--source_compound_file must be provided. It should be \
            the tab-delineated .smi file used to seed generation zero of the \
            AutoGrow run. This is a mandatory file.")

    inputs["source_compound_file"] = os.path.abspath(inputs["source_compound_file"])
    if os.path.exists(inputs["source_compound_file"]) is False:
        raise Exception("source_compound_file could not be found \
            at: {}".format(inputs["source_compound_file"]))
    if inputs["source_compound_file"].split(".")[-1] != "smi":
        raise Exception("--source_compound_file must be provided. It should be \
        the tab-delineated .smi file used to seed generation zero of the \
        AutoGrow run. This is a mandatory file.")

    # assign the destination for our pickle files (may already exist)
    inputs["ranked_mol_dict_pickle"] = inputs["input_dir"] + "ranked_mol_dict_pickle"
    inputs["comp_dict_pickle"] = inputs["input_dir"] + "comp_dict_pickle"
    inputs["master_mol_dict_pickle"] = inputs["input_dir"] + "master_mol_dict_pickle"
    inputs["master_shortname_mol_dict_pickle"] = inputs["input_dir"] \
        + "master_shortname_mol_dict_pickle"

    # handle singles image folder
    inputs["single_image_folder"] = inputs["output_dir"] + "single_image_folder" + os.sep
    if os.path.exists(inputs["single_image_folder"]) is False:
        os.mkdir(inputs["single_image_folder"])
    if "mol_name" not in inputs.keys():
        inputs["mol_name"] = None
    if inputs["pre_run"] is False:
        inputs["ancestry_image_folder"] = inputs["output_dir"]  \
            + "ancestry_"+ inputs["mol_name"]  + os.sep
    # Will wait to create this folder until its needed

    # Handle the cleanup variable purge_previous_pickled_files
    if "purge_previous_pickled_files" in inputs.keys():
        if inputs["purge_previous_pickled_files"] in [True, "true", "True"]:
            # We will delete files
            inputs["purge_previous_pickled_files"] = True
        elif inputs["purge_previous_pickled_files"] in [False, "false", "False"]:
            # We will not delete files
            inputs["purge_previous_pickled_files"] = False
        else:
            # Can not understand the input option
            raise Exception("Please check the --purge_previous_pickled_files setting provided." \
                " --purge_previous_pickled_files should be True or False.")
    else:
        inputs["purge_previous_pickled_files"] = False
    # If true delete files and terminate program
    if inputs["purge_previous_pickled_files"] is True:
        run_purge_previous_pickled_files(inputs)

    return inputs
#
def run_everything(vars):
    """
    This script runs everything
    Inputs:
    :params dict INPUTS: dictionary of argparse parameters
    """
    master_mol_dict, master_shortname_mol_dict = get_mol_dict(vars)

    if vars["pre_run"] is True or vars["mol_name"] in [None, "None", ""]:
        print("pre-run completed")
        sys.exit(0)
    mol_name = get_full_length_mol_name(vars, master_mol_dict, master_shortname_mol_dict)

    print("The full-length name of the ligand is: ", mol_name)
    print("")

    lineage_dict = get_all_ancestors(mol_name, master_shortname_mol_dict)

    # make a simplified master_mol_dict
    mol_dict = {}
    for gen_list in lineage_dict.keys():
        for lig_name in lineage_dict[gen_list]:
            if lig_name is not None:
                mol_dict[lig_name] = master_mol_dict[lig_name]
    del master_mol_dict
    del master_shortname_mol_dict

    # Write all information of lieage to .smi file
    lineage_smi = vars["output_dir"] + \
            str(mol_name) + "_lineage.smi"
    lineage_list = []
    for gen_num in lineage_dict.keys():
        lineage_list.extend(lineage_dict[gen_num])
    printout = ""
    for lig_name in lineage_list:
        if lig_name is None:
            continue
        temp = copy.deepcopy(mol_dict[lig_name])
        del temp[-1] # remove last item which is rdkit mol
        temp = "\t".join([str(x) for x in temp]) + "\n"
        printout = printout + temp

    with open(lineage_smi, 'w') as f:
        f.write(printout)

    # generate images
    make_image_files(vars, lineage_dict, mol_dict)

######################################
######################################
######################################
PARSER = argparse.ArgumentParser()

# Get needed info
PARSER.add_argument(
    "--output_dir",
    "-o",
    metavar="param.output_dir",
    required=True,
    help="Path to folder to output files. will be created if does not exist",
)
PARSER.add_argument(
    "--input_dir",
    "-i",
    metavar="param.input_dir",
    required=True,
    help="Path to input folder containing the AutoGrow run. This should be the \
        top folder which contains the vars.json file.",
)
PARSER.add_argument(
    "--complementary_mol_directory",
    metavar="param.complementary_mol_directory",
    required=False,
    default="",
    help="If using a custom complementary molecule library for mutations this \
    path is required. If not the script will try to autodetect the location of \
    the predefined complementary_mol_directory. Many molecules generated by \
    mutation will required the complementary molecule that helped spawn them.",
)
PARSER.add_argument(
    "--source_compound_file",
    metavar="param.source_compound_file",
    required=True,
    default="",
    help="This is the source .smi file used to seed generation zero of the \
    AutoGrow run. This is an essential file.",
)
PARSER.add_argument(
    "--pre_run",
    metavar="param.pre_run",
    default=False,
    help="If True this will compile the necessary dictions/picklefiles and then \
    terminate. These pickle files are stored in the input folder containing the \
    vars.json file from the AutoGrow run.",
)
PARSER.add_argument(
    "--mol_name",
    metavar="param.mol_name",
    default=None,
    help="This is the name of the molecule whose lineage will be traced back. \
    If not provided or None, the script will simply compile the necessary \
    dictions/picklefiles and then terminate. These pickle files are stored \
    in the input folder containing the vars.json file from the AutoGrow run.\
    example mol_name: Gen_5_Cross_203131 or Gen_4_Mutant_7_802531 \
        can also be provided as full-name ie: \
            (Gen_2_Mutant_7_97143)Gen_4_Mutant_7_802531",
)
PARSER.add_argument(
    "--use_docked_source_compounds",
    metavar="param.use_docked_source_compounds",
    choices=[True, False, "True", "False", "true", "false"],
    default=None,
    help="If True source ligands were docked prior to seeding generation 1. \
    If True and the source_compound file may already have the docking/fitness \
    metric score in -2 column of .smi file.\
    If False, generation 1 was randomly seeded by the source compounds with \
    no preference and there was no generation 0 testing. \
    If not provided this script will autodetect it from the vars.json \
    file if possible.",
)
PARSER.add_argument(
    "--purge_previous_pickled_files",
    metavar="param.purge_previous_pickled_files",
    choices=[True, False, "True", "False", "true", "false"],
    default=False,
    help="If True the script will delete the four pickled files previously \
    created by this script: `comp_dict_pickle`, `master_mol_dict_pickle`, \
    `master_shortname_mol_dict_pickle`, and `ranked_mol_dict_pickle`. \
    These files save time when you are tracing the lineage of multiple \
    compounds, however purging these files may be helpful for space saving \
    or if it had been previously run with an invalid input variable. \
    This does not affect the lineage files located in `output_dir`. \
    Program will terminate once these files are deleted.",
)


ARGSDICT = vars(PARSER.parse_args())

# copying ARGSDICT so we can delete out of while iterating through the
# original ARGSDICT
INPUTS = copy.deepcopy(ARGSDICT)

for k, v in ARGSDICT.items():
    if v is None:
        del INPUTS[k]

VARS = process_inputs(INPUTS)

run_everything(VARS)
