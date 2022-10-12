"""
This script will test a complementary molecule library to ensure all compounds
react in all reactions they may be used in.

Example submit:

python autogrow4/accessory_scripts/test_complementary_mol_library.py \
--rxn_library_file \
autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/ClickChem_rxn_library.json \
--function_group_library \
autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/ClickChem_functional_groups.json \
--complementary_mol_directory \
autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/complementary_mol_dir \
--output_folder autogrow4/accessory_scripts/output/
"""
import __future__

import os
import json
import copy
import argparse

import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import AllChem

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")


import support_scripts.Multiprocess as mp
import support_scripts.mol_object_handling as MOH

class SmilesClickChem():
    """
    This class will take a molecule and Mutate it by reacting it.

    This is modified from the AutoGrow source code file:
        /autogrow4/autogrow/operators/mutation/smiles_click_chem/smiles_click_chem.py
    Unused sections were removed for simplicity.
    """

    def __init__(self, rxn_library_variables, list_of_already_made_smiles):
        """
        init for SmilesClickChem. This will set up all the reaction and
        functional dictionaries required to Mutate a molecular

        Inputs:
        :param list rxn_library_variables: a list of user variables which
            define the rxn_library, rxn_library_file,
            complementary_mol_directory, and function_group_library. ie.
            rxn_library_variables = [vars['rxn_library'],
            vars['rxn_library_file'],
            vars['function_group_library'],vars['complementary_mol_directory']]
        :param list list_of_already_made_smiles: a list of lists. Each
            sublist contains info about a smiles made in this generation via
            mutation ie.[['O=C([O-])',
            '(Gen_3_Mutant_37_747+ZINC51)Gen_4_Mutant_15_52']]
        """

        # Unpackage the rxn_library_variables

        rxn_library = rxn_library_variables[0]
        rxn_library_file = rxn_library_variables[1]
        function_group_library = rxn_library_variables[2]
        complementary_mol_dir = rxn_library_variables[3]
        self.reaction_dict = self.retrieve_reaction_dict(
            rxn_library, rxn_library_file
        )
        # Retrieve the dictionary containing
        # all the possible ClickChem Reactions
        self.list_of_reaction_names = list(self.reaction_dict.keys())

        self.functional_group_dict = self.retrieve_functional_group_dict(
            rxn_library, function_group_library
        )
        self.complementary_mol_dict = self.retrieve_complementary_dictionary(
            rxn_library, complementary_mol_dir
        )

        # List of already predicted smiles
        self.list_of_already_made_smiles = [x[0] for x in list_of_already_made_smiles]

    def rxn_lib_format_json_dict_of_dict(self, old_dict):
        """
        json dictionaries  import as type unicode. This script converts all
        the keys and items to strings, with a few specific exceptions. It
        takes both the functional group dictionary and the reaction library.

        The reaction library is a dictionary of dictionary and has a few
        exceptions which are not intended to be strings. ie. the num_reactants
        which converts to interger and functional_groups which convert to a
        list of strings.

        The functional_group_dictionary is simply a dictionary with all items
        and keys needing to be strings.

        Inputs:
        :param dic old_dict: a dictionary of the the reaction library or
            functional groups. This is what is importanted from the .json file.

        Returns:
        :returns: dic new_dict: a dictionary of the the reaction library or
            functional groups where the unicode type items have been replaced with
            the proper python data types.
        """
        new_dict = {}
        for rxn_key in old_dict.keys():
            rxn_dic_old = old_dict[rxn_key]
            key_str = str(rxn_key)

            # For reaction libraries
            if type(rxn_dic_old) == dict:
                new_sub_dict = {}
                for key in rxn_dic_old.keys():
                    sub_key_str = str(key)
                    item = rxn_dic_old[key]

                    if sub_key_str == "num_reactants":
                        item = int(item)
                    elif sub_key_str == "functional_groups":
                        new_list = []
                        for i in item:
                            i_str = str(i)
                            new_list.append(i_str)

                        item = new_list
                    else:
                        item = str(item)

                    new_sub_dict[sub_key_str] = item
                new_dict[key_str] = new_sub_dict

            # For functional groups
            else:
                item = old_dict[rxn_key]
                new_dict[key_str] = str(item)

        return new_dict

    def retrieve_reaction_dict(self, rxn_library, rxn_library_file):
        """
        This is where all the chemical reactions for SmartClickChem are
        retrieved. If you want to add more just add a Custom set of reactions
        please add a folder to
        PATH/autogrow/operators/mutation/smiles_click_chem/Reaction_libraries/.
        They should be formatted as a dictionary of dictionary using the same
        format as :
        os.path.join(pwd,"reaction_libraries",
                    "click_chem_rxns","ClickChem_rxn_library.json")

        The reactions are written as SMARTS-reaction strings.

        This dictionary uses the reaction name as the key and the Reaction
        Smarts as the value.

        Inputs:
        :param str rxn_library: A string defining the choice of the reaction
            library. ClickChem uses the set of reactions from Autogrow 3.1.2.
            Custom means you've defined a path to a Custom library in
            vars['rxn_library_file']
        :param str rxn_library_file: a PATH to a Custom reaction library file
            formatted in a dictionary of dictionaries. in a .json file. This will
            be a blank string if one choses a predefined rxn_library option.

        Returns:
        :returns: dict reaction_dict: A dictionary containing all the
            reactions for ClickChemistry and all the information required to run
            the reaction
        """
        # Get the JSON file to import the proper reaction library
        pwd = os.path.dirname(__file__)
        if rxn_library_file == "":

            if rxn_library == "click_chem_rxns":
                rxn_library_file = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "click_chem_rxns",
                    "ClickChem_rxn_library.json"
                )
            elif rxn_library == "robust_rxns":
                rxn_library_file = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "robust_rxns",
                    "Robust_Rxns_rxn_library.json"
                )
            elif rxn_library == "all_rxns":
                rxn_library_file = os.path.join(
                    pwd,
                    "reaction_libraries",
                    "all_rxns",
                    "All_Rxns_rxn_library.json"
                    )
            elif rxn_library == "Custom":
                if os.path.exists(rxn_library_file) is False:
                    raise Exception(
                        "Custom rxn_library_file cannot be found. "
                        + "Please check the path: ",
                        rxn_library_file,
                    )
            else:
                raise Exception(
                    "rxn_library is not incorporated into smiles_click_chem.py"
                )

            # Import the proper reaction library JSON file
            try:
                with open(rxn_library_file, "r") as rxn_file:
                    reaction_dict_raw = json.load(rxn_file)
            except:
                raise Exception(
                    "rxn_library_file json file not able to be imported."
                    + " Check that the rxn_library is formatted correctly"
                )

        elif type(rxn_library_file) == str:
            if os.path.exists(rxn_library_file) is False:
                raise Exception(
                    "Custom specified rxn_library_file directory can not be found"
                )

            if os.path.isfile(rxn_library_file) is False:
                raise Exception(
                    "Custom specified rxn_library_file is not a file"
                )

            try:
                extension = os.path.splitext(rxn_library_file)[1]
            except:
                raise Exception(
                    "Custom specified rxn_library_file is not .json file."
                    + " It must be a .json dictionary"
                )

            if extension != ".json":
                raise Exception(
                    "Custom specified rxn_library_file is not .json file."
                    + " It must be a .json dictionary"
                )

            # Import the proper reaction library JSON file
            try:
                with open(rxn_library_file, "r") as rxn_file:
                    reaction_dict_raw = json.load(rxn_file)
            except:
                raise Exception(
                    "Custom specified rxn_library_file json file not able to "
                    + "be imported. Check that the rxn_library is "
                    + "formatted correctly"
                )

        else:
            raise Exception(
                "Custom specified rxn_library_file directory can not be found"
            )

        # Convert the reaction_dict_raw from unicode to the proper
        reaction_dict = self.rxn_lib_format_json_dict_of_dict(reaction_dict_raw)

        return reaction_dict

    def retrieve_functional_group_dict(self, rxn_library, function_group_library):
        """
        This retrieves a dictionary of all functional groups required for the
        respective reactions. This dictionary will be used to identify
        possible reactions.

        This is where all the functional groups which will be used in the
        SmartClickChem reactions are retrieved. If you want to add more just
        add a Custom set of reactions please add a folder to
        PATH/autogrow/operators/mutation/smiles_click_chem/Reaction_libraries/.
        They should be formatted as a dictionary of dictionary using the same
        format as :
        os.path.join(pwd,"reaction_libraries","click_chem_rxns",
                     "ClickChem_functional_groups.json")

        IF YOU CHOSE TO DO A Custom REACTION SET YOU MUST PROVIDE A DICTIONARY
        OF ALL FUNCTIONAL GROUPS IT WILL REACT. IF YOU FORGET TO ADD A
        FUNCTIONAL GROUP TO YOUR Custom DICTIONARY, THE REACTION MAY NEVER BE
        UTILIZED.

        Please note if your functional groups involve stereochemistry
            notations such as '\' please replace with '\\' (all functional
            groups should be formatted as SMARTS)

        Inputs:
        :param str rxn_library: A string defining the choice of the reaction
            library. ClickChem uses the set of reactions from Autogrow 3.1.2.
            Custom means you've defined a path to a Custom library in
            vars['function_group_library']
        :param str function_group_library: a PATH to a Custom functional group
            dictionary in a .json file. This will be a blank string if one choses
            a predefined functional groups option.

        Returns:
        :returns: dict functional_group_dict: A dictionary containing all
            SMARTS for identifying the functional groups for ClickChemistry
        """

        # Get the JSON file to import the proper reaction library
        pwd = os.path.dirname(__file__)

        if function_group_library == "":

            if rxn_library == "click_chem_rxns":
                function_group_library = os.path.join(
                    pwd, "reaction_libraries",
                    "click_chem_rxns",
                    "ClickChem_functional_groups.json",
                )
            elif rxn_library == "robust_rxns":
                function_group_library = os.path.join(
                    pwd, "reaction_libraries",
                    "robust_rxns",
                    "Robust_Rxns_functional_groups.json",
                )
            elif rxn_library == "all_rxns":
                function_group_library = os.path.join(
                    pwd, "reaction_libraries",
                    "all_rxns", "All_Rxns_functional_groups.json",
                )
            elif rxn_library == "Custom":
                if os.path.exists(function_group_library) is False:
                    raise Exception(
                        "Custom function_group_library cannot be found. "
                        + "Please check the path: ",
                        function_group_library,
                    )
            else:
                raise Exception(
                    "rxn_library is not incorporated into smiles_click_chem.py"
                )

            # Import the proper function_group_library JSON file
            try:
                with open(function_group_library, "r") as func_dict_file:
                    functional_group_dict_raw = json.load(func_dict_file)
            except:
                raise Exception(
                    "function_group_library json file not able to be imported. "
                    + "Check that the rxn_library is formatted correctly"
                )

        elif type(function_group_library) == str:
            if os.path.exists(function_group_library) is False:
                raise Exception(
                    "Custom specified function_group_library directory can not be found"
                )

            if os.path.isfile(function_group_library) is False:
                raise Exception("Custom specified function_group_library is not a file")

            try:
                extension = os.path.splitext(function_group_library)[1]
            except:
                raise Exception(
                    "Custom specified function_group_library is not .json "
                    + "file. It must be a .json dictionary"
                )

            if extension != ".json":
                raise Exception(
                    "Custom specified function_group_library is not .json "
                    + "file. It must be a .json dictionary"
                )

            # Import the proper function_group_library JSON file
            try:
                with open(function_group_library, "r") as func_dict_file:
                    functional_group_dict_raw = json.load(func_dict_file)
            except:
                raise Exception(
                    "function_group_library json file not able to be imported."
                    + " Check that the rxn_library is formatted correctly"
                )
        else:
            raise Exception(
                "Custom specified function_group_library directory can not be found"
            )

        # Convert the reaction_dict_raw from unicode to the proper
        functional_group_dict = self.rxn_lib_format_json_dict_of_dict(
            functional_group_dict_raw
        )

        return functional_group_dict

    def retrieve_complementary_dictionary(self, rxn_library, complementary_mol_dir):
        """
        Based on user controlled variables, this definition will retrieve a
        dictionary of molecules separated into classes by their functional
        groups. The sorting of a .smi file into this should be handled in the
        user parameter testing when autogrow is initially started.

        Inputs:
        :param str rxn_library: A string defining the choice of the reaction
            library. ClickChem uses the set of reactions from Autogrow 3.1.2.
            Custom means you've defined a path to a Custom library in
            vars['complementary_mol_dir']
        :param dict complementary_mol_dir: the path to the
            complementary_mol_dir directory. It may be an empty string in which
            case the complementary_mol_dir directory will default to those of the
            rxn_library

        Returns:
        :returns: dict complementary_mols_dict: a dictionary of complementary molecules
        """
        script_dir = os.path.dirname(os.path.realpath(__file__))

        if complementary_mol_dir == "":
            if rxn_library == "click_chem_rxns":
                complementary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "click_chem_rxns",
                    "complementary_mol_dir",
                )
            elif rxn_library == "robust_rxns":
                complementary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "robust_rxns",
                    "complementary_mol_dir",
                )
            elif rxn_library == "all_rxns":
                complementary_mol_dir = os.path.join(
                    script_dir,
                    "reaction_libraries",
                    "all_rxns",
                    "complementary_mol_dir",
                )
            elif rxn_library == "Custom":
                if os.path.isdir(complementary_mol_dir) is False:
                    raise Exception(
                        "Custom complementary_mol_dir cannot be found. "
                        + "Please check the path: ",
                        complementary_mol_dir,
                    )
            else:
                raise Exception(
                    "rxn_library is not incorporated into smiles_click_chem.py"
                )

        else:
            if os.path.isdir(complementary_mol_dir) is False:
                raise Exception(
                    "complementary_mol_dir is not a directory. It must be a \
                    directory with .smi files containing SMILES specified by \
                    functional groups.These .smi files must be named the same \
                    as the files in the complementary_mol_dir."
                )

        # Make a list of all the functional groups. These will be the name of
        # the .smi folders already separated by group.
        functional_groups = self.functional_group_dict.keys()

        missing_smi_files = []
        complementary_mols_dict = {}
        for group in functional_groups:
            filepath = "{}{}{}.smi".format(complementary_mol_dir, os.sep, group)

            if os.path.isfile(filepath) is True:
                complementary_mols_dict[group] = filepath

            else:
                missing_smi_files.append(filepath)
                print(
                    "Could not find the following .smi file for complementary "
                    + " molecules for Mutation: {}".format(filepath)
                )

        if len(missing_smi_files) != 0:
            raise Exception(
                "The following .smi file for complementary molecules "
                + "for Mutation is missing: ",
                missing_smi_files,
            )

        return complementary_mols_dict
#
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
#
def react_with_multiple_reactants(mol_tuple, mol_name, rxn_obj):
    """
    This will run a single molecule through a 1-reactant reaction.

    If it fails it will return the name of mol (mol_info[1])
    If it passes it will return None
    Inputs:
    :param tuple mol_tuple: a tuple of all mols to react
    :param str mol_name: name of the molecule being tested
    :param rdkit.Chem.rdChemReactions.ChemicalReaction rxn_obj: the reaction object to use

    Returns:
    :returns: str mol_name: returns the mol_name if it fails to react;
        returns None if it passes reaction
    """
    try:
        # if reaction works keep it
        reaction_products_list = [
                        x[0] for x in rxn_obj.RunReactants(mol_tuple)
        ]
    except:
        return mol_name

    if len(reaction_products_list) == 0:
        return mol_name
    # created a new compound so it passes
    return None
#
def run_a_single_reactant_reaction(mol_info, rxn_obj):
    """
    This will run a single molecule through a 1-reactant reaction.

    If it fails it will return the name of mol (mol_info[1])
    If it passes it will return None
    Inputs:
    :param list mol_info: list of mol info
        mol_info[0] is the SMILES,
        mol_info[1] is the name,
        mol_info[-1] is the rdkit mol obj,
    :param rdkit.Chem.rdChemReactions.ChemicalReaction rxn_obj: the reaction object to use

    Returns:
    :returns: str mol_name: returns the mol_name if it fails to react;
        returns None if it passes reaction
    """
    mol_name = mol_info[1]
    mol_1 = mol_info[-1]

    try:
        # if reaction works keep it
        reaction_products_list = rxn_obj.RunReactants((mol_1,))
    except:
        return mol_name
    if len(reaction_products_list) == 0:
        return mol_name
    # created a new compound so it passes
    return None
#
def get_rxn_and_examples(current_rxn_dict):
    """
    get the example reaction molecules from current_rxn_dict, create the rxn_obj,
    and test examples in the rxn.

    Inputs:
    :param dict current_rxn_dict: a dictionary of information about a reaction

    Returns:
    :returns: tuple example_rxn_reactants: a tuple of rdkit
            mol objects that are example compounds
    :returns: rdkit.Chem.rdChemReactions.ChemicalReaction rxn_obj: the
        reaction object to use
    """
    rxn_name = current_rxn_dict["reaction_name"]
    # Test example reactants
    example_smiles_rxn_reactants = current_rxn_dict["example_rxn_reactants"]
    example_smiles_rxn_reactants = example_smiles_rxn_reactants.replace("['", "").replace("']", "")
    example_smiles_rxn_reactants = example_smiles_rxn_reactants.replace(" ", "").replace('"', "")
    example_smiles_rxn_reactants = example_smiles_rxn_reactants.split("','")

    example_rxn_reactants = []
    for smile_str in example_smiles_rxn_reactants:
        smile_str = smile_str.replace("'", "").replace('"', "")
        smile_str = smile_str.replace(" ", "")

        example_mol = Chem.MolFromSmiles(smile_str)

        example_mol = MOH.check_sanitization(example_mol)
        if example_mol is None:
            print(smile_str)
            printout = "example mol from rxn: {}".format(rxn_name)
            printout = printout + " failed to sanitize in RDKit"
            print(printout)
            raise Exception(printout)
        example_rxn_reactants.append(example_mol)

    # convert example_rxn_reactants to tuple
    example_rxn_reactants = tuple(example_rxn_reactants)
    reaction_string = current_rxn_dict["reaction_string"]
    try:
        rxn_obj = AllChem.ReactionFromSmarts(reaction_string)
        rxn_obj.Initialize()
    except:
        printout = "rxn {} failed to be created.".format(rxn_name)
        printout = printout + "Rxn SMART is flawed"
        print(printout)
        raise Exception(printout)


    # Demo on example reactants
    example_results = react_with_multiple_reactants(example_rxn_reactants,
                                                    "test_reactions", rxn_obj)
    if example_results is not None:
        printout = "rxn {} failed to run on example compounds.".format(rxn_name)
        printout = printout + "\nPlease check example compounds"
        print(printout)
        raise Exception(printout)

    return example_rxn_reactants, rxn_obj
#
def run_all_for_fun_group(vars, fun_group, rxns_by_fun_group, a_smiles_click_object):
    """
    This runs the all testing for a single functional group.

    This will also write the compounds which pass to a .smi file.

    Inputs:
    :param dict vars: Dictionary of User variables
    :param str fun_group: functional group name
    :param dict rxns_by_fun_group: Dictionary of rxns names organized by
        functional groups
    :param obj a_smiles_click_object: a a_smiles_click_object class object.
        This provides useful pathing information.

    Returns:
    :returns: list failed_to_react: a list of mol names which failed to react
    :returns: list failed_to_sanitize: a list of mol names which failed to sanitize
    """
    # unpack variables
    complementary_mol_dict = a_smiles_click_object.complementary_mol_dict
    reaction_dict = a_smiles_click_object.reaction_dict
    number_of_processors = vars["number_of_processors"]
    output_folder = vars["output_folder"]

    smi_comp_file = complementary_mol_dict[fun_group]
    fun_group_list = get_usable_format(smi_comp_file)
    fun_group_mol_list = []
    failed_to_sanitize = []
    for info in fun_group_list:
        mol = Chem.MolFromSmiles(info[0])
        mol = MOH.check_sanitization(mol)
        if mol is None:
            failed_to_sanitize.append(info)
            continue
        temp = copy.deepcopy(info)
        temp.append(mol)
        fun_group_mol_list.append(temp)

    # print info about failures
    if len(failed_to_sanitize) != 0:
        printout = "{} compounds ".format(len(failed_to_sanitize))
        printout = printout + "failed to sanitize from: {}".format(fun_group)
        print(printout)

    failed_to_react = []
    for rxn_name in rxns_by_fun_group[fun_group]:

        current_rxn_dict = reaction_dict[rxn_name]
        example_reactants, rxn_obj = get_rxn_and_examples(current_rxn_dict)

        list_of_reactants = []
        functional_groups_rxn = current_rxn_dict["functional_groups"]
        i_count_to_use = None
        for i_count in range(len(functional_groups_rxn)):
            f_group = functional_groups_rxn[i_count]

            if fun_group == f_group:
                i_count_to_use = i_count
            else:
                continue
        if i_count_to_use is None:
            raise Exception("This is a code error.")

        list_of_reactants = []
        for mol_info in fun_group_mol_list:
            mol_tuple_temp = []
            for i_count in range(len(functional_groups_rxn)):
                if i_count == i_count_to_use:
                    mol_tuple_temp.append(mol_info[-1])
                else:
                    mol_tuple_temp.append(example_reactants[i_count])

            list_of_reactants.append(tuple([tuple(mol_tuple_temp), mol_info[1], rxn_obj]))

        output = mp.multi_threading(list_of_reactants, number_of_processors,
                                    react_with_multiple_reactants)
        output = [x for x in output if x is not None]
        failed_to_react.append([rxn_name, output])

        # print info about failures
        if len(output) != 0:
            printout = "{} compounds failed to react from ".format(len(output))
            printout = printout + "react from {} ".format(fun_group)
            printout = printout + "in rxn: {}".format(rxn_name)
            print(printout)

    master_failed_to_react = []
    master_passes_reactions = []
    for fail_mol_list in failed_to_react:
        master_failed_to_react.extend(fail_mol_list[1])
    for mol_info in fun_group_list:
        if mol_info[1] in master_failed_to_react:
            continue
        master_passes_reactions.append("    ".join(mol_info))
    # write to output .smi file
    with open(output_folder + fun_group + ".smi", "w") as f:
        f.write("\n".join(master_passes_reactions))

    return failed_to_react, failed_to_sanitize
#
def run_main(vars):
    """
    This runs the main testing.

    Inputs:
    :param dict vars: Dictionary of User variables
    """

    # Force rxn_library to be custom because why else run this
    rxn_library = "Custom"

    output_folder = vars["output_folder"]
    rxn_library_file = vars["rxn_library_file"]
    function_group_library = vars["function_group_library"]
    complementary_mol_dir = vars["complementary_mol_directory"]

    rxn_library_variables = [
        rxn_library,
        rxn_library_file,
        function_group_library,
        complementary_mol_dir
    ]
    new_mutation_smiles_list = []

    a_smiles_click_chem_object = SmilesClickChem(
        rxn_library_variables, new_mutation_smiles_list
    )

    list_of_reaction_names = a_smiles_click_chem_object.list_of_reaction_names
    functional_group_dict = a_smiles_click_chem_object.functional_group_dict
    reaction_dict = a_smiles_click_chem_object.reaction_dict

    rxns_by_fun_group = {}
    for fun_group in functional_group_dict.keys():
        rxns_by_fun_group[fun_group] = []
    for rxn_name in list_of_reaction_names:
        current_rxn_dict = reaction_dict[rxn_name]
        for fun_group in current_rxn_dict["functional_groups"]:
            temp_list = rxns_by_fun_group[fun_group]
            temp_list.append(rxn_name)
            rxns_by_fun_group[fun_group] = temp_list

    failed_to_sanitize_by_fun_group = {}
    failed_to_react_by_fun_group = {}

    for fun_group in rxns_by_fun_group.keys():
        failed_to_react, failed_to_sanitize = run_all_for_fun_group(vars, fun_group,
                                                                    rxns_by_fun_group,
                                                                    a_smiles_click_chem_object)
        failed_to_react_by_fun_group[fun_group] = failed_to_react
        failed_to_sanitize_by_fun_group[fun_group] = failed_to_sanitize

    # Handle saving log
    with open(output_folder + "failed_to_sanitize_mol_by_fun_group.json", "w") as fp:
        json.dump(failed_to_sanitize_by_fun_group, fp, indent=4)

    with open(output_folder + "failed_to_react_by_fun_group.json", "w") as fp:
        json.dump(failed_to_react_by_fun_group, fp, indent=4)

    master_failed_list = []
    for fun_group in failed_to_react_by_fun_group.keys():
        temp = [x[1] for x in failed_to_react_by_fun_group[fun_group]]
        for x in temp:
            master_failed_list.extend(x)
    master_failed_list = list(set(master_failed_list))
    if len(master_failed_list) == 0:
        print("All compounds passed!")
    else:
        print("{} compounds failed. Please check logs".format(len(master_failed_list)))
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
    if args_dict["rxn_library_file"] == "" or args_dict["function_group_library"] == "":
        raise ValueError(
            "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY \
                THE PATH TO THE REACTION LIBRARY USING INPUT PARAMETER rxn_library"
        )
    if os.path.exists(args_dict["rxn_library_file"]) is False:
        raise ValueError(
            "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY \
            THE PATH TO THE REACTION LIBRARY USING INPUT PARAMETER rxn_library"
        )

    if args_dict["complementary_mol_directory"] == "":
        raise ValueError(
            "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY THE PATH \
            TO THE REACTION LIBRARY USING INPUT PARAMETER function_group_library"
        )
    if os.path.isdir(args_dict["complementary_mol_directory"]) is False:
        raise ValueError(
            "TO USE Custom REACTION LIBRARY OPTION, ONE MUST SPECIFY THE PATH \
            TO THE REACTION LIBRARY USING INPUT PARAMETER complementary_mol_directory"
        )

    if "number_of_processors" not in args_dict.keys():
        args_dict["number_of_processors"] = -1
    try:
        args_dict["number_of_processors"] = int(args_dict["number_of_processors"])
    except:
        raise ValueError(
            "number_of_processors must be an int. \
            To use all processors set to -1.")

    if "output_folder" not in args_dict.keys():
        printout = "output_folder is a required variable. it is the PATH to where " + \
            "filtered .smi file and log files will be placed. Will save a file " + \
            "in this directory for mols which failed sanitization, mols which " + \
            "failed to react in specific reactions, and .smi files that contain " + \
            "all mols that reacted properly."
        raise ValueError(printout)

    if type(args_dict["output_folder"]) != str or args_dict["output_folder"] == "":
        printout = "output_folder is a required variable. it is the PATH to where " + \
            "filtered .smi file and log files will be placed. Will save a file " + \
            "in this directory for mols which failed sanitization, mols which " + \
            "failed to react in specific reactions, and .smi files that contain " + \
            "all mols that reacted properly."
        raise ValueError(printout)

    args_dict["output_folder"] = os.path.abspath(args_dict["output_folder"]) + os.sep
    if os.path.exists(args_dict["output_folder"]) is True:
        if os.path.isdir(args_dict["output_folder"]) is False:
            print(args_dict["output_folder"])
            printout = "output_folder must be a directory. Please check input arguments"
            raise ValueError(printout)
    else:
        try:
            os.mkdir(args_dict["output_folder"])
        except:
            pass
        if os.path.exists(args_dict["output_folder"]) is False:
            raise Exception("output_folder could not be made or found.")

    return args_dict
#


# Argument parsing
PARSER = argparse.ArgumentParser()
# Mutation Settings
PARSER.add_argument(
    "--rxn_library_file",
    type=str,
    default="",
    required=True,
    help="This PATH to a Custom json file of SMARTS reactions to use for Mutation."
)
PARSER.add_argument(
    "--function_group_library",
    type=str,
    default="",
    required=True,
    help="This PATH for a dictionary of functional groups to be used for Mutation.",
)
PARSER.add_argument(
    "--complementary_mol_directory",
    type=str,
    default="",
    required=True,
    help="This PATH to the directory containing all the molecules being used \
    to react with. The directory should contain .smi files contain SMILES of \
    molecules containing the functional group represented by that file. Each file \
    should be named with the same title as the functional groups described in \
    rxn_library_file & function_group_library +.smi \
    All Functional groups specified function_group_library must have its \
    own .smi file. We recommend you filter these dictionaries prior to Autogrow \
    for the Drug-likeliness and size filters you will Run Autogrow with.",
)
PARSER.add_argument(
    "--output_folder",
    type=str,
    default="",
    required=True,
    help="This PATH to where filtered .smi file and log files will be placed. \
        Will save a file in this directory for mols which failed sanitization, \
        mols which failed to react in specific reactions, and .smi files \
        that contain all mols that reacted properly.",
)
# processors and multithread mode
PARSER.add_argument(
    "--number_of_processors",
    "-p",
    type=int,
    default=-1,
    help="Number of processors to use for parallel calculations. \
    Set to -1 for all available CPUs.",
)



ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)
run_main(ARGS_DICT)
print("done")
