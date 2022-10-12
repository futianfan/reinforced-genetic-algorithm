import argparse

PARSER = argparse.ArgumentParser()

# Allows the run commands to be submitted via a .json file.
PARSER.add_argument(
    "--json",
    "-j",
    metavar="param.json",
    help="Name of a json file containing all parameters. \
    Overrides other arguments.",
)

# Allows the run in debug mode. Doesn't delete temp files.
PARSER.add_argument(
    "--debug_mode",
    "-d",
    action="store_true",
    default=False,
    help="Run Autogrow in Debug mode. This keeps all \
    temporary files and adds extra print statements.",
)

# receptor information
PARSER.add_argument(
    "--filename_of_receptor",
    "-r",
    metavar="receptor.pdb",
    default='./tutorial/PARP/4r6eA_PARP1_prepared.pdb', 
    help="The path to the receptor file. Should be .pdb file.",
)
PARSER.add_argument(
    "--center_x",
    "-x",
    type=float,
    default=-70.76,
    help="x-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
)
PARSER.add_argument(
    "--center_y",
    "-y",
    type=float,
    default=21.82,
    help="y-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
)
PARSER.add_argument(
    "--center_z",
    "-z",
    type=float,
    default=28.33,
    help="z-coordinate for the center of the pocket to be tested by docking. (Angstrom)",
)

PARSER.add_argument(
    "--size_x",
    type=float,
    default=25.0,
    help="dimension of box to dock into in the x-axis (Angstrom)",
)
PARSER.add_argument(
    "--size_y",
    type=float,
    default=20.0,
    help="dimension of box to dock into in the y-axis (Angstrom)",
)
PARSER.add_argument(
    "--size_z",
    type=float,
    default=25.0,
    help="dimension of box to dock into in the z-axis (Angstrom)",
)


# Input/Output directories
PARSER.add_argument(
    "--root_output_folder",
    "-o",
    type=str,
    help="The Path to the folder which all output files will be placed.",
)
PARSER.add_argument(
    "--source_compound_file",
    "-s",
    type=str,
    default='./source_compounds/naphthalene_smiles.smi',
    help="PATH to the file containing the source compounds. It must be \
    tab-delineated .smi file. These ligands will seed the first generation.",
)
PARSER.add_argument(
    "--filter_source_compounds",
    choices=[True, False, "True", "False", "true", "false"],
    default=True,
    help="If True source ligands from source_compound_file will be \
    filter using the user defined filter choices prior to the 1st generation being \
    created. If False, ligands which would fail the ligand filters could seed \
    the 1st generation. Default is True.",
)
PARSER.add_argument(
    "--use_docked_source_compounds",
    choices=[True, False, "True", "False", "true", "false"],
    default=False,
    help="If True source ligands will be docked prior to seeding generation 1. \
    If True and the source_compound file already has docking/fitness metric score \
    in -2 column of .smi file, it will not redock but reuse the scores from \
    the source_compound_file.\
    If True and no fitness metric score in -2 column of .smi file, it will \
    dock each ligand from the source_compound_file and displayed as generation 0.\
    If False, generation 1 will be randomly seeded by the source compounds with \
    no preference and there will be no generation 0. \
    If performing multiple simulations using same source compounds and protein, \
    we recommend running once this and using the generation 0 ranked file as the \
    source_compound_file for future simulations. \
    Default is True.",
)
PARSER.add_argument(
    "--start_a_new_run",
    action="store_true",
    default=False,
    help="If False make a new folder and start a fresh simulation with Generation 0.  \
    If True find the last generation in the root_output_folder and continue to fill.\
    Default is False.",
)



# SmilesMerge Settings
PARSER.add_argument(
    "--max_time_MCS_prescreen",
    type=int,
    default=1,
    help="amount time the pre-screen MCS times out. Time out doesnt prevent \
    mcs matching just takes what it has up to that point",
)
PARSER.add_argument(
    "--max_time_MCS_thorough",
    type=int,
    default=1,
    help="amount time the thorough MCS times out. Time out doesnt prevent \
    mcs matching just takes what it has up to that point",
)
PARSER.add_argument(
    "--min_atom_match_MCS",
    type=int,
    default=4,
    help="Determines the minimum number of atoms in common for a substructurematch. \
    The higher the more restrictive, but the more likely for two ligands not to match",
)
PARSER.add_argument(
    "--protanate_step",
    action="store_true",
    default=False,
    help="Indicates if Smilesmerge uses protanated mols (if true) or deprot \
    (if False) SmilesMerge is 10x faster when deprotanated",
)


# Mutation Settings
PARSER.add_argument(
    "--rxn_library",
    choices=["click_chem_rxns", "robust_rxns", "all_rxns", "Custom"],
    default="all_rxns",
    help="This set of reactions to be used in Mutation. \
    If Custom, one must also provide rxn_file Path and function_group_library path",
)
PARSER.add_argument(
    "--rxn_library_file",
    type=str,
    default="",
    help="This PATH to a Custom json file of SMARTS reactions to use for Mutation. \
    Only provide if using the Custom option for rxn_library.",
)
PARSER.add_argument(
    "--function_group_library",
    type=str,
    default="",
    help="This PATH for a dictionary of functional groups to be used for Mutation. \
    Only provide if using the Custom option for rxn_library.",
)
PARSER.add_argument(
    "--complementary_mol_directory",
    type=str,
    default="",
    help="This PATH to the directory containing all the molecules being used \
    to react with. The directory should contain .smi files contain SMILES of \
    molecules containing the functional group represented by that file. Each file \
    should be named with the same title as the functional groups described in \
    rxn_library_file & function_group_library +.smi \
    All Functional groups specified function_group_library must have its \
    own .smi file. We recommend you filter these dictionaries prior to Autogrow \
    for the Drug-likeliness and size filters you will Run Autogrow with.",
)


# processors and multithread mode
PARSER.add_argument(
    "--number_of_processors",
    "-p",
    type=int,
    metavar="N",
    default=1,
    help="Number of processors to use for parallel calculations. Set to -1 for all available CPUs.",
)
PARSER.add_argument(
    "--multithread_mode",
    default="multithreading",
    choices=["mpi", "multithreading", "serial"],
    help="Determine what style \
    multithreading: mpi, multithreading, or serial. serial will override \
    number_of_processors and force it to be on a single processor.",
)

# Genetic Algorithm Options
PARSER.add_argument(
    "--selector_choice",
    choices=["Roulette_Selector", "Rank_Selector", "Tournament_Selector"],
    default="Roulette_Selector",
    help="This determines whether the fitness criteria are chosen by a Weighted Roulette, \
    Ranked, or Tournament style Selector. The Rank option is a non-redundant selector.\
    Roulette and Tournament chose without replacement and are stoichastic options. \
    Warning do not use Rank_Selector for small runs as there is potential that \
    the number of desired ligands exceed the number of ligands to chose from.",
)
PARSER.add_argument(
    "--tourn_size",
    type=float,
    default=0.1,
    help="If using the Tournament_Selector this determines the size of each \
    tournament. The number of ligands used for each tournament will the \
    tourn_size * the number of considered ligands.",
)

# Seeding next gen and diversity
PARSER.add_argument(
    "--top_mols_to_seed_next_generation_first_generation",
    type=int,
    help="Number of mols that seed next generation, for the first generation.\
    Should be less than number_of_crossovers_first_generation + number_of_mutations_first_generation\
    If not defined it will default to top_mols_to_seed_next_generation",
)
PARSER.add_argument(
    "--top_mols_to_seed_next_generation",
    type=int,
    default=10,
    help="Number of mols that seed next generation, for all generations after the first.\
    Should be less than number_of_crossovers_first_generation \
    + number_of_mutations_first_generation",
)
PARSER.add_argument(
    "--diversity_mols_to_seed_first_generation",
    type=int,
    default=10,
    help="Should be less than number_of_crossovers_first_generation \
    + number_of_mutations_first_generation",
)
PARSER.add_argument(
    "--diversity_seed_depreciation_per_gen",
    type=int,
    default=2,
    help="Each gen diversity_mols_to_seed_first_generation will decrease this amount",
)

# Populations settings
PARSER.add_argument(
    "--num_generations",
    type=int,
    default=10,
    help="The number of generations to be created.",
)
PARSER.add_argument(
    "--number_of_crossovers_first_generation",
    type=int,
    help="The number of ligands which will be created via crossovers in the \
    first generation. If not defined it will default to number_of_crossovers",
)
PARSER.add_argument(
    "--number_of_mutants_first_generation",
    type=int,
    help="The number of ligands which will be created via mutation in \
    the first generation. If not defined it will default to number_of_mutants",
)
PARSER.add_argument(
    "--number_elitism_advance_from_previous_gen_first_generation",
    type=int,
    help="The number of ligands chosen for elitism for the first generation \
    These will advance from the previous generation directly into the next \
    generation.  This is purely advancing based on Docking/Rescore fitness. \
    This does not select for diversity. If not defined it will default to \
    number_elitism_advance_from_previous_gen",
)
PARSER.add_argument(
    "--number_of_crossovers",
    type=int,
    default=10,
    help="The number of ligands which will be created via crossover in each \
    generation besides the first",
)
PARSER.add_argument(
    "--number_of_mutants",
    type=int,
    default=10,
    help="The number of ligands which will be created via mutation in each \
    generation besides the first.",
)
PARSER.add_argument(
    "--number_elitism_advance_from_previous_gen",
    type=int,
    default=10,
    help="The number of ligands chosen for elitism. These will advance from \
    the previous generation directly into the next generation. \
    This is purely advancing based on Docking/Rescore \
    fitness. This does not select for diversity.",
)
PARSER.add_argument(
    "--redock_elite_from_previous_gen",
    choices=[True, False, "True", "False", "true", "false"],
    default=False,
    help="If True than ligands chosen via Elitism (ie advanced from last generation) \
    will be passed through Gypsum and docked again. This provides a better exploration of conformer space \
    but also requires more computation time. If False, advancing ligands are simply carried forward by \
    copying the PDBQT files.",
)

####### FILTER VARIABLES
PARSER.add_argument(
    "--LipinskiStrictFilter",
    action="store_true",
    default=False,
    help="Lipinski filters for orally available drugs following Lipinski rule of fives. \
    Filters by molecular weight, logP and number of hydrogen bond donors and acceptors. \
    Strict implementation means a ligand must pass all requirements.",
)
PARSER.add_argument(
    "--LipinskiLenientFilter",
    action="store_true",
    default=False,
    help="Lipinski filters for orally available drugs following Lipinski rule of fives. \
    Filters by molecular weight, logP and number of hydrogen bond donors and acceptors. \
    Lenient implementation means a ligand may fail all but one requirement and still passes.",
)
PARSER.add_argument(
    "--GhoseFilter",
    action="store_true",
    default=False,
    help="Ghose filters for drug-likeliness; filters by molecular weight,\
    logP and number of atoms.",
)
PARSER.add_argument(
    "--GhoseModifiedFilter",
    action="store_true",
    default=False,
    help="Ghose filters for drug-likeliness; filters by molecular weight,\
    logP and number of atoms. This is the same as the GhoseFilter, but \
    the upper-bound of the molecular weight restrict is loosened from \
    480Da to 500Da. This is intended to be run with Lipinski Filter and \
    to match AutoGrow 3's Ghose Filter.",
)
PARSER.add_argument(
    "--MozziconacciFilter",
    action="store_true",
    default=False,
    help="Mozziconacci filters for drug-likeliness; filters by the number of \
    rotatable bonds, rings, oxygens, and halogens.",
)
PARSER.add_argument(
    "--VandeWaterbeemdFilter",
    action="store_true",
    default=False,
    help="VandeWaterbeemd filters for drug likely to be blood brain barrier permeable. \
    Filters by the number of molecular weight and Polar Sureface Area (PSA).",
)
PARSER.add_argument(
    "--PAINSFilter",
    action="store_true",
    default=False,
    help="PAINS filters against Pan Assay Interference Compounds using \
    substructure a search.",
)
PARSER.add_argument(
    "--NIHFilter",
    action="store_true",
    default=False,
    help="NIH filters against molecules with undersirable functional groups \
    using substructure a search.",
)
PARSER.add_argument(
    "--BRENKFilter",
    action="store_true",
    default=False,
    help="BRENK filter for lead-likeliness, by matching common false positive \
    molecules to the current mol.",
)
PARSER.add_argument(
    "--No_Filters",
    action="store_true",
    default=False,
    help="No filters will be applied to compounds.",
)
PARSER.add_argument(
    "--alternative_filter",
    action="append",
    help="If you want to add Custom filters to the filter child classes \
    Must be a list of lists \
    [[name_filter1, Path/to/name_filter1.py],[name_filter2, Path/to/name_filter2.py]]",
)

# dependency variables
# DOCUMENT THE file conversion for docking inputs
PARSER.add_argument(
    "--conversion_choice",
    choices=["MGLToolsConversion", "ObabelConversion", "Custom"],
    default="MGLToolsConversion",
    help="Determines how .pdb files will be converted \
    to the final format for docking. For Autodock Vina and QuickVina style docking software, \
    files must be in .pdbqt format. MGLToolsConversion: uses MGLTools and is the \
    recommended converter. MGLTools conversion is required for NNScore1/2 rescoring. \
    ObabelConversion: uses commandline obabel. Easier to install but Vina docking has \
    been optimized with MGLTools conversion.",
)
PARSER.add_argument(
    "--custom_conversion_script",
    metavar="custom_conversion_script",
    default="",
    help="The path to a python script for which is used to convert \
    ligands. This is required for custom conversion_choice choices. \
    Must be a list of strings \
    [name_custom_conversion_class, Path/to/name_custom_conversion_class.py]",
)
PARSER.add_argument(
    "--mgltools_directory",
    metavar="mgltools_directory",
    help="Required if using MGLTools conversion option \
    (conversion_choice=MGLToolsConversion) \
    Path may look like: /home/user/MGLTools-1.5.6/",
)
PARSER.add_argument(
    "--mgl_python",
    metavar="mgl_python",
    required=False,
    help="/home/user/MGLTools-1.5.4/bin/pythonsh",
)
PARSER.add_argument(
    "--prepare_ligand4.py",
    metavar="prepare_ligand4.py",
    required=False,
    help="/home/user/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py",
)
PARSER.add_argument(
    "--prepare_receptor4.py",
    metavar="prepare_receptor4.py",
    required=False,
    help="/home/userMGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py",
)
PARSER.add_argument(
    "--obabel_path",
    help="required if using obabel conversion \
    option (conversion_choice=ObabelConversion).\
    Path may look like PATH/envs/py37/bin/obabel; \
    may be found on Linux by running: which obabel",
)


###################################
######### docking #################
###################################
PARSER.add_argument(
    "--dock_choice",
    metavar="dock_choice",
    default="QuickVina2Docking",
    choices=["VinaDocking", "QuickVina2Docking", "Custom"],
    help="dock_choice assigns which docking software module to use.",
)
PARSER.add_argument(
    "--docking_executable",
    metavar="docking_executable",
    default=None,
    help="path to the docking_executable",
)
PARSER.add_argument(
    "--docking_exhaustiveness",
    metavar="docking_exhaustiveness",
    default=None,
    help="exhaustiveness of the global search (roughly proportional to time. \
    see docking software for settings. Unless specified Autogrow uses the \
    docking softwares default setting. For AutoDock Vina 1.1.2 that is 8",
)
PARSER.add_argument(
    "--docking_num_modes",
    metavar="docking_num_modes",
    default=None,
    help=" maximum number of binding modes to generate in docking. \
    See docking software for settings. Unless specified Autogrow uses the \
    docking softwares default setting. For AutoDock Vina 1.1.2 that is 9",
)
PARSER.add_argument(
    "--docking_timeout_limit",
    type=float,
    default=120,
    help="The maximum amount of time allowed to dock a single ligand into a \
    pocket in seconds. Many factors influence the time required to dock, such as: \
    processor speed, the docking software, rotatable bonds, exhaustiveness docking,\
    and number of docking modes... \
    The default docking_timeout_limit is 120 seconds, which is excess for most \
    docking events using QuickVina2Docking under default settings. If run with \
    more exhaustive settings or with highly flexible ligands, consider increasing \
    docking_timeout_limit to accommodate. Default docking_timeout_limit is 120 seconds",
)
PARSER.add_argument(
    "--custom_docking_script",
    metavar="custom_docking_script",
    default="",
    help="The name and path to a python script for which is used to \
    dock ligands. This is required for Custom docking choices Must be a list of \
    strings [name_custom_conversion_class, Path/to/name_custom_conversion_class.py]",
)

# scoring
PARSER.add_argument(
    "--scoring_choice",
    metavar="scoring_choice",
    choices=["VINA", "NN1", "NN2", "Custom"],
    default="VINA",
    help="The scoring_choice to use to assess the ligands docking fitness. \
    Default is using Vina/QuickVina2 ligand affinity while NN1/NN2 use a Neural Network \
    to assess the docking pose. Custom requires providing a file path for a Custom \
    scoring function. If Custom scoring function, confirm it selects properly, \
    Autogrow is largely set to select for a more negative score.",
)
PARSER.add_argument(
    "--rescore_lig_efficiency",
    action="store_true",
    default=False,
    help="This will divide the final scoring_choice output by the number of \
    non-Hydrogen atoms in the ligand. This adjusted ligand efficiency score will \
    override the scoring_choice value. This is compatible with all scoring_choice options.",
)
PARSER.add_argument(
    "--custom_scoring_script",
    metavar="custom_scoring_script",
    type=str,
    default="",
    help="The path to a python script for which is used to \
    assess the ligands docking fitness. Autogrow is largely set to select for a most \
    negative scores (ie binding affinity the more negative is best). Must be a list of \
    strings [name_custom_conversion_class, Path/to/name_custom_conversion_class.py]",
)

# gypsum # max variance is the number of conformers made per ligand
PARSER.add_argument(
    "--max_variants_per_compound",
    type=int,
    default=3,
    help="number of conformers made per ligand. \
    See Gypsum-DL publication for details",
)
PARSER.add_argument(
    "--gypsum_thoroughness",
    "-t",
    type=int,
    default = 3, 
    help="How widely Gypsum-DL will search for \
    low-energy conformers. Larger values increase \
    run times but can produce better results. \
    See Gypsum-DL publication for details",
)
PARSER.add_argument(
    "--min_ph",
    metavar="MIN",
    type=float,
    default=6.4,
    help="Minimum pH to consider.See Gypsum-DL \
    and Dimorphite-D publication for details.",
)
PARSER.add_argument(
    "--max_ph",
    metavar="MAX",
    type=float,
    default=8.4,
    help="Maximum pH to consider.See Gypsum-DL \
    and Dimorphite-D publication for details.",
)
PARSER.add_argument(
    "--pka_precision",
    metavar="D",
    type=float,
    default=1.0,
    help="Size of pH substructure ranges. See Dimorphite-DL \
    publication for details.",
)
PARSER.add_argument(
    "--gypsum_timeout_limit",
    type=float,
    default=15,
    help="Maximum time gypsum is allowed to run for a given ligand in seconds. \
    On average Gypsum-DL takes on several seconds to run for a given ligand, but \
    factors such as mol size, rotatable bonds, processor speed, and gypsum \
    settings (ie gypsum_thoroughness or max_variants_per_compound) will change \
    how long it takes to run. If increasing gypsum settings it is best to increase \
    the gypsum_timeout_limit. Default gypsum_timeout_limit is 15 seconds",
)

# Reduce files down. This compiles and compresses the files in the PDBs folder
# (contains docking outputs, pdb, pdbqt...). This reduces the data size and
# makes data transfer quicker, but requires running the
# file_concatenation_and_compression.py in the Utility script folder to
# separate these files out for readability.
PARSER.add_argument(
    "--reduce_files_sizes",
    choices=[True, False, "True", "False", "true", "false"],
    default=True,
    help="Run this combines all files in the PDBs folder into a \
    single text file. Useful when data needs to be transferred.",
)

# Make a line plot of the simulation at the end of the run.
PARSER.add_argument(
    "--generate_plot",
    choices=[True, False, "True", "False", "true", "false"],
    default=True,
    help="Make a line plot of the simulation at the end of the run.",
)

# mpi mode pre-Run so there are python cache files without EOF Errors
PARSER.add_argument(
    "--cache_prerun",
    "-c",
    action="store_true",
    help="Run this before running gypsum in mpi-mode.",
)


args_dict = vars(PARSER.parse_args())
from autogrow.user_vars import multiprocess_handling, define_defaults, determine_bash_timeout_vs_gtimeout
# args_dict = define_defaults()


import copy 
INPUTS = copy.deepcopy(args_dict)

for k, v in args_dict.items():
    if v is None:
        del INPUTS[k]

if args_dict["cache_prerun"] is False:
    # load the commandline parameters
    from autogrow.user_vars import load_in_commandline_parameters
    args_dict, printout = load_in_commandline_parameters(INPUTS)



args_dict = multiprocess_handling(args_dict)

timeout_option = determine_bash_timeout_vs_gtimeout()
if timeout_option in  ["timeout", "gtimeout"]:
    args_dict["timeout_vs_gtimeout"] = timeout_option
else:
    raise Exception("Something is very wrong. This OS may not be supported by \
                    Autogrow or you may need to execute through Bash.")



import numpy as np
import os, json, pickle, time, sys  


file = 'source_compounds/naphthalene_smiles.smi'
with open(file, 'r') as fin:
    smiles_list = fin.readlines() 
smiles_list = [smiles.split()[0] for smiles in smiles_list]
vars = args_dict


##########################################################
# A. mutate 
##########################################################
"""
import autogrow.operators.mutation.smiles_click_chem.smiles_click_chem as SmileClickClass

rxn_library_variables = [
    vars["rxn_library"],
    vars["rxn_library_file"],
    vars["function_group_library"],
    vars["complementary_mol_directory"],
]    # Package user vars specifying the Reaction library to use for mutation 

new_mutation_smiles_list = [] # List of SMILES from mutation
a_smiles_click_chem_object = SmileClickClass.SmilesClickChem(rxn_library_variables, new_mutation_smiles_list, vars["filter_object_dict"])

for smiles in smiles_list:
    result_of_run = a_smiles_click_chem_object.run_smiles_click2(smiles)
    print('-------------', result_of_run)
"""

##########################################################
### B. crossover 
##########################################################
## crossover between 2 ligands 
''' crossover: find common structure '''
import autogrow.operators.crossover.smiles_merge.smiles_merge as smiles_merge 
import autogrow.operators.crossover.execute_crossover as execute_crossover
import autogrow.operators.filter.execute_filters as Filter
# execute_crossover.test_for_mcs
# test_for_mcs(vars, mol_1, mol_2)
# execute_crossover.convert_mol_from_smiles 
for smiles in smiles_list:
    mol = execute_crossover.convert_mol_from_smiles(smiles)
    for s2 in smiles_list:
        if smiles == s2: 
            continue 
        mol2 = execute_crossover.convert_mol_from_smiles(s2)
        if execute_crossover.test_for_mcs(vars, mol, mol2) is None:
            continue 
        for i in range(3):
            ligand_new_smiles = smiles_merge.run_main_smiles_merge(vars, smiles, s2)
            if ligand_new_smiles is not None:
                break 
        if ligand_new_smiles is None:
            print('====fail to merge')
        else:
            pass_or_not = Filter.run_filter_on_just_smiles(ligand_new_smiles, vars["filter_object_dict"])  #### True, False
            if pass_or_not:
                print('-----success----', ligand_new_smiles)


'''
python demo_GAoperation.py \
    --filename_of_receptor ./tutorial/PARP/4r6eA_PARP1_prepared.pdb \
    --center_x -70.76 --center_y  21.82 --center_z 28.33 \
    --size_x 25.0 --size_y 16.0 --size_z 25.0 \
    --source_compound_file ./source_compounds/naphthalene_smiles.smi \
    --root_output_folder ./output \
    --number_of_mutants_first_generation 50 \
    --number_of_crossovers_first_generation 50 \
    --number_of_mutants 50 \
    --number_of_crossovers 50 \
    --top_mols_to_seed_next_generation 50 \
    --number_elitism_advance_from_previous_gen 50 \
    --number_elitism_advance_from_previous_gen_first_generation 10 \
    --diversity_mols_to_seed_first_generation 10 \
    --diversity_seed_depreciation_per_gen 10 \
    --num_generations 5 \
    --mgltools_directory ./mgltools_x86_64Linux2_1.5.6/ \
    --number_of_processors -1 \
    --scoring_choice VINA \
    --LipinskiLenientFilter \
    --start_a_new_run \
    --rxn_library click_chem_rxns \
    --selector_choice Rank_Selector \
    --dock_choice VinaDocking \
    --max_variants_per_compound 5 \
    --redock_elite_from_previous_gen False \
    --generate_plot True \
    --reduce_files_sizes True \
    --use_docked_source_compounds True  
'''











