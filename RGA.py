'''
- import and config 
- policy network 
- for i in 1,...,generation
	- crossover
		- use policy network to select ligand *** 
		- crossover
		- docking 
		- update policy network *** 
	- mutation 
		- use policy network to select ligand *** 
		- mutation 
		- docking 
		- update policy network ***
	- elites 

'''

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

import numpy as np
from tqdm import tqdm 
from collections import defaultdict 
import os, json, pickle, time, sys, copy, random  
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
if timeout_option in ["timeout", "gtimeout"]:
    args_dict["timeout_vs_gtimeout"] = timeout_option
else:
    raise Exception("Something is very wrong. This OS may not be supported by Autogrow or you may need to execute through Bash.")

vars = args_dict
topk = 10

############## canonical ##############
from rdkit import Chem 
def canonicalize(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if mol is not None:
    return Chem.MolToSmiles(mol, isomericSmiles=True)
  else:
    return None 
##########################################################
### A. mutate 
import autogrow.operators.mutation.smiles_click_chem.smiles_click_chem as SmileClickClass
rxn_library_variables = [vars["rxn_library"], vars["rxn_library_file"], vars["function_group_library"], vars["complementary_mol_directory"]] # Package user vars specifying Reaction library for mutation 
new_mutation_smiles_list = [] # List of SMILES from mutation
a_smiles_click_chem_object = SmileClickClass.SmilesClickChem(rxn_library_variables, new_mutation_smiles_list, vars["filter_object_dict"])
##########################################################
##########################################################
### B. crossover 
## crossover between 2 ligands: find common structure 
import autogrow.operators.crossover.smiles_merge.smiles_merge as smiles_merge 
import autogrow.operators.crossover.execute_crossover as execute_crossover
import autogrow.operators.filter.execute_filters as Filter
##########################################################
##########################################################
# C. smiles2docking
#	- smiles2pdbqt: smiles -> sdf -> pdb -> pdbqt 
#	- docking pdbqt
###########  smiles -> sdf -> pdb -> pdbqt -> pdbqt.vina  
import autogrow.operators.convert_files.conversion_to_3d as conversion_to_3d
# conversion_to_3d.convert_smi_to_sdfs_with_gypsum 
# conversion_to_3d.convert_sdf_to_pdbs
# convert_sdf_to_pdbs(vars, gen_folder_path, sdfs_folder_path)
# conversion_to_3d.convert_single_sdf_to_pdb
# convert_ligand_pdb_file_to_pdbqt #### in run_docking_common lig_convert_multithread 

def smiles_to_sdfs(vars, gen_smiles_file, smile_file_directory):
    # adapted from conversion_to_3d.convert_smi_to_sdfs_with_gypsum 
    max_variants_per_compound = vars["max_variants_per_compound"]
    gypsum_thoroughness = vars["gypsum_thoroughness"]
    min_ph = vars["min_ph"]
    max_ph = vars["max_ph"]
    pka_precision = vars["pka_precision"]
    gypsum_timeout_limit = vars["gypsum_timeout_limit"]

    # Make a new folder to put gypsum .smi's and json. Name folder gypsum_submission_files.
    folder_path = "{}gypsum_submission_files{}".format(smile_file_directory, os.sep)
    if os.path.exists(folder_path) is False:
        os.makedirs(folder_path)

    # Make Output for Gypsum folder (where .sdf's go)
    gypsum_output_folder_path = "{}_SDF{}".format(smile_file_directory, os.sep)
    if os.path.exists(gypsum_output_folder_path) is False:
        os.makedirs(gypsum_output_folder_path)

    # Make a folder to put the log files into within the 3D_SDFs folder
    gypsum_log_path = "{}log{}".format(gypsum_output_folder_path, os.sep)
    if os.path.exists(gypsum_log_path) is False:
        os.makedirs(gypsum_log_path)

    # Make All of the json files to submit to gypsum
    list_of_gypsum_params = conversion_to_3d.make_smi_and_gyspum_params(
        gen_smiles_file,
        folder_path,
        gypsum_output_folder_path,
        max_variants_per_compound, gypsum_thoroughness,
        min_ph, max_ph, pka_precision, )

    # create a the job_inputs to run gypsum in multithread
    job_input = tuple([(gypsum_log_path, gypsum_params, gypsum_timeout_limit) for gypsum_params in list_of_gypsum_params])

    sys.stdout.flush()
    failed_to_convert = vars["parallelizer"].run(job_input, conversion_to_3d.run_gypsum_multiprocessing)
    sys.stdout.flush()

    ###    fail: return smiles 
    ###    success: return None     
    lig_failed_to_convert = [x for x in failed_to_convert if x is not None]
    lig_failed_to_convert = list(set(lig_failed_to_convert))
    if len(lig_failed_to_convert) > 0:
        print("The Following ligands Failed to convert in Gypsum")
        print("Likely due to a Timeout")
        print(lig_failed_to_convert)
    sys.stdout.flush()
    return gypsum_output_folder_path


from autogrow.docking.execute_docking import pick_run_conversion_class_dict, pick_docking_class_dict, lig_convert_multithread
def pdb_to_pdbqt(vars, pdb_dir):
    ### adapted from run_docking_common
    dock_choice = vars["dock_choice"]
    conversion_choice = vars["conversion_choice"]
    receptor = vars["filename_of_receptor"]

    # Use a temp vars dict so you don't put mpi multiprocess info through itself...
    temp_vars = {}
    for key in list(vars.keys()):
        if key == "parallelizer":
            continue
        temp_vars[key] = vars[key]

    file_conversion_class_object = pick_run_conversion_class_dict(conversion_choice)
    file_conversion_class_object = file_conversion_class_object(temp_vars, receptor, test_boot=False)

    dock_class = pick_docking_class_dict(dock_choice)
    docking_object = dock_class(temp_vars, receptor, file_conversion_class_object, test_boot=False)

    if vars["docking_executable"] is None:
        docking_executable = docking_object.get_docking_executable_file(temp_vars)
        vars["docking_executable"] = docking_executable
    ##### vina or Qvina 

    # Find PDB's
    pdbs_in_folder = docking_object.find_pdb_ligands(pdb_dir)
    print('    pdb files:', pdbs_in_folder[:2], pdb_dir, len(pdbs_in_folder))
    job_input_convert_lig = tuple([(docking_object, pdb) for pdb in pdbs_in_folder])

    # print("    Convert Ligand from PDB to PDBQT format")
    smiles_names_failed_to_convert = vars["parallelizer"].run(job_input_convert_lig, lig_convert_multithread)

    pdbqts_in_folder = docking_object.find_converted_ligands(pdb_dir)
    print('    pdbqt file: ', len(pdbqts_in_folder), pdbqts_in_folder[:2])
    return docking_object

from autogrow.docking.execute_docking import run_dock_multithread, run_docking_common
import autogrow.docking.scoring.execute_scoring_mol as Scoring
import autogrow.docking.ranking.ranking_mol as Ranking


def docking_pdbqt(vars, docking_object, pdbqt_folder, full_smiles_file):
    pdbqts_in_folder = docking_object.find_converted_ligands(pdbqt_folder)
    job_input_dock_lig = tuple([tuple([docking_object, pdbqt]) for pdbqt in pdbqts_in_folder])
    smiles_names_failed_to_dock = vars["parallelizer"].run(job_input_dock_lig, run_dock_multithread)  
    ### main docking, (including delete failed docking file)

    deleted_smiles_names_list_dock = [x for x in smiles_names_failed_to_dock if x is not None]
    deleted_smiles_names_list_dock = list(set(deleted_smiles_names_list_dock))
    print("THE FOLLOWING LIGANDS WHICH FAILED TO DOCK:", deleted_smiles_names_list_dock)
    # print("#################### \n Begin Ranking and Saving results")
    # folder_with_pdbqts = current_generation_dir + "PDBs" + os.sep
    # Run any compatible Scoring Function
    print(full_smiles_file, pdbqt_folder)
    smiles_list = Scoring.run_scoring_common(vars, full_smiles_file, pdbqt_folder)
    print('---------', smiles_list[:3], 'smiles_list[:3] --------------')

    # Output format of the .smi file will be: SMILES    Full_lig_name
    # shorthandname   ...AnyCustominfo... Fitness_metric  diversity
    # Normally the docking score is the fitness metric but if we use a
    # Custom metric than dock score gets moved to index -3 and the new
    # fitness metric gets -2

    # sort list by the affinity of each sublist (which is the last index of sublist)
    smiles_list.sort(key=lambda x: float(x[-1]), reverse=False)
    # ["[N-]=[NH+]/N=C/c1[nH+]nc(-c2cccc3ccccc23)o1", "naphthalene_35", "naphthalene_35", "naphthalene_35__3", -9.2]
    # score the diversity of each ligand compared to the rest of the ligands in the group this adds on a float in the last column for the
    # sum of pairwise comparisons the lower the diversity score the more unique a molecule is from the other mols in the same generation
    smiles_list = Ranking.score_and_append_diversity_scores(smiles_list)
    # ["[N-]=[NH+]/N=C/c1[nH+]nc(-c2cccc3ccccc23)o1", "naphthalene_35", "naphthalene_35", "naphthalene_35__3", -9.2, 40.14 (diversity)]
    pdbqts_in_folder = [pdbqt + '.vina' for pdbqt in pdbqts_in_folder if os.path.exists(pdbqt + '.vina')]
    print('pdbqts [:4]', pdbqts_in_folder[:3], len(pdbqts_in_folder))

    id2pdbqt = defaultdict(lambda:[])
    for pdbqt in pdbqts_in_folder:
        smiles_id = pdbqt.split('/')[-1].split('__')[0]
        id2pdbqt[smiles_id].append(pdbqt)
    for idx,ss in enumerate(smiles_list):
        smiles_id = ss[1]
        smiles_list[idx].append(id2pdbqt[smiles_id])
    # ["[N-]=[NH+]/N=C/c1[nH+]nc(-c2cccc3ccccc23)o1", "naphthalene_35", "naphthalene_35", "naphthalene_35__3", 
    # '-9.2', 40.14 (diversity), ['results_xxxx/xxxx__1.pdbqt', 'results_xxxx/xxxxx__2.pdbqt']]
    return smiles_list

def docking(smiles_folder, smiles_file, args_dict):
    sdfs_folder_path = smiles_folder.strip('/') + '_SDF/'
    pdb_dir = smiles_folder.strip('/') + '_PDB/'
    smiles_to_sdfs(args_dict, gen_smiles_file=os.path.join(smiles_folder,smiles_file), smile_file_directory=smiles_folder)
    conversion_to_3d.convert_sdf_to_pdbs(args_dict, gen_folder_path=smiles_folder, sdfs_folder_path=sdfs_folder_path)
    docking_object = pdb_to_pdbqt(vars = args_dict, pdb_dir = pdb_dir)
    smiles_list = docking_pdbqt(args_dict, docking_object, pdb_dir, os.path.join(smiles_folder, smiles_file))
    return smiles_list 
##########################################################
# ['N=[N+]=[N+]=C(Cc1ccc2ccccc2c1)[N+](=O)[O-]', 'N=[N+]=Nc1c(N=[N+]=N)c(N=[N+]=N)c2ccccc2c1O', ...] 
############# receptor ############## 
receptor_info_list = [
    ('4r6e', './pdb/4r6e.pdb', -70.76, 21.82, 28.33, 15.0, 15.0, 15.0), 
    ('3pbl', './pdb/3pbl.pdb', 9, 22.5, 26, 15, 15, 15), 
    ('1iep', './pdb/1iep.pdb', 15.6138918, 53.38013513, 15.454837, 15, 15, 15), ]
    # ('2rgp', './pdb/2rgp.pdb', 16.29212, 34.870818, 92.0353, 15, 15, 15),
    # ('3eml', './pdb/3eml.pdb', -9.06363, -7.1446, 55.86259999, 15, 15, 15),
    # ('3ny8', './pdb/3ny8.pdb', 2.2488, 4.68495, 51.39820000000001, 15, 15, 15),
    # ('4rlu', './pdb/4rlu.pdb', -0.73599, 22.75547, -31.23689, 15, 15, 15),
    # ('4unn', './pdb/4unn.pdb', 5.684346153, 18.1917, -7.3715, 15, 15, 15),
    # ('5mo4', './pdb/5mo4.pdb', -44.901, 20.490354, 8.48335, 15, 15, 15),
    # ('7l11', './pdb/7l11.pdb', -21.81481, -4.21606, -27.98378, 15, 15, 15), ]


def update_receptor_info(vars, receptor_info):
    name_of_receptor, filename_of_receptor, center_x, center_y, center_z, size_x, size_y, size_z = receptor_info
    vars['name_of_receptor'] = name_of_receptor 
    vars['filename_of_receptor'] = filename_of_receptor 
    vars['center_x'] = center_x 
    vars['center_y'] = center_y 
    vars['center_z'] = center_z 
    vars['size_x'] = size_x 
    vars['size_y'] = size_y 
    vars['size_z'] = size_z  
    return vars 

smiles2info = defaultdict(lambda: dict())
id2smiles = dict() 
def random_generate_id(id2smiles):
    while True:
        smiles_id = str(random.randint(1000000, 9999999))
        if smiles_id not in id2smiles:
            return smiles_id 
##########################################################
################ initialize population ################
source_compound_file = args_dict['source_compound_file']
smiles_file = 'smiles.txt'
with open(source_compound_file, 'r') as fin:
    smiles_list = fin.readlines() 
initial_smiles_list = [smiles.split()[0] for smiles in smiles_list]
#######  docking initial smiles list
for receptor_info in receptor_info_list:
    vars = update_receptor_info(vars, receptor_info) 
    name_of_receptor = receptor_info[0]
    print("---------- 0.1. save smiles ----------")    ######    new_smiles_set -> smiles_file
    meta_result_folder = './results_' + name_of_receptor + '_'
    results_folder = meta_result_folder + "000"  ### 'results_4r6e_000' 
    new_smiles_list = initial_smiles_list[:30] #### debug 
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    full_smiles_file = os.path.join(results_folder, smiles_file)
    with open(full_smiles_file, 'w') as fout:
        for smiles in new_smiles_list:
            smiles_id = random_generate_id(id2smiles)
            fout.write(smiles + '\t' + smiles_id + '\n')
            id2smiles[smiles_id] = smiles 

    print("---------- 0.2. docking ----------")
    smiles_list = docking(smiles_folder = results_folder, smiles_file = smiles_file, args_dict = vars)
    # ["[N-]=[NH+]/N=C/c1[nH+]nc(-c2cccc3ccccc23)o1", "naphthalene_35", "naphthalene_35", "naphthalene_35__3", 
    # '-9.2', 40.14 (diversity), ['results_xxxx/xxxx__1.pdbqt', 'results_xxxx/xxxxx__2.pdbqt']]
    for info in smiles_list:
        smiles, smiles_id, binding_score, pdbqt_list = info[0], info[1], float(info[-3]), info[-1]
        smiles2info[name_of_receptor][smiles] = [smiles_id, binding_score, pdbqt_list]

    print('------ 0.3. top-K smiles for next generation -------')
    new_smiles_list = [(smiles, smiles2info[name_of_receptor][smiles][1], smiles2info[name_of_receptor][smiles][2]) \
                            for smiles in new_smiles_list if smiles in smiles2info[name_of_receptor]]
    new_smiles_list.sort(key=lambda x:x[1])
    new_smiles_list = new_smiles_list[:topk]
    smiles_info_list = [(smiles,pdbqt_list) for smiles,binding_score,pdbqt_list in new_smiles_list]
    smiles2info[name_of_receptor]['smiles_info_list'] = copy.deepcopy(smiles_info_list)
 
##########################################################
################# model ################
import torch 
from model import Ligand2D, Ligand2D_product, ENN, featurize_receptor_and_ligand
# featurize_receptor_and_ligand(pdbfile, centers, pocket_size, pdbqt_file,) 
# pdbtofeature(pdbfile, centers, pocket_size) & pdbqtvina2feature(pdbqt_file)
# crossover_policy_net_1 = Ligand2D() ##### TODO pocket & center 
# crossover_policy_net_2 = Ligand2D_product()
# crossover_policy_net_1 = ENN()
# crossover_policy_net_2 = ENN()
crossover_policy_net_1 = torch.load('save_model/crossover_policy_net_1.ckpt')
crossover_policy_net_2 = torch.load('save_model/crossover_policy_net_2.ckpt')
crossover_budget = 20 
crossover_optimizer = torch.optim.Adam(list(crossover_policy_net_1.parameters()) + list(crossover_policy_net_2.parameters()), lr=1e-3)


# mutation_policy_net_1 = ENN() 
# mutation_policy_net_2 = Ligand2D_product()
mutation_policy_net_1 = torch.load('save_model/mutation_policy_net_1.ckpt')
mutation_policy_net_2 = torch.load('save_model/mutation_policy_net_2.ckpt')
mutation_budget = 20
mutation_optimizer = torch.optim.Adam(list(mutation_policy_net_1.parameters()) + list(mutation_policy_net_2.parameters()), lr=1e-3)
################# model ################
crossover_train_data = defaultdict(lambda: dict())
mutation_train_data = defaultdict(lambda: dict())
canonicalize = lambda x:x 

crossover_done_set = defaultdict(lambda: set()) 
mutation_done_set = defaultdict(lambda: set())
########################## main loop ############################ 
for num_gen in tqdm(range(args_dict['num_generations'])):
    for receptor_info in receptor_info_list:
        ##### input: smiles_list (including pdbqtvina, from previous-generation) & receptor  
        vars = update_receptor_info(vars, receptor_info) 
        name_of_receptor = vars['name_of_receptor']
        smiles_info_list = copy.deepcopy(smiles2info[name_of_receptor]['smiles_info_list']) ###### [(smiles_1, pdbqtvina_list_1), (smiles_2, pdbqtvina_list_2), ...]
        print('===== 1. beginning of the generation: smiles_info_list =====', len(smiles_info_list), smiles_info_list[:5], )
        for info in smiles_info_list:
            pdbqtvina = info[1][0]
            assert os.path.exists(pdbqtvina)
        smiles_list = [smiles for smiles, pdbqtvina_list in smiles_info_list]
        smiles2pdbqtvina_local = {smiles:pdbqtvina_list[0] for smiles,pdbqtvina_list in smiles_info_list} 
        if len(smiles_list) <=3:
            continue 
        new_smiles_set = set()
        print("---------- 2. RGA: crossover ----------")
        pdbqtvina_list = [smiles2pdbqtvina_local[smiles] for smiles in smiles_list]
        print('length of ligand:', len(pdbqtvina_list))
        if num_gen > -1:
            ##### evaluate probability distribution in RGA 
            _, crossover_sample_probability_list = crossover_policy_net_1.forward_ligand_list(
                                                    name_of_receptor = vars['name_of_receptor'], 
                                                    pdbqtvina_list = pdbqtvina_list)

        else:
            crossover_sample_probability_list = [1.0/len(pdbqtvina_list) for i in pdbqtvina_list]
        sampled_idx = random.choices(list(range(len(crossover_sample_probability_list))), 
                                     weights = crossover_sample_probability_list, 
                                     k = crossover_budget)
        ####### RGA outer loop, first ligand ########
        for idx in tqdm(sampled_idx):
            selected_smiles_1 = smiles_list[idx]
            mol = execute_crossover.convert_mol_from_smiles(selected_smiles_1) 
            holdout_smiles_list = [smiles for smiles in smiles_list if smiles!=selected_smiles_1]
            if holdout_smiles_list == []:
                continue 
            holdout_pdbqtvina_list = [smiles2pdbqtvina_local[smiles] for smiles in holdout_smiles_list]
            if num_gen > -1:
                ##### evaluate probability distribution in RGA 
                _, crossover_sample_probability_list_2 = crossover_policy_net_2.forward_ligand_list(
                                                            name_of_receptor = vars['name_of_receptor'], 
                                                            pdbqtvina_list = holdout_pdbqtvina_list)
            else:
                crossover_sample_probability_list_2 = [1.0/len(holdout_pdbqtvina_list) for i in holdout_pdbqtvina_list]
            sampled_idx_2 = random.choices(list(range(len(holdout_pdbqtvina_list))), 
                                           weights = crossover_sample_probability_list_2, k = 10)
            ########## RGA inner loop, second ligand #########
            for idx2 in sampled_idx_2:
                smiles_2 = holdout_smiles_list[idx2]
                if (selected_smiles_1, smiles_2) not in crossover_done_set[name_of_receptor]:
                    crossover_done_set[name_of_receptor].add((selected_smiles_1, smiles_2))
                    crossover_done_set[name_of_receptor].add((smiles_2, selected_smiles_1))
                else:
                    continue 
                mol2 = execute_crossover.convert_mol_from_smiles(smiles_2) 
                if execute_crossover.test_for_mcs(vars, mol, mol2) is None:
                    continue 
                for i in range(3):
                    ligand_new_smiles = smiles_merge.run_main_smiles_merge(vars, selected_smiles_1, smiles_2)
                    if ligand_new_smiles is not None:
                        break 
                if ligand_new_smiles is not None:
                    ligand_new_smiles = canonicalize(ligand_new_smiles)
                    pass_or_not = Filter.run_filter_on_just_smiles(ligand_new_smiles, vars["filter_object_dict"])  #### True, False
                    if pass_or_not:
                        new_smiles_set.add(ligand_new_smiles)
                        crossover_train_data[name_of_receptor][ligand_new_smiles] = [smiles_list, selected_smiles_1, pdbqtvina_list, 
                                                                                     holdout_smiles_list, smiles_2, holdout_pdbqtvina_list,]
                    else:
                        # print("        >>> not pass filter") 
                        pass 
                else:
                    # print('        >>> '+ selected_smiles_1 + ' ' + smiles_2 +' merge fail ')
                    pass 
        print("     >>>>>> number of smiles generated by crossover", len(new_smiles_set), ' >>>')
        print("---------- 3. RGA: mutation ----------")
        pdbqtvina_list = [smiles2pdbqtvina_local[smiles] for smiles in smiles_list]
        if num_gen > -1:
            ##### evaluate probability distribution for first action in mutation in RGA 
            _, mutation_sample_probability_list = mutation_policy_net_1.forward_ligand_list(
                                                            name_of_receptor = vars['name_of_receptor'], 
                                                            pdbqtvina_list = pdbqtvina_list)
        else:
            mutation_sample_probability_list = [1/len(pdbqtvina_list) for i in pdbqtvina_list]
        sampled_idx = random.choices(list(range(len(mutation_sample_probability_list))), 
                                     weights = mutation_sample_probability_list, 
                                     k = mutation_budget)
        for idx in sampled_idx:
            selected_smiles_1 = smiles_list[idx]
            if selected_smiles_1 in mutation_done_set[name_of_receptor]:
                continue 
            else:
                mutation_done_set[name_of_receptor].add(selected_smiles_1)
            new_mutation_smiles_list = []  
            a_smiles_click_chem_object = SmileClickClass.SmilesClickChem(rxn_library_variables=rxn_library_variables, 
    																 list_of_already_made_smiles=new_mutation_smiles_list, 
    																 filter_object_dict=vars["filter_object_dict"])
            s_list = a_smiles_click_chem_object.run_smiles_click2(selected_smiles_1) ### smiles list
            if s_list is None or len(s_list)<=1:
                continue 
            s_list = [canonicalize(s) for s in s_list] 
            ##### evaluate probability distribution for first action in mutation in RGA 
            _, mutation_sample_probability_list_2 = mutation_policy_net_2(selected_smiles_1, s_list)
            sampled_idx_2 = random.choices(list(range(len(mutation_sample_probability_list_2))), 
                                           weights = mutation_sample_probability_list_2, 
                                           k = 2)
            for idx2 in sampled_idx_2:
                smiles_2 = s_list[idx2]
                mutation_train_data[name_of_receptor][smiles_2] = [smiles_list, selected_smiles_1, pdbqtvina_list, s_list]  
                new_smiles_set.add(smiles_2)
        print(">>>>>> number of smiles generated by crossover and mutation", len(new_smiles_set))
        print("---------- 4. elite from previous generation ----------")
        print("---------- 5. save smiles ----------")    ######    new_smiles_set -> smiles_file
        meta_result_folder = './results_' + name_of_receptor + '_'
        results_folder = meta_result_folder + str(num_gen)  ### 'results_0', 'results_1', ...
        new_smiles_list = list(new_smiles_set)
        new_smiles_list = new_smiles_list[:30] #### debug 
        if not os.path.exists(results_folder):
            os.makedirs(results_folder)
        full_smiles_file = os.path.join(results_folder, smiles_file)
        with open(full_smiles_file, 'w') as fout:
            for smiles in new_smiles_list:
                smiles_id = random_generate_id(id2smiles)
                fout.write(smiles + '\t' + smiles_id + '\n')
                id2smiles[smiles_id] = smiles 
        print('new_smiles_list', len(new_smiles_list))
        print("---------- 6. docking ----------")
        smiles_list = docking(smiles_folder = results_folder, smiles_file = smiles_file, args_dict = vars) #### receptor is in "vars"  
        # ["[N-]=[NH+]/N=C/c1[nH+]nc(-c2cccc3ccccc23)o1", "naphthalene_35", "naphthalene_35", "naphthalene_35__3", '-9.2', 40.14 (diversity), ['results_xxxx/xxxx__1.pdbqt', 'results_xxxx/xxxxx__2.pdbqt']]
        for info in smiles_list:
            smiles, smiles_id, binding_score, pdbqtvina_list = info[0], info[1], float(info[-3]), info[-1]
            smiles2info[name_of_receptor][smiles] = [smiles_id, binding_score, pdbqtvina_list]
        ############# policy gradient ##############
        print('------ 7. RGA: crossover policy network optimization -------')
        for ligand_new_smiles, (smiles_list, selected_smiles_1, pdbqtvina_list, holdout_smiles_list, smiles_2, holdout_pdbqtvina_list,) \
             in tqdm(crossover_train_data[name_of_receptor].items()):
            if not (selected_smiles_1 in smiles2info[name_of_receptor] \
                and smiles_2 in smiles2info[name_of_receptor] \
                and ligand_new_smiles in smiles2info[name_of_receptor]):
                continue 
            ##### log_likelihood
            log_crossover_sample_probability, _ = crossover_policy_net_1.forward_ligand_list(
                                                    name_of_receptor = vars['name_of_receptor'], 
                                                    pdbqtvina_list = pdbqtvina_list)
            idx = smiles_list.index(selected_smiles_1)
            log_likelihood_1 = log_crossover_sample_probability[idx]
            log_crossover_sample_probability_2, _ = crossover_policy_net_2.forward_ligand_list(
                                                    name_of_receptor = vars['name_of_receptor'], 
                                                    pdbqtvina_list = holdout_pdbqtvina_list)
            idx2 = holdout_smiles_list.index(smiles_2)
            log_likelihood_2 = log_crossover_sample_probability_2[idx2]
            ####### reward 
            v1 = smiles2info[name_of_receptor][selected_smiles_1][1]
            v2 = smiles2info[name_of_receptor][smiles_2][1]
            v3 = smiles2info[name_of_receptor][ligand_new_smiles][1]
            reward = -v3 - max(-v1, -v2) #### max is better 
            ####### policy gradient 
            log_likelihood = - (log_likelihood_1 + log_likelihood_2) * reward 
            crossover_optimizer.zero_grad() 
            log_likelihood.backward() 
            crossover_optimizer.step() 

        print('------ 8. RGA: mutation policy network optimization -------')
        for smiles_2, (smiles_list, selected_smiles_1, pdbqtvina_list, s_list) in tqdm(mutation_train_data[name_of_receptor].items()): 
            if not (selected_smiles_1 in smiles2info[name_of_receptor] and smiles_2 in smiles2info[name_of_receptor]):
                continue
            ##### log_likelihood
            log_mutation_sample_probability, _ = mutation_policy_net_1.forward_ligand_list(
                                                    name_of_receptor = vars['name_of_receptor'], 
                                                    pdbqtvina_list = pdbqtvina_list)
            idx = smiles_list.index(selected_smiles_1)
            log_likelihood_1 = log_mutation_sample_probability[idx]
            log_mutation_sample_probability_2, _ = mutation_policy_net_2(selected_smiles_1, s_list)
            idx2 = s_list.index(smiles_2)
            log_likelihood_2 = log_mutation_sample_probability_2[idx2]
            ####### reward 
            v1 = smiles2info[name_of_receptor][selected_smiles_1][1]
            v2 = smiles2info[name_of_receptor][smiles_2][1]
            reward = -v2 - (-v1) #### max is better 
            ####### policy gradient 
            log_likelihood = - (log_likelihood_1 + log_likelihood_2) * reward 
            mutation_optimizer.zero_grad() 
            log_likelihood.backward() 
            mutation_optimizer.step() 

        print('------ 9. top-K smiles for next generation -------')
        new_smiles_list = [(smiles, smiles2info[name_of_receptor][smiles][1], smiles2info[name_of_receptor][smiles][2]) for smiles in new_smiles_list if smiles in smiles2info[name_of_receptor]]
        new_smiles_list.sort(key=lambda x:x[1])
        new_smiles_list = new_smiles_list[:topk]
        print('new_smiles_list', len(new_smiles_list))
        smiles_info_list = [(smiles,pdbqt_list) for smiles,binding_score,pdbqt_list in new_smiles_list]
        smiles2info[name_of_receptor]['smiles_info_list'] = copy.deepcopy(smiles_info_list)


    for name_of_receptor, smiles_info_list in smiles2info.items():
        print('-------- Evaluating ' + name_of_receptor + ' ----------')
        # print(smiles2info[name_of_receptor])
        results = [[k,v[0],v[1],v[2]] for k,v in smiles2info[name_of_receptor].items() if k!='smiles_info_list']
        results.sort(key=lambda x:x[2])
        with open('result_'+name_of_receptor + '.txt', 'w') as fo:
            for result in results:
                fo.write('\t'.join([str(i) for i in result]) + '\n')


torch.save(crossover_policy_net_1, 'save_model/crossover_policy_net_1.ckpt')
torch.save(crossover_policy_net_2, 'save_model/crossover_policy_net_2.ckpt')
torch.save(mutation_policy_net_1, 'save_model/mutation_policy_net_1.ckpt')
torch.save(mutation_policy_net_2, 'save_model/mutation_policy_net_2.ckpt')




'''

rm -rf results_* 



python RGA.py \
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
    --num_generations 10 \
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



























