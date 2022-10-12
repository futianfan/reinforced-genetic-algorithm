"""
convert pdbs into smiles

This script will take a folder and convert all pdb files into a single texted file.
The text file will contain smiles strings of the respective pdb and the name of the file.

Run example:

output example:
python convert_ligands_pdb_to_smi.py \
  --source_folder $PATH/OF/PDBS/ \
  --output_folder $PATH/TO/OUTPUT/ \
  --number_of_processors -1

This will convert all .pdb files within $PATH/OF/PDBS/ into a single
.smi file at $PATH/TO/OUTPUT/PDBS.smi
using all available processors

CC1COC(=O)OC(=O)O1    ZINC60039447
O = C1OC(=O)N2CCC12    ZINC59199492
O = C1CC(C(=O)O)C(=O)O1    ZINC59901386
"""
import __future__

import glob
import os
import argparse

from rdkit import Chem

import support_scripts.mol_object_handling as MOH
import support_scripts.Multiprocess as mp

def run_convert_on_single_pdb(pdb):

    """
    This function converts a ligand into SMILES
    and returns the list with the smiles with a name.
    The names are the basename of the file minus the .pdb

    Inputs:
    :param str pdb: path to the folder to a pdb file
    Returns:
    :returns: list output_data: A list containing all SMILES info from the file
    """

    try:
        mol = Chem.MolFromPDBFile(pdb)

        mol_sanitized = MOH.check_sanitization(mol)
        if mol_sanitized is not None:
            smiles = Chem.MolToSmiles(mol_sanitized, isomericSmiles=True)
            file_name = os.path.basename(pdb)
            file_stripped = file_name.replace(".pdb", "")
            output_data = smiles + "\t" + file_stripped
    except:
        pass
    return output_data
#

def make_smile_list(sub_folder):
    """
    This function converts every ligand within a folder into SMILES
    and returns the list of smiles with a name.
    The names are the basename of the file minus the .pdb

    Inputs:
    :param str sub_folder: path to the folder to search for pdb files
    Returns:
    :returns: list smiles_list: A list of lists containing all SMILES from
        the .pdb files and their respective name
    """
    sub_folder = sub_folder+os.sep
    smiles_list = []
    pdb_list = glob.glob(os.sep + sub_folder+"*.pdb")
    pdb_list.extend(glob.glob(os.sep + sub_folder+"*.PDB"))
    pdb_list = tuple([tuple([pdb]) for pdb in pdb_list])

    # run convert in multithread
    smiles_list = mp.multi_threading(pdb_list, -1, run_convert_on_single_pdb)

    return smiles_list
#
def start_run_main(vars):
    """
    This will run the main arguments for the script.

    Inputs:
    :param dict vars: dictionary of user variables.
    """
    # Running converter
    smiles_list = make_smile_list(vars["source_folder"])
    name = [x for x in vars["source_folder"].split(os.sep)if x != ""][-1]
    output_file = vars["output_folder"] + os.sep + name + ".smi"
    with open(output_file, "w") as f:
        f.write("\n".join(smiles_list))

    print("Converted ligands to .smi located:\n\t{}".format(output_file))

def get_arguments_from_argparse(args_dict):
    """
    This function handles the arg parser arguments for the script.

    Inputs:
    :param dict args_dict: dictionary of parameters
    Returns:
    :returns: dict args_dict: dictionary of parameters
    """


    # Argument handling
    if type(args_dict["source_folder"]) != str:
        raise Exception("provided source folder must be a directory.")
    if type(args_dict["output_folder"]) != str:
        raise Exception("provided output_folder must be a directory.")

    #  argument_handling
    if os.path.exists(args_dict["source_folder"]) is False or \
        os.path.isdir(args_dict["source_folder"]) is False:
        raise Exception("provided source folder can not be found or is not a directory.")
    args_dict["source_folder"] = os.path.abspath(args_dict["source_folder"]) + os.sep

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

    return args_dict
#

# Argument parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument(
    '--source_folder', '-s', required=True, default=None,
    help='Path to folder containing .pdb files to convert. \
    File must contain a single small molecules. Without protein. \
    Files must end with either .pdb or .PDB'
)
PARSER.add_argument(
    '--output_folder', '-o', required=True, default=None,
    help='Path to folder where we will output a .smi file of converted .pdb files.'
)
# processors and multithread mode
PARSER.add_argument(
    '--number_of_processors', '-p', type=int, metavar='N', default=1,
    help='Number of processors to use for parallel calculations. \
    This script is not MPI enable but is able to multithread using SMP architecture. \
    Set to -1 for all available CPUs.'
)

ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)

# Running converter
start_run_main(ARGS_DICT)
