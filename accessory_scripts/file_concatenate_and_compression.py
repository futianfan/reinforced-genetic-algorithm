"""
This script is used to decompress or recompress AutoGrow data.

If you use the reduce_files_sizes option AutoGrow will convert concatenate and compress
all files in the PDBs directory of each generation. This is useful when doing larger runs as
data transfer is faster and data storage is reduced when files are merged and compressed.
    -The concatenation script that is run in AutoGrow 4 can be found at:
            autogrow4/autogrow/docking/concatenate_files.py
This script will either:
    1) Return the files back to their original uncompressed and deconcatenated formatting
                or
    2) concatenate and then compress the files into a single file.


The formatting of the concatenation is:
    "\n##############################File_name: {}\n".format(os.path.basename(file_name_1))
    ... Content of the 1st file...
    "\n##############################$$END_FILE$$ {}".format(os.path.basename(file_name_1))
    "\n##############################File_name: {}\n".format(os.path.basename(file_name_2))
    ... Content of the 2nd file...
    "\n##############################$$END_FILE$$ {}".format(os.path.basename(file_name_2))

Example decompression:
    python autogrow4/accessory_scripts/file_concatenation_and_compression.py \
    --compress_or_decompress decompress \
    --input_folder_or_file PATH_TO_RUN/Run_0/generation_1/PDBs/compressed_PDBS.txt.gz
Example compression:
    python autogrow4/accessory_scripts/file_concatenation_and_compression.py \
    --compress_or_decompress compress \
    --input_folder_or_file PATH_TO_RUN/Run_0/generation_1/PDBs/

This concatenated file is tar.gz compressed.
"""
import __future__

import glob
import os
import gzip
import shutil
import argparse

import support_scripts.Multiprocess as mp


def compress_file(file_name):
    """
    Compress the concatenated file

    Inputs:
    :param str file_name: the path to the file to compress.
    """

    with open(file_name, "r") as f:
        printout = f.read()
    printout = printout.encode("utf-8")
    with gzip.open(file_name + ".gz", "wb") as f:
        f.write(printout)


#######
def decompress_file(decompressed_file):
    """
    Decompress a file. Not used in running the program but is the counter of
    def compress_file(file_name)

    Inputs:
    :param str decompressed_file: the path to the file to decompress.

    Returns:
    :returns: str decompressed_file: the path to the file to decompress.
    """

    out_file = decompressed_file.replace(".gz", "")
    with gzip.open(decompressed_file, "rb") as f_comp:
        with open(out_file, "wb") as f_decomp:
            shutil.copyfileobj(f_comp, f_decomp)
    return out_file


#######
def separate_files(compressed_file):
    """
    Separate a concatenated file. Not used in running the program but is the
    counter of def compress_file(file_name)

    Inputs:
    :param str compressed_file: the path to the file to separate/decompress.
    """

    directory = (
        os.path.abspath(compressed_file.split(os.path.basename(compressed_file))[0])
        + os.sep
    )
    compressed_file = os.path.abspath(compressed_file)

    decompressed_file = decompress_file(compressed_file)
    if os.path.exists(decompressed_file) is False:
        raise Exception("Failed to decompress the file")

    printout = ""
    list_of_new_files = []
    out_file = None
    with open(decompressed_file, "r") as f:
        for line in f.readlines():
            if "$$END_FILE$$" in line:
                if out_file is not None and os.path.exists(out_file) is False:
                    with open(out_file, "w") as f:
                        f.write(printout + "\n")
                out_file = None
                printout = ""
                continue
            if "File_name:" in line:

                printout = ""

                # Split the line up and grab the relative file path convert to
                # absolute path
                out_file = (
                    directory
                    + os.sep
                    + line.split("##############################File_name: ")[
                        1
                    ].replace("\n", "")
                )
                out_file = os.path.abspath(out_file)
                list_of_new_files.append(out_file)
                continue

            printout = printout + line
            continue

    all_are_made = True
    for f in list_of_new_files:
        if os.path.exists(f) is False:
            print("file failed to decompress: {}".format(f))
            all_are_made = False
    if all_are_made is True:
        to_run = "rm {}".format(decompressed_file)
        os.system(to_run)


#######
def get_file_info(file_name):
    """
    Used for concatenating files together. This function appends a seperator
    and the filename of a file before and after the text of the file
    file_name. It returns it as a string

    Inputs:
    :param str file_name: the path to the file to compress.

    Returns:
    :returns: str concat: the text of the file file_name with a seperator and
        label before and after the file text.
    """
    file_name_insert = "\n##############################File_name: {}\n".format(
        os.path.basename(file_name)
    )
    file_termination_insert = "\n##############################$$END_FILE$$ {}".format(
        os.path.basename(file_name)
    )
    concat = file_name_insert + open(file_name).read() + file_termination_insert
    return concat


#######
def del_files(file_name):
    """
    This function deletes a given file file_name.

    Inputs:
    :param str file_name: the path to delete.
    """

    if os.path.exists(file_name):
        try:
            os.system("rm {}".format(file_name))
        except:
            print("couldn't delete file: {}".format(file_name))


#######
def run_concatenation(directory):
    """
    This function concatenates and compresses every file in a directory. This
    makes data transfer easier later on.

    To decompress the folder please use script in
    $PATH/autogrow4/accessory_scripts/file_concatenation_and_compression.py

    Inputs:
    :param str directory: the path to the folder which will be compiled and compressed.
    """

    concat_file = directory + os.sep + "compressed_PDBS.txt"
    print(
        "Start Concatenation: To separate files use the \
        file_concatenation_and_compression.py in the Utility script folder."
    )
    file_list = glob.glob(directory + os.sep + "*")
    file_list = [os.path.abspath(x) for x in file_list]

    with open(concat_file, "a+") as f:
        for file_name in file_list:
            f.write(get_file_info(file_name))

    job_list = tuple([(file_path,) for file_path in file_list])
    print("\tFinish Concatenation")
    print("\tRemoving files that were concatenated")
    mp.multi_threading(job_list, -1, del_files)
    print("\tCompressing file")
    compress_file(concat_file)
    if os.path.exists(concat_file + ".gz"):
        del_files(concat_file)
    print("Finished Compression")
######
def run_main(vars):
    """
    This function runs the functions for compression or decompression.
    Inputs:
    :param dict vars: dictionary of user variables.
    """

    if vars["compress_or_decompress"] == "compress":

        print("BEFORE")
        input_folder = vars["input_folder_or_file"]
        print(os.path.getsize(input_folder))

        run_concatenation(input_folder)
        print("FINISH CONCATENATE")
        print("After concatenate")
        print(os.path.getsize(input_folder))

    elif vars["compress_or_decompress"] == "decompress":
        compressed_file = vars["input_folder_or_file"]

        if os.path.exists(compressed_file) is False:
            raise Exception("File to Decompress doesn't exist")

        input_folder = os.path.abspath(compressed_file.split(os.path.basename(compressed_file))[0]) + os.sep

        print("BEFORE")
        print(os.path.getsize(input_folder))

        separate_files(compressed_file)
        print("After deconcatenate")
        print(os.path.getsize(input_folder))

        del_files(compressed_file)
        print("After deconcatenate")
        print(os.path.getsize(input_folder))

#######



def get_arguments_from_argparse(arg_dict):
    """
    This function handles the arg parser arguments for the script.

    Inputs:
    :param dict arg_dict: dictionary of parameters
    Returns:
    :returns: dict arg_dict: dictionary of parameters
    """

    # Argument handling
    if type(arg_dict["compress_or_decompress"]) != str:
        raise Exception("Must chose between compress or decompress")
    if arg_dict["compress_or_decompress"].lower() not in ["compress", "decompress"]:
        raise Exception("Must chose between compress or decompress")

    # set to lower case to prevent issues
    arg_dict["compress_or_decompress"] = arg_dict["compress_or_decompress"].lower()

    #  argument_handling
    if type(arg_dict["input_folder_or_file"]) is not str:
        raise Exception("--input_folder_or_file required: Path to directory to \
                        compress or decompress.")
    if os.path.exists(arg_dict["input_folder_or_file"]) is False:
        raise Exception("--input_folder_or_file could not be found: \
                        {}".format(arg_dict["input_folder_or_file"]))

    # Make sure variable is full path and add os.sep if directory
    arg_dict["input_folder_or_file"] = os.path.abspath(arg_dict["input_folder_or_file"])
    if os.path.isdir(arg_dict["input_folder_or_file"]):
        arg_dict["input_folder_or_file"] = arg_dict["input_folder_or_file"] + os.sep

    return arg_dict
#

# Argument parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument('--compress_or_decompress', type=str, required=True,
                    choices=["compress", "decompress"],
                    help='Chose whether to compress or decompress a directory')
PARSER.add_argument('--input_folder_or_file', '-i', type=str,
                    required=True, default=None,
                    help='Path to directory/file to compress or decompress.')

ARGS_DICT = vars(PARSER.parse_args())
ARGS_DICT = get_arguments_from_argparse(ARGS_DICT)
run_main(ARGS_DICT)
