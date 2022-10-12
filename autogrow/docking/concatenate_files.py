"""
This script compresses files which makes it easier to transfer data To
    decompress the files use the script in
    $PATH/autogrow4/accessory_scripts/file_concatenation_and_compression.py .
"""
import __future__

import glob
import os
import gzip
import shutil


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
    separate a concatenated file. Not used in running the program but is the
    counter of def compress_file(file_name)

    Inputs:
    :param str compressed_file: the path to the file to separate/decompress.
    """

    directory = (
        os.path.abspath(compressed_file.split(os.path.basename(compressed_file))[0])
        + os.sep
    )
    compressed_file = os.path.abspath(compressed_file)

    decompressed_file = decompress_file(directory, compressed_file)
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
            elif "File_name:" in line:

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
            else:
                printout = printout + line
                continue

    all_are_made = True
    for f in list_of_new_files:
        if os.path.exists(f) is False:
            print("file failed to decompress: {}".format(f))
            all_are_made = False
    if all_are_made is True:
        torun = "rm {}".format(decompressed_file)
        os.system(torun)


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
def run_concatenation(parallelizer_object, directory):
    """
    This function concatenates and compresses every file in a directory. This
    makes data transfer easier later on.

    To decompress the folder please use script in
    $PATH/autogrow4/file_concatenation_and_compression.py

    parallelizer_object type is <class
    'autogrow.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer.Parallelizer'>

    Inputs:
    :param Parallelizer_obj parallelizer_object: a paralellizer object used to
        multiprocess. initialized from
        autogrow/operators/convert_files/gypsum_dl/gypsum_dl/Parallelizer.py.
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
    parallelizer_object.run(job_list, del_files)
    print("\tCompressing file")
    compress_file(concat_file)
    if os.path.exists(concat_file + ".gz"):
        del_files(concat_file)
    print("Finished Compression")
