"""
delete_failed_mol.py deletes all files associated with a ligand
This is done to clean up failed attempts to dock.
"""
#######################
# Delete extra ligands #
########################

import glob
import os


def delete_all_associated_files(pdb_filename):
    """Delete files associated with a compound

    Inputs:
    :param str pdb_filename: the filename of the compounds.
    """

    toremove = [pdb_filename]
    toremove.extend(glob.glob(pdb_filename[:-3] + "*"))
    toremove.extend(
        glob.glob(
            os.path.dirname(pdb_filename)
            + os.sep
            + "support"
            + os.sep
            + os.path.basename(pdb_filename)[:-3]
            + "*"
        )
    )

    # Remove any redundancy
    toremove = list(set(toremove))

    # print("DELETING FOLLOWING!:", toremove)
    for todel in toremove:
        if os.path.exists(todel):
            os.remove(todel)


#
