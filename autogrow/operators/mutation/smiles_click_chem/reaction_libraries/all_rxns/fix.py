"""
Developers script for merging robust_rxns and click_chem_rxns
into all_rxns... This script renumbers a rxn_libraries numbering
"""

import sys




def renumber_file(old_path, new_path, new_rxn_num):
    """
    This is a developers tool to renumber a rxn_library .json file
        Must provide the following:
        1) old_path to original .json rxn_library file
        2) new_path to output renumbered .json rxn_library file
        3) new_rxn_num: number to reindex the 1st reaction to:
            ie if new_rxn_num=37 then rxn_num=1 becomes rxn_num=37
    Inputs:
    :param str old_path: path to original .json rxn_library
        file to be renumbered
    :param str new_path: path to output renumbered .json rxn_library
    :param int new_rxn_num: number to reindex the 1st reaction to:
            ie if new_rxn_num=37 then rxn_num=1 becomes rxn_num=37
    """
    printout = ""


    original_rxn_num = 1 # index of the 1st reaction to adjust to rxn_num
    with open(old_path, "r") as f:
        for line in f.readlines():
            if '": {' in line or "reaction_name" in line:
                line = line.replace("{}_".format(
                    original_rxn_num), "{}_".format(new_rxn_num))
            elif '"RXN_NUM": ' in line:
                line = line.split(": ")[0] + ": {}\n".format(new_rxn_num)
                new_rxn_num = new_rxn_num + 1
                original_rxn_num = original_rxn_num + 1

            printout = printout + line

    with open(new_path, "w") as f:
        f.write(printout)

if __name__ == "__main__":
    # OLD_PATH ="$PATH/Input_rxn_library/temp.json"
    # new_path="$PATH/output_renumbered_rxn_library/temp_new.json"
    # NEW_RXN_NUM = 37
    try:
        OLD_PATH = sys.argv[1]
        NEW_PATH = sys.argv[2]
        NEW_RXN_NUM = sys.argv[3]
    except:
        print("This is a developers tool to renumber a rxn_library .json file \
                Must provide the following: \
                1) OLD_PATH to original .json rxn_library file \
                2) NEW_PATH to output renumbered .json rxn_library file \
                3) NEW_RXN_NUM: number to reindex the 1st reaction to: \
                    ie if NEW_RXN_NUM=37 then rxn_num=1 becomes rxn_num=37")

    renumber_file(OLD_PATH, NEW_PATH, NEW_RXN_NUM)
#
