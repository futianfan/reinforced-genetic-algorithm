"""
This script holds the parent class for docking.
This is used as the basis for all docking classes.
"""
import __future__

class ParentDocking(object):
    """
    Docking

    Inputs:
    :param class object: a class to initialize on
    """

    def __init__(self, vars=None, receptor_file=None,
                 file_conversion_class_object=None, test_boot=True):
        """
        Require to initialize any pdbqt conversion class.

        Inputs:
        :param dict vars: Dictionary of User variables
        :param str receptor_file: the path for the receptor pdb
        :param obj file_conversion_class_object: object which is used to
            convert files from pdb to pdbqt
        :param bool test_boot: used to initialize class without objects for
            testing purpose
        """

        pass

    def get_name(self):
        """
        Returns the current class name.

        Returns:
        :returns: str self.__class__.__name__: the current class name.
        """
        return self.__class__.__name__

    def run_dock(self, pdbqt_filename):
        """
        run_dock is needs to be implemented in each class.

        Inputs:
        :param str pdbqt_filename: a string for docking process raise exception
            if missing
        """

        raise NotImplementedError("run_dock() not implemented")

    def rank_and_save_output_smi(self, vars, current_generation_dir,
                                 current_gen_int, smile_file,
                                 deleted_smiles_names_list):
        """
        rank_and_save_output_smi is needs to be implemented in each class.
        raise exception if missing

        Given a folder with PDBQT's, rank all the SMILES based on docking
        score (High to low). Then format it into a .smi file. Then save the
        file.

        Inputs:
        :param dict vars: vars needs to be threaded here because it has the
            paralizer object which is needed within Scoring.run_scoring_common
        :param str current_generation_dir: path of directory of current
            generation
        :param int current_gen_int: the interger of the current generation
            indexed to zero
        :param str smile_file: File path for the file with the ligands for the
            generation which will be a .smi file
        :param list deleted_smiles_names_list: list of SMILES which may have
            failed the conversion process

        Returns:
        :returns: str output_ranked_smile_file: the path of the output ranked
            .smi file
        """

        raise NotImplementedError("rank_and_save_output_smi() not implemented")
