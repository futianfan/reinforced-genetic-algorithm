"""
This script holds the parent class for file conversion for docking.
This is used as the basis for all file conversion classes.
"""
import __future__

class ParentPDBQTConverter(object):
    """
    Docking

    Inputs:
    :param class object: a class to initialize on
    """

    def __init__(self, vars=None, receptor_file=None, test_boot=True):
        """
        Require to initialize any pdbqt conversion class.

        Inputs:
        :param dict vars: Dictionary of User variables
        :param str receptor_file: the path for the receptor pdb
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

    def convert_receptor_pdb_files_to_pdbqt(self, receptor_file, mgl_python,
                                            receptor_template,
                                            number_of_processors):
        """
        Make sure a PDB file is properly formatted for conversion to pdbqt

        Inputs:
        :param str receptor_file:  the file path of the receptor
        :param str mgl_python: file path of the pythonsh file of mgl tools
        :param str receptor_template: the receptor4.py file path from mgl
            tools.
        :param int number_of_processors: number of processors to multithread
        """

        raise NotImplementedError(
            "convert_receptor_pdb_files_to_pdbqt() not implemented"
        )

    def convert_ligand_pdb_file_to_pdbqt(self, pdb_file):
        """
        Convert the ligands of a given directory from pdb to pdbqt format

        Inputs:
        :param str pdb_file: the file name, a string.
        Returns:
        :returns: bool bool: True if it worked; False if its the gypsum param
            file or if it failed to make PDBQT
        :returns: str smile_name: name of the SMILES string from a pdb file
            None if its the param file
        """

        raise NotImplementedError("rank_and_save_output_smi() not implemented")
