"""Mozziconacci Filter
This runs a Mozziconacci filter. Mozziconacci filter is a filter for
Drug-likeliness which filters molecules by the number of: rotatable bonds,
rings, oxygens, and halogens.

To pass the filter a molecule must be:
    # of Rotatable bonds: Max 15
    # of Rings: Max 6
    # of Oxygens: Min 1
    # of Nitrogens: Min 1
    # of Halogens: Max 7

If you use the Mozziconacci Filter please cite: Mozziconacci, J. C. et al.
Preparation of a Molecular Database from a Set of 2 Million Compounds for
Virtual Screening Applications: Gathering, Structural Analysis and Filtering.
9th Electronic Computational Chemistry Conference, World Wide Web, March
(2003).
"""

import __future__

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.Lipinski as Lipinski

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

from autogrow.operators.filter.filter_classes.parent_filter_class import ParentFilter


class MozziconacciFilter(ParentFilter):
    """
    This runs a Mozziconacci filter. Mozziconacci filter is a filter for
    Drug-likeliness which filters molecules by the number of:

    To pass the filter a molecule must be:
        # of Rotatable bonds: Max 15
        # of Rings: Max 6
        # of Oxygens: Min 1
        # of Nitrogens: Min 1
        # of Halogens: Max 7

    If you use the Mozziconacci Filter please cite: Mozziconacci, J. C. et al.
    Preparation of a Molecular Database from a Set of 2 Million Compounds for
    Virtual Screening Applications: Gathering, Structural Analysis and
    Filtering. 9th Electronic Computational Chemistry Conference, World Wide
    Web, March (2003).

    Inputs:
    :param class ParentFilter: a parent class to initialize off
    """

    def run_filter(self, mol):
        """
        This runs a Mozziconacci filter. Mozziconacci filter is a filter for
        Drug-likeliness which filters molecules by the number of:

        To pass the filter a molecule must be:
            # of Rotatable bonds: Max 15
            # of Rings: Max 6
            # of Oxygens: Min 1
            # of Nitrogens: Min 1
            # of Halogens: Max 7

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be
            tested if it passes the filters

        Returns:
        :returns: bool bool: True if the mol passes the filter; False if it
            fails the filter
        """

        halogen = Chem.MolFromSmarts("[*;#9,#17,#35,#53,#85]")
        number_of_halogens = len(mol.GetSubstructMatches(halogen, maxMatches=8))
        if number_of_halogens > 7:
            return False

        oxygen = Chem.MolFromSmarts("[#8]")
        number_of_oxygens = len(mol.GetSubstructMatches(oxygen, maxMatches=2))
        if number_of_oxygens < 1:
            return False

        nitrogen = Chem.MolFromSmarts("[#7]")
        number_of_nitrogen = len(mol.GetSubstructMatches(nitrogen, maxMatches=2))
        if number_of_nitrogen < 1:
            return False

        num_rotatable_bonds = Lipinski.NumRotatableBonds(mol)
        if num_rotatable_bonds > 15:
            return False

        ring_count = Chem.rdmolops.GetSSSR(mol)
        if ring_count > 6:
            return False

        # Passes everything
        return True
