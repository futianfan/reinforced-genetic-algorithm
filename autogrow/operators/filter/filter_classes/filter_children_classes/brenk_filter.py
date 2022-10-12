"""#BRENK filter
This will filter a ligand using the BRENK filter for lead-likeliness, by
matching common false positive molecules to the current mol..

This script relies on the RDKit predefined FilterCatalog. FilterCatalog is
maintained by RDKit.

If using the BRENK filter please cite: Brenk R et al. Lessons Learnt from
Assembling Screening Libraries for Drug Discovery for Neglected Diseases.
ChemMedChem 3 (2008) 435-444. doi:10.1002/cmdc.200700139.
"""

import __future__

from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from autogrow.operators.filter.filter_classes.parent_filter_class import ParentFilter


class BRENKFilter(ParentFilter):
    """
    This will filter a ligand using a BRENK screening filter for
    lead-likeliness, by matching common false positive molecules to the
    current mol.

    This script relies on the RDKit predefined FilterCatalog. FilterCatalog is
        maintained by RDKit.

    If using the BRENK filter please cite: Brenk R et al. Lessons Learnt from
    Assembling Screening Libraries for Drug Discovery for Neglected Diseases.
    ChemMedChem 3 (2008) 435-444. doi:10.1002/cmdc.200700139.

    Inputs:
    :param class ParentFilter: a parent class to initialize off of.
    """

    def __init__(self):
        """
        This loads in the filters which will be used.
        """

        self.filters = self.get_filters()

    def get_filters(self):
        """
        This loads in the filters which will be used.

        Returns:
        :returns: rdkit.Chem.rdfiltercatalog.FilterCatalog filters: A set of
            RDKit Filters
        """

        # Make a list of the BRENK filter.
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        # This is our set of all the BRENK filters
        filters = FilterCatalog.FilterCatalog(params)
        return filters

    def run_filter(self, mol):
        """
        Runs a BRENK filter by matching common false positive molecules to the
        current mol. Filters for for lead-likeliness.

        Based on the PAINS filter implementation in RDKit described in
        http://rdkit.blogspot.com/2016/04/changes-in-201603-release-filtercatalog.html

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be
            tested if it passes the filters

        Returns:
        :returns: bool bool: True if the mol passes the filter; False if it
            fails the filter
        """

        # If the mol matches a mol in the filter list. we return a False (as
        # it failed the filter).
        if self.filters.HasMatch(mol) is True:
            return False

        # if No matches are found to filter list this will return a True
        # as it Passed the filter.
        return True
