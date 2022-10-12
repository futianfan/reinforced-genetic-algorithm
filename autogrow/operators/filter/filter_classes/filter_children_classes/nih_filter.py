"""#NIH Filter
This will filter a ligand using the NIH filter to remove ligands with
undersirable functional groups.

This script relies on the RDKit predefined FilterCatalog. FilterCatalog is
maintained by RDKit.

If using the NIH filter please cite: Jadhav A, et al. Quantitative Analyses of
Aggregation, Autofluorescence, and Reactivity Artifacts in a Screen for
Inhibitors of a Thiol Protease. J Med Chem 53 (2009) 37D51.
doi:10.1021/jm901070c.
"""
import __future__

import rdkit
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from autogrow.operators.filter.filter_classes.parent_filter_class import ParentFilter


class NIHFilter(ParentFilter):
    """
    This will filter a ligand using a NIH screening filter. This script relies
    on the RDKit predefined FilterCatalog. FilterCatalog is maintained by
    RDKit.

    The NIH filter is used to eliminate ligands with undersirable functional
    groups

    If using the NIH filter please cite: Jadhav A, et al. Quantitative
    Analyses of Aggregation, Autofluorescence, and Reactivity Artifacts in a
    Screen for Inhibitors of a Thiol Protease. J Med Chem 53 (2009) 37D51.
    doi:10.1021/jm901070c.

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

        # Make a list of the NIH filter.
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        # This is our set of all the NIH filters
        filters = FilterCatalog.FilterCatalog(params)
        return filters

    def run_filter(self, mol):
        """
        Runs a NIH filter by matching common false positive molecules to the
        current mol.

        The NIH filter is used to eliminate ligands with undersirable
        functional groups

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
        # it failed the filter)
        if self.filters.HasMatch(mol) is True:
            return False

        # if No matches are found to filter list this will return a True
        # as it Passed the filter.
        return True
