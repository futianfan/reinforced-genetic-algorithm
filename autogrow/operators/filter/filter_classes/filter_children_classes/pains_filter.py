"""#PAINS Filter
This will filter a ligand using a PAINS filter. PAINS eliminates of Pan Assay
Interference Compounds using substructure a search.

This will include PAINS_A, PAINS_B, and PAINS_C filtering.

This script relies on the RDKit predefined FilterCatalog. FilterCatalog is
    maintained by RDKit.

If using the PAINS filter please cite: Baell JB, Holloway GA. New Substructure
Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening
Libraries and for Their Exclusion in Bioassays. J Med Chem 53 (2010) 2719D40.
doi:10.1021/jm901137j.
"""

import __future__

from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from autogrow.operators.filter.filter_classes.parent_filter_class import ParentFilter


class PAINSFilter(ParentFilter):
    """
    This will filter a ligand using a PAINS filter. PAINS eliminates of Pan
    Assay Interference Compounds using substructure a search.

    This will include PAINS_A, PAINS_B, and PAINS_C filtering.

    This script relies on the RDKit predefined FilterCatalog. FilterCatalog is
        maintained by RDKit.

    If using the PAINS filter please cite: Baell JB, Holloway GA. New
    Substructure Filters for Removal of Pan Assay Interference Compounds
    (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. J
    Med Chem 53 (2010) 2719D40. doi:10.1021/jm901137j.

    Inputs:
    :param class ParentFilter: a parent class to initialize off
    """

    def __init__(self):
        """
        This loads in the filters which will be used.
        """

        self.filters_list = self.get_filters_list()

    def get_filters_list(self):
        """
        This loads in the filters which will be used.

        Returns:
        :returns: rdkit.Chem.rdfiltercatalog.FilterCatalog filters: A set of
            RDKit Filters
        """

        # Make a list of all the different PAINS Filters. PAINS should include
        # PAINS_A,PAINS_B, and PAINS_C, but because RDKit documentation
        # doesn't specify this explicitly we have included all 4 of the PAINS
        # FilterCatalogs for precaution.
        params_PAINS_A = FilterCatalogParams()
        params_PAINS_A.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params_PAINS_B = FilterCatalogParams()
        params_PAINS_B.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params_PAINS_C = FilterCatalogParams()
        params_PAINS_C.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        params_PAINS = FilterCatalogParams()
        params_PAINS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)

        params_list = [params_PAINS_A, params_PAINS_B, params_PAINS_C, params_PAINS]
        filters_list = []
        for param in params_list:
            filter = FilterCatalog.FilterCatalog(param)
            filters_list.append(filter)

        return filters_list

    def run_filter(self, mol):
        """
        Runs a PAINS filter by matching common false positive molecules to the
        current mol.

        This will filter a ligand using a PAINS filter. PAINS eliminates of
        Pan Assay Interference Compounds using substructure a search.

        Based on the PAINS filter implementation in RDKit described in
        http://rdkit.blogspot.com/2016/04/changes-in-201603-release-filtercatalog.html

        Inputs:
        :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be
            tested if it passes the filters Returns:
        Returns:
        :returns: bool bool: True if the mol passes the filter;
            False if it fails the filter
        """

        # This is our set of all the PAINS filters
        for filters in self.filters_list:

            # If the mol matches a mol in the filter list we return a False
            # (as it failed the filter)
            if filters.HasMatch(mol) is True:
                return False

        # if No matches are found to filter list this will return a True as it
        # Passed the filter.
        return True
