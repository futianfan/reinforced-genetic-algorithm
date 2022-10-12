"""
This script holds the parent class for filtering.
This is used as the basis for all filter classes.
"""
import __future__
class ParentFilter(object):
    """
    This is a script containing all of the filters for drug likeliness

    Filters for orally bio-available drugs:
        1) Lipinski

    Filters for for lead-likeness:
        1) GhoseFilter
        2) GhoseModifiedFilter
        3) MozziconacciFilter

    Filters for CNS/Blood Brain Barrier Permeable:
        1) VandeWaterbeemdFilter

    False-Positive/Metabolite substructure searches:
        1) PAINSFilter
        2) NIHFilter
        3) BRENKFilter
    """
    def get_name(self):
        """
        Returns the current class name.

        Returns:
        :returns: str self.__class__.__name__: the current class name.
        """

        return self.__class__.__name__

    def run_filter(self, input_string):
        """
        run_filter is needs to be implemented in each class.

        Inputs:
        :param str input_string:  A string to raise an exception
        """

        raise NotImplementedError("run_filter() not implemented")
