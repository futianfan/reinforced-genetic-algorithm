"""
This script holds the parent class for scoring/rescoring.
This is used as the basis for all scoring/rescoring classes.
"""
import __future__

class ParentScoring(object):
    """
    This is a script containing all of the scoring functions.
    """

    def get_name(self):
        """
        Returns the current class name.

        Returns:
        :returns: str self.__class__.__name__: the current class name.
        """
        return self.__class__.__name__

    def run_scoring(self, input_string):
        """
        run_scoring is needs to be implemented in each class.

        Inputs:
        :param str input_string:  A string to raise an exception
        """

        raise NotImplementedError("run_scoring() not implemented")
