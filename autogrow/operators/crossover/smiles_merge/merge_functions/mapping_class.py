"""
This script holds the Mapping class.
This is used when mapping most common substructure (MCS)
    to combine two molecules.
"""
import __future__

import random
import copy


class Mapping(object):
    """
    # Notes on terminology:
        -most common substructure (MCS): The substructure shared between
            the two parent ligands
        -node/anchor/I: an atom in the MCS which has 1 or more atom(s)
            connected to it which are not part of the MCS
                -the anchors are labeled by their Isotope numbers as
                    those do not get modified, where as atom Idx are modified
                    by many Rdkit functions.  Anchors have Isotope labels of
                    10,000 or higher and that label is applied to the MCS and
                    ligands 1 and 2 so everything is trackable.
        -R-group: a chain of one or more atoms connected to a single node
                -if an anchor has dimethyl's which are not part of MCS then
                    each methyl is considered its own R-group
        -B-group: 1 or more R-groups which branch off a single node.
                -if an anchor has dimethyl's which are not part of MCS then
                    the combination of both methyls is considered a single
                    B-group.

                    B-group Naming scheme: '{first_number}B{second_number}'
                        -first_number: the number before the B corresponds to
                            the parent ligand from which the B-group is
                            derived.
                        -second_number: the number which follows the B is the
                            order for which that B-group was determined when
                            condensing the R-groups into B-groups. Numbering
                            is indexed to 1. So the 1st three B groups for
                            parent ligand 1 are: 1B1,1B2,1B3
                        ie) 1B1 is the 1st B-group from parent ligand 1
                            1B2 is the second B-group from parent ligand 1
                            2B1 is the 1st B-group from parent ligand 2
                            2B2 is the second B-group from parent ligand 2


    This class handles mapping for Bs and Is to chose B-groups which will
    later be used to make a child molecule. All the choices for B-groups are
    handled here.

    This is important because if a B-group connects to more than one anchor
    atom, the selection of that B-group determines the selection of both
    anchor atoms.

    ie) if 1B1 connects to anchor atom 10003 and 10004; and 2B1 connects to
        10003 and 10005 then the decision for which B-group is chosen for
        anchor 10003 determines the options which will be viable for anchor
        atoms 10003,10004, and 10005.

    These type of decisions are handled by this class.
    """

    def __init__(self, b_to_is, i_to_bs):
        """
        When a Mapping object is initialized, it imports 2 input dictionaries,
        which can be referenced throughout the class.

        Inputs:
        :param dict b_to_is: Dictionary converting B-groups to anchor/node/I
            atoms. This contains groups from both parent molecules. This is the
            inverse of i_to_bs. keys are B-groups; items are the anchor atoms
            isotope label. ie) {'1B1': [10003], '2B3': [10005], '2B2': [10003],
            '2B1': [10002]}
        :param dict i_to_bs: Dictionary converting Anchor/node/I atoms to
            corresponding B-groups. This contains groups from both parent
            molecules. This is the inverse of b_to_is. keys are the anchor
            atoms
            isotope labels; items are B-groups ie) {10002: ['2B1'], 10003: ['1B1',
            '2B2'], 10005: ['2B3']}
        """

        self.b_to_is = copy.deepcopy(
            b_to_is
        )  # B-I mapping dictionary from outside class
        self.i_to_bs = copy.deepcopy(
            i_to_bs
        )  # I-B mapping dictionary from outside class

    def locate_b(self, i):
        """
        Given a specified anchor/I return a list of all the B-groups from both
        parent ligands, bound to that anchor

        Inputs:
        :param int i: the isolabel of an anchor atom which will be used to
            search for B-groups within self.i_to_bs.

        Returns
        :returns: list self.i_to_bs[i]: A list of all the B-groups, from both
            parent ligands which are bound to that anchor. ie) ['1B1','2B1']
        """
        return self.i_to_bs[i]
    #

    def locate_i(self, b):
        """
        Given a specified B-group return the anchor/I/node it connects to.

        Inputs:
        :param str b: the name of a B-groups within self.b_to_is.

        Returns
        :returns: list self.b_to_is[b]: A list of the anchor the given
            B-groups is bound. ie) [10001]
        """
        return self.b_to_is[b]
    #

    def delete_b(self, b):
        """
        Removes the b from b_to_is and all references to b in i_to_bs.
            b is a Key in b_to_is. B is one or more items in i_to_bs.

        Inputs:
        :param str b: A B-group to be removed from the b_to_is and B in
            i_to_bs dicts.
        """

        i_list_to_modify = self.locate_i(b)
        for i in i_list_to_modify:
            blank = self.i_to_bs[i].remove(b)
        del self.b_to_is[b]
    #

    def delete_i(self, i):
        """
        Removes the i from i_to_bs and all references to i in b_to_is. i is
            a Key in i_to_bs. i is one or more items in b_to_is.

        Inputs:
        :param int i: An interger representing the isolabel for an
            anchor/node/i atom to be removed from the b_to_is and b in i_to_bs
            dicts.
        """

        bs_to_modify = self.locate_b(i)
        for b in bs_to_modify:
            self.b_to_is[b].remove(i)
        del self.i_to_bs[i]
    #

    def chose_b_from_i(self, i):
        """
        Chose your B from a given i. This makes the decision which B-group
        will be chosen for a specific i.

        Current implementation is that there are no null choice options. ie.
        if an anchor has only 1 option to pick from then it must pick that
        B-group. It can not chose nothing or to leave it blank, even if
        choosing that B-group forces the future decision because of it's
        connections to anchors which have yet to have B-group decisions.

        this has bearings on B-groups which connect to multiple anchors as
        well as on anchors which have B-groups from only one parent ligand,
        but the other parent has nothing connected to that anchor.

        in the case of one parent having a B-group attached to an anchor but
        nothing attached to the anchor for the other parent, this
        implementation will always chose to keep the B-group and never can
        leave it blank.

        ie (if 1B1 is connected to multiple anchors)

        Lack of an B-groups bound to an anchor is not considered a B-group

        Inputs:
        :param int i:   An interger representing the isolabel for an
            anchor/node/i atom. This function choses which B-group will be bound
            to this anchor in the child molecule.

        Returns:
        :returns: str b_x: A string of the name of a chosen B-group; None if
            not in the dictionary or if there is no available choices
        """

        # Developers Notes

        # Current implementation has no Null/None as choices this means that
        # if one parent has a B-group bound to anchor and the other has
        # nothing bound to that anchor, then the program will always chose to
        # add the B-group, resulting in a larger child molecule.

        # This biases the output. It also results in B-groups with multiple
        # connections are weighted against because 1 decision on 1 node will
        # determine if they can't be chosen...

        # Two alternative implementations which could be added are listed
        # below, both with advantages and disadvantages:

        # 1) Dominant Nulls: (Null is an option which cannot be override)

        #  When we reach an anchor which only has only one B-group choice then
        #  we add a Null B-group. This Null group means nothing can be added
        #  to that anchor.

        #  This could be implemented in this step in mapping_class.py or this
        #  could be implemented at the R-groups to B-group consolidation step
        #  -implemented it at the R-group to B-group consolidation may be a
        #  better option because it will simplify issues of when is it
        #  appropriate to add Null.
        #    ie) if parent lig_1 has 1B1 attached to anchors 10000,10001,10002
        #                            1B2 attached to anchors 10003
        #        if parent lig_2 has 2B1 at 10000 and 10003
        #                            2B at 10002
        #        if 1B1 is chosen 1st using anchor 10000, then anchors 10000,10001,10002
        #        are determined

        #            THIS ALSO MEANS THAT 2B1 is impossible eliminating 2B1
        #            from an option for anchor 10003

        #            When the program needs to make a decision for anchor
        #            10003 what should its option be????:
        #                - It should have to chose 1B2
        #                   IF WE implement THE NULLS IN THIS PART OF THE CODE WE
        #                   WOULD HAVE TO CODE IN THAT AS A CONSIDERATION IF WE
        #                   implement NULLS AT THE R- TO B-GROUP CONSOLIDATION PHASE
        #                   WE WOULDN'T NEED TO ADD EXTRA CODE HERE TO PREVENT A NULL
        #                   FROM BEING ADDED

        #  If a dominant null group (which is only added when an anchor has no
        #  R-groups attached for 1 parent but some on the other) then when the
        #  decision for B-group is occuring here; if a null is chosen then
        #  nothing can be added

        #  Effects:
        #  1) easier to implement
        #  2) affects the weighting of multi connection point B-groups -A
        #     single decision more dramatically impacts chains with many
        #     connections
        #  3) a Null is permanent so its easier to code and process

        # 2) Recessive Nulls: (A Null can be chosen but it can be overriden)
        #     (soft Null). If a recessive null is chosen instead of a B-group
        #     with multiple connections to the MCS then the B-group which
        #     wasn't chosen does not get removed from the dictionaries.
        #  -Then the next anchor that the not chosen B-group is connected
        #   to is assessed, there is still the option to chose that group.
        #    -If that group is chosen then we write over the Null option.
        #    -If that group is not chosen then the null remains as the choice
        #     for the 1st anchor

        #  Effects:
        #  -Recessive Nulls favors the selection of B-groups with multiple
        #   connections, but still allows for a null to be chosen.
        #  -A more balanced option between No-Nulls (the current
        #   implementation) and Dominant Nulls. But this does bias the
        #   statistics of choices
        #  -this also makes the decision tree more complicated and makes
        #   coding this more difficult

        #  There's no right answer to this but there are certain pro's and
        #  con's to each approach. The current approach is justified as the
        #  most code and computational efficient method, with no distict
        #  preference for multi-chain B-groups, but certainly with some biases
        #  against shrinking the child molecule

        # Select an B to keep
        if i in list(self.i_to_bs.keys()):

            options = self.locate_b(i)
            if len(options) > 1:
                b_x = random.choice(options)
            elif len(options) == 1:
                b_x = options[0]
            else:
                return  "None"

            list_is = self.locate_i(b_x)
            list_bs = []
            for x in list_is:
                list_bs.append(self.locate_b(x))

            flattened = [val for sublist in list_bs for val in sublist]
            unique_bs = list(
                set(flattened)
            )  # convert list to set to list to remove redundant B's
            # delete the B's and I's
            for b in unique_bs:
                self.delete_b(b)

            for x in list_is:
                self.delete_i(x)
            return b_x

        # the i is not in list(self.i_to_bs.keys())
        # return the string "None"
        return "None"
#
    def testing_function_return_self_dicts(self):
        """
        Return the properties: self.b_to_is and self.i_to_bs

        Returns:
        :returns: dict b_to_is: Dictionary converting B-groups to
            anchor/node/I atoms. This contains groups from both parent molecules.
            This is the inverse of i_to_bs. keys are B-groups; items are the
            anchor atoms isotope label. ie) {'1B1': [10003], '2B3': [10005],
            '2B2': [10003], '2B1': [10002]}
        :returns: dict i_to_bs: Dictionary converting Anchor/node/I atoms to
            corresponding B-groups. This contains groups from both parent
            molecules. This is the inverse of b_to_is. keys are the anchor atoms
            isotope labels; items are B-groups. ie) {10002: ['2B1'], 10003:
            ['1B1', '2B2'], 10005: ['2B3']}
        """
        return self.b_to_is, self.i_to_bs

# i_dict = {10000: ['1B1', '2B1'], 10004: ['2B2'], 10005: ['2B3'], 10006: \
#       ['2B4'], 10007: ['1B3'], 10008: ['1B2']}
# b_dict = {'1B1': [10000], '1B2': [10008], '1B3': [10007], '2B4': [10006], \
# '2B3': [10005], '2B2': [10004], '2B1': [10000]}
def run_mapping(b_dict, i_dict):
    """
    This runs the mapping class which can determine which B-groups/R-groups we
    will append in SmileMerge.

    Inputs:
    :param dict b_dict: Dictionary converting B-groups to anchor/node/I
        atoms. This contains groups from both parent molecules. This is the
        inverse of i_to_bs. keys are B-groups; items are the anchor atoms isotope
        label. ie) {'1B1': [10003], '2B3': [10005], '2B2': [10003], '2B1':
        [10002]}
    :param dict i_dict: Dictionary converting Anchor/node/I atoms to
        corresponding B-groups. This contains groups from both parent molecules.
        This is the inverse of b_to_is. keys are the anchor atoms isotope labels;
        items are B-groups. ie) {10002: ['2B1'], 10003: ['1B1', '2B2'], 10005:
        ['2B3']}

    Returns:
    :returns: list bs_chosen: A list of all the chosen B-groups to be used to
        generate a child molecule later.
    """
    a_mapping_object = Mapping(b_dict, i_dict)
    bs_chosen = []
    for i in i_dict:
        b_choice = a_mapping_object.chose_b_from_i(i)
        bs_chosen.append(b_choice)

    bs_chosen = list(set(bs_chosen))

    for i in bs_chosen:
        if i == "None":
            bs_chosen.remove(i)

    return bs_chosen
    #
#
