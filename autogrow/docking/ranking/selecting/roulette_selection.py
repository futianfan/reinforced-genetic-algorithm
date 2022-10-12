"""
This script is use to select molecules using a roulette selector
"""
import __future__

import numpy.random as rn


def spin_roulette_selector(usable_list_of_smiles, number_to_chose,
                           docking_or_diversity):
    """
    make a list of ligands chosen by a random weighted roulette selection,
    without replacement, weighted by its docking score

    Inputs:
    :param list usable_list_of_smiles: a list with all the information of all
        the mols in the previous generation
    :param int number_to_chose: the number of molecules to chose based on
        docking score
    :param str docking_or_diversity: an string describing either "docking" or
        "diversity" this tells the function how to adjust the weighted scores

    Returns:
    :returns: list top_choice_smile_order: list of ligands chosen by a random
        weighted selection, without replacement, -weighted by its docking score
    """

    if type(usable_list_of_smiles) is not type([]):
        raise Exception("usable_list_of_smiles Must be a list, wrong data type")

    num_ligands = len(usable_list_of_smiles)
    if num_ligands == 0:
        raise Exception(
            "usable_list_of_smiles is an empty list. There is nothing to chose from."
        )

    if number_to_chose <= 0:
        top_choice_smile_order = []
        return top_choice_smile_order

    adjusted = adjust_scores(usable_list_of_smiles, docking_or_diversity)

    total = sum(adjusted)
    probability = [x / total for x in adjusted]
    smiles_list = [x[0] for x in usable_list_of_smiles]

    top_choice_smile_order = rn.choice(
        smiles_list, size=number_to_chose, replace=False, p=probability
    )

    return top_choice_smile_order


def adjust_scores(usable_list_of_smiles, docking_or_diversity):
    """
    This function adjusts the scores appropriately. This is where we weight
    the scores so smaller differences are more pronounced and where we adjust
    for the fact that docking score is better with a more negative number
    while diversity score is the smallest positive number is the most unique.

    Inputs:
    :param list usable_list_of_smiles: a list with all the information of all
        the mols in the previous generation
    :param str docking_or_diversity: an string describing either "docking"
        or "diversity" this tells the function how to adjust the weighted
        scores

    Returns:
    :returns: list adjusted: list of ligand scores which have been weighted
        and adjusted
    """

    if docking_or_diversity == "docking":
        weight_scores = [float(x[-1]) for x in usable_list_of_smiles]
        # minimum is the most positive value from usable_list_of_smiles the
        # more negative the docking score the better the dock
        minimum = max(weight_scores) + 0.1
        if minimum < 0:
            minimum = 0

        adjusted = [(x ** 10) + minimum for x in weight_scores]

    elif docking_or_diversity == "diversity":

        weight_scores = [float(x[-1]) for x in usable_list_of_smiles]

        # adjust by squaring the number to make the discrpency larger and
        # invert by dividing 1/x^2 (because the more diverse a mol is the
        # smaller the number)
        adjusted = [(x ** -2) for x in weight_scores]

    else:
        raise Exception("docking_or_diversity choice not an option")

    return adjusted
