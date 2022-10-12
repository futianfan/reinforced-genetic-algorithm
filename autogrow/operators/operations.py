"""
Populates an AutoGrow generation via mutation, crossover, and elitism.
Also filters and converts SMILES to 3d SDFS.
"""
import __future__

import os
import random
import copy
import sys

import rdkit
import rdkit.Chem as Chem

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.filter.execute_filters as Filter
import autogrow.docking.ranking.ranking_mol as Ranking
import autogrow.operators.mutation.execute_mutations as Mutation
import autogrow.operators.crossover.execute_crossover as execute_crossover
import autogrow.operators.convert_files.conversion_to_3d as conversion_to_3d
import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH



#############
# Main run Autogrow operators to make a generation
#############
def populate_generation(vars, generation_num):
    """
    This will run all of the mutations, crossovers, and filters for a single
        generation. Populates a new generation of ligands.

    Inputs:
    :param dict vars: a dictionary of all user variables
    :param int generation_num: the generation number, 1,2,3,4,...,

    Returns:
    :returns: str full_generation_smiles_file: the name of the .smi file
        containing the new population
    :returns: list full_generation_smiles_list: list with the new population
        of ligands
    :returns: bool None: returns None twice if any step failed. This will
        result in the program ending
    """
    number_of_processors = int(vars["number_of_processors"])

    # Determine which generation it is and how many mutations and crossovers
    # to make
    if generation_num == 1:
        # If 1st generation
        num_crossovers = vars["number_of_crossovers_first_generation"]
        num_mutations = vars["number_of_mutants_first_generation"]

        # How many advance from previous generation to the next generation
        # directly This will be done later but we are unpacking vars here
        num_elite_to_advance_from_previous_gen = vars["number_elitism_advance_from_previous_gen_first_generation"]
    else:
        # Later generations
        num_crossovers = vars["number_of_crossovers"]
        num_mutations = vars["number_of_mutants"]
        num_elite_to_advance_from_previous_gen = vars["number_elitism_advance_from_previous_gen"]

    # Get the Source compound list. This list is the full population from
    # either the previous generations or if its Generation 1 than the its the
    # entire User specified Source compound list If either has a SMILES that
    # does not sanitize in RDKit it will be excluded and a printout of its
    # Name and SMILES string will be printed.
    source_compounds_list = get_complete_list_prev_gen_or_source_compounds(
        vars, generation_num
    )

    num_seed_diversity, num_seed_dock_fitness = determine_seed_population_sizes(
        vars, generation_num
    )

    # Total Population size of this generation
    total_num_desired_new_ligands = (
        num_crossovers + num_mutations + num_elite_to_advance_from_previous_gen
    )

    ###################################################################
    ###################################################################
    # A. Making Mutations
    ###################################################################
    ###################################################################
    # Get starting compounds for Mutations
    seed_list_mutations = make_seed_list(
        vars,
        source_compounds_list,
        generation_num,
        num_seed_diversity,
        num_seed_dock_fitness,
    )

    # Save seed list for Mutations
    save_ligand_list(
        vars["output_directory"],
        generation_num,
        seed_list_mutations,
        "Mutation_Seed_List",
    )
    sys.stdout.flush()
    print("MAKE MUTATIONS")


    # Package user vars specifying the Reaction library to use for mutation
    rxn_library_variables = [
        vars["rxn_library"],
        vars["rxn_library_file"],
        vars["function_group_library"],
        vars["complementary_mol_directory"],
    ]

    new_mutation_smiles_list = [] # List of SMILES from mutation

    # Make all the required ligands by mutations
    while len(new_mutation_smiles_list) < num_mutations:
        sys.stdout.flush()

        num_mutants_to_make = num_mutations - len(new_mutation_smiles_list)

        ##################################
        # Make all mutants -------- main 
        ##################################
        new_mutants = Mutation.make_mutants(
            vars,
            generation_num,
            number_of_processors,
            num_mutants_to_make,
            seed_list_mutations,
            new_mutation_smiles_list,
            rxn_library_variables,
        )
        if new_mutants is None:
            # try once more
            new_mutants = Mutation.make_mutants(
                vars,
                generation_num,
                number_of_processors,
                num_mutants_to_make,
                seed_list_mutations,
                new_mutation_smiles_list,
                rxn_library_variables,
            )

        if new_mutants is None:
            break

        # Remove Nones:
        new_mutants = [x for x in new_mutants if x is not None]

        for i in new_mutants:
            new_mutation_smiles_list.append(i)
            if len(new_mutation_smiles_list) == num_mutations:
                break
    sys.stdout.flush()

    # save new_mutation_smiles_list
    save_ligand_list(
        vars["output_directory"],
        generation_num,
        new_mutation_smiles_list,
        "Chosen_Mutants",
    )

    if (
            new_mutation_smiles_list is None
            or len(new_mutation_smiles_list) < num_mutations
    ):
        print("")
        print("")
        print("We needed to make {} ligands through Mutation".format(num_mutations))
        print(
            "We only made {} ligands through Mutation".format(
                len(new_mutation_smiles_list)
            )
        )
        print("")
        print("")
        raise Exception("Mutation failed to make enough new ligands.")

    print("FINISHED MAKING MUTATIONS")


    ###################################################################
    ###################################################################
    # B. Making Crossovers
    # List of smiles from crossover
    ###################################################################
    ###################################################################

    # Get starting compounds to seed Crossovers
    seed_list_crossovers = make_seed_list(
        vars,
        source_compounds_list,
        generation_num,
        num_seed_diversity,
        num_seed_dock_fitness,
    )

    # Save seed list for Crossovers
    save_ligand_list(
        vars["output_directory"],
        generation_num,
        seed_list_crossovers,
        "Crossover_Seed_List",
    )

    print("MAKE CROSSOVERS")
    sys.stdout.flush()

    new_crossover_smiles_list = []

    # Make all the required ligands by Crossover
    while len(new_crossover_smiles_list) < num_crossovers:
        sys.stdout.flush()
        num_crossovers_to_make = num_crossovers - len(new_crossover_smiles_list)

        # Make all crossovers
        new_crossovers = execute_crossover.make_crossovers(
            vars,
            generation_num,
            number_of_processors,
            num_crossovers_to_make,
            seed_list_crossovers,
            new_crossover_smiles_list,
        )
        if new_crossovers is None:
            # try once more
            new_crossovers = execute_crossover.make_crossovers(
                vars,
                generation_num,
                number_of_processors,
                num_crossovers_to_make,
                seed_list_crossovers,
                new_crossover_smiles_list,
            )
        if new_crossovers is None:
            break

        # Remove Nones:
        new_crossovers = [x for x in new_crossovers if x is not None]

        # append those which passed the filter
        for i in new_crossovers:
            new_crossover_smiles_list.append(i)
            if len(new_crossover_smiles_list) == num_crossovers:
                break

    # save new_crossover_smiles_list
    save_ligand_list(
        vars["output_directory"],
        generation_num,
        new_crossover_smiles_list,
        "Chosen_Crossovers",
    )

    if (
            new_crossover_smiles_list is None
            or len(new_crossover_smiles_list) < num_crossovers
    ):
        print("")
        print("")
        print("We needed to make {} ligands through Crossover".format(num_crossovers))
        print(
            "We only made {} ligands through Crossover".format(
                len(new_crossover_smiles_list)
            )
        )
        print("")
        print("")
        raise Exception("Crossover failed to make enough new ligands.")

    print("FINISHED MAKING CROSSOVERS")

    # Get unaltered samples from the previous generation
    print("GET SOME LIGANDS FROM THE LAST GENERATION")
    sys.stdout.flush()

    # Make a list of the ligands chosen to pass through to the next generation
    # via Elitism This handles creating a seed list and defining the advance
    # to next generation final selection

    chosen_mol_to_pass_through_list = make_pass_through_list(
        vars,
        source_compounds_list,
        num_elite_to_advance_from_previous_gen,
        generation_num,
    )

    if type(chosen_mol_to_pass_through_list) == str:
        printout = (
            chosen_mol_to_pass_through_list
            + "\nIf this is the 1st generation, it may be due to the starting "
            + "library has SMILES which could not be converted to sanitizable "
            + "RDKit Molecules"
        )

        raise Exception(printout)
    sys.stdout.flush()

    # save chosen_mol_to_pass_through_list
    save_ligand_list(
        vars["output_directory"],
        generation_num,
        chosen_mol_to_pass_through_list,
        "Chosen_Elite_To_advance",
    )

    print("GOT LIGANDS FROM THE LAST GENERATION")

    # make a list of all the ligands from mutations, crossovers, and from the
    # last generation
    new_generation_smiles_list = []
    full_generation_smiles_list = []
    for i in new_mutation_smiles_list:
        new_generation_smiles_list.append(i)
        full_generation_smiles_list.append(i)

    for i in new_crossover_smiles_list:
        new_generation_smiles_list.append(i)
        full_generation_smiles_list.append(i)

    if vars["redock_elite_from_previous_gen"] is False and generation_num != 1:
        for i in chosen_mol_to_pass_through_list:
            # Doesn't append to the new_generation_smiles_list
            full_generation_smiles_list.append(i)

    # Generation 0 pass through gets added to the convert and dock list
    # because it has no docking score to compare with This is independent of
    # the vars['redock_elite_from_previous_gen']
    else:
        for i in chosen_mol_to_pass_through_list:
            new_generation_smiles_list.append(i)
            full_generation_smiles_list.append(i)

    if len(full_generation_smiles_list) < total_num_desired_new_ligands:
        print("We needed ", total_num_desired_new_ligands)
        print("We made ", len(full_generation_smiles_list))
        print(
            "population failed to make enough mutants or crossovers... \
            Errors could include not enough diversity, too few seeds to \
            the generation, the seed mols are unable to cross-over due \
            to lack of similariy, or all of the seed lack functional groups \
            for performing reactions"
        )
        return None, None, None

    # Save the Full Generation
    full_generation_smiles_file, new_gen_folder_path = save_generation_smi(
        vars["output_directory"], generation_num, full_generation_smiles_list, None
    )

    # Save the File to convert to 3d
    smiles_to_convert_file, new_gen_folder_path = save_generation_smi(
        vars["output_directory"],
        generation_num,
        new_generation_smiles_list,
        "_to_convert",
    )

    sys.stdout.flush()
    # CONVERT SMILES TO .sdf USING GYPSUM and convert .sdf to .pdb with rdkit
    # This will output sdf files into a folder. The .smi.0.sdf file is not a
    # valid mol, but all the others will be valid the 1st Smiles in the
    # original .smi file is saved as .smi.1.sdf and 2nd file is saved as
    # .smi.2.sdf
    conversion_to_3d.convert_to_3d(vars, smiles_to_convert_file, new_gen_folder_path)
    sys.stdout.flush()

    return full_generation_smiles_file, full_generation_smiles_list


def populate_generation_zero(vars, generation_num=0):
    """
    This will handle all that is required for generation 0redock and handle
    the generation 0

    Inputs:
    :param dict vars: a dictionary of all user variables
    :param int generation_num: the generation number

    Returns:
    :returns: str full_generation_smiles_file: the name of the .smi file
        containing the new population
    :returns: list full_generation_smiles_list: list with the new population
        of ligands.
    :returns: bool already_docked: if true we won't redock the source ligands.
        If False we will dock the source ligands.
    """
    number_of_processors = int(vars["number_of_processors"])

    num_crossovers = 0
    num_mutations = 0

    # Get the Source compound list. This list is the full population from
    # either the previous generations or if its Generation 1 than the its the
    # entire User specified Source compound list If either has a SMILES that
    # does not sanitize in RDKit it will be excluded and a printout of its
    # Name and SMILES string will be printed.
    source_compounds_list = get_complete_list_prev_gen_or_source_compounds(
        vars, generation_num
    )
    num_elite_to_advance_from_previous_gen = len(source_compounds_list)

    num_seed_diversity, num_seed_dock_fitness = determine_seed_population_sizes(
        vars, generation_num
    )

    # Total Population size of this generation
    total_num_desired_new_ligands = num_crossovers + num_mutations + 1

    # Get unaltered samples from the previous generation
    print("GET SOME LIGANDS FROM THE LAST GENERATION")

    # Make a list of the ligands chosen to pass through to the next generation
    # This handles creating a seed list and defining the advance to next
    # generation final selection
    chosen_mol_to_pass_through_list = make_pass_through_list(
        vars, source_compounds_list, 1, 0
    )

    if type(chosen_mol_to_pass_through_list) == str:
        printout = (
            chosen_mol_to_pass_through_list
            + "\nIf this is the 1st generation, it may be due to the starting "
            + "library has SMILES which could not be converted to "
            + "sanitizable RDKit Molecules"
        )

        raise Exception(printout)

    # save chosen_mol_to_pass_through_list
    save_ligand_list(
        vars["output_directory"],
        generation_num,
        chosen_mol_to_pass_through_list,
        "Chosen_Elite_To_advance",
    )

    print("GOT LIGANDS FROM THE LAST GENERATION")

    # make a list of all the ligands from mutations, crossovers, and from the
    # last generation
    new_generation_smiles_list = []
    full_generation_smiles_list = []

    # These will be docked and scored for generation 0
    for i in chosen_mol_to_pass_through_list:
        new_generation_smiles_list.append(i)
        full_generation_smiles_list.append(i)

    if len(full_generation_smiles_list) == 0:
        print(
            "population failed to import any molecules from the source_compounds_list."
        )
        return None, None, None

    # Save the Full Generation
    full_generation_smiles_file, new_gen_folder_path = save_generation_smi(
        vars["output_directory"], generation_num, full_generation_smiles_list, None
    )

    # Save the File to convert to 3d
    smiles_to_convert_file, new_gen_folder_path = save_generation_smi(
        vars["output_directory"],
        generation_num,
        new_generation_smiles_list,
        "_to_convert",
    )

    # order files by -2 of each lig
    try:
        full_generation_smiles_list.sort(key=lambda x: float(x[-2]), reverse=False)
        full_generation_smiles_list_printout = [
            "\t".join(x) for x in full_generation_smiles_list
        ]
        already_docked = True
    except:
        print(
            "Not all ligands in source compound list are scored. "
            + "We will convert and redock them all."
        )
        already_docked = False

    if already_docked is True:
        # Write all the ligands to a ranked file
        full_generation_smiles_list_printout = "\n".join(
            full_generation_smiles_list_printout
        )
        ranked_file = new_gen_folder_path + os.sep + "generation_0_ranked.smi"
        with open(ranked_file, "w") as f:
            f.write(full_generation_smiles_list_printout)
        return already_docked, full_generation_smiles_file, full_generation_smiles_list

    # If you are to redock and convert the generation zero you will also need
    # to do the following:

    # CONVERT SMILES TO .sdf USING GYPSUM and convert .sdf to .pdb with
    # rdkit This will output sdf files into a folder. The .smi.0.sdf file
    # is not a valid mol, but all the others will be valid the 1st Smiles
    # in the original .smi file is saved as .smi.1.sdf and 2nd file is
    # saved as .smi.2.sdf
    conversion_to_3d.convert_to_3d(
        vars, smiles_to_convert_file, new_gen_folder_path
    )

    return already_docked, full_generation_smiles_file, full_generation_smiles_list


#############
# Get seeds
#############
def test_source_smiles_convert(smile_info):
    """
    This attempts to convert a SMILES string to an rdkit.Chem.rdchem.Mol
    object
        -done in a try statement so that it is some bad SMILES string which is
        incapable of being converted
        - it also checks that the SMILES string is able to be sanitized

    Inputs:
    :param list smile_info: a list containing the SMILES of a ligand, its ID
        and potentially additional information about the ligand

    Returns:
    :returns: list smile_info: If it passed the test, it returns the list
        containing the SMILES of a ligand, its ID and potentially additional
        information about the ligand
    :returns: str printout: If it failed to convert it returns the error
        message. This passess out to prevent MPI print issues
    """
    if smile_info is None or len(smile_info) == 0:
        printout = (
            "REMOVING SMILES FROM SOURCE LIST: Blank "
            + "entry in source compound list.\n"
        )
        printout = printout + "\tRemoving: {}".format(smile_info)
        return printout

    if len(smile_info) == 1:
        printout = "REMOVING SMILES FROM SOURCE LIST: Unformatted or blank "
        printout = printout + "entry in source compound list.\n"
        printout = printout + "\tRemoving: {}".format(smile_info)
        return printout

    # separate out SMILES str and ID
    smile_str = smile_info[0]
    smile_id = str(smile_info[1])

    if type(smile_str) is not type(""):
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string is not a "
        printout = printout + "String. Check for formatting errors. \n"
        printout = printout + "\tIgnored SMILES is: {}".format(smile_str)
        return printout

    # Try importing it into RDKit with Sanitization off. Tests for errors in
    # having the wrong data type
    try:
        mol = Chem.MolFromSmiles(str(smile_str), sanitize=False)
    except:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string failed "
        printout = printout + "to import into RDKit.\n\t "
        printout = printout + "Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # This will fail if there are valence errors. We won't try to correct
    # someones source compound list Although the MOH.check_sanitization will
    # do that. try sanitizing, which is necessary later
    try:
        Chem.SanitizeMol(mol)
    except:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES "
        printout = printout + "string failed to Sanitize in RDKit.\n"
        printout = printout + "\t Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # Make the mol again fresh and try running it through MOH.handleHs() This
    # will try protanating and Deprotanating the mol. If it can't handle that
    # We reject it as many functions will require this sort of manipulation.
    # More advanced sanitization issues will also be removed in this step
    mol = Chem.MolFromSmiles(str(smile_str), sanitize=False)
    mol = MOH.handleHs(mol, True)

    if mol is None:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string failed \
                    to be protanated or deprotanated.\n"
        printout = (
            printout
            + "\t This is often an issue with valence and sanitization "
            + "issues with the SMILES string."
        )
        printout = printout + "\t Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # Check there are no * which are atoms with atomic number=0
    mol = MOH.check_for_unassigned_atom(mol)
    if mol is None:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string contained "
        printout = printout + "an unassigned atom type labeled as *.\n"
        printout = printout + "\t Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # Check for fragments.
    if len(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)) != 1:

        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string was fragmented.\n"
        printout = printout + "\t Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # the ligand is good enough to use throughout the program!
    return smile_info


def get_complete_list_prev_gen_or_source_compounds(vars, generation_num):
    """
    Get the source compounds list from the previous generation of the source
    compound list

    This also filters the list to ensure mols can be imported to RDKit and
    that they pass the drug-likliness filters.

    If generation = 1 use the User specified starting compounds If generation
    is >1 than use the previous generations top ligands. This takes an .smi
    file

    Inputs:
    :param dict vars: a dictionary of all user variables
    :param int generation_num: the interger of the current generation

    Returns:
    :returns: list usable_list_of_smiles: a list with SMILES strings, names,
        and information about the smiles from the previous generation or the
        source compound list
    """
    source_file_gen_0 = vars[
        "output_directory"
    ] + "generation_{}{}generation_{}_ranked.smi".format(0, os.sep, 0)
    if generation_num == 0:
        # This will be the full length list of starting molecules as the seed
        source_file = str(vars["source_compound_file"])
        usable_list_of_smiles = Ranking.get_usable_format(source_file)

        if len(usable_list_of_smiles) == 0:
            print(
                "\nThere were no available ligands in source compound. Check formatting\n"
            )
            raise Exception(
                "There were no available ligands in source compound. Check formatting"
            )

    elif generation_num == 1 and os.path.exists(source_file_gen_0) is False:
        # This will be the full length list of starting molecules as the seed
        source_file = str(vars["source_compound_file"])
        usable_list_of_smiles = Ranking.get_usable_format(source_file)

        if len(usable_list_of_smiles) == 0:
            print(
                "\nThere were no available ligands in source compound. Check formatting\n"
            )
            raise Exception(
                "There were no available ligands in source compound. Check formatting"
            )

    else:
        source_file = vars[
            "output_directory"
        ] + "generation_{}{}generation_{}_ranked.smi".format(
            generation_num - 1, os.sep, generation_num - 1
        )
        if os.path.exists(source_file) is False:
            printout = (
                "\n"
                + "There were no available ligands in previous"
                + " generation ranked ligand file.\n"
            )
            printout = printout + "\tCheck formatting or if file has been moved.\n"
            print(printout)
            raise Exception(printout)

        usable_list_of_smiles = Ranking.get_usable_format(source_file)

        if len(usable_list_of_smiles) == 0:
            printout = (
                "\n"
                + "There were no available ligands in previous"
                + " generation ranked ligand file.\n"
            )
            printout = printout + "\tCheck formatting or if file has been moved. \n"
            print(printout)
            raise Exception(printout)

    # Test that every SMILES in the usable_list_of_smiles is a valid SMILES
    # which will import and Sanitize in RDKit. SMILES will be excluded if they
    # are fragmented, contain atoms with no atomic number (*), or do not
    # sanitize
    job_input = tuple([tuple([i]) for i in usable_list_of_smiles])

    usable_list_of_smiles = vars["parallelizer"].run(
        job_input, test_source_smiles_convert
    )
    usable_list_of_smiles = [x for x in usable_list_of_smiles if x is not None]
    print_errors = [x for x in usable_list_of_smiles if type(x) is str]
    usable_list_of_smiles = [x for x in usable_list_of_smiles if type(x) is list]
    for x in print_errors:
        print(x)

    if len(usable_list_of_smiles) == 0:
        printout = "\nThere were no ligands in source compound or previous \
            generation which could sanitize.\n"
        print(printout)
        raise Exception(printout)

    if vars["filter_source_compounds"] is True:

        prefilter_list = copy.deepcopy(usable_list_of_smiles)
        print("")
        print("Running Filter on the Compounds from last generation/Source")
        usable_list_of_smiles = Filter.run_filter(vars, usable_list_of_smiles)

        # Remove Nones:
        usable_list_of_smiles = [x for x in usable_list_of_smiles if x is not None]

        if len(usable_list_of_smiles) == 0:
            printout = "\nThere were no ligands in source compound which \
                        passed the User-selected Filters.\n"
            print(printout)
            raise Exception(printout)

        for lig in usable_list_of_smiles:
            failed_filter_list = []
            if lig not in prefilter_list:
                failed_filter_list.append(lig[1])

        if len(failed_filter_list) != 0:
            printout = "\n THE FOLLOWING LIGANDS WERE REMOVED FROM THE\
                        SOURCE LIST: Failed the User-selected Filters\n"
            printout = printout + "\t{}".format(failed_filter_list)
            print(printout)

    random.shuffle(usable_list_of_smiles)

    return usable_list_of_smiles


def make_seed_list(vars, source_compounds_list, generation_num, num_seed_diversity,
                   num_seed_dock_fitness):
    """
    Get the starting compound list for running the Mutation and Crossovers

    If generation = 0 use the User specified starting compounds If generation
    is >0 than use the previous generations top ligands. This takes an .smi
    file

    Inputs:
    :param dict vars: a dictionary of all user variables
    :param list source_compounds_list: a list with SMILES strings, names, and
        information about the smiles from either the previous generation or the
        source compound list
    :param int generation_num: the interger of the current generation
    :param int num_seed_diversity: the number of seed molecules which come
        from diversity selection
    :param int num_seed_dock_fitness: the number of seed molecules which come
        from eite selection by docking score

    Returns:
    :returns: list usable_list_of_smiles: a list with SMILES strings, names,
        and information about the smiles which will be used to seed the next
        generation
    """
    usable_list_of_smiles = copy.deepcopy(source_compounds_list)

    full_length = False
    if generation_num == 0:
        # Get starting compounds for Mutations
        full_length = True
    elif generation_num == 1:
        if vars["use_docked_source_compounds"] is False:
            # Get starting compounds for Mutations
            full_length = True
        else:
            source_file_gen_0 = vars[
                "output_directory"
            ] + "generation_{}{}generation_{}_ranked.smi".format(0, os.sep, 0)
            if os.path.exists(source_file_gen_0) is False:
                full_length = True

            else:
                # generation_num 1 may run into problems if the source
                # compounds are smaller than the seed pool required to seed
                # generation 1. Because the seeding options are connected to
                # the generation number (due to the depreciation of diversity
                # option) Because of this we may need to ignore the ranking
                # for the seeds of generation 1 to accomidate the smaller
                # source size. This is especially important with
                # lead-optimization in which the source pool may be much
                # smaller For this reason we will override the seeding of
                # generation 1 if the number to seed is greater than exists
                # but will provide a warning message.
                if (
                        len(usable_list_of_smiles) < num_seed_diversity
                        or len(usable_list_of_smiles) < num_seed_diversity
                ):
                    # This is problematic so just use what is available
                    printout = "\n\nNot enough ligands in source compound \
                        list to seed generation 1. We will use the entire \
                        list of every ligand in the source compound list \
                        to seed generation 1. This means there is no \
                        selection in generation 1's seeding process.\n\n"
                    print(printout)
                    full_length = True
                else:
                    full_length = False
    else:
        full_length = False

    if full_length is True or generation_num == 0:
        # This will be the full length list of starting molecules as the seed
        random.shuffle(usable_list_of_smiles)

    else:
        selector_choice = vars["selector_choice"]
        tourn_size = vars["tourn_size"]
        # Get subset of the source_file based on diversity scores and docking
        # scores
        usable_list_of_smiles = Ranking.create_seed_list(
            usable_list_of_smiles,
            num_seed_diversity,
            num_seed_dock_fitness,
            selector_choice,
            tourn_size,
        )

    random.shuffle(usable_list_of_smiles)

    return usable_list_of_smiles


def determine_seed_population_sizes(vars, generation_num):
    """
    This function determines how many molecules will be chosen to seed a
    generation because of their docking score and how many will assess because
    of their diversity score.

    Inputs:
    :param dict vars: a dictionary of all user variables
    :param int generation_num: the interger of the current generation

    Returns:
    :returns:  int num_seed_diversity: the number of seed molecules which come
        from diversity selection
    :returns:  int num_seed_dock_fitness: the number of seed molecules which
        come from eite selection by docking score

    """
    # How many fewer seed mols are chosen from diversity compared to the 1st
    # generation This is also how many more get chosen from elitist selection
    diversity_depreciation = (
        int(generation_num - 1) * vars["diversity_seed_depreciation_per_gen"]
    )

    if generation_num == 1:
        top_mols_to_seed_next_generation = vars[
            "top_mols_to_seed_next_generation_first_generation"
        ]

    else:
        top_mols_to_seed_next_generation = vars["top_mols_to_seed_next_generation"]

    # Number of mols chosen because of their diversity score
    num_seed_diversity = (
        vars["diversity_mols_to_seed_first_generation"] - diversity_depreciation
    )

    # Number of mols chosen because of their docking score. Elitist style
    # selection but this will be chosen using a weight roulette selector later
    # on.
    if num_seed_diversity <= 0:
        num_seed_dock_fitness = (
            top_mols_to_seed_next_generation
            + vars["diversity_mols_to_seed_first_generation"]
        )
        num_seed_diversity = 0
    elif num_seed_diversity > 0:
        num_seed_dock_fitness = (
            top_mols_to_seed_next_generation + diversity_depreciation
        )

    return num_seed_diversity, num_seed_dock_fitness


def make_pass_through_list(vars, smiles_from_previous_gen_list,
                           num_elite_to_advance_from_previous_gen,
                           generation_num):
    """
    This function determines the molecules which elite ligands will advance
    from the previous generation without being altered into the next
    generation.

    Inputs:
    :param dict vars: a dictionary of all user variables
    :param list smiles_from_previous_gen_list: List of SMILES from the last
        generation chosen to seed the list of molecules to advance to the next
        generation without modification via elitism.
    :param int num_elite_to_advance_from_previous_gen: the number of molecules
        to advance from the last generation without modifications.
    :param int generation_num: the interger of the current generation

    Returns:
    :returns: list list_of_ligands_to_advance: a list of ligands which should
        advance into the new generation without modifications, via elitism from
        the last generation. Returns a printout of why it failed if it fails
    """
    # this will be a list of lists. Each sublist will be  [SMILES_string, ID]
    list_of_ligands_to_advance = []

    # If not enough of your previous generation sanitize to make the list
    # Return None and trigger an Error
    if (
            generation_num != 0
            and len(smiles_from_previous_gen_list) < num_elite_to_advance_from_previous_gen
    ):
        printout = "Not enough ligands in initial list the filter to progress"
        printout = (
            printout
            + "\n len(smiles_from_previous_gen_list): {} ; \
                num_elite_to_advance_from_previous_gen: {}".format(
                    len(smiles_from_previous_gen_list),
                    num_elite_to_advance_from_previous_gen)
        )
        return printout

    smiles_from_previous_gen_list = [
        x for x in smiles_from_previous_gen_list if type(x) == list
    ]

    if generation_num == 0 and vars["filter_source_compounds"] is True:
        # Run Filters on ligand list
        ligands_which_passed_filters = Filter.run_filter(
            vars, smiles_from_previous_gen_list
        )
        # Remove Nones:
        ligands_which_passed_filters = [
            x for x in ligands_which_passed_filters if x is not None
        ]
    else:
        ligands_which_passed_filters = [
            x for x in smiles_from_previous_gen_list if x is not None
        ]
    # If not enough of your previous generation sanitize to make the list
    # Return None and trigger an Error
    if (
            generation_num != 0
            and len(ligands_which_passed_filters) < num_elite_to_advance_from_previous_gen
    ):
        printout = "Not enough ligands passed the filter to progress"
        return printout

    # Save seed list of all ligands which passed which will serve as the seed
    # list.
    save_ligand_list(
        vars["output_directory"],
        generation_num,
        ligands_which_passed_filters,
        "Previous_Gen_Elite_Seed_List",
    )

    # check if ligands_which_passed_filters has docking scores
    has_dock_score = False
    try:
        temp = [float(x[-2]) for x in ligands_which_passed_filters]
        has_dock_score = True
    except:
        has_dock_score = False

    if generation_num == 0 and has_dock_score is False:
        # Take the 1st num_elite_to_advance_from_previous_gen number of
        # molecules from ligands_which_passed_filters
        random.shuffle(ligands_which_passed_filters)
        list_of_ligands_to_advance = []
        for x in range(0, len(ligands_which_passed_filters)):
            selected_mol = ligands_which_passed_filters[x]
            list_of_ligands_to_advance.append(selected_mol)
    elif generation_num == 0 and has_dock_score is True:
        # Use the make_seed_list function to select the list to advance.
        # This list will be chosen strictly by
        list_of_ligands_to_advance = make_seed_list(
            vars,
            ligands_which_passed_filters,
            generation_num,
            len(ligands_which_passed_filters),
            num_elite_to_advance_from_previous_gen,
        )

    elif generation_num != 0 and has_dock_score is False:
        # Take the 1st num_elite_to_advance_from_previous_gen number of
        # molecules from ligands_which_passed_filters
        random.shuffle(ligands_which_passed_filters)
        list_of_ligands_to_advance = []
        for x in range(0, num_elite_to_advance_from_previous_gen):
            selected_mol = ligands_which_passed_filters[x]
            list_of_ligands_to_advance.append(selected_mol)

    elif generation_num != 0 and has_dock_score is True:
        # Use the make_seed_list function to select the list to advance. This
        # list will be chosen strictly by
        list_of_ligands_to_advance = make_seed_list(
            vars,
            ligands_which_passed_filters,
            generation_num,
            0,
            num_elite_to_advance_from_previous_gen,
        )

    if generation_num == 0:
        return list_of_ligands_to_advance
    elif len(list_of_ligands_to_advance) >= num_elite_to_advance_from_previous_gen:
        return list_of_ligands_to_advance

    printout = "Not enough ligands were chosen to advance to the next generation."
    return printout


#############
# Saving Output files for generations and seeds
#############
def save_generation_smi(output_directory, generation_num,
                        formatted_smile_list, nomenclature_tag):
    """"
    This function saves a list of newly generated population of ligands as an
    .smi file. .smi file column 1 is the SMILES string and column 2 is its
    smile ID


    Inputs:
    :param dict output_directory: the directory of the run to save the
        generation
    :param int generation_num: the interger of the current generation
    :param list formatted_smile_list: list of the newly generated population
        of ligands
    :param str nomenclature_tag: The str describing the ligand list if None
        than don't add tag. It is the full list of all ligands for the
        generation
        If it says '_to_convert' its the list of ligands which will need to be
        converted to 3D -this may or may not have the ligands which pass
        through from the last gen.

    Returns:
    :returns: str output_file_name: name of the output file
    :returns: str new_gen_folder_path: the path to the folder containing all
        that will be in this generation
    """
    # folder for this new generation
    new_gen_folder_path = output_directory + "generation_{}{}".format(
        generation_num, os.sep
    )
    if nomenclature_tag is None:
        # make the name for the new generation file
        output_file_name = new_gen_folder_path + "generation_{}.smi".format(
            generation_num
        )
    else:
        # make the name for the new generation file
        output_file_name = new_gen_folder_path + "generation_{}{}.smi".format(
            generation_num, nomenclature_tag
        )

    # write as a tab delineated .smi file
    with open(output_file_name, "w") as f:
        for smile in formatted_smile_list:
            # smile_string = smile[0]
            # smile_id = smile[1]
            x = str(smile[0] + "\t" + str(smile[1]) + "\n")
            f.write(x)

    sys.stdout.flush()
    return output_file_name, new_gen_folder_path


def save_ligand_list(output_directory, generation_num,
                     list_of_chosen_ligands, nomenclature_tag):
    """
    Save the list of ligands. nomenclature_tag is a string such as "Mutation"
    or "Crossover" or "Previous_Gen_choice" describing what this data is used
    for. If it says seeding it is the chosen mols from the previous generation
    being used to seed the next generation

    Inputs:
    :param dict output_directory: the directory of the run to save the
        generation
    :param int generation_num: The generation number
    :param list list_of_chosen_ligands: The formatted list of ligands to seed
        a generation
    :param str nomenclature_tag: The str describing the ligand list
        -ie seeding_mutations is the list that seeded the mutations while
            mutations would be the list of mutations generated from the
            seeding_mutations list
        -ie. mutation, crossover, Previous_Gen_choice
    """
    # make a folder for the new generation
    new_gen_folder_path = output_directory + "generation_{}{}".format(
        generation_num, os.sep
    )

    # make a folder for the Seed files
    seed_folder_path = new_gen_folder_path + "SeedFolder" + os.sep

    # check if folders exist, if not make them
    if not os.path.isdir(new_gen_folder_path):
        os.makedirs(new_gen_folder_path)
    if not os.path.isdir(seed_folder_path):
        os.makedirs(seed_folder_path)

    output_file_name = "{}{}_Gen_{}.smi".format(
        seed_folder_path, nomenclature_tag, generation_num
    )

    # save to a new output smiles file. ie. save to ranked_smiles_file
    with open(output_file_name, "w") as output:
        for line in list_of_chosen_ligands:
            output_line = "\t".join(line) + "\n"
            output.write(output_line)

    sys.stdout.flush()
