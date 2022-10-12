# AutoGrow4

AutoGrow4 is an evolutionary algorithm that optimizes candidate ligands for
predicted binding affinity and other drug-like properties. Though no
substitute for the medicinal chemist, AutoGrow4 attempts to introduce some
chemical intuition into the automated optimization process.

AutoGrow4 takes advantage of recent advancements in multiprocessing,
cheminformatics, and docking. It handles all ligand manipulations using SMILES
strings via the RDKit API. Docking programs require 3D small-molecule models,
so AutoGrow4 use Gypsum-DL to convert the SMILES to 3D.

AutoGrow4 was designed with modularity in mind. New options allow users to
control drug-likeness, access new reaction libraries, use additional docking
programs, rank AutoGrow4-generated compounds using alternate schemes, and
advance candidate compounds via several selection algorithms. The codebase has
been redesigned so users can easily add their own custom options as well.
Please see the tutorial (`/autogrow4/tutorial/TUTORIAL.md`), which describes
how to expand these various options.

The default fragment libraries included with AutoGrow4 were derived from a
subset of the ZINC database (https://zinc.docking.org/). We thank ZINC for
allowing us to distribute these fragment libraries to AutoGrow4 users.

## Getting Started with AutoGrow4

Use the `RunAutogrow.py` script to run AutoGrow4. The script supports both
command-line parameters and a JSON parameter file. A detailed tutorial
describing AutoGrow4's dependencies and general use is located at
`/autogrow4/tutorial/TUTORIAL.md`. Additional descriptions of all AutoGrow4
parameters can be obtained by running: `python RunAutogrow.py -h`.

We strongly recommend running AutoGrow4 via Docker using
`/autogrow4/docker/autogrow_in_docker.py`. See the tutorial at
`/autogrow4/tutorial/TUTORIAL.md` for more details.

## Dependencies Notes

AutoGrow4 version 4.0.1 has been tested with the following dependencies:

```python
>>> rdkit.__version__
'2020.03.1'
>>> numpy.__version__
'1.18.1'
>>> scipy.__version__
'1.4.1'
>>> matplotlib.__version__
'3.2.1'
>>> func_timeout.__version__
'4.3.5'
```

If you are unable to run AutoGrow4, please try running AutoGrow4 in a python
environment with these specific dependencies. Alternatively, the Docker
version of AutoGrow4 automatically installs dependencies that are verified to
work with AutoGrow4. If you discover that AutoGrow4 is no longer compatible
with current library releases, please contact us, and we will attempt to
correct the code.

## Developer Note

Dimorphite is Version 1.2.3 with the citation print statements silenced.
Please remember to cite Dimorphite-DL: Ropp PJ, Kaminsky JC,
Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program for
enumerating the ionization states of drug-like small molecules. J Cheminform
11:14. doi:10.1186/s13321-019-0336-9.

As of May 2020, AutoGrow4 works in the provided Docker container (with the
specific, hard-coded versions of the dependencies). These dependencies are
known to be compatible with AutoGrow4. When developing future versions of
AutoGrow4, be sure to test the current versions of these dependencies and to
update `$PATH/docker/Dockerfile` appropriately.












## conda 

```bash
source activate autogrow 
```


## data 








## run 

see `tutorial/TUTORIAL.md`. 

```
python RunAutogrow.py \
    --filename_of_receptor ./tutorial/PARP/4r6eA_PARP1_prepared.pdb \
    --center_x -70.76 --center_y  21.82 --center_z 28.33 \
    --size_x 25.0 --size_y 16.0 --size_z 25.0 \
    --source_compound_file ./source_compounds/naphthalene_smiles.smi \
    --root_output_folder ./output \
    --number_of_mutants_first_generation 50 \
    --number_of_crossovers_first_generation 50 \
    --number_of_mutants 50 \
    --number_of_crossovers 50 \
    --top_mols_to_seed_next_generation 50 \
    --number_elitism_advance_from_previous_gen 50 \
    --number_elitism_advance_from_previous_gen_first_generation 10 \
    --diversity_mols_to_seed_first_generation 10 \
    --diversity_seed_depreciation_per_gen 10 \
    --num_generations 5 \
    --mgltools_directory ./mgltools_x86_64Linux2_1.5.6/ \
    --number_of_processors -1 \
    --scoring_choice VINA \
    --LipinskiLenientFilter \
    --start_a_new_run \
    --rxn_library click_chem_rxns \
    --selector_choice Rank_Selector \
    --dock_choice VinaDocking \
    --max_variants_per_compound 5 \
    --redock_elite_from_previous_gen False \
    --generate_plot True \
    --reduce_files_sizes True \
    --use_docked_source_compounds True  
```


- `./RunAutogrow.py`: ```AutogrowMainExecute.main_execute(vars)```
- `./autogrow/autogrow_main_execute.py`: ```operations.populate_generation```
- `./autogrow/operators/oprations.py`: ```execute_crossover.make_crossovers``` 







