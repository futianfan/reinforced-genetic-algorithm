# Welcome to AutoGrow4

This document will break down how to run AutoGrow4. It will also cover what
dependencies are required.

Please note that there are several paths used in this tutorial that need to be
replaced with the proper paths on one's own system. These paths will be
indicated by a string of ALL_CAPS. For example (from Section: Installing
AutoGrow):

`cd /PATH_TO/DESIRED_DIR/`

You must replace `/PATH_TO/DESIRED_DIR/` with the path to the `autogrow4`
directory on your own system. On a Ubuntu OS, this may look like:

`cd /home/jacob/Documents/`

For brevity, we simplify the main `autogrow4` directory to `/autogrow4/`
throughout this tutorial. You may need to supply the `/PATH_TO/DESIRED_DIR/`
described above before `/autogrow4/`.

## Computer Requirements

AutoGrow4 has been tested on Ubuntu 16.04 and higher, as well as MacOS 10.13
High Sierra. It has been verified to work on an HPC cluster using SMP
multithreading (RedHat Enterprise Server release 7.3 Maipo).

AutoGrow4 has not been configured for Windows OS, but a script capable of
running AutoGrow4 within a docker container on Windows can be found:

`/autogrow4/docker/autogrow_in_docker.py`

This script should run on any docker-enabled machine, and should be capable of
multithreading. Details on running AutoGrow4 within a docker container can be
found below, in the Section: Docker submission.

## Installing AutoGrow4

Download a copy of AutoGrow4 from the Durrant-lab website at
[http://git.durrantlab.com/jdurrant/autogrow4](http://git.durrantlab.com/jdurrant/autogrow4).

You can also install AutoGrow4 using the git clone command:

```bash
cd /PATH_TO/DESIRED_DIR/
git clone https://git.durrantlab.pitt.edu/jdurrant/autogrow4
```

## Dependencies

AutoGrow4 has several dependencies that may need to be installed separately.

### Bash (Required)

A modern installation of bash is required to run AutoGrow4. AutoGrow4 has been
tested using GNU bash, version 4.4.19. macOS and Linux come with Bash
preinstalled.

### Coreutils (Required For macOS)

Most Linux operating systems include the `timeout` tool (part of the
`coreutils` package) that AutoGrow4 requires. Use on macOS requires the
separate installation of the `coreutils` package, available through
`homebrew`, which provides the equivalent `gtimeout` binary.

```bash
sudo brew install coreutils
```

### Python Installation (Required)

AutoGrow4 is primarily written in python. A modern version of python can be
installed using `conda`:

- [https://docs.conda.io/projects/conda/en/latest/user-guide/install/](https://docs.conda.io/projects/conda/en/latest/user-guide/install/),
  or
- [http://www.python.org/getit/](http://www.python.org/getit/).

AutoGrow4 has been tested with python 2.7, 3.6, and 3.7. Future support and
updates will focus on 3.7. We recommend using the most current version of
python available, 3.7 or newer.

### MGLTools

MGLTools is written by the creators of Autodock Vina. It is used by AutoGrow4
to convert .pdb files to the .pdbqt format. The .pdbqt format is required by
Vina-type docking programs, including Autodock Vina and QuickVina2.

Morris, G. M., Huey, R., Lindstrom, W., Sanner, M. F., Belew, R. K., Goodsell,
D. S. and Olson, A. J. (2009) Autodock4 and AutoDockTools4: automated docking
with selective receptor flexibility. J. Computational Chemistry 2009, 16:
2785-91

If you prefer not to use MGLTools to perform this conversion, you may also use
Open Babel (`obabel`) or custom file-converting/docking software.

#### Installation

WARNING: MGLTools installation can be tricky! We recommend that you DO NOT
`pip` or `conda` install MGLTools, as those package managers provide an
outdated python package and creates issues with environments.

The best way to install MGLTools is to download the latest release of the
command-line version (NOT THE GUI VERSION) from
[http://mgltools.scripps.edu/downloads](http://mgltools.scripps.edu/downloads)

Once the command-line version of the MGLTools package has been downloaded,
follow this example installation (Linux, MGLTools 1.5.6):

1. To extract files, unzip/untar the package: `tar -xvf
   /PATH_TO/mgltools_x86_64Linux2_1.5.6.tar.gz`
2. Go to the extract folder: `cd  /PATH_TO/mgltools_x86_64Linux2_1.5.6`
3. Run the installation script and make sure MGLToolsPckgs is unpacked:
   - If `/PATH_TO/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/` is a folder:
     `bash install.sh`
   - If `/PATH_TO/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/` is not a folder,
     you must manually unzip/untar `MGLToolsPckgs.tar.gz`: `tar -xvf
     /PATH_TO/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs.tar.gz`
4. Click 'OK' to the licensing agreement. Please note MGLTools is free for
   academic use but may require a license for commercial usage. This should
   open automatically.
5. Find the proper pathing for the AutoGrow4 variable (see next section,
   Additional Pathing Instructions).

#### Additional Pathing Instructions

To use MGLTools to convert files, AutoGrow4 must know the path to the MGLTools
directory. The path can be found by:

1. Going to the extract folder: `cd /PATH_TO/mgltools_x86_64Linux2_1.5.6`
2. Using the `pwd` command in bash to get the absolute path to the MGLTools
   directory.

When using AutoGrow, specify the path to this MGLTools directory using the
`--mgltools_directory` parameter: `python RunAutogrow.py ..
--mgltools_directory /PATH_TO/mgltools_x86_64Linux2_1.5.6 ...`

On Linux and macOS machines, AutoGrow4 will auto-locate three important file
paths based on `mgltools_directory`:

1. prepare_ligand4.py: `mgltools_directory` +
   `/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py`
2. prepare_receptor4.py: `mgltools_directory` +
   `/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py`
3. mgl_python: `mgltools_directory` + `/bin/pythonsh`

If running on Windows OS (limited support) please provide those paths to
AutoGrow4 explicitly:

```bash
python RunAutogrow.py .. \
    --mgltools_directory /PATH_TO/mgltools_win32_1.5.6\ \
    --prepare_ligand4.py /PATH_TO/mgltools_win32_1.5.6\ \
    --prepare_receptor4.py /PATH_TO/mgltools_win32_1.5.6\ \
    --mgl_python /PATH_TO/mgltools_win32_1.5.6\ + /PATH_TO/pythonsh ...
```

Specifying the `--prepare_ligand4.py`, `--prepare_receptor4.py`, and
`--mgl_python` parameters explicitly will override the paths that AutoGrow
determines using the `--mgltools_directory` parameter alone.

### obabel

AutoDock can also perform the .pdb to .pdbqt conversion using Open Babel if
the user prefers not to use MGLTools or custom file-conversion/docking
software.

N M O'Boyle, M Banck, C A James, C Morley, T Vandermeersch, and G R Hutchison.
"Open Babel: An open chemical toolbox." J. Cheminf. (2011), 3, 33.
DOI:10.1186/1758-2946-3-33

Open Babel includes the `obabel` command-line tool for cheminformatic file
conversion. It is used by AutoGrow4 to convert .pdb files to .pdbqt format, as
required for docking with Vina-type programs such as Autodock Vina and
QuickVina2.

#### Installation

An easy installation on Linux/macOS machines is:

```bash
sudo apt-get install openbabel
sudo apt-get update
```

Full instructions for `obabel` installation can be found on their site:
https://openbabel.org/docs/dev/Installation/install.html

AutoGrow4 has been tested with `obabel` version 2.4.1

#### Additional Pathing Instructions

To use `obabel` to convert files, AutoGrow4 requires the path to the `obabel`
executable. Once Open Babel is installed, the path to the `obabel` executable
can be found by running `which obabel`

This path should be provided to AutoGrow4 using the `--obabel_path` variable:

```bash
python RunAutogrow.py .. --obabel_path /PATH_TO/obabel ...
```

### Python APIs (Required)

AutoGrow4 also uses several python API libraries beyond those found in the
standard library. These must be installed separately. Most can be installed
via `conda` or `pip`.

#### Mandatory Installations

RDKit: Cheminformatic library. RDKit can be downloaded via `conda`/`pip`. To
install using `conda` use the command:

```bash
conda install -c rdkit rdkit
```

We use the following RDKit sub-libraries in AutoGrow4:

```python
import rdkit
from rdkit import RDLogger, Chem, DataStructs
from rdkit.Chem import MolSurf, Crippen, rdFMCS, Descriptors, AllChem, FilterCatalog,  Lipinski, rdDepictor
from rdkit.Chem.Draw import PrepareMolForDrawing, rdMolDraw2D
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem.rdchem import BondStereo
```

NumPy (mathematical functions) can be downloaded via `conda`/`pip`. It can be
`conda` installed using the command `conda install -c anaconda numpy`.
AutoGrow4 has been tested using `numpy` version 1.15.0.

SciPy (mathematical functions) can be downloaded via `conda`/`pip`. It can be
`conda` installed using the command `conda install -c anaconda scipy`.
AutoGrow4 has been tested using `scipy` version 1.1.0.

Matplotlib (python graphing tool) can be downloaded via `conda`/`pip`. It can
be `conda` installed using the command `conda install matplotlib`. AutoGrow4
has been tested using `matplotlib` version 3.0.2.

func_timeout (pythonic timeout tool) can be downloaded via `pip`. It can be
`pip` installed using the command `pip install func-timeout`. AutoGrow4 has
been tested using `func_timeout` version 4.3.5.

#### Optional Installations

mpi4py (MPI multithreading python library) is required for MPI multithreading.
It can be downloaded via `conda`/`pip`. It can be `conda` installed using the
command `conda install -c anaconda mpi4py`. AutoGrow4 has been tested using
`mpi4py` version 3.0.1. This may require a preinstallation of `mpich`: `sudo
apt install mpich`

AutoGrow4 requires `mpi4py` version 2.1.0 and higher. To check the version:

1. open a python window.
2. enter into the window:

```python
>>> import mpi4py
>>> mpi4py.__version__
    3.0.1
```

MPI mode also requires an MPI-enabled computer environment. The authors use
OpenMPI. OpenMPI installation instructions can be found:
http://lsi.ugr.es/jmantas/pdp/ayuda/datos/instalaciones/Install_OpenMPI_en.pdf

Quick OpenMPI installation is possible in a bash terminal:

```bash
sudo apt-get install openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.3 libopenmpi-dbg libopenmpi-dev
```

Establishing a fully MPI-enabled computer network is complicated and should
only be attempted by qualified technicians. The authors used an Intel’s
Omni-Path communication architecture that was established by experts at the
University of Pittsburgh’s Center for Research Computing. The authors DO NOT
RECOMMEND ATTEMPTING THIS ON YOUR OWN.

### Pre-Installed Python and Binary Dependencies

AutoGrow4 comes with several dependencies preinstalled, requiring no
additional effort by the user. These packages have licenses that allow them to
be freely redistributed. If a dependency updates, please feel free to contact
us, and we will do our best to make our code future-compatible.

#### Docking Programs

AutoGrow4 comes preinstalled with two docking programs:

- Autodock Vina 1.1.2 (packaged with executables for Linux, macOS, and
  Windows)
  - Version: 1.1.2
  - Location: `/autogrow4/autogrow/docking/docking_executables/vina/`
  - Citation: Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the
    speed and accuracy of docking with a new scoring function, efficient
    optimization, and multithreading. Journal of computational chemistry,
    31(2), 455–461. doi:10.1002/jcc.21334
  - License: Apache version 2

- QuickVina2.1 (compatible with Linux OS and macOS)
  - Version: 2.1
  - Location: `/autogrow4/autogrow/docking/docking_executables/q_vina_2/`
  - Citation: Amr Alhossary, Stephanus Daniel Handoko, Yuguang Mu, and
    Chee-Keong Kwoh. Bioinformatics (2015) 31 (13): 2214-2216.
    [DOI:10.1093/bioinformatics/btv082](https://doi.org/10.1093/bioinformatics/btv082)
  - License: Apache version 2

These programs can be found within the directory
`/autogrow4/autogrow/docking/docking_executables/`

AutoGrow4 allows users to provide custom docking software. This could be as
simple as using a different version of Autodock Vina:

```bash
python RunAutogrow.py ... --docking_executable /PATH_TO/Autodock_Vina_version_X_executable
```

More advanced use allows users to provide a custom docking program. Details
re. how to use custom docking suites are provided below in the section
"Providing Custom Options".

#### Scoring/Rescoring Programs

NNScore 1 and NNScore 2 are free and open-source programs that are distributed
with AutoGrow4. Both NNScore1 and NNScore2 reassess ligand docking. They were
trained using Autodock Vina 1.1.2, so they can only be used when Vina 1.1.2 is
the user-specified AutoGrow4 docking engine.

AutoGrow4 allows users to provide custom Scoring/Rescoring software. Details
for custom Scoring/Rescoring suites are provided below in the section
"Providing Custom Options."

- NNScore 1:
  - Version: 1.1
  - Location: `/autogrow4/autogrow/docking/scoring/nn_score_exe/nnscore1/`
  - Citation: NScore: A Neural-Network-Based Scoring Function for the
    Characterization of Protein-Ligand Complexes. Jacob D. Durrant, J. Andrew
    McCammon. Journal of Chemical Information and Modeling, 2010, 50 (10),
    pp865-1871.
  - License: GNU General Public version 3

- NNScore 2:
  - Version: 2.02
  - Location: `/autogrow4/autogrow/docking/scoring/nn_score_exe/nnscore2/`
  - Citation: NNScore 2.0: A Neural-Network Receptor–Ligand Scoring Function.
    Jacob D. Durrant, Andrew McCammon. Journal of Chemical Information and
    Modeling, 2011, 51 (11), pp 2897-2903.
  - License: GNU General Public version 3

#### SMILES Conversion to 3D and Protonation Adjustments

AutoGrow4 performs most of its ligand handling using 2D SMILES. AutoGrow4 uses
the free and open-source program Gypsum-DL to convert from SMILES to 3D SDF
format. Gypsum-DL is prepackaged in AutoGrow4. Gypsum-DL itself also includes
the MolVS and Dimorphite-DL packages.

- Gypsum-DL:
  - Version: 1.1.2
  - Location: `/autogrow4/autogrow/operators/convert_files/gypsum_dl/`
  - Citation: Ropp PJ, Spiegel JO, Walker JL, Green H, Morales GA, Milliken
    KA, Ringe JJ, Durrant JD. Gypsum-DL: An Open-Source Program for Preparing
    Small-Molecule Libraries for Structure-Based Virtual Screening. J
    Cheminform. 11(1):34, 2019. [PMID: 31127411] [doi:
    10.1186/s13321-019-0358-3]
  - License: Apache version 2.0

- Dimorphite-DL:
  - Version: 1.2.2
  - Location:
    `/autogrow4/autogrow/operators/convert_files/gypsum_dl/gypsum_dl/Steps/SMILES/dimorphite_dl`
  - Citation: Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019)
    Dimorphite-DL: An open-source program for enumerating the ionization
    states of drug-like small molecules. J Cheminform 11:14.
    doi:10.1186/s13321-019-0336-9.
  - License: Apache version 2.0

- MolVS:
  - Version: v0.1.1 2019 release
  - Location:
    `/autogrow4/autogrow/operators/convert_files/gypsum_dl/gypsum_dl/molvs`
  - Citation: https://molvs.readthedocs.io; Take from
    https://github.com/mcs07/MolVS
  - License: MIT License

## Running AutoGrow4

To run AutoGrow4, use the python script `RunAutogrow.py`, located in the top
AutoGrow4 directory, from the command line. AutoGrow4 accepts user input via
two methods:

1. Command-line submission: executing directly from the command line.

```bash
cd /PATH_TO/autogrow4/

python RunAutogrow.py \
    --filename_of_receptor /autogrow4/autogrow/tutorial/PARP/4r6eA_PARP1_prepared.pdb \
    --center_x -70.76 --center_y  21.82 --center_z 28.33 \
    --size_x 25.0 --size_y 16.0 --size_z 25.0 \
    --source_compound_file /autogrow4/autogrow/source_compounds/naphthalene_smiles.smi \
    --root_output_folder /PATH_TO/output_directory/ \
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
    --mgltools_directory /PATH_TO/mgltools_x86_64Linux2_1.5.6/ \
    --number_of_processors -1 \
    --scoring_choice VINA \
    --LipinskiLenientFilter \
    --start_a_new_run \
    --rxn_library ClickChem \
    --selector_choice Rank_Selector \
    --dock_choice VinaDocking \
    --max_variants_per_compound 5 \
    --redock_elite_from_previous_gen False \
    --generate_plot True \
    --reduce_files_sizes True \
    --use_docked_source_compounds True \
    >  /PATH_TO/OUTPUT/text_file.txt 2>  /PATH_TO/OUTPUT/text_errormessage_file.txt
```

2. json file submission: store AutoGrow4 parameters in a .json file

```bash
cd /PATH_TO/autogrow4/
python RunAutogrow.py -j /PATH_TO/json_file_with_variable.json
```

Examples of the json files can be found in the folder
`/autogrow4/sample_sub_scripts/`.

## Understanding AutoGrow4 Parameters

An explanation of every parameter can be retrieved by running:

```bash
python /autogrow4/RunAutogrow.py --help
```

Custom options such as custom filters, docking software, reaction libraries,
etc., are described in other parts of the tutorial. Additionally, a
description of how to prepare the receptor file is also provided in the
section "Prepping Receptor". Details for preparing source compound files are
provided directly below.

### Source Compound Files

Source compound files simply tab-delineated SMILES files (.SMI). Specify the
path using the parameter `--source_compound_file`.

Examples of source compound files can be found at
`/autogrow4/source_compounds/`

A detail log of how the examples files were prepared is located at
`/autogrow4/source_compounds/Example_source_compound_notes.txt`

An accessory script that converts a folder of PDB files to a tab-delineated
.SMI file is provided at
`/autogrow4/accessory_scripts/convert_directory_ligands_pdb_to_smi.py`

Details for using this accessory script are provided near the bottom of this
document, in the section
"/autogrow4/accessory_scripts/convert_directory_ligands_pdb_to_smi.py".

#### Seeding an AutoGrow4 Run with Already Docked Compounds

AutoGrow4 can be set to assess the fitness of source compounds, in addition to
any new evolving molecules. This will dock the source compounds prior to
seeding the first generation. It creates a generation 0, consisting only of
the compounds in the source compound file.

To use this option, set `use_docked_source_compounds` as `True`.

If one is running multiple independent runs using the same source compounds,
it may be worth it to test all source compounds first and seed all runs with
the same scores for generation 0. This is accomplished by providing a source
compound file with a float value in the second to last tab-delineated slot.
Such a SMI file would be formatted like this:

```text
SMILE   Name    ...any info...  primary_fitness diversity_fitness
```

Please remember that docking scores are relative to a specific protein and
pocket. Changing the coordinates, protein, docking software, or (re)scoring
method will invalidate any information.

More information is provided in the `use_docked_source_compounds` section of
the `RunAutoGrow.py` help menu, which can be displayed by running:

```bash python RunAutoGrow.py --help```

## Docker Submission

The `/autogrow4/docker/` directory contains the scripts to run AutoGrow4
within a docker container. These scripts are useful when using an OS that is
not compatible with AutoGrow4 or its dependencies, such as Windows.

Prior to running these scripts, please install the docker software. Also, be
sure to ***also always run these scripts with sudo (linux/macOS) or
administrator privileges (Windows).***.

Running AutoGrow4 via docker the first time will take a few minutes longer
because docker must install the dependencies. The same is true if docker
images have been purged.

Depending on the AutoGrow4 settings, processor speed/count, etc., AutoGrow4
may complete within minutes or may take as long as multiple days. Please make
sure to use settings that are appropriate for your system. Using `nohup` may
be a useful wrapper for longer runs or when running jobs remotely (i.e., over
ssh).

More details are provided directly below and in the
`/autogrow4/docker/README.md` section.

### How to Setup AutoGrow4 in Docker

Dockerized AutoGrow4 requires the user to specify parameters via a JSON file
(not the command line).

To run the `autogrow_in_docker.py` script:

Linux/MacOS:

1. Change into the `/autogrow4/docker/` directory in a bash terminal: `cd
   /autogrow4/docker/`
2. Run `autogrow_in_docker.py` with `sudo` and supply a json file using the
   normal pathing of your system. Please note that the docker downloads its
   own copy of `obabel` and `MGLTools`, so you do not need to specify those
   paths.
3. Execute `autogrow_in_docker.py` with `sudo` privileges, providing it with a
   JSON file (MUST EXECUTE FROM `/autogrow4/docker/`): `sudo python
   autogrow_in_docker.py -j ./examples/sample_autogrow_docker_json.json`
4. Results will appear in the output directory specified by the
   `--root_output_folder` parameter.

Windows OS:

1. Open a docker-enabled and bash-enabled terminal with administrator
   privileges.
2. Change into the `/autogrow4/docker/` directory in a bash enabled terminal:
   `cd /autogrow4/docker/`
3. Execute `autogrow_in_docker.py` with `sudo` privileges, providing it with a
   JSON file (MUST EXECUTE FROM `/autogrow4/docker/`): `python
   autogrow_in_docker.py -j ./examples/sample_autogrow_docker_json.json`
4. Results will appear in the output directory specified by the
   `--root_output_folder` parameter.

## Providing Custom Plugins

AutoGrow4 was designed to be modular. This allows for the easy swapping of
code. AutoGrow4 is intended to be a living codebase. If you have added good
custom code and would like to make it open source, please contact the authors
so that we can grow the user options.

Many of the AutoGrow4 functions can be supplemented with custom options. These
functions include:

1. Custom Ligand Filters ***
2. Custom Docking Code ***
3. Custom ligand conversion code from PDB to dockable format (ie PDBQT) ***
4. Custom Scoring/Rescoring code ***
5. Custom Reaction libraries
6. Custom complementary molecules libraries

*** Indicates that when using this feature, the code is automatically copied
into the appropriate AutoGrow4 directory. This is only done once, so please
unittest the code prior to incorporating it into AutoGrow4. A print message
will indicate where the file has been copied. That file can be manually
deleted or overwritten by the user. Restart AutoGrow4 after the custom files
have been automatically copied into the proper locations. After that the new
script should be integrated into AutoGrow4.

AutoGrow4 ASSUMES ALL CUSTOM CODE HAS BEEN TESTED AND FUNCTIONS WITH SPECIFIED
I/O. For example, it assumes that scoring favors the most negative docking
score.

- AutoGrow4 will continue to assume all custom Scoring scripts set the most
  fit score to the most negative for all metrics besides diversity.
- It also assumes in ranked .smi files that the last column is the diversity
  fitness and assumes the second to last column is the metric for
  "docking/rescored" fitness.
- If a custom script scores ligands such that the most fit ligand has the
  highest score, AutoGrow4 may inadvertently be favoring ligands that are
  least fit. In some cases, users may need to multiple custom docking scores
  by -1 so that the best ligands have the most negative scores.

### 1. Custom Ligand Filters ***

This feature allows the user to incorporate custom python scripts for
filtering ligands. These filters are applied to ligands after they are created
by mutation/crossover but before Gypsum-DL conversion to 3D.

This custom code will be copied to the directory:
`/autogrow4/autogrow/operators/filter/filter_classes/filter_children_classes/`

#### Script Formatting

These filters use a class-based inheritance architecture with filter classes
that must

1. Inherit ParentFilterClass located at
   `/autogrow4/autogrow/operators/filter/filter_classes/parent_filter_class.py`
2. Have a unique name: `class unique_name(ParentFilter)` (`unique_name` cannot
   match one of the predefined filters)
3. Have at least one function called `run_filter` (`run_filter` takes a single
   variable which must be an rdkit molecule object).

#### Running Custom Filters

Because parameters can be specified via command-line or JSON file, we provide
an example of each when submitting custom filters.

1. Command-line submission:
   - Custom file is located at `/PATH_TO/custom_filter_1.py`
   - Unique class name is `custom_filter_1` (this will be how it is called in
     future versions)
   - To run multiple custom filters replace
       `[["custom_filter_1","/PATH_TO/custom_filter_1.py"]]` with:
       `[["custom_filter_1","/PATH_TO/custom_filter_1.py"],["custom_filter_2","/PATH_TO/custom_filter_2.py"]]`

```bash
python RunAutogrow.py \
    ... \
    --alternative_filter [["custom_filter_1","/PATH_TO/custom_filter_1.py"]]
```

2. JSON file submission:
    - Where the custom file is located at `/PATH_TO/custom_filter_1.py`
    - Unique class name is `custom_filter_1` (this will be how it is called in
      future submissions)
    - To run multiple files  Replace
      `[["custom_filter_1","/PATH_TO/custom_filter_1.py"]]` with:
      `[["custom_filter_1","/PATH_TO/custom_filter_1.py"],["custom_filter_2","/PATH_TO/custom_filter_2.py"]]`

```json
{
    ... ,
    "alternative_filter": [["custom_filter_1","/PATH_TO/custom_filter_1.py"]]
}
```

Submit in terminal: `python RunAutogrow.py -j
/PATH_TO/json_file_with_variable.json`

### 2. Custom Docking Code ***

This feature allows the user to incorporate custom python scripts for docking
ligands. Currently AutoGrow4 is configured to dock using Autodock Vina and
QuickVina2, but AutoGrow4 is not limited to these docking programs. A custom
script can be added to run docking using virtually any software.

This custom code will be copied to the directory:
`/autogrow4/autogrow/docking/docking_class/docking_class_children/`

#### Script Formatting

These docking scripts use a class-based inheritance architecture which
require:

1. Docking class must inherit ParentDocking:
   `/autogrow4/autogrow/docking/docking_class/parent_dock_class.py`
2. Must have a unique name: `class unique_name(ParentDocking)` (`unique_name`
   cannot be one of the predefined docking scripts, currently just VinaDocking
   and QuickVina2Docking)
3. Must have at least have three functions following the below formatting:

```python
def __init__(self, vars=None, receptor_file=None, test_boot=True):
    """
    get the specifications for ligand assessment/docking from vars
    load them into the self variables we will need
    and convert the receptor to the proper file format (ie pdb-> pdbqt)

    Inputs:
    :param dict vars: Dictionary of User variables
    :param str receptor_file: the path for the receptor pdb
    :param bool test_boot: used to initialize class without objects for testing purpose
    """

def run_dock(self, pdbqt_filename):
    """
    this function runs the docking. Returns None if it worked and the name if it failed to dock.

    Inputs:
    :param str pdbqt_filename: the pdbqt file of a ligand to dock and score
                if using a docking software that use a file format other than pdbqt please substitute that file here
    Returns:
    :returns: str smile_name: name of smiles if it failed to dock
                            returns None if it docked properly
    """

def rank_and_save_output_smi(self, vars, current_generation_dir, current_gen_int, smile_file, deleted_smiles_names_list):
    """
    Given a folder with PDBQT's, rank all the SMILES based on docking score (High to low).
    Then format it into a .smi file.
    Then save the file.

    Inputs:
    :param dict vars: vars needs to be threaded here because it has the paralizer object which is needed within Scoring.run_scoring_common
    :param str current_generation_dir: path of directory of current generation
    :param int current_gen_int: the interger of the current generation indexed to zero
    :param str smile_file:  File path for the file with the ligands for the generation which will be a .smi file
    :param list deleted_smiles_names_list: list of SMILES which may have failed the conversion process

    Returns:
    :returns: str output_ranked_smile_file: the path of the output ranked .smi file
    """
```

#### Running Custom Docking scripts

Please note integrating a new docking engine into AutoGrow4 will likely
require corresponding custom conversion and scoring scripts. Documentation for
these is provided in the next two subsections. The example below ignores these
extras.

AutoGrow4 will need to be restarted after this has been incorporated into the
code base.

Submission through .json format:

- Where JSON is located at: `/PATH_TO/To/json_file_with_variable.json`
- Where docking software executable is located at:
  `/PATH_TO/EXECUTABLE_FOR_CUSTOM_DOCKING/custom_docking`
- Where python script for running docking is located at:
  `/PATH_TO/CLASS_OBJECT_FOR_CUSTOM_DOCKING/custom_docking.py`
- Where name of custom docking class is: custom_docking

```json
{
    ...
    "docking_executable": "/PATH_TO/EXECUTABLE_FOR_CUSTOM_DOCKING/custom_docking",
    "dock_choice": "Custom",
    "custom_docking_script": ["custom_docking", "/PATH_TO/CLASS_OBJECT_FOR_CUSTOM_DOCKING/custom_docking.py"]
}
```

Submit via terminal: `python RunAutogrow.py -j
/PATH_TO/To/json_file_with_variable.json`

Command-line submission format:

- Where docking software executable is located at:
  `/PATH_TO/EXECUTABLE_FOR_CUSTOM_DOCKING/custom_docking`
- Where python script for running docking is located at:
  `/PATH_TO/CLASS_OBJECT_FOR_CUSTOM_DOCKING/custom_docking.py`
- Where name of custom docking class is: `custom_docking`

```bash
python RunAutogrow.py \
    ... \
    --docking_executable "/PATH_TO/EXECUTABLE_FOR_CUSTOM_DOCKING/custom_docking" \
    --dock_choice Custom \
    --alternative_filter ["custom_docking", "/PATH_TO/CLASS_OBJECT_FOR_CUSTOM_DOCKING/custom_docking.py"]
```

### 3. Custom Ligand Conversion Code From PDB To Dockable Format (i.e., PDBQT)

If using a docking software other than VINA/QuickVina, you may need to convert
the pdb-formatted ligands into a different format. In this case, you must
provide a custom script to convert the ligands.

This custom code will be copied to the directory:
`/autogrow4/autogrow/docking/docking_class/docking_class_children/`

#### Script Formatting

These conversion scripts use a class-based inheritance architecture:

1. Conversion class object must inherit ParentPDBQTConverter:
   `/autogrow4/autogrow/docking/docking_class/parent_pdbqt_converter.py`
2. Have a unique name: `class unique_name(ParentPDBQTConverter)`
   (`unique_name` can not be one of the predefined docking scripts)
   - Currently files named: `convert_with_mgltools.py` and
    `convert_with_obabel.py`
   - Class names already in use are: `MGLToolsConversion` and `ObabelConversion`
3. Must have at least two functions following the below formatting:

```python
def convert_receptor_pdb_files_to_pdbqt(self, receptor_file, mgl_python, receptor_template, number_of_processors):
    """
    Make sure a PDB file is properly formatted for conversion to pdbqt

    Inputs:
    :param str receptor_file:  the file path of the receptor
    :param str mgl_python: file path of the pythonsh file of mgl tools
    :param str receptor_template: the receptor4.py file path from mgl tools.
    :param int number_of_processors: number of processors to multithread
    """
    raise NotImplementedError("convert_receptor_pdb_files_to_pdbqt() not implemented")

def convert_ligand_pdb_file_to_pdbqt(self, pdb_file):
    """
    Convert the ligands of a given directory from pdb to pdbqt format

    Inputs:
    :param str pdb_file: the file name, a string.
    Returns:
    :returns: bool bool: True if it worked;
                        False if its the gypsum param file or if it failed to make PDBQT
    :returns: str smile_name: name of the SMILES string from a pdb file
                                None if its the param file
    """
    raise NotImplementedError("rank_and_save_output_smi() not implemented")
```

#### Running Custom Conversion Scripts

AutoGrow4 will need to be restarted once after this has been incorporated into
the code base.

Submission through .json format:

- Where JSON is located at: `/PATH_TO/To/json_file_with_variable.json`
- Where python conversionscript is located at:
  `/PATH_TO/CLASS_OBJECT_FOR/custom_conversion.py`
- Where name of custom conversion class is: `custom_conversion`

```json
{
    ...
    "conversion_choice": "Custom",
    "custom_conversion_script": ["custom_conversion", "/PATH_TO/CLASS_OBJECT_FOR/custom_conversion.py"]
}
```

Submit via terminal: `python RunAutogrow.py -j
/PATH_TO/JSON_FILE/json_file_with_variable.json`

Command-line submission format:

- Where python conversionscript is located at:
  `/PATH_TO/CLASS_OBJECT_FOR/custom_conversion.py`
- Where name of custom conversion class is: `custom_conversion`

```bash
python RunAutogrow.py \
    ... \
    --conversion_choice Custom \
    --custom_conversion_script ["custom_conversion", "/PATH_TO/CLASS_OBJECT_FOR/for/custom_conversion.py"]
```

### 4. Custom Scoring/Rescoring Code

This feature allows the user to incorporate custom python scripts for scoring
and rescoring ligands.

Currently AutoGrow4 is configured to dock using Autodock Vina and QuickVina2.
There are also two options to rescore a ligand using either NNScore 1 or
NNScore 2. Additionally, ligand efficiency (dividing the score/rescore value
by the number of non-hydrogen atoms) can be applied with any float-based
scoring value.

Users can incorporate custom scoring and rescoring options into AutoGrow4.

This custom code will be copied to the directory:
`/autogrow4/autogrow/docking/scoring/scoring_classes/`

#### Script Formatting

These (re)scoring scripts use a class-based inheritance architecture:

1. Scoring class object must inherit parent_scoring_class:
   `/autogrow4/autogrow/docking/scoring/scoring_classes/parent_scoring_class.py`
2. Have a unique name: `class unique_name(parent_scoring_class)`
   - `unique_name` can not be one of the predefined docking scripts.
   - Currently files named: `vina.py`, `nn1.py`, `nn2.py`, and
     `lig_efficiency.py`
   - Class names already in use are: `VINA`, `NN1`, `NN2`, and `LigEfficiency`
3. Must have at least have two functions following the below formatting:

```python
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
```

#### Running Custom Scoring/Rescoring Scripts

AutoGrow4 will need to be restarted once after this has been incorporated into
the code base.

Submission through .json format:

- Where JSON is located at: `/PATH_TO/To/json_file_with_variable.json`
- Where python scoring script is located at:
  `/PATH_TO/CLASS_OBJECT_FOR/custom_scoring.py`
- Where name of custom scoring class is: `custom_scoring_name`

```json
{
    ...
    "scoring_choice": "Custom",
    "custom_scoring_script": ["custom_scoring_name", "/PATH_TO/CLASS_OBJECT_FOR/custom_scoring.py"]
}
```

Submit via terminal: `python RunAutogrow.py -j
/PATH_TO/JSON_FILE/json_file_with_variable.json`

Command-line submission format:

- Where python scoring script is located at:
  `/PATH_TO/CLASS_OBJECT_FOR/custom_scoring.py`
- Where name of custom scoring class is: `custom_scoring_name`

```bash
python RunAutogrow.py \
    ... \
    --conversion_choice Custom \
    --custom_conversion_script ["custom_scoring_name", "/PATH_TO/CLASS_OBJECT_FOR/custom_scoring.py"]
```

### 5. Custom Reaction libraries

AutoGrow4 assumes all custom scripts have been unittested. Please ensure all
reactions and libraries are accurate before using this option.

Unlike the other custom options, reaction libraries are stored as
human-readable json dictionaries. In contrast, all other custom options use
inherited class scripts. These json files do not need to be incorporated into
AutoGrow4 and thus require no restarting or copying of files.

Reaction Libraries are stored as .json files and are dictionaries of
dictionaries. The outer dictionary uses the reaction's name as the key and the
sub-dictionary containing all information about the reaction as the item.

We provide a script to check complementary and Reaction libraries molecule
libraries at: `/autogrow4/accessory_scripts/test_complementary_mol_library.py`

A tutorial is provided in the Accessory Scripts section of this document

#### Three Requirements For Custom Reaction Libraries

Custom reaction libraries require three pieces of information, each explained
below.

##### Reaction Library .json File Contains Reactions and All Reaction Information

Each sub-dictionary must contain the following information:

- "reaction_name": "Name of the reaction",
- "example_rxn_product": "SMILES of Product using example
  example_rxn_reactants",
- "example_rxn_reactants": ["SMILES of example reactant_1"],
  - if two or more reactants in reaction ["SMILES of example
    reactant_1","SMILES of example reactant_2",...]
- "functional_groups": ["functional group name reactant_1"],
  - if two or more reactants in reaction ["functional group name
    reactant_1","functional group name reactant_2",...]
- "group_smarts": ["functional_group SMARTS reactant_1"],
  - if two or more reactants in reaction ["functional_group SMARTS
    reactant_1","functional_group SMARTS reactant_2",...]
- "num_reactants": 1,
  - (int) if 2 or more reactants change accordingly
- "reaction_string": "reaction string ie
  reactant_1_smart.reactant_2_smart>>product_SMART",
  - This uses Daylights SMARTS reaction notation
- "RXN_NUM": 3
  - (int) a unique reaction number. This is used in naming products of
    mutation. For example, a ligand named Gen_1_Mutant_72_867328 is a ligand
    from generation 1 created by the 72 reaction in a reaction library

Simplified Example of a Reaction library (from click_chem_rxns_library.json):

```json
{
    "1_Epoxide_and_Alcohol":     {
        "reaction_name": "1_Epoxide_and_Alcohol",
        "example_rxn_product": "CC(C)(C)C(O)(OCCCF)C(C)(C)O",
        "example_rxn_reactants": ["CC(C)(C)C1(O)OC1(C)C", "FCCC(O)"],
        "functional_groups": ["Epoxide_clickchem", "Alcohol_clickchem"],
        "group_smarts": ["[CR1;H2,H1X4,H0X4]1O[CR1;H2,H1X4,H0X4]1", "[#6&$([CR0,R1X3,R1X4])&!$([#6](=,-[OR0,SR0])[OR0])]-[OR0;H1,-]"],
        "num_reactants": 2,
        "reaction_string": "[CR1;H2,H1X4,H0X4:1]1O[CR1;H2,H1X4,H0X4:2]1.[#6&$([CR0,R1X3,R1X4])&!$([#6](=,-[OR0,SR0])[OR0]):3]-[OR0;H1,-]>>O[C:1][C:2]O-[#6:3]",
        "RXN_NUM": 1
        },
    "2_Epoxide_and_Thiol":     {
        "reaction_name": "2_Epoxide_and_Thiol",
        "example_rxn_product": "CC(C)(C)C(O)(SCCC(=O)OC(=O)[O-])C(C)(C)O",
        "example_rxn_reactants": ["CC(C)(C)C1(O)OC1(C)C", "O=C([O-])OC(=O)CCS"],
        "functional_groups": ["Epoxide_clickchem", "Thiol_1R_clickchem"],
        "group_smarts": ["[CR1;H2,H1X4,H0X4]1O[CR1;H2,H1X4,H0X4]1", "[#6&$([CR0,R1X3,R1X4])&!$([#6](=,-[OR0,SR0])[SR0])]-[SR0;H1,-]"],
        "num_reactants": 2,
        "reaction_string": "[CR1;H2,H1X4,H0X4:1]1O[CR1;H2,H1X4,H0X4:2]1.[#6&$([CR0,R1X3,R1X4])&!$([#6](=,-[OR0,SR0])[SR0]):3]-[SR0;H1,-]>>O[C:1][C:2]S-[#6:3]",
        "RXN_NUM": 2
        },
    "3_Alkene_Oxidized_To_Epoxide":     {
        "reaction_name": "3_Alkene_Oxidized_To_Epoxide",
        "example_rxn_product": "CNC1(C)OC1(O)Br",
        "example_rxn_reactants": ["BrC(O)=C(C)NC"],
        "functional_groups": ["Alkene_clickchem"],
        "group_smarts": ["[CR0;X3,X2H1,X1H2]=[CR0;X3,X2H1,X1H2]"],
        "num_reactants": 1,
        "reaction_string": "[CR0;X3,X2H1,X1H2:1]=[CR0;X3,X2H1,X1H2:2]>>[C:1]1O[C:2]1",
        "RXN_NUM": 3
        },
    ...,
}
```

PLEASE SEE THE EXAMPLE REACTION LIBRARIES FOUND AT:

- `/autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/All_Rxns_rxn_library.json`
- `/autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/click_chem_rxns_library.json`
- `/autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/robust_rxns/Robust_Rxns_rxn_library.json`

Reaction libraries identify ligands capable of participating in a given
reaction using the information found in the sub-dictionary's items
"functional_groups" and "group_smarts".

##### Functional Group Library .json File Simple JSON Dictionary Containing Each Functional Group and Its Smarts Definition

Functional group libraries are simple dictionaries of the functional groups
used by a reaction library. Every moiety used by the reaction library must
have an entry in the functional group library.

Functional group libraries are formatted as such: (From
`click_chem_functional_groups.json`):

```json
{
    "Acid_Anhydride_Noncyclic_clickchem": "[*]C(=O)-[O;R0]-C(=O)[*]",
    "Alcohol_clickchem": "[#6&$([CR0,R1X3,R1X4])&!$([#6](=,-[OR0,SR0])[OR0])]-[OR0;H1,-]",
    "Alkene_clickchem": "[CR0;X3,X2H1,X1H2]=[CR0;X3,X2H1,X1H2]",
    "Alkyne_clickchem": "[CR0;X2,X1H1]#[CR0;X2,X1H1]",
    "Amine_2R_clickchem":  "[#7;$([#7;H3+,H2R0X1]-[#6]),$([#7&!H3;H1R1X3](:,-[#6R1]):,-[#6R1,#7R1]),$([#7&!H3;H2]-[#6]),$([#7&!H3;H0R1X2](:,-[#6R1;X3H1]):,-[#6R1X3H1]),$([#7&!H3;H0R1X2](:,-[#6R1;X3]):,-[#7R1X3]),$([#7&!H3;H1R0X3](-[#6])-[#6R0])]",
    "Azide_1R_clickchem": "[*;#6]-[$(N=[N+]=[N-]),$([N-][N+]#N)]",
    "Carbonochloridate_clickchem": "Cl[C;X3](=O)-O[*]",
    "Carboxylate_clickchem": "[*;!O]-[$([CR0;X3](=[OR0&D1])[OR0&H1]),$([CR0;X3](=[OR0&D1])[OR0-])]",
    "Epoxide_clickchem": "[CR1;H2,H1X4,H0X4]1O[CR1;H2,H1X4,H0X4]1",
    "Ester_clickchem": "[*;#6]C(=O)-O[*]",
    "Halide_clickchem": "[Cl,Br,I][$([CX4,c]),$([#6X3]=[O,S])]",
    "Isocyanate_clickchem": "[#6]N=C=O",
    "Isothiocyanate_clickchem": "[#6]N=C=S",
    "Primary_Amine_1R_clickchem": "[#7;$([H3+]),$([H2R0;!+])]-[#6]",
    "Sulfonyl_Azide_clickchem": "[*]S(=O)(=O)-[$(N=[N+]=[N-]),$([N-][N+]#N)]",
    "Thio_Acid_clickchem": "[C]-[$([CX3R0]([S;H1,X1])=[OX1]),$([CX3R0]([O;H1,X1])=[SX1])]",
    "Thiol_1R_clickchem": "[#6&$([CR0,R1X3,R1X4])&!$([#6](=,-[OR0,SR0])[SR0])]-[SR0;H1,-]"
}
```

Examples can be found here:

- `/autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/All_Rxns_functional_groups.json`
- `/autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/click_chem_functional_groups.json`
- `/autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/robust_rxns/Robust_Rxns_functional_groups.json`

The SMARTS strings provided in this file should also be present in each
sub-dictionary of the Reaction library .json file that references that
functional group, placing the name of the group in the list of functional
group names of reactants found under sub-dictionary key "functional_groups" and
placing the SMARTS string of the group in the list of functional group SMARTS
of reactants found under sub-dictionary key "group_smarts"

##### Directory of Complementary Molecule Libraries, Directory of .smi Files

Any reaction containing more than one reactant will require a complementary
molecule to supplement the reaction.

For this reason we require a directory populated with .smi files containing
small molecules that match each functional group.

The name of each .smi file should be the name of the functional group (the
keys of the functional-group-library .json file) +.smi

Example: The .smi file for the functional group
"Acid_Anhydride_Noncyclic_clickchem" should be
`/PATH_TO/complementary_mol_directory/Acid_Anhydride_Noncyclic_clickchem.smi`

THERE MUST BE ONE ENTRY PER FUNCTIONAL GROUP. BECAUSE NAMES ARE CAP SENSITIVE
IN SOME OS AND NOT IN OTHERS, PLEASE CHECK THAT YOUR NAME IS UNIQUE
INDEPENDENT OF CAPS.

#### Important Formatting Notes About the .smi File for complementary_mol_directory

1. No headers are allowed in the file.
2. .smi files can be either tab or 4-space delineated.
3. The only columns are the 1st two columns.
   - Column 1: SMILES string
   - Column 2: ligand name/identifier (1 WORD, NO SPACES)

We strongly recommend thoroughly checking that each molecule in each library
matches the intended functional group. If a ligand does not match the intended
functional group the reaction will fail and it will slow the process of mutant
creation.

#### Running Custom Reactions

Running a custom reaction library requires 4 parameters to be set:

1. `rxn_library`
2. `rxn_library_file`
3. `function_group_library`
4. `complementary_mol_directory`

Submission through .json format:

- Where JSON is located at `/PATH_TO/To/json_file_with_variable.json`
- Where reaction library JSON file is located at
  `/PATH_TO/rxn_library_file.json`
- Where function group JSON file is located at
  `/PATH_TO/function_group_library.json`
- Where directory of SMI for complementary libraries is located at
  `/PATH_TO/complementary_mol_directory/`

```json
{
    ...
    "rxn_library": "Custom",
    "rxn_library_file": "/PATH_TO/rxn_library_file.json",
    "function_group_library": "/PATH_TO/function_group_library.json",
    "complementary_mol_directory": "/PATH_TO/complementary_mol_directory/",
}
```

Submit via terminal: `python RunAutogrow.py -j
/PATH_TO/json_file_with_variable.json`

Command-line submission format:

- Where reaction library JSON file is located at:
  `/PATH_TO/rxn_library_file.json`
- Where function group JSON file is located at:
  `/PATH_TO/function_group_library.json`
- Where directory of SMI for complementary libraries is located at:
  `/PATH_TO/complementary_mol_directory/`

```bash
python RunAutogrow.py \
    ... \
    --rxn_library Custom \
    --rxn_library_file /PATH_TO/rxn_library_file.json \
    --function_group_library /PATH_TO/function_group_library.json \
    --complementary_mol_directory /PATH_TO/complementary_mol_directory/
```

### 6. Custom Complementary Molecule libraries

One can provide custom libraries of molecules (to supplement reactions) using
the `--complementary_mol_directory` option.

This can be used in conjunction with any of the predefined reactions sets
(i.e., `click_chem_rxns`, `robust_rxns`, `all_rxns`), but this requires that
all functional groups used by those reaction libraries have a corresponding
.smi file in the custom `complementary_mol_directory`

We strongly recommend thoroughly checking that each molecule in each library
matches the intended functional group. If a ligand does not match the intended
functional group the reaction will fail and it will slow the process of mutant
creation.

We provide a script to check complementary molecule libraries at
`/autogrow4/accessory_scripts/test_complementary_mol_library.py`

A tutorial is provided in the Accessory Scripts section of this document.

THERE MUST BE ONE ENTRY PER FUNCTIONAL GROUP. BECAUSE NAMES ARE CAP SENSITIVE
IN SOME OS'S AND NOT IN OTHERS, PLEASE CHECK THAT YOUR NAME IS UNIQUE
INDEPENDENT OF CAPS.

#### Important Formatting Notes About the .smi File for complementary_mol_directory

1. No headers are allowed in the file.
2. .smi files can be either tab or 4-space delineated.
3. The only columns are the 1st two columns.
   - Column 1: SMILES string
   - Column 2: ligand name/identifier (1 WORD, NO SPACES)

#### Running Custom Reactions

Submission through .json format:

- Where JSON is located at: `/PATH_TO/To/json_file_with_variable.json`
- Where directory of SMI for complementary libraries is located at:
  `/PATH_TO/complementary_mol_directory/`

```json
{
    ...
    "complementary_mol_directory": "/PATH_TO/complementary_mol_directory/",
}
```

Submit via terminal: `python RunAutogrow.py -j
/PATH_TO/json_file_with_variable.json`

Command- line submission format, where directory of SMI for complementary
libraries is located at: `/PATH_TO/complementary_mol_directory/`

```bash
python RunAutogrow.py \
    ... \
    --complementary_mol_directory /PATH_TO/complementary_mol_directory/
```

## Preparing the Receptor

AutoGrow4 takes a single .pdb file for the receptor. Although not required, we
recommend carefully preparing the receptor file prior to submitting to
AutoGrow4.

1. Remove all ligands, water, or non-protein atoms. This can be done in a PDB
   viewer such as Pymol or VMD.
   - If a ligand is already bound to target pocket, you may want to use that
     ligand to identify the pocket location prior to removing it.
2. Remove chains not being tested.
   - i.e., many protein structures contains multiple protein chains and even
     multiple proteins. We recommend removing all chains you are not
     explicitly testing. This can be done in a PDB viewer such as Pymol or
     VMD.
3. Adjust the protonation of the receptor to the appropriate pH. Crystal
   structures rarely include hydrogen atoms.
    - More accurate scoring requires proper protonation. This can be done
      using the program PDB2PQR, available via the webserver
      http://nbcr-222.ucsd.edu/pdb2pqr_2.0.0/
      - If you use the PDB2PQR to protonate the receptor, you will need to
        convert it back to pdb.
        - To convert back we recommend obabel. Installation instructions for
          obabel are provided in the Dependencies section.
        - `obabel -ipqr /PATH_TO/PQR_FILE.pqr -opdb -O
            /PATH_TO/PDB_OUTPUT_FILE.pdb`
4. Determine and define the binding pocket:
    - Docking software such as Vina and QuickVina require 6 float parameters
      to define a binding pocket:
      - Coordinates: The center of the pocket location in x,y,z axis:
        `center_x`, `center_y`, `center_z`
      - Dimensions: The distance from the center of the pocket which will be
        considered part of the pocket in x,y,z axis: `size_x`, `size_y`,
        `size_z`
    - AutoGrow4 requires all 6 parameters to run the docking portion of the
      code.
    - To determine these we recommend using the python API library scoria:
      Citation Scoria: Ropp, P., Friedman, A., & Durrant, J. D. (2017).
      Scoria: A Python module for manipulating 3D molecular data. Journal of
      cheminformatics, 9(1), 52. doi:10.1186/s13321-017-0237-8
      - Installation of Scoria:
        - Scoria can be installed either by pip installation or manual
          download.
        - We recommend pip installation: `pip install scoria`
        - Download scoria from https://durrantlab.pitt.edu/scoria/
      - Once Scoria is installed:
        1. Manually inspect the pocket of your protein in a protein
           visualizer such as Pymol, Chimera, or VMD.
           - Pick out 3 to 6 residues which will be used to define the
             protein pocket.
           - For the AutoGrow4 publication we used Chain A of the PARP-1
             catalytic domain X-ray structure 4r6e. The selected residues used
             to define the pocket were: 763, 872, 888, 907, 988
        2. Determine the geometric center of the pocket with Scoria's
           `get_geometric_center` function in python.
           - In a python terminal or in a jupyter notebook environment:

```python
# Import the scoria API
>> import scoria
>>

# define your protein pdb file
# The protein pdb file used for the publication can be found at: /autogrow4/autogrow/tutorial/PARP/4r6eA_PARP1_prepared.pdb
>> pdb_file = "/PATH_TO/OF/PDB_FILE.pdb"

# create a scoria mol object from the protein pdb file
>> mol = scoria.Molecule(pdb_file)

# select which residues are going to be used to define pocket with resseq (the residue number)
>> sel = mol.select_atoms({"resseq":[763, 872, 888, 907, 988]})

# get geometric center of the protein pocket
>> geometric_center = mol.get_geometric_center(sel)
>> print(geometric_center)
array([-70.75619481,  21.815     ,  28.32835065])
```

From this you can set: `"center_x" = -70.756,"center_y" =21.815 ,"center_z"=
28.328`

5. Determine the dimensions of the pocket with Scoria's `bounding_box`
   function in python.

```python
# Import the scoria API
>> import scoria
>>
# define the protein molecule from the PDB file
>> mol = scoria.Molecule("/PATH_TO/OF/PDB_FILE.pdb")

# select which residues are going to be used to define pocket with resseq (the residue number)
>> sel = mol.select_atoms({"resseq":[763, 872, 888, 907, 988]})

# get the dimensions of the box that encompasses the protein pocket
>> bounding_box = mol.get_bounding_box(sel)
>> mol.get_bounding_box(sel)
array([[-83.764,  15.015,  15.305],
    [-60.814,  29.578,  36.727]])
```

From this we need to take the difference from the 1st and 2nd coordinate for
x,y,z:

1. 1st box coordinate: `x_1st = -83.764,  y_1st = 15.015, z_1st = 15.305`
2. 2nd box coordinate: `x_2nd = -60.814,  y_2nd = 29.578, z_2nd = 36.727`
3. Absolute value of diff from 1st and 2nd: `"size_x" = 22.950,"size_y" =
   14.563,"size_z"= 21.422`

We suggest rounding these up to ensure the entire pocket is included:
`"size_x" = 25.00,"size_y" = 16.00,"size_z"= 25.00

## Other Factors for Consideration Prior to Running AutoGrow4

### Processors and Multiprocessing Style

AutoGrow4 is recommended to be run on a larger computer or a cluster but it
can be run on a local computer such as a laptop or PC.

#### If Running on a Laptop or PC

We recommend lowering some AutoGrow4 parameters to reduce the computational
overhead for smaller machines.

- Lower the population size and number of generations. This will mean a less
  intense search of chemistry space but will make run times more reasonable.
- Lower the `max_variation` to 1. This means for every ligand created by
  AutoGrow4, we will only create 1 conformer and thus only dock once per
  ligand. This of course means a trade-off of getting more useful information
  for each ligand for computational efficiency.

We also recommend considering how long you can allow the computer to run. If
you need to continually use the computer while running AutoGrow4 then you want
to fix the `number_of_processors` to leave several available to perform other
activities.

If you can leave the computer to run undisturbed for an extended period we
recommend setting `number_of_processors = -1`, which will use all available
processors.

#### If Running On A Larger Super Computer

We recommend fixing the `number_of_processors` to however many processors you
will be dedicating to AutoGrow4. If `number_of_processors = -1` than all
available processors will be use to run AutoGrow4.

#### If Running On A Cluster

We recommend setting the `number_of_processors = -1` and defining the number
of processors in an SBATCH-type submission script.

## Multiprocessing/MPI/Parallelization/Parallelizer

AutoGrow4 uses the `Parallelizer.py` script from Gypsum-DL
(`/autogrow4/autogrow/operators/convert_files/gypsum_dl/gypsum_dl/Parallelizer.py`).

This script creates a Parallelizer class object which can divide jobs in three
manners:

1. Serial: run all jobs 1 at a time
2. Multiprocessing: dynamically allocated distribution of jobs across multiple
   cpus on the same device
3. MPI: static allocation of jobs across many cpus across multiple machines.

### Important Notes when Running on Clusters Using SLURM

1. Multiprocessing: When running AutoGrow4 in **Multiprocessing mode** using
   SLURM, one should:
   1. 1st run the `cache_prerun` option on a single processor. `srun -n 1
      python RunAutogrow.py -c`
      - USE `srun` or `mpirun` for the `cache_prerun`. This limits the
        `prerun` to a single processor thus preventing errors caused by race
        conditions when creating pycache files.
   2. Then run AutoGrow4 as intended. `python RunAutogrow.py -j
      custom_parameters.json`
      - Do not use `srun` or `mpirun` for the production run. cpu/job
        distribution is handled internally. Using `srun` or `mpirun` can cause
        errors with the `mpi4py` universe.
2. MPI: When running AutoGrow4 in **MPI mode** using SLURM, one should:
    1. 1st run the `cache_prerun` option on a single processor. `srun -n 1
       python RunAutogrow.py -c`
       - USE `srun` or `mpirun` for the `cache_prerun`. This limits the prerun
         to a single processor thus preventing errors caused by race
         conditions when creating pycache files.
    2. Then run the simulation as intended.
        - `mpirun -n num_processors python -m mpi4py RunAutogrow.py -j
          custom_parameters.json`
        - Make sure to provide the `-m mpi4py` before `RunAutoGrow.py`. This
          tells python how to handle Exceptions.

## Accessory Scripts

AutoGrow4 provides several accessory scripts for preparing files, processing
data, and analyzing data.

These files can be found within the `/autogrow4/accessory_scripts/` folder.

### Preparation Scripts Pre-Run

#### /autogrow4/accessory_scripts/remove_duplicates_from_smi.sh

This script accepts a file path to a tab-delineated .smi file. It then filters
the file for redundancies in the 1st and 2nd columns of the file.

The output file is the input file + '_no_dup.smi'

This script uses Bash rather than Python because it is less memory intensive
when dealing with large .smi files in the millions-of-compounds range. This is
important when filtering through large databases such as ZINC15.

This script takes one input variable (`filename` str: Required). This is the
path to the tab-delineated .smi file to remove any redundancies.

Example submit:

```bash
bash /autogrow4/accessory_scripts/remove_duplicates_from_smi.sh \
    /PATH_TO/TO/SMILES.smi
```

#### /autogrow4/accessory_scripts/convert_directory_ligands_pdb_to_smi.py

This script converts a directory of pdb files (small molecules only, not
proteins) to SMILES and creates a single .smi file with all SMILES.

This script takes 3 input arguments:

1. `--source_folder` str (-s) Required. Path to folder containing .pdb files to
    convert. File must contain a single small molecules. Without proteins.
    Files must end with either .pdb or .PDB'
2. `--output_folder` str (-o) Required. Path to folder where we will output a
    .smi file of converted .pdb files.
3. `--number_of_processors` int (-p). Number of processors to use for parallel
    calculations. This script is not MPI enable but is able to multithread
    using SMP architecture. Set to -1 for all available CPUs.

Example run:

```bash
python /autogrow4/accessory_scripts/convert_directory_ligands_pdb_to_smi.py \
    --source_folder /PATH_TO/OF/PDBS/ \
    --output_folder /PATH_TO/TO/OUTPUT/ \
    --number_of_processors -1
```

#### /autogrow4/accessory_scripts/fragmenter_of_smi_mol.py

This script will fragment compounds from a .smi file. It is useful for lead
optimization. This script was used for the PARPi lead-optimization runs in the
AutoGrow4 paper.

This can fragment compounds in two manners:

1. BRICS decomposition: This fragments along synthesizable bonds
2. Fragment rotatable bonds: This breaks compounds along rotatable bonds.
   There is an option to skip carbon-carbon single bonds.

For each molecule, all permutation of fragments are calculated. For example,
fragment rotatable bonds `C-O-C1CCCC1` could produce the following fragments:

- `C-O-C1CCCC1`, not breaking any bonds
- `C` and `O-C1CCCC1`, breaking the 1st bond
- `C-O` and `C1CCCC1`, breaking the 2nd bond
- `C` and `O` and `C1CCCC1`, breaking the 1st bond and 2nd bond

A limit on maximum number of fragments per compound and a minimum number of
atoms per fragment can be set.

This script takes seven input arguments:

1. `--smi_file` str Required. Path to tab-delineated .smi file to fragment
2. `--output_smi_file` str (-o). Path to output tab-delineated .smi file of
   fragments. If not provided it will play a file in the same directory as
   smi_file titled smi_file + _Fragmented.smi
3. `--frags_per_seed_lig` int. Number of fragments to create per input SMILES.
   default is -1 which mean all possible fragments.
4. `--run_brics` bool. Whether to fragment ligands using BRICS fragmentation.
   This fragments along synthesizable bonds. Default is True.
5. `--run_frag` bool. Whether to fragment ligands over all rotatable bonds.
   Default is True.
6. `--c_c_bonds_off` bool. Whether to exclude fragmenting carbon-carbon single
   bonds. Default is True. If True it will ignore fragments on C-C bonds; if
   False it will fragment.
7. `--number_of_processors` int (-p). Number of processors to use for parallel
   calculations. This script is not MPI enable but is able to multithread
   using SMP architecture. Set to -1 for all available CPUs.

Example run:

```bash
python /autogrow4/accessory_scripts/fragmenter_of_smi_mol.py \
    -smi_file /PATH_TO/OF/SMILES.smi
```

### Preparing Custom Reaction Libraries Pre-Run

#### /autogrow4/accessory_scripts/test_complementary_mol_library.py

This script will test a complementary molecule library to ensure all compounds
react in all reactions they may be used in.

We recommend running this test if creating custom complementary libraries or
reaction libraries.

This script takes 5 input arguments:

1. `--rxn_library_file` str: Required. This PATH to a Custom json file of SMARTS
   reactions to use for Mutation.
2. `--function_group_library` str: Required. type=strThis PATH for a dictionary
   of functional groups to be used for Mutation.
3. `--complementary_mol_directory` str: Required.  This PATH to the directory
   containing all the molecules being used to react with. The directory should
   contain .smi files contain SMILES of molecules containing the functional
   group represented by that file. Each file should be named with the same
   title as the functional groups described in rxn_library_file &
   function_group_library +.smi All functional groups specified
   function_group_library must have its own .smi file. We recommend you filter
   these dictionaries prior to AutoGrow4 for the Drug-likeliness and size
   filters you will Run AutoGrow4 with.
4. `--output_folder` str: Required. This PATH to where filtered .smi file and
   log files will be placed. Will save a file in this directory for mols which
   failed sanitization, mols which failed to react in specific reactions, and
   .smi files that contain all mols that reacted properly.
5. `--number_of_processors` int (-p). Number of processors to use for parallel
   calculations. This script is not MPI enable but is able to multithread
   using SMP architecture. Set to -1 for all available CPUs.

Example submit:

```bash
python /autogrow4/accessory_scripts/test_complementary_mol_library.py \
    --rxn_library_file /autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/ClickChem_rxn_library.json \
    --function_group_library /autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/ClickChem_functional_groups.json \
    --complementary_mol_directory /autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/complementary_mol_dir \
    --output_folder /autogrow4/accessory_scripts/output/
```

### File Handling Post-Run

#### /autogrow4/accessory_scripts/convert_single_ligand_pdbqt_to_pdb.py

This script will convert a pdbqt file into a .pdb file. This is done by
removing a column of the PDB file. This script takes 2 input arguments:

1. `--pdbqt_file` str (-f) Required. Path to .pdbqt file to convert to a .pdb
   file. This must be a single ligand and must end with .pdbqt
2. `--output_folder` str (-o) Required. Path to file where we will output .pdb
   file. If not provide the output .pdb will be the same as the input
   pdbqt_file but ending with .pdb instead of .pdbqt.

Example run:

```bash
python /autogrow4/accessory_scripts/convert_single_ligand_pdbqt_to_pdb.py \
    --pdbqt_file /PATH_TO/OF/PDBQT_file.pdbqt \
    --output_folder /PATH_TO/TO/OUTPUT/FOLDER/
```

#### /autogrow4/accessory_scripts/convert_vina_docked_pdbqt_to_pdbs.py

This script will convert a docked .pdbqt.vina file into separate .pdb file.
This is done by splitting up a single .pdbqt.vina into separate .pdbqt files
for each docked pose. Then it removes a column of the .pdbqt and saves as a
.pdb file.

- If parameter `--max_num_of_poses` is not set it will convert all poses to
  .pdb
- If `--max_num_of_poses == 1` it will only convert the top docked pose to
  .pdb
- If `--max_num_of_poses == 2` it will only convert the top 2 docked poses to
  .pdb
- If `--max_num_of_poses == 10` but there only 8 poses it will convert the 8
  poses and stop
- If `--max_docking_score` is not set it will convert all poses to .pdb.

- If `--max_docking_score == -10.0` it will only convert poses with docking
  scores less than or equal to -10.0
- If both `--max_docking_score` and `--max_num_of_poses` are set they work as
  AND type operators.
- If `--max_docking_score == -11.4` and `--max_num_of_poses == 5` it will take
  the top 5 poses as long as they also have docking scores <=-11.4

Remember docking scores are better when more negative.

This script takes 6 input arguments:

1. `--vina_docked_pdbqt_file` str (-f): Required. Path to .pdbqt.vina file to
   split into 1 .pdb file per pose that matches all criteria. If this is a
   directory it will convert all of the files with the extension .pdbqt.vina
2. `--output_folder` str (-o). Path to folder where the .pdb files will be
   placed. Files will be the basename of the docked file with
   _pose_{pose_number}.pdb replacing the extension .pdbqt.vina.
3. `--max_num_of_poses` int. Each docked file will have 1 or more poses of the
   ligand. This setting controls how many are converted. default is -1 which
   means all poses possible. max_num_of_poses=1 means only the best docked
   pose will be converted. If additional criteria like max_docking_score is
   applied a pose must meet both criteria to be converted. ie) if
   max_num_of_poses= 5 and max_docking_score=-13.0 for a pose to be converted
   it must be between the 1st and 5th pose in the file and must have docked
   with a score less than or equal to -13.0.
4. `--max_docking_score` float. The most positive docking score to be converted.
   (More negative scores     are predicted to bind better). If additional
   criteria such as max_num_of_poses is applied a pose must meet both
   criteria to be converted. ie) if max_num_of_poses= 5 and
   max_docking_score=-13.0 for a pose to be converted it must be between the
   1st and 5th pose in the file and must have docked with a score less than or
   equal to -13.0.
5. `--min_docking_score` float. The most negative docking score to be converted.
   (More negative scores are predicted to bind better). If additional criteria
   such as max_num_of_poses is applied a pose must meet both criteria to be
   converted. ie) if min_docking_score= -15.0 and max_docking_score=-13.0 for
   a pose to be converted it must: -13.0. <= docking score <= -15.0
6. `--number_of_processors` int (-p). Number of processors to use for parallel
   calculations. This script is not MPI enable but is able to multithread
   using SMP architecture. Set to -1 for all available CPUs.

Example submit:

```bash
python /autogrow4/accessory_scripts/convert_vina_docked_pdbqt_to_pdbs.py \
    --vina_docked_pdbqt_file /PATH_TO/Run_1/Run_0/generation_30/PDBs/Gen_30_Cross_313228__1.pdbqt.vina \
    --output_folder /PATH_TO/outfolder/ \
    --max_num_of_poses 1 --number_of_processors -1
```

#### /autogrow4/accessory_scripts/convert_single_ligand_pdbqt_to_pdb.py

This script is used to decompress or recompress AutoGrow4 data.

If you use the `reduce_files_sizes`, option AutoGrow4 will convert concatenate
and compress all files in the PDBs directory of each generation. This is
useful when doing larger runs because data transfer will be faster and data
storage is reduced when files are merged and compressed.

The concatenation script that is run in AutoGrow4 4 can be found at
`/autogrow4/autogrow/docking/concatenate_files.py`

This script will either:

1. Return the files back to their original uncompressed and deconcatenated
   formatting
2. Concatenate and then compress the files into a single file.

The formatting of the concatenation is:

```python
"\n##############################File_name: {}\n".format(os.path.basename(file_name_1))
... Content of the 1st file...
"\n##############################$$END_FILE$$ {}".format(os.path.basename(file_name_1))
"\n##############################File_name: {}\n".format(os.path.basename(file_name_2))
... Content of the 2nd file...
"\n##############################$$END_FILE$$ {}".format(os.path.basename(file_name_2))
```

This concatenated file is then tar.gz compressed.

This script takes 2 input arguments:

1. `--compress_or_decompress` str (-s) Required. choices=["compress",
    "decompress"]. Chose whether to compress or decompress a directory
2. `--input_folder_or_file` str (-i) Required. Path to directory/file to
    compress or decompress.

Example decompression:

```bash
python /autogrow4/accessory_scripts/file_concatenation_and_compression.py \
    --compress_or_decompress decompress \
    --input_folder_or_file PATH_TO_RUN/Run_0/generation_1/PDBs/compressed_PDBS.txt.gz
```

Example compression:

```bash
python /autogrow4/accessory_scripts/file_concatenation_and_compression.py \
    --compress_or_decompress compress \
    --input_folder_or_file PATH_TO_RUN/Run_0/generation_1/PDBs/
```

### Graph Generation For Post-Run Analysis

#### /autogrow4/accessory_scripts/plot_autogrow_run.py

This script will create a line plot of the average score for each AutoGrow4
generation. This is the same type of figure as the `--generate_plot` option
that AutoGrow4 already provides, but it also allows plotting of reference
lines.

This script takes 4 input arguments:

1. `--infolder` str (-i): Required. Path to input folder containing the
    AutoGrow4 run. This should be the top folder which contains the vars.json
    file.
2. `--outfile` str (-o). Path to folder to output files. It will be created if
    does not exist. If not provide it will be placed in the
    infolder/data_line_plot.svg
3. `--outfile_format` str. choices=[ svg  png  jpg  pdf ]. The type of file for
    figure to be exported as default is .svg file.
4. `--plot_reference_lines` list. This will be a list of lists, with each
    sublist being a different dotted-line reference to plot. For each sublist
    the order of information should be: [name, value, matplotlib_color]

For example, `[['Olaparib score',-12.8,'y'],['Niraparib score',-10.7,'k']]`
will add a horizontal dotted lines at -12.8 (yellow) and -10.7 (black) with
Olaparib and Niraparib added to the legend. Spaces must be within quotes and
not be between variables. `matplotlib` colors can be found using
`mcolors.get_named_colors_mapping().keys()`

Example submit:

```bash
    python /autogrow4/accessory_scripts/plot_autogrow_run.py \
        -i /PATH_TO/Run_1/Run_0/ \
        --plot_reference_lines [['Olaparib Score',-12.8,'y'],['Niraparib',-10.7,'k'],['NAD/NADH',-10.3,'purple'],['ADP-ribose',-9.3,'maroon']]
```

#### /autogrow4/accessory_scripts/make_lineage_figures.py

This script creates figures that list all ligands that parented a given
ligand.

All compounds for the entire AutoGrow4 run will be compiled into a dictionary
that is used to trace lineages. We pickle these dictionaries so these
dictionaries do not need to be recreated every time the script is run. For
this reason the 1st time running this script will take longer than future
runs. A pre-run option will compile these data sets without generating
figures.

1. `--output_dir` str (-o): Required. Path to folder to output files. will be
    created if does not exist
2. `--input_dir` str (-i): Required. Path to input folder containing the
    AutoGrow4 run. This should be the top folder which contains the vars.json
    file.
3. `--mol_name` str: Required unless prerun. This is the name of the molecule
    whose lineage will be traced back. If not provided or None, the script
    will simply compile the necessary dictionaries/picklefiles and then
    terminate. These pickle files are stored in the input folder containing
    the vars.json file from the AutoGrow4 run. Example mol_name:
    `Gen_5_Cross_203131 or Gen_4_Mutant_7_802531`. Can also be provided as
    full-name ie: (`Gen_2_Mutant_7_97143`)`Gen_4_Mutant_7_802531`
4. `--complementary_mol_directory` str. If using a custom complementary molecule
    library for mutations this path is required. If not the script will try to
    autodetect the location of the predefined complementary_mol_directory.
    Many molecules generated by mutation will required the complementary
    molecule that helped spawn them.
5. `--source_compound_file` str: Required. This is the source .smi file used to
    seed generation zero of the AutoGrow4 run. This is an essential file.
6. `--pre_run` bool. If True this will compile the necessary
    dictionaries/picklefiles and then terminate. These pickle files are stored
    in the input folder containing the vars.json file from the AutoGrow4 run.

Example submit:

```bash
python /autogrow4/accessory_scripts/make_lineage_figures.py \
    -o /PATH_TO/output_folder/ \
    -i /PATH_TO/INPUT_RUN/Run_3/ \
    --source_compound_file /autogrow4/source_compounds/PARPI_BRICS_frags.smi \
    --mol_name Gen_17_Cross_727024
```
