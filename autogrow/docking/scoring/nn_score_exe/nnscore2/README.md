# NNScore 1 and NNScore 2 #

## Home ##

As high-throughput biochemical screens are both expensive and labor intensive,
researchers in academia and industry are turning increasingly to
virtual-screening methodologies. Virtual screening relies on scoring functions
to quickly assess ligand potency. Unfortunately, these scoring functions
generally have many false positives and negatives; indeed, a properly trained
human being can often assess ligand potency by visual inspection with greater
accuracy.

Given the success of the human mind at receptor-ligand complex
characterization, we here present two scoring functions based on neural
networks, computational models that simulate the microscopic organization of
the brain. Computer-aided drug design depends on fast and accurate scoring
functions to aid in the identification of small-molecule ligands. The scoring
functions presented here, either on their own or used in conjunction with
other more traditional functions, may prove useful in drug-discovery projects.
Additional information about NNScore 1 can be found [in the original
paper](http://pubs.acs.org/doi/full/10.1021/ci100244v). NNScore 2 is described
in [a separate publication](https://pubs.acs.org/doi/abs/10.1021/ci2003889).

Note that NNScore 2 is not necessarily superior to NNScore 1. The best scoring
function to use is highly system dependent. Including positive controls (known
inhibitors) in virtual screens is a useful way to identify which scoring
function is best suited to your needs.

If you use NNScore in your research, **please cite the appropriate
reference**:

NNScore: A Neural-Network-Based Scoring Function for the Characterization of
Protein-Ligand Complexes. Jacob D. Durrant, J. Andrew McCammon. Journal of
Chemical Information and Modeling, 2010, 50 (10), pp 1865-1871.

NNScore 2.0: A Neural-Network Receptor–Ligand Scoring Function. Jacob D.
Durrant, Andrew McCammon. Journal of Chemical Information and Modeling, 2011,
51 (11), pp 2897-2903.

## Usage for Version 1.X ##

NNScore 1 has been implemented as a [python](http://www.python.org/) script.
The program accepts the following parameters:

```text
-receptor <pdbqt filename>
-ligand <pdbqt filename>
-network <network filename>
-networks_dir <directory>
```

Note: It is best to use multiple neural networks to judge ligand binding by
consensus. Command-line parameters can be used to add neural-network files to
the list of those that will be used. To add a single neural network to the
list, use the -network parameter to specify a single network file. To add
multiple networks to the list, create a directory containing only network
files and specify the path to that directory using the -networks_dir
parameter.

Note: Only PDBQT files of the receptor and ligand are accepted. Scripts to
convert from PDB to PDBQT are included in the [AutoDockTools
package](http://autodock.scripps.edu/resources/adt). **Be sure to use
AutoDockTools to convert from PDB to PDBQT, not Open Babel.** These two
programs do not assign the same partial atomic charges, and NNScore was
trained using AutoDockTools-assigned charges.

Examples:

```bash
python NeuroScore.py -receptor neuraminidase.pdbqt
    -ligand oseltamivir.pdbqt
    -network ./networks/top_3_networks/12.net
```

```bash
python NeuroScore.py -receptor integrase.pdbqt
    -ligand raltegravir.pdbqt
    -networks_dir ./networks/top_3_networks/
```

```bash
python NeuroScore.py -receptor protease.pdbqt
    -ligand tipranavir.pdbqt
    -networks_dir ./networks/top_24_networks/
    -network ./networks/top_3_networks/16.net
```

## Usage for Version 2.X ##

NNScore 2s has also been implemented as a [python](http://www.python.org/)
script. As demonstrated in [our
paper](https://pubs.acs.org/doi/abs/10.1021/ci2003889), NNScore 2 is not
necessarily superior to NNScore 1. The best scoring function to use is highly
system dependent. Including positive controls (known inhibitors) in virtual
screens is a useful way to identify which scoring function is best suited to
your needs.

### Requirements ###

* Python3: A copy of the Python interpreter can be downloaded from
  [http://www.python.org/getit/](http://www.python.org/getit/).
* AutoDock Vina 1.1.2: NNScore 2 uses AutoDock Vina 1.1.2 to obtain some
  information about the receptor-ligand complex. Note that previous versions
  of AutoDock Vina are not suitable. AutoDock Vina 1.1.2 can be downloaded
  from
  [http://vina.scripps.edu/download.html](http://vina.scripps.edu/download.html).
* MGLTools: As receptor and ligand inputs, NNScore 2s accepts models in the
  PDBQT format. Files in the more common PDB format can be converted to the
  PDBQT format using scripts included in MGLTools (`prepare_receptor4.py` and
  `prepare_ligand4.py`). MGLTools can be obtained from
  [http://mgltools.scripps.edu/downloads](http://mgltools.scripps.edu/downloads).
  **Be sure to use MGLTools to convert from PDB to PDBQT, not Open Babel.**
  These two programs do not assign the same partial atomic charges, and
  NNScore was trained using MGLTools-assigned charges.

### Command-Line Parameters ###

`-receptor`: File name of the receptor PDBQT file.

`-ligand`: File name of the ligand PDBQT file. AutoDock Vina output files,
typically containing multiple poses, are also permitted.

`-vina_executable`: The location of the AutoDock Vina 1.1.2 executable. If you
don't wish to specify the location of this file every time you use NNScore 2,
simply edit the `vina_executable` variable defined near the beginning of the
NNScore2.py script.

### Program Output ###

NNScore 2 evaluates each of the ligand poses contained in the file specified
by the -ligand tag using 20 distinct neural-network scoring functions. The
program then seeks to identify which of the poses has the highest predicted
affinity using several metrics:

1) Each of the 20 networks are considered separately. The poses are ranked in
   20 different ways by the scores assigned by each network.
2) The poses are ranked by the best score given by any of the 20 networks.
3) The poses are ranked by the average of the scores given by the 20 networks.
   **This is the recommended metric.**

### Example of Usage ###

```bash
python NNScore2.py -receptor myreceptor.pdbqt -ligand myligand.pdbqt -vina_executable /PATH/TO/VINA/1.1.2/vina
```

## Download ##

All versions of NNScore are released under the [GNU General Public
License](https://git.durrantlab.pitt.edu/jdurrant/nnscore2/blob/master/gpl-3.0.txt).
Your use of NNScore implies acceptance of the terms stipulated in that
license.

* Visit
  [https://git.durrantlab.pitt.edu/jdurrant/nnscore1](https://git.durrantlab.pitt.edu/jdurrant/nnscore1)
  to download the latest version of NNScore 1.
* Visit
  [https://git.durrantlab.pitt.edu/jdurrant/nnscore2](https://git.durrantlab.pitt.edu/jdurrant/nnscore2)
  to download the latest version of NNScore 2.

## Contact ##

If you have any questions, comments, or suggestions, please don't hesitate to
contact me, [Jacob Durrant](http://durrantlab.com), at durrantj [at] pitt
[dot] edu. I'd be happy to help. :)
