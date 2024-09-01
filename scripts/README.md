## System requirements

### Software requirements

These scripts are supported on Linux and have been tested on the following distributions:
  - Debian GNU/Linux 12
  - Ubuntu 22.04.3 LTS

Dependencies:
  - The `run_ligprep.sh` script requires [Schrödinger](https://www.schrodinger.com/)'s LigPrep to be installed.
  - All other scripts rely solely on Python dependencies. Ensure that the necessary Python packages are installed as specified in the installation guide.

### Hardware Requirements

- Minimum: A standard computer with at least 2GB of RAM to handle in-memory operations.
- Recommended: For handling large virtual screens or intensive computations, it is recommended to run these scripts on a High-Performance Computing (HPC) cluster.

## Installation guide




## Demo

### script 1: run ligprep
- input: ligand name and smiles
- output: 3D representation of ligand
- out format: mol2

### script 2: dock
- input: receptor structure, ligand structure
- output: docked ligand pose, flexible residues pose, docking score
- out format: pdbqt, log

#### script 3: repair smina output
- input: receptor structure, smiles (?), smina output 
- output: un-protonated (?) ligand in docking pose inside receptor w/ adjusted flexible residues 
- out format: pdb

#### script 4: docking md
- input: smiles, output of script 3 (pdb-file) 
- output: a) protonated ligand in docking pose inside receptor w/ adjusted flexible residues; b) trajectory of ligand-receptor system in solvent 
- out format: pdb, xtc

#### script 5.1: pre-compute observables for scoring
- input: output of script 4 (b) 
- output: scoring-relevant observables sorted by ligand 
- out format: pickle files containing dictionaries

#### script 5.2: compute scores
- input: observables (output of 5.1) 
- output: sorted list of scored ligands 
- out format: csv (pandas) 

#### utilities:
- periodic.py: provides a function to map ligand back into the binding pocket when broken across periodic boundaries.


## Instructions for use

1. Generating ligand structures
`./run_ligprep.sh data/ligands data/smiles.smi`

2. Preparing receptor and ligand structures
`conda activate vs_prep`
`cp ../receptor_ensemble/* data/receptors`
`./prepare_receptors.sh data/receptors data/receptors.list`
`./prepare_ligands.sh data/ligands data/ligands.list`

3. Running docking simulations


4. Running molecular dynamics simulations


5. Calculating custom scores



## License
Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany

These scripts are free software; you can redistribute them and/or modify
them under the terms of the MIT license. See LICENSE for details.
