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




## Instructions

### script 1: `run_ligprep.sh`
- input: ligand name and smiles
- output: 3D representation of ligand
- out format: mol2

### script 2: `prepare_receptors.sh`
- input: receptor pdb file
- output: prepared receptor
- out format: pdbqt

### script 3: `prepare_ligands.sh`
- input: ligand mol2 file
- output: prepared ligand
- out format: pdbqt

### script 4: `dock.sh`
- input: receptor structure, ligand structure
- output: docked ligand pose, flexible residues pose, docking score
- out format: pdbqt, log

### script 5: `smina_atomassign.py`
- input: receptor structure, smiles, output of script 4 
- output: unprotonated ligand docked to receptor w/ adjusted flexible residues 
- out format: pdb

### script 6: `docking_md.py`
- input: smiles, output of script 5 
- output: a) protonated ligand docked to receptor w/ adjusted flexible residues; b) trajectory of ligand-receptor system in solvent 
- out format: pdb, xtc

### script 7: `observables.py`
- input: output of script 6 (b)
- output: scoring-relevant observables sorted by ligand 
- out format: pickle files containing dictionaries

### script 8: `compute_mdscores.py`
- input: output of script 7
- output: sorted list of scored ligands 
- out format: csv

### utilities:
- periodic.py: provides a function to map ligand back into the binding pocket when broken across periodic boundaries.


## Demo

1. Generating ligand structures\
`./run_ligprep.sh data/ligands data/smiles.smi`

2. Preparing receptor and ligand structures\
`conda activate vs_prep`\
`mkdir data/receptors`\
`cp ../receptor_ensemble/* data/receptors`\
`./prepare_receptors.sh data/receptors data/receptors.list`\
`./prepare_ligands.sh data/ligands data/ligands.list`

3. Running docking simulations\
`conda activate vs_dock`\
`mkdir data/dockings`\
`./dock.sh data/ligands data/receptors data/dockings drugbank_ZINC000003874467 apo_hmm0_rndm_frm130`

4. Repairing docking output\


5. Running molecular dynamics simulations\


6. Calculating custom scores



## License
Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany

These scripts are free software; you can redistribute them and/or modify
them under the terms of the MIT license. See LICENSE for details.
