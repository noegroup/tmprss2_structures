## System requirements

### Software requirements

These scripts are supported on Linux and have been tested on the following distributions:
  - Debian GNU/Linux 12
  - Ubuntu 22.04.3 LTS

Dependencies:
  - The `run_ligprep.sh` script requires [Schrödinger](https://www.schrodinger.com/)'s LigPrep (Schrödinger Suite 2020-2) to be installed.
  - All other scripts rely solely on Python dependencies. Ensure that the necessary Python packages are installed as specified in the installation guide.

### Hardware Requirements

- Minimum: A standard computer with at least 2GB of RAM to handle in-memory operations.
- Recommended: For handling large virtual screens or intensive computations, it is recommended to run these scripts on a High-Performance Computing (HPC) cluster.


## Installation guide

The installation process involves setting up three environments: preparation, docking and molecular dynamics (MD). The typical install time on a standard desktop computer with a stable internet connection is around 15-30 minutes, accounting for both package downloads and the time it takes for the package manager to solve the environment.

### Create preparation environment
`conda create --name vs_prep`\
`conda activate vs_prep`\
`conda install -c bioconda mgltools=1.5.6`

### Create docking environment
`conda create --name vs_dock`\
`conda activate vs_dock`\
`conda install -c bioconda smina=2017.11.9`\
`conda install -c openbabel openbabel=2.4.1`\
`conda install -c conda-forge mdtraj=1.9.4`\
`conda install -c omnia openmm=7.4.2`

### Create MD environment
`conda create --name vs_md`\
`conda activate vs_md`\
`conda install -c conda-forge python=3.7`\
`conda install -c conda-forge mdtraj=1.9.4`\
`conda install -c omnia openmm=7.4.2`\
`conda install -c conda-forge -c omnia openmmforcefields=0.8.0`\
`conda install -c omnia openforcefield=0.7.1`\
`conda install -c conda-forge pathos`\
`conda install -c conda-forge tqdm`

## Details

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

### script 5: `run_obabel.sh`
- input: docked ligand pose, flexible residues pose
- output: docked ligand pose, flexible residues pose
- out format: pdb

### script 6: `smina_atomassign.py`
- input: receptor structure, output of script 5 
- output: unprotonated ligand docked to receptor w/ adjusted flexible residues 
- out format: pdb

### script 7: `docking_md.py`
- input: smiles, output of script 6
- output: a) protonated ligand docked to receptor w/ adjusted flexible residues; b) trajectory of ligand-receptor system in solvent 
- out format: pdb, xtc

### script 8: `observables.py`
- input: output of script 7 (b)
- output: scoring-relevant observables sorted by ligand 
- out format: pickle files containing dictionaries

### script 9: `compute_mdscores.py`
- input: output of script 8
- output: sorted list of scored ligands 
- out format: csv

### utilities:
- periodic.py: provides a function to map ligand back into the binding pocket when broken across periodic boundaries.


## Demo

1. Generating ligand structures (already provided for test dataset)\
`./run_ligprep.sh data/ligands data/ligands.smi`\
Expected run time: ~1 second per ligand.

2. Preparing receptor and ligand structures\
`conda activate vs_prep`\
`mkdir data/receptors`\
`cp ../receptor_ensemble/* data/receptors`\
`./prepare_receptors.sh data/receptors data/receptors.list`\
`./prepare_ligands.sh data/ligands data/ligands.list`\
Expected run time: ~20 seconds for 20 receptors + ~1 second per ligand.

3. Running docking simulations\
`conda activate vs_dock`\
`mkdir data/dockings`\
`./dock.sh data/ligands data/receptors data/dockings drugbank_ZINC000003874467 apo_hmm0_rndm_frm130`\
Expected run time: ~20 minutes per receptor-ligand pair on a single CPU core.

4. Converting docking output\
`conda activate vs_dock`\
`./run_obabel.sh data/dockings drugbank_ZINC000003874467`\
`./smina_atomassign.py --dir data/dockings --ligand drugbank_ZINC000003874467 --receptors_dir data/receptors --receptors_list data/receptors.list`\
Expected run time: ~30 seconds per receptor-ligand pair.

5. Running molecular dynamics simulations\
`conda activate vs_md`\
`./docking_md.py --dir data/dockings --ligand drugbank_ZINC000003874467 --receptor apo_hmm0_rndm_frm130 --smiles_file data/ligands.smi --platform CPU`\
Expected run time: ~50 minutes per receptor-ligand pair on a single GPU core.

6. Calculating target-specific scores
`conda activate vs_md`\
`./observables.py --dir data/dockings --pairs_list data/pairs.list --smiles_file data/ligands.smi --n_proc 4`\
`./compute_mdscores.py --dir data/dockings --pairs_list data/pairs.list --pre_computed data/dockings/all_results_merged.pickle`\
Expected run time: ~2 minutes per receptor-ligand pair.

In total, for a single receptor-ligand pair, the demo is expected to take approximately 70-75 minutes, depending on the performance of the CPU and GPU cores used.


## License
Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany

These scripts are free software; you can redistribute them and/or modify
them under the terms of the MIT license. See LICENSE for details.
