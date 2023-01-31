#### script 1: run ligprep
- input: ligand name and smiles
- output: 3D representation of ligand
- out format: mol2

#### script 2: dock
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


#### License
Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany

These scripts are free software; you can redistribute them and/or modify
them under the terms of the MIT license. See LICENSE for details.
