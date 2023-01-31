#!/usr/bin/env python

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.


import numpy as np
import matplotlib.pyplot as plt
import mdtraj
import scipy.stats
from glob import glob
import pickle
from pathos.multiprocessing import Pool
from contextlib import closing
import itertools
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

import sys

JOBID = sys.argv[1]
print('jobid', JOBID)

n_proc = 24
traj_length_str = '100ns'
total_n = 2215  # number of computations to distribute

# devide by number of compute nodes (20); TODO: replace manual parallelization
_lig_rec_pairs = np.loadtxt(f'../enamine_MD_list_long.txt',
                            dtype=str, comments=None)[
                 int(total_n / 20) * int(JOBID):int(total_n / 20) * (int(JOBID) + 1)]

residue_groups = {'S1': np.concatenate([np.arange(180, 187), np.arange(204, 210)]) - 1,
                  'hydroph': np.arange(24, 27) - 1,
                  'flex': np.arange(41, 46) - 1}

sasa_residues = np.concatenate(list(residue_groups.values()))

smiles_path = 'smiles_final_enamine_fixed.csv'
smiles = pd.read_csv(smiles_path, sep=' ', header=None, index_col=0).to_dict()[1]

covalent_warhead_smarts = {
    'ester': '[*]-C(=O)O-[*]',  # check C (atom in pos 1)
    'PMSF': 'O=S(=O)(F)[*]',  # check S (atom in pos 1)
    'CMK': '[*]-C(=O)CCl',  # check C (atom in pos 1)
    'aldehyde': '[*]-C-[CH]=O',  # check second C (atom in pos 2)
    'TFK': '[*]-C(=O)C(F)(F)F'  # check first C (atom in pos 1)
}

_dir = '../../'


def compute_observables(d_run, lig, rec):
    results_list_dict = {}
    _files = glob(f'{_dir}/dockings_{d_run}/{lig}/{rec}_*_{traj_length_str}.pdb')
    trajfiles = [f.replace(f'_{traj_length_str}.pdb', '-protein.xtc').replace(f'{rec}_', f'{rec}-') for f in _files]
    for trajfile in trajfiles:
        name = trajfile.split('/')[-1]
        results = {}
        try:
            traj = mdtraj.load(trajfile,
                               top=f'{_dir}/dockings_{d_run}/{lig}/{rec}_docked2_equil.pdb',
                               stride=10)

            # unfortunate combination of 2 problems: a) MD SASA crashes if two atoms overlap;
            # b) forcefield (?) fails to resolve hydrogen-hydrogen repulsion at least once
            # remove these trajectories for now
            for n in range(traj.n_frames):
                d = np.linalg.norm(traj.xyz[n] - traj.xyz[n, :, None], axis=2)
                d[np.diag_indices(d.shape[0])] = 999

                if d.min() < 0.025:
                    print('atom overlap for ', d_run, lig, rec)
                    return

            results['n_heavy'] = traj.top.select('not protein and mass > 2').shape[0]

            pairs = list(itertools.product(range(traj.n_residues), [traj.n_residues - 1]))
            all_distances = mdtraj.compute_contacts(traj, contacts=pairs)[0]

            for k, positions in residue_groups.items():
                distances = all_distances[:, positions]
                num_contacts_ligand = (distances <= 0.35).sum(axis=1).astype(np.float32)
                results['contacts_' + k] = num_contacts_ligand

                if k == 'S1':
                    results['dist_asp180'] = distances[:, 0]

            sasa_complex = mdtraj.shrake_rupley(traj, mode='residue')
            sasa_protein = mdtraj.shrake_rupley(traj.atom_slice(traj.top.select('protein')),
                                                mode='residue')

            results['dsasa'] = (sasa_complex[:, :-1] - sasa_protein)[:, sasa_residues]

            try:
                mol = Chem.MolFromPDBFile(f'{_dir}/dockings_{d_run}/{lig}/{rec}_docked2_equil.pdb')
                mol = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles[lig]), mol)
                coord_OG = np.asarray(mol.GetConformer().GetAtomPosition(1439))  # OG atom coord of SER186
                min_dist_RC = np.inf  # if there are multiple RCs takes the one closest to OG of SER186
                min_dist_arr_RC = np.zeros((traj.n_frames)) + np.inf

                for warhead in covalent_warhead_smarts:
                    pattern = Chem.MolFromSmarts(covalent_warhead_smarts[warhead])
                    if mol.HasSubstructMatch(pattern):
                        matches = mol.GetSubstructMatches(pattern)  # returns a tuple of tuples
                        for match in matches:
                            idx = match[2 if warhead == 'aldehyde' else 1]  # position in the pattern
                            coord_RC = np.asarray(mol.GetConformer().GetAtomPosition(idx))  # reactive center coord
                            dist = np.linalg.norm(coord_OG - coord_RC)  # should be <=3.5

                            # assume that the atomic order did not get messed up again
                            # (only validation: element order is o.k., no atom names found)
                            # keeping it separate for that resason
                            idx_withH = traj.top.select('mass > 2')[idx]
                            serOG_withH = traj.top.select('residue 186 and name OG')[0]
                            min_dist_arr = mdtraj.compute_distances(traj, [[idx_withH, serOG_withH]]).squeeze()

                            if dist < min_dist_RC:
                                min_dist_RC = dist
                                min_dist_arr_RC = min_dist_arr
                results['dist_reactive'] = min_dist_RC
                results['dist_reactive_arr'] = min_dist_arr_RC

            except ValueError as e:  # rdkit fails to load the smiles
                print('failed loc 1')
                print(e)
                results['dist_reactive'] = np.NaN

            results_list_dict[name] = results

        except Exception as e:  # was: OSError; there were 2 failures in 20,000 that crhased the whole thing
            print('failed loc 2')
            print(e)
            results_list_dict[name] = None
    return results_list_dict


distance_asp180 = {}
contacts = {k: {} for k in residue_groups.keys()}
dsasas = {}
n_heavy = {}

args = [(l[0], l[1], l[2]) for l in _lig_rec_pairs]

pbar = tqdm(total=len(args))
pool = Pool(processes=n_proc)

all_results = {lig: {} for lig in np.unique(_lig_rec_pairs[:, 1])}

with closing(pool):
    for a in args:
        x = pool.apply_async(compute_observables, a,
                             callback=lambda _: pbar.update(1))
        all_results[a[1]][a[2]] = x.get()

pbar.close()

pickle.dump(all_results, open('all_results_' + str(JOBID) + '.pickle', 'wb'))
