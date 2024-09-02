#!/usr/bin/env python
# coding: utf-8

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.


import argparse
import os
from copy import deepcopy

import mdtraj
import numpy as np
import simtk.openmm.app as app
import simtk.unit as u
from scipy.optimize import linear_sum_assignment as lsa

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir', type=str,
                    help='docking directory')
parser.add_argument('-l', '--ligand', type=str,
                    help='ligand name')
parser.add_argument('-r', '--receptors_dir', type=str,
                    help='path to receptors files')
parser.add_argument('-x', '--receptors_list', type=str,
                    help='path to receptors list')
args = parser.parse_args()

# set up basic receptor to get started with atomic indexing
receptors_list = np.loadtxt(args.receptors_list, dtype=str)
receptor_base = mdtraj.load_pdb(f'{args.receptors_dir}/{receptors_list[0]}.pdb')
receptor_base_noh = receptor_base.atom_slice(receptor_base.top.select('mass > 2'))


# katarina's functions, minor TH edit
def generate_splitted_flex_files(flex_path, num_flex=5):
    handles = []
    fnames = []
    for file_num in range(num_flex):
        fname = flex_path.split('/')[-2] + '_' + flex_path.split('/')[-1][:-10] + str(file_num + 1) + '.pdb'
        handles.append(open(fname, 'w'))
        fnames.append(fname)
    with open(flex_path) as f:
        file_num = 0
        for line in f:
            handles[file_num % num_flex].write(line)
            if line[0:6] == 'ENDMDL':
                file_num += 1
            if file_num == num_flex:  # first frame only
                break

    for file_num in range(num_flex):
        handles[file_num].close()

    return fnames


def load_flex_structure(flex_path, num_flex=5):
    flex_residues = []
    for file_num in range(num_flex):
        flex_residues.append(mdtraj.load_pdb(flex_path.split('/')[-2] + '_' +
                                             flex_path.split('/')[-1][:-10] +
                                             str(file_num + 1) + '.pdb'))

    flex_structure = flex_residues[0]
    for file_num in range(1, num_flex):
        flex_structure = flex_structure.stack(flex_residues[file_num])

    return flex_structure


def remove_splitted_flex_files(flex_path, num_flex=5):
    for file_num in range(num_flex):
        os.remove(flex_path.split('/')[-2] + '_' +
                  flex_path.split('/')[-1][:-10] +
                  str(file_num + 1) + '.pdb')


# Simon's function, TH edited
def match_atom_assignment(mol1, mol2, enforce_atype_match=True,
                          penalty_weight=1000.):
    """
    Arguments
    --------
    mol1 : mdtraj.core.trajectory.Trajectory - reference molecule
    mol2 : mdtraj.core.trajectory.Trajectory -
    enforce_atype_match : bool
    penalty_weight : float

    Returns
    -------
    (reference_permuation, minimal_error_permutation)
        
    """
    # Compute distance matrices
    Dmat1 = np.sum(((mol1.xyz[0, ..., None] - mol1.xyz[0, ..., None].T)) ** 2., axis=1) ** (0.5)
    Dmat2 = np.sum(((mol2.xyz[0, ..., None] - mol2.xyz[0, ..., None].T)) ** 2., axis=1) ** (0.5)

    # Get Atom types
    atypes1 = np.array([a.element.symbol for a in mol1.top.atoms])
    atypes2 = np.array([a.element.symbol for a in mol2.top.atoms])

    # compute atom type enforcement penalty
    if enforce_atype_match:
        atom_type_penalty = penalty_weight * (1 - (atypes1.reshape(-1, 1) == atypes2.reshape(1, -1)))
    else:
        atom_type_penalty = 0

    # compute assignment cost matrix
    C = np.sum(np.abs(Dmat1[..., None] - Dmat2[None, ...]), axis=1)

    ref_perm, target_perm = lsa(C + atom_type_penalty)
    return target_perm


# smina selection of "sidechains" for the selected residue
# it generally does not have hydrogens except for 
# trp HE1
# lys hz1 hz2 hz3
# gln he21 he22
# it further includes CA and C, even though these are not moved (according to some
# post by a person involved in programming this software)

# selection string A or (B or C or D) causes problems -> converted to A or B or C or D
flexres_selection = '(sidechain or name C CA) and (mass > 2 or (resname TRP and name HE1) or (resname LYS and name HZ1 HZ2 HZ3) or (resname GLN and name HE21 HE22)) and residue 44 45 180 183 206'
flexres_atomidx = receptor_base.top.select(flexres_selection)
flexres_elems = {
    'C': receptor_base.top.select(flexres_selection + ' and element C'),
    'O': receptor_base.top.select(flexres_selection + ' and element O'),
    'N': receptor_base.top.select(flexres_selection + ' and element N'),
    'H': receptor_base.top.select(flexres_selection + ' and element H')}

elements = ['C', 'O', 'N', 'H']
# additionally, smina changes atomic order.
# manual assignment works only sometimes

all_but_flexres_missing_hydros_selstr = '(mass < 2 and not (resname TRP and name HE1) and not (resname LYS and name HZ1 HZ2 HZ3) and not (resname GLN and name HE21 HE22)) and residue 44 45 180 183 206 and not name H HA'

# this is updated in the for loop again (apo vs. drugbound structures are different)
all_but_flexres_missing_hydros = np.setdiff1d(np.arange(receptor_base.n_atoms),
                                              receptor_base.top.select(all_but_flexres_missing_hydros_selstr))

n_missing_hydros = len(receptor_base.top.select(all_but_flexres_missing_hydros_selstr))

# these hydrogens will be deleted and re-added by openMM
for n in range(receptor_base.n_atoms):
    if n not in all_but_flexres_missing_hydros:
        print(receptor_base.top.atom(n))


malfunctions = []

for receptor_name in receptors_list:
    receptor_base = mdtraj.load_pdb(f'{args.receptors_dir}/{receptor_name}.pdb')

    all_but_flexres_missing_hydros = np.setdiff1d(np.arange(receptor_base.n_atoms),
                                                    receptor_base.top.select(all_but_flexres_missing_hydros_selstr))

    # for openMM hydrogen adding: define residue variants (only cystein bridges necessary)
    cbridge_resids = np.array([155, 171, 182, 210, 26, 42]) - 1
    res_variants = ['CYX' if (res.index in cbridge_resids) else None for res in receptor_base.top.residues]
    
    if os.path.exists(f'{args.dir}/{args.ligand}/{receptor_name}.pdbqt.pdb') and not os.path.exists(
            f'{args.dir}/{args.ligand}/{receptor_name}_docked2.pdb'):
        try:
            pose = mdtraj.load_pdb(f'{args.dir}/{args.ligand}/{receptor_name}.pdbqt.pdb',
                                    frame=0)

            flex_path = f'{args.dir}/{args.ligand}/{receptor_name}_flex.pdbqt.pdb'

            fnames = generate_splitted_flex_files(flex_path)
            flex_structure = load_flex_structure(flex_path)

            smina_elems = []
            for fname in fnames:
                with open(fname) as f:
                    for line in f.readlines():
                        if line.startswith("ATOM"):
                            smina_elems.append(line[:-1].split()[-1])
            smina_elems = np.array(smina_elems)
            assert set(smina_elems) == set(elements)
            assert len(smina_elems) == flex_structure.n_atoms

            remove_splitted_flex_files(flex_path)

            assert flexres_atomidx.shape[0] == flex_structure.n_atoms

            receptor = deepcopy(receptor_base)
            n_atoms_orig = receptor.n_atoms

            # overwrite flexible residue coordinates
            for e in elements:
                # exclude atoms that have same coordinates, such as CA and C
                # otherwise VMD graphics are messed up (e.g. if CA is at CD position)

                # identify where the coordinates are the same
                same_rec, same_flex = np.where((receptor.xyz[0, flexres_elems[e]][..., np.newaxis] ==
                                                flex_structure.xyz[0, smina_elems == e][..., np.newaxis].T).all(
                    1))
                # only overwrite the ones that are not the same within the same element
                receptor.xyz[0, flexres_elems[e][np.setdiff1d(np.arange(len(flexres_elems[e])), same_rec)]] = \
                flex_structure.xyz[0, smina_elems == e][
                    np.setdiff1d(np.arange(len(flexres_elems[e])), same_flex)]

            # hydrogen free version
            receptor_noh = receptor.atom_slice(receptor.top.select('mass > 2'))

            # fix atom assignments
            # all flexible residue atoms (noh)
            flexres_atomidx_noh = receptor_noh.top.select('residue 44 45 180 183 206')

            # find permutation by residue wise optimization
            # otherwise only minor flips can be recognized
            perm = []
            for resid in [44, 45, 180, 183, 206]:
                _flexres_atomidx_noh = receptor_noh.top.select(f'residue {resid}')
                mol1 = receptor_noh.atom_slice(_flexres_atomidx_noh)
                mol2 = receptor_base_noh.atom_slice(_flexres_atomidx_noh)
                p = match_atom_assignment(mol2, mol1, enforce_atype_match=True)

                perm += list(_flexres_atomidx_noh[p])
            perm = np.array(perm)

            # overwrite positions using found permutation
            receptor_noh.xyz[0, flexres_atomidx_noh] = receptor_noh.xyz[0, perm]

            # write into all-atom version:
            receptor.xyz[0, receptor.top.select('mass > 2')] = receptor_noh.xyz[0, :]

            # slice out missing sidechain hydrogens
            receptor_noflexH = receptor.atom_slice(all_but_flexres_missing_hydros)

            # fix missing sidechain hydrogens with openMM
            modeller = app.Modeller(receptor_noflexH.top.to_openmm(),
                                    receptor_noflexH.xyz[0] * u.nanometer)

            # make sure cystein bridges stay intact
            modeller.addHydrogens(variants=res_variants)
            assert modeller.getTopology().getNumAtoms() == n_atoms_orig

            # write to disk
            # port back to openMM. overwrite existing positions because
            # openMM changes backbone hydrogen positions
            new_receptor = mdtraj.Trajectory(modeller.getPositions() / u.nanometer,
                                                mdtraj.Topology.from_openmm(modeller.getTopology()))
            new_receptor.xyz[0, all_but_flexres_missing_hydros] = receptor.xyz[
                0, all_but_flexres_missing_hydros]
            new_receptor.save_pdb(f'{args.dir}/{args.ligand}/{receptor_name}_receptor.pdb')

            # stack with top scoring docked pose
            receptor_ligand = new_receptor.stack(pose)
            receptor_ligand.save_pdb(f'{args.dir}/{args.ligand}/{receptor_name}_docked2.pdb')

            # save mdtraj converted pose (only pose)
            pose.save_pdb(f'{args.dir}/{args.ligand}/{receptor_name}_bestpose.pdb')

        except IndexError:
            malfunctions.append(f'{args.dir}/{args.ligand}/{receptor_name}.pdbqt.pdb')

np.savetxt(f'{args.dir}/{args.ligand}/failed.log', malfunctions, fmt='%s')
