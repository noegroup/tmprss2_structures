#!/usr/bin/env python
# coding: utf-8

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.


import mdtraj
import numpy as np
import os, copy
import bz2
import pickle

from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as u
import openmmforcefields.generators
import openforcefield.topology

import networkx as nx
import networkx.algorithms.isomorphism as iso

from mdtraj.reporters import DCDReporter
import time, datetime

from periodic import *
import argparse

# adjustables
boxsize = 7.2 * u.nanometer  # x, y, z extend (cubic box)
integrator_timestep_ps_equil = 0.002  # picoseconds, pre-equilibration without heavy hydrogens
integrator_timestep_ps = 0.004  # picoseconds
simulation_time_ns = 100  # for production run
simulation_time_ns_equil_nvt = .1  # for each NVT equilibration (w/ constraints)
simulation_time_ns_equil_npt = .9  # for NPT equil (w/ constraints)
temperature_K = 310  # kelvin

# saving interval (production run only)
save_traj_ps = 100
save_rmsd_ps = 10

_dir = f'/data/scratch/kelez/virtual_screening/'
cache_dir = f'{_dir}/openff_cache/'
smiles_file = f'{_dir}/smiles_final_enamine_fixed.csv'

# time stamp for output files
stamp = datetime.datetime.now().strftime('%Y%b%d-%H%M')

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--ligand', type=str)
parser.add_argument('-r', '--receptor', type=str)
parser.add_argument('-d', '--run', type=str)

args = parser.parse_args()
ligand = args.ligand
receptor_name = args.receptor
d_run = args.run

print('\nMD docking refinement with')
print(f'ligand: {ligand}\nreceptor: {receptor_name} \nrun {d_run}')

# load smiles from here:
smiles = np.loadtxt(smiles_file, dtype=str, comments=None)  # ensure `#` is not interpreted as comment

# write smiles to disc
with open(f'/tmp/sm-{os.getpid()}.smi', 'w') as f:
    f.write([s[1] for s in smiles if s[0] == ligand][0] + ' ' + ligand)

# load topology from smiles
molecule = openforcefield.topology.Molecule(f'/tmp/sm-{os.getpid()}.smi')
molecule.generate_unique_atom_names()  # necessary, otherwise empty strings

# generate conformer to start MD from
molecule.generate_conformers(n_conformers=1)

# load docked pose
pose = mdtraj.load_pdb(_dir + f'dockings_{d_run}/{ligand}/{receptor_name}_bestpose.pdb',
                       frame=0)

positions = pose.xyz[0] * u.nanometer
topology = pose.top.to_openmm()

# template generation is slow and native openforcefields impl doesn't work
# DIY cache

if os.path.exists(f'{cache_dir}/{molecule.name}/drug_system.p'):
    print('found drug template')
    system = pickle.load(open(f'{cache_dir}/{molecule.name}/drug_system.p', 'rb'))
else:
    print('creating fresh drug template')
    os.makedirs(f'{cache_dir}/{molecule.name}', exist_ok=True)
    tmplt_gen = openmmforcefields.generators.SMIRNOFFTemplateGenerator(
        forcefield='openff-1.1.0')

    tmplt = tmplt_gen.generate_residue_template(molecule)
    system = tmplt_gen.get_openmm_system(molecule)

    pickle.dump(system, open(f'{cache_dir}/{molecule.name}/drug_system.p', 'wb'))

    with open(f'{cache_dir}/{molecule.name}/drug_ff.xml', 'w') as f:
        f.write(tmplt)

# convert to openMM topology
t = openforcefield.topology.Topology()
t.add_molecule(molecule)
topology = t.to_openmm()

# extract graph from openMM topology
# only use heavy atoms because hydrogens are messed up in docking pose
graph_molecule = nx.Graph()

for bond in topology.bonds():
    a, b = bond.atom1, bond.atom2
    if 'H' in (a.element.symbol, b.element.symbol):
        continue

    graph_molecule.add_node(a.index, type=a.element.symbol)
    graph_molecule.add_node(b.index, type=b.element.symbol)
    graph_molecule.add_edge(a.index, b.index)

# graph for docking pose
# infer bond graph from distance matrix
# good god...
Dmat1 = np.sum(((pose.xyz[0, ..., None] - pose.xyz[0, ..., None].T)) ** 2., axis=1) ** (0.5)
Dmat1[np.diag_indices(Dmat1.shape[0])] = 99


# match molecular graphs on heavy atoms
# Moritz' code
class MolGraphMatcher(iso.GraphMatcher):
    def semantic_feasibility(self, G1_node, G2_node):
        return G1attr[G1_node] == G2attr[G2_node]


G1attr = nx.get_node_attributes(graph_molecule, 'type')

# distance cutoffs strictly need to take into account atom types
# mostly works with 0.2 nm though; brute force others by increasing cutoff
# this is safe because graph isomorphism is checked later
bondgraph_distcutoff = [0.2, 0.205, 0.21, 0.215, 0.22, 0.225]
for cutoff in bondgraph_distcutoff:
    bondsx = np.where(Dmat1 < cutoff)
    graph_pose = nx.Graph()

    for a_idx, b_idx in zip(*bondsx):
        a, b = pose.top.atom(a_idx), pose.top.atom(b_idx)
        if 'H' in (a.element.symbol, b.element.symbol):
            continue
        graph_pose.add_node(a.index, type=a.element.symbol)
        graph_pose.add_node(b.index, type=b.element.symbol)
        graph_pose.add_edge(a.index, b.index)

    G2attr = nx.get_node_attributes(graph_pose, 'type')
    GM = MolGraphMatcher(graph_molecule, graph_pose)
    if GM.is_isomorphic():
        print(f'used {cutoff} as dist cutoff for inferring bonds from docking pose')
        break

# check if graphs are isomorphic
assert GM.is_isomorphic(), 'molecular graphs not isomorphic'
GM.match()

# check if elements match
for i, j in GM.mapping.items():
    assert list(topology.atoms())[i].element.symbol == pose.top.atom(j).element.symbol, \
        'sorted atom elements don`t match'

# add constraint force to put molecule at docking pose position
force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
force.addGlobalParameter("k", 500000.0)
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")

for i, j in GM.mapping.items():
    force.addParticle(int(i), mm.Vec3(*pose.xyz[0, j]))
system.addForce(force)


def get_rmsd(pos=None):
    if pos is None:
        state = context.getState(getPositions=True)

        pos = state.getPositions(asNumpy=True)
    rmsd = np.sqrt(np.square(np.linalg.norm((pos / u.nanometer)[list(GM.mapping.keys())] - \
                                            pose.xyz[0, list(GM.mapping.values())], axis=1)).sum() / molecule.n_atoms)
    return rmsd


print('starting drug rmsd:', get_rmsd(pos=molecule.conformers[0]))

# Create OpenM Context
platform = mm.Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'single')
integrator = mm.LangevinIntegrator(temperature_K * u.kelvin,
                                   1.0 / u.picoseconds, .5 * u.femtosecond)

context = mm.Context(system, integrator, platform)
context.setPositions(molecule.conformers[0])

# get 'molecule' to mirror docking pose
# local energy minimization is almost enough to fix it
mm.LocalEnergyMinimizer.minimize(context, tolerance=1e-16)
integrator.step(100)

state = context.getState(getPositions=True)
drug_pos = state.getPositions()


def get_intramolecule_mindist(drug_pos):
    test_pos = np.array(drug_pos.value_in_unit(u.nanometer))
    Dmat2 = np.sum((test_pos[..., None] - test_pos[..., None].T) ** 2., axis=1) ** (0.5)
    Dmat2[np.diag_indices(Dmat2.shape[0])] = 99
    return Dmat2.min()


# check for atom overlap here; can be an artifact of a too strong constraint force
if get_intramolecule_mindist(drug_pos) <= 0.02:
    print('ligand conformation damaged. Retry cautiously.')

    # re-set to original position
    context.setPositions(molecule.conformers[0])

    # increase force constant slowly, do not go as high as before because that may break it again
    for _k in [.1, 1, 10, 1000, 10000]:
        context.setParameter('k', _k)
        mm.LocalEnergyMinimizer.minimize(context, tolerance=1e-16)
        integrator.step(1000)

    state = context.getState(getPositions=True)
    drug_pos = state.getPositions()

if get_intramolecule_mindist(drug_pos) <= 0.02:
    raise RuntimeError('System could not be equilibrated without atom overlap. Abort.')

print('optimized rmsd:', get_rmsd())

#### PHASE 2
# add receptor and general MD of docked pose
# load repaired receptor structure with flexible residues
receptor = app.PDBFile(_dir + f'dockings_{d_run}/{ligand}/{receptor_name}_receptor.pdb')

# for openMM: define residue variants (only cystein bridges necessary)
cbridge_resids = np.array([155, 171, 182, 210, 26, 42]) - 1
res_variants = ['CYX' if (res.index in cbridge_resids) else None for res in receptor.topology.residues()]

m = app.Modeller(receptor.topology, receptor.positions)

# explicitly create cystein bridges
m.addHydrogens(variants=res_variants)

# add drug
m.add(topology, drug_pos)

# check if disulfid bridges are intact
ss_bonds = [frozenset([25, 41]), frozenset([181, 209]), frozenset([154, 170])]

cross_bonds = set()
for b in m.topology.bonds():
    if abs(b.atom1.residue.index - b.atom2.residue.index) > 2:
        cross_bonds.add(frozenset([b.atom1.residue.index, b.atom2.residue.index]))

for ss_bond in ss_bonds:
    assert (ss_bond in cross_bonds)

# force field consists of amber14 and previously written small molecule template
forcefield = app.ForceField('amber14/protein.ff14SB.xml',
                            'amber14/tip3p.xml',
                            f'{cache_dir}/{molecule.name}/drug_ff.xml')

# add solvent
m.addSolvent(forcefield, boxSize=[boxsize / u.nanometer for _ in range(3)],
             positiveIon='Na+', negativeIon='Cl-',
             ionicStrength=0.1 * u.molar)

# atomids for saving etc
ntop_mdtraj = mdtraj.Topology.from_openmm(m.topology)

lig_rec_atomids = ntop_mdtraj.select(f'chainid 0 1')
lig_atomids = ntop_mdtraj.select(f'chainid 1')


def rmsd_to_reference(pos, ref, atomsel):
    """
    compute root mean square deviation to a reference
    
    :param pos, ref: numpy array (in units nanometer)
    :param atomsel: list-like of atom indices to consider
    return: float
    """

    rmsd = np.sqrt(np.square(np.linalg.norm(pos[atomsel] - \
                                            ref[atomsel], axis=1)).sum() / len(atomsel))

    return rmsd


outpath = _dir + f'dockings_{d_run}/{ligand}/'

# pre-equilibration: normal hydrogen mass. local energy minimizer + NVT (2fs) + NPT (2fs)
# with position restraint force on heavy atoms to mirror docking pose

system = forcefield.createSystem(m.topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.0 * u.nanometers, rigidWater=True,
                                 ewaldErrorTolerance=0.0005, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature_K * u.kelvin, 1.0 / u.picoseconds,
                                   integrator_timestep_ps_equil * u.picosecond)
integrator.setConstraintTolerance(0.00001)

# add restraint force to keep molecule at docking pose position
force = mm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
force.addGlobalParameter("k", 5000.0)
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")

for i in range(m.topology.getNumAtoms()):
    if system.getParticleMass(i) / u.dalton > 2:
        force.addParticle(i, m.positions[i])

cforce_id = system.addForce(force)
app.PDBFile.writeFile(m.topology, m.positions,
                      file=open(
                          os.path.join(outpath, f'{receptor_name}_init.pdb'), 'w'))

simulation = app.Simulation(m.topology, system, integrator, platform)
simulation.context.setPositions(m.positions)

print('mini')
simulation.minimizeEnergy()

# NVT equil
print('constraint nvt')
simulation.step(int(simulation_time_ns_equil_nvt / (integrator_timestep_ps_equil * 1e-3)))

# add barostat
barostat = mm.MonteCarloBarostat(1. * u.bar, temperature_K * u.kelvin)
system.addForce(barostat)
simulation.context.reinitialize(preserveState=True)

print('constraint npt')
simulation.step(int(simulation_time_ns_equil_npt / (integrator_timestep_ps_equil * 1e-3)))

state = simulation.context.getState(getPositions=True, getVelocities=True)
newpositions = state.getPositions()
newvelocities = state.getVelocities()
pbox_after_nptequil = state.getPeriodicBoxVectors()

del system, integrator, simulation

# save solvated pose for later
m.topology.setPeriodicBoxVectors(pbox_after_nptequil)
app.PDBFile.writeFile(m.topology, newpositions,
                      file=open(
                          os.path.join(outpath, f'{receptor_name}_docked2_equil_solv.pdb'), 'w'))

_tmp = mdtraj.load(os.path.join(outpath, f'{receptor_name}_docked2_equil_solv.pdb')).atom_slice(lig_rec_atomids)
make_protein_complex_whole(_tmp)
_tmp.save_pdb(os.path.join(outpath, f'{receptor_name}_docked2_equil.pdb'))
del _tmp

# production run with heavy hydrogens
system = forcefield.createSystem(m.topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.0 * u.nanometers, rigidWater=True,
                                 ewaldErrorTolerance=0.0005, constraints=app.HBonds,
                                 hydrogenMass=4. * u.amu)

barostat = mm.MonteCarloBarostat(1. * u.bar, temperature_K * u.kelvin)
system.addForce(barostat)

# save topology and system in compressed format
# not sure if there are complex problems inferring with creating the system. 
# it should be the same for every pose.
# this could be done more efficiently
with bz2.BZ2File(f'{cache_dir}/{molecule.name}/{receptor_name}_system.xml.bz2', 'w') as f:
    xml = mm.XmlSerializer.serialize(system)
    f.write(xml.encode('ascii'))

# simulation setup
integrator = mm.LangevinIntegrator(temperature_K * u.kelvin, 1.0 / u.picoseconds,
                                   integrator_timestep_ps * u.picosecond)
integrator.setConstraintTolerance(0.00001)

simulation = app.Simulation(m.topology, system, integrator,
                            platform)

# re-start from previous position
simulation.context.setPositions(newpositions)
simulation.context.setVelocities(newvelocities)
simulation.context.setPeriodicBoxVectors(*pbox_after_nptequil)
simulation.context.applyConstraints(1e-12)  # HBond constraints

simulation.minimizeEnergy()

# output reporters
ener_outfile = os.path.join(outpath, f'{receptor_name}-{stamp}-energy.csv')
simulation.reporters.append(app.StateDataReporter(ener_outfile,
                                                  int(float(save_traj_ps) / integrator_timestep_ps),
                                                  time=True, potentialEnergy=True, kineticEnergy=True,
                                                  totalEnergy=True, temperature=True))

prot_outfile = os.path.join(outpath, f'{receptor_name}-{stamp}-protein.dcd')
simulation.reporters.append(DCDReporter(prot_outfile,
                                        int(float(save_traj_ps) / integrator_timestep_ps),
                                        atomSubset=lig_rec_atomids))

ref_pos = np.array(m.positions / u.nanometer)

# trj object for mirroring drug into binding pocket
trj = mdtraj.Trajectory(ref_pos[None], ntop_mdtraj,
                        unitcell_lengths=np.array(pbox_after_nptequil / u.nanometer)[np.diag_indices(3)],
                        unitcell_angles=[90 for _ in range(3)])
ref = copy.deepcopy(trj)

# time for benchmarking
starting_time = time.time()
print('production MD')

# run simulation
simulation.step(int(simulation_time_ns / (integrator_timestep_ps * 1e-3)))

# benchmarking time
time_diff = time.time() - starting_time

# print benchmark
ns_per_day = float(simulation_time_ns) / time_diff * 60. * 60. * 24.
print('ns per day: ', ns_per_day)

state = simulation.context.getState(getPositions=True)
pbox = state.getPeriodicBoxVectors()
m.topology.setPeriodicBoxVectors(pbox)

outname_solv = os.path.join(outpath, f'{receptor_name}-{stamp}_{simulation_time_ns}ns_solv.pdb')
app.PDBFile.writeFile(m.topology, state.getPositions(),
                      file=open(outname_solv, 'w'))

# save protein ligand coordinates only; map to closest ligand image
out_struct = mdtraj.load(outname_solv).atom_slice(lig_rec_atomids)
make_protein_complex_whole(out_struct)
out_struct.save_pdb(os.path.join(outpath, f'{receptor_name}_{stamp}_{simulation_time_ns}ns.pdb'))

# compress solvated frame
with open(outname_solv, 'r') as infile, bz2.BZ2File(outname_solv + '.bz2', 'w') as outfile:
    outfile.write(infile.read().encode('ascii'))
os.unlink(outname_solv)

# post-processing
trj = mdtraj.load(prot_outfile, top=ntop_mdtraj.subset(lig_rec_atomids))
make_protein_complex_whole(trj)
trj.save_xtc(prot_outfile.replace('.dcd', '.xtc'))

os.unlink(prot_outfile)
