#!/usr/bin/env python
# coding: utf-8

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.

# parse arguments
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_list', type=str,
                    help='text file that contains lines with \"d_run ligand receptor\"')
parser.add_argument('-p', '--pre_computed', type=str,
                    help='pickle file with pre-computed sasa, reactive distances etc; all_results_merged.pickle')
args = parser.parse_args()

print(' - COMPUTE REACTIVE MD-SCORES -')

import numpy as np
import pickle
import pandas as pd

lig_rec_pairs = np.loadtxt(args.input_list, dtype=str, comments=None)

lig_rec_pairs_pd = pd.DataFrame(lig_rec_pairs, columns=['d_run', 'lig', 'rec'])


residue_groups = {'S1':np.concatenate([np.arange(180, 187), np.arange(204, 210)])-1,
                  'hydroph':np.arange(24, 27)-1,
                  'flex':np.arange(41, 46)-1}
sasa_residues = np.concatenate(list(residue_groups.values()))

# results have been computed on allegro in 'virtual_screening/md_refinement/analysis_XXX/run.py'
# load subset corresponding to ligands in selection
_l, _n = np.unique(lig_rec_pairs[:, 1], return_index=True)
ligands_sorted = _l[_n.argsort()]


all_results = {lig:{} for lig in ligands_sorted}
n_results = 0
failed_cmpnds = []

result_part = pickle.load(open(args.pre_computed, 'rb'))

for lig, ligdict in result_part.items():
    if ligdict is None or ligdict == {}:
        continue

    if lig not in ligands_sorted:
        continue

    incoming_results_n = sum([sum([len(__v) if __v is not None else 0 for __v in _v.values()]) for _v in ligdict.values() if _v is not None])
    if incoming_results_n == 0:
        continue

    for rec, recdict in ligdict.items():
        if recdict is None:
            continue

        if rec not in lig_rec_pairs[lig_rec_pairs[:, 1] == lig][:, 2]:
            continue


        for run, rundict in recdict.items():
            if rundict is None:
                continue


            # exclude double entries...
            # deactivated dec 15, 2022
            #if run.split('-')[0] in [s.split('_202')[0] for s in all_results[lig].keys()]:
            #    continue

            newlabel = run.split('-')[0]+'_'+'-'.join(run.split('-')[1:3])
            all_results[lig][newlabel] = rundict

            n_results += 1        

# gather some stats about input data (failed compounds, missing lig/rec pairs)
print(f'loaded {n_results} of {len(ligands_sorted)*10} lig/rec pairs')
failed_cmpnds = []
for lig, ligdict in all_results.items():
    incoming_results_n = len([a for a in ligdict if ligdict[a] is not None])
    if incoming_results_n == 0:
        failed_cmpnds.append(lig)

failed_cmpnds = np.unique(failed_cmpnds)

if len(failed_cmpnds) > 0:
    print(f'{len(failed_cmpnds)} compounds are complete failures; writing to failed_cmpnds.txt.')
    np.savetxt('failed_cmpnds.txt', failed_cmpnds, fmt='%s')


n_md_per_compound = []
something_missing = []
for lig, ligdict in all_results.items():
    n_md_per_compound.append(len(ligdict))
    if len(ligdict) < 10 and len(ligdict) > 0:
        something_missing.append(lig)


bc = np.bincount(n_md_per_compound)
un = np.unique(n_md_per_compound)
print('\nTRAJ STATS:')
for _length, _count in zip(un, bc[bc!=0]):
    print(f'{_count} \tcompounds with {_length} trajectories')

# check for missing compounds + write into file for follow-up submission
something_missing_writeout = []
for lig in something_missing:
    ligdict = result_part[lig]
    for rec, recdict in ligdict.items():
        if recdict == {} or recdict is None or \
                recdict[list(recdict.keys())[0]] is None or \
                'n_heavy' not in recdict[list(recdict.keys())[0]].keys():
            something_missing_writeout.append([lig_rec_pairs_pd[lig_rec_pairs_pd.lig == lig]['d_run'].iloc[0],
                                               lig,
                                               rec])

if len(failed_cmpnds) > 0:
    with open(args.input_list, 'r') as infile, open('missing_cmpnds.txt', 'w') as outfile:
        for line in infile.readlines():
            if line.split()[1] in failed_cmpnds:
                outfile.write(line)

with open('missing_cmpnds.txt', 'a') as outfile:
    for line in something_missing_writeout:
        outfile.write(' '.join(line) + '\n')



# re-arrange data into other structure...
# this is an atavism of previous code, could be more efficient
def add_observable(observable_dict, lig, observable):
    observable_lig = observable_dict.get(lig, [])
    observable_lig.append(observable)
    observable_dict[lig] = observable_lig

distance_asp180 = {}
distance_reactive = {}
contacts = {k:{} for k in residue_groups.keys()}
dsasas = {}
n_heavy = {}
has_reactive_group = {}
rec_sorting = {}

for lig in ligands_sorted:
    
    lig_dict = all_results[lig]
    rec_sorting[lig] = ['_'.join(k.split('_')[:4]) for k in lig_dict.keys() if lig_dict[k] is not None]
    
    for rec, res in lig_dict.items():
        if res is not None:
            add_observable(n_heavy, lig, res['n_heavy'])
            add_observable(dsasas, lig, res['dsasa'])
            add_observable(distance_asp180, lig, res['dist_asp180'])
            if res['dist_reactive_arr'].shape == ( ):
                res['dist_reactive_arr'] = res['dist_reactive_arr'][None]
            add_observable(distance_reactive, lig, res['dist_reactive_arr'])
            add_observable(has_reactive_group, lig, np.isfinite(res['dist_reactive']))
            
            
            for name, resgroup in residue_groups.items():
                add_observable(contacts[name], lig, res[f'contacts_{name}'])
    
    if lig in n_heavy.keys():
        n_heavy[lig] = np.unique(n_heavy[lig])[0]
        if not len(np.unique(has_reactive_group[lig])) == 1:
            print(f'failed to label reactive/nonreactive {lig} -> non-reactive')
            print('... maybe problematic:', np.array(rec_sorting[lig])[~np.array(has_reactive_group[lig])])
            has_reactive_group[lig] = False
            
            distance_reactive[lig] = [np.zeros_like(a) + np.inf for a in distance_asp180[lig]]
        else:
            has_reactive_group[lig] = np.unique(has_reactive_group[lig])[0]
            
        if not has_reactive_group[lig]:
            distance_reactive[lig] = [np.zeros_like(a) + np.inf for a in distance_asp180[lig]]
            
        for a, b in zip(distance_asp180[lig], distance_reactive[lig]):
            assert a.shape == b.shape


regionalsasa = {}
# inline for-loops to fix different length sasa...

for k, reg in residue_groups.items():
    ind = np.where(np.in1d(sasa_residues, reg))[0]
    #regionalsasa[k] = {k:np.asarray(v)[:, :, ind].sum(axis=2) for k, v in dsasas.items()}
    regionalsasa[k] = {k:np.array([np.array(_v)[:, ind].sum(axis=1) for _v in v]) for k, v in dsasas.items()}

#regionalsasa['all'] = {k:np.asarray(v).sum(axis=2) for k, v in dsasas.items()}
regionalsasa['all'] = {k:np.array([np.array(_v).sum(axis=1) for _v in v]) for k, v in dsasas.items()}

#regionalsasa['flexhydro'] = {k:np.asarray(v)[:, :, 13:].sum(axis=2) for k, v in dsasas.items()}
regionalsasa['flexhydro'] = {k:np.array([np.array(_v)[:, 13:].sum(axis=1) for _v in v]) for k, v in dsasas.items()}

# ensure that scores are comparable between different datasets
mean_reactive_distance = 0.827 # previously computed average

def get_sorted_dataframes(scores, return_topn=None):
    """
    helper function to sort ligands along scores
    """
    score_means = []


    for lig in all_results.keys():
        if lig in failed_cmpnds.tolist():
            continue
        n_rec = 0
        for rec in all_results[lig].keys():
            d_run = lig_rec_pairs_pd[lig_rec_pairs_pd.lig == lig]['d_run'].to_numpy()[0]
            try:
                scores[lig][n_rec].mean()
            except:
                print(lig, n_rec)
            score_means.append([f'{d_run}/{lig}', 
                                rec, 
                                scores[lig][n_rec].mean(),
                               regionalsasa['S1'][lig][n_rec].mean(),
                               regionalsasa['flexhydro'][lig][n_rec].mean(),
                               distance_asp180[lig][n_rec].mean(),
                               distance_reactive[lig][n_rec].mean()
                               ])
            n_rec += 1

    score_means = pd.DataFrame(score_means, 
                               columns=['ligand', 'receptor', 'score', 
                                        'sasa_s1', 'sasa_flexhydro', 'dist_asp180', 
                                        'dist_react'
                                       ])
    
    top3 = pd.DataFrame()
    topn = pd.DataFrame()
    scored_md = []
    for ligand in np.unique(score_means.ligand.to_numpy()):

        ligand_frame = score_means[score_means.ligand == ligand]
        if return_topn is not None:
            lf2 = ligand_frame.copy()
            lf2 = lf2.sort_values('score', ascending=False)[:return_topn]
            topn = topn.append(lf2)
            
        ligand_frame = ligand_frame.sort_values('score', ascending=False)[:3]
        
        top3 = top3.append(ligand_frame)            

        scored_md.append([ligand_frame['ligand'].iloc[0], 
                          ligand_frame['score'].mean(),
                         has_reactive_group[ligand.split('/')[1]]])

    top3.sort_index(inplace=True)
    if return_topn is not None:
        topn.sort_index(inplace=True)
        
    scored_md = pd.DataFrame(scored_md, columns=['ligand', 'mscore', 'isReact'])

    scored_md.sort_values('mscore', ascending=False, inplace=True)
    scored_md.reset_index(inplace=True, drop=True)
    
    if return_topn:
        return top3, scored_md, topn
    else:
        return top3, scored_md


# score taking into account reactivity
# if a drug has no reactive center, the reactive distance is 
# defined as the mean reactive distance of all reactive drugs
# to keep not bias it, still probably take the other score if 
# there's both, covalent and non-covalent

print('computing scores')
scores_nodock = {}
for lig in ligands_sorted:
    if lig in failed_cmpnds.tolist():
        continue
    scores_nodock[lig] = []
    assert len(regionalsasa['S1'][lig]) == len(regionalsasa['flexhydro'][lig])
    assert len(regionalsasa['S1'][lig]) == len(distance_asp180[lig])
    assert len(regionalsasa['S1'][lig]) == len(distance_reactive[lig])
    assert len(regionalsasa['S1'][lig]) == len(rec_sorting[lig])
    
    for x, y, z1, _z2, rec in zip(regionalsasa['S1'][lig], 
                       regionalsasa['flexhydro'][lig], 
                       distance_asp180[lig],
                        distance_reactive[lig], 
                        rec_sorting[lig]):
        
        is_reactive = np.all(np.isfinite(_z2))
        z2 = _z2 if is_reactive else mean_reactive_distance
        
        scores_nodock[lig].append(x * y / (n_heavy[lig]**(4/3) * z1**2 * z2**2))

print('sorting to mean score')
top3_nodock, scored_md_nodock, top5_nodock = get_sorted_dataframes(scores_nodock, return_topn=5)

print('saving files scoredMD.csv, top3.csv, top5.csv')
scored_md_nodock.to_csv('scoredMD.csv')
top3_nodock.to_csv('top3.csv')
top5_nodock.to_csv('top5.csv')

print('Finished; exiting')
