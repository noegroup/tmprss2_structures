#!/usr/bin/env python

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.


import numpy as np


def make_protein_complex_whole(trj):
    """
    For systems with periodic boundary conditions, finds the image of a molecule
    that is closest to the 0th chain. Currently set up for 2 chains.
    Caution: operates inplace

    :param trj: MDTraj.Trajectory that contains multiple chains
    :return: MDTraj.Trajectory
    """
    chains = list(trj.top.chains)[:2]
    n_chains = len(chains)

    def get_coms(frame, top):
        coms = np.zeros((n_chains, 3))
        for chain_n in range(n_chains):
            coms[chain_n] = frame[[a.index for a in chains[chain_n].atoms], :].mean(axis=0)
        return coms

    for frame_n in range(trj.n_frames):
        pbox = trj.unitcell_vectors[frame_n]

        coms = get_coms(trj.xyz[frame_n], trj.top)
        diffvecs = np.array([coms[1] - coms[0]])#, coms[2] - coms[0]])
        shift_chain_n = 1# np.linalg.norm(diffvecs, axis=1).argmax() + 1

        longest_diffvec = diffvecs[shift_chain_n-1]

        fuckitcounter = 0
        cutoff = 1.
        while (np.linalg.norm(longest_diffvec) > cutoff and fuckitcounter < 100):
            component = abs(longest_diffvec).argmax()
            sign = np.sign(longest_diffvec[component])
            shift = sign * pbox[component]

            trj.xyz[frame_n, [a.index for a in chains[shift_chain_n].atoms], :] -= shift

            coms = get_coms(trj.xyz[frame_n], trj.top)
            diffvecs = np.array([coms[1] - coms[0]])#, coms[2] - coms[0]])
            shift_chain_n = 1#np.linalg.norm(diffvecs, axis=1).argmax() + 1

            longest_diffvec = diffvecs[shift_chain_n - 1]

            fuckitcounter += 1

            if fuckitcounter > 10:
                cutoff += .5
        if fuckitcounter > 99:
            print('something didn`t work with centering the complex. procede.')
    return trj
