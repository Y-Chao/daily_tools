#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"


import os
import sys

import numpy as np
from ase.geometry import get_distances
from ase.io import read, write
from ase.mep import NEB

if len(sys.argv) < 4:
    raise RuntimeError(
        "Usage: python run_neb.py init.traj final.traj num_images [skip_atom_index]"
    )

init_atoms = read(sys.argv[1])
final_atoms = read(sys.argv[2])
num = int(sys.argv[3])

if len(sys.argv) == 5:
    skip_atoms = [int(sys.argv[4])]
else:
    skip_atoms = []


def pos_swap(atoms, ind1, ind2):
    # swap position of two atoms
    x1, y1, z1 = [atoms[ind1].x, atoms[ind1].y, atoms[ind1].z]
    x2, y2, z2 = [atoms[ind2].x, atoms[ind2].y, atoms[ind2].z]

    atoms[ind1].x = x2
    atoms[ind1].y = y2
    atoms[ind1].z = z2

    atoms[ind2].x = x1
    atoms[ind2].y = y1
    atoms[ind2].z = z1

    return atoms


def get_closest(ref_atoms, atoms, ind, mic=True):
    # find the index of the closest atom between two states
    # making sure that the symbol is the same
    if not mic:
        pbc = (False, False, False)
    else:
        pbc = (True, True, True)
    dists = []
    for atom in atoms:
        if atom.symbol != ref_atoms[ind].symbol:
            continue
        dist = get_distances(
            p1=(atom.x, atom.y, atom.z),
            p2=(ref_atoms[ind].x, ref_atoms[ind].y, ref_atoms[ind].z),
            cell=atoms.cell,
            pbc=pbc,
        )[1][0][0]
        dists.append((atom.index, dist))
    dists.sort(key=lambda x: x[-1])
    return dists[0][0]


def reindex_atoms(ref_atoms, reindex_atoms, manual_skip_atoms=[]):
    # used to reindex atoms in reindex_atoms to match those in
    # ref_atoms. Necessary for e.g. NEB interpolation.
    for atom in reindex_atoms:
        if atom.index in manual_skip_atoms:
            continue
        closest_ind = get_closest(ref_atoms, reindex_atoms, atom.index)
        if atom.index == closest_ind:
            continue
        else:
            pos_swap(reindex_atoms, closest_ind, atom.index)
    return reindex_atoms


if len(skip_atoms) > 0:
    final_reindex_atoms = reindex_atoms(
        init_atoms, final_atoms, manual_skip_atoms=[skip_atoms]
    )
    final_atoms = final_reindex_atoms

_, dist = get_distances(
    init_atoms.positions,
    final_atoms.positions,
    cell=init_atoms.cell,
    pbc=init_atoms.pbc,
)

rmse_dist = np.sqrt(np.sum([dist[i][i] ** 2 for i in range(len(init_atoms))], axis=0))

if num == 0:
    num = int(rmse_dist / 0.8)

print(f"The rmse_dist is {rmse_dist}")

images = [init_atoms] + [init_atoms.copy() for i in range(num)] + [final_atoms]

neb = NEB(images)

neb.interpolate(method="idpp", mic=True)


write("neb.traj", images, format="traj")

for i, atoms in enumerate(images):
    os.mkdir(f"{i:02d}")
    write(f"{i:02d}/POSCAR", atoms)
