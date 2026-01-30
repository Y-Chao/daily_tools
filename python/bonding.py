#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"


"""
Define the bonding related functions.

Electronic negativity strategy:
    Refer to Matterviz developed by Janosh.
    * Electronegativity-based bonding with chemical preferences.
    * This algorithm considers electronegativity differences between atoms, metal/nonmetal
    * properties, and distance to determine bond strength. Bonds are only created if the
    * computed strength exceeds the strength_threshold parameter (default: 0.3).

"""

import math as m
from typing import Optional, Union

import numpy as np
from ase.data import covalent_radii
from pymatgen.core import Molecule, Site, Structure


def build_spatial_grid(
    sites: list[Site], cell_size: float
) -> dict[tuple[int, int, int], list[int]]:
    """
    Build a spatial grid for efficient neighbor searching.

    Parameters:
    -----------
    sites: list of Site objects
        The atomic sites in the structure.
    cell_size: float
        The size of each grid cell.
    Returns:
    --------
    grid: dict [tuple[int, int, int], list[int]]
        A dictionary mapping grid cell indices to lists of site indices.
    """
    grid = {}
    for i, site in enumerate(sites):
        x, y, z = site.coords
        ix = int(x // cell_size)
        iy = int(y // cell_size)
        iz = int(z // cell_size)
        key = (ix, iy, iz)
        if key not in grid:
            grid[key] = []
        grid[key].append(i)
    return grid


def setup_spatial_grid(
    sites: list[Site], cutoff: float
) -> tuple[Union[dict[tuple[int, int, int], list[int]], None], Union[float, None]]:
    """
    Setup the spatial grid for neighbor searching.

    Parameters:
    -----------
    sites: list of Site objects
        The atomic sites in the structure.
    cutoff: float
        The maximum cutoff distance for bonding.
    Returns:
    --------
    spatial_grid: dict [tuple[int, int, int], list[int]]
        The spatial grid for neighbor searching.
    cell_size: float or None
        The cell size used for the grid, or None if grid is not used.
    """
    if len(sites) > 50:
        return None, None
    return build_spatial_grid(sites, cutoff), cutoff


def get_neighbors_from_grid(
    pos: np.ndarray,
    grid: dict[tuple[int, int, int], list[int]],
    cell_size: float,
) -> list[int]:
    """
    Get neighboring site indices from the spatial grid.
    Parameters:
    -----------
    pos: tuple[float, float, float]
        The position to search around.
    grid: dict[tuple[int, int, int], list[int]]
        The spatial grid.
    cell_size: float
        The size of each grid cell.
    Returns:
    --------
    neighbors: list[int]
        List of neighboring site indices.
    """
    [cx, cy, cz] = [
        m.floor(pos[0] / cell_size),
        m.floor(pos[1] / cell_size),
        m.floor(pos[2] / cell_size),
    ]
    neighbors = []
    # Iterate over neighboring cells
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                cell_key = (cx + dx, cy + dy, cz + dz)
                if cell_key in grid:
                    neighbors.extend(grid[cell_key])
    return neighbors


def get_neighboring_candidates(
    pos: np.ndarray,
    sites: list[Site],
    spatial_grid: Union[dict[tuple[int, int, int], list[int]], None],
    cell_size: Union[float, None],
) -> list[Site]:
    """
    Get neighboring candidate sites from the spatial grid.
    """
    if spatial_grid is None or cell_size is None:
        return [i for i in range(len(sites))]
    return get_neighbors_from_grid(pos, spatial_grid, cell_size)


# For bond/cylinder rendering between two atoms
def compute_bond_transform(pos1: np.ndarray, pos2: np.ndarray) -> np.ndarray:
    """
    Compute the transformation matrix for the bond between two positions.
    """
    dx, dy, dz = pos2[0] - pos1[0], pos2[1] - pos1[1], pos2[2] - pos1[2]
    height = m.hypot(dx, dy, dz)  # Euclidean distance

    if height < 1e-10:
        return np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=np.float32
        )

    dir_x, dir_y, dir_z = dx / height, dy / height, dz / height
    [m00, m01, m02, m10, m11, m12, m20, m21, m22] = [0, 0, 0, 0, 0, 0, 0, 0, 0]

    # Special case: bond pointing straight up (+Y)
    if abs(dir_y - 1.0) < 1e-10:
        [m00, m01, m02, m10, m11, m12, m20, m21, m22] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    elif abs(dir_y + 1.0) < 1e-10:
        [m00, m01, m02, m10, m11, m12, m20, m21, m22] = [1, 0, 0, 0, -1, 0, 0, 0, 1]
    else:
        # General case: construct orthonormal basis (right, dir, up)
        # Right vector: perpendicular to dir in XZ plane
        [rx, rz] = [-dir_z, dir_x]
        right_length = m.hypot(rx, rz)
        [right_x, right_z] = [rx / right_length, rz / right_length]
        # Up vector: cross product of dir and right
        [up_x, up_y, up_z] = [
            dir_y * right_z,
            dir_z * right_x - dir_x * right_z,
            -dir_y * right_x,
        ]
        [m00, m01, m02, m10, m11, m12, m20, m21, m22] = [
            right_x,
            dir_x,
            up_x,
            0,
            dir_y,
            up_y,
            right_z,
            dir_z,
            up_z,
        ]

    # Position at midpoint between pos1 and pos2
    [px, py, pz] = [
        0.5 * (pos1[0] + pos2[0]),
        0.5 * (pos1[1] + pos2[1]),
        0.5 * (pos1[2] + pos2[2]),
    ]
    return np.array(
        [
            [m00, m10, m20, 0],
            [m01 * height, m11 * height, m21 * height, 0],
            [m02, m12, m22, 0],
            [px, py, pz, 1],
        ],
        dtype=np.float32,
    )


def electroneg_ratio(
    structure: Structure,
    electroneg_threshold: float = 1.7,  # Maximum electronegativity difference for bonding
    max_distance_ratio: float = 2.0,  # Factor for the sum of covalent radii
    min_bond_distance: float = 0.4,  # Minimum bond distance in Angstrom
    metal_metal_penalty: float = 0.5,  # Strength penalty for metal-metal bonding
    metal_nonmetal_bonus: float = 1.5,  # Strength bonus for metal-nonmetal bonding
    similarity_x_bonus: float = 1.2,  # Bonus for similar electronegativity
    same_element_penalty: float = 0.5,  # Penalty for same element bonding
    strength_threshold: float = 0.3,  # Minimum bonding strength to include in results
):
    """
    Calculate the bonding status based on the electronegativity ratio.
    """

    min_dist_sq = min_bond_distance**2
    sites = structure.sites
    num_sites = len(sites)
    if num_sites < 2:
        return []

    prop = []
    for site in sites:
        prop.append(
            {
                "element": site.specie.symbol,
                "electroneg": site.specie.X,
                "is_metal": site.specie.is_metal,
                "is_nonmetal": not site.specie.is_metal,
                "covalent_radius": covalent_radii[site.specie.Z],
            }
        )

    max_radius = max(covalent_radii)
    max_cutoff = (
        max_radius * 2 * max_distance_ratio
    )  # Maximum distance to consider bonding
    spatial, cell_size = setup_spatial_grid(sites, max_cutoff)
    bonds = []
    closet = {}

    for i in range(num_sites):
        xi, yi, zi = sites[i].coords
        pi = prop[i]

        site_bonds = []
        for j in get_neighboring_candidates(sites[i].coords, sites, spatial, cell_size):
            if j == i:
                continue
            xj, yj, zj = sites[j].coords
            pj = prop[j]
            dx, dy, dz = xi - xj, yi - yj, zi - zj
            dist_sq = dx * dx + dy * dy + dz * dz
            dist = m.sqrt(dist_sq)
            if (
                dist_sq < min_dist_sq
                or not pi["covalent_radius"]
                or not pj["covalent_radius"]
            ):
                continue
            expected = pi["covalent_radius"] + pj["covalent_radius"]

            if dist > expected * max_distance_ratio:
                continue

            en_diff = abs(pi["electroneg"] - pj["electroneg"])
            en_ratio = en_diff / (pi["electroneg"] + pj["electroneg"])

            bond_strength = 1.0
            if pi["is_metal"] and pj["is_metal"]:
                bond_strength *= metal_metal_penalty
            elif (pi["is_metal"] and pj["is_nonmetal"]) or (
                pi["is_nonmetal"] and pj["is_metal"]
            ):
                bond_strength *= metal_nonmetal_bonus
            if en_diff > electroneg_threshold:
                bond_strength *= 1.3
            elif en_diff < 0.5:
                bond_strength *= similarity_x_bonus

            dist_weight = m.exp(-((dist / expected - 1) ** 2) / 0.18)
            en_weight = 1.0 - 0.3 * en_ratio
            strength = bond_strength * dist_weight * en_weight

            if pi["element"] == pj["element"]:
                strength *= same_element_penalty

            ca = closet.get(i, float("inf"))
            cb = closet.get(j, float("inf"))

            if dist > ca:
                strength *= m.exp(-(dist / ca - 1) / 0.5)
            if dist > cb:
                strength *= m.exp(-(dist / cb - 1) / 0.5)

            if strength > strength_threshold:
                site_bonds.append(
                    {
                        "elem_1": sites[i].specie.symbol,
                        "elem_2": sites[j].specie.symbol,
                        "pos_1": sites[i].coords,
                        "pos_2": sites[j].coords,
                        "index_1": i,
                        "index_2": j,
                        "bond_length": dist,
                        "strength": strength,
                        "transform_matrix": compute_bond_transform(
                            sites[i].coords, sites[j].coords
                        ),
                    }
                )

        bonds.append(site_bonds)
    return bonds


def get_1st_shell(sites: list[Site], index: int, bonds: list[list[dict]]) -> list[int]:
    """
    Get the first shell neighbors of a given site index.
    """
    first_shell = []
    for bond in bonds[index]:
        first_shell.append(bond["index_2"])
    return first_shell


def get_2nd_shell(
    sites: list[Site], index: int, bonds: list[list[dict]], include: bool = True
) -> list[int]:
    """
    Get the second shell neighbors of a given site index.
    """
    second_shell = set()
    first_shell = get_1st_shell(sites, index, bonds)
    for neighbor_index in first_shell:
        neighbor_first_shell = get_1st_shell(sites, neighbor_index, bonds)
        for nn_index in neighbor_first_shell:
            if nn_index not in first_shell and nn_index != index:
                second_shell.add(nn_index)
    if include:
        second_shell.update(first_shell)
    return list(second_shell)


def get_3rd_shell(
    sites: list[Site], index: int, bonds: list[list[dict]], include: bool = True
) -> list[int]:
    """
    Get the third shell neighbors of a given site index.
    """
    third_shell = set()
    second_shell = get_2nd_shell(sites, index, bonds, include)
    for neighbor_index in second_shell:
        neighbor_first_shell = get_1st_shell(sites, neighbor_index, bonds)
        for nn_index in neighbor_first_shell:
            if nn_index not in second_shell and nn_index != index:
                third_shell.add(nn_index)
    if include:
        third_shell.update(second_shell)
    return list(third_shell)


def get_4th_shell(
    sites: list[Site], index: int, bonds: list[list[dict]], include: bool = True
) -> list[int]:
    """
    Get the fourth shell neighbors of a given site index.
    """
    fourth_shell = set()
    third_shell = get_3rd_shell(sites, index, bonds, include)
    for neighbor_index in third_shell:
        neighbor_first_shell = get_1st_shell(sites, neighbor_index, bonds)
        for nn_index in neighbor_first_shell:
            if nn_index not in third_shell and nn_index != index:
                fourth_shell.add(nn_index)
    if include:
        fourth_shell.update(third_shell)
    return list(fourth_shell)


def get_cn_structure(
    structure: Structure,
    index: int,
    bonds: Optional[list[list[dict]]] = None,
    order: int = 1,
) -> list[Site]:
    """
    Get the coordination number structure (first shell neighbors) of a given site index.
    Parameters:
    -----------
    structure: Structure
        The pymatgen Structure object.
    index: int
        The index of the site to get neighbors for.
    bonds: list of list of dict
        The bonding status as returned by electroneg_ratio function.
    order: int
        The order of the shell to retrieve (1 for first shell, 2 for second shell).
    Returns:
    --------
    cn_structure: Molecule
        A pymatgen Molecule object containing the neighboring sites.
    """
    if bonds is None:
        bonds = electroneg_ratio(structure)

    if order == 1:
        first_shell_indices = get_1st_shell(structure.sites, index, bonds)
        return Molecule.from_sites(
            [structure.sites[i] for i in first_shell_indices + [index]]
        )
    elif order == 2:
        second_shell_indices = get_2nd_shell(structure.sites, index, bonds)
        return Molecule.from_sites(
            [structure.sites[i] for i in second_shell_indices + [index]]
        )
    elif order == 3:
        third_shell_indices = get_3rd_shell(structure.sites, index, bonds)
        return Molecule.from_sites(
            [structure.sites[i] for i in third_shell_indices + [index]]
        )
    else:
        raise ValueError("Order must be 1, 2, or 3.")


def main():
    Structure.from_file(
        "working/surface_model/surface/miller_010_layers2_minlat10.0_symTrue/slab_1.cif"
    )


if __name__ == "__main__":
    main()

