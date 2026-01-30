
#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

from typing import Union

import numpy as np
from ase import Atom, Atoms
from ase.data import covalent_radii
from ase.data.colors import cpk_colors, jmol_colors
from matplotlib import colors as mc
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, PathPatch
from matplotlib.path import Path
from pymatgen.io.ase import AseAtomsAdaptor

from bonding import electroneg_ratio

"""
Plot atoms.

For bond plot implementation, still need to work on it.
"""


def adjust_color(color: Union[tuple, str], adjustment: float) -> tuple:
    """Adjust the brightness of a color.

    Parameters:
        color (tuple): RGB color as a tuple of three floats (0 to 1).
        adjustment (float): Amount to adjust the brightness (-1 to 1).

    Returns:
        tuple: Adjusted RGB color.
    """
    color = mc.to_rgb(color)
    adjusted_color = tuple(max(0, min(1, c + adjustment)) for c in color)
    return adjusted_color


def rot_vector(vector: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rotate a 2D vector by a given angle in degrees.

    Parameters:
        vector (np.ndarray): 2D vector to be rotated.
        angle_deg (float): Angle in degrees to rotate the vector.

    Returns:
        np.ndarray: Rotated 2D vector.
    """
    angle_rad = np.radians(angle_deg)
    rotation_matrix = np.array(
        [
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)],
        ]
    )
    rotated_vector = rotation_matrix.dot(vector[:2])
    return rotated_vector


class PlotAtoms(Atoms):
    def __init__(self, numbers, positions, cell, pbc, style="cpk", **kwargs):
        super().__init__(numbers=numbers, positions=positions, cell=cell, pbc=pbc)
        self._plot_atoms = [PlotAtom.from_ase_atom(atom, **kwargs) for atom in self]
        self.style = style if style in ["cpk", "vdw"] else "cpk"
        if self.style == "cpk":
            self.initialize_bonds(**kwargs)
        else:
            self.bonds = []

    def initialize_bonds(self, **kwargs):
        # print("Initializing bonds...")
        self.bonds = self.identify_bonds(**kwargs)

    def plot(self, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()

        for i, atom in enumerate(self._plot_atoms):
            if self.bonds is not None and len(self.bonds) > 0:
                bonds = self.bonds[i]
                for bond in bonds:
                    bond_patch = bond.bond
                    if isinstance(bond_patch, list):
                        for patch in bond_patch:
                            ax.add_artist(patch)
                    else:
                        ax.add_artist(bond_patch)

            ax.add_artist(atom.circle)
            if atom.shine is not None:
                ax.add_artist(atom.shine)

        ax.set_aspect("equal")

        if "outer" in kwargs.keys():
            outer = kwargs["outer"]
        else:
            outer = 1
        self.plot_cell(outer=outer, ax=ax)
        return ax

    def plot_cell(self, outer=1, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
	
        cell = kwargs.get("cell", None)

        if cell is None:
            a1 = self.cell.array[0]
            a2 = self.cell.array[1]
        else:
            a1 = cell[0]
            a2 = cell[1]

        corners = np.array(
            [
                [0, 0],
                a1[:2],
                (a1 + a2)[:2],
                a2[:2],
                [0, 0],
            ]
        )

        ax.plot(corners[:, 0], corners[:, 1], color="black", linestyle="--")

        xmax = max(a1[0], a2[0]) + outer
        ymax = max(a1[1], a2[1]) + outer
        ax.set_xlim(0 - outer, xmax)
        ax.set_ylim(0 - outer, ymax)
        return ax

    def identify_bonds(self, **kwargs) -> list[PlotBond]:
        er_bonds = electroneg_ratio(AseAtomsAdaptor.get_structure(self), metal_nonmetal_bonus=0.62, min_bond_distance=0.8, max_distance_ratio=1.5, strength_threshold=0.53)
        # Default bond plotting parameters
        bond_radius = kwargs.get("bondradius", 0.1)
        bcolor = kwargs.get("bcolor", "black")
        bcolor_adj = kwargs.get("bcolor_adj", 0)
        becolor_adj = kwargs.get("becolor_adj", -1)
        half = kwargs.get("half", True)
        edgewidth = kwargs.get("edgewidth", 0.01)
        becolor_adj = kwargs.get("becolor_adj", -1)
        alpha = kwargs.get("alpha", 1.0)

        # print(
        #     "Bond radius:",
        #     bond_radius,
        #     "bcolor:",
        #     bcolor,
        #     "bcolor_adj:",
        #     bcolor_adj,
        #     "becolor_adj:",
        #     becolor_adj,
        #     "half:",
        #     half,
        #     "edgewidth:",
        #     edgewidth,
        #     "alpha:",
        #     alpha,
        # )

        new_bonds = []
        for site in er_bonds:
            site_bonds = []
            for bond in site:
                atom_1 = self._plot_atoms[bond["index_1"]]
                atom_2 = self._plot_atoms[bond["index_2"]]
                site_bonds.append(
                    PlotBond(
                        atom1=atom_1,
                        atom2=atom_2,
                        bondradius=bond_radius,
                        bcolor=bcolor,
                        bcolor_adj=bcolor_adj,
                        becolor_adj=becolor_adj,
                        half=half,
                        edgewidth=edgewidth,
                        alpha=alpha,
                    )
                )
            new_bonds.append([b for b in site_bonds])
        return new_bonds

    @classmethod
    def resort_atoms_by_z(cls, positions):
        """
        Resort atoms by their z-coordinate for proper layering in 2D plots.
        Parameters:
            positions (list of tuple): List of atomic positions.
        Returns:
            list of int: Indices that would sort the atoms by z-coordinate.
        """
        if positions.ndim != 2 and positions.shape[1] != 3:
            raise ValueError("Positions should be a 2D array with shape (N, 3).")
        sorted_indices = np.argsort([position[2] for position in positions])
        return sorted_indices

    @classmethod
    def from_ase_atoms(cls, ase_atoms: Atoms, linewidth=1, **kwargs) -> PlotAtoms:

        resort_indices = cls.resort_atoms_by_z(ase_atoms.positions)

        return cls(
            numbers=[ase_atoms.numbers[i] for i in resort_indices],
            positions=[ase_atoms.positions[i] for i in resort_indices],
            cell=ase_atoms.cell,
            pbc=ase_atoms.pbc,
            linewidth=linewidth,
            **kwargs,
        )

    @classmethod
    def from_pymatgen_structure(cls, structure, linewidth=1, **kwargs) -> PlotAtoms:

        resort_indices = cls.resort_atoms_by_z(structure.cart_coords)

        return cls(
            numbers=[structure.atomic_numbers[i] for i in resort_indices],
            positions=[structure.cart_coords[i] for i in resort_indices],
            cell=structure.lattice.matrix,
            pbc=structure.pbc,
            linewidth=linewidth,
            **kwargs,
        )


class PlotAtom(Atom):
    def __init__(
        self,
        symbol="X",
        position=(0, 0, 0),
        color=None,
        color_adj=0,
        color_scheme="jmol",
        scale=0.5,
        linewidth=0.5,
        alpha=1.0,
        shine_shape="ellipse",
        shine_scale=0.1,
        shine_shift=0.3,
        shine_alpha=0.3,
        shine_angle=45,
        **kwargs,
    ):
        super().__init__(symbol=symbol, position=position)

        self.color_adj = color_adj if (color_adj > -1) and (color_adj < 1) else 0.0

        if color is None and color_scheme is not None:
            if color_scheme == "jmol":
                self.color = adjust_color(jmol_colors[self.number], self.color_adj)
            elif color_scheme == "cpk":
                self.color = adjust_color(cpk_colors[self.number], self.color_adj)
            else:
                raise ValueError("Unsupported color scheme.")
        elif color is not None:
            self.color = adjust_color(color, self.color_adj)
        else:
            self.color = adjust_color("gray", self.color_adj)

        self.scale = scale
        self.linewidth = linewidth
        self.alpha = alpha

        self.radius = covalent_radii[self.number]
        self._circle = {
            "scale": scale,
            "linewidth": linewidth,
            "alpha": alpha,
            "color": self.color,
            "color_adj": self.color_adj,
        }

        self._shine = {
            "shape": shine_shape,
            "scale": shine_scale,
            "shift": shine_shift,
            "alpha": shine_alpha,
            "angle": shine_angle,
        }

    @property
    def circle(self):
        radius = self._circle["scale"] * self.radius
        linewidth = self._circle["linewidth"]
        alpha = self._circle["alpha"]
        color = self._circle["color"]
        color_adj = self._circle["color_adj"]

        return Circle(
            self.position[:2],
            radius,
            fc=adjust_color(color, color_adj),
            ec=adjust_color(color, color_adj - 0.3),
            linewidth=linewidth,
            alpha=alpha,
        )

    @circle.setter
    def circle(self, params: dict):
        if "linewidth" in params:
            self._circle["linewidth"] = params["linewidth"]
        if "alpha" in params:
            self._circle["alpha"] = params["alpha"]
        if "color" in params:
            self._circle["color"] = adjust_color(params["color"], self.color_adj)
        if "color_adj" in params:
            self._circle["color_adj"] = params["color_adj"]
        if "shine_params" in params:
            self.shine = params["shine_params"]

    @property
    def shine(self):
        if self._shine["shape"] is None:
            return None
        else:
            if self._shine["shape"] == "circle":
                shift = self._shine["shift"]
                scale = self._shine["scale"]
                alpha = self._shine["alpha"]
                shine_circle = Circle(
                    self.position[:2] + np.array([-shift, +shift]) * self.radius,
                    self.radius * scale,
                    fc="white",
                    alpha=alpha,
                )
                return shine_circle
            elif self._shine["shape"] == "ellipse":
                shift = self._shine["shift"]
                scale = self._shine["scale"]
                alpha = self._shine["alpha"]
                angle = self._shine["angle"]
                shine_ellipse = Ellipse(
                    self.position[:2] + np.array([-shift, +shift]) * self.radius,
                    width=self.radius * scale,
                    height=self.radius * scale / 2,
                    angle=angle,
                    fc="white",
                    alpha=alpha,
                )
                return shine_ellipse
            else:
                raise ValueError("Unsupported shine shape.")

    @shine.setter
    def shine(
        self, shine_shape="ellipse", shine_scale=0.3, shine_shift=0.1, shine_alpha=0.6
    ):
        self._shine = {
            "shape": shine_shape,
            "scale": shine_scale,
            "shift": shine_shift,
            "alpha": shine_alpha,
        }

    @classmethod
    def from_ase_atom(
        cls, ase_atom: Atom, shine_params: Union[dict, None] = None, **kwargs
    ) -> PlotAtom:
        if shine_params is None:
            shine_params = {}
        return cls(
            symbol=ase_atom.symbol,
            position=ase_atom.position,
            **kwargs,
        )


class PlotBond:
    def __init__(
        self,
        atom1: PlotAtom,
        atom2: PlotAtom,
        bondradius=0.1,
        bcolor="black",
        bcolor_adj=0,
        becolor_adj=-1,
        half=True,
        edgewidth=0.0,
        alpha=1.0,
        **kwargs,
        # direction="high",
    ):
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_length = np.linalg.norm(self.atom2.position - self.atom1.position)
        self.bond_vector = (
            self.atom2.position - self.atom1.position
        ) / self.bond_length

        if not half:
            self.bond_color = [adjust_color(bcolor, bcolor_adj)]
        else:
            self.bond_color = [
                adjust_color(atom1.color, bcolor_adj),
                adjust_color(atom2.color, bcolor_adj),
            ]

        self.bond_edge_color_adj = becolor_adj
        self.bond_radius = bondradius
        self.edge_width = edgewidth
        self.alpha = alpha

    @property
    def bond(self):
        """
        Return a PathPatch representing the bond between atom1 and atom2.
        """
        if len(self.bond_color) == 2:
            midpoint = (self.atom1.position + self.atom2.position) / 2 - (
                self.atom2.radius * self.atom2.scale
                - self.atom1.radius * self.atom1.scale
            ) * self.bond_vector
            patch = []
            for i in range(2):
                if i == 0:
                    prefix = 1
                    start = (
                        self.atom1.position
                        + self.atom1.radius
                        * self.atom1.scale
                        * self.bond_vector  # * 0.5
                    )
                    end = (
                        midpoint
                        # * (
                        #     self.atom2.radius * self.atom2.scale
                        #     - self.atom1.radius * self.atom1.scale
                        # )
                        # * self.bond_vector
                        # * 0.5
                    )
                    bcolor = self.bond_color[i]

                else:
                    prefix = -1
                    start = (
                        self.atom2.position
                        - self.atom2.radius
                        * self.atom2.scale
                        * self.bond_vector  # * 0.5
                    )
                    end = (
                        midpoint
                        # * (
                        #     self.atom2.radius * self.atom2.scale
                        #     - self.atom1.radius * self.atom1.scale
                        # )
                        # * self.bond_vector
                    )
                    bcolor = self.bond_color[i]

                path_data = [
                    (
                        Path.MOVETO,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 90) * prefix,
                    ),
                    (
                        Path.CURVE4,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 120) * prefix,
                    ),
                    (
                        Path.CURVE4,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 150) * prefix,
                    ),
                    (
                        Path.CURVE4,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 180) * prefix,
                    ),
                    (
                        Path.LINETO,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 190) * prefix,
                    ),
                    (
                        Path.CURVE4,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 200) * prefix,
                    ),
                    (
                        Path.CURVE4,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 230) * prefix,
                    ),
                    (
                        Path.CURVE4,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 260) * prefix,
                    ),
                    (
                        Path.LINETO,
                        end[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 260) * prefix,
                    ),
                    (
                        Path.LINETO,
                        end[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 90) * prefix,
                    ),
                    (
                        Path.CLOSEPOLY,
                        start[:2]
                        + self.bond_radius * rot_vector(self.bond_vector, 90) * prefix,
                    ),
                ]
                codes, verts = zip(*path_data)
                path = Path(verts, codes)
                patch.append(
                    PathPatch(
                        path,
                        fc=bcolor,
                        ec=adjust_color(bcolor, -0.5),
                        linewidth=self.edge_width,
                        alpha=self.alpha,
                    )
                )
            return patch

        elif len(self.bond_color) == 1:
            path = Line2D(
                [self.atom1.position[0], self.atom2.position[0]],
                [self.atom1.position[1], self.atom2.position[1]],
                color=self.bond_color[0],
                linewidth=self.bond_radius * 2,
                alpha=self.alpha,
            )
            return path

        else:
            raise ValueError("bond_color should be a list of length 1 or 2.")

