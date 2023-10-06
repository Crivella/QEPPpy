import importlib
import json
from importlib import resources
from typing import Any

import numpy as np
from attrs import define, field
from numpy import typing as npt
from scipy.spatial import KDTree

from ..graphics import mpl_graphics as mplg
from ..validators import check_allowed, check_shape, converter_none

periodic_table = {}
try:
    importlib.import_module('qepppy')
except ImportError:
    pass
else:
    periodic_table = json.load(resources.files('qepppy.data').joinpath('periodic_table.json').open(encoding='utf-8'))

def cart_to_cryst(cls: 'AtomsList', coord: np.ndarray) -> np.ndarray:
    direct = cls.direct
    if len(direct) == 0:
        return []
    return coord.dot(np.linalg.inv(direct))

def cryst_to_cart(cls: 'AtomsList', coord: np.ndarray) -> np.ndarray:
    direct = cls.direct
    if len(direct) == 0:
        return []
    return coord.dot(direct)

def undo_unique(cls: 'AtomsList', lst: list) -> list[Any]:
    """Get the full list of properties from the unique one.
    Uses the knowledge of the full and unique list of atoms names.
    """
    res = []
    names     = cls.atoms_typ
    all_names = list(cls.unique_atoms_typ)
    for n in names:
        res.append(lst[all_names.index(n)])

    return res

def get_unique(lst: list) -> list[Any]:
    """Get unique elements of a list."""
    res = []
    for name in lst:
        if not name in res:
            res.append(name)

    return res

def get_unique_atm_prop(cls: 'AtomsList', lst: list) -> list[Any]:
    """Get unique list of atomic properties from a list of properties of the same shape as cls.atoms_typ."""
    res = []

    names = list(cls.atoms_typ)
    unq_names = cls.unique_atoms_typ
    for n in unq_names:
        res.append(lst[names.index(n)])

    return res

def split_atom_list_by_name(atom_coord: np.ndarray, atom_names: list[str]) -> tuple[list, np.ndarray, list[float]]:
    """Split a list of atomic coordinates by atom name."""
    trees  = []
    rad    = []
    names  = []

    atom_names = np.array(atom_names)
    for n in atom_names:
        if n in names:
            continue
        coord = atom_coord[atom_names == n,:]

        names.append(n)
        trees.append(KDTree(coord))
        rad.append(periodic_table[n]['radius'])
    return trees, np.array(names), rad

@define(slots=False)
class AtomsList():
    atoms_coord_cart: npt.ArrayLike = field(
        validator=check_shape((-1,3)),
        converter=converter_none(lambda x: np.array(x, dtype=float).reshape(-1,3)),
        default=None
    )

    atoms_forces: npt.ArrayLike = field(
        validator=check_shape((-1,3)),
        converter=converter_none(lambda x: np.array(x, dtype=float).reshape(-1,3)),
        default=None
    )

    atoms_velocities: npt.ArrayLike = field(
        validator=check_shape((-1,3)),
        converter=converter_none(lambda x: np.array(x, dtype=float).reshape(-1,3)),
        default=None
    )

    atoms_typ: npt.ArrayLike = field(
        validator=check_allowed(list(periodic_table.keys())),
        converter=converter_none(lambda x: np.array(x, dtype=str)),
        default=None
    )

    atoms_mass: npt.ArrayLike = field(
        validator=check_shape((-1,)),
        converter=converter_none(lambda x: np.array(x, dtype=float)),
        default=None
    )

    atoms_pseudo: npt.ArrayLike = field(
        validator=check_shape((-1,)),
        converter=converter_none(lambda x: np.array(x, dtype=str)),
        default=None
    )

    @property
    def atoms_coord_cryst(self):
        if self.atoms_coord_cart is None or getattr(self, 'direct', None) is None:
            return
        return cart_to_cryst(self, self.atoms_coord_cart)
    @atoms_coord_cryst.setter
    def atoms_coord_cryst(self, value):
        if getattr(self, 'direct', None) is None:
            return
        self.atoms_coord_cart = cryst_to_cart(self, value)

    @property
    def unique_atoms_typ(self):
        return get_unique(self.atoms_typ)

    @property
    def unique_atoms_mass(self):
        return get_unique_atm_prop(self, self.atoms_mass)
    @unique_atoms_mass.setter
    def unique_atoms_mass(self, value):
        self.atoms_mass = undo_unique(self, value)

    @property
    def unique_atoms_pseudo(self):
        return get_unique_atm_prop(self, self.atoms_pseudo)
    @unique_atoms_pseudo.setter
    def unique_atoms_pseudo(self, value):
        self.atoms_pseudo = undo_unique(self, value)

    @property
    def n_atoms(self):
        res = None
        if not self.atoms_coord_cart is None:
            res = len(self.atoms_coord_cart)

        if hasattr(self, '_n_atoms'):
            try:
                res = self._n_atoms
            except AttributeError:
                self._n_atoms = res
        return res

    @n_atoms.setter
    def n_atoms(self, value):
        self._n_atoms = value

    @property
    def n_types(self):
        res = None
        if not self.unique_atoms_typ is None:
            res = len(self.unique_atoms_typ)

        if hasattr(self, '_n_types'):
            try:
                res = self._n_types
            except AttributeError:
                self._n_types = res
        return res

    @n_types.setter
    def n_types(self, value):
        self._n_types = value

    @property
    def atoms_group_coord_cart(self):
        """List of atomic coordinate in cartesian basis grouped by atom name
        in a dictionary."""
        return {a:np.array(self.atoms_coord_cart[np.array(self.atoms_typ) == a]) for a in self.unique_atoms_typ}

    @property
    def atoms_group_coord_cryst(self):
        """List of atomic coordinate in crystal basis grouped by atom name
        in a dictionary."""
        return {a:np.array(self.atoms_coord_cryst[np.array(self.atoms_typ) == a]) for a in self.unique_atoms_typ}

    def draw_atoms(self, ax, atom_coord=None, atom_names=None, **kwargs):
        """
        Draw atoms onto a matplotlib axis object.
        Params:
         - ax:         Matplotlib axis object
         - atom_coord: List of atomic coordinates (not necessarily the original
                       one in order to use a supercell).
         - atom_namse: List of atom names (same shape as atom_coord). Used to
                       plot the proper color and atomic radius.
        """
        if atom_coord is None:
            atom_coord = self.atoms_coord_cart
        if atom_names is None:
            atom_names = self.atoms_typ

        trees, names, rad = split_atom_list_by_name(atom_coord, atom_names)

        for tree,n,r in zip(trees,names,rad):
            X,Y,Z = tree.data.T
            color = periodic_table[n].get('color', 'k')
            mplg.draw_atom(ax, X,Y,Z, color=color, name=n, radius=r, **kwargs)


    def draw_bonds(self, ax, atom_coord=None, atom_names=None, **kwargs):
        """
        Draw atomic bonds onto a matplotlib axis object.
        Params:
         - ax:         Matplotlib axis object
         - atom_coord: List of atomic coordinates (not necessarily the original
                       one in order to use a supercell).
         - atom_namse: List of atom names (same shape as atom_coord). Used to
                       plot the proper color and atomic radius.
        """
        from itertools import combinations_with_replacement as cwr

        if atom_coord is None:
            atom_coord = self.atoms_coord_cart
        if atom_names is None:
            atom_names = self.atoms_typ

        trees, names, rad = split_atom_list_by_name(atom_coord, atom_names)

        for rad_t, (tree1,tree2), (name1,name2) in zip(cwr(rad,2), cwr(trees,2), cwr(names,2)):
            bonds = tree1.query_ball_tree(tree2, sum(rad_t))
            c1 = periodic_table[name1]['color']
            c2 = periodic_table[name2]['color']
            for i1,b in enumerate(bonds):
                for i2 in b:
                    if i1 == i2 and tree1 == tree2:
                        continue
                    start = tree1.data[i1]
                    end   = tree2.data[i2]
                    mplg.draw_bond(ax, start, end, c1, c2, **kwargs)
