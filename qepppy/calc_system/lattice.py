from typing import Iterable

import numpy as np
from attrs import define, field

from .. import utils
from ..graphics import mpl_graphics as mplg
from ..validators import check_shape, converter_none


@define(slots=False)
class lattice():
    """Class describing a bravais lattice."""

    direct: Iterable = field(
        validator=[
            check_shape((3,3)),
            # mutually_exclusive('recipr'),
            # post_setter('recipr', lambda x: 2*np.pi*np.linalg.inv(x).T),
            ],
        converter=converter_none(lambda x: np.array(x, dtype=float).reshape(3,3)),
        default=None
    )

    @property
    def recipr(self):
        if self.direct is None:
            return
        return utils.recipr_base(self.direct)
    @recipr.setter
    def recipr(self, value):
        self.direct = utils.recipr_base(value)

    def draw_direct_cell(self, ax, center=[0,0,0]):
        """
        Draw a cell centered on 'center'.
        Params:
         - ax: matplotlib 3D axis object
         - center: center for plotting the cell
        """
        mplg.draw_cell(ax, self.direct, center=center)

    def draw_Wigner_Seitz(self, ax, draw_corners=True):
        """
        Draw the Wigner Seitz cell for the lattice.
        Params:
         - ax: matplotlib 3D axis object
         - draw_corners: Draw the cell corners as green dots
        """
        mplg.draw_Wigner_Seitz(ax, self.recipr, draw_corners=draw_corners)

    @staticmethod
    def _transalte_points(
        base, coord,
        mode='cryst',
        num=5,
        ):
        from itertools import product

        from scipy.spatial import KDTree

        base_inv = np.linalg.inv(base)

        if mode == 'cryst':
            coord_cryst = coord
            coord_cart  = coord.dot(base)

        elif mode == 'cart':
            coord_cart  = coord
            coord_cryst = coord.dot(base_inv)
        else:
            raise ValueError("Mode must be either 'cart' or 'cryst'.")

        l       = range(-num,num+1)
        L_cryst = np.array(list(product(l, repeat=3)))
        L_cart  = L_cryst.dot(base)
        L_tree  = KDTree(L_cart)

        G_l   = []
        for c in coord_cart:
            d,i = L_tree.query(c, k=L_cart.shape[0])
            w   = i[np.where(np.abs(d - d[0]) < 1E-6)]
            i0  = w[np.argmin(np.linalg.norm(L_cryst[w], axis=1))]

            G_l.append(i0)

        if mode == 'cryst':
            G_res = L_cryst[G_l]
            C_res = coord_cryst - G_res
        elif mode == 'cart':
            G_res = L_cart[G_l]
            C_res = coord_cart - G_res

        return C_res, G_res

    def translate_coord_into_PC(self, coord, **kwargs):
        """
        Translate a list of atoms coordinates into the Primitive Cell.
        Params:
         - coord: np.array of shape (-1,3) containing the coordinates
                  of the atoms to translate.
         - mode: 'crystal'/'cart' Specify the basis for the coordinates.
         - num:  integer that specify the maximum number of cells on
                 which to test.
        """
        return self._transalte_points(self.direct, coord, **kwargs)

    def translate_coord_into_FBZ(self, coord, **kwargs):
        """
        Translate a list of k-points coordinates into the First BZ.
        Params:
         - coord: np.array of shape (-1,3) containing the coordinates
                  of the k-points to translate.
         - mode: 'cryst'/'cart' Specify the basis for the coordinates.
         - num:  integer that specify the maximum number of cells on
                 which to test.
        """
        return self._transalte_points(self.recipr, coord, **kwargs)
