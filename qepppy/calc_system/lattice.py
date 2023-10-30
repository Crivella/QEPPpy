from typing import Iterable

import numpy as np
from attrs import define, field

from .. import types, utils
from ..graphics import mpl_graphics as mplg
from ..validators import check_shape, converter_none


@define(slots=False,)
class Lattice():
    """Class describing a bravais lattice."""

    direct: Iterable = field(
        validator=check_shape((3,3)),
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

    def draw_direct_cell(
            self,
            ax:
            mplg.Axes,
            center: types.LArray3 = np.zeros(3),
            ):
        """Draw a cell centered on 'center'.

        Args:
            ax (mplg.Axes): matplotlib 3D axis object
            center (types.LArray3, optional): center for plotting the cell. Defaults to np.zeros(3).
        """
        mplg.draw_cell(ax, self.direct, center=center)


    def draw_Wigner_Seitz(
            self,
            ax: mplg.Axes,
            draw_corners: bool = True
            ):
        """Draw the Wigner Seitz cell for the lattice.

        Args:
            ax (mplg.Axes): matplotlib 3D axis object
            draw_corners (bool, optional): Whether to draw the corners as dots. Defaults to True.
        """
        mplg.draw_Wigner_Seitz(ax, self.recipr, draw_corners=draw_corners)

    @staticmethod
    def _transalte_points(
            base: types.LArray3_3,
            coord: types.LArrayN_3,
            mode: str = 'cryst',
            num: int = 5,
            ) -> tuple[types.LArrayN_3, types.LArrayN_3]:
        """Translate a list of points into the Voronoi cell.

        Args:
            base (types.LArray3_3): The basis vectors of the cell as rows of a matrix.
            coord (types.LArrayN_3): The coordinates of the points to translate.
            mode (str, optional): Whether input/output coordinates are in 'cartesian' or 'crystal'. Defaults to 'cryst'.
            num (int, optional): The maximum value of i,j,k where the lattice vector G is G = iGx + jGy + kGz
               Defaults to 5.

        Raises:
            ValueError: If mode is not 'cryst' or 'cart'.

        Returns:
            tuple[types.LArrayN_3, types.LArrayN_3]: The translated coordinates and the translation vectors.
        """
        from itertools import product

        from scipy.spatial import KDTree

        base_inv = np.linalg.inv(base)

        # Operate in cartesian coordinates
        if mode == 'cryst':
            coord_cryst = coord
            coord_cart  = coord.dot(base)
        elif mode == 'cart':
            coord_cart  = coord
            coord_cryst = coord.dot(base_inv)
        else:
            raise ValueError("Mode must be either 'cart' or 'cryst'.")

        # Generate a repetition grid with all possible translation vetors
        l       = range(-num,num+1)
        L_cryst = np.array(list(product(l, repeat=3)))
        L_cart  = L_cryst.dot(base)
        L_tree  = KDTree(L_cart)

        # Find the closest translation vector for each point
        G_l   = []
        for c in coord_cart:
            d,i = L_tree.query(c, k=L_cart.shape[0])
            # There can be multiple vectors with the same distance
            w   = i[np.where(np.abs(d - d[0]) < 1E-6)]
            # Choose the one with the smallest norm in crystal coordinates
            i0  = w[np.argmin(np.linalg.norm(L_cryst[w], axis=1))]

            G_l.append(i0)

        if mode == 'cryst':
            G_res = L_cryst[G_l]
            C_res = coord_cryst - G_res
        elif mode == 'cart':
            G_res = L_cart[G_l]
            C_res = coord_cart - G_res

        return C_res, G_res

    def translate_coord_into_PC(
            self,
            coord: types.LArrayN_3,
            mode: str = 'cryst',
            num: int = 5
            ) -> tuple[types.LArrayN_3, types.LArrayN_3]:
        """Translate a list of atoms coordinates into the Primitive Cell.

        Args:
            coord (types.LArrayN_3): np.array of shape (-1,3) containing the coordinates of the atoms to translate.
            mode (str, optional): Whether input/output coordinates are in 'cartesian' or 'crystal'. Defaults to 'cryst'.
            num (int, optional): The maximum value of i,j,k where the lattice vector G is G = iGx + jGy + kGz
               Defaults to 5.

        Returns:
            tuple[types.LArrayN_3, types.LArrayN_3]: The translated coordinates and the translation vectors.
        """
        return self._transalte_points(self.direct, coord, mode=mode, num=num)

    def translate_coord_into_FBZ(
            self,
            coord: types.LArrayN_3,
            mode: str = 'cryst',
            num: int = 5
            ) -> tuple[types.LArrayN_3, types.LArrayN_3]:
        """Translate a list of reciprocal coordinates into the First Brillouin Zone.

        Args:
            coord (types.LArrayN_3): np.array of shape (-1,3) containing the coordinates of the points to translate.
            mode (str, optional): Whether input/output coordinates are in 'cartesian' or 'crystal'. Defaults to 'cryst'.
            num (int, optional): The maximum value of i,j,k where the lattice vector G is G = iGx + jGy + kGz
                Defaults to 5.
        Returns:
            tuple[types.LArrayN_3, types.LArrayN_3]: The translated coordinates and the translation vectors.
        """
        return self._transalte_points(self.recipr, coord, mode=mode, num=num)
