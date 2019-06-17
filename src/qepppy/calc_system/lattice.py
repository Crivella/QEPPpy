import numpy as np
# from ._decorators import store_property
from ..meta     import PropertyCreator
from ..graphics import mpl_graphics as cg
from .. import utils

class lattice(metaclass=PropertyCreator):
	direct={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape':(3,3),
		'conv_func':lambda x: np.array(x, dtype=np.float).reshape(3,3),
		'post_set_name':'_recipr',
		'post_set_func':lambda x,y: utils.recipr_base(y),
		'doc':"""Matrix of direct basis vector (as rows)."""
		}
	recipr={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape':(3,3),
		'conv_func':lambda x: np.array(x, dtype=np.float).reshape(3,3),
		'post_set_name':'_direct',
		'post_set_func':lambda x,y: utils.recipr_base(y),
		'doc':"""Matrix of reciprocal basis vector (as rows)."""
		}

	def __init__(self, *args, 
		direct=None,
		recipr=None,
		**kwargs
		):
		if not direct is None:
			if not recipr is None:
				raise ValueError("Cannot set direct and recipr at the same time.")
			self.direct = direct
		if not recipr is None:
			self.recipr = recipr
		super().__init__(*args, **kwargs)

	def draw_direct_cell(self, ax, center=[0,0,0]):
		"""
		Draw a cell centered on 'center'.
		Params:
		 - ax: matplotlib 3D axis object
		 - center: center for plotting the cell
		"""
		cg.draw_cell(ax, self.direct, center=center)

	def draw_Wigner_Seitz(self, ax):
		"""
		Draw the Wigner Seitz cell for the lattice.
		Params:
		 - ax: matplotlib 3D axis object
		"""
		cg.draw_Wigner_Seitz(ax, self.recipr)

	@staticmethod
	def _transalte_points(
		base, coord, num=5, mode='cryst',
		# verbose=True
		):
		from scipy.spatial import KDTree
		from ..utils import generate_repetition_grid

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
		L_cryst = generate_repetition_grid(l,l,l)
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

	def translate_atoms_into_cell(self, coord, **kwargs):
		return self._transalte_points(self.direct, coord, **kwargs)

	def translate_kpt_into_fbz(self, coord, **kwargs):
		return self._transalte_points(self.recipr, coord, **kwargs)
