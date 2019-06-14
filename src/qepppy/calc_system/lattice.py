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
		recipr=None,
		**kwargs
		):
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
