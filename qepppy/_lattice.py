import numpy as np
# from ._decorators import store_property
from .meta import PropertyCreator
from . import utils

class _lattice(metaclass=PropertyCreator):
	# direct={
	# 	'typ':(list, np.ndarray,),
	# 	'subtyp':(float, np.float,),
	# 	'func':np.array,
	# 	'excp_val':'_recipr',
	# 	'excp_func':utils.recipr_base,
	# 	'doc':"""Matrix of direct lattice vectors as rows."""
	# 	}
	# recipr={
	# 	'typ':(list, np.ndarray,),
	# 	'subtyp':(float, np.float,),
	# 	'func':np.array,
	# 	'excp_val':'_direct',
	# 	'excp_func':utils.recipr_base,
	# 	'doc':"""Matrix of reciprocal lattice vectors as rows."""
	# 	}
	@property
	# @store_property
	def direct(self):
		try:
			return np.array(self._direct)
		except:
			return utils.recipr_base(self._recipr)

	@property
	# @store_property
	def recipr(self):
		try:
			return np.array(self._recipr)
		except:
			return utils.recipr_base(self._direct)

	def draw_direct_cell(self, ax, center=[0,0,0]):
		"""
		Draw a cell centered on 'center'.
		Params:
		 - ax: matplotlib 3D axis object
		 - center: center for plotting the cell
		"""
		V = self.direct
		for n1 in range(3):
			orig = np.array(center)
			v0 = V[n1]
			for n2 in range(4):
				v = np.array(np.vstack((orig, orig + v0)))
				if n2 == n1:
					orig = V[(n2+1)%3] + V[(n2+2)%3]
				else:
					orig = V[n2%3]
				ax.plot(v[:,0], v[:,1], v[:,2], color="black", linewidth=0.5)
			ax.plot(v[:,0], v[:,1], v[:,2], color="black", linewidth=0.5)

	def draw_Wigner_Seitz(self, ax):
		"""
		Draw the Wigner Seitz cell for the lattice.
		Params:
		 - ax: matplotlib 3D axis object
		"""
		from scipy.spatial import Voronoi
		L = utils.generate_repetition_grid([-1,0,1],[-1,0,1],[-1,0,1], self.recipr)

		vor = Voronoi(L)
		P = vor.vertices
		R = vor.ridge_vertices

		rad     = max(np.linalg.norm(self.recipr, axis=1)) * np.sqrt(2)/2
		cond    = np.where(np.linalg.norm(P, axis=1) > rad)[0]
		P[cond] = np.zeros(3)

		for i1, e in enumerate(R):
			for i2, r in enumerate(e):
				if r in cond:
					R[i1][i2] = -1

		X = P[:,0]
		Y = P[:,1]
		Z = P[:,2]
		ax.scatter(X,Y,Z, color='green')

		for vert in R:
			vert.append(vert[0])
			v = np.asarray(vert)
			if np.all(v >= 0):
				ax.plot(P[v, 0], P[v, 1], P[v, 2], color='k')