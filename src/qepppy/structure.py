import numpy as np
from .atoms_list  import atoms_list   as atm
from ._lattice    import _lattice     as latt
from . import utils

def cart_to_cryst(cls, coord):
	return coord.dot(np.linalg.inv(cls.direct))

def cryst_to_cart(cls, coord):
	return coord.dot(cls.direct)

class structure(atm, latt):
	atoms_coord_cart={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,3),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'post_set_name':'_atoms_coord_cryst',
		'post_set_func':cart_to_cryst,
		'doc':"""List of atomic coordinate in CARTESIAN basis."""
		}

	atoms_coord_cryst={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,3),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'post_set_name':'_atoms_coord_cart',
		'post_set_func':cryst_to_cart,
		'doc':"""List of atomic coordinate in CRYSTAL basis."""
		}

	def make_supercell(self,repX,repY,repZ):
		"""
		Build a supercell using the repJ repetitions of the primitive cell along
		the J-th direction.

		Params:
		 - repX: repetitions along X
		 - repY: repetitions along Y
		 - repZ: repetitions along Z

		Return:
		 - res: np.array of shape (n_atoms*repX*repY*repZ) containing the 
		        positions of the atoms in the supercell.
		 - typ: list with len = (n_atoms*repX*repY*repZ) containing the 
		        name/type of all the atoms. (1 to 1 correspondence with 'res')
		"""
		from functools import reduce
		assert(isinstance(repX,(int,range,list,tuple)))
		assert(isinstance(repY,(int,range,list,tuple)))
		assert(isinstance(repZ,(int,range,list,tuple)))

		rep    = [range(a) if isinstance(a,int) else a for a in [repX,repY,repZ]]
		R_vec  = utils.generate_repetition_grid(*rep, self.direct)
 
		n_cell = [max(a)-min(a)+1 for a in rep]
		typ    = self.atoms_typ * (reduce(lambda x,y: x*y, n_cell))
 
		res    = np.empty(shape=(0,3))
		coord  = self.atoms_coord_cart
		for vec in R_vec:
			res = np.vstack((res, coord+vec))

		return res, typ

	def nearest_neighbour(self, max_shell=5):
		"""
		Find the nearest neighbour up to a max_shell.

		Params:
		 - max_shell: last shell of nearest neighour returned

		Return:
		 - dist: Array of shape (max_shell,) containing the radius of every shell
		         of nearest neighours.
		 - counts: Array of shape (max_shell,) containing the number of atoms for
		           every shell of nearest neighours.
		"""
		assert(isinstance(max_shell, int))

		m        = max_shell
		l        = range(-m, m+1)
		coord, _ = self.make_supercell(l,l,l)
		norm     = np.linalg.norm(coord, axis=1)

		dist, counts = np.unique(norm, return_counts=True)

		return dist[1:max_shell], counts[1:max_shell]

	def _plot(
		self, ax,
		repX=1, repY=1, repZ=1, 
		cell=False, 
		bonds=True,
		recipr=False,
		graph_lvl=1,):

		if recipr:
			ax.set_xlabel(r"$x (Bohr^{-1})$")
			ax.set_ylabel(r"$y (Bohr^{-1})$")
			ax.set_zlabel(r"$z (Bohr^{-1})$")
			self.draw_Wigner_Seitz(ax)
			# plt.show()
			return
		else:
			ax.set_xlabel("x (Bohr)")
			ax.set_ylabel("y (Bohr)")
			ax.set_zlabel("z (Bohr)")

		L, typ = self.make_supercell(repX,repY,repZ)

		self.draw_atoms(ax, L, typ, graph_lvl=graph_lvl)
		if cell:
			self.draw_direct_cell(ax)
		if bonds:
			self.draw_bonds(ax, L, typ, graph_lvl=graph_lvl)
		ax.legend()

	def plot(
		self, 
		*args, **kwargs,
		# repX=1, repY=1, repZ=1, 
		# cell=False, 
		# bonds=True,
		# recipr=False,
		# graph_lvl=1,
		):
		"""
		Plot the crystal cell structure.
		Args:
		 - reprX/Y/Z=1/1/1[or any positive integer]:
		        repetitions of the cell along X/Y/Z (basis vector not Cartesian!!!).
		        NOTE: They are 3 separate arguments repX, repY, repZ
		 - cell=False[or True]: plot the contour of the cell.
		 - bonds=True[or False]: plot the chemical bonds between atoms.
		 - recip=False[or True]: plot the Brilloiun Zone instead of the real cell.
		        If True all other flags are ignored.
		 - graph_lvl=1[or 0/2/3]:
		   - 0: Basic plot with circle dots as atoms and black lines as bonds.
		   - 1: Colored line as bonds (color of the nearest atom).
		   - 2: Use 3d spheres for atoms and bonds as in 1.
		   - 3: Use 3d spheres for atoms and cylinders for bonds.
		"""
		# from . import cell_graphic as cg
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		self._plot(ax, *args, **kwargs)

		plt.show()

