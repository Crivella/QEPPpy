import numpy as np
from ._atoms   import _atoms   as atm
from ._lattice import _lattice as latt
from . import cell_graphic as cg


class _cell(atm, latt):
	def make_supercell(self,repX,repY,repZ):
		rep = [range(a) for a in [repX,repY,repZ]]
		shift = cg.generate_repetition_grid(*rep, self.direct)
		typ = self.atoms_typ * (repX*repY*repZ)

		res = np.empty(shape=(0,3))
		for vec in shift:
			res = np.vstack((res, self.atoms_coord_cart+vec))

		return res, typ

	def plot(
		self, 
		repX=1, repY=1, repZ=1, 
		cell=False, 
		bonds=True,
		recipr=False,
		graph_lvl=1,
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


		if recipr:
			ax.set_xlabel(r"$x (Bohr^{-1})$")
			ax.set_ylabel(r"$y (Bohr^{-1})$")
			ax.set_zlabel(r"$z (Bohr^{-1})$")
			self.draw_Wigner_Seitz(ax)
			plt.show()
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
		plt.show()

