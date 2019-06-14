import numpy as np
import scipy.fftpack
from .. import utils

class FFTgrid():
	def make_density_grid(self, bnd_list=[1]):
		rho = 0
		for bnd in bnd_list:
			rho += np.abs(self.make_psi_grid(bnd))**2

		return rho

	def make_psi_grid(self, bnd, shape=None):
		bnd -= 1
		GRID = self._generate_g_grid_(bnd, shape)

		return scipy.fftpack.ifftn(GRID)
	

	def test_norm(self):
		if not self.binary:
			raise Exception("Must first read a wavefunction file.")
		for val in self.C_kn:
			norm = np.linalg.norm(val)
			if np.abs(norm - 1) > 1.E-7:
				raise Exception("Wavefunction not normalized.")

	def _generate_g_grid_(self,bnd,shape=None):
		if shape is None:
			shape = np.max(self.gvect, axis=0) - np.min(self.gvect, axis=0)+1
			shape = [int(self.nspin)] + list(shape)
			
		GRID = np.zeros(
			shape,
			dtype=np.complex
			)

		GRID[:,self.gvect[:,0],self.gvect[:,1],self.gvect[:,2]] = self.C_kn[bnd]

		return GRID

	def plot_density_zslice(
		self,
		rep=1,
		bnd_list=[1],
		z_slice=[0],
		cmap='inferno'
		):

		rho = self.make_density_grid(bnd_list=bnd_list)
		X,Y,Z = utils.xyz_mesh(
			rho.shape,
			base=self.direct,
			rep=rep,
			)

		import matplotlib.pyplot as plt
		for zs in z_slice:
			fig = plt.figure()
			ax1 = fig.add_subplot(111)

			x,y, index = utils.remap_plane(
				self.recipr/(2*np.pi),
				(X.min(),X.max()),
				(Y.min(),Y.max()),
				(zs,zs),
				np.array(rho.shape[1:]),
				(rep,)*3
				)

			ax1.contourf(
				x, y, rho[(0,*index)].reshape(x.shape) + rho[(1,*index)].reshape(x.shape),
				30, cmap=cmap
				)

			ax1.set_title('Slice {}'.format(zs))
			ax1.set_xlabel('X')
			ax1.set_ylabel('Y')

			plt.show()