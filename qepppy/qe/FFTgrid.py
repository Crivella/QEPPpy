import numpy as np
import scipy.fftpack
from .. import utils

class FFTgrid():
	def make_density_grid(self, bnd_list=[1]):
		rho = None
		for bnd in bnd_list:
			if rho is None:
				rho = np.abs(self.make_psi_grid(bnd))**2
			else:
				rho += np.abs(self.make_psi_grid(bnd))**2

		return rho

	def make_psi_grid(self, bnd, shape=None):
		bnd -= 1
		GRID = self._generate_g_grid_(bnd, shape)

		return scipy.fftpack.ifftn(GRID)
	

	def test_norm(self):
		if not self.binary:
			raise Exception("Must first read a wavefunction file.")
		for val in self.val:
			norm = np.linalg.norm(val)
			if np.abs(norm - 1) > 1.E-7:
				raise Exception("Wavefunction not normalized.")

	def _generate_g_grid_(self,bnd,shape=None):
		if shape is None:
			shape = np.max(self.gvect, axis=0) - np.min(self.gvect, axis=0)+1
		GRID = np.zeros(
			shape,
			dtype=np.complex
			)
		GRID[tuple(self.gvect.T)] = self.val[bnd]
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

			s1, s2, _ = rho.shape
			a = np.linspace(X.min(),X.max(), s1 * rep)
			b = np.linspace(Y.min(),Y.max(), s2 * rep)

			y,x = np.meshgrid(b,a)
			z   = np.ones(x.shape) * zs

			rec = np.round(
					np.dot(
					self.direct.T.I,
					[x.flatten(),y.flatten(),z.flatten()]
					),
					decimals=8
				) % 1
			rec = np.array(rec) * np.array(rho.shape).reshape(3,1)
			rec = rec.astype(dtype='int')
			i,j,k = rec

			ax1.contourf(
				x, y, rho[tuple((i,j,k))].reshape(x.shape),
				100,cmap=cmap
				)

			ax1.set_title('Slice {}'.format(zs))
			ax1.set_xlabel('X')
			ax1.set_ylabel('Y')

			plt.show()