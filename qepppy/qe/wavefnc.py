import numpy as np
import scipy.fftpack
from .parser.binary_io import binary_io as bin_io
from .._decorators import store_property

class wavefnc(bin_io):
	binary_format =[
		[
			{'type':'i4', 'shape':(1,), 'name':'kpt_num'},
			{'type':'f8', 'shape':(3,), 'name':'kpt'},
			{'type':'i4', 'shape':(1,), 'name':'ispin'},
			{'type':'i4', 'shape':(1,), 'name':'gamma_only'},
			{'type':'f8', 'shape':(1,), 'name':'scale_factor'},
		],
		[
			{'type':'i4', 'shape':(1,), 'name':'max_index'},
			{'type':'i4', 'shape':(1,), 'name':'igwx'},
			{'type':'i4', 'shape':(1,), 'name':'nspin'},
			{'type':'i4', 'shape':(1,), 'name':'nbnd'},
		],
		[
			{'type':'f8', 'shape':(3,3), 'name':'recipr'},
		],
		[
			{'type':'i4', 'shape':('igwx',3,), 'name':'gvect'},
		],
		([
			{'type':'c16', 'shape':('igwx',), 'name':'val'},
		], 'nbnd'),
	]
	def __init__(self, src=""):
		super().__init__()
		self.rep = 1
		self.src = src
		if src:
			self.read_binary(self.src)

	@property
	@store_property
	def direct(self):
		b1,b2,b3 = self.recipr
		vol = np.linalg.norm(np.dot(b1,np.cross(b2,b3)))
		a1 = np.cross(b2,b3) / vol
		a2 = np.cross(b3,b1) / vol
		a3 = np.cross(b1,b2) / vol

		return np.array([a1,a2,a3]) * 2 * np.pi

	def make_density_grid(self, bnd_list=[1]):
		rho = None
		for band in bnd_list:
			band -= 1
			GRID = self._generate_g_grid_(band)
			shape = GRID.shape
			if rho is None:
				rho = scipy.fftpack.fftn(GRID,shape)
				rho = rho.real**2 + rho.imag**2
			else:
				app = scipy.fftpack.fftn(GRID,shape)
				app = app.real**2 + app.imag**2
				rho += app

		self.dgrid = rho
		return rho
	

	def test_norm(self):
		if not self.binary:
			raise Exception("Must first read a wavefunction file.")
		for val in self.val:
			norm = np.linalg.norm(val)
			if np.abs(norm - 1) > 1.E-7:
				raise Exception("Wavefunction not normalized.")

	def _generate_g_grid_(self,band):
		GRID = np.zeros(
			np.max(self.gvect, axis=0) - np.min(self.gvect, axis=0)+1,
			dtype=np.complex
			)
		GRID[tuple(self.gvect.T)] = self.val[band]
		return GRID

	def _plot_grid_slice(
		self,
		X,Y,Z,grid,
		slice_index=0,
		cmap='inferno'
		):

		# from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
		import scipy.interpolate
		fig = plt.figure()
		# ax1 = fig.add_subplot(131,
		# 	projection='3d'
		# 	)
		# ax2 = fig.add_subplot(132,
		# 	projection='3d'
		# 	)

		# ax3 = fig.add_subplot(133,
		# 	# projection='3d'
		# 	)
		ax1 = fig.add_subplot(111)

		s1, s2, _ = self.dgrid.shape
		a = np.linspace(X.min(),X.max(), s1 * self.rep)
		b = np.linspace(Y.min(),Y.max(), s2 * self.rep)

		y,x = np.meshgrid(b,a)
		z   = np.ones(x.size) * slice_index
		z = scipy.interpolate.griddata(
			(X,Y,Z), grid,
			(x.flatten(),y.flatten(),z),
			'linear',
			fill_value = 0
			).reshape(x.shape)

		ax1.contourf(x,y,z,100,cmap=cmap)
		# ax1.plot_surface(x,y,z,cmap=cmap)
		# ax1.scatter(X,Y,grid,cmap=cmap)
		# ax2.scatter(X,Y,grid,cmap=cmap)
		# ax3.contourf(x,y,z,100,cmap=cmap)
		ax1.set_title('Slice {}'.format(slice_index))
		ax1.set_xlabel('X')
		ax1.set_ylabel('Y')

		# ax1.view_init(90,0)
		plt.show()

	# @property
	def _xyz_mesh_(self):
		n1,n2,n3 = self.dgrid.shape
		rep = self.rep
		a = np.linspace(0, rep, n1*rep + 1)[:-1]
		b = np.linspace(0, rep, n2*rep + 1)[:-1]
		c = np.linspace(0, rep, n3*rep + 1)[:-1]

		# Specific order to obtain the array with shape (n1,n2,n3) as the data grid
		# The 'b,a,c' order is because for a 3d meshgrid the resulting shape is (1,2,3) --> (2,1,3)
		# The 'y,x,z' order is because of how the 3d meshgrid output behaves:
		#    x,y,z=np.meshgrid(1,2,3) 
		#       will cause the x to change value along axis=1
		#					   y to change value along axis=0
		#					   z to change value along axis=2
		# Since the FFT grid has the axis=0,1,2 corresponding to x,y,z i need to do the proper remapping
		y,x,z = np.meshgrid(b,a,c)
		XYZ  = np.dot(
			self.direct.T,
			[x.flatten(),y.flatten(),z.flatten()]
			).reshape(3,*x.shape)

		return np.round(XYZ, decimals=5)

	def plot_density(
		self,
		rep=1,
		bnd_list=[1],
		z_slice=[0]
		):
		arep = 3*rep+2
		rep+=1
		self.rep = arep 
		rho = self.make_density_grid(bnd_list=bnd_list)

		l_slice  = (arep-1)//2
		r_slice  = arep//2
		n1,n2,n3 = rho.shape
		rho = np.pad(rho, 
			(
				(n1 * l_slice, n1 * r_slice),
				(n2 * l_slice, n2 * r_slice),
				(n3 * l_slice, n3 * r_slice)
			), 'wrap')
		X,Y,Z = self._xyz_mesh_()

		xs = np.unique(X)
		ys = np.unique(Y)
		zs = np.unique(Z)
		s1 = (xs.size+1)//arep
		s2 = (ys.size+1)//arep
		s3 = (zs.size+1)//arep
		c_x = ((xs[s1*rep-1] <= X) & (X <= xs[-s1*rep+1]))
		c_y = ((ys[s2*rep-1] <= Y) & (Y <= ys[-s2*rep+1]))
		c_z = ((zs[s3*rep-1] <= Z) & (Z <= zs[-s3*rep+1]))
		w = np.where(c_x & c_y & c_z)
		# print(w)
		# XYZ = XYZ[:,w] # - np.array([xs[s1*rep-1], ys[s2*rep-1], zs[s3*rep-1]]).reshape(3,1,1,1)
		for z in z_slice:
			self._plot_grid_slice(
				X[w]-xs[s1*rep-1], Y[w]-ys[s2*rep-1], Z[w]-zs[s3*rep-1], rho[w],
				z,
				)
		return
		# for i in zs[s3*rep:-s3*rep+2]:
		# 	c_z = ((-1E-4 < Z - i) & (Z - i < 1E-4))
		# 	w = np.where(c_x & c_y & c_z)

		# 	self._plot_grid_slice(
		# 		X[w]-xs[s1*rep-1], Y[w]-ys[s2*rep-1], rho[w],
		# 		i-zs[s3*rep],
		# 		)



