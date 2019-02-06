import numpy as np
from .parser.binary_io import binary_io as bin_io
from ..logger import logger, warning

# @logger()
class wavefnc( bin_io):
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

	def test_norm( self):
		if not self.binary:
			raise Exception( "Must first read a wavefunction file.")
		for i in range( self.nbnd):
			norm = np.linalg.norm( self.val[i])
			if np.abs( norm - 1) > 1.E-7:
				raise Exception( "Wavefunction not normalized.")

	@staticmethod
	def _get_recipr_basis(a1,a2,a3):
		b1 = np.cross(a2,a3)
		b2 = np.cross(a3,a1)
		b3 = np.cross(a1,a2)

		return b1,b2,b3

	@staticmethod
	def _interpolate_2d_grid(grid,n1,n2):
		from scipy.interpolate import griddata
		old_n1,old_n2 = grid.shape
		old_X,old_Y = np.meshgrid(
			np.linspace(0,1,old_n1),
			np.linspace(0,1,old_n2)
			)

		X,Y = np.meshgrid(
			np.linspace(0,1,n1),
			np.linspace(0,1,n2)
			)

		new = griddata(
			(old_X.reshape(old_X.size),old_Y.reshape(old_Y.size)),
			grid.reshape(grid.size),
			(X.reshape(X.size),Y.reshape(Y.size))
			)

		return X,Y,new.reshape(X.shape)

	def _generate_g_grid(self,band):
		C = self.val[band]
		GRID = np.zeros(
			(
				np.max(self.gvect[:,0])*2+1,
				np.max(self.gvect[:,1])*2+1,
				np.max(self.gvect[:,2])*2+1
				),
			dtype=np.complex
			)
		for g,c in zip(self.gvect,C):
			i1,i2,i3 = g
			GRID[i1,i2,i3] = c

		return GRID

	@staticmethod
	def _plot_grid_slice(
		X,Y,grid,
		slice_index=0,
		# interpolate_shape = None,
		cmap='inferno'
		):

		# from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax1 = fig.add_subplot(111,
			# projection='3d'
			)
		# ax2 = fig.add_subplot(122,
		# 	# projection='3d'
		# 	)
		z_slice = grid[slice_index,:,:]


		values = np.real(z_slice*np.conjugate(z_slice))

		# ax1.imshow(values,cmap=cmap)
		ax1.contourf(Y,-X,values,100,cmap=cmap)
		ax1.set_title('Slice {}'.format(slice_index))
		ax1.set_xlabel('X')
		ax1.set_ylabel('Y')
		# if not interpolate_shape is None:
		# 	n1,n2 = interpolate_shape
		# 	X,Y,values = wavefnc._interpolate_2d_grid(values,n1,n2)
		# 	# ax2.imshow(values,cmap=cmap)
		# 	ax2.contourf(Y,-X,values,100,cmap=cmap)
		plt.show()

	def _make_xy_mesh(self,a1,a2,grid):
		n1,n2,n3 = grid.shape
		A,B = np.meshgrid(
			np.arange(grid.shape[2]),
			np.arange(grid.shape[1])
			)
		X = A*a1[0]/n3 + B*a2[0]/n2
		Y = A*a1[1]/n3 + B*a2[1]/n2

		return X,Y

	def charge_density( self,
		bnd_list=[1],
		plot=True,
		# interpolate_shape=None,
		):

		b1,b2,b3 = self.recipr
		a1,a2,a3 = self._get_recipr_basis(b1,b2,b3)
		# print(a1,a2,a3)

		rho = None
		from scipy.fftpack import fftn
		for band in bnd_list:
			band -= 1
			GRID = self._generate_g_grid(band)
			shape = GRID.shape
			if rho is None:
				rho = fftn(GRID,shape)
				rho = rho.real**2 + rho.imag**2
			else:
				app = fftn(GRID,shape)
				app = app.real**2 + app.imag**2
				rho += app

		if plot:
			X,Y = self._make_xy_mesh(a1,a2,rho)

			for i in range(rho.shape[0]):
				self._plot_grid_slice(
					X,Y,rho,
					i,
					# interpolate_shape=interpolate_shape
					)
			return

		return rho

	def spatial_charge_density( self, 
		x=7, y=7, z=7,
		bnd_list=[],
		plot=True, 
		pfile=False, out="charge_density"
		):

		vol = np.dot( self.recipr[0], np.cross( self.recipr[1], self.recipr[2]))
		a1 = np.cross( self.recipr[1], self.recipr[2])
		a2 = np.cross( self.recipr[2], self.recipr[0])
		a3 = np.cross( self.recipr[0], self.recipr[1])
		cell = np.array( (a1,a2,a3)) * 2*np.pi / vol
		#cell = np.array([[1,0,0],[0,1,0],[0,0,1]]) * 10.2 

		div = np.array( (x-1,y-1,z-1))
		delta  = cell / div[:,None]
		delta = delta.transpose()
		recipr = self.recipr.transpose()

		try:
			weight = self.weight
		except Exception as e:
			#logger.warning( "Weight data for kpt #{} missing... settting it to 1".format( self.kpt_num))
			warning.print( "Weight data for kpt #{} missing... setting it to 1".format( self.kpt_num))
			weight = 1

		gvect_cart = np.array( [np.dot(recipr, G) for G in self.gvect])

		grid = np.empty( (x,y,z,3))
		dens = np.empty( (x,y,z))
		coord = np.empty( 3)
		for b in range( self.nbnd):
			if bnd_list:
				if not b+1 in bnd_list:
					continue
			coeff = self.val[b]
			for n1 in range( x):
				for n2 in range( y):
					for n3 in range( z):
						coord = np.dot( delta, (n1,n2,n3))
						grid[n1][n2][n3] = coord
						cc = 0
						for ng, G in enumerate( gvect_cart):
							exp = np.dot( G+self.kpt, coord)
							cc += coeff[ng] * np.exp( 1j * exp)
						dens[n1][n2][n3] = (cc.real**2 + cc.imag**2)

		dens *= weight / vol
		#print( dens)

		"""
		from mayavi import mlab
		src = mlab.pipeline.scalar_field( dens)
		print( src)
		mlab.pipeline.iso_surface( src, contours=[dens.max()/3, ], opacity=0.3)
		mlab.show()
		#"""

		if plot:
			import matplotlib.pyplot as plt
			from mpl_toolkits.mplot3d import Axes3D

			fig = plt.figure()
			for i in range( z):
				ax = fig.add_subplot( 1,z,i+1 , projection='3d')
				ax.set_title( "z slice #{}".format( i))
				ax.view_init( 90, 0)
				ax.margins( tight=True)
				ax.set_zticks([])

				ax.xaxis.pane.fill = False
				ax.yaxis.pane.fill = False
				ax.zaxis.pane.fill = False
				ax.grid( False)
				ax.plot_surface( grid[:,:,i,0], grid[:,:,i,1], dens[:,:,i], cmap='inferno')

			plt.show()
		return dens





