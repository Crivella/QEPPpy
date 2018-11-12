import numpy as np
from .qe_binary import qe_binary_reader as qbr
from .logger import *

@logger( msg_lvl='error')
class qe_wfc( qbr):
	wfc_format =[
		[
			{'t':4, 's':(1,), 'n':'kpt_num'},
			{'t':8, 's':(3,), 'n':'kpt'},
			{'t':4, 's':(1,), 'n':'ispin'},
			{'t':4, 's':(1,), 'n':'gamma_only'},
			{'t':8, 's':(1,), 'n':'scale_factor'},
		],
		[
			{'t':4, 's':(1,), 'n':'max_index'},
			{'t':4, 's':(1,), 'n':'igwx'},
			{'t':4, 's':(1,), 'n':'nspin'},
			{'t':4, 's':(1,), 'n':'nbnd'},
		],
		[
			{'t':8, 's':(3,3), 'n':'recipr'},
		],
		[
			{'t':4, 's':('igwx',3,), 'n':'gvect'},
		],
		([
			{'t':16, 's':('igwx',), 'n':'val'},
		], 'nbnd'),
	]

	def test_norm( self):
		if not self.binary:
			raise Exception( "Must first read a wavefunction file.")
		for i in range( self.nbnd):
			norm = np.linalg.norm( self.val[i])
			if np.abs( norm - 1) > 1.E-7:
				raise Exception( "Wavefunction not normalized.")

	def spatial_charge_density( self, 
		x=7, y=7, z=7,
		bnd_list=[],
		plot=True, 
		pfile=False, out="charge_density"):
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
			raise warning( "Weight data for kpt #{} missing... setting it to 1".format( self.kpt_num))
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





