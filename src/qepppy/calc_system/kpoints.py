import numpy as np
from .lattice  import lattice
from .._decorators import set_self
# from ..meta import PropertyCreator

def u(r, q):
	return (2.*r - q - 1.), (2. * q)

def cart_to_cryst(cls, coord):
	recipr = cls.recipr
	if recipr == []:
		return
	return coord.dot(np.linalg.inv(recipr))

def cryst_to_cart(cls, coord):
	recipr = cls.recipr
	if recipr == []:
		return
	return coord.dot(recipr)


class kpoints(lattice):
	kpt_cart={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,3),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'post_set_name':'_kpt_cryst',
		'post_set_func':cart_to_cryst,
		'doc':"""List of k-points coordinate in cartesian basis (k_i, i=x,y,z in units of 2pi/a)."""
		}
	kpt_cryst={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,3),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'post_set_name':'_kpt_cart',
		'post_set_func':cryst_to_cart,
		'doc':"""List of k-points coordinate in cartesian basis (k_i, i=x,y,z in units of 2pi/a)."""
		}
	weight={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape':('n_kpt',),
		'conv_func':lambda x: np.array(x, dtype=np.float)/np.array(x, dtype=np.float).sum(),
		'doc':"""List of k-points weights."""
		}
	kpt_mesh={
		'typ':(tuple,),
		'sub_typ':(int, np.integer),
		'shape': (3,),
		'conv_func':lambda x: np.array(x, dtype=np.int),
		'doc':"""Density of the Monkhorst-Pack mesh grid."""
		}
	kpt_shift={
		'typ':(tuple,),
		'sub_typ':(int, np.integer),
		'shape': (3,),
		'default':(0,0,0),
		'allowed':(0,1),
		'conv_func':lambda x: np.array(x, dtype=np.int),
		'doc':"""Shift of the Monkhorst-Pack mesh grid."""
		}
	kpt_edges={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,4,),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'doc':"""k-path where every i-th row is the 3 coordinates + the number of point
		between the i-th and i-th + 1  k-points."""
		}
	mode={
		'typ':(str,),
		'allowed':['cart','cryst', 'crystal', 'cartesian'],
		'doc':"""Set mode for the k-path."""
		}

	def __init__(self, *args, **kwargs):
		if not hasattr(self, 'symmetries'):
			from .symmetry import symmetries
			self.symmetries = symmetries()

		if self.kpt_mesh != ():
			self.generate_monkhorst_pack_grid()
		if self.kpt_edges != []:
			self.generate_kpath()

		super().__init__(*args, **kwargs)


	@property
	def n_kpt(self):
		return len(self.kpt_cart)

	def generate_kpath(
		self, 
		edges=None,
		mode='crystal'
		):
		"""
		Generate a k-point path.
		Params:
		 - edges: list or np.array of shape (N+1, 4), where N is the number of
		          k-points lines.
		          The first 3 columns have to contain the coordinates of the 
		          k-point.
		          The 4th columns has to contain an integer > 0 that indicates
		          the number of k-points between the current and the next edge.
		 - mode: Referse to the basis set for the k-point coordinates.
		         - 'crystal': coordinates given in b1,b2,b3 units
		         - 'cart':    coordinates given in kx,ky,kz units (2pi/a).
		Return:
		 np.array of shape (edges[:-1,-1].sum(), 3), where every row is the 3D
		 coordinate of a k-point.
		"""
		self.mode = mode
		if edges is None:
			edges = self.kpt_edges
		else:
			self.kpt_edges = edges
		self.kpt_mesh  = None
		self.kpt_shift = (0,0,0)

		n_pt  = np.array(edges)[:,3].astype(dtype=int) + 1
		edges = np.array(edges)[:,:3]

		path  = np.empty((0,3))
		for n,i in enumerate(n_pt[:-1]):
			new  = np.vstack((
				np.linspace(edges[n,0], edges[n+1,0], i),
				np.linspace(edges[n,1], edges[n+1,1], i),
				np.linspace(edges[n,2], edges[n+1,2], i)
				)).T
			new  = new[(n>0):,:]
			path = np.vstack((path, new))

		if   'cryst' in mode:
			self.kpt_cryst = path 
		elif 'cart' in mode:
			self.kpt_cart  = path

		return path

	def generate_unfolding_path(self, SC_rec, mode='cryst', return_all=False):
		from .symmetry import symmetries

		PC_rec     = self.recipr
		SC_rec_inv = np.linalg.inv(SC_rec)
		M          = SC_rec_inv.dot(PC_rec).T

		K_G     = self.kpt_cryst.dot(M.T)
		
		max_G_i = int(np.ceil(((np.linalg.det(M) ** (1/3)-1)/2))) + 3
		new, _  = self._transalte_points(SC_rec, K_G, mode='cryst', num=max_G_i)

		if mode == 'cart':
			new = new.dot(SC_rec)

		full_kpts = new

		symm = symmetries()
		symm.append(np.diag([1]*3))
		symm.append(np.diag([-1]*3))

		ind, red_kpts = symm.reduce(new)

		if return_all:
			return ind, full_kpts, red_kpts

		return red_kpts


	@set_self('kpt_cryst,weight')
	def generate_monkhorst_pack_grid(self, mesh=None, shift=None):
		"""
		Generate a Monkhorst-Pack grid of k-point.
		Params:
		 -mesh:  tuple of 3 ints > 0
		 -shift: tuple of 3 ints that can be either 0 or 1
		"""
		from itertools import product

		if mesh is None:
			mesh = self.kpt_mesh
		if shift is None:
			shift = self.kpt_shift

		self.kpt_mesh  = mesh
		self.kpt_shift = shift
		self.kpt_edges = self.mode = None

		s1,s2,s3   = shift

		l1,l2,l3 = [
			(
				(n+shift[i])/d for n,d in 
					(
						u(r+1,q) for r in range(q)
					)
			) for i,q in enumerate(mesh)
			]


		kpts = np.array(list(product(l1,l2,l3)))

		# from ase.spacegroup import Spacegroup
		# print(Spacegroup(227).unique_sites(kpts))

		self.full_kpt_cryst = kpts

		ind, kpts = self.symmetries.reduce(kpts)
		self.irrep_mapping  = ind

		unique, counts = np.unique(ind, return_counts=True)
		weight = counts[np.argsort(unique)]

		return kpts, weight


	@set_self('kpt_cart,weight')
	def kpt_crop(
		self, 
		center=(0,0,0), radius=np.inf, 
		verbose=True
		):
		"""
		Crop the k-points in a sphere of radius 'radius' with center 'center':
		Params:
		 - center: tuple of 3 floats containing the coordinate of the center
		           of the crop sphere. default = (0,0,0)
		 - radius: Radius of the crop sphere. Defaulr = np.inf
		 - verbose: Print information about the cropping. Default = True
		"""
		if radius < 0:
			raise ValueError("Radius must be greather than 0.")
		self.mesh  = self.shift = self.edges = None

		center = np.array(center).reshape(3)
		norms  = np.linalg.norm(self.kpt_cart - center, axis=1)
		w      = np.where(norms <= radius)[0]

		if verbose:
			print(f"# Cropping k-points around {center} with radius {radius}")
			print(f"# Cropped {len(w)} k-points out of {self.n_kpt}")

		return self.kpt_cart[w], self.weight[w]/np.sum(self.weight[w])

	def _kpt_plot(self, ax, edges_name=[]):
		self.draw_Wigner_Seitz(ax)
		ax.scatter(*self.kpt_cart.T)
		if not self.kpt_edges is None:
			ax.scatter(
				*self.kpt_edges[:,:3].dot(self.recipr).T, 
				s=10, color='r'
				)

			for name,edge in zip(edges_name, self.kpt_edges):
				e = edge[:3]
				if self.mode == 'crystal':
					e = e.dot(self.recipr)
				ax.text(*e, name, size=15)


	def kpt_plot(self, **kwargs):
		"""
		Plot the FBZ and the with the selected k-points inside.
		Params:
		 - edges_name: List containing the names of the highsymm points.
		               Used if the kpt are generated using generate_kpath.
		"""
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		self._kpt_plot(ax, **kwargs)

		ax.set_xlabel(r'$k_x$')
		ax.set_ylabel(r'$k_y$')
		ax.set_zlabel(r'$k_z$')
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_zticklabels([])

		plt.show()









