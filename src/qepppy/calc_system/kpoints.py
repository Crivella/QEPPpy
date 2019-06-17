import numpy as np
from .symmetry import symmetries
from .lattice  import lattice
# from ..meta import PropertyCreator

def u(r, q):
	return (2.*r - q - 1.), (2. * q)

def cart_to_cryst(cls, coord):
	return coord.dot(np.linalg.inv(cls.recipr))

def cryst_to_cart(cls, coord):
	return coord.dot(cls.recipr)


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

	def __init__(
		self, *args, 
		kpoint_list_cart=None,
		kpoint_list_cryst=None,
		kpoint_mesh=None, kpoint_shift=None,
		**kwargs
		):
		from itertools import combinations

		cond = [not a is None for a in [kpoint_list_cart, kpoint_list_cryst, kpoint_mesh]]
		cond = combinations(cond, 2)
		if any(all(a) for a in cond):
			raise ValueError("Cannot set k-points in multiple different ways at the same time!!!")
		if not kpoint_list_cart is None:
			self.kpt_cart = kpoint_list_cart
		if not kpoint_list_cryst is None:
			self.kpt_cryst = kpoint_list_cryst

		if not kpoint_mesh is None:
			shift = (0,0,0)
			if not kpoint_shift is None:
				shift = kpoint_shift
			self.generate_monkhorst_pack_grid(kpoint_mesh, shift)

		self.symmetries = symmetries()

		super().__init__(*args, **kwargs)

	@property
	def n_kpt(self):
		return len(self.kpt_cart)

	def generate_kpath(
		self, 
		edges,
		# edges_name=[],
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
		self.edges      = edges
		# self.edges_name = edges_name
		self.mode       = mode
		self.mesh       = self.shift = None

		n_pt  = np.array(edges)[:,3].astype(dtype=int)
		edges = np.array(edges)[:,0:3]

		path  = np.empty((0,3))
		for n,i in enumerate(n_pt[:-1]):
			new  = np.vstack((
				np.linspace(edges[n,0], edges[n+1,0], i),
				np.linspace(edges[n,1], edges[n+1,1], i),
				np.linspace(edges[n,2], edges[n+1,2], i)
				)).T
			new  = new[(n>0):,:]
			path = np.vstack((path, new))

		if   mode == 'crystal':
			self.kpt_cryst = path 
		elif mode == 'cart':
			self.kpt_cart  = path

		return path

	def generate_monkhorst_pack_grid(
		self, 
		shape, shift=(0,0,0), 
		set_self=True
		):
		"""
		Generate a Monkhorst-Pack grid of k-point.
		Params:
		 -shape: tuple of 3 ints > 0
		 -shift: tuple of 3 ints that can be either 0 or 1
		 -set_self: If True set the resulting k_point to self.kpt_cryst
		            If False return the generated kpt_list
		"""
		from itertools import product

		assert all(isinstance(a, int) for a in shape)
		assert all(isinstance(a, int) for a in shift)
		assert all(a == 0 or a == 1   for a in shift)
		self.mesh  = shape
		self.shift = shift
		self.edges = self.mode = None
		# self.edges_name = None


		s1,s2,s3   = shift

		l1,l2,l3 = [
			(
				(n+shift[i])/d for n,d in (
					u(r+1,q) for r in range(q)
					)
				) for i,q in enumerate(shape)
			]

		res    = np.array(list(product(l1,l2,l3)))
		_, res = self.symmetries.reduce(res)

		if not set_self:
			return res
		self.kpt_cryst = res


	def kpt_crop(
		self, 
		center=(0,0,0), radius=np.inf, 
		verbose=True, set_self=True
		):
		"""
		Crop the k-points in a sphere of radius 'radius' with center 'center':
		Params:
		 - center: tuple of 3 floats containing the coordinate of the center
		           of the crop sphere. default = (0,0,0)
		 - radius: Radius of the crop sphere. Defaulr = np.inf
		 - verbose: Print information about the cropping. Default = True
		"""
		self.mesh  = self.shift = self.edges = None

		center = np.array(center).reshape(3)
		norms  = np.linalg.norm(self.kpt_cart - center, axis=1)
		w      = np.where(norms <= radius)[0]

		if verbose:
			print(f"# Cropping k-points around {center} with radius {radius}")
			print(f"# Cropped {len(w)} k-points out of {self.n_kpt}")

		if not set_self:
			return self.kpt_cart[w]
		self.kpt_cart = self.kpt_cart[w]
		self.weight   = self.weight[w]

	def _kpt_plot(self, ax, edges_name=[]):
		self.draw_Wigner_Seitz(ax)
		ax.scatter(*self.kpt_cart.T)
		if not self.edges is None:
			ax.scatter(
				*self.edges[:,:3].dot(self.recipr).T, 
				s=10, color='r'
				)

			for name,edge in zip(edges_name, self.edges):
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









