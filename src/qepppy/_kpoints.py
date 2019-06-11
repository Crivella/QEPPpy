import numpy as np
from .meta import PropertyCreator
from .utils import _cart_to_cryst_, _cryst_to_cart_

def u(r, q):
	return (2.*r - q - 1.), (2. * q)

class _kpoints(metaclass=PropertyCreator):
	kpt_cart={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,3),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'post_set_name':'_kpt_cryst',
		'post_set_func':_cart_to_cryst_,
		'doc':"""List of k-points coordinate in cartesian basis (k_i, i=x,y,z in units of 2pi/a)."""
		}
	kpt_cryst={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,3),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'post_set_name':'_kpt_cart',
		'post_set_func':_cryst_to_cart_,
		'doc':"""List of k-points coordinate in cartesian basis (k_i, i=x,y,z in units of 2pi/a)."""
		}
	recipr={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		# 'size':9,
		'shape':(3,3),
		'conv_func':lambda x: np.array(x, dtype=np.float).reshape(3,3),
		'doc':"""Matrix of reciprocal basis vector (as rows)."""
		}
	weight={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape':('n_kpt',),
		'conv_func':lambda x: np.array(x, dtype=np.float)/np.array(x, dtype=np.float).sum(),
		'doc':"""List of k-points weights."""
		}

	def __init__(self):
		pass

	@property
	def n_kpt(self):
		return len(self.kpt_cart)

	def generate_kpath(self, edges, mode='crystal'):
		"""
		Generate a k-point path.
		Params:
		 - edges: list or np.array of shape (N+1, 4), where N is the number of
		          k-points lines.
		          The first 3 columns have to contain the coordinates of the 
		          k-point.
		          The 4th columns has to contain and integer > 0 that indicates
		          the number of k-points between the current and the next edge.
		 - mode: Referse to the basis set for the k-point coordinates.
		         - 'crystal': coordinates given in b1,b2,b3 units
		         - 'cart':    coordinates given in kx,ky,kz units (2pi/a).
		Return:
		 np.array of shape (edges[:-1,-1].sum(), 3), where every row is the 3D
		 coordinate of a k-point.
		"""
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

	def generate_monkhorst_pack_grid(self, shape, shift=(0,0,0), set_self=True):
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
		s1,s2,s3 = shift

		l1,l2,l3 = [
			(
				(n+shift[i])/d for n,d in (
					u(r+1,q) for r in range(q)
					)
				) for i,q in enumerate(shape)
			]
		res = np.array(list(product(l1,l2,l3)))
		if not set_self:
			return res
		self.kpt_cryst = res


	def reduce_by_space_symmetry(self):
		pass

	def reduce_by_time_reversal(self):
		pass

	def crop(self, center=(0,0,0), radius=np.inf, verbose=True):
		center = np.array(center).reshape(3)
		norms  = np.linalg.norm(self.kpt_cart - center, axis=1)
		w      = np.where(norms <= radius)[0]

		if verbose:
			print(f"# Cropping k-points around {center} with radius {radius}")
			print(f"# Cropped {len(w)} k-points out of {self.n_kpt}")
		self.kpt_cart = self.kpt_cart[w]
		self.weight   = self.weight[w]









