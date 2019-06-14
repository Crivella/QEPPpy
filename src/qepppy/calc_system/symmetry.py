import numpy as np

from ..meta import PropertyCreator

def get_symmetry_order(symm_matrix):
	from scipy.linalg import eigvals

	egv   = eigvals(symm_matrix)
	check = [1] * len(egv)

	order = 1
	while not np.allclose(egv**order, check):
		order += 1
		if order > 100:
			raise ValueError('Order > 100, stopping iteration')

	return order

def _iteration(k, sm, thr):
	from scipy.spatial import KDTree

	n_pt  = k.shape[0]
	res_i = np.arange(n_pt, dtype=int)
	res_p = []

	tree  = KDTree(k)
	itree = KDTree(k.dot(sm))
	res  = tree.query_ball_tree(itree, thr)
	avoid = []
	n1    = 0
	for n,l in enumerate(res):
		if n in avoid:
			continue
		res_p.append(tree.data[n])
		res_i[n] = n1
		for i in l:
			res_i[i] = n1
			avoid.append(i)
		n1 += 1

	return res_i, np.array(res_p)

class symmetry(metaclass=PropertyCreator):
	matrix={
		'typ':(np.ndarray,),
		'sub_typ':(np.number,),
		'shape':(3,3),
		'default':np.diag([1.]*3),
		'doc':"""Matrix representation of the symmetry."""
		}

	def __init__(self, matrix=None):
		if not matrix is None:
			self.matrix = np.array(matrix)

	@property
	def order(self):
		"""Order (N) of the symmetry (A) so that A^N = I."""
		return get_symmetry_order(self.matrix)

	def reduce(self, coord, thr=1E-5):
		"""
		Reduce a list of points using the symmetry.
		Params:
		 - coord: np.array of shape(-1,3) containing the list of point to reduce

		Return: (new_i,new_p)
		 - new_i: List of indices that map the not-reduced set of points to the
		          reduced one.
		          new_p[new_i] will give a list of points that is equivalent to
		          the not-reduced one by symmetry.
		 - new_p: The list of points that are not equivalent by symmetry.
		"""
		assert(isinstance(coord, np.ndarray))
		assert(coord.shape[1] == 3)

		symm_matrix = self.matrix

		test     = symm_matrix.copy()
		new_p    = coord.copy()
		new_i    = np.arange(coord.shape[0])
		order    = self.order - 1
		order    = max(order, 1)
		for _ in range(order):
			old_i        = new_i
			new_i, new_p = _iteration(new_p, test, thr)
			new_i        = new_i[old_i]
			test         = test.dot(test)

		return new_i, new_p

class symmetries(list):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
	def append(self, value):
		if isinstance(value, symmetry):
			new = value
		else:
			new = symmetry(value)

		super().append(new)


	def reduce(self, coord, thr=1E-5):
		"""
		Reduce a list of points using all the symmetries.
		Params:
		 - coord: np.array of shape(-1,3) containing the list of point to reduce

		Return: (new_i,new_p)
		 - new_i: List of indices that map the not-reduced set of points to the
		          reduced one.
		          new_p[new_i] will give a list of points that is equivalent to
		          the not-reduced one by symmetry.
		 - new_p: The list of points that are not equivalent by symmetry.
		"""
		assert(isinstance(coord, np.ndarray))
		assert(coord.shape[1] == 3)

		res_i = np.arange(coord.shape[0])
		res_p = coord.copy()
		for matrix in self:
			new_i, new_p = matrix.reduce(res_p, thr=thr)
			res_p = new_p
			res_i = new_i[res_i]

		return res_i, res_p

	

