import numpy as np
from ..meta.property_creator import PropertyCreator

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

def inequivalent_iteration(k, sm, thr):
	"""
	Find symmetry inequivalent points starting from a set of point,
	a symm_matrix and a precision threshold.
	Params:
	  - k:   np.ndarray, starting point coordinates
	  - sm:  np.ndarray, symmetry matrix (3x3)
	  - thr: float (def=1e-5), precision threshold for two point being equivalend
	"""
	from scipy.spatial import KDTree

	n_pt  = k.shape[0]
	res_i = np.arange(n_pt, dtype=int)
	res_p = np.empty((0,3))

	tree  = KDTree(k)
	itree = KDTree(k.dot(sm))
	res  = tree.query_ball_tree(itree, thr)

	avoid = []
	n1    = 0
	for n,l in enumerate(res):
		if n in avoid:
			continue
		res_p = np.vstack((res_p, tree.data[n]))
		res_i[n] = n1
		for i in l:
			res_i[i] = n1
			avoid.append(i)
		n1 += 1

	return res_i, res_p # np.array(res_p)

class symmetry(metaclass=PropertyCreator):
	rotation={
		'typ':(np.ndarray,),
		'sub_typ':(np.number,),
		'shape':(3,3),
		'default':np.diag([1.]*3),
		'doc':"""Matrix representation of the symmetry."""
		}

	translation={
		'typ':(np.ndarray,),
		'sub_typ':(np.number,),
		'shape':(3,),
		'default':np.zeros((3,)),
		'doc':"""Translation vector."""
		}

	def __init__(self, *args, **kwargs):
		args = list(args)
		if len(args) > 0:
			self.rotation = args.pop(0)
		if len(args) > 0:
			self.translation = args.pop(0)

		super().__init__(*args, **kwargs)

	def __str__(self):
		return f'R={self.rotation}\nT={self.translation}'

	@property
	def order(self):
		"""Order (N) of the symmetry (A) so that A^N = I."""
		return get_symmetry_order(self.rotation)

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

		symm_matrix = self.rotation

		test     = symm_matrix.copy()
		new_p    = coord.copy() - self.translation
		new_i    = np.arange(coord.shape[0])
		order    = self.order - 1
		order    = max(order, 1)
		for _ in range(order):
			old_i        = new_i
			new_i, new_p = inequivalent_iteration(new_p, test, thr)
			new_i        = new_i[old_i]
			test         = test.dot(test)

		return new_i, new_p

class symmetries(list):
	# def __init__(self, *args, **kwargs):
	# 	super().__init__(*args, **kwargs)
		
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
		          e.g.: [0,1,2,3,2,1,0]
		                means that kpt
		                  4 was remapped onto 2
		                  5 was remapped onto 1
		                  6 was remapped onto 0
		                getting a list of point equivalent by symetry to the original one
		                is done by 'new_p[new_i]'
		 - new_p: The list of points that are not equivalent by symmetry.
		 - symm_remap: List of the indices of the symmetry that was used to remap
		               every k-point. If some indexes are -1, it means that that
		               k-point was not remapped onto any other k-point using a symmetry
		               (even if other points might have been remapped onto it).
		               eg: [-1 -1  0  0  1]
		                   means that kpt
		                     0,1 were not remapped
		                     2,3 were remapped using the symmetry with index 0
		                     4   was  remapped using the symmetry with index 1
		"""
		assert(isinstance(coord, np.ndarray))
		assert(coord.shape[1] == 3)

		res_i = np.arange(coord.shape[0])
		res_p = coord.copy()

		symm_remap = np.ones((coord.shape[0],), dtype=int) * -1

		# print()
		# print(coord)
		# print(coord.shape[0])
		for n,matrix in enumerate(self):
			new_i, new_p = matrix.reduce(res_p, thr=thr)
			res_p = new_p
			res_i = new_i[res_i]
			# print(len(res_i))
			# raise NotImplementedError("NEED FIX. should use res_i nad not new_i !!!!!!!!!!!!!!!")
			np.set_printoptions(linewidth=2000, formatter={'int':lambda x:f'{x:3d}'})
			# print()
			# print('-'*50)
			# print(n, symm_remap)
			# print(matrix)
			uniq, count = np.unique(new_i, return_counts=True)
			for u,c in zip(uniq, count):
				if c == 1:
					continue
				w = np.where(new_i == u)
				# print(w[0][0], end=', ')
				w = w[0][1:]

				# if np.any(symm_remap[w] != -1):
				# 	raise NotImplementedError()
				symm_remap[w] = n

		# print(symm_remap)
		return res_i, res_p, symm_remap

	

