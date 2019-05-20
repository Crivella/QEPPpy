import numpy as np
from .meta import PropertyCreator

class _kpoints(metaclass=PropertyCreator):
	n_kpt={
		'typ':(int,),
		'default':0,
		'doc':"""Number of k-points."""
		}
	kpt_cart={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'size':'n_kpt * 3',
		'usize':True,
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'doc':"""List of k-points coordinate in cartesian basis (k_i, i=x,y,z in units of 2pi/a)."""
		}
	weight={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'size':'n_kpt',
		'usize':True,
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'doc':"""List of k-points weights."""
		}

	def generate_kpath(self, edges, type='crystal'):
		"""
		Generate a k-point path
		"""
		pass

