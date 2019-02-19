import numpy as np
from .FFTgrid import FFTgrid
from .parser.binary_io import binary_io as bin_io
from .._decorators import store_property

class charge_density(bin_io, FFTgrid):
	binary_format =[
		[
			{'type':'i4', 'shape':(1,), 'name':'????'},
			{'type':'i4', 'shape':(1,), 'name':'igwx'},
			{'type':'i4', 'shape':(1,), 'name':'ispin'},
		],
		[
			{'type':'f8', 'shape':(3,3), 'name':'recipr'},
		],
		[
			{'type':'i4', 'shape':('igwx',3,), 'name':'gvect'},
		],
		([
			{'type':'c16', 'shape':('igwx',), 'name':'val'},
		], 'ispin'),
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

		return np.mat([a1,a2,a3]) * 2 * np.pi