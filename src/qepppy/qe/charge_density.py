from .FFTgrid import FFTgrid
from ..parsers import fortran_binary_io as bin_io
from ..utils import recipr_base

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
			{'type':'c16', 'shape':('igwx',), 'name':'C_kn'},
		], 'ispin'),
	]
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		# self.rep = 1

	@property
	def direct(self):
		return recipr_base(self.recipr)
