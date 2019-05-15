from .FFTgrid import FFTgrid
from .parser.binary_io import binary_io as bin_io
# from .._decorators import store_property
from .. import utils

class wavefnc(bin_io, FFTgrid):
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
			{'type':'c16', 'shape':('nspin','igwx',), 'name':'C_kn'},
		], 'nbnd'),
	]
	def __init__(self, src=""):
		super().__init__()
		self.rep = 1
		self.src = src
		if src:
			self.read_binary(self.src)

	@property
	# @store_property
	def direct(self):
		return utils.recipr_base(self.recipr)





