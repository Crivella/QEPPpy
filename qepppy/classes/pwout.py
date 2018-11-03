from .bands     import bands     as bands
from .structure import structure as structure

class pw_out( bands, structure):
	"""
	Instance used to handle QE outputs (by parsing the "data-file*.xml" file")
	fname: name of the "data-file*.xml" to parse
	"""
	__name__ = "pw_out"
	def __init__( self, fname="", **kwargs):
		super().__init__( fname=fname, **kwargs)
		self.validate()

 
















