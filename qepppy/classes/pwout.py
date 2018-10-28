from .bands     import bands     as bands
from .structure import structure as structure

class pw_out( bands, structure):
	__name__ = "pwout"
	def __init__( self, **kwargs):
		super().__init__( **kwargs)
		self.validate()

 
















