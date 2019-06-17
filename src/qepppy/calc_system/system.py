from .structure import structure
from .bands     import bands
from .._decorators import set_self

class system(structure, bands):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
		self.symmetries = self.get_symmetries()

	@set_self('atoms_coord_cryst')
	def translate_into_PC(self):
		"""
		Translate all the atoms into the Primitive Cell.
		"""
		res, _ = self.translate_coord_into_PC(self.atoms_coord_cryst)
		return res

	@set_self('kpt_cryst')
	def translate_into_FBZ(self):
		"""
		Translate all the k-points into the First BZ.
		"""
		res, _ = self.translate_coord_into_FBZ(self.kpt_cryst)
		return res

	@set_self('symmetries')
	def get_symmetries(self):
		pass
