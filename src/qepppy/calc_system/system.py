from .structure import structure
from .bands     import bands
# from .symmetry  import symmetries
# from ..meta     import PropertyCreator


class system(structure, bands):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
		self.symmetries = self.get_symmetries()

	def get_symmetries(self):
		pass