from .structure import structure
from .bands     import bands

class system(structure, bands):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
		self.symmetries = self.get_symmetries()

	def translate_into_cell(self):
		self.atoms_coord_cryst, _ = self.translate_atoms_into_cell(self.atoms_coord_cryst)

	def translate_into_fbz(self):
		self.kpt_cryst, _ = self.translate_atoms_into_cell(self.kpt_cryst)

	def get_symmetries(self):
		pass