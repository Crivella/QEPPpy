from .structure import structure
from .parser.input_files import input_files as inp_f

class pw_in(inp_f, structure):
	"""
	Instance used to handle QE pw.x input files
	kwargs:
	 - inp = Name of the file to parse
	"""
	templ_file = "INPUT_PW.templ"

	@property
	def _atoms(self):
		name, x, y, z = self.find("X", "x", "y", "z", up="ATOMIC_POSITIONS")
		coord = [(a,b,c) for a,b,c in zip(x, y, z)]
		return [{'name':n, 'coord':c} for n,c in zip(name, coord)]

	@property
	def _atom_spec(self):
		name, mass, pfile = self.find("X", "Mass_X", "PseudoPot_X", up="ATOMIC_SPECIES")
		return [{'name':n, 'mass':m, 'pseudo_file':p} for n,m,p in zip(name, mass, pfile)]

	@property
	def _cell(self):
		a1, a2, a3 = self.find("v1", "v2", "v3")
		return [{'a1':a1, 'a2':a2, 'a3':a3}]

	@property
	def celldm(self):
		return self.find("celldm")

	@property
	def alat(self):
		return self.find("celldm(1)")

	@property
	def ibrav(self):
		return self.find("ibrav")

	@property
	def _atom_p(self):
		return self.find("ATOMIC_POSITIONS")

	@property
	def _cell_p(self):
		return self.find("CELL_PARAMETERS")

	@property
	def n_atoms(self):
		return self.find("nat")

	@property
	def n_types(self):
		return self.find("ntyp")

	def __str__(self):
		msg = inp_f.__str__(self)
		if not 'ATOMIC_POSITIONS' in self.card:
			msg += structure.__str__(self)
		return msg



