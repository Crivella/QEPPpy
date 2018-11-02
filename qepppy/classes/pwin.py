from .structure import structure as structure
from .qe_in import qe_in

class pw_in( qe_in):
	"""
	Instance used to handle QE pw.x input files
	"""
	templ_file = "INPUT_PW.templ"
	def __iadd__(self, other):
		if isinstance( other, structure):
			self._add_stc_( other)

		return self

	def _add_stc_( self, stc):
		if isinstance( stc.ibrav, int):
			self.set_nl( nl="SYSTEM", k="ibrav", v=stc.ibrav)
		else:
			raise Exception( "Must pass a valid cell structure")

		self.set_nl( nl="SYSTEM", k="celldm(1)", v=stc.lp)
		self.set_nl( nl="SYSTEM", k="ntyp", v=len( stc.atom_spec))
		self.set_nl( nl="SYSTEM", k="nat", v=len( stc.atoms))

		return


