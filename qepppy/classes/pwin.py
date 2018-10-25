import numpy as np

from qepppy.classes.structure import structure as structure#, bravais_index as bi
from qepppy.classes.qe_in import qe_in

class pw_in( qe_in):
	templ_file = "INPUT_PW.def"
	def __iadd__(self, other):
		if isinstance( other, structure):
			self._add_stc_( other)

		return self

	def _add_stc_( self, stc):
		if isinstance( stc.bravais_n, int):
			self.set( nl="SYSTEM", k="ibrav", v=stc.bravais_n)
		else:
			raise Exception( "Must pass a valid cell structure")

		self.set( nl="SYSTEM", k="celldm(1)", v=stc.lp)
		#self._d["SYSTEM"]["celldm(1)"]['v'] = stc.lp
		if self.get( nl="SYSTEM", k="ibrav") == 0:
			if not isinstance( stc.a, np.ndarray):
				raise Exception( "Basis vector must be set with ibrav = 0")

		self.set( nl="SYSTEM", k="ntyp", v=len( stc.atom_spec))
		self.set( nl="SYSTEM", k="nat", v=len( stc.atoms))
		#self._d["SYSTEM"]["ntyp"]['v'] = len( stc.atom_spec)
		#self._d["SYSTEM"]["nat"]['v'] = len( stc.atoms)

		return


