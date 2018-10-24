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
			self._d["SYSTEM"]["ibrav"]['v'] = stc.bravais_n
		else:
			raise Exception( "Must pass a valid cell structure")

		self._d["SYSTEM"]["celldm(1)"]['v'] = stc.lp
		if self._d["SYSTEM"]["ibrav"]['v'] == 0:
			if not isinstance( stc.a, np.ndarray):
				raise Exception( "Basis vector must be set with ibrav = 0")

		self._d["SYSTEM"]["ntyp"]['v'] = len( stc.atom_spec)
		self._d["SYSTEM"]["nat"]['v'] = len( stc.atoms)

		return


