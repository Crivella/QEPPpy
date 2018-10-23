import sys
import numpy as np

from qepppy.classes.structure import structure as structure#, bravais_index as bi
from qepppy.classes.namelists import pw_nl, namelist_handler

class pwin( namelist_handler):
	def __init__( self, stc=None, fname="", **kwargs):
		self._d = pw_nl.copy()
		self.stc = None
		if stc:
			if isinstance( stc, structure):
				self._add_stc_( stc)
			else:
				raise Exception( "stc is not an instance of structure")
		if fname:
			self.pw_read( fname)

		super().__init__( **kwargs)


	def __str__( self):
		if not self.stc:
			raise Exception("No valid cell structure for the input file")
		try:
			c = super().__str__()
		except Exception as e:
			return str( e)
		c += "\n"
		try:
			c += self.stc.__str__()
		except Exception as e:
			return str( e)

		return c

	def __iadd__(self, other):
		if isinstance( other, structure):
			self._add_stc_( other)

		return self

	def _add_stc_( self, stc):
		self.stc = stc
		if isinstance( stc.bravais_n, int):
			self._d["SYSTEM"]["ibrav"] = stc.bravais_n
		else:
			#if not stc.a:
			raise Exception( "Must pass a valid cell structure")
			#self._d["SYSTEM"]["ibrav"] = 0

		self._d["SYSTEM"]["celldm(1)"] = stc.lp
		if self._d["SYSTEM"]["ibrav"] == 0:
			if not isinstance( stc.a, np.ndarray):
				raise Exception( "Basis vector must be set with ibrav = 0")

		#if stc.atom_spec_n != len( stc.atom_spec):
		#	raise Exception( "Invalide structure data, ntyp does not match")
		self._d["SYSTEM"]["ntyp"] = len( stc.atom_spec)
		self._d["SYSTEM"]["nat"] = len( stc.atoms)

		return

	def pw_read( self, fname=""):
		try:
			self.namelist_read( fname)
		except Exception as e:
			print( e)
		self.stc = structure( fname)

		return


