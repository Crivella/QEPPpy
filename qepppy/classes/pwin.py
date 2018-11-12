from .structure import structure as structure
from .qe_in import qe_in
from .logger import *

@logger( )
class pw_in( qe_in, structure):
	"""
	Instance used to handle QE pw.x input files
	kwargs:
	 - inp = Name of the file to parse
	"""
	templ_file = "INPUT_PW.templ"
	def __init__( self, **kwargs):
		super().__init__( **kwargs)
		if 'parse' in kwargs:
			self.pwin_read( inp=self)
		return

	def __str__( self):
		msg = self.convert()
		if not self.check_used( 'ATOMIC_POSITIONS'):
			msg += structure.__str__( self)
		return msg
	
	def __iadd__( self, other):
		if isinstance( other, structure):
			self._add_stc_( other)

		return self

	def _add_stc_( self, stc):
		if isinstance( stc.ibrav, int):
			self.set_nl( nl="SYSTEM", k="ibrav", v=stc.ibrav)
		else:
			raise Exception( "Must pass a valid cell structure")

		self.set_nl( nl="SYSTEM", k="celldm(1)", v=stc.alat)
		self.set_nl( nl="SYSTEM", k="ntyp", v=len( stc.atom_spec))
		self.set_nl( nl="SYSTEM", k="nat", v=len( stc.atoms))

		for k, v in stc.__dict__.items():
			self.__dict__[k] = v

		return
	


