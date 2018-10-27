from .data_file_parser import data_file_parser as dfp


class bands( dfp):
	i=0
	e_units = 27.21138602
	__name__ = "bands"
	__toinit__ = ['kpt', 'egv', 'weight', 'occ', 'fermi', 'n_kpt', 'n_bnd', 'fermi_up=None; fermi_down', 'homo', 'lsda', 'noncolin', 'n_spin']
	data={
		'n_kpt':{'x':'text', 'f':'output//nk', 'n':None, 't':int},
		'n_bnd':{'x':'attr', 'f':'output//ks_energies/eigenvalues', 'n':'size', 't':int},
		'n_el':{'x':'text', 'f':'output//nelec', 'n':None, 't':float},
		'fermi':{'x':'text', 'f':'output//fermi_energy', 'n':None, 't':float},
		'fermi_s':{'x':'nodelist', 'f':'output//two_fermi_energies', 'n':None, 't':list},
		'homo':{'x':'text', 'f':'output//highestOccupiedLevel', 'n':None, 't':float},
		'lsda':{'x':'text', 'f':'output//lsda', 'n':None, 't':bool},
		'noncolin':{'x':'text', 'f':'output//noncolin', 'n':None, 't':bool},
		'kpt':{'x':'nodelist', 'f':'output//ks_energies/k_point', 'n':None, 't':list},
		'egv':{'x':'nodelist', 'f':'output//ks_energies/eigenvalues', 'n':None, 't':list},
		'occ':{'x':'nodelist', 'f':'output//ks_energies/occupations', 'n':None, 't':list},
		
		}
	"""
	def __init__( self, fname=""):
		self.i = 0
		if not fname:
			raise Exception( "Must initialize class giving the name of the .xml file")
		if not os.path.isfile( fname):
			raise IOError( "File '{}' does not exist".format( fname))
		self.fname = fname
		self._parse_xml_()
	"""

	def __str__( self):
		bnd = len(self.egv[0]['eigenvalues'])-1
		kpt_fmt = "\nkpt(#{{:5d}}):  " + "{:8.4f}"*3 + " [2pi/alat]"
		egv_fmt = "\nEigenvalues( eV):\n"+("  "+"{:12.6f}"*8+"\n")*int(bnd/8)
		egv_fmt += "  "+"{:12.6f}"*(bnd%8+1)+"\n"
		msg = ""
		for i in range( self.n_kpt):
			msg += kpt_fmt.format( *self.kpt[i]['k_point']).format( i)
			msg += egv_fmt.format( *self.egv[i]['eigenvalues'])
		return msg

	def __getitem__( self, key):
		if( isinstance( key, int)):
			if( 0<=key<self.n_kpt):
				return { 'kpt':self.kpt[key], 'egv':self.egv[key], 'occ':self.occ[key]} 
			else:
				raise Exception( "Index '{}' out of range".format( key))
		val = self.__dict__.get( key)
		if val: return val
		else: 
			raise KeyError( "'{}' object does not support key '{}'".format( self.__name__, key))

	def __iter__( self):
		return self

	def __next__( self):
		if( self.i < self.n_kpt):
			i = self.i
			self.i += 1
			return self[ i]
		else:
			self.i = 0
			raise StopIteration()

	def validate( self):
		if( self.n_kpt <= 0):
			raise Exception( "No kpt read from file '{}'.".format( self.fname))
		if( self.n_bnd <= 0):
			raise Exception( "No band read from file '{}'.".format( self.fname))
		if( not self.n_kpt == len( self.egv) == len( self.occ)):
			raise Exception( "Corrupted file. Number of kpoints does not match number egv or occ")
		return True









