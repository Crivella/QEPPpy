import os.path
import xml.etree.ElementTree as ET

class bands:
	__name__ = "bands"
	def __init__( self, fname=""):
		self.i = 0
		if not fname:
			raise Exception( "Must initialize class giving the name of the .xml file")
		if not os.path.isfile( fname):
			raise IOError( "File '{}' does not exist".format( fname))
		self.fname = fname
		self._parse_xml_()

	def __getattr__( self, key):
		if key in self.__dict__:
			return self.__dict__[key]
		else:
			raise AttributeError( "'{}' object has no attribute '{}'".format( self.__name__, key))

	def __str__( self):
		bnd = len(self.egv[0])-1
		kpt_fmt = "\nkpt(#{{:5d}}):  " + "{:8.4f}"*3 + " [2pi/alat]"
		egv_fmt = "\nEigenvalues( eV):\n"+("  "+"{:12.6f}"*8+"\n")*int(bnd/8)
		egv_fmt += "  "+"{:12.6f}"*(bnd%8+1)+"\n"
		msg = ""
		for i in range( self.n_kpt):
			msg += kpt_fmt.format( *self.kpt[i]).format( i)
			msg += egv_fmt.format( *self.egv[i])
		return msg

	def __getitem__( self, key):
		if( isinstance( key, int)):
			if( 0<=key<self.n_kpt):
				return { 'kpt':self.kpt[key], 'egv':self.egv[key], 'occ':self.occ[key]} 
			else:
				raise Exception( "Index '{}' out of range".format( key))
		if( isinstance( key, str)):
			if key in self.__dict__:
				return self.__dict__[ key]
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

	def __enter__( self):
		return self

	def __exit__( self, *args):
		del self

	def _parse_xml_( self):
		root = ET.parse( self.fname).getroot()
		self.kpt = list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.findall("output//ks_energies/k_point"))))
		self.weight = list( map( float, map( lambda x: x.get('weight'), \
			root.findall("output//ks_energies/k_point"))))
		self.e_units = 27.21138602
		self.egv = list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.findall("output//ks_energies/eigenvalues"))))
		self.egv = list( map( lambda y: list(map( lambda x: x*self.e_units, y)), self.egv))
		self.occ = list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.findall("output//ks_energies/occupations"))))
		self.n_el    = float( root.find("output//nelec").text)

		#Fermi energy
		node         = root.find("output//fermi_energy")
		if node is None:
			self.fermi = float('NaN')
		else:
			self.fermi = float( node.text)

		#Fermi energy lsda
		node         = root.find("output//two_fermi_energies")
		if node is None:
			self.fermi = float('NaN')
		else:
			app = node.text.split( " ")
			self.fermi_up   = float( app[0])
			self.fermi_down = float( app[1])
			self.fermi_up   *= self.e_units
			self.fermi_down *= self.e_units

		#HOMO
		node         = root.find("output//highestOccupiedLevel")
		if node is None:
			self.homo = float('NaN')
		else:
			self.homo = float( node.text)



		self.n_kpt   = len( self.kpt)
		if( self.n_kpt <= 0):
			raise Exception( "No kpt read from file '{}'.".format( self.fname))
		self.n_bnd   = len( self.egv[0])
		if( self.n_bnd <= 0):
			raise Exception( "No band read from file '{}'.".format( self.fname))
		if( not self.n_kpt == len( self.egv) == len( self.occ)):
			raise Exception( "Corrupted file. Number of kpoints does not match number egv or occ")

		self.lsda    = root.find("output//lsda").text
		self.noncolin = root.find("output//noncolin").text
		if( self.lsda == 'true'):
			self.n_spin = 2
		if( self.noncolin == 'true'):
			self.n_spin = 4








