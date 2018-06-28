import os.path
import xml.etree.ElementTree as ET

from .bands import bands as bnd
#from .eigenv  import egv

class pwout( bnd):
	__name__ = "pwout"
	def __init__( self, fname=""):
		self.i = 0
		if not fname:
			raise Exception( "Must initialize class giving the name of the .xml file")
		if not os.path.isfile( fname):
			raise IOError( "File '{}' does not exist".format( fname))
		self.fname = fname
		self._parse_xml_()
		#bnd._parse_xml_( self)
		#egv._parse_xml( self)

	#def __getattr__( self, key):
	#	raise AttributeError( "The attribute '{}' does not exist in the class '{}'".format( key, __name__))


	def __str__( self):
		#bnd = len(self.egv[0])-1
		#kpt_fmt = "\nkpt:\t" + "{:8.4f}"*3
		#egv_fmt = "\nEigenvalues( eV):\n"+("\t"+"{:10.6f}"*8+"\n")*int(bnd/8)
		#egv_fmt += "\t"+"{:10.6f}"*(bnd%8+1)+"\n"
		msg = ""
		#for i in range( self.n_kpt):
		#	msg += kpt_fmt.format( *self.kpt[i])
		#	msg += egv_fmt.format( *self.egv[i])
		msg += bnd.__str__( self)
		return msg

#########################################################33 boh
	def __getitem__( self, key):
		if( isinstance( key, int)):
			return bnd.__getitem__( self, key)
		else:
			raise TypeError( "'{}' object does not support indexing for key '{}'".format( self.__name__, key))

	def __iter__( self):
		return bnd.__iter__( self)

	def __enter__( self):
		return self

	def __exit__( self, *args):
		del self

	def _parse_xml_( self):
		root = ET.parse( self.fname).getroot()

		bnd._parse_xml_( self)
		self.bravais = root.find("output//atomic_structure").get("bravais_index")
		self.lp      = float(root.find("output//atomic_structure").get("alat"))
		self.a= list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.find("output//cell").getchildren())))
		self.b= list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.find("output//reciprocal_lattice").getchildren())))

		#self.n_kpt   = int( root.find("output//nk").text)
		#self.n_bnd   = int( root.find("output//nbnd").text)
		self.n_kpt   = len( self.kpt)
		if( self.n_kpt <= 0):
			raise Exception( "No kpt read from file '{}'.".format( self.fname))
		self.n_bnd   = len( self.egv[0])
		if( self.n_bnd <= 0):
			raise Exception( "No band read from file '{}'.".format( self.fname))
		if( not self.n_kpt == len( self.egv) == len( self.occ)):
			raise Exception( "Corrupted file. Number of kpoints does not match number egv or occ")








