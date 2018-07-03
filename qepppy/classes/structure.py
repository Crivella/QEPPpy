import os.path
import re
import numpy as np
import xml.etree.ElementTree as ET

bravais_index={	'0':'free', '1':'simple cubic (sc)', '2':'face-centered cubic (fcc)', '3':'body-centered cubic (bcc)', \
		'-3':'bcc more symm. axis', '4':'hexagonal', '5':'trigonal', '-5':'trigonal <111>', '6':'simple tetragonal (st)', \
		'7':'body-centered tetragonal (bct)', '8':'orthorombic P', '9':'base-centered orthorombic (bco)', '-9':'as 9 different axis', \
		'10':'face-centered orthorombic', '11':'body-centered orthorombic', '12':'monoclinic P', '-12':'as 12 unique axis', \
		'13':'base-centered monoclinic', '14':'triclinic'}

class structure:
	__name__ = "structure"
	def __init__( self, fname=""):
		if not fname:
			raise Exception( "Must initialize class giving the name of the .xml file")
		if not os.path.isfile( fname):
			raise IOError( "File '{}' does not exist".format( fname))
		self.fname = fname
		self._parse_xml_()

	def __str__( self):
		return ''

	def __getattr__( self, key):
		if key in self.__dict__:
			return self.__dict__[key]
		else:
			raise AttributeError( "'{}' object has no attribute '{}'".format( self.__name__, key))

	def __getitem__( self, key):
		if( isinstance( key, str)):
			if key in self.__dict__:
				return self.__dict__[ key]
		raise KeyError( "'{}' object does not support key '{}'".format( self.__name__, key))

	def __enter__( self):
		return self

	def __exit__( self, *args):
		del self

	def _parse_xml_( self):
		root = ET.parse( self.fname).getroot()

		self.bravais = root.find("output//atomic_structure").get("bravais_index")
		if self.bravais in bravais_index:
			self.bravais = bravais_index[ self.bravais]
		self.lp      = float(root.find("output//atomic_structure").get("alat"))
		self.a= list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.find("output//cell").getchildren())))
		self.b= list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.find("output//reciprocal_lattice").getchildren())))

		
		self.symm = []
		for node in root.findall( "output//symmetry"):
			mnode   = node.find( "rotation")
			m       = np.array( list( map( float, filter( None, re.split( "\n +| ", mnode.text)))))
			m.shape = (3,3)
			n       = node.find( "info").get( "name")
			cl      = node.find( "info").get( "class")
			rk      = node.find( "rotation").get( "rank")
			self.symm.append( {'m':m, 'name':n, 'class':cl, 'rank':rk})









