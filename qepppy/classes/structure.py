import os.path
import re
import numpy as np
import xml.etree.ElementTree as ET

bravais_index={	'0':'free', '1':'simple cubic (sc)', '2':'face-centered cubic (fcc)', '3':'body-centered cubic (bcc)',
	'-3':'bcc more symm. axis', '4':'hexagonal', '5':'trigonal', '-5':'trigonal <111>', '6':'simple tetragonal (st)',
	'7':'body-centered tetragonal (bct)', '8':'orthorombic P', '9':'base-centered orthorombic (bco)', 
	'-9':'as 9 different axis', '10':'face-centered orthorombic', '11':'body-centered orthorombic', '12':'monoclinic P',
	'-12':'as 12 unique axis', '13':'base-centered monoclinic', '14':'triclinic'}

class structure:
	__name__ = "structure"
	def __init__( self, fname="", **kwargs):
		#if not fname:
		#	raise Exception( "Must initialize class giving the name of the .xml file")
		self.bravais_n = None
		self.bravais = None
		self.lp = None
		self.a = None
		self.b = None
		self.atoms = None
		self.atom_spec = None
		self.atom_spec_n= None
		self.symm = None

		if fname:
			if not os.path.isfile( fname):
				raise IOError( "File '{}' does not exist".format( fname))
			self.fname = fname
			self._parse_xml_()
		else:
			for k, v in kwargs.items():
				#print( k, v)
				if k in self.__dict__:
					self.__dict__[k] = v
				else:
					raise Exception(" unrecognized keyword argument {}".format( k))


		return		

	def __str__( self, info=0):
		t=""
		if self.bravais_n == 0 or info > 0:
			t += "CELL_PARAMETERS"
			for l in self.a:
				for e in l:
					t += "{}".format(e)
				t += "\n"
			t += "\n\n"

		if self.atom_spec:
			t += "ATOMIC_SPECIES\n"
			for s in self.atom_spec:
				t += "{:6}{:12.4f}  {}".format(s['name'],s['mass'],s['pfile'])
			t += "\n\n"

		if self.atoms:
			t += "ATOMIC_POSITIONS\n"
			for a in self.atoms:
				t += "{:4}  ".format(a['name'])
				for c in a['coord']:
					t += "{:10.5f}".format(c)
				t += "\n"
			t += "\n"
		else:
			raise Exception( "No atom specified")

		if info == 1:
			t += "dbg info"

		return t

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

		self.bravais_n = root.find("output//atomic_structure").get("bravais_index")
		if self.bravais_n in bravais_index:
			self.bravais = bravais_index[ self.bravais_n]
		self.bravais_n = int( self.bravais_n)
		self.lp      = float(root.find("output//atomic_structure").get("alat"))
		self.a = list( map( lambda y: list(map( float, y)),
			map( lambda x: x.text.split(" "), root.find("output//cell").getchildren())))
		self.b = list( map( lambda y: list(map( float, y)),
			map( lambda x: x.text.split(" "), root.find("output//reciprocal_lattice").getchildren())))
		self.atoms = list( map( lambda x: {
			'name':x.get("name"), 
			'i':int(x.get("index")),
			'coord':list(map( float, x.text.split(" ")))
			}, root.findall("output//atom")
			))
		node = root.find("input/atomic_species")
		self.atom_spec_n = int( node.get("ntyp"))
		self.atom_spec = list( map( lambda x: {
			'mass':float(x.find("mass").text),
			'name':x.get("name"),
			'pfile':x.find("pseudo_file").text,
			'sm':float(x.find("starting_magnetization").text)
			},node.getchildren()
			))
		
		self.symm = []
		for node in root.findall( "output//symmetry"):
			mnode   = node.find( "rotation")
			m       = np.array( list( map( float, filter( None, re.split( "\n +| ", mnode.text)))))
			m.shape = (3,3)
			n       = node.find( "info").get( "name")
			cl      = node.find( "info").get( "class")
			rk      = node.find( "rotation").get( "rank")
			self.symm.append( {'m':m, 'name':n, 'class':cl, 'rank':rk})









