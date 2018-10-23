import os.path
import re
import numpy as np
import xml.etree.ElementTree as ET

bravais_index={	'0':'free', '1':'simple cubic (sc)', '2':'face-centered cubic (fcc)', '3':'body-centered cubic (bcc)',
	'-3':'bcc more symm. axis', '4':'hexagonal', '5':'trigonal', '-5':'trigonal <111>', '6':'simple tetragonal (st)',
	'7':'body-centered tetragonal (bct)', '8':'orthorombic P', '9':'base-centered orthorombic (bco)', 
	'-9':'as 9 different axis', '10':'face-centered orthorombic', '11':'body-centered orthorombic', '12':'monoclinic P',
	'-12':'as 12 unique axis', '13':'base-centered monoclinic', '14':'triclinic'}

class structure():
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
		#self.atom_spec_n= None
		self.symm = None

		if fname:
			if not os.path.isfile( fname):
				raise IOError( "File '{}' does not exist".format( fname))
			self.fname = fname
			if ".xml" in fname:
				self._parse_xml_()
			else:
				self.pw_read( fname)
		else:
			for k, v in kwargs.items():
				#print( k, v)
				if k in self.__dict__:
					self.__dict__[k] = v
				else:
					raise Exception( "Unrecognized keyword argument '{}'".format( k))
		
		if not self.validate():
			raise Exception( "Failed to initialize object '{}'.".format( self.__name__))

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

		t += "ATOMIC_SPECIES\n"
		for s in self.atom_spec:
			t += "{:6}{:12.4f}  {}".format(s['name'],s['mass'],s['pfile'])
		t += "\n\n"

		t += "ATOMIC_POSITIONS\n"
		for a in self.atoms:
			t += "{:4}  ".format(a['name'])
			for c in a['coord']:
				t += "{:10.5f}".format(c)
			t += "\n"
		t += "\n"

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

	def validate( self):
		if not self.bravais_n:														return False
		if not self.atom_spec:														return False
		#if self.atom_spec_n != len( self.atom_spec):			return False
		if not self.atoms:																return False

		for a in self.atoms:
			if not any( a['name'] == s['name'] for s in self.atom_spec): return False

		return True

	def pw_read( self, fname=""):
		with open(fname) as f:
			content = f.readlines()

		self.atoms=[]
		self.atom_spec=[]
		self.a=[]
		toup = None
		for l in content:
			lu = l.strip().upper()
			#print( lu, "ATOMIC_POSITIONS" in lu, toup)
			if "ATOMIC_POSITIONS" in lu:
				toup = self.atoms
				continue
			if "ATOMIC_SPECIES" in lu:
				toup = self.atom_spec
				continue
			if "CELL_PARAMETERS" in lu:
				toup = self.a
				continue
			if "K_POINTS" in lu:
				toup = None
				continue
			if "ibrav" in l:
				l = l[l.find("ibrav"):]
				self.bravais_n = int( l.split( "=")[1].split( ",")[0])
				continue
			if "celldm" in l:
				continue
			if toup != None:
				#print( l.split( " "))
				toup.append( list( filter( None, l.split( " "))))
				#print( toup)


		#print( self.atoms, self.atom_spec, self.a)
		self.atoms = list( map( lambda x: {
			'name':x[0],
			'coord':list(map(lambda y: float( y), x[1:4]))
			}, self.atoms))
		self.atom_spec = list( map( lambda x: {
			'name':x[0],
			'mass':float(x[1]),
			'pfile':x[2]
			}, self.atom_spec))
		if self.a:
			self.a = np.array( self.a)
			self.a.shape = (3,3)
		#print( self.atoms, self.atom_spec, self.a)



		return

	def _parse_xml_( self):
		#Read from data-file-schema.xml >6.2
		root = ET.parse( self.fname).getroot()

		#Read bravais lattice number
		self.bravais_n = root.find("output//atomic_structure").get("bravais_index")
		#If possible associate bravais lattice number to description
		if self.bravais_n in bravais_index:
			self.bravais = bravais_index[ self.bravais_n]
		self.bravais_n = int( self.bravais_n)
		#Read alat
		self.lp      = float(root.find("output//atomic_structure").get("alat"))
		#Read direct lattice vectors
		self.a = list( map( lambda y: list(map( float, y)),
			map( lambda x: x.text.split(" "), root.find("output//cell").getchildren())))
		#Read reciprocal lattice vectors
		self.b = list( map( lambda y: list(map( float, y)),
			map( lambda x: x.text.split(" "), root.find("output//reciprocal_lattice").getchildren())))
		#Read list of atom coordinates
		self.atoms = list( map( lambda x: {
			'name':x.get("name"), 
			'i':int(x.get("index")),
			'coord':list(map( float, x.text.split(" ")))
			}, root.findall("output//atom")
			))
		#Read list of atom types
		node = root.find("input/atomic_species")
		#self.atom_spec_n = int( node.get("ntyp"))
		self.atom_spec = list( map( lambda x: {
			'mass':float(x.find("mass").text),
			'name':x.get("name"),
			'pfile':x.find("pseudo_file").text,
			'sm':float(x.find("starting_magnetization").text)
			},node.getchildren()
			))
		
		#Read symmetries of the system
		self.symm = []
		for node in root.findall( "output//symmetry"):
			mnode   = node.find( "rotation")
			m       = np.array( list( map( float, filter( None, re.split( "\n +| ", mnode.text)))))
			m.shape = (3,3)
			n       = node.find( "info").get( "name")
			cl      = node.find( "info").get( "class")
			rk      = node.find( "rotation").get( "rank")
			self.symm.append( {'m':m, 'name':n, 'class':cl, 'rank':rk})









