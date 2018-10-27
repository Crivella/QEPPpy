import os.path
import re
import numpy as np

from .data_file_parser import data_file_parser as dfp

bravais_index={	'0':'free', '1':'simple cubic (sc)', '2':'face-centered cubic (fcc)', '3':'body-centered cubic (bcc)',
	'-3':'bcc more symm. axis', '4':'hexagonal', '5':'trigonal', '-5':'trigonal <111>', '6':'simple tetragonal (st)',
	'7':'body-centered tetragonal (bct)', '8':'orthorombic P', '9':'base-centered orthorombic (bco)', 
	'-9':'as 9 different axis', '10':'face-centered orthorombic', '11':'body-centered orthorombic', '12':'monoclinic P',
	'-12':'as 12 unique axis', '13':'base-centered monoclinic', '14':'triclinic'}

class structure( dfp):
	__name__ = "structure";
	__toinit__ = [ 'bravais_n', 'bravais', 'lp', 'a', 'b', 'atoms', 'atom_spec', 'symm'] 
	data={
		'bravais_n':{'x':'attr', 'f':'output//atomic_structure', 'n':'bravais_index', 't':int},
		'lp':{'x':'attr', 'f':'output//atomic_structure', 'n':'alat', 't':float},
		'a':{'x':'allchild', 'f':'output//cell', 'n':None, 't':np.array},
		'b':{'x':'allchild', 'f':'output//reciprocal_lattice', 'n':None, 't':np.array},
		'atoms':{'x':'allchild', 'f':'output//atomic_positions', 'n':'coord', 't':list},
		'atom_spec':{'x':'allchild', 'f':'input/atomic_species', 'n':'coord', 't':list},
		'symm':{'x':'nodelistnested', 'f':'output//symmetry', 'n':'rotation', 't':list}
		}
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
			t += "{:6}{:12.4f}  {}".format(s['name'],s['mass'],s['pseudo_file'])
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



	def __getitem__( self, key):
		val = self.__dict__.get( key)
		if val: return val
		else: 
			raise KeyError( "'{}' object does not support key '{}'".format( self.__name__, key))


	def validate( self):
		if not self.bravais_n:  return False
		if not self.atom_spec:  return False
		#if self.atom_spec_n != len( self.atom_spec):			return False
		if not self.atoms:      return False

		for a in self.atoms:
			if not any( a['name'] == s['name'] for s in self.atom_spec): return False

		if self.bravais_n == 0:
			if not isinstance( self.a, np.ndarray):
				return False
				#raise Exception( "Basis vector must be set with ibrav = 0")

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









