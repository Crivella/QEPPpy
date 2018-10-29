import os.path
import re
import numpy as np

from .data_file_parser import data_file_parser as dfp

bravais_index={	'0':'free', '1':'simple cubic (sc)', '2':'face-centered cubic (fcc)', '3':'body-centered cubic (bcc)',
	'-3':'bcc more symm. axis', '4':'hexagonal', '5':'trigonal', '-5':'trigonal <111>', '6':'simple tetragonal (st)',
	'7':'body-centered tetragonal (bct)', '8':'orthorombic P', '9':'base-centered orthorombic (bco)', 
	'-9':'as 9 different axis', '10':'face-centered orthorombic', '11':'body-centered orthorombic', '12':'monoclinic P',
	'-12':'as 12 unique axis', '13':'base-centered monoclinic', '14':'triclinic'}

data={
	'bravais_n':{'x':'attr', 'f':'output//atomic_structure', 'n':'bravais_index', 't':int},
	'lp':{'x':'attr', 'f':'output//atomic_structure', 'n':'alat', 't':float},
	'a':{'x':'nodelist', 'f':'output//cell', 'n':None, 't':list},
	'b':{'x':'nodelist', 'f':'output//reciprocal_lattice', 'n':None, 't':list},
	'atoms':{'x':'nodelist', 'f':'output//atom', 'n':'coord', 't':list},
	'atom_spec':{'x':'nodelist', 'f':'input//species', 'n':None, 't':list},
	'symm':{'x':'nodelist', 'f':'output//symmetry', 'n':None, 't':list}
	}

class structure( dfp):
	__name__ = "structure";
	def __init__( self, d={}, **kwargs):
		d.update( data)
		super().__init__( d=d, **kwargs)
		return

	def __str__( self, info=0):
		msg = super().__str__()						
		if self.bravais_n == 0 or info > 0:
			msg += "CELL_PARAMETERS"
			for l in self.a:
				for e in l:
					msg += "{}".format(e)
				msg += "\n"
			msg += "\n\n"

		msg += "ATOMIC_SPECIES\n"
		for s in self.atom_spec:
			msg += "{:6}{:12.4f}  {}".format(s['name'],s['mass'],s['pseudo_file'])
		msg += "\n\n"

		msg += "ATOMIC_POSITIONS\n"
		for a in self.atoms:
			msg += "{:4}  ".format(a['name'])
			for c in a['coord']:
				msg += "{:10.5f}".format(c)
			msg += "\n"
		msg += "\n"

		if info == 1:
			msg += "dbg info"

		return msg

	def pwin_read( self, fname=""):
		from .pwin import pw_in
		inp = pw_in( fname)
		inp.validate()

		atm = inp.find( "X", "x", "y", "z", up="ATOMIC_POSITIONS")
		self.atoms = []
		for n in range( len( atm[0])):
			self.atoms.append( {
				'name':atm[0][n],
				'coord':[atm[1][n],atm[2][n],atm[3][n]]
				})

		spc = inp.find( "X", "Mass_X", "PseudoPot_X", up="ATOMIC_SPECIES")
		self.atom_spec = []
		for n in range( len( spc[0])):
			self.atom_spec.append( {
				'name':spc[0][n],
				'mass':spc[1][n],
				'pseudo_file':spc[2][n]
				})

		self.a = np.array( inp.find( "v1", "v2", "v3"))
		return

	def validate( self):
		l = [ self.bravais_n, self.atom_spec, self.atoms]
		if any( a == None for a in l): return False
		#if self.atom_spec_n != len( self.atom_spec): return False

		for a in self.atoms:
			if not any( a['name'] == s['name'] for s in self.atom_spec): return False

		if self.bravais_n == 0:
			if not isinstance( self.a, np.ndarray):
				return False
		return True and super().validate()












