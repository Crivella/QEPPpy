import numpy as np
from .data_file_parser import data_file_parser as dfp

import logging
logger = logging.getLogger( __name__)
logging.basicConfig( format='%(levelname)s: %(name)s\n%(message)s\n')


bravais_index={	'0':'free', '1':'simple cubic (sc)', '2':'face-centered cubic (fcc)', '3':'body-centered cubic (bcc)',
	'-3':'bcc more symm. axis', '4':'hexagonal', '5':'trigonal', '-5':'trigonal <111>', '6':'simple tetragonal (st)',
	'7':'body-centered tetragonal (bct)', '8':'orthorombic P', '9':'base-centered orthorombic (bco)', 
	'-9':'as 9 different axis', '10':'face-centered orthorombic', '11':'body-centered orthorombic', '12':'monoclinic P',
	'-12':'as 12 unique axis', '13':'base-centered monoclinic', '14':'triclinic'}

data={
	'ibrav':{'x':'attr', 'f':'output//atomic_structure', 'n':'bravais_index', 't':int},
	'alat':{'x':'attr', 'f':'output//atomic_structure', 'n':'alat', 't':float},
	'cell':{'x':'nodelist', 'f':'output//cell', 'n':None, 't':list},
	'recip':{'x':'nodelist', 'f':'output//reciprocal_lattice', 'n':None, 't':list},
	'atoms':{'x':'nodelist', 'f':'output//atom', 'n':'coord', 't':list},
	'atom_spec':{'x':'nodelist', 'f':'input//species', 'n':None, 't':list},
	'symm':{'x':'nodelist', 'f':'output//symmetry', 'n':None, 't':list}
	}




class structure( dfp):
	__name__ = "structure";
	def __init__( self, d={}, **kwargs):
		self.atom_p = 'bohr'
		self.cell_p = 'bohr'

		d.update( data)
		super().__init__( d=d, **kwargs)
		return

	def __str__( self, info=0):
		msg = super().__str__()
		fact = self.alat if self.cell_p == 'alat' else 1				
		if not self.ibrav or info > 0:
			msg += "CELL_PARAMETERS\n"
			for l in self.cell[0].values():
				for e in l:
					msg += "{:9.4f}".format( e * fact)
				msg += "\n"
			msg += "\n"

		msg += "ATOMIC_SPECIES\n"
		for s in self.atom_spec:
			msg += "{:6}{:12.4f}  {}".format( s['name'], s['mass'], s['pseudo_file'])
		msg += "\n\n"

		fact = self.alat if self.atom_p == 'alat' else 1
		msg += "ATOMIC_POSITIONS\n"
		for a in self.atoms:
			msg += "{:4}  ".format( a['name'])
			for c in a['coord']:
				msg += "{:10.5f}".format( c * fact)
			msg += "\n"
		msg += "\n"

		if info == 1:
			msg += "dbg info"

		return msg

	def pwin_read( self, fname=""):
		from .pwin import pw_in
		inp = pw_in( fname)
		inp.validate()

		name, x, y, z = inp.find( "X", "x", "y", "z", up="ATOMIC_POSITIONS")
		coord = [(a,b,c) for a,b,c in zip( x, y, z)]
		self.atoms = [{'name':n, 'coord':c} for n,c in zip( name, coord)]

		name, mass, pfile = inp.find( "X", "Mass_X", "PseudoPot_X", up="ATOMIC_SPECIES")
		self.atom_spec = [{'name':n, 'mass':m, 'pseudo_file':p} for n,m,p in zip( name, mass, pfile)]

		a1, a2, a3 = inp.find( "v1", "v2", "v3")
		self.cell = [{'a1':a1, 'a2':a2, 'a3':a3}]

		self.alat  = inp.find( "celldm(1)")
		self.ibrav = inp.find( "ibrav")
		self.atom_p = inp.find( "ATOMIC_POSITIONS")
		self.cell_p = inp.find( "CELL_PARAMETERS")
		#print( self.atom_p, self.cell_p)
		return

	def validate( self):
		l = [ self.ibrav, self.atom_spec, self.atoms]
		if any( a == None for a in l): return False
		#if self.atom_spec_n != len( self.atom_spec): return False

		for a in self.atoms:
			if not any( a['name'] == s['name'] for s in self.atom_spec):
				logger.warning( "Atoms in ATOMIC_POSITION do not match the type in ATOMIC_SPECIES")
				return False

		if self.ibrav == 0:
			if not isinstance( self.a, np.ndarray):
				return False
		return True and super().validate()

	def plot( self, 
		repX=1, repY=1, repZ=1, 
		cell=False, 
		bonds=True
		):
		"""
		Plot the crystal cell structure.
		reprX/Y/Z: repetitions of the cell along X/Y/Z (basis vector not Cartesian!!!)
		cell=True/False: plot the contour of the cell
		bonds=True/False: plot the chemical bonds between atoms
		"""
		import qepppy.classes.cell_graphic as cg
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		typ = [a['name'] for a in self.atoms]

		fact = self.alat if self.cell_p == 'alat' else 1
		v1 = np.array( self.cell[0]['a1']) * fact
		v2 = np.array( self.cell[0]['a2']) * fact
		v3 = np.array( self.cell[0]['a3']) * fact

		fact = self.alat if self.atom_p == 'alat' else 1
		if self.atom_p != 'crystal':
			L = np.array( [np.array(a['coord'])*fact for a in self.atoms])
		else:
			U = np.array([v1,v2,v3])
			L = np.array([U.dot( a['coord']) for a in self.atoms])

		typ = typ * repX
		typ = typ * repY
		typ = typ * repZ

		L = cg.cell_repetitions( L, v1, repX)
		L = cg.cell_repetitions( L, v2, repY)
		L = cg.cell_repetitions( L, v3, repZ)

		atom_list = list( zip( typ, L))

		for t in self.atom_spec:
			cg.draw_atoms( ax, atom_list, t['name'])

		if cell:
			cg.draw_cell( ax, v1, v2, v3)

		if bonds:
			cg.draw_bonds( ax, atom_list)

		#print( atom_list)

		plt.show()



