import numpy as np
from .data_file_parser import data_file_parser as dfp


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
		d.update( data)
		super().__init__( d=d, **kwargs)
		return

	def __str__( self, info=0):
		msg = super().__str__()						
		if not self.ibrav or info > 0:
			msg += "CELL_PARAMETERS\n"
			for l in self.cell[0].values():
				for e in l:
					msg += "{:9.4f}".format( e * self.alat)
				msg += "\n"
			msg += "\n"

		msg += "ATOMIC_SPECIES\n"
		for s in self.atom_spec:
			msg += "{:6}{:12.4f}  {}".format( s['name'], s['mass'], s['pseudo_file'])
		msg += "\n\n"

		msg += "ATOMIC_POSITIONS\n"
		for a in self.atoms:
			msg += "{:4}  ".format( a['name'])
			for c in a['coord']:
				msg += "{:10.5f}".format( c)
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
		return

	def validate( self):
		l = [ self.ibrav, self.atom_spec, self.atoms]
		if any( a == None for a in l): return False
		#if self.atom_spec_n != len( self.atom_spec): return False

		for a in self.atoms:
			if not any( a['name'] == s['name'] for s in self.atom_spec): return False

		if self.ibrav == 0:
			if not isinstance( self.a, np.ndarray):
				return False
		return True and super().validate()

	def plot( self, repX=1, repY=1, repZ=1, cell=False):
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		#U = np.array([a for a in self.cell[0].values()])
		#print( U)
		#L0 = np.array([U.dot( a['coord']) for a in self.atoms])
		L = np.array( [a['coord'] for a in self.atoms])

		v1 = self.cell[0]['a1']
		v2 = self.cell[0]['a2']
		v3 = self.cell[0]['a3']

		L0 = L.copy()
		T = v1
		for n in range( 1, repX):
			L = np.vstack( ( L, L0 + T*n))

		L0 = L.copy()
		T = v2
		for n in range ( 1, repY):
			L = np.vstack( ( L, L0 + T*n))

		L0 = L.copy()
		T = v3
		for n in range( 1, repZ):
			L = np.vstack( ( L, L0 + T*n))

		X = L[:,0]
		Y = L[:,1]
		Z = L[:,2]
		ax.scatter( 
			X, Y, Z, 
			s=50,
			marker="o",
			depthshade=False,
			#c="black"
			)

		if cell:
			V = [v1,v2,v3]
			for n1 in range( 3):
				orig = np.array([0,0,0])
				v0 = V[n1]
				for n2 in range( 4):
					print( orig, v0)
					v = np.vstack( (orig, orig + v0))
					if n2 == n1:
						orig = np.array(V[(n2+1)%3]) + np.array(V[(n2+2)%3])
					else:
						orig = np.array( V[n2%3])
					ax.plot( v[:,0], v[:,1], v[:,2])
				ax.plot( v[:,0], v[:,1], v[:,2], c="black")
				
		#ax.plot( X, Y, Z)
		plt.show()












