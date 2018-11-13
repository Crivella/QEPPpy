import numpy as np
from .parser.data_file_parser import data_file_parser as dfp
from ..logger import *

bravais_index={	'0':'free', '1':'simple cubic (sc)', '2':'face-centered cubic (fcc)', '3':'body-centered cubic (bcc)',
	'-3':'bcc more symm. axis', '4':'hexagonal', '5':'trigonal', '-5':'trigonal <111>', '6':'simple tetragonal (st)',
	'7':'body-centered tetragonal (bct)', '8':'orthorombic P', '9':'base-centered orthorombic (bco)', '-9':'as 9 different axis',
	'91':'one-face base-centered orthorombic', '10':'face-centered orthorombic', '11':'body-centered orthorombic',
	'12':'monoclinic P', '-12':'as 12 unique axis', '13':'base-centered monoclinic', '14':'triclinic'}

data={
	'ibrav':{'x':'attr', 'f':'output//atomic_structure', 'n':'bravais_index', 't':int},
	'alat':{'x':'attr', 'f':'output//atomic_structure', 'n':'alat', 't':float},
	'cell':{'x':'nodelist', 'f':'output//cell', 'n':None, 't':list},
	'recip':{'x':'nodelist', 'f':'output//reciprocal_lattice', 'n':None, 't':list},
	'atoms':{'x':'nodelist', 'f':'output//atom', 'n':'coord', 't':list},
	'atom_spec':{'x':'nodelist', 'f':'input//species', 'n':None, 't':list},
	'symm':{'x':'nodelist', 'f':'output//symmetry', 'n':None, 't':list}
	}



@logger()
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

	def pwin_read( self, parse="", inp=None):
		if inp == None:
			if parse:
				from .pwin import pw_in
				inp = pw_in( parse=parse)
				inp.validate()
			else:
				raise error( "Must give a file name or a pw_in instance as args")

		name, x, y, z = inp.find( "X", "x", "y", "z", up="ATOMIC_POSITIONS")
		coord = [(a,b,c) for a,b,c in zip( x, y, z)]
		self.atoms = [{'name':n, 'coord':c} for n,c in zip( name, coord)]

		name, mass, pfile = inp.find( "X", "Mass_X", "PseudoPot_X", up="ATOMIC_SPECIES")
		self.atom_spec = [{'name':n, 'mass':m, 'pseudo_file':p} for n,m,p in zip( name, mass, pfile)]

		a1, a2, a3 = inp.find( "v1", "v2", "v3")
		self.cell = [{'a1':a1, 'a2':a2, 'a3':a3}]

		self.celldm = inp.find( "celldm")
		self.alat  = inp.find( "celldm(1)")
		self.ibrav = inp.find( "ibrav")
		self.atom_p = inp.find( "ATOMIC_POSITIONS")
		self.cell_p = inp.find( "CELL_PARAMETERS")
		#print( self.atom_p, self.cell_p)
		return

	def validate( self):
		ret = True
		if self.ibrav == None:
			warning.print( "ibrav is not set.")
			ret = False
		if self.atom_spec == None:
			warning.print( "List of atom types is not set.")
			ret = False
		if self.atoms == None:
			warning.print( "List of atomic positions is not set.")
			ret = False

		for a in self.atoms:
			if not any( a['name'] == s['name'] for s in self.atom_spec):
				warning.print( "Atoms in ATOMIC_POSITION do not match the type in ATOMIC_SPECIES")
				ret = False

		if self.ibrav == 0:
			if self.cell is None:
				warning.print( "Cell structure is not set with 'ibrav = 0'.")
				ret = False
		return ret and super().validate()

	def plot( self, 
		repX=1, repY=1, repZ=1, 
		cell=False, 
		bonds=True,
		recip=False,
		graph_lvl=1,
		):
		"""
		Plot the crystal cell structure.
		Args:
		 - reprX/Y/Z=1/1/1[or any positive integer]:
		        repetitions of the cell along X/Y/Z (basis vector not Cartesian!!!).
		        NOTE: They are 3 separate arguments repX, repY, repZ
		 - cell=False[or True]: plot the contour of the cell.
		 - bonds=True[or False]: plot the chemical bonds between atoms.
		 - recip=False[or True]: plot the Brilloiun Zone instead of the real cell.
		        If True all other flags are ignored.
		 - graph_lvl=1[or 0/2/3]:
		   - 0: Basic plot with circle dots as atoms and black lines as bonds.
		   - 1: Colored line as bonds (color of the nearest atom).
		   - 2: Use 3d spheres for atoms and bonds as in 1.
		   - 3: Use 3d spheres for atoms and cylinders for bonds.
		"""
		from ..cell_graphic import cell_graphic as cg
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		typ = [a['name'] for a in self.atoms]

		#print( self.cell)
		if self.cell == None:
			self._ibrav_to_cell_()
		else:
			if len( self.cell[0]['a1']) != 3:
				self._ibrav_to_cell_()

		if recip:
			try:
				test = isinstance( self.recip, list) and len( self.recip[0]['b1'])
			except:
				test = False
			if not test:
				self._cell_to_recip_()
			CELL = self.recip[0]
			fact = 1
			ax.set_xlabel("x (Bohr^-1)")
			ax.set_ylabel("y (Bohr^-1)")
			ax.set_zlabel("z (Bohr^-1)")
		else:
			CELL = self.cell[0]
			fact = self.alat if self.cell_p == 'alat' else 1
			ax.set_xlabel("x (Bohr)")
			ax.set_ylabel("y (Bohr)")
			ax.set_zlabel("z (Bohr)")

		v1 = np.array( list(CELL.values())[0]) * fact
		v2 = np.array( list(CELL.values())[1]) * fact
		v3 = np.array( list(CELL.values())[2]) * fact

		if recip:
			cg.draw_Wigner_Seitz( ax, v1, v2, v3)
			plt.show()
			return

		fact = self.alat if self.atom_p == 'alat' else 1
		if self.atom_p != 'crystal':
			L = np.array( [np.array(a['coord'])*fact for a in self.atoms])
		else:
			U = np.array([v1,v2,v3])
			L = np.array([U.dot( a['coord']) for a in self.atoms])

		for rep, v in zip( [repX,repY,repZ],[v1,v2,v3]):
			typ = typ * rep
			L = cg.cell_repetitions( L, v, rep)

		atom_list = list( zip( typ, L))

		if cell:
			cg.draw_cell( ax, v1, v2, v3)
		if bonds:
			cg.draw_bonds( ax, atom_list, graph_lvl=graph_lvl)
		for t in self.atom_spec:
			cg.draw_atoms( ax, atom_list, t['name'], graph_lvl=graph_lvl)

		ax.legend()
		plt.show()


	def _cell_to_recip_( self):
		fact = self.alat if self.cell_p == 'alat' else 1
		CELL = self.cell[0]
		v1 = np.array( CELL['a1']) * fact
		v2 = np.array( CELL['a2']) * fact
		v3 = np.array( CELL['a3']) * fact

		vol = v1.dot( np.cross(v2, v3))
		c = 1 / vol
		b1 = c * np.cross( v2, v3)
		b2 = c * np.cross( v3, v1)
		b3 = c * np.cross( v1, v2)
		self.recip= [{'b1':b1, 'b2':b2, 'b3':b3}]
		return


	def _ibrav_to_cell_( self):
		if self.ibrav == None:
			raise error( "Failed to generate cell structure from self.ibrav: self.ibrav not set.")
		"""
			v1 = np.array( [,,]) * lp
			v2 = np.array( [,,]) * lp
			v3 = np.array( [,,]) * lp
		"""

		lp = self.alat

		if   self.ibrav ==  1:
			v1 = np.array( [1,0,0]) * lp
			v2 = np.array( [0,1,0]) * lp
			v3 = np.array( [0,0,1]) * lp
		elif self.ibrav ==  2:
			v1 = np.array( [-1,0,1]) * lp/2
			v2 = np.array( [0,1,1]) * lp/2
			v3 = np.array( [-1,1,0]) * lp/2
		elif self.ibrav ==  3:
			v1 = np.array( [1,1,1]) * lp/2
			v2 = np.array( [-1,1,1]) * lp/2
			v3 = np.array( [-1,-1,1]) * lp/2
		elif self.ibrav == -3:
			v1 = np.array( [-1,1,1]) * lp/2
			v2 = np.array( [1,-1,1]) * lp/2
			v3 = np.array( [1,1,-1]) * lp/2
		elif self.ibrav ==  4:
			c = self.celldm[2]
			v1 = np.array( [1,0,0]) * lp
			v2 = np.array( [-1,np.sqrt(3),0]) * lp/2
			v3 = np.array( [0,0,c]) * lp
		elif self.ibrav ==  5:
			c = self.celldm[3]
			tx = (1-c)/2
			ty = (1-c)/6
			tz = (1+2*c)/3
			v1 = np.array( [tx,-ty,tz]) * lp
			v2 = np.array( [0,2*ty,tz]) * lp
			v3 = np.array( [-tx,-ty,tz]) * lp
		elif self.ibrav ==  -5:
			c = self.celldm[3]
			ty = (1-c)/6
			tz = (1+2*c)/3
			u = tz - 2*np.sqrt(2)*ty
			v = tz +np.sqrt(2)*ty
			v1 = np.array( [u,v,v]) * lp/np.sqrt(3)
			v2 = np.array( [v,u,v]) * lp/np.sqrt(3)
			v3 = np.array( [v,v,u]) * lp/np.sqrt(3)
		elif self.ibrav ==  6:
			c = self.celldm[2]
			v1 = np.array( [1,0,0]) * lp
			v2 = np.array( [0,1,0]) * lp
			v3 = np.array( [0,0,c]) * lp
		elif self.ibrav ==  7:
			c = self.celldm[2]
			v1 = np.array( [1,-1,c]) * lp/2
			v2 = np.array( [1,1,c]) * lp/2
			v3 = np.array( [-1,-1,c]) * lp/2
		elif self.ibrav ==  8:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array( [1,0,0]) * lp
			v2 = np.array( [0,b,0]) * lp
			v3 = np.array( [0,0,c]) * lp
		elif self.ibrav ==  9:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array( [1,b,0]) * lp/2
			v2 = np.array( [-1,b,0]) * lp/2
			v3 = np.array( [0,0,c]) * lp
		elif self.ibrav == -9:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array( [1,-b,0]) * lp/2
			v2 = np.array( [1,b,0]) * lp/2
			v3 = np.array( [0,0,c]) * lp
		elif self.ibrav == 91:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array( [1,0,0]) * lp
			v2 = np.array( [0,b,-c]) * lp/2
			v3 = np.array( [0,b,c]) * lp/2
		elif self.ibrav ==  10:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array( [1,0,c]) * lp/2
			v2 = np.array( [1,b,0]) * lp/2
			v3 = np.array( [0,b,c]) * lp/2
		elif self.ibrav ==  11:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array( [1,b,c]) * lp/2
			v2 = np.array( [-1,b,b]) * lp/2
			v3 = np.array( [-1,-b,c]) * lp/2
		elif self.ibrav ==  12:
			b = self.celldm[1]
			c = self.celldm[2]
			cab = self.celldm[3]
			sab = np.sqrt(1 - cab**2)
			v1 = np.array( [1,0,0]) * lp
			v2 = np.array( [b*cab,b*sab,0]) * lp
			v3 = np.array( [0,0,c]) * lp
		elif self.ibrav == -12:
			b = self.celldm[1]
			c = self.celldm[2]
			cac = self.celldm[4]
			sac = np.sqrt(1 - cac**2)
			v1 = np.array( [1,0,0]) * lp
			v2 = np.array( [0,b,0]) * lp
			v3 = np.array( [c*cac,0,c*sac]) * lp
		elif self.ibrav ==  13:
			b = self.celldm[1]
			c = self.celldm[2]
			cg = self.celldm[3]
			sg = np.sqrt(1 - cg**2)
			v1 = np.array( [1,0,-c]) * lp/2
			v2 = np.array( [b*cg,b*sg,0]) * lp/2
			v3 = np.array( [1,0,c]) * lp
		elif self.ibrav ==  -13:
			b = self.celldm[1]
			c = self.celldm[2]
			cb = self.celldm[4]
			sb = np.sqrt(1 - cb**2)
			v1 = np.array( [1,-b,0]) * lp/2
			v2 = np.array( [1,b,0]) * lp/2
			v3 = np.array( [c*cb,0,c*sb]) * lp
		elif self.ibrav ==  14:
			b = self.celldm[1]
			c = self.celldm[2]
			cbc = self.celldm[3]
			cac = self.celldm[4]
			cab = self.celldm[5]
			cg = cab
			sg = np.sqrt(1 - cg**2)
			v1 = np.array( [1,0,0]) * lp
			v2 = np.array( [b*cg,b*sg,0]) * lp
			v3 = np.array( [c*cac,
				c*(cbc-cac*cg)/sg,
				c*np.sqrt(1+2*cbc*cac*cg-cbc**2-cac**2-cg**2)/sg]) * lp

		self.cell_p = 'bohr'
		self.cell=[{'a1':v1, 'a2':v2, 'a3':v3}]


