import numpy as np
from .parser.data_file_parser import data_file_parser as dfp
from ..logger import logger, warning, error
from .._decorators import store_property

bravais_index={	
	'0':'free',
	'1':'simple cubic (sc)',
	'2':'face-centered cubic (fcc)',
	'3':'body-centered cubic (bcc)',
	'-3':'bcc more symm. axis',
	'4':'hexagonal',
	'5':'trigonal',
	'-5':'trigonal <111>',
	'6':'simple tetragonal (st)',
	'7':'body-centered tetragonal (bct)',
	'8':'orthorombic P',
	'9':'base-centered orthorombic (bco)',
	'-9':'as 9 different axis',
	'91':'one-face base-centered orthorombic',
	'10':'face-centered orthorombic',
	'11':'body-centered orthorombic',
	'12':'monoclinic P',
	'-12':'as 12 unique axis',
	'13':'base-centered monoclinic',
	'14':'triclinic'
	}

data={
	'n_atoms':{
		'xml_ptype':'attr', 
		'xml_search_string':'output//atomic_structure', 
		'extra_name':'nat', 
		'res_type':int,
		'outfile_regex':r'number of atoms/cell\s*=\s*'
		},
	'n_types':{
		'xml_ptype':'attr', 
		'xml_search_string':'output//atomic_species', 
		'extra_name':'ntyp', 
		'res_type':int,
		'outfile_regex':r'number of atomic types\s*=\s*'
		},
	'ibrav':{
		'xml_ptype':'attr', 
		'xml_search_string':'output//atomic_structure', 
		'extra_name':'bravais_index', 
		'res_type':int,
		'outfile_regex':r'bravais-lattice index\s*='
		},
	'alat':{
		'xml_ptype':'attr', 
		'xml_search_string':'output//atomic_structure', 
		'extra_name':'alat', 
		'res_type':float,
		'outfile_regex':r'lattice parameter \(alat\)\s*='
		},
	'cell_p':{
		'res_type':str,
		'outfile_regex':r'cart\. coord\. in units of (?P<flag>.*)\)'
		},
	'_cell':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//cell', 
		'extra_name':None, 
		'res_type':list,
		'outfile_regex':
			r'\s*a\(1\) = \((?P<a1>[\s\d.\-]*)\)\s*\n' + 
			r'\s*a\(2\) = \((?P<a2>[\s\d.\-]*)\)\s*\n' + 
			r'\s*a\(3\) = \((?P<a3>[\s\d.\-]*)\)\s*\n'
		},
	'_recip':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//reciprocal_lattice', 
		'extra_name':None, 
		'res_type':list,
		'outfile_regex':
			r'\s*b\(1\) = \((?P<b1>[\s\d.\-]*)\)\s*\n' + 
			r'\s*b\(2\) = \((?P<b2>[\s\d.\-]*)\)\s*\n' + 
			r'\s*b\(3\) = \((?P<b3>[\s\d.\-]*)\)\s*\n'
		},
	'atom_p':{
		'res_type':str,
		'outfile_regex':r'positions \((?P<flag>.*) units\)'
		},
	'_atoms':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//atom', 
		'extra_name':'coord', 
		'res_type':list,
		'outfile_regex':r'\d[\t ]+(?P<name>[\w]+).*\((?P<index>[ \d]+)\) = \((?P<coord>[ \d\.\-]+)\)'
		},
	'_atom_spec':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'input//species', 
		'extra_name':None, 
		'res_type':list,
		'outfile_regex':r'\s*(?P<name>\w+)\s+(?P<valence>[\d\.]+)\s+(?P<mass>[\d\.]+)\s+(?P<pseudo_file>\w+\s*\([ \d\.]+\))'
		},
	'_symm':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//symmetry', 
		'extra_name':None, 
		'res_type':list,
		'outfile_regex':
			r'isym =\s*\d{1,2}\s*(?P<name>[\S ]*)\n\s*' +
			r'cryst.\s*s\([\s\d]{2}\) = ' +
			r'(?P<rotation>(\(.*\)\s*){3})'
		}
	}

def _recip_space_(v1, v2, v3):
	vol = v1.dot(np.cross(v2, v3))
	c   = 1. / vol
	b1  = c * np.cross(v2, v3)
	b2  = c * np.cross(v3, v1)
	b3  = c * np.cross(v1, v2)
	return np.array([b1,b2,b3])

# @logger()
class structure(dfp):
	__name__ = "structure";
	def __init__(self, d={}, **kwargs):
		self.atom_p = 'bohr'
		self.cell_p = 'bohr'

		d.update(data)
		super().__init__(d=d, **kwargs)
		return

	def _format_cell_(self, info):
		if self.ibrav and info == 0:
			return ""
		msg = "CELL_PARAMETERS\n"
		for l in self.cell:
			for e in l:
				msg += "{:9.4f}".format(e)
			msg += "\n"
		msg += "\n"
		return msg

	def _format_atom_spec_(self):
		msg = "ATOMIC_SPECIES\n"
		for s in self._atom_spec:
			msg += "{:6}{:12.4f}  {}".format(s['name'], s['mass'], s['pseudo_file'])
			msg += "\n"
		msg += "\n\n"
		return msg

	def _format_atom_pos_(self):
		msg = "ATOMIC_POSITIONS\n"
		for coord,name in zip(self.atoms_coord_cart, self.atoms_typ):
			msg += "{:4}  ".format(name)
			for c in coord:
				msg += "{:10.5f}".format(c)
			msg += "\n"
		msg += "\n"
		return msg

	def __str__(self, info=0):
		msg = super().__str__()
		msg += self._format_cell_(info)
		msg += self._format_atom_spec_()
		msg += self._format_atom_pos_()

		return msg

	@property
	@store_property
	def atoms_coord_cart(self):
		fact = self.alat if self.atom_p == 'alat' else 1
		res = np.array([a['coord'] for a in self._atoms]) * fact
		n = self.n_atoms
		if n and len(res) == n*2:
			res = res[0:n,:]
		if self.atom_p == 'crystal':
			res = np.dot(self.cell.T, res.T).T
		return res

	@property
	@store_property
	def atoms_coord_cryst(self):
		fact = self.alat if self.atom_p == 'alat' else 1
		res = np.array([a['coord'] for a in self._atoms]) * fact
		n = self.n_atoms
		if len(res) == n*2:
				res = res[n:n*2,:]
		else:
			if self.atom_p != 'crystal':
				# import scipy.linalg
				cell = np.mat(self.cell.T)
				res = np.dot(cell.I, res.T).T
				# raise NotImplementedError()
		return res

	@property
	@store_property
	def atoms_group_coord_cryst(self):
		return list([np.array(self.atoms_coord_cryst[np.where(np.array(self.atoms_typ) == a)[0]]) for a in self.all_atoms_typ])

	@property
	@store_property
	def atoms_typ(self):
		res = list([a['name'] for a in self._atoms])
		# if len(res) > self.n_types:
		# 	raise ValueError("Found {} types vs ntyp = {}".format(res, self.n_types))
		return res

	@property
	@store_property
	def atoms_mass(self):
		return np.array([a['mass'] for a in self._atom_spec])

	@property
	@store_property
	def atoms_pseudo(self):
		return np.array([a['pseudo_file'] for a in self._atom_spec])

	@property
	@store_property
	def all_atoms_typ(self):
		return set(self.atoms_typ)

	@property
	@store_property
	def cell(self):
		res = None
		res =  np.array(list(self._cell[0].values()))
		if res.size != 9:
			return self._ibrav_to_cell_()
		fact = self.alat if self.cell_p == 'alat' else 1
		res *= fact
		return res

	@property
	@store_property
	def recip(self):
		try:
			res =  np.array(list(self._recip[0].values()))
		except Exception as e:
			res = self._cell_to_recip_()
		return res
	
	@property
	@store_property
	def symm_matrix(self):
		t = type(self._symm[0]['rotation'])
		if t == np.ndarray:
			res = np.array([a['rotation'].reshape(3,3) for a in self._symm])
		elif t == str:
			g = lambda x: x.replace('(',' ').replace(')',' ').replace('\n',' ')
			res = np.array([np.fromstring(g(a['rotation']), sep=' ').reshape(3,3) for a in self._symm])
		else:
			raise NotImplementedError()
		return res

	@property
	@store_property
	def symm_name(self):
		return list([a['name'] for a in self._symm])

	def pwin_read(self, parse="", inp=None):
		if inp == None:
			if parse:
				from .pwin import pw_in
				inp = pw_in(parse=parse)
				inp.validate()
			else:
				raise error("Must give a file name or a pw_in instance as args")

		name, x, y, z = inp.find("X", "x", "y", "z", up="ATOMIC_POSITIONS")
		coord = [(a,b,c) for a,b,c in zip(x, y, z)]
		self._atoms = [{'name':n, 'coord':c} for n,c in zip(name, coord)]

		name, mass, pfile = inp.find("X", "Mass_X", "PseudoPot_X", up="ATOMIC_SPECIES")
		self._atom_spec = [{'name':n, 'mass':m, 'pseudo_file':p} for n,m,p in zip(name, mass, pfile)]

		a1, a2, a3 = inp.find("v1", "v2", "v3")
		self._cell = [{'a1':a1, 'a2':a2, 'a3':a3}]

		self.celldm  = inp.find("celldm")
		self.alat    = inp.find("celldm(1)")
		self.ibrav   = inp.find("ibrav")
		self.atom_p  = inp.find("ATOMIC_POSITIONS")
		self.cell_p  = inp.find("CELL_PARAMETERS")
		self.n_atoms = inp.find("nat")
		self.n_types = inp.find("ntyp")

	def validate(self):
		ret = True
		if self.ibrav == None:
			warning.print("ibrav is not set.")
			ret = False
		if self._atom_spec == None:
			warning.print("List of atom types is not set.")
			ret = False
		if self._atoms == None:
			warning.print("List of atomic positions is not set.")
			ret = False

		for typ in self.atoms_typ:
			if not typ in self.all_atoms_typ:
				warning.print("Atoms in ATOMIC_POSITION do not match the type in ATOMIC_SPECIES")
				ret = False

		if self.ibrav == 0:
			if self.cell is None:
				warning.print("Cell structure is not set with 'ibrav = 0'.")
				ret = False
		return ret and super().validate()

	def plot(
		self, 
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
		from .. import cell_graphic as cg
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		typ = self.atoms_typ

		if recip:
			ax.set_xlabel(r"$x (Bohr^{-1})$")
			ax.set_ylabel(r"$y (Bohr^{-1})$")
			ax.set_zlabel(r"$z (Bohr^{-1})$")
			cg.draw_Wigner_Seitz(ax, self.recip)
			plt.show()
			return
		else:
			ax.set_xlabel("x (Bohr)")
			ax.set_ylabel("y (Bohr)")
			ax.set_zlabel("z (Bohr)")


		L = self.atoms_coord_cart
		for rep, v in zip([repX,repY,repZ],self.cell):
			typ = typ * rep
			L = cg.cell_repetitions(L, v, rep)

		cg.draw_atoms(ax, L, typ, graph_lvl=graph_lvl)
		if cell:
			cg.draw_cell(ax, *self.cell)
		if bonds:
			cg.draw_bonds(ax, L, typ, graph_lvl=graph_lvl)
		ax.legend()
		plt.show()

	def _cell_to_recip_(self):
		return _recip_space_(*self.cell)

	def _ibrav_to_cell_(self):
		if self.ibrav == None:
			raise error("Failed to generate cell structure from self.ibrav: self.ibrav not set.")
		"""
			v1 = np.array([,,]) * lp
			v2 = np.array([,,]) * lp
			v3 = np.array([,,]) * lp
		"""

		lp = self.alat

		if   self.ibrav ==  1:
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,1,0]) * lp
			v3 = np.array([0,0,1]) * lp
		elif self.ibrav ==  2:
			v1 = np.array([-1,0,1]) * lp/2
			v2 = np.array([0,1,1]) * lp/2
			v3 = np.array([-1,1,0]) * lp/2
		elif self.ibrav ==  3:
			v1 = np.array([1,1,1]) * lp/2
			v2 = np.array([-1,1,1]) * lp/2
			v3 = np.array([-1,-1,1]) * lp/2
		elif self.ibrav == -3:
			v1 = np.array([-1,1,1]) * lp/2
			v2 = np.array([1,-1,1]) * lp/2
			v3 = np.array([1,1,-1]) * lp/2
		elif self.ibrav ==  4:
			c = self.celldm[2]
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([-1,np.sqrt(3),0]) * lp/2
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav ==  5:
			c = self.celldm[3]
			tx = (1-c)/2
			ty = (1-c)/6
			tz = (1+2*c)/3
			v1 = np.array([tx,-ty,tz]) * lp
			v2 = np.array([0,2*ty,tz]) * lp
			v3 = np.array([-tx,-ty,tz]) * lp
		elif self.ibrav ==  -5:
			c = self.celldm[3]
			ty = (1-c)/6
			tz = (1+2*c)/3
			u = tz - 2*np.sqrt(2)*ty
			v = tz +np.sqrt(2)*ty
			v1 = np.array([u,v,v]) * lp/np.sqrt(3)
			v2 = np.array([v,u,v]) * lp/np.sqrt(3)
			v3 = np.array([v,v,u]) * lp/np.sqrt(3)
		elif self.ibrav ==  6:
			c = self.celldm[2]
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,1,0]) * lp
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav ==  7:
			c = self.celldm[2]
			v1 = np.array([1,-1,c]) * lp/2
			v2 = np.array([1,1,c]) * lp/2
			v3 = np.array([-1,-1,c]) * lp/2
		elif self.ibrav ==  8:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,b,0]) * lp
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav ==  9:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array([1,b,0]) * lp/2
			v2 = np.array([-1,b,0]) * lp/2
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav == -9:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array([1,-b,0]) * lp/2
			v2 = np.array([1,b,0]) * lp/2
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav == 91:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,b,-c]) * lp/2
			v3 = np.array([0,b,c]) * lp/2
		elif self.ibrav ==  10:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array([1,0,c]) * lp/2
			v2 = np.array([1,b,0]) * lp/2
			v3 = np.array([0,b,c]) * lp/2
		elif self.ibrav ==  11:
			b = self.celldm[1]
			c = self.celldm[2]
			v1 = np.array([1,b,c]) * lp/2
			v2 = np.array([-1,b,b]) * lp/2
			v3 = np.array([-1,-b,c]) * lp/2
		elif self.ibrav ==  12:
			b = self.celldm[1]
			c = self.celldm[2]
			cab = self.celldm[3]
			sab = np.sqrt(1 - cab**2)
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([b*cab,b*sab,0]) * lp
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav == -12:
			b = self.celldm[1]
			c = self.celldm[2]
			cac = self.celldm[4]
			sac = np.sqrt(1 - cac**2)
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,b,0]) * lp
			v3 = np.array([c*cac,0,c*sac]) * lp
		elif self.ibrav ==  13:
			b = self.celldm[1]
			c = self.celldm[2]
			cg = self.celldm[3]
			sg = np.sqrt(1 - cg**2)
			v1 = np.array([1,0,-c]) * lp/2
			v2 = np.array([b*cg,b*sg,0]) * lp/2
			v3 = np.array([1,0,c]) * lp
		elif self.ibrav ==  -13:
			b = self.celldm[1]
			c = self.celldm[2]
			cb = self.celldm[4]
			sb = np.sqrt(1 - cb**2)
			v1 = np.array([1,-b,0]) * lp/2
			v2 = np.array([1,b,0]) * lp/2
			v3 = np.array([c*cb,0,c*sb]) * lp
		elif self.ibrav ==  14:
			b = self.celldm[1]
			c = self.celldm[2]
			cbc = self.celldm[3]
			cac = self.celldm[4]
			cab = self.celldm[5]
			cg = cab
			sg = np.sqrt(1 - cg**2)
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([b*cg,b*sg,0]) * lp
			v3 = np.array([c*cac,
				c*(cbc-cac*cg)/sg,
				c*np.sqrt(1+2*cbc*cac*cg-cbc**2-cac**2-cg**2)/sg]) * lp

		self.cell_p = 'bohr'
		return np.array([v1,v2,v3])



