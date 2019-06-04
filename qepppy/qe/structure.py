import numpy as np
from .parser.data_file_parser import data_file_parser as dfp
from ..errors import ValidateError
from .._decorators import store_property
from .._cell import _cell as cell
from .. import utils

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
	'_n_atoms':{
		'xml_ptype':'attr', 
		'xml_search_string':'output//atomic_structure', 
		'extra_name':'nat', 
		'res_type':int,
		'outfile_regex':r'number of atoms/cell\s*=\s*'
		},
	'_n_types':{
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
	'_cell_p':{
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
	'_atom_p':{
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
		},
	'_fft_grid':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//fft_smooth', 
		'extra_name':None, 
		'res_type':list,
		'outfile_regex':
			r'FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*)'
		}
	}

# @logger()
class structure(dfp, cell):
	__name__ = "structure";
	def __init__(self, d={}, **kwargs):
		self._atom_p = 'bohr'
		self._cell_p = 'bohr'

		d.update(data)
		super().__init__(d=d, **kwargs)

	def _format_cell_(self, info):
		if self.ibrav and info == 0:
			return ""
		msg = "CELL_PARAMETERS\n"
		for l in self.direct:
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
	def _atoms_coord_cart(self):
		res = np.array([a['coord'] for a in self._atoms]) * self.atom_p
		n = self._n_atoms
		if n and len(res) == n*2:
			res = res[0:n,:]
		if self._atom_p == 'crystal':
			res = np.dot(self.direct.T, res.T).T
		return res

	@property
	def _atoms_coord_cryst(self):
		res = np.array([a['coord'] for a in self._atoms]) * self.atom_p
		n = self._n_atoms
		if len(res) == n*2:
				res = res[n:n*2,:]
		else:
			if self._atom_p != 'crystal':
				res = np.dot(self.recipr/(2*np.pi), res.T).T
		return res

	@property
	def _atoms_typ(self):
		res = list([a['name'] for a in self._atoms])
		n = self.n_atoms

		# Cut the crystal coordinates that are taken using regex on output files
		if len(res) == n*2:
			res = res[:n]
		return res

	@property
	def _atoms_mass(self):
		return np.array([a['mass'] for a in self._atom_spec])

	@property
	def _atoms_pseudo(self):
		return np.array([a['pseudo_file'] for a in self._atom_spec])

	@property
	def _all_atoms_typ(self):
		return list([a['name'] for a in self._atom_spec])

	@property
	def atom_p(self):
		"""Conversion factor for atom coordinates to atomic units"""
		if self._atom_p == 'alat':
			return self.alat
		if self._atom_p == 'angstrom':
			return 1/0.529177
		return 1

	@property
	def cell_p(self):
		"""Conversion factor for cell vetors to atomic units"""
		if self._cell_p == 'alat':
			return self.alat
		if self._cell_p == 'angstrom':
			return 1/0.529177
		return 1

	@cell_p.setter
	def cell_p(self, value):
		self._cell_p = value

	@property
	def _direct(self):
		res =  np.array(list(self._cell[0].values()))
		if res.size != 9:
			return self._ibrav_to_cell_()
		res *= self.cell_p
		return res

	@property
	def _recipr(self):
		return utils.recipr_base(self._direct)
	
	@property
	@store_property
	def symm_matrix(self):
		"""List of symmetry operation matrices"""
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
		"""List of symmetry operation names"""
		return list([a['name'] for a in self._symm])

	@property
	@store_property
	def fft_dense_grid(self):
		"""FFT Grid shape: np.ndarray of shape (3,)"""
		return np.array([self._fft_grid[0]['nr1'], self._fft_grid[0]['nr2'], self._fft_grid[0]['nr3']], dtype='int')

	def validate(self):
		# if self.ibrav == None:
		# 	raise ValidateError("ibrav is not set.")
		if self._atom_spec == None:
			raise ValidateError("List of atom types is not set.")
		if self._atoms == None:
			raise ValidateError("List of atomic positions is not set.")

		for typ in self.atoms_typ:
			if not typ in self.all_atoms_typ:
				raise ValidateError("Atoms in ATOMIC_POSITION do not match the type in ATOMIC_SPECIES")

		if self.ibrav == 0:
			if self.direct is None:
				raise ValidateError("Cell structure is not set with 'ibrav = 0'.")

		super().validate()


	def _ibrav_to_cell_(self):
		if self.ibrav == None:
			raise error("Failed to generate cell structure from self.ibrav: self.ibrav not set.")

		lp = self.alat

		if len(self.celldm) > 1:
			b = self.celldm[1]
		if len(self.celldm) > 2:
			c = self.celldm[2]
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
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,1,0]) * lp
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav ==  7:
			v1 = np.array([1,-1,c]) * lp/2
			v2 = np.array([1,1,c]) * lp/2
			v3 = np.array([-1,-1,c]) * lp/2
		elif self.ibrav ==  8:
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,b,0]) * lp
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav ==  9:
			v1 = np.array([1,b,0]) * lp/2
			v2 = np.array([-1,b,0]) * lp/2
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav == -9:
			v1 = np.array([1,-b,0]) * lp/2
			v2 = np.array([1,b,0]) * lp/2
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav == 91:
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,b,-c]) * lp/2
			v3 = np.array([0,b,c]) * lp/2
		elif self.ibrav ==  10:
			v1 = np.array([1,0,c]) * lp/2
			v2 = np.array([1,b,0]) * lp/2
			v3 = np.array([0,b,c]) * lp/2
		elif self.ibrav ==  11:
			v1 = np.array([1,b,c]) * lp/2
			v2 = np.array([-1,b,b]) * lp/2
			v3 = np.array([-1,-b,c]) * lp/2
		elif self.ibrav ==  12:
			cab = self.celldm[3]
			sab = np.sqrt(1 - cab**2)
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([b*cab,b*sab,0]) * lp
			v3 = np.array([0,0,c]) * lp
		elif self.ibrav == -12:
			cac = self.celldm[4]
			sac = np.sqrt(1 - cac**2)
			v1 = np.array([1,0,0]) * lp
			v2 = np.array([0,b,0]) * lp
			v3 = np.array([c*cac,0,c*sac]) * lp
		elif self.ibrav ==  13:
			cg = self.celldm[3]
			sg = np.sqrt(1 - cg**2)
			v1 = np.array([1,0,-c]) * lp/2
			v2 = np.array([b*cg,b*sg,0]) * lp/2
			v3 = np.array([1,0,c]) * lp
		elif self.ibrav ==  -13:
			cb = self.celldm[4]
			sb = np.sqrt(1 - cb**2)
			v1 = np.array([1,-b,0]) * lp/2
			v2 = np.array([1,b,0]) * lp/2
			v3 = np.array([c*cb,0,c*sb]) * lp
		elif self.ibrav ==  14:
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

		self._cell_p = 'bohr'
		return np.array([v1,v2,v3])



