import numpy as np
from ..parsers import Parser_xml, Parser_regex
# from .parser.data_file_parser import data_file_parser as dfp
from ..errors import ValidateError
# from ..calc_system.structure import structure as structure
# from ..calc_system import system
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

# _data={
# 	'_n_atoms':{
# 		'xml_ptype':'attr', 
# 		'xml_search_string':'output//atomic_structure', 
# 		'extra_name':'nat', 
# 		'res_type':int,
# 		'rstring':r'number of atoms/cell\s*=\s*'
# 		},
# 	'_n_types':{
# 		'xml_ptype':'attr', 
# 		'xml_search_string':'output//atomic_species', 
# 		'extra_name':'ntyp', 
# 		'res_type':int,
# 		'rstring':r'number of atomic types\s*=\s*'
# 		},
# 	'ibrav':{
# 		'xml_ptype':'attr', 
# 		'xml_search_string':'output//atomic_structure', 
# 		'extra_name':'bravais_index', 
# 		'res_type':int,
# 		'rstring':r'bravais-lattice index\s*='
# 		},
# 	'alat':{
# 		'xml_ptype':'attr', 
# 		'xml_search_string':'output//atomic_structure', 
# 		'extra_name':'alat', 
# 		'res_type':float,
# 		'rstring':r'lattice parameter \(alat\)\s*='
# 		},
# 	'_app_cell_p':{
# 		'res_type':str,
# 		'rstring':r'cart\. coord\. in units of (?P<flag>.*)\)'
# 		},
# 	'_cell':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//cell', 
# 		'extra_name':None, 
# 		'res_type':list,
# 		'rstring':
# 			r'\s*a\(1\) = \((?P<a1>[\s\d.\-]*)\)\s*\n' + 
# 			r'\s*a\(2\) = \((?P<a2>[\s\d.\-]*)\)\s*\n' + 
# 			r'\s*a\(3\) = \((?P<a3>[\s\d.\-]*)\)\s*\n'
# 		},
# 	'_recip':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//reciprocal_lattice', 
# 		'extra_name':None, 
# 		'res_type':list,
# 		'rstring':
# 			r'\s*b\(1\) = \((?P<b1>[\s\d.\-]*)\)\s*\n' + 
# 			r'\s*b\(2\) = \((?P<b2>[\s\d.\-]*)\)\s*\n' + 
# 			r'\s*b\(3\) = \((?P<b3>[\s\d.\-]*)\)\s*\n'
# 		},
# 	'_app_atom_p':{
# 		'res_type':str,
# 		'rstring':r'positions \((?P<flag>.*) units\)'
# 		},
# 	'_atoms':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'input//atom', 
# 		'extra_name':'coord', 
# 		'res_type':list,
# 		'rstring':r'\d[\t ]+(?P<name>[\w]+).*\((?P<index>[ \d]+)\) = \((?P<coord>[ \d\.\-]+)\)'
# 		},
# 	'_atom_spec':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'input//species', 
# 		'extra_name':None, 
# 		'res_type':list,
# 		'rstring':r'\s*(?P<name>\w+)\s+(?P<valence>[\d\.]+)\s+(?P<mass>[\d\.]+)\s+(?P<pseudo_file>\w+\s*\([ \d\.]+\))'
# 		# 'rstring':
# 		# 	r'PseudoPot. \#.*\s+(.*/)*(?P<pseudo_file>.+(\.UPF|\.upf))' +
# 		# 	r'(.*\n)+\s*atomic species.*' +
# 		# 	r'\s*(?P<name>\w+)\s+(?P<valence>[\d\.]+)\s+(?P<mass>[\d\.]+)'
# 		},
# 	'_symm':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//symmetry', 
# 		'extra_name':None, 
# 		'res_type':list,
# 		'rstring':
# 			r'isym =\s*\d{1,2}\s*(?P<name>[\S ]*)\n\s*' +
# 			r'cryst.\s*s\([\s\d]{2}\) = ' +
# 			r'(?P<rotation>(\(.*\)\s*){3})'
# 		},
# 	'_fft_dense_grid':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//fft_grid', 
# 		'extra_name':None, 
# 		'res_type':list,
# 		'rstring':
# 			r'Dense.*FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*\)'
# 		},
# 	'_fft_smooth_grid':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//fft_smooth', 
# 		'extra_name':None, 
# 		'res_type':list,
# 		'rstring':
# 			r'Smooth.*FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*\)'
# 		}
# 	}

data_regex={
	'_n_atoms':{
		'rstring':r'number of atoms/cell\s*=\s*',
		'typ':int
		},
	'_n_types':{
		'rstring':r'number of atomic types\s*=\s*',
		'typ':int
		},
	'ibrav':{
		'rstring':r'bravais-lattice index\s*=',
		'typ':int
		},
	'alat':{
		'rstring':r'lattice parameter \(alat\)\s*=',
		'typ':float
		},
	'_app_cell_p':{
		'rstring':r'cart\. coord\. in units of (?P<flag>.*)\)',
		'typ':str
		},
	'direct':{
		'rstring':
			r'\s*a\(1\) = \((?P<a1>[\s\d.\-]*)\)\s*\n' + 
			r'\s*a\(2\) = \((?P<a2>[\s\d.\-]*)\)\s*\n' + 
			r'\s*a\(3\) = \((?P<a3>[\s\d.\-]*)\)\s*\n',
		'typ':np.ndarray,
		'mode':'get_all',
		'scale_fact':'_cell_p'
		},
	'_app_atom_p':{
		'rstring':r'positions \((?P<flag>.*) units\)',
		'typ':str
		},
	'atoms_typ,_,atoms_coord_cart':{
		'rstring':r'\d[\t ]+(?P<name>[\w]+).*\((?P<index>[ \d]+)\) = \((?P<coord>[ \d\.\-]+)\)',
		'typ':np.ndarray,
		'max_num':'_n_atoms',
		'scale_fact':'_atom_p'
		},
	'unique_atoms_typ,_,unique_atoms_mass,unique_atoms_pseudo':{
		'rstring':r'\s*(?P<name>\w+)\s+(?P<valence>[\d\.]+)\s+(?P<mass>[\d\.]+)\s+(?P<pseudo_file>\w+\s*\([ \d\.]+\))',
		'typ':np.ndarray
		},
	'_symm':{
		'rstring':
			r'isym =\s*\d{1,2}\s*(?P<name>[\S ]*)\n\s*' +
			r'cryst.\s*s\([\s\d]{2}\) = ' +
			r'(?P<rotation>(\(.*\)\s*){3})',
		'typ':np.ndarray
		},
	'fft_dense_grid':{
		'rstring':
			r'Dense.*FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*\)',
		'typ':np.ndarray,
		'mode':'get_all'
		},
	'fft_smooth_grid':{
		'rstring':
			r'Smooth.*FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*\)',
		'typ':np.ndarray,
		'mode':'get_all'
		}
	}

data_xml={
	'_n_atoms':{
		'xml_search_string':'output//atomic_structure',
		'mode':'attr=nat',
		'typ':int
		},
	'_n_types':{
		'xml_search_string':'output//atomic_species', 
		'mode':'attr=ntyp',
		'typ':int
		},
	'ibrav':{
		'xml_search_string':'output//atomic_structure', 
		'mode':'attr=bravais_index',
		'typ':int
		},
	'alat':{
		'xml_search_string':'output//atomic_structure', 
		'mode':'attr=alat',
		'typ':float
		},
	# '_app_cell_p':{
	# 	'rstring':r'cart\. coord\. in units of (?P<flag>.*)\)'
	# 	},
	'direct':{
		'xml_search_string':'output//cell',
		'typ':np.ndarray
		},
	# '_app_atom_p':{
	# 	'res_type':str,
	# 	'rstring':r'positions \((?P<flag>.*) units\)'
	# 	},
	'atoms_coord_cart':{
		'xml_search_string':'input//atomic_positions/atom', 
		'typ':np.ndarray
		},
	'atoms_typ':{
		'xml_search_string':'input//atomic_positions/atom',
		'mode':'attr=name',
		'typ':np.ndarray
		},
	# '_atom_spec':{
	# 	'xml_search_string':'input//species', 
	# 	},
	'unique_atoms_mass':{
		'xml_search_string':'input//atomic_species//mass',
		'typ':np.ndarray
		},
	'unique_atoms_pseudo':{
		'xml_search_string':'input//atomic_species//pseudo_file',
		'typ':np.ndarray
		},
	# '_symm':{
	# 	'xml_search_string':'output//symmetry', 
	# 	},
	'fft_dense_grid':{
		'xml_search_string':'output//fft_grid',
		'mode':r'attr=nr\d',
		'typ':np.ndarray
		},
	'fft_smooth_grid':{
		'xml_search_string':'output//fft_smooth', 
		'mode':r'attr=nr\d',
		'typ':np.ndarray
		}
	}

class qe_structure(Parser_xml, Parser_regex):
	__name__ = "qe_structure";
	def __init__(self, xml_data={}, regex_data={}, **kwargs):
		if not hasattr(self,'_app_atom_p') or not self._app_atom_p:
			self._app_atom_p = 'bohr'
		if not hasattr(self,'_app_cell_p') or  not self._app_cell_p:
			self._app_cell_p = 'bohr'
		if self._direct is None or self._direct == []:
			self._direct = np.diag([1]*3)

		xml_data.update(data_xml)
		regex_data.update(data_regex)
		super().__init__(xml_data=xml_data, regex_data=regex_data, **kwargs)

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
	def _atom_p(self):
		"""Conversion factor for atom coordinates to atomic units"""
		cmp = self._app_atom_p #.strip()
		if cmp == 'alat':
			return self.alat
		if cmp == 'angstrom':
			return 1./0.529177
		return 1.

	@property
	def _cell_p(self):
		"""Conversion factor for cell vetors to atomic units"""
		cmp =  self._app_cell_p
		if cmp == 'alat':
			return self.alat
		if cmp == 'angstrom':
			return 1./0.529177
		return 1.

	@_cell_p.setter
	def _cell_p(self, value):
		self._app_cell_p = value
	
	@property
	def symm_matrix(self):
		"""List of symmetry operation matrices"""
		t = type(self._symm[0]['rotation'])
		if t == np.ndarray:
			res = np.array([a['rotation'].reshape(3,3) for a in self._symm])
		elif t == str:
			f = lambda x: ' '.join([a.split('(')[1] for a in x.split('\n')])
			g = lambda x: f(x).replace('(',' ').replace(')',' ').replace('\n',' ').replace('f =', ' ')
			res = np.array([np.fromstring(g(a['rotation']), sep=' ').reshape(3,3) for a in self._symm])
		else:
			raise NotImplementedError()
		return res

	@property
	def symm_name(self):
		"""List of symmetry operation names"""
		return list([a['name'] for a in self._symm])

	def validate(self):
		if self.ibrav == None:
			raise ValidateError("ibrav is not set.")
		# if self._atom_spec == None:
		if self._atoms_typ is None:
			raise ValidateError("List of atom types is not set.")
		# if self._atoms == None:
		if self._atoms_coord_cart is None:
			raise ValidateError("List of atomic positions is not set.")

		for typ in self.atoms_typ:
			if not typ in self.unique_atoms_typ:
				raise ValidateError("Atoms in ATOMIC_POSITION do not match the type in ATOMIC_SPECIES")

		if self.ibrav == 0:
			if self.direct is None:
				raise ValidateError("Cell structure is not set with 'ibrav = 0'.")

		try:
			super().validate()
		except AttributeError:
			pass


	def _ibrav_to_cell_(self):
		if self.ibrav == None:
			raise ValueError("Failed to generate cell structure from self.ibrav: self.ibrav not set.")

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

		self._app_cell_p = 'bohr'
		return np.array([v1,v2,v3])



