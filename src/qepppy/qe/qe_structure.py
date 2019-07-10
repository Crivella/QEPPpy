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

data ={
	'_n_atoms':{
		'xml_search_string':'output//atomic_structure',
		'rstring':r'number of atoms/cell\s*=\s*',
		'mode':'attr=nat',
		'typ':int
		},
	'_n_types':{
		'xml_search_string':'output//atomic_species', 
		'rstring':r'number of atomic types\s*=\s*',
		'mode':'attr=ntyp',
		'typ':int
		},
	'ibrav':{
		'xml_search_string':'output//atomic_structure', 
		'rstring':r'bravais-lattice index\s*=',
		'mode':'attr=bravais_index',
		'typ':int
		},
	'alat':{
		'xml_search_string':'output//atomic_structure', 
		'rstring':r'lattice parameter \(alat\)\s*=',
		'mode':'attr=alat',
		'typ':float,
		're_scale_fact':'_cell_p'
		},
	'_app_cell_p':{
		'rstring':r'cart\. coord\. in units of (?P<flag>.*)\)',
		'typ':str
		},
	'direct':{
		'xml_search_string':'output//cell',
		'rstring':
			r'\s*a\(1\) = \((?P<a1>[\s\d.\-]*)\)\s*\n' + 
			r'\s*a\(2\) = \((?P<a2>[\s\d.\-]*)\)\s*\n' + 
			r'\s*a\(3\) = \((?P<a3>[\s\d.\-]*)\)\s*\n',
		'typ':np.ndarray,
		're_scale_fact':'_cell_p'
		},
	'_app_atom_p':{
		'rstring':r'positions \((?P<flag>.*) units\)',
		'typ':str
		},
	'atoms_typ,atoms_coord_cart':{
		'xml_search_string':'output//atomic_positions/atom', 
		'mode':'attr=name,value',
		'typ':np.ndarray,
		'rstring':r'\d[\t ]+(?P<name>[\w]+).*\(([ \d]+)\) = \((?P<coord>[ \d\.\-]+)\)',
		'max_num':'_n_atoms',
		're_scale_fact':'_atom_p'
		},
	'unique_atoms_typ,_,unique_atoms_mass,unique_atoms_pseudo':{
		'rstring':r'\s*(?P<name>\w+)\s+(?P<valence>[\d\.]+)\s+(?P<mass>[\d\.]+)\s+(?P<pseudo_file>\w+\s*\([ \d\.]+\))',
		'typ':np.ndarray
		},
	'unique_atoms_mass':{
		'xml_search_string':'input//atomic_species//mass',
		'typ':np.ndarray
		},
	'unique_atoms_pseudo':{
		'xml_search_string':'input//atomic_species//pseudo_file',
		'typ':np.ndarray
		},
	'symm_name,symm_matrix':{
		'rstring':
			r'isym =\s*\d{1,2}\s*(?P<name>[\S ]*)\n\s*' +
			r'cryst.\s*s\([\s\d]{2}\) = ' +
			r'(?P<rotation>(\(.*\)\s*){3})',
		'typ':np.ndarray
		},
	'symm_matrix':{
		'xml_search_string':'output//symmetry/rotation',
		'typ':np.ndarray
		},
	'symm_name':{
		'xml_search_string':'output//symmetry/info', 
		'mode':'attr=name'
		},
	'fft_dense_grid':{
		'xml_search_string':'output//fft_grid',
		'rstring':
			r'Dense.*FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*\)',
		'mode':r'attr=nr\d',
		'typ':np.ndarray
		},
	'fft_smooth_grid':{
		'xml_search_string':'output//fft_smooth', 
		'rstring':
			r'Smooth.*FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*\)',
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

		xml_data.update(data)
		regex_data.update(data)
		super().__init__(xml_data=xml_data, regex_data=regex_data, **kwargs)

	def _format_cell_(self, info):
		if self.ibrav and info == 0:
			return ""
		msg = "CELL_PARAMETERS\n"
		for l in self.direct:
			for e in l:
				msg += " {:12.7f}".format(e)
			msg += "\n"
		msg += "\n"
		return msg

	def _format_atom_spec_(self):
		msg = "ATOMIC_SPECIES\n"
		for t,m,p in zip(self.unique_atoms_typ, self.unique_atoms_mass, self.unique_atoms_pseudo):
			msg += "  {:6}{:12.4f}  {}".format(t, m, p)
			msg += "\n"
		msg += "\n\n"
		return msg

	def _format_atom_pos_(self):
		msg = "ATOMIC_POSITIONS {crystal}\n"
		for coord,name in zip(self.atoms_coord_cryst, self.atoms_typ):
			msg += "  {}  ".format(name)
			for c in coord:
				msg += "{:14.9f}".format(c)
			msg += "\n"
		msg += "\n"
		return msg

	def __str__(self, info=0):
		# msg = super().__str__()
		msg  = ''
		msg += self._format_cell_(info)
		msg += self._format_atom_spec_()
		msg += self._format_atom_pos_()

		return msg

	@property
	def _atom_p(self):
		"""Conversion factor for atom coordinates to atomic units"""
		if self._app_atom_p == 'alat':
			return self.alat
		if self._app_atom_p == 'angstrom':
			return 1./0.529177
		return 1.

	@property
	def _cell_p(self):
		"""Conversion factor for cell vetors to atomic units"""
		if self._app_cell_p == 'alat':
			return self.alat
		if self._app_cell_p == 'angstrom':
			return 1./0.529177
		return 1.

	# @_cell_p.setter
	# def _cell_p(self, value):
	# 	self._app_cell_p = value
	
	# @property
	# def symm_matrix(self):
	# 	"""List of symmetry operation matrices"""
	# 	t = type(self._symm[0]['rotation'])
	# 	if t == np.ndarray:
	# 		res = np.array([a['rotation'].reshape(3,3) for a in self._symm])
	# 	elif t == str:
	# 		f = lambda x: ' '.join([a.split('(')[1] for a in x.split('\n')])
	# 		g = lambda x: f(x).replace('(',' ').replace(')',' ').replace('\n',' ').replace('f =', ' ')
	# 		res = np.array([np.fromstring(g(a['rotation']), sep=' ').reshape(3,3) for a in self._symm])
	# 	else:
	# 		raise NotImplementedError()
	# 	return res

	# @property
	# def symm_name(self):
	# 	"""List of symmetry operation names"""
	# 	return list([a['name'] for a in self._symm])

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

		celldm = self['celldm']
		if len(celldm) > 1:
			b = celldm[1]
		if len(celldm) > 2:
			c = celldm[2]
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

		# self._app_cell_p = 'bohr'
		return np.array([v1,v2,v3])



