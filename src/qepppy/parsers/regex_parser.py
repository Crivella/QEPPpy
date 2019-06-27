import re
import numpy as np

matches = {
	int:   r'[\s]*(?P<flag>[\d\-]+)',
	float: r'[\s]*(?P<flag>[\d\.\-EedD]+)',
	str:   r'[\s]*(?P<flag>.*)',
	bool:  r'[\s]*(?P<flag>.)',
}

def dbg1(*args, **kwargs):
	return
	print(*args, **kwargs)

class Parser_regex():
	def __init__(self, *args, file=None, regex_data={}, **kwargs):
		super().__init__(*args, **kwargs)

		self.regex_data = regex_data

		if file:
			self.parse_file(file)

	@property
	def regex_data(self):
		return self._regex_data

	@regex_data.setter
	def regex_data(self, value):
		if not isinstance(value, dict):
			raise ValueError("Value must be a dictionary.")

		self._regex_data = value

	@staticmethod
	def get_single_val(content, rstring, typ):
		res = None

		if not 'flag' in rstring:
			rstring += matches[typ]
		val = re.search(rstring, content)
		if not val is None:
			res = typ(val.group('flag'))

		return res

	@staticmethod
	def get_list_val(content, rstring, typ):
		res = re.finditer(rstring, content)
		res = [x.groupdict() for x in res]
		for n,e in enumerate(res):
			for k,v in e.items():
				b = np.fromstring(v, sep=' ')
				if len(b) == 0 or re.findall(r'\s[a-zA-Z]', v):
					b = str(v).strip()
				elif len(b) == 1:
					b = b[0]
				res[n][k] = b

		return res

	# @staticmethod
	def max_num_truncate(self, max_num, val):
		if not max_num is None:
			if isinstance(max_num, str):
				max_num = getattr(self, max_num)
			val = val[:max_num]

		return val

	def scale(self, val, fact):
		if not fact is None:
			if isinstance(fact, str):
				fact = getattr(self, fact)
			if len(val)>0 and isinstance(val.flatten()[0], (int,float,np.number)):
				val *= fact

		return val

	

	def parse_file(self, file):
		# res = {}

		with open(file, 'r') as f:
			content = f.read()

		for k,v in self.regex_data.items():
			if not 'rstring' in v:
				continue
			dbg1('-'*40)
			dbg1('Key:', k)
			dbg1('val:', v)

			rstring    = v.get('rstring')
			typ        = v.get('typ')
			max_num    = v.get('max_num', None)
			scale_fact = v.get('re_scale_fact', None)

			if typ in matches:
				val = self.get_single_val(content, rstring, typ)
				setattr(self, k, val)
				dbg1('found_single:', val)
			elif typ in (list,np.ndarray):
				app = self.get_list_val(content, rstring, typ)
				params = k.split(',')
				dbg1(params)
				dbg1(app)
				for num,name in enumerate(params):
					if name == '_':
						continue
					# n1  = num
					# n2 = slice(None)
					n_elem  = len(app)
					n_group = len(app[0].values()) if n_elem > 0 else 0
					if len(params) == 1 and not (n_elem > 0 and n_group == 1):
						# n1  = slice(None)
						# n2 = slice(0,n_group)
						app2 = np.array([list(a.values()) for a in app])
					else:
						app2 = np.array([list(a.values())[num] for a in app])

					# dbg1('n1,n2:', n1,n2)



					val = self.max_num_truncate(max_num, app2)
					val = self.scale(val, scale_fact)

					dbg1(f'found_mul({name}):', val)
					setattr(self, name, val)



# data_regex={
# 	# '_n_kpt':{
# 	# 	'rstring':r'number of k points[\s]*=',
# 	# 	'typ':int
# 	# 	},
# 	# '_n_bnd':{
# 	# 	'rstring':r'number of Kohn-Sham states[\s]*=',
# 	# 	'typ':int
# 	# 	},
# 	# '_n_el':{
# 	# 	'rstring':r'number of electrons[\s]*=',
# 	# 	'typ':float
# 	# 	},
# 	# 'fermi':{
# 	# 	'rstring':r'the Fermi energy is',
# 	# 	'typ':float
# 	# 	},
# 	# 'kpt,weight':{
# 	# 	'rstring':r'[\s]{4,}k\([ \d]+\) = \((?P<kpt>[ \d\.\-]+)\).*wk = (?P<weight>[ \d\.]+)',
# 	# 	'typ':np.ndarray,
# 	# 	'max_num':'_n_kpt'
# 	# 	},


# 	'_n_atoms':{
# 		'rstring':r'number of atoms/cell\s*=\s*',
# 		'typ':int
# 		},
# 	'_n_types':{
# 		'rstring':r'number of atomic types\s*=\s*',
# 		'typ':int
# 		},
# 	'ibrav':{
# 		'rstring':r'bravais-lattice index\s*=',
# 		'typ':int
# 		},
# 	'alat':{
# 		'rstring':r'lattice parameter \(alat\)\s*=',
# 		'typ':float
# 		},
# 	'_app_cell_p':{
# 		'rstring':r'cart\. coord\. in units of (?P<flag>.*)\)',
# 		'typ':str
# 		},
# 	'direct':{
# 		'rstring':
# 			r'\s*a\(1\) = \((?P<a1>[\s\d.\-]*)\)\s*\n' + 
# 			r'\s*a\(2\) = \((?P<a2>[\s\d.\-]*)\)\s*\n' + 
# 			r'\s*a\(3\) = \((?P<a3>[\s\d.\-]*)\)\s*\n',
# 		'typ':np.ndarray,
# 		'mode':'get_all'
# 		},
# 	'_app_atom_p':{
# 		'rstring':r'positions \((?P<flag>.*) units\)',
# 		'typ':str
# 		},
# 	'atoms_typ,_,atoms_coord_cart':{
# 		'rstring':r'\d[\t ]+(?P<name>[\w]+).*\((?P<index>[ \d]+)\) = \((?P<coord>[ \d\.\-]+)\)',
# 		'typ':np.ndarray,
# 		'max_num':'_n_atoms'
# 		},
# 	'unique_atoms_typ,unique_atoms_mass,unique_atoms_pseudo':{
# 		'rstring':r'\s*(?P<name>\w+)\s+(?P<valence>[\d\.]+)\s+(?P<mass>[\d\.]+)\s+(?P<pseudo_file>\w+\s*\([ \d\.]+\))',
# 		},
# 	'_symm':{
# 		'rstring':
# 			r'isym =\s*\d{1,2}\s*(?P<name>[\S ]*)\n\s*' +
# 			r'cryst.\s*s\([\s\d]{2}\) = ' +
# 			r'(?P<rotation>(\(.*\)\s*){3})',
# 		'typ':int
# 		},
# 	'fft_dense_grid':{
# 		'rstring':
# 			r'Dense.*FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*\)',
# 		'typ':np.ndarray,
# 		'mode':'get_all'
# 		},
# 	'fft_smooth_grid':{
# 		'rstring':
# 			r'Smooth.*FFT dimensions:\s*\(\s*(?P<nr1>\d*),\s*(?P<nr2>\d*),\s*(?P<nr3>\d*)\s*\)',
# 		'typ':np.ndarray,
# 		'mode':'get_all'
# 		}

# }

# a = Parser_regex(regex_data=data_regex)
# a.parse_file('/home/crivella/app/1_bands.out')
# print(b)
# print(b['kpt'].shape)
# print(b['weight'].shape)
# print(b['weight'])
# print(a.kpt.shape)
# print(a.weight.shape)