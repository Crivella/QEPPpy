from ...parsers     import fortran_namelist_collection as fnc
from ...errors      import ParseInputError, ValidateError

tf90_to_py = {
	'INTEGER': int,
	'REAL': float,
	'LOGICAL': bool,
	'CHARACTER': str
	}

tf90_to_np = {
	'INTEGER': 'i4',
	'REAL': 'f8',
	'LOGICAL': 'bool',
	'CHARACTER': 'U64'
	}

class qe_namelists(fnc):
	def __init__(self, tpl=None, **kwargs):
		if tpl:
			self.load_templ(tpl)
		super().__init__(**kwargs)

	@property
	def templ(self):
		return self._templ

	def load_templ(self, tpl=""):
		"""
		Load a QE template from a specified file (tpl) or the internal data
		"""
		import os
		if os.path.isfile(tpl):
			with open(tpl) as f:
				file = f.read()
		else:
			from pkg_resources import resource_string, resource_listdir
			if tpl in resource_listdir('qepppy.qe.parser.data', ''):
				file = resource_string('qepppy.qe.parser.data', tpl).decode('utf-8')

		import ast
		self._templ = ast.literal_eval(file)

	# def validate(self):
	# 	for name,nl in self.items():
	# 		name = name.upper()
	# 		if not name in self.templ:
	# 			raise ValidateError("Invalid namelist '{}'.".format(name))
	# 		for k,v in nl.items():
	# 			if not k in self.templ[name]:
	# 				raise ValidateError("\n\tIn namelist '{}' invalid param '{}'.".format(name,k))
	# 			self._validate_tmpl_item(name, k, v)

	# 	self._validate_default_()

	# def _validate_default_(self):
	# 	for name in self.templ['nl']:
	# 		for elem,v in self.templ[name].items():
	# 			if v['v'] == '***':
	# 				raise ValidateError("\n\tMandatory param '{}/{}' is not set.".format(name, elem))

	# def _validate_vec_(self, value, typ, lim):
	# 	inf = lim[0]
	# 	sup = lim[1]
	# 	if isinstance(sup, str):
	# 		sup = self.deep_find(sup)

	# 	if isinstance(value, (int,float)):
	# 		value = [value,]
	# 	if len(value) <= sup-inf+1:
	# 		return [tf90_to_py[typ](a) if not a is None else None for a in value]
	# 	else:
	# 		raise ValidateError("Too many elements.")

	# def _validate_tmpl_item(self, namelist, param, value):
	# 	if namelist is None:
	# 		for n in self.templ['namelist']:
	# 			if param in self.templ[n]:
	# 				namelist = n
	# 				break

	# 	v = self.templ[namelist][param]

	# 	typ    = v['t']
	# 	possib = v['c']
	# 	vec    = v['vec']

	# 	if not vec is None:
	# 		try:
	# 			v['v'] = self._validate_vec_(value, typ, vec)
	# 		except Exception as e:
	# 			raise ValidateError("\n\t{}/{}:  invalid vec value '{}'. ".format(namelist, param, value) + str(e))
	# 	else:
	# 		value = tf90_to_py[typ](value)
	# 		if possib and all(not value in a for a in possib) and value != '':
	# 			raise ValidateError("\n\t{}/{}: '{}' not among possibilities {}".format(namelist, param, value, possib))
	# 		v['v'] = value
	