from ...parsers     import fortran_namelist_collection as fnc
from ...errors      import ParseInputError, ValidateError

tf90_to_py = {
	'INTEGER': int,
	'REAL': float,
	'LOGICAL': bool,
	'CHARACTER': str
	}

# tf90_to_np = {
# 	'INTEGER': 'i4',
# 	'REAL': 'f8',
# 	'LOGICAL': 'bool',
# 	'CHARACTER': 'U64'
# 	}

def key_getter(key, item=None, attr=None):
	def getter(cls):
		if not item is None:
			return cls[item]
		elif not attr is None:
			return getattr(cls, attr)
		else:
			raise ValueError(f"Cannot retrieve {key}")

	return getter

def key_setter(key, item=None, attr=None):
	def setter(cls, key, value):
		if not item is None:
			cls[item] = value
		elif not attr is None:
			setattr(cls, attr, value)
		else:
			raise ValueError(f"Cannot set {key}")

	return setter

class VariableLinker(type):
	def __new__(cls, clsname, base, dct):
		new_dct = {}
		for k,v in dct.items():
			new_dct[k] = v
			if not isinstance(k, str) or not k.startswith("_link__") or not isinstance(v, dict):
				continue

			k       = k[7:]
			hk      = "_" + k
			attrib  = v.pop('attrib',    None)
			item    = v.pop('item',      None)

			# if not (attrib or item):
			# 	raise ValueError("Must link the variable to either an attibute or an item.")
			if attrib and item:
				raise ValueError("Must link the variable to only an attibute or an item, not both.")

			getter = key_getter(hk, item, attrib)
			setter = key_setter(hk, item, attrib)


			new_dct[k] = property(getter, setter)

		res = super().__new__(cls, clsname, base, new_dct)

		return res

class qe_namelists(fnc, metaclass=VariableLinker):
	_link__n_atoms={'item':'SYSTEM/nat'}

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

from ...calc_system import system
from ..qe_structure import qe_structure as structure
from ...meta.property_creator import PropertyCreator
class meta_app(PropertyCreator, VariableLinker):
	pass
class pw_in(qe_namelists, structure, system, metaclass=meta_app):
	pass