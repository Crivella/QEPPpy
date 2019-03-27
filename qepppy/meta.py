def key_getter(key):
	def getter(self):
		return getattr(self, key)
		
	return getter

def key_setter(key, typ=None):
	def setter(self, value):
		if not typ is None:
			if isinstance(typ, tuple):
				for t in typ:
					if isinstance(value, typ):
						return setattr(self, key, value)
				raise TypeError("'{}' must be of type '{}'.".format(key[1:],typ))
			elif not isinstance(value, typ):
				raise TypeError("'{}' must be of type '{}'.".format(key[1:],typ))
		setattr(self, key, value)

	return setter

def set_typ(typ, value):
	if typ is None:
		return value

	if value is None:
		return typ()

	if isinstance(typ, tuple):
		for t in typ:
			if isinstance(value, t):
				return value
		raise TypeError("Invalid type value '{}' for type '{}'.".format(type(value), typ))

	return typ(value)


class PropertyCreator(type):
	def __new__(cls, clsname, base, dct):
		# print("PropertyCreator: __new__", cls, clsname, base, dct)
		new_dct = {}
		for k,v in dct.items():
			new_dct[k] = v
			if not isinstance(k, str) or k.startswith("_") or not isinstance(v, dict):
				continue

			hk      = "_" + k
			typ     = v.get('typ',None)
			default = v.get('default',None)
			doc     = v.get('doc','')

			doc = ('type          = {}.\n'
				   'default_value = {}.\n\n'
				   '{}').format(typ, default, doc.lstrip().replace('\t', '')) 

			new_dct[hk] = set_typ(typ, default)

			getter = key_getter(hk)
			setter = key_setter(hk, typ)

			getter.__doc__ = doc


			new_dct[k] = property(getter, setter)

		return super().__new__(cls, clsname, base, new_dct)