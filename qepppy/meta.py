def flatten_iter(iterable):
	try:
		iter(iterable)
	except:
		return [iterable,]

	if isinstance(iterable, str):
		return [iterable,]

	res = []
	for i in iterable:
		res += flatten_iter(i)
	return res

def key_getter(key):
	def getter(self):
		return getattr(self, key)

	return getter

def check_type(typ, value):
	if typ is None or any(isinstance(value,t) for t in flatten_iter(typ)):
		return True


def check_sub_type(sub_typ, value):
	if sub_typ is None or all(any(isinstance(v,t) for t in sub_typ) for v in flatten_iter(value)):
		return True

def convert_var(self, value, typ=int):
	if not str in flatten_iter(typ) and isinstance(value,str):
		return getattr(self, value)
	if isinstance(value, typ):
		return value


def key_setter(key, typ=None, sub_typ=None, size=None):
	def setter(self, value):
		if not check_type(typ, value):
			raise TypeError("'{}' must be of type '{}'.".format(key[1:],typ))

		if not check_sub_type(sub_typ, value):
			raise TypeError("Elements of '{}' must be of type '{}'.".format(value,sub_typ))

		if not size is None:
			l = len(flatten_iter(value))
			if isinstance(size, int):
				s = size
			if isinstance(size, str):
				app = size.replace(' ', '')
				name,mul = (app.split('*') + [1,])[:2]
				name = convert_var(self, name)
				mul  = convert_var(self, int(mul))
				s = name * mul
			if l != s:
				raise TypeError("'{}' must be of size {}.".format(value, s))
		setattr(self, key, value)

	return setter

def set_typ(typ, value):
	if typ is None:
		return value

	if value is None:
		return flatten_iter(typ)[0]()

	if isinstance(typ, tuple):
		for t in typ:
			if isinstance(value, t):
				return value
		raise TypeError("Invalid type value '{}' for type '{}'.".format(type(value), typ))

	return typ(value)

class PropertyCreator(type):
	def __new__(cls, clsname, base, dct):
		new_dct = {}
		for k,v in dct.items():
			new_dct[k] = v
			if not isinstance(k, str) or k.startswith("_") or not isinstance(v, dict):
				continue

			hk      = "_" + k
			typ     = v.get('typ', None)
			sub_typ = v.get('sub_typ', None)
			size    = v.get('size', None)
			default = v.get('default', None)
			doc     = v.get('doc','')

			doc = ('type          = {}.\n'
				   'sub_type      = {}.\n'
				   'default_value = {}.\n\n'
				   '{}').format(typ, sub_typ, default, doc.lstrip().replace('\t', '')) 

			new_dct[hk] = set_typ(typ, default)

			getter = key_getter(hk)
			setter = key_setter(hk, typ, sub_typ, size)

			getter.__doc__ = doc


			new_dct[k] = property(getter, setter)

		return super().__new__(cls, clsname, base, new_dct)



