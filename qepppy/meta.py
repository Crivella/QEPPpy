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

# def key_getter(key, func=None, excp_val=None, excp_func=None):
def key_getter(key):
	def getter(self):
		val = getattr(self, key)
		# if val:
		# 	if not func is None:
		# 		val = func(val)
		# else:
		# 	# print(excp_val)
		# 	if isinstance(excp_val, str):
		# 		val = getattr(self, excp_val)
		# 	# print(excp_func)
		# 	if not excp_func is None :
		# 		val = excp_func(val)
		return val

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


def key_setter(key, 
	typ=None, sub_typ=None, 
	size=None, usize=None, 
	conv_func=lambda x:x,
	set_other_name=None, set_other_func=lambda x: x
	):
	def setter(self, value, size=size):
		# print(f'Setting attribute {key} with value {value}')
		err_header = f"While assigning '{key[1:]}':\n" + "*"*8
		if not check_type(typ, value):
			raise TypeError(
				err_header + 
				f" value='{value}' must be of type '{typ}'."
				)

		if not check_sub_type(sub_typ, value):
			raise TypeError(
				err_header + 
				f"elements of value='{value}' must be of type '{sub_typ}'."
				)

		if size:
			l = len(flatten_iter(value))
			if isinstance(size, int):
				s = size
			if isinstance(size, str):
				app = size.replace(' ', '')
				name,mul = (app.split('*') + [1,])[:2]
				name_v = convert_var(self, name)
				mul_v  = convert_var(self, int(mul))
				s = name_v * mul_v
				if l != s and usize and l % mul_v == 0:
					s = l
					setattr(self, name, l // mul_v)
			if l != s:
				raise TypeError(err_header + f"value='{value}' must be of size '{name}'={s}.")
		try:
			value = conv_func(value)
		except Exception as e:
			print(
				err_header + 
				f"Failed to run conversion function '{conv_func}' on value='{value}'"
				)
			raise e
		setattr(self, key, conv_func(value))
		if not set_other_name is None:
			setattr(self, set_other_name, set_other_func(self, value))

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
	"""
	Metaclass to creaty property object following rules specified in a 
	dictionary.
	The name of the dict being transformed into a property object must be a
	string and not begin with '_'.
	Rules:
	 - typ: 
	    Instance or tuple of instances to be used for check with 
	    'isinstance' for the property (used when calling the setter).
	    If the property is not of any of the specified type, a 'TypeError' will
	    be raised during assignment (cls.prop = var).
	    Examples: 'typ':int,    'typ':(int,float,),
	 - sub_type: 
	    Instance or tuple of instances to be used for check with 
	    'isinstance' for all sub-elements of an iterable (used when calling the 
	    setter).	    
	    If the property is not of any of the specified type, a 'TypeError' will
	    be raised during assignment (cls.prop = var).
	    Examples: 'typ':int,    'typ':(int,float,),
	 - conv_func:
	    Function to be applied to the input data before it is set.
	 - set_other_name:
	    When setting this property, will also set another property using the
	    same value. (Look 'set_other_func' for more details).
	 - set_other_func:
	    Before assgning the value to 'set_other_name', pass the value through
	    this functiond and assign the return value instead.
	    The funcion must accept 2 positional parameters:
	     - cls: instance of the method containing the property.
	     - value: value over which to run the function.
	 - size:
	    Size specification for an iterable object.
	    It can be an 'int' or a 'str' pointing to the name of another method of
	    the object.
	    The 'str' can be 'method_name' or 'method_name * multiplier', where
	    multiplier can either be a number or another method_name.
	    If the lenght of the flattened array is not equal to 'size'  and 'usize'
	    is set to 'False', a 'TypeError' will be raised.
	    Examples: 'size':9,  'size':'n_atoms * 3', 'size':'n_atoms * n_types'.
	 - usize:
	    For size specified using a method name, setting 'usize':True, allows for
	    the change in value of 'method_name' if 'size % multiplier == 0'.
	    For example it allows to set the atoms coordinate without having to
	    manually set the number of atoms.
	 - default:
	    Set a default value for the property.
	    If no value is given the default will be set as typ() where type is the
	    property type (if only one is given), or the first type specified in the
	    tuple.
	 - doc:
	 	String to be used as a docstring for the property.
	"""
	def __new__(cls, clsname, base, dct):
		new_dct = {}
		for k,v in dct.items():
			new_dct[k] = v
			if not isinstance(k, str) or k.startswith("_") or not isinstance(v, dict):
				continue

			hk        = "_" + k
			typ       = v.get('typ',       None)
			sub_typ   = v.get('sub_typ',   None)
			size      = v.get('size',      0 )
			usize     = v.get('usize',     False)
			default   = v.get('default',   None)
			# func      = v.get('func',      None)
			# excp_val  = v.get('excp_val',  None)
			# excp_func = v.get('excp_func', None)
			conv_func = v.get('conv_func', lambda x: x)
			set_other_name = v.get('set_other_name', None)
			set_other_func = v.get('set_other_func', lambda x: x)
			doc       = v.get('doc',       '')

			doc = ('type          = {}.\n'
				   'sub_type      = {}.\n'
				   'default_value = {}.\n\n'
				   '{}').format(typ, sub_typ, default, doc.lstrip().replace('\t', '')) 

			new_dct[hk] = set_typ(typ, default)

			# getter = key_getter(hk, func, excp_val, excp_func)
			getter = key_getter(hk)
			setter = key_setter(hk, typ, sub_typ, size, usize, conv_func,
				set_other_name, set_other_func)

			getter.__doc__ = doc


			new_dct[k] = property(getter, setter)

		return super().__new__(cls, clsname, base, new_dct)



