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
	def getter(cls):
		val = getattr(cls, key)
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

def set_other(cls, other_name, other_func, value):
	if other_name is None:
		return

	if isinstance(other_name, tuple):
		l_name = other_name
		l_func = other_func
	else:
		l_name = (other_name,)
		l_func = (other_func,)

	for name,func in zip(l_name, l_func):
		setattr(cls, name, func(cls, value))


def key_setter(key, 
	typ=None, sub_typ=None, 
	shape=None,
	conv_func=lambda x:x,
	pre_set_name=None, pre_set_func=None,
	post_set_name=None, post_set_func=None,
	):
	def setter(cls, value):
		# print(f'Setting attribute {key} with value {value}')
		err_header = f"While assigning '{key[1:]}':\n" + "*"*8 + " "

		set_other(cls, pre_set_name, pre_set_func, value)

		if not check_type(typ, value):
			raise TypeError(
				err_header + 
				f"value='{value}' must be of type '{typ}'."
				)

		if not check_sub_type(sub_typ, value):
			raise TypeError(
				err_header + 
				f"elements of value='{value}' must be of type '{sub_typ}'."
				)

		if shape:
			if hasattr(value, 'shape'):
				v_shape = value.shape
				if len(v_shape) != len(shape):
					raise ValueError(err_header + f"Mismatch in number of dimension required value='{v_shape}' vs required='{shape}'.")
				for ve, ce in zip(v_shape, shape):
					app = convert_var(cls, ce)
					if ve != app and app != -1:
						raise ValueError(err_header + f"Shape mismatch '{ce}'={app} not equal to '{v_shape}'")
			else:
				l   = len(value)
				if len(shape) > 1:
					raise ValueError(err_header + f"Can't confront shape of value='{l}' with '{shape}'")
				app = convert_var(cls, shape[0])
				if l != app:
					raise ValueError(err_header + f"Shape mismatch '{shape[0]}'={app} not equal to {l}")

		try:
			value = conv_func(value)
		except Exception as e:
			print(
				err_header + 
				f"Failed to run conversion function '{conv_func}' on value='{value}'"
				)
			raise e

		set_other(cls, post_set_name, post_set_func ,value)
		setattr(cls, key, value)

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
	 - settable:
	    True/False
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
	 - pre_set_name:
	    When setting this property, will also set another property using the
	    same value BEFORE all the checks are done. (Look 'pre_set_func' for more details).
	 - pre_set_func:
	    Before assgning the value to 'pre_set_name', pass the value through
	    this functiond and assign the return value instead.
	    The funcion must accept 2 positional parameters:
	     - cls: instance of the method containing the property.
	     - value: value over which to run the function.
	 - post_set_name:
	    When setting this property, will also set another property using the
	    same value AFTER all the checks are done. (Look 'post_set_func' for more details).
	 - post_set_func:
	    Before assgning the value to 'post_set_name', pass the value through
	    this functiond and assign the return value instead.
	    The funcion must accept 2 positional parameters:
	     - cls: instance of the method containing the property.
	     - value: value over which to run the function.
	 - shape:
	    Shape of the object. Must be a tuple.
	    The element of the tuple must be an 'int' or a 'str' pointing to the name 
	    of another attribute of the object.
	    The assigned value must have the 'shape' attribute, unless the tuple only has
	    one element. In this case,  shape[0] will be compared with len(value).
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
			settable  = v.pop('settable',  True)
			typ       = v.pop('typ',       None)
			sub_typ   = v.pop('sub_typ',   None)
			shape     = v.pop('shape',     0)
			default   = v.pop('default',   None)
			conv_func = v.pop('conv_func', lambda x: x)

			pre_set_name = v.pop('pre_set_name', None)
			pre_set_func = v.pop('pre_set_func', lambda x: x)

			post_set_name = v.pop('post_set_name', None)
			post_set_func = v.pop('post_set_func', lambda x: x)

			doc       = v.pop('doc',       '')

			if v:
				raise ValueError(f"Invalid keywords '{list(v.keys())}' for metaclass 'PropertyCreator'.")

			doc = ('type          = {}.\n'
				   'sub_type      = {}.\n'
				   'default_value = {}.\n\n'
				   '{}').format(typ, sub_typ, default, doc.lstrip().replace('\t', '')) 

			new_dct[hk] = set_typ(typ, default)

			# getter = key_getter(hk, func, excp_val, excp_func)
			getter = key_getter(hk)
			if settable:
				setter = key_setter(
					hk, typ, sub_typ, 
					shape,
					conv_func,
					pre_set_name, pre_set_func,
					post_set_name, post_set_func)
			else:
				setter = None

			new_dct[k] = property(getter, setter, doc=doc)

		return super().__new__(cls, clsname, base, new_dct)



