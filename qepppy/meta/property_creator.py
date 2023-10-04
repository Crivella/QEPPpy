import numpy as np

err_header = ''

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
    def getter(cls):
        val = getattr(cls, key)
        return val

    return getter

def check_type(typ, value):
    if typ is None or isinstance(value, typ):
        return True

    raise TypeError(
        err_header + 
        f"value='{value}' must be of type '{typ}'."
        )


def check_sub_type(sub_typ, value):
    if sub_typ is None or all(isinstance(v,sub_typ) for v in flatten_iter(value)):
        return True

    raise TypeError(
        err_header + 
        f"elements of value='{value}' must be of type '{sub_typ}'."
        )

def _check_shape_np(cls, value, shape):
    l = len(shape)
    if hasattr(value, 'shape'):
        v_shape = value.shape
        lv = len(v_shape)
        if lv == 1 and v_shape[0] == 0:
            return
        if lv != l:
            raise ValueError(err_header + f"Mismatch in number of dimensions: value='{v_shape}' vs required='{shape}'.")
        for ve, ce in zip(v_shape, shape):
            try:
                app = convert_var(cls, ce)
            except KeyError:
                raise ValueError(err_header + f"Attribute '{ce}' is not present in '{cls}'.")
            if ve != app and app != -1:
                raise ValueError(err_header + f"Shape mismatch '{ce}'={app} not equal to '{v_shape}'.")

def _check_value_list(cls, value, shape):
    l  = len(shape)
    lv = len(value)
    if lv == 0:
        return
    if l > 1:
        raise ValueError(err_header + f"Can't confront shape of value='{lv}' with '{shape}'")
    try:
        app = convert_var(cls, shape[0])
    except KeyError:
        raise ValueError(err_header + f"Attribute '{shape[0]}' is not present in '{cls}'.")

    if lv != app and app != -1:
        raise ValueError(err_header + f"Shape mismatch '{shape[0]}'={app} not equal to {lv}")

def check_shape(cls, value, shape):
    if not shape:
        return
    if hasattr(value, 'shape'):
        _check_shape_np(cls, value, shape)
    else:
        _check_value_list(cls, value, shape)

def check_allowed(value, allowed, nonetyp=None):
    if len(allowed) == 0:
        return
    if not nonetyp is None and type(value) == type(nonetyp) and value == nonetyp:
        return
    for v in flatten_iter(value):
        if not v in allowed:
            raise ValueError(err_header + f"Value '{v}' is not allowed.")

def convert_var(cls, value, typ=int):
    if not str in flatten_iter(typ) and isinstance(value,str):
        return getattr(cls, value)
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

def set_typ(typ, value):
    if typ is None:
        return value

    if value is None:
        app = flatten_iter(typ)[0]
        if app == np.ndarray:
            return np.empty((0,))
        return app()

    if isinstance(typ, tuple):
        for t in typ:
            if isinstance(value, t):
                return value
        raise TypeError("Invalid type value '{}' for type '{}'.".format(type(value), typ))

    return typ(value)

def key_setter(key, 
    typ=None, sub_typ=None, 
    shape=None,
    conv_func=lambda x:x,
    allowed=[],
    pre_set_name=None, pre_set_func=None,
    post_set_name=None, post_set_func=None,
    ):
    def setter(cls, value):
        # print(f'Setting attribute {key} with value {value}')
        global err_header
        err_header = f"While assigning '{key[1:]}':\n" + "*"*8 + " "

        nonetyp = None
        if value is None and not None in typ:
            value = set_typ(typ, value)
            nonetyp = value 

        set_other(cls, pre_set_name, pre_set_func, value)

        check_type(typ, value)
        check_sub_type(sub_typ, value)

        try:
            value = conv_func(value)
        except Exception as e:
            raise Exception(
                err_header + 
                f"Failed to run conversion function '{conv_func}' on value='{value}'"
                ) from e

        check_allowed(value, allowed, nonetyp)
        check_shape(cls, value, shape)

        set_other(cls, post_set_name, post_set_func ,value)
        setattr(cls, key, value)

    return setter

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
     - allowed:
        List or Tuple of allowed values (if not iterable) or sub_values for the property.
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
        lprop   = []
        for k,v in dct.items():
            new_dct[k] = v
            if not isinstance(k, str) or k.startswith("_") or not isinstance(v, dict):
                continue

            lprop.append(k)
            hk        = "_" + k
            settable  = v.pop('settable',  True)
            typ       = v.pop('typ',       None)
            sub_typ   = v.pop('sub_typ',   None)
            shape     = v.pop('shape',     0)
            default   = v.pop('default',   None)
            conv_func = v.pop('conv_func', lambda x: x)
            allowed   = v.pop('allowed',   [])

            pre_set_name = v.pop('pre_set_name', None)
            pre_set_func = v.pop('pre_set_func', lambda cls,x: x)

            post_set_name = v.pop('post_set_name', None)
            post_set_func = v.pop('post_set_func', lambda cls,x: x)

            doc       = v.pop('doc',       '')

            if v:
                raise ValueError(f"Invalid keywords '{list(v.keys())}' for metaclass 'PropertyCreator'.")

            doc = ('type          = {}.\n'
                   'sub_type      = {}.\n'
                   'allowed       = {}.\n'
                   'default_value = {}.\n\n'
                   '{}').format(typ, sub_typ, allowed, default, doc.lstrip().replace('\t', '')) 

            new_dct[hk] = set_typ(typ, default)

            getter = key_getter(hk)
            setter = None
            if settable:
                setter = key_setter(
                    hk, typ, sub_typ, 
                    shape,
                    conv_func,
                    allowed,
                    pre_set_name, pre_set_func,
                    post_set_name, post_set_func)

            new_dct[k] = property(getter, setter, doc=doc)

        res = super().__new__(cls, clsname, base, new_dct)

        old_init = new_dct.get('__init__', None)
        def __init__(self, *args, **kwargs):
            for p in lprop:
                app = kwargs.pop(p, None)
                if not app is None:
                    setattr(self, p, app)

            if not old_init is None:
                old_init(self, *args, **kwargs)
            else:
                super(res, self).__init__(*args, **kwargs)

        if not old_init is None:
            __init__.__doc__ = old_init.__doc__
        new_dct['__init__'] = __init__

        res.__init__ = __init__

        return res



