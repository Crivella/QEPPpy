import pytest
import re
import numpy as np
from qepppy.meta.property_creator import PropertyCreator
# import qepppy.meta.Pr as meta

initl_corr={
	'empty':{},
	'partial':{
		'p1':123,
		'p2':'123',
		},
	'fix_dim_shape':{
		'p4':np.random.rand(10,25),
		},
	'arg_shape':{
		'p1':3,
		'p5':np.random.rand(24,44,3),
		},
	'all':{
		'p1':3,
		'p2':'123',
		'p3':np.arange(3),
		'p4':np.random.rand(10,25),
		'p5':np.random.rand(24,44,3),
		'p6':np.arange(9).reshape(3,3)
		},
}

class A(metaclass=PropertyCreator):
	p1={
		}

	p2={
		'typ':(str,),
		'doc':"""C'era una volta tanto tempo fa..."""
		}

	p3={
		'typ':(np.ndarray,),
		}
	p4={
		'typ':(list,tuple,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape':(-1,-1),
		}
	p5={
		'typ':(list,tuple,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape':(-1,-1, 'p1'),
		}

	p6={
		'typ':(list,tuple,np.ndarray),
		'sub_typ':(int, np.integer),
		'shape':(-1,3),
		'allowed':range(10),
		'conv_func':lambda x: np.array(x, dtype=int)
		}

	p7={
		'pre_set_name':'_p1',
		}

	p8={
		'post_set_name':'_p1',
		'post_set_func':lambda cls, x: x - cls.p1,
		}

	default={
		'default':1
		}

class B(A):
	pass

@pytest.fixture(
	scope='session',
	params=list(initl_corr.keys())
	)
def init_kwargs_corr(request):
	return initl_corr[request.param]

@pytest.mark.parametrize(
	'val',
	[1, 1., '1', [1,2,3], (1,2,3), np.arange(10), None]
	)
def test_meta_default(val):
	class app(metaclass=PropertyCreator):
		test_def = {'default':val}
		
	new = app()
	get = getattr(new, 'test_def')

	check = (get == val)
	if isinstance(get, np.ndarray):
		check = np.allclose(get, val)

	assert check, f"Failed to set default value: set'{get}' vs to_be'{val}'." 


def test_meta_init_corr(init_kwargs_corr):
	A(**init_kwargs_corr)

def test_meta_init_inheritance(init_kwargs_corr):
	B(**init_kwargs_corr)

def test_meta_doc():
	assert "C'era una volta tanto tempo fa..." in A.p2.__doc__, (
		"Documentation was not assigned correctly to property."
		) 


@pytest.mark.parametrize('num', range(5))
def test_meta_conv_func(num):
	new  = A()
	new.p6 = [[1,2,3]]*num

	assert isinstance(new.p6, np.ndarray), "Failed to convert list to np.ndarray"
	assert new.p6.shape[0] == num, f"Failed conversion: wrong shape '{A.p6.shape}'"

def test_meta_allowed_wrong():
	new    = A()
	try:
		new.p6 = np.arange(12).reshape(-1,3)
	except ValueError as e:
		assert not re.search(r"Value '10' is not allowed", str(e)) is None, "Wrong error msg: " + str(e)
		return

@pytest.mark.parametrize(
	'name_val',
	[
		('p2','123'),
		('p4',np.random.rand(3,5)),
		('p6',np.random.randint(0,9,(10,3))),
	]
	)
def test_meta_reset_attr(name_val):
	name, val = name_val
	new  = A()

	assert not getattr(new, name), f"Wrong initialization of attr."
	setattr(new, name, val)
	v = getattr(new, name)
	if isinstance(v, np.ndarray) and np.shape != (0,):
		v = True
	assert v, "Failed to assign attr."
	setattr(new, name, None)
	assert not getattr(new, name), "Failed to reset attr."

def test_meta_pre_set():
	new    = A()
	new.p7 = 12345

	assert new.p1 == 12345, "Failed pres_set."

def test_meta_post_set_corr():
	new    = A()
	new.p1 = 444
	new.p8 = 12345

	assert new.p1 == 12345-444, "Failed post_set."

def test_meta_post_set_wrong():
	new    = A()
	try:
		new.p8 = 12345
	except:
		return

	raise AssertionError("post_set that should fail did not fail.")



