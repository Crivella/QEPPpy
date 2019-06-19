import pytest
import re
import numpy as np
import qepppy.meta as meta

class appdict(dict):
	def __getattr__(self, key):
		return self[key]

	def __setattr__(self, key, value):
		self[key] = value

def wrong(func, exc, *args):
	args = list(args)
	err = args.pop()
	try:
		func(*args)
	except exc as e:
		assert not re.search(err, str(e)) is None, (
			'Exception msg did not match the expected one\n' + 
			'expected:  ' + err + '\n' +
			'obtained:  ' + str(e)
			)
		return
	except Exception as e:
		raise Exception("Wrong exception obtained:  \n" + str(e))

	raise Exception("Check_type that should fail did not fail: " + str(args))


iterables = [
		None, '', (), [], 
		1, 1.5, 
		[1,2,[3,(4, [12])]],
		(1,2,[3,(4,(5,))]),
		np.arange(16).reshape(2,2,2,2)
	]
@pytest.mark.parametrize(
	'iterable',
	iterables,
	ids=['iterable' + str(n) for n in range(len(iterables))]
	)

def test_flatten_iter(iterable):
	meta.flatten_iter(iterable)

@pytest.mark.parametrize(
	'typ_value',
	[
		(
			None,
			1
		),
		(
			int,
			1
		),
		(
			float,
			1.
		),
		(
			(int,float,np.number),
			2
		),
		(
			(int,float,np.number),
			2.3
		),
		(
			(int,float,np.number),
			np.int(2)
		),
		(
			(str,),
			'boh66554'
		),
		(
			(list,),
			[1,2]
		),
		(
			(tuple,),
			(1,2)
		),
		(
			(list,tuple,),
			[1,2]
		),
		(
			(list,tuple,),
			(1,2)
		),
		(
			(np.ndarray,),
			np.arange(3)
		),
	]
	)
def test_check_type_corr(typ_value):
	typ,value = typ_value
	meta.check_type(typ, value)

err = r"value='.*' must be of type '.*'"
@pytest.mark.parametrize(
	'typ_value',
	[
		(
			int,
			1.,
			err
		),
		(
			float,
			1,
			err
		),
		(
			(int,float,np.number),
			'a',
			err
		),
		(
			(int,float,np.number),
			range(1,2),
			err
		),
		(
			(int,float,np.number),
			[1,3],
			err
		),
		(
			(list,),
			(1,2),
			err
		),
		(
			(tuple,),
			[1,2],
			err
		),
		(
			(list,tuple,),
			"[1,2]",
			err
		),
		(
			(list,tuple,),
			range(1,2),
			err
		),
		(
			(np.ndarray,),
			[1,3],
			err
		),
	]
	)
def test_check_type_wrong(typ_value):
	wrong(meta.check_type, TypeError, *typ_value)

@pytest.mark.parametrize(
	'subtyp_value',
	[
		(
			None,
			[1,1.,('a', range(5))]
			),
		(
			int,
			[1,2,3,4,5]
		),
		(
			float,
			[1.,2.,3.,4.,5.]
		),
		(
			(int,float,np.number),
			np.arange(5)
		),
		(
			(int,float,np.number),
			np.linspace(1,10,10)
		),
		(
			(int,float,np.number),
			[1,3, (3., np.int(5))]
		),
		(
			(str,),
			['a','b',('c',('d',)), 'e']
		),
	]
	)
def test_check_sub_type_corr(subtyp_value):
	subtyp,value = subtyp_value
	meta.check_sub_type(subtyp, value)

err = r"elements of value='.*' must be of type '.*'\."
@pytest.mark.parametrize(
	'subtyp_value',
	[
		(
			int,
			[1,2,3,4,5.],
			err
		),
		(
			int,
			['a',2,3,4,5.],
			err
		),
		(
			float,
			[1.,2.,3.,4.,5],
			err
		),
		(
			(int,float,np.number),
			['a',2,3,4,5.],
			err
		),
		(
			(str,),
			['a','b',('c',('d',)), 'e', 2],
			err
		),
	]
	)
def test_check_sub_type_wrong(subtyp_value):
	wrong(meta.check_sub_type, TypeError, *subtyp_value)

@pytest.mark.parametrize(
	'cls_value_shape',
	[
		(
			appdict([('a',1)]),
			np.empty((3,2,1)),
			(3,2,'a')
		),
		(
			appdict([('a',1),('b',2),('c',3),]),
			np.empty((3,2,1)),
			('c','b','a')
		),
		(
			appdict(),
			np.empty((3,2,1)),
			(-1,2,-1)
		),
		(
			appdict(),
			np.empty((3,2,1)),
			None
		),		
		(
			appdict(),
			[1,2,3,4],
			(-1,)
		),
		(
			appdict(),
			[1,2,3,4],
			(4,)
		),
		(
			appdict(),
			[1,2,3,4,5,6,7,8,9],
			(9,)
		),
	]
	)
def test_check_shape_corr(cls_value_shape):
	cls, value, shape = cls_value_shape
	meta.check_shape(cls, value, shape)

@pytest.mark.parametrize(
	'cls_value_shape',
	[
		(
			appdict([('a',1)]),
			np.empty((3,2,2)),
			(3,2,'a'),
			r"Shape mismatch '.*'=.* not equal to '.*'"
		),
		(
			appdict(),
			np.empty((3,2,1)),
			(-1,3,-1),
			r"Shape mismatch '.*'=.* not equal to '.*'"
		),
		(
			appdict(),
			np.empty((3,2,1)),
			('a',2,-1),
			r"Attribute '.*' is not present in '.*'"
		),
		(
			appdict([('a',2)]),
			np.empty((3,2,1)),
			('a',2,-1),
			r"Shape mismatch '.*'=.* not equal to '.*'"
		),
		(
			appdict([('a',2)]),
			np.empty((3,2,1)),
			('a',2,),
			r"Mismatch in number of dimensions: value='.*' vs required='.*'"
		),
		(
			appdict(),
			[1,2,3],
			(4,),
			r"Shape mismatch '.*'=.* not equal to .*"
		),
		(
			appdict(),
			[1,2,3],
			(3,1),
			r"Can't confront shape of value='.*' with '.*'"
		),
		(
			appdict(),
			[1,2,3],
			('a',),
			r"Attribute '.*' is not present in '.*'"
		),
	]
	)
def test_check_shape_wrong(cls_value_shape):
	wrong(meta.check_shape, ValueError, *cls_value_shape)

@pytest.mark.parametrize(
	'val_al_non',
	[
		(
			[1,2,3,4],
			range(10),
			None
		),
		(
			(0,1,0),
			[0,1],
			None
		),
		(
			'Y',
			['Y','N'],
			None
		),
		(
			'',
			['Y','N'],
			''
		),
	]
	)
def test_check_allowed_corr(val_al_non):
	val, al, non = val_al_non
	meta.check_allowed(val, al, non)

@pytest.mark.parametrize(
	'val_al_non',
	[
		(
			[1,2,3,4,11],
			range(10),
			None,
			r"Value '.*' is not allowed"
		),
		(
			(0,1,0,2),
			[0,1],
			None,
			r"Value '.*' is not allowed"
		),
		(
			'Yes',
			['Y','N'],
			None,
			r"Value '.*' is not allowed"
		),
		(
			'',
			['Y','N'],
			None,
			r"Value '.*' is not allowed"
		),
	]
	)
def test_check_allowed_wrong(val_al_non):
	wrong(meta.check_allowed, ValueError, *val_al_non)