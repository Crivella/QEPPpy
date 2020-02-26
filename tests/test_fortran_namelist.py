import pytest
from qepppy.parsers import fortran_namelist_collection as fnc
from qepppy.parsers.fortran_namelist import fortran_namelist as fn

string = """
&NAME1
  p1 = '123',
  p2 = "321",
  p3 = 1,
  p4 = 1.15
  p5 = 1, 2, 3
  p6(1) = 1
  p6(3) = 3
  p7 = 1, p8 = 2, ! p9 = 3
  ! p10 = 1
  p12 = .TRUE., p13 = .FaLsE.
  p14 = 1e-05, p15=1e3
/
&NAME2
  p1  = 'testo'
  p11 = 'only me'
/

"""

@pytest.fixture(scope='module')
def load_str():
	from io import StringIO
	buffer = StringIO(string)
	return fnc(src=buffer)

def test_from_str(load_str):
	pass

def test_single_quote(load_str):
	assert load_str['NAME1']['p1'] == '123', 'Failed to parse the single quoted element'

def test_double_quote(load_str):
	assert load_str['NAME1']['p2'] == '321', 'Failed to parse the double quoted element'

def test_int(load_str):
	assert load_str['NAME1']['p3'] == 1, 'Failed to parse integer value'

def test_float(load_str):
	assert load_str['NAME1']['p4'] == 1.15, 'Failed to parse float value'

def test_list(load_str):
	assert load_str['NAME1']['p5'] == [1,2,3], 'Failed to parse list value'

def test_case_bool(load_str):
	assert load_str['NAME1']['p12'] == True, 'Failed to parse case insensitive bool'
	assert load_str['NAME1']['p13'] == False, 'Failed to parse case insensitive bool'

def test_exp_notation(load_str):
	assert load_str['NAME1']['p14'] == 1e-5, 'Failed to parse exponential notation'
	assert load_str['NAME1']['p15'] == 1e+3, 'Failed to parse exponential notation'

def test_list_assigned_1by1(load_str):
	assert load_str['NAME1']['p6'] == [1,None,3], 'Failed to parse list value'

def test_list_retrive_single_value(load_str):
	assert  load_str['NAME1']['p5(2)'] == 2, 'Failed to retrieve single value from list'

def test_deep_find_multiple(load_str):
	assert load_str['p1'] == '123', 'Failed to retrieve first value from deep_find'

def test_deep_find_single(load_str):
	assert load_str['p11'] == 'only me', 'Failed to retrieve value from deep_find'

def test_deep_find_list(load_str):
	assert load_str['p6(3)'] == 3, 'Failed to retrieve list value from deep_find'

def test_multiple_same_line(load_str):
	assert load_str['p7'] == 1, 'Failed to assign multiple params on one line'
	assert load_str['p8'] == 2, 'Failed to assign multiple params on one line'

def test_inline_comment(load_str):
	assert load_str['p9'] is None, 'Not ignored an inline comment'

def test_normal_comment(load_str):
	assert load_str['p10'] is None, 'Not ignored an normal comment'

def test_tokenize_pattern(load_str):
	assert load_str['NAME1/p1'] == '123', 'Failed to use tokenized pattern'

def test_namelist_pattern(load_str):
	assert isinstance(load_str['NAME2'], fn), 'Failed to return nemelist object'

def test_case_insensitive(load_str):
	assert load_str['NAME1/p1'] == '123'
	assert load_str['NaMe1/P1'] == '123'
	assert isinstance(load_str['NAME1'], fn)

def test_manual_set(load_str):
	load_str['NAME1']['a1'] = 2
	assert load_str['a1'] == 2, 'Failed to assign value going through namelist'

def test_pattern_set(load_str):
	load_str['NAME1/a2'] = 3
	assert load_str['a2'] == 3, 'Failed to assign value going through pattern'


def test_pattern_set_new_nl(load_str):
	load_str['NAME3/a3'] = 1
	assert load_str['a3'] == 1, 'Failed to assign value in new namelist'



