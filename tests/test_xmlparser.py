import pytest
import os
from qepppy.parsers.xml_parser import Parser_xml

cwd = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture(scope='module')
def cls():
	xml = os.path.join(cwd, 'qe/test_files/1_bands.xml')
	sch = os.path.join(cwd, 'qe/test_files/qes-1.0.xsd')
	new = Parser_xml(xml=xml, schema=sch)

	return new

@pytest.fixture(scope='module')
def cls1():
	xml = os.path.join(cwd, 'test.xml')
	new = Parser_xml(xml=xml)

	return new


def test_root(cls1):
	res = cls1.find('//')
	assert res['@version'] == '1.0'

def test_full_path(cls1):
	res = cls1.find('t1/t3')
	assert isinstance(res, dict)

def test_xpath_single(cls1):
	res = cls1.find('//t2')
	assert isinstance(res, dict)

def test_xpath_multiple_crossnode(cls1):
	res = cls1.find('//t3')
	assert isinstance(res, list)
	assert len(res) == 2

def test_xpath_multiple_same_node(cls1):
	res = cls1.find('//t6')
	assert isinstance(res, list)
	assert len(res) == 9

def test_xpath_multiple_cross_same_node(cls1):
	res = cls1.find('//t10')
	assert isinstance(res, list)
	assert len(res) == 6


# def test_validate_schema(cls):
# 	cls.validate_schema()

@pytest.mark.parametrize(
	'xpath',
	[
		'parallel_info',
		'parallel_info/nthreads',
		'//nthreads',
		'//cell//a1',
		'//a1',
		'input//a1',
		'//eigenvalues',
		'//rotation',
	]
	)
def test_xmlparser_xpath_corr(cls, xpath):
	if not cls.find(xpath):
		raise AssertionError("Failed to properly parse")


@pytest.mark.parametrize(
	'xpath',
	[
		'parallel_info/boh',
		'nthreads',
		'cell//a1',
		'input//a4',
		'input//eigenvalues',
		'/rotation',
	]
	)
def test_xmlparser_xpath_wrong(cls, xpath):
	if cls.find(xpath):
		raise AssertionError("Failed to properly parse")

