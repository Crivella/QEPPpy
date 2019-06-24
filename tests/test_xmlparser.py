import pytest
import os
from qepppy.parsers.xmlschema_parser import Parser_xmlschema

cwd = os.path.dirname(os.path.realpath(__file__))

@pytest.fixture(scope='module')
def cls():
	xml = os.path.join(cwd, 'qe/test_files/1_bands.xml')
	sch = os.path.join(cwd, 'qe/test_files/qes-1.0.xsd')
	new = Parser_xmlschema(xml=xml, schema=sch)

	return new

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
def test_xmlschemaparser_xpath_corr(cls, xpath):
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
def test_xmlschemaparser_xpath_wrong(cls, xpath):
	if cls.find(xpath):
		raise AssertionError("Failed to properly parse")

