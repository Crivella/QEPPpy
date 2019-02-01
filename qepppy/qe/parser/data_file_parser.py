import numpy as np
from ...logger import logger, warning


def _format_(var):
	if isinstance(var, np.ndarray): return var
	if var == None: return None
	try: v = int(var)
	except:
		try: v = float(var)
		except:
			if var == "true": v = True
			elif var == "false": v = False
			else: v = var
	return v

matches = {
	int:   r'[\s]*(?P<flag>[\d\-]+)',
	float: r'[\s]*(?P<flag>[\d\.\-EedD]+)',
	str:   r'[\s]*(?P<flag>.*)',
	bool:  r'[\s]*(?P<flag>.)',
}

def _get_value_(f, search_data, dtype=str):
	import re

	string = search_data['outfile_regex']
	m = search_data.get('m', 1)

	a = None
	try:
		if dtype == list:
			a = re.finditer(string, f)
			a = [ x.groupdict() for x in a]
			for n,e in enumerate(a):
				for k,v in e.items():
					b = np.fromstring(v, sep=' ')
					if len(b) == 0:
						b = str(v).strip()
					elif len(b) == 1:
						b = b[0]
					a[n][k] = b*m
		else:
			a = re.search(string + matches[dtype], f).group('flag')
	except Exception as e:
		#print(e)
		pass

	val = None
	try:
		val = dtype(a)
	except:
		warning.print( "Failed to convert '{}'(from '{}') to dtype '{}'".format(a, string, dtype))
	return val

def _xml_attr_(node, f="", n=""):
	if f:
		node = node.find(f)
	return _format_(node.get(n))

def _xml_text_(node, f="", n=""):
	if f:
		node = node.find(f)
	return _format_(node.text)

def _xml_node_list_(node, f="", n=""):
	if f: node=node.findall(f)
	ret = []
	for c in node:
		d=c.attrib
		d = {k:_format_(v) for k,v in d.items()}
		add = c.text.strip().replace( "\n", " ")
		if add:
			add = list(filter(None, add.split( " ")))
			if len(add) == 1:
				add = add[0]
			if isinstance(add, list):
				add = np.array(add, dtype=float)
			if n: 
				tag = n
			else: 
				tag = c.tag
			d[tag] = _format_(add)
		else:
			for e in  _xml_node_list_(c.getchildren()):
				d.update(e)
		ret.append(d)

	return ret

@logger()
class data_file_parser(object):
	"""
	Parser for QE data'file'schema.xml (QE>=6.2)

	- d: rule dictionary defining how the parser will operate.
	- schema: Name of the "data-file*.xml" to parse (Parsing will run if the schema is set).
	- **kwargs: Overwrite parsed variables with user's one given as a var_name=value keyword arg.
	            The var_name must already exist after the parsing (cannot set unrecognized variables).
	This parser accept a rule dict passed as a keyword argument 'd' to the __init__ method.
	The rule dict has to follow the format:
	{
	'varname':{
	    'res_type':type, (int/float/ndarray/...)
	    'xml_search_string':element to find, (string)
	    'xml_ptype':xmlacquisition rule, ('attr'/text/nodelist)
	    'extra_name':name
	    }
	'varname1':{...}
	}
	The parsing will generate internal variables using as name the keys of the rule dict.
	The rule are as follow:
	- attr: Get the value of the attribute of name "n" of the node given by root.find(f)
	- text: Get the text of the node given by root.find(f)
	- nodelist: Find a list of nodes using root.findall(f) and analyze them.
	            The resulting variable will be a list dictionary (one element for every node found).
	            "n" is used as a key name for the text value of the node, otherwise the node.tag is used.
	            text that contains arrays of number are automatically converted into np.ndarray objects
	            All the node attributes are saved are saved as key:values in the dict.
	            If the node contains children instead of text, explore them recursively:
	            - saves all the attributes found in the dict as key:value pairs.
	            - saves all the the text values found in the dict as node.tag:text pairs.
	            NOTE: no conflict resolution is present for attributes with the same name across children.
	                  The attribute are updated at every step, so the value of the last one will be stored.
	"""
	__name__ = "data_file_parser"
	def __init__(self, schema="", outfile="", d={}, **kwargs):
		self._data_ = d
		for i in d:
			self.__dict__[i] = None
		if schema:
			self.schema = schema
			self._parse_xml_()
		elif outfile:
			self.outfile = outfile
			self._parse_outfile_( )
		if kwargs:
			for k, v in kwargs.items():
				if not k in self.__dict__:
					continue
				self.__dict__[k] = v
		return

	def __getitem__(self, key):
		return self.__dict__.get(key)

	def __str__(self):
		return ""

	def _parse_outfile_(self):
		with open(self.outfile, "r") as f:
			content = f.read()

		for k, v in self._data_.items():
			t = v['res_type']
			search = v.get( 'outfile_regex', None)
			if search is None:
				continue
			#print(k,v)
			#print(search)

			val = None
			try:
				val = _get_value_(content, v, dtype=t)
			except Exception as e:
				print( "ERROR: ", e)
			#print(val)
			self.__dict__[k] = val
		return

	def _parse_xml_(self):
		import xml.etree.ElementTree as ET
		root = ET.parse(self.schema).getroot()

		xml_acq_rule={
			'attr':_xml_attr_,
			'text':_xml_text_,
			'nodelist':_xml_node_list_,
		}

		for k, v in self._data_.items():
			res = None
			t = v['res_type']
			n = v['extra_name']
			f = v['xml_search_string']
			func = xml_acq_rule[v['xml_ptype']]
			try: 
				res = func(root, f, n)
			except Exception as e: 
				continue # raise Exception( "{}, {}".format(k,e));print(k, e);
			#print(k, t, res)
			#if res == None: continue
			if t == np.array: 
				self.__dict__[k] = t(res, dtype=float)
			else: 
				self.__dict__[k] = t(res)
		return

	def validate(self):
		try:
			return True and super().validate()
		except:
			return True


