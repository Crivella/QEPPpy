import numpy as np
from ...logger import logger


def _format_( var):
	if isinstance( var, np.ndarray): return var
	if var == None: return None
	try: v = int( var)
	except:
		try: v = float( var)
		except:
			if var == "true": v = True
			elif var == "false": v = False
			else: v = var
	return v


def _get_value_( f, string, delim='=', dtype=str):
	app = f.find( string)
	val = f[app:].split( '\n')[0].split( delim)[1]
	try:
		val = dtype( list( filter( None, val.split( " ") ))[0])
	except:
		raise Exception( "Failed to convert '{}' to dtype '{}'".format( val, dtype))
	return val

def _xml_attr_( node, f="", n=""):
	if f: node = node.find( f)
	return _format_( node.get( n))

def _xml_text_( node, f="", n=""):
	if f: node = node.find( f)
	return _format_( node.text)

def _xml_node_list_( node, f="", n=""):
	if f: node=node.findall( f)
	ret = []
	for c in node:
		d=c.attrib
		#add = c.text.strip().split( "\n")
		add = c.text.strip().replace( "\n", " ")
		if add:
			add = list( filter( None, add.split( " ")))
			#add = list( map( lambda x: list( filter( None, x.split( " "))), add))
			if len( add) == 1: add = add[0]
			#if len( add) == 1: add = add[0]
			if isinstance( add, list):
				add = np.array( add, dtype=float)
			if n: tag = n
			else: tag = c.tag
			d[tag] = _format_( add)
		else:
			for e in  _xml_node_list_( c.getchildren()):
				d.update( e)
		ret.append( d)

	return ret

@logger( )
class data_file_parser( object):
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
	    't':type, (int/float/ndarray/...)
	    'f':element to find, (string)
	    'x':xmlacquisition rule, ('attr'/text/nodelist)
	    'n':name
	    }
	'varname1':{...}
	}
	The parsing will generate internal variables using as name the keys of the rule dict.
	The rule are as follow:
	- attr: Get the value of the attribute of name "n" of the node given by root.find( f)
	- text: Get the text of the node given by root.find( f)
	- nodelist: Find a list of nodes using root.findall( f) and analyze them.
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
	def __init__( self, schema="", outfile="", d={}, **kwargs):
		self._data_ = d
		for i in d:
			self.__dict__[i] = None
		if schema:
			self.schema = schema
			self.parse_xml( schema)
		elif outfile:
			self.outfile = outfile
			self.parse_outfile( )
		if kwargs:
			for k, v in kwargs.items():
				if not k in self.__dict__:
					continue
				self.__dict__[k] = v
		return

	def __getitem__( self, key):
		return self.__dict__.get( key)

	def __str__( self):
		return ""

	def parse_outfile( self):
		with open( self.outfile, "r") as f:
			content = f.read()

		return

	def parse_xml( self):
		import xml.etree.ElementTree as ET
		root = ET.parse( self.schema).getroot()

		

		xml_acq_rule={
			'attr':_xml_attr_,
			'text':_xml_text_,
			'nodelist':_xml_node_list_,
		}

		for k, v in self._data_.items():
			res = None
			t = v['t']
			n = v['n']
			f = v['f']
			func = xml_acq_rule[v['x']]
			try: res = func( root, f, n)
			except Exception as e: continue # raise Exception( "{}, {}".format(k,e));print( k, e);
			#print( k, t, res)
			#if res == None: continue
			if t == np.array: self.__dict__[k] = t( res, dtype=float)
			else: self.__dict__[k] = t( res)
		return

	def validate( self):
		try:
			return True and super().validate()
		except:
			return True


