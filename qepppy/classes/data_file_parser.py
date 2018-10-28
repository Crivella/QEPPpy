import os
import numpy as np


"""
Parser for QE data'file'schema.xml (QE>=6.2)
This parser will generate internal variables.
The name of this var and the rules for the acquisition are dictated by the dictionary data.
This dicitonary is supposed to be set in the child class
Example of an element in data:
'varname':{
	't':type, (int/float/ndarray/...)
	'f':element to find, (string)
	'x':xmlacquisition rule, ('attr'/text/nodelist)
	'n':name 
}
"""

class data_file_parser( object):
	data={}
	__name__ = "data_file_parser"
	def __init__( self, fname="", d={}, **kwargs):
		#print( d)
		self.data = d
		for i in d:
			self.__dict__[i] = None
		if fname:
			self.fname = fname
			self.parse_xml( fname)
		if kwargs:
			for k, v in kwargs.items():
				if k in self.__dict__:
					self.__dict__[k] = v
				else:
					raise Exception( "{}: Unrecognized keyword argument '{}'.\n".format( self.__name__, k))
		else:
			pass
			#raise Exception( "{}: Failed to initialize object.\n".format( self.__name__))
		return

	def __getitem__( self, key):
		return self.__dict__.get( key)

	def __str__( self):
		return ""

	def parse_xml( self, fname=""):
		import xml.etree.ElementTree as ET
		root = ET.parse( fname).getroot()

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

		xml_acq_rule={
			'attr':_xml_attr_,
			'text':_xml_text_,
			'nodelist':_xml_node_list_,
		}

		for k, v in self.data.items():
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
		return True


