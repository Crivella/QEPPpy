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
	'x':xmlacquisition rule, ('attr'/text/allchild/nodelist/nodelistnested)
	'n':name 
}
"""
class data_file_parser():
	data={}
	__name__ = "data_file_parser"
	__toinit__ = []
	def __init__( self, fname="", **kwargs):
		for i in self.__toinit__:
			self.__dict__[i] = None
		if fname:
			#if not os.path.isfile( fname):
			#	raise IOError( "File '{}' does not exist".format( fname))
			self.fname = fname
			self.parse_xml( fname)
		elif kwargs:
			for k, v in kwargs.items():
				if k in self.__dict__:
					self.__dict__[k] = v
				else:
					#pass
					raise Exception( "{}: Unrecognized keyword argument '{}'.\n".format( self.__name__, k))
		else:
			pass
			#raise Exception( "{}: Failed to initialize object.\n".format( self.__name__))
		return

	def __getattr__( self, key):
		if key in self.__dict__:
			return self.__dict__[key]
		else:
			raise AttributeError( "'{}' object has no attribute '{}'".format( self.__name__, key))

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

		def _xml_attr_( node, f, n):
			if f: node = node.find( f)
			return _format_( node.get( n))

		def _xml_text_( node, f, n):
			if f: node = node.find( f)
			return _format_( node.text)

		def _xml_all_child_( node, f, n):
			if f: node = node.find( f)
			ret = []
			for child in node.getchildren():
				text = child.text.strip()
				#print( child.attrib, ", '{}'".format( text))
				if text:
					l = text.split( " ")
					l = list( map( _format_, l))
					if len( l) == 1: l = l[0]
					if n:
						d = child.attrib.copy()
						for k, v in d.items():
							d[k] = _format_( v)
						d[n] = l
					else: 
						d = l
					ret.append( d)
				else:
					d = child.attrib
					#print( "TEST, ", child.getchildren())
					for c in child.getchildren():
						d[c.tag] = _format_( c.text)
					ret.append( d)
			return ret

		def _xml_node_list_( node, f, n):
			if f: node=node.findall( f)
			ret = []
			for c in node:
				d={}
				#add = c.text.strip().split( "\n")
				add = c.text.strip().replace( "\n", " ")
				add = list( filter( None, add.split( " ")))
				#add = list( map( lambda x: list( filter( None, x.split( " "))), add))
				if len( add) == 1: add = add[0]
				#if len( add) == 1: add = add[0]
				if isinstance( add, list):
					add = np.array( add, dtype=float)
				d[c.tag] = _format_( add)
				d.update( c.attrib)
				ret.append( d)
			return ret

		def _xml_node_list_nest_( node, f, n):
			if f: node=node.findall( f)
			ret = []
			for child in node:
				d={}
				for c in child:
					#add = c.text.strip().split( "\n")
					add = c.text.strip().replace( "\n", " ")
					add = list( filter( None, add.split( " ")))
					#add = list( map( lambda x: list( filter( None, x.split( " "))), add))
					if len( add) == 1: add = add[0]
					#if len( add) == 1: add = add[0]
					if isinstance( add, list):
						add = np.array( add, dtype=float)
					d[c.tag] = _format_( add)
					d.update( c.attrib)
				ret.append( d)
			return ret

		xml_acq_rule={
			'attr':_xml_attr_,
			'text':_xml_text_,
			'allchild':_xml_all_child_,
			'nodelist':_xml_node_list_,
			'nodelistnested':_xml_node_list_nest_
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


