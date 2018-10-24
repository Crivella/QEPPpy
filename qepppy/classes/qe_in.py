import sys

"""
Import a class capable of generating a namelist dict template by parsing a file
This class should provide the following methods:
	-parse(): Function to read the file and generate a namelist template
	-validate(): Function to check namelist template after reading
	-convert(): Function to set a value in the namelist template
	-set:() Function to convert the list in a string QE input file
"""
from .qe_doc import qe_doc_parser as parser

def trim_ws( str):
	#Trim all withspace not included in a string
	ws=[ " ", "\t", "\n"]
	l = str.strip()
	new = ""
	check_str_1 = False
	check_str_2 = False
	for e in l:
		if not e in ws or check_str_1 or check_str_2:
			new += e
		if e == "\""  and not check_str_2:
			check_str_1 = not check_str_1
		if e == "'" and not check_str_1:
			check_str_2 = not check_str_2
	return new


pw_nl={}#qe_doc_read( fname="INPUT_PW.def")

class qe_in( parser):
	def __init__( self, fname="" , templ_file="", **kwargs):
		#A namelist template '_d' must be declared in the child class!!!!!!
		if templ_file:
			self.templ_file = templ_file
		if not self.templ_file: raise Exception( "Must give a template file.\n")
		self.parse( self.templ_file)

		if fname: self.namelist_read( fname=fname)
		#Check if initialization keyword arguments are compliant with the given namelist template
		for nl, v in kwargs.items():
			if not isinstance( v, dict): raise Exception( "Invalid kwargs.\n{}".format( kwargs))
			for k, v1 in v.items():
				self.set( nl, k, v1)

		return

	def __str__( self):
		self.validate()
		return self.convert()

	def fprint( self, fname=""):
		"""
		Print this or a child class __str__() to a file or to the stdout
		"""
		if fname: f = open( fname, "w+")
		else: f = sys.stdout

		f.write( self.__str__())

		if fname: f.close()

		return


	def namelist_read( self, fname=""):
		"""
		Read a the namelists of an input file
		"""
		if not fname:
			raise Exception( "Must pass a filename to open")

		#Read all the file content into 'content'
		with open(fname) as f:
			content = f.readlines()
		
		nl = None
		for l in content:
			#Ignore comments
			ls = l.strip().split( "!")[0]
			if not ls: continue
			#CASE: Namelist name
			if '&' == ls[0]:
				nl = ls[1:].upper()
				if not nl in self._d['nl']:
					raise Exception( "Reading unrecognized namelist '{}'".format( nl))
			#CASE other
			else:
				#(not '/' in l or '=' in l) => Recognize namelist field from otehr fields or namelist end '/'
				#nl => Recognize if reading field outside of namelist
				if (not '/' in ls or '=' in ls) and nl:
					if not nl: raise Exception( "Corrupted input file:\n{}".format( l))
					#Read namelist fields separated by endline ('\n') or by commas (',')
					for e in filter( None, ls.split( ",")):
						l1 = trim_ws(e).split( "=")
						if len(l1) != 2: raise Exception( "Corrupt input file:\n{}".format( l))
						v = l1[1].replace("\"", "").replace("'", "")
						
						#Check if the field/parameter name is present in the namelist template
						try: self.set( nl=nl, k=l1[0], v=v)
						except NameError as e: print( e)
				else:
					nl = None

		return






