import sys

"""
Import a class capable of generating a namelist dict template by parsing a file
This class should provide the following methods:
	-parse(): Read the file and generate an internal namelist template
	-validate(): Check namelist template after reading
	-convert():  Convert the internal dict in a string QE input file
	-check_nl( nl="namelist"): Check if nl is valid (present in the itnernal namelist)
	-set_nl:( nl="namelist", k="param", v="value to set") Set a namelist value in the namelist template
	-set_card: ( card="", v="", el=[]) Set a card value in the namelist template
		if v is set, set the card main value
		if el is set set a card list value
			el must be an entire lie to parse
	-get:( nl="namelist", k="param") Retrieve the value of a parameter
	-find: (name) Find a variable with name=name in the namelist template
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

class qe_in( parser):
	def __init__( self, fname="" , templ_file="", **kwargs):
		#A namelist template '_d' must be declared in the child class!!!!!!
		if templ_file:
			self.templ_file = templ_file
		if not self.templ_file: raise Exception( "Must give a template file.\n")
		self.parse( self.templ_file)

		if fname: self.in_parse( fname=fname)
		#Check if initialization keyword arguments are compliant with the given namelist template
		for nl, v in kwargs.items():
			if not isinstance( v, dict): raise Exception( "Invalid kwargs.\n{}".format( kwargs))
			for k, v1 in v.items():
				self.set_nl( nl, k, v1)

		return

	def __str__( self):
		return self.convert()

	def fprint( self, fname=""):
		"""
		Print this or a child class __str__() to a file or to the stdout
		"""
		if fname: f = open( fname, "w+")
		else: f = sys.stdout

		self.validate()
		f.write( self.__str__())

		if fname: f.close()

		return


	def in_parse( self, fname=""):
		"""
		Read a the namelists of an input file
		"""
		if not fname:
			raise Exception( "Must pass a filename to open")

		#Read all the file content into 'content'
		with open(fname) as f:
			content = f.readlines()
		
		nl = None
		card = None
		for l in content:
			#Ignore comments
			ls = l.strip().split( "!")[0]
			if not ls: continue
			#CASE: Namelist name
			if '&' == ls[0]:
				nl = ls[1:].upper()
				if not self.check_nl( nl):
					raise Exception( "Reading unrecognized namelist '{}'".format( nl))
			#CASE other
			else:
				#nl => Recognize if reading field outside of namelist
				if nl:
					#(not '/' in l or '=' in l) => Recognize namelist field from otehr fields or namelist end '/'
					if not '/' in ls or '=' in ls:
						#if not nl: raise Exception( "Corrupted input file:\n{}".format( l))
						#Read namelist fields separated by endline ('\n') or by commas (',')
						for e in filter( None, ls.split( ",")):
							l1 = trim_ws(e).split( "=")
							if len(l1) != 2: raise Exception( "Corrupt input file:\n{}".format( l))
							v = l1[1].replace("\"", "").replace("'", "")
							
							#Check if the field/parameter name is present in the namelist template
							try: self.set_nl( nl=nl, k=l1[0], v=v)
							except NameError as e: print( e)
					else:
						nl = None
				#Case reading a card
				else:
					lt = ' '.join( filter( None, ls.replace( "{", "").replace( "}", "").split( " ")))
					lt = lt.split( " ")
					if self.check_card( lt[0]):
						try: v = lt[1]
						except: v = None
						card=lt[0]
						self.set_card( card=card, v=v)
						#print( card, v)
					else:
						self.set_card( card=card, el=ls)
		return






