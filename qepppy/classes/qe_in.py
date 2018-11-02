import sys


from .qe_templ import qe_templ as templ

templ
def trim_ws( str):
	"""
	Trim all withspace not included in a string
	"""
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

	trim_ws

class qe_in( templ):
	"""
	Class to handle any QE input (after loading the proper template).
	Supposed to be used as a parent class for child specific to the qe file.
	"""
	def __init__( self, fname="" , templ_file="", **kwargs):
		#A namelist template '_d' must be declared in the child class!!!!!!
		if templ_file:
			self.templ_file = templ_file
		if not self.templ_file: raise Exception( "Must give a template file.\n")
		self.load_templ( self.templ_file)

		if fname: self.parse_input( fname=fname)
		#Check if initialization keyword arguments are compliant with the given namelist template
		for nl, v in kwargs.items():
			if not isinstance( v, dict):
				raise Exception( "Invalid kwargs.\n{}".format( kwargs))
			for k, v1 in v.items():
				self.set_nl( nl, k, v1)

		return

	def __str__( self):
		return self.convert()

	def print_input( self, fname="", check=True):
		"""
		Print this or a child class __str__() to a file or to the stdout
		If check, validate the file before printing it
		"""
		if fname: f = open( fname, "w+")
		else: f = sys.stdout

		if check: self.validate()
		f.write( self.__str__())

		if fname: f.close()

		return


	def parse_input( self, fname=""):
		"""
		Read an input file and load it into the template values.
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
							self.set_nl( nl=nl, k=l1[0], v=v)
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






