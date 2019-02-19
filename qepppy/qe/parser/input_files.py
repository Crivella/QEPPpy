import sys
from .templates import template_parser as templ
from ...logger import logger, error

#import logging
#logger = logging.getLogger(__name__)
#logging.basicConfig(format='%(levelname)s: %(name)s\n%(message)s\n')

def trim_ws(str):
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

# logger()()
class input_files(templ):
	"""
	Class to handle any QE input (after loading the proper template).
	Supposed to be used as a parent class for child specific to the qe file.
	kwargs:
	 - parse = Name of the file to parse
	 - templ_name = Name of the template file to use
	"""
	def __init__(self, templ_file="", parse="", **kwargs):
		try:
			self.templ_file
		except AttributeError as e:
			self.templ_file = None
		if templ_file:
			self.templ_file = templ_file
		if not self.templ_file:
			raise error("Must give a template file.\n")
		self.load_templ(self.templ_file)

		#parse = kwargs.get('parse')
		if parse:
			self.parse_input(parse=parse)
		else:
			pass
			#logger.warning("Creating instance of 'qe_in' without parsing any input file.")

		#Check if initialization keyword arguments are compliant with the given namelist template
		for nl, v in kwargs.items():
			if not isinstance(v, dict):
				continue
				#raise Exception("Invalid kwargs.\n{}".format(kwargs))
			for k, v1 in v.items():
				self.set_nl(nl, k, v1)

		super().__init__(**kwargs)
		return

	def __str__(self):
		return self.convert()

	def print_input(self, out="", check=True):
		"""
		Print this or a child class __str__() to a file or to the stdout
		If check, validate the file before printing it
		"""
		if out:
			f = open(out, "w+")
		else:
			f = sys.stdout

		if check:
			self.validate()
		f.write(self.__str__())

		if out: 
			f.close()

		return


	def parse_input(self, parse):
		"""
		Read an input file and load it into the template values.
		"""
		with open(parse) as f:
			content = f.readlines()
		
		nl = None
		card = None
		for l in content:
			ls = l.strip()
			if not ls: 
				continue
			#CASE: Namelist name
			if '&' == ls[0]:
				nl = ls[1:].upper()
				self.check_nl(nl)
				continue
			if '/' == ls:
				if not nl:
					raise Exception("Terminating namelist before opening it ?!")
				nl = None
				continue

			if nl:
				ls = ls.split("!")[0].strip()
				if not ls:
					continue
				if not '=' in ls:
					raise Exception("Invalid syntax '{}'".format(ls))
					#if not nl: raise Exception("Corrupted input file:\n{}".format(l))
				#Read namelist fields separated by endline ('\n') or by commas (',')
				for e in filter(None, ls.split(",")):
					l1 = trim_ws(e).split("=")
					if len(l1) != 2: 
						raise Exception("Invalid syntax '{}'".format(ls))
					v = l1[1].replace("\"", "").replace("'", "")
					
					#Check if the field/parameter name is present in the namelist template
					self.set_nl(nl=nl, k=l1[0], v=v)
				continue
			#Case reading a card
			else:
				# lt = ' '.join(filter(None, ls.replace("{", "").replace("}", "").split(" ")))
				# lt = lt.split(" ")
				ls = ls.split("#")[0].strip()
				lt = list(filter(None, ls.replace("{", "").replace("}", "").split(" ")))
				if self.check_card(lt[0]):
					card=lt[0]
					try: 
						v = lt[1]
					except: 
						v = None
					self.set_card(card=card, v=v)
					#print(card, v)
				else:
					self.set_card(card=card, el=ls)
				continue

			# raise Exception("Invalid syntax '{}'".format(l))






