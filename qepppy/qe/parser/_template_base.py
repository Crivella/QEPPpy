debug_templ = False #Enable convert of an empty template (Print all empty values. Put the max value for array var to 7 if not defined)
from ...logger import *


#import logging
#logger_msg = logging.getLogger( __name__)
#logging.basicConfig( format='%(levelname)s: %(name)s\n%(message)s\n')

@logger()
class templ_base( object):
	def convert( self):
		return ""

	def validate( self):
		try:
			return True and super().validate()
		except:
			return True

	@staticmethod
	def _check_type_( v, type):
		type_check = { "CHARACTER":str, "INTEGER":int, "LOGICAL":bool, "REAL":float, "STRING":str}
		ty = type_check.get( str( type).upper())
		if ty == float or ty == int: 
			if isinstance( v, str):
				v=v.replace( "D", "e").replace( "d", "e")
		try: val = ty( v)
		except:
			val = '***'
			warning.print( "Input: {}\nExpec: {}".format( v, ty))
		return val

	@staticmethod
	def _format_( v, f=0, a="'"):
		c = ""
		if isinstance( v, str):
			if f: a = "{2}{0:>{1}}{2}".format( v, f, a)
			else: a = "{1}{0}{1}".format( v, a)
			c += a
		if isinstance( v, bool):
			if v: c += ".TRUE."
			else: c += ".FALSE."
		else:
			if isinstance( v, (int, float)):
				if f: a = "{0:>{1}}".format( v, f)
				else: a = "{}".format( v)
				c += a.replace( "e", "D").replace( "E", "D")
		return c

	def _get_arr_ext_( self, sa, ea, n=-1):
		try: st= int( sa)
		except: st = self.find( sa)
		try: et= int( ea)
		except: et = self.find( ea)
		if not isinstance( et, int) and debug_templ: et = 7
		if any( not isinstance( a, int) for a in [st,et]):
			warning.print( "Failed to find boundary '{}-{}'.".format( sa, ea))
			return
			#raise Exception( "Failed to find boundary '{}-{}'.\n".format( sa, ea))
		if n >= 0:
			if not st <= n <= et:
				warning.print( "'{}' out of array range '{}-{}'".format( n, st, et))
				return
				#raise Exception( "'{}' out of array range '{}-{}'".format( n, st, et))
		return ( st, et)



@logger()
class namelist( templ_base):
	def check_nl( self, nl):
		"""
		Check if NAMELIST name actually exist in template
		"""
		if nl in self._templ_['nl']: return True
		return False

	def convert( self):
		"""
		Convert the filled template in a QE input file and return it.
		Return instance: str
		"""
		content = super().convert()
		nl = self._templ_['nl'].copy()
		#Check for unused namelist (does not print it)
		for namelist in self._templ_['nl']:
			if not any( v['v'] for v in self._templ_[namelist].values()) and not debug_templ:
				nl.pop( nl.index(namelist))
		#Write all the used namlists/parameters
		longest = self._maxl_()
		for namelist in nl:
			content += "&{}\n".format(namelist)
			for el, v in self._templ_[namelist].items():
				if v['v'] != None and v['v'] != '' or debug_templ:
					if v['vec']:
						if debug_templ:
							try: end = self._get_arr_ext_( v['vec'][0], v['vec'][1])[1]
							except: end = 7
							if not v['v']: v['v'] = ['']*end
						for n, val in enumerate( v['v']):
							if not val and not debug_templ: continue
							app = el + "({})".format( n+1)
							content += "{0:>{1}} = ".format( app, longest,)
							content += self._format_( val)
							content += " ,\n"
					else:
						content += "{0:>{1}} = ".format( el, longest)
						content += self._format_( v['v'])
						content += " ,\n"
			content += "/\n\n"
		return content

	def set_nl( self, nl, k, v):
		"""
		Set the value "v" for param "k" in the NAMELIST "nl"
		"""
		n=None
		#print( nl, k, v)
		#If array case
		if '(' in k:
			n=int( k.split('(')[1].split(')')[0])
			k=str( k.split( '(')[0])

		#Check if k is present in preset namelist
		#print ( self._templ_[nl])
		if not k in self._templ_[nl]: 
			warning.print( "Ignored unrecognized parameter '{}'".format( k))
			return
			#raise NameError( "Ignored unrecognized parameter '{}'\n".format( k))
		ptr = self._templ_[nl][k]
		v = self._check_type_( v, ptr['t'])

		#Check value agains possible values
		if ptr['c']:
			if not any( v == opt for opt in ptr['c']):
				warning.print( "Parameter '{}/{}' = '{}' is not within range of possible values: \n{}".format( nl, k, v, ptr['c']))
				#return
		#If array case
		if n:
			if not isinstance( ptr['v'], list):
				warning.print( "'{}' from namelist '{}' is not an array variable.".format( k, nl))
				return
				#raise Exception( "'{}' from namelist '{}' is not an array variable.\n".format( k, nl))
			while len(ptr['v']) < n: ptr['v'].append( '')
			ptr['v'][n-1] = v
		else: ptr['v'] = v

		return

	def validate( self):
		"""
		Validate namelist after a read call. Checks:
			-all REQUIRED var have been set
			-all var type are compliant with their type
		"""
		ret = True
		for nl in self._templ_['nl']:
			for el, v in self._templ_[nl].items():
				val = v['v']
				if val == "***":
					warning.print( "Required input parameter '{}' in namelist '{}' not set.\n".format( el, nl))
					ret = False
				elif val:
					if v['c']:
						if not any( val == opt for opt in v['c']):
							warning.print( "Parameter '{}/{}' = '{}' is not within range of possible values: \n{}\n".format( nl, el, val, v['c']))
							ret = False
					if v['vec']:
						n = len( val)
						if not self._get_arr_ext_( v['vec'][0], v['vec'][1], n): 
							warning.print( "Parameter: {}({})".format( el, n))
							ret = False
		return ret and super().validate()




	def _maxl_( self):
		#Get the longest element ot print among all namelists
		#Uset to align all the '=' in the printed QE input file
		longest = 0
		for n in self._templ_['nl']:
			for el, v in self._templ_[n].items():
				if v['v'] != None and v['v'] != '' or debug_templ:
					app = len( el)
					if v['vec']: app += 4
					if longest < app: longest = app
		longest += 2
		return longest




@logger()
class card( templ_base):
	def check_card( self, card):
		"""
		Check if a CARD name actually exist in the template
		"""
		if card in self._templ_['card']: return True
		return False

	def check_used( self, card):
		if self.check_card(  card):
			return self._templ_[card]['u']
		else:
			return False

	def convert( self):
		"""
		Convert the filled template in a QE input file and return it.
		Return instance: str
		"""
		content = super().convert()

		def _conv_syntax_( l, lvl=0, arr_lvl=0):
			#Convert the syntax dictionary of the internal template
			c = ""
			if not isinstance( l, list): return ""
			if lvl == 2:
				if not any( a['v'] for a in l): return ""
			for e in l:
				#print( e)
				if isinstance( e, dict):
					#print( e['n'], " =?= ", name)
					v = e['v']
					if isinstance( v, list):
						try: v = v[arr_lvl]
						except: v = 'empty1'
					c += " "
					c += self._format_( v, 10, a="")
				if isinstance( e, list): 
					c += _conv_syntax_( e, lvl+1, arr_lvl)
				if isinstance( e, tuple):
					ext = self._get_arr_ext_( e[1], e[2])
					num = ext[1]-ext[0]+1
					if e[3] == 'rows':
						for n in range( num):
							c +=  _conv_syntax_( e[0], lvl+1, n)
					elif e[3] == 'cols': 
						for r in e[0]:
							for n in range( num):
								try: v = r['v'][n]
								except: v = 'empty2'
								c += self._format_( v, 10, a="")
							c += "\n"
			if lvl == 1: c += "\n"
			return c
		#Convert the cards
		cards = self._templ_['card'].copy()
		#Check for unused cards (does not print it)
		for card in self._templ_['card']:
			if not self._templ_[card]['u'] and not debug_templ:
				cards.pop( cards.index(card))

		for card in cards:
			ptr = self._templ_[card]
			content += card
			if ptr['v']: content += " {{{}}}".format( ptr['v'])
			content += "\n"

			l = self._get_syntax_( ptr)
			content += _conv_syntax_( l)
			content += "\n"
		return content

	def set_card( self, card, v="", el=[]):
		"""
		Set the values inside a CARD.
		"v" sets the value supposed to be next the card name.
		"el" sets the card values line by line.
		"""
		ptr = self._templ_[card]
		synt= self._get_syntax_( ptr)
		ptr['u'] = True

		if v:
			if not v in ptr['c'] and ptr['c']:
				warning.print( "Invalid value '{}' for card '{}' ({}).\n".format( v, card, ptr['c']))
				return
			ptr['v'] = v
			return
		if el:
			lt = list( filter( None, el.strip().split( " ")))
			self._set_line_( lt, s=synt, card=card)
		return

	def validate( self):
		"""
		Validate card after a read call. Checks:
			-Card value is valid
			-all lines are compliant with the syntax
		Return True if valid
		"""
		ret = True
		def _validate_syntax_( l, lvl=0, arr=0):
			ret1 = True
			#print( "Validate: ", l, lvl, arr)
			if not isinstance( l, list): return True
			#Optional values can be either all set or unset
			if lvl == 2:
				#Check all set.... avoid because of lists :(
				#cat = not any( not a['v'] for a in l)
				#Check all unset
				cau = not any( a['v'] for a in l)
				if cau: return cau
			for e in l:
				if isinstance( e, dict):
					val = e['v']
					if not val: return False
					if arr:
						#print( val, arr)
						if len( val) != arr or any( a==None or a =='' for a in val):
							warning.print( "Param '{}: {}': Number of lines does not match specified value '{}'.".format( card, e['n'], arr))
							ret1 = False
				if isinstance( e, list): return _validate_syntax_( e, lvl+1, arr)
				if isinstance( e, tuple):
					ext = self._get_arr_ext_( e[1], e[2])
					return _validate_syntax_( e[0], lvl+1, ext[1]-ext[0]+1)

				"""
					if e['n'] == name:return e['v']
				if isinstance( e, tuple): return _sfind_( e[0])
				if isinstance( e, list): return _sfind_( e)
				"""
			return ret1

		for card in self._templ_['card']:
			c = self._templ_[card]
			if c['u']:
				v = c['v']
				if not v:
					if c['c']:
						for opt in c['c']:
							if opt in c['d']:
								v = opt
								break
						if not v:
							warning.print( "No option for card '{}' and no default value either.".format( card))
							ret = False

				l = self._get_syntax_( c)
				if not l:
					warning.print( "Cannot find a parsed syntax for card: '{} {}'.".format( card, v))
					ret = False
				if not _validate_syntax_( l):
					warning.print( "Syntax in card '{}' is invalid.".format( card))
					ret = False
			else:
				if c['r']:
					warning.print( "Mandatory card '{}' is not set.".format( card))
					ret = False
		return ret and super().validate()



	@staticmethod
	def _get_syntax_( card):
		v = card['v']
		for k1, v1 in card.items():
			if not isinstance( v1, dict): continue
			if not 'syntax' in k1: continue
			if isinstance( v1['cond'], str):
				if not v in v1['cond']: continue
			return v1['l']
		return None

	def _set_line_( self, line=[], s=[], col="", card=""):
		#Parse a line inside a card and set it on the first usnet syntax line
		if not isinstance( s, list): return
		#print( line)
		for e in s:
			if isinstance( e, dict):
				#Look for the first unset line in the syntax
				if e['v']: 
					if not isinstance( e['v'], list): return
				#Col case
				if col == 'cols':
					if e['v']: continue
					for v in line:
						val = self._check_type_( v, e['t'])

						e['v'].append( val)
					return True
				#Set all elements of the line
				app = line.copy()
				app.reverse()
				sprint = list( a['n'] for a in s if isinstance(a, dict))
				emsg = "Syntax error in card:\n{}.\nInput: {}\nExpec: {}"
				for el in s:
					if isinstance( el, dict):
						try: v = app.pop()
						except: #v = None
							warning.print( emsg.format( card, line, sprint))
							return False
					#Cicle through optional agruments
					if isinstance( el, list):
						if not app: break #No optional arguments present
						sprint += [list( a['n'] for a in el)] #Set for error msg
						for el1 in el:
							try: v = app.pop()
							except:
								warning.print( emsg.format( card, line, sprint))
								return False
							val = self._check_type_( v, el['t'])
							if isinstance( el1['v'], list): el1['v'].append( val)
							else: el1['v'] = val
					elif isinstance( el, dict):
						val = self._check_type_( v, el['t'])
						if isinstance( el['v'], list): el['v'].append( val)
						else: el['v'] = val
					else:
						warning.print( "Unrecognized syntax, expected dict or list.")
						return False
				if app:
					warning.print( emsg.format( card, line, sprint))
					return False
				return True
			if isinstance( e, list): 
				if self._set_line_( line, s=e, card=card): return
			if isinstance( e, tuple): 
				if self._set_line_( line, s=e[0], col=e[3], card=card): return
		return




