"""
Template format:
{
	nl:   ['", ... (List of namelists name)]
	card: ['", ... (List of cards name)]
	NAMELIST_NAME1: {
		VAR_NAME: {
			v: (Value of the parameter)
			t: (Type of the parameter)
			d: (Default value)
			c: (List of possible acceppted value for the parameter)
			vec:None/(start,end) (Info for array like variables)
			}
		...
		}
	NAMELIST_NAME2={...}
	...

	CARD_NAME1:{
		v: (Value associated with the card)
		c: (List of possible acceppted value for the card)
		d: (Default value for v)
		r: True/False (is card REQUIRED?)
		u: True/False (True if any values are set in syntax)
		syntax:{
			cond: "..." (Condition on card value)
			l:[ [{n: varname, v: value, t: TYPE}, ..., [{...}, ...]], [...], ([{...}, {...}, ..., ], s, e, kw), ...]
				(Every element of the list represent a line
				A Tuple represent a repeating line ( [line], start, end, keyword)
				A List within a list marks optional arguments)
		}
		syntax1:{...} (if multiple syntaxes are provided)
	}
	CARD_NAME2:{...}
	...
}
"""

debug = False #Enable convert of an empty template (Print all empty values. Put the max value for array var to 7 if not defined)

class namelist( object):
	pass

class card( object):
	pass

class qe_templ( namelist, card):
	def check_nl( self, nl):
		#Check if namelist exist in template
		if nl in self._templ_['nl']: return True
		return False

	def check_card( self, card):
		#Check if card exist in template
		if card in self._templ_['card']: return True
		return False

	def convert( self):
		#Convert the template into a string containing a QE input file
		content = ""
		#Convert the namelists
		nl = self._templ_['nl'].copy()
		#Check for unused namelist (does not print it)
		for namelist in self._templ_['nl']:
			if not any( v['v'] for v in self._templ_[namelist].values()) and not debug:
				nl.pop( nl.index(namelist))
		#Write all the used namlists/parameters
		longest = self._maxl_()
		for namelist in nl:
			content += "&{}\n".format(namelist)
			for el, v in self._templ_[namelist].items():
				if v['v'] != None and v['v'] != '' or debug:
					if v['vec']:
						if debug:
							try: end = self._get_arr_ext_( v['vec'][0], v['vec'][1])[1]
							except: end = 7
							if not v['v']: v['v'] = ['']*end
						for n, val in enumerate( v['v']):
							if not val and not debug: continue
							app = el + "({})".format( n+1)
							content += "{0:>{1}} = ".format( app, longest,)
							content += self._format_( val)
							content += " ,\n"
					else:
						content += "{0:>{1}} = ".format( el, longest)
						content += self._format_( v['v'])
						content += " ,\n"
			content += "/\n\n"

		#Convert the cards
		cards = self._templ_['card'].copy()
		#Check for unused cards (does not print it)
		for card in self._templ_['card']:
			if not self._templ_[card]['u'] and not debug:
				cards.pop( cards.index(card))

		for card in cards:
			ptr = self._templ_[card]
			content += card
			if ptr['v']: content += " {{{}}}".format( ptr['v'])
			content += "\n"

			l = self._get_syntax_( ptr)
			content += self._conv_syntax_( l)
			content += "\n"
		return content

	def find( self, *args, up=""):
		"""
		Find a var in all possible namelists and cards.
		Return its value if found otherwise return None
		"""
		def _syntax_find_( el, tof):
			#Recursive find to descend into syntax elements
			if not isinstance( el, list): return None
			for e in el:
				f = None
				if isinstance( e, dict):
					if e['n'] == tof: f = e['v']
				if isinstance( e, list):
					f = _syntax_find_( e, tof=tof)
				if isinstance( e, tuple): 
					f = _syntax_find_( e[0], tof=tof)
				if f != None: return f
			return None

		l=[]
		for name in args:
			n = None
			tof = name
			if "(" in name:
				tof = name.split( "(")[0]
				n = int( name.split( "(")[1].split( ")")[0])

			ret = None
			for nl in self._templ_['nl']:
				if up:
					if up != nl: continue
				f = self._templ_[nl].get( tof)
				if f: 
					if n: 
						try: ret = f['v'][n-1]
						except: pass
					else: ret = f['v']
					break
			for card in self._templ_['card']:
				if up:
					if up != card: continue
				synt = self._get_syntax_( self._templ_[card])
				f = _syntax_find_( synt, tof=tof)
				if f != None: 
					if n: 
						try: ret = f[n-1]
						except: pass
					else: ret = f
					break
			l.append( ret)

		if len( l) == 1: l = l[0]
		else: l = tuple( l)
		return l

	def load_templ( self, fname=""):
		import os
		if os.path.isfile( fname):
			with open( fname) as f:
				file = f.read()
		else:
			from pkg_resources import resource_string, resource_listdir
			#print( resource_listdir('qepppy.data', ''))
			if fname in resource_listdir('qepppy.data', ''):
				file = resource_string( 'qepppy.data', fname).decode('utf-8')

		import ast
		self._templ_ = ast.literal_eval( file)
		return

	def set_card( self, card, v="", el=[]):
		ptr = self._templ_[card]
		synt= self._get_syntax_( ptr)
		ptr['u'] = True

		if v:
			if not v in ptr['c'] and ptr['c']:
				raise Exception( "Invalid value '{}' for card '{}' ({}).\n".format( v, card, ptr['c']))
			ptr['v'] = v
			return
		if el:
			lt = list( filter( None, el.strip().split( " ")))
			self._set_line_( lt, s=synt, card=card)
		return


	def set_nl( self, nl, k, v):
		n=None
		#print( nl, k, v)
		#If array case
		if '(' in k:
			n=int( k.split('(')[1].split(')')[0])
			k=str( k.split( '(')[0])

		#Check if k is present in preset namelist
		#print ( self._templ_[nl])
		if not k in self._templ_[nl]: raise NameError( "Ignored unrecognized parameter '{}'\n".format( k))
		ptr = self._templ_[nl][k]
		v = self._check_type_( v, ptr['t'])
		#If array case
		if n:
			if not isinstance( ptr['v'], list):
				raise Exception( "'{}' from namelist '{}' is not an array variable.\n".format( k, nl))
			while len(ptr['v']) < n: ptr['v'].append( '')
			ptr['v'][n-1] = v
		else: ptr['v'] = v
			
		#Check value agains possible values
		if ptr['c']:
			if not any( v == opt for opt in ptr['c']):
				raise Exception( 
					"Parameter '{}/{}' = '{}' is not within range of possible values: \n{}".format( 
						nl, k, v, ptr['c']))
		return

	def validate( self):
		"""
		Validate all namelists and cards after a read call. Checks:
		Raise Exception if invalid
		Return True if valid
		"""
		"""
		from pprint import pprint as pp	
		with open("report3", "w") as f:
			pp( self._templ_, stream=f, indent=2)
			#f.write( json.dumps( nl, indent=2))
		#"""
		for n in self._templ_['nl']:
			self._validate_namelist_( n)
		for c in self._templ_['card']:
			self._validate_card_( c)
		return True








	#Internal non public methods
	def _check_type_( self, v, type):
		type_check = { "CHARACTER":str, "INTEGER":int, "LOGICAL":bool, "REAL":float, "STRING":str}
		ty = type_check.get( str( type).upper())
		if ty == float or ty == int: 
			if isinstance( v, str):
				v=v.replace( "D", "e").replace( "d", "e")
		try: val = ty( v)
		except: raise TypeError( "\ninput: {}\nexpec: {}\n".format( v, ty))
		return val

	def _conv_syntax_( self, l, lvl=0, arr_lvl=0):
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
				c += self._conv_syntax_( e, lvl+1, arr_lvl)
			if isinstance( e, tuple):
				ext = self._get_arr_ext_( e[1], e[2])
				num = ext[1]-ext[0]+1
				if e[3] == 'rows':
					for n in range( num):
						c +=  self._conv_syntax_( e[0], lvl+1, n)
				elif e[3] == 'cols': 
					for r in e[0]:
						for n in range( num):
							try: v = r['v'][n]
							except: v = 'empty2'
							c += self._format_( v, 10, a="")
						c += "\n"
		if lvl == 1: c += "\n"
		return c

	def _format_( self, v, f=0, a="'"):
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

	def _get_syntax_( self, card):
		v = card['v']
		for k1, v1 in card.items():
			if not isinstance( v1, dict): continue
			if not 'syntax' in k1: continue
			if isinstance( v1['cond'], str):
				if not v in v1['cond']: continue
			return v1['l']
		return None

	def _get_arr_ext_( self, sa, ea, n=-1):
		try: st= int( sa)
		except: st = self.find( sa)
		try: et= int( ea)
		except: et = self.find( ea)
		if not isinstance( et, int) and debug: et = 7
		if any( not isinstance( a, int) for a in [st,et]):
			raise Exception( "Failed to find boundary '{}-{}'.\n".format( sa, ea))
		if n >= 0:
			if not st <= n <= et:
				raise Exception( "'{}' out of array range '{}-{}'".format( n, st, et))
		return ( st, et)

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
				for el in s:
					if isinstance( el, dict):
						try: v = app.pop()
						except: #v = None
							raise Exception( "\nSyntax error in card:\n{}.\ninput: {}\nexpec: {}".format( 
								card, line, sprint))
					#Cicle through optional agruments
					if isinstance( el, list):
						if not app: break #No optional arguments present
						sprint += [list( a['n'] for a in el)] #Set for error msg
						for el1 in el:
							try: v = app.pop()
							except:
								raise Exception( "\nSyntax error in card:\n{}.\ninput: {}\nexpec: {}".format( 
									card, line, sprint))
							val = self._check_type_( v, el['t'])
							if isinstance( el1['v'], list): el1['v'].append( val)
							else: el1['v'] = val
					elif isinstance( el, dict):
						val = self._check_type_( v, el['t'])
						if isinstance( el['v'], list): el['v'].append( val)
						else: el['v'] = val
					else:
						raise Exception( "\nUnrecognized syntax, expected dict or list.\n")
				if app:
					raise Exception( "\nSyntax error in card:\n{}.\ninput: {}\nexpec: {}".format( 
						card, line, sprint))
				return True
			if isinstance( e, list): 
				if self._set_line_( line, s=e, card=card): return
			if isinstance( e, tuple): 
				if self._set_line_( line, s=e[0], col=e[3], card=card): return
		return

	def _validate_namelist_( self, nl):
		"""
		Validate namelist after a read call. Checks:
			-all REQUIRED var have been set
			-all var type are compliant with their type
		"""
		check_mand = False
		err = ""
		for el, v in self._templ_[nl].items():
			if v['v'] == "***":
				check_mand = True
				err += "ERROR: Mandatory input parameter '{}' in namelist '{}' not set.\n".format( el, nl)
			elif v['v']:
				if v['c']:
					if not any( v == opt for opt in v['c']):
						err += \
							"Parameter '{}/{}' = '{}' is not within range of possible values: \n{}\n".format( 
								nl, el, v['v'], v['c'])
				if v['vec']:
					n = len( v['v'])
					self._get_arr_ext_( v['vec'][0], v['vec'][1], n)
		if check_mand:
			raise Exception( err)

		return True

	def _validate_card_( self, card):
		"""
		Validate card after a read call. Checks:
			-Card value is valid
			-all lines are compliant with the syntax
		Return True if valid otherwise raise Exception
		"""
		def _validate_syntax_( l, lvl=0, arr=0):
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
							raise Exception(
								"Param '{}/{}': Number of lines does not match specified value '{}'.\n".format( 
									card, e['n'], arr))
				if isinstance( e, list): return _validate_syntax_( e, lvl+1, arr)
				if isinstance( e, tuple):
					ext = self._get_arr_ext_( e[1], e[2])
					return _validate_syntax_( e[0], lvl+1, ext[1]-ext[0]+1)

				"""
					if e['n'] == name:return e['v']
				if isinstance( e, tuple): return _sfind_( e[0])
				if isinstance( e, list): return _sfind_( e)
				"""
			return True

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
						raise Exception( "No option for card '{}' and no default value either.\n".format( card))

			l = self._get_syntax_( c)
			if _validate_syntax_( l): return True
			else: raise Exception( "Syntax in card '{}' is invalid.".format( card))
			raise Exception ( "Cannot find a parsed syntax for card: '{} {}'.\n".format( card, v))
		else:
			if c['r']:
				raise Exception( "Mandatory card '{}' is not set.\n".format( card))

		return True

	def _maxl_( self):
		#Get the longest element ot print among all namelists
		#Uset to align all the '=' in the printed QE input file
		longest = 0
		for n in self._templ_['nl']:
			for el, v in self._templ_[n].items():
				if v['v'] != None and v['v'] != '' or debug:
					app = len( el)
					if v['vec']: app += 4
					if longest < app: longest = app
		longest += 2
		return longest