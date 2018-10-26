"""
Importable a class capable of generating a namelist dict template by parsing a file
This class provides the following public methods:
	-parse(): Read the file and generate an internal namelist template
	-validate(): Check set values against the namelist template (eg: use after reading)
	-convert():  Convert the internal dict in a string QE input file
	-check_nl( nl="namelist"): Check if nl is valid (present in the internal namelist)
	-check_card( card=name): Check if card is valid (present in the internal cardlist)
	-set_nl:( nl="namelist", k="param", v="value to set") Set a namelist value in the namelist template
	-set_card: ( card="", v="", el=[]) Set a card value in the namelist template
		if v is set, set the card main value
		if el is set set a card list value
			el must be an entire lie to parse
	-get:( nl="namelist", k="param") Retrieve the value of a namelist parameter
	-find(name): Find a variable with name=name in the namelist/card template
"""

"""
Parse a .def file from the Quantum ESPRESSO documentation.
Produce a  template dictionary with the following structure:
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
class qe_doc_parser():
	type_check = { "CHARACTER":str, "INTEGER":int, "LOGICAL":bool, "REAL":float, "STRING":str}
	
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
			if not any( v['v'] for v in self._templ_[namelist].values()):
				nl.pop( nl.index(namelist))
		#Write all the used namlists/parameters
		longest = self._maxl_()
		for namelist in nl:
			content += "&{}\n".format(namelist)
			for el, v in self._templ_[namelist].items():
				if v['v'] != None and v['v'] != '':
					if v['vec']:
						for n, val in enumerate( v['v']):
							if not val: continue
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
			if not self._templ_[card]['u']:
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

	def find( self, name):
		"""
		Find a var in all possible namelists and cards.
		Return its value if found otherwise return None
		"""
		n = None
		tof = name
		if "(" in name:
			tof = name.split( "(")[0]
			n = int( name.split( "(")[1].split( ")")[0])

		for nl in self._templ_['nl']:
			f = self._templ_[nl].get( tof)
			if f: 
				if n: 
					try: ret = f['v'][n-1]
					except: ret = None
				else: ret = f['v']
				return ret
		for card in self._templ_['card']:
			for k, v in self._templ_[card].items():
				if not 'syntax' in k: continue
				f = self._syntax_find_( v['l'], tof=tof)
				if f != None: 
					if n: 
						try: ret = f[n-1]
						except: ret = None
					else: ret = f
					return ret
		return None

	def get( self, nl, k):
		#Get a value fron the internal namelist
		ptr = self._templ_.get( nl)
		ptr = ptr.get( k)
		v = ptr.get( 'v')
		return v

	def parse( self, fname=""):
		"""
		Parse a .def file from the Quantum ESPRESSO documentation.
		Produce a dictionary with the following structure: see above
		"""
		self._templ_={
			'nl':[],
			'card':[]
			}
		intermediate = self._intermediate_parse_( fname)	
		#Cicle through nested dict and parse it into the final dict template
		#Cicle through all namelists and cards					
		for k1, v1 in intermediate['input_description'].items():
			if "unnamed" in k1: continue
			if isinstance( v1, dict):
				self._templ_[self.gkw( v1).replace("namelist", "nl")].append( k1)
				if self.gkw( v1) == "namelist":
					self._templ_[k1]={}
					#Cicle through all var/vargroup/group
					for k2, v2 in v1.items():
						self._parse_nl_var_( namelist=k1, k=k2, v=v2)
				elif self.gkw( v1) == "card":
					self._parse_card_( name=k1, card=v1)
			
		"""
		from pprint import pprint as pp	
		with open("report2", "w") as f:
			pp( self._templ_, stream=f, indent=2)
			#f.write( json.dumps( nl, indent=2))
		#"""
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
		t = ptr['t']
		if t == float or t == int: 
			if isinstance( v, str):
				v=v.replace( "D", "e").replace( "d", "e")
		try: v = t( v)
		except: raise TypeError( "Parameter '{}'' must be of type '{}': value '{}' is invalid.\n".format( k, t, v))
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
		#if not isinstance( et, int): et = 11
		if any( not isinstance( a, int) for a in [st,et]):
			raise Exception( "Failed to find boundary '{}-{}'.\n".format( sa, ea))
		if n >= 0:
			if not st <= n <= et:
				raise Exception( "'{}' out of array range '{}-{}'".format( n, st, et))
		return ( st, et)

	def gkw( self, d):
		if isinstance( d, dict):
			return d.get( 'keyword')
		return ""

	def _intermediate_parse_( self, fname=""):
		with open(fname) as f:
			content = f.read()

		intermediate = {}
		ptr=intermediate
		
		ptr_list=[] #Stack of ptr storing the namelist lvl
		#Keywords(KWs) that trigger the parser
		#Every { found after a KW generate a new dictionary
		nl_keywords = [ '', 'info', 'namelist', 'var', 'dimension', 'group', 
			'vargroup', 'status', 'options', 'default', 'opt', 'input_description',]
		cd_keywords = [ 'card', 'syntax', 'table', 'rows', 'colgroup', 'optional', 'flag', 'enum', 'label',
			'choose', 'when', 'otherwise', 'rowgroup', 'row', 'col', 'cols', 'elsewhen', 'line', 'conditional']
		keywords = nl_keywords + cd_keywords

		#endline as interrupt:
		eflag_l = [ 'var', 'row', 'col']
		#overwrite mode: overwrite entire dictionary associated to this keyword with a string
		oflag_l = [ 'default', 'status', 'enum', 'label']
		# { found if one of this flag are triggered are ignored (considered as text)
		#string-like mode flag
		sflag_l = [ 'info', 'opt', '']
		sflag_l += oflag_l

		gflag_l= [ 'vargroup', 'colgroup', 'rowgroup']
		
		last=[] #Stack of keyword that where opened by {
		
		#Initialize flags from keyword
		flag = {}
		for e in keywords:
			flag[e] = False
		parse=""; kw = ""; name = ""

		#Argument mode initialized by -mname -mvalue
		#mlist is also used to print the generating keyword inside a dictionary
		check_m = False; mname=""; mvalue="";	mlist=[];	mbcheck=False

		str_c1=False;	str_c2=False #String mode initialized by ' or "
		comm = False #Mark seection as comment
		lvl = 0 #Handle paranthesis inside text

		#Loop over file letter by letter
		#Generate a nested dictionary from the .def file to be further parsed in the final dictionary
		for c in content:
			#If endline in vargroup (every line is a var without {})
			if c == "\n" and any( flag[gf] for gf in gflag_l):
				if kw in eflag_l:
					#Check against nameless var
					if not name:
						if not parse: raise Exception( "Corrupt .def file\n")
						name = parse
					#Create empty var in vargroup dictionary
					ptr[name]=""
					kw = ""; name = ""; parse = ""
			#If endline switch out of comment mode and skip  the rest
			elif c == "\n": comm = False; continue
			#If comment mode skip te rest
			if comm: continue
			#If string mode read c into parse and skip the rest
			if (str_c1 or str_c2) and c!="'" and c!='"': parse += c; continue
			#Set sflag if any of the KW flags in sflag_l is true
			#If sflag is set, treat al subsequent text as string until } that reset the flag is met
			sflag = any( flag[cf] for cf in sflag_l)
			if not sflag:
				#String mode handler (if not already in sflag mode)
				if c == "'" and not str_c2: str_c1 = not str_c1; continue
				if c == '"' and not str_c1: str_c2 = not str_c2; continue
			#Whitespace/{}/# handler
			if c == " " or c == "\t" or c == "{" or c == "}" or c == "#":
				#In sflag mode handle internal brackets as string
				if sflag:
					if c=='{': lvl += 1
					if c=='}': lvl -= 1
					if lvl >= 0: parse+=c; continue
					else: lvl = 0
				# # enables comment mode
				if c == "#": comm = True
				#If -arg mode, read arg name and arg value, store it as a tuple in mlist
				if check_m:
					#If arg name has not been read set it
					if not mname:
						if c != " ": raise Exception( "Corrupt .def file\n")
						mname = parse
					#If arg name is set, read the arg value
					else:
						mvalue+=parse
						if c=="{":
							mbcheck=True
						#Treat whitespace as break if not in { mode (mbcheck)
						if c == " " and mbcheck:
							mvalue += " "
						#Exit from -arg mode
						if c == "}" or (c == " " and mvalue and not mbcheck):
							if kw == "opt":
								name = mvalue
							else:
								mlist.append((mname,mvalue))
							mvalue = ""; mname = ""; mbcheck = False; check_m = False
					parse = ""
					continue
				#Set name associated to keyword
				if kw and not name: name = parse
				#Close dictionary and go back one level
				if c == "}":
					l = last.pop()
					flag[l] = False
					ptr = ptr_list.pop()
					#If KW is default or status, delete the dictionary and save the internal string
					if l in oflag_l: ptr[l]=' '.join( parse.split())#; print( last, l, " overwriteing ", parse)
					mname = ""
				#Check if KW and set it
				if parse in keywords and not kw: kw = parse
				#Create new dictionary associated to KW
				if c == "{":
					flag[kw] = True
					last.append( kw)
					#Handle unnamed dictionaries
					if not name:
						if kw: 
							name = kw
						else:
							u = 0
							while True:
								u += 1
								name = "unnamed{}".format( u)
								if not name in ptr: break
					ptr_list.append( ptr)
					if name:
						#Handle and merge repeated dictionaries by adding a number
						aname=name
						num = 0
						while aname in ptr:
							num += 1
							aname = name + str( num)

						ptr[aname]={}
						ptr=ptr[aname]
						#Store -arg and KW inside new dictionary
						mlist.append(("keyword", kw))
						while mlist:
							a = mlist.pop()
							ptr[a[0]] = a[1]

					name =  ""; kw = ""
				parse = ""
				continue
			#Enter -arg mode
			if c == "-" and kw: check_m = True; continue

			parse += c

		"""
		import json
		with open("report", "w") as f:
			f.write( json.dumps( intermediate, indent=2))
		#"""
		return intermediate

	def _parse_card_ ( self, name, card):
		self._templ_[name]={
			'v':"",
			'c':[],
			'd':"",
			'r':True,
			'u':False
		}
		ptr = self._templ_[name]
		def _parse_group_( card, tab=False):
			#Internal function to parse group elements
			ta = card.get( 'type', None)
			t = self.type_check.get( str( ta).upper())
			l = []
			for k, v in card.items():
				if v == "":
					if tab: l.append( {'n':k, 'v':[], 't':t})
					else: l.append( {'n':k, 'v':'', 't':t})
			return l
		def _parse_table_elements_( card):
			#Internal function to parse the elements of a table
			l=[]		
			for k1, v1 in card.items():
				if not isinstance( v1, dict): continue
				kw = self.gkw( v1)
				if kw== 'col' or kw == 'row':
					ta = v1.get( 'type')
					t = self.type_check.get( str( ta).upper())
					l.append( {'n':k1, 'v':[], 't':t})
				elif 'group' in kw:
					l += _parse_group_( v1, tab=True)
				elif kw == "optional" or kw == "conditional":
					l += [ _parse_table_elements_( v1)]
				else: raise Exception( "Unexpected '{}' in _parse_table_elements_.\n".format( kw))
			return l
		def _parse_table_( card):
			#Internal function to parse the table elements in syntax
			l = []
			s = None
			e = None
			#print( "Parsing: ", card)
			for v in card.values():
				if not isinstance( v, dict): continue
				kw = self.gkw( v)
				if kw == 'rows' or kw == 'cols':
					e = v.get( 'end')
					s = v.get( 'start')
					k = kw
					l = _parse_table_elements_( v)
					break
			#print( "Table: ", (l, s, e))
			return (l, s, e, k)
		def _parse_syntax_( card):
			"""
			Internal function to parse all syntax elements in a card
			Made to check recursively all the subelements of the card
			"""
			if self.gkw( card) == 'syntax':
				aname="syntax"
				num = 0
				while aname in ptr:
					num += 1
					aname = "syntax" + str(num)
				ptr[aname]={
					'cond':"",
					'l':[],
				}
				ptr[aname]['cond'] = card.get( 'flag')
				for k, v in card.items():
					if not isinstance( v, dict): continue
					if self.gkw( v) == 'line':
						new=[]
						for k1, v1 in v.items():
							if not isinstance( v1, dict): continue
							kw = self.gkw( v1)
							if kw == 'var':
								ta = v1.get( 'type', None)
								t = self.type_check.get( str( ta).upper())
								new.append( {'n':k1, 'v':'', 't':t})
								#ptr[aname]['l'].append( t)
							if 'group' in kw:
								new += _parse_group_( v1)
						ptr[aname]['l'].append( new)
					if self.gkw( v) == 'table':
						ptr[aname]['l'].append(  _parse_table_( v))
				#print( "SYNTAX: ", name)
			else:
				for v1 in card.values():
					if isinstance( v1, dict):
						_parse_syntax_( v1)

		a = str( card.get( 'label'))
		ptr['r'] = not "optional card" in a.lower()
		for k, v in card.items():
			if isinstance( v, dict):
				if self.gkw( v) == "flag":
					ptr['c'] = v.get( 'enum').split( " | ")
					ptr['d'] = v.get( 'default')
		_parse_syntax_( card)
		return

	def _parse_nl_var_( self, namelist="", k="", v={}):
		#Function to set a var in the final namelist parsing the temporary nested dict
		t = None #Handle variable type
		s = None #Handle array var start
		e = None #Handle array var end

		if "unnamed" in k: return #Skip unnamed dicitonaries
		if "info" in k: return #Skip info dicitonaries
		if not isinstance( v , dict): return #Check if element is a dictionary

		kw = self.gkw(v)
		#Case vargroup: read all variable inside
		if "vargroup" in kw:
			#print ("Vargroup found: ", namelist, k, v)
			if not isinstance( v, dict):
				raise Exception( "The keyword '{}' in namelist '{}' has not been parsed as a dict...\n".format( 
					k, namelist))
			ta = v.get( 'type', None)
			t = self.type_check.get( str( ta).upper())
			for k2, v2 in v.items():
				if v2 == "":
					self._templ_[namelist][k2] ={
						'v':"",
						't':t,
						'd':None,
						'c':[],
						'vec':None
						}
			return

		#Case group: read all variable inside (recursive)
		if "group" == kw:
			#print ("Group found: ", namelist, k, v)
			for k2, v2 in v.items():
				self._parse_nl_var_( namelist=namelist, k=k2, v=v2)
			return

		if "choose" == kw:
			for k2, v2 in v.items():
				if self.gkw( v2) in ['when','elsewhen','otherwise']:
					for k3, v3 in v2.items(): self._parse_nl_var_( namelist=namelist, k=k3, v=v3)
			return

		#Case var
		self._templ_[namelist][k] = {
			'v':"",
			't':None,
			'd':None,
			'c':[],
			'vec':None
			}
		ptr = self._templ_[namelist][k]
		#Check if var is array
		if "start" in v:
			s = v['start']
			try: s = int( s)
			except : s = str( s)
			if "end" in v:
				e = v['end']
				try: e = int( e)
				except : e = str( e)
			else:
				raise Exception( "No END value for the array.\n")
			ptr['vec'] = (s,e)
			ptr['v'] = []
		#Get var type
		ta = v.get( 'type', None)
		t = self.type_check.get( str( ta).upper())
		ptr['t']=t

		#Get default value
		a = str( v.get( 'default', '')).replace("'", "").replace("D", "E")
		try: ptr['d'] = t( a)
		except: ptr['d'] = a
		#Check if status is set to REQUIRED/MANDATORY
		reqstat =str( v.get( 'status')).strip().upper()
		if reqstat == "REQUIRED" or reqstat == "MANDATORY":
			 ptr['v'] = "***"
		#Check if possible option list is present
		for k3, v3 in v.get( 'options', {}).items():
			if isinstance( v3, dict):
				if "info" in k3: continue
				ptr['c'].append( k3)

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
						try: val = e['t']( v)
						except: raise Exception( "\nSyntax error in card:\n{}: {}.\ninput: {}\nexpec: {}".format( 
							card, e['n'], v, e['t']))
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
							except: #v = None
								raise Exception( "\nSyntax error in card:\n{}.\ninput: {}\nexpec: {}".format( 
									card, line, sprint))
							try: val = el1['t']( v)
							except: raise TypeError( "'{}' is not of type '{}'.\n".format( v, el1['t']))
							if isinstance( el1['v'], list): el1['v'].append( val)
							else: el1['v'] = val
					elif isinstance( el, dict):
						val = None
						try: val = el['t']( v)
						except: raise TypeError( "'{}' is not of type '{}'.\n".format( v, el['t']))
						#print ( el)
						if isinstance( el['v'], list): el['v'].append( val)
						else: el['v'] = val
						#print ( el)
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

	def _syntax_find_( self, el, tof):
		#Recursive find to descend into syntax elements
		if not isinstance( el, list): return None
		for e in el:
			f = None
			if isinstance( e, dict):
				if e['n'] == tof: f = e['v']
			if isinstance( e, list):
				f = self._syntax_find_( e, tof=tof)
			if isinstance( e, tuple): 
				f = self._syntax_find_( e[0], tof=tof)
			if f != None: return f
		return None


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
				if v['v']:
					app = len( el)
					if v['vec']: app += 4
					if longest < app: longest = app
		longest += 2
		return longest





