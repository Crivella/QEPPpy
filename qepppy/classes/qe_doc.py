
def gkw( d):
	if isinstance( d, dict):
		return d.get( 'keyword')
	return ""

class qe_doc_parser():
	def parse( self, fname=""):
		"""
		Parse a .def file from the Quantum ESPRESSO documentation.
		Produce a dictionary with the following structure:
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
				l: (List of card elements)
				syntax:{
					cond: "..." (Condition on card value)
					l:[ [{n: varname, t: TYPE}, ..., [{...}, ...]], [...], ([{...}, {...}, ..., ], s, e), ...]
						(Every element of the list represent a line
						A Tuple represent a repeating line
						A List within a list marks optional arguments)
				}
				syntax1:{...} (if multiple syntaxes are provided)
			}
			CARD_NAME2:{...}
			...
		}
		"""
		type_check = { "CHARACTER":str, "INTEGER":int, "LOGICAL":bool, "REAL":float}
		nl={
			'nl':[],
			'card':[]
			}

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

		#"""
		import json
		with open("report", "w") as f:
			f.write( json.dumps( intermediate, indent=2))
		#"""

		def _parse_nl_var_( k="", v={}, namelist=""):
			#Function to set a var in the final namelist parsing the temporary nested dict
			t = None #Handle variable type
			s = None #Handle array var start
			e = None #Handle array var end

			if "unnamed" in k: return #Skip unnamed dicitonaries
			if "info" in k: return #Skip info dicitonaries
			if not isinstance( v , dict): return #Check if element is a dictionary

			#Case vargroup: read all variable inside
			if "vargroup" in gkw(v):
				#print ("Vargroup found: ", namelist, k, v)
				if not isinstance( v, dict):
					raise Exception( "The keyword '{}' in namelist '{}' has not been parsed as a dict...\n".format( 
						k, namelist))
				ta = v.get( 'type', None)
				t = type_check.get( str( ta).upper())
				for k2, v2 in v.items():
					if v2 == "":
						nl[namelist][k2] ={
							'v':"",
							't':t,
							'd':None,
							'c':[],
							'vec':None
							}
				return

			#Case group: read all variable inside (recursive)
			if "group" ==gkw(v):
				#print ("Group found: ", namelist, k, v)
				for k2, v2 in v.items():
					_parse_nl_var_( namelist=namelist, k=k2, v=v2)
				return

			#Case var
			nl[namelist][k] = {
				'v':"",
				't':None,
				'd':None,
				'c':[],
				'vec':None
				}
			ptr = nl[namelist][k]
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
					raise Exception( "No top value for the array.\n")
				ptr['vec'] = (s,e)
				ptr['v'] = []
			#Get var type
			ta = v.get( 'type', None)
			t = type_check.get( str( ta).upper())
			ptr['t']=t

			#Get default value
			a = str( v.get( 'default', '')).replace("'", "").replace("D", "E")
			try: ptr['d'] = t( a)
			except: ptr['d'] = a
			#Check if status is set to REQUIRED
			if str( v.get( 'status')).strip().upper() == "REQUIRED":
				 ptr['v'] = "***"
			#Check if possible option list is present
			for k3, v3 in v.get( 'options', {}).items():
				if isinstance( v3, dict):
					if "info" in k3: continue
					ptr['c'].append( k3)


		def _parse_card_ ( name, card):
			nl[name]={
				'v':"",
				'c':[],
				'd':"",
				'r':True,
				'l':[]
			}
			ptr = nl[name]
			def _parse_group_( card):
				ta = card.get( 'type', None)
				t = type_check.get( str( ta).upper())
				l = []
				for k, v in card.items():
					if v == "":
						l.append( {'n':k, 't':t})
				return l
			def _parse_table_elements_( card):
				l=[]		
				for k1, v1 in card.items():
					if not isinstance( v1, dict): continue
					kw = gkw( v1)
					if kw== 'col' or kw == 'row':
						l.append( {'n':k1, 't':v1.get( 'type')})
					elif 'group' in kw:
						l += _parse_group_( v1)
					elif kw == "optional" or kw == "conditional":
						l += [ _parse_table_elements_( v1)]
					else: raise Exception( "Unexpected '{}' in _parse_table_elements_.\n".format( kw))
				return l

			def _parse_table_( card):
				l = []
				s = None
				e = None
				#print( "Parsing: ", card)
				for v in card.values():
					if not isinstance( v, dict): continue
					kw = gkw( v)
					if kw == 'rows' or kw == 'cols':
						e = v.get( 'end')
						s = v.get( 'start')
						l = _parse_table_elements_( v)
						break

				print( "Table: ", (l, s, e))
				return (l, s, e)

			def _parse_syntax_( card):
				if gkw( card) == 'syntax':
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
						if gkw( v) == 'line':
							new=[]
							for k1, v1 in v.items():
								if not isinstance( v1, dict): continue
								kw = gkw( v1)
								if kw == 'var':
									ta = v1.get( 'type', None)
									t = type_check.get( str( ta).upper())
									new.append( {'n':k1, 't':t})
									#ptr[aname]['l'].append( t)
								if 'group' in kw:
									new += _parse_group_( v1)
							ptr[aname]['l'].append( new)
						if gkw( v) == 'table':
							ptr[aname]['l'].append(  _parse_table_( v))

					print( "SYNTAX: ", name)
				else:
					for v1 in card.values():
						if isinstance( v1, dict):
							_parse_syntax_( v1)

			a = str( card.get( 'label'))
			ptr['r'] = not "optional card" in a.lower()
			for k, v in card.items():
				if isinstance( v, dict):
					if gkw( v) == "flag":
						ptr['c'] = v.get( 'enum').split( " | ")
						ptr['d'] = v.get( 'default')
			_parse_syntax_( card)


			return

		#Cicle through nested dict and parse it into the final dict
		#Cicle through all namelists and cards					
		for k1, v1 in intermediate['input_description'].items():
			if "unnamed" in k1: continue
			if isinstance( v1, dict):
				nl[gkw( v1).replace("namelist", "nl")].append( k1)
				#Cicle through all var/vargroup/group
				if gkw( v1) == "namelist":
					nl[k1]={}
					for k2, v2 in v1.items():
						_parse_nl_var_( namelist=k1, k=k2, v=v2)
				elif gkw( v1) == "card":
					_parse_card_( name=k1, card=v1)
			
		#"""
		from pprint import pprint as pp	
		with open("report2", "w") as f:
			pp( nl, stream=f, indent=2)
			#f.write( json.dumps( nl, indent=2))
		#"""

		self._templ_ = nl
		return

	def validate( self):
		"""
		Validate a namelist after a read call. Checks:
			-all REQUIRED var have been set
			-all var type are compliant with the real type
		Raise Exception if invalid
		Return True if valid
		"""
		check_mand = False
		err = ""
		for n in self._templ_['nl']:
			#Check for unused namelist (does not print it)
			for el, v in self._templ_[n].items():
				if v['v'] == "***":
					check_mand = True
					err += "ERROR: Mandatory input parameter {} in namelist {} not set.\n".format( el, n)
				elif v['v']:
					if v['c']:
						if not any( v == opt for opt in v['c']):
							err += \
								"Parameter '{}/{}' = '{}' is not within range of possible values: \n{}\n".format( 
									n, el, v['v'], v['c'])
		if check_mand:
			raise Exception( err)

		return True

	def _maxl_( self):
		longest = 0
		for n in self._templ_['nl']:
			for el, v in self._templ_[n].items():
				if v['v']:
					app = len( el)
					if v['vec']: app += 4
					if longest < app: longest = app
		longest += 2
		return longest

	def convert( self):
		def _format_( v):
			c = ""
			if isinstance( v, str):
				c += "'{0}'".format( v)
			if isinstance( v, bool):
				if v: c += ".TRUE."
				else: c += ".FALSE."
			else:
				if isinstance( v, (int, float)):
					c += "{}".format( v).replace( "e", "D").replace( "E", "D")
			return c

		nl = self._templ_['nl'].copy()
		for namelist in self._templ_['nl']:
			#Check for unused namelist (does not print it)
			if not any( v['v'] for v in self._templ_[namelist].values()):
				nl.pop( nl.index(namelist))

		#Write all the used namlists/parameters
		content = ""
		longest = self._maxl_()
		for namelist in nl:
			content += "&{}\n".format(namelist)
			for el, v in self._templ_[namelist].items():
				if v['v']:
					if v['vec']:
						for n, val in enumerate( v['v']):
							if not val: continue
							app = el + "({})".format( n+1)
							content += "{0:>{1}} = ".format( app, longest,)
							content += _format_( val)
							content += " ,\n"
					else:
						content += "{0:>{1}} = ".format( el, longest)
						content += _format_( v['v'])
						content += " ,\n"
			content += "/\n\n"
		return content

	def check_nl( self, nl):
		if nl in self._templ_['nl']: return True
		return False

	def get( self, nl, k):
		return self._templ_[nl][k]['v']

	def set( self, nl, k, v):
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
		#If array case
		if n:
			if not isinstance( ptr['v'], list):
				raise Exception( "'{}' from namelist '{}' is not an array variable.\n".format( k, nl))
			#Check if array index is within maximum possible value
			s = ptr['vec'][0]
			e = ptr['vec'][1]
			if s in self._templ_[nl]:
				if self._templ_[nl][s]['v']:
					s=self._templ_[nl][s]['v']
			if e in self._templ_[nl]:
				if self._templ_[nl][e]['v']:
					e=self._templ_[nl][e]['v']
			if isinstance( s, int) and isinstance( e, int):
				if not s <= n <= e: 
					raise Exception( "'{}({})' out of array range '{}-{}'".format( k, n, s, e))
			while len(ptr['v']) < n: ptr['v'].append( '')

		if t == float or t == int: 
			if isinstance( v, str):
				v=v.replace( "D", "e").replace( "d", "e")
		try: v = t( v)
		except: raise TypeError( "Parameter '{}'' must be of type '{}': value '{}' is invalid.\n".format( k, t, v))

		if n: ptr['v'][n-1] = v
		else: ptr['v'] = v
			
		#Check value agains possible values
		if ptr['c']:
			if not any( v == opt for opt in ptr['c']):
				raise Exception( 
					"Parameter '{}/{}' = '{}' is not within range of possible values: \n{}".format( 
						nl, k, v, ptr['c']))





