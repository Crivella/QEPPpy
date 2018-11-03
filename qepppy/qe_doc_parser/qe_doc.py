class qe_doc_parser():
	"""
	Instance used to parse a .def file from the Quantum ESPRESSO documentation.
	Produce a  template dictionary with the following structure:
	{ (Template dictionary)
	    nl:   ["...", ... (List of namelists name)]
	    card: ["...", ... (List of cards name)]
	    NAMELIST_NAME1: { (dictionary containing all the namelist parameters)
	        VAR_NAME: { (dictionary containing the details of the parameter)
	            v: (Value of the parameter)
	            t: (Type of the parameter as a string)
	            d: (Default value)
	            c: (List of possible accepted value for the parameter)
	            vec:None/(start,end) (Info for array like variables e.g. celldm(1/2/3/4/5/6))
	            }
	        ...
	        }
	    NAMELIST_NAME2={...}
	    ...
	    CARD_NAME1:{ (Dictionary containing the data and syntax of a card)
	        v: (Value associated with the card)
	        c: (List of possible accepted value for the card)
	        d: (Default value for v)
	        r: True/False (is card REQUIRED?)
	        u: True/False (True if any values are set in syntax)
	        syntax:{ (Dictionary defining the syntax that the card should follow)
	            cond: "..." (Condition on card value)
	            l:[ [{n: varname, v: value, t: TYPE}, ..., [{...}, ...]], [...], ([{...}, {...}, ..., ], s, e, kw), ...]
	                Every element of the list represent a line
	                A Tuple represent a repeating line ( [line], start, end, keyword)
	                A List within a list marks optional arguments
	        }
	        syntax1:{...} (if multiple syntaxes are provided)
	    }
	    CARD_NAME2:{...}
	    ...
	}
	"""
	type_check = { "CHARACTER":str, "INTEGER":int, "LOGICAL":bool, "REAL":float, "STRING":str}

	def __init__( self, fname="", out=""):
		if fname:
			self.parse( fname, out)
		else:
			import os
			for e in list(os.walk("."))[0][2]:
				if ".def" in e:
					print( e)
					self.parse( e, "{}.json".format( e[:-4]))
					#input()
		return

	def parse( self, fname="", out="templ"):
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
				self._templ_[self._gkw_( v1).replace("namelist", "nl")].append( k1)
				if self._gkw_( v1) == "namelist":
					self._templ_[k1]={}
					#Cicle through all var/vargroup/group
					for k2, v2 in v1.items():
						self._parse_nl_var_( namelist=k1, k=k2, v=v2)
				elif self._gkw_( v1) == "card":
					self._parse_card_( name=k1, card=v1)
			
		#"""
		from pprint import pprint as pp
		#import json
		with open( out, "w") as f:
			pp( self._templ_, stream=f, indent=2)
			#f.write( json.dumps( self._templ_, indent=2))
		#"""
		return








	#Internal non public methods
	def _gkw_( self, d):
		if isinstance( d, dict):
			return d.get( 'keyword')
		return ""

	def _intermediate_parse_( self, fname=""):
		"""
		Parser for the .def file that produced an intermediate dictionary to be
		further processed to get to the final template structure.
		"""
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
					if not mlist:
						ptr[name]=""
					else:
						ptr[name]={}
						mlist.append(("keyword", kw))
						while mlist:
							a = mlist.pop()
							ptr[name][a[0]] = a[1]

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
			"""
			Internal function to parse group elements
			"""
			ta = card.get( 'type', None)
			#t = self.type_check.get( str( ta).upper())
			l = []
			for k, v in card.items():
				if v == "":
					if tab: l.append( {'n':k, 'v':[], 't':ta})
					else: l.append( {'n':k, 'v':'', 't':ta})
			return l
		def _parse_table_elements_( card):
			#Internal function to parse the elements of a table
			l=[]		
			for k1, v1 in card.items():
				if not isinstance( v1, dict): continue
				kw = self._gkw_( v1)
				if kw== 'col' or kw == 'row':
					ta = v1.get( 'type', None)
					#t = self.type_check.get( str( ta).upper())
					l.append( {'n':k1, 'v':[], 't':ta})
				elif 'group' in kw:
					l += _parse_group_( v1, tab=True)
				elif kw == "optional" or kw == "conditional":
					l += [ _parse_table_elements_( v1)]
				else: raise Exception( "Unexpected '{}' in _parse_table_elements_.\n".format( kw))
			return l
		def _parse_table_( card):
			"""
			Internal function to parse the table elements in syntax
			"""
			l = []
			s = None
			e = None
			#print( "Parsing: ", card)
			for v in card.values():
				if not isinstance( v, dict): continue
				kw = self._gkw_( v)
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
			if self._gkw_( card) == 'syntax':
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
					if self._gkw_( v) == 'line':
						new=[]
						for k1, v1 in v.items():
							if not isinstance( v1, dict): continue
							kw = self._gkw_( v1)
							if kw == 'var':
								ta = v1.get( 'type', None)
								#t = self.type_check.get( str( ta).upper())
								new.append( {'n':k1, 'v':'', 't':ta})
								#ptr[aname]['l'].append( t)
							if 'group' in kw:
								new += _parse_group_( v1)
						ptr[aname]['l'].append( new)
					if self._gkw_( v) == 'table':
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
				if self._gkw_( v) == "flag":
					ptr['c'] = v.get( 'enum').split( " | ")
					ptr['d'] = v.get( 'default')
		_parse_syntax_( card)
		return

	def _parse_nl_var_( self, namelist="", k="", v={}):
		"""
		Function to set a var in the final namelist parsing the temporary nested dict
		"""
		t = None #Handle variable type
		s = None #Handle array var start
		e = None #Handle array var end

		if "unnamed" in k: return #Skip unnamed dicitonaries
		if "info" in k: return #Skip info dicitonaries
		if not isinstance( v , dict): return #Check if element is a dictionary

		kw = self._gkw_(v)
		#Case vargroup: read all variable inside
		if "vargroup" in kw:
			#print ("Vargroup found: ", namelist, k, v)
			if not isinstance( v, dict):
				raise Exception( "The keyword '{}' in namelist '{}' has not been parsed as a dict...\n".format( 
					k, namelist))
			ta = v.get( 'type', None)
			#t = self.type_check.get( str( ta).upper())
			for k2, v2 in v.items():
				if v2 == "":
					self._templ_[namelist][k2] ={
						'v':"",
						't':ta,
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
				if self._gkw_( v2) in ['when','elsewhen','otherwise']:
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
		ptr['t']=ta

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







