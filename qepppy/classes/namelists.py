import sys

def trim_ws( str):
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

def make_namelist( fname=""):
	"""
	Parse a .def file from the Quantum ESPRESSO documentation.
	Produce a dictionary with the following structure:
	{
	nl:   ['", ... (List of namelists name)]
	card: ['", ... (List of cards name)]
	NAMELIST_NAME: {
		VAR_NAME: {
			v: (Value of the parameter)
			t: (Type of the parameter)
			d: (Default value)
			c: (List of possible acceppted value for the parameter)
			vec:None/(start,end) (Info for array like variables)
			}
		...
		}
	}
	"""
	type_check = { "CHARACTER":str, "INTEGER":int, "LOGICAL":bool, "REAL":float}
	nl={
		'nl':[],
		'card':[]
		}

	with open(fname) as f:
		content = f.read()

	app = {}
	ptr=app
	
	ptr_list=[] #Stack of ptr storing the namelist lvl
	#Keywords(KWs) that trigger the parser
	#Every { found after a KW generate a new dictionary
	keywords = [ "", "info", "namelist", "card", "var", "dimension", "group", 
		"vargroup", "status", "options", "default", "opt", "input_description",]

	# { found if one of this flag are triggered are ignored (considered as text)
	cflag_l = [ 'info', 'default', 'opt', 'status', '']
	
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
		if c == "\n" and flag['vargroup']:
			if kw == "var":
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
		#Set cflag if any of the KW flags in cflag_l is true
		#If cflag is set, treat al subsequent text as string until } that reset the flag is met
		cflag = any( flag[cf] for cf in cflag_l)
		if not cflag:
			#String mode handler (if not already in cflag mode)
			if c == "'" and not str_c2: str_c1 = not str_c1; continue
			if c == '"' and not str_c1: str_c2 = not str_c2; continue
		#Whitespace/{}/# handler
		if c == " " or c == "{" or c == "}" or c == "#":
			#In cflag mode handle internal brackets as string
			if cflag:
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
				if l == "default" or l == "status": ptr[l]=' '.join( parse.split())
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
					#Handle and merge repeated dictionaries (eg: {..., group:{a,b}, group:{c,d}} => {..., group{a,b,c,d}})
					if name in ptr:
						if not isinstance( ptr, dict):
							ptr[name]={}
					#Create new dictionary
					else:
						ptr[name]={}
					ptr=ptr[name]
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
		f.write( json.dumps( app, indent=2))
	#"""

	def set_var( k="", v={}, namelist=""):
		#Function to set a var in the final namelist parsing the temporary nested dict
		t = None #Handle variable type
		s = None #Handle array var start/end
		e = None

		if "unnamed" in k: return #Skip unnamed dicitonaries
		if "info" in k: return #Skip info dicitonaries
		if not isinstance( v , dict): return #Check if element is a dictionary

		#Case vargroup: read all variable inside
		if "vargroup" in k:
			if not isinstance( v, dict):
				raise Exception( "The keyword '{}' in namelist '{}' has not been parsed as a dict...\n".format( 
					k, namelist))
			if "type" in v:
				if v['type'] in type_check:
					t = type_check[v['type']]
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
		if "group" in k:
			if not isinstance( v, dict):
				raise Exception( "The keyword '{}' in namelist '{}' has not been parsed as a dict...\n".format( 
					k, namelist))
				for k2, v2 in v.items():
					set_var( namelist=namelist, k=k2, v=v2)
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
		#Check if valid type and assigni it
		if "type" in v:
			t = v['type'].upper()
			if t in type_check:
				t = type_check[t]
				ptr['t'] = t
			else: raise Exception( "'{}': Unrecognized type '{}'.\n".format( k, v['type']))
		#Check if default value is present
		if "default" in v:
			if ptr['t']:
				if isinstance( v['default'], t):
					ptr['d'] = t( v['default'].replace("'", "").replace("D", "E"))
				else:
					ptr['d'] = str( v['default'].replace("'", "").replace("D", "E"))
		#Check if status is set to REQUIRED
		if "status" in v:
			if v['status'].strip().upper() == "REQUIRED":
				 ptr['v'] = "***"
		#Check if possible option list is present
		if "options" in v:
			for k3, v3 in v['options'].items():
				if isinstance( v3, dict):
					if "info" in k3: continue
					ptr['c'].append( k3)

	#Cicle through nested dict and parse it into the final dict
	#Cicle through all namelists and cards					
	for k1, v1 in app['input_description'].items():
		if "unnamed" in k1: continue
		if isinstance( v1, dict):
			nl[v1['keyword'].replace("namelist", "nl")].append( k1)
			nl[k1]={}
			#Cicle through all var/vargroup/group
			for k2, v2 in v1.items():
				set_var( namelist=k1, k=k2, v=v2)
		
	"""
	from pprint import pprint as pp	
	with open("report2", "w") as f:
		pp( nl, stream=f, indent=2)
		#f.write( json.dumps( nl, indent=2))
	#"""

	return nl

class namelist_handler():
	def __init__( self, **kwargs):
		#A namelist template '_d' must be declared in the child class!!!!!!
		if not '_d' in self.__dict__:
			raise Exception( "Must first associate the parent class with a namelist dictionary")
		#Check if initialization keyword arguments are compliant with the given namelist template
		for k, v in kwargs.items():
			check_kw = False
			for n in self._d['nl']:
				if k in self._d[n]:
					check_kw = True
					self._d[n][k]['v'] = v
			if not check_kw:
				raise Exception( "Invalid keyword argument '{}'.".format( k))

		return

	def __str__( self):
		content = ""

		longest = 0
		check_mand = False
		err = ""
		nl = self._d['nl'].copy()
		for namelist in self._d['nl']:
			#Check for unused namelist (does not print it)
			check_used = False
			for el, v in self._d[namelist].items():
				if v['v']:
					check_used = True
					app = len(el)
					if v['vec']: app += 4
					#Adjust lenght of written parameters to have the = all in column
					if longest < app: longest = app
					#Check for mandatory parameters
					if v['v'] == "***":
						check_mand = True
						err += "ERROR: Mandatory input parameter {} in namelist {} not set.\n".format( el, namelist)
			if not check_used:
				nl.pop( nl.index(namelist))
		if check_mand:
			raise Exception( err)
		longest += 2

		def _format( v):
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


		#Write all the used namlists/parameters
		for namelist in nl:
			content += "&{}\n".format(namelist)
			for el, v in self._d[namelist].items():
				if v['v']:
					if v['vec']:
						for n, val in enumerate( v['v']):
							app = el + "({})".format( n+1)
							content += "{0:>{1}} = ".format( app, longest,)
							content += _format( val)
					else:
						content += "{0:>{1}} = ".format( el, longest)
						content += _format( v['v'])
					content += " ,\n"
			content += "/\n\n"

		return content

	def _namelist_arg_set_( self, nl, k, v):
		"""
		Function to set a value in the namelist template
		"""
		n=None

		#If array case
		if '(' in k:
			n=int( k.split('(')[1].split(')')[0])
			k=str( k.split( '(')[0])

		#Check if k is present in preset namelist
		if not k in self._d[nl]: raise NameError( "Ignored unrecognized parameter '{}'\n".format( k))
		ptr = self._d[nl][k]
		t = ptr['t']
		#If array case
		if n:
			if not isinstance( ptr['v'], list):
				raise Exception( "'{}' from namelist '{}' is not an array variable.\n".format( k, nl))
			#Check if array index is within maximum possible value
			s = ptr['vec'][0]
			e = ptr['vec'][1]
			if s in self._d[nl]:
				if self._d[nl][s]['v']:
					s=self._d[nl][s]['v']
			if e in self._d[nl]:
				if self._d[nl][e]['v']:
					e=self._d[nl][e]['v']
			if isinstance( s, int) and isinstance( e, int):
				if not s <= n <= e: 
					raise Exception( "'{}({})' out of array range '{}-{}'".format( k, n, s, e))
			while len(ptr['v']) < n: ptr['v'].append( '')

		if t == float or t == int: v=v.replace( "D", "e").replace( "d", "e")
		try: v = t( v)
		except: raise TypeError( "Parameter '{}'' must be '{}': value '{}' is invalid.\n".format( k, t, v))

		if n: ptr['v'][n-1] = v
		else: ptr['v'] = v
			
		#Check value agains possible values
		if ptr['c']:
			if not any( v == opt for opt in ptr['c']):
				raise Exception( 
					"Parameter '{}/{}' = '{}' does not respect possible values {}.\n".format( 
						nl, k, v, ptr['c']))

		return

	def fprint( self, fname=""):
		"""
		Print this or a child class __str__() to a file or to the stdout
		"""
		if fname:
			f = open( fname, "w+")
		else:
			f = sys.stdout

		content = self.__str__()
		f.write( content)

		if fname:
			f.close()

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
		err = ""
		for l in content:
			#Ignore comments
			l = l.strip().split( "!")[0]
			if not l: continue
			#CASE: Namelist name
			if '&' == l[0]:
				nl = l[1:].upper()
				if not nl in self._d['nl']:
					raise Exception( "Reading unrecognized namelist '{}'".format( nl))
			#CASE other
			else:
				#(not '/' in l or '=' in l) => Recognize namelist field from otehr fields or namelist end '/'
				#nl => Recognize if reading field outside of namelist
				if (not '/' in l or '=' in l) and nl and l:
					if not nl: raise Exception( "Corrupted input file at line '{}'".format( l))
					#Read namelist fields separated by endline ('\n') or by commas (',')
					for e in filter( None, l.split( ",")):
						l1 = trim_ws(e).split( "=")
						v = l1[1].replace("\"", "").replace("'", "")
						
						#Check if the field/parameter name is present in the namelist template
						try: self._namelist_arg_set_( nl=nl, k=l1[0], v=v)
						except NameError as e: print( e)
				else:
					nl = None

		if err:
			raise Exception( err)
		return

pw_nl=make_namelist( fname="INPUT_PW.def")
"""
 pw_nl = {
	"nl" : ['CONTROL", "SYSTEM", "ELECTRONS", "IONS", "CELL'],
	"CONTROL":{
		"calculation" : "scf",
		"title" : "",
		"verbosity" : "low",
		"restart_mode" : "",
		"wf_collect":True,
		"nstep" : "",
		"iprint" : "",
		"tstress" : "",
		"tprnfor" : "",
		"dt" : "",
		"outdir" : "./tmp",
		"wfcdir" : "",
		"prefix" : "***",
		"lkpoint_dir" : "",
		"max_seconds" : "",
		"etot_conv_thr" : "",
		"forc_conv_thr" : "",
		"disk_io" : "",
		"pseudo_dir" : "***",
		"tefield" : "",
		"dipfield" : "",
		"lelfield" : "",
		"nberrycyc" : "",
		"lorbm" : "",
		"lberry" : "",
		"gdir" : "",
		"nppstr" : "",
		"lfcpopt" : "",
		"gate" : ""
	},
	"SYSTEM":{
		"ibrav" : "***",
		"celldm(1)" : "***",
		"celldm(2)" : "",
		"celldm(3)" : "",
		"celldm(4)" : "",
		"celldm(5)" : "",
		"celldm(6)" : "",
		"A" : "",
		"B" : "",
		"C" : "",
		"cosAB" : "",
		"cosAC" : "",
		"cosBC" : "",
		"nat" : "***",
		"ntyp" : "***",
		"nbnd" : "",
		"tot_charge" : "",
		"starting_charge" : "",
		"tot_magnetization" : "",
		"starting_magnetization" : "",
		"ecutwfc" : "***",
		"ecutrho" : "",
		"ecutfock" : "",
		"nr1" : "",
		"nr2" : "",
		"nr3" : "",
		"nr1s" : "",
		"nr2s" : "",
		"nr3s" : "",
		"nosym" : "",
		"nosym_evc" : "",
		"noinv" : "",
		"no_t_rev" : "",
		"force_symmorphic" : "",
		"use_all_frac" : "",
		"occupations" : "",
		"one_atom_occupations" : "",
		"starting_spin_angle" : "",
		"degauss" : "",
		"smearing" : "",
		"nspin" : "",
		"noncolin" : "",
		"ecfixed" : "",
		"qcutz" : "",
		"q2sigma" : "",
		"input_dft" : "",
		"exx_fraction" : "",
		"screening_parameter" : "",
		"exxdiv_treatment" : "",
		"x_gamma_extrapolation" : "",
		"ecutvcut" : "",
		"nqx1" : "",
		"nqx2" : "",
		"nqx3" : "",
		"lda_plus_u" : "",
		"lda_plus_u_kind" : "",
		"Hubbard_U" : "",
		"Hubbard_J0" : "",
		"Hubbard_alpha" : "",
		"Hubbard_beta" : "",
		"Hubbard_J" : "",
		"starting_ns_eigenvalue" : "",
		"U_projection_type" : "",
		"edir" : "",
		"emaxpos" : "",
		"eopreg" : "",
		"eamp" : "",
		"angle1" : "",
		"angle2" : "",
		"lforcet" : "",
		"constrained_magnetization" : "",
		"fixed_magnetization" : "",
		"lambda" : "",
		"report" : "",
		"lspinorb" : "",
		"assume_isolated" : "",
		"esm_bc" : "",
		"esm_w" : "",
		"esm_efield" : "",
		"esm_nfit" : "",
		"fcp_mu" : "",
		"vdw_corr" : "",
		"london" : "",
		"london_s6" : "",
		"london_c6" : "",
		"london_rvdw" : "",
		"london_rcut" : "",
		"dftd3_version" : "",
		"dftd3_threebody" : "",
		"ts_vdw_econv_thr" : "",
		"ts_vdw_isolated" : "",
		"xdm" : "",
		"xdm_a1" : "",
		"xdm_a2" : "",
		"space_group" : "",
		"uniqueb" : "",
		"origin_choice" : "",
		"rhombohedral" : "",
		"zgate" : "",
		"relaxz" : "",
		"block" : "",
		"block_1" : "",
		"block_2" : "",
		"block_height" : ""
	},
	"ELECTRONS":{
		"electron_maxstep" : "",
		"scf_must_converge" : "",
		"conv_thr" : "",
		"adaptive_thr" : "",
		"conv_thr_init" : "",
		"conv_thr_multi" : "",
		"mixing_mode" : "",
		"mixing_beta" : "",
		"mixing_ndim" : "",
		"mixing_fixed_ns" : "",
		"diagonalization" : "",
		"ortho_para" : "",
		"diago_thr_init" : "",
		"diago_cg_maxiter" : "",
		"diago_david_ndim" : "",
		"diago_full_acc" : "",
		"efield" : "",
		"efield_cart" : "",
		"efield_phase" : "",
		"startingpot" : "",
		"startingwfc" : "",
		"tqr" : ""
	},
	"IONS":{
		"ion_dynamics" : "",
		"ion_positions" : "",
		"pot_extrapolation" : "",
		"wfc_extrapolation" : "",
		"remove_rigid_rot" : "",
		"ion_temperature" : "",
		"tempw" : "",
		"tolp" : "",
		"delta_t" : "",
		"nraise" : "",
		"refold_pos" : "",
		"upscale" : "",
		"bfgs_ndim" : "",
		"trust_radius_max" : "",
		"trust_radius_min" : "",
		"trust_radius_ini" : "",
		"w_1" : "",
		"w_2" : ""
	},
	"CELL":{
		"cell_dynamics" : "",
		"press" : "",
		"wmass" : "",
		"cell_factor" : "",
		"press_conv_thr" : "",
		"cell_dofree" : ""
	}
}
"""






