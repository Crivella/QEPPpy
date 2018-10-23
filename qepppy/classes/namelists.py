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
					self._d[n][k] = v
			if not check_kw:
				raise Exception( "Invalid keyword argument '{}'.".format( k))

		return

	def __str__( self):
		content = ""

		longest = 0
		check_mand = False
		err = ""
		nl = self._d["nl"].copy()
		for namelist in self._d["nl"]:
			#Check for unused namelist (does not print it)
			check_used = False
			for el, v in self._d[namelist].items():
				if v != "":
					check_used = True
					app = len(el)
					#Adjust lenght of written parameters to have the = all in column
					if longest < app: longest = app
					#Check for mandatory parameters
					if v == "***":
						check_mand = True
						err += "ERROR: Mandatory input parameter {} in namelist {} not set.\n".format( el, namelist)
			if not check_used:
				nl.pop( nl.index(namelist))
		if check_mand:
			raise Exception( err)
		longest += 2

		#Write all the used namlists/parameters
		for namelist in nl:
			content += "&{}\n".format(namelist)
			for el, v in self._d[namelist].items():
				if v != "":
					content += "{0:>{1}} = ".format( el, longest)
					if isinstance( v, str):
						content += "'{0}'".format( v)
					if isinstance( v, bool):
						if v: content += ".TRUE."
						else: content += ".FALSE."
					else:
						if isinstance( v, (int, float)):
							content += "{}".format( v).replace( "e", "D").replace( "E", "D")
					content += " ,\n"
			content += "/\n\n"

		return content

	def namelist_arg_set( self, nl, k, v):
		"""
		Function to set a value in the namelist template
		"""
		if nl in self._d:
			if k in self._d[nl]:
				self._d[nl][k] = v
				return 

		raise Exception( "Parameter '{}' not found in namelist '{}'.".format( k, nl))

	def print( self, fname=""):
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
			#CASE: Namelist name
			if '&' in l and l[1:].upper() in self._d['nl']:
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
						#Is the field value an integer?
						try:
							v = int(v)
						except:
							#Is the field value a float?
							#Also replace the fortran exponential notation 'D/d' with the python/c notation 'E/e'
							try:
								v = float(v.replace( "D", "e").replace( "d", "e"))
							except:
								#Is the field value a boolean?
								#Also replace the fortran notation '.true./.false.' with the Python 'True/False'
								if v.upper() == '.TRUE.':
									v = True
								elif v.upper() == '.FALSE.':
										v = False
								pass
							pass
						
						#Check if the field/parameter name is present in the namelist template
						if l1[0] in self._d[nl]:
							self._d[nl][l1[0]] = v
						else:
							err += "Ignored unrecognized parameter '{}'\n".format( l1[0])
				else:
					nl = None

		if err:
			raise Exception( err)
		return

pw_nl = {
	"nl" : ["CONTROL", "SYSTEM", "ELECTRONS", "IONS", "CELL"],
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