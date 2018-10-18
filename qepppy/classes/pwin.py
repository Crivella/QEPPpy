import sys
import numpy as np

from qepppy.classes.structure import structure as structure#, bravais_index as bi

def trim_ws( str):
	l = str.strip().split(" ")
	new = ""
	check_str_1 = False
	check_str_2 = False
	for e in l:
		if check_str_1 or check_str_2:
			new += " "
		new += e
		if "\"" in l:
			check_str_1 = not check_str_1
		if "'" in l:
			check_str_2 = not check_str_2

	return new

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

class qe_in():
	def __init__( self, **kwargs):
		if not '_d' in self.__dict__:
			raise Exception( "Must first associate the parent class with a namelist dictionary")
		for k, v in kwargs.items():
			check_kw = False
			for n in self._d['nl']:
				if k in self._d[n]:
					check_kw = True
					self._d[n][k] = v
			if not check_kw:
				raise Exception( "Invalid keyword argument '{}'.".format( k))

		return

	def nl_set( self, nl, k, v):
		if nl in self._d:
			if k in self._d[nl]:
				self._d[nl][k] = v
				return True

		raise Exception( "Namelist not found")

	def qe_print( self, f=None):
		#Checks
		longest = 0
		check_mand = False
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
						sys.stderr.write( "Mandatory input parameter {} in namelist {} not set.\n".format( el, namelist))
			if not check_used:
				nl.pop( nl.index(namelist))
		if check_mand:
			return False
		longest += 2

		#Write all the used namlists/parameters
		for namelist in nl:
			f.write( "&{}\n".format(namelist))
			for el, v in self._d[namelist].items():
				if v != "":
					f.write( "{0:>{1}} = ".format( el, longest))
					if isinstance( v, str):
						f.write( "'{0}'".format( v))
					if isinstance( v, bool):
						if v: f.write( ".TRUE.")
						else: f.write( ".FALSE.")
					else:
						if isinstance( v, (int, float)):
							f.write( "{}".format( v).replace( "e", "D").replace( "E", "D"))
					f.write( " ,\n")
			f.write("/\n\n")
		return True

	def qe_read( self, fname=""):
		if not fname:
			raise Exception( "Must pass a filename to open")

		with open(fname) as f:
			content = f.readlines()
		
		nl = None
		for l in content:
			l = l.strip().split( "!")[0]
			if '&' in l and l[1:].upper() in self._d['nl']:
				nl = l[1:].upper()
				if not nl in self._d['nl']:
					raise Exception( "Reading unrecognized namelist '{}'".format( nl))
			else: 
				if (not '/' in l or '=' in l) and nl and l:
					if not nl: raise Exception( "Corrupted input file at line '{}'".format( l))
					for e in filter( None, l.split( ",")):
						l1 = trim_ws(e).split( "=")
						v = l1[1].replace("\"", "").replace("'", "")
						try:
							v = int(v)
						except:
							try:
								v = float(v.replace( "D", "e").replace( "d", "e"))
							except:
								if v.upper() == '.TRUE.':
									v = True
								elif v.upper() == '.FALSE.':
										v = False
								pass
							pass
						
						#print( nl, e, l1, v)
						if l1[0] in self._d[nl]:
							self._d[nl][l1[0]] = v
						else:
							sys.stderr.write( "Ignored unrecognized parameter '{}'\n".format( l1[0]))
				else:
					nl = None

		return




class pwin( qe_in):
	def __init__( self, stc=None, fname="", **kwargs):
		self._d = pw_nl.copy()
		self.stc = None
		if stc:
			if isinstance( stc, structure):
				self._add_stc_( stc)
			else:
				raise Exception( "stc is not an instance of structure")
		if fname:
			self.qe_read( fname)

		super().__init__( **kwargs)


	def __str__( self):
		self.qprint()
		return ""

	def __iadd__(self, other):
		if isinstance( other, structure):
			self._add_stc_( other)

		return self

	def _add_stc_( self, stc):
		self.stc = stc
		if isinstance( stc.bravais_n, int):
			self._d["SYSTEM"]["ibrav"] = stc.bravais_n
		else:
			#if not stc.a:
			raise Exception( "Must pass a valid cell structure")
			#self._d["SYSTEM"]["ibrav"] = 0

		self._d["SYSTEM"]["celldm(1)"] = stc.lp
		if self._d["SYSTEM"]["ibrav"] == 0:
			if not isinstance( stc.a, np.ndarray):
				raise Exception( "Basis vector must be set with ibrav = 0")

		#if stc.atom_spec_n != len( stc.atom_spec):
		#	raise Exception( "Invalide structure data, ntyp does not match")
		self._d["SYSTEM"]["ntyp"] = len( stc.atom_spec)
		self._d["SYSTEM"]["nat"] = len( stc.atoms)

		return

	def qe_print( self, fname=""):
		if not self.stc:
			raise Exception("No valid cell structure for the input file")

		#Define output 
		if fname:
			f = open( fname, "w+")
		else:
			f = sys.stdout

		if not super().qe_print( f):
			return

		#Write the structure data
		f.write( self.stc.__str__())

		if fname:
			f.close()

		return

	def qe_read( self, fname=""):
		super().qe_read( fname)
		self.stc = structure( fname)

		return


