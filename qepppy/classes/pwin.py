import sys
import numpy as np

from qepppy.classes.structure import structure as structure#, bravais_index as bi

class pwin( ):
	def __init__( self, stc=None, **kwargs):
		#if not isinstance( stc, structure):
		#	raise Exception( "Must pass a valid strucutre")
		self._d = {
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
		self.stc = None
		if stc:
			if isinstance( stc, structure):
				self._add_stc_( stc)
			else:
				raise Exception( "stc is not an instance of structure")

		for k, v in kwargs.items():
			check_kw = False
			for n in self._d['nl']:
				if k in self._d[n]:
					check_kw = True
					self._d[n][k] = v
			if not check_kw:
				raise Exception( "Invalid keyword argument '{}'.".format( k))


	def __str__( self):
		self.print_in()
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
			if not stc.a:
				raise Exception("Must pass a valid cel structure")
			self._d["SYSTEM"]["ibrav"] = 0

		self._d["SYSTEM"]["celldm(1)"] = stc.lp
		if self._d["SYSTEM"]["ibrav"] == 0:
			if not isinstance( stc.a, np.ndarray):
				raise Exception( "Basis vector must be set with ibrav = 0")

		if stc.atom_spec_n != len( stc.atom_spec):
			raise Exception( "Invalide structure data, ntyp does not match")
		self._d["SYSTEM"]["ntyp"] = stc.atom_spec_n
		self._d["SYSTEM"]["nat"] = len( stc.atoms)

		return

	def print_in( self, fname=""):
		if not self.stc:
			raise Exception("No valid cell structure for the input file")

		#Define output 
		if fname:
			f = open( fname, "w+")
		else:
			f = sys.stdout

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
			return
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
							f.write( "{}".format( v))
					f.write( " ,\n")
			f.write("/\n\n")

		#Write the structure data
		f.write( self.stc.__str__())

		if fname:
			f.close()

		return


