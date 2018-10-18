import sys

from qepppy.classes.structure import structure as structure, bravais_index as bi

class pwin( ):
	def __init__( self, stc=""):
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
				"prefix" : "",
				"lkpoint_dir" : "",
				"max_seconds" : "",
				"etot_conv_thr" : "",
				"forc_conv_thr" : "",
				"disk_io" : "",
				"pseudo_dir" : "",
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
				"ibrav" : "",
				"celldm" : "",
				"A" : "",
				"B" : "",
				"C" : "",
				"cosAB" : "",
				"cosAC" : "",
				"cosBC" : "",
				"nat" : "",
				"ntyp" : "",
				"nbnd" : "",
				"tot_charge" : "",
				"starting_charge" : "",
				"tot_magnetization" : "",
				"starting_magnetization" : "",
				"ecutwfc" : "",
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

	def __str__( self):
		self.print()
		return ""

	def print( self, fname=""):
		if fname:
			f = open( fname, "w+")
			pfunc = f.write
		else:
			pfunc = sys.stdout.write
		longest = 0
		for namelist in self._d["nl"]:
			for el, v in self._d[namelist].items():
				if v != "":
					app = len(el)
					if longest < app: longest = app
		longest += 2

		for namelist in self._d["nl"]:
			pfunc( "&{}\n".format(namelist))
			for el, v in self._d[namelist].items():
				if v != "":
					pfunc( "{0:>{1}} = ".format( el, longest))
					if isinstance( v, str):
						pfunc( "'{0}'".format( v))
					if isinstance( v, bool):
						if v: pfunc( ".TRUE.")
						else: pfunc( ".FALSE.")
					else:
						if isinstance( v, (int, float)):
							pfunc( v)
					pfunc( " ,\n")
			pfunc("/\n")

		return

"""
if __name__ == "__main__":
	a = pwin()
	a.print()
"""

