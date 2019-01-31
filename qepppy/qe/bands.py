import numpy as np
from .parser.data_file_parser import data_file_parser as dfp
from ..logger import logger, warning
from .._decorators import save_opt, plot_opt


data={
	'n_kpt':{'x':'text', 'f':'output//nk', 'n':None, 't':int, 
		'bu':r'number of k points[\s]*='},
	'n_bnd':{'x':'attr', 'f':'output//ks_energies/eigenvalues', 'n':'size', 't':int, 
		'bu':r'number of Kohn-Sham states[\s]*='},
	'n_el':{'x':'text', 'f':'output//nelec', 'n':None, 't':float, 
		'bu':r'number of electrons[\s]*='},
	'fermi':{'x':'text', 'f':'output//fermi_energy', 'n':None, 't':float, 
		'bu':r'the Fermi energy is'},
	'fermi_s':{'x':'nodelist', 'f':'output//two_fermi_energies', 'n':'fermi', 't':list},
	'homo':{'x':'text', 'f':'output//highestOccupiedLevel', 'n':None, 't':float},
	'lsda':{'x':'text', 'f':'output//lsda', 'n':None, 't':bool},
	'noncolin':{'x':'text', 'f':'output//noncolin', 'n':None, 't':bool, 
		'bu':r'spin'},
	'kpt':{'x':'nodelist', 'f':'output//ks_energies/k_point', 'n':'kpt', 't':list, 
		'bu':r'[\s]{4,}k\([ \d]+\) = \((?P<kpt>[ \d\.\-]+)\).*wk = (?P<weight>[ \d\.]+)'},
	'egv':{'x':'nodelist', 'f':'output//ks_energies/eigenvalues', 'n':'egv', 't':list, 
		'bu':r'bands \(ev\):(?P<egv>[\s\d\.\-]+)', 'm':1/27.21138602},
	'occ':{'x':'nodelist', 'f':'output//ks_energies/occupations', 'n':'occ', 't':list, 
		'bu':r'occupation numbers(?P<occ>[\s\d\.]+)'},
	}

@logger()
class bands( dfp):
	"""
	Instance used for QE eigenvalues/vector(k-points) and occupations numbers.
	Uses the internal "data_file_parser" to read from a "data-file*.xml" input.
	Can be printed as a string.
	Each k-point and its info can be called as a dictionary value using its number as the key.
	Provide the following PostProcessing methods:
	- band_structure(): Plot/print_to_file the band structure.
	- smallest_gap(): Print an analysis of the band gap.
	"""
	__name__ = "bands"
	e_units = 27.21138602
	def __init__( self, d={}, **kwargs):
		d.update( data)
		super().__init__( d=d, **kwargs)
		return

	def __str__( self):
		msg = super().__str__()
		bnd = len(self.egv[0]['egv'])-1
		kpt_fmt = "\nkpt(#{{:5d}}):  " + "{:8.4f}"*3 + " [2pi/alat]"
		egv_fmt = "\nEigenvalues( eV):\n"+("  "+"{:12.6f}"*8+"\n")*int(bnd/8)
		egv_fmt += "  "+"{:12.6f}"*(bnd%8+1)+"\n"
		for i in range( self.n_kpt):
			msg += kpt_fmt.format( *self.kpt[i]['kpt']).format( i)
			msg += egv_fmt.format( *self.egv[i]['egv'])
		return msg

	def __getitem__( self, key):
		if( isinstance( key, int)):
			if( 0 <= key < self.n_kpt):
				return { 'kpt':self.kpt[key], 'egv':self.egv[key], 'occ':self.occ[key]} 
			else:
				raise Exception( "Index '{}' out of range {}-{}".format( key, 0, self.n_kpt - 1))
		return super().__getitem__( key)

	@plot_opt
	@save_opt(_fname="plotted.dat")
	def band_structure( 
		self, *args, 
		# fname="plotted.dat", 
		fmt="", 
		ylab="Energy (eV)",
		**kwargs
		):
		"""
		Plot/print_to_file the band structure of.
		Use plot=True to plot the band structure using matplotlib
		Use pfile=True to print the band structure data to a file "fname"
		"""
		kpt = np.array([a['kpt'] for a in self.kpt])
		kpt = kpt[:self.n_kpt,:]
		egv = np.array( [a['egv'] for a in self.egv]) * self.e_units

		x = np.linalg.norm(kpt, axis=1)
		res = np.column_stack( ( x, egv))

		return res

	@plot_opt
	@save_opt(_fname="dos.dat")
	def density_of_states( 
		self, *args, 
		emin=-20, emax=20, deltaE=0.001, deg=0.00, 
		# fname="dos.dat",
		**kwargs
		):
		"""
		Holy ship
		"""
		x = np.linspace( emin, emax, (emax-emin)/deltaE+1)
		y = np.zeros( x.size)

		for n,egv in enumerate( self.egv):
			for e in egv['egv']:
				index = int((e*self.e_units - emin) / deltaE)
				if 0 <= index < x.size:
					y[index] += self.kpt[n]['weight']

		y /= deltaE

		data = np.vstack((x,y))
		if deg > 0:
			from ..tools.broad import broad
			data = broad( data, t='gauss', deg=deg, axis=1)

		return data.T

	def smallest_gap( self, radius=0., comp_point=(0.,0.,0.)):
		"""
		Can focus on only a small portion of k-points by defining:
		- radius: radius of the crop sphere centered around comp_point
		- comp_point: tuple of coordinates for the center of the crop sphere
		Print to screen the following information concerning the band gap:
		- Fermi energy
		- Direct gap
		- Min/Max of conduction/valence band (Indirect gap)
		- Optical gap (condition: valence < Fermi < conduction)
		"""
		print( "\nSMALLEST_GAP: radius={}, comp_point={}\n".format( radius, comp_point))

		if( radius < 0):
			raise ValueError( "Radius can't be negative")
		cp = np.array(comp_point, dtype='float')
		if( cp.shape != (3,)):
			raise ValueError("'comp_point' should be an [x,y,z] vector {}".format( cp))

		kpt   = np.array( [a['kpt'] for a in self.kpt])
		egv   = np.array([ a['egv'] for a in self.egv])*self.e_units
		n_el  = self.n_el
		ef    = self.fermi
		kpt   = kpt[:self.n_kpt,:]

		if ef == None: 
			ef = np.nan
		print( "E_fermi(from file):\t{:f} eV".format( ef))

		if self.lsda:
			n_el /= 2
		if self.noncolin:
			vb = (n_el - 1)
			print("spin-orbit correction detected")
		else:
			vb = (n_el/2 - 1)
			print( "No spin-orbit correction")
		vb = int(vb)
		cb = vb + 1
		print("vb = {}, cb = {}".format(vb+1, cb+1))

		mod  = np.linalg.norm(cp-kpt,axis=1) - radius
		in_range = np.where(mod >= 0)[0]

		num = np.arange(self.n_kpt)
		num = num[in_range]
		kpt = kpt[in_range,:]
		egv = egv[in_range,:]
		found = kpt.shape[0]
		if found <= 0:
			raise Exception( "No k-point found for the given criteria")
		print( "\nFound {} points with the given criteria.".format( found))

		mg1 = np.argmax(egv[:,vb])
		top_valence = egv[mg1,vb]
		print( "\nMax_vb_energy: vb= {:f} eV".format( top_valence))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( kpt[mg1], num[mg1]+1))

		mg2 = np.argmin(egv[:,cb])
		bot_conduction = egv[mg2,cb]
		print( "\nMin_cb_energy: cb = {:f} eV".format( bot_conduction))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( kpt[mg2], num[mg2]+1))

		gap = bot_conduction - top_valence
		if mg1 == mg2:
			print( "DIRECT GAP {:.4f} eV".format( gap))
		else:
			if( gap < 0):
				print( "METALLIC")
			else:
				print( "INDIRECT GAP {:.5f} eV".format(gap))

		
		if not np.isnan(ef):
			w = np.where((egv[:,vb] < ef) & (egv[:,cb] > ef))
			res = np.argmin(egv[w,cb] - egv[w,vb])
			opt_gap = egv[res,cb] - egv[res,vb]
			print( "\nMin_opt_gap: {:f} eV".format( opt_gap))
			print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1})".format(kpt[res], num[res]+1))
			print("\t{} -> {}   Ef: {} eV".format(egv[res,vb], egv[res,cb], ef))		
		else:
			print( "\nCannot calculate min_opt_gap with invalid fermi energy")

		res = np.argmin(egv[:,cb] - egv[:,vb])
		print( "\nMin_gap_energy (vb->cb): {:f} eV".format( egv[res,cb] - egv[res,vb]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( kpt[res], num[res]+1))
		print("\t{} -> {}   Ef: {} eV".format(egv[res,vb], egv[res,cb], ef))

		res = np.argmin(egv[:,cb+1] - egv[:,vb])
		print( "\nMin_gap_energy (vb->cb+1): {:f} eV".format( egv[res,cb+1] - egv[res,vb]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( kpt[res], num[res]+1))
		print("\t{} -> {}   Ef: {} eV".format(egv[res,vb], egv[res,cb+1], ef))

		res = np.argmin(egv[:,cb+1] - egv[:,cb])
		print( "\nMin_gap_energy (cb->cb+1): {:f} eV".format( egv[res,cb+1] - egv[res,cb]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( kpt[res], num[res]+1))
		print("\t{} -> {}   Ef: {} eV".format(egv[res,cb], egv[res,cb+1], ef))

	def validate( self):
		ret = True
		if self.n_kpt <= 0:
			warning.print( "Failed to read nkpt from file '{}'.".format( self.schema))
			ret = False
			#raise Exception( "No kpt read from file '{}'.".format( self.fname))
		if self.n_bnd <= 0:
			warning.print( "Failed to read nbnd from file '{}'.".format( self.schema))
			ret = False
			#raise Exception( "No band read from file '{}'.".format( self.fname))
		legv = len( self.egv)
		if self.occ:
			locc = len( self.occ)
		else:
			locc = legv
		if not self.n_kpt == legv == locc:
			warning.print( "Corrupted file. Number of kpoints does not match number egv or occ")
			ret = False
			#raise Exception( "Corrupted file. Number of kpoints does not match number egv or occ")
		return ret and super().validate()









