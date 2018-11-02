import numpy as np
from .data_file_parser import data_file_parser as dfp

data={
	'n_kpt':{'x':'text', 'f':'output//nk', 'n':None, 't':int},
	'n_bnd':{'x':'attr', 'f':'output//ks_energies/eigenvalues', 'n':'size', 't':int},
	'n_el':{'x':'text', 'f':'output//nelec', 'n':None, 't':float},
	'fermi':{'x':'text', 'f':'output//fermi_energy', 'n':None, 't':float},
	'fermi_s':{'x':'nodelist', 'f':'output//two_fermi_energies', 'n':'fermi', 't':list},
	'homo':{'x':'text', 'f':'output//highestOccupiedLevel', 'n':None, 't':float},
	'lsda':{'x':'text', 'f':'output//lsda', 'n':None, 't':bool},
	'noncolin':{'x':'text', 'f':'output//noncolin', 'n':None, 't':bool},
	'kpt':{'x':'nodelist', 'f':'output//ks_energies/k_point', 'n':'kpt', 't':list},
	'egv':{'x':'nodelist', 'f':'output//ks_energies/eigenvalues', 'n':'egv', 't':list},
	'occ':{'x':'nodelist', 'f':'output//ks_energies/occupations', 'n':'occ', 't':list},
	}

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

	def band_structure( self, fname="plotted.dat", plot=True, pfile=True):
		"""
		Plot/print_to_file the band structure of.
		Use plot=True to plot the band structure using matplotlib
		Use pfile=True to print the band structure data to a file "fname"
		"""
		n_kpt = self.n_kpt
		kpt = [a['kpt'] for a in self.kpt]
		egv = [a['egv'] for a in self.egv]
		n_bnd = self.n_bnd
		x = [0]*n_kpt
		for i in range( 1, n_kpt):
			x[i] = x[i-1] + np.linalg.norm( kpt[i] - kpt[i-1])
		x = np.array( x)
		res = np.column_stack( ( x, egv))
		#print( res[:,0], res[:,1:])

		if pfile:
			np.savetxt( fname=fname, X=res, fmt="%13.8f"+"%11.6f"*n_bnd)

		if plot:
			import matplotlib.pyplot as plt
			from matplotlib.ticker import AutoMinorLocator as AML
			fig, ax = plt.subplots()
			plt.plot( res[:,0], res[:,1:])
			plt.ylabel( "Energy( eV)")
			plt.xlabel( "")
			ml1 = AML(5)
			ax.yaxis.set_minor_locator(ml1)
			ax.yaxis.set_tick_params(which='both', right = True)
			plt.legend()
			plt.show()
		return 0

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
		print( "SMALLEST_GAP: radius={}, comp_point={}".format( radius, comp_point))

		if( radius < 0):
			print( "Invalid negative 'radius'")
			return 1
		if( len(comp_point) != 3):
			print( "'comp_point' should be an [x,y,z] vector {}".format( comp_point))
			return 1
		for x in comp_point:
			try:
				x+=1
			except TypeError:
				print( "Invalid type %s for 'comp_point'" % type(x))
				return 1

		kpt   = np.array( [a['kpt'] for a in self.kpt])
		egv = [ a['egv'] for a in self.egv]
		n_kpt = self.n_kpt
		n_el  = self.n_el
		ef    = self.fermi
		if ef == None: ef = float('Nan')

		print( "E_fermi(from file):\t{:f} eV".format( ef))
		if ( self.lsda):
			n_el /= 2
		if ( self.noncolin):
			vb = (n_el - 1)
			print( "spin-orbit correction detected");
		else:
			vb = (n_el/2 - 1)
			print( "No spin-orbit correction")
		vb = int(vb)
		cb = vb + 1
		print( "vb = %d, cb = %d" % (vb+1, cb+1));

		base = np.array( egv)
		cp = np.array( comp_point)
		if( radius > 0):
			mod = map( lambda x: radius - np.linalg.norm(cp-x), self.kpt)
		else:
			mod = map( lambda x: 1, self.kpt)
		lcb  = base[:,cb]
		lvb  = base[:,vb]
		lcb1 = base[:,cb+1]
		#lg1  = lcb - lvb

		lll  = list( filter( lambda x: x[1]>= 0, zip( range( n_kpt), mod, kpt, lvb, lcb, lcb1)))
		found = len( lll)
		if( found <= 0):
			raise Exception( "No k-point found for the given criteria")
		print( "\nFound {} points with the given criteria.".format( found))

		res = mg1 = max( lll, key = lambda a: a[3])
		print( "\nMax_vb_energy: vb= {:f} eV".format( res[3]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))

		res = mg2 = min( lll, key = lambda a: a[4])
		print( "\nMin_cb_energy: cb= {:f} eV".format( res[4]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))

		if( mg1[0] == mg2[0]):
			print( "DIRECT GAP {:.4f} eV".format( mg2[4] - mg1[3]))
		else:
			if( mg2[4] < mg1[3]):
				print( "METALLIC")
			else:
				print( "INDIRECT GAP {:.5f} eV".format(mg2[4] - mg1[3]))

		if( ef == ef):
			res = min( (i for i in lll if i[3] < ef < i[4]), key = lambda a: a[4] - a[3])
			#mog = l[ lg1.index( m)]
			print( "\nMin_opt_gap: {:f} eV".format( res[4] - res[3]))
			print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1})".format(res[2], res[0]+1))
			print("\t%lf -> %lf   Ef: %lf eV" % (res[3], res[4], ef))		
		else:
			print( "\nCannot calculate min_opt_gap with invalid fermi energy")

		res = min( lll, key = lambda a: a[4] - a[3])
		print( "\nMin_gap_energy (vb->cb): {:f} eV".format( res[4] - res[3]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))
		print("\t%lf -> %lf   Ef: %lf eV" % (res[3], res[4], ef))

		res = min( lll, key = lambda a: a[5] - a[3])
		print( "\nMin_gap_energy (vb->cb+1): {:f} eV".format( res[5] - res[3]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))
		print("\t%lf -> %lf   Ef: %lf eV" % (res[3], res[5], ef))	

		res = min( lll, key = lambda a: a[5] - a[4])
		print( "\nMin_gap_energy (cb->cb+1): {:f} eV".format( res[5] - res[4]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))
		print("\t%lf -> %lf   Ef: %lf eV" % (res[4], res[5], ef))	

		return

	def validate( self):
		if( self.n_kpt <= 0):
			raise Exception( "No kpt read from file '{}'.".format( self.fname))
		if( self.n_bnd <= 0):
			raise Exception( "No band read from file '{}'.".format( self.fname))
		if( not self.n_kpt == len( self.egv) == len( self.occ)):
			raise Exception( "Corrupted file. Number of kpoints does not match number egv or occ")
		return True and super().validate()









