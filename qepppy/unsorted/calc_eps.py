import numpy as np


HARTREE = 27.2113


def calc_eps_dft( w, mel, en, fact, FAQ, weight):
	"""
	Calculate the macroscopic dielectric function for a precise value of energy
	(w) in the spectra.
	 -w     = spectral energy (as a complex number) in atomic units (HARTREE)
			  w = _w + i \eta where \eta is the lorentzian broadening.
	 -mel   = numpy array of shape (pol, nt) 
			  contains the square modulus of the matrix element for every
			  polarization and transition
			  columns 4,5,6(xx,yy,zz) of 'matrixelements' (output of pw2gw)
	 -en    = numpy array of shape (nt,)
			  contains the transition energy for every transition
			  column 7 of 'matrixelements' (output of pw2gw) (in eV)
			  !!!Must be passed in atomic units (HARTREE)
	 -fact  = numpy array of shape (nt,)
			  contains the fact = (f_v - f_c) for every transition
			  column 8 of 'matrixelements' (output of pw2gw)     		  
	 -FAQ   = Multiplicative factor to calculate dielectric function
			  8 * pi / vol
	 -weight= numpy array of shape (nt,)
			  contains the weight of the k_point in the BZ for every transition
			  Elemnts of 4th columns of 'k.dat' taken using the 1st colum of the
			  'matrixelements' as indexes
	Return a numpy array of shape (nqpol,) containing the complex dielectric
	function for every polarization.
	Implement the formula:
	  \varepsilon(\omega) = 1 + FAQ * 
		sum_{t}[fact_t* rhotw_t*conjg(rhotw_t) * (1/(en_t - w) - 1/(-en_t - w))]
	"""
	return 1 + FAQ * np.sum( mel / en**2 * fact *( 1./ (en - w) - 1./ (-en - w)) * weight, axis=1)

def calc_eps_pw2gw( 
	vol,
	Emin=0, Emax=30, deltaE=0.05, 
	band_low=1, band_high=np.inf,
	broad=(5E-2, ),
	mel_name='matrixelements', k_name='k.dat', 
	save_name='eps_dft.dat'):

	mel = np.loadtxt( mel_name)
	mel = np.transpose( mel)
	mel = mel[:,np.where((band_low <= mel[1]) & (mel[2] <= band_high))[0]]
	mel = mel[:,np.where((Emin-1 <= mel[6]) & (mel[6] <= Emax+1))[0]]

	kpt_loc = np.array(mel[0], dtype=int)-1
	en = mel[6]/HARTREE
	fact = mel[7]
	mel = mel[(3,4,5),:]

	kpt = np.loadtxt( k_name)
	kpt = np.transpose( kpt)
	weight = kpt[3, kpt_loc]
	del( kpt)

	#vol = 532.6224
	FAQ = 8 * np.pi / vol
	for n,broad in enumerate((5E-2,)):#enumerate((2.5E-1, 1E-1, 5E-2, 1E-2, 5E-3)):
		l1 = Emin -15*broad
		if l1<0:
			l1 = 0
		l2 = Emax + 15*broad
		omega = np.arange( l1, l2, deltaE) + 1j * broad
		omega /= HARTREE
		eps = np.empty((3, omega.shape[0]),dtype=complex)

		for n1, w in enumerate( omega):
			eps[:,n1] = calc_eps_dft( w, mel, en, fact, FAQ, weight)

		omega *= HARTREE
		omega = omega[np.where((Emin<=omega) & (omega<=Emax))]
		eps = eps[:,np.where((Emin<=omega) & (omega<=Emax))[0]]
		eps = np.transpose(eps)

		out = np.empty((omega.shape[0],7))
		out[:,0] = np.real(omega)
		out[:,1::2] = np.real(eps)
		out[:,2::2] = np.imag(eps)

		if save_name:
			header = '{:>11s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}'.format('omega (eV)', 'realX', 'imagX', 'realY', 'imagY', 'realZ', 'imagZ')
			fmt = ['%13.7f'] + ['%15.5e']*6
			np.savetxt(save_name, out, fmt=fmt, header=header)

	return out

