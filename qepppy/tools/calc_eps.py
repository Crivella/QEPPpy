import numpy as np
from .._decorators import numpy_save_opt, numpy_plot_opt


HARTREE = 27.2113


def calc_eps_dft(w, mel, en, fact, FAQ, weight):
	"""
	Calculate the macroscopic dielectric function for a precise value of energy
	(w) in the spectra.
	Params:
	 -w     = Spectral energy (as a complex number) in atomic units (HARTREE)
			  w = _w + i \eta where \eta is the lorentzian broadening.
	 -mel   = Numpy array of shape (pol, nt) 
			  contains the square modulus of the matrix element for every
			  polarization and transition
			  columns 4,5,6(xx,yy,zz) of 'matrixelements' (output of pw2gw)
	 -en    = Numpy array of shape (nt,)
			  contains the transition energy for every transition
			  column 7 of 'matrixelements' (output of pw2gw) (in eV)
			  !!!Must be passed in atomic units (HARTREE)
	 -fact  = Numpy array of shape (nt,)
			  contains the fact = (f_v - f_c) for every transition
			  column 8 of 'matrixelements' (output of pw2gw)     		  
	 -FAQ   = Multiplicative factor to calculate dielectric function
			  8 * pi / vol
	 -weight= Numpy array of shape (nt,)
			  contains the weight of the k_point in the BZ for every transition
			  Elemnts of 4th columns of 'k.dat' taken using the 1st colum of the
			  'matrixelements' as indexes
	Return:
	A numpy array of shape (nqpol,) containing the complex dielectric
	function for every polarization.
	Implement the formula:
	  \varepsilon(\omega) = 1 + FAQ * 
		sum_{t}[fact_t* rhotw_t*conjg(rhotw_t) * (1/(en_t - w) - 1/(-en_t - w))]
	"""
	return 1 + FAQ * np.sum(mel / en**2 * fact *(1./ (en - w) - 1./ (-en - w)) * weight, axis=1)

@numpy_plot_opt(
	_xlab=r"$\hbar\omega (eV)$", 
	_ylab=r"$\varepsilon(\omega) (arb. units)$",
	_labels=['real_x','imag_x','real_y','imag_y','real_z','imag_z',]
	)
@numpy_save_opt(
	_fname="eps_pw2gw.dat",
	_fmt=['%13.7f'] + ['%15.5e']*6,
	_header='{:>11s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}'.format(
		'omega (eV)', 'realX', 'imagX', 'realY', 'imagY', 'realZ', 'imagZ'
		)
	)
def calc_eps_pw2gw(
	vol,
	Emin=0, Emax=30, deltaE=0.05, 
	band_low=1, band_high=np.inf,
	deg=5E-2,
	mel='matrixelements', weight='k.dat', 
	):

	if isinstance(mel, str):
		mel = np.loadtxt(mel).T
	mel = mel[:,np.where((band_low <= mel[1]) & (mel[2] <= band_high))[0]]
	mel = mel[:,np.where((Emin-1 <= mel[6]) & (mel[6] <= Emax+1))[0]]

	en   = mel[6]/HARTREE
	fact = mel[7]
	mel  = mel[(3,4,5),:]

	if isinstance(weight, str):
		if weight == 'k.dat':
			kpt_loc = np.array(mel[0], dtype=int)-1
			kpt = np.loadtxt(weight).T
			weight = kpt[3,kpt_loc]
			del(kpt)
			del(kpt_loc)
		else:
			raise NotImplementedError()

	FAQ = 8 * np.pi / vol
	l1 = Emin -15*deg
	if l1<0:
		l1 = 0
	l2 = Emax + 15*deg
	omega = np.arange(l1, l2, deltaE) + 1j * deg
	omega /= HARTREE
	eps = np.empty((3, omega.shape[0]),dtype=complex)

	for n1, w in enumerate(omega):
		eps[:,n1] = calc_eps_dft(w, mel, en, fact, FAQ, weight)

	omega *= HARTREE
	w     = np.where((Emin<=omega) & (omega<=Emax))[0]
	omega = omega[w]
	eps   = eps[:,w].T

	res = np.empty((omega.shape[0],7))
	res[:,0]    = np.real(omega)
	res[:,1::2] = np.real(eps)
	res[:,2::2] = np.imag(eps)

	return res

