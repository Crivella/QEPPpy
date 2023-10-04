import re

import numpy as np

from .._decorators import join_doc, numpy_plot_opt, numpy_save_opt
from ..parsers import fortran_binary_io

HARTREE = 27.2113

class dpforexc_rhotw(fortran_binary_io):
    binary_format =[
        [
            {'type':'i4', 'shape':(1,), 'name':'nt'},
            {'type':'i4', 'shape':(1,), 'name':'nqpol'},
        ],
        ([
            {'type':'c8', 'shape':('nt',), 'name':'rhotw'},
        ], 'nqpol'),
        [
            {'type':'f4', 'shape':('nqpol',), 'name':'vcol'},
        ],
        [
            {'type':'f4', 'shape':(1,), 'name':'FAQ'},
        ],
    ]
    def __init__(self, src='out.rhotw'):
        self.read_binary(src=src)


def dpforexc_read_trans(fname='exc.out'):
    """
    Read from the dpforexc main output.
    Use regex syntax to extract the information about the transitions.
    Return a tuple containing the following variables (in order):
     en   = numpy array of shape (nt,) with real elements
            contains the transition energy in eV for every transition
     fact = numpy array of shape (nt,) with real elements
            contains the fact (f_i - f_j) for every transition
    """
    r = re.compile(
        r'\sis,ikp,iv,ik,ic,it,fact,gwten =' + 
        r'\s+(?P<is>[\d]+)\s+(?P<kpt>[\d]+)\s+(?P<iv>[\d]+)\s+(?P<ik>[\d]+)\s+' + 
        r'(?P<ic>[\d]+)\s+(?P<num>[\d]+)\s+(?P<fact>[\d\.]+)\s+(?P<en>[\d\.]+)'
        )
    with open('exc.out', 'r') as f:
        trans = [a.groupdict() for a in  r.finditer(f.read())]
    iv   = np.array([a['iv'] for a in trans], dtype='int')
    ib   = np.array([a['ib'] for a in trans], dtype='int')
    en   = np.array([a['en']   for a in trans], dtype='float')/HARTREE
    fact = np.array([a['fact'] for a in trans], dtype='float')
    
    return iv, ib, en, fact


def _calc_eps_dft_(w, mel, en, fact, FAQ, weight):
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
        sum_{t}[fact_t* mel/(en^2) * (1/(en_t - w) - 1/(-en_t - w))]
    """
    return 1 + FAQ * np.sum(mel / en**2 * fact *(1./ (en - w) - 1./ (-en - w)) * weight, axis=1)

def _get_dpforexc_data_(
    Emin=0, Emax=30, deltaE=0.05, 
    band_low=1, band_high=np.inf,
    deg=5E-2,
    rhotw='out.rhotw', exc='exc.out',
    **kwargs
    ):
    """
    dpforexc mode:
     - band_low:   Exclude all the transition originating from bands lower than 
                   bands_low (default=1).
     - bands_high: Exclude all the transitions arriving at bands higher than 
                   bands_high (default=np.inf).
     - rhotw:      Path/name of the out.rhotw file from a dpforexc calculation
                   (default='out.rhotw')
     - exc:        Path/name of the exc.out file from a dpforexc calculation
                   (default='exc.out')"""
    data   = dpforexc_rhotw(rhotw)
    iv, ic, en, fact = dpforexc_read_trans(exc)

    mel    = data.rhotw * np.conj(data.rhotw) * en**2
    weight = 1
    FAQ    = data.FAQ * data.vcol

    w      = np.where(
        (band_low <= iv)    & 
        (ic <= band_high)   & 
        (Emin-20*deg <= en) & 
        (en <= Emax+20*deg) & 
        (en != 0)
        )[0]
    en     = en[w]
    fact   = fact[w]
    mel    = mel[:,w]
    return mel, weight, en, fact, FAQ

def _get_pw2gw_data(
    Emin=0, Emax=30, deltaE=0.05, 
    band_low=1, band_high=np.inf,
    deg=5E-2,
    mel='matrixelements', kdat='k.dat',
    vol=1,
    ):
    """
    pw2gw mode:
     - band_low:   Exclude all the transition originating from bands lower than 
                   bands_low (default=1).
     - bands_high: Exclude all the transitions arriving at bands higher than 
                   bands_high (default=np.inf).
     - mel:        Path/name of the matrixelements file from a pw2gw calculation
                   (default='matrixelements')
     - kdat:       Path/name of the k.dat file from a pw2gw calculation
                   (default='k.dat')
     - vol:        Volume in atomic units of the cell to be used for the 
                   normalization constant."""
    data   = np.loadtxt(mel,  comments='#').T
    weight = np.loadtxt(kdat, comments='#')[:,3]

    w      = np.where(
        (band_low <= data[1])    &
        (data[2] <= band_high)   &
        (Emin-20*deg <= data[6]) & 
        (data[6] <= Emax+20*deg) & 
        (data[6] != 0.0)
        )[0]
    data   = data[:,w]
    weight = weight[data[0].astype(dtype='int')-1]
    mel    = data[(3,4,5),:]
    en     = data[6]/HARTREE
    fact   = data[7]
    FAQ    = 2 * np.pi / vol
    return mel, weight, en, fact, FAQ

@numpy_plot_opt(
    _plot=False,
    _xlab=r"$\hbar\omega$ (eV)", 
    _ylab=r"$\varepsilon(\omega)$ (arb. units)",
    )
@numpy_save_opt(
    _fname="eps.dat",
    _fmt="%12.4E",
    _header=("{:>10s}" + "{:>12s}"*2).format("En (eV)", 'real', 'imag'),
    _delimiter='',
    )
def calc_eps(
    mode="pw2gw",
    Emin=0, Emax=30, deltaE=0.05, 
    deg=5E-2,
    **kwargs
    ):
    r"""
    Compute the complex dielectric function using the formula:
      \varepsilon(\omega) = 1 + FAQ * 
        sum_{t}[fact_t* mel/(en^2) * (1/(en_t - w) - 1/(-en_t - w))]

    Params:
     - mode:   A value from ['pw2gw', 'dpforexc']. Select the kind of output
               files to be used for the calculation.
     - Emin:   Low limit of energy range.
     - Emax:   High limit of energy range.
     - deltaE: Tick size.
     - deg:    Lorentzian broadening (must be > 0) used in 'w = r(w) + 1j * deg'"""
    modes = ['pw2gw', 'dpforexc']
    if not mode in modes:
        raise NotImplementedError("Mode {} is not implemented.".format(mode))
    if mode == 'pw2gw':
        mel, weight, en,fact, FAQ = _get_pw2gw_data(Emin=Emin, Emax=Emax, deltaE=deltaE, deg=deg, **kwargs)
    elif mode == 'dpforexc':
        mel, weight, en,fact, FAQ = _get_dpforexc_data_(Emin=Emin, Emax=Emax, deltaE=deltaE, deg=deg, **kwargs)

    npol = mel.shape[0]

    omega = np.arange(Emin, Emax+deltaE, deltaE, dtype='complex') + 1j * deg
    omega /= HARTREE

    eps = np.empty((npol, omega.size), dtype='complex')

    for n1, w in enumerate(omega):
        eps[:,n1] = _calc_eps_dft_(w, mel, en, fact, FAQ, weight)

    omega *= HARTREE
    w     = np.where((Emin<=omega) & (omega<=Emax))[0]
    omega = omega[w]
    eps   = eps[:,w].T

    res = np.empty((omega.size, 2*npol+1))

    res[:,0]    = np.real(omega)
    res[:,1::2] = np.real(eps)
    res[:,2::2] = np.imag(eps)

    res = (res[:-1] + res[1:])/2

    return res


@numpy_plot_opt(
    _plot=False,
    _xlab=r"$\hbar\omega$ (eV)", 
    _ylab=r"$\varepsilon_2(\omega)$ (arb. units)",
    _labels=['eps2 X','eps2 Y','eps2 Z','eps2 AVG',]
    )
@numpy_save_opt(
    _fname="eps.dat",
    _fmt="%9.4f" + "%12.4E"*4,
    _header=("{:>7s}" + "{:>12s}"*4).format("En (eV)", 'eps2 X','eps2 Y','eps2 Z','eps2 AVG'),
    _delimiter='',
    )
def calc_eps_pw2gw_light(
    mel="matrixelements", kdat="k.dat",
    Emin=0, Emax=30, deltaE=0.05, 
    band_low=1, band_high=np.inf,
    deg=0, deg_type='gauss',
    **kwargs
    ):
    vol = kwargs.get('vol', 1)

    x = np.arange(Emin, Emax+deltaE, deltaE)
    res = np.zeros((x.size, 4))
    res[:,0] = x

    weight = np.loadtxt(kdat, comments="#")
    weight = weight[:,3]

    with open(mel) as f:
        for line in f:
            k, v, c, px, py, pz, en, fact = np.fromstring(line, sep=' ')
            if en < Emin or Emax < en or en == 0:
                continue
            if v < band_low or c > band_high:
                continue
            p = np.array([px,py,pz])
            index = int((en-Emin)/deltaE + 0.5)
            res[index,1:] +=  p* fact  *weight[int(k)-1] / en**2

    res[:,1:] *= 2 * np.pi**2 * HARTREE**3 / vol / deltaE

    res = (res[:-1] + res[1:])/2

    res = np.column_stack((res,np.sum(res[:,1:], axis=1)/3))

    if deg > 0:
        from .broad import broad
        res = broad(res, deg_type, deg)

    return res


join_doc(calc_eps, _get_pw2gw_data.__doc__)
join_doc(calc_eps, _get_dpforexc_data_.__doc__)
