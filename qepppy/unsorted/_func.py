import re

import numpy as np

HARTREE = 27.2113

def fortran_bin_read( f, dtype, num):
    """
    Read the result of a fortran WRITE statement.
     -f: file object opened in 'rb' mode to read from
     -dtype: numpy dtype object specifing the type of data to be read
     -num: number of objects to be read
    The function will chec if dtype.size * num != number_of_writte_bytes.
    It will still works, but print a WARING message and skip the reading of the
    final 4byte in the tail of the fortran WRITE statement containing a
    duplicate of the WRITE size.
    """ 
    app1 = np.fromfile( f, np.dtype('<i4'),1)[0]
    to_read = dtype.itemsize*num
    if to_read != app1:
        print( "WARNING: reading only {} bytes out of {} from WRITE statement".format(to_read, app1))
    res = np.fromfile( f, dtype, num)
    if to_read == app1:
        app2 = np.fromfile( f, np.dtype('<i4'),1)[0]
        if app2 != app1:
            print( "WARNING: mismatch between WRITE header and tail")

    return res

def read_rhotw( fname='out.rhotw'):
    """
    Read from the dpforexc output 'out.rhotw' containing the rho_tilde elements.
     -fname: name of the file (default=out.rhotw)

    Return a tuple containing the following variables (in order):
     nt    = number of transitions
     nqpol = number of polarizations
     rhotw = numpy array of shape (nqpol, nt) with complex elements
             contains the elements of rho_tilde for every polarization and
             transition
     vcol  = numpy array of shape (nqpol,) with real elements
             contains the values of the coulomb potential for every polarization
     FAQ   = Multiplicative factor used to calculate the Polarizability (chi)
    """
    with open( 'out.rhotw', 'rb') as f:
        nt, nqpol = fortran_bin_read( f, np.dtype('<i4'), 2)
        print( "nt = ", nt, "   nqpol =", nqpol)

        rhotw=np.empty( (nqpol, nt), dtype=complex)
        for i in range(6):
            rhotw[i,:] = fortran_bin_read( f, np.dtype('<c8'), nt)
        rhotw = np.array( rhotw)

        vcol = fortran_bin_read( f, np.dtype('<f4'), nqpol)
        print( "Vc(ipol) = ", vcol)

        FAQ, = fortran_bin_read( f, np.dtype('<f4'), 1)
        print( "FAQ = ", FAQ)
    print( '\n\n')

    return nt, nqpol, rhotw, vcol, FAQ

def read_trans(fname='exc.out'):
    """
    Read from the dpforexc main output.
    Use regex syntax to extract the information about the transitions.
    Return a tuple containing the following variables (in order):
     en   = numpy array of shape (nt,) with real elements
            contains the transition energy in eV for every transition
     fact = numpy array of shape (nt,) with real elements
            contains the fact (f_i - f_j) for every transition
    """
    r = re.compile( r'\sis,ikp,iv,ik,ic,it,fact,gwten =\s+(?P<is>[\d]+)\s+(?P<kpt>[\d]+)\s+(?P<iv>[\d]+)\s+(?P<ik>[\d]+)\s+(?P<ic>[\d]+)\s+(?P<num>[\d]+)\s+(?P<fact>[\d\.]+)\s+(?P<en>[\d\.]+)')
    with open(fname, 'r') as f:
        trans = [a.groupdict() for a in  r.finditer( f.read())]
    en   = np.array([a['en']   for a in trans], dtype=float)/HARTREE
    fact = np.array([a['fact'] for a in trans], dtype=float)
    
    return en, fact


def _calc_eps_w( w, rhotw, en, fact, vcol, FAQ):
    """
    Calculate the macroscopic dielectric function for a precise value of energy
    (w) in the spectra.
     -w     = spectral energy (as a complex number) in atomic units (HARTREE)
              w = _w + i \eta where \eta is the lorentzian broadening.
     -rhotw = numpy array of shape (nqpol, nt) (from output of read_rhotw)
     -en    = numpy array of shape (nt,) (from output of read_trans)
     -fact  = numpy array of shape (nt,) (from output of read_trans)
     -FAQ   = real number (from output of read_trans)
    Return a numpy array of shape (nqpol,) containing the complex dielectric
    function for every polarization.
    Implement the formula:
      \varepsilon(\omega) = 1 + FAQ * vc * 
        sum_{t}[fact_t* rhotw_t*conjg(rhotw_t) * (1/(en_t - w) - 1/(-en_t - w))]
    """
    res = 1. + FAQ * vcol * np.sum( fact * np.absolute(rhotw)**2 * 
            (
                  1. / (en - w) 
                - 1. / (-en - w)
            ),
            axis=1)
    return res

def _calc_eps_load( broad=(5E-3,)):
    """
    Load the complex dielectric function from file.
    """
    out_res=[]
    for b in broad:
        out_res.append( np.loadtxt('eps_test_b{}.dat'.format(b)))

    return out_res

def _calc_eps_full( Emin=0, Emax=10, deltaE=1E-2, broad=(5E-3,), rhotw_name='out.rhotw', exc_name='exc.out', save=False):
    """
    Calculate the complex dielectric function starting from dpforexc output
     -Emin   = Lowest energy in spectra
     -Emax   = Highest energy in spectra
     -deltaE = Tick size in spectra
     -broad  = Lorentzian broadening used to compute spectra
               !!!!!!! Must use broad > deltaE/4 otherwise the calculation will
               still work, but the contribution of some unlucky transition will
               'magically' disappear giving a wrong result!!!!!!
     -rhotw_name = Name of the dpforexc output containing the rho_tilde elements
     -exc_name = Name of the dpforexc main output
     -save = True/False (default=False)
             Save the output in a file named eps_test_b{boad}.dat
             Format comparable to outrpanlf from dpforexc
    """
    _, nqpol, rhotw, vcol, FAQ = read_rhotw(rhotw_name)
    en, fact = read_trans(exc_name)

    out_res=[]
    for b in broad:
        omega = np.arange( Emin, Emax, deltaE) + 1j * b
        omega /= HARTREE
        eps = np.zeros((nqpol, omega.size),dtype=complex)

        for n1, w in enumerate( omega):
            eps[:,n1] = _calc_eps_w( w, rhotw, en, fact, vcol, FAQ)

        omega *= HARTREE
        res = np.zeros( (omega.shape[0], 17))
        avg1= (eps[0] + eps[1] + eps[2]) / 3
        avg2= (eps[3] + eps[4] + eps[5]) / 3
        res[:,0] = np.real(omega)
        res[:,1] = np.real(avg1)
        res[:,2] = np.imag(avg1)
        res[:,3] = np.real(avg2)
        res[:,4] = np.imag(avg2)
        res[:,5::2] = np.real(np.transpose( eps))
        res[:,6::2] = np.imag(np.transpose( eps))

        out_res.append( res)
        if save:
            np.savetxt( 'eps_test_b{}.dat'.format(b), res, fmt='%9.4f' + '%13.4e'*16)

    return out_res

def calc_eps( load=False, save=False, broad=(5E-3,), **kwargs):
    """
    Wrapper to calculate/save or load the complex dielectric function
    """
    try:
        broad = tuple(broad)
    except:
        broad = tuple( (broad,))
    if load:
        if save:
            print( "WARNING, can't load and save at same time... setting save to False")
        out_res = _calc_eps_load( broad=broad, **kwargs)
    else:
        out_res = _calc_eps_full( save=save, broad=broad, **kwargs)

    return np.array( out_res)


