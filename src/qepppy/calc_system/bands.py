import numpy as np
from .kpoints import kpoints
from .._decorators import numpy_plot_opt, numpy_save_opt
# from ..meta import PropertyCreator

class bands(kpoints):
	egv = {
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'doc':"""Array of shape (n_kpt, n_bnd,) where the Nth row represents the
		energy ordered eigenvalues for the Nth kpoint."""
		}

	occ = {
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'doc':"""Array of shape (n_kpt, n_bnd,) where the Nth row and Jth
		column represent the occupation for the Jth band of the Nth kpoint."""
		}


	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	@property
	def n_bnd(self):
		return self.egv.shape[1]

	@numpy_plot_opt(_ylab="Energy (eV)")
	@numpy_save_opt(_fname="plotted.dat",_fmt="")
	def band_structure(self, thr=np.inf):
		"""
		Compute the band structure.
		Params:
		  - thr: Threshold of delta_K to be considered a cut in the band plot
		Return:
		  numpy array with shape(n_kpt,n_bnd+1).
		  The first column is the coordinates of |dK| to be used as X axis for 
		  a band plot.
		  The other column are the ordered band eigenvalue for the 
		  corresponding kpt.
		"""
		kpt = np.array(self.kpt_cart)
		kpt[1:] -= kpt[:-1]
		kpt[0]  -= kpt[0]
		norm = np.linalg.norm(kpt, axis=1)

		norm[norm > thr] = 0
		
		x = [norm[:i+1].sum() for i in range(len(norm))]
		egv = self.egv

		res = np.column_stack((x, egv))

		return res

	@numpy_plot_opt(_xlab="Energy (eV)",_ylab="DOS (arb. units)")
	@numpy_save_opt(_fname="dos.dat")
	def density_of_states(self, 
		Emin=-20, Emax=20, deltaE=0.001, deg=0.00
		):
		"""
		Compute the DOS.
		  DOS(E) = sum_{n,K} [delta(E - E_{n}(K)) * weight(K)]
		Params:
		  - Emin:   Starting energy for the DOS
		  - Emax:   Final energy for the DOS
		  - deltaE: Tick separation on the X axis
		  - deg:    Sigma to be used for a gaussian broadening.
		            Default = 0.0: Does not apply any broadening.

		Return:
		  numpy array of shape ((Emax-Emin)/(deltaE)+1,2)
		  The first column is the value of the energies generated using 
		   np.linspace(Emin,Emax,(Emax-Emin)/(deltaE)+1)
		  The second column is the value of the DOS
		"""
		res = np.linspace(Emin, Emax, (Emax-Emin)/deltaE+1).reshape(1,-1)
		res = np.pad(res, ((0,1),(0,0)), 'constant')

		for n,egv in enumerate(self.egv):
			i = np.floor((egv - Emin) / deltaE +0.5).astype(dtype='int')
			i = i[np.where( (0 <= i) & (i < res[0].size))]
			res[1,i] += self.weight[n]

		res[1:] /= deltaE

		if deg > 0:
			from ..tools.broad import broad
			res = broad(res, t='gauss', deg=deg, axis=1)

		return res.T


	