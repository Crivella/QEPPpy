import numpy as np
from .kpoints import kpoints
from .._decorators import numpy_plot_opt, numpy_save_opt

def cumul_norm(
	path,
	base=np.diag([1,1,1]),
	thr=None
	):
	app       = path.dot(base)
	diff_path = app - np.vstack((app[0], app[:-1]))
	norm      = np.linalg.norm(diff_path, axis=1)

	if thr is None:
		thr = 3 * np.average(norm)
	norm[norm > thr] = 0

	res   = norm.copy()
	for i in range(1, res.shape[0]):
		res[i] += res[i-1]
	res   = np.array(res)

	return res

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

	@property
	def n_bnd(self):
		return self.egv.shape[1]

	def band_unfolding_noproject(self, SC_rec):
		PC_path  = self.kpt_cryst
		unf      = self.generate_unfolding_path(SC_rec, return_all=True)
		ri, SC_path, SC_path_red = unf

		res = np.empty((PC_path.shape[0], self.n_bnd))
		for unf_i,fol_i in enumerate(ri):
			res[unf_i] = self.egv[fol_i]

		x_coord = cumul_norm(self.kpt_cart)
		return x_coord, res

	def band_unfolding(
		self, SC_rec, wavefunctions,
		verbose = True,
		**kwargs):
		from itertools     import product
		from scipy.spatial import KDTree

		if verbose:
			import sys
			np.set_printoptions(precision=2, formatter={'float':lambda x: f' {x:+14.10f}'})

		Emin   = kwargs.pop('Emin', self.egv.min())
		Emax   = kwargs.pop('Emax', self.egv.max())
		deltaE = kwargs.pop('deltaE', 0.05)
		N      = int((Emax - Emin)/deltaE + 1)


		PC_rec   = self.recipr
		PC_path  = self.kpt_cryst
		unf      = self.generate_unfolding_path(SC_rec, return_all=True)
		ri, SC_path, SC_path_red = unf

		max_g       = np.max(np.abs(wavefunctions[0].gvect)) // 2 + 3
		l           = range(-max_g, max_g+1)
		g_list      = np.array(list(product(l,l,l)))
		g_list_cart = np.dot(g_list, PC_rec)

		A = np.zeros((PC_path.shape[0], N))
		done = []
		for unf_i,fol_i in enumerate(ri):
			check_symm = fol_i in done
			done.append(fol_i)

			psi_Km = wavefunctions[fol_i]
			c_Kmg  = psi_Km.C_kn
			G_tree = KDTree(psi_Km.gvect.dot(SC_rec))

			k0     = PC_path[unf_i].dot(PC_rec)
			K      = SC_path[unf_i].dot(SC_rec)
			dk     = k0 - K

			star = -1. if np.allclose(-SC_path[unf_i], SC_path_red[fol_i], rtol=1e-4) else 1.
			if verbose:
				print(f'KPT ({unf_i+1:3d})' + '-'*40)
				print(f'{"Unfolding vector K = ":>30s} {SC_path[unf_i]}  (K  crystal-SC) = {K} (2pi)\n'
					  f'{"on k0 = ":>30s} {PC_path[unf_i]}  (k0 crystal-PC) = {k0} (2pi)')
				print(f'{"dk = ":>30s} {dk}')
				if check_symm:
					print(f'{"":15s}' + '*'*15) 
					print(f'{"":15s}{ SC_path[unf_i]} is already calculated using symmetry.') 
					print(f'{"":15s}{ SC_path_red[fol_i]}\'s coefficients will be used (KPT #{fol_i+1:3d} in SC calclulation).')

				sys.stdout.flush()

			l = G_tree.query_ball_point(star*(g_list_cart+dk), 1E-4)
			list_rec = []
			for a in l:
				list_rec += a

			for ne,e in enumerate(self.egv[fol_i]):
				i = int(np.round((e - Emin) / deltaE))
				if 0 <= i < N:
					W = (np.abs(c_Kmg[ne, :, list_rec])**2).sum()
					A[unf_i,i] += W

		E       = np.linspace(Emin, Emax, N)
		x_coord = cumul_norm(self.kpt_cart)
		X,Y     = np.meshgrid(x_coord, E)
		Z       = A.T

		return X,Y,Z

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
		# kpt = np.array(self.kpt_cart)
		# kpt[1:] -= kpt[:-1]
		# kpt[0]  -= kpt[0]
		# norm = np.linalg.norm(kpt, axis=1)

		# norm[norm > thr] = 0
		
		# x = [norm[:i+1].sum() for i in range(len(norm))]

		x   = cumul_norm(self.kpt_cart)
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


	