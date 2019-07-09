import numpy as np
from .kpoints import kpoints
from .._decorators import numpy_plot_opt, numpy_save_opt, IO_stdout_redirect, set_self

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
	n_el = {
		'typ':(int, np.integer),
		'doc':"""Number of electrons in the system."""
		}
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
	def vb(self):
		"""Valence band index"""
		return self.n_el - 1

	@property
	def cb(self):
		"""Conduction band index"""
		return self.vb + 1

	@property
	def n_bnd(self):
		try:
			res = self.egv.shape[1]
		except AttributeError:
			res = None
		if hasattr(self, '_n_bnd'):
			try:
				res = self._n_bnd
			except:
				self._n_bnd = res
		return res

	@n_bnd.setter
	def n_bnd(self, value):
		self._n_bnd = value

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

		# calc_kpt = np.array([a.kpt for a in wavefunctions])
		# print(SC_path.shape, SC_path_red.shape)
		# for n in range(len(calc_kpt)):
		# 	print(calc_kpt[n]*20.42/6.2831, '   ', SC_path_red[n], '  ->  ', calc_kpt[n]*20.42/6.2831 - SC_path_red[n])
		# exit()

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
			res[1,i] += self.kpt_weight[n]

		res[1:] /= deltaE

		if deg > 0:
			from ..tools.broad import broad
			res = broad(res, t='gauss', deg=deg, axis=1)

		return res.T

	def bands_extrema(self, bnd):
		app = self.egv[:, bnd-1]
		return {
			'min':(np.argmin(app),app.min()),
			'max':(np.argmax(app),app.max())
			}

	@set_self('egv,occ')
	def kpt_crop(
		self,
		center=(0,0,0), radius=np.inf, **kwargs
		):
		center = np.array(center).reshape(3)
		norms  = np.linalg.norm(self.kpt_cart - center, axis=1)
		w      = np.where(norms <= radius)[0]

		super().kpt_crop(center=center, radius=radius, **kwargs)

		return self.egv[w], self.occ[w]

	@IO_stdout_redirect()
	def smallest_gap(self,
		center=(0.,0.,0.), radius=np.inf, 
		verbose=True, **kwargs
		):
		"""
		Print to screen the following information concerning the band gap:
		  - Fermi energy
		  - Direct gap
		  - Min/Max of conduction/valence band (Indirect gap)
		  - Optical gap (condition: valence < Fermi < conduction)
		Can focus on only a small portion of k-points by defining:
		Params:
		  - center: tuple of coordinates for the center of the crop sphere
		  - radius: radius of the crop sphere centered around comp_point
		"""
		def _print(*args, **kwargs):
			if verbose:
				print(*args, **kwargs)

		_print("SMALLEST_GAP: radius={}, comp_point={}\n".format(radius, center))

		center = np.array(center, dtype=np.float).reshape(3)

		norms  = np.linalg.norm(self.kpt_cart - center, axis=1)
		w      = np.where(norms <= radius)[0]

		num    = np.arange(self.n_kpt)[w]
		kpt    = self.kpt_cart[w]
		egv    = self.egv[w]
		# n_el        = self.n_el

		if len(kpt) == 0:
			_print("No k-point found for the given criteria.")
			return

		ef = np.nan
		# ef    = self.fermi
		# if ef == None: 
		# 	ef = np.nan
		# _print("E_fermi(from file):\t{:f} eV".format(ef))

		vb = self.vb
		cb = self.cb
		_print("vb = {}, cb = {}".format(vb+1, cb+1))


		_print("\nFound {} points with the given criteria.".format(kpt.shape[0]))

		mg1 = np.argmax(egv[:,vb])
		top_valence = egv[mg1,vb]
		_print("\nMax_vb_energy: vb= {:f} eV".format(top_valence))
		_print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[mg1], num[mg1]+1))

		mg2 = np.argmin(egv[:,cb])
		bot_conduction = egv[mg2,cb]
		_print("\nMin_cb_energy: cb = {:f} eV".format(bot_conduction))
		_print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[mg2], num[mg2]+1))

		gap = bot_conduction - top_valence
		if mg1 == mg2:
			_print("DIRECT GAP {:.5f} eV".format(gap))
		else:
			if(gap < 0):
				_print("METALLIC")
			else:
				_print("INDIRECT GAP {:.5f} eV".format(gap))

		
		if not np.isnan(ef):
			w = np.where((egv[:,vb] < ef) & (egv[:,cb] > ef))[0]
			app_egv = egv[w,:]
			app_kpt = kpt[w,:]
			app_num = num[w]
			res = np.argmin(app_egv[:,cb] - app_egv[:,vb])
			opt_gap = app_egv[res,cb] - app_egv[res,vb]
			_print("\nMin_opt_gap: {:f} eV".format(opt_gap))
			_print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1})".format(app_kpt[res], app_num[res]+1))
			_print("\t{} -> {}   Ef: {} eV".format(app_egv[res,vb], app_egv[res,cb], ef))		
		else:
			_print("\nCannot calculate min_opt_gap with invalid fermi energy")

		res = np.argmin(egv[:,cb] - egv[:,vb])
		_print("\nMin_gap_energy (vb->cb): {:f} eV".format(egv[res,cb] - egv[res,vb]))
		_print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[res], num[res]+1))
		_print("\t{} -> {}   Ef: {} eV".format(egv[res,vb], egv[res,cb], ef))

		res = np.argmin(egv[:,cb+1] - egv[:,vb])
		_print("\nMin_gap_energy (vb->cb+1): {:f} eV".format(egv[res,cb+1] - egv[res,vb]))
		_print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[res], num[res]+1))
		_print("\t{} -> {}   Ef: {} eV".format(egv[res,vb], egv[res,cb+1], ef))

		res = np.argmin(egv[:,cb+1] - egv[:,cb])
		_print("\nMin_gap_energy (cb->cb+1): {:f} eV".format(egv[res,cb+1] - egv[res,cb]))
		_print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[res], num[res]+1))
		_print("\t{} -> {}   Ef: {} eV".format(egv[res,cb], egv[res,cb+1], ef))


	# def fit_analysis(self, n_pt=5):
	# 	"""
	# 	Fit analysis of the valnece/conduction band extrema.
	# 	Params:
	# 	 -n_pt = 5: Number of points around the band extrema to use for the fit.

	# 	Return :
	# 	 Dictionary with the following structure
	# 	 {
	# 	 	'linear':{      #Linear fit result
	# 	 		'vb':(...), #Return of scipy.optimize.curve_fit for valence bnd
	# 	 		'cb':(...)  #Return of scipy.optimize.curve_fit for conduct bnd
	# 	 	},
	# 	 	'quadratic':{   #Quadratic fit result
	# 	 		'vb':(...), #Return of scipy.optimize.curve_fit for valence bnd
	# 	 		'cb':(...)  #Return of scipy.optimize.curve_fit for conduct bnd
	# 	 	}
	# 	 }
	# 	"""
	# 	import scipy
	# 	def linear(x, a, b):
	# 		return a*x + b
	# 	def quadratic(x, a, b, c):
	# 		return a*x**2 + b*x + c

	# 	bands = self.band_structure(pFile=False, plot=False)
	# 	ext_v = self.bands_extrema(self.vb)
	# 	ext_c = self.bands_extrema(self.cb)

	# 	x = bands[:,0]

	# 	res = {}

	# 	print("LINEAR:")
	# 	ptr = res['linear'] = {}
	# 	i   = ext_v['max'][0]
	# 	sl  = slice(i-n_pt,i+1)
	# 	print(sl)
	# 	app = scipy.optimize.curve_fit(linear, x[sl], bands[sl,self.vb])
	# 	ptr['vb'] = app
	# 	print("--Valence    band (left fit):")
	# 	print(" "*4 + "Vf = {:>12.4E} eV/(2pi/alat)".format(app[0][0]))
	# 	print(" "*9 + "{:>12.4E} m/s".format(app[0][0] * self.alat * 1.28E4))

	# 	i   = ext_c['min'][0]
	# 	sl  = slice(i-n_pt,i+1)
	# 	app = scipy.optimize.curve_fit(linear, x[sl], bands[sl,self.cb])
	# 	ptr['cb'] = app
	# 	print("--Conduction band (left fit):")
	# 	print(" "*4 + "Vf = {:>12.4E} eV/(2pi/alat)".format(app[0][0]))
	# 	print(" "*9 + "{:>12.4E} m/s".format(app[0][0] * self.alat * 1.28E4))

	# 	print("QUADRATIC:")
	# 	ptr = res['quadratic'] = {}
	# 	i   = ext_v['max'][0]
	# 	sl  = slice(i-n_pt,i+n_pt)
	# 	app = scipy.optimize.curve_fit(quadratic, x[sl], bands[sl,self.vb])
	# 	ptr['vb'] = app
	# 	print("--Valence band")
	# 	print("    f(x) = {:>12.4E} * x^2 + {:>12.4E} * x + {:>12.4E}".format(*app[0]))

	# 	i   = ext_c['min'][0]
	# 	sl  = slice(i-n_pt,i+n_pt)
	# 	app = scipy.optimize.curve_fit(quadratic, x[sl], bands[sl,self.cb])
	# 	ptr['cb'] = app
	# 	print("--Conduction band")
	# 	print("    f(x) = {:>12.4E} * x^2 + {:>12.4E} * x + {:>12.4E}".format(*app[0]))

	# 	return res


	