import numpy as np
from .parser.data_file_parser import data_file_parser as dfp
# from ..logger import logger, warning
from ..errors import ValidateError
from .._decorators import numpy_save_opt, numpy_plot_opt, store_property, IO_stdout_redirect

HA_to_eV = 27.21138602


data={
	'n_kpt':{
		'xml_ptype':'text', 
		'xml_search_string':'output//nks', 
		'extra_name':None, 
		'res_type':int,
		'outfile_regex':r'number of k points[\s]*='
		},
	'n_bnd':{
		'xml_ptype':'attr', 
		'xml_search_string':'output//ks_energies/eigenvalues', 
		'extra_name':'size', 
		'res_type':int,
		'outfile_regex':r'number of Kohn-Sham states[\s]*='
		},
	'n_el':{
		'xml_ptype':'text', 
		'xml_search_string':'output//nelec', 
		'extra_name':None, 
		'res_type':float,
		'outfile_regex':r'number of electrons[\s]*='
		},
	'fermi':{
		'xml_ptype':'text', 
		'xml_search_string':'output//fermi_energy', 
		'extra_name':None, 
		'res_type':float,
		'outfile_regex':r'the Fermi energy is'
		},
	'fermi_s':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//two_fermi_energies', 
		'extra_name':'fermi', 
		'res_type':list
		},
	'homo':{
		'xml_ptype':'text', 
		'xml_search_string':'output//highestOccupiedLevel', 
		'extra_name':None, 
		'res_type':float
		},
	'lsda':{
		'xml_ptype':'text', 
		'xml_search_string':'output//lsda', 
		'extra_name':None, 
		'res_type':bool
		},
	'noncolin':{
		'xml_ptype':'text', 
		'xml_search_string':'output//noncolin', 
		'extra_name':None, 
		'res_type':bool,
		'outfile_regex':r'spin'
		},
	'_kpt':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//ks_energies/k_point', 
		'extra_name':'kpt', 
		'res_type':list,
		'outfile_regex':r'[\s]{4,}k\([ \d]+\) = \((?P<kpt>[ \d\.\-]+)\).*wk = (?P<weight>[ \d\.]+)'
		},
	'_egv':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//ks_energies/eigenvalues', 
		'extra_name':'egv', 
		'res_type':list,
		'outfile_regex':r'bands \(ev\):(?P<egv>[\s\d\.\-]+)', 
		'modifier':1/HA_to_eV
		},
	'_occ':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'output//ks_energies/occupations', 
		'extra_name':'occ', 
		'res_type':list,
		'outfile_regex':r'occupation numbers(?P<occ>[\s\d\.]+)'
		},
	}

# @logger()
class bands(dfp):
	"""
	Instance used for QE eigenvalues/vector(k-points) and occupations numbers.
	Uses the internal "data_file_parser" to read from a "data-file-schema.xml"
	or from a pw.x output file.
	Can be printed as a string.
	Each k-point and its info can be called as a dictionary value using its number as the key.
	Provide the following PostProcessing methods:
	- band_structure(): Plot/print_to_file the band structure.
	- smallest_gap(): Print an analysis of the band gap.
	"""
	__name__ = "bands"
	e_units = HA_to_eV
	def __init__(self, d={}, **kwargs):
		d.update(data)
		super().__init__(d=d, **kwargs)
		return

	def __str__(self):
		msg = super().__str__()
		bnd = self.n_bnd
		kpt_fmt = "\nkpt(#{{:5d}}):  " + "{:8.4f}"*3 + " [2pi/alat]"
		egv_fmt = "\nEigenvalues(eV):\n" + ("  "+"{:12.6f}"*8+"\n")*(bnd//8)
		egv_fmt += "  " + "{:12.6f}"*(bnd%8) + "\n"
		for i in range(self.n_kpt):
			msg += kpt_fmt.format(*self.kpt_cart[i]).format(i)
			msg += egv_fmt.format(*self.egv[i])
		return msg

	def __getitem__(self, key):
		if(isinstance(key, int)):
			if(0 <= key < self.n_kpt):
				return {'kpt':self.kpt_cart[key], 'egv':self.egv[key], 'occ':self.occ[key]} 
			else:
				raise KeyError("Index '{}' out of range {}-{}".format(key, 0, self.n_kpt - 1))
		return super().__getitem__(key)

	@property
	@store_property
	def kpt_cart(self):
		"""
		np.ndarray of shape(nkpt, 3).
		Contains the coordinates in cartesian units of all k-points.
		"""
		kpt = np.array([a['kpt'] for a in self._kpt])
		kpt = kpt[:self.n_kpt,:]
		return kpt

	@property
	@store_property
	def kpt_cryst(self):
		"""
		np.ndarray of shape(nkpt, 3).
		Contains the coordinates in crystal units of all k-points.
		"""
		n = self.n_kpt
		kpt = np.array([a['kpt'] for a in self._kpt])
		if kpt.shape[0] > n:
			kpt = kpt[n:2*n,:]
		else:
			raise NotImplementedError()
		return kpt

	@property
	@store_property
	def weight(self):
		"""
		np.ndarray of shape (nkpt,)
		Contains the weights for all k-points.
		"""
		n = self.n_kpt
		occ = np.array([a['weight'] for a in self._kpt])
		if occ.shape[0] > n:
			occ = occ[:n]
		return occ

	@property
	@store_property
	def egv(self):
		"""
		np.ndarray of shape (nkpt,nbnd).
		Contains all the eigenvalues for the bands.
		"""
		return np.array([a['egv'] for a in self._egv]) * self.e_units

	@property
	@store_property
	def occ(self):
		"""
		np.ndarray of shape (nkpt,nbnd).
		Contains the occupation number of all the bands.
		"""
		return np.array([a['occ'] for a in self._occ])

	@numpy_save_opt(_fname='kpt.dat', _fmt="%14.6f")
	def kpt_crop(self, center, radius, mode='cart'):
		"""
		Crop all k-point from the grid, in a sphere around a center.
		Params:
		 - center: Iterable containing the X,Y,Z coordinates of the center of 
		           the sphere.
		 - radius: Radius of the sphere
		 - mode:   Crop the k-points from cartesian('cart') or crystal('cryst')
		           coordinates.
		"""
		center = np.array(center)
		if center.size != 3:
			raise ValueError("First positional argument must be a 3D coordinate for the center")

		if mode == 'cart':
			kpt = self.kpt_cart
		elif mode == 'cryst':
			kpt = self.kpt_cryst

		norm = np.linalg.norm(kpt - center, axis=1)
		w = np.where(norm <= radius)[0]

		print("Cropping {} points out of {}.".format(w.size, norm.size))
		w_red = np.sum(self.weight[w])
		w_tot = np.sum(self.weight)
		print("Total weight {:.5f} / {:.5f} = {:.5f} = Renormalization factor".format(
			w_red, w_tot, w_red/w_tot
			))

		return np.column_stack((kpt[w,:], self.weight[w]))
	
	

	@numpy_plot_opt(_ylab="Energy (eV)")
	@numpy_save_opt(_fname="plotted.dat",_fmt="")
	def band_structure(self):
		"""
		Compute the band structure.
		Params:
		  -
		Return:
		  numpy array with shape(n_kpt,n_bnd+1).
		  The first column is the coordinates of |dK| to be used as X axis for 
		  a band plot.
		  The other column are the ordered band eigenvalue for the 
		  corresponding kpt.
		"""
		kpt = self.kpt_cart
		# kpt = kpt[:self.n_kpt,:]
		egv = self.egv

		x = np.linalg.norm(kpt, axis=1)
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
		  The first column is the value of the energies generated using np.linspace(Emin,Emax,(Emax-Emin)/(deltaE)+1)
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

	@IO_stdout_redirect()
	def smallest_gap(self,
		radius=0., comp_point=(0.,0.,0.), **kwargs
		):
		"""
		Print to screen the following information concerning the band gap:
		  - Fermi energy
		  - Direct gap
		  - Min/Max of conduction/valence band (Indirect gap)
		  - Optical gap (condition: valence < Fermi < conduction)
		Can focus on only a small portion of k-points by defining:
		Params:
		  - radius: radius of the crop sphere centered around comp_point
		  - comp_point: tuple of coordinates for the center of the crop sphere
		"""
		print("SMALLEST_GAP: radius={}, comp_point={}\n".format(radius, comp_point))

		if(radius < 0):
			raise ValueError("Radius can't be negative")
		cp = np.array(comp_point, dtype='float')
		if(cp.shape != (3,)):
			raise ValueError("'comp_point' should be an [x,y,z] vector {}".format(cp))

		kpt   = self.kpt_cart
		egv   = self.egv
		n_el  = self.n_el
		ef    = self.fermi

		if ef == None: 
			ef = np.nan
		print("E_fermi(from file):\t{:f} eV".format(ef))

		if self.lsda:
			n_el /= 2
		if self.noncolin:
			vb = (n_el - 1)
			print("spin-orbit correction detected")
		else:
			vb = (n_el/2 - 1)
			print("No spin-orbit correction")
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
			raise Exception("No k-point found for the given criteria")
		print("\nFound {} points with the given criteria.".format(found))

		mg1 = np.argmax(egv[:,vb])
		top_valence = egv[mg1,vb]
		print("\nMax_vb_energy: vb= {:f} eV".format(top_valence))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[mg1], num[mg1]+1))

		mg2 = np.argmin(egv[:,cb])
		bot_conduction = egv[mg2,cb]
		print("\nMin_cb_energy: cb = {:f} eV".format(bot_conduction))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[mg2], num[mg2]+1))

		gap = bot_conduction - top_valence
		if mg1 == mg2:
			print("DIRECT GAP {:.4f} eV".format(gap))
		else:
			if(gap < 0):
				print("METALLIC")
			else:
				print("INDIRECT GAP {:.5f} eV".format(gap))

		
		if not np.isnan(ef):
			w = np.where((egv[:,vb] < ef) & (egv[:,cb] > ef))[0]
			app_egv = egv[w,:]
			app_kpt = kpt[w,:]
			app_num = num[w]
			res = np.argmin(app_egv[:,cb] - app_egv[:,vb])
			opt_gap = app_egv[res,cb] - app_egv[res,vb]
			print("\nMin_opt_gap: {:f} eV".format(opt_gap))
			print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1})".format(app_kpt[res], app_num[res]+1))
			print("\t{} -> {}   Ef: {} eV".format(app_egv[res,vb], app_egv[res,cb], ef))		
		else:
			print("\nCannot calculate min_opt_gap with invalid fermi energy")

		res = np.argmin(egv[:,cb] - egv[:,vb])
		print("\nMin_gap_energy (vb->cb): {:f} eV".format(egv[res,cb] - egv[res,vb]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[res], num[res]+1))
		print("\t{} -> {}   Ef: {} eV".format(egv[res,vb], egv[res,cb], ef))

		res = np.argmin(egv[:,cb+1] - egv[:,vb])
		print("\nMin_gap_energy (vb->cb+1): {:f} eV".format(egv[res,cb+1] - egv[res,vb]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[res], num[res]+1))
		print("\t{} -> {}   Ef: {} eV".format(egv[res,vb], egv[res,cb+1], ef))

		res = np.argmin(egv[:,cb+1] - egv[:,cb])
		print("\nMin_gap_energy (cb->cb+1): {:f} eV".format(egv[res,cb+1] - egv[res,cb]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format(kpt[res], num[res]+1))
		print("\t{} -> {}   Ef: {} eV".format(egv[res,cb], egv[res,cb+1], ef))

	def validate(self):
		if self.n_kpt <= 0:
			raise ValidateError("Failed to read nkpt from file '{}'.".format(self.xml))
		if self.n_bnd <= 0:
			raise ValidateError("Failed to read nbnd from file '{}'.".format(self.xml))
		legv = self.egv.shape[0]
		if self.occ.size:
			locc = self.occ.shape[0]
		else:
			locc = legv
		if not self.n_kpt == legv == locc:
			raise ValidateError("Corrupted file. Number of k-points does not match number egv or occ")
		super().validate()









