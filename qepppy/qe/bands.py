import numpy as np
from .parser.data_file_parser import data_file_parser as dfp
from ..logger import logger, warning
from .._decorators import numpy_save_opt, numpy_plot_opt, store_property, IO_stdout_redirect

HA_to_eV = 27.21138602


data={
	'n_kpt':{
		'xml_ptype':'text', 
		'xml_search_string':'output//nk', 
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

@logger()
class bands(dfp):
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
		kpt = np.array([a['kpt'] for a in self._kpt])
		kpt = kpt[:self.n_kpt,:]
		return kpt

	@property
	@store_property
	def kpt_cryst(self):
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
		return np.array([a['weight'] for a in self._kpt])

	@property
	@store_property
	def egv(self):
		return np.array([a['egv'] for a in self._egv]) * self.e_units

	@property
	@store_property
	def occ(self):
		return np.array([a['occ'] for a in self._occ])
	
	

	@numpy_plot_opt(_ylab="Energy (eV)")
	@numpy_save_opt(_fname="plotted.dat",_fmt="")
	def band_structure(
		self, *args,
		**kwargs
		):
		"""
		Compute the band structure.
		Params:
		  -
		Return:
		  numpy array with shape(n_kpt,n_bnd+1).
		  The first column is the coordinates of |dK| to be used as X axis for a band plot.
		  The other column are the ordered band eigenvalue for the corresponding kpt.

		numpy_save_opt specific params:
		  - pFile: (True/False) Enable/disable save functionality (default = True)
		  - fname: Output file name (must be present)
		  - fmt:   Format string to pass to np.savetxt

		numpy_plot_opt specific params:
		  - plot:      (True/False) Enable/disable plot functionality (default = True)
		  - xlab:      String to use as x label
		  - ylab:      String to use as y label
		  - start:     First column of the Y axis data to be plotted
		  - end:       Last column of the Y axis data to be plotted
		  - colors:    List of matplotlib color string to be used.
		               % is used to loop if (end-start > len(colors))
		  - labels:    List of strings to be used as labes.
		               No label is set if (end-start > len(labels))
		  - dash_list: List of tuples of dashes option to be used.
		               % is used to loop if (end-start > len(colors))
		               If no dash_list is specified, the lines will switch from nodash to dash=(8,2)
		               for every loop of the colors
		"""
		kpt = self.kpt_cart
		# kpt = kpt[:self.n_kpt,:]
		egv = self.egv

		x = np.linalg.norm(kpt, axis=1)
		res = np.column_stack((x, egv))

		return res

	@numpy_plot_opt(_xlab="Energy (eV)",_ylab="DOS (arb. units)")
	@numpy_save_opt(_fname="dos.dat")
	def density_of_states(
		self, *args, 
		emin=-20, emax=20, deltaE=0.001, deg=0.00, 
		**kwargs
		):
		"""
		Compute the DOS.
		  DOS(E) = sum_{n,K} [delta(E - E_{n}(K)) * weight(K)]
		Params:
		  - emin:   Starting energy for the DOS
		  - emax:   Final energy for the DOS
		  - deltaE: Tick separation on the X axis
		  - deg:    Sigma to be used for a gaussian broadening.
		            Default = 0.0: Does not apply any broadening.

		Return:
		  numpy array of shape ((emax-emin)/(deltaE)+1,2)
		  The first column is the value of the energies generated using np.linspace(emin,emax,(emax-emin)/(deltaE)+1)
		  The second column is the value of the DOS

		numpy_save_opt specific params:
		  - pFile: (True/False) Enable/disable save functionality (default = True)
		  - fname: Output file name (must be present)
		  - fmt:   Format string to pass to np.savetxt

		numpy_plot_opt specific params:
		  - plot:      (True/False) Enable/disable plot functionality (default = True)
		  - xlab:      String to use as x label
		  - ylab:      String to use as y label
		  - start:     First column of the Y axis data to be plotted
		  - end:       Last column of the Y axis data to be plotted
		  - colors:    List of matplotlib color string to be used.
		               % is used to loop if (end-start > len(colors))
		  - labels:    List of strings to be used as labes.
		               No label is set if (end-start > len(labels))
		  - dash_list: List of tuples of dashes option to be used.
		               % is used to loop if (end-start > len(colors))
		               If no dash_list is specified, the lines will switch from nodash to dash=(8,2)
		               for every loop of the colors
		"""
		x = np.linspace(emin, emax, (emax-emin)/deltaE+1)
		y = np.zeros(x.size)

		for n,egv in enumerate(self.egv):
			for e in egv:
				index = int((e - emin) / deltaE)
				if 0 <= index < x.size:
					y[index] += self.weight[n]

		y /= deltaE

		data = np.vstack((x,y))
		if deg > 0:
			from ..tools.broad import broad
			data = broad(data, t='gauss', deg=deg, axis=1)

		return data.T

	@IO_stdout_redirect()
	def smallest_gap(self, radius=0., comp_point=(0.,0.,0.), **kwargs):
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

		# kpt   = np.array([a['kpt'] for a in self.kpt])
		kpt   = self.kpt_cart
		# egv   = np.array([ a['egv'] for a in self.egv])*self.e_units
		egv   = self.egv
		n_el  = self.n_el
		ef    = self.fermi
		# kpt   = kpt[:self.n_kpt,:]

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
			w = np.where((egv[:,vb] < ef) & (egv[:,cb] > ef))
			app_egv = egv[w,:]
			app_kpt = kpt[w,:]
			app_num = num[w,:]
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
		ret = True
		if self.n_kpt <= 0:
			warning.print("Failed to read nkpt from file '{}'.".format(self.schema))
			ret = False
			#raise Exception("No kpt read from file '{}'.".format(self.fname))
		if self.n_bnd <= 0:
			warning.print("Failed to read nbnd from file '{}'.".format(self.schema))
			ret = False
			#raise Exception("No band read from file '{}'.".format(self.fname))
		legv = len(self.egv)
		if self.occ:
			locc = len(self.occ)
		else:
			locc = legv
		if not self.n_kpt == legv == locc:
			warning.print("Corrupted file. Number of kpoints does not match number egv or occ")
			ret = False
			#raise Exception("Corrupted file. Number of kpoints does not match number egv or occ")
		return ret and super().validate()









