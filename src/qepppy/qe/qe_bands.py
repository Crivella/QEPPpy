import numpy as np
from ..errors import ValidateError
# from .parser.data_file_parser import data_file_parser as dfp
from ..parsers import Parser_xml

HA_to_eV = 27.21138602


# data={
# 	'_n_kpt':{
# 		'xml_ptype':'text', 
# 		'xml_search_string':'output//nks', 
# 		'extra_name':None, 
# 		'res_type':int,
# 		'outfile_regex':r'number of k points[\s]*='
# 		},
# 	'_n_bnd':{
# 		'xml_ptype':'attr', 
# 		'xml_search_string':'output//ks_energies/eigenvalues', 
# 		'extra_name':'size', 
# 		'res_type':int,
# 		'outfile_regex':r'number of Kohn-Sham states[\s]*='
# 		},
# 	'_n_el':{
# 		'xml_ptype':'text', 
# 		'xml_search_string':'output//nelec', 
# 		'extra_name':None, 
# 		'res_type':float,
# 		'outfile_regex':r'number of electrons[\s]*='
# 		},
# 	'fermi':{
# 		'xml_ptype':'text', 
# 		'xml_search_string':'output//fermi_energy', 
# 		'extra_name':None, 
# 		'res_type':float,
# 		'outfile_regex':r'the Fermi energy is'
# 		},
# 	'fermi_s':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//two_fermi_energies', 
# 		'extra_name':'fermi', 
# 		'res_type':list
# 		},
# 	'homo':{
# 		'xml_ptype':'text', 
# 		'xml_search_string':'output//highestOccupiedLevel', 
# 		'extra_name':None, 
# 		'res_type':float
# 		},
# 	'lsda':{
# 		'xml_ptype':'text', 
# 		'xml_search_string':'output//lsda', 
# 		'extra_name':None, 
# 		'res_type':bool
# 		},
# 	'noncolin':{
# 		'xml_ptype':'text', 
# 		'xml_search_string':'output//noncolin', 
# 		'extra_name':None, 
# 		'res_type':bool,
# 		'outfile_regex':r'spin'
# 		},
# 	'_kpt':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//ks_energies/k_point', 
# 		'extra_name':'kpt', 
# 		'res_type':list,
# 		'outfile_regex':r'[\s]{4,}k\([ \d]+\) = \((?P<kpt>[ \d\.\-]+)\).*wk = (?P<weight>[ \d\.]+)'
# 		},
# 	'_app_egv':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//ks_energies/eigenvalues', 
# 		'extra_name':'egv', 
# 		'res_type':list,
# 		'outfile_regex':r'bands \(ev\):(?P<egv>[\s\d\.\-]+)', 
# 		'modifier':1/HA_to_eV
# 		},
# 	'_app_occ':{
# 		'xml_ptype':'nodelist', 
# 		'xml_search_string':'output//ks_energies/occupations', 
# 		'extra_name':'occ', 
# 		'res_type':list,
# 		'outfile_regex':r'occupation numbers(?P<occ>[\s\d\.]+)'
# 		},
# 	'_E_tot':{
# 		'xml_ptype':'text', 
# 		'xml_search_string':'output//total_energy/etot', 
# 		'extra_name':None, 
# 		'res_type':float,
# 		'outfile_regex':r'\!\s*total energy\s*='	
# 		}
# 	}

_data={
	'_n_kpt':{
		'xml_search_string':'output//nks',
		'typ':int,
		},
	'_n_bnd':{
		'xml_search_string':'output//nbnd',
		'typ':int,
		},
	'_n_el':{
		'xml_search_string':'output//nelec',
		'typ':float,
		},
	'fermi':{
		'xml_search_string':'output//fermi_energy',
		'typ':float,
		},
	'fermi_s':{ 
		'xml_search_string':'output//two_fermi_energies',
		'typ':np.ndarray,
		},
	'homo':{
		'xml_search_string':'output//highestOccupiedLevel',
		'typ':float,
		},
	'lsda':{
		'xml_search_string':'output//magnetization/lsda',
		'typ':bool,
		},
	'noncolin':{
		'xml_search_string':'output//magnetization/noncolin',
		'typ':bool,
		},
	'_weight':{
		'xml_search_string':'output//ks_energies/k_point',
		'mode':'attr=weight',
		'typ':np.ndarray
		},
	'kpt_cart':{
		'xml_search_string':'output//ks_energies/k_point', 
		'typ':np.ndarray,
		},
	'egv':{
		'xml_search_string':'output//ks_energies/eigenvalues', 
		'typ':np.ndarray,
		'modifier':lambda x: x*HA_to_eV
		},
	'occ':{
		'xml_search_string':'output//ks_energies/occupations', 
		'typ':np.ndarray
		},
	'_E_tot':{
		'xml_search_string':'output//total_energy/etot',
		'typ':float,
		}
	}

# @logger()
# class qe_bands(dfp):
class qe_bands(Parser_xml):
	"""
	Instance used for QE eigenvalues/vector(k-points) and occupations numbers.
	Uses the internal "data_file_parser" to read from a "data-file-schema.xml"
	or from a pw.x output file.
	Can be printed as a string.
	Each k-point and its info can be called as a dictionary value using its
	 number as the key.
	Provide the following PostProcessing methods:
	- band_structure(): Plot/print_to_file the band structure.
	- smallest_gap(): Print an analysis of the band gap.
	"""
	__name__ = "bands"
	e_units = HA_to_eV
	def __init__(self, data={}, **kwargs):
		data.update(_data)
		super().__init__(data=data, **kwargs)

	def __str__(self):
		msg = super().__str__()
		bnd = self._n_bnd
		kpt_fmt = "\nkpt(#{{:5d}}):  " + "{:8.4f}"*3 + " [2pi/alat]"
		egv_fmt = "\nEigenvalues(eV):\n" + ("  "+"{:12.6f}"*8+"\n")*(bnd//8)
		egv_fmt += "  " + "{:12.6f}"*(bnd%8) + "\n"
		for i in range(self._n_kpt):
			msg += kpt_fmt.format(*self.kpt_cart[i]).format(i)
			msg += egv_fmt.format(*self.egv[i])
		return msg

	def __getitem__(self, key):
		if(isinstance(key, int)):
			if(0 <= key < self._n_kpt):
				return {'kpt':self.kpt_cart[key], 'egv':self.egv[key], 'occ':self.occ[key]} 
			else:
				raise KeyError("Index '{}' out of range {}-{}".format(key, 0, self._n_kpt - 1))
		return super().__getitem__(key)

	@property
	def E_tot(self):
		"""Total energy"""
		if self._E_tot == 0.0:
			return None
		return self._E_tot

	@property
	def vb(self):
		"""Valence band index"""
		if self.noncolin:
			return int(self.n_el) - 1
		return int(self.n_el//2) - 1

	def validate(self):
		if self._n_kpt <= 0:
			raise ValidateError("Failed to read nkpt from file '{}'.".format(self.xml))
		if self._n_bnd <= 0:
			raise ValidateError("Failed to read nbnd from file '{}'.".format(self.xml))
		legv = self._egv.shape[0]
		if self.occ.size:
			locc = self._occ.shape[0]
		else:
			locc = legv
		# if not self.n_kpt == legv == locc:
		# 	raise ValidateError("Corrupted file. Number of k-points does not match number egv or occ {}/{}/{}".format(
		# 		self.n_kpt, legv, locc))
		super().validate()









