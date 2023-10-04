import numpy as np

from .._decorators import store_property
from ..parsers import Parser_xml

data={
	'header':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'PP_HEADER',
		'typ':list,
		},
	'mesh':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_MESH/PP_R',
		'typ':np.ndarray,
		},
	'rab':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_MESH/PP_RAB',
		'typ':np.ndarray,
		},
	'nlcc':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_NLCC',
		'typ':np.ndarray,
		},
	'nlcc':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_NLCC',
		'typ':np.ndarray,
		},
	'_semilocal':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'PP_SEMILOCAL', 
		'extra_name':'vnl', 
		'typ':list,
		},
	'local':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_LOCAL',
		'typ':np.ndarray,
		},
	'_nonloc':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'PP_NONLOCAL',
		'typ':list,
		},
	'_pswfc':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'PP_PSWFC//',
		'extra_name':'PP_CHI',
		'typ':list,
		},
	'rho_atom':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_RHOATOM',
		'typ':np.ndarray,
		},
	}

class UPF(Parser_xml):
	def __init__(self, xml_data={}, **kwargs):
		xml_data.update(data)
		super().__init__(xml_data=xml_data, **kwargs)

	@property
	@store_property
	def l_max(self):
		return self.header[0]['l_max']

	@property
	@store_property
	def l_local(self):
		# return getattr(self, 'local')
		return self.header[0]['l_local']

	@property
	@store_property
	def n_wfc(self):
		return self.header[0]['number_of_wfc']

	@property
	@store_property
	def n_proj(self):
		return self.header[0]['number_of_proj']

	@property
	@store_property
	def semilocal(self):
		try:
			return np.array([v for k,v in self._semilocal[0].items() if 'PP_VNL.' in k])
		except:
			return None

	@property
	@store_property
	def beta(self):
		try:
			return np.array([v for k,v in self._nonloc[0].items() if 'PP_BETA.' in k])
		except:
			return None

	@property
	@store_property
	def dij(self):
		try:
			res = self._nonloc[0].get('PP_DIJ', None)
			if not res is None:
				res = res.reshape(self.n_proj,self.n_proj)
			return res
		except:
			return None

	@property
	@store_property
	def pswfc(self):
		try:
			return np.array([np.fromstring(a, sep=' ') for a in self._pswfc])
		except:
			return None

	@property
	@store_property
	def pswfc_l(self):
		try:
			return np.array([a['l'] for a in self._pswfc])
		except Exception as e:
			print(e)
			return None

