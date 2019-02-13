import numpy as np
from .parser.data_file_parser import data_file_parser as dfp
from .._decorators import store_property

data={
	'header':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'PP_HEADER',
		'res_type':list,
		},
	'mesh':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_MESH/PP_R',
		'res_type':np.array,
		},
	'rab':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_MESH/PP_RAB',
		'res_type':np.array,
		},
	'nlcc':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_NLCC',
		'res_type':np.array,
		},
	'nlcc':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_NLCC',
		'res_type':np.array,
		},
	'_semilocal':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'PP_SEMILOCAL', 
		'extra_name':'vnl', 
		'res_type':list,
		},
	'local':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_LOCAL',
		'res_type':np.array,
		},
	'_nonloc':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'PP_NONLOCAL',
		'res_type':list,
		},
	'_pswfc':{
		'xml_ptype':'nodelist', 
		'xml_search_string':'PP_PSWFC',
		'res_type':list,
		},
	'rho_atom':{
		'xml_ptype':'text', 
		'xml_search_string':'PP_RHOATOM',
		'res_type':np.array,
		},
	}

class UPF(dfp):
	def __init__(self, d={}, **kwargs):
		d.update(data)
		super().__init__(d=d, **kwargs)

	@property
	@store_property
	def mesh_size(self):
		return self.header[0]['mesh_size']

	@property
	@store_property
	def l_max(self):
		return self.header[0]['l_max']

	@property
	@store_property
	def l_local(self):
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
			return np.array([v for k,v in self._pswfc[0].items() if 'PP_CHI.' in k])
		except:
			return None

