import re
import numpy as np
from ..parsers import Parser_regex
# from ..logger import logger, warning
from .._decorators import numpy_save_opt, numpy_plot_opt, store_property, IO_stdout_redirect


data={
	'n_states':{
		'typ':int,
		'rstring':r'\s*natomwfc\s*='
		},
	'_n_bnd':{
		'typ':int,
		'rstring':r'\s*nbnd\s*='
		},
	'_states':{
		'typ':list,
		# state #   1: atom   1 (C  ), wfc  1 (l=0 m= 1)
		'rstring':
			r'state #\s*(?P<state_num>\d+)\s*\:\s*' + 
			r'atom\s*(?P<atom_num>\d+)\s*'          + 
			r'\(\s*(?P<atom_name>\S+)\s*\)\s*,'     +
			r'\s*wfc\s*(?P<wfc_num>\d+)\s*'         +
			# r'(\(l=\s*(?P<l>\d+)\s+m=\s*(?P<m>\d+)\s*\))' + 
			# r'|' + 
			r'(\(l=\s*(?P<lj>\d+)\s+j=\s*(?P<j>\S+)\s+m_j=\s*(?P<mj>\S+)\s*\))'
		},
	'_kpt':{
		'typ':list,
		'rstring':r'\s*k =\s*(?P<kpt>[\d.\- ]+)'
		},
	'_egv':{
		'typ':list,
		'rstring':
			r'\s*=*\s*e(' +
			r'( \= \s*)' + 
			r'|' + 
			r'(\(\s*(?P<egv_num>\d+)\s*\)\s*=\s*))' + 
			r'(?P<egv>[\d\-\.]+)\s*(?P<unit>\S+)\s*=*\s*' +
			r'psi = (?P<components>[\s\d\.\+\*#\[\]]+)' +
			r'\|psi\|\^2 = (?P<sum2>[\d\.]+)\s*'
			# r'\s*e(( \= \s*)|(\((?P<egv_num>\s*\d+)\)\s*=\s*))(?P<egv>[\d\-\.]+)\s*(?P<unit>\S+)\s*' +
			# r'psi = (?P<components>[\s\d\.\+\*#\[\]]+)\|psi\|\^2 = (?P<sum2>[\d\.]+)\s*'
		},
	}

class pdos(Parser_regex):
	__name__ = "pdos"
	def __init__(self, regex_data={}, **kwargs):
		regex_data.update(data)
		super().__init__(regex_data=regex_data, **kwargs)

	@property
	@store_property
	def n_kpt(self):
		return len(self._kpt)

	@property
	@store_property
	def n_bnd(self):
		if self._n_bnd != len(self._egv)//self.n_kpt:
			raise NotImplementedError()
		return self._n_bnd

	@property
	# @store_property
	def kpt(self):
		return self._kpt

	@property
	@store_property
	def egv(self):
		conversion = {
			'eV':1,
			}
		res = np.empty(self.n_kpt*self.n_bnd)
		for n,e in enumerate(self._egv):
			_, egv, unit, _, _ = e
			egv = float(egv)
			if not unit in conversion:
				raise NotImplementedError("Unit {} is not implemented".format(unit))
			res[n] = egv * conversion[unit]

		return res.reshape(self.n_kpt,self.n_bnd)

	@property
	@store_property
	def components(self):
		r = re.compile(r'\+?([\d.]+)\*\[\#\s*(\d+)\]')
		res = np.zeros((self.n_kpt*self.n_bnd, self.n_states))
		for n,e in enumerate(self._egv):
			_, _, _, components, _ = e
			if components.strip():
				comp = np.array([(float(a.group(1)),int(a.group(2))) for a in r.finditer(components)])
				res[n,comp[:,1].astype(dtype=int)-1] = comp[:,0]
		return res.reshape(self.n_kpt,self.n_bnd,self.n_states)

	@property
	@store_property
	def states(self):
		res = []
		m = max([len(a['atom_name']) for a in self._states])
		for e in self._states:
			msg = ""
			msg += "atom (#{:3d}): ".format(int(e['atom_num']))
			msg += "{1:>.{0}s} ".format(m,e['atom_name'])
			# msg += "(#{:5d}) ".format(int(e['atom_num']))
			if not e['l'] is None:
				msg += "(l = {}".format(e['l']) + "  "
				msg += "m = {})".format(e['m'])

			res.append(msg)
		if len(res) != self.n_states:
			raise NotImplementedError()
		return res

	@IO_stdout_redirect()
	def pdos_char(self, kpt_list=[], bnd_list=[], thr=1E-2, **kwargs):
		for k in kpt_list:
			print(("KPT(#{:5d}): " + "{:9.4f}"*3).format(k, *self.kpt[k-1]))
			for b in bnd_list:
				print("\tBND (#{:3d}): {} eV".format(b, self.egv[k-1][b-1]))
				for p in np.where(self.components[k-1,b-1,:] >= thr)[0]:
					print("\t\t{}: {:8.3f}%".format(self.states[p], self.components[k-1,b-1,p]*100))
		print()

	@numpy_plot_opt()
	@numpy_save_opt()
	def sum_pdos(
		self, *args,
		emin=-20, emax=20, deltaE=0.001, deg=0.00,
		weight=None,
		**kwargs
		):
		if weight is None:
			weight = np.ones(self.n_kpt) / self.n_kpt
		res = np.linspace(emin, emax, (emax-emin)/deltaE+1).reshape(1,-1)
		res = np.pad(res, ((0,self.n_states),(0,0)), 'constant')

		for k,egv in enumerate(self.egv):
			i = np.floor((egv - emin) / deltaE +0.5).astype(dtype='int')
			w = np.where( (0 <= i) & (i < res[0].size))[0]
			i = i[w]
			res[1:,i] += self.components[k,w,:].T * weight[k]

		res[1:,:] /= deltaE

		if deg > 0:
			from ..tools.broad import broad
			res = broad(res, t='gauss', deg=deg, axis=1)

		return res.T

