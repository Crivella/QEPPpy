import re
import numpy as np
from ..parsers import Parser_regex
from ..calc_system import bands
# from ..logger import logger, warning
from .._decorators import numpy_save_opt, numpy_plot_opt, store_property, IO_stdout_redirect

from dataclasses import dataclass

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
	'kpt_cart':{
		'typ':list,
		'rstring':r'\s*k =\s*(?P<kpt>[\d.\- ]+)'
		},
	'_raw_egv':{
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

@dataclass(frozen=True)
class state:
	n: int
	atm: int
	atm_str: str
	wfc: int
	l: int
	j: float
	m: float


class pdos(Parser_regex, bands):
	__name__ = "pdos"
	def __init__(self, regex_data={}, **kwargs):
		regex_data.update(data)
		super().__init__(regex_data=regex_data, **kwargs)

	# @property
	# @store_property
	# def n_kpt(self):
	# 	return len(self._kpt)

	# @property
	# @store_property
	# def n_bnd(self):
	# 	# if self._n_bnd != len(self._egv)//self.n_kpt:
	# 	# 	raise NotImplementedError()
	# 	return self._n_bnd

	# @property
	# # @store_property
	# def kpt(self):
	# 	return self._kpt

	@property
	@store_property
	def _egv(self):
		conversion = {
			'eV':1,
			}
		res = np.empty((self.n_kpt,self.n_bnd))
		old = -1
		nkpt = 0
		for n,e in enumerate(self._raw_egv):
			nbnd, egv, unit, _, _ = e
			nbnd = int(float(nbnd))-1
			egv = float(egv)
			if not unit in conversion:
				raise NotImplementedError("Unit {} is not implemented".format(unit))
			if nbnd < old:
				nkpt +=1
			old = nbnd
			res[nkpt,nbnd] = egv * conversion[unit]

		return res.reshape(self.n_kpt,self.n_bnd)

	@property
	@store_property
	def components(self):
		r = re.compile(r'\+?([\d.]+)\*\[\#\s*(\d+)\]')
		res = np.zeros((self.n_kpt,self.n_bnd, self.n_states))
		old = -1
		nkpt = 0
		for n,e in enumerate(self._raw_egv):
			nbnd, _, _, components, _ = e
			nbnd = int(float(nbnd))-1
			if nbnd < old:
				nkpt +=1
			old = nbnd
			if components.strip():
				comp = np.array([(float(a.group(1)),int(a.group(2))) for a in r.finditer(components)])
				res[nkpt,nbnd,comp[:,1].astype(dtype=int)-1] = comp[:,0]
		return res.reshape(self.n_kpt,self.n_bnd,self.n_states)

	@property
	@store_property
	def states(self):
		res = []
		for e in self._states:
			n, a, astr, wfc, l, j, m = e

			new = state(
				int(float(n)),
				int(float(a)),
				astr,
				int(float(wfc)),
				int(float(l)),
				float(j),
				float(m),
			)
			res.append(new)
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

	@numpy_plot_opt(_xlab="Energy (eV)", _ylab="PDOS (arb. units)")
	@numpy_save_opt(_fname="pdos.dat")
	def sum_pdos(
		self, *args,
		emin=-20, emax=20, deltaE=0.001, deg=0.00,
		weight=None,
		**kwargs
		):
		if weight is None:
			weight = np.ones(self.n_kpt) / self.n_kpt
		res = np.linspace(emin, emax, int((emax-emin)/deltaE)+1).reshape(1,-1)
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

