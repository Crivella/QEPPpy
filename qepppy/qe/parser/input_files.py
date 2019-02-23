import re
import os
import numpy as np
from collections import OrderedDict

from ...errors import ParseInputError, ValidateError
from ... import fortran_namelist as f90nml

tf90_to_py = {
	'INTEGER': int,
	'REAL': float,
	'LOGICAL': bool,
	'CHARACTER': str
	}

tf90_to_np = {
	'INTEGER': 'i4',
	'REAL': 'f8',
	'LOGICAL': 'bool',
	'CHARACTER': 'U64'
	}

class qe_namelist(f90nml.fortran_namelist):
	def load_templ(self, tpl=""):
		"""
		Load a QE template from a specified file (tpl) or the internal data
		"""
		import os
		if os.path.isfile(tpl):
			with open(tpl) as f:
				file = f.read()
		else:
			from pkg_resources import resource_string, resource_listdir
			if tpl in resource_listdir('qepppy.qe.parser.data', ''):
				file = resource_string('qepppy.qe.parser.data', tpl).decode('utf-8')

		import ast
		self._templ_ = ast.literal_eval(file)

	def validate(self):
		for name in self._templ_['nl']:
			for elem,v in self._templ_[name].items():
				try:
					comp = self.deep_find("{0}/{1}".format(name,elem))
				except:
					continue
				typ    = v['t']
				possib = v['c']
				vec    = v['vec']

				if not vec is None:
					try:
						v['v'] = self._validate_vec_(comp, typ, vec)
					except Exception as e:
						raise ValidateError("\n\t{}/{}:  invalid vec value '{}'. ".format(name,elem,comp) + str(e))
				else:
					comp = tf90_to_py[typ](comp)
					if possib and all(not comp in a for a in possib) and comp != '':
						raise ValidateError("\n\t{}/{}: '{}' not among possibilities {}".format(name,elem,comp,possib))
					v['v'] = comp

	def _validate_vec_(self, value, typ, lim):
		inf = lim[0]
		sup = lim[1]
		if isinstance(lim, str):
			sup = self.deep_find(sup)

		if isinstance(value, (int,float)):
			value = [value,]
		if len(value) <= sup-inf+1:
			return [tf90_to_py[typ](a) for a in value]
		else:
			raise ValidateError("Too many elements.")

class qe_card(OrderedDict):
	_namelist = qe_namelist()
	def __init__(self, nl, **kwargs):
		self.namelist = nl
		super().__init__(**kwargs)

	def __str__(self):
		return self.format_output()

	def __getitem__(self, value):
		try:
			return super().__getitem__(value)
		except KeyError:
			return self.deep_find(value)

	@property
	def namelist(self):
		return self._namelist

	@namelist.setter
	def namelist(self, value):
		if not isinstance(value, qe_namelist):
			raise ValueError("Namelist element must be instance of qe_namelist.")
		self._namelist = value

	@property
	def _templ_(self):
		return self.namelist._templ_

	def format_output(self):
		return '\n\n'.join(
			card + " " +( "{{{}}}".format(value) if not value is None else '') + '\n  ' + '\n  '.join(body) 
			for card,(value,body) in self.items()
			)

	def deep_find(self, pattern, up=None):
		try:
			return self.namelist.deep_find(pattern, up)
		except:
			tof_card, tof_param, n  = self.namelist._tokenize_pattern_(pattern, up)

			if tof_param in self._templ_['card']:
				return self._templ_[tof_param]['v']
			
			res = self._deep_find_cards_(tof_card, tof_param)

			for i in n:
				res = res[int(i)-1]
			return res

	def _deep_find_cards_(self, tof_card, tof_param):
		for card in self._templ_['card']:
			if tof_card and tof_card != card:
				continue
			synt = self._get_syntax_(card)
			try:
				return self._deep_find_syntax_(synt, tof_param)
			except:
				pass
		raise


	def _deep_find_syntax_(self, synt, tof_param):
		for e in synt:
			if isinstance(e, dict):
				if tof_param == e['n']:
					return e['v']
			if isinstance(e, list):
				return self._deep_find_syntax_(e, tof_param)
			if isinstance(e, tuple):
				return self._deep_find_syntax_(e[0], tof_param)
		raise


	def parse(self, src):
		if isinstance(src, str):
			with open(src) as f:
				content = f.read()
		else:
			content = src.read()

		card_split_re = '|'.join([r'({}.*\n)'.format(a) for a in self._templ_['card']])
		mid = list(filter(None,re.split(card_split_re,content)))[1:]

		for name,body in zip(mid[::2],mid[1::2]):
			name,value = re.match(r'(?P<name>\S+)[ \t]*(?P<value>((\{[\w_]+\})|([\w_]+)))?', name).groupdict().values()
			if not value is None:
				value = value.replace('{', '').replace('}', '')

			self[name] = (value, body.replace('\t', ' ').rstrip().split('\n'))

	def validate(self):
		# from io import StringIO
		for name in self._templ_['card']:
			ptr = self._templ_[name]
			if not name in self:
				if ptr['r']:
					raise ValidateError("Mandatory CARD '{}' not set.".format(name))
				continue
			val, body = self[name]

			ptr['u'] = True
			possib = ptr['c']
			if possib and all(not val in a for a in possib) and val != '':
				raise ValidateError("\n\t{}: '{}' not among possibilities {}".format(name,val,possib))
			ptr['v'] = val

			synt = self._get_syntax_(ptr, val)
			if synt is None:
				continue

			with open("tmp-card.123", "w") as f:
				f.write('\n'.join(body))
			with open("tmp-card.123") as f:
				try:
					n = self._validate_syntax_(synt, f)
				except ValidateError as e:
					raise ValidateError("{}: ".format(name) + str(e))

				remain = f.read()
				if remain:
					print("Ignoring:\n'''\n{}\n'''".format(remain))
					n = remain.count('\n') + 1
					for i in range(n):
						self[name][1].pop()

			os.remove("tmp-card.123")

	def _validate_syntax_(self, synt, data):
		for elem in synt:
			if isinstance(elem, dict):
				typ = tf90_to_np[elem['t']]
				extract = np.fromfile(data, dtype=typ, count=1, sep=' ')[0]
				if isinstance(elem['v'], list):
					elem['v'].append(extract)
				else:
					elem['v']=extract
			if isinstance(elem, list):
				self._validate_syntax_(elem, data)
			if isinstance(elem, tuple):
				try:
					max_lines = int(elem[2])
				except:
					max_lines = self.deep_find(elem[2])
				lines_tok = []
				for i in range(max_lines):
					r = data.readline()
					if not r:
						raise ValidateError("Too few lines.")
					lines_tok.append(list(filter(None, r.rstrip().split(" "))))


				extract = np.array(lines_tok)
				if elem[3] == 'cols':
					extract = extract.T
					if extract.shape[0] > max_lines:
						raise ValidateError("Too many elements on line. Expected {}".format(max_lines))
				self._validate_shape_(extract, elem)
				self._assign_tuple_(extract, elem)


	def _validate_shape_(self, extract, elem):
		types, opt_types = self._get_tuple_types_(elem)
		if extract.shape[1] > len(types):
			types += opt_types
			if extract.shape[1] > len(types):
				raise ValidateError("Too many elements on line. Expected {}".format(len(types)))

	@staticmethod
	def _get_tuple_types_(elem):
		types = [tf90_to_py[a['t']] for a in elem[0] if isinstance(a, dict)]
		opt_types = []
		for e in elem[0]:
			if isinstance(e, list):
				opt_types = [tf90_to_py[a['t']] for a in e]

		return types, opt_types

	def _assign_tuple_(self, extract, elem):
		for n1,e1 in enumerate(elem[0]):
			if isinstance(e1, dict):
				e1['v'] = extract[:,n1].astype(tf90_to_np[e1['t']])
			if isinstance(e1, list):
				if n1 == extract.shape[1]:
					return
				for n2,e2 in enumerate(e1):
					e2['v'] = extract[:,n1+n2].astype(tf90_to_np[e2['t']])


	def _get_syntax_(self, card, val=None):
		if isinstance(card, str):
			card = self._templ_[card]
		if not val:
			val = card['v']
		for k,v in card.items():
			if not 'syntax' in k:
				continue
			if isinstance(v['cond'], str) and not val in v['cond']:
				continue
			return v['l']

class input_files():
	"""
	Class to handle any QE input (after loading the proper template).
	Supposed to be used as a parent class for child specific to the qe file.
	kwargs:
	 - parse      = Name of the file to parse
	 - templ_file = Name of the template file to use
	"""
	templ_file = None
	def __init__(self, src):
		if not self.templ_file:
			raise ParseInputError("Must give a template file.\n")
	
		self.parse_input(src)

	def __getitem__(self, key):
		return self.card.__getitem__(key)

	def __str__(self):
		return str(self.namelist) + str(self.card)

	def parse_input(self, src):
		nl   = qe_namelist()
		with open(src) as f:
			nl.parse(f)
			nl.load_templ(self.templ_file)
		nl.validate()
		self.namelist = nl
		
		card = qe_card(nl)
		with open(src) as f:
			card.parse(f)
		card.validate()
		self.card = card

	def find(self, *args, up=None):
		res = []
		for pattern in args:
			res.append(self.card.deep_find(pattern, up))
		if len(res) == 1:
			res = res[0]
		return res







