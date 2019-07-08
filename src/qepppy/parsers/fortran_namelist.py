import re
from collections import OrderedDict
from .._decorators import file_name_handle

nl_param = r'(?P<param>[^\s\(\!]+)(\((?P<vec>[\+\d, ]+)\))?'
nl_value = r'(?P<value>('                                                             + \
		       r'((\'(?P<strip_single>[^\']*)\')|(\"(?P<strip_double>[^\"]*)\"))|'    + \
		       r'(\.true\.)|(\.false\.)|'                                             + \
		       r'([\+\-\.\deEdD \t,]+)'                                                   + \
		     r'))'
nl_term  = r'(?P<term>[ \t]*,?[ \t]*(?P<inline_comment>!.*\n)?)'
nl_comm  = r'(?P<comment>\s*!.*\n)'
nl_body_line = r'('                   + \
			     r'(\s*'              + \
			       nl_param           + \
			       r'\s*=\s*'         + \
			       nl_value + nl_term + \
			     r')|'                + \
			       nl_comm            + \
			   r')'
namelist_re = r'\s*&(?P<name>\S+).*\n' + r'(?P<body>' + nl_body_line + r'*)' + r'\s*/'

def format_f90_to_py(val, strip_s=False, strip_d=False):
	val = val.replace("\n", "")
	if strip_s:
		val = val.replace("'", "")
	if strip_d:
		val = val.replace('"', '')
	if strip_s or strip_d:
		return val
	if ',' in val:
		res = re.findall(r'[\+\-\.\d]+', val)
		res = [format_f90_to_py(a, strip_s, strip_d) for a in res]
		if len(res) == 1:
			res = res[0]
		return res
	try:
		val = int(val)
	except:
		try:
			val = float(val)
		except:
			if val.lower() == '.false.':
				val = False
			elif val.lower() == '.true.':
				val = True
	return val

def _tokenize_pattern_(pattern, up=None):
	if not isinstance(pattern, list):
		pattern = ([None]*2 +  pattern.split("/"))[-2:]
	tof_major, tof_minor  = pattern
	if up:
		tof_major = up

	n = []
	if '(' in tof_minor:
		n = tof_minor.split('(')[1].split(')')[0]
		n = n.replace(" ", "").split(",")
		tof_minor = tof_minor.split('(')[0]

	return tof_major, tof_minor, n

class fortran_namelist(OrderedDict):
	def __init__(self, d=None, name=''):
		self.name = name.lower()
		if d is None:
			return
			
		if not isinstance(d, dict):
			raise TypeError("Must initialize fortran_namelist using a dictionary")

		for k,v in d.items():
			self.set_item(k,v)

	def __getitem__(self, key):
		if not isinstance(key, str):
			raise ValueError("Key for {} must be of string type.".format(self))
		return self.deep_find(key.lower())

	def __setitem__(self, key, value):
		if not isinstance(key, str):
			raise ValueError("Key for {} must be of string type.".format(self))
		return super().__setitem__(key.lower(), value)

	def __str__(self):
		return self.format_output()

	def  format_output(self, align_p=">=", align_v="<", **kwargs):
		d = {
			'al_p':0,
			'al_v':0,
			'al_1':'>',
			'al_2':'<'
			}

		if '=' in align_p:
			d['al_p'] = self.max_length_param()

		if '<' in align_p:
			d['al_1'] = '<'
		if '>' in align_v:
			d['al_v'] = self.max_length_value()
			d['al_2'] = '>'

		kwargs.update(d)

		# if self.name is None:
		# 	return '\n\n'.join(a._format_output_(**kwargs) for a in self.values())  + '\n\n'
		return self._format_output_(**kwargs) + '\n\n'

	def _format_output_(self, tabs="  ", comma=",", upper_nl=True, al_p=0, al_v=0, al_1='>', al_2='<'):
		res = '&{}\n'.format(self.name.upper() if upper_nl else self.name.lower())
		for k,v in self.items():
			if isinstance(v, list):
				for i,sv in enumerate(v):
					if sv is None:
						continue
					res += '{0}{1:{4}{6}s} = {2:{5}{7}}{3}\n'.format(tabs, "{}({})".format(k,i+1), self._value_repr_(v[i]), 
						comma, al_1, al_2, al_p, al_v)
			else:
				res += '{0}{1:{4}{6}s} = {2:{5}{7}}{3}\n'.format(tabs, k, self._value_repr_(v), comma, al_1, al_2, al_p, al_v)
		res += '/'

		return res

	def max_length_param(self):
		if len(self):
			return max(len(a) + len(str(len(v))) + 2 if isinstance(v,list) else len(a) for a,v in self.items())
		return 0

	def max_length_value(self):
		if len(self):
			return max(len(self._value_repr_(a)) for a in self.values())
		return 0

	def parse(self, src):
		if isinstance(src, dict):
			app = src
		else:
			if isinstance(src, str):
				with open(src) as f:
					content = f.read()
			else:
				content = src.read()

			app = [a.groupdict() for a in re.finditer(namelist_re, content, re.IGNORECASE)][0]

		r = [a.groupdict() for a in re.finditer(nl_body_line, app['body'], re.IGNORECASE)]
		self.name = app['name'].lower()
		for a in r:
			p = a['param']
			v = a['value']
			i = a['vec'] 
			if v is None:
				continue
			v = format_f90_to_py(v, not a['strip_single'] is None, not a['strip_double'] is None)
			self.set_item(p, v, i)


	def deep_find(self, pattern, up=None):
		if pattern in self:
			return super().__getitem__(pattern)
		tof_nl, tof_param, n  = _tokenize_pattern_(pattern, up)
		if tof_nl is None or tof_nl.lower() == self.name:
			res =  self[tof_param]
			for i in n:
				res = res[int(i)-1]
			return res

		raise

	def set_item(self, key, value, i=None):
		tof_nl, tof_param, n = _tokenize_pattern_(key)
		if not tof_nl is None:
			app = self[tof_nl]
		else:
			app = self

		if not n:
			n = i

		if not n:
			app[tof_param] = value
		else:
			app._set_vec_value_(tof_param.lower(), value, n)

	def _set_vec_value_(nl, param, value, index):
		index = str(index)
		index = re.findall(r'[\+\d]+', index)

		if not param in nl or not isinstance(nl[param], list):
			nl[param] = []
		ptr = nl[param]
		for n,i in enumerate(index):
			i = int(i)-1
			while len(ptr) < i+1:
				ptr.append(None)
			if n == len(index)-1:
				ptr[i] = value
			else:
				if not isinstance(ptr[i], list):
					ptr[i] = []
				ptr = ptr[i]

	@staticmethod
	def _value_repr_(value):
		if isinstance(value, str):
			return "'{}'".format(value)
		if isinstance(value, bool):
			if value:
				return '.TRUE.'
			return '.FALSE.'
		if isinstance(value, list):
			return ', '.join(fortran_namelist._value_repr_(a) for a in value)
		if value is None:
			return '0'
		return str(value)


class fortran_namelist_collection(OrderedDict):
	def __init__(self, src=None, input_data={}):
		if not src is None:
			self.parse(src)
		# super().__init__(**kwargs)
		for k,v in input_data.items():
			new = fortran_namelist(v, name=k)
			self[k] = new

	def __getitem__(self, key):
		if not isinstance(key, str):
			raise ValueError("Key for {} must be of string type.".format(self))
		return self.deep_find(key.lower())
		# return super().__getitem__(key.lower())

	def __setitem__(self, key, value):
		if not isinstance(value, fortran_namelist):
			raise ValueError("Value in {} must be of type {}.".format(type(self), type(fortran_namelist)))
		super().__setitem__(key.lower(), value)

	def __contains__(self, key):
		if not isinstance(key, str):
			raise ValueError("Key for {} must be of string type.".format(repr(self)))
		return super().__contains__(key.lower())

	def __str__(self):
		return self.format_output()

	def format_output(self, align_p=">=", align_v="<", **kwargs):
		d = {
			'al_p':0,
			'al_v':0,
			'al_1':'>',
			'al_2':'<'
			}

		if '=' in align_p:
			d['al_p'] = self.max_length_param()

		if '<' in align_p:
			d['al_1'] = '<'
		if '>' in align_v:
			d['al_v'] = self.max_length_value()
			d['al_2'] = '>'

		kwargs.update(d)

		return '\n\n'.join(a._format_output_(**kwargs) for a in self.values())  + '\n\n'

	def max_length_param(self):
		return max(a.max_length_param() for a in self.values())

	def max_length_value(self):
		return max(a.max_length_value() for a in self.values())

	@file_name_handle('r')
	def parse(self, src):
		content = src.read()

		mid = [a.groupdict() for a in re.finditer(namelist_re, content, re.IGNORECASE)]

		# print('CONTENT:', content)
		# print('-------------------\n', mid, '\n---------------------')
		for elem in mid:
			# print(elem)
			new =  fortran_namelist()
			new.parse(elem)

			self[elem['name']] = new

	def deep_find(self, pattern, up=None):
		# print("SEARCHING FOR:  ", pattern, up)
		if pattern in self:
			return super().__getitem__(pattern.lower())
		for elem in self.values():
			try:
				return elem.deep_find(pattern, up)
			except:
				pass
		raise



