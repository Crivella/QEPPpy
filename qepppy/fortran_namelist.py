import re
from collections import OrderedDict

nl_param = r'(?P<param>[^\s\(\!]+)(\((?P<vec>[\+\d, ]+)\))?'
nl_value = r'(?P<value>('                                                             + \
		       r'((\'(?P<strip_single>[^\']*)\')|(\"(?P<strip_double>[^\"]*)\"))|'    + \
		       r'(\.true\.)|(\.false\.)|'                                             + \
		       r'([\+\-\.\d \t,]+)'                                                   + \
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

class fortran_namelist(OrderedDict):
	name = None
	def __init__(self, src=None, **kwargs):
		if isinstance(src, str):
			self.parse(src)
		super().__init__(**kwargs)
		
	def __getitem__(self, key):
		if isinstance(key, str):
			key = key.lower()
		return super().__getitem__(key)

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

		if self.name is None:
			return '\n\n'.join(a._format_output_(**kwargs) for a in self.values())  + '\n\n'
		return self._format_output_(**kwargs) + '\n\n'

	def _format_output_(self, tabs="  ", comma=",", upper_nl=True, al_p=0, al_v=0, al_1='>', al_2='<'):
		res = '&{}\n'.format(self.name.upper() if upper_nl else self.name.lower())
		for k,v in self.items():
			if isinstance(v, list):
				for i,sv in enumerate(v):
					res += '{0}{1:{4}{6}s} = {2:{5}{7}}{3}\n'.format(tabs, "{}({})".format(k,i+1), self._value_repr_(v[i]), 
						comma, al_1, al_2, al_p, al_v)
			else:
				res += '{0}{1:{4}{6}s} = {2:{5}{7}}{3}\n'.format(tabs, k, self._value_repr_(v), comma, al_1, al_2, al_p, al_v)
		res += '/'

		return res

	def max_length_param(self):
		if self.name is None:
			return max(a.max_length_param() for a in self.values()) #self.max_length_nl_param()
		if len(self):
			return max(len(a) + len(str(len(v))) + 2 if isinstance(v,list) else len(a) for a,v in self.items())
		return 0

	def max_length_value(self):
		if self.name is None:
			return max(a.max_length_value() for a in self.values()) #self.max_length_nl_value()
		if len(self):
			return max(len(self._value_repr_(a)) for a in self.values())
		return 0

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
		
	def parse(self, src):
		if isinstance(src, str):
			with open(src) as f:
				content = f.read()
		else:
			content = src.read()

		mid = [a.groupdict() for a in re.finditer(namelist_re, content)]

		for elem in mid:
			r = [a.groupdict() for a in re.finditer(nl_body_line, elem['body'])]
			new =  fortran_namelist()
			new.name = elem['name'].lower()
			for a in r:
				p = a['param']
				v = a['value']
				i = a['vec'] 
				if v is None:
					continue
				v = format_f90_to_py(v, not a['strip_single'] is None, not a['strip_double'] is None)
				if not i is None:
					self._set_vec_value_(new, p.lower(), v, i)
				else:
					new[p.lower()] = v
			self[elem['name'].lower()] = new

	def deep_find(self, pattern, up=None):
		if self.name is None:
			for elem in self.values():
				try:
					return elem.deep_find(pattern, up)
				except:
					pass
			raise

		tof_nl, tof_param, n  = self._tokenize_pattern_(pattern, up)
		if tof_nl is None or tof_nl == self.name:
			res =  self[tof_param]
			for i in n:
				res = res[i]
			return res

		raise

	@staticmethod
	def _tokenize_pattern_(pattern, up=None):
		if not isinstance(pattern, list):
			pattern = pattern.split("/")
			pattern = [None]*2 + pattern
			pattern = pattern[-2:]
		tom_major, tof_minor  = pattern
		if up:
			tom_major = up

		n = []
		if '(' in tof_minor:
			n = tof_minor.split('(')[1].split(')')[0]
			n = n.replace(" ", "").split(",")
			tof_minor = tof_minor.split('(')[0]

		return tom_major, tof_minor, n

	@staticmethod
	def _set_vec_value_(nl, param, value, index):
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
