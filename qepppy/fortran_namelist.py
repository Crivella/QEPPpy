import re
from collections import OrderedDict

nl_param = r'(?P<param>[^\s\(\!]+)(\((?P<vec>[\+\d, ]+)\))?'
nl_value = r'(?P<value>('                                                             + \
		       r'((\'(?P<strip_single>[^\']*)\')|(\"(?P<strip_double>[^\"]*)\"))|'    + \
		       r'(\.true\.)|(\.false\.)|'                                             + \
		       r'([\+\-\.\d \t,]+)'                                                   + \
		     r'))'
nl_term  = r'(?P<term>[ \t]*,?[ \t]*(!.*\n)?)'
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

def inp_fmt(val, strip_s=False, strip_d=False):
	val = val.replace("\n", "")
	if strip_s:
		val = val.replace("'", "")
	if strip_d:
		val = val.replace('"', '')
	if strip_s or strip_d:
		return val
	if ',' in val:
		res = re.findall(r'[\+\-\.\d]+', val)
		res = [inp_fmt(a, strip_s, strip_d) for a in res]
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

def read(src):
	with open(src) as f:
		res = fortran_namelist()
		res.parse(f)
	return res


class fortran_namelist(OrderedDict):
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
			for a in r:
				p = a['param']
				v = a['value']
				i = a['vec'] 
				if v is None:
					continue
				v = inp_fmt(v, not a['strip_single'] is None, not a['strip_double'] is None)
				if not i is None:
					self._set_vec_value_(new, p.lower(), v, i)
				else:
					new[p.lower()] = v
			self[elem['name'].lower()] = new

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

	def deep_find(self, pattern, up=None):
		tof_nl, tof_param, n  = self._tokenize_pattern_(pattern, up)

		for nl,elem in self.items():
			if tof_nl is None or nl == tof_nl.lower():
				try:
					res = elem[tof_param]
					for i in n:
						res = res[i]
					return res
				except:
					pass

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
