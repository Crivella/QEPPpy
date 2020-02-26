import re
import numpy as np

matches = {
	int:   r'[\s]*(?P<flag>[\d\-]+)',
	float: r'[\s]*(?P<flag>[\d\.\-EedD]+)',
	str:   r'[\s]*(?P<flag>.*)',
	bool:  r'[\s]*(?P<flag>.)',
}

# def dbg1(*args, **kwargs):
# 	# return
# 	print(*args, **kwargs)

class Parser_regex():
	def __init__(self, *args, file=None, regex_data={}, **kwargs):
		super().__init__(*args, **kwargs)

		self.regex_data = regex_data

		if file:
			self.parse_file(file)

	@property
	def regex_data(self):
		return self._regex_data

	@regex_data.setter
	def regex_data(self, value):
		if not isinstance(value, dict):
			raise ValueError("Value must be a dictionary.")

		self._regex_data = value

	@staticmethod
	def get_single_val(content, rstring, typ):
		res = None

		if not 'flag' in rstring:
			rstring += matches[typ]
		val = re.search(rstring, content)
		if not val is None:
			res = typ(val.group('flag'))

		return res

	@staticmethod
	def get_list_val(content, rstring, typ):
		res = re.finditer(rstring, content)
		res = [x.groupdict() for x in res]
		for n,e in enumerate(res):
			for k,v in e.items():
				b = np.fromstring(v, sep=' ')
				if len(b) == 0 or re.findall(r'\s[a-zA-Z]', v):
					b = str(v).strip()
				elif len(b) == 1:
					b = b[0]
				res[n][k] = b

		return res

	# @staticmethod
	def max_num_truncate(self, max_num, val):
		fact = 1
		if not max_num is None:
			if isinstance(max_num, str):
				if max_num.startswith('-'):
					fact = -1
					max_num = max_num[1:]
				max_num = getattr(self, max_num.strip()) * fact
			if max_num > 0:
				val = val[:max_num]
			else:
				val = val[max_num:]

		return val

	def scale(self, val, fact):
		if not fact is None:
			if isinstance(fact, str):
				fact = getattr(self, fact)
			if len(val)>0 and isinstance(val.flatten()[0], (int,float,np.number)):
				val *= fact

		return val

	

	def parse_file(self, file):
		# res = {}

		with open(file, 'r') as f:
			content = f.read()

		for k,v in self.regex_data.items():
			if not 'rstring' in v:
				continue
			# dbg1('-'*40)
			# dbg1('Key:', k)
			# dbg1('val:', v)

			rstring    = v.get('rstring')
			typ        = v.get('typ')
			max_num    = v.get('max_num', None)
			scale_fact = v.get('re_scale_fact', None)

			if typ in matches:
				val = self.get_single_val(content, rstring, typ)
				setattr(self, k, val)
				# dbg1('found_single:', val)
			elif typ in (list,np.ndarray):
				app = self.get_list_val(content, rstring, typ)
				params = k.split(',')
				# dbg1(params)
				# dbg1(app)
				for num,name in enumerate(params):
					if name == '_':
						continue
					n_elem  = len(app)
					n_group = len(app[0].values()) if n_elem > 0 else 0
					if len(params) == 1 and not (n_elem > 0 and n_group == 1):
						app2 = np.array([list(a.values()) for a in app])
					else:
						app2 = np.array([list(a.values())[num] for a in app])




					val = self.max_num_truncate(max_num, app2)
					val = self.scale(val, scale_fact)

					# dbg1(f'found_mul({name}):', val)
					setattr(self, name, val)
			else:
				raise NotImplementedError()

