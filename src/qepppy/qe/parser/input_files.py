import re
import itertools
from collections import OrderedDict

from ...errors  import ParseInputError, ValidateError
from ...parsers import fortran_namelist as f90nml

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

class qe_namelist_collection(f90nml.fortran_namelist_collection):
	def __init__(self, tpl=None, **kwargs):
		if tpl:
			self.load_templ(tpl)
		super().__init__(**kwargs)

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
		for name,nl in self.items():
			name = name.upper()
			if not name in self._templ_:
				raise ValidateError("Invalid namelist '{}'.".format(name))
			for k,v in nl.items():
				if not k in self._templ_[name]:
					raise ValidateError("\n\tIn namelist '{}' invalid param '{}'.".format(name,k))
				self._validate_tmpl_item(name, k, v)

		self._validate_default_()

	def _validate_default_(self):
		for name in self._templ_['nl']:
			for elem,v in self._templ_[name].items():
				if v['v'] == '***':
					raise ValidateError("\n\tMandatory param '{}/{}' is not set.".format(name, elem))

	def _validate_vec_(self, value, typ, lim):
		inf = lim[0]
		sup = lim[1]
		if isinstance(sup, str):
			sup = self.deep_find(sup)

		if isinstance(value, (int,float)):
			value = [value,]
		if len(value) <= sup-inf+1:
			return [tf90_to_py[typ](a) if not a is None else None for a in value]
		else:
			raise ValidateError("Too many elements.")

	def _validate_tmpl_item(self, namelist, param, value):
		if namelist is None:
			for n in self._templ_['namelist']:
				if param in self._templ_[n]:
					namelist = n
					break

		v = self._templ_[namelist][param]

		typ    = v['t']
		possib = v['c']
		vec    = v['vec']

		if not vec is None:
			try:
				v['v'] = self._validate_vec_(value, typ, vec)
			except Exception as e:
				raise ValidateError("\n\t{}/{}:  invalid vec value '{}'. ".format(namelist, param, value) + str(e))
		else:
			value = tf90_to_py[typ](value)
			if possib and all(not value in a for a in possib) and value != '':
				raise ValidateError("\n\t{}/{}: '{}' not among possibilities {}".format(namelist, param, value, possib))
			v['v'] = value

def conv_typ(typ, value):
	try:
		return typ(value)
	except:
		return None

class qe_card_syntax():
	data  = []
	nameless = False

	def __init__(self, nl, synt, nameless=False):
		self.namelist = nl
		self.nameless = nameless
		self.initialize(synt)

	def convert_text_block(self, body):
		res = OrderedDict()

		if not self.nameless:
			body = list(body)
		data = list(self.data)

		while data:
			l_s = data.pop(0)
			ml = 1
			mode = ''
			if isinstance(l_s[0], dict):
				mode    = l_s[0]['mode']
				hgh_lim = l_s[0]['hgh_lim']
				l_s     = l_s[0]['def_el'] + l_s[0]['opt_el']
				ml      = self.max_lines(res, hgh_lim)

			if 'row' in mode:
				for i in range(ml):
					l_e = body.pop(0) if body else ''
					l_e = list(filter(None,l_e.split(' ')))
					for app,value in itertools.zip_longest(l_s, l_e):
						if app is None:
							res[None] = 0
							continue
						name, typ = app
						if not name in res:
							res[name] = []
						res[name].append(conv_typ(typ, value))

			elif 'col' in mode:
				for (name,typ) in l_s:
					l_e = body.pop(0) if body else ''
					l_e = list(filter(None,l_e.split(' ')))
					if not name in res:
						res[name] = []
					res[name] = [conv_typ(typ, a) for a in l_e]
			else:
				l_e = body.pop(0) if body else ''
				l_e = list(filter(None,l_e.split(' ')))
				for (name,typ), value in itertools.zip_longest(l_s, l_e):
					res[name] = conv_typ(typ, value)
				
		if body and not self.nameless:
			print("WARNING: Ignoring lines:\n{}".format(body))

		return res

	def convert_dict(self, dic):
		res = ""

		for d in self.data:
			if isinstance(d[0], dict):
				tpl = d[0]
				v = tpl['def_el'] #+ tpl['opt_el']
				to_print = [dic[a[0]] for a in v if not all(b is None for b in dic[a[0]])]
				v = tpl['opt_el']
				if any(a[0] in dic for a in v):
					to_print += [dic[a[0]] for a in v if not all(b is None for b in dic[a[0]])]
				mode = tpl['mode']
				if 'row' in mode:
					app1 = to_print
					app2 = zip(*to_print)
					fmt = "  {{:{}s}}" * len(to_print)
				else:
					app1 = zip(*to_print)
					app2 = to_print
					fmt = "  {{:{}s}}" * len(to_print[0])

				max_l = [max(len(str(a)) for a in v) for v in app1]
				for v in app2:
					res += fmt.format(*max_l).format(*[str(a) for a in v]) + '\n'
			else:
				for name,typ in d:
					res += ' {}'.format(dic[name])
				res += '\n'

		return res


	def convert_item(self, key, value):
		for d in self.data:
			if isinstance(d[0], dict):
				for name,typ in d[0]['def_el'] + d[0]['opt_el']:
					if name == key:
						return [conv_typ(typ, a) for a in value]
			else:
				for name,typ in d:
					if name == key:
						try:
							return conv_typ(typ, value)
						except:
							return []

	def max_lines(self, dic, hgh_lim):
		lim = hgh_lim
		try:
			return int(lim)
		except:
			try:
				return dic[lim]
			except:
				return self.namelist.deep_find(lim)



	def initialize(self, synt):
		self.data  = []
		for elem in synt:
			if isinstance(elem, list):
				self.data.append([])
				ptr = self.data[-1]
				for v in elem:
					ptr.append((
						v['n'],
						tf90_to_py[v['t']],
						))

			if isinstance(elem, tuple):
				self.data.append([])
				ptr = self.data[-1]
				low_lim, hgh_lim, mode = elem[1:]

				def_el     = [
					(a['n'], tf90_to_py[a['t']])
					for a in elem[0] if isinstance(a, dict)
					]
				# def_el_num = len(self.def_el)
				opt_el = []
				for v in elem[0]:
					if isinstance(v, list):
						opt_el = [
							(a['n'], tf90_to_py[a['t']])
							for a in v if isinstance(a, dict)
							]
						# self.opt_el_num = len(self.opt_el)
				ptr.append({
					'low_lim':low_lim,
					'hgh_lim':hgh_lim,
					'def_el':def_el,
					'opt_el':opt_el,
					'mode':mode,
					})

	def validate(self, dic):
		if None in dic:
			raise ValidateError("Too many element in line.")

		if None in dic.values():
			raise ValidateError("Too few elements in line.")

		for d in self.data:
			if isinstance(d[0], dict):
				tpl = d[0]
				hgh_lim = tpl['hgh_lim']

				len_old = None
				for (name,typ) in tpl['def_el']:
					app = len(dic[name])
					if not len_old is None and app != len_old:
						raise ValidateError("Lenght mismatch of card parameters '{}'.".format(name))
					len_old = app
					ml = self.max_lines(dic, hgh_lim)
					if app != ml:
						raise ValidateError("Mismatch between number of params '{}': '{}' vs required '{}={}'".format(
							name, app, hgh_lim, ml))

					if None in dic[name]:
						raise ValidateError("Missing element in param '{}'.".format(name))

				app = []
				len_old = None
				if any(a[0] in dic for a in tpl['opt_el']):
					for (name,typ) in tpl['opt_el']:
						if not name in dic or (not len_old is None and len(dic[name]) != len_old):
							raise ValidateError("Lenght mismatch of card parameters.")

						len_old = len(dic[name])
						app += dic[name]

					if None in app and any(not a is None for a in app):
						raise ValidateError("Some optional parameters are set, but not all of them.")

				break


class qe_card_collection(OrderedDict):
	def __init__(self, nl, **kwargs):
		self._namelist = nl
		# super().__init__(**kwargs)
		for k,v in kwargs.items():
			if not k in self._templ_['card']:
				continue
			new = qe_card(self,v)
			new.name = k
			self[k]  = new

	def __str__(self):
		return '\n\n'.join(str(a) for a in self.values())

	def __getitem__(self, key):
		if not isinstance(key, str):
			raise ValueError("Key for {} must be of string type.".format(self))
		key = key.lower()
		try:
			return super().__getitem__(key)
		except KeyError:
			return self.deep_find(key)

	def __setitem__(self, key, value):
		if not isinstance(key, str):
			raise ValueError("Key for {} must be of string type.".format(self))
		super().__setitem__(key.lower(), value)

	def __contains__(self, key):
		if not isinstance(key, str):
			raise ValueError("Key for {} must be of string type.".format(repr(self)))
		return super().__contains__(key.lower())

	@property
	def _templ_(self):
		return self.namelist._templ_

	@property
	def namelist(self):
		return self._namelist

	@namelist.setter
	def namelist(self, value):
		if not isinstance(value, qe_namelist_collection):
			raise ValueError("Namelist element must be instance of qe_namelist_collection.")
		self._namelist = value

	def deep_find(self, pattern, up=None):
		try:
			return self.namelist.deep_find(pattern, up)
		except:
			tof_card, tof_param, n  = f90nml._tokenize_pattern_(pattern, up)
			if tof_param in self:
				return self[tof_param].value
			if tof_card:
				res = self[tof_card][tof_param]
			else:
				res = None
				for card in self.values():
					if pattern in card:
						res = card[tof_param]

			for i in n:
				res = res[int(i) - 1]

			return res

		raise

	def parse(self, src):
		if isinstance(src, str):
			with open(src) as f:
				content = f.read()
		else:
			content = src.read()

		if any('nameless' in self._templ_[c] for c in self._templ_['card']):
			app  = content.split('\n')
			app  = [a.split("#")[0] for a in app]
			app  = '\n'.join(app)
			body = app.split('/')[-1].strip().split('\n')
			for name in self._templ_['card']:
				new      = qe_card(self,nameless=True)
				new.name = name
				new.body = body
				new.assimilate()

				self[name] = new
				if not body:
					break
		else:
			card_split_re = '|'.join([r'({}.*\n)'.format(a) for a in self._templ_['card']])
			mid = list(filter(None,re.split(card_split_re,content)))[1:]

			for name,body in zip(mid[::2],mid[1::2]):
				app   = re.match(r'(?P<name>\S+)[ \t]*(?P<value>((\{[\w_]+\})|([\w_]+)))?', name).groupdict()
				name  = app['name']
				value = app['value']
				if not value is None:
					value = value.replace('{', '').replace('}', '')
					
				new       = qe_card(self)
				new.name  = name #.lower()
				new.value = value
				new.body  = body.replace('\t', ' ').rstrip().split('\n')
				new.assimilate()

				self[name] = new

	def validate(self):
		for card in self.values():
			card.validate()



class qe_card(OrderedDict):
	_parent  = None
	nameless = False
	value    = None

	def __init__(self, parent=None, nameless=False, **kwargs):
		self._parent = parent
		self.nameless = nameless
		super().__init__(**kwargs)

	def __str__(self):
		return self.format_output()

	def __setitem__(self, key, value):
		value = self.syntax.convert_item(key, value)
		super().__setitem__(key, value)

	@property
	def namelist(self):
		return self._parent.namelist

	@namelist.setter
	def namelist(self, value):
		self.parent.namelist = value

	@property
	def _templ_(self):
		return self.namelist._templ_

	@property
	def syntax(self):
		if not hasattr(self, '_syntax'):
			# card = self._templ_[self.name.upper()]
			card = self._templ_[self.name]
			val  = self.value
			if not val:
				val = card['d']

			for k,v in card.items():
				if not 'syntax' in k:
					continue
				if isinstance(v['cond'], str) and not val in v['cond']:
					continue
				self._syntax = qe_card_syntax(self.namelist, v['l'], nameless=self.nameless)
				break

		return self._syntax

	def assimilate(self):
		synt = self.syntax
		self.update(synt.convert_text_block(self.body))

	def format_output(self):
		res = ""
		if not self.nameless:
			res += self.name + (" {{{}}}".format(self.value) if not self.value is None  else '') + '\n'
		return res + self.syntax.convert_dict(self)

	def validate(self):
		try:
			self.syntax.validate(self)
		except ValidateError as e:
			raise ValidateError("ERROR in card '{}':\n".format(self.name), str(e))

class input_files():
	"""
	Class to handle any QE input (after loading the proper template).
	Supposed to be used as a parent class for child specific to the qe file.
	kwargs:
	 - input_file        = Name of the file to parse
	"""
	templ_file = None
	def __init__(self, input_file=None, **kwargs):
		if not self.templ_file:
			raise ParseInputError("Must give a template file.\n")
	
		self.namelist_c = qe_namelist_collection(tpl=self.templ_file, **kwargs)
		self.card_c     = qe_card_collection(self.namelist_c, **kwargs)
		if input_file:
			self.parse_input(input_file)

	def __getitem__(self, key):
		return self.card_c.__getitem__(key)

	def __setitem__(self, key, value):
		nl, name = ([None]*2 + key.split('/'))[-2:]
		key, i = (name.split('(') + [None]*2)[:2]
		if isinstance(i, str):
			i = int(i.replace(')', ''))
		try:
			if not nl in self.card_c:
				if nl.upper() in self.namelist_c._templ_['card']:
					new = qe_card(self.card_c)
					new.name = nl
					self.card_c[nl] = new
			if not i is None:
				self.card_c[nl][key+'({})'.format(i)] = value
			else:
				if key:
					self.card_c[nl][key] = value
				else:
					self.card_c[nl].value = value
		except:
			if not nl in self.namelist_c:
				if nl in self.namelist_c._templ_['nl']:
					new = f90nml.fortran_namelist()
					new.name = nl
					self.namelist_c[nl.lower()] = new
			self.namelist_c[nl].set_item(key, value, i)

	def __str__(self):
		return str(self.namelist_c) + str(self.card_c)

	def parse_input(self, input_file):
		with open(input_file) as f:
			self.namelist_c.parse(f)
		
		with open(input_file) as f:
			self.card_c.parse(f)

	def _find(self, *args, up=None):
		res = []
		for pattern in args:
			res.append(self.card_c.deep_find(pattern, up))
		if len(res) == 1:
			res = res[0]
		return res

	def validate(self):
		self.namelist_c.validate()
		self.card_c.validate()
		try:
			super().validate()
		except:
			pass







