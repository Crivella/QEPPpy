import re
import xmlschema
import numpy as np

def dbg(*args, **kwargs):
	return
	print(*args, **kwargs)

class Parser_xmlschema():
	def __init__(self, *args, xml=None, schema=None, data={}, **kwargs):
		self.xml    = xml
		self.schema = schema
		self.data   = data

		if xml:
			self.xml_parse()
			if data:
				self.load_data()

		super().__init__(*args, **kwargs)

	# def xml_parse(self, path):
	def xml_parse(self):
		self.xml_dict = xmlschema.to_dict(self.xml, schema=self.schema)
		# return xmlschema.to_dict(self.xml, schema=self.schema, path=path)

	# Much slower as xmlschema reparse the XMLTree for every find call
	# def find_native_xpath(self, path, mode='all'):
	# 	res =  xmlschema.to_dict(self.xml, schema=self.schema, path=path)

	# 	if mode.startswith('attr'):
	# 		attr = mode.split('=')[1]
	# 		res = self.get_attr(res, attr)
	# 	elif mode.startswith('value'):
	# 		res = self.get_value(res)

	# 	return res

	def _find_xpath(self, to_find, dct={}, deep=False):
		if to_find == '':
			return dct
		if not isinstance(dct, list):
			dct = [dct,]
		res = []
		for elem in dct:
			if not isinstance(elem, dict):
				continue
			if to_find in elem:
				app = elem[to_find]
				if isinstance(app, list):
					res += app
				else:
					res.append(elem[to_find])

			if not deep:
				continue
			for k,v in elem.items():
				if k.startswith('@') or k == '$':
					continue
				if isinstance(v, dict):
					res += self._find_xpath(to_find, v, deep)
				if isinstance(v, list):
					for e in v:
						res += self._find_xpath(to_find, e, deep)


		return res


	def find(self, path):
		"""Xpath find function"""
		dct       = self.xml_dict
		find_list = [a.split('/') for a in path.split('//')]

		ptr = dct
		deep = False
		for f1 in find_list:
			for f2 in f1:
				ptr = self._find_xpath(f2, ptr, deep)
				deep = False

			deep = True

		if len(ptr) == 1:
			ptr = ptr[0]

		return ptr

	@staticmethod
	def _get_dict_attr(dct, attr):
		res = []
		for k,v in dct.items():
			if not k.startswith('@'):
				continue
			if re.match('@' + attr, k):
				res.append(v)

		return res

	@staticmethod
	def get_attr(l, attr):
		# attr = attr.split(',')
		if not isinstance(l, list):
			if isinstance(l, dict):
				res = Parser_xmlschema._get_dict_attr(l, attr)
			else:
				raise ValueError("Unexpected behavior!!")
		else:
			res = []
			for e in l:
				if isinstance(e, dict):
					res += Parser_xmlschema._get_dict_attr(e, attr)

		if len(res) == 1:
			res = res[0]

		return res

	@staticmethod
	def _get_dict_value(dct):
		if '$' in dct:
			return [dct['$'],]
		res = []
		for k,v in dct.items():
			if isinstance(v, dict):
				continue
			if isinstance(v, list) and isinstance(v[0], dict):
				continue
			res.append(v)
		return res


	@staticmethod
	def get_value(l):
		if not isinstance(l, list):
			if isinstance(l, dict):
				return Parser_xmlschema._get_dict_value(l)
			else:
				raise ValueError("Unexpected behavior!!")

		res = []
		for e in l:
			res += Parser_xmlschema._get_dict_value(e)
			# res.append(e['$'])

		return res

	def load_data(self):
		for k,v in self.data.items():
			dbg("-"*40)
			dbg(f"Setting attribute '{k}'")
			mode    = v.get('mode',     'all')
			typ     = v.get('typ',      None)
			modifier = v.get('modifier', lambda x: x)

			if not any(mode.startswith(a) for a in ['all','attr','value']):
				raise ValueError(f"Mode = {mode}  is not a valid value.")

			try:
				val = self.find(v['xml_search_string'])
			except Exception as e:
				raise type(e)(f'While finding {k}:' + str(e))
			dbg(f'Found value: {val}')

			if mode.startswith('attr'):
				attr = mode.split('=')[1]
				val  = self.get_attr(val, attr)
			elif mode.startswith('value'):
				val  = self.get_value(val)

			if not typ is None:
				if typ == np.ndarray:
					val = np.array(val)
					if val.shape == ():
						val = val.reshape(-1)
				else:
					val = typ(val)
			val = modifier(val)
			dbg(f'Assigning:  {val}')
			# if hasattr(self, k):
			# 	print(f'WARNING self already has attribute "{k}" of type "{type(getattr(self,k))}')
			setattr(self, k, val)

	def validate(self):
		xmlschema.validate(self.xml, schema=self.schema)








