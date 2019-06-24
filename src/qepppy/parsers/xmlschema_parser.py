import xmlschema

class Parser_xmlschema():
	def __init__(self, xml=None, schema=None, data={}):
		self.xml    = xml
		self.schema = schema
		self.data   = data

		if xml:
			self.xml_parse()
		if data:
			self.load_data()

	# def xml_parse(self, path):
	def xml_parse(self):
		self.xml_dict = xmlschema.to_dict(self.xml, schema=self.schema)
		# return xmlschema.to_dict(self.xml, schema=self.schema, path=path)

	# def find(self, path, mode='all'):
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


	def find(self, path, mode='all'):
		"""Xpath find function"""
		if not any(mode.startswith(a) for a in ['all','attr','value']):
			raise ValueError(f"Mode = {mode}  is not a valid value.")
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

		if mode.startswith('attr'):
			attr = mode.split('=')[1]
			ptr = self.get_attr(ptr, attr)
		elif mode.startswith('value'):
			ptr = self.get_value(ptr)

		return ptr

	@staticmethod
	def get_attr(l, attr):
		res = []
		for e in l:
			res.append(e['@' + attr])

		return res

	@staticmethod
	def get_value(l):
		res = []
		for e in l:
			res.append(e['$'])

		return res

	def load_data(self):
		for k,v in self.data.items():
			mode = v.get('mode', 'all')
			try:
				val = self.find(v['xml_search_string'], mode=mode)
			except Exception as e:
				raise type(e)(f'While finding {k}:' + str(e))
			setattr(self, k, val)








