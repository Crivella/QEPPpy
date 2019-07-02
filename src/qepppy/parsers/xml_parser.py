import re
# import xmltodict
# import xmlschema
import numpy as np

# def dbg(*args, **kwargs):
# 	return
# 	print(*args, **kwargs)
# def dbg2(*args, **kwargs):
# 	return
# 	print(*args, **kwargs)

class xmltodict():
	@staticmethod
	def parse(file, attrib=True):
		from xml.etree import ElementTree as ET
		tree = ET.parse(file)
		root = tree.getroot()

		res = xmltodict._parse_childrens(root, attrib=attrib)
		if attrib:
			for k,v in root.attrib.items():
				res['@' + k] = v

		return res

	@staticmethod
	def _parse_childrens(node, attrib=True):
		res = {}
		for child in node:
			new = {}

			if attrib:
				for k,v in child.attrib.items():
					new['@' + k] = v

			text = child.text
			if isinstance(text, str):
				text = text.strip()
			if text:
				new['$'] = text

			new.update(xmltodict._parse_childrens(child, attrib=attrib))

			tag = child.tag
			if tag in res:
				if isinstance(res[tag], list):
					res[tag].append(new)
				else:
					res[tag] = [res[tag], new]
			else:
				res[tag] = new

		return res

class Parser_xml():
	def __init__(self, *args, xml=None, schema=None, xml_data={}, **kwargs):
		self.xml      = xml
		self.schema   = schema
		self.xml_data = xml_data

		if xml:
			self.xml_parse()

		super().__init__(*args, **kwargs)

	def xml_parse(self, xml=None):
		if not xml is None:
			self.xml = xml

		with open(self.xml, 'rb') as f:
			self.xml_dict = xmltodict.parse(f)

		if self.xml_data:
			self.load_data()

	def _find_xpath(self, to_find, dct={}, deep=False):
		# dbg2("    ", to_find, str(dct)[:80])
		if to_find == '':
			return dct
		if not isinstance(dct, list):
			dct = [dct,]
		res = []
		for elem in dct:
			# dbg2('       deep:', deep)
			# dbg2('       elem: ', str(elem)[:80])
			if not isinstance(elem, dict):
				# dbg2(str(elem)[:80])
				raise ValueError("Unexpected behavior. All elements should be dicts.")

			# dbg2('       keys: ', elem.keys())
			if to_find in elem:
				# dbg2("         FOUND!!!!!", str(elem[to_find])[:80])
				app = elem[to_find]
				if isinstance(app, list):
					res += app
				else:
					res.append(app)

			if not deep:
				continue
			for k,v in elem.items():
				if k.startswith('@') or k == '$':
					continue
				if isinstance(v, (dict,list)):
					res += self._find_xpath(to_find, v, deep)

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
		if not isinstance(l, list):
			if isinstance(l, dict):
				res = Parser_xml._get_dict_attr(l, attr)
			else:
				raise ValueError("Unexpected behavior!!")
		else:
			res = []
			for e in l:
				if isinstance(e, dict):
					res += Parser_xml._get_dict_attr(e, attr)

		if len(res) == 1:
			res = res[0]

		return res

	@staticmethod
	def _get_dict_value(dct):
		if '$' in dct:
			return [dct['$'],]
		if not isinstance(dct, dict):
			return [dct,]
		res = []
		for k,v in dct.items():
			if isinstance(v, dict):
				if '$' in v:
					res.append(v['$'])
				continue
			if isinstance(v, list) and isinstance(v[0], dict):
				continue
			res.append(v)

		return res


	@staticmethod
	def get_value(l):
		if not isinstance(l, list):
			if isinstance(l, dict):
				res = Parser_xml._get_dict_value(l)
			else:
				raise ValueError("Unexpected behavior!!")
		else:
			res = []
			for e in l:
				res += Parser_xml._get_dict_value(e)

		if len(res) == 1:
			res = res[0]

		return res

	@staticmethod
	def typ_conversion(val, typ):
		if typ == np.ndarray:
			if isinstance(val, str):
				val = list(filter(None, re.split(r'\s+', val)))
			elif isinstance(val, list):
				if len(val) > 0:
					val = [list(filter(None, re.split(r'\s+', a))) for a in val]
					val = [a[0] if len(a) == 1 else a for a in val]
			try:
				val = np.array(val, dtype=np.float)
			except ValueError:
				val = np.array(val)

			if val.shape == ():
				val = val.reshape(-1)
		elif typ in (int,float):
			if not val == []:
				val = typ(val)
		elif typ in (bool,):
			if val.upper() == 'FALSE':
				val = False
			elif val.upper() == 'TRUE':
				val = True
			else:
				val = bool(val)
		else:
			val = typ(val)

		return val

	def scale(self, val, fact):
		if not fact is None:
			if isinstance(fact, str):
				fact = getattr(self, fact)
			if len(val)>0 and isinstance(val.flatten()[0], (int,float,np.number)):
				val *= fact

		return val

	def load_data(self):
		for k,v in self.xml_data.items():
			if not 'xml_search_string' in v:
				continue
			# dbg("-"*40)
			# dbg(f"Setting attribute '{k}'")

			xml_string = v.get('xml_search_string')
			modes      = v.get('mode',     'value')
			typ        = v.get('typ',      None)
			# modifier = v.get('modifier', lambda x: x)
			scale_fact = v.get('xml_scale_fact', None)

			params = k.split(',')
			lmodes = modes.split(',')

			if len(params) != len(lmodes):
				raise ValueError("Must give the same number of assign names and modes!!")

			if not all(any(mode.startswith(a) for a in ['all','attr','value']) for mode in lmodes):
				raise ValueError(f"Mode = {modes}  is not a valid value.")

			try:
				val = self.find(xml_string)
			except Exception as e:
				raise type(e)(f'While finding {k}:' + str(e))
			# dbg(f'Found value: {val}')

			for num,(name,mode) in enumerate(zip(params,lmodes)):
				if mode.startswith('attr'):
					attr = mode.split('=')[1]
					app  = self.get_attr(val, attr)
				elif mode.startswith('value'):
					app  = self.get_value(val)

				# dbg(f'pretyp:', app)
				if not typ is None:
					app = self.typ_conversion(app, typ)

				# val = modifier(val)
				app = self.scale(app, scale_fact)
				# dbg(f'Assigning:  {app}')

				setattr(self, name, app)

	# def validate_schema(self):
	# 	xmlschema.validate(self.xml, schema=self.schema)






