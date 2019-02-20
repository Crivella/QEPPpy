from ...errors import ValidateError, ParseInputError


class templ_base(object):
	debug_templ = False
	def convert(self):
		return ""

	def validate(self):
		pass

	@staticmethod
	def _check_type_(v, type):
		type_check = { "CHARACTER":str, "INTEGER":int, "LOGICAL":bool, "REAL":float, "STRING":str}
		ty = type_check.get(str(type).upper())
		if ty == float or ty == int: 
			if isinstance(v, str):
				v=v.replace("D", "e").replace("d", "e")
		try:
			val = ty(v)
		except:
		 	# val = '***'
		 	raise ParseInputError("\n\tInvalid type:\n\t\tInput: {}\n\t\tExpec: {}".format(v, ty))
		 	# warning.print("Input: {}\nExpec: {}".format(v, ty))
		return val

	@property
	def formatter(self):
		return self._formatter

	@formatter.setter
	def formatter(self, fmt):
		from inspect import signature
		old_sig = signature(self._formatter)
		new_sig = signature(fmt)
		if old_sig == new_sig:
			self._formatter = fmt
		else:
			raise ValueError("Invalid formatter signature '{}'".format(new_sig))
	
	@staticmethod
	def _formatter(v, f=0, a="'"):
		if isinstance(v, str):
			return "{2}{0:>{1}}{2}".format(v, f, a)
		elif isinstance(v, bool):
			return ".{}.".format(v).upper()
		elif isinstance(v, (int, float)):
			return "{0:>{1}}".format(v, f).replace("e", "D").replace("E", "D")
		elif isinstance(v, list):
			return " ".join(["{:>10}".format(str(a)) for a in v])
		raise ParseInputError("Cannot apply formater to '{}'.".format(v))

	def _get_arr_ext_(self, sa, ea, n=-1):
		try:
			st= int(sa)
		except:
			st = self.find(sa)
		try:
			et= int(ea)
		except:
			et = self.find(ea)
		if not isinstance(et, int) and self.debug_templ:
			et = 7
		if any(not isinstance(a, int) for a in [st,et]):
			raise ParseInputError("\n\tFailed to find boundary '{}-{}'.".format(sa, ea))
			# warning.print("Failed to find boundary '{}-{}'.".format(sa, ea))
			# return
			#raise Exception("Failed to find boundary '{}-{}'.\n".format(sa, ea))
		if n >= 0:
			if not st <= n <= et:
				raise ParseInputError("\n\t'{}' out of array range '{}-{}'".format(n, st, et))
				# warning.print("'{}' out of array range '{}-{}'".format(n, st, et))
				# return
				#raise Exception("'{}' out of array range '{}-{}'".format(n, st, et))
		return (st, et)



class namelist(templ_base):
	def check_nl(self, nl):
		"""
		Check if NAMELIST name actually exist in template
		"""
		if nl in self._templ_['nl']:
			return True
		raise ParseInputError("\n\tInvalid namelist: '{}'".format(nl))

	def _get_used_namelists_(self):
		nl = self._templ_['nl'].copy()
		for namelist in self._templ_['nl']:
			if not any(v['v'] for v in self._templ_[namelist].values()) and not self.debug_templ:
				nl.pop(nl.index(namelist))
		return nl

	def convert(self):
		"""
		Convert the filled template in a QE input file and return it.
		Return instance: str
		"""
		res = super().convert()

		nl = self._get_used_namelists_()
		longest = self._maxl_()

		for namelist in nl:
			res += self._convert_namelist_(namelist, longest)
		return res

	def set_nl(self, nl, k, v):
		"""
		Set the value "v" for param "k" in the NAMELIST "nl"
		"""
		n = None
		#If array case
		if '(' in k:
			n = int(k.split('(')[1].split(')')[0])
			k = str(k.split('(')[0])

		ptr = self._templ_[nl][k]
		v   = self._check_type_(v, ptr['t'])

		#Check value agains possible values
		if ptr['c'] and not any(v in opt for opt in ptr['c']):
			raise ParseInputError("\n\tParameter '{}/{}' = '{}' is not within range of possible values: \n\t\t{}".format(nl, k, v, ptr['c']))
		#If array case
		if n:
			if not isinstance(ptr['v'], list):
				raise ParseInputError("\n\t'{}' from namelist '{}' is not an array variable.".format(k, nl))
			while len(ptr['v']) < n:
				ptr['v'].append('')
			ptr['v'][n-1] = v
		else:
			ptr['v'] = v

	def validate(self):
		"""
		Validate namelist after a read call. Checks:
			-all REQUIRED var have been set
			-all var type are compliant with their type
		"""
		for nl in self._templ_['nl']:
			self._validate_namelist_(nl)
		super().validate()

	def _convert_namelist_(self, nl, longest=0):
		res = "&{}\n".format(nl)
		for el, v in self._templ_[nl].items():
			if (v['v'] != None and v['v'] != '') or self.debug_templ:
				if v['vec']:
					res += self._convert_vec_(el, v, longest)
				else:
					res += "{0:>{1}} = ".format(el, longest)
					res += self.formatter(v['v'])
					res += " ,\n"
		return res + "/\n\n"

	def _convert_vec_(self, el, v, longest=0):
		res = ''
		if self.debug_templ:
			try:
				end = self._get_arr_ext_(v['vec'][0], v['vec'][1])[1]
			except:
				end = 7
			if not v['v']:
				v['v'] = ['']*end
		for n, val in enumerate(v['v']):
			if not val and not self.debug_templ:
				continue
			app = el + "({})".format(n+1)
			res += "{0:>{1}} = ".format(app, longest,)
			res += self.formatter(val)
			res += " ,\n"
		return res

	def _maxl_(self):
		#Get the longest element to print among all namelists
		#Uset to align all the '=' in the printed QE input file
		longest = 0
		for n in self._templ_['nl']:
			for el, v in self._templ_[n].items():
				if (v['v'] != None and v['v'] != '') or self.debug_templ:
					app = len(el)
					if v['vec']:
						app += 4
					longest = max((longest, app))
		return longest + 2

	def _validate_namelist_(self, nl):
		for el,v in self._templ_[nl].items():
			val = v['v']
			if val == "***":
				raise ValidateError("\n\tRequired input parameter '{}' in namelist '{}' not set.\n".format(el, nl))
			elif val:
				if v['c'] and not any(val in opt for opt in v['c']):
					raise ValidateError("Parameter '{}/{}' = '{}' is not within range of possible values: \n{}\n".format(nl, el, val, v['c']))
				if v['vec']:
					self._validate_vec_(el, v)

	def _validate_vec_(self, el, v):
		n = len(v['v'])
		if not self._get_arr_ext_(v['vec'][0], v['vec'][1], n): 
			raise ValidateError("Parameter: {}({})".format(el, n))



class card(templ_base):
	def check_card(self, card):
		"""
		Check if a CARD name actually exist in the template
		"""
		if card in self._templ_['card']:
			return True

	def check_used(self, card):
		if self.check_card(card):
			return self._templ_[card]['u']

	def convert(self):
		"""
		Convert the filled template in a QE input file and return it.
		Return instance: str
		"""
		res = super().convert()

		#Convert the cards
		cards = self._get_used_cards_()

		for card in cards:
			ptr = self._templ_[card]
			res += card
			if ptr['v']:
				res += " {{{}}}".format(ptr['v'])
			res += "\n"

			synt = self._get_syntax_(ptr)
			# res += _conv_syntax_(synt)
			res += self._convert_syntax_line_(synt) + "\n"
		return res

	def set_card(self, card, value="", line=""):
		"""
		Set the values inside a CARD.
		"v" sets the value supposed to be next the card name.
		"el" sets the card values line by line.
		"""
		ptr = self._templ_[card]
		synt= self._get_syntax_(ptr)
		ptr['u'] = True

		if value:
			if ptr['c'] and not value in ptr['c']:
				raise ParseInputError("\n\tInvalid value '{}' for card '{}' ({}).\n".format(value, card, ptr['c']))
			ptr['v'] = value
		if line:
			split_line = list(filter(None, line.strip().replace("\t", " ").split(" ")))
			self._set_syntax_line_(split_line, synt)



	def validate(self):
		"""
		Validate card after a read call. Checks:
			-Card value is valid
			-all lines are compliant with the syntax
		Return True if valid
		"""
		for card in self._templ_['card']:
			c = self._templ_[card]
			if c['u']:
				self._card_get_value_(c)

				synt = self._get_syntax_(c)
				self._validate_syntax_line_(synt)
			elif c['r']:
				raise ValidateError("\n\tMandatory card '{}' is not set.".format(card))
		super().validate()

	@staticmethod
	def _get_card_default_value_(card_elem):
		for opt in card_elem['c']:
			if opt in card_elem['d']:
				return opt
		raise ParseInputError("\nNo option for card '{}' and no default value either.".format(card_elem))

	@staticmethod
	def _card_get_value_(card_elem):
		val = card_elem['v']
		if not val:
			if card_elem['c']:
				val = card._get_card_default_value_(card_elem)
		return val

	# @staticmethod
	def _convert_syntax_line_(self, synt, endl=''):
		res = ""
		for elem in synt:
			if isinstance(elem,dict):
				res += '  ' + self.formatter(elem['v'], a='') + endl
			elif isinstance(elem,list):
				res += self._convert_syntax_line_(elem)
			elif isinstance(elem,tuple):
				if elem[3] == 'rows':
					ext = self._get_arr_ext_(elem[1], elem[2])
					res += self._convert_syntax_rows_(elem[0], elem[3], ext[1]-ext[0]+1) + '\n'
				else:
					res += self._convert_syntax_line_(elem[0], endl='\n')

		return res + '\n'

	def _convert_syntax_rows_(self, synt, mode, arr):
		res = ''
		for i in range(arr):
			for elem in synt:
				if isinstance(elem, list):
					if all(e['v'] == [] for e in elem):
						continue
					for e in elem:
						m = max(len(str(a)) for a in e['v'])
						res += self.formatter(e['v'][i], m+4, '')
				else:
					m = max(len(str(a)) for a in elem['v'])
					res += self.formatter(elem['v'][i], m+4, '')
			res += '\n'
		return res

	@staticmethod
	def _get_syntax_(card):
		v = card['v']
		for k1, v1 in card.items():
			if not isinstance(v1, dict):
				continue
			if not 'syntax' in k1:
				continue
			if isinstance(v1['cond'], str) and not v in v1['cond']:
				continue
			return v1['l']
		raise ParseInputError("\n\tCannot find a parsed syntax for card: '{} {}'.".format(card, v))

	def _get_used_cards_(self):
		cards = self._templ_['card'].copy()
		for card in self._templ_['card']:
			if not self._templ_[card]['u'] and not self.debug_templ:
				cards.pop(cards.index(card))
		return cards

	@staticmethod
	def _set_syntax_line_(line, synt):
		for n,(val,elem) in enumerate(zip(line,synt)):
			if isinstance(elem,dict):
				card._set_sytanx_elem_(elem, val)
			elif isinstance(elem,list):
				return card._set_syntax_line_(line[n:], elem)
			elif isinstance(elem,tuple):
				if elem[3] == 'cols':
					return card._set_syntax_cols_(line,elem[0])
				return card._set_syntax_line_(line, elem[0])
		if len(line[n+1:]) > 0:
			raise ParseInputError("\n\tLine containes too many elements: {}".format(line))

	@staticmethod
	def _set_sytanx_elem_(elem, val):
		try:
			elem['v'].append(val)
		except AttributeError:
			elem['v'] = val

	@staticmethod
	def _set_syntax_cols_(line, synt):
		for elem in synt:
			if elem['v']:
				continue
			for val in line:
				val = card._check_type_(val, elem['t'])
				elem['v'].append(val)
			return
		raise ParseInputError("\n\tToo many lines")

	@staticmethod
	def _validate_syntax_element_(elem, arr=None):
		val = elem['v']
		if val is None or val == '':
			raise ValidateError("\n\tUnset value for {}.".format(elem))
		if not arr is None and (len(val) != arr or any(a==None or a =='' for a in val)):
			raise ValidateError(
				"\n\tNumber of lines/elements in '{}' does not match specified value '{}'.".format(elem['n'], arr)
				)

	# @staticmethod
	def _validate_syntax_optional_(self, elem, arr):
		if not any(len(a['v']) for a in elem):
			return
		if not any(len(a['v']) != arr for a in elem):
			return
		self._validate_syntax_line_(elem, arr)

	# @staticmethod
	def _validate_syntax_line_(self, synt, arr=None):
		for line in synt:
			if isinstance(line, dict):
				card._validate_syntax_element_(line, arr)
			if isinstance(line, list):
				self._validate_syntax_optional_(line, arr)
			if isinstance(line, tuple):
				ext = self._get_arr_ext_(line[1], line[2])
				self._validate_syntax_line_(line[0], ext[1]-ext[0]+1)






