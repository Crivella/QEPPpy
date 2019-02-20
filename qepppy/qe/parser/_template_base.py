debug_templ = False #Enable convert of an empty template (Print all empty values. Put the max value for array var to 7 if not defined)
# from ...logger import logger, warning
from ...errors import ValidateError, ParseInputError

class templ_base(object):
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

	@staticmethod
	def _format_(v, f=0, a="'"):
		c = ""
		if isinstance(v, str):
			if f:
				a = "{2}{0:>{1}}{2}".format(v, f, a)
			else:
				a = "{1}{0}{1}".format(v, a)
			c += a
		if isinstance(v, bool):
			if v:
				c += ".TRUE."
			else:
				c += ".FALSE."
		else:
			if isinstance(v, (int, float)):
				if f:
					a = "{0:>{1}}".format(v, f)
				else:
					a = "{}".format(v)
				c += a.replace("e", "D").replace("E", "D")
		return c

	def _get_arr_ext_(self, sa, ea, n=-1):
		try:
			st= int(sa)
		except:
			st = self.find(sa)
		try:
			et= int(ea)
		except:
			et = self.find(ea)
		if not isinstance(et, int) and debug_templ:
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
			return
		raise ParseInputError("\n\tInvalid namelist: '{}'".format(nl))

	def convert(self):
		"""
		Convert the filled template in a QE input file and return it.
		Return instance: str
		"""
		content = super().convert()
		nl = self._templ_['nl'].copy()
		#Check for unused namelist (does not print it)
		for namelist in self._templ_['nl']:
			if not any(v['v'] for v in self._templ_[namelist].values()) and not debug_templ:
				nl.pop(nl.index(namelist))
		#Write all the used namlists/parameters
		longest = self._maxl_()
		for namelist in nl:
			content += "&{}\n".format(namelist)
			for el, v in self._templ_[namelist].items():
				if v['v'] != None and v['v'] != '' or debug_templ:
					if v['vec']:
						if debug_templ:
							try:
								end = self._get_arr_ext_(v['vec'][0], v['vec'][1])[1]
							except:
								end = 7
							if not v['v']:
								v['v'] = ['']*end
						for n, val in enumerate(v['v']):
							if not val and not debug_templ:
								continue
							app = el + "({})".format(n+1)
							content += "{0:>{1}} = ".format(app, longest,)
							content += self._format_(val)
							content += " ,\n"
					else:
						content += "{0:>{1}} = ".format(el, longest)
						content += self._format_(v['v'])
						content += " ,\n"
			content += "/\n\n"
		return content

	def set_nl(self, nl, k, v):
		"""
		Set the value "v" for param "k" in the NAMELIST "nl"
		"""
		n=None
		#print(nl, k, v)
		#If array case
		if '(' in k:
			n=int(k.split('(')[1].split(')')[0])
			k=str(k.split('(')[0])

		#Check if k is present in preset namelist
		#print (self._templ_[nl])
		# if not k in self._templ_[nl]: 
		# 	warning.print("Ignored unrecognized parameter '{}'".format(k))
		# 	return
		# 	#raise NameError("Ignored unrecognized parameter '{}'\n".format(k))
		ptr = self._templ_[nl][k]
		v = self._check_type_(v, ptr['t'])

		#Check value agains possible values
		if ptr['c']:
			if not any(v in opt for opt in ptr['c']):
				raise ParseInputError("\n\tParameter '{}/{}' = '{}' is not within range of possible values: \n\t\t{}".format(nl, k, v, ptr['c']))
				# warning.print("Parameter '{}/{}' = '{}' is not within range of possible values: \n{}".format(nl, k, v, ptr['c']))
				#return
		#If array case
		if n:
			if not isinstance(ptr['v'], list):
				raise ParseInputError("\n\t'{}' from namelist '{}' is not an array variable.".format(k, nl))
				# warning.print("'{}' from namelist '{}' is not an array variable.".format(k, nl))
				# return
				#raise Exception("'{}' from namelist '{}' is not an array variable.\n".format(k, nl))
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
		# ret = True
		for nl in self._templ_['nl']:
			for el,v in self._templ_[nl].items():
				val = v['v']
				if val == "***":
					raise ValidateError("\n\tRequired input parameter '{}' in namelist '{}' not set.\n".format(el, nl))
					# warning.print("Required input parameter '{}' in namelist '{}' not set.\n".format(el, nl))
					# ret = False
				elif val:
					if v['c']:
						if not any(val in opt for opt in v['c']):
							raise ValidateError("Parameter '{}/{}' = '{}' is not within range of possible values: \n{}\n".format(nl, el, val, v['c']))
							# warning.print("Parameter '{}/{}' = '{}' is not within range of possible values: \n{}\n".format(nl, el, val, v['c']))
							# ret = False
					if v['vec']:
						n = len(val)
						if not self._get_arr_ext_(v['vec'][0], v['vec'][1], n): 
							raise ValidateError("Parameter: {}({})".format(el, n))
							# warning.print("Parameter: {}({})".format(el, n))
							# ret = False
		super().validate()
		# return ret and super().validate()




	def _maxl_(self):
		#Get the longest element to print among all namelists
		#Uset to align all the '=' in the printed QE input file
		longest = 0
		for n in self._templ_['nl']:
			for el, v in self._templ_[n].items():
				if v['v'] != None and v['v'] != '' or debug_templ:
					app = len(el)
					if v['vec']:
						app += 4
					if longest < app:
						longest = app
		longest += 2
		return longest




# logger()()
class card(templ_base):
	def check_card(self, card):
		"""
		Check if a CARD name actually exist in the template
		"""
		if card in self._templ_['card']:
			return True
		return False

	def check_used(self, card):
		if self.check_card( card):
			return self._templ_[card]['u']
		else:
			return False

	def convert(self):
		"""
		Convert the filled template in a QE input file and return it.
		Return instance: str
		"""
		content = super().convert()

		def _conv_syntax_(l, lvl=0, arr_lvl=0):
			#Convert the syntax dictionary of the internal template
			c = ""
			if not isinstance(l, list):
				return ""
			if lvl == 2:
				if not any(a['v'] for a in l):
					return ""
			for e in l:
				#print(e)
				if isinstance(e, dict):
					#print(e['n'], " =?= ", name)
					v = e['v']
					if isinstance(v, list):
						try:
							v = v[arr_lvl]
						except:
							v = 'empty1'
					c += " "
					c += self._format_(v, 10, a="")
				if isinstance(e, list): 
					c += _conv_syntax_(e, lvl+1, arr_lvl)
				if isinstance(e, tuple):
					ext = self._get_arr_ext_(e[1], e[2])
					num = ext[1]-ext[0]+1
					if e[3] == 'rows':
						for n in range(num):
							c +=  _conv_syntax_(e[0], lvl+1, n)
					elif e[3] == 'cols': 
						for r in e[0]:
							for n in range(num):
								try:
									v = r['v'][n]
								except:
									v = 'empty2'
								c += self._format_(v, 10, a="")
							c += "\n"
			if lvl == 1:
				c += "\n"
			return c
		#Convert the cards
		cards = self._templ_['card'].copy()
		#Check for unused cards (does not print it)
		for card in self._templ_['card']:
			if not self._templ_[card]['u'] and not debug_templ:
				cards.pop(cards.index(card))

		for card in cards:
			ptr = self._templ_[card]
			content += card
			if ptr['v']:
				content += " {{{}}}".format(ptr['v'])
			content += "\n"

			l = self._get_syntax_(ptr)
			content += _conv_syntax_(l)
			content += "\n"
		return content

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
				for line in synt:
					self._validate_syntax_line_(line)
			elif c['r']:
				raise ValidateError("\n\tMandatory card '{}' is not set.".format(card))
		super().validate()

	@staticmethod
	def _validate_syntax_element_(elem, arr=None):
		val = elem['v']
		if val is None or val == '':
			raise ValidateError("\n\tUnset value for {}.".format(elem))
		if not arr is None and (len(val) != arr or any(a==None or a =='' for a in val)):
			raise ValidateError(
				"\n\tParam '{}: {}': Number of lines does not match specified value '{}'.".format(card, elem['n'], arr)
				)

	@staticmethod
	def _validate_syntax_optional_(elem):
		print(elem)
		# print([a['v'] for a in elem])
		if not any(a['v'] for a in elem):
			return
		card._validate_syntax_line_(elem)

	@staticmethod
	def _validate_syntax_line_(line, arr=None):
		for elem in line:
			if isinstance(elem, dict):
				card._validate_syntax_element_(elem)
			if isinstance(elem, list):
				card._validate_syntax_optional_(elem)
			if isinstance(elem, tuple):
				ext = card._get_arr_ext_(elem[1], elem[2])
				card._validate_syntax_line_(elem[0], ext[1]-ext[0]+1)

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

	@staticmethod
	def _get_syntax_(card):
		v = card['v']
		for k1, v1 in card.items():
			if not isinstance(v1, dict):
				continue
			if not 'syntax' in k1:
				continue
			if isinstance(v1['cond'], str):
				if not v in v1['cond']:
					continue
			return v1['l']
		raise ParseInputError("\n\tCannot find a parsed syntax for card: '{} {}'.".format(card, v))

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






