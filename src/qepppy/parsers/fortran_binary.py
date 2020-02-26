import sys
import numpy as np
import scipy.io


endianess = sys.byteorder
_endian = {
	'little':'<',
	'big':'>'
}

class fortran_binary_io():
	def __init__(self, src="", endian=endianess, **kwargs):
		self.binary = False
		
		self.src = src
		if src:
			self.read_binary(src=src, endian=endian)

		super().__init__(**kwargs)
		
	def read_binary(self, src="", endian=endianess):
		endian = _endian[endian]
		with scipy.io.FortranFile(src) as f:
			for record in self.binary_format:
				rep = 1
				c = False
				if isinstance(record, tuple):
					rep = self._conv_val_(record[1])
					record = record[0]
					c = True
				names = [a['name'] for a in record]
				self.__dict__.update({n:[] for n in names})
				rec   = [
					endian + 
					str(self._convert_shape_(a['shape'])) + 
					a['type'] for a in record
					]
				for n in range(rep):
					res = f.read_record(*rec)
					if not isinstance(res, tuple):
						res = tuple((res,))
					for name,val in zip(names,res):
						self.__dict__[name].append(val)
				for n in names:
					self.__dict__[n] = np.array(self.__dict__[n])
					if not c:
						self.__dict__[n] = self.__dict__[n].squeeze()

		self.binary = True

	def _convert_shape_(self, tupl):
		l = []

		for e in tupl:
			l.append(self._conv_val_(e))

		return tuple(l)

	def _conv_val_(self, val):
		if isinstance(val, int):
			return val
		elif isinstance(val, str):
			return int(self.__dict__[val])
		raise ValueError("{} must be either int or str... Correct your code!!!".format(val))

