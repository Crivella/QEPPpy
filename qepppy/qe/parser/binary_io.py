import sys
import numpy as np
from functools import reduce


endianess = sys.byteorder
_endian = {
	'little':'<',
	'big':'>'
}

class binary_io():
	def __init__(self):
		self.binary = False
		return
		
	def read_binary(self, src="", endian=endianess):
		endian = _endian[endian]
		with open(src, "rb") as file:
			for vect in self.binary_format:
				rep = 1
				if isinstance(vect, tuple):
					rep = self._conv_val_(vect[1])
					vect = vect[0]
				for n in range(rep):
					size = np.fromfile(file, '{}i4'.format(endian), 1) # _int(file.read(4))
					if not size:
						break
					for s in vect:
						name = s['name']
						if rep > 1 and n==0:
							self.__dict__[name] = ['']*rep
						t = s['type']
						shape = self._convert_shape_(s['shape'])
						chunk = reduce(lambda x,y: x*y, shape)
						print('{}{}'.format(endian, t))
						res = np.fromfile(file, '{}{}'.format(endian, t), chunk)
						if res.size == 1:
							res = res[0]
						if not shape == (1,):
							res = np.array(res).reshape(shape)

						if rep == 1:
							self.__dict__[name] = res
						else:
							self.__dict__[name][n] = res
					size = np.fromfile(file, '{}i4'.format(endian), 1)
		self.binary = True
		return



	def _convert_shape_(self, tupl):
		n = len(tupl)
		l = [''] * n

		for i, e in enumerate(tupl):
			l[i] = self._conv_val_(e)

		return tuple(l)

	def _conv_val_(self, val):
		if isinstance(val, int):
			return val
		elif isinstance(val, str):
			a = self.__dict__.get(val, None)
			if a is None:
				raise ValueError("Failed to find preassigned variable {}.".format(val))
			return a
		raise NotImplementedError("{} must be either int or str... Correct your code!!!".format(val))

