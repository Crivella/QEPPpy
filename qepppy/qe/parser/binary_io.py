import sys
import struct
import numpy as np
from ...logger import logger


endianess = sys.byteorder
_endian = {
	'little':'<',
	'big':'>'
}

def _int( binary, n=1, endian="little"):
	res = struct.unpack( '{}{}l'.format( _endian[endianess], n), binary)
	if n == 1:
		return res[0]
	return res

def _dbl( binary, n=1, endian="little"):
	res = struct.unpack( '{}{}d'.format( _endian[endianess], n), binary)
	if n == 1:
		return res[0]
	return res

def _cpl( binary, n=1, endian="little"):
	res = struct.unpack( '{}{}d'.format( _endian[endianess], 2*n), binary)
	res = np.array( res)
	res.shape = (n, 2)
	res = res[:,0] + 1j * res[:,1]
	return res

wrap = {
	4:_int,
	8:_dbl,
	16:_cpl,
	}

@logger()
class binary_io():
	def __init__( self):
		self.binary = False
		return
		
	def read_binary( self, parse="", endian=endianess):
		with open( parse, "rb") as file:
			for vect in self.binary_format:
				rep = 1
				if isinstance( vect, tuple):
					rep = self._conv_val_( vect[1])
					vect = vect[0]
				for n in range( rep):
					size = _int( file.read(4))
					if not size:
						break
					content = file.read( size)
					start = 0
					end = 0
					for s in vect:
						name = s['n']
						if rep > 1 and n==0:
							self.__dict__[ name] = ['']*rep
						t = s['t']

						shape = self._convert_shape_( s['s'])
						chunk = self._get_chunk_( shape)

						end += chunk*t
						if end > size:
							raise Exception( "Binary file does not match the specified format.")
						splice = content[start:end]
						start = end

						res = wrap[t]( splice, n=chunk, endian=endianess)
						if not shape == (1,):
							res = np.array( res)
							res.shape = shape

						if rep == 1:
							self.__dict__[ name] = res
						else:
							self.__dict__[name][n] = res
					size = _int( file.read(4))
		self.binary = True
		return



	def _convert_shape_( self, tupl):
		n = len( tupl)
		l = [''] * n

		for i, e in enumerate( tupl):
			l[i] = self._conv_val_( e)

		return tuple( l)

	def _conv_val_( self, val):
		if isinstance( val, int):
			return val
		elif isinstance( val, str):
			a = self.__dict__.get( val, None)
			if a is None:
				raise Exception( "Failed to find preassigned variable {}.".format( val))
			return a
		raise Exception( "{} must be either int or str... Correct your code!!!".format( val))

	def _get_chunk_( self, shape):
		chunk = 1
		for d in shape:
			chunk *= d
		return chunk

