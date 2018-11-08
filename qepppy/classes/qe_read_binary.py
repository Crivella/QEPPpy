import struct
import numpy as np


endianess = 'little'
_endian = {
	'little':'<',
	'big':'>'
}

def _int( binary, n=1):
	#print( '{}{}l'.format( _endian[endianess], n))
	res = struct.unpack( '{}{}l'.format( _endian[endianess], n), binary)
	if n == 1:
		return res[0]
	return res

def _dbl( binary, n=1):
	#print( '{}{}d'.format( _endian[endianess], n))
	res = struct.unpack( '{}{}d'.format( _endian[endianess], n), binary)
	if n == 1:
		return res[0]
	return res

wrap = {
	4:_int,
	8:_dbl,
	}

wfc_format =[
		[
			{'t':4, 's':(1,), 'n':'nkpt'},
			{'t':8, 's':(3,), 'n':'kpt'},
			{'t':4, 's':(1,), 'n':'ispin'},
			{'t':4, 's':(1,), 'n':'gamma_only'},
			{'t':8, 's':(1,), 'n':'scale_factor'},
		],
		[
			{'t':4, 's':(1,), 'n':'max_index'},
			{'t':4, 's':(1,), 'n':'igwx'},
			{'t':4, 's':(1,), 'n':'nspin'},
			{'t':4, 's':(1,), 'n':'nbnd'},
		],
		[
			{'t':8, 's':(3,3), 'n':'recipr'},
		],
		[
			{'t':4, 's':('igwx',3,), 'n':'gvect'},
		],
		([
			{'t':8, 's':('igwx',2,), 'n':'val'},
		], 'nbnd'),
	]

class qe_binary_reader():
	def __init__( self):
		return

	def binary_read( self, parse=""):
		with open( parse, "rb") as file:
			for vect in wfc_format:
				rep = 1
				if isinstance( vect, tuple):
					rep = self._conv_val_( vect[1])
					vect = vect[0]
				#print( rep)
				for n in range( rep):
					size = _int( file.read(4))
					if not size:
						break
					content = file.read( size)
					start = 0
					end = 0
					for s in vect:
						name = s['n']
						t = s['t']
						shape = self._convert_shape_( s['s'])

						chunk = 1
						for d in shape:
							chunk *= d

						end += chunk*t
						if end > size:
							raise Exception( "Binary file does not match the specified format.")

						splice = content[start:end]
						res = wrap[t]( splice, n=chunk)
						start = end
						if not shape == (1,):
							res = np.array( res)
							res.shape = shape

						if not name in self.__dict__:
							self.__dict__[ name] = res
						else:
							self.__dict__[name] = np.vstack( ( np.array( self.__dict__[name]), np.array( res)))
					size = _int( file.read(4))

		print( "\n\n\n\n", self.val)
		print( self.val.shape)

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




a = qe_binary_reader( )
a.binary_read( parse="wfc1.dat")