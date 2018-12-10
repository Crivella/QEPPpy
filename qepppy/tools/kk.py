#!/usr/bin/env python3

import numpy as np
from scipy.integrate import trapz


def kk_im( data, off=1, broad=1E-4):
	data = np.array(data)
	res = np.zeros( data.shape)
	res[:,0] = data[:,0]

	for n1 in range( 1, data.shape[1]):
		f = lambda o: trapz( data[:,n1] * data[:,0] / (data[:,0]**2 - o**2 + 1j*broad)*2/np.pi, data[:,0])
		res[:,n1] = np.array( [np.real( f(o)) + off for o in data[:,0]])
	return res

def kk_real( data, off=1, broad=1E-4):
	data = np.array(data)
	res = np.zeros( data.shape)
	res[:,0] = data[:,0]

	for n1 in range( 1, data.shape[1]):
		f = lambda o: trapz( (data[:,n1] - off) / (data[:,0]**2 - o**2 + 1j*broad)* -2/np.pi*o, data[:,0])
		res[:,n1] = np.array( [np.real( f(o)) for o in data[:,0]])
	return res

def kk( data, sf='imag', off=1, broad=1E-4):
	switch = {
		'imag':kk_im,
		'real':kk_real,
	}
	if not sf in switch:
		print( "Wrong operation type. Must be:", switch.keys)
		return

	return switch[sf]( data, off=off, broad=broad)


if __name__ == "__main__":
	import sys
	argc = len( sys.argv)
	if( not 2<=argc<=5):
		print("Incorrect use. Pleas pass arguments:"
			"\n\t'data_file'\t (),"
			"\n\t'start_from\t(optional) (imag/real, default=imag)',"
			"\n\t'off (add to final result in imag->real or subtract to starting data in real->imag)\t(optional) (default=1)'"
			"\n\t'broad\t(optional) (default=1E-4)'")
		exit()
	data = np.loadtxt( sys.argv[1])
	if( argc==2):
		res = kk( data)
	if( argc==3):
		res = kk( data, str( sys.argv[2]))
	if( argc==4):
		res = kk( data, str( sys.argv[2]), int( sys.argv[3]))
	if( argc==5):
		res = kk( data, str( sys.argv[2]), int( sys.argv[3]), float( sys.argv[4]))
	np.savetxt( 'kk.dat', res)


