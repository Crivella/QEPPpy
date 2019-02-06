#!/usr/bin/env python3
	
import numpy as np
from scipy.signal import hilbert
from qepppy._decorators import numpy_save_opt, numpy_plot_opt

def _kk_(data, im=1, start_offset=0, end_offset=0):
	n_pt = data.shape[0]
	data = data.reshape(n_pt,-1)
	if data[0,0] == 0.0:
		n_pt -= 1
		data = data[1:,:]
	if data.shape[1] > 1:
		x    = data[:,0]
		y    = data[:,1:] + start_offset
	else:
		x = np.arange(n_pt)
		y = data[:,0:] + start_offset

	s1 = 0
	if np.all(y[0] >= 0):
		s1 = n_pt
		y = np.vstack((im*y[::-1],y))
	y = np.pad(y, ((0,n_pt), (0,0)), 'edge')
	res  = im*np.imag(hilbert(y, axis=0)) + end_offset
	res  = res[s1:s1+n_pt]
	return np.column_stack((x,res))

@numpy_plot_opt(_plot=False)
@numpy_save_opt(_fname='kk_imag.dat')
def kk_eps_real2imag(data):
	"""
	Apply the Hilbert transformation (Kramers-Kronig) to the real part of the 
	dielectric function.
	1 is subtracted to the starting data to the result to get the imaginary part
	of the dielectric function.
	Params:
	 - data: An array where the first column is the x-axis data and all the 
	         other column are y-axis to which the Hilbert transform is applied.
	         If a an array of shape (1,n_pt) or (n_pt,) is given, it is treated 
	         as a y-axis and the x-axis is generated using np.arange(n_pt).
	""" 
	return _kk_(data, im=1, start_offset=-1)

@numpy_plot_opt(_plot=False)
@numpy_save_opt(_fname='kk_real.dat')
def kk_eps_imag2real(data):
	"""
	Apply the Hilbert transformation (Kramers-Kronig) to the imaginary part of 
	the dielectric function.
	1 is added to the result to get the real part of the dielectric function.
	Params:
	 - data: An array where the first column is the x-axis data and all the 
	         other column are y-axis to which the Hilbert transform is applied.
	         If a an array of shape (1,n_pt) or (n_pt,) is given, it is treated 
	         as a y-axis and the x-axis is generated using np.arange(n_pt).
	"""
	return _kk_(data, im=-1, end_offset=1)


def kk(data, sf='imag'):
	switch = {
		'imag':kk_eps_imag2real,
		'real':kk_eps_real2imag,
	}

	return switch[sf](data)


if __name__ == "__main__":
	import sys
	argc = len(sys.argv)
	if not 2<= argc <=5:
		print(
			"Incorrect use. Pleas pass arguments:"
			"\n\t'data_file'\t (),"
			"\n\t'start_from\t(optional) (imag/real, default=imag)',"
			)
		exit()
	data = np.loadtxt(sys.argv[1])
	if argc == 2:
		res = kk(data)
	if argc == 3:
		res = kk(data, str(sys.argv[2]))


