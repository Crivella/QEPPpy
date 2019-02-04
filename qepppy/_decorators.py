import re
import functools
import numpy as np
from decorator import decorator

def join_doc(func, add):
	tabs='\t'
	if not func.__doc__:
		func.__doc__ = ""
	for line in add.split("\n"):
		func.__doc__ += tabs + re.sub(r"^\t*", "", line) + "\n"

save_opt_doc = """
	numpy_save_opt specific params:
	  - pFile: (True/False) Enable/disable save functionality (default = True)
	  - fname: Output file name (must be present)
	  - fmt:   Format string to pass to np.savetxt"""

plot_opt_doc = """
	numpy_plot_opt specific params:
	  - plot:      Enable/disable plot functionality (default = True)
	  - xlab:      String to use as x label
	  - ylab:      String to use as y label
	  - start:     First column of the Y axis data to be plotted
	  - end:       Last column of the Y axis data to be plotted
	  - colors:    List of matplotlib color string to be used.
	               % is used to loop if (end-start > len(colors))
	  - labels:    List of strings to be used as labes.
	               No label is set if (end-start > len(labels))
	  - dash_list: List of tuples of dashes option to be used.
	               % is used to loop if (end-start > len(colors))
	               If no dash_list is specified, the lines will switch from 
	               nodash to dash=(8,2) for every loop of the colors"""

def numpy_save_opt(_fname='',_fmt=''):
	"""
	Decorator factory to add functionality to save return value to file 
	using np.savetxt
	Params:
	  - _fname: Default save_file name to be used if not specified
	  - _fmt: Default format string to be used if not specified
	"""
	def decorator(func):
		@functools.wraps(func)
		def wrapped(*args, pFile=True, fname=_fname, fmt=_fmt, **kwargs):
			res = func(*args, **kwargs)
			if not pFile:
				return res

			if not fname:
				raise ValueError("Must pass valid file name to arg 'fname'")
			if fmt:
				np.savetxt( fname=fname, X=res, fmt=fmt)
			else:
				np.savetxt( fname=fname, X=res)
			return res

		join_doc(wrapped, save_opt_doc)
		return wrapped
	return decorator

def numpy_plot_opt(_xlab='',_ylab='', _plot=True):
	"""
	Decorator factory to add functionality to plot return value to file 
	using matplotlib.
	Params:
	  - _plot: Default enable/disable plot (default = True)
	  - _xlab: Default label for the x axis
	  - _ylab: Default label for the x axis
	"""
	def decorator(func):
		@functools.wraps(func)
		def wrapped(	
			*args,
			ax=None,
			plot=_plot,
			start=1, end=None,
			xlab=_xlab, ylab=_ylab,
			colors=['k','r','b','g','c','m'],
			labels=[''],
			dash_list=[],
			**kwargs
			):
			res = func(*args, **kwargs)
			if not plot:
				return res

			import matplotlib.pyplot as plt
			from matplotlib.ticker import AutoMinorLocator as AML
			to_plot = False
			if ax is None:
				offset = 0
				to_plot = True
				fig, ax = plt.subplots()
			else:
				offset = len(ax.get_lines())
			cl = len(colors)
			X = res[:,0]

			if res.shape[1] == 2:
				y_data = res[:,1].reshape(1,X.size)
			else:
				if not end is None:
					end += 1
				y_data = res[:,start:end].T

			for i,Y in enumerate(y_data):
				if dash_list:
					dash = dash_list[(i+offset)%len(dash_list)]
				else:
					dash = (8,2*(((i+offset)//cl)%2))
				ax.plot( 
					X, Y, 
					color=colors[(i+offset)%cl], 
					label=labels[i] if i<len(labels) else '',
					dashes=dash,
					)
			ax.set_xlabel(xlab)
			ax.set_ylabel(ylab)
			ml1 = AML(5)
			ax.yaxis.set_minor_locator(ml1)
			ax.yaxis.set_tick_params(which='both', right = True)
			if to_plot:
				fig.legend()
				plt.show()

			return res

		join_doc(wrapped, plot_opt_doc)
		return wrapped
	return decorator

@decorator
def store_property(func, *args, **kwargs):
	"""
	The first time the property value is accessed, generate it and store the 
	result in the object's __dict__ for future calls.
	"""
	name = func.__name__
	cls=args[0]
	try:
		res = cls.__dict__[name]
	except:
		res = func(*args, **kwargs)
		cls.__dict__[name] = res
	return res

@decorator
def IO_stdout_redirect(
	func,
	_outfile=None,
	*args,
	outfile=None,
	**kwargs
	):
	"""
	Decorator factory to catch and redirect stdout from a function call.
	If not output file is given, the stdout will not be redirected
	params:
	  - _outfile: Default output filename
	"""
	import sys
	from contextlib import redirect_stdout

	if outfile is None:
		outfile = _outfile

	f = None
	if isinstance(outfile, str):
		f = open(outfile, "w")
	elif not outfile is None and hasattr(outfile, 'close') and not outfile is sys.stdout:
		f = outfile

	if f:
		with redirect_stdout(f):
			res = func(*args, **kwargs)
	else:
		res = func(*args, **kwargs)

	if f:
		f.close()

	return res