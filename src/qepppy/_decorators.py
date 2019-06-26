import re
import functools
import numpy as np

def join_doc(func, add):
	tabs='\t'
	if not func.__doc__:
		func.__doc__ = ""

	app = func.__doc__
	func.__doc__ = ""
	for line in app.split('\n'):
		func.__doc__ += tabs + re.sub(r"^\t*", "", line) + "\n"

	for line in add.split("\n"):
		func.__doc__ += tabs + re.sub(r"^\t*", "", line) + "\n"

save_opt_doc = """
	numpy_save_opt specific params:
	  - pFile:     Enable/disable save functionality (default = True)
	  - fname:     Output file name (must be present)
	  - fmt:       Format string to pass to np.savetxt
	  - header:    Header for np.savetxt
	  - delimiter: Delimiter between coloumns"""

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

set_self_doc = """
	set_self specific params:
	 - set_self: If True set the variable as an attribute of the
	             calling class.
	             If False return the value."""

def numpy_save_opt(_fname='',_fmt='', _header='', _delimiter=' '):
	"""
	Decorator factory to add functionality to save return value to file 
	using np.savetxt
	Params:
	  - _fname:     Default save_file name to be used if not specified
	  - _fmt:       Default format string to be used if not specified
	  - _ header:   Default header string to be used if not specified
	  - _delimiter: Default delimiter to be used to separate coloumns
	"""
	def decorator(func):
		@functools.wraps(func)
		def wrapped(*args, pFile=True, fname=_fname, fmt=_fmt, header=_header, delimiter=_delimiter, **kwargs):
			res = func(*args, **kwargs)
			if not pFile:
				return res

			if not fname:
				raise ValueError("Must pass valid file name to arg 'fname'")
			save_args = {}
			if fmt:
				save_args['fmt'] = fmt
			if header:
				save_args['header'] = header
			save_args['delimiter'] = delimiter
			np.savetxt( fname=fname, X=res, **save_args)
			return res

		join_doc(wrapped, save_opt_doc)
		return wrapped
	return decorator

def numpy_plot_opt(_xlab='',_ylab='', _plot=True, _labels=['']):
	"""
	Decorator factory to add functionality to plot return value to file 
	using matplotlib.
	Params:
	  - _plot:   Default enable/disable plot (default = True)
	  - _xlab:   Default label for the x axis
	  - _ylab:   Default label for the x axis
	  - _labels: Default labels
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
			labels=_labels,
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

def set_self(_name, _default=True):
	"""
	Decorator factory to add functionality to return value, or set
	it as a class attribute.
	Params:
	 - _name:    Name of the attribute to set.
	 - _default: [True/False] Default behavior of set_self.
	"""
	def decorator(func):
		@functools.wraps(func)
		def wrapped(cls, *args, set_self=_default, **kwargs):
			res = func(cls, *args, **kwargs)
			if not set_self:
				return res

			names = _name.split(',')
			if len(names) == 1:
				setattr(cls, names[0], res)
			else:
				for name,val in zip(names,res):
					setattr(cls, name, val)

		join_doc(wrapped, set_self_doc)
		return wrapped
	return decorator

def store_property(func):	
	"""
	The first time the property value is accessed, generate it and store the 
	result in the object's __dict__ for future calls.
	"""
	@functools.wraps(func)
	def wrapped(*args, **kwargs):
		name = func.__name__
		cls=args[0]
		if name in cls.__dict__:
			return cls.__dict__[name]

		res = func(*args, **kwargs)
		cls.__dict__[name] = res
		return res
	return wrapped

def IO_stdout_redirect(_outfile=None):
	"""
	Decorator factory to catch and redirect stdout from a function call.
	If not output file is given, the stdout will not be redirected
	params:
	  - _outfile: Default output filename
	"""
	def decorator(func):	
		@functools.wraps(func)
		def wrapped(*args, outfile=_outfile, **kwargs):
			import sys
			from contextlib import redirect_stdout

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

		return wrapped
	return decorator

def IO_stderr_redirect(_errfile=None):
	"""
	Decorator factory to catch and redirect stderr from a function call.
	If not output file is given, the stdout will not be redirected
	params:
	  - _errfile: Default output filename
	"""
	def decorator(func):
		@functools.wraps(func)
		def wrapped(*args, errfile=_errfile, **kwargs):
			import sys
			from contextlib import redirect_stderr

			f = None
			if isinstance(errfile, str):
				f = open(errfile, "w")
			elif not errfile is None and hasattr(errfile, 'close') and not errfile is sys.stdout:
				f = errfile

			if f:
				with redirect_stderr(f):
					res = func(*args, **kwargs)
			else:
				res = func(*args, **kwargs)

			if f:
				f.close()

			return res

		return wrapped
	return decorator



def file_name_handle(
	_mode,
	):
	"""
	Decorator factory.
	Make it so that the first arg of a function or bound method can either be
	a file name or file handle.
	The decorated function should accept a file handle as its first argument.
	Params:
	 _mode: open mode for the file if the name is passed ['r','w',...]
	"""
	def decorator(func):
		@functools.wraps(func)
		def wrapped(*args, **kwargs):
			args = list(args)

			app = func
			if isinstance(args[0], object):
				cls  = args.pop(0)
				app = functools.partial(app, cls)

			file = args.pop(0)

			if hasattr(file, 'close'):
				return app(file, *args, **kwargs)

			f   = open(file, _mode)
			res = app(f, *args, **kwargs)
			f.close()

			return res

		return wrapped
	return decorator