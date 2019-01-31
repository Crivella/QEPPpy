import numpy as np
from decorator import decorator


@decorator
def save_opt(
	func, _fname='', *args, 
	pFile=True,
	fname="", fmt="", 
	**kwargs
	):
	"""
	Added functionality: save return value to file using np.savetxt
	params:
	  - pFile: (True/False) Enable/disable save functionality (default = True)
	  - fname: output file name (must be present)
	  - fmt: format string to pass to np.savetxt
	"""
	# try:
	# 	func.doc_save_opt
	# except:
	# 	if func.__doc__ is None:
	# 		func.__doc__ = ''
	# 	func.__doc__ += '\n' + save_opt.__doc__
	# 	func.doc_save_opt = True

	if not fname:
		fname = _fname

	res = func(*args, **kwargs)
	if pFile:
		if not fname:
			raise ValueError("Must pass valid file name to arg 'fname'")
		if fmt:
			np.savetxt( fname=fname, X=res, fmt=fmt)
		else:
			np.savetxt( fname=fname, X=res)
	return res

@decorator
def plot_opt(
	func, *args, 
	plot=True, 
	xlab="", ylab="", 
	**kwargs
	):
	"""
	Added functionality: plot the TRANSPOSE of the return value
	Uses first row as X value and the other rows as y values
	params:
	  - plot: (True/False) Enable/disable plot functionality (default = True)
	  - xlab: string to use as x label
	  - ylab: string to use as y label
	"""
	# try:
	# 	func.doc_plot_opt
	# except:
	# 	if func.__doc__ is None:
	# 		func.__doc__ = ''
	# 	func.__doc__ += '\n' + plot_opt.__doc__
	# 	func.doc_plot_opt = True

	res = func(*args, **kwargs)
	if plot:
		import matplotlib.pyplot as plt
		from matplotlib.ticker import AutoMinorLocator as AML
		fig, ax = plt.subplots()
		plt.plot( res[:,0], res[:,1:])
		plt.xlabel( xlab)
		plt.ylabel( ylab)
		ml1 = AML(5)
		ax.yaxis.set_minor_locator(ml1)
		ax.yaxis.set_tick_params(which='both', right = True)
		plt.legend()
		plt.show()

	return res