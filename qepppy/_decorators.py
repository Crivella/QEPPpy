import numpy as np
from decorator import decorator


@decorator
def save_opt(
	func, 
	_fname='', _fmt="",
	*args, 
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
	if not fmt:
		fmt= _fmt

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

@decorator
def plot_opt(
	func,
	_xlab="", _ylab="",
	*args, 
	plot=True,
	start=1, end=-1,
	xlab="", ylab="",
	colors=['k','r','b','g','c','m'],
	labels=[''],
	dash_list=[],
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

	if not xlab:
		xlab = _xlab
	if not ylab:
		ylab = _ylab

	res = func(*args, **kwargs)
	if not plot:
		return res

	import matplotlib.pyplot as plt
	from matplotlib.ticker import AutoMinorLocator as AML
	fig, ax = plt.subplots()
	cl = len(colors)
	X = res[:,0]

	if res.shape[1] == 2:
		y_data = res[:,1].reshape(1,X.size)
	else:
		y_data = res[:,start:end].T

	for i,Y in enumerate(y_data):
		if dash_list:
			dash = dash_list[i%len(dash_list)]
		else:
			dash = (8,2*((i//cl)%2))
		ax.plot( 
			X, Y, 
			color=colors[i%cl], 
			label=labels[i%len(labels)],
			dashes=dash,
			)
	ax.set_xlabel(xlab)
	ax.set_ylabel(ylab)
	ml1 = AML(5)
	ax.yaxis.set_minor_locator(ml1)
	ax.yaxis.set_tick_params(which='both', right = True)
	fig.legend()
	plt.show()

	return res