import numpy as np
from decorator import decorator

@decorator
def save_opt(
	func, cls, *args, 
	pFile=True,
	fname="", fmt="", 
	**kwargs
	):
	res = func(cls, *args, **kwargs)
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
	func, cls, *args, 
	plot=True, 
	xlab="", ylab="", 
	**kwargs
	):
	res = func(cls, *args, **kwargs)

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