import numpy as np

def _gaussian( x, broad):
	dx = x[1] - x[0]
	sigma = broad
	twoSigmaSq = 2.0 * sigma**2
	norm = np.sqrt(2 * np.pi) * sigma / dx 
	#Divided for dx to take into account implementation as sum instead of integral (volume element)

	g = np.arange( -10*sigma, 10*sigma, dx)
	g = np.exp( -(g**2) / twoSigmaSq)
	g /= norm
	return g

def _lorentz( x, broad):
	dx = x[1] - x[0]
	gamma = broad
	gammaSq = gamma**2
	norm = np.pi / gamma / dx 
	#Divided for dx to take into account implementation as sum instead of integral (volume element)

	g = np.arange( -15*gamma, 15*gamma, dx)
	g = 1/(g**2 + gammaSq)
	g /= norm
	return g

sw = {
	'gauss':_gaussian,
	'lorentz':_lorentz,
}
def broad( data, t="gauss", deg=0.1, axis=0):
	if not t in sw:
		raise Exception( "Invalide type of broadening '{}'.`n".format( t))

	if axis:
		data = np.transpose( data)
	res = np.zeros( data.shape)

	x = res[:,0] = data[:,0]
	y = data[:,1:]

	conv = sw[t]( x, deg)
	lgh = int( len( conv)/2)
	for n,yi in enumerate( y.transpose()):
		new = np.convolve( yi, conv, mode='full')
		d = len( x) - len( new) + lgh*2
		new = new[ lgh:-lgh+d]

		res[:,n+1] = np.array( new)

	if axis:
		res = np.transpose( res)

	return res


if __name__ == "__main__":
	import sys
	argc = len( sys.argv)
	if not 2<=argc<=5 or sys.argv[1] == 'help':
		print("Incorrect use. Pleas pass arguments:"
			"\n\t'fname'\t (comma separated),"
			"\n\t'type\t(gauss/lorentz, default=gauss)',"
			"\n\t'broad\t (default=0.1)'"
			"\n\t'oname\t(optional) (comma separated)'")
		exit()

	t = 'gauss'
	if argc >= 3:
		t = str( sys.argv[2])
	deg = 0.1
	if argc >= 4:
		deg = float( sys.argv[3])
	oname = None
	if argc >= 5:
		oname = str( sys.argv[4]).split( ',')
	for n, name in enumerate(str( sys.argv[1]).split( ',')):
		print("Applying '{}' broadening of '{}' to '{}'".format( t.upper(), deg, name))
		data = np.loadtxt( name, usecols=None)

		b = broad( data, t, deg)

		if oname:
			n_name = oname[n]
		else:
			n_name = "{}_{}".format( name, broad)
		print( "Saving broadened data to '{}'".format( n_name))
		np.savetxt( n_name, b)

		
