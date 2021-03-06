import numpy as np

class TypedList(list):
	def __init__(self, typ):
		self.typ = typ

	def append(self, value):
		if not isinstance(value, self.typ):
			raise TypeError(f"List '{repr(self)}', only accept type '{self.typ}'.")
		super().append(value)

def recipr_base(base):
	return np.linalg.inv(base).T * 2 * np.pi

def generate_repetition_grid(r1,r2,r3, vect_matrix=None):
	from itertools import product
	res = np.array(list(product(r1,r2,r3)))
	if not vect_matrix is None:
		res = res.dot(vect_matrix)

	return res

def xyz_mesh(shape, base=None, rep=1, reverse=False):
	nspin, n1,n2,n3 = shape
	try:
		r1,r2,r3 = rep
	except:
		r1 = r2 = r3 = rep
	a = np.linspace(0, r1, n1*r1 + 1)[1:] #[:-1] #+ .5/n1
	b = np.linspace(0, r2, n2*r2 + 1)[1:] #[:-1] #+ .5/n2
	c = np.linspace(0, r3, n3*r3 + 1)[1:] #[:-1] #+ .5/n3

	if reverse:
		a,b,c = c,b,a

	# Specific order to obtain the array with shape (n1,n2,n3) as the data grid
	# The 'b,a,c' order is because for a 3d meshgrid the resulting shape is (1,2,3) --> (2,1,3)
	# The 'y,x,z' order is because of how the 3d meshgrid output behaves:
	#    x,y,z=np.meshgrid(1,2,3) 
	#       will cause the x to change value along axis=1
	#					   y to change value along axis=0
	#					   z to change value along axis=2
	# Since the FFT grid has the axis=0,1,2 corresponding to x,y,z i need to do the proper remapping
	y,x,z = np.meshgrid(b,a,c)

	if reverse:
		x,y,z = z,y,x

	if not base is None:
		XYZ  = np.dot(
			base.T,
			[x.flatten(),y.flatten(),z.flatten()],
			)
	else:
		XYZ = [x,y,z]

	return np.array(XYZ).reshape(3,*x.shape)

def lowdin_ortho(base):
	"""
	O = (M . M^T)^{-1/2} . M
	Params:
	 - base: np.array with states along axis=0. If the states are 
	         multidimensional (eg: FFT grids), flatten them before and after the
	         orthonormalization, restore the shape
	"""
	from scipy.linalg import sqrtm, inv

	shape   = base[0].shape
	base    = np.array([a.flatten() for a in base])
	overlap = np.dot(base.conj(), base.T)
	print("OVERLAP: \n", overlap)

	base = np.dot(
		inv(sqrtm(overlap)).conj(),
		base
		)
	
	overlap = np.dot(base.conj(), base.T)
	print("OVERLAP: \n", overlap)

	return base.reshape(-1,*shape)

def remap_plane(
	invT, 
	Xlim, Ylim, Zlim, 
	shape, rep, fixaxis=2
	):
	lim = (Xlim, Ylim, Zlim)
	shape = np.array(shape)

	xmin, xmax = Xlim
	ymin, ymax = Ylim
	zmin, zmax = Zlim

	m1 = np.linspace(*lim[(fixaxis+1)%3], shape[0] * rep[0])
	m2 = np.linspace(*lim[(fixaxis+2)%3], shape[1] * rep[1])
	a,b = np.meshgrid(m1,m2)

	fix = np.ones(a.shape) * (sum(lim[fixaxis]))/2
	rec = np.round(
			np.dot(
			invT,
			[a.flatten(),b.flatten(),fix.flatten()]
			),
			decimals=8
		) % 1

	rec = np.array(rec) * shape.reshape(3,1)
	rec = rec.astype(dtype='int')
	i,j,k = rec
	
	return a, b, tuple((i,j,k))


