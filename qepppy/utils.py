import numpy as np

def xyz_mesh(shape, base=None, rep=1):
	n1,n2,n3 = shape
	a = np.linspace(0, rep, n1*rep + 1)[:-1] + .5/n1
	b = np.linspace(0, rep, n2*rep + 1)[:-1] + .5/n2
	c = np.linspace(0, rep, n3*rep + 1)[:-1] + .5/n3

	# Specific order to obtain the array with shape (n1,n2,n3) as the data grid
	# The 'b,a,c' order is because for a 3d meshgrid the resulting shape is (1,2,3) --> (2,1,3)
	# The 'y,x,z' order is because of how the 3d meshgrid output behaves:
	#    x,y,z=np.meshgrid(1,2,3) 
	#       will cause the x to change value along axis=1
	#					   y to change value along axis=0
	#					   z to change value along axis=2
	# Since the FFT grid has the axis=0,1,2 corresponding to x,y,z i need to do the proper remapping
	y,x,z = np.meshgrid(b,a,c)

	XYZ = np.array([x,y,z])

	if not base is None:
		XYZ  = np.dot(
			base.T,
			XYZ.reshape(3,-1)
			).reshape(3,*x.shape)

	return XYZ