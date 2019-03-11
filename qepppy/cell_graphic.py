import numpy as np

def generate_repetition_grid(r1,r2,r3, vect_matrix):
	from itertools import product
	res = np.array(list(product(r3,r2,r1)))
	res = np.dot(vect_matrix.T, res.T).T

	return res

def draw_sphere(ax, radius=1, center=[0,0,0], color="b"):
	x0, y0, z0 = center
	u = np.linspace(0, 2 * np.pi, 10)
	v = np.linspace(0, np.pi, 7)
	x = radius * np.outer(np.cos(u), np.sin(v)) + x0
	y = radius * np.outer(np.sin(u), np.sin(v)) + y0
	z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + z0

	ax.plot_surface(x, y, z, color=color)

def draw_cylinder(ax, radius=1, axis=[0,0,1], start=[0,0,0], color="b"):
	x0, y0, z0 = start

	norm   = np.linalg.norm(axis)
	axis  /= norm
	c_teta = axis[2]
	s_teta = np.sqrt(1 - c_teta**2)
	c_phi  = axis[0]/s_teta if s_teta else 1
	s_phi  = axis[1]/s_teta if s_teta else 0

	u = np.linspace(0, 2*np.pi, 8)
	v = np.linspace(0, 1*norm, 2)

	x = radius * np.outer(np.cos(u)*s_phi + np.sin(u)*c_phi*c_teta, np.ones(np.size(v)))   + axis[0] * v + x0
	y = radius * np.outer(-np.cos(u)*c_phi + np.sin(u)*s_phi*c_teta, np.ones(np.size(v)) ) + axis[1] * v + y0
	z = radius * np.outer(-np.sin(u)*s_teta, np.ones(np.size(v)))                          + axis[2] * v + z0

	ax.plot_surface(x, y, z, color=color)

def draw_bond(ax, start, end, color1, color2, graph_lvl=0):
	if graph_lvl == 0:
		v = np.vstack((start, end))
		ax.plot(v[:,0], v[:,1], v[:,2], color="black", linewidth=1.5)
	elif graph_lvl == 1 or graph_lvl == 2:
		mid = (start + end) / 2
		v = np.vstack((start, mid))
		ax.plot(v[:,0], v[:,1], v[:,2], color=color1, linewidth=3.5)
		v = np.vstack((mid, end))
		ax.plot(v[:,0], v[:,1], v[:,2], color=color2, linewidth=3.5)
	elif graph_lvl == 3:
		mid = (start + end) / 2
		if color1 != color2:
			draw_cylinder(ax, radius=0.15, axis=(end-start)/2, start=start,  color=color1)
			draw_cylinder(ax, radius=0.15, axis=(end-start)/2, start=mid, color=color2)
		else:
			draw_cylinder(ax, radius=0.15, axis=end-start, start=start,  color=color1)
	else:
		raise ValueError("arg 'graph_lvl' must be <= 3")

def draw_atom(ax, X,Y,Z, color='k', name='None', radius=1, graph_lvl=0):
	if graph_lvl == 0 or graph_lvl == 1:
		ax.scatter(
			X, Y, Z, 
			s=80*radius,
			marker="o",
			depthshade=False,
			c=color,
			label=name
			)
	elif graph_lvl == 2 or graph_lvl == 3:
		ax.scatter(
			X, Y, Z, 
			s=10,
			marker="o",
			depthshade=False,
			c=color,
			label=name
			)
		for x,y,z in zip(X,Y,Z):
			draw_sphere(ax, radius=radius*0.3, center=[x,y,z], color=color)
	else:
		raise ValueError("arg 'graph_lvl' must be <= 3")



