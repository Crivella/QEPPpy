import numpy as np
from .logger import logger, error

import json
from pkg_resources import resource_string
periodic_table = json.loads(resource_string('qepppy.qe.parser.data', 'periodic_table.json').decode('utf-8'))


def generate_repetition_grid(reps, vect_matrix):
	from itertools import product
	res = np.array(list(product(reps,reps,reps)))
	res = np.dot(vect_matrix.T, res.T).T

	return res


@logger()
def draw_sphere(ax, radius=1, center=[0,0,0], color="b"):
	x0, y0, z0 = center
	u = np.linspace(0, 2 * np.pi, 10)
	v = np.linspace(0, np.pi, 7)
	x = radius * np.outer(np.cos(u), np.sin(v)) + x0
	y = radius * np.outer(np.sin(u), np.sin(v)) + y0
	z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + z0

	ax.plot_surface(x, y, z, color=color)

@logger()
def draw_cylinder(ax, radius=1, axis=[0,0,1], start=[0,0,0], color="b"):
	x0, y0, z0 = start

	norm = np.linalg.norm(axis)
	axis = axis / norm
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

@logger()
def cell_repetitions(base, vect, num):
	"""
	Replicate a list of atoms along a vector for num times:
	base: list of atom coords []
	"""
	L0 = base.copy()
	for n in range(1, num):
		base = np.vstack((base, L0 + vect*n))
	return base

def split_atom_list_by_name(atom_coord, atom_names):
	from scipy.spatial import KDTree
	coords = []
	trees  = []
	rad    = []
	names  = []

	atom_names = np.array(atom_names)
	for n in set(atom_names):
		coord = atom_coord[np.where(atom_names == n)[0],:]

		coords.append(coord)
		names.append(n)
		trees.append(KDTree(coord))
		rad.append(periodic_table[n]['radius'])
	return np.array(coords), np.array(names), rad, trees


def draw_atoms(ax, atom_coord, atom_names, graph_lvl=0):
	coords, names, rad, _ = split_atom_list_by_name(atom_coord, atom_names)

	for c,n,r in zip(coords,names,rad):
		X,Y,Z = c.T
		color = periodic_table[n]['color']
		if graph_lvl == 0 or graph_lvl == 1:
			ax.scatter(
				X, Y, Z, 
				s=80*r,
				marker="o",
				depthshade=False,
				c=color,
				label=n
				)
		elif graph_lvl == 2 or graph_lvl == 3:
			ax.scatter(
				X, Y, Z, 
				s=10,
				marker="o",
				depthshade=False,
				c=color,
				label=n
				)
			for x,y,z in zip(X,Y,Z):
				draw_sphere(ax, radius=r*0.3, center=[x,y,z], color=color)
		else:
			raise ValueError("arg 'graph_lvl' must be <= 3")


@logger()
def draw_cell(ax, v1=[1,0,0], v2=[0,1,0], v3=[0,0,1], center=[0,0,0]):
	V = [v1,v2,v3]
	for n1 in range(3):
		orig = np.array(center)
		v0 = V[n1]
		for n2 in range(4):
			# print(orig, v0)
			v = np.vstack((orig, orig + v0))
			if n2 == n1:
				orig = V[(n2+1)%3] + V[(n2+2)%3]
			else:
				orig = V[n2%3]
			ax.plot(v[:,0], v[:,1], v[:,2], color="black", linewidth=0.5)
		ax.plot(v[:,0], v[:,1], v[:,2], color="black", linewidth=0.5)

@logger()
def draw_Wigner_Seitz(ax, recip):
	try:
		from scipy.spatial import Voronoi
	except:
		raise ImportError("Scipy module must be installed to print Wigner-Seitz cell.")
	L = generate_repetition_grid([-1,0,1], recip)

	vor = Voronoi(L)
	P = vor.vertices
	R = vor.ridge_vertices

	rad     = max(np.linalg.norm(recip, axis=1)) * np.sqrt(2)/2
	cond    = np.where(np.linalg.norm(P, axis=1) > rad)[0]
	P[cond] = np.zeros(3)

	for i1, e in enumerate(R):
		for i2, r in enumerate(e):
			if r in cond:
				R[i1][i2] = -1

	X = P[:,0]
	Y = P[:,1]
	Z = P[:,2]
	ax.scatter(X,Y,Z, color='green')

	for vert in R:
		vert.append(vert[0])
		v = np.asarray(vert)
		if np.all(v >= 0):
			ax.plot(P[v, 0], P[v, 1], P[v, 2], color='k')


def _draw_bond_(ax, start, end, color1, color2, graph_lvl=0):
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


def draw_bonds(ax, atom_coord, atom_names, **kwargs):
	from itertools import combinations_with_replacement as cwr

	coords, names, rad, trees = split_atom_list_by_name(atom_coord, atom_names)

	for tree1, tree2 in cwr(trees,2):
		n1 = trees.index(tree1)
		n2 = trees.index(tree2)
		
		bonds = tree1.query_ball_tree(tree2, rad[n1] + rad[n2])
		c1 = periodic_table[names[n1]]['color']
		c2 = periodic_table[names[n2]]['color']
		for i1,b in enumerate(bonds):
			for i2 in b:
				if i1 == i2 and tree1 == tree2:
					continue
				start = coords[n1][i1]
				end   = coords[n2][i2]
				_draw_bond_(ax, start, end, c1, c2, **kwargs)

