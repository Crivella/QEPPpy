import numpy as np

import logging
logger = logging.getLogger( __name__)
logging.basicConfig( format='%(levelname)s: %(name)s\n%(message)s\n')

import json
from pkg_resources import resource_string
periodic_table = json.loads( resource_string( 'qepppy.data', 'periodic_table.json').decode('utf-8'))


def draw_sphere( ax, radius=1, center=[0,0,0], color="b"):
	x0, y0, z0 = center
	u = np.linspace(0, 2 * np.pi, 10)
	v = np.linspace(0, np.pi, 7)
	x = radius * np.outer(np.cos(u), np.sin(v)) + x0
	y = radius * np.outer(np.sin(u), np.sin(v)) + y0
	z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + z0

	ax.plot_surface( x, y, z, color=color)
	return

def draw_cylinder( ax, radius=1, axis=[0,0,1], start=[0,0,0], color="b"):
	x0, y0, z0 = start

	norm = np.linalg.norm( axis)
	axis = axis / norm
	c_teta = axis[2]
	s_teta = np.sqrt( 1 - c_teta**2)
	c_phi  = axis[0]/s_teta if s_teta else 1
	s_phi  = axis[1]/s_teta if s_teta else 0

	u = np.linspace(0, 2 * np.pi, 8)
	v = np.linspace(0, 1*norm, 2)

	x = radius * np.outer(np.cos(u)*s_phi + np.sin(u)*c_phi*c_teta, np.ones(np.size(v)))   + axis[0] * v + x0
	y = radius * np.outer(-np.cos(u)*c_phi + np.sin(u)*s_phi*c_teta, np.ones(np.size(v)) ) + axis[1] * v + y0
	z = radius * np.outer(-np.sin(u)*s_teta, np.ones(np.size(v)))                          + axis[2] * v + z0

	#print( "ct={},\nst={},\ncp={},\nsp={}".format( c_teta, s_teta, c_phi, s_phi))

	ax.plot_surface( x, y, z, color=color)
	return

def cell_repetitions( base, vect, num):
	L0 = base.copy()
	for n in range( 1, num):
		base = np.vstack( ( base, L0 + vect*n))
	return base

def draw_atoms( ax, atom_list, name, graph_lvl=0):
	LP = np.array( [a[1] for a in filter( lambda x: x[0] == name, atom_list)])
	if len( LP) == 0: return
	X = LP[:,0]
	Y = LP[:,1]
	Z = LP[:,2]
	radius = periodic_table[name]['radius']
	color = periodic_table[name]['color']
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
		for x,y,z in zip( X,Y,Z):
			draw_sphere( ax, radius=radius*0.3, center=[x,y,z], color=color)
	else:
		logger.error( "arg 'graph_lvl' must be <= 2")
	return

def draw_cell( ax, v1, v2, v3):
	V = [v1,v2,v3]
	for n1 in range( 3):
		orig = np.array([0,0,0])
		v0 = V[n1]
		for n2 in range( 4):
			#print( orig, v0)
			v = np.vstack( (orig, orig + v0))
			if n2 == n1:
				orig = V[(n2+1)%3] + V[(n2+2)%3]
			else:
				orig = V[n2%3]
			ax.plot( v[:,0], v[:,1], v[:,2], color="black", linewidth=0.5)
		ax.plot( v[:,0], v[:,1], v[:,2], color="black", linewidth=0.5)
	return



def draw_bonds( ax, atom_list, graph_lvl=0):
	for n1, a1 in enumerate( atom_list):
		name1 = a1[0]
		v1 = a1[1]
		r1 = periodic_table[name1]['radius']
		for n2, a2 in enumerate( atom_list[n1+1:]):
			name2 = a2[0]
			v2 = a2[1]
			r2 = periodic_table[name2]['radius']
			delta = np.linalg.norm( v1 - v2)
			max_l = r1+r2
			#print( name1, a1, name2, a2, delta, max_l)
			if delta > max_l: continue
			if graph_lvl == 0:
				v = np.vstack( (v1, v2))
				ax.plot( v[:,0], v[:,1], v[:,2], color="black", linewidth=1.5)
			elif graph_lvl == 1 or graph_lvl == 2:
				mid = (v1 + v2) / 2
				v = np.vstack( (v1, mid))
				ax.plot( v[:,0], v[:,1], v[:,2], color=periodic_table[name1]['color'], linewidth=3.5)
				v = np.vstack( (mid, v2))
				ax.plot( v[:,0], v[:,1], v[:,2], color=periodic_table[name2]['color'], linewidth=3.5)
			elif graph_lvl == 3:
				mid = (v1+v2)/2
				if name1 != name2:
					draw_cylinder( ax, radius=0.15, axis=(v2-v1)/2, start=v1,  color=periodic_table[name1]['color'])
					draw_cylinder( ax, radius=0.15, axis=(v2-v1)/2, start=mid, color=periodic_table[name2]['color'])
				else:
					draw_cylinder( ax, radius=0.15, axis=v2-v1, start=v1,  color=periodic_table[name1]['color'])
			else:
				logger.error( "arg 'graph_lvl' must be <= 2")

	return