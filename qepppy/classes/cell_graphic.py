import numpy as np

import logging
logger = logging.getLogger( __name__)
logging.basicConfig( format='%(levelname)s: %(name)s\n%(message)s\n')

periodic_table={
	'H':{'color':'','radius':1.001511716,},
	'He':{'color':'','radius':0.585789872,},
	'Li':{'color':'','radius':3.155706727,},
	'Be':{'color':'','radius':2.116402116,},
	'B':{'color':'','radius':1.64399093,},
	'C':{'color':'','radius':1.26606198,},
	'N':{'color':'','radius':1.058201058,},
	'O':{'color':'','radius':0.907029478,},
	'F':{'color':'','radius':0.793650794,},
	'Ne':{'color':'','radius':0.718065004,},
	'Na':{'color':'','radius':3.590325019,},
	'Mg':{'color':'','radius':2.739984883,},
	'Al':{'color':'','radius':2.229780801,},
	'Si':{'color':'blue','radius':2.210884354,},
	'P':{'color':'','radius':1.851851852,},
	'S':{'color':'','radius':1.64399093,},
	'Cl':{'color':'','radius':1.49281935,},
	'Ar':{'color':'','radius':1.34164777,},
	'K':{'color':'','radius':4.591836735,},
	'Ca':{'color':'','radius':3.665910809,},
	'Sc':{'color':'','radius':3.476946334,},
	'Ti':{'color':'','radius':3.325774754,},
	'V':{'color':'','radius':3.231292517,},
	'Cr':{'color':'','radius':3.13681028,},
	'Mn':{'color':'','radius':3.042328042,},
	'Fe':{'color':'','radius':2.947845805,},
	'Co':{'color':'','radius':2.872260015,},
	'Ni':{'color':'','radius':2.815570673,},
	'Cu':{'color':'','radius':2.739984883,},
	'Zn':{'color':'','radius':2.68329554,},
	'Ga':{'color':'yellow','radius':2.569916856,},
	'Ge':{'color':'','radius':2.362055933,},
	'As':{'color':'','radius':2.154195011,},
	'Se':{'color':'','radius':1.946334089,},
	'Br':{'color':'','radius':1.776266062,},
	'Kr':{'color':'','radius':1.64399093,},
	'Rb':{'color':'','radius':5.007558579,},
	'Sr':{'color':'','radius':4.138321995,},
	'Y':{'color':'','radius':4.006046863,},
	'Zr':{'color':'','radius':3.892668178,},
	'Nb':{'color':'','radius':3.741496599,},
	'Mo':{'color':'','radius':3.590325019,},
	'Tc':{'color':'','radius':3.458049887,},
	'Ru':{'color':'','radius':3.363567649,},
	'Rh':{'color':'','radius':3.269085412,},
	'Pd':{'color':'','radius':3.193499622,},
	'Ag':{'color':'','radius':3.117913832,},
	'Cd':{'color':'','radius':3.042328042,},
	'In':{'color':'','radius':2.947845805,},
	'Sn':{'color':'','radius':2.739984883,},
	'Sb':{'color':'','radius':2.513227513,},
	'Te':{'color':'','radius':2.324263039,},
	'I':{'color':'','radius':2.173091459,},
	'Xe':{'color':'','radius':2.040816327,},
	'Cs':{'color':'','radius':5.631141345,},
	'Ba':{'color':'','radius':4.780801209,},
	'La':{'color':'','radius':0,},
	'Ce':{'color':'','radius':0,},
	'Pr':{'color':'','radius':4.667422525,},
	'Nd':{'color':'','radius':3.892668178,},
	'Pm':{'color':'','radius':3.873771731,},
	'Sm':{'color':'','radius':4.497354497,},
	'Eu':{'color':'','radius':4.365079365,},
	'Gd':{'color':'','radius':4.40287226,},
	'Tb':{'color':'','radius':4.25170068,},
	'Dy':{'color':'','radius':4.308390023,},
	'Ho':{'color':'','radius':4.270597128,},
	'Er':{'color':'','radius':4.270597128,},
	'Tm':{'color':'','radius':4.195011338,},
	'Yb':{'color':'','radius':4.195011338,},
	'Lu':{'color':'','radius':4.100529101,},
	'Hf':{'color':'','radius':3.930461073,},
	'Ta':{'color':'','radius':3.779289494,},
	'W':{'color':'','radius':3.647014361,},
	'Re':{'color':'','radius':3.552532124,},
	'Os':{'color':'','radius':3.495842782,},
	'Ir':{'color':'','radius':3.401360544,},
	'Pt':{'color':'','radius':3.344671202,},
	'Au':{'color':'','radius':3.287981859,},
	'Hg':{'color':'','radius':3.231292517,},
	'Tl':{'color':'','radius':2.947845805,},
	'Pb':{'color':'','radius':2.91005291,},
	'Bi':{'color':'','radius':2.702191988,},
	'Po':{'color':'','radius':2.551020408,},
	'At':{'color':'','radius':2.399848828,},
	'Rn':{'color':'','radius':2.267573696,},
	'Fr':{'color':'','radius':0,},
	'Ra':{'color':'','radius':0,},
	'Ac':{'color':'','radius':0,},
	'Th':{'color':'','radius':0,},
	'Pa':{'color':'','radius':0,},
	'U':{'color':'','radius':0},
	'Np':{'color':'','radius':0,},
	'Pu':{'color':'','radius':0,},
	'Am':{'color':'','radius':0,},
	'Cm':{'color':'','radius':0,},
	'Bk':{'color':'','radius':0,},
	'Cf':{'color':'','radius':0,},
	'Es':{'color':'','radius':0,},
	'Fm':{'color':'','radius':0,},
	'Md':{'color':'','radius':0,},
	'No':{'color':'','radius':0,},
	'Lr':{'color':'','radius':0,},
	'Rf':{'color':'','radius':0,},
	'Db':{'color':'','radius':0,},
	'Sg':{'color':'','radius':0,},
	'Bh':{'color':'','radius':0,},
	'Hs':{'color':'','radius':0,},
	'Mt':{'color':'','radius':0,},
	}

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
	if graph_lvl == 0 or graph_lvl == 1:
		ax.scatter( 
			X, Y, Z, 
			s=50,
			marker="o",
			depthshade=False,
			c=periodic_table[name]['color']
			)
	elif graph_lvl == 2 or graph_lvl == 3:
		radius = periodic_table[name]['radius']
		color = periodic_table[name]['color']
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
				ax.plot( v[:,0], v[:,1], v[:,2], color="black", linewidth=0.5)
			elif graph_lvl == 1 or graph_lvl == 2:
				mid = (v1 + v2) / 2
				v = np.vstack( (v1, mid))
				ax.plot( v[:,0], v[:,1], v[:,2], color=periodic_table[name1]['color'], linewidth=0.5)
				v = np.vstack( (mid, v2))
				ax.plot( v[:,0], v[:,1], v[:,2], color=periodic_table[name2]['color'], linewidth=0.5)
			elif graph_lvl == 3:
				mid = (v1+v2)/2
				draw_cylinder( ax, radius=0.15, axis=(v2-v1)/2, start=v1,  color=periodic_table[name1]['color'])
				draw_cylinder( ax, radius=0.15, axis=(v2-v1)/2, start=mid, color=periodic_table[name2]['color'])
			else:
				logger.error( "arg 'graph_lvl' must be <= 2")

	return