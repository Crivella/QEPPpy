import numpy as np
from .meta import PropertyCreator
from . import cell_graphic as cg


import json
from pkg_resources import resource_string
periodic_table = json.loads(resource_string('qepppy.data', 'periodic_table.json').decode('utf-8'))

def get_all_atoms_typ(cls, value):
	res = []
	for name in value:
		if not name in res:
			res.append(name)

	return res

def split_atom_list_by_name(atom_coord, atom_names):
	from scipy.spatial import KDTree
	trees  = []
	rad    = []
	names  = []

	atom_names = np.array(atom_names)
	for n in set(atom_names):
		coord = atom_coord[atom_names == n,:]

		names.append(n)
		trees.append(KDTree(coord))
		rad.append(periodic_table[n]['radius'])
	return trees, np.array(names), rad

class _atoms(metaclass=PropertyCreator):
	atoms_coord_cart={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,3),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		# 'post_set_name':'_atoms_coord_cryst',
		# 'post_set_func':_cart_to_cryst_,
		'doc':"""List of atomic coordinate in CARTESIAN basis."""
		}

	atoms_coord_cryst={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape': (-1,3),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		# 'post_set_name':'_atoms_coord_cart',
		# 'post_set_func':_cryst_to_cart_,
		'doc':"""List of atomic coordinate in CRYSTAL basis."""
		}

	atoms_typ={
		'typ':(list,),
		'sub_typ':(str,np.ndarray,),
		'shape':('n_atoms',),
		'post_set_name':'_all_atoms_typ',
		'post_set_func':get_all_atoms_typ,
		'doc':"""List of atom names (same order as the list of coordinates)."""
		}

	atoms_mass={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'shape':('n_types',),
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'doc':"""List of atomic masses."""
		}

	atoms_pseudo={
		'typ':(list,np.ndarray,),
		'shape':('n_types',),
		'doc':"""List of atomic pseudopotential files."""
		}

	all_atoms_typ={
		'typ':(list,np.ndarray,),
		'sub_typ':(str,),
		'doc':"""List of atom names (same order as list of masses)."""
		}

	@property
	def n_atoms(self):
		return len(self.atoms_coord_cart)

	@property
	def n_types(self):
		return len(self.all_atoms_typ)

	@property
	def atoms_group_coord_cart(self):
		"""List of atomic coordinate in cartesian basis grouped by atom name
		in a dictionary."""
		return {a:np.array(self.atoms_coord_cart[np.array(self.atoms_typ) == a]) for a in self.all_atoms_typ}

	@property
	def atoms_group_coord_cryst(self):
		"""List of atomic coordinate in crystal basis grouped by atom name
		in a dictionary."""
		return {a:np.array(self.atoms_coord_cryst[np.array(self.atoms_typ) == a]) for a in self.all_atoms_typ}

	def draw_atoms(self, ax, atom_coord=None, atom_names=None, **kwargs):
		"""
		Draw atoms onto a matplotlib axis object.
		Params:
		 - ax:         Matplotlib axis object
		 - atom_coord: List of atomic coordinates (not necessarily the original 
		               one in order to use a supercell).
		 - atom_namse: List of atom names (same shape as atom_coord). Used to
		               plot the proper color and atomic radius.
		"""
		if atom_coord is None:
			atom_coord = self.atoms_coord_cart
		if atom_names is None:
			atom_names = self.atoms_typ

		trees, names, rad = split_atom_list_by_name(atom_coord, atom_names)

		for tree,n,r in zip(trees,names,rad):
			X,Y,Z = tree.data.T
			color = periodic_table[n].get('color', 'k')
			cg.draw_atom(ax, X,Y,Z, color=color, name=n, radius=r, **kwargs)


	def draw_bonds(self, ax, atom_coord=None, atom_names=None, **kwargs):
		"""
		Draw atomic bonds onto a matplotlib axis object.
		Params:
		 - ax:         Matplotlib axis object
		 - atom_coord: List of atomic coordinates (not necessarily the original 
		               one in order to use a supercell).
		 - atom_namse: List of atom names (same shape as atom_coord). Used to
		               plot the proper color and atomic radius.
		"""
		from itertools import combinations_with_replacement as cwr

		if atom_coord is None:
			atom_coord = self.atoms_coord_cart
		if atom_names is None:
			atom_names = self.atoms_typ

		trees, names, rad = split_atom_list_by_name(atom_coord, atom_names)

		for rad_t, (tree1,tree2), (name1,name2) in zip(cwr(rad,2), cwr(trees,2), cwr(names,2)):
			bonds = tree1.query_ball_tree(tree2, sum(rad_t))
			c1 = periodic_table[name1]['color']
			c2 = periodic_table[name2]['color']
			for i1,b in enumerate(bonds):
				for i2 in b:
					if i1 == i2 and tree1 == tree2:
						continue
					start = tree1.data[i1]
					end   = tree2.data[i2]
					cg.draw_bond(ax, start, end, c1, c2, **kwargs)





