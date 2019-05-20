import numpy as np
# from ._decorators import store_property
from .meta import PropertyCreator
from . import cell_graphic as cg
from .utils import _cart_to_cryst_, _cryst_to_cart_

import json
from pkg_resources import resource_string
periodic_table = json.loads(resource_string('qepppy.qe.parser.data', 'periodic_table.json').decode('utf-8'))

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
	n_atoms={
		'typ':(int,np.int),
		'default':0,
		'doc':"""Number of atoms."""
		}
		
	n_types={
		'typ':(int,np.int),
		'default':0,
		'doc':"""Number of atomic types."""
		}

	atoms_coord_cart={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'size':'n_atoms * 3',
		'usize':True,
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'set_other_name':'_atoms_coord_cryst',
		'set_other_func':_cart_to_cryst_,
		'doc':"""List of atomic coordinate in cartesian basis."""
		}

	atoms_coord_cryst={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'size':'n_atoms * 3',
		'usize':True,
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'set_other_name':'_atoms_coord_cart',
		'set_other_func':_cryst_to_cart_,
		'doc':"""List of atomic coordinate in crystal basis."""
		}

	atoms_typ={
		'typ':(list,),
		'sub_typ':(str,np.ndarray,),
		'size':'n_atoms',
		'doc':"""List of atom names (same order as the list of coordinates)."""
		}

	atoms_mass={
		'typ':(list,np.ndarray),
		'sub_typ':(int,float,np.number),
		'size':'n_types',
		'conv_func':lambda x: np.array(x, dtype=np.float),
		'doc':"""List of atomic masses."""
		}

	atoms_pseudo={
		'typ':(list,np.ndarray,),
		'size':'n_types',
		'doc':"""List of atomic pseudopotential files."""
		}

	all_atoms_typ={
		'typ':(list,np.ndarray,),
		'sub_typ':(str,),
		'size':'n_types',
		'usize':True,
		'doc':"""List of atom names (same order as list of masses)."""
		}


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

	def draw_atoms(self, ax, atom_coord, atom_names, **kwargs):
		"""
		Draw atoms onto a matplotlib axis object.
		Params:
		 - ax:         Matplotlib axis object
		 - atom_coord: List of atomic coordinates (not necessarily the original 
		               one in order to use a supercell).
		 - atom_namse: List of atom names (same shape as atom_coord). Used to
		               plot the proper color and atomic radius.
		"""
		trees, names, rad = split_atom_list_by_name(atom_coord, atom_names)

		for tree,n,r in zip(trees,names,rad):
			X,Y,Z = tree.data.T
			color = periodic_table[n].get('color', 'k')
			cg.draw_atom(ax, X,Y,Z, color=color, name=n, radius=r, **kwargs)


	def draw_bonds(self, ax, atom_coord, atom_names, **kwargs):
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

		trees, names, rad = split_atom_list_by_name( atom_coord, atom_names)

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





