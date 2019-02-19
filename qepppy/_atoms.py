import numpy as np
from ._decorators import store_property
from . import cell_graphic as cg

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

class _atoms():
	@property
	@store_property
	def atoms_group_coord_cart(self):
		return {a:np.array(self.atoms_coord_cart[np.array(self.atoms_typ) == a]) for a in self.all_atoms_typ}

	@property
	@store_property
	def atoms_group_coord_cryst(self):
		return {a:np.array(self.atoms_coord_cryst[np.array(self.atoms_typ) == a]) for a in self.all_atoms_typ}

	@property
	@store_property
	def atoms_typ(self):
		return self._atoms_typ

	@property
	@store_property
	def atoms_mass(self):
		return self._atoms_mass

	@property
	@store_property
	def atoms_pseudo(self):
		return self._atoms_pseudo

	@property
	@store_property
	def all_atoms_typ(self):
		return self._all_atoms_typ


	def draw_atoms(self, ax, atom_coord, atom_names, **kwargs):
		trees, names, rad = split_atom_list_by_name(atom_coord, atom_names)

		for tree,n,r in zip(trees,names,rad):
			X,Y,Z = tree.data.T
			color = periodic_table[n].get('color', 'k')
			cg.draw_atom(ax, X,Y,Z, color=color, name=n, radius=r, **kwargs)


	def draw_bonds(self, ax, atom_coord, atom_names, **kwargs):
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