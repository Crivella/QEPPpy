import numpy as np
from ..calc_system import system
from .qe_structure import qe_structure as structure
from .parser.input_files import input_files as inp_f

class pw_in(inp_f, structure, system):
	"""
	Instance used to handle QE pw.x input files.
	Provides an interface for the 'structure' output file methods, in order to 
	call methods related to the atomic structure.
	kwargs:
	 - input_file = Name of the file to parse
	 - input_data = Dictionary used to initialize namelists of the input.

	Parse parameters can be accessed using the following syntax:
	 - obj.param_name
	 - obj['param_name']

	Namelist parameters can be set using the syntax:
	 - obj['NAMELIST/param_name']    = value
	 - obj['NAMELIST/param_name(n)'] = value   (For vector values)
	   param_name are the same as in the QE documentation (online or *.def)

	Card parameters can be set using the syntax:
	 - obj['CARD/param_name']    = value
	   param_name are the same as in the QE documentation (online or *.def)

	"""
	templ_file = "INPUT_PW.templ"

	def _set_atoms_coord(self, value):
		value = np.array(value)
		if value.shape[1] != 3:
			raise ValueError("Coordinates must be a 3 element list/array.")
		self['ATOMIC_POSITIONS/x'] = value[:,0]
		self['ATOMIC_POSITIONS/y'] = value[:,1]
		self['ATOMIC_POSITIONS/z'] = value[:,2]

	def _check_atoms_len(self, value):
		if len(value) != self._n_atoms:
			raise ValueError("Number of assigned atoms does not match n_atoms ({} != {}).".format(
				len(value), self._n_atoms))

	@property
	def _atoms(self):
		name, x, y, z = self._find("X", "x", "y", "z", up="ATOMIC_POSITIONS")
		coord = list(zip(x, y, z))
		return [{'name':n, 'coord':c} for n,c in zip(name, coord)]

	@structure._atoms_coord_cart.setter
	def _atoms_coord_cart(self, value):
		value = np.array(value)

		self._check_atoms_len(value)
		self._set_atoms_coord(value)

		self._app_atom_p = 'alat'

	@structure._atoms_coord_cryst.setter
	def _atoms_coord_cryst(self, value):
		value = np.array(value)

		self._check_atoms_len(value)
		self._set_atoms_coord(value)
		
		self._app_atom_p = 'crystal'

	@structure._atoms_typ.setter
	def _atoms_typ(self, value):
		self._check_atoms_len(value)
		self['ATOMIC_POSITIONS/X'] = value

	def _check_atoms_spec_len(self, value):
		if len(value) != self.n_types:
			raise ValueError("Number of assigned atoms does not match n_types ({} != {}).".format(
				len(value), self.n_types))

	def _check_atoms_spec_type(self, value, t, msg="..."):
		if any(not isinstance(a, t) for a in value):
			raise TypeError("{} types must all be of type {}.".format(msg, t))

	@property
	def _atom_spec(self):
		# name, mass, pfile = self._find("X", "Mass_X", "PseudoPot_X", up="ATOMIC_SPECIES")
		names = self._find("X", up="ATOMIC_SPECIES")
		try:
			mass = self._find("Mass_X", up="ATOMIC_SPECIES")
		except KeyError:
			mass = [0] * len(names)
		try:
			pfile = self._find("PseudoPot_X", up="ATOMIC_SPECIES")
		except KeyError:
			pfile =[n + '.UPF' for n in names]
		return [{'name':n, 'mass':m, 'pseudo_file':p} for n,m,p in zip(names, mass, pfile)]

	@structure._unique_atoms_typ.setter
	def _unique_atoms_typ(self, value):
		self._check_atoms_spec_len(value)
		self._check_atoms_spec_type(value, str, 'Atoms')

		self['ATOMIC_SPECIES/X'] = value

	@structure._unique_atoms_mass.setter
	def _unique_atoms_mass(self, value):
		self._check_atoms_spec_len(value)
		self._check_atoms_spec_type(value, (int, float), 'Mass')

		self['ATOMIC_SPECIES/Mass_X'] = value

	@structure._unique_atoms_pseudo.setter
	def _unique_atoms_pseudo(self, value):
		self._check_atoms_spec_len(value)
		self._check_atoms_spec_type(value, str, 'Pseudo')

		self['ATOMIC_SPECIES/PseudoPot_X'] = value

	@property
	def _cell(self):
		a1, a2, a3 = self._find("v1", "v2", "v3")
		return [{'a1':a1, 'a2':a2, 'a3':a3}]

	@structure._direct.setter
	def _direct(self, value):
		# if self.ibrav != 0:
		# 	raise ValueError("CELL can be set only if ibrav == 0.")
		a1, a2, a3 = np.array(value)
		self['CELL_PARAMETERS/v1'] = a1
		self['CELL_PARAMETERS/v2'] = a2
		self['CELL_PARAMETERS/v3'] = a3

	@property
	def celldm(self):
		return self._find("celldm")

	@property
	def alat(self):
		return self._find("celldm(1)")

	@alat.setter
	def alat(self, value):
		if not isinstance(value, (int,float)):
			raise TypeError("Value for alat must be of type {} not {}".format(
				float, type(value)))
		self['SYSTEM/celldm(1)'] = float(value)

	@property
	def ibrav(self):
		return self._find("ibrav")

	@ibrav.setter
	def ibrav(self, value):
		self['SYSTEM/ibrav'] = value

	@property
	def _app_atom_p(self):
		res = self._find("ATOMIC_POSITIONS")
		if res is None:
			res = 'alat'
		return res

	@_app_atom_p.setter
	def _app_atom_p(self, value):
		possib = self.namelist_c._templ_['ATOMIC_POSITIONS']['c']
		if possib and not value in possib:
			raise ValueError("ATOMIC_POSITIONS card value must be one of {} not '{}'.".format(
				possib, value))
		self["ATOMIC_POSITIONS"].value = value

	@property
	def _app_cell_p(self):
		res = self._find("CELL_PARAMETERS")
		if res is None:
			res = 'alat'
		return res

	@_app_cell_p.setter
	def _app_cell_p(self, value):
		possib = self.namelist_c._templ_['CELL_PARAMETERS']['c']
		if possib and not value in possib:
			raise ValueError("CELL_PARAMETERS card value must be one of {} not '{}'.".format(
				possib, value))
		self['CELL_PARAMETERS/'] = value

	@property
	def _n_atoms(self):
		res = self._find("nat")
		if res is None:
			return 0
		return res

	@_n_atoms.setter
	def _n_atoms(self, value):
		if not isinstance(value, int):
			raise TypeError("Number of atoms must be integer.")

		self['SYSTEM/nat'] = value

	@property
	def _n_types(self):
		res = self._find("ntyp")
		if res is None:
			return 0
		return res

	@_n_types.setter
	def _n_types(self, value):
		if not isinstance(value, int):
			raise TypeError("Number of types must be integer.")

		self['SYSTEM/ntyp'] = value


	def __getattr__(self, attr):
		try:
			return super().__getitem__(attr)
		except:
			return object.__getattribute__(self, attr)

	def __str__(self):
		msg = inp_f.__str__(self)
		# if not 'ATOMIC_POSITIONS' in self.card_c:
		# 	msg += structure.__str__(self)
		return msg

	def validate(self):
		super().validate()



