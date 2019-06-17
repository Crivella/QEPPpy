import pytest
import numpy as np
from .conftest import BaseTest
from qepppy.calc_system.atoms_list import atoms_list

class Test_atoms(BaseTest):
	cls_typ = atoms_list

	@pytest.fixture(scope='class')
	def cls_wcc(self, cls, coord):
		cls.atoms_coord_cryst = coord

		assert np.allclose(cls.atoms_coord_cryst, coord), "Failed to assign coordinates"

		return cls

	def test_atoms_coord_cryst(self, cls_wcc):
		pass

	def test_atoms_typ(self, cls_wcc):
		n_atoms = cls_wcc.n_atoms
		cls_wcc.atoms_typ = (['Si', 'Ge'] * (n_atoms))[:n_atoms]
		assert cls_wcc.all_atoms_typ == ['Si', 'Ge'][:min(2,n_atoms)]

	def test_atoms_mass(self, cls_wcc):
		n_types = cls_wcc.n_types
		cls_wcc.atoms_mass = np.random.rand(n_types)

	def test_atoms_pseudo(self, cls_wcc):
		app = [a + '.UPF' for a in cls_wcc.all_atoms_typ]
		cls_wcc.atoms_pseudo = app





