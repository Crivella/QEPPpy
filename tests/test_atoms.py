import pytest
import numpy as np

plist   = [1,5,15,40]

@pytest.fixture(
	scope='class'
	)
def cls_typ():
	from qepppy.atoms_list import atoms_list
	return atoms_list

@pytest.fixture(
	scope='class'
	)
def cls(cls_typ):
	res = cls_typ()
	assert isinstance(res, cls_typ), "Failed to initialize empty instance of " + repr(cls_typ)
	return res

@pytest.fixture(
	scope='class',
	params=[1,2]
	)
def cls_wcc(cls, request):
	cls.atoms_coord_cryst = np.array([[0,0,0],[.25,.25,.25]])[:request.param]
	return cls

class Test_atoms():
	def test_atoms_coord_cryst(self, cls_wcc):
		pass

	def test_atoms_typ(self, cls_wcc):
		n_atoms = cls_wcc.n_atoms
		# print(n_atoms)
		cls_wcc.atoms_typ = (['Si', 'Ge'] * (n_atoms))[:n_atoms]
		assert cls_wcc.all_atoms_typ == ['Si', 'Ge'][:min(2,n_atoms)]

	def test_atoms_mass(self, cls_wcc):
		n_types = cls_wcc.n_types
		cls_wcc.atoms_mass = np.random.rand(n_types)

	def test_atoms_pseudo(self, cls_wcc):
		app = [a + '.UPF' for a in cls_wcc.all_atoms_typ]
		cls_wcc.atoms_pseudo = app





