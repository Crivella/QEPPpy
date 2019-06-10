import pytest
import numpy as np
from qepppy._atoms import _atoms

cls_typ = _atoms
plist   = [1,5,15,40]

# def random_elements(num):
# 	import random
# 	import string

# 	l = string.ascii_lowercase + ' '
# 	L = string.ascii_lowercase

# 	res = [random.choice(L) + random.choice(l) for _ in range(num)]

# 	return res


# class Test_atoms(class_tester):
# 	typ            = _atoms
# 	run_empty_init = True
# 	test_recipr={}
# 	test_atoms_coord_cart={
# 		'value':np.random.rand(15,3)*0.7
# 		}
# 	test_atoms_coord_cryst={
# 		'value':np.random.rand(15,3)*0.7
# 		}
# 	test_atoms_typ={
# 		'value':random_elements(15),
# 		'prev_call':'test_atoms_coord_cryst',
# 	}

class Test_atoms():
	cls_typ = _atoms

	@pytest.fixture(scope='module')
	def cls(self):
		res = self.cls_typ()
		assert isinstance(res, self.cls_typ), "Failed to initialize empty instance of " + repr(self.cls_typ)
		return res

	@staticmethod
	@pytest.mark.parametrize('nat', [1,5,15,40])
	def test_atoms_coord_cart(cls, nat):
		cls.atoms_coord_cart = np.random.rand(nat,3)

		return cls

	@staticmethod
	def test_atoms_typ(cls):
		n_atoms = cls.n_atoms
		cls.atoms_typ = (['Si', 'Ge'] * (n_atoms))[:n_atoms]
		assert cls.all_atoms_typ == ['Si', 'Ge'][:min(2,n_atoms)]

	@staticmethod
	def test_atoms_mass(cls):
		n_types = cls.n_types
		cls.atoms_mass = np.random.rand(n_types)

	@staticmethod
	def test_atoms_pseudo(cls):
		app = [a + '.UPF' for a in cls.all_atoms_typ]
		cls.atoms_pseudo = app






