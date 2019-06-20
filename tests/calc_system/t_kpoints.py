import pytest
import numpy as np
from .t_lattice import Test_lattice
from qepppy.calc_system.kpoints import kpoints

class Test_kpoints(Test_lattice):
	cls_typ = kpoints

	@pytest.fixture(
		params=[1,5,10]
		)
	def point_list_3d(self, request):
		return np.random.randint(0,1,(request.param,3))

	def test_kpt_cart(self, cls_rec, point_list_3d):
		cls_rec.kpt_cart = point_list_3d
		assert np.all(cls_rec.kpt_cart == cls_rec.kpt_cryst)

	def test_kpt_cryst(self, cls_rec, point_list_3d):
		cls_rec.kpt_cryst = point_list_3d
		assert np.all(cls_rec.kpt_cart == cls_rec.kpt_cryst)

	@pytest.mark.parametrize('shift',[(0,0,0),(1,1,1)])
	@pytest.mark.parametrize('shape',[(1,1,1),(5,5,5),(10,10,7)])
	def test_kpt_mesh(self, cls_rec, shape, shift):
		cls_rec.generate_monkhorst_pack_grid(shape, shift, set_self=False)