import pytest
import numpy as np

@pytest.fixture
def cls():
	from qepppy._kpoints import _kpoints
	return _kpoints

@pytest.fixture
def cls_rec(cls):
	res = cls()
	res.recipr = np.diag([1]*3)

	return res

@pytest.fixture(
	params=[1,5,10]
	)
def point_list_3d(request):
	return np.random.randint(0,1,(request.param,3))

def test_kpt_cart(cls_rec, point_list_3d):
	cls_rec.kpt_cart = point_list_3d
	assert np.all(cls_rec.kpt_cart == cls_rec.kpt_cryst)

def test_kpt_cryst(cls_rec, point_list_3d):
	cls_rec.kpt_cryst = point_list_3d
	assert np.all(cls_rec.kpt_cart == cls_rec.kpt_cryst)

@pytest.mark.parametrize('shift',[(0,0,0),(1,1,1)])
@pytest.mark.parametrize('shape',[(1,1,1),(5,5,5),(10,10,7)])
def test_kpt_mesh(cls_rec, shape, shift):
	cls_rec.generate_monkhorst_pack_grid(shape, shift)
