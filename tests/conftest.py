import pytest
import numpy as np

@pytest.fixture
def tmpfile(tmpdir):
	import tempfile
	with tempfile.NamedTemporaryFile(dir=tmpdir) as f:
		yield f

@pytest.fixture(
	scope='class'
	)
def cls_typ():
	return object

# @pytest.fixture(
# 	scope='class'
# 	)
# def cls(cls_typ):
# 	res = cls_typ()
# 	assert isinstance(res, cls_typ), "Failed to initialize empty instance of " + repr(cls_typ)
# 	return res

# @pytest.fixture(
# 	scope='class',
# 	params=[1,5,15,40]
# 	)
# def cls_wcc(cls, request):
# 	cls.atoms_coord_cart = np.random.rand(request.param, 3)
# 	return cls



