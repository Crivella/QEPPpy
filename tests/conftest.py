import pytest
import itertools
import numpy as np

basis={
	'sc':np.diag([1]*3)*5.1,
	'fcc':np.array(list(set(itertools.permutations([0,.5,.5]))))*10.2,
}
coords=np.array([
	[0. , 0. , 0. ],
	[.25, .25, .25],
	])


l = list(basis.keys())
@pytest.fixture(
	scope='module',
	params=l,
	ids=['LATT:' + str(a) for a in l]
	)
def base(request):
	return basis[request.param]

l = list(range(1, len(coords)+1))
@pytest.fixture(
	scope='module',
	params=l,
	ids=['#ATM:' + str(a) for a in l]
	)
def coord(request):
	return coords[:request.param]


# clist=[
# 	np.array([
# 			[.25,.25,.25],
# 			[1,1,1],
# 			[1,.25,1]
# 			]),
# 	]

# l = range(len(clist))
# @pytest.fixture(
# 	scope='module',
# 	params=l,
# 	ids=['coord_set:' + str(a+1) for a in l]
# 	)
# def coord2(request):
# 	return clist[request.param]


@pytest.fixture
def tmpfile(tmpdir):
	import tempfile
	with tempfile.NamedTemporaryFile(dir=tmpdir) as f:
		yield f


class BaseTest():
	cls_typ = object

	@pytest.fixture(scope='class')
	def cls(self):
		cls_typ = self.cls_typ
		res = cls_typ()

		assert isinstance(res, cls_typ), "Failed to initialize empty instance of " + repr(cls_typ)

		return res
