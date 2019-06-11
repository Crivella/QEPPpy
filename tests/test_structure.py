import pytest
import itertools
import numpy as np
from .test_atoms import cls_wcc, Test_atoms

Test_atoms.__test__ = False

basis={
	'sc':np.diag([1]*3)*5.1,
	'fcc':np.array(list(set(itertools.permutations([0,.5,.5]))))*10.2,
}

@pytest.fixture(
	scope='class'
	)
def cls_typ():
	from qepppy._structure import _structure
	return _structure

@pytest.fixture(
	scope='class',
	params=['sc', 'fcc',]
	)
def cls(cls_typ, request):
	res = cls_typ()
	assert isinstance(res, cls_typ), "Failed to initialize empty instance of " + repr(cls_typ)

	res.direct = basis[request.param]
	return res


class Test_cell(Test_atoms):
	__test__ = True
	
	@pytest.mark.mpl_image_compare
	@pytest.mark.parametrize('rep', [1,2,3])
	def test_plot_cell(self, cls_wcc, rep):
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		cls_wcc._plot(ax, rep, rep, rep, cell=True, graph_lvl=1)

		return fig




