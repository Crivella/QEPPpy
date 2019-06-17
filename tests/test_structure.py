import pytest
from .test_atoms   import Test_atoms
from .test_lattice import Test_lattice
from qepppy.calc_system.structure import structure

Test_atoms.__test__   = False
Test_lattice.__test__ = False



class Test_structure(Test_atoms, Test_lattice):
	__test__ = True
	cls_typ  = structure

	@pytest.fixture(scope='class')
	def cls_wcc(self, cls_rec, coord):
		cls_rec.atoms_coord_cryst = coord
		return cls_rec
	
	@pytest.mark.mpl_image_compare
	@pytest.mark.parametrize('rep', [1,2,3], ids=['rep:' + str(a) for a in [1,2,3]])
	def test_plot_cell(self, cls_wcc, rep):
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		cls_wcc._plot(ax, rep, rep, rep, cell=True, graph_lvl=1)

		return fig

	@pytest.mark.mpl_image_compare
	def test_plot_recipr(self, cls_wcc):
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		cls_wcc._plot(ax, recipr=True)

		return fig




