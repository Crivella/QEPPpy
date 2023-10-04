import numpy as np
import pytest

from qepppy.calc_system.structure import structure

from .t_atoms_list import Test_atoms
from .t_lattice import Test_lattice


class Test_structure(Test_atoms, Test_lattice):
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

        cls_wcc._plot_mpl(ax, rep, rep, rep, cell=True, graph_lvl=1)

        return fig

    @pytest.mark.mpl_image_compare
    def test_plot_recipr(self, cls_wcc):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        cls_wcc._plot_mpl(ax, recipr=True)

        return fig

    def test_save_xyz(self, new_cls, cls_wcc, tmpfile):
        cls_wcc.atoms_typ    = ['Si']*cls_wcc.n_atoms
        cls_wcc.atoms_forces = np.random.rand(cls_wcc.n_atoms,3)
        cls_wcc.save_xyz(tmpfile.file)

        tmpfile.file.seek(0)

        new_cls.load_xyz(tmpfile.file)

        assert np.allclose(cls_wcc.direct, new_cls.direct)
        assert np.allclose(cls_wcc.atoms_coord_cryst, new_cls.atoms_coord_cryst)
        assert np.allclose(cls_wcc.atoms_forces, new_cls.atoms_forces)
        assert cls_wcc.atoms_typ == new_cls.atoms_typ