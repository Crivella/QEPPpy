import os

import numpy as np
import pytest

from qepppy.calc_system.kpoints import Kpoints

from .t_lattice import Test_lattice

cwd = os.path.dirname(os.path.realpath(__file__))

edges = np.array([
    [0.500, 0.500, 0.500,   30],
    [0.000, 0.000, 0.000,   30],
    [0.500, 0.000, 0.500,   15],
    [0.625, 0.250, 0.625,   1],
    [0.375, 0.375, 0.750,   30],
    [0.000, 0.000, 0.000,   1],
    ])

class Test_kpoints(Test_lattice):
    cls_typ = Kpoints

    @pytest.fixture(
        params=[1,5,10]
        )
    def point_list_3d(self, request):
        return np.random.randint(0,1,(request.param,3))

    @pytest.fixture(scope='class')
    def cls_edges(self, cls):
        cls.direct = np.array([
            [-.5,0,.5],
            [0,.5,.5],
            [-.5,.5,0],
            ]) * 10.21
        cls.generate_kpath(edges=edges)

        return cls


    def test_kpt_cart(self, cls_rec, point_list_3d):
        cls_rec.kpt_cart = point_list_3d
        assert np.all(cls_rec.kpt_cart == cls_rec.kpt_cryst)

    def test_kpt_cryst(self, cls_rec, point_list_3d):
        cls_rec.kpt_cryst = point_list_3d
        assert np.all(cls_rec.kpt_cart == cls_rec.kpt_cryst)

    def test_kpt_path(self, cls_rec):
        cls_rec.generate_kpath(edges=edges, mode='crystal')

    def test_kpt_unfold_path(self, cls_edges):
        file = os.path.join(cwd, '../qe/test_files', 'unf.pwscf')
        unf = cls_edges.generate_unfolding_path(cls_edges.recipr/3)
        test = np.loadtxt(file)[:,:3]

        assert np.allclose(unf, test), 'Incorrect unfold path.'

    @pytest.mark.parametrize('shift',[(0,0,0),(1,1,1)])
    @pytest.mark.parametrize('shape',[(1,1,1),(5,5,5),(10,10,7)])
    def test_kpt_mesh(self, cls_rec, shape, shift):
        cls_rec.generate_monkhorst_pack_grid(shape, shift, do_set_self=False)
