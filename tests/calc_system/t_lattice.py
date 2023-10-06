import numpy as np
import pytest

from qepppy.calc_system.lattice import Lattice

from .conftest import BaseTest


class Test_lattice(BaseTest):
    cls_typ = Lattice

    @pytest.fixture(scope='class')
    def cls_rec(self, cls, base):
        cls.direct = base

        return cls

    def test_translate_atoms_cryst(self, cls_rec):
        coord = np.array([
                [.25,.25,.25],
                [1,1,1],
                [1,.25,1]
                ])
        C,G = cls_rec.translate_coord_into_PC(coord)
        res_C = np.array([
            [.25,.25,.25],
            [0,0,0],
            [0,.25,0]
            ])
        assert np.all(res_C == C), 'Failed to bring atoms into cell correctly.'
