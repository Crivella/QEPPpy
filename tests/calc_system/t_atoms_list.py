import numpy as np
import pytest

from qepppy.calc_system.atoms_list import atoms_list

from .conftest import BaseTest


class Test_atoms(BaseTest):
    cls_typ = atoms_list

    @pytest.fixture(scope='class')
    def cls_wcc(self, cls, coord):
        cls.atoms_coord_cryst = coord

        assert np.allclose(cls.atoms_coord_cryst, coord), "Failed to assign coordinates"

        return cls

    def test_atoms_coord_cryst(self, cls_wcc):
        pass

    def test_set_atoms_typ(self, cls_wcc):
        n_atoms = cls_wcc.n_atoms
        cls_wcc.atoms_typ = (['Si', 'Ge'] * (n_atoms))[:n_atoms]
        assert cls_wcc.unique_atoms_typ == ['Si', 'Ge'][:min(2,n_atoms)]

    def test_set_atoms_mass(self, cls_wcc):
        n_atoms = cls_wcc.n_atoms
        n_types = cls_wcc.n_types
        masses  = np.arange(n_atoms)
        cls_wcc.atoms_mass = masses

        assert np.allclose(cls_wcc.atoms_mass, masses) , "Failed to assing atom_mass"
        assert cls_wcc.unique_atoms_mass == list(masses)[:n_types], "Failed to assign all_atom_mass from atom_mass"

    def test_set_unique_atoms_mass(self, cls_wcc):
        n_atoms = cls_wcc.n_atoms
        n_types = cls_wcc.n_types
        masses  = np.arange(n_types)
        cls_wcc.unique_atoms_mass = masses

        assert np.allclose(cls_wcc.unique_atoms_mass, masses), "Failed to assign unique_atoms_mass"
        assert np.allclose(
            (list(masses)*n_atoms)[:n_atoms],
            cls_wcc.atoms_mass
            ), "Failed to assign atoms_mass from unique_atoms_mass"

    def test_set_atoms_pseudo(self, cls_wcc):
        app  = [a + '.UPF' for a in cls_wcc.atoms_typ]
        app2 = [a + '.UPF' for a in cls_wcc.unique_atoms_typ]
        cls_wcc.atoms_pseudo = app

        assert list(cls_wcc.atoms_pseudo) == app, "Failed to assign atoms_pseudo"
        assert list(cls_wcc.unique_atoms_pseudo) == app2, "Failed to assign unique_atoms_pseudo form atoms_pseudo"

    def test_set_unique_atoms_pseudo(self, cls_wcc):
        app  = [a + '.UPF' for a in cls_wcc.atoms_typ]
        app2 = [a + '.UPF' for a in cls_wcc.unique_atoms_typ]
        cls_wcc.unique_atoms_pseudo = app2

        assert list(cls_wcc.atoms_pseudo) == app, "Failed to assign atoms_pseudo from unique_atoms_pseudo"
        assert list(cls_wcc.unique_atoms_pseudo) == app2, "Failed to assign unique_atoms_pseudo"