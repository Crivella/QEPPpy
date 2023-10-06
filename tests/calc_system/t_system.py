import numpy as np

from qepppy.calc_system import System

from .t_structure import Test_structure


class Test_system(Test_structure):
    cls_typ = System

    def test_save_steps(self, new_cls, cls_wcc, tmpfile):
        from qepppy.calc_system.structure import Structure
        cls_wcc.atoms_typ    = ['Si']*cls_wcc.n_atoms
        cls_wcc.atoms_forces = np.random.rand(cls_wcc.n_atoms,3)

        steps = []
        for i in np.linspace(0,1,11):
            new = Structure()
            new.direct           = cls_wcc.direct
            new.atoms_coord_cart = cls_wcc.atoms_coord_cart + i
            new.atoms_typ        = cls_wcc.atoms_typ
            new.atoms_forces     = cls_wcc.atoms_forces
            # new.atoms_velocities = cls_wcc.atoms_velocities
            steps.append(new)

        cls_wcc.steps = steps
        cls_wcc.save_step_xyz(tmpfile.file)

        tmpfile.file.seek(0)

        new_cls.load_step_xyz(tmpfile.file)

        for sn, so in zip(new_cls.steps, cls_wcc.steps):
            assert np.allclose(sn.direct, so.direct)
            assert np.allclose(sn.atoms_coord_cart, so.atoms_coord_cart)
            assert np.allclose(sn.atoms_forces, so.atoms_forces)
            # assert np.allclose(sn.atoms_velocities == so.atoms_velocities)
            assert np.all(sn.atoms_typ == so.atoms_typ)
