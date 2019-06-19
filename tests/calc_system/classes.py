import pytest
import numpy as np
from qepppy.calc_system.atoms_list import atoms_list
from qepppy.calc_system.lattice    import lattice
from qepppy.calc_system.structure  import structure
from qepppy.calc_system.kpoints    import kpoints
from qepppy.calc_system            import system

class BaseTest():
	cls_typ = object

	@pytest.fixture
	def new_cls(self):
		cls_typ = self.cls_typ
		res = cls_typ()

		assert isinstance(res, cls_typ), "Failed to initialize empty instance of " + repr(cls_typ)

		return res

	@pytest.fixture(scope='class')
	def cls(self):
		cls_typ = self.cls_typ
		res = cls_typ()

		assert isinstance(res, cls_typ), "Failed to initialize empty instance of " + repr(cls_typ)

		return res

class Test_atoms(BaseTest):
	cls_typ = atoms_list

	@pytest.fixture(scope='class')
	def cls_wcc(self, cls, coord):
		cls.atoms_coord_cryst = coord

		assert np.allclose(cls.atoms_coord_cryst, coord), "Failed to assign coordinates"

		return cls

	def test_atoms_coord_cryst(self, cls_wcc):
		pass

	def test_atoms_typ(self, cls_wcc):
		n_atoms = cls_wcc.n_atoms
		cls_wcc.atoms_typ = (['Si', 'Ge'] * (n_atoms))[:n_atoms]
		assert cls_wcc.all_atoms_typ == ['Si', 'Ge'][:min(2,n_atoms)]

	def test_atoms_mass(self, cls_wcc):
		n_types = cls_wcc.n_types
		cls_wcc.atoms_mass = np.random.rand(n_types)

	def test_atoms_pseudo(self, cls_wcc):
		app = [a + '.UPF' for a in cls_wcc.all_atoms_typ]
		cls_wcc.atoms_pseudo = app


class Test_lattice(BaseTest):
	cls_typ = lattice

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
		assert np.all(res_C == C), "Failed to bring atoms nito cell correctly."


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


class Test_kpoints(Test_lattice):
	cls_typ = kpoints

	@pytest.fixture(
		params=[1,5,10]
		)
	def point_list_3d(self, request):
		return np.random.randint(0,1,(request.param,3))

	def test_kpt_cart(self, cls_rec, point_list_3d):
		cls_rec.kpt_cart = point_list_3d
		assert np.all(cls_rec.kpt_cart == cls_rec.kpt_cryst)

	def test_kpt_cryst(self, cls_rec, point_list_3d):
		cls_rec.kpt_cryst = point_list_3d
		assert np.all(cls_rec.kpt_cart == cls_rec.kpt_cryst)

	@pytest.mark.parametrize('shift',[(0,0,0),(1,1,1)])
	@pytest.mark.parametrize('shape',[(1,1,1),(5,5,5),(10,10,7)])
	def test_kpt_mesh(self, cls_rec, shape, shift):
		cls_rec._generate_monkhorst_pack_grid(shape, shift)



class Test_system(Test_structure):
	cls_typ = system

	def test_save_steps(self, new_cls, cls_wcc, tmpfile):
		from qepppy.calc_system.structure import structure
		cls_wcc.atoms_typ    = ['Si']*cls_wcc.n_atoms
		cls_wcc.atoms_forces = np.random.rand(cls_wcc.n_atoms,3)

		steps = []
		for i in np.linspace(0,1,11):
			new = structure()
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
			assert sn.atoms_typ == so.atoms_typ
