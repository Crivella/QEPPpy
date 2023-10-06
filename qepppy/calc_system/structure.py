import re
from functools import reduce

import numpy as np
from attrs import define

from .. import utils
from .._decorators import file_name_handle, set_self
from .atoms_list import AtomsList
from .lattice import Lattice
from .symmetry import Symmetries

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

try:
    import ase
    import ase.visualize
except ImportError:
    ase = None


@define(slots=False)
class Structure(AtomsList, Lattice):
    symmetries: Symmetries = None

    @file_name_handle('w')
    def save_xyz(self, file):
        latt_str = (
            'Lattice="' +
            ' '.join(str(a) for a in self.direct.flatten()) +
            '"'
            )
        prop_str = (
            'Properties=' +
            ':'.join([
                'species:S:1',
                'pos:R:3'
                ])
            )

        data_str = self.atoms_coord_cart

        if not self.atoms_velocities is None:
            prop_str += ':vel:R:3'
            data_str = np.hstack((data_str, self.atoms_velocities))

        if not self.atoms_forces is None:
            prop_str += ':forces:R:3'
            data_str = np.hstack((data_str, self.atoms_forces))

        file.write(str(self.n_atoms) + '\n')
        file.write(latt_str + ' ' + prop_str + '\n')
        ml = max(len(a) for a in self.atoms_typ) + 1

        for i in range(self.n_atoms):
            file.write(f'{{:{ml}s}}'.format(self.atoms_typ[i]))
            np.savetxt(file, data_str[i], fmt='%14.7e', newline=' ')
            file.write('\n')

    @file_name_handle('r')
    def load_xyz(self, file):
        conv={
            'R': float,
            'I': bool,
            'S':'|U30',
            }

        test    = file.readline()
        if test == '':
            return 1
        n_atoms = int(test)
        prop    = file.readline()

        latt  = re.search(r'Lattice="(.*)"', prop).group(1)
        latt  = np.fromstring(latt, sep=' ').reshape(3,3)
        dprop = re.search(r'Properties=(.*)', prop).group(1).split(':')

        typ = {
            'names':dprop[::3],
            'formats':[(conv[a], int(n)) for n,a in zip(dprop[2::3], dprop[1::3])],
            }

        data = np.loadtxt(file, dtype=typ, max_rows=n_atoms)

        self.direct = latt
        self.atoms_coord_cart = data['pos'].reshape(-1,3)
        if 'vel' in data.dtype.names:
            self.atoms_velocities = data['vel'].reshape(-1,3)
        if 'forces' in data.dtype.names:
            self.atoms_forces = data['forces'].reshape(-1,3)

        self.atoms_typ = list(data['species'].reshape(n_atoms,))


    def make_supercell(self,repX,repY,repZ):
        """
        Build a supercell using the repJ repetitions of the primitive cell along
        the J-th direction.

        Params:
         - repX: repetitions along X
         - repY: repetitions along Y
         - repZ: repetitions along Z

        Return:
         - res: np.array of shape (n_atoms*repX*repY*repZ) containing the
                positions of the atoms in the supercell.
         - typ: list with len = (n_atoms*repX*repY*repZ) containing the
                name/type of all the atoms. (1 to 1 correspondence with 'res')
        """
        assert(isinstance(repX,(int,range,list,tuple)))
        assert(isinstance(repY,(int,range,list,tuple)))
        assert(isinstance(repZ,(int,range,list,tuple)))

        rep    = [range(a) if isinstance(a,int) else a for a in [repX,repY,repZ]]
        R_vec  = utils.generate_repetition_grid(rep, vect_matrix=self.direct)

        n_cell = [max(a)-min(a)+1 for a in rep]
        typ    = np.array(list(self.atoms_typ) * (reduce(lambda x,y: x*y, n_cell)))

        res    = np.empty(shape=(0,3))
        coord  = self.atoms_coord_cart
        for vec in R_vec:
            res = np.vstack((res, coord+vec))

        return res, typ

    def nearest_neighbour(self, max_shell=5):
        """
        Find the nearest neighbour up to a max_shell.

        Params:
         - max_shell: last shell of nearest neighour returned

        Return:
         - dist: Array of shape (max_shell,) containing the radius of every shell
                 of nearest neighours.
         - counts: Array of shape (max_shell,) containing the number of atoms for
                   every shell of nearest neighours.
        """
        assert isinstance(max_shell, int)

        m        = max_shell
        l        = range(-m, m+1)
        coord, _ = self.make_supercell(l,l,l)
        norm     = np.linalg.norm(coord, axis=1)

        dist, counts = np.unique(norm, return_counts=True)

        return dist[1:max_shell], counts[1:max_shell]

    def _plot_mpl(
        self, ax,
        repX=1, repY=1, repZ=1,
        cell=False,
        bonds=True,
        recipr=False,
        graph_lvl=1,):

        if recipr:
            ax.set_xlabel(r'$x (Bohr^{-1})$')
            ax.set_ylabel(r'$y (Bohr^{-1})$')
            ax.set_zlabel(r'$z (Bohr^{-1})$')
            self.draw_Wigner_Seitz(ax)
            # plt.show()
            return
        else:
            ax.set_xlabel('x (Bohr)')
            ax.set_ylabel('y (Bohr)')
            ax.set_zlabel('z (Bohr)')

        L, typ = self.make_supercell(repX,repY,repZ)

        self.draw_atoms(ax, L, typ, graph_lvl=graph_lvl)
        if cell:
            self.draw_direct_cell(ax)
        if bonds:
            self.draw_bonds(ax, L, typ, graph_lvl=graph_lvl)
        ax.legend()

    def plot_mpl(
        self,
        *args, **kwargs
        ):
        """
        Plot the crystal cell structure.
        Args:
         - reprX/Y/Z=1/1/1[or any positive integer]:
                repetitions of the cell along X/Y/Z (basis vector not Cartesian!!!).
                NOTE: They are 3 separate arguments repX, repY, repZ
         - cell=False[or True]: plot the contour of the cell.
         - bonds=True[or False]: plot the chemical bonds between atoms.
         - recip=False[or True]: plot the Brilloiun Zone instead of the real cell.
                If True all other flags are ignored.
         - graph_lvl=1[or 0/2/3]:
           - 0: Basic plot with circle dots as atoms and black lines as bonds.
           - 1: Colored line as bonds (color of the nearest atom).
           - 2: Use 3d spheres for atoms and bonds as in 1.
           - 3: Use 3d spheres for atoms and cylinders for bonds.
        """
        if plt is None:
            raise ImportError('matplotlib must be installed and accessible to the PYTHONPATH.')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        self._plot_mpl(ax, *args, **kwargs)

        plt.show()

    def _make_ase_atoms(self):
        if ase is None:
            raise ImportError('ase must be installed and accessible to the PYTHONPATH.')

        new = ase.Atoms(
            symbols          = self.atoms_typ,
            scaled_positions = self.atoms_coord_cryst,
            cell             = self.direct * ase.units.Bohr,
            masses           = self.atoms_mass
            )

        return new

    def plot_ase(self, *args, **kwargs): # pylint: disable=unused-argument
        if ase is None:
            raise ImportError('ase must be installed and accessible to the PYTHONPATH.')

        ase.visualize.view(self._make_ase_atoms())

    def plot(self, *args, mode='mpl', **kwargs):
        if mode == 'mpl':
            self.plot_mpl(*args, **kwargs)
        elif mode == 'ase':
            self.plot_ase(*args, **kwargs)
        else:
            raise NotImplementedError(f'Plot mode {mode} not implemented.')


    @set_self('symmetries')
    def get_symmetries(self):
        # import ase.spacegroup
        # from .symmetry import Symmetries, Symmetry

        # sg = ase.spacegroup.get_spacegroup(self._make_ase_atoms())

        new = Symmetries()

        # for rot in sg.rotations:
        #     new.append(symmetry(rotation=rot))

        return new

    # @set_self('symmetries')
    # def get_symmetries_spglib(self):
    #     import spglib
    #     from .symmetry import symmetries, symmetry

    #     nums = [self.unique_atoms_typ.index(a) for a in self.atoms_typ]
    #     cell = (
    #         self.direct,
    #         self.atoms_coord_cryst,
    #         nums
    #         )
    #     dataset = spglib.get_symmetry(cell)

    #     new = symmetries()
    #     for rot in dataset['rotations']:
    #         new.append(symmetry(rotation=rot))

    #     return new
