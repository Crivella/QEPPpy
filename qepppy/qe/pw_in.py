import numpy as np

from ..calc_system import system
from ..parsers import fortran_namelist_collection as fnc
from .parser.input_files import qe_input
from .qe_structure import qe_structure as structure


class pw_in(qe_input, structure, system):
    cards = ['CELL_PARAMETERS', 'ATOMIC_SPECIES','ATOMIC_POSITIONS','K_POINTS',
            'CONSTRAINTS','OCCUPATIONS','ATOMIC_FORCES']

    _link__n_atoms={'item':'SYSTEM/nat'}
    _link__n_types={'item':'SYSTEM/ntyp'}
    _link__n_bnd={'item':'SYSTEM/nbnd'}
    _link_ibrav={'item':'SYSTEM/ibrav'}
    _link_alat={'item':'SYSTEM/celldm(1)'}
    _link_test={'item':'SYSTEM/nat, SYSTEM/ntyp'}

    def __init__(self, *args, input_file=None, **kwargs):
        input_file = kwargs.pop('src', input_file)

        fnc.__init__(self, *args, src=input_file, **kwargs)
        kwargs.pop('input_data', None)
        structure.__init__(self, *args, **kwargs)

    def __str__(self):
        res  = fnc.__str__(self)
        res += structure.__str__(self)

        res += '\n\nK_POINTS'
        if len(self.kpt_mesh) > 0:
            res += ' {automatic}\n'
            res += '   {0[0]:d} {0[1]:d} {0[2]:d}    {1[0]:d} {1[1]:d} {1[2]:d}'.format(
                self.kpt_mesh, self.kpt_shift)
        elif len(self.kpt_edges) > 0:
            if 'cryst' in self.kpt_mode:
                res += ' {crystal_b}\n'
            else:
                res += ' {tpiba_b}\n'
            res += f'  {len(self.kpt_edges)}\n'
            for e in self.kpt_edges:
                res += '    {0[0]:11.8f} {0[1]:11.8f} {0[2]:11.8f}    {0[3]:.0f}\n'.format(e)
        else:
            res += f' {{crystal}}\n  {self.n_kpt:d}\n'
            for k,w in zip(self.kpt_cryst, self.kpt_weight):
                res += '  {0[0]:11.8f} {0[1]:11.8f} {0[2]:11.8f}    {1:9.6f}\n'.format(k,w)


        return res

    def __getattr__(self, key):
        try:
            return super().__getattribute__(key)
        except AttributeError:
            res = super().__getitem__(key)
            if res is None:
                raise
            return res

    def read_atomic_species(self, content):
        typs   = []
        mass   = []
        pseudo = []

        for n,l in enumerate(content.strip().split('\n')[1:]):
            if n >= self.n_types:
                break
            t,m,p = list(filter(None, l.split(' ')))
            typs.append(t)
            mass.append(float(m))
            pseudo.append(p)

        self.unique_atoms_typ    = typs
        self.unique_atoms_mass   = mass
        self.unique_atoms_pseudo = pseudo

    def read_atomic_positions(self, content):
        typ    = []
        x_l    = []
        y_l    = []
        z_l    = []
        # if_x   = []
        # if_y   = []
        # if_z   = []

        mode = self.get_mode(content)
        if not mode is None:
            mode = mode.lower()

        if self.ibrav != 0:
            self.direct = self._ibrav_to_cell_()

        ls_old = None
        for n,l in enumerate(content.strip().split('\n')[1:]):
            if n >= self.n_atoms:
                break
            s = list(filter(None, l.split(' ')))

            ls = len(s)
            if   ls == 4 and (ls == ls_old or ls_old is None):
                t,x,y,z = s
            elif ls == 7 and (ls == ls_old or ls_old is None):
                t,x,y,z,ix,iy,iz = s
                raise NotImplementedError()
            else:
                raise ValidateError('Invalid format for ATOMIC_POSITIONS card.')

            typ.append(t)
            x_l.append(float(x))
            y_l.append(float(y))
            z_l.append(float(z))

            ls_old = ls

        res = np.array([x_l,y_l,z_l], dtype=float).T.reshape(-1,3)

        self.atoms_typ = typ
        if mode is None or mode == 'alat':
            res *= self.alat
            self.atoms_coord_cart = res
        elif mode == 'bohr':
            self.atoms_coord_cart = res
        elif mode == 'angstrom':
            from ..units import BOHR_to_A
            res /= BOHR_to_A
            self.atoms_coord_cart = res
        elif mode == 'crystal':
            self.atoms_coord_cryst = res
        elif mode == 'crystal_sg':
            raise NotImplementedError()
        else:
            raise NotImplementedError()

    def read_k_points(self, content):
        nks    = None
        xk     = []
        yk     = []
        zk     = []
        wk     = []

        mode = self.get_mode(content)
        if not mode is None:
            mode = mode.lower()

        if mode == 'automatic':
            l = content.strip().split('\n')[1]
            s = list(filter(None, l.split(' ')))
            mx,my,mz,sx,sy,sz = s
            mesh  = tuple(int(a) for a in [mx,my,mz])
            shift = tuple(int(a) for a in [sx,sy,sz])

            self.kpt_mesh = mesh
            self.kpt_shift = shift
            return
        elif mode == 'gamma':
            raise NotImplementedError()
            return

        content = content.strip().split('\n')[1:]
        nks     = int(content[0].strip())

        for n,l in enumerate(content[1:]):
            if n >= nks:
                break
            s = list(filter(None, l.split(' ')))

            x,y,z,w = s

            xk.append(float(x))
            yk.append(float(y))
            zk.append(float(z))
            wk.append(float(w))

        self.n_kpt = nks
        res_c = np.array([xk,yk,zk], dtype=float).T
        res_w = np.array(wk, dtype=float)
        if mode is None or mode == 'tpiba':
            self.kpt_cart   = res_c
            self.kpt_weight = res_w
        elif mode == 'crystal':
            self.kpt_cryst  = res_c
            self.kpt_weight = res_w
        elif mode == 'tpiba_b':
            self.kpt_mode = 'cart'
            res = np.hstack((res_c,res_w.reshape(-1,1)))
            self.kpt_edges = res
        elif mode == 'crystal_b':
            self.kpt_mode = 'cryst'
            res = np.hstack((res_c,res_w.reshape(-1,1)))
            self.kpt_edges = res
        elif mode == 'tpiba_c':
            raise NotImplementedError()
        elif mode == 'crystal_c':
            raise NotImplementedError()
        else:
            raise NotImplementedError()

    def read_cell_parameters(self, content):
        res = [[]]*3

        if self.ibrav != 0:
            return
        mode = self.get_mode(content)
        if not mode is None:
            mode = mode.lower()

        for n,l in enumerate(content.strip().split('\n')[1:]):
            if n >= 3:
                break
            s = list(filter(None, l.split(' ')))
            res[n] = s

        res = np.array(res, dtype=float)
        if mode is None or mode == 'alat':
            res *= self.alat
        elif mode == 'bohr':
            # self.alat = res[0,0]
            pass
        elif mode == 'angstrom':
            from ..units import BOHR_to_A
            res /= BOHR_to_A
            # self.alat = res[0,0]

        self.direct = res

    def read_constraints(self, content):
        pass

    def read_occupations(self, content):
        pass

    def read_atomic_forces(self, content):
        pass
