import numpy as np

from .FFTgrid import FFTgrid
from ..parsers import fortran_binary_io as bin_io
from .. import utils

class wavefnc(bin_io, FFTgrid):
    binary_format =[
        [
            {'type':'i4', 'shape':(1,), 'name':'kpt_num'},
            {'type':'f8', 'shape':(3,), 'name':'kpt'},
            {'type':'i4', 'shape':(1,), 'name':'ispin'},
            {'type':'i4', 'shape':(1,), 'name':'gamma_only'},
            {'type':'f8', 'shape':(1,), 'name':'scale_factor'},
        ],
        [
            {'type':'i4', 'shape':(1,), 'name':'max_index'},
            {'type':'i4', 'shape':(1,), 'name':'igwx'},
            {'type':'i4', 'shape':(1,), 'name':'nspin'},
            {'type':'i4', 'shape':(1,), 'name':'nbnd'},
        ],
        [
            {'type':'f8', 'shape':(3,3), 'name':'recipr'},
        ],
        [
            {'type':'i4', 'shape':('igwx',3,), 'name':'gvect'},
        ],
        ([
            {'type':'c16', 'shape':('nspin','igwx',), 'name':'C_kn'},
        ], 'nbnd'),
    ]
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # self.rep = 1

    @property
    def direct(self):
        return utils.recipr_base(self.recipr)
    
    def get_parity_bnd(self, bnd:int, thr: float = 1e-4) -> int:
        """Get the parity of the wavefunction 'wfc' for band  number 'bnd' (indexing from 0).
        !!! This works only for k-points where  2k = G  (TRIM points  k + G = -G)  
            where G is any reciprocal lattice vector.

        Args:
            bnd (int): Index of the band starting from 0
            thr (float, optional): Threshold for checking when numbers are too differeng from 1. Defaults to 1e-4.

        Returns:
            int: The eigenvalue of the parity operator (+1 or -1) or 0 if the parity is not well defined
        """
        # Not al wfcs are normalized
        norm = np.sum(self.C_kn[bnd] * np.conj(self.C_kn[bnd]))

        # Get 2*k-point in crystal coordinates
        kc2 = np.dot(self.direct, self.kpt)/(np.pi)
        # print(f'KPT_CRYS: {kc2}')
        kc2r = np.round(kc2, decimals=0)
        if np.any(np.abs(kc2r - kc2) > thr):
            print(f"{bnd = :3d}:  KPT is not a TRIM point   {kc2 = }")
            return 0
        kc2 = (kc2r).astype(int)
        # print(kc2)

        grid = self._generate_g_grid_(bnd)

        # Parity at k-point:
        # <psi| P |psi> = <psi(r) | psi(-r)> = ... planw-wave expansion  = 
        # = sum{G} (C_knG * C_knG'*)
        # with  G'  such that   G + G' + 2k = 0   (only for TRIM points)
        i1 = np.arange(grid.shape[1]) + kc2[0]
        i2 = np.arange(grid.shape[2]) + kc2[1]
        i3 = np.arange(grid.shape[3]) + kc2[2]

        ri1, ri2, ri3 = np.meshgrid(-i1, -i2, -i3, indexing='ij')
        ri1 = ri1 % grid.shape[1]
        ri2 = ri2 % grid.shape[2]
        ri3 = ri3 % grid.shape[3]
        # grid on G'
        rev_grid = grid[(..., ri1, ri2, ri3)]

        res = np.sum(grid * np.conj(rev_grid)) / norm

        if (diff := np.abs(np.linalg.norm(res) - 1)) > thr:
            print(f"{bnd=:4d}:  EGV_NORM != 1   {diff = }")
            return 0
        else:
            if (diff := np.abs(np.abs(np.real(res)) - 1)) > thr:
                print(f"{bnd = :3d}:  EGV_REAL != 1   {diff = }")
                return 0
    
        return int(np.sign(np.real(res)))


    def get_parity_manyfold(self, n_occ: int) -> int:
        """Get the parity of the wavefunction 'wfc' for the manyfold of occupied states 'n_occ'.

        Args:
            n_occ (int): Number of occupied states

        Returns:
            int: Product of the parities of all the occupied states
        """
        res = 1
        for bnd in range(n_occ):
            app = self.get_parity_bnd(bnd)
            res *= app
        return res


