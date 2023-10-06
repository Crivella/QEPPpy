import numpy as np
from attrs import Factory, define

from .._decorators import file_name_handle, set_self
from .bands import Bands
from .kpoints import cart_to_cryst, cryst_to_cart
from .structure import Structure

try:
    import ase
    import ase.visualize
except ImportError:
    ase = None

@define(slots=False)
class System(Structure, Bands):
    steps: list = Factory(list)
    irrep_mapping: np.ndarray = None

    @file_name_handle('w')
    def save_step_xyz(self, file):
        for step in self.steps:
            step.save_xyz(file)

    @file_name_handle('r')
    def load_step_xyz(self, file):
        steps = []
        while True:
            new = Structure()
            if new.load_xyz(file):
                break
            steps.append(new)

        self.steps = steps

    def _make_ase_steps(self):
        return [a._make_ase_atoms() for a in self.steps] # pylint: disable=protected-access

    def plot_step_ase(self):
        ase.visualize.view(self._make_ase_steps())

    @set_self('atoms_coord_cryst')
    def translate_into_PC(self):
        """
        Translate all the atoms into the Primitive Cell.
        """
        res, _ = self.translate_coord_into_PC(self.atoms_coord_cryst)
        return res

    @set_self('kpt_cryst')
    def translate_into_FBZ(self):
        """
        Translate all the k-points into the First BZ.
        """
        res, _ = self.translate_coord_into_FBZ(self.kpt_cryst)
        return res

    @set_self('kpt_cryst,kpt_weight')
    def generate_monkhorst_pack_grid(
            self,
            mesh: tuple[int] = None, shift: tuple[int] = None,
            use_symmeties: bool = True, symm_thr: float = 1E-5
        ) -> tuple[np.ndarray, np.ndarray]:
        """
        Generate a Monkhorst-Pack grid of k-point.
        Params:
         -mesh:  tuple of 3 ints > 0
         -shift: tuple of 3 ints that can be either 0 or 1
         -use_symmeties: if True, use symmetries to reduce the number of k-points
         -symm_thr: threshold for symmetries
        """
        kpts, weight = super().generate_monkhorst_pack_grid(mesh, shift)

        if use_symmeties:
            # Must check symmetries on point in cartesian coordinates
            # Checking on crystal only works if [M,B] = 0
            kpts = cryst_to_cart(self, kpts)
            if getattr(self, 'symmetries', None) is None:
                self.get_symmetries()
            ind, kpts, _ = self.symmetries.reduce(kpts, thr=symm_thr)
            # ind, kpts = self.symmetries.reduce2(kpts, thr=symm_thr)
            self.irrep_mapping  = ind
            kpts = cart_to_cryst(self, kpts)

            unique, counts = np.unique(ind, return_counts=True)
            weight = counts[np.argsort(unique)]

            weight = self._normalize_weight(weight, 1)

        return kpts, weight
