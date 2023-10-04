from .._decorators import file_name_handle, set_self
from .bands import bands
from .structure import structure


class system(structure, bands):
    steps={
        'typ':(list,),
        'sub_typ':(structure,),
        'doc':"""List of configurations for the atoms during time/relaxation steps."""
        }

    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)

    #     if self.direct != [] and self.atoms_coord_cryst != []:
    #         self.get_symmetries()

    @file_name_handle('w')
    def save_step_xyz(self, file):
        for step in self.steps:
            step.save_xyz(file)

    @file_name_handle('r')
    def load_step_xyz(self, file):
        steps = []
        while True:
            new = structure()
            if new.load_xyz(file):
                break
            steps.append(new)

        self.steps = steps

    def _make_ase_steps(self):
        return [a._make_ase_atoms() for a in self.steps]

    def plot_step_ase(self):
        from ase.visualize import view

        view(self._make_ase_steps())

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
