from ..parsers import FortranBinaryIO as FBIO
from ..utils import recipr_base
from .FFTgrid import FFTgrid


class charge_density(FBIO, FFTgrid):
    binary_format =[
        [
            {'type':'i4', 'shape':(1,), 'name':'unknown'},
            {'type':'i4', 'shape':(1,), 'name':'igwx'},
            {'type':'i4', 'shape':(1,), 'name':'ispin'},
        ],
        [
            {'type':'f8', 'shape':(3,3), 'name':'recipr'},
        ],
        [
            {'type':'i4', 'shape':('igwx',3,), 'name':'gvect'},
        ],
        ([
            {'type':'c16', 'shape':('igwx',), 'name':'C_kn'},
        ], 'ispin'),
    ]
    def __init__(self, **kwargs):
        self.igwx = None
        self.ispin = None
        self.recipr = None
        self.gvect = None
        self.C_kn = None

        super().__init__(**kwargs)

    @property
    def direct(self):
        return recipr_base(self.recipr)
