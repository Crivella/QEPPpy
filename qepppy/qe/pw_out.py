import os

import numpy as np

from ..calc_system import system
from .qe_bands import qe_bands as bands
from .qe_structure import qe_structure as structure
from .tmp import tmp
from .UPF import UPF


# from ..utils       import xyz_mesh
# from .._decorators import store_property
# from ..logger import logger
# @logger()
# class pw_out(bands, structure):
class pw_out(structure, bands, system):
    """
    Instance used to handle QE outputs (by parsing the "data-file*.xml" file")
    fname: name of the "data-file*.xml" to parse
    kwargs:
     - outfile = Name of the pw outfile to parse
     - xml     = Name of the data-file*.xml to parse
    """
    __name__ = "pw_out"
    def __init__(self, *args, xml=None, **kwargs):
        self.set_data_file(xml)
        super().__init__(*args, xml=self.xml, **kwargs)
        self.validate()

    def set_data_file(self, xml):
        if xml is None:
            self.xml = xml
            return
        if os.path.isfile(xml):
            self.xml = xml
            self.data_path = os.path.dirname(os.path.realpath(xml))
        elif os.path.isdir(xml):
            file = os.path.join(xml, "data-file-schema.xml")
            if not os.path.isfile(file):
                raise FileNotFoundError(file)
            self.data_path = xml
            self.xml       = file
        else:
            raise FileNotFoundError("The file/folder {} does not exist".format(xml))
        if "data-file-schema.xml" in self.xml:
            self.prefix    = '.'.join(os.path.basename(self.data_path).split('.')[:-1])
            self.data_path = os.path.abspath( os.path.join(self.data_path, os.path.pardir))

    @property
    # @store_property
    def tmp(self):
        """Iterable. Every element is of type <class wavefnc> and contains the 
        data from the wafeunction binaries."""
        if not hasattr(self, 'prefix'):
            return []
        return tmp(self.prefix, path=self.data_path)

    @property
    # @store_property
    def pseudo(self):
        """Iterable. Every element is of type <class UPF> and contains the data
        of the atoms pseudopotentials.
        """
        pseudo = []
        if not hasattr(self, 'xml'):
            return pseudo
        path = os.path.dirname(self.xml)
        for pp in self.unique_atoms_pseudo:
            file = os.path.join(path,pp)
            pseudo.append(UPF(xml=file))
        return pseudo

    def calc_matrixelements(self, bnd_low=1, bnd_high=np.inf):
        fname = base = "matrixelements"
        i = 1
        while os.path.exists(fname):
            fname = base + '_' + str(i)
            i += 1
        bnd_low -= 1
        bnd_high = min(bnd_high, self.n_bnd)
        
        f = open(fname, "a")
        for k,psi in enumerate(self.tmp):
            egv = self.egv[k][bnd_low:bnd_high]
            occ = self.occ[k][bnd_low:bnd_high]*2
            kpt = self.kpt_cart[k].reshape(3,1)
            G   = np.dot(psi.recipr.T, psi.gvect.T)
            kG  = G + kpt
            for v in range(bnd_low, bnd_high):
                if occ[v-bnd_low] < 1E-4:
                    continue
                c  = np.where((2 - occ) > 1E-4)[0]
                c  = c[c>v]
                dE = egv[c] - egv[v-bnd_low]

                pp = np.sum(
                    np.conj(psi.C_kn[v]) * 
                    psi.C_kn[c + bnd_low].reshape(len(c), 1, psi.nspin, psi.igwx) * 
                    kG.reshape(1,3,1,psi.igwx), 
                    axis=(2,3))
                pp = np.real(np.conj(pp) * pp)

                res = np.column_stack((c+1 + bnd_low, pp, dE, occ[v-bnd_low]-occ[c]))

                fmt="{:5d}{:5d}".format(k+1,v+1) + "%5d" + "%16.8E"*3 + "%8.4f"*2
                np.savetxt(f, res, fmt=fmt)
                f.flush()

        f.close()

 
















