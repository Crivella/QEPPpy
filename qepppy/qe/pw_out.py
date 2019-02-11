import os
import numpy as np
from .bands     import bands     as bands
from .structure import structure as structure
from .tmp       import tmp
# from ..logger import logger

# @logger()
class pw_out( bands, structure):
	"""
	Instance used to handle QE outputs (by parsing the "data-file*.xml" file")
	fname: name of the "data-file*.xml" to parse
	kwargs:
	 - schema = Name of the data-file*.xml to parse
	"""
	__name__ = "pw_out"
	def __init__( self, **kwargs):
		super().__init__( **kwargs)
		self.validate()

	def load_tmp(self):
		try:
			self.tmp
		except AttributeError:
			self.tmp = tmp(self.prefix, path=self.data_path)

	def calc_matrixelements(self, bnd_low=1, bnd_high=None):
		self.load_tmp()
		fname = base = "matrixelements"
		i = 1
		while os.path.exists(fname):
			fname = base + '_' + str(i)
			i += 1
		bnd_low -= 1
		if bnd_high is None:
			bnd_high = self.n_bnd
		f = open(fname, "a")
		for k,psi in enumerate(self.tmp):
			egv = self.egv[k][bnd_low:bnd_high]
			occ = self.occ[k][bnd_low:bnd_high]*2
			kpt = self.kpt_cart[k].reshape(3,1)
			G   = np.dot(psi.recipr.T, psi.gvect.T)
			for v in range(bnd_low, bnd_high):
				if occ[v-bnd_low] < 1E-4:
					continue
				c  = np.where((2 - occ) > 1E-4)[0]
				dE = egv[c] - egv[v-bnd_low]

				pp = np.sum(np.conj(psi.val[v]) * psi.val[c + bnd_low].reshape(len(c),1,psi.igwx) * (G + kpt), axis=-1)
				pp = np.real(np.conj(pp) * pp)

				res = np.column_stack((c+1 + bnd_low, pp, dE, occ[v-bnd_low]-occ[c]))

				fmt="{:5d}{:5d}".format(k+1,v+1) + "%5d" + "%16.8E"*3 + "%8.4f"*2
				np.savetxt(f, res, fmt=fmt)
				f.flush()

			del(egv,occ,kpt,G,psi)

		f.close()

 
















