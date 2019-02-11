import os
import numpy as np
from .pw_out  import pw_out
from .wavefnc import wavefnc

class tmp():
	def __init__(self, prefix, path="."):
		self.current   = 0
		self.prefix    = prefix
		self.path      = path
		self.data_path = os.path.join(self.path, "{}.save".format(self.prefix))
		self.data      = pw_out(schema=os.path.join(self.data_path, 'data-file-schema.xml'))

	def __iter__(self):
		return self

	def __next__(self):
		if self.current < self.data.n_kpt:
			self.current += 1
			return self._get_wfc_num_(self.current)
		self.current = 0
		raise StopIteration

	def _get_wfc_num_(self, n):
		file = os.path.join(self.data_path, "wfc{}.dat".format(n))
		return wavefnc(src=file)

	def calc_matrixelements(self, bnd_low=1, bnd_high=None):
		import os
		fname = base = "matrixelements"
		i = 1
		while os.path.exists(fname):
			fname = base + '_' + str(i)
			i += 1
		bnd_low -= 1
		if bnd_high is None:
			bnd_high = self.data.n_bnd
		f = open(fname, "a")
		for k,psi in enumerate(self):
			egv = self.data.egv[k][bnd_low:bnd_high]
			occ = self.data.occ[k][bnd_low:bnd_high]*2
			kpt = self.data.kpt_cart[k].reshape(3,1)
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


