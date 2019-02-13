import os
import numpy as np
from .bands        import bands     as bands
from .structure    import structure as structure
from .tmp          import tmp
from .UPF          import UPF
from .._decorators import store_property
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

	@property
	@store_property
	def tmp(self):
		return tmp(self.prefix, path=self.data_path)

	@property
	@store_property
	def pseudo(self):
		pseudo = []
		path = os.path.dirname(self.schema)
		for pp in self.atoms_pseudo:
			file = os.path.join(path,pp)
			pseudo.append(UPF(schema=file))
		return pseudo

	def test_pseudo(self):
		import scipy
		for psi in self.tmp:
			psi.make_density_grid(bnd_list=range(1,int(self.n_el/2)+1))
			XYZ = psi._xyz_mesh_()
			print(XYZ.shape)
			for nt,coord,pp in zip(range(len(self.pseudo)), self.atoms_group_coord_cryst, self.pseudo):
				coord = coord.reshape(-1,3)

				val = np.zeros((pp.n_wfc, *XYZ[0].shape))
				rab = np.zeros(XYZ[0].shape)
				f_rab = scipy.interpolate.interp1d(pp.mesh[:pp.rab.size], pp.rab, 'nearest', bounds_error=False, fill_value=0)
				for c in coord:
					c = c.reshape(3,1,1,1)
					norm = np.linalg.norm(XYZ - c, axis=0)
					rab += f_rab(norm.flatten()).reshape(XYZ[0].shape)

					for nc,chi in enumerate(pp.pswfc):
						f_val = scipy.interpolate.interp1d(pp.mesh[:chi.size], chi, 'nearest', bounds_error=False, fill_value=0)
						val[nc] += f_val(norm.flatten()).reshape(XYZ[0].shape)
						# rab += f_rab(norm.flatten()).reshape(XYZ[0].shape)
						# print(pp.mesh[:chi.size])
						# print(norm.min(), norm.max())
						# print()
						# print(pp.mesh.min(), pp.mesh.max())
						# val += func(norm.flatten()).reshape(XYZ[0].shape)

				for i in range(pp.n_wfc):
					print("-------> nt: {}  nc: {}   sum: {}".format(nt+1, i+1, np.sum(psi.dgrid**.5 * val[i] * rab)))
				# print(pp.pswfc.shape)
				# print(coord)


			break
		pass

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

				pp = np.sum(np.conj(psi.val[v]) * psi.val[c + bnd_low].reshape(len(c),1,psi.igwx) * kG, axis=-1)
				pp = np.real(np.conj(pp) * pp)

				res = np.column_stack((c+1 + bnd_low, pp, dE, occ[v-bnd_low]-occ[c]))

				fmt="{:5d}{:5d}".format(k+1,v+1) + "%5d" + "%16.8E"*3 + "%8.4f"*2
				np.savetxt(f, res, fmt=fmt)
				f.flush()

		f.close()

 
















