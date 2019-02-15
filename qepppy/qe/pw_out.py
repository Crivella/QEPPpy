import os
import scipy
import scipy.interpolate
import numpy as np
from .bands        import bands     as bands
from .structure    import structure as structure
from .tmp          import tmp
from .UPF          import UPF
# from ..utils       import xyz_mesh
from .. import utils
from .._decorators import store_property
# from ..logger import logger

# @logger()
class pw_out(bands, structure):
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

	def test_pdos_orthonormalized_states(self):
		names, states = self.test_pdos_states()

		# Generate overlap matrix O_ij = <phi_i | phi_j>
		states_p = [a.flat/np.linalg.norm(a) for a in states]
		O = np.empty((len(states),)*2, dtype='complex')
		for i in range(O.shape[0]):
			for j in range(i,O.shape[1]):
				O[i,j] = np.vdot(states_p[i], states_p[j])
			# if i != j:
			# 	O[j,i] = np.conj(O[i,j])
		del(states_p)
		np.set_printoptions(
			linewidth=2000,
			formatter={
				'complex_kind':lambda x: "{:7.4f} {:+7.4f}i".format(x.real,x.imag) if np.abs(x) > 1E-4 else "#"*16
				}
				)
		print(O)

		egval,egvec = np.linalg.eigh(O, 'U')
		egvec = np.mat(egvec)

		O_egvec_base = np.diag(1/egval**.5)
		new_O = np.dot(egvec, np.dot(O_egvec_base, egvec.I))

		states = np.tensordot(new_O, states, axes=1)
		return names, states


	def test_pdos_states(self):
		states_mesh = []
		states_name = []
		XYZ = utils.xyz_mesh(
			np.array(self.fft_dense_grid)//2,
			base = self.cell
			)
		"""
		Instead of generating everythin at the coordinates of the atom, in order to
		take into account the periodicity of the system, generate everything at the center
		and than roll the grid to move the center on the atom position
		"""
		center = np.dot(self.cell.T, [.5,.5,.5])
		dXYZ   = np.array((
			XYZ[0,1,0,0] - XYZ[0,0,0,0],
			XYZ[1,0,1,0] - XYZ[1,0,0,0],
			XYZ[2,0,0,1] - XYZ[2,0,0,0],
			))
		grid_shape = XYZ[0].shape
		na = 0
		for nt,(name,pp) in enumerate(zip(self.atoms_typ, self.pseudo)):
			coord = self.atoms_group_coord_cart[name]
			coord = coord.reshape(-1,3)
			for c in coord:
				na   += 1
				delta = (c - center)/dXYZ
				delta = delta.astype(dtype='int')

				cXYZ  = np.array(XYZ - center.reshape(3,1,1,1))
				norm  = np.linalg.norm(cXYZ, axis=0)
				theta = np.arctan2(cXYZ[1],cXYZ[0])
				theta[theta<0] += 2*np.pi
				phi   = np.arccos(cXYZ[2]/norm)
				for nc,(l,chi) in enumerate(zip(pp.pswfc_l,pp.pswfc)):
					f_val = scipy.interpolate.interp1d(pp.mesh[:chi.size], chi,
						'linear',
						bounds_error=False,
						fill_value=0
						)
					val = f_val(norm.flat).reshape(grid_shape)
					val = val.astype(dtype='complex')
					# val = scipy.fftpack.ifftn(val)
					for m in range(-l,l+1):
						harm = scipy.special.sph_harm(m,l,theta,phi)
						# rrrr = self.fft_dense_grid[2]//4
						# print("rrrr:", rrrr)
						# import matplotlib.pyplot as plt
						# fig, ax = plt.subplots(2,2)
						# p = ax[0,0].contourf(XYZ[0,:,:,rrrr],XYZ[1,:,:,rrrr],harm[:,:,rrrr].real,100)
						# ax[0,1].contourf(XYZ[0,:,:,rrrr],XYZ[1,:,:,rrrr],val[:,:,rrrr].real)
						# ax[1,0].contourf(XYZ[0,:,:,rrrr],XYZ[1,:,:,rrrr],(val * harm)[:,:,rrrr].real)

						# ax[1,1].contourf(XYZ[0,:,:,0],XYZ[1,:,:,0],np.roll(res, delta, axis=(0,1,2))[:,:,0].real)
						# ax[1,1].scatter(c[0],c[1],color='r')
						# # plt.colorbar(p)
						# plt.suptitle("l={} m={}  z={}".format(l,m,rrrr))
						# plt.show()
						states_mesh.append(
							np.roll(val * harm, delta, axis=(0,1,2))
							)
						states_name.append("atom{:>4d} ({:<3s}), wfc{:>3d} (l={} m={:2d})".format(na,name,nc,l,m))
					# states_name.append("atom{:>4d} ({:<3s}), wfc{:>3d}".format(na,name,nc))

		return states_name,states_mesh
		
	def test_pdos(self):
		# nlist, slist = self.test_pdos_orthonormalized_states()
		nlist, slist = self.test_pdos_states()
		XYZ = utils.xyz_mesh(
			self.fft_dense_grid//2,
			base = self.cell
			)
		grid_shape = slist[0].shape
		for nk,psi in enumerate(self.tmp):
			print("KPT (#{:>4d}): {}".format(nk+1, self.kpt_cart[nk]))
			for nb in range(self.n_bnd):	
				print("  BND (#{:>3d}): {:14.6f} eV".format(nb+1,self.egv[nk,nb]))
				dgrid = psi._generate_g_grid_(band=nb, shape=grid_shape)
				# dgrid = np.abs(scipy.fftpack.fftn(dgrid))
				# dgrid = scipy.fftpack.fftn(dgrid).real
				dgrid = scipy.fftpack.fftn(dgrid)

				for name,val in zip(nlist, slist):
					comp = np.abs(np.vdot(dgrid.flat, val.flat) / (np.linalg.norm(dgrid) * np.linalg.norm(val)))**2 
					print("    {}: {:8.5f}%".format(name, comp * 100))
					import matplotlib.pyplot as plt
					for i in range(-5,min(XYZ.shape[-1],6),5):
						x = XYZ[0][:,:,i]
						y = XYZ[1][:,:,i]
						# z = np.abs(val[:,:,i])
						z = val[:,:,i].real
						# fig, ax = plt.subplots(1,2)
						# p = ax[0].contourf(x,y,z,100,cmap='inferno')
						# ax[1].contourf(x,y,np.abs(dgrid[:,:,i]),100,cmap='inferno')
						# ax[0].set_title(name)
						# ax[1].set_title('dgrid')
						# fig.suptitle(str(i))
						# fig.colorbar(p)
						# plt.show()
			break


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

 
















