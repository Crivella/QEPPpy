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
		path = os.path.dirname(self.xml)
		for pp in self.atoms_pseudo:
			file = os.path.join(path,pp)
			pseudo.append(UPF(xml=file))
		return pseudo

	def test_pdos_orthonormalized_states(self):
		names, states = self.test_pdos_states()

		# Orthonormalize the base using formula:
		# M_orth = (M^T . M)^{-1/2} . M
		# Where M is a matrix with base vectors on the rows
		shape  = states[0].shape
		states = np.array([a.flat/np.linalg.norm(a) for a in states])
		states = np.dot(
			scipy.linalg.inv(scipy.linalg.sqrtm(np.dot(np.conj(states), states.T))),
			states
			)
		return names, states.reshape(-1,*shape)


	def test_pdos_states(self):
		states_mesh = []
		states_name = []
		print(self.cell)
		XYZ = utils.xyz_mesh(
			self.fft_dense_grid//2,
			base = self.cell
			)
		"""
		Instead of generating everything at the coordinates of the atom, in order to
		take into account the periodicity of the system, generate everything at the center
		and than roll the grid to move the center on the atom position
		"""
		center = np.dot(self.cell.T, [.5,.5,.5])
		# dXYZ   = np.array((
		# 	XYZ[0,1,0,0] - XYZ[0,0,0,0],
		# 	XYZ[1,0,1,0] - XYZ[1,0,0,0],
		# 	XYZ[2,0,0,1] - XYZ[2,0,0,0],
		# 	))
		# print("dXYZ: ", dXYZ)
		grid_shape = XYZ[0].shape

		cXYZ  = np.array(XYZ - center.reshape(3,1,1,1))
		norm  = np.linalg.norm(cXYZ, axis=0)+1E-16
		max_norm = np.min(
			np.hstack((
			norm[0,:,:],
			norm[:,:,0],
			norm[:,0,:],
			# norm[-1,-1,:],
			# norm[-1,:,-1],
			# norm[:,-1,-1],
			))
			)
		theta = np.arctan2(cXYZ[1],cXYZ[0])
		theta[theta<0] += 2*np.pi
		phi   = np.arccos(cXYZ[2]/norm)

		na = 0
		for nt,(name,pp) in enumerate(zip(self.atoms_typ, self.pseudo)):
			print(name, "<-->", pp.schema)
			coord = self.atoms_group_coord_cart[name]
			coord = coord.reshape(-1,3)
			print(coord)
			for c in coord:
				na   += 1
				delta = np.dot(
					scipy.linalg.inv(self.cell.T/(self.fft_dense_grid//2)),
					c-center
					)
				delta = delta.astype(dtype='int')
				print("DELTA", delta)

				test_harm = []
				for nc,(l,chi) in enumerate(zip(pp.pswfc_l,pp.pswfc)):
					f_val = scipy.interpolate.interp1d(pp.mesh[:chi.size], chi,
						'linear',
						bounds_error=False,
						fill_value=0
						)
					val = f_val(norm.flat).reshape(grid_shape)
					val[norm > max_norm] = 0
					val = val.astype(dtype='complex')
					# val = scipy.fftpack.ifftn(val)
					for m in range(-l,l+1):
						harm = scipy.special.sph_harm(m,l,theta,phi) # * 1j**l
						harm /= np.linalg.norm(harm)
						harm[norm > max_norm] = 0
						test_harm.append(harm.flatten())

						###########################################################################
						# Test plot sph_harm and states
						# rrrr = self.fft_dense_grid[2]//4
						# import matplotlib.pyplot as plt
						# fig, ax = plt.subplots(2,2, figsize=(15,10))
						# ax[0,0].contourf(XYZ[0,:,:,rrrr],XYZ[1,:,:,rrrr],harm[:,:,rrrr].real,100,cmap='seismic')
						# ax[0,1].contourf(XYZ[0,:,:,rrrr],XYZ[1,:,:,rrrr],val[:,:,rrrr].real,cmap='seismic')
						# ax[1,0].contourf(XYZ[0,:,:,rrrr],XYZ[1,:,:,rrrr],(val * harm)[:,:,rrrr].real,cmap='seismic')

						# toplot = np.roll(val * harm, delta, axis=(0,1,2))[:,:,(rrrr+delta[2])%val.shape[2]].real
						# p = ax[1,1].contourf(XYZ[0,:,:,0],XYZ[1,:,:,0],toplot,
						# 	cmap='seismic',
						# 	vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max()
						# 	)
						# ax[1,1].scatter(c[0],c[1],color='r')
						# for i in range(XYZ[0].shape[0]):
						# 	ax[1,1].plot([XYZ[0,i,0,0],XYZ[0,i,-1,0]], [XYZ[1,i,0,0],XYZ[1,i,-1,0]], color='k')
						# 	ax[1,1].plot([XYZ[0,0,i,0],XYZ[0,-1,i,0]], [XYZ[1,0,i,0],XYZ[1,-1,i,0]], color='k')
						# ax[0,0].set_title("l={} m={}".format(l,m,rrrr))
						# ax[0,1].set_title("CHI")
						# ax[1,0].set_title("CHI * sph_harm")
						# ax[1,1].set_title("After ROLL")
						# fig.colorbar(p)
						# plt.show()
						###########################################################################

						states_mesh.append(
							np.roll(val * harm, delta, axis=(0,1,2))
							)
						states_name.append("atom{:>4d} ({:<3s}), wfc{:>3d} (l={} m={:2d})".format(na,name,nc,l,m))
					# states_name.append("atom{:>4d} ({:<3s}), wfc{:>3d}".format(na,name,nc))
				th = np.array(test_harm)
				th[:,np.where(norm > max_norm)[0]] = 0
				np.set_printoptions(
					linewidth=2000,
					formatter={
						'complex_kind':lambda x: "{:5.2f} {:+5.2f}i".format(x.real,x.imag) if np.abs(x) > 1E-2 else ":"*12
						}
					)
				print("TEST sph_harm overlap: \n", np.dot(np.conj(th),th.T))
		return states_name,states_mesh
		
	def test_pdos(self):
		# nlist, slist = self.test_pdos_orthonormalized_states()
		nlist, slist = self.test_pdos_states()
		# print(slist.shape)
		XYZ = utils.xyz_mesh(
			self.fft_dense_grid//2,
			base = self.cell
			)
		grid_shape = slist[0].shape
		for nk,psi in enumerate(self.tmp):
			print("KPT (#{:>4d}): {}".format(nk+1, self.kpt_cart[nk]))
			for nb in range(self.n_bnd):	
				print("  BND (#{:>3d}): {:14.6f} eV".format(nb+1,self.egv[nk,nb]))
				dgrid = psi.make_psi_grid(bnd=nb+1, shape=grid_shape)
				# dgrid = np.abs(scipy.fftpack.fftn(dgrid))
				# dgrid = scipy.fftpack.fftn(dgrid).real
				# dgrid = scipy.fftpack.fftn(dgrid)

				for name,val in zip(nlist, slist):
					comp = np.abs(np.vdot(dgrid.flat, val.flat) / (np.linalg.norm(dgrid) * np.linalg.norm(val)))**2 
					print("    {}: {:8.5f}%".format(name, comp * 100))

					###########################################################################
					# Test plot
					import matplotlib.pyplot as plt
					fig,ax = plt.subplots(3,3,figsize=(20, 10))
					ai = 0
					for i in range(-5,min(XYZ.shape[-1],6),5):
						x = XYZ[0][:,:,i]
						y = XYZ[1][:,:,i]
						# z = np.abs(val[:,:,i])
						z = val[:,:,i]

						toplot = z.real
						p1 = ax[0,ai].contourf(x,y,toplot,100,
							cmap='seismic',
							vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max()
							)
						cb = fig.colorbar(p1, ax=ax[0,ai])
						cb.set_clim(vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max())

						toplot = dgrid[:,:,i].real
						p2 = ax[1,ai].contourf(x,y,toplot,100,
							cmap='seismic',
							vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max()
							)
						# ax[2,ai].contourf(x,y,(np.conj(dgrid[:,:,i])*z).real,100,cmap='seismic')
						fig.colorbar(p2, ax=ax[1,ai])

						toplot = (np.conj(dgrid[:,:,:])*val).real.sum(axis=2) / (np.linalg.norm(dgrid) * np.linalg.norm(val))
						p3 = ax[2,ai].contourf(x,y,toplot,100,
							cmap='seismic',
							vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max()
							)
						fig.colorbar(p3, ax=ax[2,ai])

						ax[0,ai].set_title(name)
						ax[1,ai].set_title('dgrid')
						ax[0,ai].set_title(name)
						ax[2,ai].set_title('prod')
						ai += 1
						# ax[0].set_title(str(i))
					plt.show()
					###########################################################################
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

 
















