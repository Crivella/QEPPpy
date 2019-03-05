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

	def test_pdos_states(self, lmax=1, nn=(1,1,1)):
		states_mesh = []
		states_name = []
		print(self.direct)

		start_shape = self.fft_dense_grid//2
		print(start_shape)
		# nn = (1,1,1)#(3,3,1)
		nn = np.array(nn)
		XYZ = utils.xyz_mesh(
			start_shape,
			rep=nn,
			base = self.direct
			)
		print(np.unique(np.round(XYZ[0],decimals=5)))
		"""
		Instead of generating everything at the coordinates of the atom, in order to
		take into account the periodicity of the system, generate everything at the center
		and than roll the grid to move the center on the atom position
		"""
		center = np.array(np.dot(self.direct.T, np.array([.5,.5,.5])*nn))
		s = tuple((slice(None),*((nn*start_shape)//2-1)))
		print("CENTER: calc {}    fromXYZ {}".format(center, XYZ[s]))
		grid_shape = XYZ[0].shape
		print(grid_shape)
		# print(XYZ.shape)
		# exit()

		from itertools import product
		ll = product(*(range(-np.ceil((nr-1)/2).astype('int'),np.ceil((nr-1)/2).astype(int)+1) for nr in nn))
		print(list(ll))
		# exit()
		# for nr in nn:
		# 	ll = np.ceil((nn-1)/2)
		# 	for lr in range(-ll,ll+1):

		cXYZ  = XYZ - center.reshape(3,1,1,1)
		norm  = np.linalg.norm(cXYZ, axis=0) + 1E-16
		max_norm = np.min(
			np.hstack((
			norm[0,:,:],
			norm[:,:,0],
			norm[:,0,:],
			))
			)*10.98

		# test_harm = np.zeros(((lmax+1)**2, *norm.shape), dtype=complex)
		# for rx,ry,rz in ll:
		# cXYZ  = XYZ - center.reshape(3,1,1,1) + np.dot(self.direct.T,[rx,ry,rz])
		theta = np.arctan2(cXYZ[1],cXYZ[0])
		theta[theta<0] += 2*np.pi
		phi   = np.arccos(cXYZ[2]/norm)

		test_harm = []
		for l in range(lmax+1):
			for m in range(-l,l+1):
				# harm = test_harm[l**2 + m + l]
				harm = scipy.special.sph_harm(m,l,theta,phi) # * 1j**l
				harm[norm > max_norm] = 0
				# print(harm[:20])
				# test_harm[l**2 + m + l] = harm
				harm /= np.linalg.norm(harm)
				test_harm.append(harm)

		# print(test_harm.shape)
		test_harm = np.array([a/np.linalg.norm(a) for a in test_harm])

		th = np.array([a.flatten() for a in test_harm])
		print("TEST sph_harm overlap: \n", np.dot(np.conj(th),th.T))
		# exit()

		na = 0
		for nt,(name,pp) in enumerate(zip(self.atoms_typ, self.pseudo)):
			print(name, "<-->", pp.xml)
			coord = self.atoms_group_coord_cart[name]
			coord = coord.reshape(-1,3)
			print(coord)


			for c in coord:
				na   += 1
				delta = np.dot(
					scipy.linalg.inv(self.direct.T),
					# self.recipr/(2*np.pi),
					(c-center)*start_shape
					)
				print("DELTA", delta)
				delta = np.around(delta).astype(int)
				print("DELTA", delta)
				abc = ((start_shape-1)//2 + delta)%(start_shape)
				abc = tuple((slice(None),*abc))
				print("Post_ROLL:", XYZ[abc])


				for nc,(l,chi) in enumerate(zip(pp.pswfc_l,pp.pswfc)):
					f_val = scipy.interpolate.interp1d(
						pp.mesh[:chi.size], chi,
						'cubic',
						bounds_error=False,
						fill_value=0
						)
					val = f_val(norm.flat).reshape(grid_shape)
					val = val.astype(dtype='complex')
					for m in range(-l,l+1):
						harm = test_harm[l**2 + m + l]

						app = val * harm
						# for d in range(3):
						# 	for i in range(nn[d]-1,-1,-1):
						# 		i *= (np.array(grid_shape)//nn)[d]
						# 		s = [slice(None)]*3
						# 		s[d] = i-1
						# 		app = np.insert(app,i,app[tuple(s)],axis=d)
						# # app = app.reshape(
						# # 	nn[0],start_shape[0]+1,
						# # 	nn[1],start_shape[1]+1,
						# # 	nn[2],start_shape[2]+1,
						# # 	).sum(axis=(0,2,4))
						# # print(app.shape)
						# app = app.reshape(
						# 	nn[0],start_shape[0]+1,
						# 	nn[1],start_shape[1]+1,
						# 	nn[2],start_shape[2]+1,
						# 	)
						app = app.reshape(
							nn[0],start_shape[0],
							nn[1],start_shape[1],
							nn[2],start_shape[2],
							)
						# print(app.shape)
						app = app.sum(axis=(0,2,4))#[1:,1:,1:]
						# print(app.shape)
						app2 = np.roll(app, delta, axis=(0,1,2))

						###########################################################################
						# Test plot sph_harm and states
						# x, y, index = utils.remap_plane(
						# 	self.recipr/ nn.reshape(3,1) / (2*np.pi),
						# 	(XYZ[0].min(), XYZ[0].max()),
						# 	(XYZ[1].min(), XYZ[1].max()),
						# 	(XYZ[2].min(), XYZ[2].max()),
						# 	XYZ.shape[1:], (1,)*3
						# 	)

						# rrrr = start_shape[2]//2
						# import matplotlib.pyplot as plt
						# aXYZ = utils.xyz_mesh(
						# 	start_shape,
						# 	rep=1,
						# 	base = self.direct
						# 	)
						# fig, ax = plt.subplots(2,2, figsize=(15,10))
						# ax[0,0].contourf(x,y,harm[index].real.reshape(x.shape),100,cmap='seismic')
						# ax[0,0].scatter(x,y,
						# 	color='k', s=1
						# 	)
						# ax[0,1].contourf(x,y,val[index].real.reshape(x.shape),cmap='seismic')
						# ax[0,1].scatter(x,y,
						# 	color='k', s=1
						# 	)

						# x, y, index = utils.remap_plane(
						# 	self.recipr/ nn.reshape(3,1) / (2*np.pi),
						# 	(aXYZ[0].min(), aXYZ[0].max()),
						# 	(aXYZ[1].min(), aXYZ[1].max()),
						# 	(aXYZ[2].min(), aXYZ[2].max()),
						# 	aXYZ.shape[1:], (1,)*3
						# 	)
						# ax[1,0].contourf(x,y,
						# 	app[index].real.reshape(x.shape),
						# 	cmap='seismic'
						# 	)
						# ax[1,0].scatter(x,y,
						# 	color='k', s=.5
						# 	)
						# z = c[2]
						# x, y, index = utils.remap_plane(
						# 	self.recipr/ nn.reshape(3,1) / (2*np.pi),
						# 	(aXYZ[0].min(), aXYZ[0].max()),
						# 	(aXYZ[1].min(), aXYZ[1].max()),
						# 	(z,z),
						# 	aXYZ.shape[1:], (1,)*3
						# 	)

						# toplot = app2[index].real.reshape(x.shape)
						# p = ax[1,1].contourf(x,y,toplot,
						# 	cmap='seismic',
						# 	vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max()
						# 	)
						# ax[1,1].scatter(x,y,
						# 	color='k', s=.5
						# 	)
						# ax[1,1].scatter(c[0],c[1],color='r')
						# # for ii in range(4):
						# # 	for i in range(XYZ[0].shape[0]):
						# # 		ax[ii//2,ii%2].plot([XYZ[0,i,0,0],XYZ[0,i,-1,0]], [XYZ[1,i,0,0],XYZ[1,i,-1,0]], color='k')
						# # 		ax[ii//2,ii%2].plot([XYZ[0,0,i,0],XYZ[0,-1,i,0]], [XYZ[1,0,i,0],XYZ[1,-1,i,0]], color='k')
						# ax[0,0].set_title("l={} m={}".format(l,m,rrrr))
						# ax[0,1].set_title("CHI")
						# ax[1,0].set_title("CHI * sph_harm + reshape_sum")
						# ax[1,1].set_title("After ROLL")
						# fig.colorbar(p)
						# plt.show()
						###########################################################################
						states_mesh.append(
							app2
							# app
							)
						states_name.append("atom{:>4d} ({:<3s}), wfc{:>3d} (l={} m={:2d})".format(na,name,nc,l,m))

		return states_name, np.array(states_mesh)
		
	def test_pdos(self, thr=0.01, **kwargs):
		nlist, slist = self.test_pdos_states(**kwargs)
		slist = utils.lowdin_ortho(slist)
		# print(slist.shape)
		XYZ = utils.xyz_mesh(
			self.fft_dense_grid//2,
			base = self.direct
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

				tot = 0
				for name,val in zip(nlist, slist):
					comp = np.abs(np.vdot(dgrid.flat, val.flat) / (np.linalg.norm(dgrid) * np.linalg.norm(val)))**2
					tot += comp
					if comp < thr:
						continue
					print("    {}: {:6.3f}%".format(name, comp * 100))

					###########################################################################
					# Test plot
					# import matplotlib.pyplot as plt
					# fig,ax = plt.subplots(3,3,figsize=(20, 10))
					# for ai,i in enumerate(range(-5,min(XYZ.shape[-1],6),5)):
					# 	x = XYZ[0][:,:,i]
					# 	y = XYZ[1][:,:,i]
					# 	# z = np.abs(val[:,:,i])
					# 	z = val[:,:,i]

					# 	toplot = z.real
					# 	p1 = ax[0,ai].contourf(x,y,toplot,100,
					# 		cmap='seismic',
					# 		vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max()
					# 		)
					# 	cb = fig.colorbar(p1, ax=ax[0,ai])
					# 	cb.set_clim(vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max())

					# 	toplot = dgrid[:,:,i].real
					# 	p2 = ax[1,ai].contourf(x,y,toplot,100,
					# 		cmap='seismic',
					# 		vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max()
					# 		)
					# 	# ax[2,ai].contourf(x,y,(np.conj(dgrid[:,:,i])*z).real,100,cmap='seismic')
					# 	fig.colorbar(p2, ax=ax[1,ai])

					# 	toplot = (np.conj(dgrid[:,:,:])*val).real.sum(axis=2) / (np.linalg.norm(dgrid) * np.linalg.norm(val))
					# 	p3 = ax[2,ai].contourf(x,y,toplot,100,
					# 		cmap='seismic',
					# 		vmin=-np.abs(toplot).max(),vmax=np.abs(toplot).max()
					# 		)
					# 	fig.colorbar(p3, ax=ax[2,ai])

					# 	for ii in range(3):
					# 		for i in range(XYZ[0].shape[0]):
					# 			ax[ii,ai].plot([XYZ[0,i,0,0],XYZ[0,i,-1,0]], [XYZ[1,i,0,0],XYZ[1,i,-1,0]], color='k')
					# 			ax[ii,ai].plot([XYZ[0,0,i,0],XYZ[0,-1,i,0]], [XYZ[1,0,i,0],XYZ[1,-1,i,0]], color='k')

					# 	ax[0,ai].set_title(name)
					# 	ax[1,ai].set_title('dgrid')
					# 	ax[0,ai].set_title(name)
					# 	ax[2,ai].set_title('prod')
					# 	ai += 1
					# 	# ax[0].set_title(str(i))
					# plt.show()
					###########################################################################
				print("{:37s}: {:6.3f}%".format("Total", tot*100))
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

 
















