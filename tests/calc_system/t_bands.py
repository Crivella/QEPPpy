import numpy as np
from .t_kpoints import Test_kpoints
from qepppy.calc_system.bands import bands

class Test_bands(Test_kpoints):
	cls_typ = bands

	def test_bands_unfolding_noproject(self, cls_edges):
		from qepppy.qe import pw_out
		cls     = cls_edges
		calc    = pw_out(outfile='bands_unf._out')
		cls.egv = calc.egv

		x,y     = cls.band_unfolding_noproject(cls.recipr/2.)
		test    = np.loadtxt('bands_unfold_noproject.dat')

		assert np.allclose(x,test[:,0]),  "Failed to produce unfolded band structure without projection."
		assert np.allclose(y,test[:,1:], atol=1e-4), "Failed to produce unfolded band structure without projection."