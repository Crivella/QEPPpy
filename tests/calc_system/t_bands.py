import os
import numpy as np
from .t_kpoints import Test_kpoints
from qepppy.calc_system.bands import bands

cwd = os.path.dirname(os.path.realpath(__file__))

class Test_bands(Test_kpoints):
	cls_typ = bands

	def test_bands_unfolding_noproject(self, cls_edges):
		from qepppy.qe import pw_out
		cls     = cls_edges
		file    = os.path.join(cwd, '../qe/test_files', 'bands_unf._out')
		calc    = pw_out(file=file)
		cls.egv = calc.egv

		x,y     = cls.band_unfolding_noproject(cls.recipr/2.)
		file    = os.path.join(cwd, '../qe/test_files', 'bands_unfold_noproject.dat')
		test    = np.loadtxt(file)

		assert np.allclose(x,test[:,0]),  "Failed to produce unfolded band structure without projection."
		assert np.allclose(y,test[:,1:], atol=1e-4), "Failed to produce unfolded band structure without projection."