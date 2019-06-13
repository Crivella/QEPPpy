import pytest
import os
import pickle
import qepppy
from .conftest import compare_std

cwd      = os.path.dirname(os.path.realpath(__file__))
test_dir = os.path.join(cwd, 'test_files')
os.chdir(test_dir)

out_files = [a[:-4] for a in os.listdir(test_dir) if a.endswith('.out')]

cmp_list = [
	'direct', 'recipr',
	'ibrav', 'alat',
	'n_atoms', 'atoms_coord_cart', 'atoms_coord_cryst', 'atoms_typ', 'atoms_mass', # 'atoms_pseudo',
	'n_kpt', 'kpt_cart', 'kpt_cryst', 'weight',
	'egv',
	# 'symm_matrix', 'symm_name',
	'E_tot', 'cb', 'vb',
	'fft_dense_grid', 'fft_smooth_grid',
	]


@pytest.fixture(
	scope='module',
	params=out_files,
	)
def out_v_pkl(request):
	name     = request.param
	out_name = name + '.out'
	xml_name = name + '.xml'
	pkl_name = name + '.pickle'

	out = qepppy.qe.pw_out(outfile=out_name)
	xml = qepppy.qe.pw_out(xml=xml_name)

	with open(pkl_name, 'rb') as f:
		pkl = pickle.load(f)

	return out, xml, pkl
	

def test_pw_out_parsing(out_v_pkl):
	out, xml, pkl = out_v_pkl

	compare_std(out, pkl, cmp_list=cmp_list)
	compare_std(xml, pkl, cmp_list=cmp_list)


if __name__ == '__main__':
	"""
	Generates the baseline files for the automated testing.
	To be run with a stable pre-tested version.
	"""
	for out_name in out_files:
		xml_name = out_name + '.xml'
		pkl_name = out_name + '.pickle'

		xml = qepppy.qe.pw_out(xml=xml_name)

		res = {}
		for k in cmp_list:
			res[k] = getattr(xml, k)

		with open(pkl_name, 'wb') as f:
			pickle.dump(res, f)

