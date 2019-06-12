import pytest
import os
import pickle
import qepppy
import numpy as np

class ElementMismatch(Exception):
	pass

def isproperty(cls, name):
	typ  = cls.__class__
	try:
		test = getattr(typ, name)
	except:
		return False
	if isinstance(test, property):
		return True

	return False

def compare_element(a,b):
	# print(a,b)
	if type(a) != type(b):
		raise ElementMismatch("Type of the two element is different")
	if isinstance(a, np.ndarray):
		if len(a) == len(b) == 0:
			return
		# if len(a) == 0 or len(b) == 0:
		# 	return
		if isinstance(a.flatten()[0], np.number):
			if not np.allclose(a,b, atol=1e-4):
				raise ElementMismatch("The two array are different (even considering noise!!!).")
		else:
			if not np.all(a == b):
				raise ElementMismatch("The two array are different.")
	elif isinstance(a, list):
		if len(a) != len(b):
			raise ElementMismatch("List of different lenghts.")
		for c1,c2 in zip(a,b):
			compare_element(c1,c2)
	elif isinstance(a, dict):
		if len(a) != len(b):
			raise ElementMismatch("List of different lenghts.")
		for k,v in a.items():
			if not k in b:
				raise ElementMismatch(f"dict is missing key '{k}'")
			compare_element(v, b[k])

cmp_list = [
	'direct', 'recipr',
	'n_atoms', 'atoms_coord_cart', 'atoms_coord_cryst', 'atoms_typ', 'atoms_mass', 'atoms_pseudo',
	'n_kpt', 'kpt_cart', 'kpt_cryst', 'weight',
	'egv',
	# 'symm_matrix', 'symm_name',
	'E_tot', 'cb', 'vb',
	'fft_dense_grid'
	]

def compare_std(cls, std):
	for name in cmp_list:
		# print("COMPARE" + "-"*20 + name)
		c1 = getattr(cls, name)
		c2 = std[name]
		compare_element(c1, c2)

cwd      = os.path.dirname(os.path.realpath(__file__))
test_dir = os.path.join(cwd, 'test_files')
os.chdir(test_dir)

out_files = [a[:-4] for a in os.listdir(test_dir) if a.endswith('.out')]

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

def test_pw_outs(out_v_pkl):
	out, xml, pkl = out_v_pkl

	compare_std(out, pkl)
	compare_std(xml, pkl)


if __name__ == '__main__':
	for out_name in out_files:
		xml_name = out_name[:-4] + '.xml'
		pkl_name = out_name[:-4] + '.pickle'

		xml = qepppy.qe.pw_out(xml=xml_name)

		res = {}
		for k in cmp_list:
			res[k] = getattr(xml, k)

		with open(pkl_name, 'wb') as f:
			pickle.dump(res, f)

