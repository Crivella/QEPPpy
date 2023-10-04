import os
import numpy as np
from re      import split
from qepppy._decorators import numpy_plot_opt, numpy_save_opt, join_doc

orb_types = ('s', 'p', 'd', 'f') 
max_orb = len(orb_types)

@numpy_plot_opt(_xlab="Energy (eV)", _ylab="PDOS (arb. units)")
@numpy_save_opt()
def _sum_pdos_(atom, path=".", deg=0.0, **kwargs):
	"""
	Use the output files from a PDOS (projwfc.x) calculation to calculated the 
	PDOS resolved per atom name and per orbital angular momentum (l).
	"""
	res = None
	print("-"*80)
	print(atom)
	for f in os.listdir(path):
		if "pdos_atm#" in f and "({})".format(atom) in f:
			name = os.path.join(path, f)
			print("READING: {}".format(name))
			data = np.loadtxt(fname=name, comments="#")

			if res is None:
				n_pt = data.shape[0]
				res = np.zeros((n_pt, max_orb+1))
				res[:,0] = data[:,0]
			l = list(filter(None, split(r".+atm#|\(.+#|\(|\)|\_j", f)))
			# wfc_n = int(l[1])-1
			try:
				orb_n = int(orb_types.index(l[2]))
			except ValueError:
				raise Exception("Non recognize orbital type {}".format(l[2]))

			res[:,orb_n+1] += data[:,1]

	if deg > 0:
		from .broad import broad
		res = broad(res, 'gauss', deg)

	return res

def sum_pdos(atoms="", path=".", deg=0.0, **kwargs):
	import matplotlib.pyplot as plt
	if not atoms:
		raise ValueError("Passing empty name for 'atoms'")

	if isinstance(atoms, str):
		atoms_l = atoms.split(",")
	elif isinstance(atoms, list):
		atoms_l = atoms
	else:
		raise NotImplementedError()

	fig, ax = plt.subplots()
	res = []
	for atom in atoms_l:
		fname  = "{}_total_pdos_orb.dat".format(atom)
		labels = [atom + '_' + a for a in orb_types]
		res.append(_sum_pdos_(atom, path=path, deg=deg, ax=ax, fname=fname, labels=labels, **kwargs))

	fig.legend()
	plt.show()

	return np.array(res)

join_doc(sum_pdos, _sum_pdos_.__doc__)

def main():
	import sys
	argc = len(sys.argv)
	if(not 2<=argc<=5):
		print("Incorrect use. Please pass arguments:"
			"\n\t'atom_list\t(comma separated e.g.: \"C,H,Cl\")',"
			"\n\t'path\t\t(default = \".\")',"
			"\n\t'plot\t\t0/1 (default = 1 = True) (plot return value using matplotlib'"
			"\n\t'pfile\t\t0/1 (default = 1 = True) (save return value output to file)")
		exit()
	if(argc==2):
		sum_pdos(sys.argv[1])
	if(argc==3):
		sum_pdos(sys.argv[1], sys.argv[2])
	if(argc==4):
		sum_pdos(sys.argv[1], sys.argv[2], bool(sys.argv[3]))
	if(argc==5):
		sum_pdos(sys.argv[1], sys.argv[2], bool(sys.argv[3]), bool(sys.argv[4]))

if __name__ == "__main__":
	main()









