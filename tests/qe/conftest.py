import pytest
import numpy as np

class ElementMismatch(Exception):
	pass

def compare_element(a,b):
	if type(a) != type(b):
		if all(isinstance(c, (list, np.ndarray)) for c in [a,b]):
			a = np.array(a)
			b = np.array(b)
		else:
			raise ElementMismatch("Type of the two elements are different")
	if isinstance(a, np.ndarray):
		if len(a) == len(b) == 0:
			return
		if isinstance(a.flatten()[0], np.number):
			if not np.allclose(a,b, atol=1e-4):
				raise ElementMismatch("The two array are different (even considering noise!!!).")
		else:
			if not np.all(a == b):
				raise ElementMismatch("The two array are different.")
	elif isinstance(a, list):
		if len(a) != len(b):
			raise ElementMismatch("Two lists of different lenghts.")
		for c1,c2 in zip(a,b):
			compare_element(c1,c2)
	elif isinstance(a, dict):
		if len(a) != len(b):
			raise ElementMismatch("Two lists of different lenghts.")
		for k,v in a.items():
			if not k in b:
				raise ElementMismatch(f"dict is missing key '{k}'")
			compare_element(v, b[k])

def compare_std(cls, std, cmp_list=[]):
	for name in cmp_list:
		# print("COMPARE" + "-"*20 + name)
		c1 = getattr(cls, name)
		c2 = std[name]
		try:
			compare_element(c1, c2)
		except ElementMismatch as e:
			raise ValueError(f"While comparing '{name}': {e}")
		except Exception as e:
			raise type(e)(f"UNEXPECTED!!!: While comparing '{name}': {e}")



