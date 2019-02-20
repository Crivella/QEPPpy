import os
# import numpy as np
# from .pw_out  import pw_out
from .wavefnc import wavefnc

class tmp():
	def __init__(self, prefix, path="."):
		self.current   = 0
		self.prefix    = prefix
		self.path      = path
		self.data_path = os.path.join(self.path, "{}.save".format(self.prefix))
		# self.data      = pw_out(xml=os.path.join(self.data_path, 'data-file-schema.xml'))

	def __iter__(self):
		return self

	def __next__(self):
		self.current += 1
		try:
			return self._get_wfc_num_(self.current)
		except FileNotFoundError:
			self.current = 0
			raise StopIteration
		raise NotImplementedError()

	def __getitem__(self,key):
		if not isinstance(key, int):
			raise NotImplementedError()
		return self._get_wfc_num_(key)

	def _get_wfc_num_(self, n):
		file = os.path.join(self.data_path, "wfc{}.dat".format(n))
		return wavefnc(src=file)



