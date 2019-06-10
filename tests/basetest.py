# import pytest
# import numpy as np

# def get_maker(
# 	name, svalue,
# 	check_equal=True,
# 	prev_call='',
# 	):
# 	def test_func(self, cls, value=svalue):
# 		if prev_call:
# 			func = getattr(self, prev_call)
# 			cls  = func(cls) 
# 		setattr(cls, name, value)
# 		check = getattr(cls, name)

# 		if check_equal:
# 			emsg = "Set and get produced different results!!"
# 			if isinstance(value, np.ndarray):
# 				assert np.all(value == check), emsg
# 			else:
# 				assert value == check, emsg

# 		return cls

# 	return test_func


# class meta_tester(type):
# 	def __new__(cls, clsname, base, dct):
# 		new_dct = {}
# 		for k,v in dct.items():
# 			if not isinstance(k, str) or not k.startswith("test_") or not isinstance(v, dict):
# 				new_dct[k] = v
# 			else:

# 				name          = '_'.join(k.split('_')[1:])

# 				f_name        = 'test_' + name

# 				doc           = v.pop('doc',         '')
# 				value         = v.pop('value',       None)
# 				check_equal   = v.pop('check_equal', True)
# 				prev_call     = v.pop('prev_call',   '')

# 				func          = get_maker(
# 					name, value,
# 					check_equal=check_equal,
# 					prev_call=prev_call,
# 					)
# 				func.__name__ = f_name
# 				func.__doc__  = doc

# 				new_dct[f_name] = func

# 		return super().__new__(cls, clsname, base, new_dct)


# class class_tester(metaclass=meta_tester):
# 	typ            = object
# 	init_args      = []
# 	init_kwargs    = {}
# 	run_empty_init = False

# 	def test_init(self):
# 		if self.run_empty_init:
# 			cls = self.typ()

# 		cls = self.typ(*self.init_args, **self.init_kwargs)
# 		assert isinstance(cls, self.typ), "Failed to create instance of" + str(self.typ)

# 		# pytest.cls = cls

# 	@pytest.fixture()
# 	def cls(self, *args, **kwargs):
# 		cls = self.typ(*args, **kwargs)

# 		return cls


