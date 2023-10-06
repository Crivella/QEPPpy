import numpy as np
import pytest


class ElementMismatch(Exception):
    pass

def compare_element(a,b):
    if type(a) != type(b):
        if any(isinstance(c, np.ndarray) for c in [a,b]) and any(isinstance(c, list) for c in [a,b]):
            a = np.array(a)
            b = np.array(b)
        else:
            raise ElementMismatch(f'Type of the two elements are different. {type(a)} <--> {type(b)}')
    if isinstance(a, np.ndarray):
        if a.shape != b.shape:
            raise ElementMismatch(f"The shape of the two array does not match: '{a.shape}' vs '{b.shape}'.")
        if len(a) == 0:
            return
        if isinstance(a.flatten()[0], np.number):
            if not np.allclose(a,b, atol=1e-4):
                raise ElementMismatch('The two array are different (even considering noise!!!).')
        else:
            if not np.all(a == b):
                raise ElementMismatch('The two array are different.')
    elif isinstance(a, list):
        if len(a) != len(b):
            raise ElementMismatch(f"Two lists of different lenghts: '{len(a)}' vs '{len(b)}'.")
        for c1,c2 in zip(a,b):
            compare_element(c1,c2)
    elif isinstance(a, dict):
        if len(a) != len(b):
            raise ElementMismatch(f"Two lists of different lenghts: '{len(a)}' vs '{len(b)}'.")
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
            raise ElementMismatch(f"While comparing '{name}': \n\t{str(e)}")
        except Exception as e:
            raise type(e)(f"UNEXPECTED!!!: While comparing '{name}': {e}")
