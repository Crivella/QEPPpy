import os

import numpy as np
import pytest

from qepppy.parsers.regex_parser import Parser_regex

cwd = os.path.dirname(os.path.realpath(__file__))

data={
    'val_int':{
        'rstring':r'val_int *=',
        'typ':int
        },
    'val_float':{
        'rstring':r'val_float *=',
        'typ':float

        },
    'val_str':{
        'rstring':r'val_str *=',
        'typ':str
        },
    'multiline':{
        'rstring':r'multiline_\d+\s*=\s*(?P<v1>\d+),\s*(?P<v2>\d+),\s*(?P<v3>\d+)\s*',
        'typ':np.ndarray
        },
    'multiline_split1,multiline_split2,multiline_split3':{
        'rstring':r'multiline_\d+\s*=\s*(?P<v1>\d+),\s*(?P<v2>\d+),\s*(?P<v3>\d+)\s*',
        'typ':np.ndarray
        },
    'multiline_num':{
        'rstring':r'multiline_\d+\s*=\s*(?P<v1>\d+),\s*(?P<v2>\d+),\s*(?P<v3>\d+)\s*',
        'typ':np.ndarray,
        'max_num':3
        },
    'multiline_scale':{
        'rstring':r'multiline_\d+\s*=\s*(?P<v1>\d+),\s*(?P<v2>\d+),\s*(?P<v3>\d+)\s*',
        'typ':np.ndarray,
        're_scale_fact':3.
        },

}

@pytest.fixture(scope='module')
def cls():
    out = os.path.join(cwd, 'test.out')
    new = Parser_regex(file=out, regex_data=data)

    return new

def test_val_int(cls):
    res = cls.val_int
    assert res == 125, "Failed to parse int"

def test_val_float(cls):
    res = cls.val_float
    assert res == 2.3, "Failed to parse float"

def test_val_str(cls):
    res = cls.val_str
    assert res == 'non saprei', "Failed to parse str"

def test_multiline(cls):
    res = cls.multiline

    assert len(res) == 6, "Failed to parse multiline num"
    assert all(all(b == n+1 for n,b in enumerate(a)) for a in res), "Failed to parse multiline vals"

def test_multiline_split(cls):
    res1 = cls.multiline_split1
    res2 = cls.multiline_split2
    res3 = cls.multiline_split3

    assert len(res1) == len(res2) & len(res1) == len(res3) & len(res1) == 6, "Failed to parse multiline num"
    assert all(np.all(a == 1.) for a in res1), "Failed to parse multiline split"
    assert all(np.all(a == 2.) for a in res2), "Failed to parse multiline split"
    assert all(np.all(a == 3.) for a in res3), "Failed to parse multiline split"

def test_multiline_num(cls):
    res  = cls.multiline_num
    assert len(res) == 3, "Failed to parse multiline max_num"

def test_multiline_scale(cls):
    res  = cls.multiline_scale
    assert len(res) == 6, "Failed to parse multiline max_num"
    assert all(all(b == (n+1)*3 for n,b in enumerate(a)) for a in res), "Failed to parse multiline scale"


