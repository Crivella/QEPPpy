import os
import pickle

import pytest

import qepppy

from .conftest import compare_std

cwd      = os.path.dirname(os.path.realpath(__file__))
test_dir = os.path.join(cwd, 'test_files')

in_files = [a[:-3] for a in os.listdir(test_dir) if a.endswith('.in')]

cmp_list = [
    'ibrav', 'alat',
    'n_atoms', 'atoms_coord_cart', 'atoms_coord_cryst', 'atoms_typ', 'unique_atoms_mass',
    ]

@pytest.fixture(
    scope='module',
    params=in_files
    )
def inputs(request):
    """
    Collection of inputs for pw.in.
    Returns:
     - inp: parsed input file
     - pkl: dictionary of expected results
    """
    name = request.param
    name = os.path.join(test_dir, name)

    in_name  = name + '.in'
    pkl_name = name + '_in.pickle'

    inp = qepppy.qe.pw_in(input_file=in_name)

    with open(pkl_name, 'rb') as f:
        pkl = pickle.load(f)

    return inp, pkl

def test_pw_in_parsing(inputs):
    inp, std = inputs

    compare_std(inp, std, cmp_list=cmp_list)
    inp.validate()

def test_pw_in_build_init(inputs):
    inp, std = inputs
    new = qepppy.qe.pw_in(
        input_data={
            'SYSTEM':{
                'nat':       inp.n_atoms,
                'ntyp':      inp.n_types,
                'ibrav':     inp.ibrav,
                'ecutwfc':   inp.ecutwfc,
                'celldm(1)': inp.alat,
                },
            'ELECTRONS':{}
            }
        )

    new.atoms_coord_cart  = inp.atoms_coord_cart
    new.atoms_typ         = inp.atoms_typ
    new.kpt_mesh          = inp.kpt_mesh
    new.kpt_cart          = inp.kpt_cart
    new.kpt_weight        = inp.kpt_weight
    new.unique_atoms_mass    = inp.unique_atoms_mass
    new.unique_atoms_pseudo  = inp.unique_atoms_pseudo

    new.validate()
    str(new)

def test_pw_in_build_scratch(inputs):
    inp, std = inputs
    new = qepppy.qe.pw_in()

    new.n_atoms           = inp.n_atoms
    new.n_types           = inp.n_types
    new.ibrav             = inp.ibrav
    new['SYSTEM/ecutwfc'] = inp.ecutwfc
    new.alat              = inp.alat
    new.atoms_coord_cart  = inp.atoms_coord_cart
    new.atoms_typ         = inp.atoms_typ
    new.kpt_mesh          = inp.kpt_mesh
    new.kpt_cart          = inp.kpt_cart
    new.kpt_weight        = inp.kpt_weight
    new.unique_atoms_mass    = inp.unique_atoms_mass
    new.unique_atoms_pseudo  = inp.unique_atoms_pseudo

    new.validate()
    str(new)

if __name__ == '__main__':
    """
    Generates the baseline files for the automated testing.
    To be run with a stable pre-tested version.
    """
    for name in in_files:
        in_name  = os.path.join(test_dir, name + '.in')
        pkl_name = os.path.join(test_dir, name + '_in.pickle')

        xml = qepppy.qe.pw_in(input_file=in_name)

        res = {}
        for k in cmp_list:
            res[k] = getattr(xml, k)

        with open(pkl_name, 'wb') as f:
            pickle.dump(res, f)
