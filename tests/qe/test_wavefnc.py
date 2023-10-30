import hashlib
import os

from qepppy.qe.wavefnc import wavefnc

cwd      = os.path.dirname(os.path.realpath(__file__))
test_dir = os.path.join(cwd, 'test_files')

file = os.path.join(test_dir, 'test-wfc1.dat')

def test_wavefnc_read():
    """Test that charge density can be read from file."""
    cd = wavefnc(src=file)

    assert cd.igwx ==  27
    assert cd.ispin == 1
    assert cd.nbnd == 2
    assert cd.recipr.shape == (3,3)
    assert cd.gvect.shape == (cd.igwx, 3)
    assert cd.C_kn.shape == (cd.nbnd, cd.ispin, cd.igwx,)

def test_wavefnc_write(tmp_path):
    """Test that charge density can be written to file."""
    cd = wavefnc(src=os.path.join(test_dir, 'test-wfc1.dat'))

    cd.write_binary(str(tmp_path / 'test.dat'))

    with open(file, 'rb') as f:
        std = f.read()
    with open(os.path.join(tmp_path, 'test.dat'), 'rb') as f:
        new = f.read()

    assert std == new
