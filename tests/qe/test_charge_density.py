import os

from qepppy.qe.charge_density import charge_density

cwd      = os.path.dirname(os.path.realpath(__file__))
test_dir = os.path.join(cwd, 'test_files')

def test_charge_density_read():
    """Test that charge density can be read from file."""
    cd = charge_density(src=os.path.join(test_dir, 'test-charge-density.dat'))

    assert cd.igwx ==  27
    assert cd.ispin == 1
    assert cd.recipr.shape == (3,3)
    assert cd.gvect.shape == (cd.igwx, 3)
    assert cd.C_kn.shape == (cd.ispin, cd.igwx,)

def test_charge_density_write(tmp_path):
    """Test that charge density can be written to file."""
    cd = charge_density(src=os.path.join(test_dir, 'test-charge-density.dat'))

    cd.write_binary(str(tmp_path / 'test.dat'))

    with open(os.path.join(test_dir, 'test-charge-density.dat'), 'rb') as f:
        std = f.read()
    with open(os.path.join(tmp_path, 'test.dat'), 'rb') as f:
        new = f.read()

    assert std == new
