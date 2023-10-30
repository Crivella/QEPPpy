import sys

import numpy as np
import pytest

from qepppy import _decorators as dec


def test_numpy_save_opt_no_fname():
    """Test that numpy_save_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_save_opt()
    def test_func():
        return np.array([1,2,3])

    with pytest.raises(ValueError):
        test_func()

def test_numpy_save_opt_pfile_false(tmp_path):
    """Test that no output is saved if pFile is set to False."""
    fname = tmp_path / 'test.dat'
    @dec.numpy_save_opt(_fname=fname)
    def test_func():
        return np.array([1,2,3])

    test_func(pFile=False)

    assert not fname.exists()

def test_numpy_save_opt_default_fname(tmp_path):
    """Test that numpy_save_opt decorator saves to default filename."""
    fname = str(tmp_path / 'test.dat')
    @dec.numpy_save_opt(_fname=fname)
    def test_func():
        return np.array([1,2,3])

    assert np.all(test_func() == [1,2,3])
    assert np.all(np.loadtxt(fname) == [1,2,3])

def test_numpy_save_opt_func_fname(tmp_path):
    """Test that numpy_save_opt decorator saves to function specified filename."""
    fname = str(tmp_path / 'test.dat')
    @dec.numpy_save_opt()
    def test_func():
        return np.array([1,2,3])

    assert np.all(test_func(fname=fname) == [1,2,3])
    assert np.all(np.loadtxt(fname) == [1,2,3])

def test_numpy_save_opt_format_int(tmp_path):
    """Test that numpy_save_opt decorator saves to function specified filename."""
    fname = str(tmp_path / 'test.dat')
    @dec.numpy_save_opt()
    def test_func():
        return np.array([[1,2,3]])

    assert np.all(test_func(fname=fname, fmt='%d') == [[1,2,3]])
    assert np.all(np.loadtxt(fname) == [[1,2,3]])

    with open(fname, 'r') as f:
        content = f.read()
    content = '\n'.join(filter(lambda x: not x.startswith('#') ,content.split('\n')))
    assert content.strip() == '1 2 3'

def test_numpy_save_opt_format_float(tmp_path):
    """Test that numpy_save_opt decorator saves to function specified filename."""
    fname = str(tmp_path / 'test.dat')
    @dec.numpy_save_opt()
    def test_func():
        return np.array([[1,2,3]])

    assert np.all(test_func(fname=fname, fmt='%.2f') == [[1,2,3]])
    assert np.all(np.loadtxt(fname) == [[1,2,3]])

    with open(fname, 'r') as f:
        content = f.read()
    content = '\n'.join(filter(lambda x: not x.startswith('#') ,content.split('\n')))
    assert content.strip() == '1.00 2.00 3.00'

def test_numpy_save_opt_header(tmp_path):
    """Test that numpy_save_opt decorator saves to function specified filename."""
    fname = str(tmp_path / 'test.dat')
    @dec.numpy_save_opt()
    def test_func():
        return np.array([[1,2,3]])

    assert np.all(test_func(fname=fname, header='random text\ntext random') == [[1,2,3]])
    assert np.all(np.loadtxt(fname) == [[1,2,3]])

    with open(fname, 'r') as f:
        content = f.read()
    assert '# random text' in content
    assert '# text random' in content

def test_numpy_save_opt_delimiter(tmp_path):
    """Test that numpy_save_opt decorator saves to function specified filename."""
    fname = str(tmp_path / 'test.dat')
    @dec.numpy_save_opt()
    def test_func():
        return np.array([[1,2,3]])

    assert np.all(test_func(fname=fname, delimiter=' -- ', fmt='%d') == [[1,2,3]])

    with open(fname, 'r') as f:
        content = f.read()
    content = '\n'.join(filter(lambda x: not x.startswith('#') ,content.split('\n')))
    assert content.strip() == '1 -- 2 -- 3'

def test_numpy_plot_plot_false(mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2,3], [2,4,9]])

    test_func(plot=False, ax=mock_nested_calls)
    assert not mock_nested_calls.called

def test_numpy_plot_plot_with_ax(mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2,3], [2,4,9]])

    test_func(ax=mock_nested_calls)
    assert mock_nested_calls.called['plot'] == 2
    assert mock_nested_calls.called['set_xlabel'] == 1
    assert mock_nested_calls.called['set_ylabel'] == 1
    assert mock_nested_calls.called['set_minor_locator'] == 1
    assert mock_nested_calls.called['set_tick_params'] == 1

def test_numpy_plot_plot_with_ax_end(mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2,3], [2,4,9]])

    test_func(ax=mock_nested_calls, end=1)
    assert mock_nested_calls.called['plot'] == 1

def test_numpy_plot_plot_with_ax_xlim(mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2,3], [2,4,9]])

    test_func(ax=mock_nested_calls, xlim=1)
    assert mock_nested_calls.called['set_xlim'] == 1

def test_numpy_plot_plot_with_ax_ylim(mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2,3], [2,4,9]])

    test_func(ax=mock_nested_calls, ylim=1)
    assert mock_nested_calls.called['set_ylim'] == 1

def test_numpy_plot_plot_with_ax_Nx2array(mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2], [2,4]])

    test_func(ax=mock_nested_calls)
    assert mock_nested_calls.called['plot'] == 1

def test_numpy_plot_plot_with_ax_no_dash_list(mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2,3], [2,4,9]])

    test_func(ax=mock_nested_calls)
    for i, kwargs in enumerate(mock_nested_calls.called_kwargs['plot']):
        assert isinstance(kwargs['dashes'], tuple)
        assert kwargs['dashes'][0] == 8

def test_numpy_plot_plot_with_ax_dash_list(mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2,3], [2,4,9]])

    dash_list = [1,2,3]
    test_func(ax=mock_nested_calls, dash_list=dash_list)
    for i, kwargs in enumerate(mock_nested_calls.called_kwargs['plot']):
        assert kwargs['dashes'] == dash_list[i % len(dash_list)]

def test_numpy_plot_plot_with_no_ax(monkeypatch, mock_nested_calls):
    """Test that numpy_plot_opt decorator raises ValueError if no fname is set."""
    @dec.numpy_plot_opt()
    def test_func():
        return np.array([[1,2], [2,4]])

    monkeypatch.setattr(dec, 'plt', mock_nested_calls)

    test_func()
    assert mock_nested_calls.called['subplots'] == 1
    assert mock_nested_calls.called['legend'] == 1
    assert mock_nested_calls.called['show'] == 1

def test_set_self_single():
    """Test that set_self decorator works with single argument."""
    class Mock():
        @dec.set_self('a')
        def mock(self):
            return 1

    m = Mock()
    assert m.mock() == 1
    assert m.a == 1

def test_set_self_multiple():
    """Test that set_self decorator works with single argument."""
    class Mock():
        @dec.set_self('a,b')
        def mock(self):
            return 1, 2

    m = Mock()
    assert m.mock() == (1, 2)
    assert m.a == 1
    assert m.b == 2

def test_set_self_disable_default():
    """Test that set_self decorator works with single argument."""
    class Mock():
        @dec.set_self('a', _default=False)
        def mock(self):
            return 1

    m = Mock()
    assert m.mock() == 1
    assert not hasattr(m, 'a')

def test_set_self_disable_arg():
    """Test that set_self decorator works with single argument."""
    class Mock():
        @dec.set_self('a')
        def mock(self):
            return 1

    m = Mock()
    assert m.mock(do_set_self=False) == 1
    assert not hasattr(m, 'a')

def test_store_property():
    """Test that store_property decorator works."""
    class Mock():
        def __init__(self):
            self.once = False

        @property
        @dec.store_property
        def mock(self):
            if not self.once:
                self.once = True
                return 1
            assert False

    m = Mock()
    assert m.mock == 1
    assert m.mock == 1

def test_IOStdoutRedirect_fname_default(tmp_path):
    """Test that IOStdoutRedirect decorator works with file name from default."""
    fname = str(tmp_path / 'test.txt')
    @dec.IOStdoutRedirect(_outfile=fname)
    def test_func():
        print('test')

    test_func()

    with open(fname, 'r') as f:
        content = f.read()
    assert content.strip() == 'test'

def test_IOStdoutRedirect_fname_arg(tmp_path):
    """Test that IOStdoutRedirect decorator works with file name from args."""
    fname = str(tmp_path / 'test.txt')
    @dec.IOStdoutRedirect()
    def test_func():
        print('test')

    test_func(outfile=fname)

    with open(fname, 'r') as f:
        content = f.read()
    assert content.strip() == 'test'

def test_IOStdoutRedirect_fobj(tmp_path):
    """Test that IOStdoutRedirect decorator works with file object."""
    fname = str(tmp_path / 'test.txt')
    @dec.IOStdoutRedirect()
    def test_func():
        print('test')

    with open(fname, 'w') as f:
        test_func(outfile=f)

    with open(fname, 'r') as f:
        content = f.read()
    assert content.strip() == 'test'

def test_IOStdoutRedirect_fobj_should_not_close(tmp_path):
    """Test that IOStdoutRedirect decorator with file pointer should not close it afterward."""
    fname = str(tmp_path / 'test.txt')
    @dec.IOStdoutRedirect()
    def test_func():
        print('test')

    with open(fname, 'w') as f:
        test_func(outfile=f)
        f.write('123')

    with open(fname, 'r') as f:
        content = f.read()
    assert content.strip() == 'test\n123'

def test_IOStdoutRedirect_no_file():
    """Test that IOStdoutRedirect decorator with no files should let stdou through."""
    from contextlib import redirect_stdout
    from io import StringIO
    @dec.IOStdoutRedirect()
    def test_func():
        print('test')

    f = StringIO()
    with redirect_stdout(f):
        test_func()

    assert f.getvalue().strip() == 'test'

################################################################################

def test_IOStderrRedirect_fname_default(tmp_path):
    """Test that IOStderrRedirect decorator works with file name from default."""
    fname = str(tmp_path / 'test.txt')
    @dec.IOStderrRedirect(_errfile=fname)
    def test_func():
        print('test', file=sys.stderr)

    test_func()

    with open(fname, 'r') as f:
        content = f.read()
    assert content.strip() == 'test'

def test_IOStderrRedirect_fname_arg(tmp_path):
    """Test that IOStderrRedirect decorator works with file name from args."""
    fname = str(tmp_path / 'test.txt')
    @dec.IOStderrRedirect()
    def test_func():
        print('test', file=sys.stderr)

    test_func(errfile=fname)

    with open(fname, 'r') as f:
        content = f.read()
    assert content.strip() == 'test'

def test_IOStderrRedirect_fobj(tmp_path):
    """Test that IOStderrRedirect decorator works with file object."""
    fname = str(tmp_path / 'test.txt')
    @dec.IOStderrRedirect()
    def test_func():
        print('test', file=sys.stderr)

    with open(fname, 'w') as f:
        test_func(errfile=f)

    with open(fname, 'r') as f:
        content = f.read()
    assert content.strip() == 'test'

def test_IOStderrRedirect_fobj_should_not_close(tmp_path):
    """Test that IOStderrRedirect decorator with file pointer should not close it afterward."""
    fname = str(tmp_path / 'test.txt')
    @dec.IOStderrRedirect()
    def test_func():
        print('test', file=sys.stderr)

    with open(fname, 'w') as f:
        test_func(errfile=f)
        f.write('123')

    with open(fname, 'r') as f:
        content = f.read()
    assert content.strip() == 'test\n123'

def test_IOStderrRedirect_no_file():
    """Test that IOStderrRedirect decorator with no files should let stderr through."""
    from contextlib import redirect_stderr
    from io import StringIO
    @dec.IOStderrRedirect()
    def test_func():
        print('test', file=sys.stderr)

    f = StringIO()
    with redirect_stderr(f):
        test_func()

    assert f.getvalue().strip() == 'test'
