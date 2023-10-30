import pytest


@pytest.fixture
def tmpfile(tmpdir):
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w+', dir=tmpdir) as f:
        yield f

@pytest.fixture(scope='function')
def mock_nested_calls():
    """Mock class with callable functions."""
    class MockNested():
        def __init__(self):
            self.name = ''
            self.called = {}
            self.called_args = {}
            self.called_kwargs = {}
        def __len__(self):
            return 0
        def __iter__(self):
            return iter([self, self])
        def __call__(self, *args, **kwargs):
            ptr = self.called_args.setdefault(self.name, [])
            ptr.append(args)
            ptr = self.called_kwargs.setdefault(self.name, [])
            ptr.append(kwargs)
            return self
        def __getattr__(self, name):
            self.name = name
            self.called.setdefault(name, 0)
            self.called[name] += 1
            return self

    return MockNested()
