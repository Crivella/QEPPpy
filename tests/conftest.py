import pytest

@pytest.fixture
def tmpfile(tmpdir):
	import tempfile
	with tempfile.NamedTemporaryFile(mode='w+', dir=tmpdir) as f:
		yield f