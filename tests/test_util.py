import io
import numpy as np
import numpy.testing as npt
from fluxpart.util import stats2


def test_stats2():
    """Test stats2 func from fluxpart.util"""

    data = "7 8 4\n6 1 3\n10 6 6\n6 7 3\n8 2 4"
    dtype = [('v0', int), ('v1', int), ('v2', int)]
    arr = np.genfromtxt(io.BytesIO(data.encode()), dtype=dtype)
    ans = stats2(arr)

    npt.assert_allclose(ans.ave_v0, 37 / 5)
    npt.assert_allclose(ans.ave_v1, 24 / 5)
    npt.assert_allclose(ans.ave_v2, 4)
    npt.assert_allclose(ans.var_v0, 14 / 5)
    npt.assert_allclose(ans.var_v1, 97 / 10)
    npt.assert_allclose(ans.var_v2, 3 / 2)
    npt.assert_allclose(ans.cov_v0_v1, 3 / 5)
    npt.assert_allclose(ans.cov_v0_v2, 2)
    npt.assert_allclose(ans.cov_v1_v0, ans.cov_v0_v1)
    npt.assert_allclose(ans.cov_v1_v2, 1)
    npt.assert_allclose(ans.cov_v2_v0, ans.cov_v0_v2)
    npt.assert_allclose(ans.cov_v2_v1, ans.cov_v1_v2)

    data = "7 8 4\n6 1 3\n10 6 6\n6 7 3\n8 2 4"
    dtype = [('v0', int), ('v1', int), ('v2', int)]
    arr = np.genfromtxt(io.BytesIO(data.encode()), dtype=dtype)
    ans = stats2(arr, names=('v0', 'v2'))

    npt.assert_allclose(ans.ave_v0, 37 / 5)
    npt.assert_allclose(ans.ave_v2, 4)
    npt.assert_allclose(ans.var_v0, 14 / 5)
    npt.assert_allclose(ans.var_v2, 3 / 2)
    npt.assert_allclose(ans.cov_v0_v2, 2)
    npt.assert_allclose(ans.cov_v2_v0, ans.cov_v0_v2)

    assert not hasattr(ans, 'ave_v1')
    assert not hasattr(ans, 'var_v1')
    assert not hasattr(ans, 'cov_v0_v1')
    assert not hasattr(ans, 'cov_v1_v0')
    assert not hasattr(ans, 'cov_v1_v2')
    assert not hasattr(ans, 'cov_v2_v1')

if __name__ == "__main__":
    test_stats2()

