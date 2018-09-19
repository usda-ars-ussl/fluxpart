import io
from unittest import TestCase

import numpy as np
import numpy.testing as npt
import pandas as pd
from pandas.util.testing import assert_frame_equal

from fluxpart.util import (
    stats2,
    multifile_read_csv,
    chunked_df,
    HFDataReadWarning,
)


def test_stats2():
    """Test stats2 func from fluxpart.util"""

    data = "7 8 4\n6 1 3\n10 6 6\n6 7 3\n8 2 4"
    dtype = [("v0", int), ("v1", int), ("v2", int)]
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
    dtype = [("v0", int), ("v1", int), ("v2", int)]
    arr = np.genfromtxt(io.BytesIO(data.encode()), dtype=dtype)
    ans = stats2(arr, names=("v0", "v2"))

    npt.assert_allclose(ans.ave_v0, 37 / 5)
    npt.assert_allclose(ans.ave_v2, 4)
    npt.assert_allclose(ans.var_v0, 14 / 5)
    npt.assert_allclose(ans.var_v2, 3 / 2)
    npt.assert_allclose(ans.cov_v0_v2, 2)
    npt.assert_allclose(ans.cov_v2_v0, ans.cov_v0_v2)

    assert not hasattr(ans, "ave_v1")
    assert not hasattr(ans, "var_v1")
    assert not hasattr(ans, "cov_v0_v1")
    assert not hasattr(ans, "cov_v1_v0")
    assert not hasattr(ans, "cov_v1_v2")
    assert not hasattr(ans, "cov_v2_v1")


def test_mulitifile_read_csv():
    file1 = io.BytesIO("1,2,3\n4,5,6\n7,8,9\n10,11,12".encode())
    file2 = io.BytesIO("13,14,15\n16,17,18\n19,20,21\n22,23,24".encode())
    reader = multifile_read_csv([file1, file2], header=None)
    dfs = [
        pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]),
        pd.DataFrame([[13, 14, 15], [16, 17, 18], [19, 20, 21], [22, 23, 24]]),
    ]
    for cnt, df in enumerate(reader):
        assert_frame_equal(df, dfs[cnt])

    file1 = io.BytesIO("1,2,3\n4,5,6\n7,8,9\n10,11,12".encode())
    file2 = io.BytesIO("13,14,15\n16,17,18\n19,20,21\n22,23,24".encode())
    reader = multifile_read_csv([file1, file2], header=None, chunksize=2)
    dfs = [
        pd.DataFrame([[1, 2, 3], [4, 5, 6]]),
        pd.DataFrame([[7, 8, 9], [10, 11, 12]], index=[2, 3]),
        pd.DataFrame([[13, 14, 15], [16, 17, 18]]),
        pd.DataFrame([[19, 20, 21], [22, 23, 24]], index=[2, 3]),
    ]
    for cnt, df in enumerate(reader):
        assert_frame_equal(df, dfs[cnt])

    file1 = io.BytesIO("1,2,3\n4,5,6\n7,8,9\n10,11,12".encode())
    file2 = io.BytesIO("13,14,15\n16,17,18\n19,20,21\n22,23,24".encode())
    reader = multifile_read_csv([file1, file2], header=None, chunksize=3)
    dfs = [
        pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
        pd.DataFrame([[10, 11, 12]], index=[3]),
        pd.DataFrame([[13, 14, 15], [16, 17, 18], [19, 20, 21]]),
        pd.DataFrame([[22, 23, 24]], index=[3]),
    ]
    for cnt, df in enumerate(reader):
        assert_frame_equal(df, dfs[cnt])


def test_chunked_df():
    df = pd.DataFrame(np.random.rand(50, 3))
    df = df.set_index(pd.date_range("2000-01-01", periods=50, freq="1h"))
    df1, df2 = df.iloc[:40], df.iloc[40:]
    cdf = chunked_df(iter([df1, df2]), time_interval="1D")
    df = next(cdf)
    assert df.index[0] == pd.to_datetime("2000-01-01 00:00:00")
    assert df.index[-1] == pd.to_datetime("2000-01-01 23:00:00")
    df = next(cdf)
    assert df.index[0] == pd.to_datetime("2000-01-02 00:00:00")
    assert df.index[-1] == pd.to_datetime("2000-01-02 23:00:00")
    df = next(cdf)
    assert df.index[0] == pd.to_datetime("2000-01-03 00:00:00")
    assert df.index[-1] == pd.to_datetime("2000-01-03 01:00:00")


def test_read_warning():
    # not working
    # see: https://bugs.python.org/issue29620
    tc = TestCase()
    with tc.assertWarns(HFDataReadWarning):
        hfread_error()


def hfread_error():
    file1 = io.BytesIO("1,2,3\n4,5,6\n7,8,9\n10,11,12".encode())
    file2 = io.BytesIO("13 14 15\n16,17,18\n19,20,21\n22,23,24".encode())
    reader = multifile_read_csv([file1, file2], header=None)
    for cnt, df in enumerate(reader):
        pass


if __name__ == "__main__":
    test_stats2()
    test_mulitifile_read_csv()
    test_chunked_df()
    # test_read_warning()
