import io
import os

import numpy as np
import numpy.testing as npt
import pandas as pd

from fluxpart.hfdata import HFData, HFDataSource
from fluxpart.fluxpart import _converter_func

TESTDIR = os.path.dirname(os.path.realpath(__file__))
DATADIR = os.path.join(TESTDIR, "data")


def test_hfdata_reader():
    toy_data = (
        "foobar baz\n"
        "asdf,0,2,3,4,5,6,7,9,0\n"
        "asdf,1,2,3,4,5,6,7,9,0\n"
        "asdf,2,2,3,4,5,6,7,9,1\n"
        "asdf,3,2,3,4,5,6,,9,0\n"
        "asdf,4,2,3,4,5,6,7,9,0\n"
        "# foo\n"
        "asdf,5,2,3,4,5,6,7,9,0\n"
        "asdf,6,2,3,4,5,6,7,xxx,0\n"
        "asdf,7,???,3,4,5,6,7,9,0\n"
        "asdf,8,2,3,4,5,6,7,9,0\n"
        "asdf,9,2,3,4,5,6,7,9,0\n"
        "asdf,10, 2,3,4,5,6,7,9,0\n"
        "asdf,11,-2,3,4,5,6,7,9,0\n"
    )

    source = HFDataSource(
        filetype="csv",
        cols=(1, 2, 3, 6, 7, 4, 5),
        comment="#",
        skiprows=1,
        na_values="???",
        converters={"q": _converter_func(10., 0)},
        flags=(9, 0),
        delimiter=",",
    )

    toy = source._readfile(io.BytesIO(toy_data.encode()), df_output=False)
    toy.cleanse(rd_tol=0.1, ad_tol=2, bounds={"v": (0, np.inf)})

    npt.assert_allclose(toy.dataframe["u"], [4, 5, 6])
    npt.assert_allclose(toy.dataframe["v"], 3 * [2])
    npt.assert_allclose(toy.dataframe["w"], 3 * [3])
    npt.assert_allclose(toy.dataframe["q"], 3 * [70])
    npt.assert_allclose(toy.dataframe["c"], 3 * [6])
    npt.assert_allclose(toy.dataframe["T"], 3 * [4])
    npt.assert_allclose(toy.dataframe["P"], 3 * [5])

    # toa5
    cols = (2, 3, 4, 5, 6, 7, 8)
    fname = os.path.join(DATADIR, "TOA5_6843.ts_Above_2012_06_07_1300.dat")

    kws = dict(
        skiprows=4,
        converters={
            "T": _converter_func(1, 273.15),
            "q": _converter_func(1e-3, 0),
            "c": _converter_func(1e-6, 0),
            "P": _converter_func(1e3, 0),
        },
        flags=[(9, 0)],
    )

    source = HFDataSource(filetype="csv", cols=cols, **kws)
    data = HFData(source._readfile(fname, df_output=True))
    assert_1300_read(data)

    kws = dict(
        skiprows=4,
        time_col=0,
        converters={
            "T": _converter_func(1, 273.15),
            "q": _converter_func(1e-3, 0),
            "c": _converter_func(1e-6, 0),
            "P": _converter_func(1e3, 0),
        },
        flags=(9, 0),
    )

    source = HFDataSource(filetype="csv", cols=cols, **kws)
    reader = source.reader(fname, interval="15min")
    time, data = next(reader)
    data = HFData(data)
    assert data.dataframe.index[0] == pd.to_datetime("2012-06-07 13:00:00.05")
    assert data.dataframe.index[-1] == pd.to_datetime("2012-06-07 13:14:59.95")
    assert source._peektime(fname) == pd.to_datetime("2012-06-07 13:00:00.05")
    assert_1300_interval_read(data)

    reader = source.reader(DATADIR, interval="10min", ext=".dat")
    for cnt, (time, df) in enumerate(reader):
        assert_10min_chunk_read(cnt, time, df)
    reader = source.reader(DATADIR, interval="15min", ext=".dat")
    for cnt, (time, df) in enumerate(reader):
        assert_15min_chunk_read(cnt, time, df)
    reader = source.reader(DATADIR, interval="20min", ext=".dat")
    for cnt, (time, df) in enumerate(reader):
        assert_20min_chunk_read(cnt, time, df)

    fnames = [
        os.path.join(DATADIR, "TOA5_6843.ts_Above_2012_06_07_1300.dat"),
        os.path.join(DATADIR, "TOA5_6843.ts_Above_2012_06_07_1245.dat"),
    ]
    reader = source.reader(fnames, interval="20min", ext=".dat")
    for cnt, (time, df) in enumerate(reader):
        assert_20min_chunk_read(cnt, time, df)

    # pd.df
    df = pd.read_csv(
        fname,
        usecols=[2, 3, 4, 5, 6, 7, 8, 9],
        names=["u", "v", "w", "c", "q", "T", "P", "ex_flag"],
        skiprows=4,
    )

    source = HFDataSource(
        filetype="pd.df",
        flags=(7, 0),
        cols=[0, 1, 2, 3, 4, 5, 6],
        converters={
            "T": _converter_func(1, 273.15),
            "q": _converter_func(1e-3, 0),
            "c": _converter_func(1e-6, 0),
            "P": _converter_func(1e3, 0),
        },
    )
    data = source._readfile(df)
    assert_1300_read(data)

    # tob
    fname = os.path.join(DATADIR, "testing.tob")
    kws = dict(
        filetype="tob1",
        cols=(3, 4, 5, 6, 7, 8, 9),
        converters={
            "T": _converter_func(1, 273.15),
            "q": _converter_func(1e-3, 0),
            "c": _converter_func(1e-6, 0),
            "P": _converter_func(1e3, 0),
        },
    )

    source = HFDataSource(**kws)
    time, data = next(source.reader(fname, interval="30min"))
    data = HFData(data)
    assert_tob_read(data)
    assert source._peektime(fname) == pd.to_datetime("2017-08-03 00:00:00.1")
    assert data.dataframe.index[0] == pd.to_datetime("2017-08-03 00:00:00.1")
    assert data.dataframe.index[-1] == pd.to_datetime("2017-08-03 00:00:14.4")

    reader = source.reader(DATADIR, interval="5s", ext=".tob")
    for cnt, (time, df) in enumerate(reader):
        assert_5S_tobchunk_read(cnt, time, df)
    reader = source.reader(DATADIR, interval="1min", ext=".tob")
    for cnt, (time, df) in enumerate(reader):
        assert_1min_tobchunk_read(cnt, time, df)


def assert_1300_read(data):
    npt.assert_allclose(data["u"].iloc[0], 0.468)
    npt.assert_allclose(data["v"].iloc[0], -0.9077501)
    npt.assert_allclose(data["w"].iloc[0], 0.1785)
    npt.assert_allclose(data["c"].iloc[0], 659.7584e-6)
    npt.assert_allclose(data["q"].iloc[0], 9.530561e-3)
    npt.assert_allclose(data["T"].iloc[0], 28.52527 + 273.15)
    npt.assert_allclose(data["P"].iloc[0], 100.1938e3)
    npt.assert_allclose(data["u"].iloc[-1], 1.3675)
    npt.assert_allclose(data["v"].iloc[-1], -0.75475)
    npt.assert_allclose(data["w"].iloc[-1], -0.1775)
    npt.assert_allclose(data["c"].iloc[-1], 658.2624e-6)
    npt.assert_allclose(data["q"].iloc[-1], 9.404386e-3)
    npt.assert_allclose(data["T"].iloc[-1], 28.35199 + 273.15)
    npt.assert_allclose(data["P"].iloc[-1], 100.1938e3)

    npt.assert_allclose(data["u"].mean(), 1.43621, atol=1e-4)
    npt.assert_allclose(data["v"].mean(), -0.634818, atol=1e-4)
    npt.assert_allclose(data["w"].mean(), 0.0619483, atol=1e-4)
    npt.assert_allclose(data["c"].mean(), 659.052e-6, atol=1e-9)
    npt.assert_allclose(data["q"].mean(), 9.56732e-3, atol=1e-7)
    npt.assert_allclose(data["T"].mean(), 28.5431 + 273.15, atol=1e-4)
    npt.assert_allclose(data["P"].mean(), 100.179e3, atol=1e0)


def assert_1300_interval_read(data):
    """last line in file is not read because belongs to next interval"""
    npt.assert_allclose(data["u"].iloc[0], 0.468)
    npt.assert_allclose(data["v"].iloc[0], -0.9077501)
    npt.assert_allclose(data["w"].iloc[0], 0.1785)
    npt.assert_allclose(data["c"].iloc[0], 659.7584e-6)
    npt.assert_allclose(data["q"].iloc[0], 9.530561e-3)
    npt.assert_allclose(data["T"].iloc[0], 28.52527 + 273.15)
    npt.assert_allclose(data["P"].iloc[0], 100.1938e3)
    npt.assert_allclose(data["u"].iloc[-1], 1.20050)
    npt.assert_allclose(data["v"].iloc[-1], -0.775)
    npt.assert_allclose(data["w"].iloc[-1], -0.1610)
    npt.assert_allclose(data["c"].iloc[-1], 658.3277e-6)
    npt.assert_allclose(data["q"].iloc[-1], 9.394239e-3)
    npt.assert_allclose(data["T"].iloc[-1], 28.35199 + 273.15)
    npt.assert_allclose(data["P"].iloc[-1], 100.1678e3)

    npt.assert_allclose(data["u"].mean(), 1.43621, atol=1e-4)
    npt.assert_allclose(data["v"].mean(), -0.634818, atol=1e-4)
    npt.assert_allclose(data["w"].mean(), 0.0619483, atol=1e-4)
    npt.assert_allclose(data["c"].mean(), 659.052e-6, atol=1e-9)
    npt.assert_allclose(data["q"].mean(), 9.56732e-3, atol=1e-7)
    npt.assert_allclose(data["T"].mean(), 28.5431 + 273.15, atol=1e-4)
    npt.assert_allclose(data["P"].mean(), 100.179e3, atol=1e0)


def assert_tob_read(data):
    npt.assert_allclose(data["u"].iloc[0], -2.57175016)
    npt.assert_allclose(data["v"].iloc[0], 1.6450001)
    npt.assert_allclose(data["w"].iloc[0], -0.12725)
    npt.assert_allclose(data["c"].iloc[0], 612.54150391e-6)
    npt.assert_allclose(data["q"].iloc[0], 13.11471748e-3)
    npt.assert_allclose(data["T"].iloc[0], 23.29580116 + 273.15)
    npt.assert_allclose(data["P"].iloc[0], 85.04070282e3)
    npt.assert_allclose(data["u"].iloc[-1], -2.4402502)
    npt.assert_allclose(data["v"].iloc[-1], 1.5402501)
    npt.assert_allclose(data["w"].iloc[-1], -0.11375)
    npt.assert_allclose(data["c"].iloc[-1], 615.627e-6)
    npt.assert_allclose(data["q"].iloc[-1], 13.200139e-3)
    npt.assert_allclose(data["T"].iloc[-1], 23.015879 + 273.15)
    npt.assert_allclose(data["P"].iloc[-1], 85.0407e3)
    assert data.dataframe.index[0] == pd.to_datetime("2017-08-03 00:00:00.1")
    assert data.dataframe.index[-1] == pd.to_datetime("2017-08-03 00:00:14.4")


def assert_10min_chunk_read(ichunk, time, df):
    times = [
        "2012-06-07 12:40:00",
        "2012-06-07 12:50:00",
        "2012-06-07 13:00:00",
        "2012-06-07 13:10:00",
    ]
    starts = [
        "2012-06-07 12:45:00.05",
        "2012-06-07 12:50:00",
        "2012-06-07 13:00:00",
        "2012-06-07 13:10:00",
    ]
    stops = [
        "2012-06-07 12:49:59.95",
        "2012-06-07 12:59:59.95",
        "2012-06-07 13:09:59.95",
        "2012-06-07 13:15:00",
    ]
    assert time == pd.to_datetime(times[ichunk])
    assert df.index[0] == pd.to_datetime(starts[ichunk])
    assert df.index[-1] == pd.to_datetime(stops[ichunk])


def assert_15min_chunk_read(ichunk, time, df):
    times = [
        "2012-06-07 12:45:00",
        "2012-06-07 13:00:00",
        "2012-06-07 13:15:00",
    ]
    starts = [
        "2012-06-07 12:45:00.05",
        "2012-06-07 13:00:00",
        "2012-06-07 13:15:00",
    ]
    stops = [
        "2012-06-07 12:59:59.95",
        "2012-06-07 13:14:59.95",
        "2012-06-07 13:15:00",
    ]
    assert time == pd.to_datetime(times[ichunk])
    assert df.index[0] == pd.to_datetime(starts[ichunk])
    assert df.index[-1] == pd.to_datetime(stops[ichunk])


def assert_20min_chunk_read(ichunk, time, df):
    times = ["2012-06-07 12:40:00", "2012-06-07 13:00:00"]
    starts = ["2012-06-07 12:45:00.05", "2012-06-07 13:00:00"]
    stops = ["2012-06-07 12:59:59.95", "2012-06-07 13:15:00"]
    assert time == pd.to_datetime(times[ichunk])
    assert df.index[0] == pd.to_datetime(starts[ichunk])
    assert df.index[-1] == pd.to_datetime(stops[ichunk])


def assert_5S_tobchunk_read(cnt, time, df):
    times = [
        "2017-08-03 00:00:00",
        "2017-08-03 00:00:05",
        "2017-08-03 00:00:10",
    ]
    starts = [
        "2017-08-03 00:00:00.1",
        "2017-08-03 00:00:05",
        "2017-08-03 00:00:10",
    ]
    stops = [
        "2017-08-03 00:00:04.9",
        "2017-08-03 00:00:09.9",
        "2017-08-03 00:00:14.4",
    ]
    assert time == pd.to_datetime(times[cnt])
    assert df.index[0] == pd.to_datetime(starts[cnt])
    assert df.index[-1] == pd.to_datetime(stops[cnt])


def assert_1min_tobchunk_read(cnt, time, df):
    times = ["2017-08-03 00:00:00"]
    starts = ["2017-08-03 00:00:00.1"]
    stops = ["2017-08-03 00:00:14.4"]
    assert time == pd.to_datetime(times[cnt])
    assert df.index[0] == pd.to_datetime(starts[cnt])
    assert df.index[-1] == pd.to_datetime(stops[cnt])


if __name__ == "__main__":
    test_hfdata_reader()
