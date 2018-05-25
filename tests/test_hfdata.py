import io
import os

import numpy as np
import numpy.testing as npt
import pandas as pd

from fluxpart.hfdata import HFData
from fluxpart.fluxpart import _converter_func, _str_converter_func

TESTDIR = os.path.dirname(os.path.realpath(__file__))
DATADIR = os.path.join(TESTDIR, 'data')


def test_hfdata_read_csv():
    cols = (2, 3, 4, 6, 5, 7, 8)
    fname = os.path.join(DATADIR, 'TOA5_6843.ts_Above_2012_06_07_1300.dat')

    kws = dict(
            converters={
                'T': _str_converter_func(1, 273.15),
                'q': _str_converter_func(1e-3, 0),
                'c': _str_converter_func(1e-6, 0),
                'P': _str_converter_func(1e3, 0)},
            delimiter=",",
            skiprows=4)

    data = HFData(cols=cols, flags={'ex_flag': (9, 0)}, **kws)
    data.read(fname)
    assert_1300_read(data)

    df = pd.read_csv(
            fname,
            usecols=[2, 3, 4, 5, 6, 7, 8, 9],
            names=['u', 'v', 'w', 'c', 'q', 'T', 'P', 'ex_flag'],
            skiprows=4)

    data = HFData(
            datasource='pd.df',
            flags={'ex_flag': (7, 0)},
            cols=[0, 1, 2, 4, 3, 5, 6],
            converters={
                'T': _converter_func(1, 273.15),
                'q': _converter_func(1e-3, 0),
                'c': _converter_func(1e-6, 0),
                'P': _converter_func(1e3, 0),
            },
           )
    data.read(df)
    assert_1300_read(data)

    fname = os.path.join(DATADIR, '2017_08_03_0000_ts_data.tob')
    kws = dict(
              datasource='tob1',
              cols=(3, 4, 5, 7, 6, 8, 9),
              converters={
                  'T': _converter_func(1, 273.15),
                  'q': _converter_func(1e-3, 0),
                  'c': _converter_func(1e-6, 0),
                  'P': _converter_func(1e3, 0),
              }
          )

    data = HFData(**kws)
    data._read_tob1(fname)
    assert_tob_read(data)

    toy_data = (
        'foobar baz\n'
        'asdf,0,2,3,4,5,6,7,9,0\n'
        'asdf,1,2,3,4,5,6,7,9,0\n'
        'asdf,2,2,3,4,5,6,7,9,1\n'
        'asdf,3,2,3,4,5,6,,9,0\n'
        'asdf,4,2,3,4,5,6,7,9,0\n'
        '# foo\n'
        'asdf,5,2,3,4,5,6,7,9,0\n'
        'asdf,6,2,3,4,5,6,7,xxx,0\n'
        'asdf,7,???,3,4,5,6,7,9,0\n'
        'asdf,8,2,3,4,5,6,7,9,0\n'
        'asdf,9,2,3,4,5,6,7,9,0\n'
        'asdf,10, 2,3,4,5,6,7,9,0\n'
        'asdf,11,-2,3,4,5,6,7,9,0\n'
    )

    toy = HFData(
        cols=(1, 2, 3, 7, 6, 4, 5),
        comment='#',
        skiprows=1,
        na_values="???",
        converters={'q': _str_converter_func(10., 0)},
        flags={'ex_flag': (9, 0)},
        delimiter=",",
        )
    toy.read(io.BytesIO(toy_data.encode()))
    toy.quality_check(
        rd_tol=0.1,
        ad_tol=2,
        bounds={'v': (0, np.inf)},
        )

    npt.assert_allclose(toy.dataframe['u'], [4, 5, 6])
    npt.assert_allclose(toy.dataframe['v'], 3 * [2, ])
    npt.assert_allclose(toy.dataframe['w'], 3 * [3, ])
    npt.assert_allclose(toy.dataframe['q'], 3 * [70, ])
    npt.assert_allclose(toy.dataframe['c'], 3 * [6, ])
    npt.assert_allclose(toy.dataframe['T'], 3 * [4, ])
    npt.assert_allclose(toy.dataframe['P'], 3 * [5, ])


def assert_1300_read(data):
    npt.assert_allclose(data['u'].iloc[0], 0.468)
    npt.assert_allclose(data['v'].iloc[0], -0.9077501)
    npt.assert_allclose(data['w'].iloc[0], 0.1785)
    npt.assert_allclose(data['c'].iloc[0], 659.7584e-6)
    npt.assert_allclose(data['q'].iloc[0], 9.530561e-3)
    npt.assert_allclose(data['T'].iloc[0], 28.52527 + 273.15)
    npt.assert_allclose(data['P'].iloc[0], 100.1938e3)
    npt.assert_allclose(data['u'].iloc[-1], 1.3675)
    npt.assert_allclose(data['v'].iloc[-1], -0.75475)
    npt.assert_allclose(data['w'].iloc[-1], -0.1775)
    npt.assert_allclose(data['c'].iloc[-1], 658.2624e-6)
    npt.assert_allclose(data['q'].iloc[-1], 9.404386e-3)
    npt.assert_allclose(data['T'].iloc[-1], 28.35199 + 273.15)
    npt.assert_allclose(data['P'].iloc[-1], 100.1938e3)

    npt.assert_allclose(data['u'].mean(), 1.43621, atol=1e-4)
    npt.assert_allclose(data['v'].mean(), -0.634818, atol=1e-4)
    npt.assert_allclose(data['w'].mean(), 0.0619483, atol=1e-4)
    npt.assert_allclose(data['c'].mean(), 659.052e-6, atol=1e-9)
    npt.assert_allclose(data['q'].mean(), 9.56732e-3, atol=1e-7)
    npt.assert_allclose(data['T'].mean(), 28.5431 + 273.15, atol=1e-4)
    npt.assert_allclose(data['P'].mean(), 100.179e3, atol=1e0)


def assert_tob_read(data):
    npt.assert_allclose(data['u'].iloc[0], -2.57175016)
    npt.assert_allclose(data['v'].iloc[0], 1.6450001)
    npt.assert_allclose(data['w'].iloc[0], -0.12725)
    npt.assert_allclose(data['c'].iloc[0], 612.54150391e-6)
    npt.assert_allclose(data['q'].iloc[0], 13.11471748e-3)
    npt.assert_allclose(data['T'].iloc[0], 23.29580116 + 273.15)
    npt.assert_allclose(data['P'].iloc[0], 85.04070282e3)
    npt.assert_allclose(data['u'].iloc[-1], -0.026)
    npt.assert_allclose(data['v'].iloc[-1], -4.57375002)
    npt.assert_allclose(data['w'].iloc[-1], -0.07250001)
    npt.assert_allclose(data['c'].iloc[-1], 620.06219482e-6)
    npt.assert_allclose(data['q'].iloc[-1], 14.72776413e-3)
    npt.assert_allclose(data['T'].iloc[-1], 22.66570282 + 273.15)
    npt.assert_allclose(data['P'].iloc[-1], 85.06277466e3)


if __name__ == '__main__':
    test_hfdata_read_csv()
