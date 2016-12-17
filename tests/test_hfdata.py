import os, io
import numpy as np
import numpy.testing as npt
from fluxpart.hfdata import HFData
from fluxpart.fluxpart import _converter_func

TESTDIR = os.path.dirname(os.path.realpath(__file__))


def test_hfdata_from_txt():
    cols = (2, 3, 4, 6, 5, 7, 8)
    fname = os.path.join(TESTDIR,
                         'data/TOA5_6843.ts_Above_2012_06_07_1300.dat')

    data = HFData(
        fname,
        cols,
        converters={
            'T': _converter_func(1, 273.15),
            'q': _converter_func(1e-3, 0),
            'c': _converter_func(1e-6, 0),
            'P': _converter_func(1e3, 0)},
        flags=(9, 0),
        delimiter=",",
        skip_header=4)

    npt.assert_allclose(data['u'][0], 0.468)
    npt.assert_allclose(data['v'][0], -0.9077501)
    npt.assert_allclose(data['w'][0], 0.1785)
    npt.assert_allclose(data['c'][0], 659.7584e-6)
    npt.assert_allclose(data['q'][0], 9.530561e-3)
    npt.assert_allclose(data['T'][0], 28.52527 + 273.15)
    npt.assert_allclose(data['P'][0], 100.1938e3)
    npt.assert_allclose(data['u'][-1], 1.3675)
    npt.assert_allclose(data['v'][-1], -0.75475)
    npt.assert_allclose(data['w'][-1], -0.1775)
    npt.assert_allclose(data['c'][-1], 658.2624e-6)
    npt.assert_allclose(data['q'][-1], 9.404386e-3)
    npt.assert_allclose(data['T'][-1], 28.35199 + 273.15)
    npt.assert_allclose(data['P'][-1], 100.1938e3)

    npt.assert_allclose(data['u'].mean(), 1.43621, atol=1e-4)
    npt.assert_allclose(data['v'].mean(), -0.634818, atol=1e-4)
    npt.assert_allclose(data['w'].mean(), 0.0619483, atol=1e-4)
    npt.assert_allclose(data['c'].mean(), 659.052e-6, atol=1e-9)
    npt.assert_allclose(data['q'].mean(), 9.56732e-3, atol=1e-7)
    npt.assert_allclose(data['T'].mean(), 28.5431 + 273.15, atol=1e-4)
    npt.assert_allclose(data['P'].mean(), 100.179e3, atol=1e0)

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
        'asdf,11,-2,3,4,5,6,7,9,0\n')

    toy = HFData(
        io.BytesIO(toy_data.encode()),
        cols=(1, 2, 3, 7, 6, 4, 5),
        comments='#',
        skip_header=1,
        missing_values="???",
        converters={'q': _converter_func(10., 0)},
        flags=(9, 0),
        delimiter=",",
        rd_tol=0.1,
        ad_tol=2,
        bounds={'v': (0, np.inf)})

    npt.assert_allclose(toy.data_table['u'], [4, 5, 6])
    npt.assert_allclose(toy.data_table['v'], 3 * [2, ])
    npt.assert_allclose(toy.data_table['w'], 3 * [3, ])
    npt.assert_allclose(toy.data_table['q'], 3 * [70, ])
    npt.assert_allclose(toy.data_table['c'], 3 * [6, ])
    npt.assert_allclose(toy.data_table['T'], 3 * [4, ])
    npt.assert_allclose(toy.data_table['P'], 3 * [5, ])

if __name__ == '__main__':
    test_hfdata_from_txt()
