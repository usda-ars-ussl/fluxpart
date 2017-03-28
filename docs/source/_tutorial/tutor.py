import glob
import os
import pickle
import textwrap
import zipfile

import fluxpart as fp
from fluxpart.plots import plot_fluxes
from fluxpart.plots import plot_hfdata
from fluxpart.hfdata import HFData
from fluxpart.fluxpart import _converter_func


DATADIR = './data'
RESDIR = './results'
FIGDIR = '../images'
TXTDIR = '../text'
ZIPFILE = '05-08June2012.zip'
TUTORPKL = 'tutor_fpout.pkl'
ZIPPKL = 'zip_fpout.pkl'

QUICKFILE = 'TOA5_6843.ts_Above_2012_06_07_1245.dat'

# BEGIN/END tags and redundant imports are for sphinx literalinclude directive


def make_quickdata():
    fname = QUICKFILE
    if not os.path.exists(os.path.join(DATADIR, fname)):
        with zipfile.ZipFile(os.path.join(DATADIR, ZIPFILE), 'r') as zarch:
            zarch.extract(fname, DATADIR)


def quickstart():
    make_quickdata()
    # BEGIN quickstart
    import fluxpart as fp
    out = fp.flux_partition(
        fname='data/TOA5_6843.ts_Above_2012_06_07_1245.dat',
        hfd_options={
            'cols': (2, 3, 4, 6, 5, 7, 8),
            'unit_convert': {'q': 1e-3, 'c': 1e-6, 'P': 1e3},
            'temper_unit': 'C',
            'delimiter': ",",
            'skip_header': 4},
        wue_options={
            'meas_ht': 7.11,
            'canopy_ht': 4.42,
            'ppath': 'C3'})
    # END quickstart
    return out


def known_wue():
    make_quickdata()
    # BEGIN known_wue
    import fluxpart as fp
    out = fp.flux_partition(
        fname='data/TOA5_6843.ts_Above_2012_06_07_1245.dat',
        hfd_options={
            'cols': (2, 3, 4, 6, 5, 7, 8),
            'unit_convert': {'q': 1e-3, 'c': 1e-6, 'P': 1e3},
            'temper_unit': 'C',
            'delimiter': ",",
            'skip_header': 4},
        meas_wue=-6.5e-3)  # kg CO2 / kg H2O
    # END known_wue
    return out


# BEGIN my_fpart
import fluxpart as fp
def my_fpart(datafile, timestamp):
    return fp.flux_partition(
        fname=datafile,
        label=timestamp,
        hfd_options={
            'cols': (2, 3, 4, 6, 5, 7, 8),
            'unit_convert': {'q': 1e-3, 'c': 1e-6, 'P': 1e3},
            'temper_unit': 'C',
            'delimiter': ",",
            'skip_header': 4,
            'flags': (9, 0)},
        wue_options={
            'meas_ht': 7.11,
            'canopy_ht': 4.42,
            'ppath': 'C3'})
# END my_fpart


def tutor_example():
    # BEGIN tutor_example
    import datetime
    from itertools import product
    import os
    DATADIR = "./data"
    fpout = []
    for day, hour, minu in product(range(5, 9), range(0, 24), range(0, 60, 15)):
        timestamp = datetime.datetime(2012, 6, day, hour, minu)
        filename = (
            "TOA5_6843.ts_Above_2012_06_{:02}_{:02}{:02}.dat".format(day, hour, minu))
        datafile = os.path.join(DATADIR, filename)
        fpout.append(my_fpart(datafile, timestamp))
    # END tutor_example
    return fpout


def tutor_fpout():
    pklout = os.path.join(RESDIR, TUTORPKL)
    if os.path.exists(pklout):
        with open(pklout, 'rb') as f:
            fpout = pickle.load(f)
    else:
        print('re-creating fp results ... this could take a few minutes')
        zip_archive = os.path.join(DATADIR, ZIPFILE)
        with zipfile.ZipFile(zip_archive, 'r') as zfile:
            zfile.extractall(DATADIR)
        fpout = tutor_example()
        make_clean_dat()
        os.makedirs(RESDIR, exist_ok=True)
        with open(pklout, 'wb') as f:
            pickle.dump(fpout, f)
    return fpout


def cat_quickstart_data():
    with zipfile.ZipFile(os.path.join(DATADIR, ZIPFILE), 'r') as zarch:
        with zarch.open(QUICKFILE, 'r') as datafile:
            with open(os.path.join(TXTDIR, 'cat_quickstart.txt'), 'w') as f:
                f.write('$ cat ' + QUICKFILE + '\n')
                for _ in range(0, 7):
                    f.write(next(datafile).decode('utf-8').strip('\r\n') + '\n')
                for line in datafile:
                    pass
                f.write('...\n')
                f.write(line.decode('utf-8').strip('\r\n'))


def quickstart_out(field):
    fname = 'quickstart_' + field + '_out.txt'
    out = quickstart()
    with open(os.path.join(TXTDIR, fname), 'w') as f:
        f.write(">>> print(out['{}'])\n".format(field) + out[field].__str__())


def quickstart_rawout(field):
    fname = 'quickstart_' + field + '_rawout.txt'
    out = quickstart()
    with open(os.path.join(TXTDIR, fname), 'w') as f:
        f.write(">>> out['{}']\n".format(field))
        for line in textwrap.wrap(out[field].__repr__()):
            f.write(line + '\n')


def hfdata_plot():
    make_quickdata()
    unit_convert = {'q': 1e-3, 'c': 1e-6, 'P': 1e3}
    converters = {
        k: _converter_func(float(v), 0.) for k, v in unit_convert.items()}
    converters['T'] = _converter_func(1., 273.15)
    with zipfile.ZipFile(os.path.join(DATADIR, ZIPFILE), 'r') as zarch:
        with zarch.open(QUICKFILE, 'r') as datafile:
            data = HFData(
                fname=datafile,
                cols=(2, 3, 4, 6, 5, 7, 8),
                converters=converters,
                delimiter=",",
                skip_header=4)
    fig = plot_hfdata(data)
    return fig


def make_clean_dat():
    for fname in glob.glob(os.path.join(DATADIR, '*.dat')):
        os.remove(fname)


def main():
    fpout = tutor_fpout()
    fig = plot_fluxes(fpout[0:96], '05 June 2012')
    fig.savefig(os.path.join(FIGDIR, 'tutor_05June2012_fluxes.png'))
    fig = plot_fluxes(fpout, '05-08 June 2012')
    fig.savefig(os.path.join(FIGDIR, 'tutor_05-08June2012_fluxes.png'))
    fig = hfdata_plot()
    fig.savefig(os.path.join(FIGDIR, 'tutor_quickstart_hfdata.png'))
    cat_quickstart_data()
    quickstart_out('fluxes')
    quickstart_rawout('fluxes')
    quickstart_out('result')
    quickstart_out('hfsummary')
    quickstart_out('wue')
    quickstart_out('numersoln')
    quickstart_out('qcdata')


if __name__ == '__main__':
    main()
