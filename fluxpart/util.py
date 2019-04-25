from collections import namedtuple
from itertools import permutations
from math import exp
import os
import warnings
import zipfile

import numpy as np
import pandas as pd
import pywt

from .constants import MOLECULAR_WEIGHT as MW  # kg/mol
from .constants import SPECIFIC_GAS_CONSTANT as Rgas  # m^3 Pa / (kg K)


NP_TYPE = {
    "IEEE4": np.float32,
    "IEEE8": np.float64,
    "LONG": np.int32,
    "ULONG": np.uint32,
}


class HFDataReadWarning(UserWarning):
    pass


def chunked_df(dataframes, time_interval):
    """Partition time-indexed dataframe sequence into time intervals."""

    if time_interval == -1:
        yield pd.concat(dataframes)
        return
    if time_interval is None:
        for df in dataframes:
            yield df
        return

    next_df = next(dataframes)
    consumed_dfs = [next_df]
    current_interval = next_df.index[0].floor(time_interval)
    for next_df in dataframes:
        consumed_dfs.append(next_df)
        if next_df.index[-1].floor(time_interval) == current_interval:
            continue

        # If we reach here, then at least one new interval has been
        # found. We yield sequentially all read intervals except the
        # last one because more of that interval could still be
        # coming

        df = pd.concat(consumed_dfs)
        # gdf = df.groupby(df.index.floor(time_interval))
        gdf = df.groupby(pd.Grouper(freq=time_interval))
        group = iter(gdf)
        for _ in range(len(gdf) - 1):
            interval, chunk = next(group)
            # len=0 can happen if a time gap > interval exists in data
            if len(chunk) == 0:
                continue
            yield chunk
        current_interval, next_df = next(group)
        consumed_dfs = [next_df]

    # Yield remaining chunks. Should always be at least one.
    # Can be more if the main loop is skipped (len(dataframes)=1)
    # and time_interval < time span of the single dataframe

    df = pd.concat(consumed_dfs)
    gdf = df.groupby(df.index.floor(time_interval))
    for interval, chunk in gdf:
        if len(chunk) <= 1:
            continue
        yield chunk


def multifile_read_csv(files, *args, **kwargs):
    """Buffered pd.read_csv of data split across multiple files."""
    for file_ in files:
        try:
            df = pd.read_csv(file_, *args, **kwargs)
        except Exception as e:
            mssg = "Skipping file " + str(file_) + " because " + e.args[0]
            warnings.warn(mssg, HFDataReadWarning)
            continue
        if isinstance(df, pd.io.parsers.TextFileReader):
            yield from df
        else:
            yield df


def multifile_read_ghg(files, *args, **kwargs):
    """Buffered pd.read_csv of data split across multiple files."""
    for file_ in files:
        with zipfile.ZipFile(file_) as z:
            fn = os.path.basename(file_)
            try:
                df = pd.read_csv(z.open(fn[:-3] + "data"), *args, **kwargs)
            except Exception as e:
                mssg = "Skipping file " + str(file_) + " because " + e.args[0]
                warnings.warn(mssg, HFDataReadWarning)
                continue
        if isinstance(df, pd.io.parsers.TextFileReader):
            yield from df
        else:
            yield df


def multifile_read_tob1(tobfiles, count=-1):
    """Buffered read of multiple tob1 files into dataframes."""
    for tobfile in tobfiles:
        try:
            df = dataframe_read_tob1(tobfile, count)
        except Exception as e:
            mssg = "Skipping file " + str(tobfile) + " because " + e.args[0]
            warnings.warn(mssg, HFDataReadWarning)
            continue
        yield df


def dataframe_read_tob1(tobfile, count=-1):
    return pd.DataFrame(ndarray_read_tob1(tobfile, count))


def ndarray_read_tob1(tobfile, count=-1):
    """Read TOB1 data file into structured numpy array."""
    with open(tobfile, "rb") as f:
        f.readline()
        names = f.readline().decode().strip().replace('"', "").split(",")
        f.readline()
        f.readline()
        types = f.readline().decode().strip().replace('"', "").split(",")
        dtype = np.dtype([(n, NP_TYPE[t]) for n, t in zip(names, types)])
        return np.fromfile(f, dtype=dtype, count=count)


def stats2(sarray, names=None):
    """Calculate means and (co)variances for structured array data."""

    if names is None:
        names = sarray.dtype.names
    nvar = len(names)
    data = tuple(sarray[name] for name in names)
    cov = np.cov(data)
    nondiag_cov = list(cov[i, j] for i, j in permutations(range(nvar), 2))

    names_ave = list("ave_" + name for name in names)
    names_var = list("var_" + name for name in names)
    names_cov = list(
        "cov_" + n1 + "_" + n2 for n1, n2 in permutations(names, 2)
    )

    out = dict(zip(names_ave, np.mean(data, axis=1)))
    out.update(zip(names_var, cov.diagonal()))
    out.update(zip(names_cov, nondiag_cov))

    NamedStats = namedtuple("Stats2", names_ave + names_var + names_cov)
    return NamedStats(**out)


def sat_vapor_press(t_kelvin):
    # e_sat
    tr = 1 - 373.15 / t_kelvin
    arg = 13.3185 * tr - 1.9760 * tr ** 2 - 0.6445 * tr ** 3 - 0.1299 * tr ** 4
    return 101325.0 * exp(arg)


def vapor_press_deficit(rho_vapor, t_kelvin):
    return sat_vapor_press(t_kelvin) - rho_vapor * Rgas.vapor * t_kelvin


def vapor_press_deficit_mass(rho_vapor, t_kelvin):
    return vapor_press_deficit(rho_vapor, t_kelvin) / Rgas.vapor / t_kelvin


def qflux_mass_to_heat(massflux, Tk):  # kg/m^2/s, K
    Lv = 2.5008e6 - 2366.8 * (Tk - 273.15)  # J/kg
    return massflux * Lv  # W/m^2


def cflux_mass_to_mol(massflux):  # kg/m^2/s
    return 1.0 / MW.co2 * massflux  # mol/m^2/s


def qflux_mass_to_mol(massflux):  # kg/m^2/s
    return 1.0 / MW.vapor * massflux  # mol/m^2/s


def progressive_lowcut_series(series):
    """Progressively remove low-frequency components of 1D series.

    Yields sequence in which the low-frequency (large-scale) components
    of `series` are progressively removed.  The sequence is obtained
    from reconstruction of a multilevel discrete haar wavelet
    decomposition of `series`.

    N.B.: The length of `series` is assumed to be a power of 2
    (does not check!)

    Parameters
    ----------
    series : array_like
        1D data series with a length that is a power of 2

    Yields
    -------
    lowcut_series : array
        Sequence of progressively lowcut filtered data `series`. The
        yielded series have the same length as `series`.

    Notes
    -----
    After an N level discrete wavelet decomposition, a data series S can
    be reconstructed in terms of wavelet 'approximations' (A) and
    'details' (D)::

        S = A(N) + D(N) + D(N-1) ... D(2) + D(1)
          = A(N-1) + D(N-1) + ... D(2) + D(1)     [A(N-1) = A(N) + D(N)]
            ...
          = A(1) + D(1)                           [(A(1) = A(2) + D(2)]

    where A(N) represents the 'lowest' level approximation. For
    the haar wavelet and data S having a length that is a power of 2,
    A(N) is equal to mean(S)]

    The sequence returned by this function is::

        S - A(N),
        S - A(N-1),
        S - A(N-2),
        ...,
        S - A(1)

    This sequence is computed by the equivalent::

        S - A(N),
        S - A(N) - D(N),
        S - A(N) - D(N) - D(N-1),
        ...,
        S - A(N) - D(N) - D(N-1) - ... - D(2),

    i.e. the details are removed in sequence::

        lowcut = S - A(N)
        for j = N to 2
            lowcut -= D(j)

    """

    series_data = np.asarray(series)
    wavelet = pywt.Wavelet("haar")
    nlevel = pywt.dwt_max_level(series_data.size, wavelet.dec_len)
    decomp_coef = pywt.wavedec(series_data, wavelet=wavelet, level=nlevel)
    cAn, cD = decomp_coef[0], decomp_coef[1:]
    lowcut_series = series_data - pywt.upcoef("a", cAn, wavelet, level=nlevel)
    yield lowcut_series
    for j, cDj in enumerate(cD[:-1]):
        Dj = pywt.upcoef("d", cDj, wavelet, level=nlevel - j)
        lowcut_series -= Dj
        yield lowcut_series


if __name__ == "__main__":
    pass
