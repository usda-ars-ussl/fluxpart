"""
High-frequency eddy covariance time series data.

The following notation is used in variable naming and documentation to
represent meteorological quantities (SI units)::

    u, v, w = wind velocities (m/s)
    q = water vapor mass concentration (kg/m^3)
    c = carbon dioxide mass concentration (kg/m^3)
    T = air temperature (K)
    P = total air pressure (Pa)

"""
from glob import glob
import os

import math
import numpy as np
import pandas as pd

from . import util
from .containers import HFSummary
from .constants import MOLECULAR_WEIGHT as MW
from .constants import SPECIFIC_HEAT_CAPACITY as CP
from .constants import SPECIFIC_GAS_CONSTANT as GC


class Error(Exception):
    pass


class HFDataReadError(Error):
    def __init__(self, message):
        self.message = message


class TooFewDataError(Error):
    def __init__(self, data_frac, rd_tol, len_max_slice, ad_tol):
        self.message = (
            "HF Data read but rejected because the longest continuous "
            "run of valid data was too short on a relative (length "
            "data / total length) = {:.4} < rd_tol = {:.4}) and/or "
            "absolute basis (data length = {} < ad_tol = {})"
            "".format(data_frac, rd_tol, len_max_slice, ad_tol)
        )


VAR_NAMES = ["u", "v", "w", "c", "q", "T", "P"]


class HFDataReader(object):
    """Reader for high-frequency eddy covariance data.

    Parameters
    ----------
    filetype : {'csv', 'tob1'}
        'csv' = delimited text file (default); 'tob1' = Campbell
        Scientific binary format file.
    cols : 7*(int,), optional
        Column indices for (u, v, w, q, c, T, P) data, in that order.
        0-based indexing. Default is (2, 3, 4, 6, 5, 7, 8).
    time_col : int, optional
        Datetime column for `csv` data. Default is None.
    converters : dict, optional
        Dictionary of functions used to convert any non-SI data to
        SI units.  Dict keys are 'u', 'v', 'w', 'q', 'c', 'T', or
        'P'. Funcs take a single argument, e.g.
        ``converters={'P': lambda arg: 1e3 * arg}``.
    flags : 2-tuple or list of 2-tuples, optional
        Specifies that one or more data columns are used to flag bad
        data records. Each tuple is of the form (col, goodval), where
        col is an int specifying the column number containing the flag
        (0-based indexing), and goodval is the value of the flag that
        indicates a good data record.
    **kwargs
        Passed to pandas.read_csv_ when filetype is csv. Should not
        include `usecols`  or `header` keywords.


    .. _pandas.read_csv:
        https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

    """

    def __init__(
        self,
        filetype="csv",
        cols=(2, 3, 4, 5, 6, 7, 8),
        time_col=None,
        converters=None,
        flags=None,
        **kwargs
    ):
        self._filetype = filetype.strip().lower()
        self._cols = cols
        self._time_col = time_col
        self._converters = {} if converters is None else converters
        if flags is None:
            self._flags = []
        elif not isinstance(flags, list):
            self._flags = [flags]
        else:
            self._flags = flags
        self._readcsv_kws = kwargs

    @property
    def _namecols(self):
        namecols = dict(zip(VAR_NAMES, self._cols))
        flags = {"flag-" + str(col): col for col, val in self._flags}
        namecols.update(flags)
        if isinstance(self._time_col, int) and self._filetype == "csv":
            namecols["Datetime"] = self._time_col
        return namecols

    @property
    def _names(self):
        # sorted because pd.read_csv sorts usecols but not names
        return sorted(self._namecols, key=self._namecols.get)

    @property
    def _usecols(self):
        return [self._namecols[k] for k in self._names]

    def peektime(self, fname):
        try:
            if self._filetype == "csv":
                df = self._read_csv(fname, nrows=1)
            elif self._filetype in ("tob", "tob1"):
                df = self._read_tob1(fname, count=5)
            elif self._filetype == "pd.df":
                df = fname
            else:
                raise HFDataReadError("Unknown file type")
        except Exception as err:
            raise HFDataReadError(err.args[0])
        return df.index[0]

    def read(self, fname, count=-1, **kwargs):
        """Read high frequency eddy covariance data into dataframe.

        Parameters
        ----------
        fname : str
            Filename.
        **kwargs
            Passed to pandas.read_csv_ when filetype is csv. Will
            override any duplicate kwargs passed to the constructor.
            Should not include `usecols`  or `header` keywords.


        .. _pandas.read_csv:
            https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

        """
        try:
            if self._filetype == "csv":
                dataframe = self._read_csv(fname, **kwargs)
            elif self._filetype in ("tob", "tob1"):
                dataframe = self._read_tob1(fname, count)
            elif self._filetype == "pd.df":
                dataframe = self._read_df(fname)
            else:
                raise HFDataReadError("Unknown file type")
        except Exception as err:
            raise HFDataReadError(err.args[0])

        dataframe = self._unit_convert(dataframe)
        dataframe = self._flagvals_to_mask(dataframe)
        return dataframe

    def _read_csv(self, csvfile, **kwargs):
        df = pd.read_csv(
            csvfile,
            usecols=self._usecols,
            header=None,
            **{**self._readcsv_kws, **kwargs}
        )
        df.columns = self._names
        if self._time_col is not None:
            df["Datetime"] = pd.to_datetime(df.iloc[:, self._time_col])
            df = df.set_index("Datetime")
        return df

    def _read_tob1(self, tobfile, count=-1):
        df = pd.DataFrame(util.tob1_to_array(tobfile, count))
        df["Datetime"] = pd.to_datetime(
            arg=df.loc[:, "SECONDS"] + 10 ** -9 * df.loc[:, "NANOSECONDS"],
            unit="s",
            origin="1990-01-01",
        )
        df = df.set_index("Datetime")
        return self._read_df(df)

    def _read_df(self, df):
        df = df.iloc[:, self._usecols]
        df.columns = self._names
        return df

    def _unit_convert(self, df):
        for var, func in self._converters.items():
            df.loc[:, var] = func(df.loc[:, var])
        return df

    def _flagvals_to_mask(self, df):
        """Mask is True if flag does not equal the good value."""
        for col, val in self._flags:
            df.loc[:, "flag-" + str(col)] = (
                df.loc[:, "flag-" + str(col)] != val
            )
        return df


class HFData(object):
    def __init__(self, dataframe):
        self.dataframe = dataframe
        self._already_corrected_external = False

    def __getitem__(self, name):
        """Column-wise get without specifying dataframe attribute"""
        return self.dataframe.loc[:, name]

    def __setitem__(self, name, value):
        """Column-wise set without specifying dataframe attribute"""
        self.dataframe.loc[:, name] = value

    def cleanse(self, bounds=None, rd_tol=0.5, ad_tol=1024):
        """Apply some data QC/QA, remove bad data.

        If problems are found, self.dataframe is modified to contain
        only the longest contiguous stretch of good data.

        Parameters
        ----------
        bounds : dict, optional
            Dictionary specifying any prescribed lower and upper bounds
            for legal data. Dict entries have the form ``varname:
            (float, float)``, where varname is one of 'u', 'v', 'w',
            'q', 'c', 'T', or 'P', and the 2-tuple holds values for the
            lower and upper bounds: ``(lower, upper)``.  Data records
            are rejected if a variable in the record is outside the
            prescribed bounds. Default is None.
        rd_tol : float, optional
            Relative tolerance for rejecting the datafile. Default is
            `rd_tol` = 0.4.  See `ad_tol` for explanation.
        ad_tol : int, optional
            Absolute tolerance for rejecting the datafile. Defaults is
            `ad_tol` = 1024. If the datafile contains bad records (not
            readable, out-of-bounds, or flagged data), the partitioning
            analysis is performed using the longest stretch of
            consecutive good data records found, unless that stretch is
            too short, in which case the analysis is aborted. The
            criteria for judging 'too short' can be specified in both
            relative and absolute terms: the datafile is rejected if the
            good stretch is a fraction of the total data that is less
            than `rd_tol`, and/or is less than `ad_tol` records long.

        """
        bounds = bounds or {}
        data = self.dataframe

        # 1D mask is True for a row if any data are nan, any flag is
        # True, or any data are out-of-bounds
        varindx = pd.Index(VAR_NAMES)
        mask = data.loc[:, varindx].isnull().any(axis=1)
        mask |= data.loc[:, data.columns.difference(varindx)].any(axis=1)
        for var, (low, high) in bounds.items():
            mask |= (data[var] < low) | (data[var] > high)

        # Find longest span of valid (unmasked) data
        marray = np.ma.array(np.zeros([data.shape[0]]), mask=mask.values)
        unmasked_slices = np.ma.clump_unmasked(marray) or [slice(0, 0)]
        max_indx = np.argmax([s.stop - s.start for s in unmasked_slices])
        max_slice = unmasked_slices[max_indx]
        len_max_slice = max_slice.stop - max_slice.start

        # verify sufficient data length
        data_frac = len_max_slice / data.shape[0]
        if data_frac < rd_tol or len_max_slice < ad_tol:
            self.dataframe = None
            raise TooFewDataError(data_frac, rd_tol, len_max_slice, ad_tol)

        self.dataframe = data.iloc[max_slice]
        return

    def correct_external(self):
        """Adjust q and c data series to correct for external effects.

        Water vapor and carbon dioxide series data in the dataframe are
        corrected for external fluctuations associated with air
        temperature and vapor density. See: [WPL80]_ and [DK07]_.

        """

        if self._already_corrected_external:
            return
        ave_vapor = self["q"].mean()
        ave_co2 = self["c"].mean()
        ave_T = self["T"].mean()
        dev_vapor = self["q"] - ave_vapor
        dev_T = self["T"] - ave_T

        Pdryair = self["P"].mean() - ave_vapor * GC.vapor * ave_T
        rho_totair = ave_vapor + Pdryair / GC.dryair / ave_T

        specific_vapor = ave_vapor / rho_totair
        specific_co2 = ave_co2 / rho_totair
        mu = MW.dryair / MW.vapor
        muq = mu * specific_vapor
        muc = mu * specific_co2

        self["q"] += muq * dev_vapor + (1 + muq) * ave_vapor * dev_T / ave_T
        self["c"] += muc * dev_vapor + (1 + muq) * ave_co2 * dev_T / ave_T
        self._already_corrected_external = True
        return

    def summarize(self):
        """Summarize high frequency dataframe statistics.

        Returns
        -------
        namedtuple
            :class:`~fluxpart.containers.HFSummary`

        """
        hfs = util.stats2(self.dataframe, VAR_NAMES)
        Pvap = hfs.ave_q * GC.vapor * hfs.ave_T
        rho_dryair = (hfs.ave_P - Pvap) / GC.dryair / hfs.ave_T
        rho_totair = rho_dryair + hfs.ave_q
        Cp = CP.dryair * (1 + 0.84 * hfs.ave_q / rho_totair)

        return HFSummary(
            T=hfs.ave_T,
            P=hfs.ave_P,
            Pvap=Pvap,
            ustar=(hfs.cov_w_u ** 2 + hfs.cov_w_v ** 2) ** 0.25,
            wind_w=hfs.ave_w,
            var_w=hfs.var_w,
            rho_vapor=hfs.ave_q,
            var_vapor=hfs.var_q,
            rho_co2=hfs.ave_c,
            var_co2=hfs.var_c,
            corr_q_c=hfs.cov_q_c / math.sqrt(hfs.var_q * hfs.var_c),
            H=rho_totair * Cp * hfs.cov_w_T,
            cov_w_q=hfs.cov_w_q,
            cov_w_c=hfs.cov_w_c,
            rho_dryair=rho_dryair,
            rho_totair=rho_totair,
            cov_w_T=hfs.cov_w_T,
            N=self.dataframe.shape[0],
        )

    def truncate_pow2(self):
        """Truncate dataframe length to largest possible power of 2."""
        truncate_len = 2 ** int(np.log2(self.dataframe.shape[0]))
        self.dataframe = self.dataframe.iloc[:truncate_len]


class HFDataSequence(object):
    """Sequence of high frequency data objects.

    Parameters
    ----------
    filepaths : str or list
        High frequency data file, data directory, or list of either.
    reader : :class:`~fluxpart.containers.HFDataReader`
    kind : {'infer', 'file', 'dir'}, optional
        Indicates whether `filepath` refers to file(s) or directory(s).
        Default ('infer') tries to determine the kind automatically.
    ext : str, optional
        When `filepath` refers to a directory, all files in the
        directory with extension `ext` will be read. Default ("") reads
        all files regardless of extension. When specifying `ext`,
        include the 'dot' where appropriate (e.g., ``ext=".dat"``)


    """

    def __init__(self, filepaths, reader, kind="infer", ext=""):
        self._ext = ext
        self._filepaths = filepaths
        self._kind = kind
        self.reader = reader

    def files(self):
        filepaths = self._filepaths
        if isinstance(filepaths, str):
            filepaths = [filepaths]
        if self._kind.lower() == "infer":
            if os.path.isfile(filepaths[0]):
                kind = "file"
            elif os.path.isdir(filepaths[0]):
                kind = "dir"
            else:
                raise HFDataReadError("Unable to infer sequence type")
        if kind.lower() == "dir":
            unsorted_files = []
            for p in filepaths:
                unsorted_files += glob(os.path.join(p, "*" + self._ext))
        elif kind.lower() == "file":
            unsorted_files = filepaths
        else:
            raise HFDataReadError("Sequence type not recognized")

        datafiles = [(self.reader.peektime(f), f) for f in unsorted_files]
        datafiles.sort(key=lambda p: p[0])
        return datafiles

    def chunk(self, interval="30min"):
        """Generator that reads/partitions data by time interval.

        Parameters
        ----------
        interval : str
            Offset alias_ specifying the chunking interval. Default is '30min'.

        Notes
        -----
        When `files` is a single file, then a rough equivalent is:

            df = HFData.read(file)
            for time, group in df.groupby(df.index.floor(interval):
                yield time, group


        .. _alias:
            http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases

        """
        timesfiles = self.files()
        dataframes = []
        current_interval = timesfiles[0][0].floor(interval)
        for time, _file in timesfiles:
            try:
                hfdat = self.reader.read(_file)
            except HFDataReadError:
                # TODO: raise Warning (or Error?)
                continue
            dataframes.append(hfdat)
            if hfdat.index[-1].floor(interval) == current_interval:
                continue

            # If we reach here, then at least one new interval has been
            # found. We yield all read intervals except the very last
            # one because more of that interval could still be in the
            # next file to be read

            df = pd.concat(dataframes)
            gdf = df.groupby(df.index.floor(interval))
            group = iter(gdf)
            for _ in range(len(gdf) - 1):
                yield next(group)
            current_interval, df = next(group)
            dataframes = [df]

        # should always be one remaining interval to clean up
        yield current_interval, pd.concat(dataframes)


def get_hfdata(fname, *args, **kws):
    """Convenience function for reading HFData."""
    return HFData(HFDataReader(*args, **kws).read(fname))
