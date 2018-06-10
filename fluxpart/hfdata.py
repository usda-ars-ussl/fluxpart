"""
Read and model eddy covariance high-frequency time series data.

The following notation is used in variable naming to represent
meteorological quantities (SI units)::

    u, v, w = wind velocities (m/s)
    q = water vapor mass concentration (kg/m^3)
    c = carbon dioxide mass concentration (kg/m^3)
    T = air temperature (K)
    P = total air pressure (Pa)

"""
from copy import deepcopy
from glob import glob
import os

import attr
import math
import numpy as np
import pandas as pd

from . import util
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
            "The longest contiguous run of valid data was too short:\n"
            "(length data / total length) = {:.4} < rd_tol = {:.4}\n"
            "and/or data length = {} < ad_tol = {}"
            "".format(data_frac, rd_tol, len_max_slice, ad_tol)
        )


VAR_NAMES = ["u", "v", "w", "c", "q", "T", "P"]


class HFData(object):
    """
    Parameters
    ----------
    dataframe
        High frequency eddy covariance dataframe. Must include columns
        for ["u", "v", "w", "c", "q", "T", "P"]. Normally, dataframe
        should have a datetime index. Dataframe may include also boolean
        mask columns.

    """

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
        """Apply some data QC, remove bad data.

        If problems are found (data not readable, out-of-bounds, or
        flagged), self.dataframe is modified to contain only the longest
        contiguous stretch of good data. An error is raised if the
        resulting data are too few. The criteria for judging 'too few'
        can be specified in both relative and absolute terms: the
        datafile is rejected if the good stretch is a fraction of the
        total data that is less than `rd_tol`, and/or is less than
        `ad_tol` records long.

        Parameters
        ----------
        bounds : dict, optional
            Dictionary specifying lower and upper bounds for legal data.
            Dict entries have the form ``varname: (lower, upper)``,
            where varname is one of 'u', 'v', 'w', 'q', 'c', 'T', or
            'P', and the 2-tuple holds values for the lower and upper
            bounds. Default is None.
        rd_tol : float, optional
            Relative tolerance for rejecting the datafile. Default is
            `rd_tol` = 0.5.
        ad_tol : int, optional
            Absolute tolerance for rejecting the datafile. Default is
            `ad_tol` = 1024.

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
        :class:`~fluxpart.hfdata.HFSummary`

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


class HFDataSource(object):
    """Reader for high-frequency eddy covariance data.

    Parameters
    ----------
    filetype : {'csv', 'tob1'}
        'csv' = delimited text file; 'tob1' = Campbell Scientific binary
        format file.
    cols : 7*(int,)
        Column indices for (u, v, w, q, c, T, P) data, in that order.
        0-based indexing.
    converters : dict, required if datafile uses any non-SI units
        Dictionary of functions used to convert any non-SI data to
        SI units.  Dict keys are 'u', 'v', 'w', 'q', 'c', 'T', or
        'P'. Functions take a single argument, e.g.
        ``converters={'P': lambda arg: 1e3 * arg}``. If all data are
        SI units, set `converters` to None (default).
    time_col : int, optional
        Datetime column for `csv` filetype. Default is None.
    flags : 2-tuple or list of 2-tuples, optional
        Indicate that one or more data columns are used to flag bad data
        records. Each tuple is of the form (col, goodval), where col is
        an int specifying the column number (0-based indexing), and
        goodval is the flag value indicating a good data record.
    **kwargs
        Passed to pandas.read_csv_ when filetype is csv. Should not
        include `usecols` or `header` keywords.


    .. _pandas.read_csv:
        https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

    """

    def __init__(
        self,
        filetype,
        cols,
        converters=None,
        time_col=None,
        flags=None,
        **kwargs,
    ):
        if flags is None:
            flags = []
        elif not isinstance(flags, list):
            flags = [flags]

        self._filetype = filetype.strip().lower()
        self._cols = cols
        self._converters = {} if converters is None else converters
        self._time_col = time_col
        self._flags = flags
        self._csvformat_kws = kwargs

    def _buffered_read(self, datafiles, **kwargs):
        """Generator providing read of data split across files.

        Yields data in file-sized chunks in dataframe format.

        Parameters
        ----------
        datafiles : list of str
            Time-sorted list of datafiles to be read.
        **kwargs
            Passed to :class:`~fluxpart.hfdata.HFDataSource._readfile`.
            'df_output' is ignored.

        """
        kws = deepcopy(kwargs)
        # Ensure we are getting dataframes and not HFData
        kws["df_output"] = True
        for dfile in datafiles:
            try:
                yield self._readfile(dfile, **kws)
            except HFDataReadError:
                # TODO: raise Warning
                pass

    def _chunk(self, buffered_io, chunk_interval):
        """Consume time-indexed data stream; yield time-chunked data.

        Incoming stream is a buffered, time-indexed dataframe. Chunks
        are yielded as a tuple: (datetime, dataframe)

        """
        df = next(buffered_io)
        dataframes = [df]
        current_interval = df.index[0].floor(chunk_interval)
        for frame in buffered_io:
            dataframes.append(frame)
            if frame.index[-1].floor(chunk_interval) == current_interval:
                continue

            # If we reach here, then at least one new interval has been
            # found. We yield sequentially all read intervals except the
            # last one because more of that interval could still be
            # coming

            df = pd.concat(dataframes)
            gdf = df.groupby(df.index.floor(chunk_interval))
            group = iter(gdf)
            for _ in range(len(gdf) - 1):
                yield next(group)
            current_interval, df = next(group)
            dataframes = [df]

        # Yield remaining intervals. Should always be at least one.
        # Can be more if buffered_io contains only one item and loop is
        # skipped with chunk_interval < time span of buffered file.

        df = pd.concat(dataframes)
        gdf = df.groupby(df.index.floor(chunk_interval))
        for interval, group in gdf:
            yield interval, group

    @property
    def _names(self):
        # sorted because pd.read_csv sorts usecols but not names
        return sorted(self._namecols, key=self._namecols.get)

    @property
    def _namecols(self):
        namecols = dict(zip(VAR_NAMES, self._cols))
        flags = {"flag-" + str(col): col for col, val in self._flags}
        namecols.update(flags)
        if isinstance(self._time_col, int) and self._filetype == "csv":
            namecols["Datetime"] = self._time_col
        return namecols

    def _peektime(self, fname):
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

    def reader(self, file_or_dir, interval, ext="", **kwargs):
        """Generator that reads high frequency eddy covariance data.

        Returns a generator that when iterated produces data chunked by
        time interval.

        Parameters
        ----------
        file_or_dir : str or list of str
            High frequency data file, data directory, or list of either.
        interval : str
            Time interval used to partition data. The interval is
            specified using the pandas string alias_ format
            (e.g., ``interval="30min"``).
        ext : str, optional
            When `file_or_dir` refers to a directory, all files in the
            directory with extension `ext` will be read. Default ("")
            reads all files regardless of extension. When specifying
            `ext`, include the 'dot' where appropriate (e.g.,
            ``ext=".dat"``)
        **kwargs
            For csv filetype, kwargs are passed to pandas.read_csv_.
            Can be used to override or add to formatting specified in
            the initializer. Should not include `usecols`  or `header`
            keywords. For tob1 filetype, kwargs is passed to
            numpy.fromfile_. In this case, the only allowable kwarg is
            'count', which can be used to limit the number of lines
            read.


        .. _numpy.fromfile
            https://docs.scipy.org/doc/numpy/reference/generated/numpy.fromfile.html

        .. _pandas.read_csv:
            https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

        .. _alias:
            http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases

        """
        if isinstance(file_or_dir, str):
            file_or_dir = [file_or_dir]
        if os.path.isfile(file_or_dir[0]):
            kind = "file"
        elif os.path.isdir(file_or_dir[0]):
            kind = "dir"
        else:
            raise HFDataReadError("Unable to infer data type")
        if kind.lower() == "dir":
            unsorted_files = []
            for p in file_or_dir:
                unsorted_files += glob(os.path.join(p, "*" + ext))
        elif kind.lower() == "file":
            unsorted_files = file_or_dir

        timesfiles = [(self._peektime(f), f) for f in unsorted_files]
        timesfiles.sort(key=lambda p: p[0])
        times, sorted_files = zip(*timesfiles)
        buff_io = self._buffered_read(sorted_files, **kwargs)
        yield from self._chunk(buff_io, interval)

    def _readfile(
        self, fname, df_output=True, units=True, flags=True, **kwargs
    ):
        """Read a single high frequency eddy covariance data file.

        Defaults: read using format specified in the initializer; unit
        conversion is applied to eddy covariance data; flag columns are
        converted to mask booleans (True => masked).

        Parameters
        ----------
        fname : str
            Filename.
        df_output : bool
            Return dataframe if True (default), HFData object if False
        units : bool
            Apply SI unit conversion to EC data if True (default).
            Ignored if `df_output` is False.
        flags : bool
            Convert flag columns to boolean mask (True => masked).
            Ignored if `df_output` is False.
        **kwargs
            For csv filetype, kwargs are passed to pandas.read_csv_.
            Can be used to override or add to formatting specified in
            the initializer. Should not include `usecols`  or `header`
            keywords. For tob1 filetype, kwargs is passed to
            numpy.fromfile_. In this case, the only allowable kwarg is
            'count', which can be used to limit the number of lines
            read.


        .. _numpy.fromfile
            https://docs.scipy.org/doc/numpy/reference/generated/numpy.fromfile.html

        .. _pandas.read_csv:
            https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

        """
        try:
            if self._filetype == "csv":
                dataframe = self._read_csv(fname, **kwargs)
            elif self._filetype in ("tob", "tob1"):
                dataframe = self._read_tob1(fname, **kwargs)
            elif self._filetype == "pd.df":
                dataframe = self._read_df(fname)
            else:
                raise HFDataReadError("Unknown file type")
        except Exception as err:
            raise HFDataReadError(err.args[0])

        if units or not df_output:
            dataframe = self._unit_convert(dataframe)
        if flags or not df_output:
            for col, val in self._flags:
                dataframe.loc[:, "flag-" + str(col)] = (
                    dataframe.loc[:, "flag-" + str(col)] != val
                )
        if df_output:
            return dataframe
        return HFData(dataframe)

    def _read_csv(self, csvfile, **kwargs):
        df = pd.read_csv(
            csvfile,
            usecols=self._usecols,
            header=None,
            **{**self._csvformat_kws, **kwargs},
        )
        df.columns = self._names
        if self._time_col is not None:
            df["Datetime"] = pd.to_datetime(df.iloc[:, self._time_col])
            df = df.set_index("Datetime")
        return df

    def _read_df(self, df):
        df = df.iloc[:, self._usecols]
        df.columns = self._names
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

    def _unit_convert(self, df):
        for var, func in self._converters.items():
            df.loc[:, var] = func(df.loc[:, var])
        return df

    @property
    def _usecols(self):
        return [self._namecols[k] for k in self._names]


@attr.s
class HFSummary(object):
    """Summary of high frequency eddy covariance data.

    Parameters
    ----------
    T : float
        Mean air temperature, K.
    P, Pvap : float
        Mean total atmospheric (`P`) and vapor (`Pvap`) pressure, Pa.
    ustar : float
        Mean friction velocity, m/s.
    wind_w : float
        Mean vertical wind velocity, m/s.
    var_w : float
        Variance of vertical wind velocity, (m/s)^2.
    rho_vapor, rho_co2 : float
        Mean H2O vapor and CO2 concentrations, kg/m^3.
    var_vapor, var_co2 : float
        Variance of H2O vapor and CO2 concentrations, (kg/m^3)^2.
    corr_q_c : float
        Correlation coefficient for H2O and CO2 concentrations
    cov_w_q, cov_w_c : float
        Covariance of vertical wind velocity (w) with water vapor (q)
        and CO2 (c) mass densities, kg/m^2/s.
    H : float
        Sensible heat flux, W/m^2.
    rho_dryair, rho_totair : float
        Dry and moist air densities, kg/m^3.
    cov_w_T : float
        Covariance of temperature and vertical wind velocity, K m/s.
    N : int
        Length of data series.

    """

    T = attr.ib(default=np.nan)
    P = attr.ib(default=np.nan)
    Pvap = attr.ib(default=np.nan)
    ustar = attr.ib(default=np.nan)
    wind_w = attr.ib(default=np.nan)
    var_w = attr.ib(default=np.nan)
    rho_vapor = attr.ib(default=np.nan)
    rho_co2 = attr.ib(default=np.nan)
    var_vapor = attr.ib(default=np.nan)
    var_co2 = attr.ib(default=np.nan)
    corr_q_c = attr.ib(default=np.nan)
    cov_w_q = attr.ib(default=np.nan)
    cov_w_c = attr.ib(default=np.nan)
    H = attr.ib(default=np.nan)
    rho_dryair = attr.ib(default=np.nan)
    rho_totair = attr.ib(default=np.nan)
    cov_w_T = attr.ib(default=np.nan)
    N = attr.ib(default=np.nan)

    def __str__(self):

        # For some fields, print common units instead of SI
        T = self.T - 273.15
        P = 1e-3 * self.P
        Pvap = 1e-3 * self.Pvap
        rho_vapor = 1e3 * self.rho_vapor
        rho_co2 = 1e6 * self.rho_co2
        var_vapor = 1e6 * self.var_vapor
        var_co2 = 1e12 * self.var_co2
        cov_w_q = 1e3 * self.cov_w_q
        cov_w_c = 1e6 * self.cov_w_c

        return (
            "---------------\n"
            "HF Data Summary\n"
            "---------------\n"
            f"T = {T:.4} C\n"
            f"P = {P:.4} kPa\n"
            f"Pvap = {Pvap:.4} kPa\n"
            f"ustar = {self.ustar:.4} m/s\n"
            f"wind_w = {self.wind_w:.4} m/s\n"
            f"var_w = {self.var_w:.4} (m/s)^2\n"
            f"rho_vapor = {rho_vapor:.4} g/m^3\n"
            f"rho_co2 = {rho_co2:.4} mg/m^3\n"
            f"var_vapor = {var_vapor:.4} (g/m^3)^2\n"
            f"var_co2 = {var_co2:.4} (mg/m^3)^2\n"
            f"corr_q_c = {self.corr_q_c:.4}\n"
            f"cov_w_q = {cov_w_q:.4} g/m^2/s\n"
            f"H = {self.H:.4} W/m^2\n"
            f"cov_w_c = {cov_w_c:.4} mg/m^2/s\n"
            f"rho_dryair = {self.rho_dryair:.4} kg/m^3\n"
            f"rho_totair = {self.rho_totair:.4} kg/m^3\n"
            f"cov_w_T = {self.cov_w_T:.4} C m/s\n"
            f"N = {self.N}"
        )
