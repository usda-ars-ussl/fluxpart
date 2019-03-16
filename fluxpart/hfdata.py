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
    pass


VAR_NAMES = ["u", "v", "w", "c", "q", "T", "P"]

_badfiletype = "File type not recognized ({})"
_toofewdata = (
    "The longest contiguous run of valid data was too short:\n"
    "(length data / total length) = {frac:.4} < rd_tol = {rd_tol:.4}\n"
    "OR data length = {N} < ad_tol = {ad_tol}"
)


class HFData(object):
    """
    Parameters
    ----------
    hf_dataframe
        High frequency eddy covariance dataframe. Must include columns
        for ["u", "v", "w", "c", "q", "T", "P"]. Normally, dataframe
        should have a datetime index. Dataframe may include also boolean
        mask columns.

    """

    def __init__(self, hf_dataframe):
        self.dataframe = hf_dataframe
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
            mssg = _toofewdata.format(
                frac=data_frac, rd_tol=rd_tol, N=len_max_slice, ad_tol=ad_tol
            )
            raise TooFewDataError(mssg)

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
    files : list or sequence of files
        Sorted sequence of data files
    filetype : {'csv', 'tob1', 'ghg'}
        'csv' = delimited text file; 'tob1' = Campbell Scientific binary
        format file; 'ghg' = LI-COR raw data format 
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
        Datetime column for `csv` and `ghg` filetype. Default is None.
    flags : 2-tuple or list of 2-tuples, optional
        Indicate that one or more data columns are used to flag bad data
        records. Each tuple is of the form (col, goodval), where col is
        an int specifying the column number (0-based indexing), and
        goodval is the flag value indicating a good data record.
    **kwargs
        Passed to pandas.read_csv_ when filetype is csv or ghg. Should
        not include `usecols` or `header` keywords.


    .. _pandas.read_csv:
        https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

    """

    def __init__(
        self,
        files,
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

        self._files = files
        self._filetype = filetype.strip().lower()
        self._cols = cols
        self._converters = {} if converters is None else converters
        self._time_col = time_col
        self._flags = flags
        self._csv_kws = kwargs
        self._dt_kws = self._csv_kws.pop("to_datetime_kws", {})

    @property
    def _names(self):
        # sorted because pd.read_csv sorts usecols but not names
        return sorted(self._namecols, key=self._namecols.get)

    @property
    def _namecols(self):
        namecols = dict(zip(VAR_NAMES, self._cols))
        flags = {"flag-" + str(col): col for col, val in self._flags}
        namecols.update(flags)
        if self._filetype in ("csv", "ghg") and self._time_col is not None:
            try:
                namecols["Date"] = self._time_col[0]
                namecols["Time"] = self._time_col[1]
            except TypeError:
                namecols["Datetime"] = self._time_col
        return namecols

    def reader(self, interval, **kwargs):
        """Consume data source and yield in chunks of the time interval.

        Data chunks are returned in dataframe format.

        Parameters
        ----------
        interval : str
            Time interval used to chunk the data. Is specified using the
            pandas string alias_ format (e.g., ``interval="30min"``).
            If set to -1, a single dataframe corresponding to a
            concatenation of all data files is returned on the first
            iteration. If set to None, dataframes corresponding to whole
            individual data files are returned with each iteration.
        **kwargs
            For csv and ghg filetypes, kwargs are passed to
            pandas.read_csv_.  Can be used to override or add to
            formatting specified in the initializer. Should not include
            `usecols`  or `header` keywords. For tob1 filetype, kwargs
            is passed to numpy.fromfile_. In this case, the only
            allowable kwarg is 'count', which can be used to limit the
            number of lines read.


        .. _numpy.fromfile:
            https://docs.scipy.org/doc/numpy/reference/generated/numpy.fromfile.html

        .. _pandas.read_csv:
            https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

        .. _alias:
            http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases

        """

        if self._filetype in ("csv", "ghg"):
            kws = {**self._csv_kws, **kwargs}
            kws["usecols"] = self._usecols
            kws["header"] = None
            if self._filetype == "csv":
                dfs = util.multifile_read_csv(self._files, **kws)
            else:
                dfs = util.multifile_read_ghg(self._files, **kws)
            indx_dfs = (self._set_indices_csv(df) for df in dfs)
        else:  # "tob1"
            count = kwargs.get("count", -1)
            dfs = util.multifile_read_tob1(self._files, count=count)
            indx_dfs = (self._set_indices_tob1(df) for df in dfs)

        si_dfs = (self._set_units(df) for df in indx_dfs)
        hf_dfs = (self._set_flags(df) for df in si_dfs)
        yield from util.chunked_df(hf_dfs, interval)

    def _set_indices_csv(self, df):
        df.columns = self._names
        if self._time_col is None:
            return df
        try:
            df["Datetime"] = pd.to_datetime(
                df.loc[:, "Datetime"], **self._dt_kws
            )
        except KeyError:
            df["Datetime"] = pd.to_datetime(
                df.loc[:, "Date"] + " " + df.loc[:, "Time"], **self._dt_kws
            )
            df = df.drop(["Date", "Time"], axis=1)
        return df.set_index("Datetime")

    def _set_indices_tob1(self, df):
        df["Datetime"] = pd.to_datetime(
            arg=df.loc[:, "SECONDS"] + 10 ** -9 * df.loc[:, "NANOSECONDS"],
            unit="s",
            origin="1990-01-01",
        )
        df = df.set_index("Datetime")
        df = df.iloc[:, self._usecols]
        df.columns = self._names
        return df

    def _set_flags(self, df):
        for col, val in self._flags:
            df.loc[:, "flag-" + str(col)] = (
                df.loc[:, "flag-" + str(col)] != val
            )
        return df

    def _set_units(self, df):
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
        # prints common units instead of SI
        return self.results_str().format(**self.common_units())

    def results_str(self):
        lab = self.common_units_labels()
        return (
            "---------------\n"
            "HF Data Summary\n"
            "---------------\n"
            "  T = {T:.4} " + lab["T"] + "\n"
            "  P = {P:.4} " + lab["P"] + "\n"
            "  Pvap = {Pvap:.4} " + lab["Pvap"] + "\n"
            "  ustar = {ustar:.4} " + lab["ustar"] + "\n"
            "  wind_w = {wind_w:.4} " + lab["wind_w"] + "\n"
            "  var_w = {var_w:.4} " + lab["var_w"] + "\n"
            "  rho_vapor = {rho_vapor:.4} " + lab["rho_vapor"] + "\n"
            "  rho_co2 = {rho_co2:.4} " + lab["rho_co2"] + "\n"
            "  var_vapor = {var_vapor:.4} " + lab["var_vapor"] + "\n"
            "  var_co2 = {var_co2:.4} " + lab["var_co2"] + "\n"
            "  corr_q_c = {corr_q_c:.4} " + lab["corr_q_c"] + "\n"
            "  cov_w_q = {cov_w_q:.4} " + lab["cov_w_q"] + "\n"
            "  H = {H:.4} " + lab["H"] + "\n"
            "  cov_w_c = {cov_w_c:.4} " + lab["cov_w_c"] + "\n"
            "  rho_dryair = {rho_dryair:.4} " + lab["rho_dryair"] + "\n"
            "  rho_totair = {rho_totair:.4} " + lab["rho_totair"] + "\n"
            "  cov_w_T = {cov_w_T:.4} " + lab["cov_w_T"] + "\n"
            "  N = {N:.0f} " + lab["N"]
        )

    def common_units(self):
        return dict(
            T=self.T - 273.15,
            P=1e-3 * self.P,
            Pvap=1e-3 * self.Pvap,
            ustar=self.ustar,
            wind_w=self.wind_w,
            var_w=self.var_w,
            rho_vapor=1e3 * self.rho_vapor,
            rho_co2=1e6 * self.rho_co2,
            var_vapor=1e6 * self.var_vapor,
            var_co2=1e12 * self.var_co2,
            corr_q_c=self.corr_q_c,
            cov_w_q=1e3 * self.cov_w_q,
            H=self.H,
            cov_w_c=1e6 * self.cov_w_c,
            rho_dryair=self.rho_dryair,
            rho_totair=self.rho_totair,
            cov_w_T=self.cov_w_T,
            N=self.N,
        )

    def common_units_labels(self):
        return dict(
            T="C",
            P="kPa",
            Pvap="kPa",
            ustar="m/s",
            wind_w="m/s",
            var_w="(m/s)^2",
            rho_vapor="g/m^3",
            rho_co2="mg/m^3",
            var_vapor="(g/m^3)^2",
            var_co2="(mg/m^3)^2",
            corr_q_c="",
            cov_w_q="g/m^2/s",
            H="W/m^2",
            cov_w_c="mg/m^2/s",
            rho_dryair="kg/m^3",
            rho_totair="kg/m^3",
            cov_w_T="C m/s",
            N="",
        )
