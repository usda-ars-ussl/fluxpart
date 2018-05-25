"""Data structure for high-frequency eddy covariance time series."""

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
            'HF Data read but rejected because the longest continuous '
            'run of valid data was too short on a relative (length '
            'data / total length) = {:.4} < dtol = {:.4}) AND/OR '
            'absolute basis (data length = {} < {})'
            ''.format(data_frac, rd_tol, len_max_slice, ad_tol))


class HFData(object):
    """High-frequency eddy covariance data.

    The following notation is used in variable naming and
    documentation to represent meteorological quantities (SI
    units)::

        u, v, w = wind velocities (m/s)
        q = water vapor mass concentration (kg/m^3)
        c = carbon dioxide mass concentration (kg/m^3)
        T = air temperature (K)
        P = total air pressure (Pa)

    """
    var_names = ['u', 'v', 'w', 'c', 'q', 'T', 'P']

    def __init__(
            self,
            datasource='csv',
            cols=(2, 3, 4, 5, 6, 7, 8),
            time_col=None,
            converters=None,
            flags=None,
            **kwargs):
        """High-frequency eddy covariance data.

        Parameters
        ----------
        datasource : {'csv', 'tob1', 'pd.df'}
            'csv' = delimited text file (default); 'tob1' = Campbell
            Scientific binary format file; 'pd.df' = pandas dataframe.
        cols : 7*(int,), optional
            Column indices for (u, v, w, q, c, T, P) data, in that
            order. 0-based indexing. Default is (2, 3, 4, 6, 5, 7, 8).
        time_col : int, optional (TODO)
            (TODO) Datetime column. Default is None.
        converters : dict, optional
            Dictionary of functions used to convert any non-SI data to
            SI units.  Dict keys are 'u', 'v', 'w', 'q', 'c', 'T', or
            'P'. With a 'csv' `datasource`, the funcs take a string
            single argument, e.g.
            ``converters={'P': lambda s: 1e3 * float(s)}``.
        flags : 2-tuple or list of 2-tuples, optional
            Specifies that one or more columns in `fname` are used to
            flag bad data records. Each tuple is of the form (col,
            goodval), where col is an int specifying the column number
            containing the flag (0-based indexing), and goodval is the
            value of the flag that indicates a good data record.
        **kwargs
            With a 'csv' `datasource`, kwargs is passed to
            pandas.read_csv_. Can be used to indicate additional
            formatting options. Should *not* include the following keys:
            `usecols`, `names`, `converters`, `dtype`, `index_col`.


        .. _pandas.read_csv:
            https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

        Attributes
        ----------
        dataframe : pandas dataframe
            Table of high frequency data series with columns 'u', 'v',
            'w', 'q', 'c', 'T', and 'P'

        """
        self._datasource = datasource
        self._cols = cols
        self._time_col = time_col
        self._converters = converters or {}
        self._flags = flags or {}
        self._read_kws = kwargs
        self._already_corrected_external = False
        self.dataframe = None

    @property
    def _dtype(self):
        dtype = dict(zip(HFData.var_names, len(HFData.var_names) * (float, )))
        dtype.update({k: type(v[1]) for k, v in self._flags.items()})
        # Avoid pandas Warning about vars being in both converters and dtype
        for k in self._converters.keys():
            dtype.pop(k, None)
        return dtype

    @property
    def _namecols(self):
        namecols = dict(zip(HFData.var_names, self._cols))
        namecols.update({k: v[0] for k, v in self._flags.items()})
        return namecols

    @property
    def _names(self):
        # sorted because pd.read_csv sorts usecols but not names
        return sorted(self._namecols, key=self._namecols.get)

    @property
    def _usecols(self):
        return [self._namecols[k] for k in self._names]

    def __getitem__(self, name):
        """Column-wise get without specifying dataframe attribute"""
        return self.dataframe.loc[:, name]

    def __setitem__(self, name, value):
        """Column-wise set without specifying dataframe attribute"""
        self.dataframe.loc[:, name] = value

    def read(self, fname, **kwargs):
        """Read high frequency eddy covariance data into dataframe.

        Parameters
        ----------
        fname : str
            Filepath/filename or dataframe.
        **kwargs
            With a 'csv' `datasource`, kwargs is passed to
            pandas.read_csv_. Will be merged with/overwrite any kwargs
            passed in the constructor. Should *not* include the
            following keys: `usecols`, `names`, `converters`, `dtype`,
            `index_col`.

        """
        try:
            if self._datasource.strip().lower() == 'csv':
                self._read_csv(fname, **kwargs)
            elif self._datasource.strip().lower() in ('tob', 'tob1'):
                self._read_tob1(fname)
            elif self._datasource.strip().lower() == 'pd.df':
                self.dataframe = fname
                self._format_df()
            else:
                raise HFDataReadError('Unknown file type')
        except Exception as err:
            raise HFDataReadError(err.args[0])

        self._already_corrected_external = False

    def _read_csv(self, fname, **kwargs):
        kws = dict(
            usecols=self._usecols,
            # TODO
            # index_col=self._time_col,
            dtype=self._dtype,
            names=self._names,
            converters=self._converters,
            **{**self._read_kws, **kwargs}
        )
        self.dataframe = pd.read_csv(fname, **kws)

    def _read_tob1(self, tobfile):
        self.dataframe = pd.DataFrame(util.tob1_to_array(tobfile))
        self._format_df()

    def _format_df(self):
        self.dataframe = self.dataframe.iloc[:, self._usecols]
        self.dataframe.columns = self._names
        for k, func in self._converters.items():
            self.dataframe.loc[:, k] = func(self.dataframe.loc[:, k])

    def quality_check(self, bounds=None, rd_tol=0.5, ad_tol=1024):
        """Apply some data QC/QA.

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
        data = self.dataframe

        # 1D bool mask is True if any `data` field is missing (=nan)
        mask = data.loc[:, HFData.var_names].isnull().any(axis=1).values

        # Also mask records with out-of-bounds or flagged data
        bounds = bounds or {}
        for var, (low, high) in bounds.items():
            mask[(data[var] < low) | (data[var] > high)] = True
        if self._flags:
            for flgname, val in self._flags.items():
                mask[data[flgname] != val[1]] = True

        # Find longest span of valid (unmasked) data
        dummy_ma = np.ma.array(np.zeros([data.shape[0], ]), mask=mask)
        unmasked_slices = np.ma.clump_unmasked(dummy_ma) or [slice(0, 0), ]
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

    def truncate(self):
        """Truncate dataframe length to largest possible power of 2."""
        truncate_len = 2 ** int(np.log2(self.dataframe.shape[0]))
        self.dataframe = self.dataframe.iloc[:truncate_len]

    def correct_external(self):
        """Adjust q and c data series to correct for external effects.

        Water vapor and carbon dioxide series data in the dataframe are
        corrected for external fluctuations associated with air
        temperature and vapor density. See: [WPL80]_ and [DK07]_.

        """

        if self._already_corrected_external:
            return
        ave_vapor = self['q'].mean()
        ave_co2 = self['c'].mean()
        ave_T = self['T'].mean()
        dev_vapor = self['q'] - ave_vapor
        dev_T = self['T'] - ave_T

        Pdryair = self['P'].mean() - ave_vapor * GC.vapor * ave_T
        rho_totair = ave_vapor + Pdryair / GC.dryair / ave_T

        specific_vapor = ave_vapor / rho_totair
        specific_co2 = ave_co2 / rho_totair
        mu = MW.dryair / MW.vapor
        muq = mu * specific_vapor
        muc = mu * specific_co2

        self['q'] += muq * dev_vapor + (1 + muq) * ave_vapor * dev_T / ave_T
        self['c'] += muc * dev_vapor + (1 + muq) * ave_co2 * dev_T / ave_T
        self._already_corrected_external = True
        return

    def summarize(self):
        """Summarize high frequency dataframe statistics.

        Returns
        -------
        namedtuple
            :class:`~fluxpart.containers.HFSummary`

        """
        hfs = util.stats2(self.dataframe, HFData.var_names)
        Pvap = hfs.ave_q * GC.vapor * hfs.ave_T
        rho_dryair = (hfs.ave_P - Pvap) / GC.dryair / hfs.ave_T
        rho_totair = rho_dryair + hfs.ave_q
        Cp = CP.dryair * (1 + 0.84 * hfs.ave_q / rho_totair)

        return HFSummary(
            T=hfs.ave_T,
            P=hfs.ave_P,
            Pvap=Pvap,
            ustar=(hfs.cov_w_u**2 + hfs.cov_w_v**2)**0.25,
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
