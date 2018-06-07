from copy import deepcopy

import attr
import numpy as np

from .__version__ import __version__
from .wue import water_use_efficiency, WUEError
from .hfdata import get_hfdata, HFDataReadError, TooFewDataError
from .partition import fvspart_progressive
from .util import vapor_press_deficit
from .containers import AllFluxes, WUE


EC_TOA5 = {
    "filetype": "csv",
    "skiprows": 4,
    "cols": (2, 3, 4, 5, 6, 7, 8),
    "time_col": 0,
    "temper_unit": 'C',
    "unit_convert": dict(q=1e-3, c=1e-6, P=1e3),
}

EC_TOB1 = {
    "filetype": "tob1",
    "cols": (3, 4, 5, 6, 7, 8, 9),
    "temper_unit": "C",
    "unit_convert": dict(q=1e-3, c=1e-6, P=1e3),
}

HFD_FORMAT = EC_TOB1

WUE_OPTIONS = {
    "ci_mod": "const_ratio",
    "ci_mod_param": None,
    "leaf_temper": None,
    "leaf_temper_corr": 0,
    "diff_ratio": 1.6,
}

HFD_OPTIONS = {
    "bounds": {"c": (0, np.inf), "q": (0, np.inf)},
    "rd_tol": 0.4,
    "ad_tol": 1024,
    "ustar_tol": 0.1,
    "correcting_external": True,
}

PART_OPTIONS = dict(adjust_fluxes=True, sun=None, interval=None)


_bad_ustar = "ustar = {:.4} is less than ustar_tol = {:.4}"
_bad_vpd = "Vapor pressure deficit({} <= 0"
_bad_qflux = "Fq = cov(w,q) = {:.4} <= 0 is incompatible with fvs partitioning"


def fvs_partition(
    file_or_dir,
    hfd_format=None,
    hfd_options=None,
    meas_wue=None,
    wue_options=None,
    part_options=None,
    label=None,
):
    """Partition CO2 & H2O fluxes into stomatal & nonstomatal components.

    Provides a full implementation of the flux variance similarity
    partitioning algorithm [SS08]_[SAAS+18]_: reads high frequency eddy
    covariance data; performs data transformations and data QA/QC;
    analyzes water vapor and carbon dioxide fluxes; and partitions the
    fluxes into stomatal (transpiration, photosynthesis) and nonstomatal
    (evaporation, respiration) components.

    The following notation is used in variable naming and documentation
    to represent meteorological quantities::

        u, v, w = wind velocities
        q = water vapor mass concentration
        c = carbon dioxide mass concentration
        T = air temperature
        P = total air pressure

    Parameters
    ----------
    file_or_dir : str or list of str
        High frequency eddy covariance data file, data directory, or
        list of either. If a directory or list of directories is passed,
        data files within the directory(s) will be analyzed.
    hfd_format : dict or {"ec-TOB1", "ec-TOA5"}, optional
        Dictionary of parameters specifying the high frequency data file
        format. See ``Other parameters`` below for an explanation of
        required and optional formatting parameters. Some pre-defined
        file formats can be specified by passing a string instead of a
        dictionary.  Currently two formats are defined: "ec-TOB1" and
        "ec-TOA5".  These formats correspond to files commonly created
        when using Campbell Scientific eddy covariance data logger
        software. See ``Notes`` below for an explanation of these
        pre-defined formats. Default is "ec-TOB1".
    hfd_options : dict, optional
        Dictionary of parameters specifying options for correcting
        high-frequency eddy covariance data and applying quality control
        measures. See ``Other parameters`` for a listing of options.
    meas_wue : float or callable, optional
        Measured (or otherwise prescribed) leaf-level water use
        efficiency (kg CO2 / kg H2O). Note that by definition,
        `meas_wue` must be a negative value (< 0). If callable, it must
        take a pandas datetime object as its sole argument and return a
        float value.
    wue_options : dict, required if `meas_wue` is not provided
        Dictionary of parameters and options used to estimate water use
        efficiency if `meas_wue` is not provided. See ``Other
        parameters`` section for a description of valid fields.  When
        passing `wue_options`, it is always required to provide entries
        'canopy_ht', 'meas_ht', and 'ppath'. Other entries are optional.
        The values for 'canopy_ht' and 'meas_ht' can be either a float
        or a callable. The callable must accept a datetime object as its
        sole argument and return a float value.
    part_options : dict, optional
        Dictionary of options for the fvs partitioning procedure. See
        ``Other parameters`` section for a listing of valid options.
    label : optional
        Optional identifier to be appended to the results.


    Returns
    -------
    :class:`~fluxpart.fluxpart.FluxpartResult` or
    :class:`~fluxpart.fluxpart.FluxpartSeriesResult`


    Other Parameters
    ----------------
    hfd_format : dict
        High frequency data file format. NOTE: When passing hfd_format,
        it is required at a minimum to provide values for 'filetype' and
        'cols' (detailed below).  'unit_convert' and 'temper_unit' are
        also required if data are not in SI units.
    hfd_format['filetype'] : {'csv', 'tob1', 'pd.df'}
        'csv' = delimited text file; 'tob1' = Campbell Scientific binary
        data table format; 'pd.df' = pandas dataframe.
    hfd_format['cols'] : 7*(int,)
        7-tuple of integers indicating the data column numbers that
        contain series data for (u, v, w, c, q, T, P), in that order.
        Uses 0-based indexing.
    hfd_format['unit_convert'] : dict
        Dictionary of multiplication factors required to convert any u,
        v, w, c, q, or P data not in SI units to SI units (m/s, kg/m^3,
        Pa). (Note T is not in that list). The dictionary keys are the
        variable names. For example, if all data are in SI units except
        P and c, which are in kPa and mg/m^3, respectively, then set:
        ``hfd_options['unit_convert'] = {'P': 1e3, 'c': 1e-6}``,
        since it is necessary to multiply the kPa pressure data by 1e3
        to obtain the SI pressure unit (Pa), and the mg/m^3 CO2 data by
        1e-6 to obtain the SI concentration unit (kg/m^3).
    hfd_format['temper_unit'] : {'K', 'C'}
        Temperature data units. Default is 'K' (Kelvin).
    hfd_format['ext'] : str, optional
        When reading multiple files from directory, all files in the
        directory with extension `ext` will be read. Default ("") reads
        all files regardless of extension. When specifying `ext`,
        include the "dot" where appropriate (e.g.,
        ``hfd_format["ext"] = ".dat"``)
    hfd_format['flags'] : 2-tuple or list of 2-tuples
        Specifies that one or more data columns are used to flag bad
        data records. Each tuple is of the form (col, goodval), where
        col is an int specifying the column number containing the flag
        (0-based indexing), and goodval is the value of the flag that
        indicates a good data record.
    hfd_format[ other keys ]
        when `hfd_format['filetype']` is 'csv', all other key:value
        pairs in `hfd_format` are passed as keyword arguments to
        pandas.read_csv_. Those keywords may be required to specify the
        details of the file formatting. Among the most commonly required
        are: 'sep', the str that is used to separate values or define
        column widths (default is sep=','); and 'skiprows', which will
        be needed if the file contains header rows. See pandas.read_csv_
        for a full description of available format options.

    hfd_options: dict
        Options for pre-processing high frequency data.
    hfd_options['correcting_external'] : bool, optional
        If True (default), the water vapor and carbon dioxide series
        data are corrected for external fluctuations associated with air
        temperature and vapor density according to [WPL80]_ and [DK07]_.
    hfd_options['bounds'] : dict, optional
        Dictionary indicating any lower and upper bounds for valid data.
        Dictionary entries have the form ``varname: (float, float)``,
        where varname is one of 'u', 'v', 'w', 'q', 'c', 'T', or 'P',
        and the 2-tuple holds values for the lower and upper bounds:
        ``(lower, upper)``.  Data records are rejected if a variable in
        the record is outside the prescribed bounds. Default is ``bounds
        = {'c': (0, np.inf), 'q': (0, np.inf)}`` such that data records
        are rejected if c or q data are not positive values.
    hfd_options['rd_tol'] : float, optional
        Relative tolerance for rejecting the datafile. Default is
        'hfd_options['rd_tol']` = 0.4. See hfd_options['ad_tol'] for
        further explanation.
    hfd_options['ad_tol'] : int, optional
        Absolute tolerance for rejecting the datafile. Default is
        `ad_tol` = 1024. If the datafile contains bad records (not
        readable, out-of-bounds, or flagged data), the partitioning
        analysis is performed using the longest stretch of consecutive
        good data records found, unless that stretch is too short, in
        which case the analysis is aborted. The criteria for judging
        'too short' can be specified in both relative and absolute
        terms: the datafile is rejected if the good stretch is a
        fraction of the total data that is less than `rd_tol`, and/or is
        less than `ad_tol` records long.
    hfd_options['ustar_tol'] : float
        If the friction velocity (m/s) determined from the high
        frequency data is less than `ustar_tol`, the
        partitioning analysis is aborted due to insufficient turbulence.
        Defalult is `hfd_options['ustar_tol']` = 0.1 (m/s).

    wue_options: dict
        Parameters for estimating water use efficiency. 'canopy_ht',
        'meas_ht', and 'ppath' are required keys.
    wue_options['canopy_ht'] : float or callable
        Vegetation canopy height (m). If callable must accept a datetime
        object as its sole argument and return a float value.
    wue_options['meas_ht'] : float or callable
        Eddy covariance measurement height (m). If callable must accept
        a datetime object as its sole argument and return a float value.
    wue_options['ppath'] : {'C3', 'C4'}
        Photosynthetic pathway.
    wue_options['ci_mod'] : str
        Valid values: 'const_ratio', 'const_ppm', 'linear', 'sqrt'.
        See: :func:`~fluxpart.wue.water_use_efficiency`.
    wue_options['ci_mod_param'] : float or (float, float)
        Paramter values to be used with `ci_mod`.
        See: :func:`~fluxpart.wue.water_use_efficiency`.
    wue_options['leaf_temper'] : float
        Canopy leaf temperature (K). If not specified, it is assumed to
        be equal to the air temperature See:
        :func:`~fluxpart.wue.water_use_efficiency`.
    wue_options['leaf_temper_corr'] : float
        Offset adjustment applied to canopy temperature (K). Default is
        zero. See: :func:`~fluxpart.wue.water_use_efficiency`.
    wue_options['diff_ratio']: float, optional
        Ratio of molecular diffusivities for water vapor and CO2.
        Default is `diff_ratio` = 1.6.
        See: :func:`~fluxpart.wue.water_use_efficiency`.

    part_options : dict
        Options for the fvs partitioning algorithm
    part_options['interval'] : str, optional
        Time interval to be used to aggregate data and partition fluxes.
        The interval is specified using the pandas string alias_ format.
        Set to None to process an entire data file without any
        consideration of time or  date. Default is None.
    part_options['adjust_fluxes'] : bool
        If True (default), the final partitioned fluxes are adjusted
        proportionally such that sum of the partitioned fluxes match
        exactly the total fluxes indicated in the original data.
    part_options['sun'] : 2-tuple = (sunrise, sunset), or callable
        A 2-tuple of times corresponding to sunrise and sunset. If
        specified, fluxes for times starting  before sunrise or after
        sunset will be taken to be all non-stomatal. If callable,
        it should take a datetime as its sole argument and return the
        2-tuple of times. Default is (None, None).


    NOTES
    -----
    Two pre-defined hfd_formats are available, 'ec.TOA5' and 'ec.TOB1'.

    'ec.TOA5'::

        hfd_format = {
            "filetype": 'csv',
            "skiprows": 4,
            "cols": (2, 3, 4, 5, 6, 7, 8),
            "time_col": 0,
            "temper_unit": 'C',
            "unit_convert": {"q": 1e-3, "c": 1e-6, "P": 1e3},

        }

    'ec.TOB1'::

        hfd_format = {
            "filetype": 'tob1',
            'cols': (3, 4, 5, 6, 7, 8, 9),
            "temper_unit": 'C',
            "unit_convert": {"q": 1e-3, "c": 1e-6, "P": 1e3},
        }


    .. _pandas.read_csv:
        https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

    .. _alias:
        http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases

    """
    if hfd_format is None:
        hfd_format = deepcopy(HFD_FORMAT)
    elif hfd_format.upper() == "EC-TOA5":
        hfd_format = deepcopy(EC_TOA5)
    elif hfd_format.upper() == "EC-TOB1":
        hfd_format = deepcopy(EC_TOB1)
    else:
        # _validate_hfd_format(hfd_format)
        pass

    hfd_options = {**HFD_OPTIONS, **(hfd_options or {})}
    wue_options = {**WUE_OPTIONS, **(wue_options or {})}
    part_options = {**PART_OPTIONS, **(part_options or {})}

    usecols = np.array(hfd_format.pop("cols"), dtype=int).reshape(7)
    unit_convert = hfd_format.pop("unit_convert", {})
    converters = {k: _converter_func(v, 0.) for k, v in unit_convert.items()}
    temper_unit = hfd_format.pop("temper_unit")
    if temper_unit.upper() == "C" or temper_unit.upper() == "CELSIUS":
        converters["T"] = _converter_func(1., 273.15)

    try:
        hfdat = get_hfdata(
            file_or_dir, cols=usecols, converters=converters, **hfd_format
        )
    except HFDataReadError as err:
        return FluxpartResult(dataread=False, mssg=err.args[0])

    # Early returns if data are not compatible with fvs partitioning
    try:
        hfdat.cleanse(
            hfd_options["bounds"], hfd_options["rd_tol"], hfd_options["ad_tol"]
        )
    except TooFewDataError as err:
        return FluxpartResult(dataread=True, mssg=err.args[0])

    hfdat.truncate_pow2()
    if hfd_options["correcting_external"]:
        hfdat.correct_external()
    hfsum = hfdat.summarize()

    if hfsum.ustar < hfd_options["ustar_tol"]:
        mssg = _bad_ustar.format(hfsum.ustar, hfd_options["ustar_tol"])
        return FluxpartResult(dataread=True, mssg=mssg, hfsummary=hfsum)

    vpd = vapor_press_deficit(hfsum.rho_vapor, hfsum.T)
    if vpd <= 0:
        return FluxpartResult(dataread=True, mssg=_bad_vpd.format(vpd))

    if hfsum.cov_w_q <= 0:
        mssg = _bad_qflux.format(hfsum.cov_w_q)
        return FluxpartResult(dataread=True, mssg=mssg, hfsummary=hfsum)

    if meas_wue:
        leaf_wue = WUE(wue=float(meas_wue))
    else:
        try:
            leaf_wue = water_use_efficiency(hfsum, **wue_options)
        except WUEError as err:
            return FluxpartResult(
                dataread=True, mssg=err.args[0], hfsummary=hfsum
            )

    # Everything seems OK so compute partitioned fluxes
    fvsp = fvspart_progressive(
        hfdat["w"].values,
        hfdat["q"].values,
        hfdat["c"].values,
        leaf_wue.wue,
        part_options["adjust_fluxes"],
    )

    if fvsp.valid_partition:
        fvsp.fluxes = AllFluxes(
            **attr.asdict(fvsp.fluxes), temper_kelvin=hfsum.T
        )
    else:
        fvsp.fluxes = None

    return FluxpartResult(
        label=label,
        dataread=True,
        attempt_partition=True,
        valid_partition=fvsp.valid_partition,
        mssg=fvsp.mssg,
        fvsp_result=fvsp,
        hfsummary=hfsum,
        wue=leaf_wue,
    )


def flux_partition(*args, **kws):
    return fvs_partition(*args, **kws)


class FluxpartResult(object):
    """Fluxpart result."""

    def __init__(
        self,
        dataread=False,
        attempt_partition=False,
        valid_partition=False,
        mssg=None,
        label=None,
        hfsummary=None,
        wue=None,
        fvsp_result=None,
    ):
        """Fluxpart result.

        Parameters
        ----------
        dataread, attempt_partition, valid_partition : bool
            Indicates success or failure in reading high frequency data,
            attempting and obtaining a valid partioning solution.
        mssg : str
            Possibly informative message if `dataread` or `valid_partition`
            are False
        label : optional
            Optional id label. Could be, e.g. a datetime object or string.
        fvsp_result : :class:`~fluxpart.containers.FVSPResult`
        wue : :class:`~fluxpart.containers.WUE`
        hfsummary : :class:`~fluxpart.containers.HFSummary`

        """
        self.version = __version__
        self.dataread = dataread
        self.attempt_partition = attempt_partition
        self.valid_partition = valid_partition
        self.mssg = mssg
        self.label = label
        self.fvsp_result = fvsp_result
        self.wue = wue
        self.hfsummary = hfsummary

    def __str__(self):
        result = (
            "===============\n"
            "Fluxpart Result\n"
            "===============\n"
            f"version = {self.version}\n"
            f"dataread = {self.dataread}\n"
            f"attempt_partition = {self.attempt_partition}\n"
            f"valid_partition = {self.valid_partition}\n"
            f"mssg = {self.mssg}\n"
        )
        if self.fvsp_result is not None:
            result += self.fvsp_result.__str__() + "\n"
        if self.wue is not None:
            result += self.wue.__str__() + "\n"
        if self.hfsummary is not None:
            result += self.hfsummary.__str__()
        return result


def _converter_func(slope, intercept):
    """Return a function for linear transform of data."""

    def func(val):
        return slope * val + intercept

    return func


if __name__ == "__main__":
    pass
"""
def _validate_hfd_format(hfd_format):
    try:
        cols = hfd_format["cols"]
    except KeyError:
        raise Error("No value for hfd_format['cols'] given.")
    if "filetype" not in hfd_format:
        raise Error("No value for hfd_format['filetype'] given.")
    if hfd_format["filetype"] not in ("csv", "tob1", "pd.df"):
        raise Error(f"Unrecognized filetype: {filetype}")
"""


