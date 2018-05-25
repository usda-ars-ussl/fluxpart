import attr
import numpy as np

import fluxpart
from .wue import water_use_efficiency, WUEError
from .hfdata import HFData, HFDataReadError
from .partition import fvspart_progressive
from .util import vapor_press_deficit
from .containers import (
    AllFluxes,
    FVSPResult,
    HFSummary,
    RootSoln,
    WQCData,
    WUE,
)

DEFAULT_WUE_OPTIONS = dict(
    ci_mod='const_ratio',
    ci_mod_param=None,
    leaf_temper=None,
    leaf_temper_corr=0,
    diff_ratio=1.6,
    )

DEFAULT_HFD_OPTIONS = dict(
    cols=(2, 3, 4, 5, 6, 7, 8),
    skiprows=4,
    sep=',',
    unit_convert={'q': 1e-3, 'c': 1e-6, 'P': 1e3},
    temper_unit='C',
    bounds={'c': (0, np.inf), 'q': (0, np.inf)},
    flags=None,
    rd_tol=0.4,
    ad_tol=1024,
    ustar_tol=0.1,
    correcting_external=True,
    )

DEFAULT_PART_OPTIONS = dict(
    adjust_fluxes=True,
    )

VERSION = fluxpart.__version__


def flux_partition(
        fname,
        meas_wue=None,
        hfd_options=None,
        wue_options=None,
        part_options=None,
        label=None,
):
    """Partition CO2 & H2O fluxes into stomatal & nonstomatal components.

    This is the primary user interface for :mod:`Fluxpart`. The function
    provides a full implementation of the flux partitioning algorithm: it
    reads high frequency eddy covariance data, performs necessary data
    transformations and QA/QC, analyzes water vapor and carbon dioxide
    fluxes, and partitions the fluxes into stomatal (transpiration,
    photosynthesis) and nonstomatal (evaporation, respiration)
    components using the flux variance similarity method of [SS08]_.

    The fluxpart submodule is imported in __init__ so this function can
    imported without referencing the submodule:
    ``import fluxpart.flux_partition``.

    The following notation is used in variable naming to represent
    meteorological quantities::

        u, v, w = wind velocities
        q = water vapor mass concentration
        c = carbon dioxide mass concentration
        T = air temperature
        P = total air pressure

    Parameters
    ----------
    fname : str
        Name of delimited file containing high-frequency eddy covariance
        time series data.
    hfd_options : dict, optional
        Dictionary of parameters specifying options for reading, quality
        control, and correcting high-frequency eddy covariance data. See
        ``Other parameters`` section for a listing of valid
        `hfd_options` fields. Note that `hfd_options` is optional only
        if the high-frequency data being read are all in SI units and
        the file is formatted according to the default options; see
        `hfd_options['unit_convert']` and `hfd_options['temper_unit']`
        for information about specifying and converting data units.
    meas_wue : float, optional
        Measured (or otherwise prescribed) leaf-level water use
        efficiency (kg CO2 / kg H2O). Note that by definition,
        `meas_wue` must be a negative value (< 0).
    wue_options : dict, required if `meas_wue` is not provided
        Dictionary of parameters and options used to estimate water use
        efficiency if `meas_wue` is not provided. See ``Other
        parameters`` section for a description of valid fields for
        `wue_options`. If specifying `wue_options`, it is always
        required to provide values for the 'canopy_ht', 'meas_ht', and
        'ppath' fields. Other entries are optional.
    part_options : dict, optional
        Dictionary of options for the flux partitioning procedure. See
        ``Other parameters`` section for a listing of valid options.
    label : optional
        Optional identification label/object for the data set. Could be
        a str, int, datetime object, etc.

    Returns
    -------
    :class:`~fluxpart.containers.FluxpartResult`

    Other Parameters
    ----------------
    hfd_options['cols'] : 7*(int,)
        7-tuple of integers indicating the column numbers of `fname`
        that contain series data for (u, v, w, c, q, T, P), in that
        order. Uses 0-based indexing. Default is (2, 3, 4, 5, 6, 7, 8)
    hfd_options['unit_convert'] : dict
        Dictionary of multiplication factors required to convert any u,
        v, w, c, q, or P data not in SI units to SI units (m/s, kg/m^3,
        Pa). (Note T is not in that list). The dictionary keys are the
        variable names. For example, if all data in `fname`
        are in SI units except P and c, which are in units of kPa and
        mg/m^3, respectively, then set:
        ``hfd_options['unit_convert'] = {'P': 1e3, 'c': 1e-6}``,
        since it is necessary to multiply the kPa pressure data by 1e3
        to obtain the SI pressure unit (Pa), and the mg/m^3 CO2 data by
        1e-6 to obtain the SI concentration unit (kg/m^3).
    hfd_options['temper_unit'] : {'K', 'C'}
        The units of the temperature data T in `fname`. Default is the
        SI unit, 'K'.
    hfd_options['bounds'] : dict
        Dictionary specifying any prescribed lower and upper bounds for
        valid data. Dictionary entries have the form
        ``varname: (float, float)``, where varname is one of 'u', 'v',
        'w', 'q', 'c', 'T', or 'P', and the 2-tuple holds values for the
        lower and upper bounds: ``(lower, upper)``.  Data records are
        rejected if a variable in the record is outside the prescribed
        bounds. Default is
        ``bounds = {'c': (0, np.inf), 'q': (0, np.inf)}`` such that data
        records are rejected if c or q data are not positive values.
    hfd_options['flags'] : dict
        Specifies that one or more columns in `fname` are used to flag
        bad data records. Dict keys are flag names (can be an arbitrary,
        unique identifier); values are 2-tuples of the form
        (col, badval), where col is an int specifying the column number
        containing the flag (0-based indexing), and badval is the value
        of the flag that indicates a bad data record. Default is None.
    hfd_options['rd_tol'] : float
        Relative tolerance for rejecting the datafile. Default is
        'hfd_options['rd_tol']` = 0.4. See
        :class:`~fluxpart.hfdata.HFData`.
    hfd_options['ad_tol'] : int
        Absolute tolerance for rejecting the datafile. Default is
        `hfd_options['ad_tol']` = 1024. See
        :class:`~fluxpart.hfdata.HFData`.
    hfd_options['ustar_tol'] : float
        If the friction velocity (m/s) determined from the high
        frequency data is less than `hfd_options['ustar_tol']`, the
        partitioning analysis is aborted due to insufficient turbulence.
        Defalult is `hfd_options['ustar_tol']` = 0.1 (m/s).
    hfd_options['correcting_external'] : bool, optional
        If True (default), the water vapor and carbon dioxide series
        data are corrected for external fluctuations associated with air
        temperature and vapor density according to [WPL80]_ and [DK07]_.
    hfd_options[ other keys ]
        All other key:value pairs in `hfd_options` are passed as keyword
        arguments to pandas.read_csv_ (where the file is read). These
        keywords are often required to specify the details of the
        formatting of the delimited datafile.  Among the most
        commonly required are: 'delimiter', a str, int, or sequence
        that is used to separate values or define column widths (default
        is that any consecutive whitespace delimits values); and
        'skiprows', an int, list of ints, or callable that specifies
        rows to skip (e.g. header rows at the beginning of the file. See
        pandas.read_csv_ for a full description of available format
        options.
    wue_options['canopy_ht'] : float
        Vegetation canopy height (m).
    wue_options['meas_ht'] : float
        Eddy covariance measurement height (m).
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
        be equal to the air temperature. See:
        :func:`~fluxpart.wue.water_use_efficiency`.
    wue_options['leaf_temper_corr'] : float
        Offset adjustment applied to canopy temperature (K).
        See: :func:`~fluxpart.wue.water_use_efficiency`.
    wue_options['diff_ratio']: float, optional
        Ratio of molecular diffusivities for water vapor and CO2.
        Default is `diff_ratio` = 1.6.
        See: :func:`~fluxpart.wue.water_use_efficiency`.
    part_options['adjust_fluxes'] : bool
        If True (default), the final partitioned fluxes are adjusted
        proportionally such that sum of the partitioned fluxes match
        exactly the total fluxes indicated in the original data.


    .. _pandas.read_csv:
        https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

    """
    def abort(dataread, mssg):
        hfs = hfsum if dataread else HFSummary()
        return FluxpartResult(
            label=label,
            version=VERSION,
            dataread=dataread,
            attempt_partition=False,
            valid_partition=False,
            mssg=mssg,
            hfsummary=hfs,
        )

    hfd_options = {**DEFAULT_HFD_OPTIONS, **(hfd_options or {})}
    wue_options = {**DEFAULT_WUE_OPTIONS, **(wue_options or {})}
    part_options = {**DEFAULT_PART_OPTIONS, **(part_options or {})}

    usecols = np.array(hfd_options.pop('cols'), dtype=int).reshape(7,)

    converters = None
    unit_convert = hfd_options.pop('unit_convert', {})
    converters = {
        k: _str_converter_func(float(v), 0.) for k, v in unit_convert.items()}

    temper_unit = hfd_options.pop('temper_unit')
    if temper_unit.upper() == 'C' or temper_unit.upper() == 'CELSIUS':
        converters = converters or {}
        converters['T'] = _str_converter_func(1., 273.15)

    correcting_external = hfd_options.pop('correcting_external')
    ustar_tol = hfd_options.pop('ustar_tol')
    bounds = hfd_options.pop('bounds')
    rd_tol = hfd_options.pop('rd_tol')
    ad_tol = hfd_options.pop('ad_tol')

    # read high frequency data
    try:
        hfdat = HFData(cols=usecols, converters=converters, **hfd_options)
        hfdat.read(fname)
        hfdat.quality_check(bounds, rd_tol, ad_tol)
    except HFDataReadError as err:
        abort(dataread=False, mssg=err.args[0])

    # preliminary data processing and analysis
    hfdat.truncate()
    if correcting_external:
        hfdat.correct_external()
    hfsum = hfdat.summarize()

    # verify data are compatible with fvs partitioning

    if hfsum.ustar < ustar_tol:
        mssg = ('ustar = {:.4} is less than ustar_tol = {:.4}'.
                format(hfsum.ustar, ustar_tol))
        abort(dataread=True, mssg=mssg)

    vpd = vapor_press_deficit(hfsum.rho_vapor, hfsum.T)
    if vpd <= 0:
        abort(dataread=True, mssg=f'Vapor pressure deficit {vpd}')

    if hfsum.cov_w_q <= 0:
        mssg = ('cov(w,q) = {:.4} <= 0 is incompatible with partitioning '
                'algorithm'.format(hfsum.cov_w_q))
        abort(dataread=True, mssg=mssg)

    if meas_wue:
        leaf_wue = WUE(wue=float(meas_wue))
    else:
        try:
            leaf_wue = water_use_efficiency(hfsum, **wue_options)
        except WUEError as err:
            abort(dataread=True, mssg=err.args[0])

    # compute partitioned fluxes

    adjust_fluxes = part_options['adjust_fluxes']
    fvsp = (
        fvspart_progressive(
            hfdat['w'].values,
            hfdat['q'].values,
            hfdat['c'].values,
            leaf_wue.wue,
            adjust_fluxes,
        )
    )

    if fvsp.valid_partition:
        fvsp.fluxes = AllFluxes(
                **attr.asdict(fvsp.fluxes), temper_kelvin=hfsum.T)
    else:
        fvsp.fluxes = AllFluxes()

    return FluxpartResult(
        label=label,
        version=VERSION,
        dataread=True,
        attempt_partition=True,
        valid_partition=fvsp.valid_partition,
        mssg=fvsp.mssg,
        fvsp_result=fvsp,
        hfsummary=hfsum,
        wue=leaf_wue,
        )


@attr.s
class FluxpartResult(object):
    """Fluxpart result.

    Parameters
    ----------
    version : str
        Fluxpart version
    dataread, attempt_partition, valid_partition : bool
        Indicates success or failure in reading high frequency data,
        attempting and obtaining a valid partioning solution.
    mssg : str
        Possibly informative message if `dataread` or `valid_partition`
        are False

    TODO

    """
    version = attr.ib()
    dataread = attr.ib()
    attempt_partition = attr.ib()
    valid_partition = attr.ib()
    mssg = attr.ib()
    fvsp_result = attr.ib(
            default=FVSPResult(WQCData(), RootSoln(), AllFluxes()))
    hfsummary = attr.ib(default=HFSummary())
    wue = attr.ib(default=WUE())
    label = attr.ib(default=None)

    # TODO
    # def __str__(self):
    #     return (
    #         'Outcome(\n'
    #         + f'    version = {self.version},\n'
    #         + f'    dataread = {self.dataread},\n'
    #         + f'    attempt_partition = {self.attempt_partition},\n'
    #         + f'    valid_partition = {self.valid_partition},\n'
    #         + f'    mssg = {self.mssg})'
    #         )


def _str_converter_func(slope, intercept):
    """Return func for linear transform of data when reading file."""
    def func(stringval):
        try:
            return slope * float(stringval) + intercept
        except ValueError:
            return np.nan
    return func


def _converter_func(slope, intercept):
    """Return a function for linear transform of data."""
    def func(val):
        return slope * val + intercept
    return func


if __name__ == "__main__":
    pass
