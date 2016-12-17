import numpy as np

import fluxpart.partition as fp
import fluxpart.wue as wue
import fluxpart.util as util
from fluxpart.hfdata import HFData
from fluxpart.containers import Fluxes, WUE, Result, NumerSoln

DEFAULT_WUE_OPTIONS = {
    'ci_mod': 'const_ratio',
    'ci_mod_param': None,
    'leaf_temper': None,
    'leaf_temper_corr': 0}

DEFAULT_HFD_OPTIONS = {
    'cols': (1, 2, 3, 4, 5, 6, 7),
    'unit_convert': None,
    'temper_unit': 'K',
    'bounds': None,
    'flags': None,
    'rd_tol': 0.4,
    'ad_tol': 1024,
    'ustar_tol': 0.1,
    'correcting_external': True}

DEFAULT_PART_OPTIONS = {
    'adjusting_fluxes': True}


def flux_partition(fname, meas_wue=None, hfd_options=None, wue_options=None,
                   part_options=None, label=None):
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
    dict
        {'result': :class:`~fluxpart.containers.Result`,
        'fluxes': :class:`~fluxpart.containers.Fluxes`,
        'datsumm': :class:`~fluxpart.containers.HFSummary`,
        'wue': :class:`~fluxpart.containers.WUE`,
        'numsoln': :class:`~fluxpart.containers.NumerSoln`,
        'label': `label`}

    Other Parameters
    ----------------
    hfd_options['cols'] : 7*(int,)
        7-tuple of integers indicating the column numbers of `fname`
        that contain series data for (u, v, w, q, c, T, P), in that
        order. Uses 0-based indexing. Default is (1, 2, 3, 4, 5, 6, 7)
        (thus the first column (=0) in the file is not read)
    hfd_options['unit_convert'] : dict
        Dictionary of multiplication factors required to convert any u,
        v, w, q, c, or P data not in SI units to SI units (m/s, kg/m^3,
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
    hfd_options['flags'] : 2-tuple or list of 2-tuples
        Specifies that one or more columns in `fname` are used to flag
        bad data records. Each tuple is of the form (col, badval),
        where col is an int specifying the column number containing the
        flag (0-based indexing), and badval is the value of the flag
        that indicates a bad data record. Default is None.
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
        arguments to numpy.genfromtxt_ (where the file is read). These
        keywords are often required to specify the details of the
        formatting of the delimited datafile.  Among the most
        commonly required are: 'delimiter', a str, int, or sequence
        that is used to separate values or define column widths (default
        is that any consecutive whitespace delimits values); and
        'skip_header', an int that specifies the number of lines to skip
        at the beginning of the file. See numpy.genfromtxt_ for a full
        description of available format options.
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
    part_options['adjusting_fluxes'] : bool
        If True (default), the final partitioned fluxes are adjusted
        proportionally such that sum of the partitioned fluxes match
        exactly the total fluxes indicated in the original data.


    .. _numpy.genfromtxt:
        http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html

    """

    # py3.5 dictionary merge/update
    hfd_options = {**DEFAULT_HFD_OPTIONS, **(hfd_options or {})}
    wue_options = {**DEFAULT_WUE_OPTIONS, **(wue_options or {})}
    part_options = {**DEFAULT_PART_OPTIONS, **(part_options or {})}

    usecols = np.array(hfd_options.pop('cols'), dtype=int).reshape(7,)

    converters = None
    unit_convert = hfd_options.pop('unit_convert')
    if unit_convert:
        converters = {
            k: _converter_func(float(v), 0.) for k, v in unit_convert.items()}

    temper_unit = hfd_options.pop('temper_unit')
    if temper_unit.upper() == 'C' or temper_unit.upper() == 'CELSIUS':
        converters = converters or {}
        converters['T'] = _converter_func(1., 273.15)

    correcting_external = hfd_options.pop('correcting_external')
    ustar_tol = hfd_options.pop('ustar_tol')

    # read high frequency data
    try:
        hfdat = HFData(fname, cols=usecols, converters=converters,
                       **hfd_options)
    except (TypeError, ValueError) as err:
        mssg = 'HFData read fail: ' + err.args[0]
        result = Result(dataread=False, attempt_partition=False,
                        valid_partition=False, mssg=mssg)
        return {'label': label,
                'result': result,
                'fluxes': Fluxes(*np.full(15, np.nan)),
                'datsumm': None,
                'wue': WUE(*np.full(11, np.nan)),
                'numsoln': NumerSoln(*np.full(10, np.nan))}

    # preliminary data processing and analysis
    hfdat.truncate()
    if correcting_external:
        hfdat.qc_correct()
    hfsum = hfdat.summarize()

    # exit if friction velocity is too low (lack of turbulence)
    if hfsum.ustar < ustar_tol:
        mssg = ('ustar = {:.4} is less than ustar_tol = {:.4}'.
                format(hfsum.ustar, ustar_tol))
        result = Result(dataread=True, attempt_partition=False,
                        valid_partition=False, mssg=mssg)
        return {'label': label,
                'result': result,
                'fluxes': Fluxes(*np.full(15, np.nan)),
                'datsumm': hfsum,
                'wue': WUE(*np.full(11, np.nan)),
                'numsoln': NumerSoln(*np.full(10, np.nan))}

    # exit if water vapor flux is downward (negative)
    if hfsum.cov_w_q <= 0:
        mssg = ('cov(w,q) = {:.4} <= 0 is incompatible with partitioning '
                'algorithm'.format(hfsum.cov_w_q))
        result = Result(dataread=True, attempt_partition=False,
                        valid_partition=False, mssg=mssg)
        return {'label': label,
                'result': result,
                'fluxes': _no_partitioned_fluxes(hfsum),
                'datsumm': hfsum,
                'wue': WUE(*np.full(11, np.nan)),
                'numsoln': NumerSoln(*np.full(10, np.nan))}

    # get or calculate water use efficiency
    if meas_wue:
        leaf_wue = WUE(float(meas_wue), *10 * (np.nan,))
    else:
        leaf_wue = wue.water_use_efficiency(hfsum, **wue_options)

    # exit if wue value is bad
    if not leaf_wue.wue < 0:
        mssg = ('wue={} must be less than zero'.format(leaf_wue.wue))
        result = Result(dataread=True, attempt_partition=False,
                        valid_partition=False, mssg=mssg)
        return {'label': label,
                'result': result,
                'fluxes': _no_partitioned_fluxes(hfsum),
                'datsumm': hfsum,
                'wue': leaf_wue,
                'numsoln': None}

    # compute partitioned fluxes
    adjusting_fluxes = part_options['adjusting_fluxes']
    pout = fp.partition_from_wqc_series(hfdat['w'], hfdat['q'], hfdat['c'],
                                        leaf_wue.wue, adjusting_fluxes)

    # collect results and return
    result = Result(dataread=True,
                    attempt_partition=True,
                    valid_partition=pout['valid_partition'],
                    mssg=pout['partmssg'])

    if pout['valid_partition']:
        fluxes = Fluxes(
            *pout['fluxcomps'],
            LE=util.qflux_mass_to_heat(pout['fluxcomps'].wq, hfsum.T),
            LEt=util.qflux_mass_to_heat(pout['fluxcomps'].wqt, hfsum.T),
            LEe=util.qflux_mass_to_heat(pout['fluxcomps'].wqe, hfsum.T),
            Fq_mol=util.qflux_mass_to_mol(pout['fluxcomps'].wq),
            Fqt_mol=util.qflux_mass_to_mol(pout['fluxcomps'].wqt),
            Fqe_mol=util.qflux_mass_to_mol(pout['fluxcomps'].wqe),
            Fc_mol=util.cflux_mass_to_mol(pout['fluxcomps'].wc),
            Fcp_mol=util.cflux_mass_to_mol(pout['fluxcomps'].wcp),
            Fcr_mol=util.cflux_mass_to_mol(pout['fluxcomps'].wcr))
    else:
        fluxes = _no_partitioned_fluxes(hfsum)

    return {'label': label,
            'result': result,
            'fluxes': fluxes,
            'datsumm': hfsum,
            'wue': leaf_wue,
            'numsoln': pout['numsoln'],
            'qcdat': pout['qcdat']}


def _converter_func(slope, intercept):
    """Return a function for linear transform of data when reading file.
    """
    def func(stringval):
        try:
            return slope * float(stringval.strip()) + intercept
        except ValueError:
            return np.nan
    return func


def _no_partitioned_fluxes(hfsum):
    return Fluxes(
        Fq=hfsum.cov_w_q,
        Fqe=np.nan,
        Fqt=np.nan,
        Fc=hfsum.cov_w_c,
        Fcr=np.nan,
        Fcp=np.nan,
        LE=util.qflux_mass_to_heat(hfsum.cov_w_q, hfsum.T),
        LEt=np.nan,
        LEe=np.nan,
        Fq_mol=util.qflux_mass_to_mol(hfsum.cov_w_q),
        Fqt_mol=np.nan,
        Fqe_mol=np.nan,
        Fc_mol=util.cflux_mass_to_mol(hfsum.cov_w_c),
        Fcp_mol=np.nan,
        Fcr_mol=np.nan)


if __name__ == "__main__":
    pass
