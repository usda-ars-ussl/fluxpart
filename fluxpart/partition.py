import math
import numpy as np
from scipy import optimize
import fluxpart.util as util
from fluxpart.containers import NumerSoln, FluxComponents, QCData


def partition_from_wqc_series(w, q, c, wue, adjust_fluxes=True):
    """Partition H2O and CO2 fluxes using series data for w, q, and c.

    If a valid partitioning solution is not found for the passed series
    data, low-frequency (large-scale) components are progressively
    removed from the data until either a valid solution is found or
    the series decomposition is exhausted.

    Parameters
    ----------
    w,q,c : array_like
        1D time series for vertical wind speed `w` (m/s), water vapor
        concentration `q` (kg/m^3), and CO2 concentration `c` (kg/m^3).
    wue : float
        leaf-level water use efficiency, `wue` < 0, kg CO2 / kg H2O.
    adjust_fluxes : bool, optional
        Indicates whether the obtained partitioned fluxes should be
        adjusted so that the totals match the original data. Default is
        `adjust_fluxes` = True.

    Returns
    -------
    dict
        {'valid_partition': bool, 'partmssg': str,
        'fluxcomps': :class:`~fluxpart.containers.FluxComponents`,
        'numsoln': :class:`~fluxpart.containers.NumerSoln`,
        'qcdat': :class:`~fluxpart.containers.QCData`}

    Notes
    -----
    If a valid partitioning is not found, the returned `numsoln` and `qcdat`
    correspond to the final iteration attempted.

    """

    max_decomp_lvl = int(np.log2(w.size))
    wq_tot = np.cov((w, q))[0, 1]
    wc_tot = np.cov((w, c))[0, 1]

    # The loop progressively filters the data until a physically valid
    # partitioning is found or the loop/filter is exhausted. The first
    # iteration of progressive_lowcut removes only the mean value
    # (q'=q-<q>, etc.), so the first iteration uses the unfiltered
    # deviations.

    for cnt, lowcut_wqc in enumerate(progressive_lowcut(w, q, c)):
        cov = np.cov(lowcut_wqc)
        qcdat = QCData(
            wq=cov[0, 1],
            wc=cov[0, 2],
            var_q=cov[1, 1],
            var_c=cov[2, 2],
            corr_qc=cov[1, 2] / math.sqrt(cov[1, 1] * cov[2, 2]),
            wave_lvl=(max_decomp_lvl - cnt, max_decomp_lvl))

        nsoln, fcomp = partition_from_qc_averages(qcdat, wue)

        if nsoln.success and nsoln.validroot:

            if adjust_fluxes:
                fcomp = adjust_partitioned_fluxes(fcomp, wue, wq_tot, wc_tot)

            valid_partition, partition_mssg = isvalid_partition(fcomp)

            if valid_partition:
                break

    if not nsoln.success:
        valid_partition = False
        partition_mssg = 'Root finding failed: ' + nsoln.mssg

    if nsoln.success and not nsoln.validroot:
        valid_partition = False
        partition_mssg = ('Root located but not physically valid:'
                          + nsoln.validmssg)
    if not valid_partition:
        fcomp = FluxComponents(*np.full(6, np.nan))

    return {'valid_partition': valid_partition,
            'partmssg': partition_mssg,
            'fluxcomps': fcomp,
            'numsoln': nsoln,
            'qcdat': qcdat}


def partition_from_qc_averages(qcdat, wue, init=None):
    """Partition H2O and CO2 fluxes using interval average q and c data.

    All arguments are passed directly to the findroot function

    Parameters
    ----------
    qcdat : QCData namedtuple or equivalent namespace
    wue : float, kg CO2 / kg H2O
        Leaf-level water use efficiency, `wue` < 0
    init : (float, float), optional
        2-Tuple is initial value for (corr_cp_cr, var_cp).  If `init`
        = None (default), an initial estimate is calculated internally.
        Note when specifying initial values: -1 < `corr_cp_cr` < 0, and
        `var_cp` has units of (kg/m^3)^2.

    Returns
    -------
    namedtuples
        :class:`~fluxpart.containers.NumerSoln`,
        :class:`~fluxpart.containers.FluxComponents`

    """

    nsoln = findroot(qcdat, wue, init)
    if nsoln.success and nsoln.validroot:
        fluxes = flux_components(nsoln.var_cp, nsoln.corr_cp_cr, qcdat, wue,
                                 nsoln.co2soln_id)
    else:
        fluxes = FluxComponents(*np.full(6, np.nan))
    return nsoln, fluxes


def findroot(qcdat, wue, init=None):
    """Solve numerically for (corr_cp_cr, var_cp).

    Parameters
    ----------
    qcdat : namedtuple or equivalent namespace
        :class:`~fluxpart.containers.QCData`
    wue : float
        Leaf-level water use efficiency, `wue` < 0, kg CO2 / kg H2O.
    init : (float, float), optional
        2-Tuple is initial value for (corr_cp_cr, var_cp).  If `init`
        = None (default), an initial estimate is calculated internally.
        Note when specifying initial values: -1 < corr_cp_cr < 0, and
        `var_cp` has units of (kg/m^3)^2.

    Returns
    -------
    namedtuple
        :class:`~fluxpart.containers.NumerSoln`

    """

    if init is None:
        init_corr_cp_cr = -0.8
        varcp_ubound0 = qcdat.var_c / (1 - init_corr_cp_cr**2)
        varcp_ubound1 = wue**2 * qcdat.var_q / (1 - init_corr_cp_cr**2)
        init = (init_corr_cp_cr, 0.5 * min(varcp_ubound0, varcp_ubound1))

    co2_ids = (1, 0) if qcdat.wc < 0 else (0, )
    for co2_id in co2_ids:
        co2soln_id = co2_id
        try:
            soln = optimize.root(residual_func, init, method='hybr',
                                 args=(qcdat, wue, co2soln_id))
        except ValueError as err:
            soln = {'success': False,
                    'message': 'optimize.root fail: ' + err.args[0],
                    'nfev': np.nan}
        if soln['success']:
            break

    if soln['success']:
        corr_cp_cr, var_cp = soln['x']
        valid_root, valid_root_mssg = isvalid_root(corr_cp_cr, var_cp)
        if var_cp > 0:
            wcr_ov_wcp = flux_ratio(var_cp, corr_cp_cr, qcdat, 'co2', co2soln_id)
            sig_cr = wcr_ov_wcp * math.sqrt(var_cp) / corr_cp_cr
        else:
            sig_cr = np.nan
    else:
        var_cp, corr_cp_cr = np.nan, np.nan
        sig_cr = np.nan
        valid_root, valid_root_mssg = False, ''

    return NumerSoln(
        corr_cp_cr=corr_cp_cr,
        var_cp=var_cp,
        sig_cr=sig_cr,
        co2soln_id=co2soln_id,
        validroot=valid_root,
        validmssg=valid_root_mssg,
        success=soln['success'],
        mssg=soln['message'],
        init=init,
        nfev=soln['nfev'])


def flux_ratio(var_cp, corr_cp_cr, qcdat, ftype, farg):
    """Compute the nonstomatal:stomatal ratio of the H2O or CO2 flux.

    The ratio (either wqe/wqt or wcr/wcp) is found by solving Eq. 13
    of [SS08]_.

    Parameters
    ---------
    qcdat : namedtuple or equivalent namespace
        :class:`~fluxpart.containers.QCData`
    ftype : {'co2', 'h2o'}
        Specifies whether the flux is CO2 or H2O
    farg : number
        If `ftype` = 'co2', then `farg` = {0 or 1} specifies the root
        of Eq. 13b to be used to calculate the CO2 flux ratio wcr/wcp:
        `farg`=1 uses the '+' solution, `farg`=0 uses the '-' solution.
        If `ftype` = 'h2o', then `farg` is a float equal to the water
        use efficiency (wue < 0), kg/kg.

    Returns
    -------
    fratio : float or np.nan
        The requested flux ratio, wqe/wqt or wcr/wcp. Set to np.nan
        if solution is not real.

    Notes
    -----
    When solving for wqe/wqt, the '-' solution of the quadratic Eq. 13a
    is not relevant because only the '+' solution yields wqe/wqt > 0,
    which is required/assumed by the physical model in [SS08]_.

    """

    if ftype == 'co2' or ftype == 'CO2':
        sign = 1 if farg == 1 else -1
        num = qcdat.var_c
    elif ftype == 'h2o' or ftype == 'H2O':
        sign = 1
        num = farg**2 * qcdat.var_q
    else:
        raise ValueError("ftype must be 'co2' or 'h2o'")
    disc = 1 - 1 / corr_cp_cr**2 + num / var_cp / corr_cp_cr**2
    if disc < 0:
        fratio = np.nan
    else:
        fratio = corr_cp_cr**2 * (sign * math.sqrt(disc) - 1)
    return fratio


def flux_components(var_cp, corr_cp_cr, qcdat, wue, co2soln_id):
    """Calculate flux components for given (var_cp, corr_cp_cr) pair."""
    wcr_ov_wcp = flux_ratio(var_cp, corr_cp_cr, qcdat, 'co2', co2soln_id)
    # TODO: handle wcr_ov_wcp ~ -1
    wcp = qcdat.wc / (wcr_ov_wcp + 1)
    wcr = qcdat.wc - wcp
    wqt = wcp / wue
    wqe = qcdat.wq - wqt
    return FluxComponents(wq=qcdat.wq, wqt=wqt, wqe=wqe,
                          wc=qcdat.wc, wcp=wcp, wcr=wcr)


def residual_func(x, qcdat, wue, co2soln_id):
    """Residual function used with root finding routine.

    The two components of the residual are for Eqs. 15 and 18 of [SS08]_.

    """

    corr_cp_cr, var_cp = x
    wcr_ov_wcp = flux_ratio(var_cp, corr_cp_cr, qcdat, 'co2', co2soln_id)
    wqe_ov_wqt = flux_ratio(var_cp, corr_cp_cr, qcdat, 'h2o', wue)

    # Eq. 3, "f1"
    lhs = wue * qcdat.wq / qcdat.wc * (wcr_ov_wcp + 1)
    rhs = wqe_ov_wqt + 1
    resid1 = lhs - rhs

    # Eq. 4, "f2"
    lhs = wue * qcdat.corr_qc * math.sqrt(qcdat.var_c * qcdat.var_q)
    rhs = var_cp * (
        1 + wqe_ov_wqt + wcr_ov_wcp + wqe_ov_wqt * wcr_ov_wcp / corr_cp_cr**2)
    resid2 = lhs - rhs
    return [resid1, resid2]


def isvalid_root(corr_cp_cr, var_cp):
    isvalid = True
    mssg = ''
    if var_cp <= 0:
        isvalid = False
        mssg += 'var_cp <= 0; '
    if not -1 < corr_cp_cr < 0:
        isvalid = False
        mssg += 'corr_cp_cr <-1 OR >0; '
    # TODO: could add other bound checks?
    return isvalid, mssg


def isvalid_partition(pf):
    """Test if partitioned flux components are valid.

    Parameters
    ----------
    pf : namedtuple
        :class:`~fluxpart.containers.FluxComponents`

    """
    isvalid = True
    mssg = ''
    if pf.wqt <= 0:
        isvalid = False
        mssg += 'Fqt <= 0; '
    if pf.wqe <= 0:
        isvalid = False
        mssg += 'Fqe <= 0; '
    if pf.wcp >= 0:
        isvalid = False
        mssg += 'Fcp >= 0; '
    if pf.wcr <= 0:
        isvalid = False
        mssg += 'Fcr <= 0; '
    return isvalid, mssg


def adjust_partitioned_fluxes(fc, wue, wq_tot, wc_tot):
    """Adjust partitioned fluxes so they match measured totals.

    If filtering has been applied to the series data, covariances in
    the filtered data may differ from those in the original data.
    Consequently, partitioned flux totals may not match exactly the
    total fluxes indicated by the original data. Here, partitioned fluxes
    are adjusted proportionally so that they match the totals in the
    original data.

    Parameters
    ----------
    fc : :class:`~fluxpart.containers.FluxComponents` or equivalent
        Container holding partitioned flux components, kg/m^2/s.
    wue : float
        Leaf-level water use efficiency (`wue` < 0), kg CO2 / kg H2O
    wq_tot, wc_tot : float
        Desired net total H2O (`wq_tot`) and CO2 (`wc_tot`) fluxes,
        kg/m^2/s.

    Returns
    -------
    namedtuple
        :class:`~fluxpart.containers.FluxComponents`

    """

    wq_diff = wq_tot - (fc.wqe + fc.wqt)
    wqe = fc.wqe + wq_diff * (fc.wqe / (fc.wqt + fc.wqe))
    wqt = wq_tot - wqe
    wcp = wue * (fc.wqt + wq_diff * (fc.wqt / (fc.wqt + fc.wqe)))
    wcr = wc_tot - wcp
    return FluxComponents(wq=wq_tot, wqt=wqt, wqe=wqe,
                          wc=wc_tot, wcp=wcp, wcr=wcr)


def progressive_lowcut(wind, vapor, co2):
    """Apply progressive lowcut filter to wind, vapor, and CO2 series.

    Uses wavelet decompostion to yield a sequence of (w, q, c) series
    in which low frequency (large scale) components are progressively
    removed from w, q, c.

    Parameters
    ----------
    wind,vapor,co2 : array-like
        1D time series of vertical wind velocity (w), water vapor
        concentration (q), and carbon dioxide concentration (c).

    Yields
    ------
    (lc_w, lc_q, lc_c) : tuple of arrays
        Arrays lc_w, lc_q, and lc_c are low cut (high pass) filtered
        versions of the passed w,q,c series.

    Notes
    -----
    The time series data are truncated at the maxium dyadic length
    (power of 2) before the filter is applied.

    """

    max_dyadic_len = 2**int(np.log2(np.asarray(co2).shape[0]))
    trunc_w = np.asarray(wind)[:max_dyadic_len]
    trunc_q = np.asarray(vapor)[:max_dyadic_len]
    trunc_c = np.asarray(co2)[:max_dyadic_len]
    lowcut_w = util.progressive_lowcut_series(trunc_w)
    lowcut_q = util.progressive_lowcut_series(trunc_q)
    lowcut_c = util.progressive_lowcut_series(trunc_c)
    for lowcut_series in zip(lowcut_w, lowcut_q, lowcut_c):
        yield lowcut_series


if __name__ == '__main__':
    pass
