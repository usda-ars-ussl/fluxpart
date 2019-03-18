"""TODO:"""
import math
from types import SimpleNamespace

import numpy as np

import fluxpart.util as util
from .containers import FVSPSolution, MassFluxes, RootSoln, WQCData


class Error(Exception):
    pass


class FVSError(Error):
    def __init__(self, message):
        self.message = message


def fvspart_progressive(w, q, c, wue, adjust_fluxes=True):
    """FVS flux partitioning using high frequency eddy covariance data.

    If a valid partitioning solution is not found for the passed series
    data, low-frequency (large-scale) components are progressively
    removed from the data until either a valid solution is found or
    the series decomposition is exhausted.

    Parameters
    ----------
    w,q,c : array_like
        1D high frequency time series for vertical wind speed `w` (m/s),
        water vapor `q` (kg/m^3), and CO2 `c` (kg/m^3).
    wue : float, `wue` < 0
        leaf-level water use efficiency, kg CO2 / kg H2O.
    adjust_fluxes : bool, optional
        Indicates whether the obtained partitioned fluxes should be
        adjusted so that the totals match the original data. Default is
        True.

    Returns
    -------
    :class:`~fluxpart.containers.FVSPSolution`

    Notes
    -----
    If a valid partitioning is not found, the returned `numersoln` and
    `wqc_data` correspond to the final iteration attempted.

    """
    max_decomp_lvl = int(np.log2(w.size))
    wq_tot = np.cov((w, q))[0, 1]
    wc_tot = np.cov((w, c))[0, 1]

    # The loop progressively filters the data until a physically valid
    # partitioning is found or the loop/filter is exhausted. The first
    # iteration of progressive_lowcut removes only the mean value
    # (q'=q-<q>, etc.), so the first iteration uses the "unfiltered"
    # deviations.

    for cnt, lowcut_wqc in enumerate(_progressive_lowcut(w, q, c)):
        wave_lvl = (max_decomp_lvl - cnt, max_decomp_lvl)

        fluxes, fvsp = fvspart_series(*lowcut_wqc, wue)
        if cnt == 0:
            if fvsp.fvsp_mssg:
                mssg_for_unfiltered_data = fvsp.fvsp_mssg
            else:
                mssg_for_unfiltered_data = fvsp.rootsoln.root_mssg

        if fvsp.rootsoln.valid_root:
            if adjust_fluxes:
                fluxes = _adjust_fluxes(fluxes, wue, wq_tot, wc_tot)
                fvsp.valid_partition, fvsp.mssg = _isvalid_partition(fluxes)
            if fvsp.valid_partition:
                break

    fvsp.wave_lvl = wave_lvl
    if not fvsp.rootsoln.valid_root:
        fvsp.valid_partition = False
        fvsp.mssg = mssg_for_unfiltered_data
    if not fvsp.valid_partition:
        fluxes = MassFluxes()
        fvsp.mssg = mssg_for_unfiltered_data
    return fluxes, fvsp


def fvspart_series(w, q, c, wue, wipe_if_invalid=False):
    """FVS partition q and c fluxes using high frequency eddy cov data.

    Parameters
    ----------
    w,q,c : array_like
        1D high frequency time series for vertical wind speed `w` (m/s),
        water vapor `q` (kg/m^3), and CO2 `c` (kg/m^3).
    wue : float, `wue` < 0
        leaf-level water use efficiency, kg CO2 / kg H2O.
    wipe_if_invalid : boolean
        If True, return default (nan) values for all mass fluxes if any
        calculated fluxes violate directional requirements. Default is
        False.

    Returns
    -------
    :class:`~fluxpart.containers.FVSPSolution`,

    """
    cov = np.cov([w, q, c])
    wqc_data = WQCData(
        wq=cov[0, 1],
        wc=cov[0, 2],
        var_q=cov[1, 1],
        var_c=cov[2, 2],
        corr_qc=cov[1, 2] / math.sqrt(cov[1, 1] * cov[2, 2]),
    )
    return fvspart_interval(wqc_data, wue)


def fvspart_interval(wqc_data, wue, wipe_if_invalid=False):
    """Partition H2O and CO2 fluxes using interval averaged data values.

    Parameters
    ----------
    wqc_data : :class:`~fluxpart.containers.WQCData`
    wue : float, kg CO2 / kg H2O
        Leaf-level water use efficiency, `wue` < 0
    wipe_if_invalid : boolean
        If True, return default (nan) values for all mass fluxes if any
        calculated fluxes violate directional requirements. Default is
        False.

    Returns
    -------
    :class:`~fluxpart.containers.FVSPSolution`

    """
    try:
        _check_fvs_assumptions(wqc_data)
    except FVSError as e:
        mass_fluxes = MassFluxes()
        fvsps = FVSPSolution(
            wqc_data=wqc_data,
            valid_partition=False,
            fvsp_mssg=e.args[0],
            rootsoln=RootSoln(valid_root=False),
        )
        return mass_fluxes, fvsps

    rootsoln = findroot(wqc_data, wue)
    if not rootsoln.valid_root:
        fvsp = FVSPSolution(
            wqc_data=wqc_data,
            rootsoln=rootsoln,
            valid_partition=False,
            fvsp_mssg=rootsoln.root_mssg,
        )
        return MassFluxes, fvsp
    mass_fluxes = _mass_fluxes(
        var_cp=rootsoln.var_cp,
        corr_cp_cr=rootsoln.corr_cp_cr,
        wqc_data=wqc_data,
        wue=wue,
        co2soln_id=rootsoln.co2soln_id,
    )
    isvalid, mssg = _isvalid_partition(mass_fluxes)
    if not isvalid and wipe_if_invalid:
        mass_fluxes = MassFluxes()
    fvsps = FVSPSolution(
        wqc_data=wqc_data,
        rootsoln=rootsoln,
        valid_partition=isvalid,
        fvsp_mssg=mssg,
    )
    return mass_fluxes, fvsps


def findroot(wqc_data, wue):
    """Calculate (corr_cp_cr, var_cp).

    Parameters
    ----------
    wqc_data : namedtuple or equivalent namespace
        :class:`~fluxpart.containers.WQCData`
    wue : float
        Leaf-level water use efficiency, `wue` < 0, kg CO2 / kg H2O.

    Returns
    -------
    namedtuple
        :class:`~fluxpart.containers.RootSoln`

    """
    co2soln_id = None
    sig_cr, var_cp, corr_cp_cr = np.nan, np.nan, np.nan
    valid_root = False
    valid_mssg = "Series data and WUE value incompatible with FVS"

    var_q, var_c = wqc_data.var_q, wqc_data.var_c
    sd_q, sd_c = math.sqrt(var_q), math.sqrt(var_c)
    wq, wc = wqc_data.wq, wqc_data.wc
    corr_qc = wqc_data.corr_qc

    wcwq = wc / wq
    scsq = sd_c / sd_q
    if corr_qc > 0 and ((wue < wcwq < 0) or (0 < wcwq < scsq * corr_qc)):
        co2soln_id = 0 # minus root
    if corr_qc < 0:
        if wcwq * corr_qc < scsq < wcwq / corr_qc:
            co2soln_id = 1 if scsq <= corr_qc * wue else 0

    if co2soln_id in (0, 1):
        valid_root = True
        valid_mssg = ""

        numer = -2 * corr_qc * sd_c * sd_q * wq * wc
        numer += var_c * wq ** 2 + var_q * wc ** 2
        numer *= -(corr_qc ** 2 - 1) * var_c * var_q * wue ** 2
        denom = -corr_qc * sd_c * sd_q * (wc + wq * wue)
        denom += var_c * wq + var_q * wc * wue
        denom = denom ** 2

        var_cp = numer / denom

        numer = -(corr_qc ** 2 - 1) * var_c * var_q * (wc - wq * wue) ** 2
        denom = -2 * corr_qc * sd_c * sd_q * wc * wq
        denom += var_c * wq ** 2 + var_q * wc ** 2
        denom *= -2 * corr_qc * sd_c * sd_q * wue + var_c + var_q * wue ** 2

        rho_sq = numer / denom
        corr_cp_cr = -math.sqrt(rho_sq)

        wcr_ov_wcp = flux_ratio(
            var_cp, corr_cp_cr, wqc_data, "co2", co2soln_id
        )
        sig_cr = wcr_ov_wcp * math.sqrt(var_cp) / corr_cp_cr

    return RootSoln(
        corr_cp_cr=corr_cp_cr,
        var_cp=var_cp,
        sig_cr=sig_cr,
        co2soln_id=co2soln_id,
        valid_root=valid_root,
        root_mssg=valid_mssg,
    )


def flux_ratio(var_cp, corr_cp_cr, wqc_data, ftype, farg):
    """Compute the nonstomatal:stomatal ratio of the H2O or CO2 flux.

    The ratio (either wqe/wqt or wcr/wcp) is found by solving Eq. 13
    of [SS08]_.

    Parameters
    ---------
    wqc_data : namedtuple or equivalent namespace
        :class:`~fluxpart.containers.WQCData`
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
    if ftype.lower() == "co2":
        sign = 1 if farg == 1 else -1
        num = wqc_data.var_c
    else:  # is "h2o"
        sign = 1
        num = farg ** 2 * wqc_data.var_q
    disc = 1 - 1 / corr_cp_cr ** 2 + num / var_cp / corr_cp_cr ** 2
    if disc < 0:
        fratio = np.nan
    else:
        fratio = corr_cp_cr ** 2 * (sign * math.sqrt(disc) - 1)
    return fratio


def _mass_fluxes(var_cp, corr_cp_cr, wqc_data, wue, co2soln_id):
    """Calculate flux components for given (var_cp, corr_cp_cr) pair."""
    wcr_ov_wcp = flux_ratio(var_cp, corr_cp_cr, wqc_data, "co2", co2soln_id)
    # TODO: handle wcr_ov_wcp ~ -1
    wcp = wqc_data.wc / (wcr_ov_wcp + 1)
    wcr = wqc_data.wc - wcp
    wqt = wcp / wue
    wqe = wqc_data.wq - wqt
    return MassFluxes(
        Fq=wqc_data.wq, Fqt=wqt, Fqe=wqe, Fc=wqc_data.wc, Fcp=wcp, Fcr=wcr
    )


def _check_fvs_assumptions(qcdat):
    pqc = qcdat.corr_qc
    wcwq = qcdat.wc / qcdat.wq
    scsq = np.sqrt(qcdat.var_c / qcdat.var_q)
    if (
        math.isclose(pqc, 0)
        or (pqc < 0 and wcwq <= scsq / pqc)
        or (wcwq >= scsq * pqc)
    ):
        raise FVSError("Series data incompatible with FVS partitioning")


def _isvalid_partition(flux_components):
    """Test if partitioned flux directions (signs) are valid."""
    fc = flux_components
    isvalid = True
    mssg = ""
    if fc.Fqt <= 0:
        isvalid = False
        mssg += "Fqt <= 0; "
    if fc.Fqe <= 0:
        isvalid = False
        mssg += "Fqe <= 0; "
    if fc.Fcp >= 0:
        isvalid = False
        mssg += "Fcp >= 0; "
    if fc.Fcr <= 0:
        isvalid = False
        mssg += "Fcr <= 0; "
    return isvalid, mssg


def _adjust_fluxes(flux_components, wue, Fq_tot, Fc_tot):
    """Adjust partitioned fluxes so they match measured totals.

    If filtering has been applied to the series data, covariances in
    the filtered data may differ from those in the original data.
    Consequently, partitioned flux totals may not match exactly the
    total fluxes indicated by the original data. Here, partitioned
    fluxes are adjusted proportionally so that they match the totals in
    the original data.

    Parameters
    ----------
    flux_components : namedtuple or equivalent
        Attributes (floats) specify the mass flux components (Fq, Fqe,
        Fq, Fc, Fcr, Fcp), kg/m^2/s.
    wue : float
        Leaf-level water use efficiency (`wue` < 0), kg CO2 / kg H2O
    Fq_tot, Fc_tot : float
        Desired net total H2O (`Fq_tot`) and CO2 (`Fc_tot`) fluxes,
        kg/m^2/s.

    Returns
    -------
    namedtuple
        :class:`~fluxpart.containers.MassFluxes`

    """
    fc = flux_components
    Fq_diff = Fq_tot - (fc.Fqe + fc.Fqt)
    Fqe = fc.Fqe + Fq_diff * (fc.Fqe / (fc.Fqt + fc.Fqe))
    Fqt = Fq_tot - Fqe
    Fcp = wue * Fqt
    Fcr = Fc_tot - Fcp
    return MassFluxes(Fq=Fq_tot, Fqt=Fqt, Fqe=Fqe, Fc=Fc_tot, Fcp=Fcp, Fcr=Fcr)


def _progressive_lowcut(wind, vapor, co2):
    """Apply progressive lowcut filter to wind, vapor, and CO2 series.

    Use wavelet decomposition to yield a sequence of (w, q, c) series
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
    Before the filter is applied, the data are truncated so that the
    length is a power of 2.

    """
    max_pow2_len = 2 ** int(np.log2(np.asarray(co2).shape[0]))
    trunc_w = np.asarray(wind)[:max_pow2_len]
    trunc_q = np.asarray(vapor)[:max_pow2_len]
    trunc_c = np.asarray(co2)[:max_pow2_len]
    lowcut_w = util.progressive_lowcut_series(trunc_w)
    lowcut_q = util.progressive_lowcut_series(trunc_q)
    lowcut_c = util.progressive_lowcut_series(trunc_c)
    for lowcut_series in zip(lowcut_w, lowcut_q, lowcut_c):
        yield lowcut_series
