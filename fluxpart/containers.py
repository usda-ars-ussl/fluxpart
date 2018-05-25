"""Classes used to group data and results."""

from collections import namedtuple

import attr
import numpy as np

from .util import qflux_mass_to_heat, qflux_mass_to_mol, cflux_mass_to_mol


@attr.s
class MassFluxes(object):
    """H2O and CO2 mass flux components.

    Attributes
    ----------
    Fq, Fqt, Fqe : float
        Total, transpiration, and evaporation H2O fluxes, kg/m^2/s
    Fc, Fcp, Fcr : float
        Total, photosynthesis, and respiration CO2 fluxes, kg/m^2/s

    """
    Fq = attr.ib(default=np.nan)
    Fqt = attr.ib(default=np.nan)
    Fqe = attr.ib(default=np.nan)
    Fc = attr.ib(default=np.nan)
    Fcp = attr.ib(default=np.nan)
    Fcr = attr.ib(default=np.nan)

    def __str__(self):
        # Print with common mass units instead of SI
        fqs = [1e3 * self.Fq, 1e3 * self.Fqt, 1e3 * self.Fqe]
        fcs = [1e6 * self.Fc, 1e6 * self.Fcp, 1e6 * self.Fcr]
        return (
            'MassFluxes(\n'
            '    Fq = {:.4} g/m^2/s,\n'
            '    Fqt = {:.4} g/m^2/s,\n'
            '    Fqe = {:.4} g/m^2/s,\n'
            '    Fc = {:.4} mg/m^2/s,\n'
            '    Fcp = {:.4} mg/m^2/s,\n'
            '    Fcr = {:.4} mg/m^2/s)'
            ''.format(*(fqs + fcs))
            )


@attr.s
class AllFluxes(object):
    """Water vapor and CO2 fluxes.

    Attributes
    ----------
    Fq, Fqt, Fqe : float
        Total, transpiration, and evaporation H2O fluxes, kg/m^2/s.
    Fc, Fcp, Fcr : float
        Total, photosynthesis, and respiration CO2 fluxes, kg/m^2/s.
    temper_kelvin : float
        Temerature, K

    """
    Fq = attr.ib(default=np.nan)
    Fqt = attr.ib(default=np.nan)
    Fqe = attr.ib(default=np.nan)
    Fc = attr.ib(default=np.nan)
    Fcp = attr.ib(default=np.nan)
    Fcr = attr.ib(default=np.nan)
    temper_kelvin = attr.ib(default=np.nan)
    
    def __attrs_post_init__(self):
        """Water vapor and CO2 fluxes.

        Attributes
        ----------
        LE, LEt, LEe : float
            Same vapor fluxes expressed as latent heat, W/m^2.
        Fq_mol, Fqt_mol, Fqe_mol : float
            Same vapor fluxes expressed as mol/m^2/s.
        Fc_mol, Fcp_mol, Fcr_mol : float
            Same CO2 fluxes expressed as mol/m^2/s.

        """
        self.LE = qflux_mass_to_heat(self.Fq, self.temper_kelvin)
        self.LEt = qflux_mass_to_heat(self.Fqt, self.temper_kelvin)
        self.LEe = qflux_mass_to_heat(self.Fqe, self.temper_kelvin)
        self.Fq_mol = qflux_mass_to_mol(self.Fq)
        self.Fqt_mol = qflux_mass_to_mol(self.Fqt)
        self.Fqe_mol = qflux_mass_to_mol(self.Fqe)
        self.Fc_mol = cflux_mass_to_mol(self.Fc)
        self.Fcp_mol = cflux_mass_to_mol(self.Fcp)
        self.Fcr_mol = cflux_mass_to_mol(self.Fcr)

    def __str__(self):
        # For some fields, print common units instead of SI
        row = "{:35}{}"
        return (
            'Fluxes(\n' +
            row.format(
                '    Fq = {:.4} g/m^2/s,'.format(1e3 * self.Fq),
                'Fc = {:.4} mg/m^2/s,\n'.format(1e6 * self.Fc)) +
            row.format(
                '    Fqt = {:.4} g/m^2/s,'.format(1e3 * self.Fqt),
                'Fcp = {:.4} mg/m^2/s,\n'.format(1e6 * self.Fcp)) +
            row.format(
                '    Fqe = {:.4} g/m^2/s,'.format(1e3 * self.Fqe),
                'Fcr = {:.4} mg/m^2/s,\n\n'.format(1e6 * self.Fcr)) +
            row.format(
                '    Fq_mol = {:.4} mmol/m^2/s,'.format(1e3 * self.Fq_mol),
                'Fc_mol = {:.4} umol/m^2/s,\n'.format(1e6 * self.Fc_mol)) +
            row.format(
                '    Fqt_mol = {:.4} mmol/m^2/s,'.format(1e3 * self.Fqt_mol),
                'Fcp_mol = {:.4} umol/m^2/s,\n'.format(1e6 * self.Fcp_mol)) +
            row.format(
                '    Fqe_mol = {:.4} mmol/m^2/s,'.format(1e3 * self.Fqe_mol),
                'Fcr_mol = {:.4} umol/m^2/s,\n\n'.format(1e6 * self.Fcr_mol)) +
            '    LE  = {:.4} W/m^2,\n'.format(self.LE) +
            '    LEt = {:.4} W/m^2,\n'.format(self.LEt) +
            '    LEe = {:.4} W/m^2)'.format(self.LEe))


@attr.s
class FVSPResult(object):
    """Result of FVS partitioning."""
    wqc_data = attr.ib()
    rootsoln = attr.ib()
    fluxes = attr.ib()
    wave_lvl = attr.ib(default=None)
    valid_partition = attr.ib(default=None)
    mssg = attr.ib(default=None)

    def __str__(self):
        return (
            'FVSPResult:\n' 
            + f'   mssg:: {self.mssg}\n'
            + f'   valid_partition: {self.valid_partition}\n'
            + f'   wave_lvl: {self.wave_lvl}\n'
            + self.wqc_data.__str__() + '\n' 
            + self.fluxes.__str__() + '\n' 
            + self.rootsoln.__str__()
            )


@attr.s
class HFSummary(object):
    """Summary of high frequency eddy covariance data.

    Attributes
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
            'HFSummary(\n'
            + f'    T = {T:.4} C,\n'
            + f'    P = {P:.4} kPa,\n'
            + f'    Pvap = {Pvap:.4} kPa,\n'
            + f'    ustar = {self.ustar:.4} m/s,\n'
            + f'    wind_w = {self.wind_w:.4} m/s,\n'
            + f'    var_w = {self.var_w:.4} (m/s)^2,\n'
            + f'    rho_vapor = {rho_vapor:.4} g/m^3,\n'
            + f'    rho_co2 = {rho_co2:.4} mg/m^3,\n'
            + f'    var_vapor = {var_vapor:.4} (g/m^3)^2,\n'
            + f'    var_co2 = {var_co2:.4} (mg/m^3)^2,\n'
            + f'    corr_q_c = {self.corr_q_c:.4},\n'
            + f'    cov_w_q = {cov_w_q:.4} g/m^2/s,\n'
            + f'    H = {self.H:.4} W/m^2,\n'
            + f'    cov_w_c = {cov_w_c:.4} mg/m^2/s,\n'
            + f'    rho_dryair = {self.rho_dryair:.4} kg/m^3,\n'
            + f'    rho_totair = {self.rho_totair:.4} kg/m^3,\n'
            + f'    cov_w_T = {self.cov_w_T:.4} C m/s,\n'
            + f'    N = {self.N})'
            )


@attr.s
class RootSoln(object):
    """Results from calcuating the root (corr_cp_cr, var_cp).

    Attributes
    ----------
    corr_cp_cr : float
        Correlation coefficient for CO2 concentrations connected with
        photosynthesis (cp) and respiration (cr).
    var_cp : float
        Variance of CO2 concentration connected with photosynthesis,
        (kg/m^3)^2.
    sig_cr : float
        Standard deviation of CO2 concentration connected with
        respiration, kg/m^3.
    co2soln_id : {0 or 1}
        Indicates the solution used for the quadratic Eq. 13b of
        [SS08]_. Equal to 1 for the '+' root, equal to 0 for the '-'
        root.
    isvalid : bool
        Indicates whether the obtained root (`corr_cp_cr`, `var_cp`) is
        physically plausible.
    mssg : str
        Possibly informative message if `isvalid` = False.

    """
    corr_cp_cr = attr.ib(default=np.nan)
    var_cp = attr.ib(default=np.nan)
    sig_cr = attr.ib(default=np.nan)
    co2soln_id = attr.ib(default=np.nan)
    isvalid = attr.ib(default=False)
    mssg = attr.ib(default="")

    def __str__(self):
        # For some fields, print common units instead of SI
        return (
            'RootSoln(\n'
            + '    corr_cp_cr = {:.4},\n'.format(self.corr_cp_cr)
            + '    var_cp = {:.4} (mg/m^3)^2,\n'.format(1e12 * self.var_cp)
            + '    sig_cr = {:.4} mg/m^3,\n'.format(1e6 * self.sig_cr)
            + '    co2soln_id = {},\n'.format(self.co2soln_id)
            + '    isvalid = {},\n'.format(self.isvalid)
            + '    mssg = {})'.format(self.mssg)
            )

#    wave_lvl : (int, int)
#        2-tuple indicating the level of filtering applied (number of
#        components removed from the series data). The second int is the
#        maximum possible wavelet decompostion level given the length of
#        the data. The first is the number of components remaining in the
#        data.  So when the first number is equal to the second, no
#        components have been removed (no filtering applied). When the
#        first number is 1, the maximum level of filtering was applied.
#   wave_lvl')):
#           + f'    wave_lvl = {})'

@attr.s
class WQCData(object):
    """Summary stats for wind, water vapor, and CO2 data.

    Attributes
    ----------
    var_q, var_c : float
        Variance of vapor (q) and CO2 (c) concentrations, (kg/m^3)^2.
    corr_qc : float
        Correlation coefficient for vapor and CO2 concentrations.
    wq, wc : float
        Mean vapor (`wq`) and CO2 (`wc`) fluxes, kg/m^2/s.

    """
    var_q = attr.ib(default=np.nan)
    var_c = attr.ib(default=np.nan)
    corr_qc = attr.ib(default=np.nan)
    wq = attr.ib(default=np.nan)
    wc = attr.ib(default=np.nan)

    def __str__(self):

        # For some fields, print common units instead of SI
        var_q=1e6 * self.var_q
        var_c=1e12 * self.var_c
        wq=1e3 * self.wq
        wc=1e6 * self.wc
        return (
            'WQCData(\n'
            + f'    var_q = {var_q:.4} (g/m^3)^2,\n'
            + f'    var_c = {var_c:.4} (mg/m^3)^2,\n'
            + f'    corr_qc = {self.corr_qc:.4},\n'
            + f'    wq = {wq:.4} g/m^2/s,\n'
            + f'    wc = {wc:.4} mg/m^2/s'
            )


@attr.s
class Outcome(object):
    """Overall outcome of partitioning.

    Attributes
    ----------
    version : str
        Fluxpart version
    dataread, attempt_partition, valid_partition : bool
        Indicates success or failure in reading high frequency data,
        attempting and obtaining a valid partioning solution.
    mssg : str
        Possibly informative message if `dataread` or `valid_partition`
        are False

    """
    version = attr.ib()
    dataread = attr.ib()
    attempt_partition = attr.ib()
    valid_partition = attr.ib()
    mssg = attr.ib()

    def __str__(self):
        return (
            'Outcome(\n'
            + f'    version = {self.version},\n'
            + f'    dataread = {self.dataread},\n'
            + f'    attempt_partition = {self.attempt_partition},\n'
            + f'    valid_partition = {self.valid_partition},\n'
            + f'    mssg = {self.mssg})'
            )
            

@attr.s
class WUE():
    """Summary of leaf-level water use efficiency calculation.

    Attributes
    ----------
    wue : float
        Leaf-level water use efficiency, kg CO2 / kg H2O.
    inter_h2o : float
        Concentration of intercellular water vapor, kg/m^3
    inter_co2 : float
        Concentrations of intercellular CO2, kg/m^3.
    ambient_h2o : float
        Concentrations of ambient atmospheric water vapor, kg/m^3.
    ambient_co2 : float
        Concentration of ambient CO2, kg/m^3.
    vpd : float
        Atmospheric vapor pressure deficit, Pa.
    ci_mod : str
        The name of the model used to estimate the intercellular CO2
        concentration.
    ci_mod_param : float or 2-tuple of floats
        Specific paramter values used with `ci_mod`.
    leaf_temper : float
    ppath : {'C3' or 'C4'}
        Photosynthetic pathway.
    meas_ht, canopy_ht : float
        Eddy covariance measument height and plant canopy height, m.

    """
    wue = attr.ib(default=np.nan)
    inter_h2o = attr.ib(default=np.nan)
    inter_co2 = attr.ib(default=np.nan)
    ambient_h2o = attr.ib(default=np.nan)
    ambient_co2 = attr.ib(default=np.nan)
    vpd = attr.ib(default=np.nan)
    ci_mod = attr.ib(default=np.nan)
    ci_mod_param = attr.ib(default=np.nan)
    leaf_temper = attr.ib(default=np.nan)
    ppath = attr.ib(default=np.nan)
    meas_ht = attr.ib(default=np.nan)
    canopy_ht = attr.ib(default=np.nan)
    diff_ratio = attr.ib(default=np.nan)

    def __str__(self):

        # For some fields, print common units instead of SI
        wue = 1e3 * self.wue
        inter_h2o = 1e3 * self.inter_h2o
        inter_co2 = 1e6 * self.inter_co2
        ambient_h2o = 1e3 * self.ambient_h2o
        ambient_co2 = 1e6 * self.ambient_co2
        vpd = 1e-3 * self.vpd
        leaf_temper = self.leaf_temper - 273.15

        return (
            'WUE(\n'
            + f'    wue = {wue:.4} mg/g,\n'
            + f'    inter_h2o = {inter_h2o:.4} g/m^3,\n'
            + f'    inter_co2 = {inter_co2:.4} mg/m^3,\n'
            + f'    ambient_h2o = {ambient_h2o:.4} g/m^3,\n'
            + f'    ambient_co2 = {ambient_co2:.4} mg/m^3,\n'
            + f'    vpd = {vpd:.4} kPa,\n'
            + f'    ci_mod = {self.ci_mod},\n'
            + f'    ci_mod_param = {self.ci_mod_param},\n'
            + f'    leaf_temper = {leaf_temper:.4} C,\n'
            + f'    ppath = {self.ppath},\n'
            + f'    meas_ht = {self.meas_ht} m,\n'
            + f'    canopy_ht = {self.canopy_ht} m,\n'
            + f'    diff_ratio = {self.diff_ratio})'
            )


@attr.s
class FluxpartResult(object):
    outcome = attr.ib()
    fvsp_result = attr.ib(
            default=FVSPResult(WQCData(), RootSoln(), AllFluxes()))
    hfsummary = attr.ib(default=HFSummary())
    wue = attr.ib(default=WUE())
    label = attr.ib(default=None)
