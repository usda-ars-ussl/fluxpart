"""Classes used to group data and results."""

import attr
import numpy as np

from .util import qflux_mass_to_heat, qflux_mass_to_mol, cflux_mass_to_mol


@attr.s
class MassFluxes(object):
    """H2O and CO2 mass flux components.

    Parameters
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
            "MassFluxes(\n"
            "    Fq = {:.4} g/m^2/s,\n"
            "    Fqt = {:.4} g/m^2/s,\n"
            "    Fqe = {:.4} g/m^2/s,\n"
            "    Fc = {:.4} mg/m^2/s,\n"
            "    Fcp = {:.4} mg/m^2/s,\n"
            "    Fcr = {:.4} mg/m^2/s)"
            "".format(*(fqs + fcs))
        )


@attr.s
class AllFluxes(object):
    """Water vapor and CO2 fluxes.

    Parameters
    ----------
    Fq, Fqt, Fqe : float
        Total, transpiration, and evaporation H2O fluxes, kg/m^2/s.
    Fc, Fcp, Fcr : float
        Total, photosynthesis, and respiration CO2 fluxes, kg/m^2/s.
    temper_kelvin : float
        Temperature, K

    """

    Fq = attr.ib(default=np.nan)
    Fqt = attr.ib(default=np.nan)
    Fqe = attr.ib(default=np.nan)
    Fc = attr.ib(default=np.nan)
    Fcp = attr.ib(default=np.nan)
    Fcr = attr.ib(default=np.nan)
    temper_kelvin = attr.ib(default=np.nan)
    LE = attr.ib(default=False)
    LEt = attr.ib(default=False)
    LEe = attr.ib(default=False)
    Fq_mol = attr.ib(default=False)
    Fqt_mol = attr.ib(default=False)
    Fqe_mol = attr.ib(default=False)
    Fc_mol = attr.ib(default=False)
    Fcp_mol = attr.ib(default=False)
    Fcr_mol = attr.ib(default=False)

    def __attrs_post_init__(self):
        """Water vapor and CO2 fluxes.

        Derived Parameters
        ------------------
        LE, LEt, LEe : float
            Water vapor fluxes expressed as latent heat, W/m^2.
        Fq_mol, Fqt_mol, Fqe_mol : float
            Water vapor fluxes expressed as mol/m^2/s.
        Fc_mol, Fcp_mol, Fcr_mol : float
            CO2 fluxes expressed as mol/m^2/s.

        """
        if self.LE is False:
            self.LE = qflux_mass_to_heat(self.Fq, self.temper_kelvin)
        if self.LEt is False:
            self.LEt = qflux_mass_to_heat(self.Fqt, self.temper_kelvin)
        if self.LEe is False:
            self.LEe = qflux_mass_to_heat(self.Fqe, self.temper_kelvin)
        if self.Fq_mol is False:
            self.Fq_mol = qflux_mass_to_mol(self.Fq)
        if self.Fqt_mol is False:
            self.Fqt_mol = qflux_mass_to_mol(self.Fqt)
        if self.Fqe_mol is False:
            self.Fqe_mol = qflux_mass_to_mol(self.Fqe)
        if self.Fc_mol is False:
            self.Fc_mol = cflux_mass_to_mol(self.Fc)
        if self.Fcp_mol is False:
            self.Fcp_mol = cflux_mass_to_mol(self.Fcp)
        if self.Fcr_mol is False:
            self.Fcr_mol = cflux_mass_to_mol(self.Fcr)

    def __str__(self):
        # print common units instead of SI
        return self.results_str().format(**self.common_units())

    def results_str(self):
        lab = self.common_units_labels()
        return (
            "------\n"
            "Fluxes\n" + "------\n"
            "  Fq = {Fq:.4} " + lab["Fq"] + "\n"
            "  Fqt = {Fqt:.4} " + lab["Fqt"] + "\n"
            "  Fqe = {Fqe:.4} " + lab["Fqe"] + "\n"
            "  Fc = {Fc:.4} " + lab["Fc"] + "\n"
            "  Fcp = {Fcp:.4} " + lab["Fcp"] + "\n"
            "  Fcr = {Fcr:.4} " + lab["Fcr"] + "\n"
            "  Fq_mol = {Fq_mol:.4} " + lab["Fq_mol"] + "\n"
            "  Fqt_mol = {Fqt_mol:.4} " + lab["Fqt_mol"] + "\n"
            "  Fqe_mol = {Fqe_mol:.4} " + lab["Fqe_mol"] + "\n"
            "  Fc_mol = {Fc_mol:.4} " + lab["Fc_mol"] + "\n"
            "  Fcp_mol = {Fcp_mol:.4} " + lab["Fcp_mol"] + "\n"
            "  Fcr_mol = {Fcr_mol:.4} " + lab["Fcr_mol"] + "\n"
            "  LE  = {LE:.4} " + lab["LE"] + "\n"
            "  LEt = {LEt:.4} " + lab["LEt"] + "\n"
            "  LEe = {LEe:.4} " + lab["LEe"]
        )

    def common_units(self):
        return dict(
            Fq=1e3 * self.Fq,
            Fqt=1e3 * self.Fqt,
            Fqe=1e3 * self.Fqe,
            Fc=1e6 * self.Fc,
            Fcp=1e6 * self.Fcp,
            Fcr=1e6 * self.Fcr,
            Fq_mol=1e3 * self.Fq_mol,
            Fqt_mol=1e3 * self.Fqt_mol,
            Fqe_mol=1e3 * self.Fqe_mol,
            Fc_mol=1e6 * self.Fc_mol,
            Fcp_mol=1e6 * self.Fcp_mol,
            Fcr_mol=1e6 * self.Fcr_mol,
            LE=self.LE,
            LEt=self.LEt,
            LEe=self.LEe,
        )

    def common_units_labels(self):
        return dict(
            Fq="g/m^2/s",
            Fqt="g/m^2/s",
            Fqe="g/m^2/s",
            Fc="mg/m^2/s",
            Fcp="mg/m^2/s",
            Fcr="mg/m^2/s",
            Fq_mol="mmol/m^2/s",
            Fqt_mol="mmol/m^2/s",
            Fqe_mol="mmol/m^2/s",
            Fc_mol="umol/m^2/s",
            Fcp_mol="umol/m^2/s",
            Fcr_mol="umol/m^2/s",
            LE="W/m^2",
            LEt="W/m^2",
            LEe="W/m^2",
        )


@attr.s
class RootSoln(object):
    """Results from calculating the root (corr_cp_cr, var_cp).

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
    valid_root : bool
        Indicates whether the obtained root (`corr_cp_cr`, `var_cp`) is
        physically plausible.
    root_mssg : str
        Possibly informative message if `isvalid` = False.

    """

    corr_cp_cr = attr.ib(default=np.nan)
    var_cp = attr.ib(default=np.nan)
    sig_cr = attr.ib(default=np.nan)
    co2soln_id = attr.ib(default=np.nan)
    valid_root = attr.ib(default=np.nan)
    root_mssg = attr.ib(default="")

    def __str__(self, head=True):
        # Print common units instead of SI
        return self.results_str(head).format(**self.common_units())

    def results_str(self, head=True):
        lab = self.common_units_labels()
        out = ""
        if head is True:
            out += " ------------\nRoot Solution\n-------------\n"
        out += (
            "  corr_cp_cr = {corr_cp_cr:.4}\n"
            "  var_cp = {var_cp:.4} " + lab["var_cp"] + "\n"
            "  sig_cr = {sig_cr:.4} " + lab["sig_cr"] + "\n"
            "  co2soln_id = {co2soln_id:.0f}\n"
            "  valid_root = {valid_root}\n"
            "  root_mssg = {root_mssg}"
        )
        return out

    def common_units(self):
        return dict(
            corr_cp_cr=self.corr_cp_cr,
            var_cp=1e12 * self.var_cp,
            sig_cr=1e6 * self.sig_cr,
            co2soln_id=self.co2soln_id,
            valid_root=self.valid_root,
            root_mssg=self.root_mssg,
        )

    def common_units_labels(self):
        return dict(
            corr_cp_cr="",
            var_cp="(mg/m^3)^2",
            sig_cr="mg/m^3",
            co2soln_id="",
            valid_root="",
            root_mssg="",
        )


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

    def __str__(self, head=True):
        # Print common units instead of SI
        return self.results_str(head).format(**self.common_units())

    def results_str(self, head=True):
        lab = self.common_units_labels()
        out = ""
        if head:
            out += "--------------\nInterval Stats\n" "--------------\n"
        out += (
            "  var_q = {var_q:.4} " + lab["var_q"] + "\n"
            "  var_c = {var_c:.4} " + lab["var_c"] + "\n"
            "  corr_qc = {corr_qc:.4}\n"
            "  cov_w_q = {wq:.4} " + lab["cov_w_q"] + "\n"
            "  cov_w_c = {wc:.4} " + lab["cov_w_c"]
        )
        return out

    def common_units(self):
        return dict(
            var_q=1e6 * self.var_q,
            var_c=1e12 * self.var_c,
            corr_qc=self.corr_qc,
            wq=1e3 * self.wq,
            wc=1e6 * self.wc,
        )

    def common_units_labels(self):
        return dict(
            var_q="(g/m^3)^2",
            var_c="(mg/m^3)^2",
            corr_qc="",
            cov_w_q="g/m^2/s",
            cov_w_c="mg/m^2/s",
        )


@attr.s
class FVSPSolution(object):
    """Result of FVS partitioning.

    Parameters
    ----------
     wave_lvl : (int, int)
         2-tuple indicating the level of filtering applied (number of
         components removed from the series data). The second int is the
         maximum possible wavelet decompostion level given the length of
         the data. The first is the number of components remaining in the
         data.  So when the first number is equal to the second, no
         components have been removed (no filtering applied). When the
         first number is 1, the maximum level of filtering was applied.

    """

    wqc_data = attr.ib(default=WQCData())
    rootsoln = attr.ib(default=RootSoln())
    wave_lvl = attr.ib(default=np.nan)
    valid_partition = attr.ib(default=np.nan)
    fvsp_mssg = attr.ib(default="")

    def __str__(self):
        return self.results_str().format(
            fvsp_mssg=self.fvsp_mssg,
            valid_partition=self.valid_partition,
            wave_lvl=self.wave_lvl,
            **self.wqc_data.common_units(),
            **self.rootsoln.common_units(),
        )

    def results_str(self):
        result = (
            "-------------\n"
            "FVSP Solution\n"
            "-------------\n"
            "  valid_partition: {valid_partition}\n"
            "  fvsp_mssg: {fvsp_mssg}\n"
        )
        result += "  Wavelet filtering level: {wave_lvl}\n"
        result += self.wqc_data.results_str(head=False) + "\n"
        result += self.rootsoln.results_str(head=False)
        return result

    def common_units(self):
        return {
            "valid_partition": self.valid_partition,
            "fvsp_mssg": self.fvsp_mssg,
            "wave_lvl": self.wave_lvl,
            **self.wqc_data.common_units(),
            **self.rootsoln.common_units(),
        }


@attr.s
class WUE:
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
        # Print common units instead of SI
        return self.results_str().format(**self.common_units())

    def results_str(self):
        lab = self.common_units_labels()
        return (
            "--------------------\n"
            "Water Use Efficiency\n"
            "--------------------\n"
            "  wue = {wue:.4} " + lab["wue"] + "\n"
            "  inter_h2o = {inter_h2o:.4} " + lab["inter_h2o"] + "\n"
            "  inter_co2 = {inter_co2:.4} " + lab["inter_co2"] + "\n"
            "  ambient_h2o = {ambient_h2o:.4} " + lab["ambient_h2o"] + "\n"
            "  ambient_co2 = {ambient_co2:.4} " + lab["ambient_co2"] + "\n"
            "  vpd = {vpd:.4} " + lab["vpd"] + "\n"
            "  ci_mod = {ci_mod} " + lab["ci_mod"] + "\n"
            "  ci_mod_param = {ci_mod_param} " + lab["ci_mod_param"] + "\n"
            "  leaf_temper = {leaf_temper:.4} " + lab["leaf_temper"] + "\n"
            "  ppath = {ppath} " + lab["ppath"] + "\n"
            "  meas_ht = {meas_ht} " + lab["meas_ht"] + "\n"
            "  canopy_ht = {canopy_ht} " + lab["canopy_ht"] + "\n"
            "  diff_ratio = {diff_ratio:.4} " + lab["diff_ratio"]
        )

    def common_units(self):
        return dict(
            wue=1e3 * self.wue,
            inter_h2o=1e3 * self.inter_h2o,
            inter_co2=1e6 * self.inter_co2,
            ambient_h2o=1e3 * self.ambient_h2o,
            ambient_co2=1e6 * self.ambient_co2,
            vpd=1e-3 * self.vpd,
            ci_mod=self.ci_mod,
            ci_mod_param=self.ci_mod_param,
            leaf_temper=self.leaf_temper - 273.15,
            ppath=self.ppath,
            meas_ht=self.meas_ht,
            canopy_ht=self.canopy_ht,
            diff_ratio=self.diff_ratio,
        )

    def common_units_labels(self):
        return dict(
            wue="mg/g",
            inter_h2o="g/m^3",
            inter_co2="mg/m^3",
            ambient_h2o="g/m^3",
            ambient_co2="mg/m^3",
            vpd="kPa",
            ci_mod="",
            ci_mod_param="",  # depends on the mod
            leaf_temper="C",
            ppath="",
            meas_ht="m",
            canopy_ht="m",
            diff_ratio="",
        )
