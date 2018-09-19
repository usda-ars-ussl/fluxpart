from types import SimpleNamespace
import numpy.testing as npt
from fluxpart.partition import fvspart_interval


def test_fvspart_interval():
    """
    References
    ----------
    [PRV14] L Palatella, G Rana, D Vitale (2014), Towards a flux-
    partitioning procedure based on the direct use of high-frequency
    eddy-covariance data.  Boundary-Layer Meteorol. 153:327--337,
    doi:10.1007/s10546-014-9947-x
    """

    # April 7 "physical" example from [PRV14] Table 1
    wue = -37.158598e-3
    qcdata = SimpleNamespace(
        var_q=0.411163e-3 ** 2,
        var_c=5.182580e-6 ** 2,
        wq=0.033140e-3,
        wc=-0.472108e-6,
        corr_qc=-0.881017,
    )
    # Results given in [PRV14]
    desired = SimpleNamespace(
        var_cp=5.230500e-6 ** 2,
        corr_cp_cr=-0.757626,
        Fqt=0.012878e-3,
        Fqe=0.020262e-3,
        Fcp=-0.476489e-6,
        Fcr=0.004381e-6,
    )
    massfluxes, fvsp = fvspart_interval(qcdata, wue)
    assert_partition(massfluxes, fvsp, desired)

    # April 5 "non-physical" example from [PRV14] Table 1
    wue = -24.558131e-3
    qcdata = SimpleNamespace(
        var_q=0.455994e-3 ** 2,
        var_c=4.544450e-6 ** 2,
        wq=0.062700e-3,
        wc=-0.712862e-6,
        corr_qc=-0.922292,
    )
    # Resutls given in [PRV14]
    desired = SimpleNamespace(
        var_cp=3.951510e-6 ** 2,
        corr_cp_cr=-0.724706,
        Fqt=0.025416e-3,
        Fqe=0.037284e-3,
        Fcp=-0.624172e-6,
        Fcr=-0.088690e-6,
    )
    mass_fluxes, fvsp = fvspart_interval(qcdata, wue)
    # ver. 0.2.0 no longer returns nonphysical soln
    if 0:
        assert_partition(mass_fluxes, fvsp, desired)
    assert not fvsp.rootsoln.valid_root

    # comparable to Ray Anderson peach data, 2012-06-07 1300
    # The valid solution in this case uses the '+' CO2 root
    wue = -7.060177e-3
    qcdata = SimpleNamespace(
        var_q=0.40639e-6,
        var_c=7.68505e-12,
        wq=0.1506337e-3,
        wc=-0.6254288e-6,
        corr_qc=-.9501656,
    )
    massfluxes, fvsp = fvspart_interval(qcdata, wue)
    assert fvsp.rootsoln.valid_root

    # TODO: Add a test using the '-' root.
    # Previous test used here was wrong


def assert_partition(fluxes, fvsp, desired):
    if fvsp.rootsoln.valid_root:
        npt.assert_allclose(fvsp.rootsoln.var_cp, desired.var_cp, atol=1e-14)
        npt.assert_allclose(
            fvsp.rootsoln.corr_cp_cr, desired.corr_cp_cr, atol=1e-2
        )
        npt.assert_allclose(fluxes.Fqt, desired.Fqt, atol=1e-7)
        npt.assert_allclose(fluxes.Fqe, desired.Fqe, atol=1e-7)
        npt.assert_allclose(fluxes.Fcp, desired.Fcp, atol=1e-9)
        npt.assert_allclose(fluxes.Fcr, desired.Fcr, atol=1e-10)
    else:
        assert False


if __name__ == "__main__":
    test_fvspart_interval()
