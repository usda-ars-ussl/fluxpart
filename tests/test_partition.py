from types import SimpleNamespace
import numpy.testing as npt
from fluxpart.partition import partition_from_qc_averages


def test_partition_from_qcdat():
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
    qcdat = SimpleNamespace(
        var_q=0.411163e-3 ** 2,
        var_c=5.182580e-6 ** 2,
        wq=0.033140e-3,
        wc=-0.472108e-6,
        corr_qc=-0.881017)
    # Results given in [PRV14]
    desired = SimpleNamespace(
        var_cp=5.230500e-6 ** 2,
        corr_cp_cr=-0.757626,
        wqt=0.012878e-3,
        wqe=0.020262e-3,
        wcp=-0.476489e-6,
        wcr=0.004381e-6)
    nsoln, fluxes = partition_from_qc_averages(qcdat, wue)
    assert_partition(nsoln, fluxes, desired)

    # April 5 "non-physical" example from [PRV14] Table 1
    wue = -24.558131e-3
    qcdat = SimpleNamespace(
        var_q=0.455994e-3 ** 2,
        var_c=4.544450e-6 ** 2,
        wq=0.062700e-3,
        wc=-0.712862e-6,
        corr_qc=-0.922292)
    # Resutls given in [PRV14]
    desired = SimpleNamespace(
        var_cp=3.951510e-6 ** 2,
        corr_cp_cr=-0.724706,
        wqt=0.025416e-3,
        wqe=0.037284e-3,
        wcp=-0.624172e-6,
        wcr=-0.088690e-6)
    nsoln, fluxes = partition_from_qc_averages(qcdat, wue)
    assert_partition(nsoln, fluxes, desired)

    # comparable to Ray Anderson peach data, 2012-06-07 1300
    # The valid solution in this case uses the '+' CO2 root
    wue = -7.060177e-3
    qcdat = SimpleNamespace(
        var_q=0.40639e-6,
        var_c=7.68505e-12,
        wq=0.1506337e-3,
        wc=-0.6254288e-6,
        corr_qc=-.9501656)
    nsoln, fluxes = partition_from_qc_averages(qcdat, wue)
    assert nsoln.success and nsoln.validroot

    # comparable to Ray Anderson peach data, 2012-06-07 0230
    # The valid solution in this case uses the '-' CO2 root
    wue = -68.77103e-3
    qcdat = SimpleNamespace(
        var_q=0.001326586e-6,
        var_c=9.948297e-12,
        wq=0.00088955655e-3,
        wc=0.07513186e-6,
        corr_qc=0.8886955)
    nsoln, fluxes = partition_from_qc_averages(qcdat, wue)
    assert nsoln.success and nsoln.validroot


def assert_partition(nsoln, fluxes, desired):
    if nsoln.success and nsoln.validroot:
        npt.assert_allclose(nsoln.var_cp, desired.var_cp, atol=1e-14)
        npt.assert_allclose(nsoln.corr_cp_cr, desired.corr_cp_cr, atol=1e-2)
        npt.assert_allclose(fluxes.wqt, desired.wqt, atol=1e-7)
        npt.assert_allclose(fluxes.wqe, desired.wqe, atol=1e-7)
        npt.assert_allclose(fluxes.wcp, desired.wcp, atol=1e-9)
        npt.assert_allclose(fluxes.wcr, desired.wcr, atol=1e-10)
    else:
        assert False


if __name__ == '__main__':
    test_partition_from_qcdat()
