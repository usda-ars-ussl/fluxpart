from types import SimpleNamespace
import numpy.testing as npt
from fluxpart.wue import water_use_efficiency


def test_water_use_efficiency():
    """ Unit test does not guarantee the wue function is correct, but if
    this fails, something has changed, and better understand why.
    """

    hf_stats = SimpleNamespace(
        rho_vapor=9.607e-3,
        rho_co2=658.8e-6,
        T=28.56 + 273.15,
        P=100.1e3,
        cov_w_q=0.1443e-3,
        cov_w_c=-1.059e-6,
        cov_w_T=0.1359,
        ustar=0.4179,
        rho_totair=1.150,
    )

    wue = water_use_efficiency(
        hf_stats,
        meas_ht=7.11,
        canopy_ht=4.42,
        ppath="C3",
        ci_mod="const_ppm",
        diff_ratio=(1 / 0.7),
    )

    npt.assert_allclose(wue.wue, -6.45e-3, atol=0.005e-3)
    npt.assert_allclose(wue.ambient_h2o, 12.4e-3, atol=0.05e-3)
    npt.assert_allclose(wue.ambient_co2, 638.e-6, atol=0.5e-6)
    npt.assert_allclose(wue.inter_h2o, 28.3e-3, atol=0.05e-3)
    npt.assert_allclose(wue.inter_co2, 492.e-6, atol=0.5e-6)


if __name__ == "__main__":
    test_water_use_efficiency()
