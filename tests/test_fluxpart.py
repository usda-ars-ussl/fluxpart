import os
from types import SimpleNamespace

import numpy.testing as npt

# from fluxpart import flux_partition
from fluxpart import fvs_partition

TESTDIR = os.path.dirname(os.path.realpath(__file__))


def test_flux_partition():
    """Integrated test of highest level flux_partition function.

    In test below, the "desired" results are not exact or even
    necessarily correct. They were obtained with an independent (matlab)
    code, and there should be reasonable agreement.
    """

    wue_data = {
        "meas_ht": 7.11,
        "canopy_ht": 4.42,
        "ppath": "C3",
        "ci_mod": "const_ppm",
        "diff_ratio": 1 / 0.7,
    }

    # soln exists for this data without any wavelet filtering
    fname = os.path.join(
        TESTDIR, "data/TOA5_6843.ts_Above_2012_06_07_1300.dat"
    )

    matlab_fluxes = SimpleNamespace(
        Fcp=-1.02807378762793,
        Fcr=0.402683101177712,
        Fqe=0.00500869036088203,
        Fqt=0.145615860044424,
        Fcp_mol=-23.3600042633021,
        Fcr_mol=9.14980916104776,
        LEe=12.1870591763138,
        LEt=354.310004313924,
    )

    fvsp = fvs_partition(fname, wue_options=wue_data, hfd_format="ec-TOA5")

    npt.assert_allclose(fvsp.df["fvsp_solution"]["var_cp"], 18.9272, atol=1)
    assert_flux_components(fvsp.df["fluxes"], matlab_fluxes)

    # test reading heights and temperature from file
    wue_data = {
        "heights": os.path.join(TESTDIR, "data/heights.csv"),
        "leaf_temper": os.path.join(TESTDIR, "data/leaf_temper.csv"),
        "ppath": "C3",
        "ci_mod": "const_ppm",
        "diff_ratio": 1 / 0.7,
    }

    fvsp = fvs_partition(fname, wue_options=wue_data, hfd_format="ec-TOA5")

    npt.assert_allclose(fvsp.df["fvsp_solution"]["var_cp"], 18.9272, atol=1)
    assert_flux_components(fvsp.df["fluxes"], matlab_fluxes)

    # soln is obtained after some wavelet filtering
    fname = os.path.join(
        TESTDIR, "data/TOA5_6843.ts_Above_2012_06_07_1245.dat"
    )

    wue_data = {
        "meas_ht": 7.11,
        "canopy_ht": 4.42,
        "ppath": "C3",
        "ci_mod": "const_ppm",
        "diff_ratio": 1 / 0.7,
    }

    matlab_fluxes = SimpleNamespace(
        Fcp=-0.866856083109642,
        Fcr=0.353428894620522,
        Fqe=0.0124697200158396,
        Fqt=0.117438136138301,
        Fcp_mol=-23.1074540435422,
        Fcr_mol=10.6590633820467,
        LEe=35.6007693518818,
        LEt=335.283229492226,
    )

    fvsp = fvs_partition(fname, wue_options=wue_data, hfd_format="ec-TOA5")

    npt.assert_allclose(
        fvsp.df.iloc[0]["fvsp_solution"]["var_cp"], 15.2944, atol=1
    )
    assert_flux_components(fvsp.df.iloc[0]["fluxes"], matlab_fluxes)


def assert_flux_components(calc, desired):
    npt.assert_allclose(calc.Fcp, desired.Fcp, atol=0.5)
    npt.assert_allclose(calc.Fcr, desired.Fcr, atol=0.1)
    npt.assert_allclose(calc.Fqe, desired.Fqe, atol=0.01)
    npt.assert_allclose(calc.Fqt, desired.Fqt, atol=0.05)
    npt.assert_allclose(calc.Fcp_mol, desired.Fcp_mol, atol=10)
    npt.assert_allclose(calc.Fcr_mol, desired.Fcr_mol, atol=10)
    npt.assert_allclose(calc.LEe, desired.LEe, atol=10)
    npt.assert_allclose(calc.LEt, desired.LEt, atol=50)


if __name__ == "__main__":
    test_flux_partition()
