from copy import deepcopy
import datetime as pydatetime
from glob import glob
import os

import attr

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .__version__ import __version__
from .wue import water_use_efficiency, WUEError
from .hfdata import (
    HFData,
    HFSummary,
    HFDataSource,
    HFDataReadError,
    TooFewDataError,
)
from .partition import fvspart_progressive, FVSPSolution
from .util import vapor_press_deficit
from .containers import AllFluxes, WUE


EC_TOA5 = {
    "filetype": "csv",
    "skiprows": 4,
    "cols": (2, 3, 4, 5, 6, 7, 8),
    "time_col": 0,
    "temper_unit": "C",
    "unit_convert": dict(q=1e-3, c=1e-6, P=1e3),
}

EC_TOB1 = {
    "filetype": "tob1",
    "cols": (3, 4, 5, 6, 7, 8, 9),
    "temper_unit": "C",
    "unit_convert": dict(q=1e-3, c=1e-6, P=1e3),
}

HFD_FORMAT = EC_TOB1

WUE_OPTIONS = {
    "ci_mod": "const_ratio",
    "ci_mod_param": None,
    "leaf_temper": None,
    "leaf_temper_corr": 0,
    "diff_ratio": 1.6,
}

HFD_OPTIONS = {
    "bounds": {"c": (0, np.inf), "q": (0, np.inf)},
    "rd_tol": 0.5,
    "ad_tol": 1024,
    "ustar_tol": 0.1,
    "correct_external": True,
}

PART_OPTIONS = dict(adjust_fluxes=True, sunrise=None, sunset=None)

_bad_ustar = "ustar = {:.4} is less than ustar_tol = {:.4}"
_bad_vpd = "Incompatible vapor pressure deficit (vpd = {} <= 0)"
_bad_qflux = "Fq = cov(w,q) = {:.4} <= 0 is incompatible with fvs partitioning"
_fp_result_str = (
    "===============\n"
    "Fluxpart Result\n"
    "===============\n"
    "fluxpart version = {version}\n"
    "time = {timenow}\n"
    "---------------\n"
    "dataread = {dataread}\n"
    "attempt_partition = {attempt_partition}\n"
    "partition_success = {partition_success}\n"
    "mssg = {mssg}\n"
    "label = {label}\n"
    + AllFluxes().results_str()
    + "\n"
    + HFSummary().results_str()
    + "\n"
    + WUE().results_str()
    + "\n"
    + FVSPSolution().results_str()
)


class Error(Exception):
    pass


class FluxpartError(Error):
    def __init__(self, message):
        self.message = message


def fvspart(
    file_or_dir,
    ext="",
    interval="30min",
    hfd_format=None,
    hfd_options=None,
    meas_wue=None,
    wue_options=None,
    part_options=None,
    label=None,
):
    """Partition CO2 & H2O fluxes into stomatal & nonstomatal components.

    Provides a full implementation of the flux variance similarity
    partitioning algorithm [SS08]_[SAAS+18]_: reads high frequency eddy
    covariance data; performs data transformations and data QA/QC;
    analyzes water vapor and carbon dioxide fluxes; and partitions the
    fluxes into stomatal (transpiration, photosynthesis) and nonstomatal
    (evaporation, respiration) components.

    The following notation is used in variable naming and documentation
    to represent meteorological quantities::

        u, v, w = wind velocities
        q = water vapor mass concentration
        c = carbon dioxide mass concentration
        T = air temperature
        P = total air pressure

    Parameters
    ----------
    For parameters explanation see: :func:`~fluxpart.api.fvs_partition`

    Returns
    -------
    :class:`~fluxpart.fluxpart.FluxpartResult`

    """
    if hfd_format is None:
        hfd_format = deepcopy(HFD_FORMAT)
    elif hfd_format.upper() == "EC-TOA5":
        hfd_format = deepcopy(EC_TOA5)
    elif hfd_format.upper() == "EC-TOB1":
        hfd_format = deepcopy(EC_TOB1)
    else:
        hfd_format = deepcopy(hfd_format)
        _validate_hfd_format(hfd_format)

    hfd_options = {**HFD_OPTIONS, **(hfd_options or {})}
    wue_options = {**WUE_OPTIONS, **(wue_options or {})}
    part_options = {**PART_OPTIONS, **(part_options or {})}

    unit_convert = hfd_format.pop("unit_convert", {})
    converters = {k: _converter_func(v, 0.) for k, v in unit_convert.items()}
    temper_unit = hfd_format.pop("temper_unit")
    if temper_unit.upper() == "C" or temper_unit.upper() == "CELSIUS":
        converters["T"] = _converter_func(1., 273.15)
    hfd_format["converters"] = converters

    files = _files(file_or_dir, ext)
    times = _peektime(files, **hfd_format)
    sorted_files = [f for t, f in sorted(zip(times, files))]

    reader = HFDataSource(sorted_files, **hfd_format).reader(interval=interval)
    results = []

    while True:
        try:
            hfdat = HFData(next(reader))
        except HFDataReadError as e:
            results.append(FluxpartResult(label=label, mssg=e.args[0]))
            continue
        except StopIteration:
            break

        datetime = hfdat.dataframe.index[0]
        date = datetime.date()
        time = datetime.time()

        try:
            hfdat, hfsum = _set_hfdata(hfdat, **hfd_options)
        except TooFewDataError as e:
            results.append(
                FluxpartResult(
                    label=datetime,
                    mssg=e.args[0],
                    dataread=True,
                    partition_success="NA",
                )
            )
            continue

        sunrise, sunset = pydatetime.time.min, pydatetime.time.max
        if part_options["sunrise"]:
            sunrise = part_options["sunrise"]
            if callable(sunrise):
                sunrise = sunrise(date)
            sunrise = pd.to_datetime(sunrise)
            # shift so that we partition if sunrise occurs during interval
            sunrise = (sunrise - pd.Timedelta(interval)).time()
        if part_options["sunset"]:
            sunset = part_options["sunset"]
            if callable(sunset):
                sunset = sunset(date)
            sunset = pd.to_datetime(sunset).time()
        if (time < sunrise or time > sunset) and hfsum.cov_w_c >= 0:
            fluxes = AllFluxes(
                temper_kelvin=hfsum.T, **_set_all_fluxes_nonstomatal(hfsum)
            )
            results.append(
                FluxpartResult(
                    label=datetime,
                    partition_success="NA",
                    mssg="Non-stomatal fluxes assumed due to time of day",
                    dataread=True,
                    fluxes=fluxes,
                    hfsummary=hfsum,
                )
            )
            continue

        try:
            if meas_wue:
                leaf_wue = WUE(wue=float(meas_wue))
            else:
                leaf_wue = water_use_efficiency(
                    hfsum, date=date, **wue_options
                )
        except WUEError as e:
            results.append(
                FluxpartResult(
                    label=datetime,
                    mssg=e.args[0],
                    dataread=True,
                    hfsummary=hfsum,
                    fluxes=AllFluxes(
                        **_set_only_total_fluxes(hfsum), temper_kelvin=hfsum.T
                    ),
                )
            )

        fluxes, fvsp = fvspart_progressive(
            hfdat["w"].values,
            hfdat["q"].values,
            hfdat["c"].values,
            leaf_wue.wue,
            part_options["adjust_fluxes"],
        )

        if fvsp.valid_partition:
            fluxes = AllFluxes(**attr.asdict(fluxes), temper_kelvin=hfsum.T)
        else:
            fluxes = AllFluxes(
                **_set_only_total_fluxes(hfsum), temper_kelvin=hfsum.T
            )

        results.append(
            FluxpartResult(
                label=datetime,
                dataread=True,
                attempt_partition=True,
                fluxes=fluxes,
                partition_success=fvsp.valid_partition,
                mssg=fvsp.mssg,
                fvsp_solution=fvsp,
                hfsummary=hfsum,
                wue=leaf_wue,
            )
        )

    return FluxpartResultSeries(results)


def _set_all_fluxes_nonstomatal(hfsum):
    return dict(
        Fq=hfsum.cov_w_q,
        Fqt=0.,
        Fqe=hfsum.cov_w_q,
        Fc=hfsum.cov_w_c,
        Fcp=0.,
        Fcr=hfsum.cov_w_c,
    )


def _set_only_total_fluxes(hfsum):
    return dict(
        Fq=hfsum.cov_w_q,
        Fqt=np.nan,
        Fqe=np.nan,
        Fc=hfsum.cov_w_c,
        Fcp=np.nan,
        Fcr=np.nan,
    )


def _set_hfdata(hfdata, bounds, rd_tol, ad_tol, correct_external, ustar_tol):
    hfdata.cleanse(bounds, rd_tol, ad_tol)
    hfdata.truncate_pow2()
    if correct_external:
        hfdata.correct_external()
    hfsum = hfdata.summarize()
    if hfsum.ustar < ustar_tol:
        raise FluxpartError(_bad_ustar.format(hfsum.ustar, ustar_tol))
    vpd = vapor_press_deficit(hfsum.rho_vapor, hfsum.T)
    if vpd <= 0:
        raise FluxpartError(_bad_vpd.format(vpd))
    if hfsum.cov_w_q <= 0:
        raise FluxpartError(_bad_qflux.format(hfsum.cov_w_q))
    return hfdata, hfsum


def flux_partition(*args, **kws):
    return fvspart(*args, **kws)


class FluxpartResult(object):
    """Fluxpart result."""

    def __init__(
        self,
        dataread=False,
        attempt_partition=False,
        partition_success=False,
        mssg=None,
        label=None,
        fluxes=AllFluxes(),
        hfsummary=HFSummary(),
        wue=WUE(),
        fvsp_solution=FVSPSolution(),
    ):
        """Fluxpart result.

        Parameters
        ----------
        dataread, attempt_partition, partition_success : bool
            Indicates success or failure in reading high frequency data,
            attempting and obtaining a valid partioning solution.
        mssg : str
            Possibly informative message if `dataread` or `partition_success`
            are False
        label : optional
            Pandas datetime.
        fluxes : :class:`~fluxpart.containers.AllFluxes`
        fvsp_solution : :class:`~fluxpart.containers.FVSPResult`
        wue : :class:`~fluxpart.containers.WUE`
        hfsummary : :class:`~fluxpart.hfdata.HFSummary`

        """
        self.version = __version__
        self.dataread = dataread
        self.attempt_partition = attempt_partition
        self.partition_success = partition_success
        self.mssg = mssg
        self.fluxes = fluxes
        self.label = label
        self.fvsp_solution = fvsp_solution
        self.wue = wue
        self.hfsummary = hfsummary

    def __str__(self):
        fluxpart = attr.asdict(self.fvsp_solution)
        wqc_data = fluxpart.pop("wqc_data")
        rootsoln = fluxpart.pop("rootsoln")
        return _fp_result_str.format(
            timenow=pydatetime.datetime.now(),
            version=self.version,
            dataread=self.dataread,
            attempt_partition=self.attempt_partition,
            partition_success=self.partition_success,
            mssg=self.mssg,
            label=self.label,
            **attr.asdict(self.fluxes),
            **attr.asdict(self.hfsummary),
            **attr.asdict(self.wue),
            **fluxpart,
            **wqc_data,
            **rootsoln,
        )


class FluxpartResultSeries(object):
    def __init__(self, results):
        index = pd.DatetimeIndex(r.label for r in results)
        df0 = pd.DataFrame(
            (r.fluxes.common_units() for r in results), index=index
        )
        df1 = pd.DataFrame(
            (r.hfsummary.common_units() for r in results), index=index
        )
        df2 = pd.DataFrame(
            (r.wue.common_units() for r in results), index=index
        )
        df3 = pd.DataFrame(
            (r.fvsp_solution.common_units() for r in results), index=index
        )
        df4 = pd.DataFrame(
            {
                "version": [r.version for r in results],
                "dataread": [r.dataread for r in results],
                "attempt_partition": [r.attempt_partition for r in results],
                "partition_success": [r.partition_success for r in results],
                "mssg": [r.mssg for r in results],
            },
            index=index,
        )
        self.df = pd.concat(
            [df0, df1, df2, df3, df4],
            axis=1,
            keys=["fluxes", "hfsummary", "wue", "fvsp_solution", "fluxpart"],
        )

    def plot_Fc(self, **kws):
        ax = self.df["fluxes"][["Fc", "Fcp", "Fcr"]].plot(**kws)
        ax.set_ylabel(r"CO2 $[\mathrm{mg/m^2/s}]$")
        return ax

    def plot_Fc_mol(self, **kws):
        ax = self.df["fluxes"][["Fc_mol", "Fcp_mol", "Fcr_mol"]].plot(**kws)
        ax.set_ylabel(r"CO2 $[\mathrm{umol/m^2/s}]$")
        ax.legend(["Fc", "Fcp", "Fcr"])
        return ax

    def plot_Fq(self, **kws):
        ax = self.df["fluxes"][["Fq", "Fqe"]].plot(**kws)
        ax.set_ylabel(r"H2O $[\mathrm{g/m^2/s}]$")
        ax.legend(["Evapotranspiration", "Evaporation"])
        return ax

    def plot_Fq_mol(self, **kws):
        ax = self.df["fluxes"][["Fq_mol", "Fqe_mol"]].plot(**kws)
        ax.set_ylabel(r"H2O $[\mathrm{mmol/m^2/s}]$")
        ax.legend(["Evapotranspiration", "Evaporation"])
        return ax

    def plot_LE(self, **kws):
        ax = self.df["fluxes"][["LE", "LEe"]].plot(**kws)
        ax.set_ylabel(r"LE $[\mathrm{W/m^2}]$")
        ax.legend(["LE", "LE-evap"])
        return ax

    def iresult(self, i):
        """Return a str represenation of the ith result"""
        return _fp_result_str.format(
            timenow=pydatetime.datetime.now(),
            label=self.df.index[i],
            **self.df.iloc[i]["fluxpart"].to_dict(),
            **self.df.iloc[i]["fluxes"].to_dict(),
            **self.df.iloc[i]["fvsp_solution"].to_dict(),
            **self.df.iloc[i]["hfsummary"].to_dict(),
            **self.df.iloc[i]["wue"].to_dict(),
        )

    def tsplot(self, indx, **kws):
        ax = self.df.loc[:, indx].plot(**kws)
        return ax


def _converter_func(slope, intercept):
    """Return a function for linear transform of data."""

    def func(val):
        return slope * val + intercept

    return func


def _files(file_or_dir, ext=""):
    if isinstance(file_or_dir, str):
        file_or_dir = [file_or_dir]
    if os.path.isfile(file_or_dir[0]):
        unsorted_files = file_or_dir
    elif os.path.isdir(file_or_dir[0]):
        unsorted_files = []
        for p in file_or_dir:
            unsorted_files += glob(os.path.join(p, "*" + ext))
    else:
        raise FluxpartError("Unable to infer data type")
    return unsorted_files


def _peektime(files, **kwargs):
    kws = deepcopy(kwargs)
    if kws["filetype"] == "csv":
        kws["nrows"] = 1
    else:  # "tob1"
        kws = {"count": 5}
    source = HFDataSource(files, **kws)
    return [df.index[0] for df in source.reader(interval=None)]


def _validate_hfd_format(hfd_format):
    if "cols" not in hfd_format:
        raise Error("No value for hfd_format['cols'] given.")
    if "filetype" not in hfd_format:
        raise Error("No value for hfd_format['filetype'] given.")
    if hfd_format["filetype"] not in ("csv", "tob1"):
        raise Error(f"Unrecognized filetype: {hfd_format['filetype']}")
