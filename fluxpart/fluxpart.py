from copy import deepcopy
import datetime as pydatetime
from functools import lru_cache
from glob import iglob
import os
import pickle
import zipfile

import attr

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
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
from .constants import MOLECULAR_WEIGHT as MW


EC_TOA5 = {
    "filetype": "csv",
    "skiprows": 4,
    "time_col": 0,
    "cols": (2, 3, 4, 5, 6, 7, 8),
    "temper_unit": "C",
    "unit_convert": dict(q=1e-3, c=1e-6, P=1e3),
    "na_values": "NAN",
}

EC_TOB1 = {
    "filetype": "tob1",
    "time_col": 0,
    "cols": (3, 4, 5, 6, 7, 8, 9),
    "temper_unit": "C",
    "unit_convert": dict(q=1e-3, c=1e-6, P=1e3),
}

EC_GHG1 = {
    "filetype": "ghg",
    "sep": "\t",
    "cols": (11, 12, 13, 7, 8, 9, 10),
    "time_col": [5, 6],
    "unit_convert": dict(q=1e-3 * MW.vapor, c=1e-3 * MW.co2, P=1e3),
    "temper_unit": "C",
    "skiprows": 8,
    "na_values": "NAN",
    "to_datetime_kws": {"format": "%Y-%m-%d %H:%M:%S:%f"},
}

HFD_FORMAT = EC_TOA5

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

PART_OPTIONS = dict(adjust_fluxes=True)

_bad_ustar = "ustar = {:.4} <= ustar_tol = {:.4}"
_bad_vpd = "vpd = {:.4} Pa <= 0"
_bad_qflux = "Fq = cov(w,q) = {:.4} <= 0"
_night_mssg = "Nighttime, fluxes all non-stomatal"
_fp_result_str = (
    "===============\n"
    "Fluxpart Result\n"
    "===============\n"
    "fluxpart version = {version}\n"
    "date = {date}\n"
    "---------------\n"
    "dataread = {dataread}\n"
    "attempt_partition = {attempt_partition}\n"
    "partition_success = {partition_success}\n"
    "mssg = {mssg}\n"
    "label = {label}\n"
    "sunrise = {sunrise}\n"
    "sunset = {sunset}\n"
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
    time_sorted=False,
    interval=None,
    hfd_format=None,
    hfd_options=None,
    meas_wue=None,
    wue_options=None,
    part_options=None,
    label=None,
    stdout=True,
    verbose=True,
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
    elif type(hfd_format) is str and hfd_format.upper() == "EC-TOA5":
        hfd_format = deepcopy(EC_TOA5)
    elif type(hfd_format) is str and hfd_format.upper() == "EC-TOB1":
        hfd_format = deepcopy(EC_TOB1)
    elif type(hfd_format) is str and hfd_format.upper() == "EC-GHG1":
        hfd_format = deepcopy(EC_GHG1)
    else:
        hfd_format = deepcopy(hfd_format)
        _validate_hfd_format(hfd_format)

    hfd_options = {**HFD_OPTIONS, **(hfd_options or {})}
    wue_options = {**WUE_OPTIONS, **(wue_options or {})}
    part_options = {**PART_OPTIONS, **(part_options or {})}

    unit_convert = hfd_format.pop("unit_convert", {})
    converters = {k: _converter_func(v, 0.0) for k, v in unit_convert.items()}
    temper_unit = hfd_format.pop("temper_unit").upper()
    if temper_unit == "C" or temper_unit == "CELSIUS":
        converters["T"] = _converter_func(1.0, 273.15)
    hfd_format["converters"] = converters

    # TODO: these should work with file objects, not just str name
    heights, leaf_temper = None, None
    canopy_ht = wue_options.pop("canopy_ht", None)
    meas_ht = wue_options.pop("meas_ht", None)
    if "heights" in wue_options:
        heights = wue_options.pop("heights")
        if not callable(heights):
            heights = _lookup(heights, 0, 1, 2)
    if "leaf_temper" in wue_options:
        leaf_temper = wue_options.pop("leaf_temper")
        if type(leaf_temper) is str:
            leaf_temper = _lookup(leaf_temper, 0, 1)
    if "daytime" in part_options:
        if type(part_options["daytime"]) is str:
            part_options["daytime"] = _lookup(part_options["daytime"], 0, 1, 2)
    if meas_wue:
        if type(meas_wue) is str:
            meas_wue = _lookup(meas_wue, 0, 1)

    if stdout:
        print("Getting filenames ...")
    files = _files(file_or_dir)
    if len(files) == 0:
        print("No files found!")
        return
    if not time_sorted:
        if stdout:
            print("Reading {} file datetimes ...".format(len(files)))
        times = _peektime(files, **hfd_format)
        if stdout:
            print("Time sorting data files ...")
        sorted_files = [f for t, f in sorted(zip(times, files))]
    if stdout:
        print("Creating data source ...")

    reader = HFDataSource(sorted_files, **hfd_format).reader(interval=interval)
    results = []
    if stdout:
        print("Beginning partitioning analyses ...")

    while True:
        try:
            hfdat = HFData(next(reader))
        except HFDataReadError as e:
            results.append(FVSResult(label=label, mssg=e.args[0]))
            continue
        except StopIteration:
            break

        datetime = hfdat.dataframe.index[0]
        date = datetime.date()
        time = datetime.time()
        if stdout:
            print("Processing {}".format(datetime))

        sunrise, sunset = None, None
        nighttime = False
        if "daytime" in part_options:
            if callable(part_options["daytime"]):
                sunrise, sunset = part_options["daytime"](date)
            else:
                sunrise, sunset = part_options["daytime"]
            sunrise = pd.to_datetime(sunrise).time()
            sunset = pd.to_datetime(sunset).time()
            # shift sunrise so we partition if rise occurs during interval
            if interval:
                sunrise = pd.Timestamp.combine(date, sunrise)
                sunrise = (sunrise - pd.Timedelta(interval)).time()
            nighttime = time < sunrise or time > sunset

        try:
            hfdat, hfsum = _set_hfdata(hfdat, **hfd_options)
        except (TooFewDataError, FluxpartError) as e:
            results.append(
                FVSResult(
                    label=datetime,
                    mssg=e.args[0],
                    dataread=True,
                    attempt_partition=False,
                    partition_success=False,
                    sunrise=sunrise,
                    sunset=sunset,
                )
            )
            if stdout and verbose:
                print(e.args[0])
            continue

        if nighttime and hfsum.cov_w_c >= 0:
            fluxes = AllFluxes(
                temper_kelvin=hfsum.T, **_set_all_fluxes_nonstomatal(hfsum)
            )
            mssg = _night_mssg
            results.append(
                FVSResult(
                    label=datetime,
                    attempt_partition=False,
                    partition_success=True,
                    mssg=mssg,
                    dataread=True,
                    fluxes=fluxes,
                    hfsummary=hfsum,
                    sunrise=sunrise,
                    sunset=sunset,
                )
            )
            if stdout and verbose:
                print(mssg)
            continue

        try:
            if meas_wue:
                if callable(meas_wue):
                    leaf_wue = WUE(wue=meas_wue(datetime))
                else:
                    leaf_wue = WUE(wue=float(meas_wue))
            else:
                if heights is not None:
                    if callable(heights):
                        canopy_ht, meas_ht = heights(date)
                    else:
                        canopy_ht, meas_ht = heights
                else:
                    if callable(canopy_ht):
                        canopy_ht = canopy_ht(date)
                    if callable(meas_ht):
                        meas_ht = meas_ht(date)
                leaf_t = None
                if leaf_temper is not None:
                    if callable(leaf_temper):
                        leaf_t = leaf_temper(datetime)
                    else:
                        leaf_t = float(leaf_temper)
                    if temper_unit == "C" or temper_unit == "CELSIUS":
                        leaf_t = leaf_t + 273.15
                leaf_wue = water_use_efficiency(
                    hfsum,
                    canopy_ht=canopy_ht,
                    meas_ht=meas_ht,
                    leaf_temper=leaf_t,
                    **wue_options,
                )

        except WUEError as e:
            results.append(
                FVSResult(
                    label=datetime,
                    mssg=e.args[0],
                    dataread=True,
                    hfsummary=hfsum,
                    sunrise=sunrise,
                    sunset=sunset,
                )
            )
            if stdout and verbose:
                print(e.args[0])
            continue

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
            fluxes = AllFluxes()

        results.append(
            FVSResult(
                label=datetime,
                sunrise=sunrise,
                sunset=sunset,
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
        if stdout and verbose:
            if fvsp.mssg:
                print(fvsp.mssg)

    return FluxpartResult(results)


def _set_all_fluxes_nonstomatal(hfsum):
    return dict(
        Fq=hfsum.cov_w_q,
        Fqt=0.0,
        Fqe=hfsum.cov_w_q,
        Fc=hfsum.cov_w_c,
        Fcp=0.0,
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


class FVSResult(object):
    """FVS partitioning result."""

    def __init__(
        self,
        dataread=False,
        attempt_partition=False,
        partition_success=False,
        mssg=None,
        label=None,
        sunrise=None,
        sunset=None,
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
        self.sunrise = sunrise
        self.sunset = sunset
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
            sunrise=self.sunrise,
            sunset=self.sunset,
            **attr.asdict(self.fluxes),
            **attr.asdict(self.hfsummary),
            **attr.asdict(self.wue),
            **fluxpart,
            **wqc_data,
            **rootsoln,
        )


class FluxpartResult(object):
    def __init__(self, fp_results):
        if type(fp_results) is str:
            with open(fp_results, "rb") as f:
                self.df = pd.read_pickle(f)
                self.meta = pickle.load(f)
            return
        index = pd.DatetimeIndex(r.label for r in fp_results)
        df0 = pd.DataFrame(
            (r.fluxes.common_units() for r in fp_results),
            index=index,
            columns=fp_results[0].fluxes.common_units().keys(),
        )
        df1 = pd.DataFrame(
            (r.hfsummary.common_units() for r in fp_results),
            index=index,
            columns=fp_results[0].hfsummary.common_units().keys(),
        )
        df2 = pd.DataFrame(
            (r.wue.common_units() for r in fp_results),
            index=index,
            columns=fp_results[0].wue.common_units().keys(),
        )
        df3 = pd.DataFrame(
            (r.fvsp_solution.common_units() for r in fp_results),
            index=index,
            columns=fp_results[0].fvsp_solution.common_units().keys(),
        )
        df4 = pd.DataFrame(
            {
                "dataread": [r.dataread for r in fp_results],
                "attempt_partition": [r.attempt_partition for r in fp_results],
                "partition_success": [r.partition_success for r in fp_results],
                "mssg": [r.mssg for r in fp_results],
                "sunrise": [r.sunrise for r in fp_results],
                "sunset": [r.sunset for r in fp_results],
            },
            index=index,
        )
        self.df = pd.concat(
            [df0, df1, df2, df3, df4],
            axis=1,
            sort=False,
            keys=["fluxes", "hfsummary", "wue", "fvsp_solution", "fluxpart"],
        )

        self.meta = {
            "version": fp_results[0].version,
            "date": str(pydatetime.datetime.now()),
        }

    def __str__(self):
        if len(self.df) == 1:
            return self.istr(0)
        else:
            return self.df.__str__()

    def __getitem__(self, item):
        return self.df[item]

    def __getattr__(self, x):
        return getattr(self.df, x)

    def plot_co2(
        self,
        start=None,
        end=None,
        units="mass",
        components=(0, 1, 2),
        ax=None,
        **kws,
    ):
        if ax is None:
            ax = plt.gca()
        if units == "mass":
            cols = ["Fc", "Fcp", "Fcr"]
            ylab = r"$\mathrm{CO_2\ Flux\ (mg\ m^{-2}\ s^{-1})}$"
        else:
            cols = ["Fc_mol", "Fcp_mol", "Fcr_mol"]
            ylab = r"$\mathrm{CO_2\ Flux\ (umol\ m^{-2}\ s^{-1})}$"
        labels = [
            r"$\mathrm{F_c}$",
            r"$\mathrm{F_{c_p}}$",
            r"$\mathrm{F_{c_r}}$",
        ]
        cols = [cols[j] for j in components]
        labels = [labels[j] for j in components]
        self.df.loc[start:end, ("fluxes", cols)].plot(ax=ax)
        ax.legend(labels)
        ax.set_ylabel(ylab)
        return ax

    def plot_h2o(
        self,
        start=None,
        end=None,
        units="mass",
        components=(0, 1, 2),
        ax=None,
        **kws,
    ):
        if ax is None:
            ax = plt.gca()
        if units == "mass":
            cols = ["Fq", "Fqt", "Fqe"]
            ylab = r"$\mathrm{H_20\ Flux\ (g\ m^{-2}\ s^{-1})}$"
        elif units == "mol":
            cols = ["Fq_mol", "Fqt_mol", "Fqe_mol"]
            ylab = r"$\mathrm{H_20\ Flux\ (mmol\ m^{-2}\ s^{-1})}$"
        else:
            cols = ["LE", "LEt", "LEe"]
            ylab = r"$\mathrm{LE\ (W\ m^{-2})}$"
        labels = [
            r"$\mathrm{F_q}$",
            r"$\mathrm{F_{q_t}}$",
            r"$\mathrm{F_{q_e}}$",
        ]

        cols = [cols[j] for j in components]
        labels = [labels[j] for j in components]

        self.df.loc[start:end, ("fluxes", cols)].plot(ax=ax)
        ax.legend(labels)
        ax.set_ylabel(ylab)
        return ax

    def istr(self, i):
        """Return a string representation of the ith result"""
        return _fp_result_str.format(
            version=self.meta["version"],
            date=self.meta["date"],
            label=self.df.index[i],
            **self.df.iloc[i]["fluxpart"].to_dict(),
            **self.df.iloc[i]["fluxes"].to_dict(),
            **self.df.iloc[i]["fvsp_solution"].to_dict(),
            **self.df.iloc[i]["hfsummary"].to_dict(),
            **self.df.iloc[i]["wue"].to_dict(),
        )

    def save(self, filename):
        with open(filename, "wb") as f:
            self.df.to_pickle(f)
            pickle.dump(self.meta, f)


def _converter_func(slope, intercept):
    """Return a function for linear transform of data."""

    def func(val):
        return slope * val + intercept

    return func


def _files(file_or_dir):
    if type(file_or_dir) is str:
        file_or_dir = [file_or_dir]
    unsorted_files = []
    for path in file_or_dir:
        if os.path.isfile(path):
            unsorted_files.append(path)
            continue
        if os.path.isdir(path):
            path = os.path.join(path, "*")
        unsorted_files += iglob(path)
    return unsorted_files


def _peektime(files, **kwargs):
    if kwargs["filetype"] == "csv" or kwargs["filetype"] == "ghg":
        dtcols = kwargs["time_col"]
        if type(dtcols) is int:
            dtcols = [dtcols]
        sep = ","
        if "delimiter" in kwargs:
            sep = kwargs["delimiter"]
        if "sep" in kwargs:
            sep = kwargs["sep"]
        datetimes = []
        to_datetime_kws = kwargs.get("to_datetime_kws", {})
    if kwargs["filetype"] == "csv":
        for file_ in files:
            with open(file_, "rt") as f:
                for _ in range(kwargs["skiprows"]):
                    f.readline()
                row = f.readline().split(sep)
                tstamp = " ".join([row[i].strip("'\"") for i in dtcols])
                datetimes.append(pd.to_datetime(tstamp, **to_datetime_kws))
    elif kwargs["filetype"] == "ghg":
        for file_ in files:
            with zipfile.ZipFile(file_) as z:
                with z.open(os.path.basename(file_)[:-3] + "data", "r") as f:
                    for _ in range(kwargs["skiprows"]):
                        f.readline()
                    row = f.readline().decode("utf-8").split(sep)
                    tstamp = " ".join([row[i].strip("'\"") for i in dtcols])
                    datetimes.append(pd.to_datetime(tstamp, **to_datetime_kws))
    else:  # "tob1"
        source = HFDataSource(files, count=5, **kwargs)
        datetimes = [df.index[0] for df in source.reader(interval=None)]
    return datetimes


def _validate_hfd_format(hfd_format):
    if "cols" not in hfd_format:
        raise Error("No value for hfd_format['cols'] given.")
    if "filetype" not in hfd_format:
        raise Error("No value for hfd_format['filetype'] given.")
    if hfd_format["filetype"] not in ("csv", "tob1", "ghg"):
        raise Error(f"Unrecognized filetype: {hfd_format['filetype']}")


def _lookup(csv_file, date_icol, icol1, icol2=None, method="ffill"):
    """Create a function for looking up data in csv file.
    date_icol, icol1, icol2 : int
        column index for the respective data
    method : str
        Interpolation method used with pandas df.index.get_loc. The
        default 'ffill' returns the PREVIOUS values if no exact date
        match is found in the lookup.

    """
    df = pd.read_csv(csv_file, index_col=date_icol, parse_dates=True)

    @lru_cache()
    def func(date):
        ix = df.index.get_loc(pd.to_datetime(date), method=method)
        if icol2 is None:
            return df.iloc[ix, icol1 - 1]
        else:
            return df.iloc[ix, icol1 - 1], df.iloc[ix, icol2 - 1]

    return func
