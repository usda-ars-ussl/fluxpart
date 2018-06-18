from copy import deepcopy
import os

import attr
import numpy as np

from .__version__ import __version__
from .wue import water_use_efficiency, WUEError
from .hfdata import HFData, HFDataSource, HFDataReadError, TooFewDataError
from .partition import fvspart_progressive
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

PART_OPTIONS = dict(adjust_fluxes=True, sun=None)

_bad_ustar = "ustar = {:.4} is less than ustar_tol = {:.4}"
_bad_vpd = "Incompatible vapor pressure deficit (vpd = {} <= 0)"
_bad_qflux = "Fq = cov(w,q) = {:.4} <= 0 is incompatible with fvs partitioning"


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
                FluxpartResult(label=label, mssg=e.args[0], dataread=True)
            )
            continue

        if part_options["sun"]:
            sunrise, sunset = part_options["sun"]
            if callable(sunrise):
                sunrise = sunrise(date)
            if callable(sunset):
                sunset = sunset(date)
            if time < sunrise or time > sunset:
                fluxes = AllFluxes(_set_all_nonstomatal_fluxes(hfsum), hfsum.T)
                results.append(
                    FluxpartResult(
                        label=label,
                        valid_partition=True,
                        mssg="Night time; fluxes assumed to be non-stomatal",
                        dataread=True,
                        fluxes=fluxes,
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
                    label=label, mssg=e.args[0], dataread=True, hfsummary=hfsum
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
            fvsp.fluxes = None

        results.append(
            FluxpartResult(
                label=label,
                dataread=True,
                attempt_partition=True,
                fluxes=fluxes,
                valid_partition=fvsp.valid_partition,
                mssg=fvsp.mssg,
                fvsp_result=fvsp,
                hfsummary=hfsum,
                wue=leaf_wue,
            )
        )

    return results


def _set_all_fluxes_nonstomatal(hfsum):
    return MassFluxes(
        Fq=hfsum.cov_w_q,
        Fqt=0.,
        Fqe=hfsum.cov_w_q,
        Fc=hfsum.cov_w_c,
        Fcp=0.,
        Fcr=hfsum.cov_w_c,
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
        valid_partition=False,
        mssg=None,
        fluxes=None,
        label=None,
        hfsummary=None,
        wue=None,
        fvsp_result=None,
    ):
        """Fluxpart result.

        Parameters
        ----------
        dataread, attempt_partition, valid_partition : bool
            Indicates success or failure in reading high frequency data,
            attempting and obtaining a valid partioning solution.
        mssg : str
            Possibly informative message if `dataread` or `valid_partition`
            are False
        label : optional
            Optional id label. Could be, e.g. a datetime object or string.
        fvsp_result : :class:`~fluxpart.containers.FVSPResult`
        wue : :class:`~fluxpart.containers.WUE`
        hfsummary : :class:`~fluxpart.hfdata.HFSummary`

        """
        self.version = __version__
        self.dataread = dataread
        self.attempt_partition = attempt_partition
        self.valid_partition = valid_partition
        self.mssg = mssg
        self.fluxes = fluxes
        self.label = label
        self.fvsp_result = fvsp_result
        self.wue = wue
        self.hfsummary = hfsummary

    def __str__(self):
        result = (
            "===============\n"
            "Fluxpart Result\n"
            "===============\n"
            f"version = {self.version}\n"
            f"dataread = {self.dataread}\n"
            f"attempt_partition = {self.attempt_partition}\n"
            f"valid_partition = {self.valid_partition}\n"
            f"mssg = {self.mssg}\n"
        )
        if self.fvsp_result is not None:
            result += self.fvsp_result.__str__() + "\n"
        if self.fluxes is not None:
            result += self.fluxes.__str__() + "\n"
        if self.wue is not None:
            result += self.wue.__str__() + "\n"
        if self.hfsummary is not None:
            result += self.hfsummary.__str__()
        return result


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
