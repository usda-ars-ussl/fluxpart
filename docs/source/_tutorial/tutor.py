"""
/media/tower is link to Beltsville tower data
/media/data is link to Beltsville processed data
/media/results is link to Beltsville processed/fluxpart/version
BEGIN/END tags are for sphinx literalinclude directive
"""
import os

FIGDIR = "../images"


def tutor_fvsp():

    # BEGIN tutorfvsp
    from fluxpart import fvs_partition

    data_files = "/media/tower/*_2014_0[4-9]_*"
    results_file = "/media/results/apr-sept_2014.pkl"
    heights_file = "/media/data/meta_heights.csv"
    daylight_file = "/media/data/daylight2014.csv"

    fvsp = fvs_partition(
        file_or_dir=data_files,
        hfd_format="EC-TOA5",
        interval="30min",
        wue_options={"heights": heights_file, "ppath": "C3"},
        part_options={"daytime": daylight_file},
    )

    fvsp.save(results_file)
    # END tutorfvsp

    return fvsp


def tutor_plot():
    # need to sync this pathname def with the one used in tutorfvsp
    results_file = "/media/results/apr-sept_2014.pkl"
    from fluxpart import fpread
    fvsp = fpread(results_file)

    # BEGIN tutorplot
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    fvsp.plot_h2o(start="2014-09-01", end="2014-09-04", ax=ax)
    fig.autofmt_xdate()
    fig.savefig("/media/results/example_h2o.png")
    # END tutorplot
    fig.savefig(os.path.join(FIGDIR, "example_h2o.png"))


def main():
    fvsp = tutor_fvsp()
    tutor_plot()


if __name__ == "__main__":
    main()
