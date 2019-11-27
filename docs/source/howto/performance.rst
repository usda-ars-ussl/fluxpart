.. _performance-howto:

Execution times and data files
------------------------------
When processing many datafiles at once, it is usually better to pass a
time-sorted list of filenames rather than relying on **fluxpart** to time-sort
the data. If files and paths are named such that alphanumeric sorting of the
filenames equates to time-sorting, then it is recommended to pre-sort the
filenames and call `fvs_partition` with `time_sorted=True`. For example, if the
datafiles are:

.. code::

    ./data/path/TOA5_2019_02_01.dat
    ./data/path/TOA5_2019_02_02.dat
    ./data/path/TOA5_2019_02_03.dat
    ./data/path/TOA5_2019_02_04.dat
    [...]

then the files can be sorted and analyzed with:

.. code:: python

    from glob import glob

    datapath = "./data/path/*"
    files = sorted(glob(datapath))
    [...]

    fvsp = fvs_partition(
        file_or_dir=files,
        time_sorted=True,
        [...]
    )

You can also manually create a list of time-sorted files.

Another advantage of using `time_sorted=True` is that `fvs_partition` can be
called with compressed data files (gzip, bz2, xz, or a single file zip).
Currently, automatic time-sorting (`time_sorted=False`) of data does not work
if datafiles are compressed. For example, the above code would work with the
following data:

.. code::

    ./data/path/TOA5_2019_02_01.dat.gz
    ./data/path/TOA5_2019_02_02.dat.gz
    ./data/path/TOA5_2019_02_03.dat.gz
    ./data/path/TOA5_2019_02_04.dat.gz
    [...]

Using compressed data files requires on-the-fly decompression, which of course
is slower than working with uncompressed data. However, when processing large
amounts of data, the lower required disk space can be worth it. We've found
that when processing, say, multiple years of high-frequency data, storing data
in gziped daily files provides a reasonable trade-off between execution time
and the number and size of datafiles that have to be stored.
