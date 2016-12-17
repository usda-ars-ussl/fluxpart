.. _ziparchive-example:

High-frequency data in ZIP archive
----------------------------------
In the Tutorial :ref:`example<tutorial-example>`, the high-frequency data files
were located in a file directory. The same analysis could be done if the files
were in a ZIP archive. Here, we assume the name of the archive is
my_archive.zip, and that it is located in ``/home/todd/data``. The wrapper
function ``my_partition`` from the Tutorial :ref:`example<tutorial-example>` is
assumed to be already defined. This code produces the same results as obtained
in the Tutorial :ref:`example<tutorial-example>`:

.. code-block:: python

    import datetime
    from itertools import product
    import os
    import zipfile

    fpout = []
    DATADIR = "/home/todd/data"
    zip_archive = 'my_archive.zip.zip'
    with zipfile.ZipFile(os.path.join(DATADIR, zip_archive), 'r') as arch:
        for dd, hh, mm in product(range(5, 9), range(0, 24), range(0, 60, 15)):
            filename = "TOA5_6843.ts_Above_2012_06_{:02}_{:02}{:02}.dat".format(dd, hh, mm)
            with arch.open(filename, 'r') as datafile:
                timestamp = datetime.datetime(2012, 6, dd, hh, mm)
                fpout.append(my_partition(datafile, timestamp))
