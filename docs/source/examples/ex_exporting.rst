.. _saving-example:

.. _pickle: https://docs.python.org/3/library/pickle.html

Exporting selected results as csv
---------------------------------
If ``fpout`` contains :func:`~fluxpart.fluxpart.flux_partition()` output (or is
list of outputs, such as in this :ref:`example<tutorial-example>`), selected
results can be exported by looping over the output and writing the data to
a file.  For example:: 

    # write partitioned vapor fluxes to comma delimited text file
    with open("qfluxes.csv", 'w') as f:
        for out in fpout:
            f.write("{},{},{}\n"
            "".format(out['label'], out['fluxes'].Fqe, out['fluxes'].Fqt))

This generates the file `qfluxes.csv`:

.. code-block:: none

    2012-06-05 00:00:00,7.08134115358e-06,5.64515166327e-06
    2012-06-05 00:15:00,2.76599159913e-06,1.10481592295e-05
      ...,
    2012-06-08 23:45:00,1.56789791135e-06,8.55723150951e-06
