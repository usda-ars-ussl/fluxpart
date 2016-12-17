.. _saving-example:

.. _pickle: https://docs.python.org/3/library/pickle.html

Saving and retrieving results
-----------------------------
The results from a **Fluxpart** analysis can be saved to a file and later
retrieved using the pickle_ module from the Python standard library. If
``fpout`` contains :func:`~fluxpart.fluxpart.flux_partition()` output (or is
list of outputs, such as in this :ref:`example<tutorial-example>`), it can be
saved by::
 
    import pickle
    with open("fpout.pkl", 'wb') as f:
        pickle.dump(fpout, f)

The results can be retrieved in a subsequent Python session with::

    import pickle
    with open("fpout.pkl", 'rb') as f:
        fpout = pickle.load(f)
