.. _qcaverage-example:

Partitioning from interval average data
---------------------------------------
The flux partitioning procedure of [SS08]_ emphasizes the use of high-frequency
eddy covariance data because it may be necessary to filter the data to satisfy
variance similarity assumptions and obtain a physically consistent
partitioning. However, once the filter has been applied (or if a filter is not
required on a particular data set), the partitioning calculation is actually
made using interval average values for Fq=<w'q'>, Fc=<w'c'>, var(q), var(c),
corr(q,c), and WUE.

In **fluxpart**, the partitioning calculation using interval average values
is made in :func:`fluxpart.partition.partition_from_qc_averages`.

**Example**: Palatella et al. [PRV14]_ (Table 1) reported the following
quantities for an hour time interval

.. code-block:: none

    stddev(q) = 0.411163 g/m^3
    stddev(c) = 5.182580 mg/m^3
    <w'q'> = 0.033140 g/m^2/s
    <w'c'> = -0.472108 mg/m^2/s
    correlation(q, c) = -0.881017
    wue = -37.158598 mg CO2 / g H2O

A flux partitioning solution for this data can found with::

    from types import SimpleNamespace
    from fluxpart.partition import partition_from_qc_averages

    >>> wue = -37.158598e-3            # kg CO2 / kg H2O
    >>> qcdat = SimpleNamespace(
    ...     var_q=0.411163e-3 ** 2,    # (kg / m^3)^2
    ...     var_c=5.182580e-6 ** 2,    # (kg / m^3)^2
    ...     wq=0.033140e-3,            #  kg / m^2 / s
    ...     wc=-0.472108e-6,           #  kg / m^2 / s
    ...     corr_qc=-0.881017)
    >>> nsoln, fluxcom = partition_from_qc_averages(qcdat, wue)
    >>> print(nsoln)
    NumerSoln(
        corr_cp_cr = -0.7585,
        var_cp = 27.36 (mg/m^3)^2,
        sig_cr = 0.06352 mg/m^3,
        co2soln_id = 1,
        validroot = True,
        validmssg = ,
        init = (-0.8, 37.3),
        success = True,
        mssg = The solution converged.,
        nfev = 10)
    >>> print(fluxcom)
    FluxComponents(
        wq = 0.03314 g/m^2/s,
        wqt = 0.01282 g/m^2/s,
        wqe = 0.02032 g/m^2/s,
        wc = -0.4721 mg/m^2/s,
        wcp = -0.4765 mg/m^2/s,
        wcr = 0.004389 mg/m^2/s)

The obtained results match closely those reported in [PRV14]_:

================= ========= ============
\                 [PRV14]_  **fluxpart**
================= ========= ============
var_cp (mg/m^3)^2 27.358130 27.36
corr_cp_cr        -0.757626 -0.7585
wqt (g/m^3/s)     0.012878  0.01282
wqe (g/m^3/s)     0.020262  0.02032
wcp (mg/m^3/s)    -0.476489 -0.4765
wcr (mg/m^3/s)    0.004381  0.004389
================= ========= ============

The above data are for the "April 7 2010" example in Table 1 of [PRV14]_.
[PRV14]_ also give input data for "April 5 2010" which they find yields a
"non-physical" partitioning result::

    >>> wue = -24.558131e-3
    >>> qcdat = SimpleNamespace(
    ...    var_q=0.455994e-3 ** 2,
    ...    var_c=4.544450e-6 ** 2,
    ...    wq=0.062700e-3,
    ...    wc=-0.712862e-6,
    ...    corr_qc=-0.922292)
    >>> nsoln, fluxcom = partition_from_qc_averages(qcdat, wue)
    >>> print(nsoln)
    NumerSoln(
        corr_cp_cr = -0.7247,
        var_cp = 15.61 (mg/m^3)^2,
        sig_cr = -0.7749 mg/m^3,
        co2soln_id = 1,
        validroot = True,
        validmssg = ,
        init = (-0.8, 28.68),
        success = True,
        mssg = The solution converged.,
        nfev = 13)
    >>> print(fluxcom)
    FluxComponents(
        wq = 0.0627 g/m^2/s,
        wqt = 0.02542 g/m^2/s,
        wqe = 0.03728 g/m^2/s,
        wc = -0.7129 mg/m^2/s,
        wcp = -0.6242 mg/m^2/s,
        wcr = -0.0887 mg/m^2/s)

The results again match those in [PRV14]_. Although the numerical solution
converges, the obtained flux components are not physically allowable because
the respiration CO2 flux is wcr < 0 (directed toward the land surface).  Note
that :func:`~fluxpart.partition.partition_from_qc_averages` returns the
obtained fluxes without checking their validity. The function
:func:`fluxpart.partition.isvalid_partition` checks the validity of a
:class:`~fluxpart.containers.FluxComponents` object and returns an explanatory
message if invalid::

   >>> from fluxpart.partition import isvalid_partition
   >>> isvalid_partition(fluxcom)
   (False, 'Fcr <= 0; ')
