.. _qcaverage-howto:

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

    TODO: Add code block

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

    TODO: add codeblock

