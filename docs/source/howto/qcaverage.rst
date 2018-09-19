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
is made in :func:`fluxpart.partition.fvspart_interval`.

**Example**: Palatella et al. [PRV14]_ (Table 1) reported the following
quantities for an hour time interval:

========== ==========================
stddev(q)  0.411163 g/m^3
stddev(c)  5.182580 mg/m^3
<w'q'>     0.033140 g/m^2/s
<w'c'>     -0.472108 mg/m^2/s
corr(q, c) -0.881017
wue        -37.158598 mg CO2 / g H2O
========== ==========================

A flux partitioning solution for this data can found with:

.. code:: python

    from fluxpart.partition import fvspart_interval
    from fluxpart.containers import WQCData

    # Note converting input to SI units

    wue = -37.158598e-3
    interval_data = WQCData(
        var_q=0.411163e-3 ** 2,                                                 
        var_c=5.182580e-6 ** 2,                                                 
        wq=0.033140e-3,                                                         
        wc=-0.472108e-6,                                                        
        corr_qc=-0.881017, 
    )
    massfluxes, fvsp = fvspart_interval(interval_data, wue, wipe_if_invalid=True) 

    print(fvsp)
    print()
    print(massfluxes)

which produces:

.. code::

    -------------
    FVSP Solution
    -------------
      valid_partition: True
      fvsp_mssg: 
      Wavelet filtering level: nan
      var_q = 0.1691 (g/m^3)^2
      var_c = 26.86 (mg/m^3)^2
      corr_qc = -0.881
      cov_w_q = 0.03314 g/m^2/s
      cov_w_c = -0.4721 mg/m^2/s
      corr_cp_cr = -0.7585
      var_cp = 27.36 (mg/m^3)^2
      sig_cr = 0.06352 mg/m^3
      co2soln_id = 1
      valid_root = True
      root_mssg = 

    MassFluxes(
        Fq = 0.03314 g/m^2/s,
        Fqt = 0.01282 g/m^2/s,
        Fqe = 0.02032 g/m^2/s,
        Fc = -0.4721 mg/m^2/s,
        Fcp = -0.4765 mg/m^2/s,
        Fcr = 0.004389 mg/m^2/s)

The obtained results match closely those reported in [PRV14]_:

================== ========= ============
\                  [PRV14]_  **fluxpart**
================== ========= ============
var(cp) (mg/m^3)^2 27.358130 27.36
corr(cp,cr)        -0.757626 -0.7585
<w'qt'> (g/m^3/s)  0.012878  0.01282
<w'qe'> (g/m^3/s)  0.020262  0.02032
<w'cp'> (mg/m^3/s) -0.476489 -0.4765
<w'cr'> (mg/m^3/s) 0.004381  0.004389
================== ========= ============
