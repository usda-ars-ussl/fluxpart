.. |H2O| replace:: H\ :sub:`2`\O
.. |CO2| replace:: CO\ :sub:`2`

.. _fluxpart-background:

==========
Background
==========

The `eddy covariance technique`__ is commonly used to measure gas fluxes
over agricultural fields.
Greater insight into the functioning of agroecosystems is possible if the
measured gas fluxes can be separated into their constitutive components:
the water vapor flux into transpiration and direct evaporation components,
and the carbon dioxide flux into photosynthesis and respiration components.
Scanlon and Sahu [SS08]_ devised a flux partitioning procedure based on flux
variance similarity relationships and correlation analyses of eddy covariance
data.
**Fluxpart** is a Python 3  module that implements the Scanlon and Sahu
[SS08]_ flux partitioning method.
Estimates of the constituent flux components are obtained from processing
high-frequency (e.g., 20 Hz) eddy covariance data.

.. _ecwiki: https://en.wikipedia.org/wiki/Eddy_covariance

__ ecwiki_

The Scanlon and Sahu [SS08]_ flux partitioning method has been
described at length in the literature ([SS08]_, [SK10]_, [PRV14]_).
We present here a brief overview.

Flux variance similarity theory implies that time series for scalars such as
water vapor (:math:`q`) and |CO2| (:math:`c`) concentrations
will exhibit perfect correlation when measured at the same point
within a homogeneous atmospheric layer.
Actual measured correlations may differ from this expectation for several
reasons, including the presence of separate, distinct source/sinks for
:math:`q` and :math:`c`:
one arising from the exchange of :math:`q` and :math:`c`
across leaf stomata during transpiration and photosynthesis,
and a second from sub-canopy (nonstomatal) direct evaporation and respiration.
If only transpiration and photosynthesis occur, similarity theory
predicts :math:`\rho_{q,c}=-1` (the correlation is negative because
transpiration acts as a :math:`q` source and photosynthesis as a :math:`c`
sink).
Conversely, if only evaporation and respiration occur, theory predicts
:math:`\rho_{q,c}=1`
(evaporation and respiration being sources for :math:`q` and :math:`c`,
respectively).
Flux variance similarity and perfect correlation may exist for
fluxes connected with one particular source/sink, but the superposition
of fluxes from both sinks degrades the :math:`q`-:math:`c` correlation.

Scanlon and Sahu [SS08]_ note that one may think of evaporation and
respiration as contaminating the transpiration and photosynthesis fluxes,
driving the :math:`q`-:math:`c` correlation away from the expected
:math:`\rho_{q,c}=-1`.
The premise of the Scanlon and Sahu [SS08]_ technique is that an analysis of
the degree of that contamination can be used to infer the relative
amounts of stomatal and nonstomatal fluxes present.

More specifically, Scanlon and Sahu [SS08]_ propose the following
partitioning analysis.
The water vapor concentration (:math:`q`), 
carbon dioxide concentration (:math:`c`), 
and vertical wind velocity (:math:`w`)
measured at a point on an eddy covariance tower can each be decomposed as

.. math::
    q = \langle q \rangle + q' \\
    c = \langle c \rangle + c' \\
    w = \langle w \rangle + w' \\

where angle brackets indicate the mean value obtained from time-averaging data
over a short interval (e.g., 15 to 60 min),
and a prime indicates the fluctuating component.
According to the conventional eddy covariance method,
the vapor and |CO2| fluxes for the interval are, respectively

.. math::
    F_q = \langle w' q' \rangle \\
    F_c = \langle w' c' \rangle \\

The gas concentrations and fluxes can be further decomposed into components
regulated by stomatal and non-stomatal controls:

.. math::
    q' = q_e' + q_t' \\
    c' = c_r' + c_p'

.. math::
    F_q = F_{q_e} + F_{q_t} = \langle w'q_e' \rangle + \langle w'q_t' \rangle \\
    F_c = F_{c_r} + F_{c_p} = \langle w'c_r' \rangle + \langle w'c_p' \rangle

where :math:`q_e'` and :math:`q_t'` 
are the vapor concentration fluctuations associated with
nonstomatal (evaporation) and
stomatal (transpiration) controls,
respectively; 
:math:`c_r'` and :math:`c_p'`
are the |CO2| concentration fluctuations associated with
nonstomatal (respiration) and
stomatal (photosynthesis) controls, respectively;
and
:math:`F_{q_e}`,
:math:`F_{q_t}`, 
:math:`F_{c_r}` and
:math:`F_{c_p}` 
are the corresponding flux components.

Scanlon and Sahu [SS08]_ assume that the transfer efficiencies
of the stomatal-controlled scalars are greater than those of the
nonstomatal scalars, which leads to the approximation

.. math::
    \rho_{q_t,q_e}
    \approx
    \frac{ \rho_{w,q_e} }{ \rho_{w,q_r} }
    =
    \frac{\langle w'q_e' \rangle}{\langle w'q_t' \rangle}
    \frac{\sigma_{q_t}}{\sigma_{q_e}}

.. math::
    \rho_{c_p,c_r}
    \approx
    \frac{ \rho_{w,c_r} }{ \rho_{w,c_p} }
    =
    \frac{\langle w'c_r' \rangle}{\langle w'c_p' \rangle}
    \frac{\sigma_{c_p}}{\sigma_{c_r}}

These definitions and approximations allow
Scanlon and Sahu [SS08]_ to derive a series of equations that can be used to
calculate the constituent flux components.
As noted by Palatella et al. [PRV14]_, the main computational task in
this procedure involves solving a system of two simultaneous
nonlinear equations.
We formulate the two-equation system as

.. math::
    :label: eq1

    W \frac{\langle w'q' \rangle}{\langle w'c' \rangle}
    \left(
    \frac{ \langle w'c_r' \rangle }{ \langle w'c_p' \rangle } + 1
    \right) 
    =
    \left(
    \frac{ \langle w'q_e' \rangle }{ \langle w'q_t' \rangle } + 1
    \right) 

.. math::
    :label: eq2

    W \rho_{q,c} \sigma_q \sigma_c \sigma_{c_p}^{-2} 
    = 1 +
    \frac{ \langle w'c_r' \rangle }{ \langle w'c_p' \rangle } + 
    \frac{ \langle w'q_e' \rangle }{ \langle w'q_t' \rangle } + 
    \rho_{c_p,c_r}^{-2}
    \frac{ \langle w'c_r' \rangle }{ \langle w'c_p' \rangle }
    \frac{ \langle w'q_e' \rangle }{ \langle w'q_t' \rangle }

where:

.. math::
    \frac{ \langle w'q_e' \rangle }{ \langle w'q_t' \rangle }
    = - \rho_{c_p,c_r}^2 + \rho_{c_p,c_r}^2
    \sqrt{1 - \rho_{c_p,c_r}^{-2}
    \left(1 - W^2 \sigma_q^2 / \sigma_{c_p}^2\right)}

.. math::
    \frac{ \langle w'c_r' \rangle }{ \langle w'c_p' \rangle }
    = - \rho_{c_p,c_r}^2 \pm \rho_{c_p,c_r}^2
    \sqrt{1 -  \rho_{c_p,c_r}^{-2}
    \left(1 - \sigma_c^2 / \sigma_{c_p}^2\right)}

This system contains three free parameters that are not known directly from the
eddy covariance data:
:math:`\sigma_{c_p}^2`,
the variance of the photosynthesis |CO2| concentration; 
:math:`\rho_{c_p,c_r}`,
the correlation coefficient for the photosynthesis and respiration |CO2|
concentrations;
and 
:math:`W`, 
the leaf-level water use efficiency, defined

.. math::
    W = \left. \langle w'c_p' \rangle \middle/ \langle w'q_t' \rangle \right.
      = \left. c_p' \middle/ q_t' \right.

A value for :math:`W` can be determined from leaf-level measurements made in
the field or estimated from concentration gradients. The system of equations
can then be solved numerically for 
:math:`\sigma_{c_p}^2`
and
:math:`\rho_{c_p,c_r}`.
The flux components are then given by,

.. math::
    F_{c_p} = 
    \left.  
    \langle w'c' \rangle 
    \middle/ 
    \left(
    \frac{ \langle w'c_r' \rangle }{ \langle w'c_p' \rangle } + 1
    \right) 
    \right.

.. math::
    F_{c_r} = F_c - F_{c_p}

.. math::
    F_{q_t} = F_{c_p} / W

.. math::
    F_{q_e} = F_q - F_{q_t}

The partitioning method is applicable only when the photosynthesis |CO2| flux
is directed downward and the other fluxes are upward,

.. math::
    F_{c_p} < 0 \\
    F_{c_r} > 0 \\
    F_{q_t} > 0 \\
    F_{q_e} > 0

Partitioning solutions that do not conform with this requirement are
considered "non-physical". If non-physical results are found, **Fluxpart**
uses wavelet filtering to progressively remove from the data low-frequency
components that do not contribute significantly to the fluxes but can
contaminant scalar correlations.
**Fluxpart** additionally includes capabilities for applying basic QA/QC to
high-frequency eddy covariance data,
for correcting high frequency data for external fluctuations associated with
air temperature and vapor density ([WPL80]_, [DK07]_),
and for estimating water use efficiency by various models. 
