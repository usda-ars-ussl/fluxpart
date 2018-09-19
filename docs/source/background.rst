.. |H2O| replace:: H\ :sub:`2`\O
.. |CO2| replace:: CO\ :sub:`2`

.. _fluxpart-background:

==========
Background
==========

(Adapted from [SAAS+18]_)

The `eddy covariance method`__ is routinely used to measure gas fluxes over agricultural fields and other landscapes.
For many applications, it is desirable to partition the measured fluxes into constitutive components:
the water vapor flux into transpiration and direct evaporation components,
and the carbon dioxide flux into photosynthesis and respiration components.
The flux variance similarity (FVS) partitioning method [SS08]_ is based on correlation analyses of eddy covariance data.
**Fluxpart** is a Python 3  module that implements the FVS method.
Estimates of constituent flux components are obtained from processing high-frequency eddy covariance data.

.. _ecwiki: https://en.wikipedia.org/wiki/Eddy_covariance

__ ecwiki_

The FVS flux partitioning method is described at length in the literature ([SS08]_, [SK10]_, [PRV14]_).
We present here a brief overview.

Monin-Obukhov similarity theory implies that high-frequency time series for scalars,
such as the water vapor (:math:`q`) and carbon dioxide (:math:`c`) concentrations,
will exhibit perfect correlation when measured at the same point within a homogeneous atmospheric layer. 
Actual measured correlations may deviate from this prediction for several reasons,
including the presence of multiple, distinct source/sinks for the scalar quantities within the layer.
Flux variance similarity and perfect correlation may exist for concentrations associated with a single source/sink,
but the superposition of fluxes from multiple source/sinks degrades the overall correlation.

In the case of :math:`q` and :math:`c`,
one source/sink arises from the exchange of :math:`q` and :math:`c` across leaf stomata during transpiration and photosynthesis,
and a second from non-stomatal direct evaporation and respiration.
If only transpiration and photosynthesis occur, similarity theory predicts :math:`\rho_{q,c}=-1`
(a negative correlation because transpiration acts as a :math:`q` source and photosynthesis as a :math:`c` sink).
Conversely, if only evaporation and respiration occur, theory predicts  :math:`\rho_{q,c}=1`
(evaporation and respiration being sources for :math:`q` and :math:`c`, respectively).

Scanlon and Sahu [SS08]_ note that one may think of evaporation and respiration as contaminating the transpiration and photosynthesis fluxes,
driving the  :math:`q`--:math:`c` correlation away from the expected :math:`\rho_{q,c}=-1`.
The premise of the FVS technique is that an analysis of the degree of that contamination can be used to infer the relative amounts of stomatal and non-stomatal fluxes present.

More specifically, Scanlon and Sahu [SS08]_ propose the following partitioning analysis.
The water vapor concentration, carbon dioxide concentration, and vertical wind velocity (:math:`w`) measured at a point on an eddy covariance tower can each be decomposed as

.. math::
    q = \langle q \rangle + q' \\
    c = \langle c \rangle + c' \\
    w = \langle w \rangle + w' \\

where angle brackets indicate the temporal mean over a short interval (e.g., 15 to 60 min) and the prime indicates the fluctuation from the mean.

According to the conventional eddy covariance method, the water vapor and |CO2| fluxes for the interval are, respectively,

.. math::
    F_q = \langle w' q' \rangle \\
    F_c = \langle w' c' \rangle \\

The gas concentrations and fluxes can be further decomposed into components regulated by stomatal and non-stomatal controls:

.. math::
    q' = q_e' + q_t' \\
    c' = c_r' + c_p'

.. math::
    F_q = F_{q_e} + F_{q_t} = \langle w'q_e' \rangle + \langle w'q_t' \rangle \\
    F_c = F_{c_r} + F_{c_p} = \langle w'c_r' \rangle + \langle w'c_p' \rangle

where
:math:`q_e'` and :math:`q_t'` are the water vapor concentration fluctuations associated with non-stomatal (evaporation) and stomatal (transpiration) controls, respectively; 
:math:`c_r'` and :math:`c_p'` are the |CO2| concentration fluctuations associated with non-stomatal (respiration) and stomatal (photosynthesis) controls, respectively;
and :math:`F_{q_e}`, :math:`F_{q_t}`, :math:`F_{c_r}`, and :math:`F_{c_p}` are the corresponding flux components.

These definitions along with some assumptions about the relative transfer efficiencies of various flux quantities allowed Scanlon and Sahu [SS08]_ to derive a system of equations
that can be used to calculate the constituent flux components.
We formulate this two-equation system as

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

Note that the system defined by Eqs. :eq:`eq1` and :eq:`eq2` has two branches due the presence of the plus-minus operator in the defintion of :math:`\langle w'c_r' \rangle / \langle w'c_p' \rangle`.

The system contains five parameters that are known directly from eddy covariance data and three unknown parameters.
The known parameters are:
:math:`F_q = \langle w'q' \rangle` and :math:`F_c = \langle w'c' \rangle`, the water vapor and |CO2| fluxes, respectively;
:math:`\sigma_{q}^2` and :math:`\sigma_{c}^2`, the variances of the water vapor and |CO2| concentrations, respectively;
and :math:`\rho_{q,c}`, the correlation coefficient for the water vapor and |CO2| concentrations.
The free parameters are:
:math:`\sigma_{c_p}^2`, the variance of the photosynthesis |CO2| concentration; 
:math:`\rho_{c_p,c_r}`, the correlation coefficient for the photosynthesis and respiration |CO2| concentrations; and 
:math:`W`, the leaf-level water use efficiency.
The latter is defined

.. math::
    W = \frac{\langle w'c_p' \rangle}{\langle w'q_t' \rangle}

If a value for :math:`W` is known from leaf-level measurements or can be otherwise estimated (see Appendix A of [SS08]_),
then Eqs. :eq:`eq1` and :eq:`eq2` can be solved for the remaining unknowns, :math:`\sigma_{c_p}^2` and :math:`\rho_{c_p,c_r}`.
A trial solution for Eqs. :eq:`eq1` and :eq:`eq2` is given by [SAAS+18]_ 

.. math::
  \sigma_{c_p}^2 = 
    \frac{(1 - \rho_{q,c}^2) \ (\sigma_{q} \sigma_{c} W)^2
            \left(\sigma_{q}^2 \langle w'c' \rangle^2 - 2 \rho_{q,c} \sigma_{q} \sigma_{c} \langle w'q' \rangle \langle w'c' \rangle + \sigma_{c}^2 \langle w'q' \rangle^2\right)}
         {\left[\sigma_{c}^2 \langle w'q' \rangle + \sigma_{q}^2 \langle w'c' \rangle W - \rho_{q,c} \sigma_{q} \sigma_{c} (\langle w'c' \rangle + \langle w'q' \rangle W) \right]^2}

.. math::
  \rho_{c_p,c_r}^2 = 
    \frac{(1 - \rho_{q,c}^2) \ \sigma_{q}^2 \sigma_{c}^2 \ \left(\langle w'c' \rangle - \langle w'q' \rangle W\right)^2}
         {\left(\sigma_{q}^2 \langle w'c' \rangle^2 - 2\rho_{q,c}\sigma_{c}\sigma_{q}\langle w'c' \rangle\langle w'q' \rangle + \sigma_{c}^2\langle w'q' \rangle^2\right)
          \left(\sigma_{c}^2 - 2\rho_{q,c}\sigma_{c}\sigma_{q}W + \sigma_{q}^2W^2\right)}

where :math:`\rho_{c_p,c_r}` is the negative square root of :math:`\rho_{c_p,c_r}^2`.
If the trial solution satisfies one of the branches of Eqs. :eq:`eq1` and :eq:`eq2`, the flux components are then

.. math::
    F_{c_p} = \left.  \langle w'c' \rangle \middle/ \left( \frac{\langle w'c_r'\rangle}{\langle w'c_p'\rangle} + 1 \right) \right.

.. math::
    F_{c_r} = F_c - F_{c_p}

.. math::
    F_{q_t} = F_{c_p} / W

.. math::
    F_{q_e} = F_q - F_{q_t}

The partitioning procedure does not always produce a result because the trial solution may not satisfy to Eqs. :eq:`eq1` and :eq:`eq2` or the solution may be non-physical.
By definition or by physical reasoning, it is required that :math:`-1 < \rho_{c_p,c_r} < 0`,  :math:`\sigma_{c_p}^2 > 0`, and :math:`W < 0`.
The FVS method also assumes/requires that the stomatal |CO2| flux is directed downward and the other fluxes are upward,

.. math::
   F_{c_p} < 0 \\
   F_{c_r},~F_{q_t},~F_{q_e} > 0

No solution for a given time interval may be the correct outcome.
For example, meteorological conditions may be incompatible with the theory or assumptions underlying the FVS method.
On the other hand, Scanlon and Sahu [SS08]_ found that, in some instances,
the root cause of failure may be the presence of large-scale eddies that affect flux variance similarity relationships but do not contribute significantly to fluxes.
**Fluxpart** retries failed analyses after filtering :math:`q`, :math:`c`, and :math:`w` time series data (using `PyWavelets`__, [LWGW+06]_).
Low-frequency components are progressively removed from the data and the partitioning procedure is applied at each stage, 
quitting when either a successful solution is found or the decomposition is exhausted.
**Fluxpart** additionally includes capabilities for applying basic QA/QC to high-frequency eddy covariance data,
for correcting high frequency data for external fluctuations associated with air temperature and vapor density ([WPL80]_, [DK07]_),
and for estimating water use efficiency by various models. 

.. _pywavelets: https://github.com/Py-Wavelets/pywt

__ pywavelets_

In sum, applying the partitioning procedure to an interval of eddy covariance data requires input values for
:math:`F_q = \langle w'q' \rangle`, :math:`F_c=\langle w'c' \rangle`, :math:`\sigma_{c}^2`, :math:`\sigma_{q}^2`, :math:`\rho_{q,c}^2`, and :math:`W`. 
Water use efficiency, :math:`W`, can be estimated by **Fluxpart** from atmospheric |CO2| and |H2O| data if direct measurements are not available.
Due to a possible need to remove low-frequency components from the data,
the partitioning algorithm is typically implemented using high-frequency time series data as input,
rather than just interval average values.
