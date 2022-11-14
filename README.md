# Inferring Extinction Times in the Fossil Record using Test-Statistic Inversion

zID: z5209633

Supervisor: Professor David Warton

Precise estimates of extinction times in the fossil record are crucial to understanding how and why extinction events occurred. However, existing methods are unable to adequately address the sampling uncertainty of the fossil record as well as uncertainty introduced by radiocarbon dating and calibration. We propose two novel approaches to estimating extinction times, both based on test-statistic inversion: the Minimum Statistic Inversion (MINMI) estimator, and the Simulated Inversion - Robbins Monro Process (SI-RM). The MINMI estimator exploits assumptions of uniformity to directly estimate the probability of a test statistic exceeding the observed test-statistic value, from which point estimates and confidence intervals can be obtained by inversion. On the other hand, the SI-RM estimator uses stochastic approximation to obtain estimates via inversion. Under uniformity assumptions, we show that the MINMI and SI-RM estimators produce confidence intervals with better coverage, more appropriate widths, and faster computation times than existing methods.

Scripts and other miscellaneous R code is in the `src` directory. Notebooks for ad-hoc figures, analysis, and post-processing are in the `notebooks` directory, albeit in a haphazard organisation.
