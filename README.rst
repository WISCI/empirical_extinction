Code for WISCI empirical extinction curves
==========================================

Calculates UV-MIR extinction curves based on pair method using stellar
atmospheres for the unreddened comparison.
Extensively uses measure_extinction repository.

In Development!
---------------

Active development.
Data still in changing.
Use at your own risk.

Contributors
------------
Karl Gordon

License
-------

This code is licensed under a 3-clause BSD style license (see the
``LICENSE`` file).

Extinction Curves
-----------------

Extinction curves created by running the `fitstars` bash script.  This fits the
STIS_Opt and UBVRIJHK (whatever exists) using `utils/fit_model.py`.  After determining
the best fit stellar parameters, the E(lambda - 0.55 um) extinction curve is created 
using those parameters.

The extinction curves are normalized to A(0.55 um) by running using `utils/fit_mir_ext_powerlaw.py`.
This fits the E(lambda - 0.55 um) curves with a a broken power + 2 modified Drudes for the 10
and 20 micron silicate features.  The 20 micron silicate feature center, width, and asymmetry 
are fixed as the MIRI IFU data beyond 20 micron have systematic issues due to the low flux of 
the stars.  A(0.55 um) is one of the parameters of these fits.  With A(0.55 um), the
A(lambda)/A(0.55 um) curve is calculated.

Figures
-------

1. UV-MIR spectra of all stars: plotting/plot_spec_stack.py

Tables
------
