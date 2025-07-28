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
STIS_Opt and JHK (whatever exists) using `utils/fit_model.py`.  After determining
the best fit stellar parameters, the E(lambda - 0.55 um) extinction curve is created 
using those parameters.

The classical A(55) values are measured by extrapolating JHK assuming the Rieke, 
Rieke, & Paul (1989) extinction shape.  This is not as accurate as what is described in the
next paragraph, but useful for UV work and as a comparison.  This is cone with the 
`utils/compute_av_JHK.py` code.

The extinction curves are normalized to A(0.55 um) by running using `utils/fit_mir_ext_powerlaw.py`.
This fits the E(lambda - 0.55 um) curves with a a broken power + 2 modified Drudes for the 10
and 20 micron silicate features.  The 20 micron silicate feature center, width, and asymmetry 
are fixed as the MIRI IFU data beyond 20 micron have systematic issues due to the low flux of 
the stars.  A(0.55 um) is one of the parameters of these fits.  With A(0.55 um), the
A(lambda)/A(0.55 um) curve is calculated.  The `fitfora55` bash script runs this for all the
stars.

The FM90 parameters are fit using the `fit_uv_ext_fm90.py` code.  Almost identical to what
was done in Gordon et al. (2024).  The bash script `fit_fm90_all` runs this for all the stars
that have STIS UV data.

Figures
-------

1. UV-MIR spectra of all stars: plotting/plot_spec_stack.py --models

2. Example fit: `meplot_model 2massj182302 --obspath ~/Python/extstar_data/MW/ --picmodname tlusty_z100_modinfo.p`

3. UV-MIR extinction curves for all stars: `plotting/plot_uvoptir_ext.py`

Tables
------

1. Photometry: utils/create_phot_table.py 

2. Fitting parameters: By hand

3. Stellar fit parameters: utils/create_param_table.py

4. Column fit parameters: utils/create_param_table.py

5. FM90 fit parameters: utils/create_param_table.py. 
