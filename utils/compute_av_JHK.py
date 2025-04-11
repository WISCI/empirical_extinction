import glob
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import PowerLaw1D, Const1D

from measure_extinction.extdata import ExtData

if __name__ == "__main__":

    fpath = "exts/"

    files = glob.glob(f"{fpath}*ext.fits")

    fontsize = 12

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    # fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=(10, 13))

    for fname in files:
        ifile = fname
        extdata = ExtData(ifile)

        # convert from E(l-V) to A(l)/A(V)
        # get E(44-55)
        bv_vals = np.interp(
            np.array([0.44, 0.55]),
            extdata.waves["STIS_Opt"].value,
            extdata.exts["STIS_Opt"],
        )
        ebv = -1.0 * np.diff(bv_vals)
        bv_vals_unc = np.interp(
            np.array([0.44, 0.55]),
            extdata.waves["STIS_Opt"].value,
            extdata.uncs["STIS_Opt"],
        )
        ebv_unc = np.sqrt(np.sum(np.square(bv_vals_unc)))
        extdata.columns["EBV"] = (ebv[0], ebv_unc)

        # get A(55)
        extdata.calc_AV_JHK()

        # get R(55)
        extdata.calc_RV()

        print(
            fname,
            extdata.columns["EBV"][0],
            extdata.columns["AV"][0],
            extdata.columns["RV"][0],
        )

        ofile = fname.replace("_ext.fits", "_ext_a55.fits")
        extdata.save(ofile)
