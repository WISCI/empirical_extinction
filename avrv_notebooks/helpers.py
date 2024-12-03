import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from dust_extinction.parameter_averages import G23
from measure_extinction.extdata import ExtData


def get_fit_params(params):
    """
    Removes the parameters that are fixed so the fitter does not fit them
    """
    fit_params = np.zeros(5)
    fit_params[0:2] = params[0:2]
    fit_params[2:4] = params[3:5]
    fit_params[4] = params[-1]

    return fit_params


def get_full_params(params):
    """
    Expands the fit parameters to the full list
    """
    full_params = np.zeros(len(params) + 2)
    full_params[0:2] = params[0:2]
    full_params[3:-1] = params[2:]
    full_params[-1] = 17.0

    return full_params

def get_fit_params_g23only(params):
    """
    Removes the parameters that are fixed so the fitter does not fit them
    """
    fit_params = np.zeros(5)
    fit_params[0:2] = params[0:2]
    fit_params[2:4] = params[3:5]
    fit_params[4] = params[-1]

    return fit_params


def get_full_params_g23only(params):
    """
    Expands the fit parameters to the full list
    """
    # dummy info needed to pass the lnprior check
    full_params = [4.2, 4.0, 0.0, 5.0, 3.1, 0.7, 3.23, 0.41, 4.59, 0.95, 20.5, 17.0,-6.0]
    full_params[0:2] = params[0:2]
    full_params[3:5] = params[2:4]
    full_params[-1]  = params[-1]

    return full_params


def plot_data_model(
    reddened_star,
    weights,
    modinfo,
    params_in,
    paramnames,
    starname,
    velocity,
    fullparamsfunc=get_full_params,
    fit_range="all",
    params_unc=None,
    prange=None,
    rjflux=False,
    savefig=False,
    figname='default',
):

    params = fullparamsfunc(params_in)

    # intrinsic sed
    modsed = modinfo.stellar_sed(params[0:3], velocity=velocity)

    # dust_extinguished sed
    ext_modsed = modinfo.dust_extinguished_sed(params[3:10], modsed, fit_range=fit_range)

    # hi_abs sed
    hi_ext_modsed = modinfo.hi_abs_sed(params[10:12], [velocity, 0.0], ext_modsed)

    # create a StarData object for the best fit SED
    # modsed_stardata = modinfo.SED_to_StarData(modsed)

    norm_model = np.average(hi_ext_modsed["BAND"])
    norm_data = np.average(reddened_star.data["BAND"].fluxes).value

    # plotting setup for easier to read plots
    fontsize = 18
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    # setup the plot
    fig, axes = plt.subplots(
        nrows=2, figsize=(13, 10), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )

    # plot the bands and all spectra for this star
    ax = axes[0]
    for cspec in modinfo.fluxes.keys():
        if cspec == "BAND":
            ptype = "o"
        else:
            ptype = "-"

        # ax.plot(reddened_star.data[cspec].waves,
        #        weights[cspec], 'k-')

        gvals = reddened_star.data[cspec].fluxes > 0.0
        modspec = hi_ext_modsed[cspec] * norm_data / norm_model
        if rjflux:
            plot_y = reddened_star.data[cspec].fluxes[gvals] * np.power(
                reddened_star.data[cspec].waves[gvals], 4.0
            )
            plot_ymod = (
                modsed[cspec][gvals]
                * (norm_data / norm_model)
                * np.power(modinfo.waves[cspec][gvals], 4.0)
            )
            plot_ymodext = (
                ext_modsed[cspec][gvals]
                * (norm_data / norm_model)
                * np.power(modinfo.waves[cspec][gvals], 4.0)
            )
            plot_ymodnotext = modspec[gvals] * np.power(
                modinfo.waves[cspec][gvals], 4.0
            )
        else:
            plot_y = reddened_star.data[cspec].fluxes[gvals]
            plot_ymod = modsed[cspec][gvals] * (norm_data / norm_model)
            plot_ymodext = ext_modsed[cspec][gvals] * (norm_data / norm_model)
            plot_ymodnoext = modspec[gvals]

        ax.plot(
            reddened_star.data[cspec].waves[gvals],
            plot_y,
            "k" + ptype,
            label="data",
            alpha=0.7,
        )

        # print(reddened_star.data[cspec].waves)
        # print(modinfo.waves[cspec])

        ax.plot(
            modinfo.waves[cspec][gvals],
            plot_ymod,
            "b" + ptype,
            label=cspec,
            alpha=0.5,
        )
        ax.plot(
            modinfo.waves[cspec][gvals],
            plot_ymodext,
            "r" + ptype,
            label=cspec,
            alpha=0.5,
        )
        ax.plot(
            modinfo.waves[cspec][gvals],
            plot_ymodnoext,
            "g" + ptype,
            label=cspec,
            alpha=0.5,
        )

        diff = 100.0 * (reddened_star.data[cspec].fluxes.value - modspec) / modspec
        axes[1].plot(reddened_star.data[cspec].waves, diff, "k" + ptype)

    # finish configuring the plot
    if prange is None:
        ax.set_ylim(4e5 * norm_data / norm_model, 1e7 * norm_data / norm_model)
    else:
        ax.set_ylim(prange)
    ax.set_yscale("log")
    ax.set_xscale("log")
    axes[1].set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
    ax.set_ylabel(r"$\lambda^4 F(\lambda)$ [RJ units]", fontsize=1.3 * fontsize)
    axes[1].set_ylabel("residuals [%]", fontsize=1.0 * fontsize)
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")
    axes[1].set_ylim(-10.0, 10.0)
    axes[1].plot([0.1, 2.5], [0.0, 0.0], "k:")

    if rjflux:
        textpos = [0.7, 0.5]
    else:
        textpos = [0.8, 0.95]
    for k, (pname, pval) in enumerate(zip(paramnames, params_in)):
        if params_unc is not None:
            ptxt = rf"{pname} = ${pval:.2f} \pm {params_unc[k]:.2f}$"
        else:
            ptxt = f"{pname} = {pval:.2f}"
        ax.text(
            textpos[0],
            textpos[1] - k * 0.04,
            ptxt,
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax.transAxes,
        )

    ax.text(0.1, 0.9, starname, transform=ax.transAxes)

    # ax.legend()

    # use the whitespace better
    fig.tight_layout()

    if savefig:
        fig.savefig(starname+'_'+figname+'.png',dpi=200)

#from measure_extinction.extdata import ExtData

def comp_ext(modinfo, fit_params, velocity, reddened_star, relband, starname):

    # intrinsic sed
    modsed = modinfo.stellar_sed(fit_params[0:3], velocity=velocity)

    # create a StarData object for the best fit SED
    modsed_stardata = modinfo.SED_to_StarData(modsed)

    # create an extincion curve and save it
    extdata = ExtData()
    extdata.calc_elx(reddened_star, modsed_stardata, rel_band=relband)

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

    if isinstance(relband, str):
        extdata.fit_band_ext()
    else:
        extdata.calc_AV_JHK()
    extdata.calc_RV()

    extdata.save(starname + "_ext.fits")

    return extdata


def plot_ext(extdata, starname, savefig=False, figname='extinction'):

    fig, ax = plt.subplots(figsize=(13, 10))
    fontsize = 18

    extdata.trans_elv_alav()

    extdata.plot(ax)  # , alax=True)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
    # ax.set_ylim(0.0, 10.0)
    ax.set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    # plot R(V) dependent relation extinction
    mod_x = np.arange(0.3, 8.7, 0.1) / u.micron
    rv = extdata.columns["RV"][0]
    G23_rv31 = G23(Rv=rv)
    ax.plot(1.0 / mod_x, G23_rv31(mod_x), "k-", label=f"G23 R(V) = {rv:.2f}")

    ax.legend()

    textpos = [0.7, 0.85]
    pvals = ["E(44-55)", "A(55)", "R(55)"]
    for k, (pname, pval) in enumerate(
        zip(pvals, extdata.columns.values())
    ):
        print(pval)
        ptxt = rf"{pname} = ${pval[0]:.2f} \pm {pval[1]:.2f}$"
        ax.text(
            textpos[0],
            textpos[1] - k * 0.04,
            ptxt,
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax.transAxes,
        )

    ax.text(0.7, 0.95, starname, transform=ax.transAxes)

    if savefig:
        fig.savefig(starname+'_'+figname+'.png',dpi=200)