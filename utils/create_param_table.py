import os
import numpy as np
from astropy.table import QTable

from measure_extinction.extdata import ExtData
from measure_extinction.extdata import conv55toAv, conv55toRv


def prettyname(name):
    return name.upper()


if __name__ == "__main__":

    for ctype in [""]:

        otab = QTable(
            # fmt: off
            names=("name", "E4455", "E4455_unc", "A55", "A55_unc", "R55", "R55_unc", 
                   "AV", "AV_unc", "RV", "RV_unc", 
                   "C1", "C1_unc", "C2", "C2_unc", "B3", "B3_unc", "C4", "C4_unc",
                   "x0", "x0_unc", "gamma", "gamma_unc"),
            dtype=("S", "f", "f", "f", "f", "f", "f", "f", "f", "f", 
                   "f", "f", "f", "f",
                   "f", "f", "f", "f", "f", "f", "f", "f", "f"),
            # fmt:on
        )

        otab_lat = QTable(
            # fmt: off
            names=("Name", 
                   r"$E(44-55)$", r"$A(55)$", r"$R(55)$", r"$A(V)$", r"$R(V)$"),
            dtype=("S", "S", "S", "S", "S", "S")
            # fmt:on
        )


        otab_lat2 = QTable(
            # fmt: off
            names=("Name", r"$C_1$", r"$C_2$", r"$B_3$", r"$C_4$", r"$x_o$", r"$\gamma$"),
            dtype=("S", "S", "S", "S", "S", "S", "S")
            # fmt:on
        )

        otab_lat3 = QTable(
            # fmt: off
            names=("Name", r"$\log(T_\mathrm{eff})$", r"$\log(g)$", r"$v_\mathrm{vturb}$", r"$vel$"),
            dtype=("S", "S", "S", "S", "S")
            # fmt:on
        )


        colnames = ["EBV", "AV", "RV", "AVd", "RVd"]
        fm90names = ["C1", "C2", "B3", "C4", "xo", "gamma"]
        stellnames = ["logTeff", "logg", "vturb", "velocity"]
        
        files = ["all"]
        tags = ["All"]
        for cfile, ctag in zip(files, tags):

            filename = "data/wisci_uv.dat"

            f = open(filename, "r")
            file_lines = list(f)
            starnames = []
            for line in file_lines:
                if (line.find("#") != 0) & (len(line) > 0):
                    name = line.rstrip()
                    starnames.append(name)
            # starnames = np.sort(starnames)

            for cname in starnames:
                ctype = ""
                cpath = "exts"
                cfile = f"{cpath}/{cname}_mefit_ext_a55.fits"
                edata = ExtData(filename=cfile)

                rdata = []
                rdata_lat = []
                rdata_lat2 = []
                rdata_lat3 = []
                rdata.append(cname)
                pcname = prettyname(cname)
                print(pcname)
                rdata_lat.append(pcname)
                rdata_lat2.append(pcname)
                rdata_lat3.append(pcname)

                for ccol in colnames:
                    if ccol == "AVd":
                        # AV, EBV in file actually A55, E4455
                        Av = conv55toAv(edata.columns["AV"], edata.columns["EBV"])
                        val = Av[0]
                        unc = Av[1]
                    elif ccol == "RVd":
                        # RV in file actually R55
                        Rv = conv55toRv(edata.columns["RV"])
                        val = Rv[0]
                        unc = Rv[1]
                    else:
                        val = edata.columns[ccol][0]
                        unc = edata.columns[ccol][1]
                    rdata.append(val)
                    rdata.append(unc)
                    if ccol == "NHI":
                        val /= 1e21
                        unc /= 1e21
                    rdata_lat.append(fr"${val:.2f} \pm {unc:.2f}$")

                cfile2 = f"{cpath}/{cname}_mefit_ext_a55_FM90.fits"
                if os.path.isfile(cfile2):
                    edata2 = ExtData(filename=cfile2)

                    fdata = edata2.fit_params["MCMC"]
                    ip50, = np.where(fdata["name"] == "p50")
                    iunc, = np.where(fdata["name"] == "unc")
                    for ccol in fm90names:
                        val = fdata[ccol][ip50[0]]
                        unc = fdata[ccol][iunc[0]]
                        rdata.append(val)
                        rdata.append(unc)
                        if ccol == "xo":
                            tstr = fr"${val:.3f} \pm {unc:.3f}$"
                        else:
                            tstr = fr"${val:.2f} \pm {unc:.2f}$"
                        rdata_lat2.append(tstr)
                    otab_lat2.add_row(rdata_lat2)
                else:
                    for ccol in fm90names:
                        rdata.append(0.0)
                        rdata.append(0.0)
                        rdata_lat2.append(" ")

                cfile3 = f"{cpath}/{cname}_mefit_ext.fits"
                edata3 = ExtData(filename=cfile3)

                fdata = edata3.fit_params["MCMC"]
                for ccol in stellnames:
                    idx, = np.where(fdata["name"] == ccol)
                    val = fdata[idx]["value"].data[0]
                    unc = fdata[idx]["unc"].data[0]
                    if ccol in ["logTeff", "logg"]:
                        tstr = fr"${val:.3f} \pm {unc:.3f}$"
                    else:
                        tstr = fr"${val:.2f} \pm {unc:.2f}$"
                    rdata_lat3.append(tstr)
                otab_lat3.add_row(rdata_lat3)

                otab.add_row(rdata)
                otab_lat.add_row(rdata_lat)

        basestr = "G25_wisci"
        otab.write(
            f"tables/{basestr}{ctype}_ensemble_params.dat", format="ascii.ipac", overwrite=True
        )

        otab_lat.write(
            f"tables/{basestr}{ctype}_ensemble_dust_params.tex",
            format="aastex",
            col_align="lccccc",
            latexdict={
                "caption": r"Extinction Parameters \label{tab_ext_col_param}",
            },
            overwrite=True,
        )

        otab_lat2.write(
            f"tables/{basestr}{ctype}_ensemble_fm90_params.tex",
            format="aastex",
            col_align="lcccccc",
            latexdict={
                "caption": r"FM90 Parameters \label{tab_ext_fm90_params}",
            },
            overwrite=True,
        )

        otab_lat3.write(
            f"tables/{basestr}{ctype}_ensemble_stell_params.tex",
            format="aastex",
            col_align="lcccc",
            latexdict={
                "caption": r"Stellar Parameters \label{tab_ext_stell_params}",
            },
            overwrite=True,
        )