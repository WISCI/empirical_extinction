import os
import numpy as np
from astropy.table import QTable

from measure_extinction.stardata import StarData


if __name__ == "__main__":

    otab_lat = QTable(
        # fmt: off
        names=("name", "MIRI Band", "MIRI Flux", "IRAC1", "IRAC2", "IRAC3", "IRAC4", 
                "WISE1", "WISE2", "WISE3"),
        dtype=("S", "S", "S", "S", "S", "S", "S",
               "S", "S", "S"),
        # fmt:on
    )

    filename = "data/wisci_all.dat"

    f = open(filename, "r")
    file_lines = list(f)
    starnames = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)

    for cname in starnames:
        # get the photometry
        sdata = StarData(f"{cname.lower()}.dat", 
                        path="/home/kgordon/Python/extstar_data/MW/",
                        photonly=True)
        bands = sdata.data["BAND"].get_band_names()

        rdata_lat = []
        rdata_lat.append(cname.upper())

        for pband in ["MIRI_FND", "MIRI_F1000W", "MIRI_F1500W"]:
            if pband in bands:
                pmag = sdata.data["BAND"].get_band_mag(pband)
                rdata_lat.append(pband.split("_")[1])
                rdata_lat.append(fr"${pmag[0]:.2f} \pm {pmag[1]:.2f}$")

        for pband in ["IRAC1", "IRAC2", "IRAC3", "IRAC4", "WISE1", "WISE2", "WISE3"]:
            if pband in bands:
                pmag = sdata.data["BAND"].get_band_mag(pband)
                rdata_lat.append(fr"${pmag[0]:.2f} \pm {pmag[1]:.2f}$")
            else:
                rdata_lat.append("\\nodata")

        otab_lat.add_row(rdata_lat)

    otab_lat.write(
        f"tables/wisci_photometry.tex",
        format="aastex",
        col_align="lccccccccc",
        latexdict={
            "caption": r"Photometry \label{tab_phot}",
        },
        overwrite=True,
    )