from astropy.table import QTable
from measure_extinction.extdata import conv55toAv, conv55toRv

if __name__ == "__main__":

    tab = QTable.read(
        "wisci_a55r55.dat", format="ascii.commented_header", header_start=-1
    )

    otab = QTable(
        names=("name", "AV", "AVunc", "RV", "RVunc"),
        dtype=("str", "float64", "float64", "float64", "float64"),
    )
    for k in range(len(tab)):
        A55 = [tab["A55"][k], tab["A55unc"][k]]
        E4455 = [tab["E4455"][k], tab["E4455unc"][k]]
        R55 = [tab["R55"][k], tab["R55unc"][k]]
        Av = conv55toAv(A55, E4455)
        Rv = conv55toRv(R55)
        otab.add_row((tab["name"][k], Av[0], Av[1], Rv[0], Rv[1]))

    otab.write("wisci_avrv.dat", format="ascii.commented_header")