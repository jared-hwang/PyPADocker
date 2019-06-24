from optparse import OptionParser
# from matplotlib import pyplot as plt
import numpy as np
from scipy import constants
from random import random
import os.path
import sys

__doc__ = "Rewritten version of the original ReadIGUN from 2010. This is now compatible with the latest " \
          "version of IGUN (2016)."
__author__ = "Daniel Winklehner, Alberto Lemut"


def read_igun(filename, npart=5000):
    """
    Legacy function for backwards compatibility
    :param filename: 
    :param npart: 
    :return: 
    """
    ir = IgunReader()
    ir.read_trj(filename)

    e_rrp, e_xxp = ir.get_emittance()
    print("rrp: {} mm-mrad".format(e_rrp))
    print("xxp: {} mm-mrad".format(e_xxp))

    return ir.generate_dist(npart=npart)

ReadIGUN = read_igun


class IgunReader(object):
    """
    This class contains the necessary functions to read IGUN .TRJ files and return a distribution of
    randomly generated particles to match the IGUN trajectories.
    """

    def __init__(self):
        """
        Constructor
        """
        self.filename = None  # path and filename of trj file
        self.run_label = None  # Unique string at beginning of trj file
        self.data = None  # Structured numpy array containing the input data
        self.ns = 0  # Number of species
        self.legacy = False  # Flag for legacy file handling

    def read_trj(self, 
                 filename=None, 
                 resolution=0.25  # polygon units --> mm
                 ):
        """
        Function that reads in the values from TRJ file
        :param filename: 
        :param resolution: 
        :return: 
        """

        if filename is None:
            return None

        rest, ext = os.path.splitext(filename)  # extract basename (incl. path) and extension

        if ext not in [".trj", ".TRJ", ".Trj"]:
            return None

        self.filename = filename

        with open(filename) as infile:

            self.run_label = infile.readline()
            raw_data = infile.readlines()

        # Flag for legacy file handling
        if "BETA" in raw_data[-1]:
            self.legacy = False
        else:
            self.legacy = True

        raw_data.pop(-1)  # Delete the last row, it contains the column headers

        mydtype = [("ray", float),  # Ray number
                   ("group", float),  # Ray group
                   ("q", float),  # Ray charge (e)
                   ("m", float),  # Ray mass (amu)
                   ("rho", float),  # Ray R (in polygon units!)
                   ("zeta", float),  # Ray Z (in polygon units!)
                   ("energy", float),  # Ray kinetic energy per charge state (i.e. source voltage eV)
                   ("atandrdz", float),  # Ray R' (rad)
                   ("i", float),  # Ray current (A)
                   ("atandtdz", float),  # Ray Theta'
                   ("phi", float),  # Ray Theta (rad)
                   ("beta", float)]  # Ray beta (relativistic parameter)

        data = []

        if self.legacy:

            for line in raw_data:
                data_ = [float(item) for item in line.split()]
                # Old IGUN .TRJ files didn't have the BETA column...
                # Calculate gamma from energy and mass (Ekin = m0c^2 * (gamma - 1)).
                # Cave: Energy is given as source voltage. i.e. needs to be multiplied with charge state.
                gamma = data_[2] * data_[6] / data_[3] / 931500000.0 + 1.0
                # Calculate beta from gamma and append to data
                data_.append(np.sqrt(1.0 - gamma ** -2.0))
                data.append(tuple(data_))

            self.data = np.array(data, dtype=mydtype)

            # noinspection PyTypeChecker
            self.data["i"] *= -6.2832e-6  # Legacy currents were given as uA/6.2832 and pos. ions had negative currents

        else:

            for line in raw_data:
                data.append(tuple([float(item) for item in line.split()]))

            self.data = np.array(data, dtype=mydtype)

            # noinspection PyTypeChecker
            self.data["i"] *= -1.0  # Positive currents are given as negative in IGUN

        self.ns = len(np.unique(self.data["group"]))

        # noinspection PyTypeChecker
        self.data["zeta"] *= resolution  # Polygon units --> mm
        # noinspection PyTypeChecker
        self.data["rho"] *= resolution   # Polygon units --> mm

        return data

    def get_emittance(self):

        groups = np.array(np.unique(self.data["group"]))

        e_rrp = []
        e_xxp = []

        for species in groups:  # process each species

            data = self.data[np.where(self.data["group"] == species)]  # Select subset of data

            # species = int(species) - 1  # So we can use it as an index

            r = data["rho"]  # (mm)
            rp = data["atandrdz"] * 1000.0  # (rad --> mrad)
            currents = np.array(data["i"])  # (A)

            currentsum = sum(currents)

            e_rrp.append(np.sqrt(
                sum(currents * r ** 2.0) * sum(currents * rp ** 2.0) - sum(currents * r * rp) ** 2.0) / currentsum)

            e_xxp.append(np.sqrt(0.5 *
                                 sum(currents * r ** 2.0) * sum(currents * rp ** 2.0) 
                                 - sum(currents * r * rp) ** 2.0) / currentsum)

        return np.array(e_rrp), np.array(e_xxp)

    def generate_dist(self, npart=5000):
        """
        Uses the loaded data to generate a random particle distribution corresponding to the trajectory info.
        :param npart: Number of particles to generate per species.
        :return:
        """

        groups = np.array(np.unique(self.data["group"]))

        x = np.zeros((self.ns, npart), 'd')
        y = np.zeros((self.ns, npart), 'd')
        xp = np.zeros((self.ns, npart), 'd')
        yp = np.zeros((self.ns, npart), 'd')
        z = np.zeros((self.ns, npart), 'd')
        vx = np.zeros((self.ns, npart), 'd')
        vy = np.zeros((self.ns, npart), 'd')
        vz = np.zeros((self.ns, npart), 'd')
        currentsum = np.zeros(self.ns, 'd')

        pps = []
        mass = []
        charge = []

        for species in groups:  # process each species

            data = self.data[np.where(self.data["group"] == species)]  # Select subset of data

            species = int(species) - 1  # So we can use it as an index

            numpart = len(data)
            pps.append(numpart)
            mass.append(data["m"][0])
            charge.append(data["q"][0])

            currentsum[species] = sum(data["i"])

            cumulative = np.zeros(numpart + 1, 'd')

            for k in range(numpart):
                cumulative[k + 1] = cumulative[k] + data["i"][k] / currentsum[species]

            # indices = []

            for k in range(npart):

                probability = random()  # get random number
                jmin = 0
                jmid = int(numpart / 2)
                jmax = numpart

                for dummy in range(200):
                    if cumulative[jmin] <= probability <= cumulative[jmid]:
                        if jmin + 1 == jmid:
                            jmid = jmin
                            break
                        jmax = jmid
                        jmid = int((jmin + jmax) / 2)
                    elif cumulative[jmid] <= probability <= cumulative[jmax]:
                        if jmid + 1 == jmax:
                            break
                        jmin = jmid
                        jmid = int((jmin + jmax) / 2)
                    else:
                        print("{}: probability {} of out boundaries cumulative[{}] = {} - cumulative[{}] = {}\n".format(
                            os.path.split(sys.argv[0])[1], probability, jmin,
                            cumulative[jmin], jmax, cumulative[jmax]))

                jmid -= 1

                theta = 2.0 * np.pi * random()
                velocity = data["beta"][jmid] * constants.c
                x[species, k] = data["rho"][jmid] * np.cos(theta)  # (mm)
                y[species, k] = data["rho"][jmid] * np.sin(theta)  # (mm)
                z[species, k] = data["zeta"][0]  # (mm)
                xp[species, k] = (data["atandrdz"][jmid] * np.cos(theta) - data["atandtdz"][jmid] * np.sin(theta))
                yp[species, k] = (data["atandrdz"][jmid] * np.sin(theta) + data["atandtdz"][jmid] * np.cos(theta))
                vz[species, k] = velocity / np.sqrt(xp[species, k] ** 2 + yp[species, k] ** 2 + 1)  # (m/s)
                vx[species, k] = xp[species, k] * vz[species, k]  # (m/s)
                vy[species, k] = yp[species, k] * vz[species, k]  # (m/s)

        # Calculate some handy additional output values
        vzmean = vz.mean(axis=1)  # Calculate mean vz for each species (m/s)
        xmax = x.max(axis=1)
        ymax = y.max(axis=1)

        xenv = np.zeros(self.ns, 'd')
        yenv = np.zeros(self.ns, 'd')

        for k in range(self.ns):

            maxid = np.where(x[k, :] == xmax[k])
            xenv[k] = xp[k, maxid[0]]  # rad
            maxid = np.where(y[k, :] == ymax[k])
            yenv[k] = yp[k, maxid[0]]  # rad

        results = {"value": 0,
                   "ns": self.ns,
                   "np": np.ones(self.ns, 'd') * npart,
                   "pps": np.array(pps),
                   "M": np.array(mass),
                   "Q": np.array(charge),
                   "totalCurrent": currentsum * 1000000.0,  # Postprocessor expects current in uA
                   "x": x,
                   "y": y,
                   "z": z,
                   "xp": xp * 1000,  # mrad
                   "yp": yp * 1000,  # mrad
                   "vx": vx,
                   "vy": vy,
                   "vz": vz,
                   "vzmean": vzmean,
                   "xmax": xmax,
                   "ymax": ymax,
                   "xenv": xenv,
                   "yenv": yenv}

        return results


# This part is only executed if ReadIGUN.py is called on it's own (e.g. from command line)
if __name__ == '__main__':
    # --- Option parser for command-line options --- #
    parser = OptionParser()
    parser.add_option("-i", "--infile", dest="ipf", help="Specify input file (*.TRJ)", metavar="IFILE")
    parser.add_option("-o", "--outfile", dest="opf", help="Specify output file", metavar="OFILE")
    parser.add_option("-p", "--particles", dest="nparticles", type="int", help="Number of particles to be generated",
                      metavar="NP")

    (options, args) = parser.parse_args()

    # Set variables according to command-line options:
    if options.ipf is not None:
        ipf = options.ipf
    else:
        ipf = None
    if options.nparticles is not None:
        nparticles = options.nparticles
    else:
        nparticles = 5000
    if options.opf is not None:
        opf = options.opf
    else:
        opf = None

    if ipf is None or opf is None:
        print("Error: Either inputfile or outputfile not specified!")
        print("Usage: 'ReadIGUN.py -i <INPUT FILE> -o <OUTPUT FILE> [-p <# of particles to calculate>]'")
        raise SystemExit

    # Call the main script
    igun_reader = IgunReader()
    igun_reader.read_trj(filename=ipf)
    res = igun_reader.generate_dist(npart=nparticles)

    # --- write results to file --- #
    for j in range(res["ns"]):
        outpath = os.path.splitext(opf)[0] + "_species" + str(j + 1) + ".dat"  # each species gets it's own file
        print("Output file {} written\n".format(outpath))
        f = open(outpath, 'w')
        f.write("Original file: %s\n" % (os.path.split(sys.argv[0])[1]))
        f.write("M = %i amu\n" % (res["M"][j]))
        f.write("Q = %i e\n" % (res["Q"][j]))
        f.write("I = %f euA\n" % (res["totalCurrent"][j]))
        f.write(
            "x (mm)		y (mm)		z (mm)		xp (mrad)	yp (mrad)	vx (m/s)	vy (m/s)	vz (m/s)\n")
        for i in range(res["np"][j]):
            f.write("%e	%e	%e	%e	%e	%e	%e	%e\n" % (
                res["x"][j][i], res["y"][j][i], res["z"][j][i],
                res["xp"][j][i] * 1000, res["yp"][j][i] * 1000, res["vx"][j][i], res["vy"][j][i],
                res["vz"][j][i]))
        f.close()
