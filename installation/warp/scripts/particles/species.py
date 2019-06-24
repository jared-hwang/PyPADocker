"""
This module contains the definition of the Species class which is used to
define the species of the simulation particles.  It also defines the particle
types that are passed into Species. These include all of the atomic elements
(Hydrogen, etc.) and the following, Electron, Positron, Proton, Neutron,
Dihydrogen, Dinitrogen, Dioxygen, Carbon_Monoxide, Carbon_Dioxide, and Water
"""
from ..warp import *


def SpRandom(loc=0., scale=1., size=None):
    if scale > 0.:
        return random.normal(loc, scale, size)
    else:
        return loc


class Particle(object):
    def __init__(self, mass=None, charge=None, Symbol=None, name=None):
        self.Symbol = Symbol
        self.name = name
        if mass is not None:
            self.mass = mass
            self.M = self.mass
        if charge is not None:
            self.charge = charge
            self.Q = self.charge


class Atom(Particle):
    def __init__(self, Symbol, A, Z, Group, Period, name=None):
        Particle.__init__(self, mass=A*amu, Symbol=Symbol, name=name)
        self.A = A
        self.Z = Z
        self.Group = Group
        self.Period = Period


class Molecule(Particle):
    def __init__(self, Symbol, composition, name=None):
        self.composition = composition
        mass = 0.
        self.A = 0
        self.Z = 0
        for c in composition:
            mass += c.mass
            self.A += c.A
            self.Z += c.Z
        Particle.__init__(self, mass=mass, Symbol=Symbol, name=name)


periodic_table = {}
periodic_table['Hydrogen'] = {'A': 1.0079400000000001, 'Symbol': 'H', 'Z': 1, 'Group': 1, 'Period': 1}
periodic_table['Deuterium'] = {'A': 2.0141017800000001, 'Symbol': 'D', 'Z': 1, 'Group': 1, 'Period': 1}
periodic_table['Tritium'] = {'A': 3.0160492000000001, 'Symbol': 'T', 'Z': 1, 'Group': 1, 'Period': 1}
periodic_table['Helium'] = {'A': 4.0026020000000004, 'Symbol': 'He', 'Z': 2, 'Group': 18, 'Period': 1}
periodic_table['Lithium'] = {'A': 6.9409999999999998, 'Symbol': 'Li', 'Z': 3, 'Group': 1, 'Period': 2}
periodic_table['Lithium6'] = {'A': 6.01512279516, 'Symbol': 'Li6', 'Z': 3, 'Group': 1, 'Period': 2}
periodic_table['Lithium7'] = {'A': 7.016004558, 'Symbol': 'Li7', 'Z': 3, 'Group': 1, 'Period': 2}
periodic_table['Beryllium'] = {'A': 9.0121819999999992, 'Symbol': 'Be', 'Z': 4, 'Group': 2, 'Period': 2}
periodic_table['Boron'] = {'A': 10.811, 'Symbol': 'B', 'Z': 5, 'Group': 13, 'Period': 2}
periodic_table['Carbon'] = {'A': 12.0107, 'Symbol': 'C', 'Z': 6, 'Group': 14, 'Period': 2}
periodic_table['Nitrogen'] = {'A': 14.0067, 'Symbol': 'N', 'Z': 7, 'Group': 15, 'Period': 2}
periodic_table['Oxygen'] = {'A': 15.9994, 'Symbol': 'O', 'Z': 8, 'Group': 16, 'Period': 2}
periodic_table['Fluorine'] = {'A': 18.998403199999998, 'Symbol': 'F', 'Z': 9, 'Group': 17, 'Period': 2}
periodic_table['Neon'] = {'A': 20.1797, 'Symbol': 'Ne', 'Z': 10, 'Group': 18, 'Period': 2}
periodic_table['Sodium'] = {'A': 22.98977, 'Symbol': 'Na', 'Z': 11, 'Group': 1, 'Period': 3}
periodic_table['Magnesium'] = {'A': 24.305, 'Symbol': 'Mg', 'Z': 12, 'Group': 2, 'Period': 3}
periodic_table['Aluminium'] = {'A': 26.981538, 'Symbol': 'Al', 'Z': 13, 'Group': 13, 'Period': 3}
periodic_table['Silicon'] = {'A': 28.0855, 'Symbol': 'Si', 'Z': 14, 'Group': 14, 'Period': 3}
periodic_table['Phosphorus'] = {'A': 30.973761, 'Symbol': 'P', 'Z': 15, 'Group': 15, 'Period': 3}
periodic_table['Sulfur'] = {'A': 32.064999999999998, 'Symbol': 'S', 'Z': 16, 'Group': 16, 'Period': 3}
periodic_table['Chlorine'] = {'A': 35.453000000000003, 'Symbol': 'Cl', 'Z': 17, 'Group': 17, 'Period': 3}
periodic_table['Argon'] = {'A': 39.948, 'Symbol': 'Ar', 'Z': 18, 'Group': 18, 'Period': 3}
periodic_table['Potassium'] = {'A': 39.098300000000002, 'Symbol': 'K', 'Z': 19, 'Group': 1, 'Period': 4}
periodic_table['Calcium'] = {'A': 40.078000000000003, 'Symbol': 'Ca', 'Z': 20, 'Group': 2, 'Period': 4}
periodic_table['Scandium'] = {'A': 44.955910000000003, 'Symbol': 'Sc', 'Z': 21, 'Group': 3, 'Period': 4}
periodic_table['Titanium'] = {'A': 47.866999999999997, 'Symbol': 'Ti', 'Z': 22, 'Group': 4, 'Period': 4}
periodic_table['Vanadium'] = {'A': 50.941499999999998, 'Symbol': 'V', 'Z': 23, 'Group': 5, 'Period': 4}
periodic_table['Chromium'] = {'A': 51.996099999999998, 'Symbol': 'Cr', 'Z': 24, 'Group': 6, 'Period': 4}
periodic_table['Manganese'] = {'A': 54.938048999999999, 'Symbol': 'Mn', 'Z': 25, 'Group': 7, 'Period': 4}
periodic_table['Iron'] = {'A': 55.844999999999999, 'Symbol': 'Fe', 'Z': 26, 'Group': 8, 'Period': 4}
periodic_table['Cobalt'] = {'A': 58.933199999999999, 'Symbol': 'Co', 'Z': 27, 'Group': 9, 'Period': 4}
periodic_table['Nickel'] = {'A': 58.693399999999997, 'Symbol': 'Ni', 'Z': 28, 'Group': 10, 'Period': 4}
periodic_table['Copper'] = {'A': 63.545999999999999, 'Symbol': 'Cu', 'Z': 29, 'Group': 11, 'Period': 4}
periodic_table['Zinc'] = {'A': 65.409000000000006, 'Symbol': 'Zn', 'Z': 30, 'Group': 12, 'Period': 4}
periodic_table['Gallium'] = {'A': 69.722999999999999, 'Symbol': 'Ga', 'Z': 31, 'Group': 13, 'Period': 4}
periodic_table['Germanium'] = {'A': 72.640000000000001, 'Symbol': 'Ge', 'Z': 32, 'Group': 14, 'Period': 4}
periodic_table['Arsenic'] = {'A': 74.921599999999998, 'Symbol': 'As', 'Z': 33, 'Group': 15, 'Period': 4}
periodic_table['Selenium'] = {'A': 78.959999999999994, 'Symbol': 'Se', 'Z': 34, 'Group': 16, 'Period': 4}
periodic_table['Bromine'] = {'A': 79.903999999999996, 'Symbol': 'Br', 'Z': 35, 'Group': 17, 'Period': 4}
periodic_table['Krypton'] = {'A': 83.798000000000002, 'Symbol': 'Kr', 'Z': 36, 'Group': 18, 'Period': 4}
periodic_table['Rubidium'] = {'A': 85.467799999999997, 'Symbol': 'Rb', 'Z': 37, 'Group': 1, 'Period': 5}
periodic_table['Strontium'] = {'A': 87.620000000000005, 'Symbol': 'Sr', 'Z': 38, 'Group': 2, 'Period': 5}
periodic_table['Yttrium'] = {'A': 88.905850000000001, 'Symbol': 'Y', 'Z': 39, 'Group': 3, 'Period': 5}
periodic_table['Zirconium'] = {'A': 91.224000000000004, 'Symbol': 'Zr', 'Z': 40, 'Group': 4, 'Period': 5}
periodic_table['Niobium'] = {'A': 92.906379999999999, 'Symbol': 'Nb', 'Z': 41, 'Group': 5, 'Period': 5}
periodic_table['Molybdenum'] = {'A': 95.939999999999998, 'Symbol': 'Mo', 'Z': 42, 'Group': 6, 'Period': 5}
periodic_table['Technetium'] = {'A': 98.0, 'Symbol': 'Tc', 'Z': 43, 'Group': 7, 'Period': 5}
periodic_table['Ruthenium'] = {'A': 101.06999999999999, 'Symbol': 'Ru', 'Z': 44, 'Group': 8, 'Period': 5}
periodic_table['Rhodium'] = {'A': 102.9055, 'Symbol': 'Rh', 'Z': 45, 'Group': 9, 'Period': 5}
periodic_table['Palladium'] = {'A': 106.42, 'Symbol': 'Pd', 'Z': 46, 'Group': 10, 'Period': 5}
periodic_table['Silver'] = {'A': 107.8682, 'Symbol': 'Ag', 'Z': 47, 'Group': 11, 'Period': 5}
periodic_table['Cadmium'] = {'A': 112.411, 'Symbol': 'Cd', 'Z': 48, 'Group': 12, 'Period': 5}
periodic_table['Indium'] = {'A': 114.818, 'Symbol': 'In', 'Z': 49, 'Group': 13, 'Period': 5}
periodic_table['Tin'] = {'A': 118.70999999999999, 'Symbol': 'Sn', 'Z': 50, 'Group': 14, 'Period': 5}
periodic_table['Antimony'] = {'A': 121.76000000000001, 'Symbol': 'Sb', 'Z': 51, 'Group': 15, 'Period': 5}
periodic_table['Tellurium'] = {'A': 127.59999999999999, 'Symbol': 'Te', 'Z': 52, 'Group': 16, 'Period': 5}
periodic_table['Iodine'] = {'A': 126.90447, 'Symbol': 'I', 'Z': 53, 'Group': 17, 'Period': 5}
periodic_table['Xenon'] = {'A': 131.29300000000001, 'Symbol': 'Xe', 'Z': 54, 'Group': 18, 'Period': 5}
periodic_table['Caesium'] = {'A': 132.90545, 'Symbol': 'Cs', 'Z': 55, 'Group': 1, 'Period': 6}
periodic_table['Barium'] = {'A': 137.327, 'Symbol': 'Ba', 'Z': 56, 'Group': 2, 'Period': 6}
periodic_table['Lutetium'] = {'A': 174.96700000000001, 'Symbol': 'Lu', 'Z': 71, 'Group': 3, 'Period': 6}
periodic_table['Hafnium'] = {'A': 178.49000000000001, 'Symbol': 'Hf', 'Z': 72, 'Group': 4, 'Period': 6}
periodic_table['Tantalum'] = {'A': 180.9479, 'Symbol': 'Ta', 'Z': 73, 'Group': 5, 'Period': 6}
periodic_table['Tungsten'] = {'A': 183.84, 'Symbol': 'W', 'Z': 74, 'Group': 6, 'Period': 6}
periodic_table['Rhenium'] = {'A': 186.20699999999999, 'Symbol': 'Re', 'Z': 75, 'Group': 7, 'Period': 6}
periodic_table['Osmium'] = {'A': 190.22999999999999, 'Symbol': 'Os', 'Z': 76, 'Group': 8, 'Period': 6}
periodic_table['Iridium'] = {'A': 192.21700000000001, 'Symbol': 'Ir', 'Z': 77, 'Group': 9, 'Period': 6}
periodic_table['Platinum'] = {'A': 195.078, 'Symbol': 'Pt', 'Z': 78, 'Group': 10, 'Period': 6}
periodic_table['Gold'] = {'A': 196.96655000000001, 'Symbol': 'Au', 'Z': 79, 'Group': 11, 'Period': 6}
periodic_table['Mercury'] = {'A': 200.59, 'Symbol': 'Hg', 'Z': 80, 'Group': 12, 'Period': 6}
periodic_table['Thallium'] = {'A': 204.38329999999999, 'Symbol': 'Tl', 'Z': 81, 'Group': 13, 'Period': 6}
periodic_table['Lead'] = {'A': 207.19999999999999, 'Symbol': 'Pb', 'Z': 82, 'Group': 14, 'Period': 6}
periodic_table['Bismuth'] = {'A': 208.98038, 'Symbol': 'Bi', 'Z': 83, 'Group': 15, 'Period': 6}
periodic_table['Polonium'] = {'A': 210.0, 'Symbol': 'Po', 'Z': 84, 'Group': 16, 'Period': 6}
periodic_table['Astatine'] = {'A': 210.0, 'Symbol': 'At', 'Z': 85, 'Group': 17, 'Period': 6}
periodic_table['Radon'] = {'A': 220.0, 'Symbol': 'Rn', 'Z': 86, 'Group': 18, 'Period': 6}
periodic_table['Francium'] = {'A': 223.0, 'Symbol': 'Fr', 'Z': 87, 'Group': 1, 'Period': 7}
periodic_table['Radium'] = {'A': 226.0, 'Symbol': 'Ra', 'Z': 88, 'Group': 2, 'Period': 7}
periodic_table['Actinium'] = {'A': 227.0, 'Symbol': 'Ac', 'Z': 89, 'Group': 3, 'Period': 7}
periodic_table['Thorium'] = {'A': 232.0381, 'Symbol': 'Th', 'Z': 90, 'Group': 102, 'Period': 7}
periodic_table['Protactinium'] = {'A': 231.0359, 'Symbol': 'Pa', 'Z': 91, 'Group': 102, 'Period': 7}
periodic_table['Uranium'] = {'A': 238.02891, 'Symbol': 'U', 'Z': 92, 'Group': 102, 'Period': 7}
periodic_table['Neptunium'] = {'A': 237.0, 'Symbol': 'Np', 'Z': 93, 'Group': 102, 'Period': 7}
periodic_table['Plutonium'] = {'A': 244.0, 'Symbol': 'Pu', 'Z': 94, 'Group': 102, 'Period': 7}
periodic_table['Americium'] = {'A': 243.0, 'Symbol': 'Am', 'Z': 95, 'Group': 102, 'Period': 7}
periodic_table['Curium'] = {'A': 247.0, 'Symbol': 'Cm', 'Z': 96, 'Group': 102, 'Period': 7}
periodic_table['Berkelium'] = {'A': 247.0, 'Symbol': 'Bk', 'Z': 97, 'Group': 102, 'Period': 7}
periodic_table['Californium'] = {'A': 251.0, 'Symbol': 'Cf', 'Z': 98, 'Group': 102, 'Period': 7}
periodic_table['Einsteinium'] = {'A': 252.0, 'Symbol': 'Es', 'Z': 99, 'Group': 102, 'Period': 7}
periodic_table['Fermium'] = {'A': 257.0, 'Symbol': 'Fm', 'Z': 100, 'Group': 102, 'Period': 7}
periodic_table['Mendelevium'] = {'A': 258.0, 'Symbol': 'Md', 'Z': 101, 'Group': 102, 'Period': 7}
periodic_table['Nobelium'] = {'A': 259.0, 'Symbol': 'No', 'Z': 102, 'Group': 102, 'Period': 7}
periodic_table['Lawrencium'] = {'A': 262.0, 'Symbol': 'Lr', 'Z': 103, 'Group': 102, 'Period': 7}
periodic_table['Rutherfordium'] = {'A': 261.0, 'Symbol': 'Rf', 'Z': 104, 'Group': 4, 'Period': 7}
periodic_table['Lawrencium'] = {'A': 262.0, 'Symbol': 'Lr', 'Z': 103, 'Group': 3, 'Period': 7}
periodic_table['Dubnium'] = {'A': 262.0, 'Symbol': 'Db', 'Z': 105, 'Group': 5, 'Period': 7}
periodic_table['Bohrium'] = {'A': 264.0, 'Symbol': 'Bh', 'Z': 107, 'Group': 7, 'Period': 7}
periodic_table['Seaborgium'] = {'A': 266.0, 'Symbol': 'Sg', 'Z': 106, 'Group': 6, 'Period': 7}
periodic_table['Meitnerium'] = {'A': 268.0, 'Symbol': 'Mt', 'Z': 109, 'Group': 9, 'Period': 7}
periodic_table['Darmstadtium'] = {'A': 271.0, 'Symbol': 'Ds', 'Z': 110, 'Group': 10, 'Period': 7}
periodic_table['Roentgenium'] = {'A': 272.0, 'Symbol': 'Rg', 'Z': 111, 'Group': 11, 'Period': 7}
periodic_table['Hassium'] = {'A': 277.0, 'Symbol': 'Hs', 'Z': 108, 'Group': 8, 'Period': 7}
periodic_table['Ununtrium'] = {'A': 284.0, 'Symbol': 'Uut', 'Z': 113, 'Group': 13, 'Period': 7}
periodic_table['Ununbium'] = {'A': 285.0, 'Symbol': 'Uub', 'Z': 112, 'Group': 12, 'Period': 7}
periodic_table['Ununpentium'] = {'A': 288.0, 'Symbol': 'Uup', 'Z': 115, 'Group': 15, 'Period': 7}
periodic_table['Ununquadium'] = {'A': 289.0, 'Symbol': 'Uuq', 'Z': 114, 'Group': 14, 'Period': 7}
periodic_table['Ununhexium'] = {'A': 292.0, 'Symbol': 'Uuh', 'Z': 116, 'Group': 16, 'Period': 7}

elements = []
for k in periodic_table:
    S = periodic_table[k]['Symbol']
    A = periodic_table[k]['A']
    Z = periodic_table[k]['Z']
    G = periodic_table[k]['Group']
    P = periodic_table[k]['Period']
    exec(k + " = Atom(S, A, Z, G, P, name='%s')"%k)
    exec('elements.append(' + k + ')')
#  exec(k + " = periodic_table['" + k + "']")
del k

Electron = Particle(charge=-echarge, mass=emass, Symbol='e-', name='Electron')
Positron = Particle(charge=+echarge, mass=emass, Symbol='e+', name='Positron')
Muon = Particle(charge=-echarge, mass=1.883531475e-28, Symbol='mu-', name='Muon')
Antimuon = Particle(charge=+echarge, mass=1.883531475e-28, Symbol='mu+', name='Antimuon')

Proton = Particle(mass=1.6726231e-27, charge=echarge, Symbol='P', name='Proton')
AntiProton = Particle(mass=1.6726231e-27, charge=-echarge, Symbol='P', name='Antiproton')
Neutron = Particle(mass=1.6749286e-27, charge=0, Symbol='N', name='Neutron')

Dihydrogen = Molecule(composition=[Hydrogen, Hydrogen], Symbol='H2', name='Dihydrogen')
Dideuterium = Molecule(composition=[Deuterium, Deuterium], Symbol='D2', name='Dideuterium')
Dinitrogen = Molecule(composition=[Nitrogen, Nitrogen], Symbol='N2', name='Dinitrogen')
Dioxygen = Molecule(composition=[Oxygen, Oxygen], Symbol='O2', name='Dioxygen')
Carbon_Monoxide = Molecule(composition=[Carbon, Oxygen], Symbol='CO', name='Carbon_Monoxide')
Carbon_Dioxide = Molecule(composition=[Carbon, Oxygen, Oxygen], Symbol='CO2', name='Carbon_Dioxide')
Water = Molecule(composition=[Hydrogen, Hydrogen, Oxygen], Symbol='H2O', name='Water')

Photon = Particle(charge=0., mass=0., Symbol='gnu', name='Photon')

for elemt in elements:
    elemt.ionization_levels = zeros(elemt.Z)
    Ionization_Levels = elemt.ionization_levels
    Z = elemt.Z
    if Z == 1:  # --- Hydrogen
        Ionization_Levels[0] = 13.6
    if Z == 2:  # --- Helium
        Ionization_Levels[0] = 24.58
        Ionization_Levels[1] = 54.4
    if Z == 3:  # --- Lithium
        Ionization_Levels[0] = 5.392
        Ionization_Levels[1] = 75.638
        Ionization_Levels[2] = 122.451
    if Z == 5:  # --- Boron
        Ionization_Levels[0] = 8.297
        Ionization_Levels[1] = 25.154
        Ionization_Levels[2] = 37.929
        Ionization_Levels[3] = 259.367
        Ionization_Levels[4] = 340.216
    if Z == 6:  # --- Carbon
        Ionization_Levels[0] = 11.261
        Ionization_Levels[1] = 24.383
        Ionization_Levels[2] = 47.887
        Ionization_Levels[3] = 64.492
        Ionization_Levels[4] = 392.081
        Ionization_Levels[5] = 489.979
    if Z == 7:  # --- N
        Ionization_Levels[0] = 14.534
        Ionization_Levels[1] = 29.601
        Ionization_Levels[2] = 47.448
        Ionization_Levels[3] = 77.472
        Ionization_Levels[4] = 97.888
        Ionization_Levels[5] = 552.057
        Ionization_Levels[6] = 667.029
    if Z == 10:  # --- Ne
        Ionization_Levels[0] = 21.564
        Ionization_Levels[1] = 40.962
        Ionization_Levels[2] = 63.45
        Ionization_Levels[3] = 97.11
        Ionization_Levels[4] = 126.21
        Ionization_Levels[5] = 157.93
        Ionization_Levels[6] = 207.27
        Ionization_Levels[7] = 239.09
        Ionization_Levels[8] = 1195.79
        Ionization_Levels[9] = 1362.164
    if Z == 13:  # --- Aluminium
        Ionization_Levels[0] = 5.968
        Ionization_Levels[1] = 18.796
        Ionization_Levels[2] = 28.399
        Ionization_Levels[3] = 119.78
        Ionization_Levels[4] = 153.561
        Ionization_Levels[5] = 190.156
        Ionization_Levels[6] = 241.34
        Ionization_Levels[7] = 284.163
        Ionization_Levels[8] = 329.564
        Ionization_Levels[9] = 398.057
        Ionization_Levels[10] = 441.232
        Ionization_Levels[11] = 2082.379
        Ionization_Levels[12] = 2300.16
    if Z == 14:  # --- Silicon:
        Ionization_Levels[0] = 7.264E0
        Ionization_Levels[1] = 1.695E1
        Ionization_Levels[2] = 3.427E1
        Ionization_Levels[3] = 4.665E1
        Ionization_Levels[4] = 1.598E2
        Ionization_Levels[5] = 2.105E2
        Ionization_Levels[6] = 2.613E2
        Ionization_Levels[7] = 3.120E2
        Ionization_Levels[8] = 3.640E2
        Ionization_Levels[9] = 4.151E2
        Ionization_Levels[10] = 5.037E2
        Ionization_Levels[11] = 5.522E2
        Ionization_Levels[12] = 2.324E3
        Ionization_Levels[13] = 2.569E3
    if Z == 18:  # --- Argon
        Ionization_Levels[0] = 15.759
        Ionization_Levels[1] = 27.6
        Ionization_Levels[2] = 40.7
        Ionization_Levels[3] = 59.81
        Ionization_Levels[4] = 75.02
        Ionization_Levels[5] = 91.007
        Ionization_Levels[6] = 124.319
        Ionization_Levels[7] = 143.456
        Ionization_Levels[8] = 422.44
        Ionization_Levels[9] = 478.68
        Ionization_Levels[10] = 538.95
        Ionization_Levels[11] = 618.24
        Ionization_Levels[12] = 686.09
        Ionization_Levels[13] = 755.73
        Ionization_Levels[14] = 854.75
        Ionization_Levels[15] = 918
        Ionization_Levels[16] = 4120.778
        Ionization_Levels[17] = 4426.114
    if Z == 29:  # --- Copper:
        Ionization_Levels[0] = 7.70
        Ionization_Levels[1] = 20.234
        Ionization_Levels[2] = 36.739
        Ionization_Levels[3] = 57.213
        Ionization_Levels[4] = 79.577
        Ionization_Levels[5] = 102.313
        Ionization_Levels[6] = 138.48
        Ionization_Levels[7] = 165.354
        Ionization_Levels[8] = 198.425
        Ionization_Levels[9] = 231.492
        Ionization_Levels[10] = 264.566
        Ionization_Levels[11] = 367.913
        Ionization_Levels[12] = 399.951
        Ionization_Levels[13] = 434.055
        Ionization_Levels[14] = 482.627
        Ionization_Levels[15] = 518.798
        Ionization_Levels[16] = 554.970
        Ionization_Levels[17] = 631.446
        Ionization_Levels[18] = 668.672
        Ionization_Levels[19] = 1691.78
        Ionization_Levels[20] = 1799.261
        Ionization_Levels[21] = 1916.0
        Ionization_Levels[22] = 2060.0
        Ionization_Levels[23] = 2182.0
        Ionization_Levels[24] = 2308.0
        Ionization_Levels[25] = 2478.0
        Ionization_Levels[26] = 2587.5
        Ionization_Levels[27] = 11062.38
        Ionization_Levels[28] = 11567.617
    if Z == 36:  # --- Krypton
        Ionization_Levels[0] = 1.3999e1
        Ionization_Levels[1] = 2.4359e1
        Ionization_Levels[2] = 3.695e1
        Ionization_Levels[3] = 5.25e1
        Ionization_Levels[4] = 6.47e1
        Ionization_Levels[5] = 7.85e1
        Ionization_Levels[6] = 1.11e2
        Ionization_Levels[7] = 1.258e2
        Ionization_Levels[8] = 2.3085e2
        Ionization_Levels[9] = 2.682e2
        Ionization_Levels[10] = 3.08e2
        Ionization_Levels[11] = 3.50e2
        Ionization_Levels[12] = 3.91e2
        Ionization_Levels[13] = 4.47e2
        Ionization_Levels[14] = 4.92e2
        Ionization_Levels[15] = 5.41e2
        Ionization_Levels[16] = 5.92e2
        Ionization_Levels[17] = 6.41e2
        Ionization_Levels[18] = 7.86e2
        Ionization_Levels[19] = 8.33e2
        Ionization_Levels[20] = 8.84e2
        Ionization_Levels[21] = 9.37e2
        Ionization_Levels[22] = 9.98e2
        Ionization_Levels[23] = 1.051e3
        Ionization_Levels[24] = 1.151e3
        Ionization_Levels[25] = 1.2053e3
        Ionization_Levels[26] = 2.928e3
        Ionization_Levels[27] = 3.07e3
        Ionization_Levels[28] = 3.227e3
        Ionization_Levels[29] = 3.381e3
        Ionization_Levels[30] = 0.0
        Ionization_Levels[31] = 0.0
        Ionization_Levels[32] = 0.0
        Ionization_Levels[33] = 0.0
        Ionization_Levels[34] = 0.0
        Ionization_Levels[35] = 0.0
    if Z == 54:  # --- Xenon
        # The data below is taken from NIST (http://physics.nist.gov/)
        Ionization_Levels[0] = 12.1298431
        Ionization_Levels[1] = 20.975
        Ionization_Levels[2] = 31.05
        Ionization_Levels[3] = 42.20
        Ionization_Levels[4] = 54.1
        Ionization_Levels[5] = 66.703
        Ionization_Levels[6] = 91.6
        Ionization_Levels[7] = 105.9778
        Ionization_Levels[8] = 179.84
        Ionization_Levels[9] = 202.0
        Ionization_Levels[10] = 229.02
        Ionization_Levels[11] = 255.0
        Ionization_Levels[12] = 280.9
        Ionization_Levels[13] = 314.1
        Ionization_Levels[14] = 342.9
        Ionization_Levels[15] = 374.1
        Ionization_Levels[16] = 403.9
        Ionization_Levels[17] = 433.9
        Ionization_Levels[18] = 549.0
        Ionization_Levels[19] = 582
        Ionization_Levels[20] = 616
        Ionization_Levels[21] = 650
        Ionization_Levels[22] = 700
        Ionization_Levels[23] = 736
        Ionization_Levels[24] = 818
        Ionization_Levels[25] = 857.0
        Ionization_Levels[26] = 1493
        Ionization_Levels[27] = 1571
        Ionization_Levels[28] = 1653
        Ionization_Levels[29] = 1742
        Ionization_Levels[30] = 1826
        Ionization_Levels[31] = 1919
        Ionization_Levels[32] = 2023
        Ionization_Levels[33] = 2113
        Ionization_Levels[34] = 2209
        Ionization_Levels[35] = 2300
        Ionization_Levels[36] = 2556
        Ionization_Levels[37] = 2637
        Ionization_Levels[38] = 2726
        Ionization_Levels[39] = 2811
        Ionization_Levels[40] = 2975
        Ionization_Levels[41] = 3068
        Ionization_Levels[42] = 3243
        Ionization_Levels[43] = 3333.8
        Ionization_Levels[44] = 7660
        Ionization_Levels[45] = 7889
        Ionization_Levels[46] = 8144
        Ionization_Levels[47] = 8382
        Ionization_Levels[48] = 8971
        Ionization_Levels[49] = 9243
        Ionization_Levels[50] = 9581
        Ionization_Levels[51] = 9810
        Ionization_Levels[52] = 40272
        Ionization_Levels[53] = 41300
    if Z == 55:  # --- Cesium
        Ionization_Levels[0] = 3.894


listofallspecies = []
class Species(object):
    """
  Creates a new species of particles. All arguments are optional.
    - js=None: global species number (usually this should not be set)
    - type=Electron: one of either an elementary particle, an element, or one of
                     the predefined molecules
    - charge=echarge: charge in Coulombs, will be obtained automatically from
                      type or charge_state
    - mass=emass: mass in kg, will be obtained automatically from type
    - charge_state=0: charge_state of the ion or molecule
    - weight=None: simulation weight, will be obtained automatically from
                   other input
    - lvariableweights=False: When true, each particle is given its own weight
                              (top.wpid is set)
    - name='': species name
    - nautodt=1: number of species to setup for automatic subcycling.
                 The slowest will have dt = 2**(nautodt-1)*top.dt.
    - efetch=1: Method to use for fetching the self fields. See documentation
                of top.efetch for more info.
    - fselfb=None: z velocity to use when applying relativistic corrections.
                 No corrections are done if it is zero.
    - limplicit=false: Flag to turn on the implicit particle advance for this
                       species.
    - color='fg', marker='\1', msize=1.0: Default values used when making particle
                                          plots of the species.
    """
    def __init__(self, js=None, pgroup=None, pgroups=None,
                 type=None, charge=echarge, mass=emass, charge_state=0,
                 weight=None, lvariableweights=False, name='', nautodt=1,
                 efetch=None, fselfb=None, limplicit=None,
                 color='fg', marker='\1', msize=1.0):
        assert type is None or isinstance(type, Particle), 'type must be one of either an elementary particle, an element, or one of the predefined molecules'
        listofallspecies.append(self)
        # --- Note that some default arguments are None in case the user had
        # --- set the values in pgroup already, in which case they should not
        # --- be overwritten here unless the inputs are explicitly set.
        self.jslist = []
        # --- Setup the pid for particles weights. Do it early in the routine so that
        # --- the pid array gets allocated appropriately.
        if lvariableweights:
            top.wpid = top.nextpid()
        # --- If pgroup is not passed in, top.pgroup is used. But, note that
        # --- a reference to top.pgroup isn't actually saved. The pgroup
        # --- property would always returns top.pgroup in that case.
        if pgroup is not None:
            self._pgroup = pgroup
            if not pgroup.allocated('sm'):
                pgroup.gchange()
        self.pgroups = pgroups
        self.type = type
        self.add_group(charge=charge, mass=mass, charge_state=charge_state, weight=weight, js=js)
        self.charge = self.pgroup.sq[self.jslist[0]]
        self.mass = self.pgroup.sm[self.jslist[0]]
        if type is None or type.__class__ is not Particle:
            self.charge_state = charge_state
        self.name = name
        self.nautodt = nautodt
        if self.nautodt > 1 and top.chdtspid == 0:
            top.chdtspid = nextpid()
        if efetch is not None:
            top.efetch[self.jslist[-1]] = efetch
        if fselfb is not None:
            self.pgroup.fselfb[self.jslist[-1]] = fselfb
        if limplicit is not None:
            self.pgroup.limplicit[self.jslist[-1]] = limplicit
        for i in range(nautodt - 1):
            self.add_group(charge=self.charge, mass=self.mass, charge_state=self.charge_state, weight=weight)
            self.pgroup.ndts[self.jslist[-1]] = 2*self.pgroup.ndts[self.jslist[-2]]
            if efetch is not None:
                top.efetch[self.jslist[-1]] = efetch
            if fselfb is not None:
                self.pgroup.fselfb[self.jslist[-1]] = fselfb
            if limplicit is not None:
                self.pgroup.limplicit[self.jslist[-1]] = limplicit
            # --- zero out sp_fract for the extra species added with larger ndts
            top.sp_fract[self.jslist[-1]] = 0.

        # --- Save the default values for plotting
        self.color = color
        self.marker = marker
        self.msize = msize

    def __setstate__(self, dict):
        self.__dict__.update(dict)
        from ..particles import species
        # --- Check if the type is already declared in this module.
        # --- If it is, replace it with the one defined here since they
        # --- are supposed to be singletons.
        try:
            self.type = species.__dict__[self.type.name]
        except (KeyError, AttributeError):
            pass
        # --- Add the species to the list if it is not already there.
        if self not in listofallspecies:
            listofallspecies.append(self)

    def add_group(self, charge, mass, charge_state, weight, js=None):
        if js is None:
            # --- If there are no species defined or if ns is 1 (the default) and
            # --- the first species has already been setup, then a new species needs
            # --- to be added.
            if self.pgroup.ns == 0 or self.pgroup.sm[0] != 0.:
                if self.pgroup.name() == 'ParticleGroup':
                    addspecies(pgroup=self.pgroup)
                else:
                    self.pgroup.addspecies()
            js = self.pgroup.ns - 1
            # --- Setup the starting index for this species in the particle
            # --- arrays.
            if js > 0:
                self.pgroup.ins[js] = self.pgroup.ins[js-1] + self.pgroup.nps[js-1]
            self.pgroup.sid[js] = js
        self.jslist.append(js)
        type = self.type
        # --- set the charge
        try:
            # --- Try type.charge first, which is the charge of the fundemental
            # --- particle, or the user set charge.
            self.pgroup.sq[js] = type.charge
        except:
            # --- If that doesn't work...
            if type is not None and type.__class__ is not Particle:
                self.pgroup.sq[js] = echarge*charge_state
            else:
                self.pgroup.sq[js] = charge
        top.zion_s[js] = nint(self.pgroup.sq[js]/echarge)
        # set mass
        try:
            self.pgroup.sm[js] = type.mass
        except:
            try:
                self.pgroup.sm[js] = type.A*amu
            except:
                self.pgroup.sm[js] = mass
        top.aion_s[js] = self.pgroup.sm[js]/amu
        if weight is not None:
            self.pgroup.sw[js] = weight
        # set atomic number, if any
        try:
            top.nion_s[js] = type.Z
        except:
            pass

    def get_density(self, xmin=None, xmax=None, nx=None, ymin=None, ymax=None, ny=None, zmin=None, zmax=None,
                    nz=None, lost=0, charge=0, dens=None, l_minmax_grid=true, l_dividebyvolume=1, l4symtry=None, l2symtry=None):
        if l_minmax_grid:
            if xmin is None: xmin = w3d.xmmin
            if xmax is None: xmax = w3d.xmmax
            if ymin is None: ymin = w3d.ymmin
            if ymax is None: ymax = w3d.ymmax
            if zmin is None: zmin = w3d.zmmin + top.zgrid
            if zmax is None: zmax = w3d.zmmax + top.zgrid
            if l4symtry is None: l4symtry = w3d.l4symtry
            if l2symtry is None: l2symtry = w3d.l2symtry
        else:
            if xmin is None: xmin = min(self.getx())
            if xmax is None: xmax = max(self.getx())
            if ymin is None: ymin = min(self.gety())
            if ymax is None: ymax = max(self.gety())
            if zmin is None: zmin = min(self.getz())
            if zmax is None: zmax = max(self.getz())
            if l4symtry is None: l4symtry = false
            if l2symtry is None: l2symtry = false
        if dens is None:
            if nx is None: nx = w3d.nx
            if ny is None: ny = w3d.ny
            if nz is None: nz = w3d.nz
            if w3d.solvergeom is w3d.XYgeom:
                density = fzeros([nx+1, ny+1], 'd')
                densityc = fzeros([nx+1, ny+1], 'd')
            else:
                if w3d.solvergeom in [w3d.Zgeom]:
                    density = fzeros([nz+1], 'd')
                    densityc = fzeros([nz+1], 'd')
                elif w3d.solvergeom in [w3d.XZgeom, w3d.RZgeom]:
                    density = fzeros([nx+1, nz+1], 'd')
                    densityc = fzeros([nx+1, nz+1], 'd')
                else:
                    density = fzeros([nx+1, ny+1, nz+1], 'd')
                    densityc = fzeros([nx+1, ny+1, nz+1], 'd')
        else:
            if w3d.solvergeom is w3d.XYgeom:
                nx = shape(dens)[0] - 1
                ny = shape(dens)[1] - 1
            else:
                if w3d.solvergeom in [w3d.Zgeom]:
                    nz = shape(dens)[0] - 1
                if w3d.solvergeom in [w3d.XZgeom, w3d.RZgeom]:
                    nx = shape(dens)[0] - 1
                    nz = shape(dens)[1] - 1
                else:
                    nx = shape(dens)[0] - 1
                    ny = shape(dens)[1] - 1
                    nz = shape(dens)[2] - 1
            density = dens
            densityc = 0.*dens
        if self.pgroups is not None:
            pgroups = self.flatten(self.pgroups)
        else:
            pgroups = [self.pgroup]
        np = 0
        for pg in pgroups:
            for js in self.jslist:
                np += pg.nps[js]
        np = globalsum(np)
        if np == 0:
            if dens is None:
                return density
            else:
                return
        for pg in pgroups:
            for js in self.jslist:
                x = getx(js=js, lost=lost, gather=0, pgroup=pg)
                y = gety(js=js, lost=lost, gather=0, pgroup=pg)
                z = getz(js=js, lost=lost, gather=0, pgroup=pg)
                if w3d.solvergeom == w3d.RZgeom:
                    x = sqrt(x*x + y*y)
                np = shape(x)[0]
                if np > 0:
                    if top.wpid == 0:
                        w = pg.sw[js]*ones(np, 'd')
                    else:
                        w = pg.sw[js]*getpid(js=js, id=top.wpid-1, gather=0, pgroup=pg)
                    if charge:
                        w *= pg.sq[js]
                    if w3d.solvergeom is w3d.Zgeom:
                        deposgrid1d(1, np, z, w, nz, density, densityc, zmin, zmax)
                    elif w3d.solvergeom is w3d.XYgeom:
                        deposgrid2d(1, np, x, y, w, nx, ny, density, densityc, xmin, xmax, ymin, ymax)
                    elif w3d.solvergeom in [w3d.XZgeom, w3d.RZgeom]:
                        deposgrid2d(1, np, x, z, w, nx, nz, density, densityc, xmin, xmax, zmin, zmax)
                    else:
                        deposgrid3d(1, np, x, y, z, w, nx, ny, nz, density, densityc, xmin, xmax, ymin, ymax, zmin, zmax)
        if w3d.solvergeom is w3d.Zgeom:
            if l_dividebyvolume:
                density *= nz/(zmax - zmin)
        elif w3d.solvergeom is w3d.XYgeom:
            if l_dividebyvolume:
                density *= nx*ny/((xmax - xmin)*(ymax - ymin))
                if l4symtry:
                    density[0,:] *= 2
                    density[:,0] *= 2
                if l2symtry:
                    density[:,0] *= 2
                if w3d.boundxy is periodic:
                    density[0,:] += density[-1,:]
                    density[-1,:] = density[0,:]
                    density[:,0] += density[:,-1]
                    density[:,-1] = density[:,0]
        else:
            if w3d.solvergeom in [w3d.XZgeom, w3d.RZgeom]:
                if l_dividebyvolume:
                    density *= nx*nz/((xmax - xmin)*(zmax - zmin))
                    if l4symtry:
                        density[0,:] *= 2
                    if w3d.boundxy is periodic:
                        density[0,:] += density[-1,:]
                        density[-1,:] = density[0,:]
                    if w3d.bound0 is periodic:
                        density[:,0] += density[:,-1]
                        density[:,-1] = density[:,0]
                    if w3d.solvergeom == w3d.RZgeom:
                        dr = (xmax - xmin)/nx
                        r = xmin + arange(nx + 1)*dr
                        if xmin == 0.:
                            density[0,:] /= pi*dr/2
                        else:
                            density[0,:] /= 2.*pi*r[0]
                        for j in range(1, nx + 1):
                            density[j,:] /= 2.*pi*r[j]
            else:
                if l_dividebyvolume:
                    density *= nx*ny*nz/((xmax - xmin)*(ymax - ymin)*(zmax - zmin))
                    if l4symtry:
                        density[0,:,:] *= 2
                        density[:,0,:] *= 2
                    if l2symtry:
                        density[:,0,:] *= 2
                    if w3d.boundxy is periodic:
                        density[0,:,:] += density[-1,:,:]
                        density[-1,:,:] = density[0,:,:]
                        density[:,0,:] += density[:,-1,:]
                        density[:,-1,:] = density[:,0,:]
                    if w3d.bound0 is periodic:
                        density[:,:,0] += density[:,:,-1]
                        density[:,:,-1] = density[:,:,0]
        density[...] = parallelsum(density)
        if dens is None:
            return density

    def fetche(self):
        """Fetches the self E field, putting it into the self.pgroup.ex etc arrays."""
        fetche3d(self.pgroup, self.ins, self.nps, self.js+1)

    def iterpgroups(self):
        if self.pgroups is None:
            yield self.pgroup
        else:
            for pgroup in self.flatten(self.pgroups):
                yield pgroup

    def xpush3d(self, dt=None):
        if dt is None:
            dt = top.dt
        for pgroup in self.iterpgroups():
            for js in self.jslist:
                i1 = pgroup.ins[js] - 1
                i2 = i1 + pgroup.nps[js]
                xpush3d(pgroup.nps[js], pgroup.xp[i1:i2], pgroup.yp[i1:i2], pgroup.zp[i1:i2],
                        pgroup.uxp[i1:i2], pgroup.uyp[i1:i2], pgroup.uzp[i1:i2], pgroup.gaminv[i1:i2], dt)

    def epush3d(self, dt=None):
        if dt is None:
            dt = top.dt
        for pgroup in self.iterpgroups():
            for js in self.jslist:
                i1 = pgroup.ins[js] - 1
                i2 = i1 + pgroup.nps[js]
                epush3d(pgroup.nps[js], pgroup.uxp[i1:i2], pgroup.uyp[i1:i2], pgroup.uzp[i1:i2],
                        pgroup.ex[i1:i2], pgroup.ey[i1:i2], pgroup.ez[i1:i2],
                        pgroup.sq[js], pgroup.sm[js], dt)
                gammaadv(pgroup.nps[js], pgroup.gaminv[i1:i2], pgroup.uxp[i1:i2], pgroup.uyp[i1:i2], pgroup.uzp[i1:i2],
                         top.gamadv, top.lrelativ)

    def bpush3d(self, dt=None):
        if dt is None:
            dt = top.dt
        for pgroup in self.iterpgroups():
            for js in self.jslist:
                i1 = pgroup.ins[js] - 1
                i2 = i1 + pgroup.nps[js]
                bpush3d(pgroup.nps[js], pgroup.uxp[i1:i2], pgroup.uyp[i1:i2], pgroup.uzp[i1:i2], pgroup.gaminv[i1:i2],
                        pgroup.bx[i1:i2], pgroup.by[i1:i2], pgroup.bz[i1:i2], pgroup.sq[js], pgroup.sm[js], dt, top.ibpush)

    def getselffields(self):
        for pgroup in self.iterpgroups():
            for js in self.jslist:
                i1 = pgroup.ins[js]
                np = pgroup.nps[js]
                fetche3d(pgroup, i1, np, js+1)

    def getappliedfields(self, dt=None, dtl=None, dtr=None):
        if dt is None:
            dt = top.dt

        for pgroup in self.iterpgroups():
            for js in self.jslist:
                i1 = pgroup.ins[js] - 1
                i2 = i1 + pgroup.nps[js]

                np = pgroup.nps[js]
                m = pgroup.sm[js]
                q = pgroup.sq[js]
                x = pgroup.xp[i1:i2]
                y = pgroup.yp[i1:i2]
                z = pgroup.zp[i1:i2]
                ux = pgroup.uxp[i1:i2]
                uy = pgroup.uyp[i1:i2]
                uz = pgroup.uzp[i1:i2]
                gi = pgroup.gaminv[i1:i2]
                ex = pgroup.ex[i1:i2]
                ey = pgroup.ey[i1:i2]
                ez = pgroup.ez[i1:i2]
                bx = pgroup.bx[i1:i2]
                by = pgroup.by[i1:i2]
                bz = pgroup.bz[i1:i2]

                # --- Handle bends
                bendres = ones(np)
                bendradi = ones(np)
                getbend(np, np, z, uz, gi, bendres, bendradi, dtl, dtr, false)
                # --- Correct Ez_self for warped mesh effect
                # --- !!! This is fragile since it depends on this being called immediately after
                # --- !!! fetching the self feilds.
                bendez3d(np, x, z, ez, bendres, bendradi, top.bends, w3d.bnezflag, top.linbend)

                # --- Add in ears and uniform focusing E field pieces
                othere3d(np, x, y, z,
                         top.zbeam, top.zimax, top.zimin, top.straight, top.ifeears, top.eears,
                         top.eearsofz, top.dzzi, top.nzzarr, top.zzmin,
                         top.dedr, top.dexdx, top.deydy, top.dbdr, top.dbxdy, top.dbydx,
                         ex, ey, ez, bx, by, bz)

                # --- Add fields from the lattice elements
                exteb3d(np, x, y, z, uz, gi, dtl, dtr,
                        bx, by, bz, ex, ey, ez, m, q,
                        bendres, bendradi, dt)

    def loadrho(self, lzero=true, lfinalize_rho=True):
        for pgroup in self.iterpgroups():
            loadrho3d(pgroup, -1, -1, -1, lzero, lfinalize_rho)
        if lfinalize_rho:
            finalizerho()

    def particleboundaries(self, lcallcontrollers=True):
        for pgroup in self.iterpgroups():
            particleboundaries3d(pgroup, -1, lcallcontrollers)

    def flatten(self, lis):
        """Given a list, possibly nested to any level, return it flattened."""
        new_lis = []
        for item in lis:
            if isinstance(item, list):
                new_lis.extend(self.flatten(item))
            else:
                new_lis.append(item)
        return new_lis

    def getappliedfieldsold(self):
        dtl = -0.5*top.dt
        dtr = +0.5*top.dt
        bendres = ones(self.nps, 'd')
        bendradi = ones(self.nps, 'd')
    # --- getbend not available in Python yet
    # getbend(self.nps, self.nps,
    #         self.zp, self.uzp, self.gaminv,
    #         bendres, bendradi, dtl, dtr, false)

        othere3d(self.nps, self.xp, self.yp, self.zp,
                 top.zbeam, top.zimax, top.zimin, top.straight, top.ifeears, top.eears,
                 top.eearsofz, top.dzzi, top.nzzarr, top.zzmin,
                 top.dedr, top.dexdx, top.deydy, top.dbdr, top.dbxdy, top.dbydx,
                 self.ex, self.ey, self.ez,
                 self.bx, self.by, self.bz)

        exteb3d(self.nps, self.xp, self.yp, self.zp,
                self.uzp, self.gaminv, dtl, dtr,
                self.bx, self.by, self.bz,
                self.ex, self.ey, self.ez,
                self.sm, self.sq, bendres, bendradi, top.dt)

    def addparticles(self, x=0., y=0., z=0., vx=0., vy=0., vz=0., gi=1., js=None, w=None, **kw):
        if js is None:
            js = self.jslist[0]

        # --- Allow added pid names to be passed in through the keyword list.
        pidpairs = kw.setdefault('pidpairs', [])
        for k, v in kw.items():
            if k in self._addedpids:
                pidpairs.append([getattr(self, k + 'pid'), v])
                kw.pop(k)
            try:
                n = self._builtinpids[k]
                pidpairs.append([getattr(n[0], n[1]), v])
                kw.pop(k)
            except KeyError:
                pass

        return addparticles(x, y, z, vx, vy, vz, gi=gi, js=js, w=w, **kw)

    addpart = addparticles

    def add_uniform_box(self, np, xmin, xmax, ymin, ymax, zmin, zmax, vthx=0.,
                        vthy=0., vthz=0., vxmean=0., vymean=0., vzmean=0., js=None,
                        lmomentum=0, spacing='random', nx=None, ny=None, nz=None,
                        lallindomain=0, solvergeom=None, **kw):
        """
     - np: number of particles to load. When using uniform spacing, the actual
           number may be lower.
     - xmin, xmax, ymin, ymax, zmin, zmax: extent of the box to fill with particles.
                                           It is OK is one or more dimensions have
                                           a length of zero.
     - vthx=0., vthy=0., vthz=0.: Thermal velocity
     - vxmean=0., vymean=0., vzmean=0.: Drift velocity
     - lmomentum=false: Set to false when velocities are input as velocities, true
                        when input as massless momentum (as WARP stores them).
                        Only used when top.lrelativ is true.
     - spacing='random': either 'random' or 'uniform' particle spacing.
     - nx=None, ny=None, nz=None: With uniform spacing, number of particles along
                                  each axis. If not given, the number along each
                                  axis will be the same, calculated to give a
                                  total number that is less then but close to np.
                                  If all three are given, then np is ignored.
                                  For zero length dimensions, the number is set
                                  to one. They will be rounded to the nearest
                                  integer.
     - lallindomain=0: If true, the code only loads particles within the domain.
                       This only matters when parallel.
     - solvergeom=w3d.solvergeom: Specifies the geometry to use
        """

        if np == 0:
            if top.debug:
                print("add_uniform_box: Warning: no particles loaded since the input np is zero")
            return

        if lallindomain:
            # --- Crop the min and max to be within the local domain
            xdmin, xdmax = top.xpminlocal, top.xpmaxlocal
            ydmin, ydmax = top.ypminlocal, top.ypmaxlocal
            zdmin, zdmax = top.zpminlocal + top.zgrid, top.zpmaxlocal + top.zgrid
        else:
            # --- Crop the min and max to be within the global domain
            xdmin, xdmax = w3d.xmmin, w3d.xmmax
            ydmin, ydmax = w3d.ymmin, w3d.ymmax
            zdmin, zdmax = w3d.zmmin + top.zgrid, w3d.zmmax + top.zgrid

        if solvergeom is None:
            solvergeom = w3d.solvergeom
        if solvergeom == w3d.XZgeom:
            ydmin, ydmax = -largepos, +largepos

        # --- The min and max in the second argument guarantee that if the min and
        # --- max are outside of the grid, then it will end up with minp = maxp.
        xminp = max(xmin, min(xmax, xdmin))
        xmaxp = min(xmax, max(xmin, xdmax))
        yminp = max(ymin, min(ymax, ydmin))
        ymaxp = min(ymax, max(ymin, ydmax))
        zminp = max(zmin, min(zmax, zdmin))
        zmaxp = min(zmax, max(zmin, zdmax))

        if spacing == 'random':
            # --- Add a random number to the number of particles so that on
            # --- average, the correct number of particles will be generated.
            np = int(np + random.random())

            # --- Adjust the number of particles to load to based on the
            # --- width of the cropped min and max and the original
            def setfac(min, max, dmin, dmax, minp, maxp):
                if min == max:
                    if min > dmax or max < dmin:
                        return None
                    return 1.
                else:
                    return (maxp - minp)/(max - min)

            xfac = setfac(xmin, xmax, xdmin, xdmax, xminp, xmaxp)
            if xfac is None:
                return
            yfac = setfac(ymin, ymax, ydmin, ydmax, yminp, ymaxp)
            if yfac is None:
                return
            zfac = setfac(zmin, zmax, zdmin, zdmax, zminp, zmaxp)
            if zfac is None:
                return

            # --- Is nint correct?
            np = nint(xfac*yfac*zfac*np)

            if np == 0:
                if top.debug:
                    print("add_uniform_box: Warning: no particles loaded for random " +
                          "spacing since all are outside of the z range of the domain " +
                          "of processor %d"%me)
                return

            x = random.random(np)
            y = random.random(np)
            z = random.random(np)

            x = xminp + (xmaxp - xminp)*x
            y = yminp + (ymaxp - yminp)*y
            z = zminp + (zmaxp - zminp)*z

        elif spacing == 'uniform':
            assert nx is None or nx > 0, "if nx is given, it must be greater than zero"
            assert ny is None or ny > 0, "if ny is given, it must be greater than zero"
            assert nz is None or nz > 0, "if nz is given, it must be greater than zero"
            if nx is None and xmax <= xmin:
                nx = 1
            if ny is None and ymax <= ymin:
                ny = 1
            if nz is None and zmax <= zmin:
                nz = 1
            dims = 0
            nptemp = np
            if nx is None:
                dims += 1
            else:
                nptemp /= nx
            if ny is None:
                dims += 1
            else:
                nptemp /= ny
            if nz is None:
                dims += 1
            else:
                nptemp /= nz
            if dims > 0:
                nptemp = nint(nptemp**(1./dims))
                if nx is None:
                    nx = nptemp
                if ny is None:
                    ny = nptemp
                if nz is None:
                    nz = nptemp

            # --- Find the range of particle locations within the cropped
            # --- min and max.
            def setnp(bmin, bmax, dmin, dmax, minp, maxp, n):
                if bmin == bmax:
                    if bmin > dmax or bmax < dmin:
                        return None, None
                    return nint(n), 0
                else:
                    d = (bmax - bmin)/n
                    iminp = nint((minp - bmin)/d)
                    imaxp = nint((maxp - bmin)/d)
                    return max(0, imaxp - iminp), iminp

            nxp, ixminp = setnp(xmin, xmax, xdmin, xdmax, xminp, xmaxp, nx)
            if nxp is None:
                return
            nyp, iyminp = setnp(ymin, ymax, ydmin, ydmax, yminp, ymaxp, ny)
            if nyp is None:
                return
            nzp, izminp = setnp(zmin, zmax, zdmin, zdmax, zminp, zmaxp, nz)
            if nzp is None:
                return

            np = nxp*nyp*nzp

            if np == 0:
                if top.debug:
                    print("add_uniform_box: Warning: no particles loaded for uniform " +
                          "spacing since all are outside of the z range of the domain " +
                          "of processor %d"%me)
                return

            x, y, z = getmesh3d((ixminp + 0.5)/nx, 1./nx, nxp - 1,
                                (iyminp + 0.5)/ny, 1./ny, nyp - 1,
                                (izminp + 0.5)/nz, 1./nz, nzp - 1)

            # --- Perform a transpose so that the data is ordered with increasing z.
            # --- The copy is needed since transposed arrays cannot be reshaped.
            x = transpose(x).copy()
            y = transpose(y).copy()
            z = transpose(z).copy()
            x.shape = (np,)
            y.shape = (np,)
            z.shape = (np,)

            x = xmin + (xmax - xmin)*x
            y = ymin + (ymax - ymin)*y
            z = zmin + (zmax - zmin)*z

        vx = SpRandom(vxmean, vthx, np)
        vy = SpRandom(vymean, vthy, np)
        vz = SpRandom(vzmean, vthz, np)
        if lmomentum:
            gi = 1./sqrt(1. + (vx*vx + vy*vy + vz*vz)/clight**2)
        else:
            gi = 1.

        kw['lallindomain'] = lallindomain
        kw['lmomentum'] = lmomentum
        return self.addparticles(x, y, z, vx, vy, vz, js=js, gi=gi, **kw)

    def add_uniform_cylinder(self, np, rmax, zmin, zmax, vthx=0., vthy=0., vthz=0.,
                             xmean=0., ymean=0, zmean=0, vxmean=0., vymean=0., vzmean=0.,
                             theta=0., phi=0., vrmax=0., vtheta=0.,
                             js=None,
                             lmomentum=0, spacing='random', nr=None, nz=None,
                             thetamin=0., thetamax=2.*pi,
                             lvariableweights=None, lallindomain=0,
                             solvergeom=None,
                             ellipticity=1., wscale=1., rmin=None,
                             **kw):
        """Creates particles, uniformly filling a cylinder.
    If top.wpid is nonzero, then the particles are uniformly spaced in radius and
    the weights are set appropriately (weight=r/rmax). Otherwise, the particles are spaced uniformly
    in radius squared.
     - np: total number of particles to load
     - rmax: radius of cylinder
     - rmin: minimum radius, below which no particles are added
     - zmin, zmax: z extent of the cylinder. Note that if one has zmin == zmax,
                  then lallindomain must be set to true to get any particles.
     - vthx, vthy, vthz: thermal velocity, defaults to 0.
     - xmean, ymean, zmean: center of the cylinder, defaults to 0.
     - vxmean, vymean, vzmean: directed velocity, defaults to 0.
     - theta, phi: angle of cylinder about the center,
                  The transformation to lab frame is a rotation about the x axis by phi,
                  followed by a rotation about the new y axis by theta.
                  Note that the position mins and maxs and the velocity averages and means
                  are relative to the rotated frame and are automatically transformed
                  into the lab frame. The position means are relative to the lab frame.
     - vrmax: diverging (or converging) velocity. A velocity of vr = r*vrmax/rmax
              is added to each particle.
     - vtheta: rotational velocity about the cylinder center [rad/s].
     - js: particle species number, don't set it unless you mean it
     - lmomentum=false: Set to false when velocities are input as velocities, true
                        when input as massless momentum (as WARP stores them).
                        Only used when top.lrelativ is true.
     - spacing='random': either 'random' or 'uniform' particle spacing. For uniform,
                         r and z are uniform, theta is still random
     - nr, nz: for 'uniform' spacing, number of particles along r and z axis
               These will be rounded to the nearest integer.
     - thetamin=0., thetamax=2.*pi: range of theta around the cylinder
     - lvariableweights: By default, if wpid is set, then the particles will be
                         weighted according to their radius, otherwise all will
                         have the same weight. Use this option to override the
                         default.
     - wscale=1.: Scale factor applied to the particle weights. Only used if
                  variable weights are being used.
     - lallindomain=0: If true, the code only loads particles within the domain.
                       This only matters when parallel.
     - solvergeom=w3d.solvergeom: Specifies the geometry to use
     - ellipticity: factor multiplying size in y, i.e. y-size is given by rmax*ellipticity
        """
        if np == 0:
            if top.debug:
                print("add_uniform_cylinder: Warning: no particles loaded since the input np is zero")
            return

        if lallindomain:
            # --- Crop the min and max to be within the local domain
            xdmin, xdmax = top.xpminlocal, top.xpmaxlocal
            ydmin, ydmax = top.ypminlocal, top.ypmaxlocal
            zdmin, zdmax = top.zpminlocal + top.zgrid, top.zpmaxlocal + top.zgrid
        else:
            # --- Crop the min and max to be within the global domain
            xdmin, xdmax = w3d.xmmin, w3d.xmmax
            ydmin, ydmax = w3d.ymmin, w3d.ymmax
            zdmin, zdmax = w3d.zmmin + top.zgrid, w3d.zmmax + top.zgrid

        if solvergeom is None:
            solvergeom = w3d.solvergeom
        if solvergeom == w3d.XZgeom:
            ydmin, ydmax = -largepos, +largepos

        # --- Precalculate the trigonometrics if needed.
        if theta != 0. or phi != 0.:
            ct = cos(theta)
            st = sin(theta)
            cp = cos(phi)
            sp = sin(phi)

        # --- Check if the loading volume is within the grid. If there is no
        # --- overlap, the just exit. This is primarily for parallel optimization.
        # --- Get the coordinates of the corners of the box that encloses
        # --- the cylinder
        xx = array([-rmax, +rmax, -rmax, +rmax, -rmax, +rmax, -rmax, +rmax])
        yy = array([-rmax, -rmax, +rmax, +rmax, -rmax, -rmax, +rmax, +rmax])
        zz = array([zmin, zmin, zmin, zmin, zmax, zmax, zmax, zmax])

        if theta != 0. or phi != 0.:
            # --- Transform the corners from the rotated frame to the lab frame.
            xx1 = +xx*ct - yy*st*sp + zz*st*cp
            yy1 =        + yy*cp    + zz*sp
            zz1 = -xx*st - yy*ct*sp + zz*ct*cp
            xx, yy, zz = xx1, yy1, zz1

        # --- Shift by the mean and handle symmetries.
        xx += xmean
        yy += ymean
        zz += zmean
        if w3d.l4symtry:
            xx = abs(xx)
        if w3d.l2symtry or w3d.l4symtry:
            yy = abs(yy)

        # --- If outside the domain, just return.
        if min(xx) > xdmax or max(xx) < xdmin:
            if top.debug:
                print("add_uniform_cylinder: Warning: no particles loaded since all " +
                      "are outside of the x range of the domain of processor %d"%me)
            return
        if min(yy) > ydmax or max(yy) < ydmin:
            if top.debug:
                print("add_uniform_cylinder: Warning: no particles loaded since all " +
                      "are outside of the y range of the domain of processor %d"%me)
            return
        if min(zz) > zdmax or max(zz) < zdmin:
            if top.debug:
                print("add_uniform_cylinder: Warning: no particles loaded since all " +
                      "are outside of the z range of the domain of processor %d"%me)
            return

        if theta == 0. and phi == 0.:
            # --- When no angle is specified, then the clipping to the domain
            # --- can be done directly on the zmin and zmax
            zminp = max(zmin + zmean, min(zmax + zmean, zdmin)) - zmean
            zmaxp = min(zmax + zmean, max(zmin + zmean, zdmax)) - zmean
        else:
            # --- When angles are specified, the clipping must be done on a particle
            # --- by particle basis. Copy the zmin and zmax directly over to start.
            zminp = zmin
            zmaxp = zmax

        if spacing == 'random':
            # --- Add a random number to the number of particles so that on
            # --- average, the correct number of particles will be generated.
            np = int(np + random.random())
            if zmax != zmin:
                # --- Adjust the number of particles to load to based on the
                # --- width of the cropped zmin and max and the original
                np = nint((zmaxp - zminp)/(zmax - zmin)*np)
                # --- If the region is zero length, zmin == zmax, then don't change np.
            if np == 0:
                if top.debug:
                    print("add_uniform_cylinder: Warning: no particles loaded for " +
                          "random spacing since all are outside of the z range of the " +
                          "domain of processor %d"%me)
                return
            r = random.random(np)
            if zmax != zmin:
                z = random.random(np)
                z = (zminp + (zmaxp - zminp)*z - zmin)/(zmax - zmin)
            else:
                z = ones(np)
        else:
            if zmax != zmin:
                if nr is None:
                    nr = np**(1./2.)
                if nz is None:
                    nz = np**(1./2.)
                nr = nint(nr)
                nz = nint(nz)

                # --- Find the range of particle z locations within the cropped
                # --- zmin and max.
                dz = (zmax - zmin)/nz
                izminp = int((zminp - zmin)/dz + 0.5)
                izmaxp = int((zmaxp - zmin)/dz + 0.5)
                nzp = max(0, izmaxp - izminp)
                np = nr*nzp
                if np == 0:
                    if top.debug:
                        print("add_uniform_cylinder: Warning: no particles loaded for" +
                              "uniform spacing since all are outside of the z range of " +
                              "the domain of processor %d"%me)
                    return

                r, z = getmesh2d(0.5/nr, 1./nr, nr-1,
                                 (izminp + 0.5)/nz, 1./nz, nzp-1)
            else:  # zmax == zmin
                if nr is None:
                    nr = np
                nr = nint(nr)
                np = nr
                r = getmesh1d(0.5/nr, 1./nr, nr-1)
                z = ones(np)

            # --- Perform a transpose so that the data is ordered with increasing z.
            # --- The copy is needed since transposed arrays cannot be reshaped.
            r = transpose(r).copy()
            z = transpose(z).copy()
            r.shape = (np,)
            z.shape = (np,)

        thetap = (thetamax - thetamin)*random.random(np) + thetamin

        if lvariableweights is None:
            lvariableweights = (top.wpid != 0)
        if lvariableweights and top.wpid == 0:
            top.wpid = nextpid()

        if not lvariableweights:
            # --- Scale the radius appriately to get a uniform distribution
            # --- in radius.
            if rmin is not None:
                r = r*(1. - (rmin/rmax)**2) + (rmin/rmax)**2
            r = sqrt(r)
        else:
            if rmin is not None:
                r = r*(1. - (rmin/rmax)) + (rmin/rmax)

        x = rmax*r*cos(thetap)
        y = rmax*r*sin(thetap)*ellipticity
        if zmin != zmax:
            z = zmin + (zmax - zmin)*z
        else:
            z = zmin*z

        if theta != 0. or phi != 0.:
            # --- Transform positions from rotated frame into the lab frame.
            x1 = +x*ct - y*st*sp + z*st*cp
            y1 =       + y*cp    + z*sp
            z1 = -x*st - y*ct*sp + z*ct*cp
            x, y, z = x1, y1, z1

        # --- Now add the means, after the rotation transform.
        x += xmean
        y += ymean
        z += zmean

        # --- Clip transversely if needed.
        if top.nxprocs > 1 or top.nyprocs > 1:
            xx = x
            yy = y
            if w3d.l4symtry:
                xx = abs(xx)
            if w3d.l2symtry or w3d.l4symtry:
                yy = abs(yy)
            indomain = logical_and(logical_and(xdmin <= xx, xx < xdmax),
                                   logical_and(ydmin <= yy, yy < ydmax))
            x = x[indomain]
            y = y[indomain]
            z = z[indomain]
            if vtheta != 0.:
                r = r[indomain]
            np = len(z)
            if np == 0:
                if top.debug:
                    print("add_uniform_cylinder: Warning: no particles loaded since " +
                          "all are outside of the transverse range of the domain of " +
                          "processor %d"%me)
                return

        if theta != 0. or phi != 0.:
            # --- When angles are specified, the clipping in z must be done on a
            # --- particle by particle basis.
            indomain = logical_and(zdmin <= z, z < zdmax)
            x = x[indomain]
            y = y[indomain]
            z = z[indomain]
            if vtheta != 0.:
                r = r[indomain]
            np = len(z)
            if np == 0:
                if top.debug:
                    print("add_uniform_cylinder: Warning: no particles loaded after " +
                          "rotation since all are outside of the z range of the " +
                          "domain of processor %d"%me)
                return

        # --- The weights can be set now if needed (after clipping the positions).
        # --- Note that rmax must be scaled out.
        if not lvariableweights:
            if top.wpid != 0:
                kw['w'] = 1.*wscale
        else:
            kw['w'] = 2*sqrt((x - xmean)**2 + (y - ymean)**2)/rmax*wscale

        # --- Now the velocities are generated (after clipping the positions).
        vx = SpRandom(vxmean, vthx, np)
        vy = SpRandom(vymean, vthy, np)
        vz = SpRandom(vzmean, vthz, np)

        if vrmax != 0.:
            vx += vrmax*x/rmax
            vy += vrmax*y/rmax
        if vtheta != 0.:
            vx -= vtheta*y
            vy += vtheta*x

        if theta != 0. or phi != 0.:
            # --- Transform velocities from rotated frame into the lab frame.
            vx1 = +vx*ct - vy*st*sp + vz*st*cp
            vy1 =        + vy*cp    + vz*sp
            vz1 = -vx*st - vy*ct*sp + vz*ct*cp
            vx, vy, vz = vx1, vy1, vz1

        if lmomentum:
            gi = 1./sqrt(1. + (vx*vx + vy*vy + vz*vz)/clight**2)
        else:
            gi = 1

        kw['lallindomain'] = lallindomain
        kw['lmomentum'] = lmomentum
        return self.addparticles(x, y, z, vx, vy, vz, js=js, gi=gi, **kw)

    def add_gaussian_dist(self, np, deltax, deltay, deltaz, vthx=0., vthy=0., vthz=0.,
                          xmean=0., ymean=0., zmean=0., vxmean=0., vymean=0., vzmean=0.,
                          vxdiv=0., vydiv=0., vzdiv=0.,
                          zdist='random', rdist='linear', nz=1000, fourfold=0, js=None, lmomentum=0, w=None, **kw):
        """
    Add particles with a Gaussian distribution.
     - np: total number of particles to load
     - deltax, deltay, deltaz: width of the distribution
     - vthx, vthy, vthz: thermal velocity, defaults to 0.
     - xmean, ymean, zmean: center of the cylinder, defaults to 0.
     - vxmean, vymean, vzmean: directed velocity, defaults to 0.
     - vxdiv, vydiv, vzdiv: Velocity divergence, defaults to 0.
     - zdist='random': type of random distribution along z, possible values are
                       'random' and 'regular'
     - rdist='linear': type of random distribution along r, possible values are
                       'linear' and 'flat'
     - nz=1000: number of data points to use along z when zdist='regular'
     - fourfold=False: whether to use four fold symmetry.
     - js: particle species number, don't set it unless you mean it
     - lmomentum=false: Set to false when velocities are input as velocities, true
                        when input as massless momentum (as WARP stores them).
                        Only used when top.lrelativ is true.
     - w=None: Particle weight (when top.wpid > 0). With rdist='flat', this will be multiplied by
               a radially dependent factor.
    Note that the lreturndata option doesn't work.
        """
        if fourfold:
            np = nint(float(np)/4)
        kw['lmomentum'] = lmomentum
        if zdist == 'random':
            if rdist == 'linear': # for uniform noise in XY.
                x = SpRandom(0., deltax, np)
                y = SpRandom(0., deltay, np)
                wp = w
            else: # rdist = flat, for uniform noise in R, need to adjust weight accordingly.
                rmax = 4.
                r = rmax*random.random(np)
                t = 2*pi*random.random(np)
                x = r*cos(t)*deltax
                y = r*sin(t)*deltay
                wp = rmax*r*exp(-0.5*r**2)
                if w is not None:
                    wp *= w
                if top.wpid==0:top.wpid=nextpid()

            z = SpRandom(0., deltaz, np)
            vx = SpRandom(0., vthx, np)
            vy = SpRandom(0., vthy, np)
            vz = SpRandom(0., vthz, np)
            if fourfold:
                xsigns = [1., -1.]
                ysigns = [1., -1.]
            else:
                xsigns = [1.]
                ysigns = [1.]
            za = z + zmean
            vza = vz + vzmean + vzdiv*z
            for ysign in ysigns:
                for xsign in xsigns:
                    xa = xsign*x + xmean
                    ya = ysign*y + ymean
                    vxa = xsign*vx + vxmean + vxdiv*xsign*x
                    vya = ysign*vy + vymean + vydiv*ysign*y
                    if lmomentum:
                        gi = 1./sqrt(1. + (vxa*vxa + vya*vya + vza*vza)/clight**2)
                    else:
                        gi = 1.
                    self.addparticles(xa, ya, za, vxa, vya, vza, gi=gi, js=js, w=wp, **kw)
        elif zdist == 'regular':
            dz = 16.*deltaz/nz
            zmin = -(float(nz/2) - 0.5)*dz
            for i in range(nz):
                zadd = zmin + i*dz
                Nadd = max(0., np*dz*exp(-0.5*(zadd/deltaz)**2)/(sqrt(2.*pi)*deltaz))
                if ranf() < (Nadd - int(Nadd)):
                    Nadd = int(Nadd) + 1
                else:
                    Nadd = int(Nadd)
                if Nadd > 0:
                    if rdist == 'linear': # for uniform noise in XY.
                        x = SpRandom(0., deltax, Nadd)
                        y = SpRandom(0., deltay, Nadd)
                        wp = w
                    else: # rdist = flat, for uniform noise in R, need to adjust weight accordingly.
                        rmax = 4.*sqrt(max(0.,(1.-(zadd/(pi*deltaz))**2)))
                        r = (arange(Nadd)+0.5)*rmax/Nadd
                        t = 2*pi*random.random(Nadd)
                        x = r*cos(t)*deltax
                        y = r*sin(t)*deltay
                        wp = rmax*r*exp(-0.5*r**2)
                        if w is not None:
                            wp *= w
                        if top.wpid==0:top.wpid=nextpid()

                    z = zadd + dz*(random.random(Nadd) - 0.5)
                    vx = SpRandom(0., vthx, Nadd)
                    vy = SpRandom(0., vthy, Nadd)
                    vz = SpRandom(0., vthz, Nadd)
                    if fourfold:
                        xsigns = [1., -1.]
                        ysigns = [1., -1.]
                    else:
                        xsigns = [1.]
                        ysigns = [1.]
                    za = z + zmean
                    vza = vz + vzmean + vzdiv*z
                    for ysign in ysigns:
                        for xsign in xsigns:
                            xa = xsign*x + xmean
                            ya = ysign*y + ymean
                            vxa = xsign*vx + vxmean + vxdiv*xsign*x
                            vya = ysign*vy + vymean + vydiv*ysign*y
                            if lmomentum:
                                gi = 1./sqrt(1. + (vxa*vxa + vya*vya + vza*vza)/clight**2)
                            else:
                                gi = 1.
                            self.addparticles(xa, ya, za, vxa, vya, vza, gi=gi, js=js, w=wp, **kw)

    def gather_zmmnts_locs(self):
        get_zmmnts_stations(len(self.jslist),
                            array(self.jslist),
                            self.pgroup,
                            self.nzmmnts_locs,
                            self.zmmnts_locs[0],
                            self.zmmnts_locs[-1],
                            self.zmmnts_locs_vfrm,
                            self.zmmnts_locs_pnum,
                            self.zmmnts_locs_xbar,
                            self.zmmnts_locs_ybar,
                            self.zmmnts_locs_xpbar,
                            self.zmmnts_locs_ypbar,
                            self.zmmnts_locs_x2,
                            self.zmmnts_locs_y2,
                            self.zmmnts_locs_xp2,
                            self.zmmnts_locs_yp2,
                            self.zmmnts_locs_xxp,
                            self.zmmnts_locs_yyp)
        self.zmmnts_gathered = false

    def gatherall_zmmnts_locs(self):
        self.zmmnts_locs_pnum = parallelsum(self.zmmnts_locs_pnum)
        self.zmmnts_locs_xbar = parallelsum(self.zmmnts_locs_xbar)
        self.zmmnts_locs_ybar = parallelsum(self.zmmnts_locs_ybar)
        self.zmmnts_locs_xpbar = parallelsum(self.zmmnts_locs_xpbar)
        self.zmmnts_locs_ypbar = parallelsum(self.zmmnts_locs_ypbar)
        self.zmmnts_locs_x2 = parallelsum(self.zmmnts_locs_x2)
        self.zmmnts_locs_y2 = parallelsum(self.zmmnts_locs_y2)
        self.zmmnts_locs_xp2 = parallelsum(self.zmmnts_locs_xp2)
        self.zmmnts_locs_yp2 = parallelsum(self.zmmnts_locs_yp2)
        self.zmmnts_locs_xxp = parallelsum(self.zmmnts_locs_xxp)
        self.zmmnts_locs_yyp = parallelsum(self.zmmnts_locs_yyp)
        if me > 0:
            self.zmmnts_locs_pnum[...] = 0.
            self.zmmnts_locs_xbar[...] = 0.
            self.zmmnts_locs_ybar[...] = 0.
            self.zmmnts_locs_xpbar[...] = 0.
            self.zmmnts_locs_ypbar[...] = 0.
            self.zmmnts_locs_x2[...] = 0.
            self.zmmnts_locs_y2[...] = 0.
            self.zmmnts_locs_xp2[...] = 0.
            self.zmmnts_locs_yp2[...] = 0.
            self.zmmnts_locs_xxp[...] = 0.
            self.zmmnts_locs_yyp[...] = 0.
        self.zmmnts_locs_gathered = true

    def set_zmmnts_locs(self, zmin, zmax, nz, vfrm=0.):
        self.zmmnts_locs = getmeshcoordinates([zmin], [(zmax - zmin)/(nz - 1)], [nz - 1])[0]
        self.nzmmnts_locs = nz
        self.zmmnts_locs_vfrm = vfrm
        self.dzmmnts_locs = (zmax - zmin)/(nz - 1)
        self.zmmnts_locs_pnum = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_xbar = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_ybar = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_zbar = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_xpbar = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_ypbar = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_x2 = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_y2 = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_z2 = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_xp2 = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_yp2 = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_xxp = zeros(self.nzmmnts_locs, 'd')
        self.zmmnts_locs_yyp = zeros(self.nzmmnts_locs, 'd')
        if top.xoldpid == 0: top.xoldpid = nextpid()
        if top.yoldpid == 0: top.yoldpid = nextpid()
        if top.zoldpid == 0: top.zoldpid = nextpid()
        if top.uxoldpid == 0: top.uxoldpid = nextpid()
        if top.uyoldpid == 0: top.uyoldpid = nextpid()
        if top.uzoldpid == 0: top.uzoldpid = nextpid()
        self.zmmnts_locs_gathered = true
        installafterstep(self.gather_zmmnts_locs)

    def set_zmmnts(self):
        self.zmmnts_pnum = AppendableArray(typecode='d')
        self.zmmnts_xbar = AppendableArray(typecode='d')
        self.zmmnts_ybar = AppendableArray(typecode='d')
        self.zmmnts_zbar = AppendableArray(typecode='d')
        self.zmmnts_xpbar = AppendableArray(typecode='d')
        self.zmmnts_ypbar = AppendableArray(typecode='d')
        self.zmmnts_xpnbar = AppendableArray(typecode='d')
        self.zmmnts_ypnbar = AppendableArray(typecode='d')
        self.zmmnts_x2 = AppendableArray(typecode='d')
        self.zmmnts_y2 = AppendableArray(typecode='d')
        self.zmmnts_z2 = AppendableArray(typecode='d')
        self.zmmnts_xp2 = AppendableArray(typecode='d')
        self.zmmnts_yp2 = AppendableArray(typecode='d')
        self.zmmnts_xxp = AppendableArray(typecode='d')
        self.zmmnts_yyp = AppendableArray(typecode='d')
        self.zmmnts_xpn2 = AppendableArray(typecode='d')
        self.zmmnts_ypn2 = AppendableArray(typecode='d')
        self.zmmnts_xxpn = AppendableArray(typecode='d')
        self.zmmnts_yypn = AppendableArray(typecode='d')
        self.zmmnts_gathered = true
        installafterstep(self.gather_zmmnts)

    def gather_zmmnts(self):
        print 'gather_zmmnts'
        self.zmmnts_pnum.append(self.getn(gather=0))
        xbar = 0.
        ybar = 0.
        zbar = 0.
        xpbar = 0.
        ypbar = 0.
        xpnbar = 0.
        ypnbar = 0.
        x2 = 0.
        y2 = 0.
        z2 = 0.
        xp2 = 0.
        yp2 = 0.
        xxpbar = 0.
        yypbar = 0.
        xpn2 = 0.
        ypn2 = 0.
        xxpnbar = 0.
        yypnbar = 0.
        pg = self.pgroup
        nparpgrp = 1024
        xpa = zeros(nparpgrp, 'd')
        ypa = zeros(nparpgrp, 'd')
        xpna = zeros(nparpgrp, 'd')
        ypna = zeros(nparpgrp, 'd')
        betaa = zeros(nparpgrp, 'd')
        for js in self.jslist:
            if pg.nps[js] == 0:
                continue
            ng = 1 + pg.nps[js]/nparpgrp
            for ig in range(ng):
                il = pg.ins[js] - 1 + nparpgrp*ig
                iu = min(il + nparpgrp, pg.ins[js] - 1 + pg.nps[js])
                np = iu - il
                x = pg.xp[il:iu]
                y = pg.yp[il:iu]
                z = pg.zp[il:iu]
                xpa[:np] = pg.uxp[il:iu]/pg.uzp[il:iu]
                ypa[:np] = pg.uyp[il:iu]/pg.uzp[il:iu]
                gaminv = pg.gaminv[il:iu]
                betaa[:np] = (1. - gaminv)*(1. + gaminv)
                xpna[:np] = xpa[:np]*betaa[:np]/gaminv
                ypna[:np] = ypa[:np]*betaa[:np]*gaminv
                xp = xpa[:np]
                yp = ypa[:np]
                xpn = xpna[:np]
                ypn = ypna[:np]
                beta = betaa[:np]
                xbar += sum(x)
                ybar += sum(y)
                zbar += sum(z)
                xpbar += sum(xp)
                ypbar += sum(yp)
                xpnbar += sum(xpn)
                ypnbar += sum(ypn)
                x2 += sum(x*x)
                y2 += sum(y*y)
                z2 += sum(z*z)
                xp2 += sum(xp*xp)
                yp2 += sum(yp*yp)
                xpn2 += sum(xpn*xpn)
                ypn2 += sum(ypn*ypn)
                xxpbar += sum(x*xp)
                yypbar += sum(y*yp)
                xxpnbar += sum(x*xpn)
                yypnbar += sum(y*ypn)
        self.zmmnts_xbar.append(xbar)
        self.zmmnts_ybar.append(ybar)
        self.zmmnts_zbar.append(zbar)
        self.zmmnts_xpbar.append(xpbar)
        self.zmmnts_ypbar.append(ypbar)
        self.zmmnts_x2.append(x2)
        self.zmmnts_y2.append(y2)
        self.zmmnts_z2.append(z2)
        self.zmmnts_xp2.append(xp2)
        self.zmmnts_yp2.append(yp2)
        self.zmmnts_xxp.append(xxpbar)
        self.zmmnts_yyp.append(yypbar)
        self.zmmnts_xpnbar.append(xpnbar)
        self.zmmnts_ypnbar.append(ypnbar)
        self.zmmnts_xpn2.append(xpn2)
        self.zmmnts_ypn2.append(ypn2)
        self.zmmnts_xxpn.append(xxpnbar)
        self.zmmnts_yypn.append(yypnbar)
        self.zmmnts_gathered = false
        print 'gather_zmmnts_done'

    def gatherall_zmmnts(self):
        self.zmmnts_pnum.data()[...] = parallelsum(self.zmmnts_pnum.data())
        self.zmmnts_xbar.data()[...] = parallelsum(self.zmmnts_xbar.data())
        self.zmmnts_ybar.data()[...] = parallelsum(self.zmmnts_ybar.data())
        self.zmmnts_xpbar.data()[...] = parallelsum(self.zmmnts_xpbar.data())
        self.zmmnts_ypbar.data()[...] = parallelsum(self.zmmnts_ypbar.data())
        self.zmmnts_xpnbar.data()[...] = parallelsum(self.zmmnts_xpnbar.data())
        self.zmmnts_ypnbar.data()[...] = parallelsum(self.zmmnts_ypnbar.data())
        self.zmmnts_x2.data()[...] = parallelsum(self.zmmnts_x2.data())
        self.zmmnts_y2.data()[...] = parallelsum(self.zmmnts_y2.data())
        self.zmmnts_xp2.data()[...] = parallelsum(self.zmmnts_xp2.data())
        self.zmmnts_yp2.data()[...] = parallelsum(self.zmmnts_yp2.data())
        self.zmmnts_xxp.data()[...] = parallelsum(self.zmmnts_xxp.data())
        self.zmmnts_yyp.data()[...] = parallelsum(self.zmmnts_yyp.data())
        self.zmmnts_xpn2.data()[...] = parallelsum(self.zmmnts_xpn2.data())
        self.zmmnts_ypn2.data()[...] = parallelsum(self.zmmnts_ypn2.data())
        self.zmmnts_xxpn.data()[...] = parallelsum(self.zmmnts_xxpn.data())
        self.zmmnts_yypn.data()[...] = parallelsum(self.zmmnts_yypn.data())
        if me > 0:
            self.zmmnts_pnum.data()[...] = 0.
            self.zmmnts_xbar.data()[...] = 0.
            self.zmmnts_ybar.data()[...] = 0.
            self.zmmnts_xpbar.data()[...] = 0.
            self.zmmnts_ypbar.data()[...] = 0.
            self.zmmnts_xpnbar.data()[...] = 0.
            self.zmmnts_ypnbar.data()[...] = 0.
            self.zmmnts_x2.data()[...] = 0.
            self.zmmnts_y2.data()[...] = 0.
            self.zmmnts_xp2.data()[...] = 0.
            self.zmmnts_yp2.data()[...] = 0.
            self.zmmnts_xxp.data()[...] = 0.
            self.zmmnts_yyp.data()[...] = 0.
            self.zmmnts_xpn2.data()[...] = 0.
            self.zmmnts_ypn2.data()[...] = 0.
            self.zmmnts_xxpn.data()[...] = 0.
            self.zmmnts_yypn.data()[...] = 0.
        self.zmmnts_gathered = true

    def getzmmnts_pnum(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            return self.zmmnts_pnum.data()
        else:
            return None

    def getzmmnts_xrms(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return sqrt(self.zmmnts_x2.data()/pnum)
        else:
            return None

    def getzmmnts_yrms(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return sqrt(self.zmmnts_y2.data()/pnum)
        else:
            return None

    def getzmmnts_xprms(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return sqrt(self.zmmnts_xp2.data()/pnum)
        else:
            return None

    def getzmmnts_yprms(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return sqrt(self.zmmnts_yp2.data()/pnum)
        else:
            return None

    def getzmmnts_xbar(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return self.zmmnts_xbar.data()/pnum
        else:
            return None

    def getzmmnts_ybar(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return self.zmmnts_ybar.data()/pnum
        else:
            return None

    def getzmmnts_xpbar(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return self.zmmnts_xpbar.data()/pnum
        else:
            return None

    def getzmmnts_ypbar(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return self.zmmnts_ypbar.data()/pnum
        else:
            return None

    def getzmmnts_xxpbar(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return self.zmmnts_xxp.data()/pnum
        else:
            return None

    def getzmmnts_yypbar(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
            return self.zmmnts_yyp.data()/pnum
        else:
            return None

    def getzmmnts_emitxrms(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me > 0:
            return None
        pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
        xbar = self.getzmmnts_xbar()
        xpbar = self.getzmmnts_xpbar()
        xxp = self.zmmnts_xxp.data()/pnum
        x2 = self.zmmnts_x2.data()/pnum
        xp2 = self.zmmnts_xp2.data()/pnum
        return sqrt((x2 - xbar*xbar)*(xp2 - xpbar*xpbar) - (xxp - xbar*xpbar)**2)

    def getzmmnts_emityrms(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me > 0:
            return None
        pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
        ybar = self.getzmmnts_ybar()
        ypbar = self.getzmmnts_ypbar()
        yyp = self.zmmnts_yyp.data()/pnum
        y2 = self.zmmnts_y2.data()/pnum
        yp2 = self.zmmnts_yp2.data()/pnum
        return sqrt((y2 - ybar*ybar)*(yp2 - ypbar*ypbar) - (yyp - ybar*ypbar)**2)

    def getzmmnts_emitxnrms(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me > 0:
            return None
        pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
        xbar = self.getzmmnts_xbar()
        xpnbar = self.zmmnts_xpnbar.data()/pnum
        xxpn = self.zmmnts_xxpn.data()/pnum
        x2 = self.zmmnts_x2.data()/pnum
        xpn2 = self.zmmnts_xpn2.data()/pnum
        return sqrt((x2 - xbar*xbar)*(xpn2 - xpnbar*xpnbar) - (xxpn - xbar*xpnbar)**2)

    def getzmmnts_emitynrms(self):
        if not self.zmmnts_gathered:
            self.gatherall_zmmnts()
        if me > 0:
            return None
        pnum = where(self.zmmnts_pnum.data() > smallpos, self.zmmnts_pnum.data(), 1.)
        ybar = self.getzmmnts_ybar()
        ypnbar = self.zmmnts_ypnbar.data()/pnum
        yypn = self.zmmnts_yypn.data()/pnum
        y2 = self.zmmnts_y2.data()/pnum
        ypn2 = self.zmmnts_ypn2.data()/pnum
        return sqrt((y2 - ybar*ybar)*(ypn2 - ypnbar*ypnbar) - (yypn - ybar*ypnbar)**2)

    def plzmmnts_data(self, x, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        if me == 0:
            zst = self.zmmnts_locs
            pla(yscale*x, xscale*(zst - xoffset), color=color, width=width)

    def plzmmnts_xbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_xbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_ybar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_ybar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_xrms(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_xrms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_yrms(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_yrms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_xprms(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_xprms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_yprms(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_yprms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_xpbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_xpbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_ypbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_ypbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_xxpbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_xxpbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_yypbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_yypbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_emitx(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_emitxrms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_emity(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_data(self.getzmmnts_emityrms(), color, width, type, xscale, yscale, xoffset)

    def getzmmnts_locs_pnum(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            return self.zmmnts_locs_pnum
        else:
            return None

    def getzmmnts_locs_xrms(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return sqrt(self.zmmnts_locs_x2/pnum)
        else:
            return None

    def getzmmnts_locs_yrms(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return sqrt(self.zmmnts_locs_y2/pnum)
        else:
            return None

    def getzmmnts_locs_xprms(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return sqrt(self.zmmnts_locs_xp2/pnum)
        else:
            return None

    def getzmmnts_locs_yprms(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return sqrt(self.zmmnts_locs_yp2/pnum)
        else:
            return None

    def getzmmnts_locs_xbar(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return self.zmmnts_locs_xbar/pnum
        else:
            return None

    def getzmmnts_locs_ybar(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return self.zmmnts_locs_ybar/pnum
        else:
            return None

    def getzmmnts_locs_xpbar(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return self.zmmnts_locs_xpbar/pnum
        else:
            return None

    def getzmmnts_locs_ypbar(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return self.zmmnts_locs_ypbar/pnum
        else:
            return None

    def getzmmnts_locs_xxpbar(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return self.zmmnts_locs_xxp/pnum
        else:
            return None

    def getzmmnts_locs_yypbar(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me == 0:
            pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
            return self.zmmnts_locs_yyp/pnum
        else:
            return None

    def getzmmnts_locs_emitxrms(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me > 0:
            return None
        pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
        xbar = self.getzmmnts_locs_xbar()
        xpbar = self.getzmmnts_locs_xpbar()
        xxp = self.zmmnts_locs_xxp/pnum
        x2 = self.zmmnts_locs_x2/pnum
        xp2 = self.zmmnts_locs_xp2/pnum
        return sqrt((x2 - xbar*xbar)*(xp2 - xpbar*xpbar) - (xxp - xbar*xpbar)**2)

    def getzmmnts_locs_emityrms(self):
        if not self.zmmnts_locs_gathered:
            self.gatherall_zmmnts()
        if me > 0:
            return None
        pnum = where(self.zmmnts_locs_pnum > smallpos, self.zmmnts_locs_pnum, 1.)
        ybar = self.getzmmnts_locs_ybar()
        ypbar = self.getzmmnts_locs_ypbar()
        yyp = self.getzmmnts_locs_yypbar()
        y2 = self.zmmnts_locs_y2/pnum
        yp2 = self.zmmnts_locs_yp2/pnum
        return sqrt((y2 - ybar*ybar)*(yp2 - ypbar*ypbar) - (yyp - ybar*ypbar)**2)

    def plzmmnts_locs_data(self, x, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        if me == 0:
            zst = self.zmmnts_locs
            pla(yscale*x, xscale*(zst - xoffset), color=color, width=width)

    def plzmmnts_locs_xbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_xbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_ybar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_ybar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_xrms(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_xrms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_yrms(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_yrms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_xprms(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_xprms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_yprms(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_yprms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_xpbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_xpbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_ypbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_ypbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_xxpbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_xxpbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_yypbar(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_yypbar(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_emitx(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_emitxrms(), color, width, type, xscale, yscale, xoffset)

    def plzmmnts_locs_emity(self, color=black, width=1, type='solid', xscale=1., yscale=1., xoffset=0.):
        self.plzmmnts_locs_data(self.getzmmnts_locs_emityrms(), color, width, type, xscale, yscale, xoffset)

    def selectparticles(self, **kw):
        # --- This adjusts the indices returned by selectparticles so that they
        # --- are relative to start of the species. The ii returned can then be
        # --- used for example like self.xp[ii]. Note that this won't work if
        # --- jslist has includes multiple species.
        ii = selectparticles(jslist=self.jslist, pgroup=self.pgroup, **kw)
        i1 = self.pgroup.ins[self.jslist[0]] - 1
        if isinstance(ii, slice):
            ii = slice(ii.start - i1, ii.stop - i1)
        else:
            ii = ii - i1
        return ii

    def getn(self, **kw):
        """Calls :py:func:`~particles.getn` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getn(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getn(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getx(self, **kw):
        """Calls :py:func:`~particles.getx` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getx(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getx(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def gety(self, **kw):
        """Calls :py:func:`~particles.gety` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return gety(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return gety(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getz(self, **kw):
        """Calls :py:func:`~particles.getz` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getz(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getz(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getr(self, **kw):
        """Calls :py:func:`~particles.getr` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getr(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getr(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def gettheta(self, **kw):
        """Calls :py:func:`~particles.gettheta` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return gettheta(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return gettheta(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getvx(self, **kw):
        """Calls :py:func:`~particles.getvx` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getvx(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getvx(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getvy(self, **kw):
        """Calls :py:func:`~particles.getvy` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getvy(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getvy(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getvz(self, **kw):
        """Calls :py:func:`~particles.getvz` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getvz(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getvz(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getvr(self, **kw):
        """Calls :py:func:`~particles.getvr` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getvr(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getvr(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getvtheta(self, **kw):
        """Calls :py:func:`~particles.getvtheta` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getvtheta(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getvtheta(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getux(self, **kw):
        """Calls :py:func:`~particles.getux` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getux(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getux(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getuy(self, **kw):
        """Calls :py:func:`~particles.getuy` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getuy(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getuy(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getuz(self, **kw):
        """Calls :py:func:`~particles.getuz` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getuz(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getuz(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getex(self, **kw):
        """Calls :py:func:`~particles.getex` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getex(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getex(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getey(self, **kw):
        """Calls :py:func:`~particles.getey` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getey(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getey(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getez(self, **kw):
        """Calls :py:func:`~particles.getez` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getez(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getez(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def geter(self, **kw):
        """Calls :py:func:`~particles.geter` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return geter(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return geter(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getetheta(self, **kw):
        """Calls :py:func:`~particles.getetheta` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getetheta(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getetheta(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getbx(self, **kw):
        """Calls :py:func:`~particles.getbx` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getbx(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getbx(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getby(self, **kw):
        """Calls :py:func:`~particles.getby` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getby(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getby(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getbz(self, **kw):
        """Calls :py:func:`~particles.getbz` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getbz(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getbz(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getbr(self, **kw):
        """Calls :py:func:`~particles.getbr` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getbr(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getbr(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getbtheta(self, **kw):
        """Calls :py:func:`~particles.getbtheta` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getbtheta(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getbtheta(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getxp(self, **kw):
        """Calls :py:func:`~particles.getxp` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getxp(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getxp(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getyp(self, **kw):
        """Calls :py:func:`~particles.getyp` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getyp(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getyp(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getrp(self, **kw):
        """Calls :py:func:`~particles.getrp` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getrp(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getrp(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def gettp(self, **kw):
        """Calls :py:func:`~particles.gettp` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return gettp(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return gettp(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getgaminv(self, **kw):
        """Calls :py:func:`~particles.getgaminv` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getgaminv(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getgaminv(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getpid(self, **kw):
        """Calls :py:func:`~particles.getpid` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getpid(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getpid(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getw(self, **kw):
        """Calls :py:func:`~particles.getw` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getw(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getw(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def getweights(self, **kw):
        """Returns the weights of the particles, the product of pid[:, wpid]*sw."""
        if top.wpid != 0:
            # If space has been allocated for the weights, use them
            if self.pgroups is not None:
                listofpgroups = self.flatten(self.pgroups)
                return getw(jslist=self.jslist, pgroups=listofpgroups, **kw)*self.sw
            else:
                return getw(jslist=self.jslist, pgroup=self.pgroup, **kw)*self.sw
        else:
            # Otherwise, return an array of uniform weights
            return ones_like(self.getx(**kw))*self.sw

    def getke(self, **kw):
        """Calls :py:func:`~particles.getke` for this species."""
        if self.pgroups is not None:
            listofpgroups = self.flatten(self.pgroups)
            return getke(jslist=self.jslist, pgroups=listofpgroups, **kw)
        else:
            return getke(jslist=self.jslist, pgroup=self.pgroup, **kw)

    def _callppfunc(self, ppfunc, pgroup=None, **kw):
        """This is an intermediary for all of the pp particle plot methods. This
        makes it easier to make changes to all of them at once, without adding alot
        of code."""
        if pgroup is None:
            if self.pgroups is not None:
                for pg in self.flatten(self.pgroups):
                    self._callppfunc(ppfunc, pg, **kw)
                return
            else:
                pgroup = self.pgroup

        kw.setdefault('color', self.color)
        kw.setdefault('marker', self.marker)
        kw.setdefault('msize', self.msize)
        return ppfunc(jslist=self.jslist, pgroup=pgroup, **kw)

    def ppxy(self, **kw):
        """Calls :py:func:`~warpplots.ppxy` for this species."""
        return self._callppfunc(ppxy, **kw)

    def ppxxp(self, **kw):
        """Calls :py:func:`~warpplots.ppxxp` for this species."""
        return self._callppfunc(ppxxp, **kw)

    def ppyyp(self, **kw):
        """Calls :py:func:`~warpplots.ppyyp` for this species."""
        return self._callppfunc(ppyyp, **kw)

    def ppxpyp(self, **kw):
        """Calls :py:func:`~warpplots.ppxpyp` for this species."""
        return self._callppfunc(ppxpyp, **kw)

    def ppxvx(self, **kw):
        """Calls :py:func:`~warpplots.ppxvx` for this species."""
        return self._callppfunc(ppxvx, **kw)

    def ppyvy(self, **kw):
        """Calls :py:func:`~warpplots.ppyvy` for this species."""
        return self._callppfunc(ppyvy, **kw)

    def ppxvz(self, **kw):
        """Calls :py:func:`~warpplots.ppxvz` for this species."""
        return self._callppfunc(ppxvz, **kw)

    def ppyvz(self, **kw):
        """Calls :py:func:`~warpplots.ppyvz` for this species."""
        return self._callppfunc(ppyvz, **kw)

    def ppzxy(self, **kw):
        """Calls :py:func:`~warpplots.ppzxy` for this species."""
        return self._callppfunc(ppzxy, **kw)

    def ppzx(self, **kw):
        """Calls :py:func:`~warpplots.ppzx` for this species."""
        return self._callppfunc(ppzx, **kw)

    def ppzy(self, **kw):
        """Calls :py:func:`~warpplots.ppzy` for this species."""
        return self._callppfunc(ppzy, **kw)

    def ppzr(self, **kw):
        """Calls :py:func:`~warpplots.ppzr` for this species."""
        return self._callppfunc(ppzr, **kw)

    def ppzxp(self, **kw):
        """Calls :py:func:`~warpplots.ppzxp` for this species."""
        return self._callppfunc(ppzxp, **kw)

    def ppzvx(self, **kw):
        """Calls :py:func:`~warpplots.ppzvx` for this species."""
        return self._callppfunc(ppzvx, **kw)

    def ppzyp(self, **kw):
        """Calls :py:func:`~warpplots.ppzyp` for this species."""
        return self._callppfunc(ppzyp, **kw)

    def ppzvy(self, **kw):
        """Calls :py:func:`~warpplots.ppzvy` for this species."""
        return self._callppfunc(ppzvy, **kw)

    def ppzvz(self, **kw):
        """Calls :py:func:`~warpplots.ppzvz` for this species."""
        return self._callppfunc(ppzvz, **kw)

    def ppxux(self, **kw):
        """Calls :py:func:`~warpplots.ppxux` for this species."""
        return self._callppfunc(ppxux, **kw)

    def ppyuy(self, **kw):
        """Calls :py:func:`~warpplots.ppyuy` for this species."""
        return self._callppfunc(ppyuy, **kw)

    def ppzux(self, **kw):
        """Calls :py:func:`~warpplots.ppzux` for this species."""
        return self._callppfunc(ppzux, **kw)

    def ppzuy(self, **kw):
        """Calls :py:func:`~warpplots.ppzuy` for this species."""
        return self._callppfunc(ppzuy, **kw)

    def ppzuz(self, **kw):
        """Calls :py:func:`~warpplots.ppzuz` for this species."""
        return self._callppfunc(ppzuz, **kw)

    def ppzrp(self, **kw):
        """Calls :py:func:`~warpplots.ppzrp` for this species."""
        return self._callppfunc(ppzrp, **kw)

    def ppzke(self, **kw):
        """Calls :py:func:`~warpplots.ppzke` for this species."""
        return self._callppfunc(ppzke, **kw)

    def ppzvr(self, **kw):
        """Calls :py:func:`~warpplots.ppzvr` for this species."""
        return self._callppfunc(ppzvr, **kw)

    def ppzvtheta(self, **kw):
        """Calls :py:func:`~warpplots.ppzvtheta` for this species."""
        return self._callppfunc(ppzvtheta, **kw)

    def ppzvperp(self, **kw):
        """Calls :py:func:`~warpplots.ppzvperp` for this species."""
        return self._callppfunc(ppzvperp, **kw)

    def ppvzvperp(self, **kw):
        """Calls :py:func:`~warpplots.ppvzvperp` for this species."""
        return self._callppfunc(ppvzvperp, **kw)

    def pptrace(self, **kw):
        """Calls :py:func:`~warpplots.pptrace` for this species."""
        return self._callppfunc(pptrace, **kw)

    def pprrp(self, **kw):
        """Calls :py:func:`~warpplots.pprrp` for this species."""
        return self._callppfunc(pprrp, **kw)

    def pprtp(self, **kw):
        """Calls :py:func:`~warpplots.pprtp` for this species."""
        return self._callppfunc(pprtp, **kw)

    def pprvr(self, **kw):
        """Calls :py:func:`~warpplots.pprvr` for this species."""
        return self._callppfunc(pprvr, **kw)

    def pprvz(self, **kw):
        """Calls :py:func:`~warpplots.pprvz` for this species."""
        return self._callppfunc(pprvz, **kw)

    def ppxex(self, **kw):
        """Calls :py:func:`~warpplots.ppxex` for this species."""
        return self._callppfunc(ppxex, **kw)

    def ppxey(self, **kw):
        """Calls :py:func:`~warpplots.ppxey` for this species."""
        return self._callppfunc(ppxey, **kw)

    def ppxez(self, **kw):
        """Calls :py:func:`~warpplots.ppxez` for this species."""
        return self._callppfunc(ppxez, **kw)

    def ppyex(self, **kw):
        """Calls :py:func:`~warpplots.ppyex` for this species."""
        return self._callppfunc(ppyex, **kw)

    def ppyey(self, **kw):
        """Calls :py:func:`~warpplots.ppyey` for this species."""
        return self._callppfunc(ppyey, **kw)

    def ppyez(self, **kw):
        """Calls :py:func:`~warpplots.ppyez` for this species."""
        return self._callppfunc(ppyez, **kw)

    def ppzex(self, **kw):
        """Calls :py:func:`~warpplots.ppzex` for this species."""
        return self._callppfunc(ppzex, **kw)

    def ppzey(self, **kw):
        """Calls :py:func:`~warpplots.ppzey` for this species."""
        return self._callppfunc(ppzey, **kw)

    def ppzez(self, **kw):
        """Calls :py:func:`~warpplots.ppzez` for this species."""
        return self._callppfunc(ppzez, **kw)

    def ppxbx(self, **kw):
        """Calls :py:func:`~warpplots.ppxbx` for this species."""
        return self._callppfunc(ppxbx, **kw)

    def ppxby(self, **kw):
        """Calls :py:func:`~warpplots.ppxby` for this species."""
        return self._callppfunc(ppxby, **kw)

    def ppxbz(self, **kw):
        """Calls :py:func:`~warpplots.ppxbz` for this species."""
        return self._callppfunc(ppxbz, **kw)

    def ppybx(self, **kw):
        """Calls :py:func:`~warpplots.ppybx` for this species."""
        return self._callppfunc(ppybx, **kw)

    def ppyby(self, **kw):
        """Calls :py:func:`~warpplots.ppyby` for this species."""
        return self._callppfunc(ppyby, **kw)

    def ppybz(self, **kw):
        """Calls :py:func:`~warpplots.ppybz` for this species."""
        return self._callppfunc(ppybz, **kw)

    def ppzbx(self, **kw):
        """Calls :py:func:`~warpplots.ppzbx` for this species."""
        return self._callppfunc(ppzbx, **kw)

    def ppzby(self, **kw):
        """Calls :py:func:`~warpplots.ppzby` for this species."""
        return self._callppfunc(ppzby, **kw)

    def ppzbz(self, **kw):
        """Calls :py:func:`~warpplots.ppzbz` for this species."""
        return self._callppfunc(ppzbz, **kw)

    def dump(self, filename='pdump.pdb'):
        if self.getn() == 0:
            return
            # x=y=z=ux=uy=uz=gi=pidNone
        else:
            x = self.getx()
            y = self.gety()
            z = self.getz()
            ux = self.getux()
            uy = self.getuy()
            uz = self.getuz()
            gi = self.getgaminv()
            if top.npid > 0:
                pid = self.getpid()
            else:
                pid = None
        if me == 0:
            from ..data_dumping import PWpickle as PW
            f = PW.PW(filename)
            f.time = top.time
            f.x = x
            f.y = y
            f.z = z
            f.ux = ux
            f.uy = uy
            f.uz = uz
            f.gi = gi
            f.pid = pid
            f.close()

    def load(self, filename='pdump.pdb'):
        from ..data_dumping import PRpickle as PR
        f = PR.PR(filename)
        self.addparticles(f.x, f.y, f.z, f.ux, f.uy, f.uz, f.gi, lmomentum=True)
        f.close()

    # --- Fancy python to provide convenient methods for getting various
    # --- species quantities. This allows something like the following, to
    # --- get and set the species weight, for example.
    # --- beam = Species(type=Potassium, charge_state=1)
    # --- beam.sw = npreal/nppart
    # --- print beam.sw
    # --- Note that the documentation can only be obtained by referencing the
    # --- class object. i.e. beam.__class__.sw.__doc__

    # --- Make pgroup a property. This is done so that when pgroup is the same
    # --- as top.pgroup, a reference to top.pgroup is not actually kept in the
    # --- self dict. This avoids problems in a dump/restart.
    def getpgroup(self):
        try:
            return self._pgroup
        except AttributeError:
            if self.pgroups is not None and len(self.pgroups) > 0:
                return self.flatten(self.pgroups)[0]
            else:
                # --- If not otherwise specified, top.pgroup is the default.
                return top.pgroup
    pgroup = property(getpgroup)

    # --- Note that this will only work reliably if there is only one js
    # --- represented.
    def getjs(self):
        return self.jslist[0]
    js = property(getjs)

    def _getpgroupattribute(name, doc=None):
        if doc is None:
            # --- Note that top is used here since self is not actually defined
            # --- at this point (since this is called as a class method).
            doc = top.pgroup.getvardoc(name)
        def fget(self):
            if len(self.jslist) == 1:
                return getattr(self.pgroup, name)[self.jslist[0]]
            else:
                return getattr(self.pgroup, name)[self.jslist]
        def fset(self, value):
            if len(self.jslist) == 1:
                getattr(self.pgroup, name)[self.jslist[0]] = value
            else:
                getattr(self.pgroup, name)[self.jslist] = value
        return fget, fset, None, doc

    sw = property(*_getpgroupattribute('sw', 'Species particle weight, real particles per simulation particle'))
    sm = property(*_getpgroupattribute('sm', 'Species particle mass (in kilograms)'))
    sq = property(*_getpgroupattribute('sq', 'Species particle charge (in coulombs)'))
    ins = property(*_getpgroupattribute('ins', 'Starting index of species'))
    nps = property(*_getpgroupattribute('nps', 'Number of particle in species'))
    sid = property(*_getpgroupattribute('sid'))
    ndts = property(*_getpgroupattribute('ndts'))
    ldts = property(*_getpgroupattribute('ldts'))
    lvdts = property(*_getpgroupattribute('lvdts'))
    iselfb = property(*_getpgroupattribute('iselfb'))
    fselfb = property(*_getpgroupattribute('fselfb'))
    l_maps = property(*_getpgroupattribute('l_maps'))
    dtscale = property(*_getpgroupattribute('dtscale'))
    limplicit = property(*_getpgroupattribute('limplicit'))
    iimplicit = property(*_getpgroupattribute('iimplicit'))
    ldoadvance = property(*_getpgroupattribute('ldoadvance'))
    lparaxial = property(*_getpgroupattribute('lparaxial'))
    zshift = property(*_getpgroupattribute('zshift'))

    # --- This handles the particles from the species.
    # --- Note that this is only for the local domain, no parallel communication
    # --- is done. If particles are needed across parallel domains, use the
    # --- get methods.
    def _getpgroupattribute(name, doc=None):
        def fget(self):
            if len(self.jslist) == 1:
                if self.pgroup.npmax > 0:
                    js = self.jslist[0]
                    i1 = self.pgroup.ins[js] - 1
                    i2 = i1 + self.pgroup.nps[js]
                    return getattr(self.pgroup, name)[i1:i2]
                else:
                    # --- Arrays are unallocated - return a zero length array.
                    # --- This avoids the need of code always having to check if
                    # --- the arrays are allocated.
                    return zeros(0)
            else:
                raise NotImplementedError('The species attributes only works with one species')
        def fset(self, value):
            if len(self.jslist) == 1:
                js = self.jslist[0]
                i1 = self.pgroup.ins[js] - 1
                i2 = i1 + self.pgroup.nps[js]
                getattr(self.pgroup, name)[i1:i2] = value
            else:
                raise NotImplementedError('The species attributes only works with one species')
        return fget, fset, None, doc

    gaminv = property(*_getpgroupattribute('gaminv', 'gamma inverse, in local domain'))
    xp = property(*_getpgroupattribute('xp', 'x position, in local domain'))
    yp = property(*_getpgroupattribute('yp', 'y position, in local domain'))
    zp = property(*_getpgroupattribute('zp', 'z position, in local domain'))
    uxp = property(*_getpgroupattribute('uxp', 'x velocity, in local domain'))
    uyp = property(*_getpgroupattribute('uyp', 'y velocity, in local domain'))
    uzp = property(*_getpgroupattribute('uzp', 'z velocity, in local domain'))
    ex = property(*_getpgroupattribute('ex', 'x E-field, in local domain'))
    ey = property(*_getpgroupattribute('ey', 'y E-field, in local domain'))
    ez = property(*_getpgroupattribute('ez', 'z E-field, in local domain'))
    bx = property(*_getpgroupattribute('bx', 'x B-field, in local domain'))
    by = property(*_getpgroupattribute('by', 'y B-field, in local domain'))
    bz = property(*_getpgroupattribute('bz', 'z B-field, in local domain'))

    # --- This handles the particles pid from the species.
    def _getpgroupattribute(name, doc=None):
        def fget(self):
            if len(self.jslist) == 1:
                if self.pgroup.npmax > 0:
                    js = self.jslist[0]
                    i1 = self.pgroup.ins[js] - 1
                    i2 = i1 + self.pgroup.nps[js]
                    return getattr(self.pgroup, name)[i1:i2,:]
                else:
                    # --- Array is unallocated - return a zero length array.
                    # --- This avoids the need of code always having to check if
                    # --- the arrays are allocated.
                    return zeros((0, self.pgroup.npid))
            else:
                raise NotImplementedError('The species attributes only works with one species')
        def fset(self, value):
            if len(self.jslist) == 1:
                js = self.jslist[0]
                i1 = self.pgroup.ins[js] - 1
                i2 = i1 + self.pgroup.nps[js]
                getattr(self.pgroup, name)[i1:i2,:] = value
            else:
                raise NotImplementedError('The species attributes only works with one species')
        return fget, fset, None, doc

    pid = property(*_getpgroupattribute('pid', 'particle ID information'))

    _builtinpids = {}

    # --- This handles the components of the pid for the particles from the species.
    def _getpgroupattribute(pkg, idname, doc=None, _builtinpids=_builtinpids):
        attrname = idname[:-3]  # strip off the pid suffix
        # --- Add attribute to dictionary, allowing it to be used as a keyword to addparticles
        _builtinpids[attrname] = [pkg, idname]
        def fget(self):
            id = getattr(pkg, idname) - 1
            if id < 0:
                raise Exception('%s not setup'%idname)
            if len(self.jslist) == 1:
                if self.pgroup.npmax > 0:
                    js = self.jslist[0]
                    i1 = self.pgroup.ins[js] - 1
                    i2 = i1 + self.pgroup.nps[js]
                    return self.pgroup.pid[i1:i2,id]
                else:
                    # --- Array is unallocated - return a zero length array.
                    # --- This avoids the need of code always having to check if
                    # --- the arrays are allocated.
                    return zeros((0))
            else:
                raise NotImplementedError('The species attributes only works with one species')
        def fset(self, value):
            id = getattr(pkg, idname) - 1
            if id < 0:
                raise Exception('%s not setup'%idname)
            if len(self.jslist) == 1:
                js = self.jslist[0]
                i1 = self.pgroup.ins[js] - 1
                i2 = i1 + self.pgroup.nps[js]
                self.pgroup.pid[i1:i2,id] = value
            else:
                raise NotImplementedError('The species attributes only works with one species')
        return fget, fset, None, doc

    def _setupgetmethod(pkg, idname):
        def getpidname(self, **kw):
            id = getattr(pkg, idname) - 1
            if id < 0:
                raise Exception('%s not setup'%idname)
            return self.getpid(id=id, **kw)
        return getpidname

    w = property(*_getpgroupattribute(top, 'wpid', 'particle variable weight'))
    # getw is already defined
    injdata = property(*_getpgroupattribute(top, 'injdatapid', 'particle injection data'))
    getinjdata = _setupgetmethod(top, 'injdatapid')
    ssn = property(*_getpgroupattribute(top, 'ssnpid', 'particle SSNs'))
    getssn = _setupgetmethod(top, 'ssnpid')
    ssnparent = property(*_getpgroupattribute(top, 'ssnparentpid', 'particle parent SSNs'))
    getssnparent = _setupgetmethod(top, 'ssnparentpid')
    rbirth = property(*_getpgroupattribute(top, 'rbirthpid', 'particle initial radius'))
    getrbirth = _setupgetmethod(top, 'rbirthpid')
    tbirth = property(*_getpgroupattribute(top, 'tbirthpid', 'particle creation time'))
    gettbirth = _setupgetmethod(top, 'tbirthpid')
    xbirth = property(*_getpgroupattribute(top, 'xbirthpid', 'particle birth x'))
    getxbirth = _setupgetmethod(top, 'xbirthpid')
    ybirth = property(*_getpgroupattribute(top, 'ybirthpid', 'particle birth y'))
    getybirth = _setupgetmethod(top, 'ybirthpid')
    zbirth = property(*_getpgroupattribute(top, 'zbirthpid', 'particle birth z'))
    getzbirth = _setupgetmethod(top, 'zbirthpid')
    zbirthlab = property(*_getpgroupattribute(top, 'zbirthlabpid', 'particle birth z in lab frame'))
    getzbirthlab = _setupgetmethod(top, 'zbirthlabpid')
    uxbirth = property(*_getpgroupattribute(top, 'uxbirthpid', 'particle birth ux'))
    getuxbirth = _setupgetmethod(top, 'uxbirthpid')
    uybirth = property(*_getpgroupattribute(top, 'uybirthpid', 'particle birth uy'))
    getuybirth = _setupgetmethod(top, 'uybirthpid')
    uzbirth = property(*_getpgroupattribute(top, 'uzbirthpid', 'particle birth uz'))
    getuzbirth = _setupgetmethod(top, 'uzbirthpid')
    xold = property(*_getpgroupattribute(top, 'xoldpid', 'particle previous position in x'))
    getxold = _setupgetmethod(top, 'xoldpid')
    yold = property(*_getpgroupattribute(top, 'yoldpid', 'particle previous position in y'))
    getyold = _setupgetmethod(top, 'yoldpid')
    zold = property(*_getpgroupattribute(top, 'zoldpid', 'particle previous position in z'))
    getzold = _setupgetmethod(top, 'zoldpid')
    uxold = property(*_getpgroupattribute(top, 'uxoldpid', 'particle previous x velocity'))
    getuxold = _setupgetmethod(top, 'uxoldpid')
    uyold = property(*_getpgroupattribute(top, 'uyoldpid', 'particle previous y velocity'))
    getuyold = _setupgetmethod(top, 'uyoldpid')
    uzold = property(*_getpgroupattribute(top, 'uzoldpid', 'particle previous z velocity'))
    getuzold = _setupgetmethod(top, 'uzoldpid')
    vdxold = property(*_getpgroupattribute(top, 'vdxoldpid', 'particle previous x drift velocity'))
    getvdxold = _setupgetmethod(top, 'vdxoldpid')
    vdyold = property(*_getpgroupattribute(top, 'vdyoldpid', 'particle previous y drift velocity'))
    getvdyold = _setupgetmethod(top, 'vdyoldpid')
    vdzold = property(*_getpgroupattribute(top, 'vdzoldpid', 'particle previous z drift velocity'))
    getvdzold = _setupgetmethod(top, 'vdzoldpid')
    dx = property(*_getpgroupattribute(top, 'dxpid', 'dx grid cell size at particle location'))
    getdx = _setupgetmethod(top, 'dxpid')
    dy = property(*_getpgroupattribute(top, 'dypid', 'dy grid cell size at particle location'))
    getdy = _setupgetmethod(top, 'dypid')
    dz = property(*_getpgroupattribute(top, 'dzpid', 'dz grid cell size at particle location'))
    getdz = _setupgetmethod(top, 'dzpid')
    bxold = property(*_getpgroupattribute(top, 'bxoldpid', 'particle old Bx field'))
    getbxold = _setupgetmethod(top, 'bxoldpid')
    byold = property(*_getpgroupattribute(top, 'byoldpid', 'particle old By field'))
    getbyold = _setupgetmethod(top, 'byoldpid')
    bzold = property(*_getpgroupattribute(top, 'bzoldpid', 'particle old Bz field'))
    getbzold = _setupgetmethod(top, 'bzoldpid')
    exold = property(*_getpgroupattribute(top, 'exoldpid', 'particle old Ex field'))
    getexold = _setupgetmethod(top, 'exoldpid')
    eyold = property(*_getpgroupattribute(top, 'eyoldpid', 'particle old Ey field'))
    geteyold = _setupgetmethod(top, 'eyoldpid')
    ezold = property(*_getpgroupattribute(top, 'ezoldpid', 'particle old Ez field'))
    getezold = _setupgetmethod(top, 'ezoldpid')
    chdts = property(*_getpgroupattribute(top, 'chdtspid', 'particle flag for dts change'))
    getchdts = _setupgetmethod(top, 'chdtspid')
    uparBo = property(*_getpgroupattribute(top, 'uparBopid', 'particle uparallel * B'))
    getuparBo = _setupgetmethod(top, 'uparBopid')
    bxpred = property(*_getpgroupattribute(top, 'bxpredpid', 'particle predicted Bx field'))
    getbxpred = _setupgetmethod(top, 'bxpredpid')
    bypred = property(*_getpgroupattribute(top, 'bypredpid', 'particle predicted By field'))
    getbypred = _setupgetmethod(top, 'bypredpid')
    bzpred = property(*_getpgroupattribute(top, 'bzpredpid', 'particle predicted Bz field'))
    getbzpred = _setupgetmethod(top, 'bzpredpid')
    uparoBpred = property(*_getpgroupattribute(top, 'uparoBpredpid', 'particle uparallel * B'))
    getuparoBpred = _setupgetmethod(top, 'uparoBpredpid')
    expred = property(*_getpgroupattribute(top, 'expredpid', 'particle predicted Ex field'))
    getexpred = _setupgetmethod(top, 'expredpid')
    eypred = property(*_getpgroupattribute(top, 'eypredpid', 'particle predicted Ey field'))
    geteypred = _setupgetmethod(top, 'eypredpid')
    ezpred = property(*_getpgroupattribute(top, 'ezpredpid', 'particle predicted Ez field'))
    getezpred = _setupgetmethod(top, 'ezpredpid')
    lmap = property(*_getpgroupattribute(top, 'lmappid', 'particle location in lmap lattice'))
    getlmap = _setupgetmethod(top, 'lmappid')
    dt = property(*_getpgroupattribute(wxy, 'dtpid', 'particle time step size for slice code'))
    getdt = _setupgetmethod(wxy, 'dtpid')

    def _gettopattribute(name, doc=None):
        if doc is None:
            doc = top.getvardoc(name)
        def fget(self):
            if len(self.jslist) == 1:
                return getattr(top, name)[self.jslist[0]]
            else:
                return getattr(top, name)[self.jslist]
        def fset(self, value):
            if len(self.jslist) == 1:
                getattr(top, name)[self.jslist[0]] = value
            else:
                getattr(top, name)[self.jslist] = value
        return fget, fset, None, doc

    # InPart
    np = property(*_gettopattribute('np_s'))
    a0 = property(*_gettopattribute('a0_s'))
    ap0 = property(*_gettopattribute('ap0_s'))
    b0 = property(*_gettopattribute('b0_s'))
    bp0 = property(*_gettopattribute('bp0_s'))
    x0 = property(*_gettopattribute('x0_s'))
    xp0 = property(*_gettopattribute('xp0_s'))
    y0 = property(*_gettopattribute('y0_s'))
    yp0 = property(*_gettopattribute('yp0_s'))
    aion = property(*_gettopattribute('aion_s'))
    ekin = property(*_gettopattribute('ekin_s'))
    emit = property(*_gettopattribute('emit_s'))
    emitx = property(*_gettopattribute('emitx_s'))
    emity = property(*_gettopattribute('emity_s'))
    emitn = property(*_gettopattribute('emitn_s'))
    emitnx = property(*_gettopattribute('emitnx_s'))
    emitny = property(*_gettopattribute('emitny_s'))
    ibeam = property(*_gettopattribute('ibeam_s'))
    zion = property(*_gettopattribute('zion_s'))
    straight = property(*_gettopattribute('straight_s'))
    vbeam = property(*_gettopattribute('vbeam_s'))
    vtilt = property(*_gettopattribute('vtilt_s'))
    vthperp = property(*_gettopattribute('vthperp_s'))
    vthz = property(*_gettopattribute('vthz_s'))
    zimin = property(*_gettopattribute('zimin_s'))
    zimax = property(*_gettopattribute('zimax_s'))
    sp_fract = property(*_gettopattribute('sp_fract'))
    xcent = property(*_gettopattribute('xcent_s'))
    ycent = property(*_gettopattribute('ycent_s'))
    xpcent = property(*_gettopattribute('xpcent_s'))
    ypcent = property(*_gettopattribute('ypcent_s'))
    efetch = property(*_gettopattribute('efetch'))

    # InjectVars
    npinject = property(*_gettopattribute('npinje_s'))
    rnpinject = property(*_gettopattribute('rnpinje_s'))
    ntinject = property(*_gettopattribute('ntinject'))

    # SelfB
    iselfb = property(*_gettopattribute('iselfb'))

    # LostParticles
    inslost = property(*_gettopattribute('inslost'))
    npslost = property(*_gettopattribute('npslost'))

    # --- This handles the lost particles from the species.
    # --- Note that this is on for the local domain, no parallel communication
    # --- is done. If lost particles are needed across parallel domains, use the
    # --- get methods.
    def _getpgroupattribute(name, doc=None):
        def fget(self):
            if len(self.jslist) == 1:
                if top.npmaxlost > 0:
                    js = self.jslist[0]
                    i1 = top.inslost[js] - 1
                    i2 = i1 + top.npslost[js]
                    return getattr(top, name)[i1:i2]
                else:
                    # --- Arrays are unallocated - return a zero length array.
                    # --- This avoids the need of code always having to check if
                    # --- the arrays are allocated.
                    return zeros(0)
            else:
                raise NotImplementedError('The species attributes only works with one species')
        def fset(self, value):
            if len(self.jslist) == 1:
                js = self.jslist[0]
                i1 = top.inslost[js] - 1
                i2 = i1 + top.npslost[js]
                getattr(top, name)[i1:i2] = value
            else:
                raise NotImplementedError('The species attributes only works with one species')
        return fget, fset, None, doc

    gaminvlost = property(*_getpgroupattribute('gaminvlost', 'gamma inverse of lost particles, in local domain'))
    xplost = property(*_getpgroupattribute('xplost', 'x position of lost particles, in local domain'))
    yplost = property(*_getpgroupattribute('yplost', 'y position of lost particles, in local domain'))
    zplost = property(*_getpgroupattribute('zplost', 'z position of lost particles, in local domain'))
    uxplost = property(*_getpgroupattribute('uxplost', 'x velocity of lost particles, in local domain'))
    uyplost = property(*_getpgroupattribute('uyplost', 'y velocity of lost particles, in local domain'))
    uzplost = property(*_getpgroupattribute('uzplost', 'z velocity of lost particles, in local domain'))
    exlost = property(*_getpgroupattribute('exlost', 'x E-field of lost particles, in local domain'))
    eylost = property(*_getpgroupattribute('eylost', 'y E-field of lost particles, in local domain'))
    ezlost = property(*_getpgroupattribute('ezlost', 'z E-field of lost particles, in local domain'))
    bxlost = property(*_getpgroupattribute('bxlost', 'x B-field of lost particles, in local domain'))
    bylost = property(*_getpgroupattribute('bylost', 'y B-field of lost particles, in local domain'))
    bzlost = property(*_getpgroupattribute('bzlost', 'z B-field of lost particles, in local domain'))
    tplost = property(*_getpgroupattribute('tplost', 'time particles were lost, in local domain'))

    def _getpgroupattribute(name, doc=None):
        def fget(self):
            if len(self.jslist) == 1:
                if top.npmaxlost > 0:
                    js = self.jslist[0]
                    i1 = top.inslost[js] - 1
                    i2 = i1 + top.npslost[js]
                    return getattr(top, name)[i1:i2,:]
                else:
                    # --- Arrays are unallocated - return a zero length array.
                    # --- This avoids the need of code always having to check if
                    # --- the arrays are allocated.
                    return zeros((0, top.npidlost))
            else:
                raise NotImplementedError('The species attributes only works with one species')
        def fset(self, value):
            if len(self.jslist) == 1:
                js = self.jslist[0]
                i1 = top.inslost[js] - 1
                i2 = i1 + top.npslost[js]
                getattr(top, name)[i1:i2,:] = value
            else:
                raise NotImplementedError('The species attributes only works with one species')
        return fget, fset, None, doc

    pidlost = property(*_getpgroupattribute('pidlost', 'particle ID information of lost particles, in local domain'))

    # --- Clean up the name space
    del _getpgroupattribute

    # --- For variables that are 2, 3, or 4-D
    def _gettopattribute(name, doc=None):
        if doc is None:
            doc = top.getvardoc(name)
        def fget(self):
            if len(self.jslist) == 1:
                return getattr(top, name)[..., self.jslist[0]]
            else:
                return getattr(top, name)[..., self.jslist]
        def fset(self, value):
            if len(self.jslist) == 1:
                getattr(top, name)[..., self.jslist[0]] = value
            else:
                getattr(top, name)[..., self.jslist] = value
        return fget, fset, None, doc

    # InPart
    depos_order = property(*_gettopattribute('depos_order'))

    # InjectVars
    vzinject = property(*_gettopattribute('vzinject'))
    finject = property(*_gettopattribute('finject'))
    winject = property(*_gettopattribute('winject'))
    npinjtmp = property(*_gettopattribute('npinjtmp'))
    ftinject = property(*_gettopattribute('ftinject'))
    wtinject = property(*_gettopattribute('wtinject'))
    tinj_npactual = property(*_gettopattribute('tinj_npactual'))
    tinjprev = property(*_gettopattribute('tinjprev'))

    # ExtPart
    nep = property(*_gettopattribute('nep'))
    tep = property(*_gettopattribute('tep'))
    xep = property(*_gettopattribute('xep'))
    yep = property(*_gettopattribute('yep'))
    uxep = property(*_gettopattribute('uxep'))
    uyep = property(*_gettopattribute('uyep'))
    uzep = property(*_gettopattribute('uzep'))
    pidep = property(*_gettopattribute('pidep'))

    # SemiTransparentDisc
    s_STdiscs = property(*_gettopattribute('s_STdiscs'))

    del _gettopattribute

    def _getw3dattribute(name, doc=None):
        if doc is None:
            doc = w3d.getvardoc(name)
        def fget(self):
            if len(self.jslist) == 1:
                return getattr(w3d, name)[self.jslist[0]]
            else:
                return getattr(w3d, name)[self.jslist]
        def fset(self, value):
            if len(self.jslist) == 1:
                getattr(w3d, name)[self.jslist[0]] = value
            else:
                getattr(w3d, name)[self.jslist[0]] = value
        return fget, fset, None, doc

    # w3d.DKInterp
    m_over_q = property(*_getw3dattribute('m_over_q'))
    qovermsq = property(*_getw3dattribute('qovermsq'))
    alpha0 = property(*_getw3dattribute('alpha0'))
    acntr = property(*_getw3dattribute('acntr'))
    usealphacalc = property(*_getw3dattribute('usealphacalc'))
    notusealphcalc = property(*_getw3dattribute('notusealphcalc'))
    interpdk = property(*_getw3dattribute('interpdk'))

    del _getw3dattribute

    # --- For variables that are 1-D, and use sid
    # --- Note that this won't handle multiple js values well. It assumes
    # --- that all of the jslist have the same jsid.
    def _gettopjsid1dattribute(name, doc=None):
        if doc is None:
            doc = top.getvardoc(name)
        def fget(self):
            jsid = self.pgroup.sid[self.jslist[0]]
            return getattr(top, name)[jsid]
        def fset(self, value):
            jsid = self.pgroup.sid[self.jslist[0]]
            getattr(top, name)[jsid] = value
        return fget, fset, None, doc

    # --- For variables that are 2, 3, or 4-D, and use sid
    # --- Note that this won't handle multiple js values well. It assumes
    # --- that all of the jslist have the same jsid.
    def _gettopjsidattribute(name, doc=None):
        if doc is None:
            doc = top.getvardoc(name)
        def fget(self):
            jsid = self.pgroup.sid[self.jslist[0]]
            return getattr(top, name)[..., jsid]
        def fset(self, value):
            jsid = self.pgroup.sid[self.jslist[0]]
            getattr(top, name)[..., jsid] = value
        return fget, fset, None, doc

    # Z_arrays
    curr = property(*_gettopjsidattribute('curr'))
    lostpars = property(*_gettopjsidattribute('lostpars'))

    # Win_Moments
    npsim = property(*_gettopjsidattribute('npsim'))
    pnum = property(*_gettopjsidattribute('pnum'))
    xbar = property(*_gettopjsidattribute('xbar'))
    ybar = property(*_gettopjsidattribute('ybar'))
    zbar = property(*_gettopjsidattribute('zbar'))
    xpbar = property(*_gettopjsidattribute('xpbar'))
    ypbar = property(*_gettopjsidattribute('ypbar'))
    vxbar = property(*_gettopjsidattribute('vxbar'))
    vybar = property(*_gettopjsidattribute('vybar'))
    vzbar = property(*_gettopjsidattribute('vzbar'))
    xybar = property(*_gettopjsidattribute('xybar'))
    xypbar = property(*_gettopjsidattribute('xypbar'))
    yxpbar = property(*_gettopjsidattribute('yxpbar'))
    xpypbar = property(*_gettopjsidattribute('xpypbar'))
    xsqbar = property(*_gettopjsidattribute('xsqbar'))
    ysqbar = property(*_gettopjsidattribute('ysqbar'))
    zsqbar = property(*_gettopjsidattribute('zsqbar'))
    xpsqbar = property(*_gettopjsidattribute('xpsqbar'))
    ypsqbar = property(*_gettopjsidattribute('ypsqbar'))
    vxsqbar = property(*_gettopjsidattribute('vxsqbar'))
    vysqbar = property(*_gettopjsidattribute('vysqbar'))
    vzsqbar = property(*_gettopjsidattribute('vzsqbar'))
    xxpbar = property(*_gettopjsidattribute('xxpbar'))
    yypbar = property(*_gettopjsidattribute('yypbar'))
    zvzbar = property(*_gettopjsidattribute('zvzbar'))
    xvzbar = property(*_gettopjsidattribute('xvzbar'))
    yvzbar = property(*_gettopjsidattribute('yvzbar'))
    vxvzbar = property(*_gettopjsidattribute('vxvzbar'))
    vyvzbar = property(*_gettopjsidattribute('vyvzbar'))
    xrms = property(*_gettopjsidattribute('xrms'))
    yrms = property(*_gettopjsidattribute('yrms'))
    zrms = property(*_gettopjsidattribute('zrms'))
    rrms = property(*_gettopjsidattribute('rrms'))
    xprms = property(*_gettopjsidattribute('xprms'))
    yprms = property(*_gettopjsidattribute('yprms'))
    epsx = property(*_gettopjsidattribute('epsx'))
    epsy = property(*_gettopjsidattribute('epsy'))
    epsz = property(*_gettopjsidattribute('epsz'))
    epsnx = property(*_gettopjsidattribute('epsnx'))
    epsny = property(*_gettopjsidattribute('epsny'))
    epsnz = property(*_gettopjsidattribute('epsnz'))
    epsr = property(*_gettopjsidattribute('epsr'))
    epsg = property(*_gettopjsidattribute('epsg'))
    epsh = property(*_gettopjsidattribute('epsh'))
    epsnr = property(*_gettopjsidattribute('epsnr'))
    epsng = property(*_gettopjsidattribute('epsng'))
    epsnh = property(*_gettopjsidattribute('epsnh'))
    vxrms = property(*_gettopjsidattribute('vxrms'))
    vyrms = property(*_gettopjsidattribute('vyrms'))
    vzrms = property(*_gettopjsidattribute('vzrms'))

    # Z_Moments
    zmmntsq = property(*_gettopjsid1dattribute('zmmntsq'))
    zmmntsm = property(*_gettopjsid1dattribute('zmmntsm'))
    zmmntsw = property(*_gettopjsid1dattribute('zmmntsw'))
    zmomentscalculated = property(*_gettopjsid1dattribute('zmomentscalculated'))
    npsimz = property(*_gettopjsidattribute('npsimz'))
    pnumz = property(*_gettopjsidattribute('pnumz'))
    xbarz = property(*_gettopjsidattribute('xbarz'))
    ybarz = property(*_gettopjsidattribute('ybarz'))
    zbarz = property(*_gettopjsidattribute('zbarz'))
    xpbarz = property(*_gettopjsidattribute('xpbarz'))
    ypbarz = property(*_gettopjsidattribute('ypbarz'))
    vxbarz = property(*_gettopjsidattribute('vxbarz'))
    vybarz = property(*_gettopjsidattribute('vybarz'))
    vzbarz = property(*_gettopjsidattribute('vzbarz'))
    xybarz = property(*_gettopjsidattribute('xybarz'))
    xypbarz = property(*_gettopjsidattribute('xypbarz'))
    yxpbarz = property(*_gettopjsidattribute('yxpbarz'))
    xpypbarz = property(*_gettopjsidattribute('xpypbarz'))
    xsqbarz = property(*_gettopjsidattribute('xsqbarz'))
    ysqbarz = property(*_gettopjsidattribute('ysqbarz'))
    zsqbarz = property(*_gettopjsidattribute('zsqbarz'))
    xpsqbarz = property(*_gettopjsidattribute('xpsqbarz'))
    ypsqbarz = property(*_gettopjsidattribute('ypsqbarz'))
    vxsqbarz = property(*_gettopjsidattribute('vxsqbarz'))
    vysqbarz = property(*_gettopjsidattribute('vysqbarz'))
    vzsqbarz = property(*_gettopjsidattribute('vzsqbarz'))
    xxpbarz = property(*_gettopjsidattribute('xxpbarz'))
    yypbarz = property(*_gettopjsidattribute('yypbarz'))
    zvzbarz = property(*_gettopjsidattribute('zvzbarz'))
    xvzbarz = property(*_gettopjsidattribute('xvzbarz'))
    yvzbarz = property(*_gettopjsidattribute('yvzbarz'))
    vxvzbarz = property(*_gettopjsidattribute('vxvzbarz'))
    vyvzbarz = property(*_gettopjsidattribute('vyvzbarz'))
    xrmsz = property(*_gettopjsidattribute('xrmsz'))
    yrmsz = property(*_gettopjsidattribute('yrmsz'))
    zrmsz = property(*_gettopjsidattribute('zrmsz'))
    rrmsz = property(*_gettopjsidattribute('rrmsz'))
    xprmsz = property(*_gettopjsidattribute('xprmsz'))
    yprmsz = property(*_gettopjsidattribute('yprmsz'))
    epsxz = property(*_gettopjsidattribute('epsxz'))
    epsyz = property(*_gettopjsidattribute('epsyz'))
    epszz = property(*_gettopjsidattribute('epszz'))
    epsnxz = property(*_gettopjsidattribute('epsnxz'))
    epsnyz = property(*_gettopjsidattribute('epsnyz'))
    epsnzz = property(*_gettopjsidattribute('epsnzz'))
    epsrz = property(*_gettopjsidattribute('epsrz'))
    epsgz = property(*_gettopjsidattribute('epsgz'))
    epshz = property(*_gettopjsidattribute('epshz'))
    epsnrz = property(*_gettopjsidattribute('epsnrz'))
    epsngz = property(*_gettopjsidattribute('epsngz'))
    epsnhz = property(*_gettopjsidattribute('epsnhz'))
    vxrmsz = property(*_gettopjsidattribute('vxrmsz'))
    vyrmsz = property(*_gettopjsidattribute('vyrmsz'))
    vzrmsz = property(*_gettopjsidattribute('vzrmsz'))
    tempmaxp = property(*_gettopjsidattribute('tempmaxp'))
    tempminp = property(*_gettopjsidattribute('tempminp'))
    tempzmmnts0 = property(*_gettopjsidattribute('tempzmmnts0'))
    tempzmmnts = property(*_gettopjsidattribute('tempzmmnts'))

    # --- For lab-window moments.
    # --- Note that this won't handle multiple js values well. It assumes
    # --- that all of the jslist have the same jsid.
    def _gettoplabmomattribute(name, doc=None):
        if doc is None:
            doc = top.getvardoc(name)
        def fget(self):
            # --- A function is needed that takes ilw as an argument
            def _getlabwindowmomentdata(ilw):
                jsid = self.pgroup.sid[self.jslist[0]]
                return getattr(top, name)[:top.ilabwn[ilw, jsid], ilw, jsid]
            return _getlabwindowmomentdata
        def fset(self, value):
            jsid = self.pgroup.sid[self.jslist[0]]
            getattr(top, name)[..., jsid] = value
        return fget, fset, None, doc

    # Lab_Moments
    ilabwn = property(*_gettopjsidattribute('ilabwn'))
    timelw = property(*_gettoplabmomattribute('timelw'))
    npsimlw = property(*_gettoplabmomattribute('npsimlw'))
    pnumlw = property(*_gettoplabmomattribute('pnumlw'))
    xbarlw = property(*_gettoplabmomattribute('xbarlw'))
    ybarlw = property(*_gettoplabmomattribute('ybarlw'))
    vzbarlw = property(*_gettoplabmomattribute('vzbarlw'))
    epsxlw = property(*_gettoplabmomattribute('epsxlw'))
    epsylw = property(*_gettoplabmomattribute('epsylw'))
    epszlw = property(*_gettoplabmomattribute('epszlw'))
    vxrmslw = property(*_gettoplabmomattribute('vxrmslw'))
    vyrmslw = property(*_gettoplabmomattribute('vyrmslw'))
    vzrmslw = property(*_gettoplabmomattribute('vzrmslw'))
    xrmslw = property(*_gettoplabmomattribute('xrmslw'))
    yrmslw = property(*_gettoplabmomattribute('yrmslw'))
    rrmslw = property(*_gettoplabmomattribute('rrmslw'))
    xxpbarlw = property(*_gettoplabmomattribute('xxpbarlw'))
    yypbarlw = property(*_gettoplabmomattribute('yypbarlw'))
    currlw = property(*_gettoplabmomattribute('currlw'))
    lostparslw = property(*_gettoplabmomattribute('lostparslw'))
    linechglw = property(*_gettoplabmomattribute('linechglw'))

    # Moments
    ek = property(*_gettopjsid1dattribute('ek'))
    ekzmbe = property(*_gettopjsid1dattribute('ekzmbe'))
    ekzbeam = property(*_gettopjsid1dattribute('ekzbeam'))
    ekperp = property(*_gettopjsid1dattribute('ekperp'))
    pz = property(*_gettopjsid1dattribute('pz'))
    xmaxp = property(*_gettopjsid1dattribute('xmaxp'))
    xminp = property(*_gettopjsid1dattribute('xminp'))
    ymaxp = property(*_gettopjsid1dattribute('ymaxp'))
    yminp = property(*_gettopjsid1dattribute('yminp'))
    zmaxp = property(*_gettopjsid1dattribute('zmaxp'))
    zminp = property(*_gettopjsid1dattribute('zminp'))
    vxmaxp = property(*_gettopjsid1dattribute('vxmaxp'))
    vxminp = property(*_gettopjsid1dattribute('vxminp'))
    vymaxp = property(*_gettopjsid1dattribute('vymaxp'))
    vyminp = property(*_gettopjsid1dattribute('vyminp'))
    vzmaxp = property(*_gettopjsid1dattribute('vzmaxp'))
    vzminp = property(*_gettopjsid1dattribute('vzminp'))

    # --- For variables that are 2, 3, or 4-D, and use sid
    # --- Note that this won't handle multiple js values well. It assumes
    # --- that all of the jslist have the same jsid.
    def _gettopjsidhistattribute(name, doc=None):
        if doc is None:
            doc = top.getvardoc(name)
        def fget(self):
            jsid = self.pgroup.sid[self.jslist[0]]
            return getattr(top, name)[..., :top.jhist+1, jsid]
        def fset(self, value):
            jsid = self.pgroup.sid[self.jslist[0]]
            getattr(top, name)[..., :top.jhist+1, jsid] = value
        return fget, fset, None, doc

    # Hist
    hbmlen = property(*_gettopjsidhistattribute('hbmlen'))
    hekzmbe = property(*_gettopjsidhistattribute('hekzmbe'))
    hekzbeam = property(*_gettopjsidhistattribute('hekzbeam'))
    hekperp = property(*_gettopjsidhistattribute('hekperp'))
    hxmaxp = property(*_gettopjsidhistattribute('hxmaxp'))
    hxminp = property(*_gettopjsidhistattribute('hxminp'))
    hymaxp = property(*_gettopjsidhistattribute('hymaxp'))
    hyminp = property(*_gettopjsidhistattribute('hyminp'))
    hzmaxp = property(*_gettopjsidhistattribute('hzmaxp'))
    hzminp = property(*_gettopjsidhistattribute('hzminp'))
    hvxmaxp = property(*_gettopjsidhistattribute('hvxmaxp'))
    hvxminp = property(*_gettopjsidhistattribute('hvxminp'))
    hvymaxp = property(*_gettopjsidhistattribute('hvymaxp'))
    hvyminp = property(*_gettopjsidhistattribute('hvyminp'))
    hvzmaxp = property(*_gettopjsidhistattribute('hvzmaxp'))
    hvzminp = property(*_gettopjsidhistattribute('hvzminp'))
    hepsx = property(*_gettopjsidhistattribute('hepsx'))
    hepsy = property(*_gettopjsidhistattribute('hepsy'))
    hepsz = property(*_gettopjsidhistattribute('hepsz'))
    hepsnx = property(*_gettopjsidhistattribute('hepsnx'))
    hepsny = property(*_gettopjsidhistattribute('hepsny'))
    hepsnz = property(*_gettopjsidhistattribute('hepsnz'))
    hepsr = property(*_gettopjsidhistattribute('hepsr'))
    hepsg = property(*_gettopjsidhistattribute('hepsg'))
    hepsh = property(*_gettopjsidhistattribute('hepsh'))
    hepsnr = property(*_gettopjsidhistattribute('hepsnr'))
    hepsng = property(*_gettopjsidhistattribute('hepsng'))
    hepsnh = property(*_gettopjsidhistattribute('hepsnh'))
    hnpsim = property(*_gettopjsidhistattribute('hnpsim'))
    hpnum = property(*_gettopjsidhistattribute('hpnum'))
    hxbar = property(*_gettopjsidhistattribute('hxbar'))
    hybar = property(*_gettopjsidhistattribute('hybar'))
    hzbar = property(*_gettopjsidhistattribute('hzbar'))
    hxybar = property(*_gettopjsidhistattribute('hxybar'))
    hxrms = property(*_gettopjsidhistattribute('hxrms'))
    hyrms = property(*_gettopjsidhistattribute('hyrms'))
    hrrms = property(*_gettopjsidhistattribute('hrrms'))
    hzrms = property(*_gettopjsidhistattribute('hzrms'))
    hxprms = property(*_gettopjsidhistattribute('hxprms'))
    hyprms = property(*_gettopjsidhistattribute('hyprms'))
    hxsqbar = property(*_gettopjsidhistattribute('hxsqbar'))
    hysqbar = property(*_gettopjsidhistattribute('hysqbar'))
    hvxbar = property(*_gettopjsidhistattribute('hvxbar'))
    hvybar = property(*_gettopjsidhistattribute('hvybar'))
    hvzbar = property(*_gettopjsidhistattribute('hvzbar'))
    hxpbar = property(*_gettopjsidhistattribute('hxpbar'))
    hypbar = property(*_gettopjsidhistattribute('hypbar'))
    hvxrms = property(*_gettopjsidhistattribute('hvxrms'))
    hvyrms = property(*_gettopjsidhistattribute('hvyrms'))
    hvzrms = property(*_gettopjsidhistattribute('hvzrms'))
    hxpsqbar = property(*_gettopjsidhistattribute('hxpsqbar'))
    hypsqbar = property(*_gettopjsidhistattribute('hypsqbar'))
    hxxpbar = property(*_gettopjsidhistattribute('hxxpbar'))
    hyypbar = property(*_gettopjsidhistattribute('hyypbar'))
    hxypbar = property(*_gettopjsidhistattribute('hxypbar'))
    hyxpbar = property(*_gettopjsidhistattribute('hyxpbar'))
    hxpypbar = property(*_gettopjsidhistattribute('hxpypbar'))
    hxvzbar = property(*_gettopjsidhistattribute('hxvzbar'))
    hyvzbar = property(*_gettopjsidhistattribute('hyvzbar'))
    hvxvzbar = property(*_gettopjsidhistattribute('hvxvzbar'))
    hvyvzbar = property(*_gettopjsidhistattribute('hvyvzbar'))
    hcurrz = property(*_gettopjsidhistattribute('hcurrz'))
    hnpsimz = property(*_gettopjsidhistattribute('hnpsimz'))
    hpnumz = property(*_gettopjsidhistattribute('hpnumz'))
    hepsxz = property(*_gettopjsidhistattribute('hepsxz'))
    hepsyz = property(*_gettopjsidhistattribute('hepsyz'))
    hepsnxz = property(*_gettopjsidhistattribute('hepsnxz'))
    hepsnyz = property(*_gettopjsidhistattribute('hepsnyz'))
    hepsrz = property(*_gettopjsidhistattribute('hepsrz'))
    hepsgz = property(*_gettopjsidhistattribute('hepsgz'))
    hepshz = property(*_gettopjsidhistattribute('hepshz'))
    hepsnrz = property(*_gettopjsidhistattribute('hepsnrz'))
    hepsngz = property(*_gettopjsidhistattribute('hepsngz'))
    hepsnhz = property(*_gettopjsidhistattribute('hepsnhz'))
    hxbarz = property(*_gettopjsidhistattribute('hxbarz'))
    hybarz = property(*_gettopjsidhistattribute('hybarz'))
    hxybarz = property(*_gettopjsidhistattribute('hxybarz'))
    hxrmsz = property(*_gettopjsidhistattribute('hxrmsz'))
    hyrmsz = property(*_gettopjsidhistattribute('hyrmsz'))
    hrrmsz = property(*_gettopjsidhistattribute('hrrmsz'))
    hxprmsz = property(*_gettopjsidhistattribute('hxprmsz'))
    hyprmsz = property(*_gettopjsidhistattribute('hyprmsz'))
    hxsqbarz = property(*_gettopjsidhistattribute('hxsqbarz'))
    hysqbarz = property(*_gettopjsidhistattribute('hysqbarz'))
    hvxbarz = property(*_gettopjsidhistattribute('hvxbarz'))
    hvybarz = property(*_gettopjsidhistattribute('hvybarz'))
    hvzbarz = property(*_gettopjsidhistattribute('hvzbarz'))
    hxpbarz = property(*_gettopjsidhistattribute('hxpbarz'))
    hypbarz = property(*_gettopjsidhistattribute('hypbarz'))
    hvxrmsz = property(*_gettopjsidhistattribute('hvxrmsz'))
    hvyrmsz = property(*_gettopjsidhistattribute('hvyrmsz'))
    hvzrmsz = property(*_gettopjsidhistattribute('hvzrmsz'))
    hxpsqbarz = property(*_gettopjsidhistattribute('hxpsqbarz'))
    hypsqbarz = property(*_gettopjsidhistattribute('hypsqbarz'))
    hxxpbarz = property(*_gettopjsidhistattribute('hxxpbarz'))
    hyypbarz = property(*_gettopjsidhistattribute('hyypbarz'))
    hxypbarz = property(*_gettopjsidhistattribute('hxypbarz'))
    hyxpbarz = property(*_gettopjsidhistattribute('hyxpbarz'))
    hxpypbarz = property(*_gettopjsidhistattribute('hxpypbarz'))
    hxvzbarz = property(*_gettopjsidhistattribute('hxvzbarz'))
    hyvzbarz = property(*_gettopjsidhistattribute('hyvzbarz'))
    hvxvzbarz = property(*_gettopjsidhistattribute('hvxvzbarz'))
    hvyvzbarz = property(*_gettopjsidhistattribute('hvyvzbarz'))

    del _gettopjsid1dattribute
    del _gettopjsidattribute
    del _gettoplabmomattribute
    del _gettopjsidhistattribute

    # --- Save a list of added pid attributes
    _addedpids = []

    @classmethod
    def addpid(cls, name, pid=None):
        """Dynamically add a new particle attribute
        - name: name of new attribute. An attribute name+'pid' will also be added to hold the pid number.
                The particle data can be access directly with the attribute or via a get method.
                For example, by calling beam.addpid('foo'), the data can be accessed using beam.foo or beam.getfoo().
        - pid=None: if specified, the pid to use, otherwise nextpid will be called
        """
        if name in cls.__dict__ or name + 'pid' in cls.__dict__:
            raise Exception('A new pid must have a unique name')
        if pid is None:
            pid = nextpid()

        def fget(self):
            if len(self.jslist) == 1:
                if self.pgroup.npmax > 0:
                    js = self.jslist[0]
                    i1 = self.pgroup.ins[js] - 1
                    i2 = i1 + self.pgroup.nps[js]
                    return self.pgroup.pid[i1:i2,pid-1]
                else:
                    # --- Arrays are unallocated - return a zero length array.
                    # --- This avoids the need of code always having to check if
                    # --- the arrays are allocated.
                    return zeros(0)
            else:
                raise NotImplementedError('The species attributes only works with one species')
        def fset(self, value):
            if len(self.jslist) == 1:
                js = self.jslist[0]
                i1 = self.pgroup.ins[js] - 1
                i2 = i1 + self.pgroup.nps[js]
                self.pgroup.pid[i1:i2,pid-1] = value
            else:
                raise NotImplementedError('The species attributes only works with one species')

        def getpidname(self, **kw):
            return self.getpid(id=pid-1, **kw)

        setattr(cls, name + 'pid', pid)
        setattr(cls, name, property(fget, fset, None, None))
        setattr(cls, 'get' + name, getpidname)
        cls._addedpids.append(name)
