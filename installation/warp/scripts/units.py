import numpy as np

# --- Set handy units conversion factors, converting from named unit to MKS.
inch = 0.0254  # inches to meter
nm = 1.e-9     # nm to meter
um = 1.e-6     # um to meter
mm = 1.e-3     # mm to meter
cm = 0.01      # cm to meter
km = 1.e3      # km to meter
ps = 1.0e-12   # picosecond to second
ns = 1.0e-9    # nanosecond to second
us = 1.0e-6    # microsecond to second
ms = 1.0e-3    # millesecond to second
kV = 1.0e3     # kV to V
keV = 1.0e3    # keV to eV
MV = 1.0e6     # MV to V
MeV = 1.0e6    # MeV to eV
GV = 1.0e9     # GV to V
GeV = 1.0e9    # GeV to eV
mA = 1.0e-3    # mA to A
uC = 1.0e-6    # micro-Coulombs to Coulombs
nC = 1.0e-9    # nano-Coulombs to Coulombs
gauss = 1.e-4  # gauss to T
deg = np.pi/180.0  # degrees to radians
