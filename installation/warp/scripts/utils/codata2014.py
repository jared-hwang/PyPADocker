from ..warp import top
# *********** Constant:
# Physical constants, in one place so they stay consistent
# Values from 2014CODATA
# http:physics.nist.govcuuConstantsbibliography.html
top.amu = 1.660539040e-27  # [kg] # Atomic Mass Unit
top.clight = 2.99792458e+8  # [ms] # Speed of light in vacuum (exact)
top.echarge = 1.6021766208e-19 # [C] # Proton charge
top.emass = 9.10938356e-31 # [kg] # Electron mass
top.boltzmann = 1.38064852e-23 # [J/K]  # Boltzmann's constant
top.avogadro = 6.022140857e23        # Avogadro's Number
top.planck = 6.626070040e-34 # [J.s] # Planck's constant

# --- Magnetic constant = 4*pi*1.e-7
top.mu0 = 4*top.pi/10000000
# --- Conversion factor from joules to eV is just echarge
top.jperev = top.echarge
# --- Epsilon_0 calculated from speed of light and mu_0
top.eps0 = 1/(top.mu0*top.clight*top.clight)
# --- Create python versions of the constants
amu       = top.amu
clight    = top.clight
echarge   = top.echarge
emass     = top.emass
eps0      = top.eps0
jperev    = top.jperev
mu0       = top.mu0
boltzmann = top.boltzmann
avogadro  = top.avogadro
planck    = top.planck

