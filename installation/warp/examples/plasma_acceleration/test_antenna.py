from warp import *
import os
from warp.data_dumping.openpmd_diag import FieldDiagnostic, ParticleDiagnostic
from scipy.constants import c, m_e, e
from warp.field_solvers.laser.laser_profiles import GaussianProfile
from mpi4py import MPI
import numpy as np
from opmd_viewer.addons import LpaDiagnostics

# ------------------------------------------------------------------------------
# Resolution and stencil
# ------------------------------------------------------------------------------

dim       = "2d"
stencil   = 0

# Flags turning off auto decomp and let user specify its decomp
top.lautodecomp   = 1 # Particles
top.lfsautodecomp = 1 # fields
top.nhist = top.nt
top.lprntpara = false
top.lpsplots = false

# Time and spatial steps
N_step = 800
dxfact = 10
dtfact = 10
dx     = 0.0035
dt     = dx

if dim=="3d":
    dxfact *= 5
    dtfact *= 5

dx *= dxfact
dy = dz = dx
dt *= dtfact
dtcoef = dt / dx

# ------------------------------------------------------------------------------
# Laser characteristics
# ------------------------------------------------------------------------------

lambda_laser = 0.8e-6
w0           = 2.*np.pi*c/lambda_laser
laser_waist  = 1.5*lambda_laser
k0           = w0/c
tau          = 4.*lambda_laser/c
ctau         = c * tau
t_peak       = 1.8*tau
a0           = 1.
Eamp         = a0*w0*m_e*c/e

# ------------------------------------------------------------------------------
# Simulation box
# ------------------------------------------------------------------------------

w3d.xmmin = 0.
w3d.xmmax = w3d.xmmin + 20.*lambda_laser
w3d.ymmin = 0.
w3d.ymmax = w3d.ymmin + 20.*lambda_laser
w3d.zmmin = 0.
w3d.zmmax = w3d.zmmin + 20.*lambda_laser

w3d.dx    = dx*lambda_laser
w3d.dy    = dy*lambda_laser
w3d.dz    = dz*lambda_laser

w3d.nx    = nint((w3d.xmmax-w3d.xmmin)/w3d.dx)
w3d.nz    = nint((w3d.zmmax-w3d.zmmin)/w3d.dz)
w3d.ny    = nint((w3d.ymmax-w3d.ymmin)/w3d.dy)

w3d.xmmax = w3d.xmmin + w3d.nx*w3d.dx
w3d.ymmax = w3d.ymmin + w3d.ny*w3d.dy
w3d.zmmax = w3d.zmmin + w3d.nz*w3d.dz

print('w3d.dx, w3d.dy, w3d.dz = %f %f %f'%(w3d.dx, w3d.dy, w3d.dz))
print('w3d.nx, w3d.ny, w3d.nz = %d %d %d'%(w3d.nx, w3d.ny, w3d.nz))
print('w3d.xmmin, w3d.ymmin, w3d.zmmin = %f %f %f'%(w3d.xmmin, w3d.ymmin, w3d.zmmin))
print('w3d.xmmax, w3d.ymmax, w3d.zmmax = %f %f %f'%(w3d.xmmax, w3d.ymmax, w3d.zmmax))

# --- sets field boundary conditions
# --- longitudinal
w3d.bound0  = w3d.boundnz = openbc
# --- transverse
w3d.boundxy = openbc

# --- sets particles boundary conditions
# --- longitudinal
top.pbound0  = absorb
top.pboundnz = absorb
# --- transverse
top.pboundxy = absorb

if dim in ["2d"]:
    w3d.ny    = 2
    w3d.ymmin = - float(w3d.ny)/2.
    w3d.ymmax = + float(w3d.ny)/2.
    w3d.dy    = 1.
    nyguard   = 0

top.depos_order    = 3
top.efetch         = 4
top.lrelativ       = true

top.fstype = -1
package('w3d');generate()

# ------------------------------------------------------------------------------
# Laser Antenna
# ------------------------------------------------------------------------------

x0 = w3d.xmmax/2.
if dim=="3d":
    y0 = w3d.ymmax/2.
else:
    y0 = 0.
z0 = 5*w3d.dz

theta           = 0.
laser_vector    = np.array([- np.sin(theta), 0., np.cos(theta)])
laser_polvector = np.array([+ np.cos(theta), 0., np.sin(theta)])
laser_spot      = np.array([x0,y0,z0])

print("n = ", laser_vector)
print("P = ", laser_polvector)
print("r = ", laser_spot)

zf = - w3d.zmmax/2.

laser_func = GaussianProfile(k0, laser_waist, tau, t_peak,
                           a0, dim, focal_length= zf )

# ------------------------------------------------------------------------------
# Creation em3dsolver
# ------------------------------------------------------------------------------

em = EM3D( laser_func      = laser_func,
           laser_vector    = laser_vector,
           laser_polvector = laser_polvector,
           laser_spot      = laser_spot,
           laser_emax      = Eamp,
           stencil         = stencil,
           l_2dxz          = dim=="2d",
           l_1dz           = dim=="1d",
           dtcoef          = dtcoef,
           l_getrho        = 0,
           l_verbose       = 0)

print('register solver')
registersolver(em)
em.finalize()
loadrho()
print('done')

# ------------------------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------------------------

# Calculate the number of steps between each output
fielddiag_period = 684

# Add diagnostic
diag_f = FieldDiagnostic( period=fielddiag_period, top=top, w3d=w3d, em=em,
                      comm_world=comm_world, fieldtypes=["E", "B"] )
installafterstep( diag_f.write )

# ------------------------------------------------------------------------------
# Run steps
# ------------------------------------------------------------------------------

if (me==0):
      print("Total number of steps: ",N_step)
      tdeb=MPI.Wtime()

em.step(N_step)

tend=MPI.Wtime()
if(me==0):
  print("Final runtime (s): "+str(tend-tdeb))

# ------------------------------------------------------------------------------
# Test Diagnostics
# ------------------------------------------------------------------------------

print("=== Test waist, ctau and a0 ===")

# Load the diags OpenPMD
ts = LpaDiagnostics('./diags/hdf5/')

# Values of waist, ctau and a0 set in the script
laser_waist_simu = ts.get_laser_waist( iteration=ts.iterations[0], pol='x' )
ctau_simu        = ts.get_ctau( iteration=ts.iterations[0], pol='x' )
a0_simu          = ts.get_a0( iteration=ts.iterations[0], pol='x' )

# Compare the values and crash the code if not equal

assert np.isclose( laser_waist, laser_waist_simu, 5.e-2 ), "waist error. "
print("waist checked. ")
assert np.isclose(    ctau    ,     ctau_simu   , 5.e-2 ), "ctau error. "
print("ctau checked. ")
assert np.isclose(     a0     ,      a0_simu    , 5.e-2 ), "a0 error. "
print("a0 checked. ")
