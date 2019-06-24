w3d
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package W3D of code WARP
# 3D - PIC package of 3d particle code
# Alex Friedman, LLNL, (510)422-0827
# David P. Grote, LLNL, (510)423-7194
{
LARGEPOS = 1.0e+36 # This must be the same as in top.v
}

*********** InPltCtl3d dump:
# Controls for when the various plots are made
# Elements 0,1,2,3,... control plots for corresponding window
icrhoxy(0:NWINDOWS) integer  # rho contours zwindows
icrhozx(0:NWINDOWS) integer  # rho contours ywindows
icrhozy(0:NWINDOWS) integer  # rho contours xwindows
icphixy(0:NWINDOWS) integer  # phi contours zwindows
icphizx(0:NWINDOWS) integer  # phi contours ywindows
icphizy(0:NWINDOWS) integer  # phi contours xwindows
icrhoxy4            integer  # rho contours 4 zwindows / page
icrhozx4            integer  # rho contours 4 ywindows / page
icrhozy4            integer  # rho contours 4 xwindows / page
icphixy4            integer  # phi contours 4 zwindows / page
icphizx4            integer  # phi contours 4 ywindows / page
icphizy4            integer  # phi contours 4 xwindows / page

*********** InGen3d dump:
# General parameters which control the mechanics of the run (input qtys)
filt(5,3)            real  [1]  /15*0./    # Filtering coeffs for fieldsolve
rwallfac             real  [1]  /1./       # factor for g factor in w3dgen
bndftol              real  [1]  /1.e-6/    # relative tol for b-s-f iteration
bndfflag             logical    /.true./   # flag, bent-self-field iteration
bnprflag             logical    /.true./   # flag, printing b-s-f iterations
bnezflag             logical    /.true./   # flag, bend Ez correction
bnjtflag             logical    /.true./   # flag, bend "jump term" correction
bndfitmx             integer    /15/       # maximum number of b-s-f iterations
pipeshpe             character*8 /" "/
      # Shape of pipe for Kz capacity matrix (default circle; "quad" hyperbola)
      # Shape of pipe for 3D capacity matrix (default circle; "hyp" hyperbola)
l2symtry             logical    /.false./  # Turns on 2-fold symmetry
l4symtry             logical    /.false./  # Turns on 4-fold symmetry
lbeforefs  logical    /.false./  # Turns on call to python function "beforefs"
lafterfs   logical    /.false./  # Turns on call to python function "afterfs"
lbeforelr logical /.false./ # Turns on call to python function "beforelr"
lafterloadrho logical /.false./ # Turns on call to python function "afterloadrho"
lothereuser logical /.false./ # Turns on call to python function "othereuser"
lbeforescraper   logical /.false./ # Turns on call to python function "beforescraper"
lafterscraper   logical /.false./ # Turns on call to python function "afterscraper"
lcallscraper   logical /.false./ # Turns on call to python function "callscraper"
lcallparticleloader logical /.false./ # Turns on call to python function "callparticleloader"
luserappliedfields logical /.false./ # Turns on call to python function "userappliedfields"
                                     # just before particle advance
solvergeom integer    /0/  # Geometry of field solver
XYZgeom    integer    /0/  # 3D-XYZ geometry will be used if solvergeom=XYZgeom
RZgeom     integer    /1/  # axisymmetric RZ geometry will be used if solvergeom=RZgeom
AMRgeom    integer    /2/  # 3D geometry using AMR if solvergeom=AMRgeom
XZgeom     integer    /3/  # 2-D sheath geometry
XYgeom     integer    /4/  # 2-D planar geometry
Zgeom      integer    /5/  # 1-D planar geometry
Rgeom      integer    /6/  # 1-D radial geometry
Ygeom      integer    /7/  # 1-D planar geometry

*********** InDiag3d dump:
lgetese3d logical /.true./ # Sets whether electrostatic-energy is calculated,
                           # the product of rho and phi.
                           # This calculation is expense so this flag was
                           # added to turn it off if desired.
lgtlchg3d logical /.true./ # Sets whether the line-charge is gathered.
                           # This calculation is expense so this flag was
                           # added to turn it off if desired.
lrhodia3d logical /.true./ # Sets whether rho diagnostics are done,
                           # calculating rhomid and rhomax along the axis.
                           # This calculation is expense so this flag was
                           # added to turn it off if desired.
lpltfld3d logical /.false./ # When true, the compiled pltfld3d calls the
                            # interpreter function 'pltfld3d' with the
                            # appropriate arguments and then returns.
lgetvzofz logical /.true./ # Sets wether vzofz is calculated
lsrhoax3d logical /.true./ # Sets wether the charge density on axis is calculated
lsphiax3d logical /.true./ # Sets wether the potential on axis is calculated
lsezax3d  logical /.true./ # Sets wether the longitudinal electric field on axis is calculated
lsetcurr             logical /.true./ # Sets wether the current on axis is calculated


*********** InPart3d dump:
# Particle input quantities (input qtys)
zjig                      real  [1]    /0./
   # Controls "jiggling" of initial positions in z, for grid loading
distrbtn                  character*20  /"none"/
   # transvese particle distribution:
   #  "none"  = no distribution (user must define in interpreter)
   #  "semigaus" or "SG" or "SemiGaussian" = semi-Gaussian
   #  "K-V" or "KV" = Kapchinskij-Vladimirskij
   #  "KV0" = alternative Kapchinskij-Vladimirskij
   #  Based on transformations of Hamiltonian defined equilibria in continuous
   #  applied focusing channels (not presently implemented, future addition):
   #    "WB" or "Waterbag" = water-bag
   #    "PA" or "Parabolic" = parabolic
   #    "TE" or "ThermalEquilibrium" = thermal equilibrium
   #  Based on zero-space-charge Courant-Snyder invariants of lattices with
   #  axially varying linear applied focusing:
   #    "WB0" or "Waterbag0" = water-bag
   #    "PA0" or "Parabolic0" = parabolic
   #    "GA0" or "Gaussian0" = Gaussian
   #  "preload" = User preloaded distribution - in this case, the species
   #              quantities are calculated from base quantities
distr_l                   character*8  /"neuffer"/
   # longitudinal velocity distribution for cigar load: either "neuffer"
   # for hard edged distribution (Vlasov equilibrium), or "gaussian" for
   # gaussian distribution with same rms variation in z.
distr_t                   character*8  /"gaussian"/
   # transverse velocity distribution: either "uniform"
   # for uniform distribution, or "gaussian" for gaussian distribution.
cigarld                   logical      /.false./
   # specifies whether or not to do a cigar load (finite beam)
xrandom                   character*8 /"pseudo"/
   # random numbers used for x,y, and z. "pseudo","fibonacc","digitrev","grid"
vtrandom                  character*8 /"digitrev"/
   # random numbers used for vx and vy. "pseudo" or "digitrev"
vzrandom                  character*8 /"digitrev"/
   # random numbers used for vz. "pseudo" or "digitrev"
ldprfile                  character*8 /"polar"/
   # load profile "polar", "streamls" or "stripes"
cylinder                  logical      /.false./
   # specifies whether or not to load a cylinder
hollow                    integer      /0/
   # specifies type of hollow beam(0:none, 1:linear in r^2, 2:n(r)~(h+(1-h)r^2))
hollow_h                  real     /.5/
   # Hollowness factor used for hollow=2, Note: cannot be 1.
nxstripe                  integer  /0/
   # Number of x stripes on which particles are loaded, for grid loading
nystripe                  integer  /0/
   # Number of y stripes on which particles are loaded, for grid loading
nzstripe                  integer  /0/
   # Number of z stripes on which particles are loaded, for grid loading
nfibgrps                  integer  /0/
   # Number used fibonacci groups
fibg1                     integer  /1958/
   # Number used for fibonacci loading (See S.K. Zaremba ref. in setptcls)
fibg2                     integer  /252/
   # Number used for fibonacci loading (See S.K. Zaremba ref. in setptcls)
fibg3                     integer  /414/
   # Number used for fibonacci loading (See S.K. Zaremba ref. in setptcls)
fibg4                     integer  /0/
   # Number used for fibonacci loading (See S.K. Zaremba ref. in setptcls)
dig1                      integer /2/
   # first base used for digit reversed loading
dig2                      integer /3/
   # second base used for digit reversed loading
dig3                      integer /5/
   # third base used for digit reversed loading
dig4                      integer /7/
   # fourth base used for digit reversed loading
dig5                      integer /11/
   # fifth base used for digit reversed loading
dig6                      integer /13/
   # sixth base used for digit reversed loading
dig7                      integer /17/
   # seventh base used for digit reversed loading
dig8                      integer /23/
   # eigth base used for digit reversed loading
dig9                      integer /29/
   # ninth base used for digit reversed loading
dig10                     integer /31/
   # tenth base used for digit reversed loading
nrdist                    integer
   # number of data points for arbitrary particle distribution in r
rdist(0:nrdist)           _real
   # User specified arbitrary particle distribution in r
   # Data is assumed to be uniformly spaced in a normalized range from 0 to 1
   # with 0 being at rmin and 1 being at rmax.
nrmrdist(0:nrdist)          _real
   # Normalized rdist. sum(nrmrdist)=1. Calculated automatically.
intrdist(0:nrdist)         _real
   # Integral of rdist.  Calculated automatically.
nvrdist                   integer
   # number of data points for arbitrary particle velocity distribution in r
vthrofr(0:nvrdist)         _real
   # User specified radial thermal velocity as a function of radius.
vrbarofr(0:nvrdist)         _real
   # User specified average radial velocity as a function of radius.
nzdist                    integer
   # number of data points for arbitrary particle distribution in z
zdist(0:nzdist)           _real
   # User specified arbitrary particle distribution in z
   # Data is assumed to be uniformly spaced in a normalized range from 0 to 1
   # with 0 being at zimin and 1 being at zimax.
lnormalize_zdist           logical /.true./
   # When true, zdist will be normalized so that the total number of particles
   # and total charge are unchanged. When false, zdist will be used as is, as
   # a direct scaling factor on the line-charge - where zdist is one, the
   # specified line-charge will be loaded. This case has the nice property
   # that if zdist is changed in one place, the number of particles loaded
   # elsewhere will not be affected.
nrmzdist(0:nzdist)          _real
   # Normalized zdist. sum(nrmzdist)=1. Calculated automatically.
intzdist(0:nzdist)         _real
   # Integral of zdist.  Calculated automatically.
nenvofz                    integer
   # Number of data points for the axially varying envelope
aofz(0:nenvofz)            _real
   # User specified axially varying X envelope.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meters.
bofz(0:nenvofz)            _real
   # User specified axially varying Y envelope.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meters.
apofz(0:nenvofz)           _real
   # User specified axially varying X' envelope.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of radians.
bpofz(0:nenvofz)           _real
   # User specified axially varying Y' envelope.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of radians.
xofz(0:nenvofz)            _real
   # User specified axially varying X centroid
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meters.
yofz(0:nenvofz)            _real
   # User specified axially varying Y centroid
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meters.
xpofz(0:nenvofz)           _real
   # User specified axially varying X' centroid
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of radians.
ypofz(0:nenvofz)           _real
   # User specified axially varying Y' centroid
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of radians.
nvbeamofz                  integer
   # Number of data points for the axially varying axial velocity
vbeamofz(0:nvbeamofz)      _real
   # User specified axially varying axial velocity
   # Data is assumed to be uniformly spaced in a normalized range from 0 to 1
   # with 0 being at zimin and 1 being at zimax.
nvthzofz                   integer
   # Number of data points for the axially varying axial velocity spread
vthzofz(0:nvthzofz)        _real
   # User specified axially varying axial velocity spread.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meter/second.
nemitofz                   integer
   # Number of data points for axially varying transverse emittances
   # emitxofz and emityofz.
emitxofz(0:nemitofz)       _real
   # User specified axially varying transverse x-emittances
   # Data is assumed to be uniformly spaced in a normalized range from 0 to 1
   # with 0 being at zimin and 1 being at zimax.
emityofz(0:nemitofz)       _real
   # User specified axially varying transverse y-emittances
   # with same mesh as emitxofz.
   #

*********** InMesh3d dump:
# Mesh specifications (input qtys)
xmmin   real  [m]  /0./             # Lower limit of mesh
xmmax   real  [m]  /0./             # Upper limit of mesh
ymmin   real  [m]  /0./             # Lower limit of mesh
ymmax   real  [m]  /0./             # Upper limit of mesh
zmmin   real  [m]  /0./             # Lower limit of mesh
zmmax   real  [m]  /0./             # Upper limit of mesh
nx      integer    /0/              # Mesh points are 0,...,nx
ny      integer    /0/              # Mesh points are 0,...,ny
nz      integer    /0/              # Full size of nz
nxlocal integer    /0/  +parallel   # Mesh points are 0,...,nxlocal
nylocal integer    /0/  +parallel   # Mesh points are 0,...,nylocal
nzlocal integer    /0/  +parallel   # Mesh points are 0,...,nzlocal
xmminlocal  real [m] /0./ +parallel # Local value of xmmin
xmmaxlocal  real [m] /0./ +parallel # Local value of xmmax
ymminlocal  real [m] /0./ +parallel # Local value of ymmin
ymmaxlocal  real [m] /0./ +parallel # Local value of ymmax
zmminlocal  real [m] /0./ +parallel # Local value of zmmin
zmmaxlocal  real [m] /0./ +parallel # Local value of zmmax
xmminglobal real [m]                # Global value of xmmin
xmmaxglobal real [m]                # Global value of xmmax
ymminglobal real [m]                # Global value of ymmin
ymmaxglobal real [m]                # Global value of ymmax
zmminglobal real [m]                # Global value of zmmin
zmmaxglobal real [m]                # Global value of zmmax
nxguardphi integer /1/ # Number of guard cells in x for phi arrays
nyguardphi integer /1/ # Number of guard cells in y for phi arrays
nzguardphi integer /1/ # Number of guard cells in z for phi arrays
nxguardrho integer /0/ # Number of guard cells in x for rho arrays
nyguardrho integer /0/ # Number of guard cells in y for rho arrays
nzguardrho integer /0/ # Number of guard cells in z for rho arrays
nxguarde   integer /0/ # Number of guard cells in x for e arrays
nyguarde   integer /0/ # Number of guard cells in y for e arrays
nzguarde   integer /0/ # Number of guard cells in z for e arrays

******* GridBoundary3d dump:
bound0    integer /0/  # Type of boundary condition at plane z=0
                       # 0 is constant potential, 1 is zero normal derivative,
                       # and 2 is periodic
boundnz   integer /0/  # Type of boundary condition at plane z=nz
                       # 0 is constant potential, 1 is zero normal derivative,
                       # and 2 is periodic
boundxy   integer /0/  # Type of boundary condition at sides
                       # 0 is constant potential, 1 is zero normal derivative,
                       # and 2 is periodic
bounds(0:5) integer /6*0/ # Type of boundaries at edge of mesh, in order of
                          # lower, upper for x, y, z.
lzerophiedge logical /.true./ # When true and when gridmode == 0, the edge of
                       # the phi array is zeroed out. This clears phi at any
                       # conductor points on the edge of the mesh.

*********** Damped_eom dump:
# All quantities associated with the damped mover algorithm
eomdamp     real  [1] /0./ # EOM Damping param and switch
itdamp      integer   /5/  # Timestep when damping starts
exeomoldpid integer        # Electric field pid index at particle one step back
eyeomoldpid integer        # Electric field pid index at particle one step back
ezeomoldpid integer        # Electric field pid index at particle one step back
exeomlagpid integer        # Electric field pid index, lag-averaged for damping
eyeomlagpid integer        # Electric field pid index, lag-averaged for damping
ezeomlagpid integer        # Electric field pid index, lag-averaged for damping

*********** Mesh3d:
# Large arrays: potential and charge density, plot scratch, etc.
xmesh(0:nx)             _real [m] +dump # X coordinates of mesh points
ymesh(0:ny)             _real [m] +dump # Y coordinates of mesh points
zmesh(0:nz)             _real [m] +dump # Z coordinates of global mesh points
xmeshlocal(0:nxlocal)   _real [m] +dump +parallel # x coordinates of local mesh points
ymeshlocal(0:nylocal)   _real [m] +dump +parallel # y coordinates of local mesh points
zmeshlocal(0:nzlocal)   _real [m] +dump +parallel # Z coordinates of local mesh points

*********** Fields3dSolver:
bndferr                 real    [1]     # Error after bent-self-field solve
bndfit                  integer [1]     # Iters used, bent-self-field solve
nmxy                    integer /0/ +dump # larger of nx, ny
nmxyz                   integer /0/ +dump +parallel # largest of nx, ny, nz
phi(:,:,:)              _real [V] +parallel # Electrostatic potential
rho(:,:,:)              _real [C/m**3] +parallel # Charge density
selfe(:,:,:,:)          _real [V/m] +parallel # Self E field
                            # calculated from phi via finite difference.
                            # Only used when top.efetch = 3
attx(0:nx-1)            _real           # Attenuation factor as fcn. of kx
atty(0:ny-1)            _real           # Attenuation factor as fcn. of ky
attz(0:nz)              _real           # Attenuation factor as fcn. of kz
kxsq(0:nx-1)            _real [1/m**2]  # Discrete analog to kx^2/4Pi
kysq(0:ny-1)            _real [1/m**2]  # Discrete analog to ky^2/4Pi
kzsq(0:nz)              _real [1/m**2]  # Discrete analog to kz^2/4Pi
rstar(-1:nzlocal+1)     _real [m]       # Radius of curv of refrnce orbit
phiprv(0:nx,0:nzlocal)  _real [V]       # Prev phi at y_mid, for error test
phisav(0:nx,-1:nzlocal) _real [V]       # Phi at current y slice (scratch)
scrtch(-1:nmxyz+1,-1:nmxy+1)  _real           # Scratch for fieldsolve, plots
xywork3d(2,0:nx,0:ny)   _real           # Work space for transverse FFTs
zwork(2,0:nx,0:nz)      _real           # Work space used to optimize vsftz

*********** Fields3dParticles parallel:
nxp  integer /0/ # Number of grid cells in x axis for phip and rhop
nyp  integer /0/ # Number of grid cells in y axis for phip and rhop
nzp  integer /0/ # Number of grid cells in z axis for phip and rhop
nxpextra  integer /0/ # Extra number of grid cells in x axis for phip and rhop
nypextra  integer /0/ # Extra number of grid cells in y axis for phip and rhop
nzpextra  integer /0/ # Extra number of grid cells in z axis for phip and rhop
xmminp real      # Lower limit of x for grid for particles
xmmaxp real      # Upper limit of x for grid for particles
ymminp real      # Lower limit of y for grid for particles
ymmaxp real      # Upper limit of y for grid for particles
zmminp real      # Lower limit of z for grid for particles
zmmaxp real      # Upper limit of z for grid for particles
nrhopndtscopies3d integer /0/ # Copy of nrhopndtscopies from top
nsndts3d          integer /0/ # Copy of nsndts from top
nsndtsphi3d       integer /0/ # Copy of nsndtsphi from top
phip(:,:,:) _real # Potential used by the particles to calculate
                  # the field from the solution of Poisson's equation.
                  # This will be pointed to each of the different ndts
                  # groups as needed.
rhop(:,:,:) _real # Charge density from the particles.
                  # This will be pointed to each of the different ndts
                  # groups as needed.
selfep(:,:,:,:) _real [V/m] # Self E field for the particles,
                            # calculated from phi via finite difference.
                            # Only used when top.efetch = 3
rhopndts(-nxguardrho:nxp+nxguardrho,
         -nyguardrho:nyp+nyguardrho,
         -nzguardrho:nzp+nzguardrho,
         0:nrhopndtscopies3d-1,
         0:nsndts3d-1) _real +fassign
                 # Charge density from the particles
                 # for groups with different time step sizes.
                 # This includes the time averaged charge density from
                 # faster particles and the old rho from the slower particles.
phipndts(-nxguardphi:nxp+nxguardphi,
         -nyguardphi:nyp+nyguardphi,
         -nzguardphi:nzp+nzguardphi,
         0:nsndtsphi3d-1) _real +fassign
                 # Potential from the particles
                 # for groups with different time step sizes.
                 # This includes the effect of time averaged charge density
                 # from faster particles and the old rho from the slower
                 # particles. With only sampling, the same phip is used
                 # for all groups, i.e. nsndtsphi==1.
lrhofinalized logical /.true./
                 # Flags whether the boundary conditions and other processing
                 # of rho has been done, so that its ready to be used in a
                 # field solve. Anytime rhop is changed, this flag is set
                 # to false. Just before a field solve, the flag is checked,
                 # and if needed, the appropriate operations done.

nsselfb3d integer /0/ # Copy of nsselfb from top
rhopselfb(-nxguardrho:nxp+nxguardrho,
          -nyguardrho:nyp+nyguardrho,
          -nzguardrho:nzp+nzguardrho,
          0:nsselfb3d-1) _real +fassign
                 # Charge density from the particles
                 # for groups which require correction for their self B.
                 # This includes the time averaged charge density from
                 # faster particles and the old rho from the slower particles.
phipselfb(-nxguardphi:nxp+nxguardphi,
          -nyguardphi:nyp+nyguardphi,
          -nzguardphi:nzp+nzguardphi,
          0:nsselfb3d-1) _real +fassign
                 # Temporary copy of the potential from the particles
                 # for groups which require correction for their self B.
                 # This includes the effect of time averaged charge density
                 # from faster particles and the old rho from the slower
                 # particles.

nxc  integer /0/ # Number of grid cells in x axis for isnearbycond
nyc  integer /0/ # Number of grid cells in y axis for isnearbycond
nzc  integer /0/ # Number of grid cells in z axis for isnearbycond
isnearbycond(0:nxc,0:nyc,0:nzc) _integer

*********** FieldSolveAPI:
jmin   integer    /0/  # index start of current particle chunk in top.pgroup
                       # (PYTHON indexed: based 0)
jmax   integer    /0/  # index end of current particle chunk in top.pgroup
                       # (PYTHON indexed: based 0)
jsfsapi     integer
js1fsapi    integer /-1/
js2fsapi    integer /-1/
ndtsfsapi   integer
npfsapi     integer
ipminfsapi  integer
lzerorhofsapi logical /.true./
lfinalizerhofsapi logical
pgroupfsapi _ParticleGroup +fassign
api_xlf2    logical /.false./
exfsapi(:) _real
eyfsapi(:) _real
ezfsapi(:) _real
bxfsapi(:) _real
byfsapi(:) _real
bzfsapi(:) _real
xfsapi(:) _real
yfsapi(:) _real
zfsapi(:) _real
phifsapi(:) _real
afsapi(:,:) _real

%%%%%%%%%%% Grid3dtype:
nx    integer /-1/
ny    integer /-1/
nz    integer /-1/
nxguard integer /0/
nyguard integer /0/
nzguard integer /0/
dx    real /LARGEPOS/
dy    real /LARGEPOS/
dz    real /LARGEPOS/
dxi   real /LARGEPOS/
dyi   real /LARGEPOS/
dzi   real /LARGEPOS/
xmin  real /+LARGEPOS/
xmax  real /-LARGEPOS/
ymin  real /+LARGEPOS/
ymax  real /-LARGEPOS/
zmin  real /+LARGEPOS/
zmax  real /-LARGEPOS/
grid(-nxguard:nx+nxguard,
     -nyguard:ny+nyguard,
     -nzguard:nz+nzguard) _real

%%%%%%%%%%% Grid2dtype:
nx    integer /-1/
ny    integer /-1/
nxguard integer /0/
nyguard integer /0/
dx    real /LARGEPOS/
dy    real /LARGEPOS/
dxi   real /LARGEPOS/
dyi   real /LARGEPOS/
xmin  real /+LARGEPOS/
xmax  real /-LARGEPOS/
ymin  real /+LARGEPOS/
ymax  real /-LARGEPOS/
grid(-nxguard:nx+nxguard,
     -nyguard:ny+nyguard) _real

*********** BoltzmannElectrons dump:
# Parameters controlling the Boltzmann-Electrons.
nberegions       integer /1/
     # Number of regions where Boltzmann-Electrons are
iondensity(nberegions)  /0./    _real
     # Base ion charge density, multiplier on Boltzmann
     # exponential for field solvers including Boltzmann electrons
electrontemperature(nberegions) _real /0./ [eV]
     # Electron temperature in units of eV
     # for field solvers including Boltzmann electrons
plasmapotential(nberegions) /0./     _real
     # Potential of plasma for field solvers including Boltzmann electrons
electrondensitymaxscale(nberegions) _real /2./
     # Limit of electron density relative to iondensity. This should be large
     # enough to not affect the solution, but small enough to prevent
     # divergences during field solve.
logelecdenmaxscale(nberegions)  _real # Log of electrondensitymaxscale,
       # calculated by code
xbemin(nberegions) _real /-LARGEPOS/ [m]
     # Minimum of the extent in x over which Boltzmann electrons are included
xbemax(nberegions) _real /+LARGEPOS/ [m]
     # Maximum of the extent in x over which Boltzmann electrons are included
ybemin(nberegions) _real /-LARGEPOS/ [m]
     # Minimum of the extent in y over which Boltzmann electrons are included
ybemax(nberegions) _real /+LARGEPOS/ [m]
     # Maximum of the extent in y over which Boltzmann electrons are included
zbemin(nberegions) _real /-LARGEPOS/ [m]
     # Minimum of the extent in z over which Boltzmann electrons are included
zbemax(nberegions) _real /+LARGEPOS/ [m]
     # Maximum of the extent in z over which Boltzmann electrons are included
lincludeiondensityinrho(nberegions) _logical /0/
luseparticledensity(nberegions) _integer /0/
                     # When nonzero, use density from the particles
                                # rather than the specified constant iondensity.
                                # For 1, uses local density; for 2, uses value
                                # averaged over z.
                                # Only applies now to the RZ solver.
lclampphitophimax(nberegions) _logical /0/
     # When true, the potential is clamped to phimax, which is calculated
     # assuming that the max electron density is the max of iondensity and
     # the density from ion particles. Only applies now to the RZ solver.
lnormtoavboltzfac   logical /.false./  # If true, divide Boltzmann factor by its z average.
                                       # Use in conjunction with luseparticledensity = 2
luselinboltz    logical    /.false./   # If 1, use linearized Boltzmann response
                                       # to get Boltzmann electron response whose z-av
                                       # value = z-av value of ion density (quasineutral)
liondensitygrid3d(nberegions)       _logical /0/
useriondensitygrid3d    Grid3dtype # Allows a spatially varying iondensity,
                                   # or is used when luseparticleldensity is
                                   # true
lelectrontemperaturegrid3d(nberegions) _logical /0/
electrontemperaturegrid3d  Grid3dtype # Allows a spatially varying electron
                                      # temperature
liondensitygrid2d(nberegions)       _logical /0/
iondensitygrid2d        Grid2dtype # Allows a spatially varying iondensity,
                                   # or is used when luseparticleldensity is
                                   # true
lelectrontemperaturegrid2d(nberegions) _logical /0/
electrontemperaturegrid2d  Grid2dtype # Electron temperature on code's 2D grid

*********** Picglb3d dump:
# Globally useful quantities for PIC simulation
dx               real   [m]  /0./      #  mesh spacing in x
dy               real   [m]  /0./      #  mesh spacing in y
dz               real   [m]  /0./      #  mesh spacing in z
nxyz             integer /0/
   # size of a field array, (nx+1)*(ny+1)*(nz+1)
ix_axis          integer [1] # x location of axis in mesh
iy_axis          integer [1] # y location of axis in mesh
iz_axis          integer [1] # z location of axis in mesh
ixlocal_axis     integer [1] # x location of axis in mesh relative to local coordinates
iylocal_axis     integer [1] # y location of axis in mesh relative to local coordinates
izlocal_axis     integer [1] # z location of axis in mesh relative to local coordinates

*********** InjectVars3d dump:
inj_ninj             integer  # Auto set to either 1 or ninject. Set to 1
                              # when no emitting source are within two grid
                              # cells of each other.
inj_ns               integer  # Auto set to either 1 or ns. Set to 1
                              # when only one species is being injected from
                              # each source.
inj_ij               integer  # Currently active injection source Useful for
                              # generateuserparticlesforinjection.
inj_js               integer  # Current speices being injected. Useful for
                              # generateuserparticlesforinjection.
inj_nx               integer  /0/  # size injection arrays in x
inj_ny               integer  /0/  # size injection arrays in y
inj_nz               integer  /0/  # size injection arrays in z
inj_dx               real [m] /0./ # mesh spacing in x for injection
inj_dy               real [m] /0./ # mesh spacing in y for injection
inj_dz               real [m] /0./ # mesh spacing in z for injection
inj_dz0              real [m] /0./ # mesh spacing in z for injection
inj_xwide            integer  /2/  # number of cells in x for each emitting pad
inj_ywide            integer  /2/  # number of cells in y for each emitting pad
inj_addfdz           real     /0./ # fraction of dz to add to z position at
                                   # injection
l_inj_rz             logical /.false./ # if true, make RZ injection with
                                       # variable weights and axisymmetric
                                       # injection grid
l_inj_rz_grid        logical /.false./ # if true, axisymmetric injection grid
                                       # is used
l_inj_regular        logical /.false./ # if true, inject one particle at each
                                       # grid node and adjust weight accordingly
l_inj_delay_temp     logical /.false./ # if true, add temperature only after
                                       # particles at distance inj_dtemp from
                                       # emitter
l_inj_addtempz_abs   logical /.false./ # if true, longitudinal thermal velocity
                                       # is positive
l_inj_zero_theta     logical /.false./ # When true, and l_inj_rz or
                                       # l_inj_rz_grid is true, then particles
                                       # are only created with the transverse
                                       # angle theta=0.
l_inj_rec_inittime   logical /.false./ # if true, time of creation is recorded
                                       # in pid
l_inj_rec_initradius logical /.false./ # if true, radius of creation is
                                       # recorded in pid
l_inj_exact          logical /.false./ # if true, position and angle of
                                       # injected particle computed
                                       # analytically rather than interpolated
l_inj_area           logical /.true./  # if false, when l_inj_rz=true, adjust
                                       # inj_dx so that inj_area is not used
                                       # (no effect if l_inj_rz=false)
l_inj_no_rho_on_emit logical /.false./ # If true, no rho deposited on emitter
l_inj_use_rho_with_mr logical /.true./ # When true, use the deposited rho
                                       # in the mr calculation of phi
                                       # in the emitting region.
l_inj_user_particles logical /.false./ # When true, user specified particles
                                       # will be injected using arrays from
                                       # Setpwork3d. Only works with
                                       # inject=1.
l_inj_user_particles_z logical /.false./ # When true, the user also provides
                                         # the z-position of specified injected
                                         # particles. They are set
                                         # automatically otherwise
                                         # from zinject and rinject (default).
l_inj_user_particles_v logical /.false./ # When true, the user also provides
                                         # the velocity of specified injected
                                         # particles. They are set
                                         # automatically otherwise (default).
l_inj_user_particles_dt logical /.false./ # When true, the user also provides
                                          # the fractional time steps
                                          # of specified injected particles.
                                          # They are set automatically
                                          # otherwise (default).
l_inj_zmminmmaxglobal   logical /.false./ # When true, sets inj_zmmin/max to
                                          # global values in parallel
inj_xmmin(inj_ninj)  _real [m] /0./ # Min x extent of injection mesh
inj_ymmin(inj_ninj)  _real [m] /0./ # Min y extent of injection mesh
inj_grid(0:inj_nx,0:inj_ny,inj_ninj) _real [m]
   # Grid giving axial field grid location of injection sources in the lab frame
inj_angl(0:inj_nx,0:inj_ny,inj_ninj) _real
   # Grid giving angle of injection sources for each transverse location
inj_throttle(0:inj_nx,0:inj_ny,inj_ninj,inj_ns) _real /1./
   # Modifies the rate of particle injection
inj_lambda(0:inj_nx,0:inj_ny,inj_ninj,inj_ns) _real /1./
  # Spatially varying coefficient for thermionic emission
inj_np(0:inj_nx,0:inj_ny,inj_ninj,inj_ns)   _real
   # Grid holding number of particles injected on the current time step.
inj_prev(0:inj_nx,0:inj_ny,inj_ninj,inj_ns) _real
   # Grid holding number of particles injected on previous time step.
inj_npactual(0:inj_nx,0:inj_ny,inj_ninj,inj_ns)   _real
   # Grid holding actual number of particles injected on the current time step.
inj_rho(0:inj_nx,0:inj_ny,inj_ninj)  _real
   # Surface charge density at the emitting surface.
inj_phi(0:inj_nx,0:inj_ny,inj_ninj)  _real
   # Electrostatic potential at the emitting surface.
inj_ex(0:inj_nx,0:inj_ny,inj_ninj)  _real
inj_ey(0:inj_nx,0:inj_ny,inj_ninj)  _real
inj_area(0:inj_nx,0:inj_ny,inj_ninj)  _real
   # Inverse of the fraction of the grid cell's area within the emitting surface
inj_q(0:inj_nx,0:inj_ny,0:inj_nz,inj_ninj)  _real
   # Charge in bins in the emitting region.
inj_phi_3d(0:inj_nx,0:inj_ny,0:inj_nz,inj_ninj)  _real
   # Electrostatic potential in the emitting region.
inj_ex_3d(0:inj_nx,0:inj_ny,0:inj_nz,inj_ninj)  _real
   # Tangential field in the emitting region.
inj_ey_3d(0:inj_nx,0:inj_ny,0:inj_nz,inj_ninj)  _real
   # Tangential field in the emitting region.
inj_ez_3d(0:inj_nx,0:inj_ny,0:inj_nz,inj_ninj)  _real
   # Normal field in the emitting region.

*********** Setpwork3d:
# Scratch arrays for subroutine setptcls
npgrp             integer /0/
npidgrp           integer /0/
indx(npgrp)      _integer
xt(npgrp)        _real
yt(npgrp)        _real
zt(npgrp)        _real
rt(npgrp)        _real
tt(npgrp)        _real
uxt(npgrp)       _real
uyt(npgrp)       _real
uzt(npgrp)       _real
perpscal(npgrp)  _real
at(npgrp)        _real
apt(npgrp)       _real
bt(npgrp)        _real
bpt(npgrp)       _real
xct(npgrp)       _real
xpct(npgrp)      _real
yct(npgrp)       _real
ypct(npgrp)      _real
pidt(npgrp,npidgrp)  _real

*********** Multipole dump:
# Electrostatic multipole moments of the electrostatic potential
nmom                  integer  /1/  # number of terms in multipole expansion
nzmom                 integer  /0/  # Number of z grid points in the
                                    # multipolar decomposition, should be the
                                    # same as w3d.nz
lazi(nmom)            _integer      # azimuthal mode numbers of moments
irpow(nmom)           _integer      # radial powers of moments
rmomcex(nmom,0:nzmom) _real    [1]  # multipoles, cos(theta) expansion terms
rmomsex(nmom,0:nzmom) _real    [1]  # multipoles, sin(theta) expansion terms

*********** Apertures dump:
# This group contains data describing the size and location of a circular
# aperture which is used to scrape particles. At the aperture, the transverse
# E field is calculated seperately, to explicitly include the conductor.
napertures integer /0/          # Number of special apertures
aper_zmax  integer /1/          # Maximum z length of the apertures
aper_zs(napertures)   _real [m] # Z start in lab frame of special apertures.
aper_ze(napertures)   _real [m] # Z end in lab frame of special apertures.
aper_rad(napertures)  _real [m] # Radius of special apertures.
aper_x(napertures)    _real [m] # X of special apertures center.
aper_y(napertures)    _real [m] # Y of special apertures center.
aper_volt(napertures) _real [V] # Voltage on special apertures.
aper_ex(0:nx,0:ny,-1:aper_zmax,napertures) _real [V/m]
                            # Calculated Ex on special apertures.
aper_ey(0:nx,0:ny,-1:aper_zmax,napertures) _real [V/m]
                            # Calculated Ey on special apertures.

*********** DKInterp dump:
# This group contains inputs and other saved data for the drift-kinetic
#  interpolation
# Note several species dependent quantities; assume max of ns is 1000
m_over_q(1000)  _real [mks] /0./ # mass over charge, calc. by code
qovermsq(1000)  _real [mks] /0./ # (charge/mass)**2, calc. by code
alpha0(1000)     _real     /0./
                # Interpolation parameter if not to be automatically set
                # Note alpha0=0 for pure drift kinetics
acntr(1000)      _real /.5/    # centering parameter for predictor-corrector
usealphacalc(1000) _real  /1./  # fraction of calculated interpolation parameter
                         # to use; will use (1-usealphacalc)*alpha0.
notusealphcalc(1000) _real /0./ # alpha0*(1-usealphacalc), calculated by code
dksmall       real  /1.e-20/ # small parameter for safe divides
igradb   integer  /2/    #  parameter to select method of calculating grad B
                         # 1 for lookup table, 2 for assumed quadrupole
                         # 3 for lookup in z, quad in x,y
interpdk(1000) _integer  /0/ # parameter specifies whether and how to do orbit
                         # interpolation: 0, full orbit. 1, interpolate
impinterp     integer   /0/  # set to 1 for implicit version of interpolated
                             # mover
alphcoef      real    /0.25/ # coefficient multiplying (omegac dt)**2
                             # in setting alpha
ipalpha       integer /1/   # power of sqrt(1+omegadt) in setting alpha.
                            # allowed values for now are  1 or 2.
                            # If any other value, then uses the arbitrary
                            # (real) power palpha in a slower calculation
                            # Note, only absolute value matters for alpha.
                            # Sign of ipalpha determines whether omegadt
                            # is multiplied by 2 (2x as long dt) for
                            # ipredcor = 0.  (Done for ipalpha>0, not
                            # for ipalpha < 0.)
palpha        real  /0.5/    # power of (1+omegadt) in setting alpha
                            # used only if ipalpha != 1 or 2.
npredcorv      integer /1/   # number of times to do predictor-corrector
                            # update of effective velocities
npcmax         integer /-1/  # maximum number of corrector steps for
                            # global predictor-corrector. If it is negative,
                            # do old-style prediction from step 0 to 1; number
                            # of iterations is then -npcmax - 1.
ipcmax         integer      # max number of predictor-corrector iterations
                            # used (0 if no interpolated species, npcmax if
                            # there IS an interpolated species and npcmax>0).
                            # Calculated by code.
ipredcor       integer /0/  # flag to designate predictor/corrector for step3d
icalceps       integer /0/  # flag to calculate polarization
                            # correction to Poisson solver if interpolating
iusing2deps    integer      # flag set by code if 2D polarization is calcluated
polfac               /1./   # multiplier on polarization term
lnonlinpol logical /.true./ # if 1, use nonlinear polarization; if 0, linear with a constant
                            # value for omega_p^2/omega_c^2 given by omegpsqomegcsq
omegpomegcsq(1000) _real /1./  # The constant value of omega_p^2/omega_c^2 per species
igetextdrift   integer /0/  # flag to invoke a python method "getextradrift"
                            # added directly to veff's after getveff.
epsilon2d       Grid2dtype # Allows a spatially varying dielectric constant
epszfac2d       Grid2dtype # spatially varying epsilon_z/epsilon_x

*********** DKInterptmp:
# This group contains temporary data for the drift-kinetic interpolation
npint          integer   /0/     # dimension for interpolation arrays
grdbsq(3,npint) _real [t**2/m]   # components of gradbsq (z,x,y)
dbsqdz(npint)   _real [T**2/m]    # dB^2/dz
bsqi(npint)     _real [1/T**2]    # 1/B^2
bsq(npint)      _real [1/T**2]    # B^2
uparoverB(npint) _real [m/s*T]    # v_parallel/B
uparsq(npint)   _real [m**2/s**2] # vparallel^2
uparsq_new(npint) _real [m**2/s**2] # vparallel^2, udated v's
                                    # for implicit, eval at predicted x's
uperpsq(npint)  _real [m**2/s**2] # vperpendicular^2
usq(npint)      _real [m**2/s**2] # v^2
usq_new(npint)  _real [m**2/s**2] # v^2, updated v's
vbx(npint)      _real [m/s]       # x magnetic drift velocity
vby(npint)      _real [m/s]       # y magnetic drift velocity
vbz(npint)      _real [m/s]       # z magnetic drift velocity
vex(npint)      _real [m/s]       # x electric drift velocity
vey(npint)      _real [m/s]       # y electric drift velocity
vez(npint)      _real [m/s]       # z electric drift velocity
vdx(npint)      _real [m/s]       # V_dx
vdy(npint)      _real [m/s]       # V_dy
vdz(npint)      _real [m/s]       # V_dz
uxeff(npint)    _real [m/s]       # interpolated ux
uyeff(npint)    _real [m/s]       # interpolated uy
uzeff(npint)    _real [m/s]       # interpolated uz
alpha(npint)    _real             # interepolation parameter
alphabar(npint) _real             # complement of interp param

*********** W3Dcollisions:
langevincollisions3d(lself:logical,np:integer,
                     ux1(np):real,uy1(np):real,uz1(np):real,
                     vxmi(np):real,vymi(np):real,vzmi(np):real,
                     density2(np):real,vthsqi2(np):real,
                     q1:real,q2:real,m1:real,m2:real,
                     vthsqinit1:real,vthsqinit2:real,
                     dt:real,loglambda(np):real,epvth:real)
                         subroutine      # Langevin collision operator for like and
                                         # unlike particle species
compton_scatter_electrons(photon_energy:real,
                          photon_density(0:nx,0:ny,0:nz):real,
                          photon_px(0:nx,0:ny,0:nz):real,
                          photon_py(0:nx,0:ny,0:nz):real,
                          photon_pz(0:nx,0:ny,0:nz):real,
                          nx:integer, ny:integer, nz:integer, dx:real, dy:real, dz:real,
                          xmmin:real, ymmin:real, zmmin:real, zgrid0:real, dt:real,
                          boost_gamma:real, time:real, weightingmethod:integer,
                          anglethreshold:real,
                          electrongroup:ParticleGroup, electrons:integer,
                          photongroup:ParticleGroup, photons:integer,
                          tbirthpid:integer, scattercountpid:integer,
                          withgammaeffective:logical, scatterelectron:logical) subroutine
          # Scatters electrons off of a background photon density

*********** AMR dump:
AMRlevels                  integer /0/   # number of mesh refinement levels
                                         # (0 = no mesh refinement)
AMRrefinement              integer /2/   # refinement ratio between levels
AMRcoalescing(AMRlevels)  _real    /0.8/ # coefficient controlling coalescence
                                         # for each refinement level;
                                         # ranges between 0 and 1;
                                         # 1=minimal coalescence;
                                         # 0=most aggressive coalescence
AMRtransit                 integer /2/   # number of transition cells around
                                         # each patch, in units of parent cells.
AMRgenerate_periodicity    integer /1/   # periodicity at which to generate new
                                         # set of AMR blocks, in units of time
                                         # steps.
AMRmaxlevel_density        integer /-1/  # maximum level for automatic
                                         # refinement based in charge density
AMRmaxlevel_gradient       integer /-1/  # maximum level for automatic
                                         # refinement based in charge density
                                         # gradient
AMRthreshold_gradient      real    /0.8/ # threshold above which to refine for
                                         # automatic refinement based on
                                         # gradient
AMRmaxsize_isolated_blocks integer /0/   # maximum size of isolated blocks.
                                         # Blocks which contain
                                         # AMRmaxsize_isolated_blocks^dim
                                         # (where dim=2 in 2D and dim=3 in 3D)
                                         # cells or less are removed. The goal
                                         # is to prevent the creation of
                                         # refinement blocks due to statistical
                                         # noise.
                                         # THIS IS EXPERIMENTAL: this might be
                                         # partially implemented or not
                                         # implemented at all. There is no
                                         # guaranty at this stage that this
                                         # would effectively remove all
                                         # isolated blocks nor that it would
                                         # not remove isolated blocks which are
                                         # not due to statistical noise.
AMRuse_inactive_regions logical /.false./# When true, the fields in inactive
                                         # regions of refined patches will be
                                         # used. An inactive region is a region
                                         # that was not flagged to be refined,
                                         # but is refined because of its
                                         # proximity to other refined areas,
                                         # for example between two patches
                                         # that are coalesced. Normally, the
                                         # fields in these regions are not
                                         # used.

*********** W3D_interpsubs:
# Subroutines in file W3d_interp.F
#setptrs(bx:real,by:real,bz:real,ex:real,ey:real,ez:real)
#            subroutine # sets database B and E arrays
#                      # equal to bx, by, bz, ex, ey, ez
mugrdbpush(pgroup:ParticleGroup,np,is,ipmin,dtb:real,needcalcgradb)
            subroutine #  does the mu grad B parallel acceleration
                       #  and corresponding change to vperp
xpush3dintrp(pgroup:ParticleGroup,np,is,ipmin,dt:real)
              subroutine #  does the interpolated x push
getvperpparsq(pgroup:ParticleGroup,np,ipmin)
              subroutine #  finds v_perp^2, v_parallel^2, vparallel/B,v^2
                         #  of particles
getveff(pgroup:ParticleGroup,np,is,ipmin,x(np):real,y(np):real,z(np):real,
        ex(np):real,ey(np):real,ez(np):real,bx(np):real,by(np):real,bz(np):real,
        predcor:string,dt:real)
    subroutine #  gets components of interpolated velocity used in x push
getvdrift(np,is,x(np):real,y(np):real,z(np):real,gaminv(np):real,
          ex(np):real,ey(np):real,ez(np):real,bx(np):real,by(np):real,bz(np):real)
    subroutine #  calculates vdrift from ExB and gradB
setfields(pgroup:ParticleGroup,np,is,ipmin,x(np):real,y(np):real,z(np):real,
          ex(np):real,ey(np):real,ez(np):real,bx(np):real,by(np):real,bz(np):real,
          dt:real)
    subroutine #  sets E,B at particle arrays; determines interpolation
               #  parameter alpha and its complement alphabar
$setfields2(pgroup:ParticleGroup,np,is,x(np):real,y(np):real,z(np):real,
$           ux(np):real,uy(np):real,uz(np):real)
$    subroutine #  sets E,B at particle arrays; determines interpolation
$               #  parameter alpha and its complement alphabar
getgradbsq(np,is,x(np):real,y(np):real,z(np):real,
           bx(np):real,by(np):real,bz(np):real)
         subroutine    #  calculates or fetches grad B^2 components
getalphas(np,is,ipmin,x(np):real,y(np):real,z(np):real,dt:real)
          subroutine   # calculate interpolation parameter alpha
geteb(pgroup:ParticleGroup,np,is,ipmin,x(np):real,y(np):real,z(np):real,dt:real)
          subroutine   #  fetch E and B fields
set_polarization(rho(0:nx,0:nz),nx,nz,dx:real,dz:real,xmin:real,zmin:real)
          subroutine   # calculates x polarization and ratio of z to x polarization
                       # x-y geometry only, for now.
$setvdrifts(pgroup:ParticleGroup,np,is,x(np):real,y(np):real,z(np):real,
$           ux(np):real,uy(np):real,uz(np):real,predcor:string)
$    subroutine #  calculates vdrifts from ExB and gradB
oldsetup() subroutine # sets up pid indices for interpolation routines

*********** W3Dsubs:
# Subroutines in package 3D
w3dgen() subroutine
w3dexe() subroutine
w3dfin() subroutine
divxy(pgroup:ParticleGroup,iz,ndiv,
      divx(0:ndiv):real,divy(0:ndiv):real,divvx(0:ndiv):real,
      divvx2(0:ndiv):real,divvy(0:ndiv):real,
      divvy2(0:ndiv):real,wnpx(0:ndiv):real,wnpy(0:ndiv):real,itask)
             subroutine # calculates RMS vx and vy versus x and y
exteb3d(np:integer,xp(np):real,yp(np):real,zp(np):real,uzp(np):real,
        gaminv(np):real,dtl:real,dtr:real,
        bx(np):real,by(np):real,bz(np):real,
        ex(np):real,ey(np):real,ez(np):real,
        m:real,q:real,bendres(np):real,bendradi(np):real,dt:real)
             subroutine # Sets external E and B fields
othere3d(np:integer,xp(np):real,yp(np):real,zp(np):real,
         zbeam:real,zimax:real,zimin:real,
         straight:real,ifeears,eears:real,eearsofz(0:nzzarr):real,
         dzzi:real,nzzarr,
         zzmin:real,dedr:real,dexdx:real,deydy:real,dbdr:real,
         dbxdy:real,dbydx:real,
         ex(np):real,ey(np):real,ez(np):real,
         bx(np):real,by(np):real,bz(np):real)
             subroutine # Sets external E field
getese3d()   subroutine # Computes electrostatic energy
getese()     subroutine # Computes electrostatic energy
gtlchg3d()   subroutine # Computes line charge density
gtlchg()     subroutine # Computes line charge density
gtlchg3dfromrho(nxlocal:integer,nylocal:integer,nzlocal:integer,
                nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                rho(-nxguardrho:nxlocal+nxguardrho,
                    -nyguardrho:nylocal+nyguardrho,
                    -nzguardrho:nzlocal+nzguardrho):real,
                dx:real,dy:real,dz:real,
                zgrid:real,zmminlocal:real,
                l2symtry:logical,l4symtry:logical,islastproc:logical)
             subroutine # Computes line charge density from the rho array in 3D
gtlchgrzfromrho(nxlocal:integer,nzlocal:integer,
                nxguardrho:integer,nzguardrho:integer,
                rho:real,dx:real,dz:real,
                zgrid:real,zmminlocal:real,islastproc:logical)
             subroutine # Computes line charge density from the rho array in RZ
getese3dfromrhophi(nx:integer,ny:integer,nz:integer,
                   nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
                   nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                   rho:real,phi:real,
                   dx:real,dy:real,dz:real,
                   l4symtry:logical,l2symtry:logical,ese:real)
             subroutine # Computes electrostatic energy from the input arrays
geteserzfromrhophi(nx:integer,nz:integer,
                   nxguardphi:integer,nzguardphi:integer,
                   nxguardrho:integer,nzguardrho:integer,
                   rho:real,phi:real,
                   dx:real,dz:real,xmmin:real,ese:real)
             subroutine # Computes electrostatic energy from the input arrays
inject3d(itask:integer,pgroup:ParticleGroup)
             subroutine # Injection routine
injctint(pgroup:ParticleGroup) subroutine # Initialization for injection
fill_inj()   subroutine # Initializes arrays describing the geometry of the
                        # emitting surface. Automatically called.
inj_sete(pgroup:ParticleGroup,ipmin:integer,np:integer,
         ex(np):real,ey(np):real,ez(np):real)
             subroutine # Calculate the E field for particles near the
                        # emitting surface.
inj_sete3d(np:integer,xp(np):real,yp(np):real,zp(np):real,pid(np):real,
           ex(np):real,ey(np):real,ez(np):real)
             subroutine # Calculate the E field for particles near the
                        # emitting surface.
inj_transform(np:integer,x(np):real,y(np):real,z(np):real,
              ni:integer,ijp(ni):integer,
              tsign:integer,lshift:logical)
             subroutine # Transforms coordinates into and out of frame
                        # of injection sources
createparticlesincells(nx:integer,ny:integer,nz:integer,
                       rnn(0:nx,0:ny,0:nz):real,
                       exgrid(0:nx+1,0:ny,0:nz):real,
                       eygrid(0:nx,0:ny+1,0:nz):real,
                       ezgrid(0:nx,0:ny,0:nz+1):real,
                       condid(0:nx,0:ny,0:nz):real,
                       lcylindrical:logical,
                       dx:real,dy:real,dz:real,
                       nn:integer,xx(nn):real,yy(nn):real,zz(nn):real,
                       ex(nn):real,ey(nn):real,ez(nn):real,pp(nn):real)
             subroutine # Creates particles at each grid cell
loadrho3d(pgroup:ParticleGroup,ins:integer,nps:integer,is:integer,
          lzero:logical,lfinalize_rho:logical)
             subroutine # Provides a simple interface to the charge density
                        # loading routine setrho3d
finalizerho() subroutine # Finish the rho accumulation, taking care of parallel
                         # communication and boundary conditions
fetchphi(n:integer,x(n):real,y(n):real,z(n):real,p(n):real)
             subroutine # Fetches the electrostatic potential at the given
                        # list of locations. It uses whatever geometry and
                        # field solver that is active.
setupFields3dParticles()
             subroutine # Sets up the Fields3dParticles group
setrhoandphiforfieldsolve(rhop:real,phip:real)
             subroutine # Copies data from rhop to rho - mainly for parallel
setrhoforfieldsolve3d(nxlocal:integer,nylocal:integer,nzlocal:integer,
                      rho(-nxguardrho:nxlocal+nxguardrho,
                          -nyguardrho:nylocal+nyguardrho,
                          -nzguardrho:nzlocal+nzguardrho):real,
                      nxp:integer,nyp:integer,nzp:integer,
                      rhop(-nxguardrho:nxp+nxguardrho,
                           -nyguardrho:nyp+nyguardrho,
                           -nzguardrho:nzp+nzguardrho):real,
                      nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                      fsdecomp:Decomposition,ppdecomp:Decomposition)
             subroutine # Copies data from rhop to rho - for parallel
getphipforparticles(indts:integer)
             subroutine # Copies data from phi to phip - mainly for parallel
getphipforparticles3d(nc:integer,nxlocal:integer,nylocal:integer,
                      nzlocal:integer,
                      nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
                      phi:real,
                      nxp:integer,nyp:integer,nzp:integer,phip:real,
                      fsdecomp:Decomposition,ppdecomp:Decomposition)
             subroutine # Calls the parallel routine for copying phi into phip
allocateselfepforparticles(lforce:logical)
             subroutine # Allocates selfep if needed
allocateselfeforfieldsolve() subroutine # Allocates selfe if needed
padvnc3d(center:string,pgroup:ParticleGroup)
             subroutine # Advances particles and rho
perphi3d()
             subroutine # Equates end slices of phi for periodicity
setboundsfromflags(bounds(0:5):integer,boundxy:integer,bound0:integer,
                   boundnz:integer,
                   l2symtry:integer,l4symtry:integer) subroutine
applyrhoboundaryconditions()
             subroutine # Applies boundary conditions to rho
applyrhoboundaryconditions3d(nc:integer,
                     nxlocal:integer,nylocal:integer,nzlocal:integer,
                     nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                     rho:real,bounds(0:5):integer,fsdecomp:Decomposition,
                     lrz:logical)
             subroutine # Applies boundary conditions to 3-D rho
prntpa3d(lprntpara:logical)
             subroutine # Prints out 3d specific stuff (like prntpara())
bendez3d(np,xp(np):real,zp(np):real,ez(np):real,
         bendres(np):real,bendradi(np):real,
         bends:logical,bnezflag:logical,linbend:logical)
             subroutine #  Corrects axial electric field for warped geometry
epush3d(np,uxp(np):real,uyp(np):real,uzp(np):real,
        ex(np):real,ey(np):real,ez(np):real,q:real,m:real,dt:real)
             subroutine # Particle velocity advance from E field
epusht3d(np,uxp(np):real,uyp(np):real,uzp(np):real,
        ex(np):real,ey(np):real,ez(np):real,q:real,m:real,dt:real,fdt:real)
             subroutine # Particle velocity advance from E field with varying dt
bpush3d(np,uxp(np):real,uyp(np):real,uzp(np):real,gaminv(np):real,
        bx(np):real,by(np):real,bz(np):real,
        q:real,m:real,dt:real,ibpush:integer)
             subroutine # Particle velocity advance from B field
bpusht3d(np,uxp(np):real,uyp(np):real,uzp(np):real,gaminv(np):real,
         bx(np):real,by(np):real,bz(np):real,
         q:real,m:real,dtp(np):real,fdt:real,ibpush:integer)
             subroutine # Particle velocity advance from B field with varying dt
xpush3d(np,xp(np):real,yp(np):real,zp(np):real,
        uxp(np):real,uyp(np):real,uzp(np):real,gaminv(np):real,dt:real)
             subroutine # Particle position advance
xpusht3d(np,xp(np):real,yp(np):real,zp(np):real,
        uxp(np):real,uyp(np):real,uzp(np):real,gaminv(np):real,fdt:real)
             subroutine # Particle position advance with varying dt
ebcancelpush3d(np,uxp(np):real,uyp(np):real,uzp(np):real,gi(np):real,
                  exp(np):real,eyp(np):real,ezp(np):real,
                  bxp(np):real,byp(np):real,bzp(np):real,
                  q:real,m:real,dt:real,which:integer) subroutine # velocity push with E+vxB cancellation
ebcancelpush3dt(np,uxp(np):real,uyp(np):real,uzp(np):real,gi(np):real,
                  exp(np):real,eyp(np):real,ezp(np):real,
                  bxp(np):real,byp(np):real,bzp(np):real,
                  q:real,m:real,dt(np):real,which:integer) subroutine # velocity push with E+vxB cancellation
seteears()  subroutine # Sets eearsofz, the axial confining field
sete3d(phi1d:real,selfe:real,np,xp(np):real,yp(np):real,zp(np):real,zgrid:real,
       xmmin:real,ymmin:real,zmmin:real,dx:real,dy:real,dz:real,
       nx:integer,ny:integer,nz:integer,
       nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
       nxguarde:integer,nyguarde:integer,nzguarde:integer,
       efetch:integer,depos_order(0:2):integer,
       ex(np):real,ey(np):real,ez(np):real,
       l2symtry:logical,l4symtry:logical,lcylindrical:logical)
             subroutine # Sets internal E field
getselfe3d(phi(-nxguardphi:nx+nxguardphi,
               -nyguardphi:ny+nyguardphi,
               -nzguardphi:nz+nzguardphi):real,
           nx:integer,ny:integer,nz:integer,
           nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
           selfe(0:2,-nxguarde:nx+nxguarde,
                     -nyguarde:ny+nyguarde,
                     -nzguarde:nz+nzguarde):real,
           nxguarde:integer,nyguarde:integer,nzguarde:integer,
           dx:real,dy:real,dz:real,
           lzero:logical)
             subroutine # Calculates the self-E via finite difference of phi
setrho3d(rho(-nxguardrho:nx+nxguardrho,
             -nyguardrho:ny+nyguardrho,
             -nzguardrho:nz+nzguardrho):real,
         np:integer,xp(np):real,yp(np):real,zp(np):real,zgrid:real,
         q:real,wght:real,depos:string,depos_order(0:2):integer,
         nx:integer,ny:integer,nz:integer,
         nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
         dx:real,dy:real,dz:real,xmmin:real,ymmin:real,zmmin:real,
         l2symtry:logical,l4symtry:logical,lcylindrical:logical)
             subroutine # Computes charge density on a 3D grid.
                        # Calls one of several different routines based on the
                        # value of depos.
setrho3dw(rho(-nxguardrho:nx+nxguardrho,
              -nyguardrho:ny+nyguardrho,
              -nzguardrho:nz+nzguardrho):real,
          np:integer,xp(np):real,yp(np):real,zp(np):real,zgrid:real,
          wfact(np):real,
          q:real,wght:real,depos:string,depos_order(0:2):integer,
          nx:integer,ny:integer,nz:integer,
          nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
          dx:real,dy:real,dz:real,xmmin:real,ymmin:real,zmmin:real,
          l2symtry:logical,l4symtry:logical,lcylindrical:logical)
             subroutine # Computes charge density on a 3D grid
                        # using variable weighted particles
setrho3ddirect2(rho(-nxguardrho:nx+nxguardrho,
                    -nyguardrho:ny+nyguardrho,
                    -nzguardrho:nz+nzguardrho):real,
                np:integer,
                xp(np):real,yp(np):real,zp(np):real,zgrid:real,q:real,wght:real,
                nx:integer,ny:integer,nz:integer,
                nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                dx:real,dy:real,dz:real,xmmin:real,ymmin:real,zmmin:real,
                l2symtry:logical,l4symtry:logical,lcylindrical:logical)
             subroutine # Computes charge density on a 2D or 3D grid, with
                        # the option of cylindrical coordinates.
sezax3d()    subroutine # Sets EZAX, Ez on axix
sezax()      subroutine # Sets EZAX, Ez on axix
sphiax3d()   subroutine # Sets PHIAX, E. S. potential on axis
sphiax()     subroutine # Sets PHIAX, E. S. potential on axis
srhoax3d()   subroutine # Sets RHOAX, charge density on axis
srhoax()     subroutine # Sets RHOAX, charge density on axis
rhodia3d()   subroutine # Sets rhomid and rhomax diagnostics
rhodia()     subroutine # Sets rhomid and rhomax diagnostics
stckxy3d(pgroup:ParticleGroup,js:integer,zbeam:real,lcountaslost:logical)
             subroutine # Enforces transverse absorbing boundaries
stptcl3d(pgroup:ParticleGroup)   subroutine # Particle initializer
setrstar(rstar(-1:nz+1):real,nz:integer,dz:real,zmmin:real,zgrid:real)
             subroutine # Loads radius of reference orbit into rstar array
fieldsol3d(iwhich) subroutine # Bent-self-field iterative solver
vp3d(iwhich) subroutine # The 3d Poisson solver

pltfld3d(fld:string,freqflag:integer)
             subroutine # Controls field plots

multpole(lmod:integer,nlmod:integer,irpowmx:integer,
         lcosex:logical,lsinex:logical,aper:real,xcen:real,ycen:real,
         nmult:integer,nres:integer,tol:real)
            subroutine # calculate the multipole moments of the potential
inj_smoother(nx:integer,ny:integer,inj_phi(0:nx,0:ny):real,
             dx:real,dy:real,xmmin:real,ymmin:real,
             x0:real,y0:real,a0:real,b0:real,
             inj_nsmooth:integer) subroutine
getinj_phi() subroutine
getinj_phi_3d() subroutine
gettinj_phi() subroutine
fetche3d(pgroup:ParticleGroup,ipmin:integer,ip:integer,is:integer) subroutine
fetchb3d(pgroup:ParticleGroup,ipmin:integer,ip:integer,is:integer) subroutine
fetche3dfrompositions(jsid:integer,indts:integer,n:integer,
                      x(n):real,y(n):real,z(n):real,
                      ex(n):real,ey(n):real,ez(n):real,
                      bx(n):real,by(n):real,bz(n):real) subroutine
particlegridboundaries3d(pgroup:ParticleGroup,js:integer) subroutine
particleboundaries3d(pgroup:ParticleGroup,js:integer,lcallcontrollers:logical) subroutine
particleboundariesxy(pgroup:ParticleGroup,js:integer,lcallcontrollers:logical) subroutine
loadperpdist0(kioff:integer,np:integer,
              x(np):real,y(np):real,xp(np):real,yp(np):real,
              rx(np):real,ry(np):real,rxp(np):real,ryp(np):real,
              epsx(np):real,epsy(np):real) subroutine
loadperpdist(kioff:integer,np:integer,
             x(np):real,y(np):real,r(np):real,t(np):real,
             xp(np):real,yp(np):real,
             rx(np):real,ry(np):real,rxp(np):real,ryp(np):real,
             epsx(np):real,epsy(np):real) subroutine
check_cc3d(pgroup:ParticleGroup,is:integer,ipmin:integer,np:integer) subroutine
set_aperture_e() subroutine
setuppadvncsubcyclingaveraging(it:integer,center:string,pgroup:ParticleGroup) subroutine
ijcoll2d(i1:integer,i2:integer,np1:integer,np2:integer,
         nxp:integer,nzp:integer,ux1:real,ux2:real,uy1:real,uy2:real,
         uz1:real,uz2:real,depvth:real,x1:real,z1:real,x2:real,z2:real,
         lx:real,lz:real,aion1:real,zion1:real,aion2:real,zion2:real,
         nion1:real,nion2:real,vthmisq1:real,vthmisq2:real,dnudt:real,
         dncint:integer,ixbc:integer,izbc:integer)
         subroutine     # Unlike particle collision routine

pdiag2d(x2:real,y2:real,ux2:real,uy2:real,uz2:real,pdmy2:real,np2:integer,
        lx:real,ly:real,nxp:integer,nyp:integer,s0:real,
        svx:real,svy:real,svz:real,se:real,ixbc:integer,izbc:integer)
         subroutine     # Computes grid deposition of moments
setupgridextent() subroutine # Sets up grid cell sizes and mesh arrays

*********** W3Dutilities:
sortparticlesbyindex1(n:integer,indx(n):integer,x(n):real,y(n):real,z(n):real,
                      uz(n):real,nw:integer,wfact(nw):real,nblocks:integer,
                      xout(n):real,yout(n):real,zout(n):real,uzout(n):real,
                      wfactout(nw):real,pcounts(0:nblocks-1):integer)
      subroutine
sortparticlesbyindex2(n:integer,indx(n):integer,x(n):real,y(n):real,z(n):real,
                      ux(n):real,uy(n):real,uz(n):real,gaminv(n):real,
                      nw:integer,wght(nw):real,
                      nblocks:integer,
                      xout(n):real,yout(n):real,zout(n):real,
                      uxout(n):real,uyout(n):real,uzout(n):real,
                      gaminvout(n):real,wghtout(nw):real,
                      pcounts(0:nblocks-1):integer)
      subroutine
sortparticlesbyindexgetisort(n:integer,indx(n):integer,
                             x(n):real,y(n):real,z(n):real,
                             nblocks:integer,
                             xout(n):real,yout(n):real,zout(n):real,
                             isort(n):integer,pcounts(0:nblocks-1):integer)
      subroutine
getichild(gridnumb:integer,np:integer,x(np):real,y(np):real,z(np):real,
          ichild(np):integer,
          nx:integer,ny:integer,nz:integer,grid(0:nx,0:ny,0:nz):integer,
          xmin:real,xmax:real,ymin:real,ymax:real,
          zmin:real,zmax:real,
          zgrid:real,l2symtry:logical,l4symtry:logical)
      subroutine # Gathers data from a 3-D grid using nearest grid point.
getichildpositiveonly(gridnumb:integer,np:integer,
                      x(np):real,y(np):real,z(np):real,
                      ichild(np):integer,
                      nx:integer,ny:integer,nz:integer,
                      grid(0:nx,0:ny,0:nz):integer,
                      xmin:real,xmax:real,ymin:real,ymax:real,
                      zmin:real,zmax:real,
                      zgrid:real,l2symtry:logical,l4symtry:logical)
      subroutine # Gathers data from a 3-D grid using nearest grid point.
                 # Only gathers positive values.
smooth121nonzero(uin:real,uout:real,nx:integer,ny:integer,nz:integer)
      subroutine # Perform a 1/4-1/2-1/4 smoothing on only the nonzero data
smooth121nonzerorz(uin:real,uout:real,nx:integer,nz:integer)
      subroutine # Perform a 1/4-1/2-1/4 smoothing on only the nonzero data for RZ array
getabsgrad3d(nx:integer,ny:integer,nz:integer,
             f(0:nx,0:ny,0:nz):real,gr(0:nx,0:ny,0:nz):real,
             dx:real,dy:real,dz:real)
      subroutine
addsortedefield(n:integer,isort(0:n-1):integer,
                tex(0:n-1):real,tey(0:n-1):real,tez(0:n-1):real,
                ex(0:n-1):real,ey(0:n-1):real,ez(0:n-1):real)
      subroutine
addsortedpotential(n:integer,isort(0:n-1):integer,
                tpotential(0:n-1):real,
                potential(0:n-1):real)
      subroutine
getextpart(pgroup:ParticleGroup)  subroutine
setupgrid3dtype(grid:Grid3dtype,check:logical) subroutine
      # Checks the consistency of Grid3dtype input and allocates the grid
      # If check is false, then the input is inconsistent.
setupgrid2dtype(grid:Grid2dtype,check:logical) subroutine
      # Checks the consistency of Grid3dtype input and allocates the grid
      # If check is false, then the input is inconsistent.
setupdecompositionw3d() subroutine
      # Setup the domain decomposition
domaindecomposefields(nz:integer,nslaves:integer,lfsautodecomp:logical,
        izfsslave(0:nslaves-1):integer,nzfsslave(0:nslaves-1):integer,
        overlap:integer) subroutine
      # Do the domain decomposition for the field solver
domaindecomposeparticles(nz:integer,nslaves:integer,
                nzguardrho:integer,zmmin:real,
                dz:real,zslave(0:nslaves-1):real,lautodecomp:logical,
                izpslave(0:nslaves-1):integer,nzpslave(0:nslaves-1):integer,
                zpslmin(0:nslaves-1):real,zpslmax(0:nslaves-1):real) subroutine
      # Do the domain decomposition for the particles
initializedecomp(decomp:Decomposition) subroutine
      # Does some setup for a decompostion object
flagparticlesbysortedssn(flag:real,np:integer,ssn(np):real,pflag(np):real,
                         nflagged:integer,ssnflagged(nflagged):real) subroutine
      # Flag all of the particles listed in the sorted ssnflagged array
transferarray1toarray23d(nx1:integer,ny1:integer,nz1:integer,
                         array1(-nxguard:nx1+nxguard,
                                -nyguard:ny1+nyguard,
                                -nzguard:nz1+nzguard):real,
                         nx2:integer,ny2:integer,nz2:integer,
                         array2(-nxguard:nx2+nxguard,
                                -nyguard:ny2+nyguard,
                                -nzguard:nz2+nzguard):real,
                         nxguard:integer,nyguard:integer,nzguard:integer,
                         lsendguards:logical,
                         decomp1:Decomposition,decomp2:Decomposition) subroutine
      # Transfers data from one decomposition to another


*********** W3Dload:
r_b      real [m]     /0./ # Code set: CFE rms equivalent beam radius
emit_b   real [m-rad] /0./ # Code set: CFE rms equivalent beam emittance
q_perv   real [1]     /0./ # Code set: CFE rms equivalent beam perveance
k_beta0  real [m^-1]  /0./ # Code set: CFE rms equivalent beam foc. wavenumber
k1re     real [1]     /0./ # Code set: CFE WB distribution k_1*r_e
f_2      real [m^-2]  /0./ # Code set: CFE WB distribution norm f_2
r_e      real [m]     /0./ # Code set: CFE WB distribution edge radius
delta    real [1]     /0./ # Code set: CFE TE distribution delta parameter
logdeltarad real [1] /3.0/   # Bracket radius in log(delta) for delta sol
                             #   in TE root find: default: logdeltarad = 3.0
logdeltatol real [1] /1.e-4/ # Tolerance to calculate log(delta) for delta sol
                             #   in TE root find: default: logdeltarad = 1.e-4
glambdad real [1]     /0./ # Code set: CFE TE distribution gamma*lambda_Debye
tscaled  real [1]     /0./ # Code set: CFE TE distribution scaled temperature
den0     real [m^-3]  /0./ # Code set: CFE TE distribuiton on-axis density
hl_rmax   real  [m] # maximum radius of continuous focusing (CF) equivalent
                    #  Hamiltonian equilibrium density and angle spreads
hl_nrdist integer   # number of radial data points of CF Hamiltonian eq
                    #  ranging from [0:hl_nrdist] -> [0,hl_rmax].
hl_rdist(0:hl_nrdist)    _real
   # Code set: particle density in r calculated for CF Hamiltonian eq load
hl_nrmrdist(0:hl_nrdist) _real
   # Code set: normalized hl_rdist with sum(hl_nrmrdist)=1.
hl_intrdist(0:hl_nrdist) _real
   # Code set: integral of r*rdist from r=0 to r where r is a scaled radius.
perp_cfe_den() subroutine # Calculates the CFE density profile
bessi0(x:real) real function # Modified Bessel function I_0(x) for real x
bessi1(x:real) real function # Modified Bessel function I_1(x) for real x
bessi(n:integer,x:real) real function # Modified Bessel function I_n(x)
                                      #  for real x and integer n >= 2

*********** W3DloadTE:
te_rhomax    real  [1] # Code set: max value of norm radius rho to use for TE
te_rhotrans  real  [1] # Code set: transition value of rho in norm density calc
te_denptrans real  [1] /1.e-3/
                       # Value of d den/d rho at transition in TE density calc
te_rhodenint real  [1] /0./ # Code set: integral of rho*den(rho) on rho grid
te_nrho integer       # number of radial grid points in rho to use for TE
te_rho(0:te_nrho)    _real
   # Code set: normalized radial coordinate grid ranging from 0 to te_rhomax
te_den(2,0:te_nrho)  _real
   # Code set: TE radial density profile on te_rho grid:
   #   te_den(1,:) = normalized density
   #   te_den(2,:) = derivative of normalized density wrt rho
te_exp_nterm integer  # max number of terms to use in series expansion for
                      #   TE density.  Should set this to at least 100.
te_exp_alpha(0:te_exp_nterm) _real
   # Code set: Coefficients of series expansion for TE density profile
te_den_exp_coeff(delta:real,nterm:integer,alpha(0:nterm):real) subroutine
   # Calculate expansion coefficients for CF TE density profile
te_den_exp(rho:real,delta:real,nterm:integer,alpha(0:nterm):real,tol:real)
   real function
   # Calculate the density of a CF TE using a series expansion
te_denp_exp(rho:real,delta:real,nterm:integer,alpha(0:nterm):real,tol:real)
   real function
   # Calculate the derivative with respect of rho of the density of a
   #   CF TE using a series expansion
te_rho_max(delta:real,frac:real) real function
   # Estimate max value of norm radius rho where norm density = frac for CF TE
te_delta_est(sc_param:real) real function
   # Estimate the value of delta needed for an rms equiv TE
te_radial_den(delta:real) subroutine
   # Calculate the normalized radial density profile of TE for dimensionless
   #   space-charge parameter delta
integrate_test(y(nvar,nstep):real,nvar:integer,nstep:integer,x1:real,x2:real)
   subroutine
   # temp routine link for debugging
te_constr_test(sc_param:real,delta:real) real function
   # temp routine link for debugging
te_delta_root_test(sc_param:real,delta_1:real,delta_2:real,tol:real)
  real function
   # temp routine link for debugging
smooth3d_121(q(0:nx,0:ny,0:nz):real,nx:integer,ny:integer,nz:integer,
               npass(3):integer,alpha(3):real) subroutine
   # three points linear smoothing
smooth3d_121withbounds(q(0:nx,0:ny,0:nz):real,nx:integer,ny:integer,nz:integer,
               npass(3):integer,alpha(3):real,bounds(0:5):integer) subroutine
   # three points linear smoothing, taking into account boundary conditions
smooth3d_121_stride(q(0:nx,0:ny,0:nz):real,nx:integer,ny:integer,nz:integer,
                    npass(3):integer,alpha(3):real,stride(3):integer) subroutine
   # three points linear smoothing with stride
smoothcirc_121_stride(q(0:nx,0:nz,0:nm):complex,nx:integer,nz:integer,
                    nm:integer,npass(3):integer,alpha(3):real,stride(3):integer) subroutine
   # three points linear smoothing with stride
smooth3d_121_stride_mask(q(0:nx,0:ny,0:nz):real,mask(0:nx,0:ny,0:nz):real,
                    nx:integer,ny:integer,nz:integer,
                    npass(3):integer,alpha(3):real,stride(3):integer) subroutine
   # three points linear smoothing with stride
smooth_poly(q(0:nx,0:ny,0:nz):real,nx:integer,ny:integer,nz:integer,sten(0:4):real) subroutine
   # performs smoothing in z with nine points stencil 
   
******** Subtimersw3d dump:
lw3dtimesubs logical /.false./
timew3dinit real /0./
timew3dvers real /0./
timew3dgen real /0./
timew3dexe real /0./
timew3dfin real /0./
timestep3d real /0./
timeexteb3d real /0./
timeothere3d real /0./
timegetese3d real /0./
timerhodia real /0./
timegtlchg real /0./
timesrhoax real /0./
timegetese real /0./
timesphiax real /0./
timesezax real /0./
timegtlchg3d real /0./
timegtlchg3dfromrho real /0./
timegtlchgrzfromrho real /0./
timegetese3dfromrhophi real /0./
timegeteserzfromrhophi real /0./
timeseteears real /0./
timepositionadvance3d real /0./
timepadvnc3d real /0./
timegetextpart real /0./
timeperphi3d real /0./
timeperb3d real /0./
timeapplyrhoboundaryconditions real /0./
timeepush real /0./
timeepush3d real /0./
timeepusht3d real /0./
timebpush real /0./
timebpush3d real /0./
timebpusht3d real /0./
timexpush3d real /0./
timexpusht3d real /0./
timeebcancelpush3d real /0./
timecheck_cc3d real /0./
timesete3d_relativity real /0./
timeedamp real /0./
timegetbend real /0./
timebendez3d real /0./
timesete3d real /0./
timesete2d real /0./
timegetselfe3d real /0./
timestptcl3d real /0./
timesetrho3d real /0./
timeloadrho3d real /0./
timefinalizerho real /0./
timesetrhoandphiforfieldsolve real /0./
timesetrhoforfieldsolve3d real /0./
timefetche3d real /0./
timefetche3dfrompositions real /0./
timefetchb3d real /0./
timeparticlegridboundaries3d real /0./
timeparticleboundaries3d real /0./
timeparticleboundariesxy real /0./
timestckxy3d real /0./
timesetrstar real /0./
timefetchphi real /0./
timeinject3d real /0./
timeinjctint real /0./
timefill_inj real /0./
timeinj_transform real /0./
timegetinj_phi real /0./
timegettinj_phi real /0./
timegetinj_phi_mr real /0./
timeinj_sete real /0./
timeinj_sete3d real /0./
timeinj_sete_3darray real /0./
timetinj_sete3d real /0./
timeinj_addtemp3d real /0./
timeinj_setrho3d real /0./
timeinj_setrhomr real /0./
timesete3d_aperture real /0./
timeset_aperture_e real /0./
timegetphipforparticles real /0./
timegetphiforfields real /0./
timefieldsol3d real /0./
timevp3d real /0./
timesortpart_byindex real /0./
timesortpart_byindexgetisort real /0./
timegetichild real /0./
timegetichildpositiveonly real /0./
timeaddrhotoowner real /0./
timegetrhofromowner real /0./
timegetabsgrad real /0./
timeaddsortedefield real /0./
timeaddsortedpotential real /0./

timeinit_w3d_parallel real /0./
timesw_globalsum real /0./
timesumsourcepondomainboundaries real /0./
timeapplyrhoboundaryconditions3d real /0./
timesetsourceforfieldsolve3d_parallel real /0./
timemakesourceperiodic_slave_work real /0./
timeperpot3d_slave real /0./
timegetphipforparticles3d_parallel real /0./
timegetphiforparticles3d real /0./
timegetphiforfields3d real /0./
timegetsourcetransfermessagedata real /0./
timetransfersourceptosource3d real /0./
timetransferarray1toarray23d real /0./

timerdata(10) real /10*0./

timemugrdbpush real /0./
timexpush3dintrp real /0./
timegetvperpparsq real /0./
timegetveff real /0./
timegetvdrift real /0./
timesetfields real /0./
timegetgradbsq real /0./
timegetalphas real /0./
timegeteb real /0./
timefetche3d1 real /0./
timeoldsetup real /0./
timeswapxarrays real /0./
timestoreoldu real /0./
timeswapxarr2 real /0./
timeset_polarization real /0./
timesetTotalE real /0./

