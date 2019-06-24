wxy
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package WXY of code WARP
# XY - PIC package of xy particle code
# Alex Friedman, LLNL, (510)422-0827
# David P. Grote, LLNL, (510)423-7194

*********** Particlesxy dump parallel:
dtpid      integer # ID in the pid array where the time step size for each
                   # particle is saved.

*********** InGenxy dump:
ds       real /0./         # Axial step size, defaults to vbeam*dt
lvzchang logical /.true./  # When true, iterative algorithm is used to
                           # estimate dt.
lexactds logical /.false./
niter_dt integer /4/       # Number of iterations in the calculation of dt
lthick   logical /.false./ # Sets whether to do thick slice model.
                           # When false, does thin slice model.
lvp3d    logical /.false./ # Sets whether or not to use the 3-D field solver
ldiag    logical /.true./  # When false, no diagnostics are done
lexbend  logical /.true./  # When true, use exact transformation in bend, when
                           # false, use same transformation as in 3-D code.
lcommonz logical /.false./ # When true, the z of all particles is set to zbeam
                           # when they are initially created. Otherwise, the
                           # z's have an artificial spread.
lwithez  logical /.false./ # When true, uses an approximate Ez calculated
                           # Using the potential from the previous step.
lzstepcorrection logical /.false./

*********** Fieldsxy:
xywork(:,:) _real # Work array for the 2D FFT solver
xyphisave(:,:) _real # Work array for the 2D FFT solver

*********** WXYsubs:
# Subroutines in package XY
wxygen() subroutine
wxyexe() subroutine
initdtp(pgroup:ParticleGroup) subroutine # Initializes dtp
initzpxy(pgroup:ParticleGroup) subroutine # Optionally sets all z values to zbeam.
extebxy(np,xp:real,yp:real,zp:real,uzp:real,gaminv:real,dtl:real,dtr:real,
        bx:real,by:real,bz:real,ex:real,ey:real,ez:real,
        m:real,q:real,bendres:real,bendradi:real,lexbend:logical,gammabar:real,
        zbeam:real,vbeam:real,dt:real,time:real)
             subroutine # Sets external E and B fields
otherexy(np,xp:real,yp:real,dedr:real,dexdx:real,deydy:real,dbdr:real,
         ex:real,ey:real,ez:real,
         bx:real,by:real,bz:real)
             subroutine # Sets external E field
padvncxy(center:string)
             subroutine # Advances particles and rho
fixrhoxy(rho:real,nx,ny,nz,periinz:logical,lthick:logical)
             subroutine # Sums end slices of rho for periodicity
bendezxy(np,xp:real,zp:real,ez:real,bendres:real,bendradi:real,
          bends:logical,bnezflag:logical,linbend:logical)
             subroutine # Corrects axial electric field for warped geometry
epushxy(np,uxp:real,uyp:real,uzp:real,ex:real,ey:real,ez:real,q:real,m:real,
        dtp:real,fdt:real)
             subroutine # Particle velocity advance from E field
bpushxy(np,uxp:real,uyp:real,uzp:real,gaminv:real,bx:real,by:real,bz:real,
        q:real,m:real,dtp:real,fdt:real,ibpush:integer)
             subroutine # Particle velocity advance from B field
xpushxy(np,xp:real,yp:real,zp:real,uxp:real,uyp:real,uzp:real,gaminv:real,
        dtp:real)
             subroutine # Particle position advance
setcurrxy(curr:real,np,zp:real,uzp:real,gaminv:real,q:real,wght:real,
          zbeam:real,dzz:real,zzmin:real,dz:real,lthick:logical)
             subroutine # Computes current
setrhoxy(rho1d:real,np:integer,xp:real,yp:real,zp:real,zgrid:real,
         uzp:real,gaminv:real,q:real,wght:real)
             subroutine # Computes charge density
fieldsolxy(iwhich:integer)
             subroutine # Complete field solve
vpxy(iwhich:integer)
             subroutine # Call field solver
loadrhoxy(pgroup:ParticleGroup,ins:integer,nps:integer,is:integer,
          lzero:logical,lfinalize_rho:logical)
             subroutine # Simple interface to setrhoxy
fetchexy(pgroup:ParticleGroup,ipmin:integer,ip:integer,is:integer,
         ex:real,ey:real,ez:real)
             subroutine # Returns electric field on particles from phi

*********** Subtimerswxy:
lwxytimesubs logical /.false./
timewxyinit real /0./
timewxyvers real /0./
timewxygen real /0./
timewxyexe real /0./
timewxyfin real /0./
timestepxy real /0./
timeextebxy real /0./
timeotherexy real /0./
timepadvncxy real /0./
timefixrhoxy real /0./
timefixcurrxy real /0./
timeepushxy real /0./
timebpushxy real /0./
timexpushxy real /0./
timeinitdtp real /0./
timesetdtp real /0./
timegetnewdtpwithe real /0./
timenextbend real /0./
timebendezxy real /0./
timeexbendcorxy real /0./
timebendcorxy real /0./
timesetrhoxy real /0./
timesetexy real /0./
timeloadrhoxy real /0./
timesetcurrxy real /0./
timebendfieldsolxy real /0./
timefieldsolxy real /0./
timevpxy real /0./
