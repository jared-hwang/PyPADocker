fxy
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package FXY of code WARPxy
# fieldsolver package and test driver
# David Grote, LLNL, (510)423-7194
{
}

******** CapMatxy dump: 
ncxymax            integer  # Maximum number of points in conductors
ncxy               integer  # Number of points within conductors
xcond(ncxymax)    _real [m] # X coordinate of points in conductors
ycond(ncxymax)    _real [m] # Y coordinate of points in conductors
vcond(ncxymax)    _real [V] # voltage of points in conductors
pcond(ncxymax)    _real [V] # actual potential on points in conductors
                            # (auto set)
qcond(ncxymax)    _real [C] # induced charge on points in conductors (auto set)
cmatxy(ncxy,ncxy) _real /0./ # Capacity matrix (auto set)
kpvtxy(ncxy)      _integer  # Pivot points for matrix solve (auto set)

*********** FXYsubs:
capmatxyf(iwhich:integer,phi:real,kxsq:real,kysq:real,attx:real,atty:real,
          filt:real,xlen:real,ylen:real,nx:integer,ny:integer,dx:real,dy:real,
          xmmin:real,ymmin:real,scrtch:real,phisave:real,xywork:real,
          l2symtry:logical,l4symtry:logical)
        subroutine
