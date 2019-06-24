"""
 _______________________________________________________________________________

 *** Copyright Notice ***

 "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
 2016, The Regents of the University of California, through Lawrence Berkeley
 National Laboratory (subject to receipt of any required approvals from the
 U.S. Dept. of Energy). All rights reserved.

 If you have questions about your rights to use or distribute this software,
 please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

 NOTICE.
 This Software was developed under funding from the U.S. Department of Energy
 and the U.S. Government consequently retains certain rights. As such, the U.S.
 Government has been granted for itself and others acting on its behalf a
 paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 reproduce, distribute copies to the public, prepare derivative works, and
 perform publicly and display publicly, and to permit other to do so.


 Class EM3DPXR for using the PICSAR subroutines in the WARP PIC loop.

 Developers:
 Henri Vincenti
 Remi Lehe
 Jean-Luc Vay
 Mathieu Lobet
 Guillaume Blaclard

 Date:
 Creation 2016

 _______________________________________________________________________________
"""

from warp.field_solvers.em3dsolverFFT import *
from warp.particles.species import *
from .laser.laser_antenna import LaserAntenna

try:
    from mpi4py import MPI
    print 'from mpi4py import MPI'
except:
    print 'Error cannot import mpi4py'

try:
    #import warp.field_solvers.GPSTD as gpstd
    import GPSTDPXR as gpstd
    print 'Import GPSTDPXR as gpstd'
except:
    #import GPSTDPXR as gpstd
    import warp.field_solvers.GPSTD as gpstd

try:
    from picsar_python import picsarpy as pxrpy
    print 'Import picsarpy as pxrpy'
    pxr = pxrpy.picsar
    l_pxr=True
except:
    l_pxr=False

try:
    import numpy as numpy
except:
    print 'Error cannot import numpy'

try:
    import os as os
    print 'Import os as os'
except:
    print 'Error cannot import os'

try:
    # Try to import fortran wrapper of FFTW
    # import pyfftw
    # fft = pyfftw.interfaces.numpy_fft
    import fastfftforpy as fftpy
    print 'Import fastfftforpy as fftpy'
    import fastfftpy as fstpy
    print 'Import fastfftpy as fstpy'
    fst=fstpy.fastfft
    fft=fftpy
    l_fftw=True
except:
    fft = np.fft
    l_fftw=False

def addparticlesPXR(self,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.,gi=1.,w=None,
                     ux=None,uy=None,uz=None,
                     lallindomain=False,
                     xmmin=None,xmmax=None,
                     ymmin=None,ymmax=None,
                     zmmin=None,zmmax=None,
                     lmomentum=False,
                     l2symtry=None,l4symtry=None,lrz=None,
                     pidpairs=None,
                     lnewparticles=True):

    """
    Adds particles to the simulation (this subroutine replaces Warp's
    particles.addparticles function, reproducing most basic functionalities).

    Parameters:
    -----------
    - x,y,z,vx,vy,vz,gi: particle coordinates and velocities.
            Can be arrays or scalars. Scalars are broadcast to all particles.
            Any that are unsupplied default to zero, except gi, which defaults
            to 1. (gi is 1/gamma, the relatistic parameter)

    - w=1.: particle weight
            this is only used if top.wpid > 0 and if lnewparticles is true.

    - ux, uy, uz: particle momentum instead of vx,vy,vz
            These can be used, passing in momentums. The flag lmomentum will be
            set to true.

    - lallindomain=false:
            Flags whether particles are within the parallel domains.
            When true, all particles are assumed to be with in the extent of
            the domain so particle scraping is not done. This is automatically
            set to true when the code is in slice mode, i.e. after
            package('wxy'). Except if the option is explicitly set. If false,
            the particles that are outside the parallel domain are not added.

    - xmmin=top.xpminlocal, xmmax=top.xpmaxlocal
            x extent of the domain, should only be set in unusual circumstances.

    - ymmin=top.ypminlocal, ymmax=top.ypmaxlocal
            y extent of the domain, should only be set in unusual circumstances.

    - zmmin=top.zpminlocal+top.zgrid, zmmax=top.zpmaxlocal+top.zgrid
            z extent of the domain, should only be set in unusual circumstances.

    - lmomentum=false:
            Flags whether momentum or velocities are being input.
            Set to false when velocities are input as velocities, true when
            input as massless momentum (as WARP stores them).
            Only used when top.lrelativ is true.

    - l2symtry, l4symtry, lrz
            System symmetries default to w3d values

    - pidpairs=None
            Allows setting specific pid columns. Argument must be a list of
            lists, each having the format [id,pidvalue].
            id is the one-based index returned by nextpid.
            The assigment pid[:,id-1] = pidvalue is done.

    - lnewparticles=true:
            Flag whether the particles are treated as new. If so, the ssn will be
            set if needed, and the position saved as the birth location.
            Set this to false if using addparticles to restore particles.

    """

    nps0 = x.size
    pids = np.zeros([nps0,pxr.npid])

    if pidpairs is not None :
        for id,pp in pidpairs:
            pids[:,id-1] = pp

    if top.wpid >0:
        if w is None:
            w=np.ones(nps0)
        pids[:,pxr.wpid-1]=w*self.sw0

    # --- Use momentum quantities if specified. These take
    # --- precedence over vx, vy, and vz.
    if ux is not None or uy is not None or uz is not None:
        lmomentum = true
        if ux is not None: vx = ux
        if uy is not None: vy = uy
        if uz is not None: vz = uz

    # --- Convert all to arrays of length maxlen, broadcasting scalars
    x = array(x)*ones(nps0,'d')
    y = array(y)*ones(nps0,'d')
    z = array(z)*ones(nps0,'d')
    gi = array(gi)*ones(nps0,'d')

    if lmomentum:
        vx = array(vx)*ones(nps0,'d')
        vy = array(vy)*ones(nps0,'d')
        vz = array(vz)*ones(nps0,'d')
    else:
        vx = array(vx)*ones(nps0,'d')/gi
        vy = array(vy)*ones(nps0,'d')/gi
        vz = array(vz)*ones(nps0,'d')/gi

    # --- Set extent of domain
    if xmmin is None: xmmin = top.xpminlocal
    if xmmax is None: xmmax = top.xpmaxlocal

    if w3d.solvergeom == w3d.XZgeom:
        if ymmin is None: ymmin = -largepos
        if ymmax is None: ymmax = +largepos
    else:
        if ymmin is None: ymmin = top.ypminlocal
        if ymmax is None: ymmax = top.ypmaxlocal

    if zmmin is None: zmmin = top.zpminlocal + top.zgrid
    if zmmax is None: zmmax = top.zpmaxlocal + top.zgrid

    if l2symtry is None: l2symtry = w3d.l2symtry
    if l4symtry is None: l4symtry = w3d.l4symtry
    if lrz is None: lrz = (w3d.solvergeom in [w3d.RZgeom,w3d.Rgeom])

    # --- Do some error checking
    if not lallindomain:
        if not lrz and ymmin == ymmax:
            print "========================================================="
            print " AddparticlesPXR: warning - no particles will be loaded. "
            print " You should either set lallindomain=true or set ymmin and"
            print " ymmax so they are different from each other.            "
            print "========================================================="

        if zmmin == zmmax:
            print "========================================================="
            print " AddparticlesPXR: warning - no particles will be loaded. "
            print " You should either set lallindomain=true or set zmmin and"
            print " zmmax so they are different from each other.            "
            print "========================================================="


    # --- if lalldomain==False, removes particles outside boundaries
    if not lallindomain:
        if lrz:
            r = sqrt(x*x+y*y)
            cond = (r>=xmmin) & (r<xmmax) \
                 & (z>=zmmin) & (z<zmmax)
        else:
            if l4symtry:
                xc = abs(x)
                yc = abs(y)
            else:
                if l2symtry:
                    xc = abs(x)
                    yc = y
                else:
                    xc = x
                    yc = y
            cond = (xc>=xmmin) & (xc<xmmax) \
                 & (yc>=ymmin) & (yc<ymmax) \
                 & (z >=zmmin) & (z <zmmax)

        if cond.size==0:return

        x = x[cond]
        y = y[cond]
        z = z[cond]
        vx = vx[cond]
        vy = vy[cond]
        vz = vz[cond]
        gi = gi[cond]
        pids = compress(cond,pids,0)

        # Update also the number of particles
        nps0 = x.size

    if lnewparticles:
        # --- Set time of creation
        if top.tbirthpid > 0:
            pids[:,top.tbirthpid-1] = top.time

        # --- Set xyz old
        if top.xoldpid > 0: pids[:,top.xoldpid-1] = x
        if top.yoldpid > 0: pids[:,top.yoldpid-1] = y
        if top.zoldpid > 0: pids[:,top.zoldpid-1] = z

        if lmomentum:
            if top.uxoldpid > 0: pids[:,top.uxoldpid-1] = vx
            if top.uyoldpid > 0: pids[:,top.uyoldpid-1] = vy
            if top.uzoldpid > 0: pids[:,top.uzoldpid-1] = vz
        else:
            if top.uxoldpid > 0: pids[:,top.uxoldpid-1] = vx/gi
            if top.uyoldpid > 0: pids[:,top.uyoldpid-1] = vy/gi
            if top.uzoldpid > 0: pids[:,top.uzoldpid-1] = vz/gi

    # --- Call to PICSAR function for adding particles
    pxr.py_add_particles_to_species(self.pxr_species_array, nps0,pxr.npid,
                                    x,
                                    y,
                                    z,
                                    vx,
                                    vy,
                                    vz,
                                    gi,
                                    pids)

    aliasparticlearrays()

def aliasparticlearrays():
    global listofallspecies
    # --- Detect if tile arrays have been reallocated in PXR
    # --- and make proper aliasing in WARP

    isrealloc=zeros((pxr.ntilex,pxr.ntiley,pxr.ntilez),dtype=dtype('i8'))
    for i,s in enumerate(listofallspecies):
        pxr.get_are_tiles_reallocated(i+1, pxr.ntilex, pxr.ntiley, pxr.ntilez,isrealloc)
        ix,iy,iz=where(isrealloc==1)
        for il in range(0,len(ix)):
            pg = s.pgroups[iz[il]][iy[il]][ix[il]]
            pxr.point_to_tile(i+1, ix[il]+1, iy[il]+1, iz[il]+1)
            pg.npmax = 0
            pxr.partnmax
            pg.ns=1
            pg.npid=pxr.npid
            pg.gchange()
            pg.sq = s.charge
            pg.sm = s.mass
            pg.sw = s.sw
            pg.npmax = pxr.partnmax
            pg.nps = pxr.partn
            pg.ins[0] = 1
            pg.sid[0]=0
            pg.xp = pxr.partx
            pg.yp = pxr.party
            pg.zp = pxr.partz
            pg.uxp = pxr.partux
            pg.uyp = pxr.partuy
            pg.uzp = pxr.partuz
            #pg.pid = fzeros([pg.npmax,top.npid])
            pg.pid = pxr.pid
            pg.gaminv = pxr.partgaminv
            pg.ex = pxr.partex
            pg.ey = pxr.partey
            pg.ez = pxr.partez
            pg.bx = pxr.partbx
            pg.by = pxr.partby
            pg.bz = pxr.partbz
        pxr.set_are_tiles_reallocated(i+1, pxr.ntilex,pxr.ntiley,pxr.ntilez,zeros((pxr.ntilex,pxr.ntiley,pxr.ntilez),dtype=dtype('i8')))

def get_quantity_pxr( self, quantity, gather=True, bcast=False, **kw ):
    """
        Rewrite the method get_quantity of the class Species when pxr is loaded.
        Return the given pxr array for a given 'quantity'.

        Parameters:
        -----------
        self: Species
            Be careful, here 'self' is not relative to EM3DPXR but is used by
            the class Species.

        quantity: String
            must be choosen as like 'x', 'ux', 'xold', 'w', 'ex' etc...

        gather: bool
            If False: this function returns the particles from the local procs
            If True: this functions returns the gathered particles on all procs

        bcast: bool
            Only used when gather is True
            If bcast=False: only proc 0 gathers and returns the particles
            If bcast=False: all proc gather and return the particles
    """
    quantity_dict = dict(x=1, y=2, z=3, ux=4, uy=5, uz=6, ex=7, ey=8,
                         ez=9, bx=10, by=11, bz=12)

    quantity_pid_dict = dict()
    if top.xoldpid is not None:
        quantity_pid_dict['xold'] = top.xoldpid
    if top.yoldpid is not None:
        quantity_pid_dict['yold'] = top.yoldpid
    if top.zoldpid is not None:
        quantity_pid_dict['zold'] = top.zoldpid
    if top.uxoldpid is not None:
        quantity_pid_dict['uxold'] = top.uxoldpid
    if top.uyoldpid is not None:
        quantity_pid_dict['uyold'] = top.uyoldpid
    if top.uzoldpid is not None:
        quantity_pid_dict['uzold'] = top.uzoldpid
    if top.ssnpid is not None:
        quantity_pid_dict['id'] = top.ssnpid
    if top.wpid is not None:
        quantity_pid_dict['w'] = top.wpid
    if(top.exoldpid) is not None:
        quantity_pid_dict['exold'] = top.exoldpid
    if(top.eyoldpid) is not None:
        quantity_pid_dict['eyold'] = top.eyoldpid
    if(top.ezoldpid) is not None:
        quantity_pid_dict['ezold'] = top.ezoldpid
    if(top.bxoldpid) is not None:
        quantity_pid_dict['bxold'] = top.bxoldpid
    if(top.byoldpid) is not None:
        quantity_pid_dict['byold'] = top.byoldpid
    if(top.bzoldpid) is not None:
        quantity_pid_dict['bzold'] = top.bzoldpid


    js = self.pxr_species_array
    nb = numpy.empty(1,dtype=numpy.int64)
    pxr.get_local_number_of_particles_from_species(js, nb )

    quantity_array = numpy.empty(nb[0], dtype=numpy.float64, order='F')

    # Usual variables such as positionn, momentum or field
    if  quantity in quantity_dict:
        pxr.getquantity(js, quantity_dict[quantity], nb,
                        quantity_array)

    # Pid variables such as old variables or weight
    elif quantity in quantity_pid_dict:
        pxr.getquantity_pid(js, quantity_pid_dict[quantity], nb,
                            quantity_array)

    else:
        return "Error in get_quantity, key '%s' undefined or top.pid=None. \
           Please choose something among 'x', 'y', 'z', 'ux', 'uy', 'uz', \
           'ex', 'ey', 'ez', 'bx', 'by', 'bz' 'w', 'id', 'xold', 'yold', \
           'zold', 'uxold', 'uyold', 'uzold' or define top.pid."%quantity

    if lparallel and gather:
        return gatherarray(quantity_array,bcast=bcast)
    else:
        return quantity_array

def getx(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('x', gather=gather, bcast=bcast, **kw)

def gety(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('y', gather=gather, bcast=bcast, **kw)

def getz(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('z', gather=gather, bcast=bcast, **kw)

def getux(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('ux', gather=gather, bcast=bcast, **kw)

def getuy(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('uy', gather=gather, bcast=bcast, **kw)

def getuz(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('uz', gather=gather, bcast=bcast, **kw)

def getxold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('xold', gather=gather, bcast=bcast, **kw)

def getyold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('yold', gather=gather, bcast=bcast, **kw)

def getzold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('zold', gather=gather, bcast=bcast, **kw)

def getuxold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('uxold', gather=gather, bcast=bcast, **kw)

def getuyold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('uyold', gather=gather, bcast=bcast, **kw)

def getuzold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('uzold', gather=gather, bcast=bcast, **kw)

def getssn(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('id', gather=gather, bcast=bcast, **kw)

def getw(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('w', gather=gather, bcast=bcast, **kw)

def getex(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('ex', gather=gather, bcast=bcast, **kw)

def getey(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('ey', gather=gather, bcast=bcast, **kw)

def getez(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('ez', gather=gather, bcast=bcast, **kw)

def getbx(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('bx', gather=gather, bcast=bcast, **kw)

def getby(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('by', gather=gather, bcast=bcast, **kw)

def getbz(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('bz', gather=gather, bcast=bcast, **kw)


def getexold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('exold', gather=gather, bcast=bcast, **kw)

def geteyold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('eyold', gather=gather, bcast=bcast, **kw)

def getezold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('ezold', gather=gather, bcast=bcast, **kw)

def getbxold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('bxold', gather=gather, bcast=bcast, **kw)

def getbyold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('byold', gather=gather, bcast=bcast, **kw)

def getbzold(self, gather=1, bcast=None, **kw ):
    return self.get_quantity_pxr('bzold', gather=gather, bcast=bcast, **kw)



def getn(self, gather=1, bcast=None, **kw ):
    js = self.pxr_species_array
    nb = numpy.empty(1,dtype=numpy.int64)
    pxr.get_local_number_of_particles_from_species(js, nb )
    return nb[0]



def initialize_virtual_particles( self, w3d ):
    """
    This function overwrites the LaserAntenna class method
    initialize_virtual_particles, when using picsar.
    It creates new antenna macroparticles in picsar.

    Initialization of the antenna particles depending on the dimension and
    the laser propagation vector.
    """

    def switch_min_max( x1, x2, u ):
        """
        Return x1 or x2 depending on the sign of u
        """
        if u >= 0 :
            return x1
        else:
            return x2

    # Shortcut definition
    x0 = self.spot[0]
    y0 = self.spot[1]
    z0 = self.spot[2]
    xmin = w3d.xmminlocal
    xmax = w3d.xmmaxlocal
    ymin = w3d.ymminlocal
    ymax = w3d.ymmaxlocal
    zmin = w3d.zmminlocal
    zmax = w3d.zmmaxlocal

    if self.dim == "1d":
        # Ux is chosen orthogonal to self.vector in the plane (x,z)
        Uy = np.array([0.,1.,0.])
        Ux = np.cross(Uy,self.vector)
        self.Ux = Ux
        self.Uy = Uy

        # 1D injection along x
        self.nn_global = 1
        self.xx_global = x0 + np.zeros(self.nn_global)
        self.yy_global = y0 + np.zeros(self.nn_global)
        self.zz_global = z0 + np.zeros(self.nn_global)

    elif self.dim == "circ":
        # 2D circ

        # Check that the normal vector is along z and that the
        # position of the antenna is on the axis
        # (Otherwise the simulation cannot be performed in cylindrical
        # coordinate)
        assert np.allclose( self.vector, np.array([0,0,1]) )
        assert np.allclose( self.spot[:1], np.array([0,0]) )
        # Get the vectors that give the coordinate system of the antenna
        Ux = self.polvector
        Uy = self.polvector_2
        self.Ux = Ux
        self.Uy = Uy

        # The points of the antenna are along a star-pattern
        imin = np.floor( xmin/w3d.dx )
        imax = np.floor( xmax/w3d.dx )
        rr = w3d.dx * np.arange( imin, imax+1 )
        self.weights_circ = 2 * np.pi * rr / w3d.dx
        self.weights_circ /= 4 * self.circ_m
        w0 = self.weights_circ.copy()
        self.xx_global = rr.copy()
        self.yy_global = np.zeros_like( self.xx_global  )
        for i in range( 1, 4*self.circ_m ):
            phase = 0.5*np.pi*float(i)/self.circ_m
            self.xx_global = np.concatenate( (self.xx_global,
                                                rr*np.cos(phase)) )
            self.yy_global = np.concatenate( (self.yy_global,
                                                rr*np.sin(phase)) )
            self.weights_circ = np.concatenate((self.weights_circ,w0))
        self.nn_global = np.shape(self.xx_global)[0]
        self.zz_global = z0 + np.zeros(self.nn_global)

    elif self.dim == "2d":
        # 2D plane

        # Ux is chosen orthogonal to self.vector in the plane (x,z)
        Uy = np.array([0.,1.,0.])
        Ux = np.cross(Uy,self.vector)
        self.Ux = Ux
        self.Uy = Uy

        # Spacing between virtual particles to ensure at least one
        # particle per cell
        # select only the Ux components different from 0
        list_Ux = []
        if not Ux[0] == 0.: list_Ux.append(w3d.dx/np.abs(Ux[0]))
        if not Ux[1] == 0.: list_Ux.append(w3d.dy/np.abs(Ux[1]))
        if not Ux[2] == 0.: list_Ux.append(w3d.dz/np.abs(Ux[2]))
        self.Sx = min(list_Ux)

        # Boundaries of the box, depending on sign of the components of Ux
        xmin_i = switch_min_max(xmin, xmax, Ux[0])
        ymin_i = switch_min_max(ymin, ymax, Ux[1])
        zmin_i = switch_min_max(zmin, zmax, Ux[2])
        xmax_i = switch_min_max(xmax, xmin, Ux[0])
        ymax_i = switch_min_max(ymax, ymin, Ux[1])
        zmax_i = switch_min_max(zmax, zmin, Ux[2])

        # Find the range of integer with which the particles will be
        # initialized
        imin = Ux[0]*(xmin_i-x0) + Ux[1]*(ymin_i-y0) + Ux[2]*(zmin_i-z0)
        imax = Ux[0]*(xmax_i-x0) + Ux[1]*(ymax_i-y0) + Ux[2]*(zmax_i-z0)
        imin = np.floor(imin/self.Sx)
        imax = np.floor(imax/self.Sx)+1
        antenna_i = np.arange(imin, imax)

        # Initialize the particle positions
        self.xx_global = x0 + self.Sx*Ux[0]*antenna_i
        self.zz_global = z0 + self.Sx*Ux[2]*antenna_i

        # Keep only the particles that are inside the global box
        is_in_global_box = (self.xx_global >= xmin) \
                            & (self.xx_global < xmax) \
                            & (self.zz_global >= zmin) \
                            & (self.zz_global < zmax)
        self.zz_global = self.zz_global[is_in_global_box]
        self.xx_global = self.xx_global[is_in_global_box]
        self.yy_global = np.zeros(len(self.xx_global))
        # Number of virtual particles
        self.nn_global = np.shape(self.xx_global)[0]

    else:
        # 3D case, Ux = polvector and Uy = polvector_2
        Ux = self.polvector
        Uy = self.polvector_2
        self.Ux = Ux
        self.Uy = Uy

        # Spacing between virtual particles to ensure at least
        # one particle per cell
        # select only the components of Ux and Uy different from 0
        list_Ux = []; list_Uy = []
        if not Ux[0] == 0.: list_Ux.append( w3d.dx/np.abs(Ux[0]) )
        if not Ux[1] == 0.: list_Ux.append( w3d.dy/np.abs(Ux[1]) )
        if not Ux[2] == 0.: list_Ux.append( w3d.dz/np.abs(Ux[2]) )
        if not Uy[0] == 0.: list_Uy.append( w3d.dx/np.abs(Uy[0]) )
        if not Uy[1] == 0.: list_Uy.append( w3d.dy/np.abs(Uy[1]) )
        if not Uy[2] == 0.: list_Uy.append( w3d.dz/np.abs(Uy[2]) )
        self.Sx = min(list_Ux)
        self.Sy = min(list_Uy)

        # Boundaries of the box, depending on sign of the components of Ux
        xmin_i = switch_min_max(xmin, xmax, Ux[0])
        ymin_i = switch_min_max(ymin, ymax, Ux[1])
        zmin_i = switch_min_max(zmin, zmax, Ux[2])
        xmax_i = switch_min_max(xmax, xmin, Ux[0])
        ymax_i = switch_min_max(ymax, ymin, Ux[1])
        zmax_i = switch_min_max(zmax, zmin, Ux[2])
        # Boundaries of the box, depending on sign of the components of Uy
        xmin_j = switch_min_max(xmin, xmax, Uy[0])
        ymin_j = switch_min_max(ymin, ymax, Uy[1])
        zmin_j = switch_min_max(zmin, zmax, Uy[2])
        xmax_j = switch_min_max(xmax, xmin, Uy[0])
        ymax_j = switch_min_max(ymax, ymin, Uy[1])
        zmax_j = switch_min_max(zmax, zmin, Uy[2])

        # Find the range of integer with which the particles will be
        # initialized
        imin = Ux[0]*(xmin_i-x0) + Ux[1]*(ymin_i-y0) + Ux[2]*(zmin_i-z0)
        imax = Ux[0]*(xmax_i-x0) + Ux[1]*(ymax_i-y0) + Ux[2]*(zmax_i-z0)
        jmin = Uy[0]*(xmin_j-x0) + Uy[1]*(ymin_j-y0) + Uy[2]*(zmin_j-z0)
        jmax = Uy[0]*(xmax_j-x0) + Uy[1]*(ymax_j-y0) + Uy[2]*(zmax_j-z0)
        imin = np.floor(imin/self.Sx)
        imax = np.floor(imax/self.Sx)+1
        jmin = np.floor(jmin/self.Sy)
        jmax = np.floor(jmax/self.Sy)+1
        array_i = np.arange(imin, imax)
        array_j = np.arange(jmin, jmax)
        antenna_i, antenna_j = np.meshgrid(array_i,array_j)

        # Initialize the particle positions
        self.xx_global = x0 + self.Sx*Ux[0]*antenna_i + self.Sy*Uy[0]*antenna_j
        self.yy_global = y0 + self.Sx*Ux[1]*antenna_i + self.Sy*Uy[1]*antenna_j
        self.zz_global = z0 + self.Sx*Ux[2]*antenna_i + self.Sy*Uy[2]*antenna_j
        self.xx_global = self.xx_global.flatten()
        self.yy_global = self.yy_global.flatten()
        self.zz_global = self.zz_global.flatten()

        # Keep only the particles that are inside the global box
        is_in_global_box = (self.xx_global >= xmin) \
                            & (self.xx_global < xmax) \
                            & (self.yy_global >= ymin) \
                            & (self.yy_global < ymax) \
                            & (self.zz_global >= zmin) \
                            & (self.zz_global < zmax)
        self.zz_global = self.zz_global[is_in_global_box]
        self.yy_global = self.yy_global[is_in_global_box]
        self.xx_global = self.xx_global[is_in_global_box]
        # Number of virtual particles
        self.nn_global = np.shape(self.xx_global)[0]

    # Set the deplacement around the initial position and normalized momenta
    # variation of each macroparticles to 0
    self.xdx_global = np.zeros(self.nn_global)
    self.ydy_global = np.zeros(self.nn_global)
    self.zdz_global = np.zeros(self.nn_global)
    self.ux_global = np.zeros(self.nn_global)
    self.uy_global = np.zeros(self.nn_global)
    self.uz_global = np.zeros(self.nn_global)
    self.gi_global = np.ones(self.nn_global)

    # Calculate the weights
    self.weights_global = np.ones(self.nn_global) * eps0*self.emax/0.01

    if self.dim == "2d":
        self.weights_global *= self.Sx
    elif self.dim == "3d" :
        self.weights_global *= self.Sy*self.Sx
    elif self.circ_m > 0 : # Circ
        # Laser initialized with particles in a star-pattern
        self.weights_global *= w3d.dx**2 * self.weights_circ


    # Create two new antenna species in PICSAR and assign their caracteristics
    js_laser_pos = numpy.empty(1,dtype=numpy.int64)
    js_laser_neg = numpy.empty(1,dtype=numpy.int64)

    pxr.init_laser_species_python(self.emax, self.spot, self.vector, Ux, Uy,
                                    1., self.weights_global, self.xx_global,
                                self.yy_global, self.zz_global, self.nn_global,
                                js_laser_pos)

    pxr.init_laser_species_python(self.emax, self.spot, self.vector, Ux, Uy,
                                -1., self.weights_global, self.xx_global,
                                self.yy_global, self.zz_global, self.nn_global,
                                js_laser_neg)
    self.js_pos = js_laser_pos
    self.js_neg = js_laser_neg

def push_virtual_particles(self, top, f, clight ):
    """
    This function overwrites the LaserAntenna class method
    push_virtual_particles, when using picsar.

    Calculate the motion parameters of the laser antenna at a given
    timestep
    """

    dt = top.dt

    # Coordinate of the antenna in the plane (Ux,Uy)
    wpid = pxr.wpid

    for js in [self.js_neg, self.js_pos]:
        nb = numpy.empty(1,dtype=numpy.int64)
        pxr.get_local_number_of_particles_from_species(js, nb )

        quantity_array = numpy.empty(nb[0], dtype=numpy.float64, order='F')
        pxr.getquantity_pid(js, wpid+1, nb, quantity_array)
        x = quantity_array

        quantity_array = numpy.empty(nb[0], dtype=numpy.float64, order='F')
        pxr.getquantity_pid(js, wpid+2, nb, quantity_array)
        y = quantity_array

        t = top.time*(1.-self.v/clight)
        amp = self.laser_func(x,y,t)

        # --- displaces fixed weight particles on "continuous" trajectories
        dispmax = 0.01*clight
        coef_ampli = dispmax * (1.-self.v/clight) / self.emax

        if isinstance(amp,list): #elliptic polarization
            amp_x = amp[0]*self.polvector[0] + amp[1]*self.polvector_2[0]
            amp_y = amp[0]*self.polvector[1] + amp[1]*self.polvector_2[1]
            amp_z = amp[0]*self.polvector[2] + amp[1]*self.polvector_2[2]

            amplitude_x = coef_ampli * amp_x
            amplitude_y = coef_ampli * amp_y
            amplitude_z = coef_ampli * amp_z

        else: #linear polarization
            amplitude_x = coef_ampli * amp * self.polvector[0]
            amplitude_y = coef_ampli * amp * self.polvector[1]
            amplitude_z = coef_ampli * amp * self.polvector[2]

        pxr.laser_pusher_profile(js, amplitude_x,amplitude_y,amplitude_z, nb,self.vx,self.vy,self.vz)

def select_particles_in_local_box(self, w3d, zgrid):
    """
    This function overwrites the LaserAntenna class method
    select_particles_in_local_box, when using picsar.

    Since the particles are exchanged via basic exchange routines,
    this function is not used anymore.

    The "self.nn = 0" command assures that the current and charge
    deposition are not done through warp.
    """
    self.nn = 0
    return

class EM3DPXR(EM3DFFT):

    __em3dpxrinputs__ = []
    __flaginputs__ = {'ntilex':1,
                      'ntiley':1,
                      'ntilez':1,
                      'listofallspecies':[],
                      'dload_balancing':0,
                      'dlb_freq':1,
                      'dlb_threshold':20,
                      'dlb_at_init':1,
                      'it_dlb_init':11,
                      'l_output_grid':0,
                      'l_output_freq':1,
                      'rhodepo':0,      # Charge deposition method
                      'currdepo':0,     # Current deposition method
                      'mpicom_curr':1,  # Com type Current deposition
                      'fieldgathe':0,   # Field gathering method
                      'partcom':0,      # Particle communication
                      'fg_p_pp_separated':0,
                      'lvec_curr_depo':8,
                      'lvec_charge_depo':64,
                      'lvec_fieldgathe':0,
                      'mpi_buf_size':2000,
                      'sorting':None,
                      'l_debug':0,
                      'l_reinject':[0, 0, 0, 0, 0, 0],
                      'offset_x_part_grid':[0.,0.],
                      'offset_y_part_grid':[0.,0.],
                      'offset_z_part_grid':[0.,0.],
                      'full_pxr': False,
                      'fftw_hybrid':False,
                      'fftw_with_mpi':False,
                      'fftw_mpi_transpose':False,
                      'p3dfft_flag':False,
                      'p3dfft_stride':False,
                      'nb_group_x':0,
                      'nb_group_y':0,
                      'nb_group_z':0,
                      'nyg_group':0,
                      'nzg_group':0,
                      'nx_pml':8,
                      'ny_pml':8,
                      'nz_pml':8,
                      'shift_x_pml_pxr':4, # number of guardcells whre fields are forced to 0
                      'shift_y_pml_pxr':4, # when using pml with full pxr mode. This parameter
                      'shift_z_pml_pxr':4, # is only
                      'absorbing_bcs_x':0,
                      'absorbing_bcs_y':0,
                      'absorbing_bcs_z':0,
                      'g_spectral':False,
                      'pxr_antenna':True,
                      }

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(EM3DPXR.__flaginputs__,kw)

        if (self.l_debug):
          print("Call __init__")
          print(' Debug prints activated')

        if(self.pxr_antenna):
          top.nextpid()
          top.nextpid()

        EM3DFFT.__init__(self,kwdict=kw)

        self.l_pxr = l_pxr
        self.l_fftw = l_fftw

        # If sorting undefined
        if self.sorting==None:
          self.sorting = Sorting([],[],activated=0,dx=1.,dy=1.,dz=1.,xshift=-0.5,yshift=-0.5,zshift=-0.5)

        if (self.l_debug): print("End __init__")

    def finalize(self,lforce=False):

        if self.finalized and not lforce: return
        if(self.l_debug): print("begin finalize")
        if self.l_pxr:
          EM3D.finalize(self)

          # Rewrite the LaserAntenna class methods for pxr
          if(self.pxr_antenna):
            LaserAntenna.initialize_virtual_particles  = \
                                                    initialize_virtual_particles
            LaserAntenna.push_virtual_particles        = push_virtual_particles
            LaserAntenna.select_particles_in_local_box = \
                                                 select_particles_in_local_box


          self.allocatefieldarraysFFT()
          self.allocatefieldarraysPXR()
          #if full_pxr == True additional computations are done in picsar.
          # This includes PSATD block initialization and fields boundaries through PMLS.
          #This mode also allows to use FFTW_MPI or P3DFFT in order to perform FFT computations.

          if(self.full_pxr):
            pxr.init_plans_blocks()
            if(pxr.absorbing_bcs):
              pxr.init_pml_arrays()
              # Creates a pointer to each pxr pml sub field
              self.exy_pxr = pxr.exy
              self.exz_pxr = pxr.exz
              self.eyx_pxr = pxr.eyx
              self.eyz_pxr = pxr.eyz
              self.ezx_pxr = pxr.ezx
              self.ezy_pxr = pxr.ezy
              self.bxy_pxr = pxr.bxy
              self.bxz_pxr = pxr.bxz
              self.byx_pxr = pxr.byx
              self.byz_pxr = pxr.byz
              self.bzx_pxr = pxr.bzx
              self.bzy_pxr = pxr.bzy


          self.lorentz_transform2d = pxr.transform_lorentz2d
          self.lorentz_transform3d = pxr.transform_lorentz3d
          self.lorentz_transform_parts_with_fields = pxr.lorentz_transform_parts_with_fields
          self.lorentz_transform_parts_without_fields = pxr.lorentz_transform_parts_without_fields


          # Rewrite the get_quantity methods to the class species for pxr
          Species.get_quantity_pxr = get_quantity_pxr
          Species.getx             = getx
          Species.gety             = gety
          Species.getz             = getz
          Species.getux            = getux
          Species.getuy            = getuy
          Species.getuz            = getuz
          Species.getxold          = getxold
          Species.getyold          = getyold
          Species.getzold          = getzold
          Species.getuxold         = getuxold
          Species.getuyold         = getuyold
          Species.getuzold         = getuzold
          Species.getssn           = getssn
          Species.getweights       = getw
          Species.getw             = getw
          Species.getex            = getex
          Species.getey            = getey
          Species.getez            = getez
          Species.getbx            = getbx
          Species.getby            = getby
          Species.getbz            = getbz
          Species.getn             = getn
          Species.getexold         = getexold
          Species.geteyold         = geteyold
          Species.getezold         = getezold
          Species.getbxold         = getbxold
          Species.getbyold         = getbyold
          Species.getbzold         = getbzold

        else:
          EM3DFFT.finalize(self)

        Species.addparticles = Species.addpart = addparticlesPXR

    def convertindtoproc(self,ix,iy,iz,nx,ny,nz):
      ixt = ix
      iyt = iy
      izt = iz

      if (ixt < 0   ): ixt = nx - 1
      if (ixt > nx-1): ixt = 0
      if (iyt < 0   ): iyt = ny - 1
      if (iyt > ny-1): iyt = 0
      if (izt < 0   ): izt = nz - 1
      if (izt > nz-1): izt = 0

      convertindextoproc = ixt + iyt*nx + izt*nx*ny

      return convertindextoproc

    def allocatefieldarraysPXR(self):

        if (self.l_debug): print("allocatefieldarraysPXR")

        # Set up case dimensionality
        if (self.l_2dxz):
          pxr.c_dim=2
        else:
          pxr.c_dim=3

        # Set up PXR MPI Data
        if (self.l_debug): print(" Setup PXR MPI Data")
        pxr.nprocx=top.fsdecomp.nxprocs
        pxr.nprocy=top.fsdecomp.nyprocs
        pxr.nprocz=top.fsdecomp.nzprocs
        ixcpu=top.fsdecomp.ixproc
        iycpu=top.fsdecomp.iyproc
        izcpu=top.fsdecomp.izproc
        for iz in range(-1,2):
            for iy in range(-1,2):
                for ix in range(-1,2):
                    indtoproc=self.convertindtoproc(ixcpu+ix,iycpu+iy,izcpu+iz,pxr.nprocx,pxr.nprocy,pxr.nprocz)
                    pxr.neighbour[ix+1,iy+1,iz+1]=indtoproc
      #  pxr.proc_x_max =  pxr.neighbour[1,0,0 ]
      #  pxr.proc_x_min =  pxr.neighbour[-1,0,0]
      #  pxr.proc_y_max =  pxr.neighbour[0,1,0 ]
      #  pxr.proc_y_min =  pxr.neighbour[0,-1,0]
      #  pxr.proc_z_max =  pxr.neighbour[0,0,1 ]
      #  pxr.proc_z_min =  pxr.neighbour[0,0,-1]

        if (ixcpu==0):
            pxr.x_min_boundary=1
        if (ixcpu==pxr.nprocx-1):
            pxr.x_max_boundary=1
        if (iycpu==0):
            pxr.y_min_boundary=1
        if (iycpu==pxr.nprocy-1):
            pxr.y_max_boundary=1
        if (izcpu==0):
            pxr.z_min_boundary=1
        if (izcpu==pxr.nprocz-1):
            pxr.z_max_boundary=1

        pxr.x_coords=ixcpu
        pxr.y_coords=iycpu
        pxr.z_coords=izcpu

        # MPI boundaries index in global array
        if (self.l_debug): print(" MPI boundaries index in global array")
        pxr.cell_x_min=top.fsdecomp.ix
        pxr.cell_y_min=top.fsdecomp.iy
        pxr.cell_z_min=top.fsdecomp.iz
        pxr.cell_x_max=pxr.cell_x_min+(top.fsdecomp.nx-1)
        pxr.cell_y_max=pxr.cell_y_min+(top.fsdecomp.ny-1)
        pxr.cell_z_max=pxr.cell_z_min+(top.fsdecomp.nz-1)
        pxr.x_grid_mins=top.fsdecomp.xmin
        pxr.x_grid_maxs=top.fsdecomp.xmax
        pxr.y_grid_mins=top.fsdecomp.ymin
        pxr.y_grid_maxs=top.fsdecomp.ymax
        pxr.z_grid_mins=top.fsdecomp.zmin
        pxr.z_grid_maxs=top.fsdecomp.zmax

        # Particle boundaries for PXR
        if (self.l_debug): print(" Setup particle boundaries for PXR")
        if (top.pbound0 == absorb):
          if (self.l_reinject[4]):
              pxr.pbound_z_min=3
          else:
              pxr.pbound_z_min=1
        elif(top.pbound0 == reflect):
          pxr.pbound_z_min=2
        else: # Default is periodic
          pxr.pbound_z_min=0

        if (top.pboundnz == absorb):
            if (self.l_reinject[5]):
                pxr.pbound_z_max=3
            else:
                pxr.pbound_z_max=1
        elif(top.pboundnz == reflect):
            pxr.pbound_z_max=2
        else: # Default is periodic
            pxr.pbound_z_max=0

        if (top.pboundxy == absorb):
            if (self.l_reinject[0]):
                pxr.pbound_x_min=3
            else:
                pxr.pbound_x_min=1
            if (self.l_reinject[1]):
                pxr.pbound_x_max=3
            else:
                pxr.pbound_x_max=1
            if (self.l_reinject[2]):
                pxr.pbound_y_min=3
            else:
                pxr.pbound_y_min=1
            if (self.l_reinject[3]):
                pxr.pbound_y_max=3
            else:
                pxr.pbound_y_max=1

        elif(top.pboundxy == reflect):
            pxr.pbound_x_min=2
            pxr.pbound_x_max=2
            pxr.pbound_y_min=2
            pxr.pbound_y_max=2
        else: # Default is periodic
            pxr.pbound_x_min=0
            pxr.pbound_x_max=0
            pxr.pbound_y_min=0
            pxr.pbound_y_max=0

        # --- number of grid cells
        if (self.l_debug): (" Setup number of grid cells for PXR")
        pxr.nx_global = w3d.nx
        pxr.ny_global = w3d.ny
        pxr.nz_global = w3d.nz
        pxr.nx_global_grid = pxr.nx_global+1
        pxr.ny_global_grid = pxr.ny_global+1
        pxr.nz_global_grid = pxr.nz_global+1

        pxr.nx = self.nxlocal
        pxr.ny = self.nylocal
        pxr.nz = self.nzlocal
        pxr.nx_grid=pxr.nx+1
        pxr.ny_grid=pxr.ny+1
        pxr.nz_grid=pxr.nz+1


        # --- number of guard cells
        if (self.l_debug): print(" Setup number of guard cells for PXR")
        pxr.nxguards = self.nxguard
        pxr.nyguards = self.nyguard
        pxr.nzguards = self.nzguard
        pxr.nxjguards = self.nxguard
        pxr.nyjguards = self.nyguard
        pxr.nzjguards = self.nzguard

        # --- Grid domain extents and dimensions
        pxr.xmin = w3d.xmmin
        pxr.ymin = w3d.ymmin
        pxr.zmin = w3d.zmmin
        pxr.xmax = w3d.xmmax
        pxr.ymax = w3d.ymmax
        pxr.zmax = w3d.zmmax
        pxr.x_grid_min=pxr.xmin
        pxr.x_grid_max=pxr.xmax
        pxr.y_grid_min=pxr.ymin
        pxr.y_grid_max=pxr.ymax
        pxr.z_grid_min=pxr.zmin
        pxr.z_grid_max=pxr.zmax

        pxr.x_min_local = self.fields.xmin
        pxr.x_max_local = self.fields.xmax
        pxr.y_min_local = self.fields.ymin
        pxr.y_max_local = self.fields.ymax
        pxr.z_min_local = self.fields.zmin
        pxr.z_max_local = self.fields.zmax
        pxr.x_grid_min_local=pxr.x_min_local
        pxr.x_grid_max_local=pxr.x_max_local
        pxr.y_grid_min_local=pxr.y_min_local
        pxr.y_grid_max_local=pxr.y_max_local
        pxr.z_grid_min_local=pxr.z_min_local
        pxr.z_grid_max_local=pxr.z_max_local
        #pxr.zgrid=top.zgrid

        pxr.length_x = pxr.xmax-pxr.xmin
        pxr.length_y = pxr.ymax-pxr.ymin
        pxr.length_z = pxr.zmax-pxr.zmin

        # --- Particle domain extents and dimensions
        # Set offset Grid/part in PXR:
        pxr.offset_grid_part_x_min = self.offset_x_part_grid[0]
        pxr.offset_grid_part_x_max = self.offset_x_part_grid[1]
        pxr.offset_grid_part_y_min = self.offset_y_part_grid[0]
        pxr.offset_grid_part_y_max = self.offset_y_part_grid[1]
        pxr.offset_grid_part_z_min = self.offset_z_part_grid[0]
        pxr.offset_grid_part_z_max = self.offset_z_part_grid[1]

        # Global part boundaries
        # Local part boundaries
        # - Xmin
        if (pxr.pbound_x_min == 3) or (pxr.pbound_x_min==1):
            pxr.xmin_part=pxr.xmin+pxr.offset_grid_part_x_min
            if (pxr.xmin_part >= pxr.x_min_local) and (pxr.xmin_part < pxr.x_max_local):
                pxr.x_min_boundary_part= 1
                pxr.x_min_local_part=pxr.xmin_part
            else:
                pxr.x_min_boundary_part= 0
                pxr.x_min_local_part = pxr.x_min_local
        else:
            pxr.xmin_part=pxr.xmin
            pxr.x_min_boundary_part=pxr.x_min_boundary
            pxr.x_min_local_part=pxr.x_min_local
        # - Xmax
        if (pxr.pbound_x_max == 3) or (pxr.pbound_x_max==1):
            pxr.xmax_part=pxr.xmax+pxr.offset_grid_part_x_max
            if (pxr.xmax_part >=  pxr.x_min_local) and (pxr.xmax_part < pxr.x_max_local):
                pxr.x_max_boundary_part= 1
                pxr.x_max_local_part=pxr.xmax_part
            else:
                pxr.x_max_boundary_part= 0
                pxr.x_max_local_part = pxr.x_max_local
        else:
            pxr.xmax_part=pxr.xmax
            pxr.x_max_boundary_part=pxr.x_max_boundary
            pxr.x_max_local_part=pxr.x_max_local
        # - Ymin
        if (pxr.pbound_y_min == 3) or (pxr.pbound_y_min==1):
            pxr.ymin_part=pxr.ymin+pxr.offset_grid_part_y_min
            if (pxr.ymin_part >=  pxr.y_min_local) and (pxr.ymin_part < pxr.y_max_local):
                pxr.y_min_boundary_part= 1
                pxr.y_min_local_part=pxr.ymin_part
            else:
                pxr.y_min_boundary_part= 0
                pxr.y_min_local_part = pxr.y_min_local
        else:
            pxr.ymin_part=pxr.ymin
            pxr.y_min_boundary_part=pxr.y_min_boundary
            pxr.y_min_local_part=pxr.y_min_local
        # - Ymax
        if (pxr.pbound_y_max == 3) or (pxr.pbound_y_max==1):
            pxr.ymax_part=pxr.ymax+pxr.offset_grid_part_y_max
            if (pxr.ymax_part >=  pxr.y_min_local) and (pxr.ymax_part < pxr.y_max_local):
                pxr.y_max_boundary_part= 1
                pxr.y_max_local_part=pxr.ymax_part
            else:
                pxr.y_max_boundary_part= 0
                pxr.y_max_local_part = pxr.y_max_local
        else:
            pxr.ymax_part=pxr.ymax
            pxr.y_max_boundary_part=pxr.y_max_boundary
            pxr.y_max_local_part=pxr.y_max_local
        # - Zmin
        if (pxr.pbound_z_min == 3) or (pxr.pbound_z_min==1):
            pxr.zmin_part=pxr.zmin+pxr.offset_grid_part_z_min
            if (pxr.zmin_part >=  pxr.z_min_local) and (pxr.zmin_part < pxr.z_max_local):
                pxr.z_min_boundary_part= 1
                pxr.z_min_local_part=pxr.zmin_part
            else:
                pxr.z_min_boundary_part= 0
                pxr.z_min_local_part = pxr.z_min_local
        else:
            pxr.zmin_part=pxr.zmin
            pxr.z_min_boundary_part=pxr.z_min_boundary
            pxr.z_min_local_part=pxr.z_min_local
        # - Zmax
        if (pxr.pbound_z_max == 3) or (pxr.pbound_z_max==1):
            pxr.zmax_part=pxr.zmax+pxr.offset_grid_part_z_max
            if (pxr.zmax_part >=  pxr.z_min_local) and (pxr.zmax_part < pxr.z_max_local):
                pxr.z_max_boundary_part= 1
                pxr.z_max_local_part=pxr.zmax_part
            else:
                pxr.z_max_boundary_part= 0
                pxr.z_max_local_part = pxr.z_max_local
        else:
            pxr.zmax_part=pxr.zmax
            pxr.z_max_boundary_part=pxr.z_max_boundary
            pxr.z_max_local_part=pxr.z_max_local

        # Particle domain extents
        pxr.length_x_part = pxr.xmax_part - pxr.xmin_part
        pxr.length_y_part = pxr.ymax_part - pxr.ymin_part
        pxr.length_z_part = pxr.zmax_part - pxr.zmin_part

        # INIT MPI_DATA FOR PICSAR
        # Init communicator variable in picsar
        if (self.l_debug): print(" Init communicator variable in PXR")
        pxr.mpi_minimal_init_python(top.fsdecomp.mpi_comm)

        # allocate grid quantities
        if(self.full_pxr):
          self.l_pxr = True
          pxr.l_spectral = True
          pxr.g_spectral = self.g_spectral
          pxr.fftw_with_mpi = self.fftw_with_mpi
          pxr.fftw_mpi_transpose = self.fftw_mpi_transpose
          pxr.fftw_hybrid = self.fftw_hybrid
          pxr.p3dfft_flag = self.p3dfft_flag
          pxr.p3dfft_stride = self.p3dfft_stride
          pxr.nxg_group = self.nxguard
          pxr.nyg_group = self.nyg_group
          pxr.nzg_group = self.nzg_group
          pxr.nb_group_z = self.nb_group_z
          pxr.nb_group_y = self.nb_group_y
          pxr.absorbing_bcs_x = self.absorbing_bcs_x
          pxr.absorbing_bcs_y = self.absorbing_bcs_y
          pxr.absorbing_bcs_z = self.absorbing_bcs_z


         # if(w3d.bound0 == openbc or w3d.boundnz==openbc):
         #   pxr.absorbing_bcs_z = True
         # if(w3d.boundxy == openbc):
         #   pxr.absorbing_bcs_x = True
         #   pxr.absorbing_bcs_y = True

          #Set aborbing_bcs flag to true if there is an absorbing bc in any direction
          #If absorbing bcs in one direction then increase the grid offset for particles
          #in order to avoid PMLS instabilities.
          if(pxr.absorbing_bcs_x or pxr.absorbing_bcs_y or pxr.absorbing_bcs_z):
            pxr.absorbing_bcs=True
            pxr.nx_pml = self.nx_pml
            pxr.ny_pml = self.ny_pml
            pxr.nz_pml = self.nz_pml

            pxr.shift_x_pml = self.shift_x_pml_pxr
            pxr.shift_y_pml = self.shift_y_pml_pxr
            pxr.shift_z_pml = self.shift_z_pml_pxr
            if(self.fftw_hybrid):
              pxr.shift_x_pml = pxr.nxguards
              pxr.shift_y_pml = pxr.nyguards
              pxr.shift_z_pml = pxr.nzguards

          self.absorbing_bcs_pxr = pxr.absorbing_bcs
          pxr.get_neighbours_python()
          if(pxr.absorbing_bcs==True):
            #if absorbing_bcs in warp then set warp bcs to periodic to avoid bugs (and useless block inits)
            pxr.g_spectral = True
            pxr.get_non_periodic_mpi_bcs()
          if(pxr.fftw_with_mpi):
            if(pxr.fftw_hybrid):
              pxr.setup_groups()
              pxr.get2d_intersection_group_mpi()
            else:
              pxr.adjust_grid_mpi_global()

        if (self.l_debug): print(" Allocate grid quantities in PXR")
        pxr.allocate_grid_quantities()
        if(self.l_debug): print("Compute simulation axis in PXR")
        pxr.compute_simulation_axis()

        # set time step
        pxr.dt = top.dt

        # --- Resolution
        if (self.l_debug): print(" Setup resolution and related variables in PXR")
        pxr.dx = self.dx
        pxr.dy = self.dy
        pxr.dz = self.dz
        pxr.dxi = 1./self.dx
        pxr.dyi = 1./self.dy
        pxr.dzi = 1./self.dz
        pxr.invvol = pxr.dxi*pxr.dyi*pxr.dzi
        pxr.dts2dx = 0.5*pxr.dt*pxr.dxi
        pxr.dts2dy = 0.5*pxr.dt*pxr.dyi
        pxr.dts2dz = 0.5*pxr.dt*pxr.dzi
        pxr.clightsq = 1.0/pxr.clight**2

        # --- Maxwell solver
        pxr.norderx = self.norderx
        pxr.nordery = self.nordery
        pxr.norderz = self.norderz

        pxr.xcoeffs = self.fields.xcoefs
        pxr.ycoeffs = self.fields.ycoefs
        pxr.zcoeffs = self.fields.zcoefs

        # Set coefficient for Maxwell solver
        if (self.l_debug): print(" Set coefficient for Maxwell solver")
        pxr.alphax = em3d.alphax
        pxr.alphay = em3d.alphay
        pxr.alphaz = em3d.alphaz
        pxr.betaxy = em3d.betaxy
        pxr.betayx = em3d.betayx
        pxr.betaxz = em3d.betaxz
        pxr.betazx = em3d.betazx
        pxr.betayz = em3d.betayz
        pxr.betazy = em3d.betazy
        pxr.gammax = em3d.gammax
        pxr.gammay = em3d.gammay
        pxr.gammaz = em3d.gammaz
        pxr.deltaz = em3d.deltaz

        pxr.ex = self.fields.Ex
        pxr.ey = self.fields.Ey
        pxr.ez = self.fields.Ez
        pxr.bx = self.fields.Bx
        pxr.by = self.fields.By
        pxr.bz = self.fields.Bz
        pxr.jx = self.fields.Jx
        pxr.jy = self.fields.Jy
        pxr.jz = self.fields.Jz
        if(self.full_pxr):
          pxr.rho = self.fields.Rho
          pxr.rhoold = self.fields.Rhoold
          if(pxr.absorbing_bcs):
            pxr.init_splitted_fields_random()


        pxr.ex_p = self.fields.Exp
        pxr.ey_p = self.fields.Eyp
        pxr.ez_p = self.fields.Ezp
        pxr.bx_p = self.fields.Bxp
        pxr.by_p = self.fields.Byp
        pxr.bz_p = self.fields.Bzp

        pxr.l_nodalgrid = self.l_nodalgrid

        pxr.nxs = 0
        pxr.nys = 0
        pxr.nzs = 0

        # Current deposition
        pxr.nox = top.depos_order[0][0]
        pxr.noy = top.depos_order[1][0]
        pxr.noz = top.depos_order[2][0]

        if (self.l_debug): print(" Set up algorithms in PXR")

        # Charge deposition algorithm
        pxr.rhodepo=self.rhodepo
        # Current deposition algorithm
        pxr.currdepo=self.currdepo
        # Tye of MPI communication for the current
        pxr.mpicom_curr=self.mpicom_curr
        # Field gathering method
        pxr.fieldgathe=self.fieldgathe
        # Particle communication
        pxr.partcom=self.partcom
        # Field gathering and PArticle pusher separated
        pxr.fg_p_pp_separated=self.fg_p_pp_separated
        # Particle pusher type
        pxr.particle_pusher = top.pgroup.lebcancel_pusher
        # lvec size for the current deposition
        pxr.lvec_curr_depo = self.lvec_curr_depo
        # lvec size for the charge deposition
        pxr.lvec_charge_depo = self.lvec_charge_depo
        # lvec size for the field gathering
        if (self.lvec_fieldgathe==0):
          if ((pxr.nox==3)and(pxr.noy==3)and(pxr.noz==3)):
            pxr.lvec_fieldgathe = 64
          else:
            pxr.lvec_fieldgathe = 512
        else:
          pxr.lvec_fieldgathe = self.lvec_fieldgathe
        # MPI buffer size for particle exchange
        pxr.mpi_buf_size = self.mpi_buf_size

        # Type of field gathering
        pxr.l4symtry=w3d.l4symtry
        pxr.l_lower_order_in_v = self.l_lower_order_in_v

        # --- Tiling parameters
        pxr.ntilex = self.ntilex
        pxr.ntiley = self.ntiley
        pxr.ntilez = self.ntilez

        # --- Sorting parameters
        if (self.l_debug): print(" Setup sorting parameters in PXR")
        pxr.sorting_activated = self.sorting.activated
        pxr.sorting_dx = self.sorting.dx*pxr.dx
        pxr.sorting_dy = self.sorting.dy*pxr.dy
        pxr.sorting_dz = self.sorting.dz*pxr.dz
        pxr.sorting_shiftx = self.sorting.xshift*pxr.dx
        pxr.sorting_shifty = self.sorting.yshift*pxr.dy
        pxr.sorting_shiftz = self.sorting.zshift*pxr.dz
        pxr.sorting_verbose = self.sorting.verbose

        # --- time statistics
        self.time_stat_loc_array = zeros([20])

        # --- allocates array of species
        if (self.l_debug): print(" Allocates array of species")
        pxr.init_species_section()

        for i,s in enumerate(self.listofallspecies):
            # Check for sorting
            if (i >= len(self.sorting.periods)):
              self.sorting.periods.append(0)
              self.sorting.starts.append(0)
            # initialize species in pxr
            pxr.set_particle_species_properties(i+1,s.name,s.mass,s.charge,0, \
                                                0.,0.,0.,0.,0.,0., \
                                                0.,0.,0.,0.,0.,0., \
                                                self.sorting.periods[i], \
                                                self.sorting.starts[i],  \
                                                s.pgroup.ldodepos[i])
            pxr.nspecies+=1
        pxr.npid=top.npid



        pxr.ssnpid=top.ssnpid
        pxr.set_tile_split()
        pxr.init_tile_arrays()

        # Add all particles of all species to PXR
        for i,s in enumerate(self.listofallspecies):
            pids=s.getpid(id=-1,bcast=0,gather=0)
            # In PXR, pid[:,wpid] is the weight of the particle
            # (but not in WARP so correct it to get good normalization)
            pids[:,top.wpid-1]*=s.sw
            s.sw0=s.sw*1.
            # Add particles of species s to PXR
            pxr.py_add_particles_to_species(i+1, s.nps,pxr.npid,
                                            s.getx(bcast=0,gather=0),
                                            s.gety(bcast=0,gather=0),
                                            s.getz(bcast=0,gather=0),
                                            s.getux(bcast=0,gather=0),
                                            s.getuy(bcast=0,gather=0),
                                            s.getuz(bcast=0,gather=0),
                                            s.getgaminv(bcast=0,gather=0),
                                            pids)
        # Removed duplicate species in WARP
        top.pgroup.npmax=0
        top.pgroup.ns=1
        top.pgroup.nps=0
        top.pgroup.gchange()

        # --- mirror PXR tile structure in Warp with list of pgroups
        if (self.l_debug): print(" Mirror PXR tile structure in Warp with list of pgroups")
        for i,s in enumerate(self.listofallspecies):
            s.pgroups = []
            s.jslist = [0]
            s.sw=1.
            s.pxr_species_array=i+1
            for iz in range(1,self.ntilez+1):
                xygroup=[]
                for iy in range(1,self.ntiley+1):
                    xgroup=[]
                    for ix in range(1,self.ntilex+1):
                        pg = ParticleGroup()
                        xgroup.append(pg)
                        pxr.point_to_tile(i+1, ix, iy, iz)
                        pg.npmax = 0
                        pxr.partnmax
                        pg.ns=1
                        pg.npid=top.npid
                        pg.gchange()
                        pg.sq = s.charge
                        pg.sm = s.mass
                        pg.sw = s.sw
                        pg.npmax = pxr.partnmax
                        pg.nps = pxr.partn
                        pg.ins[0] = 1
                        pg.sid[0]=0
                        pg.xp = pxr.partx
                        pg.yp = pxr.party
                        pg.zp = pxr.partz
                        pg.uxp = pxr.partux
                        pg.uyp = pxr.partuy
                        pg.uzp = pxr.partuz
                        #pg.pid = fzeros([pg.npmax,top.npid])
                        pg.pid = pxr.pid
                        pg.gaminv = pxr.partgaminv
                        pg.ex = pxr.partex
                        pg.ey = pxr.partey
                        pg.ez = pxr.partez
                        pg.bx = pxr.partbx
                        pg.by = pxr.partby
                        pg.bz = pxr.partbz
                        pg.lebcancel_pusher=top.pgroup.lebcancel_pusher
                    xygroup.append(xgroup)
                s.pgroups.append(xygroup)
            pxr.set_are_tiles_reallocated(i+1, self.ntilex,self.ntiley,self.ntilez,zeros((self.ntilex,self.ntiley,self.ntilez),dtype=dtype('i8')))
#        for i,s in enumerate(self.listofallspecies):
#            def ppzx(self,**kw):
#                for pg in self.pgroups:
#                   self._callppfunc(ppzx,pgroup=pg,**kw)

        if (self.l_debug): print(" Allocates antenna")
        for i in range(len(self.laser_antenna)):
            self.laser_antenna[i].initialize_virtual_particles(w3d)

        if (self.l_debug): print("End allocatefieldarraysPXR")

#            s.ppzx = ppzx

    def print_nptiles(self,ispecies):
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(ispecies, ix, iy, iz)
                    print ix,iy,iz,pxr.partn[0], pxr.partnmax
    def print_nptiles_sp0(self):
        s=self.listofallspecies[0]
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(1, ix, iy, iz)
                    print ix,iy,iz,pxr.partn, pxr.partnmax
                    print ix,iy,iz,s.pgroups[iz-1][iy-1][ix-1].nps, s.pgroups[iz-1][iy-1][ix-1].npmax
    def ppzx_ptiles(self,ispecies,ppg,colors=['black','blue','red','green'],msize=2):
        ncolor = len(colors)
        ic=0
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(ispecies, ix, iy, iz)
                    ppg(pxr.partx[:pxr.partn[0]],pxr.partz[:pxr.partn[0]],color=colors[ic%ncolor],msize=msize)
                    ic+=1

    def ppzx_ptiles_v2(self,ispecies,ppg,**kw):
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(ispecies, ix, iy, iz)
                    ppg(pxr.partx[:pxr.partn[0]],pxr.partz[:pxr.partn[0]],kwdict=kw)

    def push_e(self,dir=1.):
        """
        Electric field Maxwell solver
        """

        tdeb=MPI.Wtime()

        dt = dir*top.dt/self.ntsub
        if self.novercycle==1:
            if dir>0.:
                doit=True
            else:
                doit=False
        else:
            if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
                doit=True
            else:
                doit=False
        if doit:
            if self.l_verbose:print 'push_e',self,dt,top.it,self.icycle

            if self.l_pxr:
                f=self.fields
                l_pushe=False
                tdebcell=MPI.Wtime()
                if self.l_2dxz:
                    if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
                        lo = [0, 0]
                        hi = [f.nx, f.nz]
                        flo = [-f.nxguard, -f.nzguard]
                        fhi = [f.nx + f.nxguard, f.nz + f.nzguard]
                        # Warp field arrays have shape (nx, 1, nz) while PICSAR
                        # function pxrpush_em2d_evec takes fields with shape (nx, nz),
                        # so field arrays have to be squeezed.
                        pxr.pxrpush_em2d_evec(lo, hi, lo, hi, lo, hi,
                                              f.Ex.squeeze(), flo, fhi,
                                              f.Ey.squeeze(), flo, fhi,
                                              f.Ez.squeeze(), flo, fhi,
                                              f.Bx.squeeze(), flo, fhi,
                                              f.By.squeeze(), flo, fhi,
                                              f.Bz.squeeze(), flo, fhi,
                                              f.Jx.squeeze(), flo, fhi,
                                              f.Jy.squeeze(), flo, fhi,
                                              f.Jz.squeeze(), flo, fhi,
                                              clight**2*mu0*dt,
                                              clight**2*dt/f.dx*f.xcoefs[0],
                                              clight**2*dt/f.dy*f.ycoefs[0],
                                              clight**2*dt/f.dz*f.zcoefs[0])
                    else:
                        pxr.pxrpush_em2d_evec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                                  f.Jx,f.Jy,f.Jz,
                                                  clight**2*mu0*dt,
                                                  clight**2*dt/f.dx*f.xcoefs,
                                                  clight**2*dt/f.dy*f.ycoefs,
                                                  clight**2*dt/f.dz*f.zcoefs,
                                                  f.nx,f.ny,f.nz,
                                                  f.norderx,f.nordery,f.norderz,
                                                  f.nxguard,f.nyguard,f.nzguard,
                                                  0,0,0,f.l_nodalgrid)
                else:
                    if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
                        lo = [0, 0, 0]
                        hi = [f.nx, f.ny, f.nz]
                        flo = [-f.nxguard, -f.nyguard, -f.nzguard]
                        fhi = [f.nx + f.nxguard, f.ny + f.nyguard, f.nz + f.nzguard]
                        pxr.pxrpush_em3d_evec(lo, hi, lo, hi, lo, hi,
                                              f.Ex, flo, fhi,
                                              f.Ey, flo, fhi,
                                              f.Ez, flo, fhi,
                                              f.Bx, flo, fhi,
                                              f.By, flo, fhi,
                                              f.Bz, flo, fhi,
                                              f.Jx, flo, fhi,
                                              f.Jy, flo, fhi,
                                              f.Jz, flo, fhi,
                                              clight**2*mu0*dt,
                                              clight**2*dt/f.dx*f.xcoefs[0],
                                              clight**2*dt/f.dy*f.ycoefs[0],
                                              clight**2*dt/f.dz*f.zcoefs[0])
                    else:
                        pxr.pxrpush_em3d_evec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              f.Jx,f.Jy,f.Jz,
                                              clight**2*mu0*dt,
                                              clight**2*dt/f.dx*f.xcoefs,
                                              clight**2*dt/f.dy*f.ycoefs,
                                              clight**2*dt/f.dz*f.zcoefs,
                                              f.nx,f.ny,f.nz,
                                              f.norderx,f.nordery,f.norderz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
                tendcell=MPI.Wtime()
                pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
            else:
                l_pushe=True
            push_em3d_eef(self.block,dt,0,self.l_pushf,self.l_pushpot,l_pushe)
        if self.refinement is not None:
            self.__class__.__bases__[1].push_e(self.field_coarse,dir)

        tend=MPI.Wtime()
        self.time_stat_loc_array[7] += (tend-tdeb)

    def push_b_part_1(self,dir=1.):
      """
      Magnetic field solver
      """

      tdeb=MPI.Wtime()

      dt = dir*top.dt/self.ntsub
      if self.novercycle==1:
        if dir>0.:
          doit=True
        else:
          doit=False
      else:
        if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
          doit=True
        else:
          doit=False
      if doit:
        if self.l_verbose:print 'push_b part 1',self,dt,top.it,self.icycle,dir
        if self.l_pxr:
          tdebcell=MPI.Wtime()
          f=self.fields
          l_pushb=False
          if self.l_2dxz:
            if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
              if (f.stencil==0): # Yee solver
                lo = [0, 0]
                hi = [f.nx, f.nz]
                flo = [-f.nxguard, -f.nzguard]
                fhi = [f.nx + f.nxguard, f.nz + f.nzguard]
                # Warp field arrays have shape (nx, 1, nz) while PICSAR
                # function pxrpush_em2d_bvec takes fields with shape (nx, nz),
                # so field arrays have to be squeezed.
                pxr.pxrpush_em2d_bvec(lo, hi, lo, hi, lo, hi,
                                      f.Ex.squeeze(), flo, fhi,
                                      f.Ey.squeeze(), flo, fhi,
                                      f.Ez.squeeze(), flo, fhi,
                                      f.Bx.squeeze(), flo, fhi,
                                      f.By.squeeze(), flo, fhi,
                                      f.Bz.squeeze(), flo, fhi,
                                      0.5*dt/f.dx*f.xcoefs[0],
                                      0.5*dt/f.dy*f.ycoefs[0],
                                      0.5*dt/f.dz*f.zcoefs[0])
              elif (f.stencil==1): # Karkkainen solver
                lo = [0, 0]
                hi = [f.nx, f.nz]
                flo = [-f.nxguard, -f.nzguard]
                fhi = [f.nx + f.nxguard, f.nz + f.nzguard]
                # Warp field arrays have shape (nx, 1, nz) while PICSAR
                # function pxrpush_em2d_bvec_ckc takes fields with shape
                # (nx, nz), so field arrays have to be squeezed.
                pxr.pxrpush_em2d_bvec_ckc(lo, hi, lo, hi, lo, hi,
                                      f.Ex.squeeze(), flo, fhi,
                                      f.Ey.squeeze(), flo, fhi,
                                      f.Ez.squeeze(), flo, fhi,
                                      f.Bx.squeeze(), flo, fhi,
                                      f.By.squeeze(), flo, fhi,
                                      f.Bz.squeeze(), flo, fhi,
                                      0.5*dt/f.dx*f.xcoefs[0],
                                      0.5*dt/f.dy*f.ycoefs[0],
                                      0.5*dt/f.dz*f.zcoefs[0])
            else: #nth order solver >  2
              pxr.pxrpush_em2d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                          0.5*dt/f.dx*f.xcoefs,
                          0.5*dt/f.dy*f.ycoefs,
                          0.5*dt/f.dz*f.zcoefs,
                          f.nx,f.ny,f.nz,
                          f.norderx,f.nordery,f.norderz,
                          f.nxguard,f.nyguard,f.nzguard,
                          0,0,0,f.l_nodalgrid)
          else:
            if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
              if (f.stencil==0): # Yee solver
                lo = [0, 0, 0]
                hi = [f.nx, f.ny, f.nz]
                flo = [-f.nxguard, -f.nyguard, -f.nzguard]
                fhi = [f.nx + f.nxguard, f.ny + f.nyguard, f.nz + f.nzguard]
                pxr.pxrpush_em3d_bvec(lo, hi, lo, hi, lo, hi,
                                      f.Ex, flo, fhi,
                                      f.Ey, flo, fhi,
                                      f.Ez, flo, fhi,
                                      f.Bx, flo, fhi,
                                      f.By, flo, fhi,
                                      f.Bz, flo, fhi,
                                      0.5*dt/f.dx*f.xcoefs[0],
                                      0.5*dt/f.dy*f.ycoefs[0],
                                      0.5*dt/f.dz*f.zcoefs[0])
              elif (f.stencil==1): # Karkkainen solver
                lo = [0, 0, 0]
                hi = [f.nx, f.ny, f.nz]
                flo = [-f.nxguard, -f.nyguard, -f.nzguard]
                fhi = [f.nx + f.nxguard, f.ny + f.nyguard, f.nz + f.nzguard]
                pxr.pxrpush_em3d_bvec_ckc(lo, hi, lo, hi, lo, hi,
                                      f.Ex, flo, fhi,
                                      f.Ey, flo, fhi,
                                      f.Ez, flo, fhi,
                                      f.Bx, flo, fhi,
                                      f.By, flo, fhi,
                                      f.Bz, flo, fhi,
                                      0.5*dt/f.dx*f.xcoefs[0],
                                      0.5*dt/f.dy*f.ycoefs[0],
                                      0.5*dt/f.dz*f.zcoefs[0])
            else: #nth order solver >  2
              pxr.pxrpush_em3d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                          0.5*dt/f.dx*f.xcoefs,
                          0.5*dt/f.dy*f.ycoefs,
                          0.5*dt/f.dz*f.zcoefs,
                          f.nx,f.ny,f.nz,
                          f.norderx,f.nordery,f.norderz,
                          f.nxguard,f.nyguard,f.nzguard,
                          0,0,0,f.l_nodalgrid)
          tendcell=MPI.Wtime()
          pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
        else:
          l_pushb=True
        push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot,l_pushb)
      if self.refinement is not None:
        self.__class__.__bases__[1].push_b_part_1(self.field_coarse,dir)

      tend=MPI.Wtime()
      self.time_stat_loc_array[5] += (tend-tdeb)

    def push_b_part_2(self):
      """
      Magnetic field solver
      """

      tdeb=MPI.Wtime()

      if top.efetch[0] != 4 and (self.refinement is None):self.node2yee3d()
      dt = top.dt/self.ntsub
      if self.ntsub<1.:
        self.novercycle = nint(1./self.ntsub)
        self.icycle = (top.it-1)%self.novercycle
      else:
        self.novercycle = 1
        self.icycle = 0
      if self.icycle==0:
        if self.l_verbose:print 'push_b part 2',self,dt,top.it,self.icycle
        if self.l_pxr:
          f=self.fields
          l_pushb=False
          tdebcell=MPI.Wtime()
          if self.l_2dxz:
            if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
              if (f.stencil==0): # Yee solver
                lo = [0, 0]
                hi = [f.nx, f.nz]
                flo = [-f.nxguard, -f.nzguard]
                fhi = [f.nx + f.nxguard, f.nz + f.nzguard]
                # Warp field arrays have shape (nx, 1, nz) while PICSAR
                # function pxrpush_em2d_bvec takes fields with shape (nx, nz),
                # so field arrays have to be squeezed.
                pxr.pxrpush_em2d_bvec(lo, hi, lo, hi, lo, hi,
                                      f.Ex.squeeze(), flo, fhi,
                                      f.Ey.squeeze(), flo, fhi,
                                      f.Ez.squeeze(), flo, fhi,
                                      f.Bx.squeeze(), flo, fhi,
                                      f.By.squeeze(), flo, fhi,
                                      f.Bz.squeeze(), flo, fhi,
                                      0.5*dt/f.dx*f.xcoefs[0],
                                      0.5*dt/f.dy*f.ycoefs[0],
                                      0.5*dt/f.dz*f.zcoefs[0])
              elif (f.stencil==1): # Karkkainen solver
                lo = [0, 0]
                hi = [f.nx, f.nz]
                flo = [-f.nxguard, -f.nzguard]
                fhi = [f.nx + f.nxguard, f.nz + f.nzguard]
                # Warp field arrays have shape (nx, 1, nz) while PICSAR
                # function pxrpush_em2d_bvec_ckc takes fields with shape
                # (nx, nz), so field arrays have to be squeezed.
                pxr.pxrpush_em2d_bvec_ckc(lo, hi, lo, hi, lo, hi,
                                      f.Ex.squeeze(), flo, fhi,
                                      f.Ey.squeeze(), flo, fhi,
                                      f.Ez.squeeze(), flo, fhi,
                                      f.Bx.squeeze(), flo, fhi,
                                      f.By.squeeze(), flo, fhi,
                                      f.Bz.squeeze(), flo, fhi,
                                      0.5*dt/f.dx*f.xcoefs[0],
                                      0.5*dt/f.dy*f.ycoefs[0],
                                      0.5*dt/f.dz*f.zcoefs[0])
            else: #nth order solver >  2
              pxr.pxrpush_em2d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                          0.5*dt/f.dx*f.xcoefs,
                          0.5*dt/f.dy*f.ycoefs,
                          0.5*dt/f.dz*f.zcoefs,
                          f.nx,f.ny,f.nz,
                          f.norderx,f.nordery,f.norderz,
                          f.nxguard,f.nyguard,f.nzguard,
                          0,0,0,f.l_nodalgrid)
          else:
            if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
              if (f.stencil==0): # Yee solver
                lo = [0, 0, 0]
                hi = [f.nx, f.ny, f.nz]
                flo = [-f.nxguard, -f.nyguard, -f.nzguard]
                fhi = [f.nx + f.nxguard, f.ny + f.nyguard, f.nz + f.nzguard]
                pxr.pxrpush_em3d_bvec(lo, hi, lo, hi, lo, hi,
                                      f.Ex, flo, fhi,
                                      f.Ey, flo, fhi,
                                      f.Ez, flo, fhi,
                                      f.Bx, flo, fhi,
                                      f.By, flo, fhi,
                                      f.Bz, flo, fhi,
                                      0.5*dt/f.dx*f.xcoefs[0],
                                      0.5*dt/f.dy*f.ycoefs[0],
                                      0.5*dt/f.dz*f.zcoefs[0])
              elif (f.stencil==1): # Karkkainen solver
                lo = [0, 0, 0]
                hi = [f.nx, f.ny, f.nz]
                flo = [-f.nxguard, -f.nyguard, -f.nzguard]
                fhi = [f.nx + f.nxguard, f.ny + f.nyguard, f.nz + f.nzguard]
                pxr.pxrpush_em3d_bvec_ckc(lo, hi, lo, hi, lo, hi,
                                      f.Ex, flo, fhi,
                                      f.Ey, flo, fhi,
                                      f.Ez, flo, fhi,
                                      f.Bx, flo, fhi,
                                      f.By, flo, fhi,
                                      f.Bz, flo, fhi,
                                      0.5*dt/f.dx*f.xcoefs[0],
                                      0.5*dt/f.dy*f.ycoefs[0],
                                      0.5*dt/f.dz*f.zcoefs[0])
            else:  #nth order solver >  2
              pxr.pxrpush_em3d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                          0.5*dt/f.dx*f.xcoefs,
                          0.5*dt/f.dy*f.ycoefs,
                          0.5*dt/f.dz*f.zcoefs,
                          f.nx,f.ny,f.nz,
                          f.norderx,f.nordery,f.norderz,
                          f.nxguard,f.nyguard,f.nzguard,
                          0,0,0,f.l_nodalgrid)
          tendcell=MPI.Wtime()
          pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
        else:
          l_pushb=True
        push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot,l_pushb)
      if self.refinement is not None:
        self.__class__.__bases__[1].push_b_part_2(self.field_coarse)

      # Time statistics
      tend=MPI.Wtime()
      self.time_stat_loc_array[5] += (tend-tdeb)

    def push_spectral_psaotd(self):
        """
        PSAOTD Maxwell solver
        """
        if self.l_pxr:
          tdebcell=MPI.Wtime()

        if(self.full_pxr):
          self.solve_maxwell_full_pxr()
          return
        if top.efetch[0] != 4 and (self.refinement is None) and not self.l_nodalgrid:self.node2yee3d()

        if self.ntsub==inf:
          self.GPSTDMaxwell.fields['rhoold']=self.fields.Rhoold
          self.fields.Rho=self.fields.Rhoarray[...,0]
          self.GPSTDMaxwell.fields['rhonew']=self.fields.Rho
        else:
          if self.l_pushf:
        #                self.fields.Rho=self.fields.Rhoarray[...,0]
            self.GPSTDMaxwell.fields['rhoold']=self.fields.Rhoold.copy()
            self.GPSTDMaxwell.fields['rhonew']=self.fields.Rho.copy()
            self.GPSTDMaxwell.fields['drho']=self.fields.Rho-self.fields.Rhoold

        self.GPSTDMaxwell.fields['jx']=self.fields.Jx
        self.GPSTDMaxwell.fields['jy']=self.fields.Jy
        self.GPSTDMaxwell.fields['jz']=self.fields.Jz

        self.GPSTDMaxwell.push_fields()

        b=self.block

        # --- sides
        if b.xlbnd==openbc:self.xlPML.push()
        if b.xrbnd==openbc:self.xrPML.push()
        if b.ylbnd==openbc:self.ylPML.push()
        if b.yrbnd==openbc:self.yrPML.push()
        if b.zlbnd==openbc:self.zlPML.push()
        if b.zrbnd==openbc:self.zrPML.push()

        # --- edges
        if(b.xlbnd==openbc and b.ylbnd==openbc):self.xlylPML.push()
        if(b.xrbnd==openbc and b.ylbnd==openbc):self.xrylPML.push()
        if(b.xlbnd==openbc and b.yrbnd==openbc):self.xlyrPML.push()
        if(b.xrbnd==openbc and b.yrbnd==openbc):self.xryrPML.push()
        if(b.xlbnd==openbc and b.zlbnd==openbc):self.xlzlPML.push()
        if(b.xrbnd==openbc and b.zlbnd==openbc):self.xrzlPML.push()
        if(b.xlbnd==openbc and b.zrbnd==openbc):self.xlzrPML.push()
        if(b.xrbnd==openbc and b.zrbnd==openbc):self.xrzrPML.push()
        if(b.ylbnd==openbc and b.zlbnd==openbc):self.ylzlPML.push()
        if(b.yrbnd==openbc and b.zlbnd==openbc):self.yrzlPML.push()
        if(b.ylbnd==openbc and b.zrbnd==openbc):self.ylzrPML.push()
        if(b.yrbnd==openbc and b.zrbnd==openbc):self.yrzrPML.push()

        # --- corners
        if(b.xlbnd==openbc and b.ylbnd==openbc and b.zlbnd==openbc):self.xlylzlPML.push()
        if(b.xrbnd==openbc and b.ylbnd==openbc and b.zlbnd==openbc):self.xrylzlPML.push()
        if(b.xlbnd==openbc and b.yrbnd==openbc and b.zlbnd==openbc):self.xlyrzlPML.push()
        if(b.xrbnd==openbc and b.yrbnd==openbc and b.zlbnd==openbc):self.xryrzlPML.push()
        if(b.xlbnd==openbc and b.ylbnd==openbc and b.zrbnd==openbc):self.xlylzrPML.push()
        if(b.xrbnd==openbc and b.ylbnd==openbc and b.zrbnd==openbc):self.xrylzrPML.push()
        if(b.xlbnd==openbc and b.yrbnd==openbc and b.zrbnd==openbc):self.xlyrzrPML.push()
        if(b.xrbnd==openbc and b.yrbnd==openbc and b.zrbnd==openbc):self.xryrzrPML.push()

        #    if em.pml_method==2:
        #      self.fields.spectral=0
        #      scale_em3d_bnd_fields(self.block,top.dt,self.l_pushf)
        #      self.fields.spectral=1

        if self.boris_cor:
          self.boris_correction()
        if self.l_pxr:
          tendcell=MPI.Wtime()
          pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
          self.time_stat_loc_array[7] += (tendcell-tdebcell)

    def current_cor_spectral(self):
        """
        Current spectral correction
        """

        if self.l_pxr:
            tdebcell=MPI.Wtime()

        j=1j      # imaginary number
        emK = self.FSpace
        em = self
        f = self.fields
        ixl,ixu,iyl,iyu,izl,izu = emK.get_ius()

        fields_shape = [ixu-ixl,iyu-iyl,izu-izl]

        if emK.planj_rfftn is None:
          emK.planj_rfftn= emK.create_plan_rfftn(np.asarray(fields_shape))
        if emK.planj_irfftn is None:
          emK.planj_irfftn= emK.create_plan_irfftn(np.asarray(fields_shape))

        self.wrap_periodic_BC([f.Rho,f.Rhoold_local,f.Jx,f.Jy,f.Jz])

        if emK.nx>1:JxF = emK.rfftn(squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)
        if emK.ny>1:JyF = emK.rfftn(squeeze(f.Jy[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)
        if emK.nz>1:JzF = emK.rfftn(squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)

        em.dRhoodtF = emK.rfftn(squeeze((f.Rho-f.Rhoold_local)[ixl:ixu,iyl:iyu,izl:izu]/top.dt),plan=emK.planj_rfftn)

        # --- get longitudinal J
        divJ = 0.
        if emK.nx>1:divJ += emK.kxmn*JxF
        if emK.ny>1:divJ += emK.kymn*JyF
        if emK.nz>1:divJ += emK.kzmn*JzF

        if emK.nx>1:
          Jxl = emK.kxpn*divJ
        if emK.ny>1:
          Jyl = emK.kypn*divJ
        if emK.nz>1:
          Jzl = emK.kzpn*divJ

        # --- get transverse J
        if emK.nx>1:
          Jxt = JxF-Jxl
        if emK.ny>1:
          Jyt = JyF-Jyl
        if emK.nz>1:
          Jzt = JzF-Jzl

        if emK.nx>1:
          Jxl = j*em.dRhoodtF*emK.kxpn/emK.kmag
        if emK.ny>1:
          Jyl = j*em.dRhoodtF*emK.kypn/emK.kmag
        if emK.nz>1:
          Jzl = j*em.dRhoodtF*emK.kzpn/emK.kmag

        if emK.nx>1:
          JxF = Jxt+Jxl
        if emK.ny>1:
          JyF = Jyt+Jyl
        if emK.nz>1:
          JzF = Jzt+Jzl


        if emK.nx>1:
          Jx = emK.irfftn(JxF, np.asarray(np.shape(squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]))), plan=emK.planj_irfftn, field_out=squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]))
          Jx.resize(fields_shape)
          f.Jx[ixl:ixu,iyl:iyu,izl:izu] = Jx.real
        if emK.ny>1:
          Jy = emK.irfftn(JyF, np.asarray(np.shape(squeeze(f.Jy[ixl:ixu,iyl:iyu,izl:izu]))), plan=emK.planj_irfftn, field_out=squeeze(f.Jy[ixl:ixu,iyl:iyu,izl:izu]))
          Jy.resize(fields_shape)
          f.Jy[ixl:ixu,iyl:iyu,izl:izu] = Jy.real
        if emK.nz>1:
          Jz = emK.irfftn(JzF, np.asarray(np.shape(squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]))), plan=emK.planj_irfftn, field_out=squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]))
          Jz.resize(fields_shape)
          f.Jz[ixl:ixu,iyl:iyu,izl:izu] = Jz.real

        # Time statistics
        if self.l_pxr:
          tendcell=MPI.Wtime()
          pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
          self.time_stat_loc_array[16] += (tendcell-tdebcell)


    def exchange_e(self,dir=1.):
        """
        Electric field boundary conditions
        """

        t0 = MPI.Wtime()
        if self.novercycle==1:
            if dir>0.:
                doit=True
            else:
                doit=False
        else:
            if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
                doit=True
            else:
                doit=False
        if doit:
            if (self.l_pxr and self.full_pxr):
                pxr.efield_bcs()
            else:
                em3d_exchange_e(self.block)
        if self.refinement is not None:
            self.__class__.__bases__[1].exchange_e(self.field_coarse)

        t1 = MPI.Wtime()
        self.time_stat_loc_array[8] += (t1-t0)

    def exchange_b(self,dir=1.):
        """
        Magnetic field boundary conditions
        """

        t0 = MPI.Wtime()

        if self.novercycle==1:
            if dir>0.:
                doit=True
            else:
                doit=False
        else:
            if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
                doit=True
            else:
                doit=False
        if doit:
            if self.l_verbose:print 'exchange_b',self,top.it,self.icycle
            if (self.l_pxr and self.full_pxr):
                pxr.bfield_bcs()
            else:
                em3d_exchange_b(self.block)
        if self.refinement is not None:
            self.__class__.__bases__[1].exchange_b(self.field_coarse,dir)

        t1 = MPI.Wtime()
        self.time_stat_loc_array[6] += (t1-t0)


    def solve_maxwell_full_pxr(self):

        """ full Maxwell push in pxr"""
        pxr.rho = self.fields.Rho
        pxr.rhoold = self.fields.Rhoold
        pxr.jx = self.fields.Jx
        pxr.jy = self.fields.Jy
        pxr.jz = self.fields.Jz

        if(self.l_debug):print("begin solve maxwell full pxr")
        if(pxr.fftw_with_mpi):
          pxr.get_ffields_mpi_lb()
        else:
          pxr.get_ffields()
        if(pxr.g_spectral):
          pxr.multiply_mat_vector(pxr.nmatrixes)
        else:
          if(pxr.c_dim == 2):
            pxr.push_psaotd_ebfielfs_2d()
          else:
            pxr.push_psaotd_ebfielfs_3d()
        if(pxr.fftw_with_mpi):
          pxr.get_fields_mpi_lb()
        else:
          pxr.get_fields()
        if(pxr.absorbing_bcs):
          pxr.field_damping_bcs()
          #pxr.merge_fields()
        if(self.l_debug):print("end solve maxwell full pxr")

    def move_cells(self,n, coord):
        # move the boundaries of the box along the coord axis
        # in case of moving window along the coordinate coord.
        #coord = 'x', 'y', 'z'
        if(self.full_pxr and self.absorbing_bcs_pxr):

          # when using full pxr mode with absorbing_bcs, the moving window needs
          # to be applied to the splitted fields of pxr

          save_ntimes = self.block.core.yf.ntimes
          if   coord=='x': shift_em3dblock_ncells_x(self.block,n)
          elif coord=='y': shift_em3dblock_ncells_y(self.block,n)
          elif coord=='z': shift_em3dblock_ncells_z(self.block,n)
          # sets ntimes to 0 to force the mv window to act only once on the currents and densities
          self.block.core.yf.ntimes = 0

          s1 = self.block.core.yf.Ex
          s2 = self.block.core.yf.Ey
          s3 = self.block.core.yf.Ez
          s4 = self.block.core.yf.Bx
          s5 = self.block.core.yf.By
          s6 = self.block.core.yf.Bz

          self.block.core.yf.Ex = self.exy_pxr
          self.block.core.yf.Ey = self.eyx_pxr
          self.block.core.yf.Ez = self.ezx_pxr
          self.block.core.yf.Bx = self.bxy_pxr
          self.block.core.yf.By = self.byx_pxr
          self.block.core.yf.Bz = self.bzx_pxr


          if   coord=='x': shift_em3dblock_ncells_x(self.block,n)
          elif coord=='y': shift_em3dblock_ncells_y(self.block,n)
          elif coord=='z': shift_em3dblock_ncells_z(self.block,n)

          self.block.core.yf.Ex = self.exz_pxr
          self.block.core.yf.Ey = self.eyz_pxr
          self.block.core.yf.Ez = self.ezy_pxr
          self.block.core.yf.Bx = self.bxz_pxr
          self.block.core.yf.By = self.byz_pxr
          self.block.core.yf.Bz = self.bzy_pxr

          if   coord=='x': shift_em3dblock_ncells_x(self.block,n)
          elif coord=='y': shift_em3dblock_ncells_y(self.block,n)
          elif coord=='z': shift_em3dblock_ncells_z(self.block,n)

          self.block.core.yf.Ex = s1
          self.block.core.yf.Ey = s2
          self.block.core.yf.Ez = s3
          self.block.core.yf.Bx = s4
          self.block.core.yf.By = s5
          self.block.core.yf.Bz = s6
          self.block.core.yf.ntimes = save_ntimes
        else:
          if   coord=='x': shift_em3dblock_ncells_x(self.block,n)
          elif coord=='y': shift_em3dblock_ncells_y(self.block,n)
          elif coord=='z': shift_em3dblock_ncells_z(self.block,n)

        listtoshift = [(self,'%s_grid' %(coord) ),
                       (self,'%smmin'  %(coord) ),
                       (self,'%smmax' %(coord) ),
                       (self,'%smminlocal'%(coord) ),
                       (self,'%smmaxlocal'%(coord) ),
                       (self.fields,'%smin'%(coord) ),
                       (self.fields,'%smax'%(coord) ),
                       (self.block,'%smin'%(coord) ),
                       (self.block,'%smax'%(coord) ),
                       (w3d,'%smmin'%(coord) ),
                       (w3d,'%smmax'%(coord) ),
                       (w3d,'%smminp'%(coord) ),
                       (w3d,'%smmaxp'%(coord) ),
                       (w3d,'%smminlocal'%(coord) ),
                       (w3d,'%smmaxlocal'%(coord) ),
                       (w3d,'%smminglobal'%(coord) ),
                       (w3d,'%smmaxglobal'%(coord) ),
                       (top,'%spmin'%(coord) ),
                       (top,'%spmax'%(coord) ),
                       (top,'%spminlocal'%(coord) ),
                       (top,'%spmaxlocal'%(coord) )]

        if   coord=='x': increment=self.dx
        elif coord=='y': increment=self.dy
        elif coord=='z': increment=self.dz

        for (coord_object,coord_attribute) in listtoshift:
            # loop equivalent to self.incrementposition(coord_object.coord_attribute, increment, n)
            # for each tupple in listtoshift
            coordtoshift=getattr(coord_object,coord_attribute)
            setattr(coord_object,coord_attribute,self.incrementposition(coordtoshift,increment,n))






    def step(self,n=1,freq_print=10,lallspecl=0,stdout_stat=10):
      """
      This function performs a range of Particle-In-Cell iterations

      Inputs:
      - n: number of iterations
      - freq_print: print frequency
      """

      if (self.l_debug): print("Call step")

      t0=MPI.Wtime()
      tdeb=MPI.Wtime()

      for i in range(n):
          if(me==0):
              if top.it%freq_print==0:print 'it = %g time = %g'%(top.it,top.time)
          if lallspecl:
              l_first=l_last=1
          else:
              if i==0:
                  l_first=1
              else:
                  l_first=0
              if i==n-1:
                  l_last=1
              else:
                  l_last=0
          self.onestep(l_first,l_last)

          if(l_pxr & (top.it%stdout_stat==0) & (pxr.rank==0)):
              tend=MPI.Wtime()
              mpi_time_per_stat=(tend-tdeb)
              tdeb=MPI.Wtime()
              print("time/stdout_stat (s)",mpi_time_per_stat)

      # Total time spend in the kernel
      tend = MPI.Wtime()
      self.total_kernel_time = (tend-t0)

      if (self.l_debug): print("End step")



    def onestep(self,l_first,l_last):
        """
        Perform a single particle-in-cell step
        """

        if (self.l_debug): print("Call onestep")

        # --- Iteration number
        pxr.it = top.it

        # --- call beforestep functions
        if (self.l_debug): print("Call beforestep functions")
        callbeforestepfuncs.callfuncsinlist()
        #top.zgrid+=top.vbeamfrm*top.dt
        #top.zbeam=top.zgrid

        # --- gather fields from grid to particles
        if (self.l_debug): print("Call Field gathering and particle push")
#        w3d.pgroupfsapi = top.pgroup
#        for js in range(top.pgroup.ns):
#          self.fetcheb(js)
        if l_pxr:
            tdebpart=0.
            tendpart=0.
            tdebfield=0.
            tendfield=0.
            tdebcell=0.
            tendcell=0.
            tdeb=MPI.Wtime()
            pxr.local_time_part=0.
            pxr.local_time_cell=0.
        # --- push
        if l_first:
            if l_pxr:
                # Particle pusher
                tdebpart=MPI.Wtime()
                if (self.l_debug): print("Call record_old_positions()")
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.record_old_positions(0,pg)
                if (self.l_debug): print("Call pxr.pxrpush_particles_part2()")
                pxr.pxrpush_particles_part2()
                tendpart=MPI.Wtime()
                pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
                self.time_stat_loc_array[0] += (tendpart-tdebpart)

                # Particle boundary consitions
                #pxr.particle_bcs_2d()
                tdebpart=MPI.Wtime()
                if (self.l_debug): print("Call pxr.particle_bcs()")
                pxr.particle_bcs()
                tendpart=MPI.Wtime()
                self.time_stat_loc_array[1] += (tendpart-tdebpart)

                #for i,s in enumerate(self.listofallspecies):
                #    for pg in s.flatten(s.pgroups):
                #        particleboundaries3d(pg,-1,False)
                #pxr.particle_bcs_tiles()
                if (self.l_debug): print("Call aliasparticlearrays()")
                aliasparticlearrays()
            else:
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.push_velocity_second_half(0,pg)
                        self.record_old_positions(0,pg)
                        self.push_positions(0,pg)
                        particleboundaries3d(pg,-1,False)
        else:
            if l_pxr:
                # Particle pusher
                if (self.l_debug): print("Call pxr.field_gathering_plus_particle_pusher()")
                tdebpart=MPI.Wtime()
                #pxr.push_particles()
                if (self.l_debug): print("Call record_old_positions()")
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.record_old_positions(0,pg)
                pxr.field_gathering_plus_particle_pusher()
                tendpart=MPI.Wtime()
                pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
                self.time_stat_loc_array[0] += (tendpart-tdebpart)

                # Particle boundary conditions
                if (self.l_debug): print("Call pxr.particle_bcs()")
                tdebpart=MPI.Wtime()
                pxr.particle_bcs()
                tendpart=MPI.Wtime()
                self.time_stat_loc_array[1] += (tendpart-tdebpart)


                #for i,s in enumerate(self.listofallspecies):
                #    for pg in s.flatten(s.pgroups):
                #        particleboundaries3d(pg,-1,False)
                #pxr.particle_bcs_tiles()
                aliasparticlearrays()

            else:
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        tendpart=MPI.Wtime()
                        self.push_velocity_full(0,pg)
                        self.record_old_positions(0,pg)
                        self.push_positions(0,pg)
                        tendpart=MPI.Wtime()
                        self.time_stat_loc_array[0] += (tendpart-tdebpart)

                        # Particle boundary conditions
                        tdebpart=MPI.Wtime()
                        particleboundaries3d(pg,-1,False)

                        tendpart=MPI.Wtime()
                        self.time_stat_loc_array[1] += (tendpart-tdebpart)


        # --- Particle sorting
        if (self.l_debug): print("Call Particle Sorting")
        if l_pxr:
          if ((self.sorting.activated)and(top.it>=0)):
            tdebpart=MPI.Wtime()
            pxr.particle_sorting_sub()
            tendpart=MPI.Wtime()
            self.time_stat_loc_array[10] += (tendpart-tdebpart)


        pgroups = []
        for i,s in enumerate(self.listofallspecies):
            pgroups+=s.flatten(s.pgroups)
        self.pgroups = pgroups
#        self.loadsource(pgroups=pgroups)
        #tdebpart=MPI.Wtime()

        inject3d(1, top.pgroup)

        # Call user-defined injection routines
        if (self.l_debug): print("Call user-defined injection routines")
        userinjection.callfuncsinlist()

        # --- call beforeloadrho functions
        if (self.l_debug): print("Call beforeloadrho functions")
        beforeloadrho.callfuncsinlist()

        xgrid=w3d.xmmin-pxr.xmin
        ygrid=w3d.ymmin-pxr.ymin
        zgrid=w3d.zmmin-pxr.zmin
        if (xgrid != 0. or ygrid!=0. or zgrid !=0.):
            pxr.pxr_move_sim_boundaries(xgrid,ygrid,zgrid)
            pxr.particle_bcs()
            aliasparticlearrays()


        if (self.l_debug): print("Call loadrho")
        self.loadrho(pgroups=pgroups)
        if (self.l_debug): print("Call loadj")
        self.loadj(pgroups=pgroups)

        # Moving window

        #tendpart=MPI.Wtime()
        #pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
#        self.solve2ndhalf()

        #tdebcell=MPI.Wtime()
        # --- dosolve
        # Current deposition + Maxwell

        if (self.l_debug): print("Call dosolve")

        self.dosolve()
        #tendcell=MPI.Wtime()
        #pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)

        tdebpart=MPI.Wtime()
        if not l_pxr or (l_pxr and pxr.fieldgathe<0):
            for i,s in enumerate(self.listofallspecies):
                for pg in s.flatten(s.pgroups):
                    w3d.pgroupfsapi = pg
                    self.fetcheb(0,pg)

        if l_last:
            if l_pxr:
                if (self.l_debug): print("Call pxr.pxrpush_particles_part1()")
                pxr.pxrpush_particles_part1()
            else:
                for pg in s.flatten(s.pgroups):
                    w3d.pgroupfsapi = pg
                    self.push_velocity_first_half(0,pg)

        tendpart=MPI.Wtime()
        if l_pxr:pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
        self.time_stat_loc_array[0] += (tendpart-tdebpart)

        # --- update time, time counter
        top.time+=top.dt
        if top.it%top.nhist==0:
#           zmmnt()
           minidiag(top.it,top.time,top.lspecial)
        top.it+=1

        #Load balance every dlb_freq time step
        if (self.l_debug): print("Call Load balance")
        if (l_pxr & (self.dload_balancing & (top.it%self.dlb_freq==0))):
            pxr.mpitime_per_it=pxr.local_time_part+pxr.local_time_cell
            pxr.get_max_time_per_it()
            pxr.get_min_time_per_it()
            ## --- Compute time per part and per cell
            pxr.compute_time_per_part()
            pxr.compute_time_per_cell()
            imbalance=(pxr.max_time_per_it-pxr.min_time_per_it)/pxr.min_time_per_it*100.
            if (imbalance>self.dlb_threshold):
              if (self.l_2dxz):
                  self.load_balance_2d(str(imbalance)+"%")
              else:
                    self.load_balance_3d(str(imbalance)+"%")
        # Try to Load balance at init
        if ((top.it==self.it_dlb_init) & self.dlb_at_init & self.dload_balancing):
          pxr.mpitime_per_it=pxr.local_time_part+pxr.local_time_cell
          pxr.get_max_time_per_it()
          pxr.get_min_time_per_it()
          ## --- Compute time per part and per cell
          pxr.compute_time_per_part()
          pxr.compute_time_per_cell()
          if (self.l_2dxz):
              self.load_balance_2d('Init')
          else:
            self.load_balance_3d('Init')

        # PXr custom outputs mpi-io
        if (self.l_debug): print("Call PXR custom outputs mpi-io")
        if(l_pxr & self.l_output_grid & (top.it % self.l_output_freq ==0)):
          self.output_pxr(top.it)

        xgrid=w3d.xmmin-pxr.xmin
        ygrid=w3d.ymmin-pxr.ymin
        zgrid=w3d.zmmin-pxr.zmin
        if (xgrid != 0. or ygrid!=0. or zgrid !=0.):
            pxr.pxr_move_sim_boundaries(xgrid,ygrid,zgrid)
            pxr.particle_bcs()
            aliasparticlearrays()

        # --- call afterstep functions
        if (self.l_debug): print("Call callafterstepfuncs.callfuncsinlist()")
        callafterstepfuncs.callfuncsinlist()



    def load_balance_3d(self,imbalance):
        """
        Load balance between MPI domains in 3D
        """
        if (l_pxr):
            tdeb = MPI.Wtime()

            ## --- Compute time per part and per cell
            pxr.compute_time_per_part()
            pxr.compute_time_per_cell()

            ## --- Compute new split along each dimension
            pxr.compute_new_split(pxr.global_time_per_part,pxr.global_time_per_cell,pxr.nx_global,pxr.ny_global,pxr.nz_global,
                              pxr.new_cell_x_min,pxr.new_cell_x_max,pxr.new_cell_y_min,pxr.new_cell_y_max,
                              pxr.new_cell_z_min,pxr.new_cell_z_max,pxr.nprocx,pxr.nprocy,pxr.nprocz)
            isnewsplit=sum(pxr.cell_x_min-pxr.new_cell_x_min)+sum(pxr.cell_x_max-pxr.new_cell_x_max)+ \
                     sum(pxr.cell_y_min-pxr.new_cell_y_min)+sum(pxr.cell_y_max-pxr.new_cell_y_max)+ \
                     sum(pxr.cell_z_min-pxr.new_cell_z_min)+sum(pxr.cell_z_max-pxr.new_cell_z_max)
            if (isnewsplit==0):
              if(pxr.rank==0):
                  print("Optimal load balancing already achieved by current implementation")
            else:
              if(pxr.rank==0):
                print("trying to load balance the simulation, imbalance=", imbalance)
                ## --- Compute limits for all procs
                ix1old=np.zeros(pxr.nproc,dtype="i8"); ix2old=np.zeros(pxr.nproc,dtype="i8")
                iy1old=np.zeros(pxr.nproc,dtype="i8"); iy2old=np.zeros(pxr.nproc,dtype="i8")
                iz1old=np.zeros(pxr.nproc,dtype="i8"); iz2old=np.zeros(pxr.nproc,dtype="i8")
                ix1new=np.zeros(pxr.nproc,dtype="i8"); ix2new=np.zeros(pxr.nproc,dtype="i8")
                iy1new=np.zeros(pxr.nproc,dtype="i8"); iy2new=np.zeros(pxr.nproc,dtype="i8")
                iz1new=np.zeros(pxr.nproc,dtype="i8"); iz2new=np.zeros(pxr.nproc,dtype="i8")

                pxr.get_1darray_proclimits(ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,
                                        pxr.cell_x_min,pxr.cell_y_min,pxr.cell_z_min,
                                        pxr.cell_x_max,pxr.cell_y_max,pxr.cell_z_max,
                                        pxr.nprocx, pxr.nprocy, pxr.nprocz, pxr.nproc,
                                        top.lcomm_cartesian)
                pxr.get_1darray_proclimits(ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,
                                        pxr.new_cell_x_min,pxr.new_cell_y_min,pxr.new_cell_z_min,
                                        pxr.new_cell_x_max,pxr.new_cell_y_max,pxr.new_cell_z_max,
                                        pxr.nprocx, pxr.nprocy, pxr.nprocz, pxr.nproc,top.lcomm_cartesian)
                ## --- Compute new sizes for grid arrays
                nx_new=pxr.new_cell_x_max[pxr.x_coords]-pxr.new_cell_x_min[pxr.x_coords]+1
                ny_new=pxr.new_cell_y_max[pxr.y_coords]-pxr.new_cell_y_min[pxr.y_coords]+1
                nz_new=pxr.new_cell_z_max[pxr.z_coords]-pxr.new_cell_z_min[pxr.z_coords]+1

                ## --- Remap field arrays
                # -- Ex
                ex_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(ex_new,nx_new,ny_new,nz_new,
                                          pxr.ex,pxr.nx,pxr.ny,pxr.nz,
                                          pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                          ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                          ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.ex=ex_new
                # -- Ey
                ey_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(ey_new,nx_new,ny_new,nz_new,
                                          pxr.ey,pxr.nx,pxr.ny,pxr.nz,
                                          pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                          ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                          ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.ey=ey_new
                # -- Ez
                ez_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(ez_new,nx_new,ny_new,nz_new,
                                          pxr.ez,pxr.nx,pxr.ny,pxr.nz,
                                          pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                          ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                          ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.ez=ez_new
                # -- Bx
                bx_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(bx_new,nx_new,ny_new,nz_new,
                                          pxr.bx,pxr.nx,pxr.ny,pxr.nz,
                                          pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                          ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                          ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.bx=bx_new
                # -- By
                by_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(by_new,nx_new,ny_new,nz_new,
                                          pxr.by,pxr.nx,pxr.ny,pxr.nz,
                                          pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                          ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                          ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.by=by_new
                # -- Bz
                bz_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(bz_new,nx_new,ny_new,nz_new,
                                          pxr.bz,pxr.nx,pxr.ny,pxr.nz,
                                          pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                          ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                          ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.bz=bz_new
                ## -- Reallocate current arrays
                # Currents are recomputed each iteration so no need to exchange them
                jx_new=zeros((nx_new+2*pxr.nxjguards+1,ny_new+2*pxr.nyjguards+1,nz_new+2*pxr.nzjguards+1),order='F')
                jy_new=zeros((nx_new+2*pxr.nxjguards+1,ny_new+2*pxr.nyjguards+1,nz_new+2*pxr.nzjguards+1),order='F')
                jz_new=zeros((nx_new+2*pxr.nxjguards+1,ny_new+2*pxr.nyjguards+1,nz_new+2*pxr.nzjguards+1),order='F')
                pxr.jx=jx_new
                pxr.jy=jy_new
                pxr.jz=jz_new

                # Update pxr new array dimensions
                pxr.nx=nx_new
                pxr.ny=ny_new
                pxr.nz=nz_new
                pxr.nx_grid=pxr.nx+1
                pxr.ny_grid=pxr.ny+1
                pxr.nz_grid=pxr.nz+1

                # Test if domain has been resized - used for particle remaping
                isnewdom=pxr.cell_x_min[pxr.x_coords]-pxr.new_cell_x_min[pxr.x_coords]+pxr.cell_x_max[pxr.x_coords]-pxr.new_cell_x_max[pxr.x_coords]+ \
                pxr.cell_y_min[pxr.y_coords]-pxr.new_cell_y_min[pxr.y_coords]+pxr.cell_y_max[pxr.y_coords]-pxr.new_cell_y_max[pxr.y_coords]+ \
                pxr.cell_z_min[pxr.z_coords]-pxr.new_cell_z_min[pxr.z_coords]+pxr.cell_z_max[pxr.z_coords]-pxr.new_cell_z_max[pxr.z_coords]

                # Update new subdomain index arrays
                pxr.cell_x_min=pxr.new_cell_x_min
                pxr.cell_x_max=pxr.new_cell_x_max
                pxr.cell_y_min=pxr.new_cell_y_min
                pxr.cell_y_max=pxr.new_cell_y_max
                pxr.cell_z_min=pxr.new_cell_z_min
                pxr.cell_z_max=pxr.new_cell_z_max
                pxr.nx_global_grid_min = pxr.cell_x_min[pxr.x_coords]
                pxr.nx_global_grid_max = pxr.cell_x_max[pxr.x_coords]+1
                pxr.ny_global_grid_min = pxr.cell_y_min[pxr.y_coords]
                pxr.ny_global_grid_max = pxr.cell_y_max[pxr.y_coords]+1
                pxr.nz_global_grid_min = pxr.cell_z_min[pxr.z_coords]
                pxr.nz_global_grid_max = pxr.cell_z_max[pxr.z_coords]+1


                # Update global simulation axis
                pxr.compute_simulation_axis()

                # Set new min and max for local domain
                pxr.x_min_local = pxr.x_grid_mins[pxr.x_coords]
                pxr.x_max_local = pxr.x_grid_maxs[pxr.x_coords]
                pxr.y_min_local = pxr.y_grid_mins[pxr.y_coords]
                pxr.y_max_local = pxr.y_grid_maxs[pxr.y_coords]
                pxr.z_min_local = pxr.z_grid_mins[pxr.z_coords]
                pxr.z_max_local = pxr.z_grid_maxs[pxr.z_coords]
                pxr.x_grid_min_local=pxr.x_min_local
                pxr.x_grid_max_local=pxr.x_max_local
                pxr.y_grid_min_local=pxr.y_min_local
                pxr.y_grid_max_local=pxr.y_max_local
                pxr.z_grid_min_local=pxr.z_min_local
                pxr.z_grid_max_local=pxr.z_max_local

                ##--- Alias WARP grid arrays on pxr new arrays
                self.nxlocal=pxr.nx
                self.nylocal=pxr.ny
                self.nzlocal=pxr.nz
                self.ymminlocal = pxr.y_min_local
                self.zmminlocal = pxr.z_min_local
                self.fields.xmin = pxr.x_min_local
                self.fields.xmax = pxr.x_max_local
                self.fields.ymin = pxr.y_min_local
                self.fields.ymax = pxr.y_max_local
                self.fields.zmin = pxr.z_min_local
                self.fields.zmax = pxr.z_max_local

                # Udpate domain decomposition in WARP
                top.fsdecomp.nx=pxr.cell_x_max-pxr.cell_x_min+1
                top.fsdecomp.ny=pxr.cell_y_max-pxr.cell_y_min+1
                top.fsdecomp.nz=pxr.cell_z_max-pxr.cell_z_min+1
                top.fsdecomp.ix=pxr.cell_x_min
                top.fsdecomp.iy=pxr.cell_y_min
                top.fsdecomp.iz=pxr.cell_z_min
                top.fsdecomp.xmin=pxr.cell_x_min*pxr.dx
                top.fsdecomp.xmax=(pxr.cell_x_max+1)*pxr.dx
                top.fsdecomp.ymin=pxr.cell_y_min*pxr.dy
                top.fsdecomp.ymax=(pxr.cell_y_max+1)*pxr.dy
                top.fsdecomp.zmin=pxr.cell_z_min*pxr.dz
                top.fsdecomp.zmax=(pxr.cell_z_max+1)*pxr.dz
                top.ppdecomp.nx=pxr.cell_x_max-pxr.cell_x_min+1
                top.ppdecomp.ny=pxr.cell_y_max-pxr.cell_y_min+1
                top.ppdecomp.nz=pxr.cell_z_max-pxr.cell_z_min+1
                top.ppdecomp.ix=pxr.cell_x_min
                top.ppdecomp.iy=pxr.cell_y_min
                top.ppdecomp.iz=pxr.cell_z_min
                top.ppdecomp.xmin=pxr.cell_x_min*pxr.dx
                top.ppdecomp.xmax=(pxr.cell_x_max+1)*pxr.dx
                top.ppdecomp.ymin=pxr.cell_y_min*pxr.dy
                top.ppdecomp.ymax=(pxr.cell_y_max+1)*pxr.dy
                top.ppdecomp.zmin=pxr.cell_z_min*pxr.dz
                top.ppdecomp.zmax=(pxr.cell_z_max+1)*pxr.dz

                # Reallocate warp arrays
                self.allocatefieldarrays()
                # Alias newly allocated arrays on WARP structure
                self.fields.Ex=pxr.ex
                self.fields.Ey=pxr.ey
                self.fields.Ez=pxr.ez
                self.fields.Bx=pxr.bx
                self.fields.By=pxr.by
                self.fields.Bz=pxr.bz
                self.fields.Exp=pxr.ex
                self.fields.Eyp=pxr.ey
                self.fields.Ezp=pxr.ez
                self.fields.Bxp=pxr.bx
                self.fields.Byp=pxr.by
                self.fields.Bzp=pxr.bz
                self.fields.Jx=pxr.jx
                self.fields.Jy=pxr.jy
                self.fields.Jz=pxr.jz


                em3d_exchange_e(self.block)
                em3d_exchange_b(self.block)

                # If domain has been resized, do a new tile split and exchange particles
                if 1:#((isnewdom != 0)):
                   # Now exchanging particles
                    pxr.create_new_tile_split()

                    self.ntilex = pxr.ntilex
                    self.ntiley = pxr.ntiley
                    self.ntilez = pxr.ntilez

                  # Alias PXR tiles to WARP pgroups
                    for i,s in enumerate(self.listofallspecies):
                      s.pgroups = []
                      s.jslist = [0]
                      for iz in range(1,self.ntilez+1):
                        xygroup=[]
                        for iy in range(1,self.ntiley+1):
                          xgroup=[]
                          for ix in range(1,self.ntilex+1):
                            pg = ParticleGroup()
                            xgroup.append(pg)
                            pxr.point_to_tile(i+1, ix, iy, iz)
                            pg.npmax = 0
                            pxr.partnmax
                            pg.ns=1
                            pg.npid=top.npid
                            pg.gchange()
                            pg.sq = s.charge
                            pg.sm = s.mass
                            pg.sw = s.sw
                            pg.npmax = pxr.partnmax
                            pg.nps = pxr.partn
                            pg.ins[0] = 1
                            pg.sid[0]=0
                            pg.xp = pxr.partx
                            pg.yp = pxr.party
                            pg.zp = pxr.partz
                            pg.uxp = pxr.partux
                            pg.uyp = pxr.partuy
                            pg.uzp = pxr.partuz
                            pg.pid = fzeros([pg.npmax,top.npid])
                            pg.pid = pxr.pid
                            pg.gaminv = pxr.partgaminv
                            pg.ex = pxr.partex
                            pg.ey = pxr.partey
                            pg.ez = pxr.partez
                            pg.bx = pxr.partbx
                            pg.by = pxr.partby
                            pg.bz = pxr.partbz
                            pg.lebcancel_pusher=top.pgroup.lebcancel_pusher
                          xygroup.append(xgroup)
                        s.pgroups.append(xygroup)
                      pxr.set_are_tiles_reallocated(i+1, self.ntilex,self.ntiley,self.ntilez,zeros((self.ntilex,self.ntiley,self.ntilez),dtype=dtype('i8')))
#                pxr.particle_bcs_mpi_blocking()
                pxr.remap_particles(ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,
                            ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,
                            pxr.cell_x_min,pxr.cell_x_max,pxr.cell_y_min,pxr.cell_y_max,
                            pxr.cell_z_min,pxr.cell_z_max,
                            pxr.rank, pxr.nproc, pxr.nprocx, pxr.nprocy,pxr.nprocz,top.lcomm_cartesian)

            # Time statistics
            tend=MPI.Wtime()
            self.time_stat_loc_array[15] += (tend-tdeb)

    def load_balance_2d(self,imbalance):
        if (l_pxr):
            ## --- Compute time per part and per cell
            pxr.compute_time_per_part()
            pxr.compute_time_per_cell()

            ## --- Compute new split along each dimension
            pxr.compute_new_split_2d(pxr.global_time_per_part,pxr.global_time_per_cell,pxr.nx_global,pxr.nz_global,
                              pxr.new_cell_x_min,pxr.new_cell_x_max,
                              pxr.new_cell_z_min,pxr.new_cell_z_max,pxr.nprocx,pxr.nprocz)
            isnewsplit=sum(pxr.cell_x_min-pxr.new_cell_x_min)+sum(pxr.cell_x_max-pxr.new_cell_x_max)+ \
                     sum(pxr.cell_z_min-pxr.new_cell_z_min)+sum(pxr.cell_z_max-pxr.new_cell_z_max)
            if (isnewsplit==0):
                if(pxr.rank==0):
                  print("Optimal load balancing already achieved by current implementation")
            else:
                if(pxr.rank==0):
                  print("trying to load balance the simulation, imbalance=", imbalance)

                ## --- Compute limits for all procs
                ix1old=np.zeros(pxr.nproc,dtype="i8"); ix2old=np.zeros(pxr.nproc,dtype="i8")
                iy1old=np.zeros(pxr.nproc,dtype="i8"); iy2old=np.zeros(pxr.nproc,dtype="i8")
                iz1old=np.zeros(pxr.nproc,dtype="i8"); iz2old=np.zeros(pxr.nproc,dtype="i8")
                ix1new=np.zeros(pxr.nproc,dtype="i8"); ix2new=np.zeros(pxr.nproc,dtype="i8")
                iy1new=np.zeros(pxr.nproc,dtype="i8"); iy2new=np.zeros(pxr.nproc,dtype="i8")
                iz1new=np.zeros(pxr.nproc,dtype="i8"); iz2new=np.zeros(pxr.nproc,dtype="i8")

                pxr.get_1darray_proclimits(ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,
                                        pxr.cell_x_min,pxr.cell_y_min,pxr.cell_z_min,
                                        pxr.cell_x_max,pxr.cell_y_max,pxr.cell_z_max,
                                        pxr.nprocx, pxr.nprocy, pxr.nprocz, pxr.nproc,
                                        top.lcomm_cartesian)
                pxr.get_1darray_proclimits(ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,
                                        pxr.new_cell_x_min,pxr.new_cell_y_min,pxr.new_cell_z_min,
                                        pxr.new_cell_x_max,pxr.new_cell_y_max,pxr.new_cell_z_max,
                                        pxr.nprocx, pxr.nprocy, pxr.nprocz, pxr.nproc,top.lcomm_cartesian)
                ## --- Compute new sizes for grid arrays
                nx_new=pxr.new_cell_x_max[pxr.x_coords]-pxr.new_cell_x_min[pxr.x_coords]+1
                nz_new=pxr.new_cell_z_max[pxr.z_coords]-pxr.new_cell_z_min[pxr.z_coords]+1

                ## --- Remap field arrays
                # -- Ex
                ex_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(ex_new,nx_new,nz_new,
                                          pxr.ex,pxr.nx,pxr.nz,
                                          pxr.nxguards,pxr.nzguards,
                                          ix1old, ix2old, iz1old, iz2old,
                                          ix1new, ix2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.ex=ex_new
                # -- Ey
                ey_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(ey_new,nx_new,nz_new,
                                          pxr.ey,pxr.nx,pxr.nz,
                                          pxr.nxguards,pxr.nzguards,
                                          ix1old, ix2old, iz1old, iz2old,
                                          ix1new, ix2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.ey=ey_new
                # -- Ez
                ez_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(ez_new,nx_new,nz_new,
                                          pxr.ez,pxr.nx,pxr.nz,
                                          pxr.nxguards,pxr.nzguards,
                                          ix1old, ix2old, iz1old, iz2old,
                                          ix1new, ix2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.ez=ez_new
                # -- Bx
                bx_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(bx_new,nx_new,nz_new,
                                          pxr.bx,pxr.nx,pxr.nz,
                                          pxr.nxguards,pxr.nzguards,
                                          ix1old, ix2old, iz1old, iz2old,
                                          ix1new, ix2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.bx=bx_new
                # -- By
                by_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(by_new,nx_new,nz_new,
                                          pxr.by,pxr.nx,pxr.nz,
                                          pxr.nxguards,pxr.nzguards,
                                          ix1old, ix2old, iz1old, iz2old,
                                          ix1new, ix2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.by=by_new
                # -- Bz
                bz_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(bz_new,nx_new,nz_new,
                                          pxr.bz,pxr.nx,pxr.nz,
                                          pxr.nxguards,pxr.nzguards,
                                          ix1old, ix2old, iz1old, iz2old,
                                          ix1new, ix2new, iz1new, iz2new,
                                          pxr.rank, pxr.nproc)
                pxr.bz=bz_new
                ## -- Reallocate current arrays
                # Currents are recomputed each iteration so no need to exchange them
                jx_new=zeros((nx_new+2*pxr.nxjguards+1,1,nz_new+2*pxr.nzjguards+1),order='F')
                jy_new=zeros((nx_new+2*pxr.nxjguards+1,1,nz_new+2*pxr.nzjguards+1),order='F')
                jz_new=zeros((nx_new+2*pxr.nxjguards+1,1,nz_new+2*pxr.nzjguards+1),order='F')
                pxr.jx=jx_new
                pxr.jy=jy_new
                pxr.jz=jz_new

                # Update pxr new array dimensions
                pxr.nx=nx_new
                pxr.nz=nz_new
                pxr.nx_grid=pxr.nx+1
                pxr.nz_grid=pxr.nz+1

                # Test if domain has been resized - used for particle remaping
                isnewdom=pxr.cell_x_min[pxr.x_coords]-pxr.new_cell_x_min[pxr.x_coords]+pxr.cell_x_max[pxr.x_coords]-pxr.new_cell_x_max[pxr.x_coords]+ \
                pxr.cell_z_min[pxr.z_coords]-pxr.new_cell_z_min[pxr.z_coords]+pxr.cell_z_max[pxr.z_coords]-pxr.new_cell_z_max[pxr.z_coords]

                # Update new subdomain index arrays
                pxr.cell_x_min=pxr.new_cell_x_min
                pxr.cell_x_max=pxr.new_cell_x_max
                pxr.cell_z_min=pxr.new_cell_z_min
                pxr.cell_z_max=pxr.new_cell_z_max
                pxr.nx_global_grid_min = pxr.cell_x_min[pxr.x_coords]
                pxr.nx_global_grid_max = pxr.cell_x_max[pxr.x_coords]+1
                pxr.nz_global_grid_min = pxr.cell_z_min[pxr.z_coords]
                pxr.nz_global_grid_max = pxr.cell_z_max[pxr.z_coords]+1


                # Update global simulation axis
                pxr.compute_simulation_axis()

                # Set new min and max for local domain
                pxr.x_min_local = pxr.x_grid_mins[pxr.x_coords]
                pxr.x_max_local = pxr.x_grid_maxs[pxr.x_coords]
                pxr.z_min_local = pxr.z_grid_mins[pxr.z_coords]
                pxr.z_max_local = pxr.z_grid_maxs[pxr.z_coords]
                pxr.x_grid_min_local=pxr.x_min_local
                pxr.x_grid_max_local=pxr.x_max_local
                pxr.z_grid_min_local=pxr.z_min_local
                pxr.z_grid_max_local=pxr.z_max_local

                ##--- Alias WARP grid arrays on pxr new arrays
                self.nxlocal=pxr.nx
                self.nzlocal=pxr.nz
                self.xmminlocal = pxr.x_min_local
                self.zmminlocal = pxr.z_min_local
                self.fields.xmin = pxr.x_min_local
                self.fields.xmax = pxr.x_max_local
                self.fields.zmin = pxr.z_min_local
                self.fields.zmax = pxr.z_max_local

                # Udpate domain decomposition in WARP
                top.fsdecomp.nx=pxr.cell_x_max-pxr.cell_x_min+1
                top.fsdecomp.nz=pxr.cell_z_max-pxr.cell_z_min+1
                top.fsdecomp.ix=pxr.cell_x_min
                top.fsdecomp.iz=pxr.cell_z_min
                top.fsdecomp.xmin=pxr.cell_x_min*pxr.dx
                top.fsdecomp.xmax=(pxr.cell_x_max+1)*pxr.dx
                top.fsdecomp.zmin=pxr.cell_z_min*pxr.dz
                top.fsdecomp.zmax=(pxr.cell_z_max+1)*pxr.dz
                top.ppdecomp.nx=pxr.cell_x_max-pxr.cell_x_min+1
                top.ppdecomp.nz=pxr.cell_z_max-pxr.cell_z_min+1
                top.ppdecomp.ix=pxr.cell_x_min
                top.ppdecomp.iz=pxr.cell_z_min
                top.ppdecomp.xmin=pxr.cell_x_min*pxr.dx
                top.ppdecomp.xmax=(pxr.cell_x_max+1)*pxr.dx
                top.ppdecomp.zmin=pxr.cell_z_min*pxr.dz
                top.ppdecomp.zmax=(pxr.cell_z_max+1)*pxr.dz

                # Reallocate warp arrays
                self.allocatefieldarrays()
                if (self.spectral==1):
                  self.allocatefieldarraysFFT()
                # Alias newly allocated arrays on WARP structure
                self.fields.Ex=pxr.ex
                self.fields.Ey=pxr.ey
                self.fields.Ez=pxr.ez
                self.fields.Bx=pxr.bx
                self.fields.By=pxr.by
                self.fields.Bz=pxr.bz
                self.fields.Exp=pxr.ex
                self.fields.Eyp=pxr.ey
                self.fields.Ezp=pxr.ez
                self.fields.Bxp=pxr.bx
                self.fields.Byp=pxr.by
                self.fields.Bzp=pxr.bz
                self.fields.Jx=pxr.jx
                self.fields.Jy=pxr.jy
                self.fields.Jz=pxr.jz


                em3d_exchange_e(self.block)
                em3d_exchange_b(self.block)
                #If domain has been resized, do a new tile split and exchange particles
                if 1:#((isnewdom != 0)):
                # Now exchanging particles
                    pxr.create_new_tile_split()
                    pxr.remap_particles_2d(ix1old,ix2old,iz1old,iz2old,
                                            ix1new,ix2new,iz1new,iz2new,
                                            pxr.cell_x_min,pxr.cell_x_max,
                                            pxr.cell_z_min,pxr.cell_z_max,
                                            pxr.rank, pxr.nproc, pxr.nprocx,pxr.nprocz,top.lcomm_cartesian)
                    self.ntilex = pxr.ntilex
                    self.ntilez = pxr.ntilez

          # Alias PXR tiles to WARP pgroups
                    for i,s in enumerate(self.listofallspecies):
                        s.pgroups = []
                        s.jslist = [0]
                        for iz in range(1,self.ntilez+1):
                          xygroup=[]
                          for iy in range(1,self.ntiley+1):
                            xgroup=[]
                            for ix in range(1,self.ntilex+1):
                              pg = ParticleGroup()
                              xgroup.append(pg)
                              pxr.point_to_tile(i+1, ix, iy, iz)
                              pg.npmax = 0
                              pxr.partnmax
                              pg.ns=1
                              pg.npid=top.npid
                              pg.gchange()
                              pg.sq = s.charge
                              pg.sm = s.mass
                              pg.sw = s.sw
                              pg.npmax = pxr.partnmax
                              pg.nps = pxr.partn
                              pg.ins[0] = 1
                              pg.sid[0]=0
                              pg.xp = pxr.partx
                              pg.yp = pxr.party
                              pg.zp = pxr.partz
                              pg.uxp = pxr.partux
                              pg.uyp = pxr.partuy
                              pg.uzp = pxr.partuz
                              pg.pid = fzeros([pg.npmax,top.npid])
                              pg.pid = pxr.pid
                              pg.gaminv = pxr.partgaminv
                              pg.ex = pxr.partex
                              pg.ey = pxr.partey
                              pg.ez = pxr.partez
                              pg.bx = pxr.partbx
                              pg.by = pxr.partby
                              pg.bz = pxr.partbz
                              pg.lebcancel_pusher=top.pgroup.lebcancel_pusher
                            xygroup.append(xgroup)
                          s.pgroups.append(xygroup)
                        pxr.set_are_tiles_reallocated(i+1, self.ntilex,self.ntiley,self.ntilez,zeros((self.ntilex,self.ntiley,self.ntilez),dtype=dtype('i8')))


    def output_pxr(self,iter):
      pxr.py_mpi_output_grid_quantity('ez',pxr.ez,pxr.nx,pxr.ny,pxr.nz,pxr.nxguards,pxr.nyguards,pxr.nzguards,iter)


    def fetcheb(self,js,pg=None):
        if self.l_verbose:print me,'enter fetcheb'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        w3d.ipminfsapi=pg.ins[js]
        w3d.npfsapi=pg.nps[js]
        pg.ex[il:iu]=0.
        pg.ey[il:iu]=0.
        pg.ez[il:iu]=0.
        pg.bx[il:iu]=0.
        pg.by[il:iu]=0.
        pg.bz[il:iu]=0.
        self.fetche()
        self.fetchb()

    def push_velocity_full(self,js,pg=None):
        if self.l_verbose:print me,'enter push_ions_velocity_full'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        if pg.lebcancel_pusher:
          ebcancelpush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                            pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                            pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                            pg.sq[js],pg.sm[js],top.dt,0)
        else:
          # --- push velocity from electric field (half step)
          epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                     pg.sq[js],pg.sm[js],0.5*top.dt)
          # --- update gamma
          self.set_gamma(js,pg)
          # --- push velocity from magnetic field
          bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                      pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                      pg.sq[js],pg.sm[js],top.dt, top.ibpush)
          # --- push velocity from electric field (half step)
          epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                     pg.sq[js],pg.sm[js],0.5*top.dt)
          # --- update gamma
          self.set_gamma(js,pg)

        if self.l_verbose:print me,'exit push_ions_velocity_first_half'

    def push_velocity_first_half(self,js,pg=None):
        if self.l_verbose:print me,'enter push_ions_velocity_first_half'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        if pg.lebcancel_pusher:
          ebcancelpush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                            pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                            pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                            pg.sq[js],pg.sm[js],top.dt,1)
        else:
          # --- push velocity from electric field (half step)
          epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                     pg.sq[js],pg.sm[js],0.5*top.dt)
          # --- update gamma
          self.set_gamma(js,pg)
          # --- push velocity from magnetic field
          bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                      pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                      pg.sq[js],pg.sm[js],0.5*top.dt, top.ibpush)

        if self.l_verbose:print me,'exit push_ions_velocity_first_half'

    def push_velocity_second_half(self,js,pg=None):
        if self.l_verbose:print me,'enter push_ions_velocity_second_half'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        if pg.lebcancel_pusher:
          ebcancelpush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                            pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                            pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                            pg.sq[js],pg.sm[js],top.dt,2)
        else:
          # --- push velocity from magnetic field
          bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                      pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                      pg.sq[js],pg.sm[js],0.5*top.dt, top.ibpush)
          # --- push velocity from electric field (half step)
          epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                     pg.sq[js],pg.sm[js],0.5*top.dt)
        # --- update gamma
        self.set_gamma(js,pg)

        if self.l_verbose:print me,'exit push_ions_velocity_second_half'

    def set_gamma(self,js,pg=None):
        if self.l_verbose:print me,'enter set_gamma'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        # --- update gamma
        gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                 top.gamadv,top.lrelativ)

        if self.l_verbose:print me,'exit push_ions_velocity_second_half'

    def push_positions(self,js,pg=None):
        if self.l_verbose:print me,'enter push_ions_positions'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                       pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                       pg.gaminv[il:iu],top.dt)

        if self.l_verbose:print me,'exit push_ions_positions'

    def loadsource(self,lzero=None,lfinalize_rho=None,pgroups=None,**kw):
        '''
        Current and charge deposition, uses particles from top directly.

        Inputs:
           - lzero
           - lfinalize_rho
           - pgroups
        '''

        # --- Note that the grid location is advanced even if no field solve
        # --- is being done.
        self.advancezgrid()
        # --- If ldosolve is false, then skip the gather of rho, unless
        # --- lzero is also false, in which case the solver is assumed to
        # --- be gathering the source (for example during an EGUN iteration).
        if not self.ldosolve and lzero: return
        if lzero is None: lzero = w3d.lzerorhofsapi
        if lfinalize_rho is None: lfinalize_rho = w3d.lfinalizerhofsapi

        self.setparticledomains()
        self.allocatedataarrays()
        if lzero: self.zerosourcep()

        if pgroups is None: pgroups = [top.pgroup]

        if  l_pxr:
            # --- PICSAR current deposition
            # --- js = 0
             f=self.fields

             for pgroup in pgroups:
                if w3d.js1fsapi >= 0: js1 = w3d.js1fsapi
                else:                 js1 = 0
                if w3d.js2fsapi >= 0: js2 = w3d.js2fsapi+1
                else:                 js2 = pgroup.ns

                jslist = kw.get('jslist',None)
                if jslist is None: jslist = range(js1,js2)

                for js in jslist:
                    n = pgroup.nps[js]
                    if n == 0: continue
                    if pgroup.ldts[js]:
                        indts = top.ndtstorho[pgroup.ndts[js]-1]
                        iselfb = pgroup.iselfb[js]
                        self.setsourcepforparticles(0,indts,iselfb)

                        if self.debug:
                            i1 = pgroup.ins[js]-1
                            i2 = pgroup.ins[js]+pgroup.nps[js]-1
                            if self.nxlocal > 0:
                                x = pgroup.xp[i1:i2]
                                if self.l4symtry: x = abs(x)
                                if self.solvergeom == w3d.RZgeom:
                                    y = pgroup.yp[i1:i2]
                                    x = sqrt(x**2 + y**2)
                                assert x.min() >= self.xmminp,\
                                       "Particles in species %d have x below the grid when depositing the source, min x = %e"%(js,x.min())
                                assert x.max() < self.xmmaxp,\
                                       "Particles in species %d have x above the grid when depositing the source, max x = %e"%(js,x.max())
                            if self.nylocal > 0:
                                y = pgroup.yp[i1:i2]
                                if self.l4symtry or self.l2symtry: y = abs(y)
                                assert y.min() >= self.ymminp,\
                                       "Particles in species %d have y below the grid when depositing the source, min y = %e"%(js,y.min())
                                assert y.max() < self.ymmaxp,\
                                       "Particles in species %d have y above the grid when depositing the source, max y = %e"%(js,y.max())
                            if self.nzlocal > 0:
                                z = pgroup.zp[i1:i2]
                                assert z.min() >= self.zmminp+self.getzgridndts()[indts],\
                                       "Particles in species %d have z below the grid when depositing the source, min z = %e"%(js,z.min())
                                assert z.max() < self.zmmaxp+self.getzgridndts()[indts],\
                                       "Particles in species %d have z above the grid when depositing the source, max z = %e"%(js,z.max())

             # ___________________________________
             # Depose currents in PXR

             if (self.l_debug): print("Call pxr.pxrdepose_currents_on_grid_jxjyjz()")

             t0 = MPI.Wtime()

             pxr.jx = self.fields.Jx
             pxr.jy = self.fields.Jy
             pxr.jz = self.fields.Jz

             if pxr.c_dim == 2:
               pxr.pxrdepose_currents_on_grid_jxjyjz_2d()

               #pxr.pxrdepose_currents_on_grid_jxjyjz_sub_openmp(f.Jx,f.Jy,f.Jz,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,
               #pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt)

             elif pxr.c_dim ==3:

               pxr.pxrdepose_currents_on_grid_jxjyjz()

             # Time statistics
             t1 = MPI.Wtime()
             self.time_stat_loc_array[2] += (t1-t0)

             # ___________________________________
             # Depose charge density in PXR if required

             if self.l_getrho : # Depose Rho in PXR

               if (self.l_debug): print("Call pxr.pxrdepose_rho_on_grid_sub_openmp()")

               t0 = MPI.Wtime()

               if pxr.c_dim == 2:

                 pxr.rho = self.fields.Rho
                 pxr.pxrdepose_rho_on_grid()
                 #pxr.pxrdepose_rho_on_grid_sub_openmp_2d(f.Rho,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt,0)

               elif pxr.c_dim ==3:

                 pxr.rho = self.fields.Rho
                 pxr.pxrdepose_rho_on_grid()

               # Time statistics
               t1 = MPI.Wtime()
               self.time_stat_loc_array[12] += (t1-t0)

             #pxr.pxrdepose_rho_on_grid_sub_openmp_3d(f.Rho,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt,0)
             if self.current_cor: # Depose Rhoold_local in PXR
                 t0 = MPI.Wtime()
                 pxr.pxrdepose_rho_on_grid_sub_openmp_3d(f.Rhoold_local,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt,1)
                 t1 = MPI.Wtime()
                 self.time_stat_loc_array[12] += (t1-t0)

        else:

            for pgroup in pgroups:

                if w3d.js1fsapi >= 0: js1 = w3d.js1fsapi
                else:                 js1 = 0
                if w3d.js2fsapi >= 0: js2 = w3d.js2fsapi+1
                else:                 js2 = pgroup.ns

                jslist = kw.get('jslist',None)
                if jslist is None: jslist = range(js1,js2)

                for js in jslist:
                    n = pgroup.nps[js]
                    if n == 0: continue
                    if pgroup.ldts[js]:
                        indts = top.ndtstorho[pgroup.ndts[js]-1]
                        iselfb = pgroup.iselfb[js]
                        self.setsourcepforparticles(0,indts,iselfb)

                        if self.debug:
                            i1 = pgroup.ins[js]-1
                            i2 = pgroup.ins[js]+pgroup.nps[js]-1
                            if self.nxlocal > 0:
                                x = pgroup.xp[i1:i2]
                                if self.l4symtry: x = abs(x)
                                if self.solvergeom == w3d.RZgeom:
                                    y = pgroup.yp[i1:i2]
                                    x = sqrt(x**2 + y**2)
                                assert x.min() >= self.xmminp,\
                                       "Particles in species %d have x below the grid when depositing the source, min x = %e"%(js,x.min())
                                assert x.max() < self.xmmaxp,\
                                       "Particles in species %d have x above the grid when depositing the source, max x = %e"%(js,x.max())
                            if self.nylocal > 0:
                                y = pgroup.yp[i1:i2]
                                if self.l4symtry or self.l2symtry: y = abs(y)
                                assert y.min() >= self.ymminp,\
                                       "Particles in species %d have y below the grid when depositing the source, min y = %e"%(js,y.min())
                                assert y.max() < self.ymmaxp,\
                                       "Particles in species %d have y above the grid when depositing the source, max y = %e"%(js,y.max())
                            if self.nzlocal > 0:
                                z = pgroup.zp[i1:i2]
                                assert z.min() >= self.zmminp+self.getzgridndts()[indts],\
                                       "Particles in species %d have z below the grid when depositing the source, min z = %e"%(js,z.min())
                                assert z.max() < self.zmmaxp+self.getzgridndts()[indts],\
                                       "Particles in species %d have z above the grid when depositing the source, max z = %e"%(js,z.max())

                        self.setsourcep(js,pgroup,self.getzgridndts()[indts])

        # --- Only finalize the source if lzero is true, which means the this
        # --- call to loadsource should be a complete operation.
        self.sourcepfinalized = False
        if lzero and lfinalize_rho: self.finalizesourcep()


    def apply_bndconditions(self,js,pg=None):
        if self.l_verbose:print me,'enter apply_ions_bndconditions'
        # --- apply boundary conditions
        if pg is None:
            pg = top.pgroup
        if pg.nps[js]==0:return
        self.apply_bnd_conditions(js,pg)
        if self.l_verbose:print me,'exit apply_ions_bndconditions'

    def apply_bnd_conditions(self,js,pg=None):
        if self.l_verbose:print me,'enter apply_bnd_conditions'
        if pg is None:
            pg = top.pgroup
        if pg.nps[js]==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        #stckxy3d(pg.nps[js],pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
        #              pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
        #              pg.zp[il:iu],w3d.zmminlocal,w3d.dz,
        #              pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
        #              top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,top.pboundxy,true)
        stckxy3d(pg,js,top.zbeam,true)
        partbndwithdata(pg.nps[js],pg.xp[il:iu],pg.uxp[il:iu],pg.gaminv[il:iu],
                        w3d.xmmaxlocal,w3d.xmminlocal,w3d.dx,0.,
                        top.pboundxy,top.pboundxy)
        partbndwithdata(pg.nps[js],pg.yp[il:iu],pg.uyp[il:iu],pg.gaminv[il:iu],
                        w3d.ymmaxlocal,w3d.ymminlocal,w3d.dy,0.,
                        top.pboundxy,top.pboundxy)
        if js==0 or js==w3d.nzp-1:
          if js==0:top.pboundnz=-1
          if js==w3d.nzp-1:top.pbound0=-1
          partbndwithdata(pg.nps[js],pg.zp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                          w3d.zmmaxlocal,w3d.zmminlocal,w3d.dz,top.zgrid,
                          top.pbound0,top.pboundnz)
          if js==0:top.pboundnz=0
          if js==w3d.nzp-1:top.pbound0=0
        if self.scraper is not None:self.scraper.scrape(js)
        processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
        if self.l_verbose:print me,'enter apply_bnd_conditions'

    def get_total_particle_number(self,**kw):
        """
        Get the total number of particles from all species

        output:
        - total number of particles
        """

        nbptot = zeros(1,dtype=numpy.int64)

        pxr.get_tot_number_of_particles(nbptot)

        return nbptot[0]

    def get_kinetic_energy(self,sp,**kw):
        """
        Get the total kinetic energy of the species sp using PICSAR fortran subroutines

        input:
        - sp: species number
        """
        total_kinetic_energy = zeros(1)
        if self.l_verbose:print me,'compute kinetic energy on species',sp
        pxr.get_kinetic_energy(sp,total_kinetic_energy)
        #print total_kinetic_energy,sp
        return total_kinetic_energy[0]

    def get_field_energy(self,field,**kw):
        """
        Get the total field energy for the given component.
        The field energy is calculated in parallel with a picsar fortran subroutine.

        input:
        - field: field component
        """
        field_energy = zeros(1)

        if pxr.c_dim==2:

          if field=='ex':
            pxr.get_field_energy_2d(self.fields.Ex,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='ey':
            pxr.get_field_energy_2d(self.fields.Ey,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='ez':
            pxr.get_field_energy_2d(self.fields.Ez,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='bx':
            pxr.get_field_energy_2d(self.fields.Bx,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='by':
            pxr.get_field_energy_2d(self.fields.By,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='bz':
            pxr.get_field_energy_2d(self.fields.Bz,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          return field_energy[0]

        else:

          if field=='ex':
            pxr.get_field_energy(self.fields.Ex,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='ey':
            pxr.get_field_energy(self.fields.Ey,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='ez':
            pxr.get_field_energy(self.fields.Ez,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='bx':
            pxr.get_field_energy(self.fields.Bx,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='by':
            pxr.get_field_energy(self.fields.By,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='bz':
            pxr.get_field_energy(self.fields.Bz,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          return field_energy[0]

    def get_normL2_divEeps0_rho(self):
        """
        Compute the L2 norm of divE*eps0 - rho
        Computation of rho has to be activated

        """
        div = zeros(1)

        pxr.calc_field_div(pxr.dive,pxr.ex, pxr.ey, pxr.ez, pxr.nx,pxr.ny,pxr.nz,pxr.nxguards,pxr.nyguards,pxr.nzguards,pxr.dx,pxr.dy,pxr.dz)

        pxr.get_norm_diverho(pxr.dive,pxr.rho,pxr.nx,pxr.ny,pxr.nz,pxr.nxguards,pxr.nyguards,pxr.nzguards,div)

        return div[0]

    def display_picsar_time_statistics(self):
        """
        Display the Picsar time statistics
        """
        pxr.time_statistics()

    def display_time_statistics(self,):
        """
        Display the time statistics
        """
        self.time_stat_ave_array = zeros([20])
        self.time_stat_min_array = zeros([20])
        self.time_stat_max_array = zeros([20])
        nproc = pxr.nprocx*pxr.nprocy*pxr.nprocz

        MPI.COMM_WORLD.Reduce([self.time_stat_loc_array,MPI.DOUBLE], [self.time_stat_ave_array,MPI.DOUBLE], op=MPI.SUM, root=0)
        MPI.COMM_WORLD.Reduce([self.time_stat_loc_array,MPI.DOUBLE], [self.time_stat_min_array,MPI.DOUBLE], op=MPI.MIN, root=0)
        MPI.COMM_WORLD.Reduce([self.time_stat_loc_array,MPI.DOUBLE], [self.time_stat_max_array,MPI.DOUBLE], op=MPI.MAX, root=0)

        self.time_stat_ave_array[:] /= nproc

        if me==0:

          print ' _______________________________________________________________________________'
          print
          print '  Time statisctics'
          print ' _______________________________________________________________________________'

          print ' Parts                              {:^8} {:^8} {:^8} {:^8}'.format('min', 'ave', 'max', '%')
          print ' -------------------------------------------------------------------------------'
          print ' Particle pusher + field gathering: {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[0],self.time_stat_ave_array[0],self.time_stat_max_array[0],self.time_stat_max_array[0]/self.total_kernel_time*100)
          print ' Particle boundary conditions:      {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[1],self.time_stat_ave_array[1],self.time_stat_max_array[1],self.time_stat_max_array[1]/self.total_kernel_time*100)
          print ' Current deposition:                {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[2],self.time_stat_ave_array[2],self.time_stat_max_array[2],self.time_stat_max_array[2]/self.total_kernel_time*100)
          print ' Current bound. cond.:              {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[3],self.time_stat_ave_array[3],self.time_stat_max_array[3],self.time_stat_max_array[3]/self.total_kernel_time*100)
          print ' Magnetic field solver:             {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[5],self.time_stat_ave_array[5],self.time_stat_max_array[5],self.time_stat_max_array[5]/self.total_kernel_time*100)
          print ' Magnetic field bound. cond.:       {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[6],self.time_stat_ave_array[6],self.time_stat_max_array[6],self.time_stat_max_array[6]/self.total_kernel_time*100)
          print ' Electric field solver:             {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[7],self.time_stat_ave_array[7],self.time_stat_max_array[7],self.time_stat_max_array[7]/self.total_kernel_time*100)
          print ' Electric field bound. cond.:       {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[8],self.time_stat_ave_array[8],self.time_stat_max_array[8],self.time_stat_max_array[8]/self.total_kernel_time*100)
          print ' Particle sorting:                  {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[10],self.time_stat_ave_array[10],self.time_stat_max_array[10],self.time_stat_max_array[10]/self.total_kernel_time*100)
          print ' Charge deposition:                 {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[12],self.time_stat_ave_array[12],self.time_stat_max_array[12],self.time_stat_max_array[12]/self.total_kernel_time*100)
          print ' Charge bound. cond.:               {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[13],self.time_stat_ave_array[13],self.time_stat_max_array[13],self.time_stat_max_array[13]/self.total_kernel_time*100)
          print ' Load balancing:                    {:8.3f} {:8.3f} {:8.3f} {:8.3}'.format(self.time_stat_min_array[15],self.time_stat_ave_array[15],self.time_stat_max_array[15],self.time_stat_max_array[15]/self.total_kernel_time*100)
          print


    def allocatefieldarraysFFT(self):
        def fc(x,norder):
            fact1 = 1
            fact2 = 1
            result = 0
            for i in range(abs(norder)/2):
              fact1 *= max(i,1)
              fact2 *= max(2*i,1)*max(2*i-1,1)
              result += x**(2*i+1)*fact2/float(2**(2*i)*fact1**2*(2*i+1))
            return result


        f=self.fields
        b=self.block
        s=self
        f.spectral = (self.spectral > 0)
        bc_periodic = [self.bounds[0]==periodic,
                       self.bounds[2]==periodic,
                       self.bounds[4]==periodic]
        if self.current_cor:
            f.nxdrho = f.nx
            f.nydrho = f.ny
            f.nzdrho = f.nz
            f.nxdrhoguard = f.nxguard
            f.nydrhoguard = f.nyguard
            f.nzdrhoguard = f.nzguard
            f.gchange()

        if self.spectral:
            kwGPSTD = {'l_staggered':s.l_spectral_staggered,\
                     'spectral':s.spectral,\
                     'norderx':s.norderx,\
                     'nordery':s.nordery,\
                     'norderz':s.norderz,\
                     'nxguard':s.nxguard,\
                     'nyguard':s.nyguard,\
                     'nzguard':s.nzguard,\
                     'dt':top.dt,\
                     'dx':w3d.dx,\
                     'dy':w3d.dy,\
                     'dz':w3d.dz,\
                     'ntsub':s.ntsub,\
                     'l_pushf':s.l_pushf,\
                     'l_pushg':s.l_pushg,\
                     'l_getrho':s.l_getrho,\
                     'clight':clight}

            if s.ntsub is np.inf:
                if not self.l_getrho:
                    self.l_getrho = True
                    f.nxr = f.nx
                    f.nyr = f.ny
                    f.nzr = f.nz
                    f.gchange()
                if(self.full_pxr == False):
                  self.GPSTDMaxwell = gpstd.PSATD_Maxwell(yf=self.fields,
                                                  eps0=eps0,
                                                  bc_periodic=bc_periodic,
                                                  **kwGPSTD)
            else:
                if self.l_pushf and not self.l_getrho:
                    self.l_getrho = True
                    f.nxr = f.nx
                    f.nyr = f.ny
                    f.nzr = f.nz
                    f.gchange()
                if(self.full_pxr==False):
                  self.GPSTDMaxwell = gpstd.GPSTD_Maxwell(yf=self.fields,
                                                  eps0=eps0,
                                                  bc_periodic=bc_periodic,
                                                  **kwGPSTD)
                if(self.full_pxr==False):
                  self.FSpace = self.GPSTDMaxwell
        else:
            kwFS = {'l_staggered':s.l_spectral_staggered,\
                     'spectral':s.spectral,\
                     'norderx':s.norderx,\
                     'nordery':s.nordery,\
                     'norderz':s.norderz,\
                     'nxguard':s.nxguard,\
                     'nyguard':s.nyguard,\
                     'nzguard':s.nzguard,\
                     'dt':top.dt,\
                     'dx':w3d.dx,\
                     'dy':w3d.dy,\
                     'nx':max([1,self.fields.nx]),\
                     'ny':max([1,self.fields.ny]),\
                     'nz':max([1,self.fields.nz]),\
                     'dz':w3d.dz}
            self.FSpace = Fourier_Space(bc_periodic=bc_periodic,**kwFS)

        # --- computes Brendan's Jz,Jx multipliers
        if self.Jmult and self.GPSTDMaxwell.nz>1:
                k = self.GPSTDMaxwell.k
                if self.GPSTDMaxwell.nx>1:kxvzdto2 = 0.5*self.GPSTDMaxwell.kx*clight*top.dt
                if self.GPSTDMaxwell.ny>1:kyvzdto2 = 0.5*self.GPSTDMaxwell.ky*clight*top.dt
                kzvzdto2 = 0.5*self.GPSTDMaxwell.kz*clight*top.dt
                sinkzvzdto2 = sin(kzvzdto2)
                coskzvzdto2 = cos(kzvzdto2)
                kdto2 = 0.5*k*clight*top.dt
                sinkdto2 = sin(kdto2)
                coskdto2 = cos(kdto2)
                numer = clight*top.dt*k*self.kz*(self.sinkzvzdto2**2-self.sinkdto2**2)
                denom = 2*sinkdto2*sinkzvzdto2 \
                      * (self.GPSTDMaxwell.kz*sinkzvzdto2*coskdto2-k*coskzvzdto2*sinkdto2)
                denomno0 = where(denom==0.,0.0001,self.denom)

                raise Exception("What is the 3-D version of Brendan's correction?")

                ktest=where((pi/2-kxvzdto2**2/(2*pi))>0,(pi/2-kxvzdto2**2/(2*pi)),0)

                Jmultiplier = where(abs(self.kzvzdto2)<ktest,numer/denomno0,0)

                self.Jmultiplier[0,:]=self.Jmultiplier[1,:]
                self.Jmultiplier[:,0]=self.Jmultiplier[:,1]

        # --- set Ex,By multipliers (ebcor=0,1,2)
        if self.l_correct_num_Cherenkov and self.spectral:
              emK = self.FSpace
#              k = emK.k
              k = sqrt(emK.kx_unmod*emK.kx_unmod+emK.ky_unmod*emK.ky_unmod+emK.kz_unmod*emK.kz_unmod)
              if top.boost_gamma==1.:
                  raise Exception('Error: l_correct_num_Cherenkov=True with top.boost_gamma=1.')

              b0 = sqrt(1.-1./top.boost_gamma**2)
              self.b0=b0
              self.ebcor = 2

              if 0:

              # --- old coefs
                  # --- set Ex,By multipliers (ebcor=0,1,2)
                  if self.ebcor==2:
                      self.kzvzdto2 = where(emK.kz_unmod==0,0.0001,0.5*emK.kz_unmod*b0*clight*top.dt)
                      self.sinkzvzdto2 = sin(self.kzvzdto2)
                      self.coskzvzdto2 = cos(self.kzvzdto2)
                      self.Exmultiplier = self.kzvzdto2*self.coskzvzdto2/self.sinkzvzdto2

                      self.kdto2 = where(k==0,0.0001,0.5*k*clight*top.dt)
                      self.sinkdto2 = sin(self.kdto2)
                      self.coskdto2 = cos(self.kdto2)
                      self.Bymultiplier = self.kdto2*self.coskdto2/self.sinkdto2

                  if self.ebcor==1:
                      self.kzvzdto2 = where(emK.kz_unmod==0,0.0001,0.5*emK.kz_unmod*b0*clight*top.dt)
                      self.sinkzvzdto2 = sin(self.kzvzdto2)
                      self.coskzvzdto2 = cos(self.kzvzdto2)
                      self.kdto2 = where(k==0,0.0001,0.5*k*clight*top.dt)
                      self.sinkdto2 = sin(self.kdto2)
                      self.coskdto2 = cos(self.kdto2)
                      self.Exmultiplier = self.kdto2*self.sinkdto2**2*self.sinkzvzdto2*self.coskzvzdto2/ \
                        (self.kzvzdto2*(self.kdto2*self.sinkdto2**2+ \
                        (self.sinkdto2*self.coskdto2-self.kdto2)*self.sinkzvzdto2**2))

              else:
              # --- new cooefs
                  if self.ebcor==2:
                      # --- set Ex multiplier
                      self.kzvzdto2 = where(emK.kz_unmod==0,0.0001,0.5*emK.kz_unmod*b0*clight*top.dt)
                      self.sinkzvzdto2 = sin(self.kzvzdto2)
                      self.coskzvzdto2 = cos(self.kzvzdto2)
                      self.Exmultiplier = self.kzvzdto2*self.coskzvzdto2/self.sinkzvzdto2
                      # --- set By multiplier
                      if self.norderx is inf:
                          self.kdto2 = where(k==0,0.0001,0.5*k*clight*top.dt)
                      else:
                          self.kdto2 = sqrt((fc(sin(emK.kx_unmod*0.5*self.dx),self.norderx)/(0.5*self.dx))**2+ \
                              (fc(sin(emK.kz_unmod*0.5*self.dz),self.norderz)/(0.5*self.dz))**2)
                          self.kdto2 = where(self.kdto2==0,0.0001,0.5*self.kdto2*clight*top.dt)
                      if 0:#self.solver==PSATD:
                          self.Bymultiplier = self.kdto2/tan(self.kdto2)
                      else:
                          self.thetadto2=self.ntsub*arcsin(self.kdto2/self.ntsub)
                          self.Bymultiplier = self.kdto2/(tan(self.thetadto2)*cos(self.thetadto2/self.ntsub))

                  if self.ebcor==1:
                      self.kzvzdto2 = where(emK.kz_unmod==0,0.0001,0.5*emK.kz_unmod*b0*clight*top.dt)
                      self.sinkzvzdto2 = sin(self.kzvzdto2)
                      self.coskzvzdto2 = cos(self.kzvzdto2)
                      if self.norderx is None:
                          self.kdto2 = where(k==0,0.0001,0.5*k*clight*top.dt)
                      else:
                          self.kdto2 = sqrt((fc(sin(emK.kx_unmod*0.5*self.dx),self.norderx)/(0.5*self.dx))**2+ \
                              (fc(sin(emK.kz_unmod*0.5*self.dz),self.norderz)/(0.5*self.dz))**2)
                          self.kdto2 = where(self.kdto2==0,0.0001,0.5*self.kdto2*clight*top.dt)
                          self.kzvzdto2 = fc(sin(emK.kz_unmod*0.5*self.dz),self.norderz)/(0.5*self.dz)
                          self.kzvzdto2 = where(self.kzvzdto2==0,0.0001,0.5*self.kzvzdto2*b0*clight*top.dt)
                      if 0:#:self.solver==PSATD:
                          self.sinkdto2 = sin(self.kdto2)
                          self.coskdto2 = cos(self.kdto2)
                          self.Exmultiplier = self.kdto2*self.sinkdto2**2*self.sinkzvzdto2*self.coskzvzdto2/ \
                           (self.kzvzdto2*(self.kdto2*self.sinkdto2**2+ \
                           (self.sinkdto2*self.coskdto2-self.kdto2)*self.sinkzvzdto2**2))
                      else:
                          self.thetadto2=self.ntsub*arcsin(self.kdto2/self.ntsub)
                          self.Exmultiplier = self.ntsub*self.sinkzvzdto2*self.coskzvzdto2*sin(self.thetadto2)**2/ \
                           (self.kzvzdto2*(self.ntsub*sin(self.thetadto2)**2-self.sinkzvzdto2**2* \
                           (self.ntsub-sin(2*self.thetadto2)/sin(2*self.thetadto2/self.ntsub))))


        if 0:#self.spectral:
                  emK = self.FSpace
                  b0 = sqrt(1.-1./top.boost_gamma**2)
                  self.cut = 0.6
                  k = sqrt(emK.kx_unmod*emK.kx_unmod+emK.kz_unmod*emK.kz_unmod)
                  self.k_source_filter = where(k*self.dz/pi>self.cut*min(1.,self.dz/(b0*clight*top.dt)),0.,1.)
                  if self.l_getrho:emK.add_Sfilter('rho',self.k_source_filter)
                  emK.add_Sfilter('jx',self.k_source_filter)
                  emK.add_Sfilter('jy',self.k_source_filter)
                  emK.add_Sfilter('jz',self.k_source_filter)
        if(self.full_pxr == False):
          if self.spectral:
              kwPML = kwGPSTD
              if s.ntsub==inf:
                  GPSTD_PML = gpstd.PSATD_Maxwell_PML
              else:
                  GPSTD_PML = gpstd.GPSTD_Maxwell_PML
              # --- sides
              if b.xlbnd==openbc: s.xlPML = GPSTD_PML(syf=b.sidexl.syf,**kwPML)
              if b.xrbnd==openbc: s.xrPML = GPSTD_PML(syf=b.sidexr.syf,**kwPML)
              if b.ylbnd==openbc: s.ylPML = GPSTD_PML(syf=b.sideyl.syf,**kwPML)
              if b.yrbnd==openbc: s.yrPML = GPSTD_PML(syf=b.sideyr.syf,**kwPML)
              if b.zlbnd==openbc: s.zlPML = GPSTD_PML(syf=b.sidezl.syf,**kwPML)
              if b.zrbnd==openbc: s.zrPML = GPSTD_PML(syf=b.sidezr.syf,**kwPML)
              # --- edges
              if(b.xlbnd==openbc and b.ylbnd==openbc): s.xlylPML = GPSTD_PML(syf=b.edgexlyl.syf,**kwPML)
              if(b.xrbnd==openbc and b.ylbnd==openbc): s.xrylPML = GPSTD_PML(syf=b.edgexryl.syf,**kwPML)
              if(b.xlbnd==openbc and b.yrbnd==openbc): s.xlyrPML = GPSTD_PML(syf=b.edgexlyr.syf,**kwPML)
              if(b.xrbnd==openbc and b.yrbnd==openbc): s.xryrPML = GPSTD_PML(syf=b.edgexryr.syf,**kwPML)
              if(b.xlbnd==openbc and b.zlbnd==openbc): s.xlzlPML = GPSTD_PML(syf=b.edgexlzl.syf,**kwPML)
              if(b.xrbnd==openbc and b.zlbnd==openbc): s.xrzlPML = GPSTD_PML(syf=b.edgexrzl.syf,**kwPML)
              if(b.xlbnd==openbc and b.zrbnd==openbc): s.xlzrPML = GPSTD_PML(syf=b.edgexlzr.syf,**kwPML)
              if(b.xrbnd==openbc and b.zrbnd==openbc): s.xrzrPML = GPSTD_PML(syf=b.edgexrzr.syf,**kwPML)
              if(b.ylbnd==openbc and b.zlbnd==openbc): s.ylzlPML = GPSTD_PML(syf=b.edgeylzl.syf,**kwPML)
              if(b.yrbnd==openbc and b.zlbnd==openbc): s.yrzlPML = GPSTD_PML(syf=b.edgeyrzl.syf,**kwPML)
              if(b.ylbnd==openbc and b.zrbnd==openbc): s.ylzrPML = GPSTD_PML(syf=b.edgeylzr.syf,**kwPML)
              if(b.yrbnd==openbc and b.zrbnd==openbc): s.yrzrPML = GPSTD_PML(syf=b.edgeyrzr.syf,**kwPML)

              # --- corners
              if(b.xlbnd==openbc and b.ylbnd==openbc and b.zlbnd==openbc): s.xlylzlPML = GPSTD_PML(syf=b.cornerxlylzl.syf,**kwPML)
              if(b.xrbnd==openbc and b.ylbnd==openbc and b.zlbnd==openbc): s.xrylzlPML = GPSTD_PML(syf=b.cornerxrylzl.syf,**kwPML)
              if(b.xlbnd==openbc and b.yrbnd==openbc and b.zlbnd==openbc): s.xlyrzlPML = GPSTD_PML(syf=b.cornerxlyrzl.syf,**kwPML)
              if(b.xrbnd==openbc and b.yrbnd==openbc and b.zlbnd==openbc): s.xryrzlPML = GPSTD_PML(syf=b.cornerxryrzl.syf,**kwPML)
              if(b.xlbnd==openbc and b.ylbnd==openbc and b.zrbnd==openbc): s.xlylzrPML = GPSTD_PML(syf=b.cornerxlylzr.syf,**kwPML)
              if(b.xrbnd==openbc and b.ylbnd==openbc and b.zrbnd==openbc): s.xrylzrPML = GPSTD_PML(syf=b.cornerxrylzr.syf,**kwPML)
              if(b.xlbnd==openbc and b.yrbnd==openbc and b.zrbnd==openbc): s.xlyrzrPML = GPSTD_PML(syf=b.cornerxlyrzr.syf,**kwPML)
              if(b.xrbnd==openbc and b.yrbnd==openbc and b.zrbnd==openbc): s.xryrzrPML = GPSTD_PML(syf=b.cornerxryrzr.syf,**kwPML)


class Sorting:
  """
    Class Sorting

    Used to setup the sorting with picsars

    - activated: >0 sorting is activated
    - periods: list containing the sorting periods for each species
    - starts: first iteration before the start of the sorting
    - dx, dy, dz: the bin size normalized to the cell size. For instance, a dx of 1 corresponds to the cell dx.
    - xshift,yshift,zshift: shift of the sorting grid. The shift is normalized to dx,dy,dz. For instance a shift of 1 corresponds of 1 space step.

  """
  def __init__(self,periods,starts,activated=1,dx=1.,dy=1.,dz=1.,xshift=0.,yshift=0,zshift=0,verbose=False):
    self.activated = activated
    self.periods = periods
    self.starts = starts
    self.dx = dx
    self.dy = dy
    self.dz = dz
    self.xshift = xshift
    self.yshift = yshift
    self.zshift = zshift
    self.verbose = verbose
