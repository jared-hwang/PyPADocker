"""
Class for doing multigrid field solve on 2-D
--------------------------------------------
"""
from ..warp import *
from find_mgparam import find_mgparam

try:
    import psyco
except ImportError:
    pass

##############################################################################
class MultiGrid2D(MultiGrid3D):
    """
  2-D field solver, axisymmetric and slab, based on the 2-D solver in f3d_mgrid.F
    """

    def __init__(self,lreducedpickle=1,**kw):
        kw['lreducedpickle'] = lreducedpickle
        self.grid_overlap = 2

        # --- Force ny (which is not used here)
        self.ny = 0

        SubcycledPoissonSolver.__init__(self,kwdict=kw)
        if (self.solvergeom != w3d.RZgeom and self.solvergeom != w3d.XZgeom):
            self.solvergeom = w3d.RZgeom
        self.ncomponents = 1
        self.nyguardphi = 0
        self.nyguardrho = 0
        self.nyguarde   = 0

        # --- Make sure that the bounds have acceptable values.
        assert 0 <= min(self.bounds) and max(self.bounds) <= 2,"The boundary conditions have an incorrect value. They must be one of dirichlet, neumann or periodic."

        # --- Kludge - make sure that the multigrid3df routines never sets up
        # --- any conductors. This is not really needed here.
        f3d.gridmode = 1

        # --- Save input parameters
        self.processdefaultsfrompackage(MultiGrid2D.__w3dinputs__,w3d,kw)
        self.processdefaultsfrompackage(MultiGrid2D.__f3dinputs__,f3d,kw)

        # --- If there are any remaning keyword arguments, raise an error.
        assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

        # --- Check for consistency
        if self.solvergeom == w3d.RZgeom:
            assert self.xmmin >= 0.,"With RZgeom, xmmin must be >= 0."
            #self.xmmin = 0.

        # --- Create conductor objects
        self.initializeconductors()

        # --- Give these variables dummy initial values.
        self.mgiters = 0
        self.mgerror = 0.

    def getrho(self):
        'Returns the rho array without the guard cells'
        return self.source[self.nxguardrho:-self.nxguardrho or None,
                           0,
                           self.nzguardrho:-self.nzguardrho or None]

    def getrhop(self):
        'Returns the rhop array without the guard cells'
        return self.sourcep[self.nxguardrho:-self.nxguardrho or None,
                            0,
                            self.nzguardrho:-self.nzguardrho or None]

    def getphi(self):
        'Returns the phi array without the guard cells'
        return self.potential[self.nxguardphi:-self.nxguardphi or None,
                              0,
                              self.nzguardphi:-self.nzguardphi or None]

    def getphip(self):
        'Returns the phip array without the guard cells'
        return self.potentialp[self.nxguardphi:-self.nxguardphi or None,
                               0,
                               self.nzguardphi:-self.nzguardphi or None]

    def getselfe(self,*args,**kw):
        return super(MultiGrid2D,self).getselfe(*args,**kw)[:,:,0,:]

    def getselfep(self,*args,**kw):
        return super(MultiGrid2D,self).getselfep(*args,**kw)[:,:,0,:]

    def fetchpotentialfrompositions(self,x,y,z,phi):
        n = len(x)
        if n == 0: return
        if isinstance(self.potentialp,float): return
        if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
        else:                           r = x
        nxp = self.nxp + 2*self.nxguardphi
        nzp = self.nzp + 2*self.nzguardphi
        xmminp = self.xmminp - self.dx*self.nxguardphi
        xmmaxp = self.xmmaxp + self.dx*self.nxguardphi
        zmminp = self.zmminp - self.dz*self.nzguardphi + self.getzgridprv()
        zmmaxp = self.zmmaxp + self.dz*self.nzguardphi + self.getzgridprv()
        getgrid2d(n,r,z,phi,nxp,nzp,self.potentialp[:,0,:],
                  xmminp,xmmaxp,zmminp,zmmaxp)

    def fetchpotentialfsfrompositions(self,x,y,z,potential):
        'Fetches potential from the field solver grid'
        n = len(x)
        if n == 0: return
        if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
        else:                           r = x
        nxlocal = self.nxlocal + 2*self.nxguardphi
        nzlocal = self.nzlocal + 2*self.nzguardphi
        xmminlocal = self.xmminlocal - self.nxguardphi*self.dx
        xmmaxlocal = self.xmmaxlocal + self.nxguardphi*self.dx
        zmminlocal = self.zmminlocal - self.nzguardphi*self.dz
        zmmaxlocal = self.zmmaxlocal + self.nzguardphi*self.dz
        getgrid2d(n,r,z,potential,nxlocal,nzlocal,self.potential[:,0,:],
                  xmminlocal,xmmaxlocal,zmminlocal,zmmaxlocal)

    def dosolve(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        self.dosolvemultigrid(iwhich,zfact,isourcepndtscopies,indts,iselfb)
        #self.dosolvesuperlu(iwhich,*args)

    def dosolvemultigrid(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        if not self.l_internal_dosolve: return
        # --- set for longitudinal relativistic contraction
        if zfact is None:
            beta = top.pgroup.fselfb[iselfb]/clight
            zfact = 1./sqrt((1.-beta)*(1.+beta))
        else:
            beta =  sqrt( (1.-1./zfact)*(1.+1./zfact) )

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source
        if isinstance(self.potential,float): return

        mgverbose = self.getmgverbose()
        mgiters = zeros(1,'l')
        mgerror = zeros(1,'d')
        # --- This takes care of clear out the conductor information if needed.
        # --- Note that f3d.gridmode is passed in below - this still allows the
        # --- user to use the addconductor method if needed.
        if self.gridmode == 0: self.clearconductors([top.pgroup.fselfb[iselfb]])
        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
        self.lbuildquads = false
        #t0 = wtime()
        multigrid2dsolve(iwhich,self.nx,self.nz,self.nxlocal,self.nzlocal,
                         self.nxguardphi,self.nzguardphi,
                         self.nxguardrho,self.nzguardrho,
                         self.dx,self.dz*zfact,
                         self._phi[:,self.nyguardphi,:],
                         self._rho[:,self.nyguardrho,:],
                         self.bounds,self.xmminlocal,
                         self.mgparam,self.mgform,mgiters,self.mgmaxiters,
                         self.mgmaxlevels,mgerror,self.mgtol,mgverbose,
                         self.downpasses,self.uppasses,
                         self.lcndbndy,self.laddconductor,self.icndbndy,
                         f3d.gridmode,conductorobject,self.solvergeom==w3d.RZgeom,
                         false,self.fsdecomp)
        #t1 = wtime()
        #print "Multigrid time = ",t1-t0

        self.mgiters = mgiters[0]
        self.mgerror = mgerror[0]

    def dosolvesuperlu(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        "Solver using the SuperLU matrix solver instead of multigrid"
        if not self.l_internal_dosolve: return

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source
        if isinstance(self.potential,float): return

        mgiters = zeros(1,'l')
        mgerror = zeros(1,'d')
        # --- This takes care of clear out the conductor information if needed.
        # --- Note that f3d.gridmode is passed in below - this still allows the
        # --- user to use the addconductor method if needed.
        if self.gridmode == 0: self.clearconductors([top.pgroup.fselfb[iselfb]])
        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
        self.lbuildquads = false

        # --- Use direct matrix solver
        t0 = wtime()
        n = self.nxlocal*self.nzlocal
        nrhs = 1
        b = -self.source[:-1,self.nyguardrho,:-1]/eps0
        phi = self.potential[1:-1,self.nyguardphi,1:-1]
        info = zeros(1,'l')

        values = fzeros((5,n),'d')
        rowind = fzeros((5,n),'l')
        colptr = arange(n+1)*5 + 1
        rowcnt = zeros(n,'l')
        rmmin = self.xmminlocal
        dr = self.dx
        dz = self.dz
        drsqi = 1./dr**2
        dzsqi = 1./dz**2
        nxlocal = self.nxlocal
        nzlocal = self.nzlocal
        coeffikm1 = dzsqi
        coeffikp1 = dzsqi
        for iz in range(0,nzlocal):
            for ix in range(0,nxlocal):
                icol = iz*nxlocal + ix
                r = rmmin + ix*dr
                if r == 0.:
                    coeffik = - 4.*drsqi - 2.*dzsqi
                    coeffim1k = 0
                    coeffip1k = 4.*drsqi
                else:
                    coeffik = - 2.*drsqi - 2.*dzsqi
                    coeffim1k = (r-0.5*dr)/r*drsqi
                    coeffip1k = (r+0.5*dr)/r*drsqi
                    if ix == nxlocal-1:
                        b[ix,iz] += -coeffip1k*phi[ix+1,iz]
                        coeffip1k = 0.

                vtemp = [coeffikm1,coeffim1k,coeffik,coeffip1k,coeffikp1]
                rtemp = [-nxlocal+icol,-1+icol,0+icol,+1+icol,+nxlocal+icol]
                if rtemp[0] < 0:
                    # --- Periodic Z boundary condition
                    rtemp = rtemp[1:] + [rtemp[0] + nzlocal*nxlocal]
                    vtemp = vtemp[1:] + [vtemp[0]]
                if rtemp[0] < 0:
                    # --- Throw away point "below" r=0 axis at iz=0.
                    del rtemp[0]
                    del vtemp[0]
                if rtemp[-1] >= n:
                    # --- Periodic Z boundary condition
                    rtemp = [rtemp[-1] - nzlocal*nxlocal] + rtemp[:-1]
                    vtemp = [vtemp[-1]]         + vtemp[:-1]
                if rtemp[-1] >= n:
                    # --- Throw away point beyond r=nxlocal, iz=nzlocal
                    del rtemp[-1]
                    del vtemp[-1]

                for i in range(len(rtemp)):
                    irow = rowcnt[rtemp[i]]
                    values[irow,rtemp[i]] = vtemp[i]
                    rowind[irow,rtemp[i]] = icol + 1
                    rowcnt[rtemp[i]] += 1

        # --- There are two values of rowind that are unset, the (i-1) term
        # --- for (ix,iz)=(0,0) and the (i+1) term for (ix,iz)=(nx,nzlocal).
        # --- Give the first one a fake value (since the coefficient is zero
        # --- anway.
        rowind[-1,0] = rowind[-2,0] + 1
        # --- The other is ignored by decrementing the last value of colptr.
        colptr[-1] -= 1

        self.values = values
        self.rowind = rowind
        self.colptr = colptr
        nnz = colptr[-1] - 1

        t1 = wtime()
        superlu_dgssv(n,nnz,nrhs,values,rowind,colptr,b,info)
        t2 = wtime()

        self.potential[1:-2,0,1:-2] = b
        self.potential[0,0,1:-2] = self.potential[2,0,1:-2]
        self.potential[-1,0,1:-2] = 2*self.potential[-2,0,1:-2]-self.potential[-3,0,1:-2]
        self.potential[:,0,-2:] = self.potential[:,0,1:3]
        self.potential[:,0,0] = self.potential[:,0,-3]
        t3 = wtime()

        print "Solve time = ",t2 - t1
        print "Total time = ",t3 - t0
        self.fstime = t2 - t1
        self.tottime = t3 - t0

    ##########################################################################
    # Define the basic plot commands
    def pfzr(self,**kw): self.genericpf(kw,pfzx)
    def pfzrg(self,**kw): self.genericpf(kw,pfzxg)

    def getresidual(self):
        res = zeros(shape(self._phi),'d')
        dxsqi  = 1./self.dx**2
        dzsqi  = 1./self.dz**2
        xminodx = self.xmminlocal/self.dx
        rho = self._rho/eps0
        conductorobject = self.getconductorobject()
        residual2d(self.nxlocal,self.nzlocal,
                   self.nxguardphi,self.nzguardphi,
                   self.nxguardrho,self.nzguardrho,
                   self.nxguardphi,self.nzguardphi,
                   dxsqi,dzsqi,xminodx,self.solvergeom==w3d.RZgeom,false,
                   self._phi[:,self.nyguardphi,:],rho[:,self.nyguardrho,:],
                   res[:,self.nyguardphi,:],0,self.bounds,
                   self.mgform,true,self.lcndbndy,self.icndbndy,conductorobject)
        return res

    def getimagecharges(self, includeboundaries=False, iselfb=0):
        """This calculates the image charges inside of any conductors.
        This is a bit of a hack. It calculates the residual, but turning off
        the zeroing out of the residual inside any conductors and on the boundaries."""
        if includeboundaries:
            # --- Normally, with Dirichlet boundaries, the phi is linearly extrapolated into
            # --- the guard cells since this gives better behavior when fetching the E fields.
            # --- However, this makes the residual zero. This call fills the guard cells with the
            # --- potential on the boundary. Also set bounds so that no boundary conditions are
            # --- applied to the residual.
            applyboundaryconditions3d(self.nxlocal,self.nylocal,self.nzlocal,
                                      self.nxguardphi,self.nyguardphi,self.nzguardphi,
                                      self._phi,1,self.bounds,false,true)
            bounds = [-1,-1,-1,-1,-1,-1]
        else:
            bounds = self.bounds

        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
        istartsave = conductorobject.interior.istart.copy()
        conductorobject.interior.istart = 1

        dxsqi  = 1./self.dx**2
        dzsqi  = 1./self.dz**2
        xminodx = self.xmminlocal/self.dx
        rho = self._rho/eps0
        result = zeros(shape(self._phi),'d')
        residual2d(self.nxlocal,self.nzlocal,
                   self.nxguardphi,self.nzguardphi,
                   self.nxguardrho,self.nzguardrho,
                   self.nxguardphi,self.nzguardphi,
                   dxsqi,dzsqi,xminodx,self.solvergeom==w3d.RZgeom,false,
                   self._phi[:,self.nyguardphi,:],rho[:,self.nyguardrho,:],
                   result[:,self.nyguardphi,:],0,bounds,
                   self.mgform,true,self.lcndbndy,self.icndbndy,conductorobject)

        conductorobject.interior.istart[:] = istartsave
        if includeboundaries:
            # --- Undo the applyboundaryconditions3d from above.
            applyboundaryconditions3d(self.nxlocal,self.nylocal,self.nzlocal,
                                      self.nxguardphi,self.nyguardphi,self.nzguardphi,
                                      self._phi,1,self.bounds,true,false)

        # --- Remove the premultiplying factor
        result *= eps0
        return result

##############################################################################
class MultiGridRZ(MultiGrid2D):
    """
  2-D axisymmetric field solver, based on the 2-D solver in f3d_mgrid.F
    """
    def __init__(self, **kw):
        self.solvergeom = w3d.RZgeom
        MultiGrid2D.__init__(self, **kw)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MultiGrid2DDielectric(MultiGrid2D):
    """
  2-D solver allowing for spatially varying dielectric
    """

    def __init__(self,epsilon=None,lreducedpickle=1,**kw):
        MultiGrid2D.__init__(self,lreducedpickle,**kw)

        if epsilon is None:
            self.epsilon = eps0*fones((self.nxlocal+2,self.nzlocal+2),'d')
        else:
            self.epsilon = epsilon

    def dosolve(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        if not self.l_internal_dosolve: return
        assert self.epsilon is not None,"epsilon must be defined"

        # --- set for longitudinal relativistic contraction
        if zfact is None:
            beta = top.pgroup.fselfb[iselfb]/clight
            zfact = 1./sqrt((1.-beta)*(1.+beta))
        else:
            beta =  sqrt( (1.-1./zfact)*(1.+1./zfact) )

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source
        if isinstance(self.potential,float): return

        mgverbose = self.getmgverbose()
        mgiters = zeros(1,'l')
        mgerror = zeros(1,'d')
        # --- This takes care of clear out the conductor information if needed.
        # --- Note that f3d.gridmode is passed in below - this still allows the
        # --- user to use the addconductor method if needed.
        if self.gridmode == 0: self.clearconductors([top.pgroup.fselfb[iselfb]])
        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
        multigrid2ddielectricsolve(iwhich,self.nx,self.nz,self.nxlocal,self.nzlocal,
                         self.nxguardphi,self.nzguardphi,
                         self.nxguardrho,self.nzguardrho,
                         self.dx,self.dz*zfact,
                         self._phi[:,self.nyguardphi,:],
                         self._rho[:,self.nyguardrho,:],
                         self.epsilon,self.bounds,
                         self.xmminlocal*zfact,
                         self.mgparam,mgiters,self.mgmaxiters,
                         self.mgmaxlevels,mgerror,self.mgtol,mgverbose,
                         self.downpasses,self.uppasses,
                         self.lcndbndy,self.laddconductor,
                         f3d.gridmode,conductorobject,self.solvergeom==w3d.RZgeom,
                         self.fsdecomp)

        self.mgiters = mgiters[0]
        self.mgerror = mgerror[0]

    def getresidual(self):
        res = zeros(shape(self._phi),'d')
        conductorobject = self.getconductorobject()
        residual2ddielectric(self.nxlocal,self.nzlocal,
                             self.nxguardphi,self.nzguardphi,
                             self.nxguardrho,self.nzguardrho,
                             self.nxguardphi,self.nzguardphi,
                             self._phi,rho,self.epsilon,res,
                             self.dx,self.dz,0,self.bounds,conductorobject)
        return res

##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MultiGridImplicit2D(MultiGrid3D):
    """
  This solves the modified Poisson equation which includes the suseptibility
  tensor that appears from the direct implicit scheme.
  It currently uses the generic sparse matrix solver SuperLU. The various
  multigrid input parameters are maintained for future use, but are ignored now.

  Initially, conductors are not implemented.
    """

    def __init__(self,lreducedpickle=1,withbadvance=true,**kw):
        kw['lreducedpickle'] = lreducedpickle
        self.withbadvance = withbadvance
        self.grid_overlap = 2

        # --- Force ny (which is not used here)
        self.ny = 0

        SubcycledPoissonSolver.__init__(self,kwdict=kw)
        if (self.solvergeom != w3d.RZgeom and self.solvergeom != w3d.XZgeom):
            self.solvergeom = w3d.RZgeom
        self.ncomponents = 1
        self.nyguardphi = 0
        self.nyguardrho = 0
        self.nyguarde   = 0

        # --- Make sure that the bounds have acceptable values.
        assert 0 <= min(self.bounds) and max(self.bounds) <= 2,"The boundary conditions have an incorrect value. They must be one of dirichlet, neumann or periodic."

        # --- Kludge - make sure that the multigrid3df routines never sets up
        # --- any conductors. This is not really needed here.
        f3d.gridmode = 1

        # --- Save input parameters
        self.processdefaultsfrompackage(MultiGrid2D.__w3dinputs__,w3d,kw)
        self.processdefaultsfrompackage(MultiGrid2D.__f3dinputs__,f3d,kw)

        # --- If there are any remaning keyword arguments, raise an error.
        assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

        # --- Create conductor objects
        self.initializeconductors()

        # --- Give these variables dummy initial values.
        self.mgiters = 0
        self.mgerror = 0.

        # --- Turn on the chi kludge, where chi is set to be an average value
        # --- of chi for grid cells where is it zero.
        self.chikludge = 1

    def __getstate__(self):
        dict = MultiGrid3D.__getstate__(self)
        if self.lreducedpickle:
            if 'chi0' in dict: del dict['chi0']
        return dict

    def getpdims(self):
        # --- This is needed to set the top.nsimplicit variable.
        setupImplicit(top.pgroup)
        dims = MultiGrid3D.getpdims(self)
        # --- The extra dimension is to hold the charge density and the chi's
        # --- for the implicit groups.
        dims = (tuple(list(dims[0])+[1+top.nsimplicit]),)+dims[1:]
        return dims

    def getdims(self):
        # --- This is needed to set the top.nsimplicit variable.
        setupImplicit(top.pgroup)
        dims = MultiGrid3D.getdims(self)
        # --- The extra dimension is to hold the charge density and the chi's
        # --- for the implicit groups.
        dims = (tuple(list(dims[0])+[1+top.nsimplicit]),)+dims[1:]
        return dims

    def getrho(self):
        return self.source[:,0,:,0]

    def getrhop(self):
        return self.sourcep[:,0,:,0]

    def getphi(self):
        'Returns the phi array without the guard cells'
        return MultiGrid3D.getphi(self)[:,0,:]

    def getphip(self):
        'Returns the phip array without the guard cells'
        return MultiGrid3D.getphip(self)[:,0,:]

    def getselfe(self,*args,**kw):
        return super(MultiGridImplicit2D,self).getselfe(*args,**kw)[:,:,0,:]

    def getselfep(self,*args,**kw):
        return super(MultiGridImplicit2D,self).getselfep(*args,**kw)[:,:,0,:]

    # --- A special version is needed since only part if source is returned.
    def _setuprhoproperty():
        doc = "Charge density array"
        def fget(self):
            return self.returnsource(0,0)[...,0]
        def fset(self,value):
            self.returnsource(0,0)[...,0] = value
        return locals()
    rho = property(**_setuprhoproperty())
    del _setuprhoproperty

    # --- A special version is needed since only part if sourcep is returned.
    def _setuprhopproperty():
        doc = "Charge density array for particles"
        def fget(self):
            return self.returnsourcep(0,0,0)[...,0]
        def fset(self,value):
            self.returnsourcep(0,0,0)[...,0] = value
        return locals()
    rhop = property(**_setuprhopproperty())
    del _setuprhopproperty

    def loadrho(self,lzero=None,lfinalize_rho=None,**kw):
        # --- top.laccumulate_rho is used as a flag by the implicit stepper.
        # --- When true, the load rho is skipped - it is not needed at some
        # --- points during a step.
        if top.laccumulate_rho: return
        MultiGrid3D.loadsource(self,lzero,lfinalize_rho,**kw)

    def fetche(self,*args,**kw):
        # --- lresetparticlee is used as a flag in the implicit stepper.
        # --- When false, skip the fetche since the field is calculated
        # --- from existing data.
        if not top.lresetparticlee: return
        MultiGrid3D.fetchfield(self,*args,**kw)

    def setsourcep(self,js,pgroup,zgrid):
        n  = pgroup.nps[js]
        if n == 0: return
        i  = pgroup.ins[js] - 1
        x  = pgroup.xp[i:i+n]
        y  = pgroup.yp[i:i+n]
        z  = pgroup.zp[i:i+n]
        ux = zeros((0,), 'd')
        uy = zeros((0,), 'd')
        uz = pgroup.uzp[i:i+n]
        gaminv = zeros((0,), 'd')
        q  = pgroup.sq[js]
        m  = pgroup.sm[js]
        w  = pgroup.sw[js]*top.pgroup.dtscale[js]
        iimp = pgroup.iimplicit[js]
        if top.wpid == 0: wght = zeros((0,), 'd')
        else:             wght = pgroup.pid[i:i+n,top.wpid-1]
        depos_order = top.depos_order[:,js]
        self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wght,zgrid,q,m,w,iimp,
                                  depos_order)

    def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wght,zgrid,q,m,w,iimp,
                             depos_order):
        n  = len(x)
        if n == 0: return
        # --- Create a temporary array to pass into setrho3d. This contributes
        # --- differently to the charge density and to chi. Also, make it a
        # --- 3-D array so it is accepted by setrho3d.
        sourcep = fzeros(self.sourcep.shape[:-1],'d')
        if top.wpid == 0:
            setrho3d(sourcep,n,x,y,z,zgrid,q,w,top.depos,depos_order,
                     self.nxp,self.nyp,self.nzp,
                     self.nxguardrho,self.nyguardrho,self.nzguardrho,
                     self.dx,1.,self.dz,
                     self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                     self.solvergeom==w3d.RZgeom)
        else:
            # --- Need top.pid(:,top.wpid)
            setrho3dw(sourcep,n,x,y,z,zgrid,wght,q,w,top.depos,depos_order,
                      self.nxp,self.nyp,self.nzp,
                      self.nxguardrho,self.nyguardrho,self.nzguardrho,
                      self.dx,1.,self.dz,
                      self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                      self.solvergeom==w3d.RZgeom)
        self.sourcep[...,0] += sourcep
        if iimp >= 0:
            # --- The extra terms convert rho to chi
            self.sourcep[...,iimp+1] += 0.5*sourcep*q/m*top.dt**2/eps0

    def setsourceforfieldsolve(self,*args):
        # --- A separate copy is needed since self.source has an extra dimension
        # --- which must be looped over.
        SubcycledPoissonSolver.setsourceforfieldsolve(self,*args)
        if self.lparallel:
            SubcycledPoissonSolver.setsourcepforparticles(self,*args)
            if isinstance(self.source,float): return
            if isinstance(self.sourcep,float): return
            for iimp in range(1+top.nsimplicit):
                setrhoforfieldsolve3d(self.nxlocal,self.nylocal,self.nzlocal,
                                      self.source[...,iimp],
                                      self.nxp,self.nyp,self.nzp,self.sourcep[...,iimp],
                                      self.nxguardrho,self.nyguardrho,self.nzguardrho,
                                      self.fsdecomp,self.ppdecomp)

    def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
        MultiGrid3D.fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js,pgroup)
        # --- Force ey to zero (is this really needed?)
        #ey[...] = 0.

    def fetchpotentialfrompositions(self,x,y,z,potential):
        n = len(x)
        if n == 0: return
        if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
        else:                           r = x
        nxp = self.nxp + 2*self.nxguardphi
        nzp = self.nzp + 2*self.nzguardphi
        xmminp = self.xmminp - self.nxguardphi*self.dx
        xmmaxp = self.xmmaxp + self.nxguardphi*self.dx
        zmminp = self.zmminp - self.nzguardphi*self.dz + self.getzgridprv()
        zmmaxp = self.zmmaxp + self.nzguardphi*self.dz + self.getzgridprv()
        getgrid2d(n,r,z,potential,nxp,nzp,self.potentialp[:,0,:],
                  xmminp,xmmaxp,zmminp,zmmaxp)

    def removedelsqphi(self):
        phi = self.potential[:,0,:]
        rho = self.source[:,0,:,0]
        if self.solvergeom==w3d.RZgeom:
            rr = self.xmminlocal + arange(self.nxlocal+1)*self.dx
            rho += eps0*(
              +(+phi[:-2,1:-1]*((rr-0.5*self.dx)/rr)[:,newaxis]
                -2.*phi[1:-1,1:-1]
                +phi[2:,1:-1]*((rr+0.5*self.dx)/rr)[:,newaxis])/self.dx**2
              +(phi[1:-1,:-2] - 2.*phi[1:-1,1:-1] + phi[1:-1,2:])/self.dz**2)
        else:
            rho += eps0*(
              +(phi[:-2,1:-1] - 2.*phi[1:-1,1:-1] + phi[2:,1:-1])/self.dx**2
              +(phi[1:-1,:-2] - 2.*phi[1:-1,1:-1] + phi[1:-1,2:])/self.dz**2)

    def dosolve(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        if not self.l_internal_dosolve: return
        # --- Do the solve, including chi
        #self.dosolvesuperlu(iwhich,*args)
        self.dosolvemg(iwhich,zfact,isourcepndtscopies,indts,iselfb)

    def dosolvemg(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        # --- set for longitudinal relativistic contraction
        if zfact is None:
            beta = top.pgroup.fselfb[iselfb]/clight
            zfact = 1./sqrt((1.-beta)*(1.+beta))
        else:
            beta =  sqrt( (1.-1./zfact)*(1.+1./zfact) )

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source[...,0]
        if isinstance(self.potential,float): return

        mgiters = zeros(1,'l')
        mgerror = zeros(1,'d')
        # --- This takes care of clear out the conductor information if needed.
        # --- Note that f3d.gridmode is passed in below - this still allows the
        # --- user to use the addconductor method if needed.
        if self.gridmode == 0: self.clearconductors([top.pgroup.fselfb[iselfb]])
        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])

        # --- Setup implicit chi
        qomdt = top.implicitfactor*top.dt # implicitfactor = q/m
        #--- chi0 = 0.5*rho*q/m*top.dt**2/eps0
        self.chi0 = self.source[...,1:]
        # --- Kludge alart!!!
        if self.chikludge:
            for js in range(self.source.shape[-1]-1):
                if maxnd(abs(self.chi0[...,js])) == 0.: continue
                avechi = sumnd(self.chi0[...,js])/sumnd(where(self.chi0[...,js] == 0.,0.,1.))
                self.chi0[...,js] = where(self.chi0[...,js]==0.,avechi,self.chi0[...,js])
        """
        # --- Test a linearly varying chi and parabolic phi
        c1 = 10.
        c2 = 2.
        alpha = 10.
        for iz in range(self.nzlocal+1):
          self.chi0[...,iz] = (c1 + c2*self.zmesh[iz])
          self.source[...,iz] = -(2.*alpha + 2.*c1*alpha + 4.*c2*alpha*w3d.zmesh[iz])*eps0
        """

        mgverbose = self.getmgverbose()
        mgsolveimplicites2d(iwhich,self.nx,self.nz,self.nxlocal,self.nzlocal,
                            self.dx,self.dz*zfact,
                            self.nxguardphi,self.nzguardphi,
                            self.nxguardrho,self.nzguardrho,
                            self.potential,self.source,
                            top.nsimplicit,qomdt,self.chi0,self.withbadvance,
                            self.bounds,self.xmminlocal,self.zmminlocal*zfact,
                            self.getzgrid()*zfact,
                            self.mgparam,mgiters,self.mgmaxiters,
                            self.mgmaxlevels,mgerror,self.mgtol,mgverbose,
                            self.downpasses,self.uppasses,
                            self.lcndbndy,self.laddconductor,self.icndbndy,
                            f3d.gridmode,conductorobject,
                            self.solvergeom==w3d.RZgeom,self.fsdecomp)

        self.mgiters = mgiters[0]
        self.mgerror = mgerror[0]

    def dosolvesuperlu(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        "Note that this does not actually include the implicit susecptibility"
        #self.grid.rho = self.source
        #self.grid.phi = self.potential
        #solve_mgridrz(self.grid,self.mgtol,false)
        #self.mgiters = nb_iters
        ##self.mgerror = mgerror[0] # not saved anywhere
        self.getconductorobject()

        # --- Use direct matrix solver
        t0 = wtime()
        n = self.nxlocal*self.nzlocal
        nrhs = 1
        b = -self.source[:-1,:-1]/eps0
        phi = self.potential[1:-1,1:-1]
        info = zeros(1,'l')

        values = fzeros((5,n),'d')
        rowind = fzeros((5,n),'l')
        colptr = arange(n+1)*5 + 1
        rowcnt = zeros(n,'l')
        rmmin = self.xmminlocal
        dr = self.dx
        dz = self.dz
        drsqi = 1./dr**2
        dzsqi = 1./dz**2
        nr = self.nxlocal
        nzlocal = self.nzlocal
        coeffikm1 = dzsqi
        coeffikp1 = dzsqi
        for iz in range(0,nzlocal):
            for ix in range(0,nr):
                icol = iz*nr + ix
                r = rmmin + ix*dr
                if r == 0.:
                    coeffik = - 4.*drsqi - 2.*dzsqi
                    coeffim1k = 0
                    coeffip1k = 4.*drsqi
                else:
                    coeffik = - 2.*drsqi - 2.*dzsqi
                    coeffim1k = (r-0.5*dr)/r*drsqi
                    coeffip1k = (r+0.5*dr)/r*drsqi
                    if ix == nr-1:
                        b[ix,iz] += -coeffip1k*phi[ix+1,iz]
                        coeffip1k = 0.

                vtemp = [coeffikm1,coeffim1k,coeffik,coeffip1k,coeffikp1]
                rtemp = [-nr+icol,-1+icol,0+icol,+1+icol,+nr+icol]
                if rtemp[0] < 0:
                    # --- Periodic Z boundary condition
                    rtemp = rtemp[1:] + [rtemp[0] + nzlocal*nr]
                    vtemp = vtemp[1:] + [vtemp[0]]
                if rtemp[0] < 0:
                    # --- Throw away point "below" r=0 axis at iz=0.
                    del rtemp[0]
                    del vtemp[0]
                if rtemp[-1] >= n:
                    # --- Periodic Z boundary condition
                    rtemp = [rtemp[-1] - nzlocal*nr] + rtemp[:-1]
                    vtemp = [vtemp[-1]]         + vtemp[:-1]
                if rtemp[-1] >= n:
                    # --- Throw away point beyond r=nr, iz=nzlocal
                    del rtemp[-1]
                    del vtemp[-1]

                for i in range(len(rtemp)):
                    irow = rowcnt[rtemp[i]]
                    values[irow,rtemp[i]] = vtemp[i]
                    rowind[irow,rtemp[i]] = icol + 1
                    rowcnt[rtemp[i]] += 1

        # --- There are two values of rowind that are unset, the (i-1) term
        # --- for (ix,iz)=(0,0) and the (i+1) term for (ix,iz)=(nx,nzlocal).
        # --- Give the first one a fake value (since the coefficient is zero
        # --- anway.
        rowind[-1,0] = rowind[-2,0] + 1
        # --- The other is ignored by decrementing the last value of colptr.
        colptr[-1] -= 1

        self.values = values
        self.rowind = rowind
        self.colptr = colptr
        nnz = colptr[-1] - 1

        t1 = wtime()
        superlu_dgssv(n,nnz,nrhs,values,rowind,colptr,b,info)
        t2 = wtime()

        self.potential[1:-2,1:-2] = b
        self.potential[0,1:-2] = self.potential[2,1:-2]
        self.potential[-1,1:-2] = 2*self.potential[-2,1:-2]-self.potential[-3,1:-2]
        self.potential[:,-2:] = self.potential[:,1:3]
        self.potential[:,0] = self.potential[:,-3]
        t3 = wtime()

        print "Solve time = ",t2 - t1
        print "Total time = ",t3 - t0
        self.fstime = t2 - t1
        self.tottime = t3 - t0

    ##########################################################################
    # Define the basic plot commands
    def pfzr(self,**kw): self.genericpf(kw,pfzx)
    def pfzrg(self,**kw): self.genericpf(kw,pfzxg)


##############################################################################
class MultiGridRZunsupported(MultiGrid3D):
    """
  2-D field solver, axisymmetric and slab, based on the 2-D solver in frz_mgrid.F
    """

    def __init__(self,lreducedpickle=1,**kw):
        kw['lreducedpickle'] = lreducedpickle
        self.grid_overlap = 2

        # --- Force ny (which is not used here)
        self.ny = 0

        SubcycledPoissonSolver.__init__(self,kwdict=kw)
        if (self.solvergeom != w3d.RZgeom and self.solvergeom != w3d.XZgeom):
            self.solvergeom = w3d.RZgeom
        self.ncomponents = 1
        self.nyguardphi = 0
        self.nyguardrho = 0
        self.nyguarde   = 0

        # --- Make sure that the bounds have acceptable values.
        assert 0 <= min(self.bounds) and max(self.bounds) <= 2,"The boundary conditions have an incorrect value. They must be one of dirichlet, neumann or periodic."

        # --- Kludge - make sure that the multigrid3df routines never sets up
        # --- any conductors. This is not really needed here.
        f3d.gridmode = 1

        # --- Save input parameters
        self.processdefaultsfrompackage(MultiGrid2D.__w3dinputs__,w3d,kw)
        self.processdefaultsfrompackage(MultiGrid2D.__f3dinputs__,f3d,kw)

        # --- If there are any remaning keyword arguments, raise an error.
        assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

        # --- Create conductor objects
        self.initializeconductors()

        # --- Give these variables dummy initial values.
        self.mgiters = 0
        self.mgerror = 0.

        self.initializegrid()

    def initializegrid(self):
        # --- Initialize the grid object
        self.grid = GRIDtype()
        init_gridrz(self.grid,self.nx,self.nz,self.dx,self.dz,self.xmmin,self.zmmin,
                    self.lparallel,self.bounds[1],self.bounds[4],self.bounds[5])

    def __getstate__(self):
        dict = MultiGrid3D.__getstate__(self)
        if self.lreducedpickle:
            if 'grid' in dict: del dict['grid']
        return dict

    def __setstate__(self,dict):
        MultiGrid3D.__setstate__(self,dict)
        if self.lreducedpickle and not self.lnorestoreonpickle:
            self.initializegrid()

    def getrho(self):
        'Returns the rho array without the guard cells'
        return self.source[self.nxguardrho:-self.nxguardrho or None,
                           0,
                           self.nzguardrho:-self.nzguardrho or None]

    def getrhop(self):
        'Returns the rhop array without the guard cells'
        return self.sourcep[self.nxguardrho:-self.nxguardrho or None,
                            0,
                            self.nzguardrho:-self.nzguardrho or None]

    def getphi(self):
        'Returns the phi array without the guard cells'
        return self.potential[self.nxguardphi:-self.nxguardphi or None,
                              0,
                              self.nzguardphi:-self.nzguardphi or None]

    def getphip(self):
        'Returns the phip array without the guard cells'
        return self.potentialp[self.nxguardphi:-self.nxguardphi or None,
                               0,
                               self.nzguardphi:-self.nzguardphi or None]

    def getselfe(self,*args,**kw):
        return super(MultiGridRZ,self).getselfe(*args,**kw)[:,:,0,:]

    def getselfep(self,*args,**kw):
        return super(MultiGridRZ,self).getselfep(*args,**kw)[:,:,0,:]

    def fetchpotentialfrompositions(self,x,y,z,phi):
        n = len(x)
        if n == 0: return
        if isinstance(self.potentialp,float): return
        if self.solvergeom==w3d.RZgeom: r = sqrt(x**2 + y**2)
        else:                           r = x
        nxp = self.nxp + 2*self.nxguardphi
        nzp = self.nzp + 2*self.nzguardphi
        xmminp = self.xmminp - self.dx*self.nxguardphi
        xmmaxp = self.xmmaxp + self.dx*self.nxguardphi
        zmminp = self.zmminp - self.dz*self.nzguardphi + self.getzgridprv()
        zmmaxp = self.zmmaxp + self.dz*self.nzguardphi + self.getzgridprv()
        getgrid2d(n,r,z,phi,nxp,nzp,self.potentialp[:,0,:],
                  xmminp,xmmaxp,zmminp,zmmaxp)

    def _installconductor(self,conductorobject,installedlist,conductordata,
                          fselfb):
        # --- Copied from MultiGrid3D, but turns on installrz.
        # --- This does that actual installation of the conductor into the
        # --- conductor object

        # --- Extract the data from conductordata (the arguments to
        # --- installconductor)
        conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill = conductordata

        if conductor in installedlist: return
        installedlist.append(conductor)

        nx,ny,nz = self.nx,self.ny,self.nz
        if fselfb == 'p':
            zscale = 1.
            nxlocal,nylocal,nzlocal = self.nxp,self.nyp,self.nzp
            mgmaxlevels = 1
            decomp = self.ppdecomp
        else:
            # --- Get relativistic longitudinal scaling factor
            # --- This is quite ready yet.
            beta = fselfb/clight
            zscale = 1./sqrt((1.-beta)*(1.+beta))
            nxlocal,nylocal,nzlocal = self.nxlocal,self.nylocal,self.nzlocal
            mgmaxlevels = None
            decomp = self.fsdecomp

        xmmin,xmmax = self.xmmin,self.xmmax
        ymmin,ymmax = self.ymmin,self.ymmax
        zmmin,zmmax = self.zmmin,self.zmmax
        installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                          self.getzgrid(),
                          nx,ny,nz,nxlocal,nylocal,nzlocal,
                          xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,
                          zscale,self.l2symtry,self.l4symtry,
                          installrz=1,
                          solvergeom=self.solvergeom,conductors=conductorobject,
                          mgmaxlevels=mgmaxlevels,
                          gridrz=self.grid)
        get_cond_rz_grid(self.grid,conductorobject)

    def setconductorvoltage(self,voltage,condid=0,discrete=false,
                            setvinject=false):
        'calls setconductorvoltage'

        # --- Set vinject first is requested.
        # --- This is copied from plot_conductors.setconductorvoltage.
        if setvinject:
            if instance(voltage,(list,tuple,ndarray)):
                # --- Set it to the voltage on the left edge
                top.vinject = voltage[0]
            elif callable(voltage):
                # --- Set it to the voltage at the source center
                top.vinject = voltage(top.xinject,top.yinject,top.zinject)
            else:
                top.vinject = voltage

        # --- Loop over all of the selfb groups to that all conductor objects
        # --- are handled.
        # --- XXX NOTE THAT SELFB IS NOT IMPLEMENTED YET FOR MultiGridRZ XXX
        for iselfb in range(top.nsselfb):
            if isinstance(voltage,(list,tuple,ndarray)):
            # --- Voltage is assumed to be the voltages are the z grid cell locations
            # --- (in the global beam frame).
                setconductorvoltagerz_grid(self.grid,voltage,self.nz,self.zmmin,
                                           self.dz,discrete,condid)
            else:
                setconductorvoltagerz_id_grid(self.grid,condid,voltage)

    def dosolve(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        if not self.l_internal_dosolve: return

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source
        if isinstance(self.potential,float): return

#   if self.izfsslave is None: self.izfsslave = top.izfsslave
#   if self.nzfsslave is None: self.nzfsslave = top.nzfsslave
        mgiters = zeros(1,'l')
        mgerror = zeros(1,'d')

        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])

        self.grid.rho = self.source[:,0,:]
        self.grid.phi = self.potential[:,0,:]
        solve_mgridrz(self.grid,self.mgtol,false)

        self.mgiters = frz.nb_iters
        self.mgerror = frz.maxerr

    ##########################################################################
    # Define the basic plot commands
    def pfzr(self,**kw): self.genericpf(kw,pfzx)
    def pfzrg(self,**kw): self.genericpf(kw,pfzxg)

# --- This can only be done after MultiGridRZ and MultiGridImplicit2D are defined.
try:
    psyco.bind(MultiGridRZ)
    psyco.bind(MultiGridImplicit2D)
except NameError:
    pass
