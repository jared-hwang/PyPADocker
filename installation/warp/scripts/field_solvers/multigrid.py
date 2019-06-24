"""
Classed for doing multigrid field solve in 3-D
----------------------------------------------
"""
# ToDo:
#  - modify setrhop to check if particles are within grid
#  - incorporate instances into the particle mover, so charge is deposited and
#    the E fields gather appropriately.
from ..warp import *
from find_mgparam import find_mgparam

try:
    import psyco
except ImportError:
    pass


##############################################################################
class MultiGrid3D(SubcycledPoissonSolver):
    """
  Creates a 3-D multigrid Poisson solver. There are many input parameters,
  most of which get there default values from one of the fortran packages.

  | Input parameters from w3d:
  |    nx,ny,nz,dx,dy,dz,nxlocal,nylocal,nzlocal,
  |    nxguardphi,nyguardphi,nzguardphi,
  |    nxguardrho,nyguardrho,nzguardrho,
  |    nxguarde,nyguarde,nzguarde,
  |    xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,
  |    xmminlocal,xmmaxlocal,
  |    ymminlocal,ymmaxlocal,
  |    zmminlocal,zmmaxlocal,
  |    bound0,boundnz,boundxy,l2symtry,l4symtry,
  |    solvergeom,
  |    iondensity,electrontemperature,plasmapotential,
  |    electrondensitymaxscale

  | Input parameters from top:
  |    pbound0,pboundnz,pboundxy,
  |    nprocs,nxprocs,nyprocs,nzprocs,
  |    lfsautodecomp,zslave,debug

  | Input parameters from f3d:
  |    gridmode,mgparam,downpasses,uppasses,
  |    mgmaxiters,mgtol,mgmaxlevels,mgform,
  |    mgverbose,mgntverbose,
  |    lcndbndy,icndbndy,laddconductor,lprecalccoeffs

  | Other input paramters:
  |    forcesymmetries=1: When true and either 2 or 4 fold symmetry are set,
  |                       the x and/or y mins are forced to zero.
  |    lreducedpickle=1: When true, when the instance is pickled, the large
  |                      arrays are not included. Experts only.
  |    lnorestoreonpickle=0: When false, when the instance is unpickled, the
  |                          large arrays are restored. Experts only.
  |    ldosolve=1: Flags sets whether a field solve is done.
  |                When true, do the solve, otherwise don't, and also don't
  |                do the charge deposition.
  |    l_internal_dosolve=1: Another flag which sets whether a field solve is done.
  |                          With this flag, if false, the charge deposition is
  |                          still done.
  |    gridvz=None: An option grid velocity, independent of top.vbeam.
  |                 Only used for special purposes.
  |    lchild=False: Internally used flag, true when the instance is a child
  |                  relative to a root grid when doing mesh refinement.
  |    userfsdecompnx=None: User specified decomposition for the field solver in x
  |    userfsdecompny=None: User specified decomposition for the field solver in y
  |    userfsdecompnz=None: User specified decomposition for the field solver in z

    """

    __w3dinputs__ = ['iondensity','electrontemperature','plasmapotential',
                     'electrondensitymaxscale']
    __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                     'mgmaxiters','mgtol','mgmaxlevels','mgform',
                     'mgverbose','mgntverbose',
                     'lcndbndy','icndbndy','laddconductor','lprecalccoeffs']

    def __init__(self,lreducedpickle=1,**kw):
        kw['lreducedpickle'] = lreducedpickle
        self.grid_overlap = 2
        SubcycledPoissonSolver.__init__(self,kwdict=kw)
        self.solvergeom = w3d.XYZgeom
        self.ncomponents = 1

        # --- Make sure that the bounds have acceptable values.
        assert 0 <= min(self.bounds) and max(self.bounds) <= 2,"The boundary conditions have an incorrect value. They must be one of dirichlet, neumann or periodic."

        # --- Kludge - make sure that the multigrid3df routines never sets up
        # --- any conductors.
        f3d.gridmode = 1

        # --- Save input parameters
        self.processdefaultsfrompackage(MultiGrid3D.__w3dinputs__,w3d,kw)
        self.processdefaultsfrompackage(MultiGrid3D.__f3dinputs__,f3d,kw)

        # --- If there are any remaning keyword arguments, raise an error.
        assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

        self.initializeconductors()

        # --- Give these variables dummy initial values.
        self.mgiters = 0
        self.mgerror = 0.

        # --- At the start, assume that there are no bends. This is corrected
        # --- in the solve method when there are bends.
        self.linbend = false

        # --- By default, the E field is not directly calculated.
        self.lwithselfe = 0
        self.lwithselfep = 0

        # --- This allows some introspection.
        self.dict_of_grids = {
            'phi':{'getter':'getphi', 'centering':'node', 'units':'V'},
            'rho':{'getter':'getrho', 'centering':'node', 'units':'C/m**3'},
            'Ex':{'getter':'getex', 'centering':'node', 'units':'V/m'},
            'Ey':{'getter':'getey', 'centering':'node', 'units':'V/m'},
            'Ez':{'getter':'getez', 'centering':'node', 'units':'V/m'},
        }

    def initializeconductors(self):
        # --- Create the attributes for holding information about conductors
        # --- and conductor objects.
        # --- Note that a conductor object will be created for each value of
        # --- fselfb. This is needed since fselfb effects how the coarsening
        # --- is done, and different conductor data sets are needed for
        # --- different coarsenings.

        # --- This stores the ConductorType objects. Note that the objects are
        # --- not actually created until getconductorobject is called.
        self.conductorobjects = {}

        # --- This stores the conductors that have been installed in each
        # --- of the conductor objects.
        self.installedconductorlists = {}

        # --- This is a list of conductors that have been added.
        # --- New conductors are not actually installed until the data is needed,
        # --- when getconductorobject is called.
        # --- Each element of this list contains all of the input to the
        # --- installconductor method.
        self.conductordatalist = []

    def __getstate__(self):
        dict = SubcycledPoissonSolver.__getstate__(self)
        if self.lreducedpickle:

            # --- Write out an empy conductorobjects since it can be big. Also,
            # --- write out an empty list of conductors so they will all be
            # --- reinstalled upon restoration.
            dict['conductorobjects'] = {}
            dict['installedconductorlists'] = {}

            if '_rho' in dict: del dict['_rho']
            if '_phi' in dict: del dict['_phi']
            if 'bfieldp' in dict: del dict['bfieldp']

        return dict

    def __setstate__(self,dict):
        SubcycledPoissonSolver.__setstate__(self,dict)

        # --- Check if an old file is being restored
        if 'conductorobjects' not in self.__dict__:

            # --- Create the appropriate attributes that are now needed.
            # --- This is not the best thing, since is replicates code in
            # --- the __init__
            self.conductorobjects = {}
            self.installedconductorlists = {}
            self.conductordatalist = []

            # --- Get the list of conductors from old formats
            if 'newconductorlist' in self.__dict__:
                conductorlist = self.newconductorlist
                del self.newconductorlist
            elif 'conductorlist' in self.__dict__:
                conductorlist = self.conductorlist
                del self.conductorlist
            else:
                conductorlist = []
            if 'lprecalccoeffs' not in self.__dict__:
                self.lprecalccoeffs = 0

            for conductor in conductorlist:
                self.installconductor(conductor)

    def getconductorobject(self,fselfb=0.):
        "Checks for and installs any conductors not yet installed before returning the object"
        # --- This is the routine that does the creation of the ConductorType
        # --- objects if needed and ensures that all conductors are installed
        # --- into it.

        # --- This method is needed during a restore from a pickle, since this
        # --- object may be restored before the conductors. This delays the
        # --- installation of the conductors until they are really needed.

        # --- There is a special case, fselfb='p', which refers to the conductor
        # --- object that has the data generated relative to the particle domain,
        # --- which can be different from the field domain, especially in parallel.
        if fselfb == 'p':
            # --- In serial, just use a reference to the conductor object for the
            # --- first iselfb group.
            if not lparallel and 'p' not in self.conductorobjects:
                self.conductorobjects['p'] = self.conductorobjects[top.fselfb[0]]
                self.installedconductorlists['p'] = self.installedconductorlists[top.fselfb[0]]
            # --- In parallel, a whole new instance is created (using the
            # --- setdefaults below).
            # --- Check to make sure that the grid the conductor uses is consistent
            # --- with the particle grid. This is needed so that the conductor
            # --- data is updated when particle load balancing is done. If the
            # --- data is not consistent, delete the conductor object so that
            # --- everything is reinstalled.
            try:
                conductorobject = self.conductorobjects['p']
                ixproc = self.ppdecomp.ixproc
                iyproc = self.ppdecomp.iyproc
                izproc = self.ppdecomp.izproc
                if (conductorobject.levelix[0] != self.ppdecomp.ix[ixproc] or
                    conductorobject.levelnx[0] != self.ppdecomp.nx[ixproc] or
                    conductorobject.leveliy[0] != self.ppdecomp.iy[iyproc] or
                    conductorobject.levelny[0] != self.ppdecomp.ny[iyproc] or
                    conductorobject.leveliz[0] != self.ppdecomp.iz[izproc] or
                    conductorobject.levelnz[0] != self.ppdecomp.nz[izproc]):
                    del self.conductorobjects['p']
                    del self.installedconductorlists['p']
            except KeyError:
                # --- 'p' object has not yet been created anyway, so do nothing.
                pass

        #conductorobject = self.conductorobjects.setdefault(fselfb,ConductorType())
        try:
            conductorobject = self.conductorobjects[fselfb]
        except KeyError:
            conductorobject = ConductorType()
            self.conductorobjects[fselfb] = conductorobject

        installedconductorlist = self.installedconductorlists.setdefault(fselfb,[])

        # --- Now, make sure that the conductors are installed into the object.
        # --- This may be somewhat inefficient, since it loops over all of the
        # --- conductors everytime. This makes the code more robust, though, since
        # --- it ensures that all conductors will be properly installed into
        # --- the conductor object.
        for conductordata in self.conductordatalist:
            self._installconductor(conductorobject,installedconductorlist,
                                   conductordata,fselfb)

        # --- Return the desired conductor object
        return conductorobject

    def setconductorvoltage(self,voltage,condid=0,discrete=false,
                            setvinject=false):
        'calls setconductorvoltage'
        # --- Loop over all of the selfb groups to that all conductor objects
        # --- are handled.
        for iselfb in range(top.nsselfb):
            conductorobject = self.getconductorobject(top.fselfb[iselfb])
            setconductorvoltage(voltage,condid,discrete,setvinject,
                                conductors=conductorobject)

    def getpdims(self):
        # --- Returns the dimensions of the arrays used by the particles

        # --- If there are any relativistic groups, then turn on the code
        # --- which uses the selfe array.
        if top.allocated('fselfb') and abs(top.fselfb).max() > 0.:
            # --- This is probably redundant, but it shouldn't hurt.
            # --- This forces all species to use the precalculated E field
            # --- if any have the B correction.
            top.efetch = 3

        try:
            self.lwithselfep
        except AttributeError:
            self.lwithselfep = 0
        self.lwithselfep = (self.lwithselfep or
                            sometrue(top.efetch == 3) or
                            maxnd(top.depos_order) > 1)

        if self.lwithselfep:
            return ((1+self.nxp+2*self.nxguardrho,
                     1+self.nyp+2*self.nyguardrho,
                     1+self.nzp+2*self.nzguardrho),
                    (3,1+self.nxp+2*self.nxguarde,
                       1+self.nyp+2*self.nyguarde,
                       1+self.nzp+2*self.nzguarde),
                    (1+self.nxp+2*self.nxguardphi,
                     1+self.nyp+2*self.nyguardphi,
                     1+self.nzp+2*self.nzguardphi))
        else:
            return ((1+self.nxp+2*self.nxguardrho,
                     1+self.nyp+2*self.nyguardrho,
                     1+self.nzp+2*self.nzguardrho),
                    (1+self.nxp+2*self.nxguardphi,
                     1+self.nyp+2*self.nyguardphi,
                     1+self.nzp+2*self.nzguardphi))

    def getdims(self):
        # --- Returns the dimensions of the arrays used by the field solver
        dims = [(1+self.nxlocal+2*self.nxguardrho,
                 1+self.nylocal+2*self.nyguardrho,
                 1+self.nzlocal+2*self.nzguardrho),
                (1+self.nxlocal+2*self.nxguardphi,
                 1+self.nylocal+2*self.nyguardphi,
                 1+self.nzlocal+2*self.nzguardphi)]
        try:
            self.lwithselfe
        except AttributeError:
            self.lwithselfe = 0
        if self.lwithselfe:
            dims[1:1] = [(3,1+self.nxlocal+2*self.nxguarde,
                            1+self.nylocal+2*self.nyguarde,
                            1+self.nzlocal+2*self.nzguarde)]
        return tuple(dims)

    def getrho(self):
        'Returns the rho array without the guard cells'
        return self.source[self.nxguardrho:-self.nxguardrho or None,
                           self.nyguardrho:-self.nyguardrho or None,
                           self.nzguardrho:-self.nzguardrho or None]

    def getrhop(self):
        'Returns the rhop array without the guard cells'
        return self.sourcep

    def getphi(self):
        'Returns the phi array without the guard cells'
        return self.potential[self.nxguardphi:-self.nxguardphi or None,
                              self.nyguardphi:-self.nyguardphi or None,
                              self.nzguardphi:-self.nzguardphi or None]

    def getphip(self):
        'Returns the phip array without the guard cells'
        return self.potentialp[self.nxguardphi:-self.nxguardphi or None,
                               self.nyguardphi:-self.nyguardphi or None,
                               self.nzguardphi:-self.nzguardphi or None]

    def getselfe(self,recalculate=None,lzero=true):
        'Returns the E field array without the guard cells'
        self.calcselfe(recalculate=recalculate,lzero=lzero)
        return self.field[:,self.nxguarde:-self.nxguarde or None,
                            self.nyguarde:-self.nyguarde or None,
                            self.nzguarde:-self.nzguarde or None]

    def getselfep(self,recalculate=None,lzero=true):
        'Returns the E fieldp array without the guard cells'
        self.calcselfep(recalculate=recalculate,lzero=lzero)
        return self.fieldp[:,self.nxguarde:-self.nxguarde or None,
                             self.nyguarde:-self.nyguarde or None,
                             self.nzguarde:-self.nzguarde or None]

    def getex(self,recalculate=None,lzero=true):
        return self.getselfe(recalculate,lzero)[0,...]

    def getey(self,recalculate=None,lzero=true):
        return self.getselfe(recalculate,lzero)[1,...]

    def getez(self,recalculate=None,lzero=true):
        return self.getselfe(recalculate,lzero)[2,...]

    def _setuprhoproperty():
        doc = "Charge density array"
        def fget(self):
            return self.returnsource(0,0)
        def fset(self,value):
            self.returnsource(0,0)[...] = value
        return locals()
    rho = property(**_setuprhoproperty())
    del _setuprhoproperty

    def _setupphiproperty():
        doc = "Electrostatic potential array, including guard cells"
        def fget(self):
            return self.returnpotential(0,0)
        def fset(self,value):
            self.returnpotential(0,0)[...] = value
        return locals()
    phi = property(**_setupphiproperty())
    del _setupphiproperty

    def _setupselfeproperty():
        doc = "Electric field array for particles"
        def fget(self):
            self.calcselfe()
            return self.returnfield(0,0)
        def fset(self,value):
            self.returnfieldp(0,0)[...] = value
        return locals()
    selfe = property(**_setupselfeproperty())
    del _setupselfeproperty

    def _setuprhopproperty():
        doc = "Charge density array for particles"
        def fget(self):
            return self.returnsourcep(0,0,0)
        def fset(self,value):
            self.returnsourcep(0,0,0)[...] = value
        return locals()
    rhop = property(**_setuprhopproperty())
    del _setuprhopproperty

    def _setupphipproperty():
        doc = "Electrostatic potential array for particles"
        def fget(self):
            return self.returnpotentialp(0,0)
        def fset(self,value):
            self.returnpotentialp(0,0)[...] = value
        return locals()
    phip = property(**_setupphipproperty())
    del _setupphipproperty

    def _setupselfepproperty():
        doc = "Electric field array for particles"
        def fget(self):
            self.calcselfep()
            return self.returnfieldp(0,0)
        def fset(self,value):
            self.returnfieldp(0,0)[...] = value
        return locals()
    selfep = property(**_setupselfepproperty())
    del _setupselfepproperty

    def loadrho(self,lzero=None,lfinalize_rho=None,pgroups=None,**kw):
        SubcycledPoissonSolver.loadsource(self,lzero,lfinalize_rho,pgroups,**kw)

    def fetche(self,*args,**kw):
        SubcycledPoissonSolver.fetchfield(self,*args,**kw)

    def fetchphi(self,*args,**kw):
        SubcycledPoissonSolver.fetchpotential(self,*args,**kw)

    def gtlchg(self):
        'Calculate the line charge, putting it into the array top.linechg'
        if self.solvergeom == w3d.XYZgeom:
            gtlchg3dfromrho(self.nxlocal,self.nylocal,self.nzlocal,
                            self.nxguardrho,self.nyguardrho,self.nzguardrho,
                            self.rho,
                            self.dx,self.dy,self.dz,
                            self.getzgrid(),self.zmminlocal,
                            self.l2symtry,self.l4symtry,
                            self.izproc==self.nzprocs-1)
        elif self.solvergeom==w3d.RZgeom:
            gtlchgrzfromrho(self.nxlocal,self.nzlocal,
                            self.nxguardrho,self.nzguardrho,
                            self.rho,
                            self.dx,self.dz,
                            self.getzgrid(),self.zmminlocal,
                            self.izproc==self.nzprocs-1)

    def getese(self):
        'Calculate the electrostatic potential energy, rho*phi, and put it in top.ese'
        # --- ese must be an array to get the returned value in the calls below
        ese = zeros(1,'d')

        if self.solvergeom == w3d.XYZgeom:
            getese3dfromrhophi(self.nxlocal,self.nylocal,self.nzlocal,
                               self.nxguardphi,self.nyguardphi,self.nzguardphi,
                               self.nxguardrho,self.nyguardrho,self.nzguardrho,
                               self.rho,self.phi,
                               self.dx,self.dy,self.dz,self.l4symtry,self.l2symtry,
                               ese)
        elif self.solvergeom==w3d.RZgeom:
            geteserzfromrhophi(self.nxlocal,self.nzlocal,
                               self.nxguardphi,self.nzguardphi,
                               self.nxguardrho,self.nzguardrho,
                               self.rho,self.phi,
                               self.dx,self.dz,self.xmminlocal,
                               ese)
        top.ese = ese[0]

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
        w  = pgroup.sw[js]*pgroup.dtscale[js]
        if top.wpid==0:
            wfact = zeros((0,), 'd')
        else:
            wfact = pgroup.pid[i:i+n,top.wpid-1]
        depos_order = top.depos_order[:,js]
        self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,w,depos_order)

    def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,w,
                             depos_order):
        n = len(x)
        if n == 0: return
        if isinstance(self.sourcep,float): return
        if top.wpid==0:
            setrho3d(self.sourcep,n,x,y,z,zgrid,q,w,top.depos,depos_order,
                     self.nxp,self.nyp,self.nzp,
                     self.nxguardrho,self.nyguardrho,self.nzguardrho,
                     self.dx,self.dy,self.dz,
                     self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                     self.solvergeom==w3d.RZgeom)
        else:
            setrho3dw(self.sourcep,n,x,y,z,zgrid,wfact,q,w,top.depos,depos_order,
                      self.nxp,self.nyp,self.nzp,
                      self.nxguardrho,self.nyguardrho,self.nzguardrho,
                      self.dx,self.dy,self.dz,
                      self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                      self.solvergeom==w3d.RZgeom)

    def returnbfieldp(self):
        try:
            self.bfieldp
        except AttributeError:
            self.bfieldp = None
        if self.bfieldp is None or self.bfieldp.shape != self.fieldp.shape:
            self.bfieldp = zeros_like(self.fieldp)
        return self.bfieldp

    def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
        # --- This is called by fetchfield from fieldsolver.py
        # --- Only sets the E field from the potential
        n = len(x)
        if n == 0: return
        if top.efetch[js] == 0: return
        if top.efetch[js] == 3 and isinstance(self.fieldp,float): return
        if top.efetch[js] != 3 and isinstance(self.potentialp,float): return
        if not f3d.lcorrectede:
            sete3d(self.potentialp,self.fieldp,n,x,y,z,self.getzgridprv(),
                   self.xmminp,self.ymminp,self.zmminp,
                   self.dx,self.dy,self.dz,self.nxp,self.nyp,self.nzp,
                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                   self.nxguarde,self.nyguarde,self.nzguarde,
                   top.efetch[js],top.depos_order[:,js],
                   ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)
        else:
            sete3dwithconductor(self.getconductorobject('p'),n,x,y,z,
                                top.efetch[js],ex,ey,ez,
                                self.getzgridprv(),
                                self.xmminp,self.ymminp,self.zmminp,
                                self.dx,self.dy,self.dz,self.nxp,self.nyp,self.nzp,
                                self.nxguardphi,self.nyguardphi,self.nzguardphi,
                                self.nxguarde,self.nyguarde,self.nzguarde,
                                self.potentialp,self.fieldp,
                                self.l2symtry,self.l4symtry,
                                self.solvergeom==w3d.RZgeom)
        if top.allocated('fselfb') and abs(top.fselfb).max() > 0.:
            # --- Note that now fetche3dfrompositions takes B field arguments so
            # --- this code should work OK now. The if statement is left just case
            # --- there is an odd situation.
            if len(bx) != n: return
            bfieldp = self.returnbfieldp()
            setb3d(bfieldp,n,x,y,z,self.getzgridprv(),bx,by,bz,
                   self.nxp,self.nyp,self.nzp,
                   self.nxguarde,self.nyguarde,self.nzguarde,
                   self.dx,self.dy,self.dz,
                   self.xmminp,self.ymminp,self.zmminp,
                   self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)

    def fetchpotentialfrompositions(self,x,y,z,phi):
        n = len(x)
        if n == 0: return
        if isinstance(self.potentialp,float): return
        nxp = self.nxp + 2*self.nxguardphi
        nyp = self.nyp + 2*self.nyguardphi
        nzp = self.nzp + 2*self.nzguardphi
        xmminp = self.xmminp - self.dx*self.nxguardphi
        xmmaxp = self.xmmaxp + self.dx*self.nxguardphi
        ymminp = self.ymminp - self.dy*self.nyguardphi
        ymmaxp = self.ymmaxp + self.dy*self.nyguardphi
        zmminp = self.zmminp - self.dz*self.nzguardphi + self.getzgridprv()
        zmmaxp = self.zmmaxp + self.dz*self.nzguardphi + self.getzgridprv()
        getgrid3d(n,x,y,z,phi,nxp,nyp,nzp,self.potentialp,
                  xmminp,xmmaxp,ymminp,ymmaxp,zmminp,zmmaxp,
                  self.l2symtry,self.l4symtry)

    def setsourceforfieldsolve(self,*args):
        SubcycledPoissonSolver.setsourceforfieldsolve(self,*args)
        if self.lparallel:
            SubcycledPoissonSolver.setsourcepforparticles(self,*args)
            if isinstance(self.source,float): return
            if isinstance(self.sourcep,float): return
            setrhoforfieldsolve3d(self.nxlocal,self.nylocal,self.nzlocal,self.source,
                                  self.nxp,self.nyp,self.nzp,self.sourcep,
                                  self.nxguardrho,self.nyguardrho,self.nzguardrho,
                                  self.fsdecomp,self.ppdecomp)

    def getpotentialpforparticles(self,*args):
        self.setpotentialpforparticles(*args)
        if not self.lparallel:
            SubcycledPoissonSolver.getpotentialpforparticles(self,*args)
        else:
            if isinstance(self.potential,float): return
            if isinstance(self.potentialp,float): return
            getphipforparticles3d(1,self.nxlocal,self.nylocal,self.nzlocal,
                                  self.nxguardphi,self.nyguardphi,self.nzguardphi,
                                  self.potential,
                                  self.nxp,self.nyp,self.nzp,self.potentialp,
                                  self.fsdecomp,self.ppdecomp)

        iselfb = args[2]
        if iselfb == 0 and f3d.lcorrectede:
            # --- This only needs to be calculated once, so is only done
            # --- when iselfb == 0.
            conductorobjectp = self.getconductorobject('p')
            # --- This sets up the icgrid
            setupconductorfielddata(self.nx,self.ny,self.nz,
                                    self.nxp,self.nyp,self.nzp,
                                    self.dx,self.dy,self.dz,conductorobjectp,
                                    self.ppdecomp)
            # --- This calculates the field
            getefieldatconductorsubgrid(conductorobjectp,
                                   self.dx,self.dy,self.dz,
                                   self.nxp,self.nyp,self.nzp,
                                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                                   self.potentialp,self.bounds)
        if sometrue(top.efetch == 3) or maxnd(top.depos_order) > 1:
            self.setfieldpforparticles(*args)
            indts = args[1]
            # --- If this is the first group, set make sure that fieldp gets
            # --- zeroed out. Otherwise, the data in fieldp is accumulated.
            # --- This coding relies on the fact that fieldsolver does the
            # --- loops in descending order.
            tmpnsndts = getnsndtsforsubcycling()
            lzero = ((indts == tmpnsndts-1) and (iselfb == top.nsselfb-1))
            #if lzero:
            #  tfieldp = transpose(self.fieldp)
            #  tfieldp[...] = 0.
            self.calcselfep(recalculate=1,lzero=lzero)
            if top.allocated('fselfb') and abs(top.fselfb[iselfb]) > 0:
                # --- If the self-B correction is nonzero, then calculate and include
                # --- the approximate correction terms A and dA/dt.
                bfieldp = self.returnbfieldp()
                bfieldp[...] = 0.
                self.getselfb(bfieldp,top.fselfb[iselfb],self.potentialp)
                self.adddadttoe(self.fieldp,top.fselfb[iselfb],self.potentialp)

            if iselfb == 0 and f3d.lcorrectede:
                # --- Now correct the E field at conductor points. This matters when
                # --- the edge of conductors are aligned with the mesh and there are no
                # --- subgrid points there.
                # --- This is done when iselfb == 0 since that will be the last
                # --- species - this routine modifies fieldp in place.
                fixefieldatconductorpoints(conductorobjectp,
                                           self.dx,self.dy,self.dz,
                                           self.nxp,self.nyp,self.nzp,
                                           self.nxguarde,self.nyguarde,self.nzguarde,
                                           self.fieldp)

    def applysourceboundaryconditions(self):
        applyrhoboundaryconditions3d(self.ncomponents,
                                 self.nxlocal,self.nylocal,self.nzlocal,
                                 self.nxguardrho,self.nyguardrho,self.nzguardrho,
                                 self.source,self.bounds,self.fsdecomp,
                                 self.solvergeom==w3d.RZgeom)

    def calcselfe(self,recalculate=None,lzero=true):
        if not self.lparallel:
            # --- If serial, then defer to calcselfep since field and fieldp would
            # --- be the same array.
            self.calcselfep(recalculate=recalculate,lzero=lzero)
            self.field = self.fieldp

        else:
            # --- The rest is the same as self.calcselfep, but using the non-p
            # --- attributes.

            # --- Since the selfe is never directly used, except for diagnostics,
            # --- it should always be recalculated by default.
            if recalculate is None: recalculate = 1

            self.lwithselfe = 1
            self.allocatedataarrays()
            self.setfieldforfieldsolve(0,0,0)
            if recalculate:
                if isinstance(self.potential,float): return
                if isinstance(self.field,float): return
                getselfe3d(self.potential,self.nxlocal,self.nylocal,self.nzlocal,
                           self.nxguardphi,self.nyguardphi,self.nzguardphi,
                           self.field,self.nxguarde,self.nyguarde,self.nzguarde,
                           self.dx,self.dy,self.dz,lzero)

    def calcselfep(self,recalculate=None,lzero=true):
        # --- Check if the E field should be recalculated.
        # --- If it had not yet been calculated at all, then definitely
        # --- calculate it now.
        if not self.lwithselfep: recalculate = 1

        # --- If the E field is not actively being used, then recalculate it,
        # --- unless recalculate is passed in by the user.
        if (alltrue(top.efetch != 3) and maxnd(top.depos_order) == 1 and
            recalculate is None):
            recalculate = 1

        self.lwithselfep = 1
        self.allocatedataarrays()
        self.setfieldpforparticles(0,0,0)
        if recalculate:
            if isinstance(self.potentialp,float): return
            if isinstance(self.fieldp,float): return
            getselfe3d(self.potentialp,self.nxp,self.nyp,self.nzp,
                       self.nxguardphi,self.nyguardphi,self.nzguardphi,
                       self.fieldp,self.nxguarde,self.nyguarde,self.nzguarde,
                       self.dx,self.dy,self.dz,lzero)

    def getslicewithguard(self,i1,i2,guard):
        if i1 is not None: i1 = i1 + guard
        if i2 is not None:
            if i2 < 0: i2 = i2 - guard
            else:      i2 = i2 + guard
        if i1 is None and guard > 0: i1 = +guard
        if i2 is None and guard > 0: i2 = -guard
        return slice(i1,i2)

    def getselfb(self,bfieldp,fselfb,potentialp):
        ix = self.getslicewithguard(None,None,self.nxguardphi)
        iy = self.getslicewithguard(None,None,self.nyguardphi)
        iz = self.getslicewithguard(None,None,self.nzguardphi)
        Az = (fselfb/clight**2)*potentialp
        if self.ny > 0:
            ysu = self.nyguardphi+1
            yeu = -self.nyguardphi+1
            if yeu==0:yeu=None
            ysl = self.nyguardphi-1
            yel = -self.nyguardphi-1
            bfieldp[0,:,:,:] += (Az[ix,ysu:yeu,iz] - Az[ix,ysl:yel,iz])/(2.*self.dy)
        if self.nx > 0:
            xsu = self.nxguardphi+1
            xeu = -self.nxguardphi+1
            if xeu==0:xeu=None
            xsl = self.nxguardphi-1
            xel = -self.nxguardphi-1
            bfieldp[1,:,:,:] -= (Az[xsu:xeu,iy,iz] - Az[xsl:xel,iy,iz])/(2.*self.dx)

    def adddadttoe(self,fieldp,fselfb,potentialp):
        """Ez = -dA/dt = -beta**2 dphi/dz"""
        ix = self.getslicewithguard(None,None,self.nxguardphi)
        iy = self.getslicewithguard(None,None,self.nyguardphi)
        # --- This assumes that nzguard is always 1
        Ez = (fselfb/clight)**2*(potentialp[ix,iy,2:]-potentialp[ix,iy,:-2])/(2.*self.dz)
        fieldp[2,:,:,:] += Ez

    def installconductor(self,conductor,
                              xmin=None,xmax=None,
                              ymin=None,ymax=None,
                              zmin=None,zmax=None,
                              dfill=None):
        'Install the conductor into the field solver'
        # --- This only adds the conductor to the list. The data is only actually
        # --- installed when it is needed, during a call to getconductorobject.
        self.conductordatalist.append((conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill))

    def uninstallconductor(self,conductor):
        """Uninstall the conductor from the field solver. Note that the conductor
    data will be cleared, and will be regenerated upon the next field solve"""
        # --- Loop over the list of conductors, removing the desired conductor.
        # --- Note that if it was installed multiple times, then every instance
        # --- will be uninstalled. Clear out the conductor data, forcing it to be
        # --- generated upon the next field solve.
        conductorfound = False
        i = 0
        while i < len(self.conductordatalist):
            if self.conductordatalist[i][0] is conductor:
                conductorfound = True
                del self.conductordatalist[i]
                self.clearconductors()
            else:
                i += 1

        if not conductorfound:
            raise Exception('conductor was not found')

    def _installconductor(self,conductorobject,installedlist,conductordata,
                          fselfb):
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
                          installrz=0,
                          solvergeom=self.solvergeom,conductors=conductorobject,
                          mgmaxlevels=mgmaxlevels,decomp=decomp)

    def hasconductors(self):
        return len(self.conductordatalist) > 0

    def clearconductors(self,fselfblist=None):
        "Clear out the conductor data"
        if fselfblist is None:
            fselfblist = self.conductorobjects.keys()
        for fselfb in fselfblist:
            if fselfb in self.conductorobjects:
                conductorobject = self.conductorobjects[fselfb]
                conductorobject.interior.n = 0
                conductorobject.evensubgrid.n = 0
                conductorobject.oddsubgrid.n = 0
                self.installedconductorlists[fselfb] = []

    def cond_potmg(self,mglevel=0,mgform=None,iselfb=0,phi=None):
        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
        if phi is None:
            phi = self.potential
        if mgform is None:
            mgform = self.mgform
        cond_potmg(conductorobject.interior,self.nxlocal,self.nylocal,self.nzlocal,
                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                   phi,mglevel,self.mgform,false)

    def find_mgparam(self,lsavephi=false,resetpasses=1):
        # --- This is a temporary kludge, the same as is done in genericpf
        self._phi = self.potential
        find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
                     solver=self,pkg3d=self)

    def getmgverbose(self):
        # --- Check if the convergence diagnostic should be printed.
        if self.mgntverbose < 0:
            mgverbose = 0
        elif self.mgntverbose > 1:
            if (top.it%self.mgntverbose) == 0:
                mgverbose = 1
            else:
                mgverbose = 0
        else:
            mgverbose = self.mgverbose
        return mgverbose

    def dosolve(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        if not self.l_internal_dosolve: return
        # --- set for longitudinal relativistic contraction
        if zfact is None:
            beta = top.pgroup.fselfb[iselfb]/clight
            zfact = 1./sqrt((1.-beta)*(1.+beta))
        else:
          beta = sqrt( (1.-1./zfact)*(1.+1./zfact) )

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source
        if isinstance(self.potential,float): return

        # --- Setup data for bends.
        rstar = zeros(3+self.nzlocal,'d')
        if top.bends:

            # --- This commented out code does the same thing as the line below
            # --- setting linbend but is a bit more complicated. It is preserved
            # --- in case of some unforeseen problem with the code below.
            #ii = (top.cbendzs <= self.zmmax+zgrid and
            #                     self.zmmin+zgrid <= top.cbendze)
            #self.linbend = sometrue(ii)

            setrstar(rstar,self.nzlocal,self.dz,self.zmminlocal,self.getzgrid())
            self.linbend = rstar.min() < largepos

        mgverbose = self.getmgverbose()
        mgiters = zeros(1,'l')
        mgerror = zeros(1,'d')
        # --- This takes care of clear out the conductor information if needed.
        # --- Note that f3d.gridmode is passed in below - this still allows the
        # --- user to use the addconductor method if needed.
        if self.gridmode == 0: self.clearconductors([top.pgroup.fselfb[iselfb]])
        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
        if self.electrontemperature == 0:
            multigrid3dsolve(iwhich,self.nx,self.ny,self.nz,
                             self.nxlocal,self.nylocal,self.nzlocal,
                             self.nxguardphi,self.nyguardphi,self.nzguardphi,
                             self.nxguardrho,self.nyguardrho,self.nzguardrho,
                             self.dx,self.dy,self.dz*zfact,self.potential,self.source,
                             rstar,self.linbend,self.bounds,
                             self.xmmin,self.ymmin,self.zmmin*zfact,
                             self.mgparam,self.mgform,mgiters,self.mgmaxiters,
                             self.mgmaxlevels,mgerror,self.mgtol,mgverbose,
                             self.downpasses,self.uppasses,
                             self.lcndbndy,self.laddconductor,self.icndbndy,
                             f3d.gridmode,conductorobject,self.lprecalccoeffs,
                             self.fsdecomp)
        else:
            self.iondensitygrid3d = Grid3dtype()
            setupiondensitygrid3d(self.xmmin,self.ymmin,self.zmmin,
                                  self.dx,self.dy,self.dz,
                                  self.nxlocal,self.nylocal,self.nzlocal,
                                  self.nxguardrho,self.nyguardrho,self.nzguardrho,
                                  self._rho,self.iondensitygrid3d)
            multigridbe3dsolve(iwhich,self.nx,self.ny,self.nz,
                               self.nxguardphi,self.nyguardphi,self.nzguardphi,
                               self.nxguardrho,self.nyguardrho,self.nzguardrho,
                               self.dx,self.dy,self.dz*zfact,
                               self.potential,self.source,
                               rstar,self.linbend,self.bounds,
                               self.xmmin,self.ymmin,self.zmmin*zfact,
                               self.mgparam,mgiters,self.mgmaxiters,
                               self.mgmaxlevels,mgerror,self.mgtol,mgverbose,
                               self.downpasses,self.uppasses,
                               self.lcndbndy,self.laddconductor,self.icndbndy,
                               f3d.gridmode,conductorobject,
                               self.iondensitygrid3d,
                               self.fsdecomp)
        self.mgiters = mgiters[0]
        self.mgerror = mgerror[0]

    ##########################################################################
    # Define the basic plot commands
    def genericpf(self,kw,pffunc):
        fselfb = kw.get('fselfb',(top.allocated('fselfb') and top.fselfb[0]) or 0.)
        if 'fselfb' in kw: del kw['fselfb']
        kw['conductors'] = self.getconductorobject(fselfb)
        kw['solver'] = self
        pffunc(**kw)
    def pfxy(self,**kw): self.genericpf(kw,pfxy)
    def pfzx(self,**kw): self.genericpf(kw,pfzx)
    def pfzy(self,**kw): self.genericpf(kw,pfzy)
    def pfxyg(self,**kw): self.genericpf(kw,pfxyg)
    def pfzxg(self,**kw): self.genericpf(kw,pfzxg)
    def pfzyg(self,**kw): self.genericpf(kw,pfzyg)
    def pcrhozx(self,*args,**kw): pcrhozx(*args,solver=self,**kw)
    def pcrhozy(self,*args,**kw): pcrhozy(*args,solver=self,**kw)
    def pcrhoxy(self,*args,**kw): pcrhoxy(*args,solver=self,**kw)
    def pcphizx(self,*args,**kw): pcphizx(*args,solver=self,**kw)
    def pcphizy(self,*args,**kw): pcphizy(*args,solver=self,**kw)
    def pcphixy(self,*args,**kw): pcphixy(*args,solver=self,**kw)
    def pcselfezx(self,*args,**kw): pcselfezx(*args,solver=self,**kw)
    def pcselfezy(self,*args,**kw): pcselfezy(*args,solver=self,**kw)
    def pcselfexy(self,*args,**kw): pcselfexy(*args,solver=self,**kw)

    def getresidual(self):
        res = zeros(shape(self._phi),'d')
        dxsqi  = 1./self.dx**2
        dysqi  = 1./self.dy**2
        dzsqi  = 1./self.dz**2
        reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
        rho = self._rho*reps0c
        conductorobject = self.getconductorobject(top.pgroup.fselfb[0])
        residual3d(self.nxlocal,self.nylocal,self.nzlocal,
                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                   self.nxguardrho,self.nyguardrho,self.nzguardrho,
                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                   dxsqi,dysqi,dzsqi,self._phi,rho,res,
                   0,self.bounds,self.mgparam,self.mgform,true,
                   self.lcndbndy,self.icndbndy,conductorobject,self.lprecalccoeffs)
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
        dysqi  = 1./self.dy**2
        dzsqi  = 1./self.dz**2
        reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
        rho = self._rho*reps0c
        result = zeros(shape(self._phi),'d')
        residual3d(self.nxlocal,self.nylocal,self.nzlocal,
                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                   self.nxguardrho,self.nyguardrho,self.nzguardrho,
                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                   dxsqi,dysqi,dzsqi,self._phi,rho,result,
                   0,bounds,self.mgparam,self.mgform,true,
                   self.lcndbndy,self.icndbndy,conductorobject,self.lprecalccoeffs)

        conductorobject.interior.istart[:] = istartsave
        if includeboundaries:
            # --- Undo the applyboundaryconditions3d from above.
            applyboundaryconditions3d(self.nxlocal,self.nylocal,self.nzlocal,
                                      self.nxguardphi,self.nyguardphi,self.nzguardphi,
                                      self._phi,1,self.bounds,true,false)

        # --- Remove the premultiplying factor
        result /= reps0c
        return result


MultiGrid = MultiGrid3D

##############################################################################
##############################################################################
class FullMultiGrid3D(MultiGrid3D):
    def __init__(self,mgmaxvcycles=2,**kw):
        self.mgmaxvcycles = mgmaxvcycles
        MultiGrid3D.__init__(self,**kw)

    def dosolve(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        if not self.l_internal_dosolve: return
        # --- set for longitudinal relativistic contraction
        if zfact is None:
            beta = top.pgroup.fselfb[iselfb]/clight
            zfact = 1./sqrt((1.-beta)*(1.+beta))
        else:
          beta = sqrt( (1.-1./zfact)*(1.+1./zfact) )

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source
        if isinstance(self.potential,float): return

        # --- Setup data for bends.
        rstar = zeros(3+self.nzlocal,'d')
        if top.bends:

            # --- This commented out code does the same thing as the line below
            # --- setting linbend but is a bit more complicated. It is preserved
            # --- in case of some unforeseen problem with the code below.
            #ii = (top.cbendzs <= self.zmmax+zgrid and
            #                     self.zmmin+zgrid <= top.cbendze)
            #self.linbend = sometrue(ii)

            setrstar(rstar,self.nzlocal,self.dz,self.zmminlocal,self.getzgrid())
            self.linbend = rstar.min() < largepos

        mgverbose = self.getmgverbose()
        mgiters = zeros(1,'l')
        mgerror = zeros(1,'d')
        # --- This takes care of clear out the conductor information if needed.
        # --- Note that f3d.gridmode is passed in below - this still allows the
        # --- user to use the addconductor method if needed.
        if self.gridmode == 0: self.clearconductors([top.pgroup.fselfb[iselfb]])
        conductorobject = self.getconductorobject(top.pgroup.fselfb[iselfb])
        if self.electrontemperature == 0:
            fullmultigrid3dsolve(iwhich,0,self.nx,self.ny,self.nz,
                             self.nxlocal,self.nylocal,self.nzlocal,
                             self.nxguardphi,self.nyguardphi,self.nzguardphi,
                             self.nxguardrho,self.nyguardrho,self.nzguardrho,
                             self.dx,self.dy,self.dz*zfact,self.potential,self.source,
                             rstar,self.linbend,self.bounds,
                             self.xmmin,self.ymmin,self.zmmin*zfact,
                             self.mgparam,mgiters,self.mgmaxiters,self.mgmaxvcycles,
                             self.mgmaxlevels,mgerror,self.mgtol,mgverbose,
                             self.downpasses,self.uppasses,
                             self.lcndbndy,self.laddconductor,self.icndbndy,
                             f3d.gridmode,conductorobject,self.lprecalccoeffs,
                             self.fsdecomp)
        else:
            # --- These are the plain multigrid solvers, since the full MG versions
            # --- have not yet be written.
            iondensitygrid3d = Grid3dtype()
            setupiondensitygrid3d(self.xmmin,self.ymmin,self.zmmin,
                                  self.dx,self.dy,self.dz,
                                  self.nxlocal,self.nylocal,self.nzlocal,
                                  self._rho,iondensitygrid3d)
            self.iondensitygrid3d = iondensitygrid3d
            multigridbe3dsolve(iwhich,self.nx,self.ny,self.nz,
                               self.nxguardphi,self.nyguardphi,self.nzguardphi,
                               self.nxguardrho,self.nyguardrho,self.nzguardrho,
                               self.dx,self.dy,self.dz*zfact,
                               self.potential,self.source,
                               rstar,self.linbend,self.bounds,
                               self.xmmin,self.ymmin,self.zmmin*zfact,
                               self.mgparam,mgiters,self.mgmaxiters,
                               self.mgmaxlevels,mgerror,self.mgtol,mgverbose,
                               self.downpasses,self.uppasses,
                               self.lcndbndy,self.laddconductor,self.icndbndy,
                               f3d.gridmode,conductorobject,
                               iondensitygrid3d,
                               self.fsdecomp)
        self.mgiters = mgiters[0]
        self.mgerror = mgerror[0]



##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MultiGridImplicit3D(MultiGrid3D):
    """
  This solves the modified Poisson equation which includes the suseptibility
  tensor that appears from the direct implicit scheme.
    """

    def __init__(self,lreducedpickle=1,withbadvance=true,**kw):
        MultiGrid3D.__init__(self,lreducedpickle,kwdict=kw)
        self.withbadvance = withbadvance
        self.solvergeom = w3d.XYZgeom
        self.ncomponents = 1

        # --- Kludge - make sure that the multigrid3df routines never sets up
        # --- any conductors. This is not really needed here.
        f3d.gridmode = 1

        # --- Save input parameters
        self.processdefaultsfrompackage(MultiGrid3D.__w3dinputs__,w3d,kw)
        self.processdefaultsfrompackage(MultiGrid3D.__f3dinputs__,f3d,kw)

        # --- If there are any remaning keyword arguments, raise an error.
        assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

        # --- Create conductor objects
        self.initializeconductors()

        # --- Give these variables dummy initial values.
        self.mgiters = 0
        self.mgerror = 0.

        # --- At the start, assume that there are no bends. This is corrected
        # --- in the solve method when there are bends.
        self.linbend = false

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
        'Returns the rho array without the guard cells'
        return MultiGrid3D.getrho(self)[:,:,:]

    def getrhop(self):
        'Returns the rhop array without the guard cells'
        return MultiGrid3D.getrhop(self)[:,:,:]

    def getphi(self):
        'Returns the phi array without the guard cells'
        return MultiGrid3D.getphi(self)[:,:,:]

    def getphip(self):
        'Returns the phip array without the guard cells'
        return MultiGrid3D.getphip(self)[:,:,:]

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
        if top.wpid == 0: wfact = zeros((0,), 'd')
        else:             wfact = pgroup.pid[i:i+n,top.wpid-1]
        depos_order = top.depos_order[:,js]
        self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,m,w,iimp,
                                  depos_order)

    def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,m,w,iimp,
                             depos_order):
        n  = len(x)
        if n == 0: return
        if iimp >= 0:
            # --- Create a temporary array to pass into setrho3d. This contributes
            # --- differently to the charge density and to chi. Also, make it a
            # --- 3-D array so it is accepted by setrho3d.
            sourcep = fzeros(self.sourcep.shape[:-1],'d')
        else:
            sourcep = self.sourcep[...,-1]
        if top.wpid == 0:
            setrho3d(sourcep,n,x,y,z,zgrid,q,w,top.depos,top.depos_order,
                     self.nxp,self.nyp,self.nzp,
                     self.nxguardrho,self.nyguardrho,self.nzguardrho,
                     self.dx,self.dy,self.dz,
                     self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                     self.solvergeom==w3d.RZgeom)
        else:
            # --- Need top.pid(:,top.wpid)
            setrho3dw(sourcep,n,x,y,z,zgrid,wfact,q,w,top.depos,top.depos_order,
                      self.nxp,self.nyp,self.nzp,
                      self.nxguardrho,self.nyguardrho,self.nzguardrho,
                      self.dx,self.dy,self.dz,
                      self.xmminp,self.ymminp,self.zmminp,self.l2symtry,self.l4symtry,
                      self.solvergeom==w3d.RZgeom)
        if iimp >= 0:
            self.sourcep[...,0] += sourcep
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
                                      self.nxp,self.nyp,self.nzp,
                                      self.sourcep[...,iimp],
                                      self.nxguardrho,self.nyguardrho,self.nzguardrho,
                                      self.fsdecomp,self.ppdecomp)

    def dosolve(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        if not self.l_internal_dosolve: return
        # --- set for longitudinal relativistic contraction
        if zfact is None:
            beta = top.pgroup.fselfb[iselfb]/clight
            zfact = 1./sqrt((1.-beta)*(1.+beta))
        else:
          beta = sqrt( (1.-1./zfact)*(1.+1./zfact) )

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source[...,0]
        if isinstance(self.potential,float): return

        # --- Setup data for bends.
        rstar = fzeros(3+self.nzlocal,'d')
        if top.bends:

            # --- This commented out code does the same thing as the line below
            # --- setting linbend but is a bit more complicated. It is preserved
            # --- in case of some unforeseen problem with the code below.
            #ii = (top.cbendzs <= self.zmmax+zgrid and
            #                     self.zmmin+zgrid <= top.cbendze)
            #self.linbend = sometrue(ii)

            setrstar(rstar,self.nzlocal,self.dz,self.zmminlocal,self.getzgrid())
            self.linbend = rstar.min() < largepos

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
        mgsolveimplicites3d(iwhich,self.nx,self.ny,self.nz,
                            self.nxlocal,self.nylocal,self.nzlocal,
                            self.dx,self.dy,self.dz*zfact,
                            self.potential,self._rho,
                            top.nsimplicit,qomdt,self.chi0,
                            rstar,self.linbend,self.withbadvance,
                            self.bounds,self.xmminlocal,self.ymminlocal,
                            self.zmminlocal*zfact,
                            self.getzgrid()*zfact,
                            self.mgparam,mgiters,self.mgmaxiters,
                            self.mgmaxlevels,mgerror,self.mgtol,mgverbose,
                            self.downpasses,self.uppasses,
                            self.lcndbndy,self.laddconductor,self.icndbndy,
                            f3d.gridmode,conductorobject,self.fsdecomp)

        self.mgiters = mgiters[0]
        self.mgerror = mgerror[0]

# --- This can only be done after MultiGrid3D is defined.
try:
    psyco.bind(MultiGrid3D)
    psyco.bind(MultiGridImplicit3D)
except NameError:
    pass
