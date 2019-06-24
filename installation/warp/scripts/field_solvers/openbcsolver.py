"""
Class for doing open boundary field solve in 3-D
----------------------------------------------
"""
from ..warp import *
import openbc_poisson

##############################################################################
class OpenBC3D(SubcycledPoissonSolver):
    """
  Creates a 3-D Poisson solver that has open boundary conditions.
  There are many input parameters, most of which get there default values from one of the fortran packages.

  | Input parameters from w3d:
  |    nx,ny,nz,dx,dy,dz,nxlocal,nylocal,nzlocal,
  |    nxguardphi,nyguardphi,nzguardphi,
  |    nxguardrho,nyguardrho,nzguardrho,
  |    nxguarde,nyguarde,nzguarde,
  |    xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,
  |    xmminlocal,xmmaxlocal,
  |    ymminlocal,ymmaxlocal,
  |    zmminlocal,zmmaxlocal,

  | Input parameters from top:
  |    nprocs,nxprocs,nyprocs,nzprocs,
  |    lfsautodecomp,zslave,debug

  | Other input paramters:
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
  |    igfflag: Whether or not to use the integrated green's function solver

    """

    def __init__(self,igfflag=True,lreducedpickle=1,**kw):
        self.igfflag = igfflag
        kw['lreducedpickle'] = lreducedpickle
        self.grid_overlap = 2

        # --- Make sure that the bounds have acceptable values.
        self.bounds = [3,3,3,3,3,3]

        SubcycledPoissonSolver.__init__(self,kwdict=kw)
        self.solvergeom = w3d.XYZgeom
        self.ncomponents = 1

        # --- If there are any remaning keyword arguments, raise an error.
        assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

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

    def __getstate__(self):
        dict = SubcycledPoissonSolver.__getstate__(self)
        if self.lreducedpickle:

            if '_rho' in dict: del dict['_rho']
            if '_phi' in dict: del dict['_phi']
            if 'bfieldp' in dict: del dict['bfieldp']

        return dict

    def getpdims(self):
        # --- Returns the dimensions of the arrays used by the particles

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
        gtlchg3dfromrho(self.nxlocal,self.nylocal,self.nzlocal,
                        self.nxguardrho,self.nyguardrho,self.nzguardrho,
                        self.rho,
                        self.dx,self.dy,self.dz,
                        self.getzgrid(),self.zmminlocal,
                        self.l2symtry,self.l4symtry,
                        self.izproc==self.nzprocs-1)

    def getese(self):
        'Calculate the electrostatic potential energy, rho*phi, and put it in top.ese'
        # --- ese must be an array to get the returned value in the calls below
        ese = zeros(1,'d')

        getese3dfromrhophi(self.nxlocal,self.nylocal,self.nzlocal,
                           self.nxguardphi,self.nyguardphi,self.nzguardphi,
                           self.nxguardrho,self.nyguardrho,self.nzguardrho,
                           self.rho,self.phi,
                           self.dx,self.dy,self.dz,self.l4symtry,self.l2symtry,
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

    def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
        # --- This is called by fetchfield from fieldsolver.py
        # --- Only sets the E field from the potential
        n = len(x)
        if n == 0: return
        if top.efetch[js] == 0: return
        if top.efetch[js] == 3 and isinstance(self.fieldp,float): return
        if top.efetch[js] != 3 and isinstance(self.potentialp,float): return
        sete3d(self.potentialp,self.fieldp,n,x,y,z,self.getzgridprv(),
               self.xmminp,self.ymminp,self.zmminp,
               self.dx,self.dy,self.dz,self.nxp,self.nyp,self.nzp,
               self.nxguardphi,self.nyguardphi,self.nzguardphi,
               self.nxguarde,self.nyguarde,self.nzguarde,
               top.efetch[js],top.depos_order[:,js],
               ex,ey,ez,self.l2symtry,self.l4symtry,self.solvergeom==w3d.RZgeom)

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

        if sometrue(top.efetch == 3) or maxnd(top.depos_order) > 1:
            self.setfieldpforparticles(*args)
            indts = args[1]
            # --- If this is the first group, set make sure that fieldp gets
            # --- zeroed out. Otherwise, the data in fieldp is accumulated.
            # --- This coding relies on the fact that fieldsolver does the
            # --- loops in descending order.
            tmpnsndts = getnsndtsforsubcycling()
            lzero = (indts == tmpnsndts-1)
            #if lzero:
            #  tfieldp = transpose(self.fieldp)
            #  tfieldp[...] = 0.
            self.calcselfep(recalculate=1,lzero=lzero)

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

    def hasconductors(self):
        return False

    def dosolve(self,iwhich=0,zfact=1.,isourcepndtscopies=None,indts=None,iselfb=None):
        if not self.l_internal_dosolve: return

        # --- This is only done for convenience.
        self._phi = self.potential
        self._rho = self.source
        if isinstance(self.potential,float): return

        # The last proc in each direction has one more cell, in order to
        # be able to calculate the finite difference for E.
        # (For the other procs, this extra cell is not included, as it
        # is obtained from a different proc)
        x_offset = int(top.ixproc == top.nxprocs-1)
        y_offset = int(top.iyproc == top.nyprocs-1)
        z_offset = int(top.izproc == top.nzprocs-1)
        # (x/y/z_offset is 1 if the extra cell is included, and 0 otherwise)
        rho = self._rho[self.nxguardrho:self.nxguardrho+self.nx+x_offset,
                        self.nyguardrho:self.nyguardrho+self.ny+y_offset,
                        self.nzguardrho:self.nzguardrho+self.nz+z_offset]
        phi = self._phi[self.nxguardphi:self.nxguardphi+self.nx+x_offset,
                        self.nyguardphi:self.nyguardphi+self.ny+y_offset,
                        self.nzguardphi:self.nzguardphi+self.nz+z_offset]

        # Local size of the arrays
        ilo, ihi = 1, self.nx + x_offset
        jlo, jhi = 1, self.ny + y_offset
        klo, khi = 1, self.nz + z_offset
        # Global size of the arrays
        ilo_rho_gbl, ihi_rho_gbl = ilo, ihi
        jlo_rho_gbl, jhi_rho_gbl = jlo, jhi
        klo_rho_gbl, khi_rho_gbl = klo, khi

        # MPI decomposition
        idecomp = -1
        # Number of processors in each dimension
        nxp = nyp = nzp = 1

        igfflag = ((self.igfflag and 1) or 0) # for normal green function (1 for IGF)
        self.ierr = zeros(1, 'l')

        # Calculate charge per cell
        charge = rho * (self.dx * self.dy * self.dz * zfact)
        # Scale the charge, to match the assumed units of `openbcpotential`
        # so that the returned potential phi is in SI units
        charge *= 1./(4.*pi*eps0)
        openbc_poisson.openbcpotential(
            charge, phi, self.dx, self.dy, self.dz*zfact,
            ilo, ihi, jlo, jhi, klo, khi,
            ilo_rho_gbl, ihi_rho_gbl, jlo_rho_gbl, jhi_rho_gbl,
            klo_rho_gbl, khi_rho_gbl, idecomp, nxp, nyp, nzp,
            igfflag, self.ierr)

        # For now, set phi in the guard cell to be equal to the value in
        # in the last physical cell (this effectively sets E to 0 in the
        # guard cells)
        # - Along x
        self._phi[:self.nxguardphi,:,:] = \
            (self._phi[self.nxguardphi,:,:])[newaxis,:,:]
        self._phi[ self.nxguardphi+self.nx+x_offset:,:,: ] = \
            (self._phi[self.nxguardphi+self.nx+x_offset-1,:,:])[newaxis,:,:]
        # - Along y
        self._phi[:,:self.nyguardphi,:] = \
            (self._phi[:,self.nyguardphi,:])[:,newaxis,:]
        self._phi[ :,self.nyguardphi+self.ny+y_offset:,: ] = \
            (self._phi[ :,self.nyguardphi+self.ny+y_offset-1,:])[:,newaxis,:]
        # - Along z
        self._phi[:,:,:self.nzguardphi] = \
            (self._phi[:,:,self.nzguardphi])[:,:,newaxis]
        self._phi[ :,:,self.nzguardphi+self.nz+z_offset:] = \
            (self._phi[ :,:,self.nzguardphi+self.nz+z_offset-1])[:,:,newaxis]

        # --- Note that the guard cells and upper edge of the domain are not set yet.
        # --- A future fix will be to pass in the entire grid into the solver, including the guard
        # --- cells. Though, this would require the rho array to be the same size as phi.
        #applyboundaryconditions3d(self.nx,self.ny,self.nz,self.nxguardphi,self.nyguardphi,self.nzguardphi,
        #                          self._phi,1,self.bounds,false,false)

    ##########################################################################
    # Define the basic plot commands
    def genericpf(self,kw,pffunc):
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
        reps0c = 1./(eps0*2.*(dxsqi+dysqi+dzsqi))
        rho = self._rho*reps0c
        conductorobject = ConductorType()
        residual3d(self.nxlocal,self.nylocal,self.nzlocal,
                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                   self.nxguardrho,self.nyguardrho,self.nzguardrho,
                   self.nxguardphi,self.nyguardphi,self.nzguardphi,
                   dxsqi,dysqi,dzsqi,self._phi,rho,res,
                   0,self.bounds,1.,1,true,
                   false,0,conductorobject,false)
        return res
