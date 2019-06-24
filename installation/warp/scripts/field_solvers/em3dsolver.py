"""Class for doing 3 D electromagnetic solver """
from ..warp import *
from ..diagnostics.palettes.mkpalette import getpalhrgb
from .laser.laser_antenna import LaserAntenna
import collections
import types
import operator

try:
    import Opyndx
    l_opyndx = 1
except:
    l_opyndx = 0

try:
    import psyco
except ImportError:
    pass

class EM3D(SubcycledPoissonSolver):

    __em3dinputs__ = []
    __flaginputs__ = {'mode':1,'l_apply_pml':true,
                      'nbndx':10,'nbndy':10,'nbndz':10,
                      'nxguard':1,'nyguard':1,'nzguard':1,
                      'norderx':2,'nordery':2,'norderz':2,
                      'l_particles_weight':false,'l_usecoeffs':false,
                      'l_pushf':false,'l_pushpot':false,'l_verbose':false,
                      'l_nodalgrid':false, # flag for nodal/staggered grid
                      'l_pushg':false,
                      'laser_func':None,
                      'laser_polangle':0.,
                      'laser_vector':array([0.,0.,1.]),
                      'laser_polvector':None,
                      'laser_spot':None,
                      'laser_source_z':None,
                      'laser_source_v':array([0., 0., 0.]),
                      'laser_emax':None,
                      'laser_depos_order_x':3,
                      'laser_depos_order_y':3,
                      'laser_depos_order_z':3,
                      'ntsub':1,
                      'l_enableovercycle':False,
                      'l_2dxz':0,'l_2drz':0,'l_1dz':0,'l_sumjx':0,
                      'l_lower_order_in_v':True,
                      'npass_smooth':array([[0],[0],[0]]),
                      'alpha_smooth':array([[0.5],[0.5],[0.5]]),
                      'stride_smooth':array([[1],[1],[1]]),
                      'mask_smooth':None,
                      'l_smooth_particle_fields':False,
                      'n_smooth_fields':None,
                      'autoset_timestep':true,'dtcoef':1.,#0.99,
                      'deposit_energy_density':false,'refinement':None,
                      'l_force_nzlocal2nz':false,'isactiveem':None,
                      'l_coarse_patch':false,
                      'stencil':false, # 0=Yee; 1=Cole-Karkkainen, 3=Lehe; if spectral=1, is used only for computing time step according to Courant condition of selected solver
                      'l_esirkepov':true, # flag for deposition using Esirkepov algorithm (real space on staggered grid with Yee EM solver)
                      'l_deposit_nodal':None,
                      'l_getrho':false,'theta_damp':0.,
                      'sigmae':0.,'sigmab':0.,
                      'colecoefs':None,'l_setcowancoefs':True,
                      'pml_method':1,
                      'vxgrid':0.,
                      'vygrid':0.,
                      'vzgrid':0.,
                      'l_correct_num_Cherenkov':False,
                      'l_fieldcenterK':False, # if using staggered grid with node-centered gather (efetch=1); centers field by shifts in k-space rather than averaging in real space
                      'V_galilean':array([0.,0.,0.]),
                      'V_pseudogalilean':array([0.,0.,0.]),
                      'circ_m':0, 'l_laser_cart':0, 'type_rz_depose':0}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.solveroff=False # flag to turn off the solver, for testing purpose
        assert top.grid_overlap<2,"The EM solver needs top.grid_overlap<2!"
        self.grid_overlap = 1
        FieldSolver.__init__(self,kwdict=kw)
#    top.allspecl = true
        top.lcallfetchb = true
        top.lgridqnt = true
        self.x_grid=0.
        self.x_gridcont=0.
        self.y_grid=0.
        self.y_gridcont=0.
        self.z_grid=0.
        self.z_gridcont=0.
        self.zgrid=top.zgrid
        self.odd=0
        # --- Save input parameters
#    self.processdefaultsfrompackage(EM3D.__w3dinputs__,w3d,kw)
#    self.processdefaultsfrompackage(EM3D.__topinputs__,top,kw)
        self.processdefaultsfrompackage(EM3D.__em3dinputs__,em3d,kw)
        self.processdefaultsfromdict(EM3D.__flaginputs__,kw)
        if self.circ_m>0:w3d.solvergeom=w3d.RZgeom
        if w3d.solvergeom==w3d.XZgeom:self.l_2dxz=True
        if w3d.solvergeom==w3d.RZgeom:self.l_2drz=True
        if w3d.solvergeom==w3d.Zgeom:
            self.l_2dxz=True
            self.l_1dz=True
        if self.l_2dxz: w3d.solvergeom=w3d.XZgeom
        if self.l_2drz:
            w3d.solvergeom=w3d.RZgeom
            self.l_2dxz=True
        if self.l_1dz:
            w3d.solvergeom=w3d.Zgeom
            self.stencil=0
        self.solvergeom = w3d.solvergeom
        particleboundaries3d(top.pgroup,-1,False)
        if self.isactiveem is not None:
            self.isactive=self.isactiveem
        else:
            try:
                self.isactive=self.isactive
            except:
                self.isactive=True
        # --- When initializing self.field_coarse, there is some inconsistency
        # --- with the way that FieldSolver.__init__ resize nzlocal in parallel.
        # --- This flag takes care of it for now.
        if self.l_force_nzlocal2nz:
            self.nxlocal=self.nx
            self.nylocal=self.ny
            self.nzlocal=self.nz

        # --- Impose that rho be calculated, if l_pushf is active
        if self.l_pushf:
            self.l_getrho = True

        # --- Set l_deposit_nodal if not set by user
        if self.l_esirkepov:
            self.l_deposit_nodal = False
        else:
            if self.l_nodalgrid:
                self.l_deposit_nodal=True
            else:
                self.l_deposit_nodal=False

        # --- Impose type_rz_depose = 0 if not in circ mode
        if self.l_2drz == False :
            self.type_rz_depose = 0

        # --- smoothing lists to arrays if needed
        self.npass_smooth  = array(self.npass_smooth)
        self.alpha_smooth  = array(self.alpha_smooth)
        self.stride_smooth  = array(self.stride_smooth)

        # --- set number of guard cells to appropriate value depending on order of
        # --- current deposition, smoothing and stencil.
        minguards = 2+aint(top.depos_order.max(1)/2)+(self.npass_smooth*self.stride_smooth).sum(1)
        # --- add guards cells for larger stencils
        if self.nxguard==1 and self.norderx is not inf:
            minguards[0] += self.norderx/2+1
            self.nxguard = minguards[0]
        if self.nyguard==1 and self.nordery is not inf:
            minguards[1] += self.nordery/2+1
            self.nyguard = minguards[1]
        if self.nzguard==1 and self.norderz is not inf:
            minguards[2] += self.norderz/2+1
            self.nzguard = minguards[2]
        # --- ensures that there is at least 2 guard cells with stencil>0
        if self.stencil>0:
            if self.nxguard<2:self.nxguard=2
            if self.nyguard<2:self.nyguard=2
            if self.nzguard<2:self.nzguard=2
        # --- sets number of cells and guards to 0 along unused dimension
        if self.l_2dxz:
            self.nyguard=self.nylocal=self.ny=self.nyp=0
        if self.l_1dz:
            self.l_2dxz = True
            self.nxguard=self.nxlocal=self.nx=self.nxp=0
            self.nyguard=self.nylocal=self.ny=self.nyp=0
        # --- enforces number of guards cells to not go over the number of local cells
        self.nxguard = min(self.nxguard,self.nxlocal)
        self.nyguard = min(self.nyguard,self.nylocal)
        self.nzguard = min(self.nzguard,self.nzlocal)

        # --- bounds is special since it will sometimes be set from the
        # --- variables bound0, boundnz, boundxy, l2symtry, and l4symtry
        if 'bounds' not in self.__dict__:
            if 'bounds' in kw:
                self.bounds = kw['bounds']
            else:
                self.bounds = zeros(6,'l')
                self.bounds[0] = self.boundxy
                self.bounds[1] = self.boundxy
                self.bounds[2] = self.boundxy
                self.bounds[3] = self.boundxy
                self.bounds[4] = self.bound0
                self.bounds[5] = self.boundnz
                if self.l2symtry:
                    self.bounds[2] = neumann
                    if self.boundxy == periodic: self.bounds[3] = neumann
                    if self.forcesymmetries: self.ymmin = 0.
                elif self.l4symtry:
                    self.bounds[0] = neumann
                    self.bounds[2] = neumann
                    if self.boundxy == periodic: self.bounds[1] = neumann
                    if self.boundxy == periodic: self.bounds[3] = neumann
                    if self.forcesymmetries: self.xmmin = 0.
                    if self.forcesymmetries: self.ymmin = 0.
        if self.l_2dxz:self.bounds[2:4] = -1
        if self.l_2drz:self.bounds[0]=-1
        if self.l_1dz:self.bounds[:2] = -1
        self.bounds = self.bounds.copy()

        # --- removes bounds from kw to prevent conflict with bounds created
        # --- by the MR class, when adding MR children.
        try:
            if 'bounds' in kw:kw.pop('bounds')
        except:
            pass

        # --- If there are any remaning keyword arguments, raise an error.
        assert len(kw.keys()) == 0,"Bad keyword arguments %s"%kw.keys()

        # --- This needs to be called again since some mesh values may have changed
        # --- since it was previously called in FieldSolver.__init__.
        # --- For instance, if using RZ geometry, ny will be set to zero.
        self.setupmeshextent()

        # --- sets coefficients of Cole solver
        # Lehe stencil (see Lehe et al., PRSTAB 16 021301 (2013))
        if self.stencil == 3 :
        # Warning : the coefficients alphaz and deltaz are calculated
        # later in the file, i.e. only once the dt has been calculated.
            if self.l_2dxz:
                if self.l_2drz:
                    # 2D cylindrical
                    em3d.betaxz = 1./4
                    em3d.betazx = 0.
                    em3d.betayx = 0.
                    em3d.betaxy = 0.
                    em3d.betayz = 1./4
                    em3d.betazy = 0.
                else:
                    # 2D Cartesian
                    em3d.betaxz = 1./8
                    em3d.betazx = self.dz**2/self.dx**2*1./8
                    em3d.betayx = 0.
                    em3d.betaxy = 0.
                    em3d.betayz = 0.
                    em3d.betazy = 0.
            else :
                # 3D Cartesian
                em3d.betaxz = 1./8
                em3d.betazx = self.dz**2/self.dx**2*1./8
                em3d.betayx = 0.
                em3d.betaxy = 0.
                em3d.betayz = 1./8
                em3d.betazy = self.dz**2/self.dy**2*1./8
        elif self.l_setcowancoefs:
            if self.l_2dxz:
                delta = min(self.dx,self.dz)
                rx = (delta/self.dx)**2
                ry = 0.
                rz = (delta/self.dz)**2
                beta = 0.125*(1.-rx*ry*rz/(ry*rz+rz*rx+rx*ry))
                em3d.betaxz = 0.125*rz
                em3d.betazx = 0.125*rx
                em3d.alphax = 1. - 2.*em3d.betaxz
                em3d.alphaz = 1. - 2.*em3d.betazx
            else:
                delta = min(self.dx,self.dy,self.dz)
                rx = (delta/self.dx)**2
                ry = (delta/self.dy)**2
                rz = (delta/self.dz)**2
                beta = 0.125*(1.-rx*ry*rz/(ry*rz+rz*rx+rx*ry))
                em3d.betaxy = ry*beta
                em3d.betaxz = rz*beta
                em3d.betayx = rx*beta
                em3d.betayz = rz*beta
                em3d.betazx = rx*beta
                em3d.betazy = ry*beta
                em3d.gammax = ry*rz*(1./16.-0.125*ry*rz/(ry*rz+rz*rx+rx*ry))
                em3d.gammay = rx*rz*(1./16.-0.125*rx*rz/(ry*rz+rz*rx+rx*ry))
                em3d.gammaz = rx*ry*(1./16.-0.125*rx*ry/(ry*rz+rz*rx+rx*ry))
                em3d.alphax = 1. - 2.*em3d.betaxy - 2.* em3d.betaxz - 4.*em3d.gammax
                em3d.alphay = 1. - 2.*em3d.betayx - 2.* em3d.betayz - 4.*em3d.gammay
                em3d.alphaz = 1. - 2.*em3d.betazx - 2.* em3d.betazy - 4.*em3d.gammaz
        else:
            if self.colecoefs is None:
                if self.l_2dxz:
                    em3d.betaxy = 1./8.
                    em3d.alphax = 1.-2*em3d.betaxy
                else:
                    em3d.alphax = 7./12.
                    em3d.betaxy  = 1./12.
                    em3d.gammax = 1./48.
            else:
                em3d.alphax = self.colecoefs[0]
                em3d.betaxy  = self.colecoefs[1]
                if not self.l_2dxz:
                    em3d.gammax = self.colecoefs[2]
            em3d.alphay = em3d.alphaz = em3d.alphax
            em3d.betaxz  = em3d.betayx = em3d.betayz = em3d.betazx = em3d.betazy  = em3d.betaxy
            em3d.gammay = em3d.gammaz = em3d.gammax

        print 'alphax,alphaz,betazx,betaxz',em3d.alphax,em3d.alphaz,em3d.betazx,em3d.betaxz
        # --- set time step as a fraction of Courant condition
        # --- also set self.ntsub if top.dt over Courant condition times dtcoef
        try:
            parentid = self.parents[0]
            parent = self.root.listofblocks[parentid]
            try:
                sibling = parent.children[0]
            except:
                sibling = None
        except:
            parent = None
            sibling = None
        if sibling is not None:            # No mesh refinement
                self.ntsub=sibling.ntsub
        else:
                if self.autoset_timestep:  # True by default (optional argument of EM3D)
#                   if self.mode==2 and self.dtcoef>0.5: self.dtcoef/=2 # Obsolete
                    if self.l_1dz:
                    ### - 1D
                        self.dtcourant=self.dz/clight
                    elif self.l_2dxz:
                    ### - 2D
                        if self.l_2drz == False : # 2D Cartesian
                            if self.stencil==0:  # Yee scheme
                                self.dtcourant=1./(clight*sqrt(1./self.dx**2+1./self.dz**2))
                            elif self.stencil in [1,2] : # Cole-Karkkainen scheme
                                self.dtcourant=min(self.dx,self.dz)/clight
                                Cx = em3d.alphax -2.*em3d.betaxz
                                Cz = em3d.alphaz -2.*em3d.betazx
                                self.dtcourant=1./(clight*sqrt(Cx/self.dx**2+Cz/self.dz**2))
                            elif self.stencil == 3 : # Lehe scheme
                                self.dtcourant = 1./clight  * min( self.dz, self.dx )
                        else :  # 2D r-z
                            if self.stencil==1:
                                raise ValueError(
                                    "The Cole-Karkkainen solver (stencil=1) cannot be used in cylindrical geometry.")
                            elif self.stencil==3: # Lehe scheme
                                self.dtcourant = 1./clight  * min( self.dz, self.dx )

                            else:  # Yee scheme
                                # In the rz case, the Courant limit has been evaluated
                                # semi-analytically by R. Lehe, and resulted in the following
                                # coefficients. For an explanation, see (not officially published)
                                # www.normalesup.org/~lehe/Disp_relation_Circ.pdf
                                # NB : Here the coefficient for m=1 as compared to this document,
                                # as it was observed in practice that this coefficient was not
                                # high enough (The simulation became unstable).
                                circ_coeffs = [ 0.2105, 1.0, 3.5234, 8.5104, 15.5059, 24.5037 ]
                                if self.circ_m < len(circ_coeffs) : # Use the table of the coefficients
                                    circ_alpha = circ_coeffs[self.circ_m]
                                else : # Use a realistic extrapolation
                                    circ_alpha = self.circ_m**2 - 0.4
                                self.dtcourant=1./(clight*sqrt((1+circ_alpha)/self.dx**2+1./self.dz**2))
                    else:
                    ### - 3D
                        if self.stencil==0:
                            self.dtcourant=1./(clight*sqrt(1./self.dx**2+1./self.dy**2+1./self.dz**2))
                        elif self.stencil in [1,2] : # Cole-Karakkainen scheme
                            self.dtcourant=min(self.dx,self.dy,self.dz)/clight
                            Cx = em3d.alphax -2.*(em3d.betaxy+em3d.betaxz)+4.*em3d.gammax
                            Cy = em3d.alphay -2.*(em3d.betayx+em3d.betayz)+4.*em3d.gammay
                            Cz = em3d.alphaz -2.*(em3d.betazx+em3d.betazy)+4.*em3d.gammaz
                            self.dtcourant=1./(clight*sqrt(Cx/self.dx**2+Cy/self.dy**2+Cz/self.dz**2))
                        elif self.stencil == 3 : # Lehe scheme
                            self.dtcourant = 1./clight * min( self.dz, 1./sqrt( 1./self.dx**2 + 1./self.dy**2 ) )
                    if self.theta_damp>0.:
                        self.dtcourant*=sqrt((2.+self.theta_damp)/(2.+3.*self.theta_damp))
                    if top.dt==0.:
                        top.dt=self.dtcourant*self.dtcoef
                    if self.ntsub is not inf and self.ntsub<2:
                       if top.dt>(self.dtcourant):
#              self.ntsub = (nint(top.dt/(self.dtcourant))+0)
                            self.ntsub = int(top.dt/(self.dtcourant)+1.)
                            print '#1', self.ntsub,top.dt,self.dtcourant
                       elif self.l_enableovercycle:
                            self.ntsub = 1./(nint((self.dtcourant)/top.dt)+0)
                            print '#2', self.ntsub,top.dt,self.dtcourant
        dtodz = clight*top.dt/self.dz
        if self.l_correct_num_Cherenkov and top.efetch[0]==4 and self.l_lower_order_in_v and dtodz>0.756 and dtodz<0.764:
            print "*** Warning: Coefficients for Galerkin algorithm are ill behaved for 0.756<c*Dt/Dz<0.764 and should not be used."
            print "*** Rescaling top.dt to 0.75*Dz/c"
            top.dt = 0.75*self.dz/clight

        self.dtinit = top.dt

        if self.stencil == 3 : # Lehe stencil, calculation of the remaining coefficients
            em3d.deltaz = -0.25*( self.dz**2/(clight*top.dt)**2 * sin( pi*clight*top.dt/(2.*self.dz) )**2 - 1 )
            em3d.alphaz = 1. - 2.*em3d.betazx - 2.* em3d.betazy - 3.*em3d.deltaz
            em3d.alphax = 1. - 2.*em3d.betaxy - 2.* em3d.betaxz
            em3d.alphay = 1. - 2.*em3d.betayx - 2.* em3d.betayz


        if  top.vbeamfrm != 0.:self.bounds[-2:]=-1

        if self.refinement is not None:
            ref=self.refinement
            self.field_coarse = self.__class__(l_force_nzlocal2nz=True,
                                     l_coarse_patch=True,
                                     nx=self.nx/ref[0],dx=self.dx*ref[0],
                                     ny=self.ny/ref[1],dy=self.dy*ref[1],
                                     nz=self.nz/ref[2],dz=self.dz*ref[2],
                                     nxlocal=self.nxlocal/ref[0],
                                     nylocal=self.nylocal/ref[1],
                                     nzlocal=self.nzlocal/ref[2],
                                     xmmin=self.xmmin,xmmax=self.xmmax,
                                     ymmin=self.ymmin,ymmax=self.ymmax,
                                     zmmin=self.zmmin,zmmax=self.zmmax,
                                     xmminlocal=self.xmminlocal,xmmaxlocal=self.xmmaxlocal,
                                     ymminlocal=self.ymminlocal,ymmaxlocal=self.ymmaxlocal,
                                     zmminlocal=self.zmminlocal,zmmaxlocal=self.zmmaxlocal,
                                     #bounds=self.bounds,
                                     isactiveem=self.isactive,
                                     ntsub=self.root.listofblocks[self.parents[0]].ntsub,
                                     lchild=True,
                                     **self.kw)

#   if not self.l_1dz:
#     if self.xmmin>w3d.xmmaxlocal or self.xmmax<w3d.xmminlocal:return
#   if not self.l_2dxz and not self.l_2drz:
#     if self.ymmin>w3d.ymmaxlocal or self.ymmax<w3d.ymminlocal:return
#   if self.zmmin>w3d.zmmaxlocal or self.zmmax<w3d.zmminlocal:return

        # --- Handle laser inputs
        self.setuplaser()

        self.initializeconductors()

        # ---- install 2nd part of field solve (only for root)
        if self.refinement is None and not self.l_coarse_patch:
            installbeforeloadrho(self.solve2ndhalf)

        self.finalized = False

        # --- This allows some introspection.
        if self.l_nodalgrid:
            self.dict_of_grids = {
              'rho':{'getter':'getrho', 'centering':'node', 'units':'C/m**3'},
              'Jx':{'getter':'getjx', 'centering':'Yee', 'units':'A/m**2'},
              'Jy':{'getter':'getjy', 'centering':'Yee', 'units':'A/m**2'},
              'Jz':{'getter':'getjz', 'centering':'Yee', 'units':'A/m**2'},
              'Ex':{'getter':'getexg', 'centering':'node', 'units':'V/m'},
              'Ey':{'getter':'geteyg', 'centering':'node', 'units':'V/m'},
              'Ez':{'getter':'getezg', 'centering':'node', 'units':'V/m'},
              'Bx':{'getter':'getbxg', 'centering':'node', 'units':'T'},
              'By':{'getter':'getbyg', 'centering':'node', 'units':'T'},
              'Bz':{'getter':'getbzg', 'centering':'node', 'units':'T'},
              'divE':{'getter':'getdive', 'centering':'node', 'units':'V/m**2'},
            }
        else:
            self.dict_of_grids = {
              'rho':{'getter':'getrho', 'centering':'node', 'units':'C/m**3'},
              'Jx':{'getter':'getjx', 'centering':'Yee', 'units':'A/m**2'},
              'Jy':{'getter':'getjy', 'centering':'Yee', 'units':'A/m**2'},
              'Jz':{'getter':'getjz', 'centering':'Yee', 'units':'A/m**2'},
              'Ex':{'getter':'getexg', 'centering':'Yee', 'units':'V/m'},
              'Ey':{'getter':'geteyg', 'centering':'Yee', 'units':'V/m'},
              'Ez':{'getter':'getezg', 'centering':'Yee', 'units':'V/m'},
              'Bx':{'getter':'getbxg', 'centering':'Yee', 'units':'T'},
              'By':{'getter':'getbyg', 'centering':'Yee', 'units':'T'},
              'Bz':{'getter':'getbzg', 'centering':'Yee', 'units':'T'},
              'divE':{'getter':'getdive', 'centering':'node', 'units':'V/m**2'},
            }

    def processdefaultsfrompackage(self,defaults,package,kw):
        for name in defaults:
            if name not in self.__dict__:
                #self.__dict__[name] = kw.pop(name,getattr(w3d,name)) # Python2.3
                self.__dict__[name] = kw.get(name,getattr(package,name))
            if name in kw: del kw[name]

    def processdefaultsfromdict(self,dict,kw):
        for name,defvalue in dict.iteritems():
            if name not in self.__dict__:
                #self.__dict__[name] = kw.pop(name,getattr(top,name)) # Python2.3
                self.__dict__[name] = kw.get(name,defvalue)
            if name in kw: del kw[name]

    def getpdims(self):
        # --- Returns the dimensions of the arrays used by the particles

        # --- If there are any relativistic groups, then turn on the code
        # --- which uses the selfe array.
        if max(abs(top.fselfb)) > 0.:
            # --- This is probably redundant, but it shouldn't hurt.
            # --- This forces all species to use the precalculated E field
            # --- if any have the B correction.
            top.efetch = 3
            # --- Number of fields (E and B)
            nfields = 2
        else:
            # --- Number of fields (E only)
            nfields = 1

        if sometrue(top.efetch == 3):
            return ((1+self.nxp,1+self.nyp,1+self.nzp),
                    (3,1+self.nxp,1+self.nyp,1+self.nzp,nfields),
                    (1+self.nxp+2*self.nxguard,1+self.nyp+2*self.nyguard,1+self.nzp+2*self.nzguard))
        else:
            return ((1+self.nxp,1+self.nyp,1+self.nzp),
                    (1+self.nxp+2*self.nxguard,1+self.nyp+2*self.nyguard,1+self.nzp+2*self.nzguard))

    def getdims(self):
        # --- Returns the dimensions of the arrays used by the field solver
        return ((1+self.nxlocal,1+self.nylocal,1+self.nzlocal),
                (1+self.nxlocal+2*self.nxguard,
                 1+self.nylocal+2*self.nyguard,
                 1+self.nzlocal+2*self.nzguard))

    def finalize(self,lforce=False):
        if self.finalized and not lforce: return
        self.setbcparallel(0) # x
        self.setbcparallel(1) # y
        self.setbcparallel(2) # z

        # --- Create field and source arrays and other arrays.
        self.allocatefieldarrays()

        w3d.lfinalizerhofsapi=True

        self.finalized = True

    def allocatefieldarrays(self):

        # --- Sets correction of Numerical Cherenkov
        if self.l_correct_num_Cherenkov:
            self.set_num_Cherenkov_cor_coefs()
            excoef,bycoef = self.get_num_Cherenkov_cor_coefs()
        else:
            excoef=bycoef=None

        self.block = pyinit_3dem_block(self.nxlocal,
                                       self.nylocal,
                                       self.nzlocal,
                                       self.nbndx,
                                       self.nbndy,
                                       self.nbndz,
                                       self.nxguard,
                                       self.nyguard,
                                       self.nzguard,
                                       self.norderx,
                                       self.nordery,
                                       self.norderz,
                                       top.dt,
                                       self.dx,
                                       self.dy,
                                       self.dz,
                                       clight,
                                       mu0,
                                       self.xmminlocal,
                                       self.ymminlocal,
                                       self.zmminlocal,
                                       int(self.bounds[0]),
                                       int(self.bounds[2]),
                                       int(self.bounds[4]),
                                       int(self.bounds[1]),
                                       int(self.bounds[3]),
                                       int(self.bounds[5]),
                                       self.deposit_energy_density,
                                       self.l_nodalgrid,
                                       self.refinement,
                                       self.l_pushf,
                                       self.l_pushg,
                                       self.l_getrho,
                                       self.stencil,
                                       self.npass_smooth,
                                       self.l_smooth_particle_fields,
                                       self.ntsub,
                                       self.l_1dz,
                                       self.l_2dxz,
                                       self.l_2drz,
                                       self.theta_damp,
                                       self.sigmae,
                                       self.sigmab,
                                       excoef,
                                       bycoef,
                                       self.pml_method,
                                       self.circ_m)
        self.fields = self.block.core.yf
        if self.l_2drz:
            self.vol = 2.*pi*(arange(self.nx+1)*self.dx+self.block.xmin)*self.dx*self.dz
            if self.block.xmin==0.:self.vol[0] = 0.25*self.dx**2*self.dz


################################################################################
#                                   LASER
################################################################################

#===============================================================================
    def setuplaser(self):
#===============================================================================
        if self.l_1dz:
            dim = "1d"
        elif self.circ_m > 0:
            dim = "circ"
        elif self.l_2dxz:
            dim = "2d"
        else:
            dim = "3d"

        # Initialize the list of antennas
        self.laser_antenna = []

        # Add a new antenna to the list if self.laser_func is not None
        if self.laser_func is not None:
            self.laser_antenna.append(LaserAntenna(
                self.laser_func, self.laser_vector,
                self.laser_polvector, self.laser_spot,
                self.laser_emax, self.laser_source_z,
                self.laser_source_v,
                self.laser_polangle,
                w3d, dim, self.circ_m))


#===============================================================================
    def add_laser(self,field, laser_antenna_i):
        """
        Calculate the parameters of the corresponding antenna and depose the
        corresponding current on the grid.

        Depose the source that generates the laser, on the grid (generation of
        the laser by an antenna).

        Notice that this current is generated by fictious macroparticles, whose
        motion is explicitly specified. Thus, this current does not correspond
        to that of the actual macroparticles of the simulation.

        Parameter :
        -----------
        field : EM3D_YEEFIELD
            The self.fields object, whose attributes are the field arrays
            Ex, Ey, etc ...

        laser_antenna_i: LaserAntenna
            ith component of the object self.laser_antenna
        """
        # If the user did not request a laser antenna,
        # this function should not be executed
        if laser_antenna_i is None:
            return

        f = self.block.core.yf

        # Push particles accordingly to laser_func
        laser_antenna_i.push_virtual_particles(top, f, clight)
        laser_antenna_i.select_particles_in_local_box( w3d, self.zgrid )
        # Skip current deposition when no laser particles are in the local box
        if laser_antenna_i.nn == 0:
            return

        # Current and charge deposition
        if top.ndts[0] != 1:
            print "Error in depose_j_laser: top.ndts[0] must be 1 if injecting\
                   a laser"
            raise
        f.Jx = self.fields.Jxarray[:,:,:,0]
        f.Jy = self.fields.Jyarray[:,:,:,0]
        f.Jz = self.fields.Jzarray[:,:,:,0]
        f.Rho = self.fields.Rhoarray[:,:,:,0]

        # q represents the sign of the charged macroparticles
        # The antenna is made of two types of fictious particles :
        # positive and negative
        l_particles_weight = True
        for q in [1.,-1.]:
            # Current deposition
            self.depose_current_density(laser_antenna_i.nn, f,
                           laser_antenna_i.xx + q*laser_antenna_i.xdx,
                           laser_antenna_i.yy + q*laser_antenna_i.ydy,
                           laser_antenna_i.zz + q*laser_antenna_i.zdz,
                           laser_antenna_i.vx + q*laser_antenna_i.ux,
                           laser_antenna_i.vy + q*laser_antenna_i.uy,
                           laser_antenna_i.vz + q*laser_antenna_i.uz,
                           laser_antenna_i.gi,
                           top.dt,
                           laser_antenna_i.weights,
                           self.zgrid, q, 1.,
                           self.laser_depos_order_x,
                           self.laser_depos_order_y,
                           self.laser_depos_order_z,
                           l_particles_weight)

            if self.l_getrho :
                # Charge deposition
                self.depose_charge_density(laser_antenna_i.nn, f,
                           laser_antenna_i.xx + q*laser_antenna_i.xdx,
                           laser_antenna_i.yy + q*laser_antenna_i.ydy,
                           laser_antenna_i.zz + q*laser_antenna_i.zdz,
                           laser_antenna_i.vx + q*laser_antenna_i.ux,
                           laser_antenna_i.vy + q*laser_antenna_i.uy,
                           laser_antenna_i.vz + q*laser_antenna_i.uz,
                           laser_antenna_i.gi,
                           top.dt,
                           laser_antenna_i.weights,
                           self.zgrid, q, 1.,
                           self.laser_depos_order_x,
                           self.laser_depos_order_y,
                           self.laser_depos_order_z,
                           l_particles_weight)

        # Call a function of the class em3dsolverFFT
        if hasattr( self, 'current_cor' ) and self.current_cor:
            self.depose_rho_old_laser(f, laser_antenna_i)


################################################################################
# FIELD FETCHING
################################################################################

#===============================================================================
    def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
#===============================================================================
        # --- This is called by fetchfield from fieldsolver.py
        if not self.finalized: return
        n = len(x)
        if n == 0: return
        nox = top.depos_order[0,w3d.jsfsapi]
        noy = top.depos_order[1,w3d.jsfsapi]
        noz = top.depos_order[2,w3d.jsfsapi]
#    if (not (self.getconductorobject(top.fselfb[iselfb]).lcorrectede or
#    else:
        f = self.block.core.yf

        # --- fetch e
        if top.efetch[w3d.jsfsapi] in [1,3,5]:
            if self.l_1dz:
                getf1dz_n(n,z,ex,ey,ez,
                         f.zmin+self.zgrid,
                         f.dz,
                         f.nz,
                         f.nzguard,
                         noz,
                         f.Exp,f.Eyp,f.Ezp)
            elif self.l_2dxz:
                if f.circ_m == 0 :
                    getf2dxz_n(n,x,y,z,ex,ey,ez,
                         f.xmin,f.zmin+self.zgrid,
                         f.dx,f.dz,
                         f.nx,f.ny,f.nz,
                         f.nxguard,f.nyguard,f.nzguard,
                         nox,noz,
                         f.Exp,f.Eyp,f.Ezp,
                         w3d.l4symtry,self.l_2drz)
                else :
                    getf2drz_circ_n(n,x,y,z,ex,ey,ez,
                             f.xmin,f.zmin+self.zgrid,
                             f.dx,f.dz,
                             f.nx, f.ny, f.nz,
                             f.nxguard, f.nyguard, f.nzguard,
                             nox,noz,
                             f.Exp, f.Eyp, f.Ezp,
                             f.Exp_circ,f.Eyp_circ,f.Ezp_circ,
                             f.circ_m)
            else:
                if nox==1 and noy==1 and noz==1 and not w3d.l4symtry:
                    getf3d_linear(n,x,y,z,ex,ey,ez,
                             f.xmin,f.ymin,f.zmin+self.zgrid,
                             f.dx,f.dy,f.dz,
                             f.nx,f.ny,f.nz,
                             f.nxguard,f.nyguard,f.nzguard,
                             f.Exp,f.Eyp,f.Ezp)
                else:
                    getf3d_n(n,x,y,z,ex,ey,ez,
                             f.xmin,f.ymin,f.zmin+self.zgrid,
                             f.dx,f.dy,f.dz,
                             f.nx,f.ny,f.nz,
                             f.nxguard,f.nyguard,f.nzguard,
                             nox,noy,noz,
                             f.Exp,f.Eyp,f.Ezp,
                             w3d.l4symtry)
        elif top.efetch[w3d.jsfsapi]==4:
            if self.l_1dz:
                gete1dz_n_energy_conserving(n,z,ex,ey,ez,
                              f.zmin+self.zgrid,
                              f.dz,
                              f.nz,
                              f.nzguard,
                              noz,
                              f.Exp,f.Eyp,f.Ezp,
                              self.l_lower_order_in_v)
            elif self.l_2dxz:
                gete2dxz_n_energy_conserving(n,x,y,z,ex,ey,ez,
                              f.xmin,f.zmin+self.zgrid,
                              f.dx,f.dz,
                              f.nx,f.nz,
                              f.nxguard,f.nzguard,
                              nox,noz,
                              f.Exp,f.Eyp,f.Ezp,
                              w3d.l4symtry,self.l_2drz,
                              self.l_lower_order_in_v)
            else:
                if 0:#nox==1 and noy==1 and noz==1 and not w3d.l4symtry:
                    gete3d_linear_energy_conserving(n,x,y,z,ex,ey,ez,
                                  f.xmin,f.ymin,f.zmin+self.zgrid,
                                  f.dx,f.dy,f.dz,
                                  f.nx,f.ny,f.nz,
                                  f.nxguard,f.nyguard,f.nzguard,
                                  f.Exp,f.Eyp,f.Ezp)
                else:
                    gete3d_n_energy_conserving(n,x,y,z,ex,ey,ez,
                                  f.xmin,f.ymin,f.zmin+self.zgrid,
                                  f.dx,f.dy,f.dz,
                                  f.nx,f.ny,f.nz,
                                  f.nxguard,f.nyguard,f.nzguard,
                                  nox,noy,noz,
                                  f.Exp,f.Eyp,f.Ezp,
                                  w3d.l4symtry,
                                  self.l_lower_order_in_v)
        # --- fetch b
        if top.efetch[w3d.jsfsapi] in [1,3,5]:
            if self.l_1dz:
                getf1dz_n(n,z,bx,by,bz,
                         f.zmin+self.zgrid,
                         f.dz,
                         f.nz,
                         f.nzguard,
                         noz,
                         f.Bxp,f.Byp,f.Bzp)
            elif self.l_2dxz :
                if f.circ_m == 0 :
                    getf2dxz_n(n,x,y,z,bx,by,bz,
                            f.xmin,f.zmin+self.zgrid,
                            f.dx,f.dz,
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            nox,noz,
                            f.Bxp,f.Byp,f.Bzp,
                            w3d.l4symtry,self.l_2drz)
                else:
                    getf2drz_circ_n(n,x,y,z,bx,by,bz,
                                f.xmin,f.zmin+self.zgrid,
                                f.dx,f.dz,
                                f.nx,f.ny,f.nz,
                                f.nxguard,f.nyguard,f.nzguard,
                                nox,noz,
                                f.Bxp,f.Byp,f.Bzp,
                                f.Bxp_circ,f.Byp_circ,f.Bzp_circ,
                                f.circ_m)
            else:
                getf3d_n(n,x,y,z,bx,by,bz,
                         f.xmin,f.ymin,f.zmin+self.zgrid,
                         f.dx,f.dy,f.dz,
                         f.nx,f.ny,f.nz,
                         f.nxguard,f.nyguard,f.nzguard,
                         nox,noy,noz,
                         f.Bxp,f.Byp,f.Bzp,
                         w3d.l4symtry)
        elif top.efetch[w3d.jsfsapi]==4:
            if self.l_1dz:
                getb1dz_n_energy_conserving(n,z,bx,by,bz,
                            f.zmin+self.zgrid,
                            f.dz,
                            f.nz,
                            f.nzguard,
                            noz,
                            f.Bxp,f.Byp,f.Bzp,
                            self.l_lower_order_in_v)
            elif self.l_2dxz:
                getb2dxz_n_energy_conserving(n,x,y,z,bx,by,bz,
                            f.xmin,f.zmin+self.zgrid,
                            f.dx,f.dz,
                            f.nx,f.nz,
                            f.nxguard,f.nzguard,
                            nox,noz,
                            f.Bxp,f.Byp,f.Bzp,
                            w3d.l4symtry,self.l_2drz,
                            self.l_lower_order_in_v)
            else:
                if 0:#nox==1 and noy==1 and noz==1 and not w3d.l4symtry:
                    getb3d_linear_energy_conserving(n,x,y,z,bx,by,bz,
                                f.xmin,f.ymin,f.zmin+self.zgrid,
                                f.dx,f.dy,f.dz,
                                f.nx,f.ny,f.nz,
                                f.nxguard,f.nyguard,f.nzguard,
                                f.Bxp,f.Byp,f.Bzp)
                else:
                    getb3d_n_energy_conserving(n,x,y,z,bx,by,bz,
                                f.xmin,f.ymin,f.zmin+self.zgrid,
                                f.dx,f.dy,f.dz,
                                f.nx,f.ny,f.nz,
                                f.nxguard,f.nyguard,f.nzguard,
                                nox,noy,noz,
                                f.Bxp,f.Byp,f.Bzp,
                                w3d.l4symtry,
                                self.l_lower_order_in_v)

    def fetchphifrompositions(self,x,z,phi):
        pass

#===============================================================================
    def getfieldsfrompositions(self,x,y,z,gather=True):
#===============================================================================
        # --- This returns e and b at positions x,y,z
        n = len(x)
        if n == 0: return None
        f = self.block.core.yf
        if self.l_1dz:
            ilist = compress((z>=(f.zmin+self.zgrid)) & (z<(f.zmax+self.zgrid)),arange(n))
        elif self.l_2dxz:
            ilist = compress((x>=f.xmin) & (x<f.xmax) &
                             (z>=(f.zmin+self.zgrid)) & (z<(f.zmax+self.zgrid)),arange(n))
        else:
            ilist = compress((x>=f.xmin) & (x<f.xmax) &
                             (y>=f.ymin) & (y<f.ymax) &
                             (z>=(f.zmin+self.zgrid)) & (z<(f.zmax+self.zgrid)),arange(n))
        nlocal = len(ilist)
        if nlocal>0:
            x = take(x,ilist)
            y = take(y,ilist)
            z = take(z,ilist)
            ex=zeros(nlocal,'d')
            ey=zeros(nlocal,'d')
            ez=zeros(nlocal,'d')
            bx=zeros(nlocal,'d')
            by=zeros(nlocal,'d')
            bz=zeros(nlocal,'d')
            self.fetchfieldfrompositions(x,y,z,ex,ey,ez,bx,by,bz)

        if not gather:
            if nlocal>0:
                return ilist,ex,ey,ez,bx,by,bz
            else:
                return None,None,None,None,None,None,None

        exp=zeros(n,'d')
        eyp=zeros(n,'d')
        ezp=zeros(n,'d')
        bxp=zeros(n,'d')
        byp=zeros(n,'d')
        bzp=zeros(n,'d')
        if me>0:
            mpisend(nlocal, dest = 0, tag = 3)
            if nlocal>0:
                mpisend((ilist,ex,ey,ez,bx,by,bz), dest = 0, tag = 3)
        else:
            if nlocal>0:
                for j,i in enumerate(ilist):
                    exp[i] = ex[j]
                    eyp[i] = ey[j]
                    ezp[i] = ez[j]
                    bxp[i] = bx[j]
                    byp[i] = by[j]
                    bzp[i] = bz[j]
            for ip in range(1,npes):
                nlocal = mpirecv(source = ip, tag = 3)
                if nlocal>0:
                    ilist,ex,ey,ez,bx,by,bz = mpirecv(source = ip, tag = 3)
                    for j,i in enumerate(ilist):
                        exp[i] = ex[j]
                        eyp[i] = ey[j]
                        ezp[i] = ez[j]
                        bxp[i] = bx[j]
                        byp[i] = by[j]
                        bzp[i] = bz[j]

        return exp,eyp,ezp,bxp,byp,bzp

################################################################################
# CHARGE/CURRENT DEPOSITION
################################################################################

#===============================================================================
    def setsourcep(self,js,pgroup,zgrid):
#===============================================================================
        if self.l_verbose:print 'setsourcep, species ',js
        if (pgroup.ldodepos[js]==False): return
        n  = pgroup.nps[js]
        if n == 0: return
        i  = pgroup.ins[js] - 1
        x  = pgroup.xp[i:i+n]
        y  = pgroup.yp[i:i+n]
        z  = pgroup.zp[i:i+n]
        ux = pgroup.uxp[i:i+n]
        uy = pgroup.uyp[i:i+n]
        uz = pgroup.uzp[i:i+n]
        gaminv = pgroup.gaminv[i:i+n]
        q  = pgroup.sq[js]
        w  = pgroup.sw[js]*pgroup.dtscale[js]
        if top.wpid==0:
            wfact = zeros((0,), 'd')
        else:
            wfact = pgroup.pid[i:i+n,top.wpid-1]
        w3d.jsfsapi=js
        self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,w)

    def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wfact,zgrid,q,w):
        n = x.shape[0]
        if n == 0: return
        # --- call routine performing current deposition
        f = self.block.core.yf
        js = w3d.jsfsapi
        nox = top.depos_order[0,js]
        noy = top.depos_order[1,js]
        noz = top.depos_order[2,js]
        dt = top.dt*top.pgroup.ndts[js]

        if top.wpid==0:
            wfact = ones((1,),'d')
            l_particles_weight = false
        else:
            l_particles_weight = true

        self.depose_current_density(n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight)

        if self.l_getrho :
          self.depose_charge_density(n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight)

    def depose_current_density(self,n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight):
        if not self.deposit_energy_density:
            if self.l_1dz:
                if 0:
                    jx = self.fields.Jx[self.fields.nxguard,self.fields.nyguard,:]*0.
                    jy = self.fields.Jy[self.fields.nxguard,self.fields.nyguard,:]*0.
                    jz = self.fields.Jz[self.fields.nxguard,self.fields.nyguard,:]*0.
                    depose_jxjyjz_esirkepov_linear_serial1d(jx,jy,jz,n,x*0.,y,z,ux,uy,uz,gaminv,
                                                      wfact*w,q/w3d.dx,
                                                      -0.5,-0.5,f.zmin+self.zgrid,
                                                      dt,
                                                      1.,1.,self.fields.dz,
                                                      self.fields.nz,
                                                      l_particles_weight)
#          for j in range(1,shape(self.fields.Jarray)[0]-1):
                    for j in range(shape(self.fields.Jxarray)[0]):
                        self.fields.Jxarray[j,self.fields.nyguard,:,0]+=jx
                        self.fields.Jyarray[j,self.fields.nyguard,:,0]+=jy
                        self.fields.Jzarray[j,self.fields.nyguard,:,0]+=jz
                else:
                    jx = self.fields.Jx[self.fields.nxguard,self.fields.nyguard,:]
                    jy = self.fields.Jy[self.fields.nxguard,self.fields.nyguard,:]
                    jz = self.fields.Jz[self.fields.nxguard,self.fields.nyguard,:]
                    depose_j_n_1dz(jx,jy,jz,n,z,ux,uy,uz,gaminv,
                                                      wfact,q*w,
                                                      f.zmin+self.zgrid,
                                                      dt,
                                                      self.fields.dz,
                                                      self.fields.nz,
                                                      f.nzguard,
                                                      noz,
                                                      l_particles_weight)
            elif self.l_2dxz:
                jx = self.fields.Jx[:,self.fields.nyguard,:]
                jy = self.fields.Jy[:,self.fields.nyguard,:]
                jz = self.fields.Jz[:,self.fields.nyguard,:]
                if self.fields.circ_m:
                    jx_circ = self.fields.Jx_circ
                    jy_circ = self.fields.Jy_circ
                    jz_circ = self.fields.Jz_circ
                if 0:
#          depose_jxjy_esirkepov_linear_serial_2d(j,n,z,x,z-gaminv*uz*top.dt,x-gaminv*ux*top.dt,
#                                                 uy,gaminv,
#                                            wfact,q*w,
#                                            f.zmin+self.zgrid,f.xmin,
#                                            dt,
#                                            f.dz,f.dx,
#                                            f.nz,f.nx,
#                                           l_particles_weight)
                    depose_jxjy_esirkepov_linear_serial_2d(jx,jy,jz,n,x,z,x-gaminv*ux*top.dt,z-gaminv*uz*top.dt,
                                                           uy,gaminv,
                                                      wfact,q*w,
                                                      f.xmin,f.zmin+self.zgrid,
                                                      dt,
                                                      f.dx,f.dz,
                                                      f.nx,f.nz,
                                                      l_particles_weight)
                else:
                    if self.l_esirkepov:
                        if self.circ_m==0:
                            depose_jxjyjz_esirkepov_n_2d(jx,jy,jz,n,
                                                x,y,z,ux,uy,uz,
                                                gaminv,wfact,q*w,
                                                f.xmin-0.5*top.dt*self.V_galilean[0],
                                                f.zmin+self.zgrid-0.5*top.dt*self.V_galilean[2],
                                                dt,
                                                f.dx,f.dz,
                                                f.nx,f.nz,
                                                f.nxguard,f.nzguard,
                                                nox,
                                                noz,
                                                l_particles_weight,
                                                w3d.l4symtry,self.l_2drz,
                                                self.type_rz_depose)
                        else:
                            depose_jxjyjz_esirkepov_n_2d_circ(jx,jy,jz,
                                                    jx_circ,jy_circ,jz_circ,
                                                    self.circ_m,n,
                                                    x,y,z,ux,uy,uz,
                                                    gaminv,wfact,q*w,
                                                    f.xmin-0.5*top.dt*self.V_galilean[0],
                                                    f.zmin+self.zgrid-0.5*top.dt*self.V_galilean[2],
                                                    dt,
                                                    f.dx,f.dz,
                                                    f.nx,f.nz,
                                                    f.nxguard,f.nzguard,
                                                    nox,
                                                    noz,
                                                    l_particles_weight,
                                                    self.type_rz_depose)
                    else:
                        depose_j_n_2dxz(jx,jy,jz,n,x,z,ux,uy,uz,
                                        gaminv,wfact,q*w,
                                        f.xmin-0.5*top.dt*self.V_galilean[0],
                                        f.zmin+self.zgrid-0.5*top.dt*self.V_galilean[2],
                                        dt,
                                        f.dx,f.dz,
                                        f.nx,f.nz,
                                        f.nxguard,f.nzguard,
                                        nox,
                                        noz,
                                        l_particles_weight,
                                        w3d.l4symtry,
                                        self.l_deposit_nodal,1,True)
            else:
                if 0:#nox==1 and noy==1 and noz==1 and not w3d.l4symtry:
                    depose_jxjyjz_esirkepov_linear_serial(self.fields.Jx,self.fields.Jy,self.fields.Jz,n,
                                                      x,y,z,ux,uy,uz,
                                                      gaminv,wfact,q*w,
                                                      f.xmin,f.ymin,f.zmin+self.zgrid,
                                                      dt,
                                                      f.dx,f.dy,f.dz,
                                                      f.nx,f.ny,f.nz,
                                                      f.nxguard,f.nyguard,f.nzguard,
                                                      l_particles_weight)
                else:
                    depose_jxjyjz_esirkepov_n(self.fields.Jx,self.fields.Jy,self.fields.Jz,n,
                                                      x,y,z,ux,uy,uz,
                                                      gaminv,wfact,q*w,
                                                      f.xmin-0.5*top.dt*self.V_galilean[0],
                                                      f.ymin-0.5*top.dt*self.V_galilean[1],
                                                      f.zmin+self.zgrid-0.5*top.dt*self.V_galilean[2],
                                                      dt,
                                                      f.dx,f.dy,f.dz,
                                                      f.nx,f.ny,f.nz,
                                                      f.nxguard,f.nyguard,f.nzguard,
                                                      nox,
                                                      noy,
                                                      noz,
                                                      l_particles_weight,w3d.l4symtry)
        else:
            depose_jxjyjz_pxpypz_esirkepov_linear_serial(self.fields.Jx,self.fields.Jy,self.fields.Jz,
                                              self.fields.Mp,n,
                                              x,y,z,ux,uy,uz,
                                              gaminv,wfact,q*w,m,
                                              f.xmin,f.ymin,f.zmin+self.zgrid,
                                              dt,
                                              f.dx,f.dy,f.dz,
                                              f.nx,f.ny,f.nz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              l_particles_weight,
                                              top.lrelativ)

    def depose_charge_density(self,n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight):

        if self.l_2dxz:
                if self.circ_m==0:
                    depose_rho_n_2dxz(self.fields.Rho,n,
                                         x,y,z,
                                         wfact,q*w,
                                         f.xmin,f.zmin+self.zgrid,
                                         f.dx,f.dz,
                                         f.nx,f.nz,
                                         f.nxguard,f.nzguard,
                                         nox,
                                         noz,
                                         l_particles_weight,w3d.l4symtry,
                                         self.l_2drz,
                                         self.type_rz_depose)
                else:
                    depose_rho_n_2d_circ(self.fields.Rho,
                                         self.fields.Rho_circ,
                                         self.fields.circ_m,
                                         n,x,y,z,
                                         wfact,q*w,
                                         f.xmin,f.zmin+self.zgrid,
                                         f.dx,f.dz,
                                         f.nx,f.nz,
                                         f.nxguard,f.nzguard,
                                         nox,
                                         noz,
                                         l_particles_weight,
                                         self.type_rz_depose)
        else:
                depose_rho_n(self.fields.Rho,n,
                                       x,y,z,
                                       wfact,q*w,
                                       f.xmin,f.ymin,f.zmin+self.zgrid,
                                       f.dx,f.dy,f.dz,
                                       f.nx,f.ny,f.nz,
                                       f.nxguard,f.nyguard,f.nzguard,
                                       nox,
                                       noy,
                                       noz,
                                       l_particles_weight,w3d.l4symtry)

    def allocatedataarrays(self):
        if self.l_verbose:print 'allocatedataarrays'
        self.finalize()
        # --- reallocate Jarray if needed
        if self.fields.ntimes != top.nsndts:
            self.fields.ntimes=top.nsndts
            self.fields.gchange()

    def zerosourcep(self):
        if self.l_verbose:print 'zerosourcep',self

        # --- copy rho to rhoold if needed
        if self.l_getrho:
            self.fields.Rhoold = self.fields.Rho.copy()
            if self.refinement is not None:
                self.field_coarse.fields.Rhoold = self.field_coarse.fields.Rho.copy()

        # --- zero proper portion of Jarray
        for indts in range(top.nsndts):
            if top.ldts[indts]:
                self.fields.Jxarray[...,indts] = 0.
                self.fields.Jyarray[...,indts] = 0.
                self.fields.Jzarray[...,indts] = 0.
                if self.fields.circ_m:
                    self.fields.Jx_circ[...] = 0.
                    self.fields.Jy_circ[...] = 0.
                    self.fields.Jz_circ[...] = 0.
#        if self.refinement is not None:
#          self.field_coarse.fields.Jarray[...,indts] = 0.
                if self.l_getrho:
                    self.fields.Rhoarray[...,indts] = 0.
                    if self.fields.circ_m:
                        self.fields.Rho_circ[...] = 0.
#          if self.refinement is not None:
#            self.field_coarse.fields.Rho[...,indts] = 0.
                if self.deposit_energy_density:
                    self.fields.Mp[...] = 0.

    def setsourcepforparticles(self,isourcepndtscopies,indts,iselfb):
        if self.l_verbose:print 'setsourcepforparticles'
        # --- point J array to proper Jarray slice
        self.fields.Jx = self.fields.Jxarray[:,:,:,indts]
        self.fields.Jy = self.fields.Jyarray[:,:,:,indts]
        self.fields.Jz = self.fields.Jzarray[:,:,:,indts]
        if self.l_getrho: self.fields.Rho = self.fields.Rhoarray[:,:,:,indts]

    def add_source_ndts_slices(self):
        # --- add slices
        if top.nsndts>1:
            if self.refinement is not None:raise Exception('Error in finalizesourcep:nsndts>1 not fully implemented yet with MR')
            for indts in range(top.nsndts-2,-1,-1):
                if top.ldts[indts]:
                    add_current_slice_3d(self.fields,indts+1)
                    if self.l_getrho:add_rho_slice_3d(self.fields,indts+1)

    def finalizesourcep(self):
        if self.sourcepfinalized: return
        self.sourcepfinalized = True
        if self.l_verbose:print 'finalizesourcep'
        # --- add slices
        self.add_source_ndts_slices()
        self.aftersetsourcep()
        # -- add laser, loop on the different antennas
        for i in range(len(self.laser_antenna)):
            self.add_laser(self.block.core.yf, self.laser_antenna[i])
        # --- smooth current density
        if any(self.npass_smooth>0):self.smoothdensity()
        if self.l_nodalgrid and not self.l_deposit_nodal:self.Jyee2node3d()
        # --- apply boundary conditions
        self.applysourceboundaryconditions()
        # --- exchange guard cells data across domains
        self.exchange_bc(self.block)
        if self.l_verbose:print 'finalizesourcep done'

    def Jyee2node3d(self):
        Jyee2node3d(self.block.core.yf)
        if self.refinement is not None:
            self.__class__.__bases__[1].Jyee2node3d(self.field_coarse)

    def apply_rho_bc(self,block):
        # --- point Rho to first slice of Rhoarray
        block.core.yf.Rho = block.core.yf.Rhoarray[:,:,:,0]
        em3d_applybc_rho(block.core.yf,
                         block.xlbnd,
                         block.xrbnd,
                         block.ylbnd,
                         block.yrbnd,
                         block.zlbnd,
                         block.zrbnd,
                         self.type_rz_depose)

    def apply_current_bc(self,block):
        # --- point J to first slice of Jarray
        block.core.yf.Jx = block.core.yf.Jxarray[:,:,:,0]
        block.core.yf.Jy = block.core.yf.Jyarray[:,:,:,0]
        block.core.yf.Jz = block.core.yf.Jzarray[:,:,:,0]
        em3d_applybc_j(block.core.yf,
                       block.xlbnd,
                       block.xrbnd,
                       block.ylbnd,
                       block.yrbnd,
                       block.zlbnd,
                       block.zrbnd,
                       self.type_rz_depose)

    def exchange_bc(self,block):
        em3d_exchange_j(block)
        if self.l_getrho:
            em3d_exchange_rho(block)
        if self.refinement is not None:
            self.exchange_bc(self.field_coarse.block)

    def applysourceboundaryconditions(self):
        # --- apply boundary condition on current
        self.apply_current_bc(self.block)
        if self.refinement is not None:
            self.apply_current_bc(self.field_coarse.block)
        if self.l_getrho:
            self.apply_rho_bc(self.block)
            if self.refinement is not None:
                self.apply_rho_bc(self.field_coarse.block)

        if self.l_sumjx:
            j1 = self.fields.Jxarray[0,:,:,0]*0.
            j2 = self.fields.Jyarray[0,:,:,0]*0.
            j3 = self.fields.Jzarray[0,:,:,0]*0.
            for i in range(shape(self.fields.Jxarray)[0]):
                j1+=self.fields.Jxarray[i,:,:,0]
                j2+=self.fields.Jyarray[i,:,:,0]
                j3+=self.fields.Jzarray[i,:,:,0]
            for i in range(shape(self.fields.Jxarray)[0]):
                self.fields.Jxarray[i,:,:,0]=j1.copy()
                self.fields.Jyarray[i,:,:,0]=j2.copy()
                self.fields.Jzarray[i,:,:,0]=j3.copy()

        # --- point J to first slice of Jarray
        self.fields.Jx = self.fields.Jxarray[:,:,:,0]
        self.fields.Jy = self.fields.Jyarray[:,:,:,0]
        self.fields.Jz = self.fields.Jzarray[:,:,:,0]
        if self.l_getrho:self.fields.Rho = self.fields.Rhoarray[:,:,:,0]

    def smoothdensity(self):
        if all(self.npass_smooth==0):return
        nx,ny,nz = shape(self.fields.Jx)
        if self.circ_m > 0 :
            nxc, nzc, circ_m = shape(self.fields.Jx_circ[:,:,:])
        nsm = shape(self.npass_smooth)[1]
        l_mask_method=1
        if self.mask_smooth is None:
            self.mask_smooth=[]
            for js in range(nsm):
                self.mask_smooth.append(None)
        if l_mask_method==2:
            Jxcopy = self.fields.Jx.copy()
            Jycopy = self.fields.Jy.copy()
            Jzcopy = self.fields.Jz.copy()
            if self.l_getrho:
                Rhocopy = self.fields.Rho.copy()
        for js in range(nsm):
            if self.mask_smooth[js] is None or l_mask_method==2:
                smooth3d_121_stride( self.fields.Jx,
                    nx-1, ny-1, nz-1, self.npass_smooth[:,js].copy(),
                    self.alpha_smooth[:,js].copy(), self.stride_smooth[:,js].copy())
                smooth3d_121_stride( self.fields.Jy,
                    nx-1, ny-1, nz-1, self.npass_smooth[:,js].copy(),
                    self.alpha_smooth[:,js].copy(), self.stride_smooth[:,js].copy())
                smooth3d_121_stride( self.fields.Jz,
                    nx-1, ny-1, nz-1, self.npass_smooth[:,js].copy(),
                    self.alpha_smooth[:,js].copy(),self.stride_smooth[:,js].copy())
                if self.circ_m > 0 :
                    smoothcirc_121_stride( self.fields.Jx_circ,
                        nxc-1, nzc-1, circ_m-1, self.npass_smooth[:,js].copy(),
                        self.alpha_smooth[:,js].copy(), self.stride_smooth[:,js].copy())
                    smoothcirc_121_stride( self.fields.Jy_circ,
                        nxc-1, nzc-1, circ_m-1, self.npass_smooth[:,js].copy(),
                        self.alpha_smooth[:,js].copy(), self.stride_smooth[:,js].copy())
                    smoothcirc_121_stride( self.fields.Jz_circ,
                        nxc-1, nzc-1, circ_m-1, self.npass_smooth[:,js].copy(),
                        self.alpha_smooth[:,js].copy(), self.stride_smooth[:,js].copy())
                if self.l_getrho:
                    smooth3d_121_stride( self.fields.Rho[...],
                        nx-1, ny-1, nz-1, self.npass_smooth[:,js].copy(),
                        self.alpha_smooth[:,js].copy(),self.stride_smooth[:,js].copy())
                    if self.circ_m > 0 :
                        smoothcirc_121_stride( self.fields.Rho_circ[...],
                        nxc-1, nzc-1, circ_m-1, self.npass_smooth[:,js].copy(),
                        self.alpha_smooth[:,js].copy(), self.stride_smooth[:,js].copy())
            else:
                smooth3d_121_stride_mask( self.fields.Jx,
                    self.mask_smooth[js], nx-1, ny-1, nz-1,
                    self.npass_smooth[:,js].copy(), self.alpha_smooth[:,js].copy(),
                    self.stride_smooth[:,js].copy())
                smooth3d_121_stride_mask( self.fields.Jy,
                    self.mask_smooth[js], nx-1, ny-1, nz-1,
                    self.npass_smooth[:,js].copy(), self.alpha_smooth[:,js].copy(),
                    self.stride_smooth[:,js].copy())
                smooth3d_121_stride_mask( self.fields.Jz,
                    self.mask_smooth[js], nx-1, ny-1, nz-1,
                    self.npass_smooth[:,js].copy(), self.alpha_smooth[:,js].copy(),
                    self.stride_smooth[:,js].copy())
                if self.l_getrho:
                    smooth3d_121_stride_mask(self.fields.Rho[...],
                        self.mask_smooth[js], nx-1, ny-1, nz-1,
                        self.npass_smooth[:,js].copy(), self.alpha_smooth[:,js].copy(),
                        self.stride_smooth[:,js].copy())
        if l_mask_method==2:
            self.fields.Jx *= self.mask_smooth
            self.fields.Jx += Jxcopy*(1.-self.mask_smooth)
            self.fields.Jy *= self.mask_smooth
            self.fields.Jy += Jycopy*(1.-self.mask_smooth)
            self.fields.Jz *= self.mask_smooth
            self.fields.Jz += Jzcopy*(1.-self.mask_smooth)
            if self.l_getrho:
                self.fields.Rho *= self.mask_smooth
                self.fields.Rho += Rhocopy*(1.-self.mask_smooth)

    def smootharray(self,f,js=None,mask=None):
        nx,ny,nz = shape(f)
        if js is None:
            jslist=arange(shape(self.npass_smooth)[1])
        else:
            jslist=[js]
        for js in jslist:
            if mask is None:
                smooth3d_121_stride(f,nx-1,ny-1,nz-1,
                                    self.npass_smooth[:,js].copy(),
                                    self.alpha_smooth[:,js].copy(),
                                    self.stride_smooth[:,js].copy())
            else:
                smooth3d_121_stride_mask(f,self.mask_smooth[js],nx-1,ny-1,nz-1,
                                         self.npass_smooth[:,js].copy(),
                                         self.alpha_smooth[:,js].copy(),
                                         self.stride_smooth[:,js].copy())

    def smoothfields(self):
        nx,ny,nz = shape(self.fields.Jx)
        nsm = shape(self.npass_smooth)[1]
        l_mask_method=1
        if self.mask_smooth is None:
            self.mask_smooth=[]
            for js in range(nsm):
                self.mask_smooth.append(None)
        if l_mask_method==2:
            Expcopy = self.fields.Exp.copy()
            Eypcopy = self.fields.Eyp.copy()
            Ezpcopy = self.fields.Ezp.copy()
            Bxpcopy = self.fields.Bxp.copy()
            Bypcopy = self.fields.Byp.copy()
            Bzpcopy = self.fields.Bzp.copy()
        for js in range(nsm):
            if self.mask_smooth[js] is None or l_mask_method==2:
                mask=None
            else:
                mask=self.mask_smooth[js]
            self.smootharray(self.fields.Exp,js,mask)
            self.smootharray(self.fields.Eyp,js,mask)
            self.smootharray(self.fields.Ezp,js,mask)
            self.smootharray(self.fields.Bxp,js,mask)
            self.smootharray(self.fields.Byp,js,mask)
            self.smootharray(self.fields.Bzp,js,mask)
        if l_mask_method==2:
            self.fields.Exp *= self.mask_smooth; self.fields.Exp += Expcopy*(1.-self.mask_smooth)
            self.fields.Eyp *= self.mask_smooth; self.fields.Eyp += Eypcopy*(1.-self.mask_smooth)
            self.fields.Ezp *= self.mask_smooth; self.fields.Ezp += Ezpcopy*(1.-self.mask_smooth)
            self.fields.Bxp *= self.mask_smooth; self.fields.Bxp += Bxpcopy*(1.-self.mask_smooth)
            self.fields.Byp *= self.mask_smooth; self.fields.Byp += Bypcopy*(1.-self.mask_smooth)
            self.fields.Bzp *= self.mask_smooth; self.fields.Bzp += Bzpcopy*(1.-self.mask_smooth)
        if self.n_smooth_fields is not None:
            if top.it%self.n_smooth_fields==0:
                nsm = shape(self.npass_smooth)[1]
                for js in range(nsm):
                    self.smootharray(self.fields.Ex,js)
                    self.smootharray(self.fields.Ey,js)
                    self.smootharray(self.fields.Ez,js)
                    self.smootharray(self.fields.Bx,js)
                    self.smootharray(self.fields.By,js)
                    self.smootharray(self.fields.Bz,js)
                    if self.l_pushf:
                        smooth3d_121(self.fields.F,nx-1,ny-1,nz-1,npass_smooth[:,js]*m,alpha_smooth[:,js].copy())

    def smoothfields_poly(self):
        nx,ny,nz = shape(self.fields.Exp)
        smooth_poly(self.fields.Exp,nx-1,ny-1,nz-1,self.fields.ex_stencil)
        smooth_poly(self.fields.Eyp,nx-1,ny-1,nz-1,self.fields.ex_stencil)
        smooth_poly(self.fields.Ezp,nx-1,ny-1,nz-1,self.fields.by_stencil)
        smooth_poly(self.fields.Bxp,nx-1,ny-1,nz-1,self.fields.by_stencil)
        smooth_poly(self.fields.Byp,nx-1,ny-1,nz-1,self.fields.by_stencil)
        smooth_poly(self.fields.Bzp,nx-1,ny-1,nz-1,self.fields.ex_stencil)

    def getsmoothx(self):
        nx = w3d.nx
        ny = 1
        nz = 1
        a = zeros([nx,ny,nz],'d')
        a[nx/2,:,:]=1.
        nsm = shape(self.npass_smooth)[1]
        for js in range(nsm):
            self.smootharray(a,js)
        return a[:,0,0]

    def getsmoothz(self):
        nz = w3d.nz
        ny = 1
        nx = 1
        a = zeros([nx,ny,nz],'d')
        a[:,:,nz/2]=1.
        nsm = shape(self.npass_smooth)[1]
        for js in range(nsm):
            self.smootharray(a,js)
        return a[0,0,:]

    def getsmoothxfftth(self):
        nx = w3d.nx
        theta=(1.+arange(nx))*pi/nx
        nsm = shape(self.npass_smooth)[1]
        f = ones(nx,'d')
        print nx,shape(theta),shape(f)
        for js in range(nsm):
            w = 0.5*(1./self.alpha_smooth[0,js]-1.)
            f *= ((1.+2.*w*cos(theta*self.stride_smooth[0,js]))/(1.+2.*w))**self.npass_smooth[0,js]
        return f,theta

    def pltsmx(self,color=black,width=1.):
        pla(self.getsmoothx(),color=color,width=width)

    def pltsmxfft(self,color=black,width=1.):
        f=abs(fft.fft(self.getsmoothx()))
        theta=(1.+arange(shape(f)[0]))*2.*pi/shape(f)[0]
        pla(f,2.*pi/theta,color=color,width=width)

    def pltsmxfftth(self,color=black,width=1.):
        f,theta = self.getsmoothxfftth()
        pla(f,2.*pi/theta,color=color,width=width)

    def pltsmz(self,color=black,width=1.):
        pla(self.getsmoothz(),color=color,width=width)

    def pltsmzfft(self,color=black,width=1.):
        f=abs(fft.fft(self.getsmoothz()))
        theta=arange(shape(f)[0])*2.*pi/shape(f)[0]
        pla(f,theta,color=color,width=width)

    def binomial_expansion(self,x,y,n):
        # --- we assume that x is real and y is an array of coefficients
        # initializes f
        f = 0.*y
        # initializes Pascal triangle
        pt = zeros(n)
        pt[0] = 1.
        for i in range(n):
            f += pt[i]*x**(n-i)*self.getpower(y,i)
            if i<n-1:pt[1:] = pt[1:]+pt[:-1]
        return f

    def getpower(self,y,n):
        # computes coefficients of y**n where y is a list of coefficients
        o = shape(y)[0]
        f = y.copy()
        for k in range(n):
            x = f.copy()
            for i in range(o):
                for j in range(o):
                    if (i+j)<o:
                        f[i+j] += y[i]*x[j]
        return f

    def fetche(self,*args,**kw):
        if self.l_verbose:print 'fetche',self
        self.fetchfield(*args,**kw)

    def loadrho(self,lzero=None,lfinalize_rho=None,pgroups=None,**kw):
        if self.l_verbose:print 'loadrho',self
        self.loadsource(lzero,lfinalize_rho,pgroups,**kw)

    def fetchphi(self):
        pass

    def fetcha(self):
        pass

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

    def installconductor(self,conductor,
                              xmin=None,xmax=None,
                              ymin=None,ymax=None,
                              zmin=None,zmax=None,
                              dfill=None):
        # --- This only adds the conductor to the list. The data is only actually
        # --- installed when it is needed, during a call to getconductorobject.
        self.conductordatalist.append((conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill))

    def init_macroscopic_coefs(self):
        if self.fields.l_macroscopic:return

        self.fields.nxs = self.fields.nx
        self.fields.nys = self.fields.ny
        self.fields.nzs = self.fields.nz
        self.fields.gchange()
        self.fields.Sigmax=0.
        self.fields.Sigmay=0.
        self.fields.Sigmaz=0.
        self.fields.Epsix=1.
        self.fields.Epsiy=1.
        self.fields.Epsiz=1.
        self.fields.Mux=1.
        self.fields.Muy=1.
        self.fields.Muz=1.
        self.fields.l_macroscopic=True

    def _installconductor(self,conductorobject,installedlist,conductordata,fselfb):
        # --- This does that actual installation of the conductor into the
        # --- conductor object

        # --- Extract the data from conductordata (the arguments to installconductor)
        conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill = conductordata

        # --- Set dfill to be a large number so that the entire interior of the conductor
        # --- gets filled in. This ensures that the field is forced to zero everywhere
        # --- inside the conductor, but does not introduce a performance penalty.
        if dfill is None: dfill = largepos

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
        # to be updated
        mgmaxlevels=1
        conductorobject.interior.n=0
        installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                          self.zgrid,
                          nx,ny,nz,
                          nxlocal,nylocal,nzlocal,
                          xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,
                          zscale,self.l2symtry,self.l4symtry,
                          installrz=0,
                          solvergeom=self.solvergeom,conductors=conductorobject,
                          mgmaxlevels=mgmaxlevels,decomp=decomp)

        self.fields.nconds = conductorobject.interior.n
        self.fields.nxcond = self.fields.nx
        self.fields.nycond = self.fields.ny
        self.fields.nzcond = self.fields.nz
        self.fields.gchange()
        self.fields.incond=False
        if self.fields.nconds>0:
            if conductor.conductivity is not None or \
               conductor.permittivity is not None or \
               conductor.permeability is not None:
                self.init_macroscopic_coefs()
                if conductor.conductivity is not None:
                    conductivity = conductor.conductivity
                else:
                    conductivity = 0.
                if conductor.permittivity is not None:
                    permittivity = conductor.permittivity
                else:
                    permittivity = 1.
                if conductor.permeability is not None:
                    permeability = conductor.permeability
                else:
                    permeability = 1.
                set_macroscopic_coefs_on_yee(self.fields, \
                                             self.fields.nconds, \
                                             aint(conductorobject.interior.indx[:,:self.fields.nconds]), \
                                             conductivity, \
                                             permittivity, \
                                             permeability)
            else:
                set_incond(self.fields, \
                           self.fields.nconds, \
                           aint(conductorobject.interior.indx[:,:self.fields.nconds]))
                if self.block.xlbnd==openbc:self.fields.incond[:3,:,:]=False
                if self.block.xrbnd==openbc:self.fields.incond[-3:,:,:]=False
                if self.block.ylbnd==openbc:self.fields.incond[:,:3,:]=False
                if self.block.yrbnd==openbc:self.fields.incond[:,-3:,:]=False
                if self.block.zlbnd==openbc:self.fields.incond[:,:,:3]=False
                if self.block.zrbnd==openbc:self.fields.incond[:,:,-3:]=False

    def hasconductors(self):
        return len(self.conductordatalist) > 0

    def clearconductors(self):
        "Clear out the conductor data"
        for fselfb in top.fselfb:
            if fselfb in self.conductorobjects:
                conductorobject = self.conductorobjects[fselfb]
                conductorobject.interior.n = 0
                conductorobject.evensubgrid.n = 0
                conductorobject.oddsubgrid.n = 0
                self.installedconductorlists[fselfb] = []

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
                if (conductorobject.leveliz[0] != self.izpslave[self.my_index] or
                    conductorobject.levelnz[0] != self.nzpslave[self.my_index]):
                    del self.conductorobjects['p']
                    del self.installedconductorlists['p']
            except KeyError:
                # --- 'p' object has not yet been created anyway, so do nothing.
                pass

        conductorobject = self.conductorobjects.setdefault(fselfb,ConductorType())
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
        return
        'calls setconductorvoltage'
        # --- Loop over all of the selfb groups to that all conductor objects
        # --- are handled.
        for iselfb in range(top.nsselfb):
            conductorobject = self.getconductorobject(top.fselfb[iselfb])
            setconductorvoltage(voltage,condid,discrete,setvinject,
                                conductors=conductorobject)

    def optimizeconvergence(self,resetpasses=1):
        pass

    def move_window_fields(self):
        # --- move window in x
        self.x_gridcont+=self.V_galilean[0]*top.dt+self.vxgrid*top.dt
        self.push_galilean('x')
        while (abs(self.x_grid-self.x_gridcont)>=0.5*self.dx):
            nxmove = int(sign(self.vxgrid))
            self.move_cells(nxmove,'x')

        # --- move window in y
        self.y_gridcont+=self.V_galilean[1]*top.dt+self.vygrid*top.dt
        self.push_galilean('y')
        while (abs(self.y_grid-self.y_gridcont)>=0.5*self.dy):
            nymove = int(sign(self.vygrid))
            self.move_cells(nymove,'y')

        # --- move window in z with vzgrid
        self.z_gridcont+=self.V_galilean[2]*top.dt+self.vzgrid*top.dt
        self.push_galilean('z')
        while (abs(self.z_grid-self.z_gridcont)>=0.5*self.dz):
            nzmove = int(sign(self.vzgrid))
            self.move_cells(nzmove,'z')

        # --- move window in z with zgrid
        if top.vbeamfrm==0.:
            top.zgridprv=top.zgrid
            top.zbeam=top.zgrid
            top.zgridndts[...]=top.zgrid
            self.zgrid=top.zgrid
        if top.vbeamfrm > 0.:
            while ((top.zgrid-self.zgrid)>=0.5*self.dz):
                self.shift_cells_z(1)
        elif top.vbeamfrm < 0.:
            while ((top.zgrid-self.zgrid)<=-0.5*self.dz):
                self.shift_cells_z(-1)

    def incrementposition(self,z,dz,n):
        # --- z is incremented using an integer to avoid round off problems.
        # --- In the old way, repeatedly doing z += dz, z will drift
        # --- because of roundoff so that z != nz*dz. Incrementing using
        # --- integers fixes this problem.
        nz = int(z/dz)
        wz = z - nz*dz
        return (nz + n)*dz + wz

    def push_galilean(self, coord):
        # move the boundaries of the box along the coord axis
        # in case of galilean frame along the coordinate coord
        #coord = 'x', 'y', 'z'
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

        for (coord_object,coord_attribute) in listtoshift:
            # loop equivalent to coord_object.coord_attribute+=self.V_galilean[..]*top.dt
            # for each tupple in listtoshift
            coordtoshift =getattr(coord_object,coord_attribute)
            if   coord=='x': coordtoshift += self.V_galilean[0]*top.dt
            elif coord=='y': coordtoshift += self.V_galilean[1]*top.dt
            elif coord=='z': coordtoshift += self.V_galilean[2]*top.dt
            setattr(coord_object,coord_attribute,coordtoshift)

    def move_cells(self,n, coord):
        # move the boundaries of the box along the coord axis
        # in case of moving window along the coordinate coord.
        #coord = 'x', 'y', 'z'
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

    def shift_cells_z(self,n):
        shift_em3dblock_ncells_z(self.block,n)
        self.zgrid = self.incrementposition(self.zgrid,self.dz,n)

    def solve2ndhalf(self):
        self.allocatedataarrays()
        if self.solveroff:return
        if self.fields.spectral:
            self.move_window_fields()
            return
        if self.mode==2:
            self.solve2ndhalfmode2()
            return
        if self.l_verbose:print 'solve 2nd half',self
        if top.dt != self.dtinit:raise Exception('Time step has been changed since initialization of EM3D.')
        self.push_b_part_2()
        if self.l_pushf:self.exchange_f()
        self.exchange_b()
        self.move_window_fields()
        if self.l_verbose:print 'solve 2nd half done'

    def dosolve(self,iwhich=0,*args):
        self.getconductorobject()
        if self.solveroff:return
        if self.ntsub<1.:
            self.novercycle = nint(1./self.ntsub)
            self.icycle = (top.it-1)%self.novercycle
        else:
            self.novercycle = 1
            self.icycle = 0
        if self.mode==2:
            self.dosolvemode2()
            return
        if any(top.fselfb != 0.):raise Exception('Error:EM solver does not work if fselfb != 0.')
        if self.l_verbose:print 'solve 1st half'
        if top.dt != self.dtinit:raise Exception('Time step has been changed since initialization of EM3D.')
        if self.fields.spectral:
            #self.move_window_fields()
            self.push_spectral_psaotd()
        else:
            self.push_e()
            self.exchange_e()
            for i in range(int(self.ntsub)-1):
                self.push_b_full()
                if self.l_pushf:self.exchange_f()
                self.exchange_b()
                self.push_e_full(i)
                self.exchange_e()
            self.push_b_part_1()
        if(hasattr(self,"full_pxr") == False):
          if self.pml_method==2:
              scale_em3d_bnd_fields(self.block,top.dt,self.l_pushf,self.l_pushg)
        elif(self.full_pxr == False):
          if self.pml_method==2:
              scale_em3d_bnd_fields(self.block,top.dt,self.l_pushf,self.l_pushg)

        if self.fields.spectral:
            self.exchange_e()
        if self.l_pushf:self.exchange_f()
        self.exchange_b()
        self.setebp()
        if not self.l_nodalgrid and top.efetch[0] != 4:
            if self.l_fieldcenterK:
                self.fieldcenterK()
            else:
                self.yee2node3d()
        if self.l_nodalgrid and top.efetch[0] == 4:
            self.fields.l_nodecentered=True
            self.node2yee3d()
        if self.l_smooth_particle_fields and any(self.npass_smooth>0):
            self.smoothfields()
        if self.l_correct_num_Cherenkov:self.smoothfields_poly()
        # --- for fields that are overcycled, they need to be pushed backward every ntsub
        self.push_e(dir=-1)
        self.exchange_e(dir=-1)
        if self.l_verbose:print 'solve 1st half done'

    def push_e(self,dir=1.,l_half=False):
        dt = dir*top.dt/self.ntsub
        if l_half:dt*=0.5
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
            push_em3d_eef(self.block,dt,0,self.l_pushf,self.l_pushpot,1)
        if self.refinement is not None:
            self.__class__.__bases__[1].push_e(self.field_coarse,dir)

    def setebp(self):
#    if (self.refinement is not None) or \
#       (self.novercycle>1) or \
#       (self.l_smooth_particle_fields and any(self.npass_smooth>0)):
        setebp(self.fields,self.icycle,self.novercycle)
        if self.refinement is not None:
            self.__class__.__bases__[1].setebp(self.field_coarse)

    def exchange_e(self,dir=1.):
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
            em3d_exchange_e(self.block)
        if self.refinement is not None:
            self.__class__.__bases__[1].exchange_e(self.field_coarse)

    def push_b_part_1(self,dir=1.):
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
            push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot,True)
        if self.refinement is not None:
            self.__class__.__bases__[1].push_b_part_1(self.field_coarse,dir)

    def push_b_part_2(self):
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
            push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot,True)
        if self.refinement is not None:
            self.__class__.__bases__[1].push_b_part_2(self.field_coarse)

    def yee2node3d(self):
        if self.l_nodalgrid:return
        yee2node3d(self.block.core.yf)
        if self.refinement is not None:
            self.__class__.__bases__[1].yee2node3d(self.field_coarse)

    def node2yee3d(self):
#        if self.l_nodalgrid:return
        node2yee3d(self.block.core.yf)
        if self.refinement is not None:
            self.__class__.__bases__[1].node2yee3d(self.field_coarse)

    def exchange_b(self,dir=1.):
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
            em3d_exchange_b(self.block)
        if self.refinement is not None:
            self.__class__.__bases__[1].exchange_b(self.field_coarse,dir)

    def exchange_f(self,dir=1.):
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
            em3d_exchange_f(self.block)
        if self.refinement is not None:
            self.__class__.__bases__[1].exchange_f(self.field_coarse,dir)

    def push_b_full(self):
        dt = top.dt/self.ntsub
        if self.l_verbose:print 'push_b full',self,dt,top.it,self.icycle
        push_em3d_bf(self.block,dt,0,self.l_pushf,self.l_pushpot,True)

    def push_e_full(self,i):
        dt = top.dt/self.ntsub
        if self.l_getrho:
            w = float(i+2)/self.ntsub
            self.fields.Rho = (1.-w)*self.fields.Rhoold + w*self.fields.Rhoarray[...,0]
        if self.l_verbose:print 'push_e full',self,dt,top.it,self.icycle
        push_em3d_eef(self.block,dt,0,self.l_pushf,self.l_pushpot,True)

    ##########################################################################
    # Define the basic plot commands
    def getdataatslice(self,data,direction=None,slice=None,l_abs=False):
        """This does the same thing that genericpfem3d does when getting a slice of the
           data, and returns that slice."""
        if direction is None and not l_opyndx:direction=2
        if self.l_2dxz:direction=1
        dataslice = None
        if self.isactive:
            f=self.block.core.yf
            if self.l_1dz:
                nxd = nyd = 1
                nzd=shape(data)[0]
            elif self.l_2dxz:
                nyd = 1
                nxd,nzd=shape(data)
            else:
                nxd,nyd,nzd=shape(data)
            if direction==0:
                if slice is None:
                    if self.l4symtry:
                        slice=0
                    else:
                        slice=self.nx/2
                xslice = w3d.xmmin+slice*w3d.dx
                selfslice = nint((xslice-self.block.xmin)/self.block.dx)
                if selfslice<0 or selfslice>nxd-1:
                    dataslice=None
                else:
                    dataslice=data[selfslice,:,:]
            if direction==1:
                if slice is None:
                    if self.l4symtry:
                        slice=0
                    else:
                        slice=self.ny/2
                yslice = w3d.ymmin+slice*w3d.dy
                selfslice = nint((yslice-self.block.ymin)/self.block.dy)
                if self.l_2dxz:slice=selfslice=0
                if selfslice<0 or selfslice>nyd-1:
                    dataslice=None
                else:
                    if not self.l_2dxz:
                        dataslice=data[:,selfslice,:]
                    else:
                        dataslice = data
            if direction==2:
                if slice is None:slice=self.nz/2
                zslice = w3d.zmmin+slice*w3d.dz
                selfslice = nint((zslice-self.block.zmin)/self.block.dz)
                if selfslice<0 or selfslice>nzd-1:
                    dataslice=None
                else:
                    dataslice=data[:,:,selfslice]
        if l_abs and dataslice is not None:dataslice=abs(dataslice)
        return slice,dataslice

    def genericpfem3d(self,data,title,titles=True,l_transpose=false,direction=None,slice=None,l_abs=False,
                      origins=None,deltas=None,display=1,scale=None,camera=None,l_box=1,
                      interactive=0,labels=['X','Y','Z'],palette=None,
                      opacitystart=None,opacityend=None,opacities=None,
                      adjust=0,labelscale=0.5,color='auto',ncolor=None,cmin=None,cmax=None,
                      niso=None,isomin=None,isomax=None,isos=None,opacity=0.5,
                      cscale=1.,l_csym=0,
                      procs=None,
                      xmmin=None,xmmax=None,
                      ymmin=None,ymmax=None,
                      zmmin=None,zmmax=None,
                      **kw):
        if direction is None and not l_opyndx:direction=2
        if self.l_2dxz:direction=1
        if 'view' in kw:
            view=kw['view']
        else:
            view=1
        if 'xscale' in kw:
            xscale=kw['xscale']
        else:
            xscale=1
        if 'yscale' in kw:
            yscale=kw['yscale']
        else:
            yscale=1
        if 'zscale' in kw:
            zscale=kw['zscale']
        else:
            zscale=1
        if 'gridscale' in kw:
            gridscale=kw['gridscale']
        else:
            gridscale=None
        if self.isactive:
            f=self.block.core.yf
            if self.l_1dz:
                nxd = nyd = 1
                nzd=shape(data)[0]
            elif self.l_2dxz:
                nyd = 1
                if len(shape(data))==3:data=data[:,0,:]
                nxd,nzd=shape(data)
            else:
                nxd,nyd,nzd=shape(data)
            if direction==0:
                if slice is None:
                    if self.l4symtry:
                        slice=0
                    else:
                        slice=self.nx/2
                xslice = w3d.xmmin+slice*w3d.dx
                selfslice = nint((xslice-self.block.xmin)/self.block.dx)
                if selfslice<0 or selfslice>nxd-1:
                    data=None
                    xmin=xmax=ymin=ymax=0.
                else:
                    if ymmin is None:
                        xmin=self.block.ymin
                    else:
                        xmin=ymmin
                    if ymmax is None:
                        xmax=self.block.ymax
                    else:
                        xmax=ymmax
                    if zmmin is None:
                        ymin=self.block.zmin+self.zgrid
                    else:
                        ymin=zmmin
                    if zmmax is None:
                        ymax=self.block.zmax+self.zgrid
                    else:
                        ymax=zmmax
                    if l_transpose:
                        xmin,ymin=ymin,xmin
                        xmax,ymax=ymax,xmax
                        data=transpose(data[selfslice,:,:])
                    else:
                        data=data[selfslice,:,:]
            if direction==1:
                if slice is None:
                    if self.l4symtry:
                        slice=0
                    else:
                        slice=self.ny/2
                yslice = w3d.ymmin+slice*w3d.dy
                selfslice = nint((yslice-self.block.ymin)/self.block.dy)
                if self.l_2dxz:slice=selfslice=0
                if selfslice<0 or selfslice>nyd-1:
                    data=None
                    xmin=xmax=ymin=ymax=0.
                else:
                    if xmmin is None:
                        xmin=self.block.xmin
                    else:
                        xmin=xmmin
                    if xmmax is None:
                        xmax=self.block.xmax
                    else:
                        xmax=xmmax
                    if zmmin is None:
                        ymin=self.block.zmin+self.zgrid
                    else:
                        ymin=zmmin
                    if zmmax is None:
                        ymax=self.block.zmax+self.zgrid
                    else:
                        ymax=zmmax
                    if l_transpose:
                        xmin,ymin=ymin,xmin
                        xmax,ymax=ymax,xmax
                    if self.l_2dxz:
                        if l_transpose and not self.l_1dz:
                            data=transpose(data)
                    else:
                        if l_transpose:
                            data=transpose(data[:,selfslice,:])
                        else:
                            data=data[:,selfslice,:]
            if direction==2:
                if slice is None:slice=self.nz/2
                zslice = w3d.zmmin+slice*w3d.dz
                selfslice = nint((zslice-self.block.zmin)/self.block.dz)
                if selfslice<0 or selfslice>nzd-1:
                    data=None
                    xmin=xmax=ymin=ymax=0.
                else:
                    if xmmin is None:
                        xmin=self.block.xmin
                    else:
                        xmin=xmmin
                    if xmmax is None:
                        xmax=self.block.xmax
                    else:
                        xmax=xmmax
                    if ymmin is None:
                        ymin=self.block.ymin
                    else:
                        ymin=ymmin
                    if ymmax is None:
                        ymax=self.block.ymax
                    else:
                        ymax=ymmax
                    if l_transpose:
                        xmin,ymin=ymin,xmin
                        xmax,ymax=ymax,xmax
                        data=transpose(data[:,:,selfslice])
                    else:
                        data=data[:,:,selfslice]
        if gridscale is not None and data is not None:
            data=data.copy()*gridscale
            kw['gridscale']=None
        if l_abs and data is not None:data=abs(data)
        if cmin is None:
            if data is None:
                mindata = 1.e36
            else:
                mindata = minnd(data)
            cmin=parallelmin(mindata)
        if cmax is None:
            if data is None:
                maxdata = -1.e36
            else:
                maxdata = maxnd(data)
            cmax=parallelmax(maxdata)
        if cscale is not None:
            cmin *= cscale
            cmax *= cscale
        if l_csym:
            cmax = max(abs(cmin),abs(cmax))
            cmin = -cmax
        if procs is None:procs=arange(npes)
        if direction in [0,1,2]:
            if not self.l_1dz:
                kw.setdefault('cmin',cmin)
                kw.setdefault('cmax',cmax)
            if me==0 and titles:
                if direction==0:
                    if l_transpose:
                        xtitle='Z';ytitle='Y'
                    else:
                        xtitle='Y';ytitle='Z'
                if direction==1:
                    if self.l_1dz:
                        xtitle='Z';ytitle=''
                    else:
                        if l_transpose:
                            xtitle='Z';ytitle='X'
                        else:
                            xtitle='X';ytitle='Z'
                if direction==2:
                    if l_transpose:
                        xtitle='Y';ytitle='X'
                    else:
                        xtitle='X';ytitle='Y'
                plsys(view)
                ptitles(title,xtitle,ytitle,'t = %gs'%(top.time))
            if isinstance(procs,types.IntType):procs=[procs]
            if me>0 and me in procs:
                mpisend(self.isactive, dest = 0, tag = 3)
                if self.isactive:
                    mpisend((xmin,xmax,ymin,ymax,data), dest = 0, tag = 3)
            else:
                lcolorbar=1
                if me in procs and data is not None:
                    if self.l_1dz:
                        plsys(view)
                        if kw.has_key('view'):kw.pop('view')
                        nz=shape(data)[0]
                        dz=(ymax-ymin)/nz
                        zmesh=ymin+arange(nz)*dz
                        pla(data,zmesh,**kw)
                    else:
                        kw.setdefault('xmin',xmin)
                        kw.setdefault('xmax',xmax)
                        kw.setdefault('ymin',ymin)
                        kw.setdefault('ymax',ymax)
                        ppgeneric(grid=data,**kw)
                        lcolorbar=0
                for i in range(1,npes):
                    isactive = mpirecv(source = i, tag = 3)
                    if isactive:
                        xminp,xmaxp,yminp,ymaxp,data=mpirecv(source = i, tag = 3)
                        if data is not None:
                            if self.l_1dz:
                                nz=shape(data)[0]
                                dz=(ymaxp-yminp)/nz
                                zmesh=yminp+arange(nz)*dz
                                pla(data,zmesh,**kw)
                            else:
                                kw['xmin']=xminp
                                kw['xmax']=xmaxp
                                kw['ymin']=yminp
                                kw['ymax']=ymaxp
                                ppgeneric(grid=data,lcolorbar=lcolorbar,**kw)
                                lcolorbar=0
        else:
            xmin=self.block.xmin
            xmax=self.block.xmax
            ymin=self.block.ymin
            ymax=self.block.ymax
            zmin=self.block.zmin
            zmax=self.block.zmax
            dx=self.block.dx
            dy=self.block.dy
            dz=self.block.dz
            if me>0 and me in procs:
                mpisend((xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz,data), dest = 0, tag = 3)
            else:
                # --- sets origins and deltas
                if origins is None:origins = [xmin*xscale,ymin*yscale,zmin*zscale]
                if deltas is None:deltas = [dx*xscale,dy*yscale,dz*zscale]
                # --- sets colormap
                if ncolor is None:ncolor = 2
                dec = (cmax-cmin)/ncolor
                colordata = cmin+arange(ncolor)*dec+dec/2
                if color=='auto':
                    if palette is not None:
                        try:
                            hp,rp,gp,bp = getpalhrgb(palette)
                            doautocolormap=False
                            color = transpose(array([rp,gp,bp])/255.)
                            ncolor = shape(color)[0]
                            colormap,opacities = Opyndx.DXColormap(data=colordata,ncolors=ncolor,
                                                            colors=color,min=cmin,max=cmax,
                                                            opacitystart=opacitystart,
                                                            opacityend=opacityend,
                                                            opacities=opacities)
                        except:
                            doautocolormap=True
                            print 'WARNING: Palette not found.'
                    else:
                        doautocolormap=True
                    if doautocolormap:
                        color = []
                        for i in range(ncolor):
                            r = array([1.,0.,0.])
                            g = array([0.,1.,0.])
                            b = array([0.,0.,1.])
                            if i<=ncolor/4:
                                color.append(b+g*(i*4./ncolor))
                            elif i<=ncolor/2:
                                color.append(b*(2.-(i*4./ncolor))+g)
                            elif i<=3*ncolor/4:
                                color.append(g+r*((i*4./ncolor)-2.))
                            else:
                                color.append(g*(4.-(i*4./ncolor))+r)
                        colormap,opacities  = Opyndx.DXColormap(data=colordata,
                                                       ncolors=ncolor,
                                                       colors=color,
                                                       opacitystart=opacitystart,
                                                       opacityend=opacityend,
                                                       opacities=opacities)
                Opyndx.DXReference(colormap)
                # --- sets isosurfaces values
                if isos is None:
                    if niso is None:niso=ncolor
                    if isomin is None:isomin = cmin+1.e-6*(cmax-cmin)
                    if isomax is None:isomax = cmax-1.e-6*(cmax-cmin)
                    dei = (isomax-isomin)/niso
                    isos = isomin+arange(niso)*dei+dei/2
                # --- renders isosurfaces
                if opacities is None:
                  opac = opacity
                else:
                  opac = opacities
                e3d,colorbar = Opyndx.viewisosurface1(data,isos,color=color,display=0,
                                origins=origins,
                                deltas=deltas,
                                opacity=opacity,
                                colormap=colormap)
#        dxob = [e3d,colorbar]
                dxob = [e3d]
                for i in procs:#range(1,npes):
                    if i != me:
                        xminp,xmaxp,dxp,yminp,ymaxp,dyp,zminp,zmaxp,dzp,data=mpirecv(source = i, tag = 3)
                        origins = [xminp*xscale,yminp*yscale,zminp*zscale]
                        deltas = [dxp*xscale,dyp*yscale,dzp*zscale]
                        Opyndx.DXReference(colormap)
                        e3d,colorbar = Opyndx.viewisosurface1(data,isos,color=color,display=0,
                                        origins=origins,
                                        deltas=deltas,
                                        opacity=opacity,
                                        colormap=colormap)
                        dxob.append(e3d)
                if l_box:
                    box = Opyndx.viewboundingbox(w3d.xmmin*xscale,
                                          w3d.xmmax*xscale,
                                          w3d.ymmin*yscale,
                                          w3d.ymmax*yscale,
                                          w3d.zmmin*zscale,
                                          w3d.zmmax*zscale,color='yellow')
                    dxob.append(box)
                if display:
                    dxob.append(colorbar)
                    dxob = Opyndx.DXCollect(dxob)
                    if camera is None:
                        camera = Opyndx.DXAutocamera(dxob,direction=[-1.,1.,-1.],width=100.,resolution=640,
                                              aspect=2.,up=[0,1,0],perspective=1,angle=60.,
                                              background='black')
                    Opyndx.DXImage(dxob,camera=camera,
                                 l_interactive=interactive,
                                 labels=labels,
                                 adjust=adjust,
                                 scale=scale,
                                 labelscale=labelscale)
        if display:
            return slice
        else:
            dxob = Opyndx.DXCollect(dxob)
            return dxob,colorbar#,slice

    ########################################################
    # Gather the requested array on processor 0.
    ########################################################
    def gatherarray(self,data,direction=None,slice=None,procs=None,guards=0) :
        """
        Gather the requested array on processor 0.
        Return only a slice of the array, if direction is in [0,1,2]

        Parameters :
        ------------
        data : ndarray
            An array of dimension ndim, containing the field on the local processor
            (without the guard cells) and which has to be combined with that of the
            other processors.

        direction : int
            A flag indicating the direction of the slice to be selected.
            If direction is not provided (slice=None), the function returns the full array.
            direction = 0 : zy slice
            direction = 1 : zx slice
            direction = 2 : xy slice

        slice : int
            Used only if direction is not None.
            The index of the slice to be gathered, in the global mesh.
            If slice is not provided (slice=None), this index is set
            to the middle of the global mesh.

        procs : ndarray or int
            An integer or array of integers that indicate the processors that
            have to send data.

        Returns :
        ---------
        A ndarray containing the values of the fields.
        If there is slicing (direction in [0,1,2]), the resulting array has dimension
        ndim-1. Otherwise, the array has dimension ndim.

        """

        if self.isactive: # If the local processor is active

            # Determine the shape of the data array
            ndim = data.ndim
            if ndim == 1 : # 1D
                nzd=shape(data)[0]
                nxd = 1
                nyd = 1
            elif ndim == 2 : # 2D Cartesian, or mode 0 of Circ
                nxd,nzd=shape(data)
                nyd = 1
            elif ( ndim == 3 and self.circ_m > 0 ) : # modes > 0 of Circ
                nxd,nzd=shape(data)[0:2]
            else : # 3D
                nxd,nyd,nzd=shape(data)

            # If direction is in [0,1,2] (i.e. there is slicing)
            # determine the index and data of the slice

            # Slicing has no meaning in the 1D case
            # Instead, return full array
            if self.l_1dz : direction = None
            # Slicing in y (xz slice) has no meaning in the 2D and Circ case
            # Instead, return full array
            if self.l_2dxz == 1 :
                if direction == 1 : direction = None

            # yz slice, at a given index and position in x
            if direction==0 :
                if slice is None :
                    if self.l4symtry:
                        slice=0
                    else:
                        slice=self.nx/2
                xslice = w3d.xmmin + slice*w3d.dx
                # selfslice : index of the slice on the local processor
                selfslice = nint((xslice-self.block.xmin)/self.block.dx)
                # Check whether the local processor is involved for this slice
                if selfslice < 0 or selfslice > nxd-1 :
                    data=None
                    xmin=xmax=ymin=ymax=0.
                else:
                    xmin=self.block.ymin
                    xmax=self.block.ymax
                    ymin=self.block.zmin
                    ymax=self.block.zmax
                    data=data[selfslice,...]

            # zx slice, at a given index and position in y
            elif direction==1 and not self.l_2dxz:
                if slice is None:
                    if self.l4symtry:
                        slice=0
                    else:
                        slice=self.ny/2
                yslice = w3d.ymmin+slice*w3d.dy
                # selfslice : index of the slice on the local processor
                selfslice = nint((yslice-self.block.ymin)/self.block.dy)
                if self.l_2dxz:slice=0
                # Check whether the local processor is involved for this slice
                if selfslice<0 or selfslice>nyd-1:
                    data=None
                    xmin=xmax=ymin=ymax=0.
                else:
                    xmin=self.block.xmin
                    xmax=self.block.xmax
                    ymin=self.block.zmin
                    ymax=self.block.zmax
                    data=data[:,selfslice,:]

            # xy slice, at a given index and position in z
            if direction==2 :
                if slice is None : slice=self.nz/2
                zslice = w3d.zmmin+slice*w3d.dz
                # selfslice : index of the slice on the local processor
                selfslice = nint((zslice-self.block.zmin)/self.block.dz)
                # Check whether the local processor is involved for this slice
                if selfslice<0 or selfslice>nzd-1:
                    data=None
                    xmin=xmax=ymin=ymax=0.
                else:
                    xmin=self.block.xmin
                    xmax=self.block.xmax
                    ymin=self.block.ymin
                    ymax=self.block.ymax
                    data=data[...,selfslice]

        if procs is None:procs=arange(npes)

        # Determine the number of guard cells
        if me==0:
            if guards:
                nxg=self.nxguard
                nyg=self.nyguard
                nzg=self.nzguard
            else:
                nxg=0
                nyg=0
                nzg=0

        # GATHER OPERATION

        # If there is slicing
        if direction in [0,1,2] and not (direction==1 and self.l_2dxz):
            if me==0:
                # Prepare the global array on the local processor
                if ndim == 2 : # 2D Cartesian, or mode 0 of Circ
                    if direction==0: datag = zeros([self.nz+1+nzg*2],'d')
                    if direction==2: datag = zeros([self.nx+1+nxg*2],'d')
                elif ( ndim == 3 and self.circ_m > 0 ) : # modes > 0 of Circ
                    if direction==0: datag = zeros([self.nz+1+nzg*2, self.circ_m],'complex')
                    if direction==2: datag = zeros([self.nx+1+nxg*2, self.circ_m],'complex')
                else: # 3D
                    if direction==0: datag = zeros([self.ny+1+nyg*2,self.nz+1+nzg*2],'d')
                    if direction==1: datag = zeros([self.nx+1+nxg*2,self.nz+1+nzg*2],'d')
                    if direction==2: datag = zeros([self.nx+1+nxg*2,self.ny+1+nyg*2],'d')
            else:
                datag=None
            if isinstance(procs,types.IntType):procs=[procs]
            # Find the processors that have to send data
            validdata = gatherlist( self.isactive and data is not None, bcast=1)
            validprocs = compress( validdata , arange(len(validdata)) )
            # Send data
            alldata = gatherlist( [xmin,ymin,data], dest=0, procs=validprocs )
            barrier() # This ensures that processor 0 will not bet overflowed with messages
            if me==0:
                # Fill the global array datag with the local arrays contained in alldata
                for i in range(len(alldata)):
                    xminp = alldata[i][0]
                    yminp = alldata[i][1]
                    data  = alldata[i][-1]
                    if data is not None:
                        if self.l_2dxz : # 2D and Circ
                            if direction==0:
                                ny = shape(data)[0]
                                iymin = nint((yminp-self.zmmin)/self.dz)+nzg
                                datag[iymin:iymin+ny] = data[...]
                            if direction==2:
                                nx = shape(data)[0]
                                ixmin = nint((xminp-self.xmmin)/self.dx)+nxg
                                datag[ixmin:ixmin+nx] = data[...]
                        else: # 3D
                            nx,ny = shape(data)
                            if direction==0:
                                ixmin = nint((xminp-self.ymmin)/self.dy)+nyg
                                iymin = nint((yminp-self.zmmin)/self.dz)+nzg
                            if direction==1:
                                ixmin = nint((xminp-self.xmmin)/self.dx)+nxg
                                iymin = nint((yminp-self.zmmin)/self.dz)+nzg
                            if direction==2:
                                ixmin = nint((xminp-self.xmmin)/self.dx)+nxg
                                iymin = nint((yminp-self.ymmin)/self.dy)+nyg
                            datag[ixmin:ixmin+nx,iymin:iymin+ny] = data[...]

        # If there is no slicing
        else:
            if me==0:
                # Prepare the global array on the local processor
                if ndim==1 : # 1D
                    datag = zeros([self.nz+1+nzg*2],'d')
                elif ndim==2 : # 2D and mode 0 of Circ
                    datag = zeros([self.nx+1+nxg*2,self.nz+1+nzg*2],'d')
                elif (ndim==3 and self.circ_m > 0) : # modes > 0 Circ
                    datag = zeros([self.nx+1+nxg*2,self.nz+1+nzg*2,self.circ_m],'complex')
                else: # 3D
                    datag = zeros([self.nx+1+nxg*2,self.ny+1+nyg*2,self.nz+1+nzg*2],'d')
            else:
                datag = None
            xmin=self.block.xmin
            xmax=self.block.xmax
            ymin=self.block.ymin
            ymax=self.block.ymax
            zmin=self.block.zmin
            zmax=self.block.zmax
            dx=self.block.dx
            dy=self.block.dy
            dz=self.block.dz
            # Determine the list of procs that have to send data
            if procs is None : procs=range(npes)
            if type(procs) is int : procs=[procs]
            # Send the data
            if me>0 and me in procs:
                mpisend((xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz,data), dest = 0, tag = 3)
            else:
                for i in procs :
                    # Receive the data on proc 0
                    if i != me:
                        xminp,xmaxp,dxp,yminp,ymaxp,dyp,zminp,zmaxp,dzp,data=mpirecv(source = i, tag = 3)
                    else:
                        xminp = xmin
                        yminp = ymin
                        zminp = zmin
                    # Fill the global array datag with the local arrays sent
                    if data is not None:
                        if self.l_1dz:
                            nz = shape(data)[0]
                            izmin = nint((zminp-self.zmmin)/self.dz)+nzg
                            datag[izmin:izmin+nz] = data[...]
                        elif self.l_2dxz: # 2D and Circ
                            nx,nz = shape(data)[0:2]
                            ixmin = nint((xminp-self.xmmin)/self.dx)+nxg
                            izmin = nint((zminp-self.zmmin)/self.dz)+nzg
                            datag[ixmin:ixmin+nx,izmin:izmin+nz] = data[...]
                        else:
                            nx,ny,nz = shape(data)
                            ixmin = nint((xminp-self.xmmin)/self.dx)+nxg
                            iymin = nint((yminp-self.ymmin)/self.dy)+nyg
                            izmin = nint((zminp-self.zmmin)/self.dz)+nzg
                            datag[ixmin:ixmin+nx,iymin:iymin+ny,izmin:izmin+nz] = data[...]
            barrier() # This ensures that processor 0 will not bet overflowed with messages

        # Return the global array
        return(datag)

    def getarray(self,g,guards=0,overlap=0):
        if guards:
            return g
        else:
            f=self.fields
            ox=oy=oz=0
            if not overlap:
                if self.block.xrbnd==em3d.otherproc:ox=1
                if self.block.yrbnd==em3d.otherproc:oy=1
                if self.block.zrbnd==em3d.otherproc:oz=1
            if self.l_1dz:
                return g[0,0,f.nzguard:-f.nzguard-oz]
            elif self.l_2dxz:
                return g[f.nxguard:-f.nxguard-ox,0,f.nzguard:-f.nzguard-oz]
            else:
                return g[f.nxguard:-f.nxguard-ox,f.nyguard:-f.nyguard-oy,f.nzguard:-f.nzguard-oz]

    def getarray_circ(self,g,guards=0,overlap=0):
        if guards:
            return g
        else:
            f=self.fields
            ox=oy=oz=0
            if not overlap:
                if self.block.xrbnd==em3d.otherproc:ox=1
                if self.block.yrbnd==em3d.otherproc:oy=1
                if self.block.zrbnd==em3d.otherproc:oz=1
            return g[f.nxguard:-f.nxguard-ox,f.nzguard:-f.nzguard-oz,:]

    def step(self,n=1,freq_print=10,lallspecl=0,l_alawarpx=False):
        for i in range(n):
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
            if l_alawarpx:
                self.onestep_alawarpx(l_first,l_last)
            else:
                self.onestep(l_first,l_last)

    def onestep(self,l_first,l_last):
        # --- call beforestep functions
        callbeforestepfuncs.callfuncsinlist()

        top.zgrid+=top.vbeamfrm*top.dt
        top.zbeam=top.zgrid
        # --- gather fields from grid to particles
#        w3d.pgroupfsapi = top.pgroup
#        for js in range(top.pgroup.ns):
#          self.fetcheb(js)

        # --- push
        if l_first:
            for js in range(top.pgroup.ns):
                self.push_velocity_second_half(js)
                self.record_old_positions(js)
                self.push_positions(js)
        else:
            for js in range(top.pgroup.ns):
                self.push_velocity_full(js)
                self.record_old_positions(js)
                self.push_positions(js)

        inject3d(1, top.pgroup)

        # --- call user-defined injection routines
        userinjection.callfuncsinlist()

        particleboundaries3d(top.pgroup,-1,False)

        # --- call beforeloadrho functions
        beforeloadrho.callfuncsinlist()

#        self.loadsource()
        self.loadrho()
        self.loadj()

#        self.solve2ndhalf()
        self.dosolve()

        w3d.pgroupfsapi = top.pgroup
        for js in range(top.pgroup.ns):
          self.fetcheb(js)

        if l_last:
            for js in range(top.pgroup.ns):
                self.push_velocity_first_half(js)

        # --- update time, time counter
        top.time+=top.dt
        if top.it%top.nhist==0:
#           zmmnt()
           minidiag(top.it,top.time,top.lspecial)
        top.it+=1

        # --- call afterstep functions
        callafterstepfuncs.callfuncsinlist()

    def onestep_alawarpx(self,l_first,l_last):
        if self.ntsub<1.:
            self.novercycle = nint(1./self.ntsub)
            self.icycle = (top.it-1)%self.novercycle
        else:
            self.novercycle = 1
            self.icycle = 0

        # --- call beforestep functions
        callbeforestepfuncs.callfuncsinlist()

        # --- push by half step if first step (where all are synchronized)
        if l_first:
            if not self.solveroff:
                self.allocatedataarrays()
                self.push_e(l_half=True)
                self.exchange_e()
            for js in range(top.pgroup.ns):
                self.push_positions(js,dtmult=0.5)

        if not self.solveroff:
            for i in range(int(self.ntsub)-1):
                self.push_b_full()
                if self.l_pushf:self.exchange_f()
                self.exchange_b()
                self.push_e_full(i)
                self.exchange_e()
            self.push_b_part_1()

            if self.pml_method==2:
                scale_em3d_bnd_fields(self.block,top.dt,self.l_pushf,self.l_pushg)

            if self.l_pushf:self.exchange_f()
            self.exchange_b()
            self.setebp()

            if not self.l_nodalgrid and top.efetch[0] != 4:
                if self.l_fieldcenterK:
                    self.fieldcenterK()
                else:
                    self.yee2node3d()
            if self.l_nodalgrid and top.efetch[0] == 4:
                self.fields.l_nodecentered=True
                self.node2yee3d()
            if self.l_smooth_particle_fields and any(self.npass_smooth>0):
                self.smoothfields()
            if self.l_correct_num_Cherenkov:self.smoothfields_poly()

        top.zgrid+=top.vbeamfrm*top.dt
        top.zbeam=top.zgrid

        w3d.pgroupfsapi = top.pgroup
        for js in range(top.pgroup.ns):
            self.fetcheb(js)
            self.push_velocity_full(js)
            self.record_old_positions(js)
            self.push_positions(js)

        inject3d(1, top.pgroup)

        # --- call user-defined injection routines
        userinjection.callfuncsinlist()

        particleboundaries3d(top.pgroup,-1,False)

        self.solve2ndhalf()

        # --- call beforeloadrho functions
        if isinstalledbeforeloadrho(self.solve2ndhalf):
            uninstallbeforeloadrho(self.solve2ndhalf)
        beforeloadrho.callfuncsinlist()

#        self.loadsource()
        self.loadrho()
        self.loadj()

        self.getconductorobject()

        if l_last:
            for js in range(top.pgroup.ns):
                self.record_old_positions(js)
                self.push_positions(js,dtmult=-0.5)
            if not self.solveroff:self.push_e(l_half=True)
        else:
            if not self.solveroff:
                self.push_e()
                self.exchange_e()

        # --- update time, time counter
        top.time+=top.dt
        if top.it%top.nhist==0:
#           zmmnt()
           minidiag(top.it,top.time,top.lspecial)
        top.it+=1

        # --- call afterstep functions
        callafterstepfuncs.callfuncsinlist()

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
        pg.ex[il:iu]=top.ex0
        pg.ey[il:iu]=top.ey0
        pg.ez[il:iu]=top.ez0
        pg.bx[il:iu]=top.bx0
        pg.by[il:iu]=top.by0
        pg.bz[il:iu]=top.bz0
        self.fetche()
        self.fetchb()

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

    def push_velocity_full(self,js,pg=None):
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

        if self.l_verbose:print me,'exit push_velocity_second_half'

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

        if self.l_verbose:print me,'exit set_gamma'

    def push_positions(self,js,pg=None,dtmult=1.):
        if self.l_verbose:print me,'enter push_positions'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        dt = top.dt*dtmult
        xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                       pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                       pg.gaminv[il:iu],dt)

        if self.l_verbose:print me,'exit push_positions'

    def record_old_positions(self,js,pg=None,dtmult=1.):
        if self.l_verbose:print me,'enter record_old_positions'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        setuppgroup(pg)
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        dt = top.dt*dtmult
        if top.xoldpid>0:pg.pid[il:iu,top.xoldpid-1] = pg.xp[il:iu].copy()
        if top.yoldpid>0:pg.pid[il:iu,top.yoldpid-1] = pg.yp[il:iu].copy()
        if top.zoldpid>0:pg.pid[il:iu,top.zoldpid-1] = pg.zp[il:iu].copy()
        if top.uxoldpid>0:pg.pid[il:iu,top.uxoldpid-1] = pg.uxp[il:iu].copy()
        if top.uyoldpid>0:pg.pid[il:iu,top.uyoldpid-1] = pg.uyp[il:iu].copy()
        if top.uzoldpid>0:pg.pid[il:iu,top.uzoldpid-1] = pg.uzp[il:iu].copy()

        if self.l_verbose:print me,'exit record_old_positions'

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

    #############################################################
    #                   Plotting methods                        #
    #############################################################

    def pfex( self, l_children=1, guards=0, direction=None, m=None,
              iz=w3d.nz/2, output=False, show=True, **kw) :
        """
        Plot the Ex field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis
            and making an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfvec_circ and genericpfem3d for further options

        """

        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            er0 = self.gatherarray( self.getarray( self.fields.Exp,guards,overlap=True) )
            etheta0 = self.gatherarray( self.getarray( self.fields.Eyp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                er = self.gatherarray( \
                            self.getarray_circ( self.fields.Exp_circ,guards,overlap=True) )
                etheta = self.gatherarray( \
                            self.getarray_circ( self.fields.Eyp_circ,guards,overlap=True ) )
            else :
                er = None
                etheta = None
            # Sum and plot them
            f = self.pfvec_circ( er0, etheta0, er, etheta, name='Ex^^',
                            component='x', m=m, direction=direction,
                            output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.Exp,guards,
                               overlap=True)
            if show :
                self.genericpfem3d( f,'E_x', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])

    def pfey( self, l_children=1, guards=0, direction=None, m=None,
                iz=w3d.nz/2, output=False, show=True, **kw) :
        """
        Plot the Ey field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots the sum
            of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfvec_circ and genericpfem3d for further options

        """

        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            er0 = self.gatherarray( self.getarray( self.fields.Exp,guards,overlap=True) )
            etheta0 = self.gatherarray( self.getarray( self.fields.Eyp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                er = self.gatherarray( self.getarray_circ( \
                                            self.fields.Exp_circ,guards,overlap=True ) )
                etheta = self.gatherarray( self.getarray_circ( \
                                            self.fields.Eyp_circ,guards,overlap=True ) )
            else :
                er = None
                etheta = None
            # Sum and plot them
            f = self.pfvec_circ( er0, etheta0, er, etheta, name='Ey^^',
                    component='y', m=m, direction=direction,
                    output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.Eyp,guards,overlap=True)
            if show :
                self.genericpfem3d(f, 'E_y', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])


    def pfez( self, l_children=1, guards=0, direction=None, m=None,
                iz=w3d.nz/2, output=False, show=True, **kw) :
        """
        Plot the Ez field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots the sum
            of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """

        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            ez0 = self.gatherarray( self.getarray( self.fields.Ezp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                ez = self.gatherarray( self.getarray_circ( \
                                    self.fields.Ezp_circ,guards,overlap=True ) )
            else :
                ez = None
            # Sum and plot them
            f = self.pfscalar_circ( ez0, ez, name='Ez^^', m=m,
                    direction=direction, output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.Ezp,guards,overlap=True)
            if show :
                self.genericpfem3d( f, 'E_z', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])

    def pfer( self, l_children=1, guards=0, direction=None, m=None,
                output=False, show=True, **kw) :
        """
        Plot the Er field (only when using azimuthal
        decomposition (l_2drz = 1))

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            er0 = self.gatherarray( self.getarray( self.fields.Exp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                er = self.gatherarray( self.getarray_circ( \
                            self.fields.Exp_circ,guards,overlap=True ) )
            else :
                er = None
            # Sum and plot them
            f = self.pfscalar_circ( er0, er, name='Er^^', m=m,
                direction=direction, output=output, show=show, **kw )

            if output and me==0 :
                return(f[::-1])

    def pfet( self, l_children=1, guards=0, direction=None, m=None,
                            output=False, show=True, **kw) :
        """
        Plot the E_theta field (only when using azimuthal decomposition
         (l_2drz = 1))

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis
            and making an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            et0 = self.gatherarray( self.getarray( self.fields.Eyp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                et = self.gatherarray( self.getarray_circ( \
                                        self.fields.Eyp_circ,guards,overlap=True ) )
            else :
                et = None
            # Sum and plot them
            f = self.pfscalar_circ( et0, et, name='Etheta^^', m=m,
                    direction=direction, output=output, show=show, **kw )

            if output and me==0 :
                return(f[::-1])


    def pfbx( self, l_children=1, guards=0, direction=None, m=None,
                            iz=w3d.nz/2, output=False, show=True, **kw) :
        """
        Plot the Bx field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis
            and making an angle theta_plot with the x axis

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfvec_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            br0 = self.gatherarray( self.getarray( self.fields.Bxp,guards,overlap=True) )
            btheta0 = self.gatherarray( self.getarray( self.fields.Byp,guards,overlap=True) )
            if self.circ_m > 0 :
                # Modes > 0
                br = self.gatherarray( self.getarray_circ( \
                                        self.fields.Bxp_circ,guards,overlap=True ) )
                btheta = self.gatherarray( self.getarray_circ( \
                                        self.fields.Byp_circ,guards,overlap=True ) )
            else :
                br = None
                btheta = None
            # Sum and plot them
            f = self.pfvec_circ( br0, btheta0, br, btheta, name='Bx^^',
                            component='x', m=m, direction=direction,
                            output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.Bxp,guards,overlap=True)
            if show :
                self.genericpfem3d( f, 'B_x', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])


    def pfby( self, l_children=1, guards=0, direction=None, m=None,
                            iz=w3d.nz/2, output=False, show=True, **kw) :
        """
        Plot the By field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis
            and making an angle theta_plot with the x axis

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfvec_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gathbr the Fourier coefficients
            # Mode 0
            br0 = self.gatherarray( self.getarray( self.fields.Bxp,guards,overlap=True) )
            btheta0 = self.gatherarray( self.getarray( self.fields.Byp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                br = self.gatherarray( self.getarray_circ( \
                                            self.fields.Bxp_circ,guards,overlap=True ) )
                btheta = self.gatherarray( self.getarray_circ( \
                                            self.fields.Byp_circ,guards,overlap=True ) )
            else :
                br = None
                btheta = None
            # Sum and plot them
            f = self.pfvec_circ( br0, btheta0, br, btheta, name='By^^',
                             component='y', m=m, direction=direction,
                             output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.Byp,guards,overlap=True)
            if show :
                self.genericpfem3d(f, 'B_y', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])


    def pfbz( self, l_children=1, guards=0, direction=None, m=None,
          iz=w3d.nz/2, output=False, show=True, **kw) :
        """
        Plot the Bz field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots the sum
            of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            bz0 = self.gatherarray( self.getarray( self.fields.Bzp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                bz = self.gatherarray( self.getarray_circ( \
                                        self.fields.Bzp_circ,guards,overlap=True ) )
            else :
                bz = None
            # Sum and plot them
            f = self.pfscalar_circ( bz0, bz, name='Bz^^', m=m,
                direction=direction, output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.Bzp,guards,overlap=True)
            if show :
                self.genericpfem3d( f, 'B_z', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])

    def pfbr( self, l_children=1, guards=0, direction=None, m=None,
                            output=False, show=True, **kw) :
        """
        Plot the Br field (only when using azimuthal
        decomposition (l_2drz = 1))

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis

        m : int
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            br0 = self.gatherarray( self.getarray( self.fields.Bxp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                br = self.gatherarray( self.getarray_circ( \
                                        self.fields.Bxp_circ,guards,overlap=True ) )
            else :
                br = None
            # Sum and plot them
            f = self.pfscalar_circ( br0, br, name='Br^^', m=m,
                direction=direction, output=output, show=show, **kw )

            if output and me==0 :
                return(f[::-1])

    def pfbt( self, l_children=1, guards=0, direction=None, m=None,
                            output=False, show=True, **kw) :
        """
        Plot the B_theta field (only when using azimuthal
        decomposition (l_2drz = 1))

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis

        m : int
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            btheta0 = self.gatherarray( self.getarray( self.fields.Byp,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                btheta = self.gatherarray( self.getarray_circ( \
                                            self.fields.Byp_circ,guards,overlap=True ) )
            else :
                btheta = None
            # Sum and plot them
            f = self.pfscalar_circ( btheta0, btheta, name='Bt^^', m=m,
                direction=direction, output=output, show=show, **kw )

            if output and me==0 :
                return(f[::-1])

    def pfjx( self, l_children=1, guards=0, direction=None, m=None,
                            output=False, show=True, **kw) :
        """
        Plot the Jx field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis
            and making an angle theta_plot with the x axis

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfvec_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            jr0 = self.gatherarray( self.getarray( self.fields.Jx,guards,overlap=True) )
            jtheta0 = self.gatherarray(self.getarray(self.fields.Jy,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                jr = self.gatherarray( self.getarray_circ( \
                                        self.fields.Jx_circ,guards,overlap=True ) )
                jtheta = self.gatherarray( self.getarray_circ( \
                                        self.fields.Jy_circ,guards,overlap=True ) )
            else :
                jr = None
                jtheta = None
            # Sum and plot them
            f = self.pfvec_circ( jr0, jtheta0, jr, jtheta, name='Jx^^',
                            component='x', m=m, direction=direction,
                            output=output, show=show, **kw )

        else: # Cartesian fields
            f =  self.getarray(self.fields.Jx,
                               guards,overlap=True)
            if show :
                self.genericpfem3d( f, 'J_x', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])

    def pfjy( self, l_children=1, guards=0, direction=None, m=None,
                            output=False, show=True, **kw) :
        """
        Plot the Jy field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and
            making an angle theta_plot with the x axis

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfvec_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            jr0 = self.gatherarray( self.getarray( self.fields.Jx,guards,overlap=True) )
            jtheta0 = self.gatherarray(self.getarray(self.fields.Jy,guards,overlap=True))
            if self.circ_m > 0 :
                jr = self.gatherarray( self.getarray_circ( \
                                        self.fields.Jx_circ,guards,overlap=True ) )
                jtheta = self.gatherarray( self.getarray_circ( \
                                        self.fields.Jy_circ,guards,overlap=True ) )
            else :
                jr = None
                jtheta = None
            # Sum and plot them
            f = self.pfvec_circ( jr0, jtheta0, jr, jtheta, name='Jy^^', \
                             component='y', m=m, direction=direction,
                             output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.Jy,guards,
                              overlap=True)
            if show :
                self.genericpfem3d( f, 'J_y', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])

    def pfjz( self, l_children=1, guards=0, direction=None, m=None,
                            output=False, show=True, **kw) :
        """
        Plot the Jz field.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            jz0 = self.gatherarray( self.getarray( self.fields.Jz,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                jz = self.gatherarray(self.getarray_circ( \
                                    self.fields.Jz_circ,guards,overlap=True ) )
            else :
                jz = None
            # Sum and plot them
            f = self.pfscalar_circ( jz0, jz, name='Jz^^', m=m,
                direction=direction, output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.Jz,guards,
                              overlap=True)
            if show :
                self.genericpfem3d( f, 'J_z', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])


    def pfjr( self, l_children=1, guards=0, direction=None, m=None,
                            output=False, show=True, **kw) :
        """
        Plot the Jr field (only when using azimuthal
        decomposition (l_2drz = 1))

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis

        m : int
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            jr0 = self.gatherarray( self.getarray( self.fields.Jx,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                jr = self.gatherarray( self.getarray_circ( \
                                    self.fields.Jx_circ,guards,overlap=True ) )
            else :
                jr = None
            # Sum and plot them
            f = self.pfscalar_circ( jr0, jr, name='Jr^^', m=m,
                direction=direction, output=output, show=show, **kw )

            if output and me==0 :
                return(f[::-1])

    def pfjt( self, l_children=1, guards=0, direction=None, m=None,
                            output=False, show=True, **kw) :
        """
        Plot the J_theta field (only when using azimuthal
        decomposition (l_2drz = 1))

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis

        m : int
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            jt0 = self.gatherarray( self.getarray( self.fields.Jy,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                jt = self.gatherarray( self.getarray_circ( \
                                    self.fields.Jy_circ,guards,overlap=True ) )
            else :
                jt = None
            # Sum and plot them
            self.pfscalar_circ( jt0, jt, name='Jt^^', m=m,
                direction=direction, output=output, show=show, **kw )

            if output and me==0 :
                return(f[::-1])


    def pfexg(self,l_children=1,guards=0,direction=None,**kw):
        self.genericpfem3d(self.getarray(self.fields.Ex,guards,overlap=True),'Eg_x',
        direction=direction,**kw)

    def pfeyg(self,l_children=1,guards=0,direction=None,**kw):
        self.genericpfem3d(self.getarray(self.fields.Ey,guards,overlap=True),'Eg_y',
        direction=direction,**kw)

    def pfezg(self,l_children=1,guards=0,direction=None,**kw):
        self.genericpfem3d(self.getarray(self.fields.Ez,guards,overlap=True),'Eg_z',
        direction=direction,**kw)

    def pfbxg(self,l_children=1,guards=0,direction=None,**kw):
        self.genericpfem3d(self.getarray(self.fields.Bx,guards,overlap=True),'Bg_x',
        direction=direction,**kw)

    def pfbyg(self,l_children=1,guards=0,direction=None,**kw):
        self.genericpfem3d(self.getarray(self.fields.By,guards,overlap=True),'Bg_y',
        direction=direction,**kw)

    def pfbzg(self,l_children=1,guards=0,direction=None,**kw):
        self.genericpfem3d(self.getarray(self.fields.Bz,guards,overlap=True),'Bg_z',
        direction=direction,**kw)

    def pff( self, l_children=1, guards=0, direction=None, m=None,
                           output=False, show=True, **kw) :
        """
        Plot the F pseudo field (dF/dt = div E - rho).

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis
            and making an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots
            the sum of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """

        if self.l_2drz == 1 : # Azimuthal decomposition
            # Gather the Fourier coefficients
            # Mode 0
            f0 = self.gatherarray( self.getarray( self.fields.F,guards,overlap=True) )
            # Modes > 0
            if self.circ_m > 0 :
                f = self.gatherarray(self.getarray_circ(self.fields.F_circ,guards,overlap=True))
            else :
                f = None
            # Sum and plot them
            f = self.pfscalar_circ( f0, f, name='F^^', m=m,
                direction=direction, output=output, show=show, **kw )

        else: # Cartesian fields
            f = self.getarray(self.fields.F,guards,overlap=True)
            self.genericpfem3d( f, 'F', direction=direction,**kw)

        if output:
            if self.l_2drz == 0 : # Cartesian fields
                f=self.gatherarray(f, direction=direction)
            if me==0:return(f[::-1])

    def pfdive(self,l_children=1,guards=0,direction=None,**kw):
        self.genericpfem3d(self.getdive(guards,overlap=True),'div(E)',
        direction=direction,**kw)

    def pfe(self,l_children=1,guards=0,direction=None,**kw):
        e = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap=True)
        self.genericpfem3d(sqrt(e),'E',
        direction=direction,**kw)

    def pfb(self,l_children=1,guards=0,direction=None,**kw):
        b = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap=True)
        self.genericpfem3d(sqrt(b),'B',
        direction=direction,**kw)

    def pfw(self,l_children=1,guards=0,direction=None,**kw):
        e = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap=True)
        b = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap=True)
        self.genericpfem3d(sqrt(e+b),'W',
        direction=direction,**kw)


    def sum_modes_circ(self, ff0, ff, modes, direction=1, theta_plot=0,
                       iz=w3d.nz/2 ) :
        """
        Calculates the real field corresponding to the Fourier coefficients
        of the specified modes, over a 2D plane.

        Note that the calculation of the sum of the contributions of
        the mode uses the same formulas as getf2dxz_n (for the mode 0)
        and getf2drz_circ_n (for the modes > 0), with the exception that
        it does not do the sum of the r and theta component of a vector
        field to obtain the x or y component (use fvec_circ for this)

        Parameters
        ----------
        ff0 : ndarray
            A 2d array of real values with shape [ w3d.nx+1, w3d.nz+1 ],
            which contains the azimuthal Fourier coefficients of the field
            for the mode m=0.

        ff : ndarray
            A 3d array of complex values with shape [ w3d.nx+1, w3d.nz+1 ,
            self.fields.circ_m ], which contains the azimuthal Fourier
            coefficients of the field for the modes m>0.

        modes : list
            A list of integers indicating the modes whose contributions are
            to be summed. The other modes are ignored.

        direction : int
            A flag that defines in which plame the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis

        theta_plot : float
            Used only if direction = -1
            The angle (in radians) between the plane of the plot and the
            zx plane.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        Returns
        -------
            A tuple of two 2D array with real values
            - the first array contains the values of fields in the plane
            defined by direction
            - the second array contains the corresponding values of theta
              (needed e.g. when doing the vector sum of two components)
            If direction is in [-1,0,1], the first axis of the returned
            arrays corresponds to z.
            If direction is 2, the first and second axis of the returned
            array correspond to x and y respectively.
        """

        # If the plane of the plot includes the z axis
        if direction in [-1,0,1] :

            # Prepare the theta array
            if direction == 0 : theta_plot = pi/2  # zy plane
            if direction == 1 : theta_plot = 0.    # zx plane
            theta = zeros([ 2*w3d.nx+1 , w3d.nz ])
            # Upper half of the plane
            theta[ w3d.nx: , : ] = theta_plot
            # Lower half of the plane
            theta[ :w3d.nx , : ] = theta_plot + pi

            # Calculate the fields in the given plane, from its Fourier
            # coefficients
            f = zeros([ 2*w3d.nx+1 , w3d.nz ])
            # Carry out the sum
            for m_mode in modes :
                if m_mode == 0 : # Mode 0
                    # Upper half of the plane
                    f[ w3d.nx: , : ] += ff0[ 0:w3d.nx+1 , :w3d.nz ]
                    # Lower half of the plane
                    f[ w3d.nx-1::-1 , : ] += ff0[ 1:w3d.nx+1 , :w3d.nz ]
                else : # Modes > 0
                    # Upper half of the plane
                    f[ w3d.nx: , : ] += \
                    cos( m_mode*theta_plot ) * ff[ 0:w3d.nx+1, :w3d.nz ,
                        m_mode-1].real + sin( m_mode*theta_plot ) * \
                         ff[ 0:w3d.nx+1, :w3d.nz , m_mode-1].imag
                    # Lower half of the plane
                    f[ w3d.nx-1::-1 , : ] += \
                    cos( m_mode*(theta_plot+pi) ) * ff[ 1:w3d.nx+1, :w3d.nz ,
                        m_mode-1].real + sin( m_mode*(theta_plot+pi) ) * \
                        ff[ 1:w3d.nx+1, :w3d.nz , m_mode-1].imag

            # Transpose the array, so that the first axis corresponds to z
            f = f.T
            theta = theta.T

        # If the plot is in the xy plane
        if direction == 2 :

            # Get the coordinates on which to calculate the field,
            # in the xy plane
            xx , yy = getmesh2d( -w3d.xmmax,w3d.dx,w3d.nx*2,-w3d.xmmax,
                                 w3d.dx,w3d.nx*2 )
            r = sqrt(xx*xx+yy*yy).flatten()
            theta=arctan2(yy,xx).flatten()

            # Calculate the fields in the given plane,
            # from its Fourier coefficients
            f = zeros( (2*w3d.nx+1)**2 )
            # Carry out the sum
            for m_mode in modes :
                if m_mode == 0 : # Mode 0
                    # Interpolate the Fourier coefficients
                    # (Arrays are flattened so as to adapt to
                    # the interface of getgrid1d)
                    f_real = zeros((2*w3d.nx+1)**2)
                    getgrid1d((2*w3d.nx+1)**2,r,f_real,w3d.nx,ff0[:,iz],
                              0.,w3d.xmmax)
                    # Sum the mode
                    f += f_real
                else : # Modes > 0
                    # Interpolate the Fourier coefficients
                    # (Arrays are flattened so as to adapt to the
                    # interface of getgrid1d)
                    f_real = zeros((2*w3d.nx+1)**2)
                    f_imag = zeros((2*w3d.nx+1)**2)
                    getgrid1d((2*w3d.nx+1)**2,r,f_real,w3d.nx,
                              ff[:,iz,m_mode-1].real,0.,w3d.xmmax)
                    getgrid1d((2*w3d.nx+1)**2,r,f_imag,w3d.nx,
                              ff[:,iz,m_mode-1].imag,0.,w3d.xmmax)
                    # Sum the mode
                    f += f_real*cos(m_mode*theta) + f_imag*sin(m_mode*theta)

            # Unflatten the arrays
            f=f.reshape([2*w3d.nx+1,2*w3d.nx+1])
            theta = theta.reshape([2*w3d.nx+1,2*w3d.nx+1])

        # Finally return the arrays
        return(f, theta)

    def pfscalar_circ(self, ff0, ff, name='', m=None, direction=1,
                      theta_plot=0, iz=w3d.nz/2, titles=1,
                      output=False, show=True, l_transpose=None, **kw):
        """
        Plot a scalar field from its azimuthal Fourier coefficients
        ff0 (mode 0) and ff (modes > 0).

        Parameters
        ----------
        ff0 : ndarray
            A 2d array of real values with shape [ w3d.nx+1, w3d.nz+1 ],
            which contains the azimuthal Fourier coefficients of the
            field for the mode m=0.

        ff : ndarray
            A 3d array of complex values with shape [ w3d.nx+1, w3d.nz+1 ,
            self.fields.circ_m ], which contains the azimuthal Fourier
            coefficients of the field for the modes m>0.

        name : str
            The name of the field considered, which is then printed on the plot.

        m : int
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots the sum
            of the contribution of all modes.

        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis

        theta_plot : float
            Used only if direction = -1
            The angle (in radians) between the plane of the plot
            and the zx plane.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        titles : int
            A flag that indicates whether the plots should be
            automatically labeled

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        """

        if me == 0 :

            # Sets default direction to 1
            if ( direction in [-1,0,1,2] ) == False : direction = 1

            # Prepare sum over the modes
            if m != None : # Single mode
                assert m<=self.fields.circ_m, "m > circ_m: %g > %g" \
                    %(m,self.fields.circ_m)
                modes = [m]
            else : # All modes
                modes = range( self.fields.circ_m+1 )

            # Sum all the selected modes, in the plane given by direction
            # (Dump the returned theta array, since it is not used here)
            f, _ = self.sum_modes_circ( ff0, ff, direction=direction,
                            modes=modes, theta_plot=theta_plot, iz=iz )

            # Set the boundaries (here X and Y represent the axes of the plot)
            # - If the plane of the plot includes the z axis
            if direction in [-1,0,1] :
                Xmin = w3d.zmmin + self.zgrid
                Xmax = w3d.zmmax + self.zgrid
                Ymin = -w3d.xmmax
                Ymax = w3d.xmmax
            # - If the plot is in the xy plane
            if direction == 2 :
                Xmin = -w3d.xmmax
                Xmax = w3d.xmmax
                Ymin = -w3d.xmmax
                Ymax = w3d.xmmax

            # Plot function : colormap of the matrix f
            if show :
                ppgeneric( f, xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax, **kw)

                # Labeling of the plot
                if titles == 1 :
                    if m!=None :
                        name = name+' (m=%g)'%m
                    if direction == -1 :
                        ptitles(name,'Z','Coordinate along theta=%f' \
                                %theta_plot)
                    elif direction == 0 :
                        ptitles(name,'Z','Y')
                    elif direction == 1 :
                        ptitles(name,'Z','X')
                    elif direction == 2 :
                        ptitles(name,'X','Y')

            if output and me==0 :
                return( f.T )

    def pfvec_circ(self, ffr0, fftheta0, ffr, fftheta, name='',
                   component='x', m=None, direction=1, theta_plot=0,
                   iz=w3d.nz/2, titles=1, output=False, show=True, l_transpose=None, **kw):
        """
        Plot the x or y component of a field from the azimuthal Fourier
        coefficients of its radial and azimuthal components ffr0, fft0
        (mode 0) and ffr, fft (modes > 0).

        Parameters
        ----------
        ffr0, fftheta0 : ndarray
            2d arrays of real values with shape [ w3d.nx+1, w3d.nz+1 ,
            self.fields.circ_m ], which contain the Fourier coefficients
            of the radial and azimuthal components of the field,
            for the mode 0.

        ffr, fftheta : ndarray
            3d arrays of complex values with shape [ w3d.nx+1, w3d.nz+1 ,
            self.fields.circ_m ], which contain the Fourier coefficients
            of the radial and azimuthal components of the field.

        name : str
            The name of the field considered, which is then printed on the plot.

        component : str
            The component of the field that is to be plotted : either 'x' or 'y'

        m : int
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots the sum o
            f the contribution of all modes.

        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis

        theta_plot : float
            Used only if direction = -1
            The angle (in radians) between the plane of the plot
            and the zx plane.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        titles : int
            A flag that indicates whether the plots should be
            automatically labeled

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        """

        if me == 0 :

            # Sets default direction to 1
            if ( direction in [-1,0,1,2] ) == False : direction = 1

            # Prepare sum over the modes
            if m is not None : # Single mode
                assert m<=self.fields.circ_m, "m > circ_m: %g > %g" \
                  %(m,self.fields.circ_m)
                modes = [m]
            else : # All modes
                modes = range( self.fields.circ_m+1 )


            # Sum all the selected modes, in the plane given by direction
            # - Radial component
            # (dump the returned theta array : it is given on the next line)
            fr, _ = self.sum_modes_circ( ffr0, ffr, direction=direction,
                                 modes=modes, theta_plot=theta_plot, iz=iz )
            # - Azimuthal component
            ftheta, theta = self.sum_modes_circ( fftheta0, fftheta,
                direction=direction, modes=modes, theta_plot=theta_plot, iz=iz)

            # Sum the radial and azimuthal components
            if component == 'x' :
                f = fr*cos(theta) - ftheta*sin(theta)
            elif component == 'y' :
                f = fr*sin(theta) + ftheta*cos(theta)

            # Set the boundaries (here X and Y represent the axes of the plot)
            # - If the plane of the plot includes the z axis
            if direction in [-1,0,1] :
                Xmin = w3d.zmmin + self.zgrid
                Xmax = w3d.zmmax + self.zgrid
                Ymin = -w3d.xmmax
                Ymax = w3d.xmmax
            # - If the plot is in the xy plane
            if direction == 2 :
                Xmin = -w3d.xmmax
                Xmax = w3d.xmmax
                Ymin = -w3d.xmmax
                Ymax = w3d.xmmax

            if show :
                # Plot function : colormap of the matrix f
                ppgeneric( f , xmin=Xmin, xmax=Xmax, ymin=Ymin,
                           ymax=Ymax, **kw )

                # Labeling of the plot
                if titles == 1 :
                    if m!=None :
                        name = name+' (m=%g)'%m
                    if direction == -1 :
                        ptitles(name,'Z','Coordinate along theta=%f' \
                                %theta_plot)
                    elif direction == 0 :
                        ptitles(name,'Z','Y')
                    elif direction == 1 :
                        ptitles(name,'Z','X')
                    elif direction == 2 :
                        ptitles(name,'X','Y')

            if output and me==0 :
                return( f.T )


    def pfrho( self, l_children=1, guards=0, direction=None, m=None,
               output=False, show=True, **kw) :
        """
        Plot the Rho charge density.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis
            (Used only with the azimuthal decomposition (l_2drz = 1))

        m : int
            Used only with the azimuthal decomposition (l_2drz = 1)
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots the sum
            of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        output, show : bool, optional
            whether to return the field matrix, and whether to
            plot it.

        Returns
        -------
        A 2d array representing the fields, if output is True

        The shape of the array is such that it can be directly
        plotted with matplotlib's imshow (non need to transpose it,
        or reverse one dimension)

        See the docstring of pfscalar_circ and genericpfem3d for further options
        """
        if self.l_getrho :
            if self.l_2drz == 1 : # Azimuthal decomposition
                # Gather the Fourier coefficients
                # Mode 0
                rho0 = self.gatherarray( self.getarray( self.fields.Rho,guards,overlap=True) )
                # Modes > 0
                if self.circ_m > 0 :
                    rho = self.gatherarray( \
                        self.getarray_circ( self.fields.Rho_circ,guards,overlap=True ) )
                else :
                    rho = None
                # Sum and plot them
                f = self.pfscalar_circ( rho0, rho, name='Rho^^', m=m,
                    direction=direction, output=output, show=show, **kw )

            else: # Cartesian fields
                f = self.getarray(self.fields.Rho,guards,
                    overlap=True)
                if show :
                    self.genericpfem3d( f,'Rho', direction=direction,**kw)

            if output and me==0 :
                return(f[::-1])

    def pfdivE( self, l_children=1, guards=0, direction=None, m=None, **kw) :
        """
        Plot the quantity epsilon_0*div(E)- rho (charge conservation)
        Only implemented in the azimuthal case.

        Main parameters :
        ------------
        direction : int
            A flag that defines in which plane the field is plotted.
            direction = 0 : zy plane
            direction = 1 : zx plane
            direction = 2 : xy plane
            direction = -1 : plane including the z axis and making
            an angle theta_plot with the x axis

        m : int
            The single mode whose contribution is to be plotted.
            If m is not provided (m=None), this function plots the sum
            of the contribution of all modes.

        iz : int
            Used only if direction = 2
            The z index of the xy slice that is to be plotted.

        See the docstring of pfscalar_circ and genericpfem3d for further options

        """
        im = 1.j  # Imaginary number im**2 = -1

        if self.l_getrho and self.l_2drz == 1 :
            # Gather the Fourier coefficients for the mode 0
            rho0 = self.gatherarray( self.getarray( self.fields.Rho,guards,overlap=True) )
            er0 = self.gatherarray( self.getarray( self.fields.Ex,guards,overlap=True ) )
            ez0 = self.gatherarray( self.getarray( self.fields.Ez,guards,overlap=True ) )
            # Calculate the divergence for the mode 0
            divE0 = zeros(rho0.shape)
            for j in range(1,len(er0)) :
                divE0[j,1:] = ( (1+0.5/j)*er0[j,1:] - \
                                (1-0.5/j)*er0[j-1,1:] )/w3d.dx + \
                                (ez0[j,1:]-ez0[j,:-1])/w3d.dz
            divE0[0,1:] = 4./w3d.dx*er0[0,1:] + (ez0[0,1:]-ez0[0,:-1])/w3d.dz
            if self.circ_m > 0 :
                # Gather the Fourier coefficients for the modes > 0
                rho=self.gatherarray(self.getarray_circ(self.fields.Rho_circ,guards,overlap=True))
                er=self.gatherarray(self.getarray_circ(self.fields.Ex_circ,guards,overlap=True))
                et=self.gatherarray(self.getarray_circ(self.fields.Ey_circ,guards,overlap=True))
                ez=self.gatherarray(self.getarray_circ(self.fields.Ez_circ,guards,overlap=True))
                divE = zeros(rho.shape, dtype='complex')
                # Calculate the divergence for the mode 0
                for M in range(self.circ_m) :
                    for j in range(1,len(rho)) :
                        divE[j,1:,M] =  ( (1+0.5/j)*er[j,1:,M]  \
                                        - (1-0.5/j)*er[j-1,1:,M] )/w3d.dx \
                                        - im*(M+1)*et[j,1:,M]/(j*w3d.dx) \
                                        + (ez[j,1:,M]-ez[j,:-1,M])/w3d.dz
                # Nothing for j=0, since the divergence of E is 0 on axis
                # for modes > 0
            else :
                divE = None
            # Sum and plot them
            self.pfscalar_circ( eps0*divE0-rho0, eps0*divE-rho,
                                name='eps0*div(E)-rho^^', m=m,
                                direction=direction, **kw )


    def sezax(self):
        pass

    def sphiax(self):
        pass

    def rhodia(self):
        pass

    def gtlchg(self):
        pass

    def srhoax(self):
        pass

    def getese(self):
        pass


    ######################################################
    #                    Getters                         #
    ######################################################


    def getjx(self,guards=0,overlap=0):
        return self.getarray(self.fields.Jx,guards,overlap)

    def getjy(self,guards=0,overlap=0):
        return self.getarray(self.fields.Jy,guards,overlap)

    def getjz(self,guards=0,overlap=0):
        return self.getarray(self.fields.Jz,guards,overlap)

    def getex(self,guards=0,overlap=0):
        return self.getarray(self.fields.Exp,guards,overlap)

    def getey(self,guards=0,overlap=0):
        return self.getarray(self.fields.Eyp,guards,overlap)

    def getez(self,guards=0,overlap=0):
        return self.getarray(self.fields.Ezp,guards,overlap)

    def getbx(self,guards=0,overlap=0):
        return self.getarray(self.fields.Bxp,guards,overlap)

    def getby(self,guards=0,overlap=0):
        return self.getarray(self.fields.Byp,guards,overlap)

    def getbz(self,guards=0,overlap=0):
        return self.getarray(self.fields.Bzp,guards,overlap)

    def getexg(self,guards=0,overlap=0):
        return self.getarray(self.fields.Ex,guards,overlap)

    def geteyg(self,guards=0,overlap=0):
        return self.getarray(self.fields.Ey,guards,overlap)

    def getezg(self,guards=0,overlap=0):
        return self.getarray(self.fields.Ez,guards,overlap)

    def getbxg(self,guards=0,overlap=0):
        return self.getarray(self.fields.Bx,guards,overlap)

    def getbyg(self,guards=0,overlap=0):
        return self.getarray(self.fields.By,guards,overlap)

    def getbzg(self,guards=0,overlap=0):
        return self.getarray(self.fields.Bz,guards,overlap)

    def getrho(self,guards=0,overlap=0):
        return self.getarray(self.fields.Rho,guards,overlap)

    def getrhop(self,guards=0,overlap=1):
        # a version of getrho with default overlap = 1, needed for particleinjection.py
        return self.getrho(guards,overlap)

    def getrhoold(self,guards=0,overlap=0):
        return self.getarray(self.fields.Rhoold,guards,overlap)

    def getf(self,guards=0,overlap=0):
        return self.getarray(self.fields.F,guards,overlap)

    def getg(self,guards=0,overlap=0):
        return self.getarray(self.fields.G,guards,overlap)

    def getex_circ(self,guards=0,overlap=0):
        return self.getarray_circ(self.fields.Ex_circ,guards,overlap)

    def getey_circ(self,guards=0,overlap=0):
        return self.getarray_circ(self.fields.Ey_circ,guards,overlap)

    def getez_circ(self,guards=0,overlap=0):
        return self.getarray_circ(self.fields.Ez_circ,guards,overlap)

    def getf_circ(self,guards=0,overlap=0):
        return self.getarray_circ(self.fields.F_circ,guards,overlap)

    def getjx_circ(self,guards=0,overlap=0):
        return self.getarray_circ(self.fields.Jx_circ,guards,overlap)

    def getjy_circ(self,guards=0,overlap=0):
        return self.getarray_circ(self.fields.Jy_circ,guards,overlap)

    def getjz_circ(self,guards=0,overlap=0):
        return self.getarray_circ(self.fields.Jz_circ,guards,overlap)

    def getrho_circ(self,guards=0,overlap=0):
        return self.getarray_circ(self.fields.Rho_circ,guards,overlap)

    def getincond(self,guards=0,overlap=0):
        return self.getarray(self.fields.incond,guards,overlap)

    def getdive(self,guards=0,overlap=0):
        dive = zeros(shape(self.fields.Ex),'d')
        f = self.fields
        if top.efetch[0] != 4:node2yee3d(f)
        getdive(f.Ex,f.Ey,f.Ez,dive,f.dx,f.dy,f.dz,
                f.nx,f.ny,f.nz,f.nxguard,f.nyguard,f.nzguard,
                f.xmin,
                self.l_2dxz,self.l_2drz,self.l_nodalgrid)
        if top.efetch[0] != 4:yee2node3d(f)
        return self.getarray(dive,guards,overlap)

    def getdiveold(self,guards=0,overlap=0):
        dive = zeros(shape(self.fields.Ex),'d')
        f = self.fields
        if top.efetch[0] != 4:node2yee3d(f)
        if self.l_2dxz:
            dive[1:-1,0,1:-1] = (f.Ex[1:-1,0,1:-1]-f.Ex[:-2,0,1:-1])/f.dx \
                             + (f.Ez[1:-1,0,1:-1]-f.Ez[1:-1,0,:-2])/f.dz
        else:
            dive[1:-1,1:-1,1:-1] = (f.Ex[1:-1,1:-1,1:-1]-f.Ex[:-2,1:-1,1:-1])/f.dx \
                                 + (f.Ey[1:-1,1:-1,1:-1]-f.Ey[1:-1,:-2,1:-1])/f.dy \
                                 + (f.Ez[1:-1,1:-1,1:-1]-f.Ez[1:-1,1:-1,:-2])/f.dz
        if top.efetch[0] != 4:yee2node3d(f)
        return dive

    def getdivj(self,guards=0,overlap=0):
        divj = zeros(shape(self.fields.Ex),'d')
        f = self.fields
        getdive(f.Jx,f.Jy,f.Jz,divj,f.dx,f.dy,f.dz,
                f.nx,f.ny,f.nz,f.nxguard,f.nyguard,f.nzguard,
                f.xmin,
                self.l_2dxz,self.l_2drz,self.l_nodalgrid)
        return self.getarray(divj,guards,overlap)

    def gete(self,guards=0,overlap=0):
        return self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap)

    def getb(self,guards=0,overlap=0):
        return self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap)

    def getw(self,guards=0,overlap=0):
        e2 = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap)
        b2 = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap)
        return sqrt(e2+clight**2*b2)

    def getw2(self,guards=0,overlap=0):
        e2 = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap)
        b2 = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap)
        return e2+clight**2*b2

    def getu(self,guards=0,overlap=0):
        e2 = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap)
        b2 = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap)
        u = e2+clight**2*b2
        if self.l_1dz:
            u*=self.dz*eps0/2
        elif self.l_2dxz:
            u*=self.dx*self.dz*eps0/2
        else:
            u*=self.dx*self.dy*self.dz*eps0/2
        return u

    def getug(self,guards=0,overlap=0):
        if not self.fields.l_nodecentered:
            self.yee2none3d()
            l_donode2yee=True
        else:
            l_donode2yee=False
        e2 = self.getarray(self.fields.Exp**2+self.fields.Eyp**2+self.fields.Ezp**2,guards,overlap)
        b2 = self.getarray(self.fields.Bxp**2+self.fields.Byp**2+self.fields.Bzp**2,guards,overlap)
        u = e2+clight**2*b2
        if self.l_1dz:
            u*=self.dz*eps0/2
        elif self.l_2dxz:
            u*=self.dx*self.dz*eps0/2
        else:
            u*=self.dx*self.dy*self.dz*eps0/2
        if l_donode2yee:
            self.node2yee3d()
        return u

    def geteg(self,guards=0,overlap=0):
        return self.getarray(self.fields.Ex**2+self.fields.Ey**2+self.fields.Ez**2,guards,overlap)

    def getbg(self,guards=0,overlap=0):
        return self.getarray(self.fields.Bx**2+self.fields.By**2+self.fields.Bz**2,guards,overlap)

    def getwg(self,guards=0,overlap=0):
        e2 = self.getarray(self.fields.Ex**2+self.fields.Ey**2+self.fields.Ez**2,guards,overlap)
        b2 = self.getarray(self.fields.Bx**2+self.fields.By**2+self.fields.Bz**2,guards,overlap)
        return sqrt(e2+clight**2*b2)

    def gets(self,guards=0,overlap=0):
        if not self.fields.l_nodecentered:
            self.yee2node3d()
            l_donode2yee=True
        else:
            l_donode2yee=False
        ex = self.getex(guards,overlap)
        ey = self.getey(guards,overlap)
        ez = self.getez(guards,overlap)
        bx = self.getbx(guards,overlap)
        by = self.getby(guards,overlap)
        bz = self.getbz(guards,overlap)
        sx = (ey*bz-ez*by)*eps0
        sy = (ez*bx-ex*bz)*eps0
        sz = (ex*by-ey*bx)*eps0
        if l_donode2yee:
            self.node2yee3d()
        return sx,sy,sz

    def getp(self,guards=0,overlap=0):
        sx,sy,sz=self.gets(guards,overlap)
        if self.l_1dz:
            sx*=self.dz
            sy*=self.dz
            sz*=self.dz
        elif self.l_2dxz:
            sx*=self.dx*self.dz
            sy*=self.dx*self.dz
            sz*=self.dx*self.dz
        else:
            sx*=self.dx*self.dy*self.dz
            sy*=self.dx*self.dy*self.dz
            sz*=self.dx*self.dy*self.dz
        return sx,sy,sz

    def getwg2(self,guards=0,overlap=0):
        e2 = self.getarray(self.fields.Ex**2+self.fields.Ey**2+self.fields.Ez**2,guards,overlap)
        b2 = self.getarray(self.fields.Bx**2+self.fields.By**2+self.fields.Bz**2,guards,overlap)
        return e2+clight**2*b2

    def gatherexg(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getexg(guards),direction=direction,**kw)

    def gathereyg(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.geteyg(guards),direction=direction,**kw)

    def gatherezg(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getezg(guards),direction=direction,**kw)

    def gatherbxg(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getbxg(guards),direction=direction,**kw)

    def gatherbyg(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getbyg(guards),direction=direction,**kw)

    def gatherbzg(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getbzg(guards),direction=direction,**kw)

    def gatherex(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getex(guards),direction=direction,**kw)

    def gatherey(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getey(guards),direction=direction,**kw)

    def gatherez(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getez(guards),direction=direction,**kw)

    def gatherbx(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getbx(guards),direction=direction,**kw)

    def gatherby(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getby(guards),direction=direction,**kw)

    def gatherbz(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getbz(guards),direction=direction,**kw)

    def gatherjx(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getjx(guards),direction=direction,**kw)

    def gatherjy(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getjy(guards),direction=direction,**kw)

    def gatherjz(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getjz(guards),direction=direction,**kw)

    def gatherrho(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getrho(guards),direction=direction,**kw)

    def gatherrhoold(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getrhoold(guards),direction=direction,**kw)

    def gatherf(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getf(guards),direction=direction,**kw)

    def gatherincond(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getincond(guards),direction=direction,**kw)

    def gatherdive(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getdive(guards),direction=direction,**kw)

    def gatheru(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getu(guards,overlap=direction is None),direction=direction,**kw)

    def gatherg(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getg(guards,overlap=direction is None),direction=direction,**kw)

    def gathers(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.gets(guards,overlap=direction is None),direction=direction,**kw)

    def gatherp(self,guards=0,direction=None,**kw):
        return self.gatherarray(self.getp(guards,overlap=direction is None),direction=direction,**kw)

    def get_tot_energy(self):
#    yee2node3d(self.fields)
        if not self.l_2drz:
            w2tot = globalsum(self.getwg2()) # must be getwg2, not getw2
            if self.l_1dz:
                w2tot*=self.dz*eps0/2
            elif self.l_2dxz:
                w2tot*=self.dx*self.dz*eps0/2
            else:
                w2tot*=self.dx*self.dy*self.dz*eps0/2
        else:
            w2 = self.getw2()
            w2=transpose(w2)*self.vol[:shape(w2)[0]]
            w2tot = globalsum(w2)*eps0/2

        return w2tot

    def sumdive(self):
        flist = [self.block.core.yf]
        if self.refinement is not None:
            flist += [self.field_coarse.block.core.yf]
        q = []
        for f in flist:
            if top.efetch[0] != 4:node2yee3d(f)
            Ex = f.Ex#self.getarray(f.Ex)
            Ey = f.Ey#self.getarray(f.Ey)
            Ez = f.Ez#self.getarray(f.Ez)
            if self.l_1dz:
                q.append( (sum(Ez[1:-1,0,1:-1]-Ez[1:-1,0,:-2])/f.dz) \
                        * eps0*(f.dx*f.dz))
            elif self.l_2dxz:
                q.append( (sum(Ex[1:-1,0,1:-1]-Ex[:-2,0,1:-1])/f.dx \
                        +  sum(Ez[1:-1,0,1:-1]-Ez[1:-1,0,:-2])/f.dz) \
                        * eps0*(f.dx*f.dz))
            else:
                q.append( (sum(Ex[1:-1,1:-1,1:-1]-Ex[:-2,1:-1,1:-1])/f.dx \
                        +  sum(Ey[1:-1,1:-1,1:-1]-Ey[1:-1,:-2,1:-1])/f.dy \
                        +  sum(Ez[1:-1,1:-1,1:-1]-Ez[1:-1,1:-1,:-2])/f.dz) \
                        * eps0*(f.dx*f.dy*f.dz))
            if top.efetch[0] != 4:yee2node3d(f)
        return q

    def sumq(self):
        flist = [self.block.core.yf]
        if self.refinement is not None:
            flist += [self.field_coarse.block.core.yf]
        q = []
        for f in flist:
            if self.l_1dz:
                q.append(sum(f.Rho)*f.dz)
            elif self.l_2dxz:
                q.append(sum(f.Rho)*(f.dx*f.dz))
            else:
                q.append(sum(f.Rho)*(f.dx*f.dy*f.dz))
        return q

    def setbcparallel(self,dir):
        if dir==0:nprocs = top.nxprocs
        if dir==1:nprocs = top.nyprocs
        if dir==2:nprocs = top.nzprocs
        ib = 2*dir
        if nprocs>1:
            down = top.procneighbors[0,dir]
            up   = top.procneighbors[1,dir]
            # --- check lower bound in z
            if self.bounds[ib+1] == periodic:
                condu = up != me
            else:
                condu = up>me
            if self.bounds[ib] == periodic:
                condd = down != me
            else:
                condd = down<me
            if condu:
                mpisend(self.isactive, dest = up)
            if condd:
                isactive = mpirecv(source = down)
                if isactive and self.isactive:self.bounds[ib]=em3d.otherproc
            # --- check upper bound in z
            if condd:
                mpisend(self.isactive, dest = down)
            if condu:
                isactive = mpirecv(source = up)
                if isactive and self.isactive:self.bounds[ib+1]=em3d.otherproc

    def connectblocks(self,blo,bup,dir):
        if dir==0:
            # --- X
            loxropen = blo.xrbnd == openbc
            loylopen = blo.ylbnd == openbc
            loyropen = blo.yrbnd == openbc
            lozlopen = blo.zlbnd == openbc
            lozropen = blo.zrbnd == openbc

            upxlopen = bup.xlbnd == openbc
            upylopen = bup.ylbnd == openbc
            upyropen = bup.yrbnd == openbc
            upzlopen = bup.zlbnd == openbc
            upzropen = bup.zrbnd == openbc

            # --- sides
            if loxropen: del blo.sidexr.syf # --- deallocate upper x PML block
            if upxlopen: del bup.sidexl.syf # --- deallocate lower x PML block
            blo.sidexr = bup.core
#      bup.sidexl = blo.core
            bup.sidexl = EM3D_FIELDtype()
            # --- edges
            if loxropen and loylopen:del blo.edgexryl.syf
            if loxropen and loyropen:del blo.edgexryr.syf
            if loxropen and lozlopen:del blo.edgexrzl.syf
            if loxropen and lozropen:del blo.edgexrzr.syf
            if upxlopen and upylopen:del bup.edgexlyl.syf
            if upxlopen and upyropen:del bup.edgexlyr.syf
            if upxlopen and upzlopen:del bup.edgexlzl.syf
            if upxlopen and upzropen:del bup.edgexlzr.syf
            if loylopen:
                blo.edgexryl = bup.sideyl
            else:
                blo.edgexryl = EM3D_FIELDtype()
            if loyropen:
                blo.edgexryr = bup.sideyr
            else:
                blo.edgexryr = EM3D_FIELDtype()
            if lozlopen:
                blo.edgexrzl = bup.sidezl
            else:
                blo.edgexrzl = EM3D_FIELDtype()
            if lozropen:
                blo.edgexrzr = bup.sidezr
            else:
                blo.edgexrzr = EM3D_FIELDtype()
#      bup.edgexlyl = blo.sideyl
#      bup.edgexlyr = blo.sideyr
#      bup.edgexlzl = blo.sidezl
#      bup.edgexlzr = blo.sidezr
            bup.edgexlyl = EM3D_FIELDtype()
            bup.edgexlyr = EM3D_FIELDtype()
            bup.edgexlzl = EM3D_FIELDtype()
            bup.edgexlzr = EM3D_FIELDtype()
            # --- corners
            if loxropen and loylopen and lozlopen:del blo.cornerxrylzl.syf
            if loxropen and loyropen and lozlopen:del blo.cornerxryrzl.syf
            if loxropen and loylopen and lozropen:del blo.cornerxrylzr.syf
            if loxropen and loyropen and lozropen:del blo.cornerxryrzr.syf
            if upxlopen and upylopen and upzlopen:del bup.cornerxlylzl.syf
            if upxlopen and upyropen and upzlopen:del bup.cornerxlyrzl.syf
            if upxlopen and upylopen and upzropen:del bup.cornerxlylzr.syf
            if upxlopen and upyropen and upzropen:del bup.cornerxlyrzr.syf
            if loylopen and lozlopen:
                blo.cornerxrylzl = bup.edgeylzl
            else:
                blo.cornerxrylzl = EM3D_FIELDtype()
            if loyropen and lozlopen:
                blo.cornerxryrzl = bup.edgeyrzl
            else:
                blo.cornerxryrzl = EM3D_FIELDtype()
            if loylopen and lozropen:
                blo.cornerxrylzr = bup.edgeylzr
            else:
                blo.cornerxrylzr = EM3D_FIELDtype()
            if loyropen and lozropen:
                blo.cornerxryrzr = bup.edgeyrzr
            else:
                blo.cornerxryrzr = EM3D_FIELDtype()
#      bup.cornerxlylzl = blo.edgeylzl
#      bup.cornerxlyrzl = blo.edgeyrzl
#      bup.cornerxlylzr = blo.edgeylzr
#      bup.cornerxlyrzr = blo.edgeyrzr
            bup.cornerxlylzl = EM3D_FIELDtype()
            bup.cornerxlyrzl = EM3D_FIELDtype()
            bup.cornerxlylzr = EM3D_FIELDtype()
            bup.cornerxlyrzr = EM3D_FIELDtype()

            blo.xrbnd = em3d.otherblock
            bup.xlbnd = em3d.otherblock

        if dir==1:
            # --- Y
            loyropen = blo.yrbnd == openbc
            lozlopen = blo.zlbnd == openbc
            lozropen = blo.zrbnd == openbc
            loxlopen = blo.xlbnd == openbc
            loxropen = blo.xrbnd == openbc

            upylopen = bup.ylbnd == openbc
            upzlopen = bup.zlbnd == openbc
            upzropen = bup.zrbnd == openbc
            upxlopen = bup.xlbnd == openbc
            upxropen = bup.xrbnd == openbc

            # --- sides
            if loyropen: del blo.sideyr.syf # --- deallocate upper y PML block
            if upylopen: del bup.sideyl.syf # --- deallocate lower y PML block
            blo.sideyr = bup.core
#      bup.sideyl = blo.core
            bup.sideyl = EM3D_FIELDtype()
            # --- edges
            if loyropen and lozlopen:del blo.edgeyrzl.syf
            if loyropen and lozropen:del blo.edgeyrzr.syf
            if loyropen and loxlopen:del blo.edgexlyr.syf
            if loyropen and loxropen:del blo.edgexryr.syf
            if upylopen and upzlopen:del bup.edgeylzl.syf
            if upylopen and upzropen:del bup.edgeylzr.syf
            if upylopen and upxlopen:del bup.edgexlyl.syf
            if upylopen and upxropen:del bup.edgexryr.syf
            if loxlopen:
                blo.edgexlyr = bup.sidexl
            else:
                blo.edgexlyr = EM3D_FIELDtype()
            if loxropen:
                blo.edgexryr = bup.sidexr
            else:
                blo.edgexryr = EM3D_FIELDtype()
            if lozlopen:
                blo.edgeyrzl = bup.sidezl
            else:
                blo.edgeyrzl = EM3D_FIELDtype()
            if lozropen:
                blo.edgeyrzr = bup.sidezr
            else:
                blo.edgeyrzr = EM3D_FIELDtype()
#      bup.edgexlyl = blo.sidexl
#      bup.edgexryl = blo.sidexr
#      bup.edgeylzl = blo.sidezl
#      bup.edgeylzr = blo.sidezr
            bup.edgexlyl = EM3D_FIELDtype()
            bup.edgexryl = EM3D_FIELDtype()
            bup.edgeylzl = EM3D_FIELDtype()
            bup.edgeylzr = EM3D_FIELDtype()
            # --- corners
            if loyropen and lozlopen and loxlopen:del blo.cornerxlyrzl.syf
            if loyropen and lozropen and loxlopen:del blo.cornerxlyrzr.syf
            if loyropen and lozlopen and loxropen:del blo.cornerxryrzl.syf
            if loyropen and lozropen and loxropen:del blo.cornerxryrzr.syf
            if upylopen and upzlopen and upxlopen:del bup.cornerxlylzl.syf
            if upylopen and upzropen and upxlopen:del bup.cornerxlylzr.syf
            if upylopen and upzlopen and upxropen:del bup.cornerxrylzl.syf
            if upylopen and upzropen and upxropen:del bup.cornerxrylzr.syf
            if loxlopen and lozlopen:
                blo.cornerxlyrzl = bup.edgexlzl
            else:
                blo.cornerxlyrzl = EM3D_FIELDtype()
            if loxropen and lozlopen:
                blo.cornerxryrzl = bup.edgexrzl
            else:
                blo.cornerxryrzl = EM3D_FIELDtype()
            if loxlopen and lozropen:
                blo.cornerxlyrzr = bup.edgexlzr
            else:
                blo.cornerxlyrzr = EM3D_FIELDtype()
            if loxropen and lozropen:
                blo.cornerxryrzr = bup.edgexrzr
            else:
                blo.cornerxryrzr = EM3D_FIELDtype()
#      bup.cornerxlylzl = blo.edgexlzl
#      bup.cornerxrylzl = blo.edgexrzl
#      bup.cornerxlylzr = blo.edgexlzr
#      bup.cornerxrylzr = blo.edgexrzr
            bup.cornerxlylzl = EM3D_FIELDtype()
            bup.cornerxrylzl = EM3D_FIELDtype()
            bup.cornerxlylzr = EM3D_FIELDtype()
            bup.cornerxrylzr = EM3D_FIELDtype()

            blo.yrbnd = em3d.otherblock
            bup.ylbnd = em3d.otherblock

        if dir==2:
            # --- Z
            lozropen = blo.zrbnd == openbc
            loxlopen = blo.xlbnd == openbc
            loxropen = blo.xrbnd == openbc
            loylopen = blo.ylbnd == openbc
            loyropen = blo.yrbnd == openbc

            upzlopen = bup.zlbnd == openbc
            upxlopen = bup.xlbnd == openbc
            upxropen = bup.xrbnd == openbc
            upylopen = bup.ylbnd == openbc
            upyropen = bup.yrbnd == openbc

            # --- sides
            if lozropen: del blo.sidezr.syf # --- deallocate upper y PML block
            if upzlopen: del bup.sidezl.syf # --- deallocate lower y PML block
            blo.sidezr = bup.core
#      bup.sidezl = blo.core
            bup.sidezl = EM3D_FIELDtype()
            # --- edges
            if lozropen and loxlopen:del blo.edgexlzr.syf
            if lozropen and loxropen:del blo.edgexrzr.syf
            if lozropen and loylopen:del blo.edgeylzr.syf
            if lozropen and loyropen:del blo.edgeyrzr.syf
            if upzlopen and upxlopen:del bup.edgexlzl.syf
            if upzlopen and upxropen:del bup.edgexrzl.syf
            if upzlopen and upylopen:del bup.edgeylzl.syf
            if upzlopen and upyropen:del bup.edgeyrzl.syf
            if loxlopen:
                blo.edgexlzr = bup.sidexl
            else:
                blo.edgexlzr = EM3D_FIELDtype()
            if loxropen:
                blo.edgexrzr = bup.sidexr
            else:
                blo.edgexrzr = EM3D_FIELDtype()
            if loylopen:
                blo.edgeylzr = bup.sideyl
            else:
                blo.edgeylzr = EM3D_FIELDtype()
            if loyropen:
                blo.edgeyrzr = bup.sideyr
            else:
                blo.edgeyrzr = EM3D_FIELDtype()
         #     bup.edgexlzl = blo.sidexl
         #     bup.edgexrzl = blo.sidexr
         #     bup.edgeylzl = blo.sideyl
         #     bup.edgeyrzl = blo.sideyr
            bup.edgexlzl = EM3D_FIELDtype()
            bup.edgexrzl = EM3D_FIELDtype()
            bup.edgeylzl = EM3D_FIELDtype()
            bup.edgeyrzl = EM3D_FIELDtype()
            # --- corners
            if lozropen and loxlopen and loylopen:del blo.cornerxlylzr.syf
            if lozropen and loxropen and loylopen:del blo.cornerxrylzr.syf
            if lozropen and loxlopen and loyropen:del blo.cornerxlyrzr.syf
            if lozropen and loxropen and loyropen:del blo.cornerxryrzr.syf
            if upzlopen and upxlopen and upylopen:del bup.cornerxlylzl.syf
            if upzlopen and upxropen and upylopen:del bup.cornerxrylzl.syf
            if upzlopen and upxlopen and upyropen:del bup.cornerxlyrzl.syf
            if upzlopen and upxropen and upyropen:del bup.cornerxryrzl.syf
            if loxlopen and loylopen:
                blo.cornerxlylzr = bup.edgexlyl
            else:
                blo.cornerxlylzr = EM3D_FIELDtype()
            if loxropen and loylopen:
                blo.cornerxrylzr = bup.edgexryl
            else:
                blo.cornerxrylzr = EM3D_FIELDtype()
            if loxlopen and loyropen:
                blo.cornerxlyrzr = bup.edgexlyr
            else:
                blo.cornerxlyrzr = EM3D_FIELDtype()
            if loxropen and loyropen:
                blo.cornerxryrzr = bup.edgexryr
            else:
                blo.cornerxryrzr = EM3D_FIELDtype()
           #   bup.cornerxlylzl = blo.edgexlyl
           #   bup.cornerxrylzl = blo.edgexryl
           #   bup.cornerxlyrzl = blo.edgexlyr
           #   bup.cornerxryrzl = blo.edgexryr
            bup.cornerxlylzl = EM3D_FIELDtype()
            bup.cornerxrylzl = EM3D_FIELDtype()
            bup.cornerxlyrzl = EM3D_FIELDtype()
            bup.cornerxryrzl = EM3D_FIELDtype()

            blo.zrbnd = em3d.otherblock
            bup.zlbnd = em3d.otherblock

    def setdtinit(self):
        self.dtinit=top.dt

    def setbcoverlaps(self):
        xlbndmatch = False
        xrbndmatch = False
        ylbndmatch = False
        yrbndmatch = False
        zlbndmatch = False
        zrbndmatch = False
        for n in self.overlapshigher:
            block = self.root.listofblocks[n].block
            blockcoarse = self.root.listofblocks[n].field_coarse.block
            lower = self.overlapshigher[n][0]
            upper = self.overlapshigher[n][1]
            # --- check grid matching in upper x
            if self.fullupper[0] == lower[0] and \
               all(self.fulllower[1:] == lower[1:]) and \
               all(self.fullupper[1:] == upper[1:]):
                if not (self.block.sidexr is block):
                    self.connectblocks(self.block,block,0)
                    self.connectblocks(self.field_coarse.block,blockcoarse,0)
                xrbndmatch = True
            # --- check grid matching in upper y
            if self.fullupper[1] == lower[1] and \
               all(self.fulllower[0::2] == lower[0::2]) and \
               all(self.fullupper[0::2] == upper[0::2]):
                if not (self.block.sideyr is block):
                    self.connectblocks(self.block,block,1)
                    self.connectblocks(self.field_coarse.block,blockcoarse,1)
                yrbndmatch = True
            # --- check grid matching in upper z
            if self.fullupper[2] == lower[2] and \
               all(self.fulllower[:2] == lower[:2]) and \
               all(self.fullupper[:2] == upper[:2]):
                if not (self.block.sidezr is block):
                    self.connectblocks(self.block,block,2)
                    self.connectblocks(self.field_coarse.block,blockcoarse,2)
                zrbndmatch = True
        for n in self.overlapslower:
            block = self.root.listofblocks[n].block
            blockcoarse = self.root.listofblocks[n].field_coarse.block
            lower = self.overlapslower[n][0]
            upper = self.overlapslower[n][1]
            # --- check grid matching in lower x
            if self.fulllower[0] == upper[0] and \
               all(self.fulllower[1:] == lower[1:]) and \
               all(self.fullupper[1:] == upper[1:]):
                if not (self.block.sidexl is block):
                    self.connectblocks(block,self.block,0)
                    self.connectblocks(blockcoarse,self.field_coarse.block,0)
                xlbndmatch = True
            # --- check grid matching in lower y
            if self.fulllower[1] == upper[1] and \
               all(self.fulllower[0::2] == lower[0::2]) and \
               all(self.fullupper[0::2] == upper[0::2]):
                if not (self.block.sideyl is block):
                    self.connectblocks(block,self.block,1)
                    self.connectblocks(blockcoarse,self.field_coarse.block,1)
                ylbndmatch = True
            # --- check grid matching in lower z
            if self.fulllower[2] == upper[2] and \
               all(self.fulllower[:2] == lower[:2]) and \
               all(self.fullupper[:2] == upper[:2]):
                if not (self.block.sidezl is block):
                    self.connectblocks(block,self.block,2)
                    self.connectblocks(blockcoarse,self.field_coarse.block,2)
                zlbndmatch = True

        if not xlbndmatch and self.block.xlbnd == em3d.otherblock:
            self.block.xlbnd = openbc
            self.field_coarse.block.xlbnd = openbc
            # --- TODO: allocate openbc
        if not xrbndmatch and self.block.xrbnd == em3d.otherblock:
            self.block.xrbnd = openbc
            self.field_coarse.block.xrbnd = openbc
            # --- TODO: allocate openbc
        if not ylbndmatch and self.block.ylbnd == em3d.otherblock:
            self.block.ylbnd = openbc
            self.field_coarse.block.ylbnd = openbc
            # --- TODO: allocate openbc
        if not yrbndmatch and self.block.yrbnd == em3d.otherblock:
            self.block.yrbnd = openbc
            self.field_coarse.block.yrbnd = openbc
            # --- TODO: allocate openbc
        if not zlbndmatch and self.block.zlbnd == em3d.otherblock:
            self.block.zlbnd = openbc
            self.field_coarse.block.zlbnd = openbc
            # --- TODO: allocate openbc
        if not zrbndmatch and self.block.zrbnd == em3d.otherblock:
            self.block.zrbnd = openbc
            self.field_coarse.block.zrbnd = openbc
            # --- TODO: allocate openbc

        self.bounds = array([self.block.xlbnd,
                             self.block.xrbnd,
                             self.block.ylbnd,
                             self.block.yrbnd,
                             self.block.zlbnd,
                             self.block.zrbnd])
        if self.refinement is not None:
            self.field_coarse.bounds = array([self.block.xlbnd,
                                 self.block.xrbnd,
                                 self.block.ylbnd,
                                 self.block.yrbnd,
                                 self.block.zlbnd,
                                 self.block.zrbnd])

    def checkconnections(self):
        yeefieldtype = -1
        if self.block.xlbnd==openbc:
            if self.block.edgexlyl.fieldtype==yeefieldtype:
                self.block.edgexlyl=EM3D_FIELDtype()
                self.block.cornerxlylzl=EM3D_FIELDtype()
                self.block.cornerxlylzr=EM3D_FIELDtype()
            if self.block.edgexlyr.fieldtype==yeefieldtype:
                self.block.edgexlyr=EM3D_FIELDtype()
                self.block.cornerxlyrzl=EM3D_FIELDtype()
                self.block.cornerxlyrzr=EM3D_FIELDtype()
            if self.block.edgexlzl.fieldtype==yeefieldtype:
                self.block.edgexlzl=EM3D_FIELDtype()
                self.block.cornerxlylzl=EM3D_FIELDtype()
                self.block.cornerxlyrzl=EM3D_FIELDtype()
            if self.block.edgexlzr.fieldtype==yeefieldtype:
                self.block.edgexlzr=EM3D_FIELDtype()
                self.block.cornerxlylzr=EM3D_FIELDtype()
                self.block.cornerxlyrzr=EM3D_FIELDtype()
        if self.block.xrbnd==openbc:
            if self.block.edgexryl.fieldtype==yeefieldtype:
                self.block.edgexryl=EM3D_FIELDtype()
                self.block.cornerxrylzl=EM3D_FIELDtype()
                self.block.cornerxrylzr=EM3D_FIELDtype()
            if self.block.edgexryr.fieldtype==yeefieldtype:
                self.block.edgexryr=EM3D_FIELDtype()
                self.block.cornerxryrzl=EM3D_FIELDtype()
                self.block.cornerxryrzr=EM3D_FIELDtype()
            if self.block.edgexrzl.fieldtype==yeefieldtype:
                self.block.edgexrzl=EM3D_FIELDtype()
                self.block.cornerxrylzl=EM3D_FIELDtype()
                self.block.cornerxryrzl=EM3D_FIELDtype()
            if self.block.edgexrzr.fieldtype==yeefieldtype:
                self.block.edgexrzr=EM3D_FIELDtype()
                self.block.cornerxrylzr=EM3D_FIELDtype()
                self.block.cornerxryrzr=EM3D_FIELDtype()

        if self.block.ylbnd==openbc:
            if self.block.edgexlyl.fieldtype==yeefieldtype:
                self.block.edgexlyl=EM3D_FIELDtype()
                self.block.cornerxlylzl=EM3D_FIELDtype()
                self.block.cornerxlylzr=EM3D_FIELDtype()
            if self.block.edgexryl.fieldtype==yeefieldtype:
                self.block.edgexryl=EM3D_FIELDtype()
                self.block.cornerxrylzl=EM3D_FIELDtype()
                self.block.cornerxrylzr=EM3D_FIELDtype()
            if self.block.edgeylzl.fieldtype==yeefieldtype:
                self.block.edgeylzl=EM3D_FIELDtype()
                self.block.cornerxlylzl=EM3D_FIELDtype()
                self.block.cornerxrylzl=EM3D_FIELDtype()
            if self.block.edgeylzr.fieldtype==yeefieldtype:
                self.block.edgeylzr=EM3D_FIELDtype()
                self.block.cornerxlylzr=EM3D_FIELDtype()
                self.block.cornerxrylzr=EM3D_FIELDtype()
        if self.block.yrbnd==openbc:
            if self.block.edgexlyr.fieldtype==yeefieldtype:
                self.block.edgexlyr=EM3D_FIELDtype()
                self.block.cornerxlyrzl=EM3D_FIELDtype()
                self.block.cornerxlyrzr=EM3D_FIELDtype()
            if self.block.edgexryr.fieldtype==yeefieldtype:
                self.block.edgexryr=EM3D_FIELDtype()
                self.block.cornerxryrzl=EM3D_FIELDtype()
                self.block.cornerxryrzr=EM3D_FIELDtype()
            if self.block.edgeyrzl.fieldtype==yeefieldtype:
                self.block.edgeyrzl=EM3D_FIELDtype()
                self.block.cornerxlyrzl=EM3D_FIELDtype()
                self.block.cornerxryrzl=EM3D_FIELDtype()
            if self.block.edgeyrzr.fieldtype==yeefieldtype:
                self.block.edgeyrzr=EM3D_FIELDtype()
                self.block.cornerxlyrzr=EM3D_FIELDtype()
                self.block.cornerxryrzr=EM3D_FIELDtype()

        if self.block.zlbnd==openbc:
            if self.block.edgexlzl.fieldtype==yeefieldtype:
                self.block.edgexlzl=EM3D_FIELDtype()
                self.block.cornerxlylzl=EM3D_FIELDtype()
                self.block.cornerxlyrzl=EM3D_FIELDtype()
            if self.block.edgexrzl.fieldtype==yeefieldtype:
                self.block.edgexrzl=EM3D_FIELDtype()
                self.block.cornerxrylzl=EM3D_FIELDtype()
                self.block.cornerxryrzl=EM3D_FIELDtype()
            if self.block.edgeylzl.fieldtype==yeefieldtype:
                self.block.edgeylzl=EM3D_FIELDtype()
                self.block.cornerxlylzl=EM3D_FIELDtype()
                self.block.cornerxrylzl=EM3D_FIELDtype()
            if self.block.edgeyrzl.fieldtype==yeefieldtype:
                self.block.edgeyrzl=EM3D_FIELDtype()
                self.block.cornerxlyrzl=EM3D_FIELDtype()
                self.block.cornerxryrzl=EM3D_FIELDtype()
        if self.block.zrbnd==openbc:
            if self.block.edgexlzr.fieldtype==yeefieldtype:
                self.block.edgexlzr=EM3D_FIELDtype()
                self.block.cornerxlylzr=EM3D_FIELDtype()
                self.block.cornerxlyrzr=EM3D_FIELDtype()
            if self.block.edgexrzr.fieldtype==yeefieldtype:
                self.block.edgexrzr=EM3D_FIELDtype()
                self.block.cornerxrylzr=EM3D_FIELDtype()
                self.block.cornerxryrzr=EM3D_FIELDtype()
            if self.block.edgeylzr.fieldtype==yeefieldtype:
                self.block.edgeylzr=EM3D_FIELDtype()
                self.block.cornerxlylzr=EM3D_FIELDtype()
                self.block.cornerxrylzr=EM3D_FIELDtype()
            if self.block.edgeyrzr.fieldtype==yeefieldtype:
                self.block.edgeyrzr=EM3D_FIELDtype()
                self.block.cornerxlyrzr=EM3D_FIELDtype()
                self.block.cornerxryrzr=EM3D_FIELDtype()

        if self.block.xlbnd==openbc and  self.block.ylbnd==openbc:
            if self.block.cornerxlylzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
            if self.block.cornerxlylzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
        if self.block.xrbnd==openbc and  self.block.ylbnd==openbc:
            if self.block.cornerxrylzl.fieldtype==yeefieldtype:self.block.cornerxrylzl=EM3D_FIELDtype()
            if self.block.cornerxrylzr.fieldtype==yeefieldtype:self.block.cornerxrylzr=EM3D_FIELDtype()
        if self.block.xlbnd==openbc and  self.block.yrbnd==openbc:
            if self.block.cornerxlyrzl.fieldtype==yeefieldtype:self.block.cornerxlyrzl=EM3D_FIELDtype()
            if self.block.cornerxlyrzr.fieldtype==yeefieldtype:self.block.cornerxlyrzr=EM3D_FIELDtype()
        if self.block.xrbnd==openbc and  self.block.yrbnd==openbc:
            if self.block.cornerxryrzl.fieldtype==yeefieldtype:self.block.cornerxryrzl=EM3D_FIELDtype()
            if self.block.cornerxryrzr.fieldtype==yeefieldtype:self.block.cornerxryrzr=EM3D_FIELDtype()

        if self.block.xlbnd==openbc and  self.block.zlbnd==openbc:
            if self.block.cornerxlylzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
            if self.block.cornerxlyrzl.fieldtype==yeefieldtype:self.block.cornerxlyrzl=EM3D_FIELDtype()
        if self.block.xrbnd==openbc and  self.block.zlbnd==openbc:
            if self.block.cornerxrylzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
            if self.block.cornerxryrzl.fieldtype==yeefieldtype:self.block.cornerxlyrzl=EM3D_FIELDtype()
        if self.block.xlbnd==openbc and  self.block.zrbnd==openbc:
            if self.block.cornerxlylzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
            if self.block.cornerxlyrzr.fieldtype==yeefieldtype:self.block.cornerxlyrzr=EM3D_FIELDtype()
        if self.block.xrbnd==openbc and  self.block.zrbnd==openbc:
            if self.block.cornerxrylzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
            if self.block.cornerxryrzr.fieldtype==yeefieldtype:self.block.cornerxlyrzr=EM3D_FIELDtype()

        if self.block.ylbnd==openbc and  self.block.zlbnd==openbc:
            if self.block.cornerxlylzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
            if self.block.cornerxrylzl.fieldtype==yeefieldtype:self.block.cornerxrylzl=EM3D_FIELDtype()
        if self.block.yrbnd==openbc and  self.block.zlbnd==openbc:
            if self.block.cornerxlyrzl.fieldtype==yeefieldtype:self.block.cornerxlylzl=EM3D_FIELDtype()
            if self.block.cornerxryrzl.fieldtype==yeefieldtype:self.block.cornerxrylzl=EM3D_FIELDtype()
        if self.block.ylbnd==openbc and  self.block.zrbnd==openbc:
            if self.block.cornerxlylzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
            if self.block.cornerxrylzr.fieldtype==yeefieldtype:self.block.cornerxrylzr=EM3D_FIELDtype()
        if self.block.yrbnd==openbc and  self.block.zrbnd==openbc:
            if self.block.cornerxlyrzr.fieldtype==yeefieldtype:self.block.cornerxlylzr=EM3D_FIELDtype()
            if self.block.cornerxryrzr.fieldtype==yeefieldtype:self.block.cornerxrylzr=EM3D_FIELDtype()

        if self.refinement is not None:
            self.field_coarse.checkconnections()

    def fillchilddomains(self):
        parent = self.root.listofblocks[self.parents[0]]
        xl,yl,zl = self.fullloweroverrefinement-parent.fulllower
#    xu,yu,zu = self.fullupperoverrefinement-parent.fulllower+[1,1,1]
        xu,yu,zu = self.fullupperoverrefinement-parent.fulllower+[0,0,0]

        if self.bounds[0]==openbc:
            xlguard = self.nguard[0]
            xlguarddepos = self.nguarddepos[0]
        else:
            xlguard = 0
            xlguarddepos = 0
        if self.bounds[1]==openbc:
            xrguard = self.nguard[0]
            xrguarddepos = self.nguarddepos[0]
        else:
            xrguard = 0
            xrguarddepos = 0

        if self.bounds[2]==openbc:
            ylguard = self.nguard[1]
            ylguarddepos = self.nguarddepos[1]
        else:
            ylguard = 0
            ylguarddepos = 0
        if self.bounds[3]==openbc:
            yrguard = self.nguard[1]
            yrguarddepos = self.nguarddepos[1]
        else:
            yrguard = 0
            yrguarddepos = 0

        if self.bounds[4]==openbc:
            zlguard = self.nguard[2]
            zlguarddepos = self.nguarddepos[2]
        else:
            zlguard = 0
            zlguarddepos = 0
        if self.bounds[5]==openbc:
            zrguard = self.nguard[2]
            zrguarddepos = self.nguarddepos[2]
        else:
            zrguard = 0
            zrguarddepos = 0

        cd = parent.childdomains
        cd[xl:xu,yl:yu,zl:zu]=parent.blocknumber
        cd[xl+xlguarddepos:xu-xrguarddepos,
           yl+ylguarddepos:yu-yrguarddepos,
           zl+zlguarddepos:zu-zrguarddepos]=-self.blocknumber
        cd[xl+xlguard:xu-xrguard,
           yl+ylguard:yu-yrguard,
           zl+zlguard:zu-zrguard]=self.blocknumber

    def fillchilddomainsold(self):
        xl,yl,zl = self.fullloweroverrefinement
        xu,yu,zu = self.fullupperoverrefinement+[1,1,1]

        if self.bounds[0]==openbc:
            xlguard = self.nguard[0]
            xlguarddepos = self.nguarddepos[0]
        else:
            xlguard = 0
            xlguarddepos = 0
        if self.bounds[1]==openbc:
            xrguard = self.nguard[0]
            xrguarddepos = self.nguarddepos[0]
        else:
            xrguard = 0
            xrguarddepos = 0

        if self.bounds[2]==openbc:
            ylguard = self.nguard[1]
            ylguarddepos = self.nguarddepos[1]
        else:
            ylguard = 0
            ylguarddepos = 0
        if self.bounds[3]==openbc:
            yrguard = self.nguard[1]
            yrguarddepos = self.nguarddepos[1]
        else:
            yrguard = 0
            yrguarddepos = 0

        if self.bounds[4]==openbc:
            zlguard = self.nguard[2]
            zlguarddepos = self.nguarddepos[2]
        else:
            zlguard = 0
            zlguarddepos = 0
        if self.bounds[5]==openbc:
            zrguard = self.nguard[2]
            zrguarddepos = self.nguarddepos[2]
        else:
            zrguard = 0
            zrguarddepos = 0

        cd = self.root.listofblocks[self.parents[0]].childdomains
        cd[xl:xu,yl:yu,zl:zu]=0
        cd[xl+xlguarddepos:xu-xrguarddepos,
           yl+ylguarddepos:yu-yrguarddepos,
           zl+zlguarddepos:zu-zrguarddepos]=-self.blocknumber
        cd[xl+xlguard:xu-xrguard,
           yl+ylguard:yu-yrguard,
           zl+zlguard:zu-zrguard]=self.blocknumber

    def initstaticfields(self, relat_species=None,
                         relat_pgroup=None, relat_jslist=None ):
        """
        Initialize the space charge fields, using a Poisson solver.

        - If no arguments are passed, this uses an electrostatic and
        magnetostatic solver, using *all* particles as the source

        - If *either* relat_species *or* relat_pgroup and relat_jslist
        are passed, this uses a relativistic Poisson solver, using *only*
        the designated relativistic particles as the source

        Parameter
        ---------
        relat_species: a Species object
            A relativistic species, whose charge density will serve
            as the source for the relativistic Poisson solver

        relat_pgroup: a ParticleGroup object
            An object containing several species, among which
            some species (which are designated by js_list, below)
            will serve as the source for the relativistic Poission solver.
            This assumes that all the particles are propagating along the z
            direction, essentially at the same speed.

        jslist: a list of integers
            Indicates which of species of the relat_pgroup will serve
            as the source.
        """
        # --- This is needed because of the order in which things
        # --- are imported in warp.py.
        # --- There, MagnetostaticMG is imported after em3dsolver.
        from magnetostaticMG import MagnetostaticMG

        # --- This is needed to create the arrays.
        # --- This should be done to make sure that the EM arrays are set up
        self.allocatedataarrays()

        # If needed, extract relat_pgroup and relat_jslist from the species
        if relat_species is not None:
            if (relat_pgroup is not None) or (relat_jslist is not None):
                raise ValueError('Both `species` and `relat_pgroup`'
                '/`relat_jslist` were passed. Please use one or the other.')
            relat_pgroup = relat_species.pgroup
            relat_jslist = relat_species.jslist

        # --- Calculate the static fields.
        top.grid_overlap = 2
        # Pick either a multigrid, or an openbc solver
        if self.solvergeom == w3d.XYZgeom:
            if w3d.boundxy == openbc:
                try:
                    # Open BC solver - Note: this requires the installation
                    # of the openbc_poisson package (separate repository)
                    from warp.field_solvers.openbcsolver import OpenBC3D as ESolver
                except:
                    ESolver = MultiGrid3D
            else:
                ESolver = MultiGrid3D
        else:
            ESolver = MultiGrid2D
        # Initialize electrostatic solver
        esolver = ESolver(nxguardphi=self.nxguard+1,
                          nyguardphi=self.nyguard+1,
                          nzguardphi=self.nzguard+1,
                          nxguarde=self.nxguard,
                          nyguarde=self.nyguard,
                          nzguarde=self.nzguard,
#                          mgtol = 1.e-20,
                          )
        esolver.conductordatalist = self.conductordatalist
        # check if used for calculation of relativistic beam
        if relat_pgroup is not None:
          # pass particle group to density deposition
          esolver.loadrho( pgroups=[relat_pgroup], jslist=relat_jslist )
          # Calculate the relativistic contraction factor along z
          gaminv = getgaminv( pgroup=relat_pgroup, jslist=relat_jslist )
          zfact = numpy.mean(1./gaminv)
          # Call the relativisitic Poisson solver
          esolver.solve(iwhich=0,zfact=zfact)
        else:
          esolver.loadrho()
          esolver.solve()

        # For non-relativistic bunch, calculate the vector potential
        # by using a magnetostatic multigrid solver
        if relat_pgroup is None:
            bsolver = MagnetostaticMG(luse2D=not(self.solvergeom == w3d.XYZgeom),
                                  nxguardphi=self.nxguard+1,
                                  nyguardphi=self.nyguard+1,
                                  nzguardphi=self.nzguard+1,
                                  nxguarde=self.nxguard,
                                  nyguarde=self.nyguard,
                                  nzguarde=self.nzguard,
                                  mgtol=1.e-20,
                                  )
            bsolver.conductordatalist = self.conductordatalist
            bsolver.loadj()
            bsolver.solve()
            if self.solvergeom == w3d.RZgeom:
                bsolver.field[0,0,...] = 0.
                bsolver.field[1,0,...] = 0.
            # Get the vector potential
            Ax = bsolver.potential[0,...]
            Ay = bsolver.potential[1,...]
            Az = bsolver.potential[2,...]
            # Get the gridstep on which A is calculated
            # (In order to calculate B via finite difference)
            zfact = 1.
            Asolver_dx = bsolver.dx
            Asolver_dy = bsolver.dy
            Asolver_dz = bsolver.dz
        # For a relativistic bunch, calculate the vector potential
        # from phi according to J.-L. Vay, Phys. Plasmas 15, 056701 (2008)
        else:
            gaminv = getgaminv( pgroup=relat_pgroup, jslist=relat_jslist )
            zfact = numpy.mean(1./gaminv)
            beta = sqrt(1.-1./zfact/zfact)
            # Get the vector potential
            Ax = zeros_like(esolver.phi)
            Ay = zeros_like(esolver.phi)
            Az = esolver.phi*beta/clight
            # Get the gridstep on which A is calculated
            # (In order to calculate B via finite difference)
            Asolver_dx = esolver.dx
            Asolver_dy = esolver.dy
            Asolver_dz = esolver.dz

        # --- Reset to the value needed by the EM solver
        top.grid_overlap = 1

        # --- This can be done since the number of guard cells is the same
        # --- for the ES and EM solvers.
        # --- These fields are all node centered.
        # --- This assumes that the domains of the EM solver are all contained within the
        # --- domains of the static solvers. This should be the case since the static solvers
        # --- have a larger grid_overlap.
        ix1 = self.fsdecomp.ix[self.fsdecomp.ixproc] - esolver.fsdecomp.ix[esolver.fsdecomp.ixproc] + 1
        ix2 = ix1 + self.fsdecomp.nx[self.fsdecomp.ixproc] + 1 + 2*self.nxguard
        iy1 = self.fsdecomp.iy[self.fsdecomp.iyproc] - esolver.fsdecomp.iy[esolver.fsdecomp.iyproc] + 1
        iy2 = iy1 + self.fsdecomp.ny[self.fsdecomp.iyproc] + 1 + 2*self.nyguard
        iz1 = self.fsdecomp.iz[self.fsdecomp.izproc] - esolver.fsdecomp.iz[esolver.fsdecomp.izproc] + 1
        iz2 = iz1 + self.fsdecomp.nz[self.fsdecomp.izproc] + 1 + 2*self.nzguard

        # --- Calculate the fields on the Yee mesh by direct finite differences
        # --- of the potential (which is on a node centered grid)
        if self.solvergeom == w3d.XYZgeom:
            Ex = (esolver.phi[ix1:ix2,iy1:iy2,iz1:iz2] - esolver.phi[ix1+1:ix2+1,iy1:iy2,iz1:iz2])/esolver.dx
            Ey = (esolver.phi[ix1:ix2,iy1:iy2,iz1:iz2] - esolver.phi[ix1:ix2,iy1+1:iy2+1,iz1:iz2])/esolver.dy
            Ez = (esolver.phi[ix1:ix2,iy1:iy2,iz1:iz2] - esolver.phi[ix1:ix2,iy1:iy2,iz1+1:iz2+1])/esolver.dz/zfact/zfact

            Bx = 0.5*(-((Ay[ix1:ix2,iy1  :iy2  ,iz1+1:iz2+1] - Ay[ix1:ix2,iy1  :iy2  ,iz1:iz2]) +
                        (Ay[ix1:ix2,iy1+1:iy2+1,iz1+1:iz2+1] - Ay[ix1:ix2,iy1+1:iy2+1,iz1:iz2]))/Asolver_dz
                      +((Az[ix1:ix2,iy1+1:iy2+1,iz1  :iz2  ] - Az[ix1:ix2,iy1:iy2,iz1  :iz2  ]) +
                        (Az[ix1:ix2,iy1+1:iy2+1,iz1+1:iz2+1] - Az[ix1:ix2,iy1:iy2,iz1+1:iz2+1]))/Asolver_dy)
            By = 0.5*(-((Az[ix1+1:ix2+1,iy1:iy2,iz1  :iz2  ] - Az[ix1:ix2,iy1:iy2,iz1  :iz2  ]) +
                        (Az[ix1+1:ix2+1,iy1:iy2,iz1+1:iz2+1] - Az[ix1:ix2,iy1:iy2,iz1+1:iz2+1]))/Asolver_dx
                      +((Ax[ix1  :ix2  ,iy1:iy2,iz1+1:iz2+1] - Ax[ix1  :ix2  ,iy1:iy2,iz1:iz2]) +
                        (Ax[ix1+1:ix2+1,iy1:iy2,iz1+1:iz2+1] - Ax[ix1+1:ix2+1,iy1:iy2,iz1:iz2]))/Asolver_dz)
            Bz = 0.5*(-((Ax[ix1  :ix2  ,iy1+1:iy2+1,iz1:iz2] - Ax[ix1  :ix2  ,iy1:iy2,iz1:iz2]) +
                        (Ax[ix1+1:ix2+1,iy1+1:iy2+1,iz1:iz2] - Ax[ix1+1:ix2+1,iy1:iy2,iz1:iz2]))/Asolver_dy
                      +((Ay[ix1+1:ix2+1,iy1  :iy2  ,iz1:iz2] - Ay[ix1:ix2,iy1  :iy2  ,iz1:iz2]) +
                        (Ay[ix1+1:ix2+1,iy1+1:iy2+1,iz1:iz2] - Ay[ix1:ix2,iy1+1:iy2+1,iz1:iz2]))/Asolver_dx)

        elif self.l_2drz:
            Ex = (esolver.phi[ix1:ix2,:,iz1:iz2] - esolver.phi[ix1+1:ix2+1,:,iz1:iz2])/esolver.dx
            Ey = zeros_like(Ex)
            Ez = (esolver.phi[ix1:ix2,:,iz1:iz2] - esolver.phi[ix1:ix2,:,iz1+1:iz2+1])/esolver.dz/zfact/zfact

            Bx =      -((Ay[ix1:ix2,:,iz1+1:iz2+1] - Ay[ix1:ix2,:,iz1:iz2]))/Asolver_dz
            By = 0.5*(-((Az[ix1+1:ix2+1,:,iz1  :iz2  ] - Az[ix1:ix2,:,iz1  :iz2  ]) +
                        (Az[ix1+1:ix2+1,:,iz1+1:iz2+1] - Az[ix1:ix2,:,iz1+1:iz2+1]))/Asolver_dx
                      +((Ax[ix1  :ix2  ,:,iz1+1:iz2+1] - Ax[ix1  :ix2  ,:,iz1:iz2]) +
                        (Ax[ix1+1:ix2+1,:,iz1+1:iz2+1] - Ax[ix1+1:ix2+1,:,iz1:iz2]))/Asolver_dz)
            ii = arange(ix1,ix2) - esolver.nxguardphi
            rr = (ii + 0.5)*Asolver_dx
            Bz = +(((ii[:,newaxis,newaxis] + 1)*Ay[ix1+1:ix2+1,:,iz1:iz2] - ii[:,newaxis,newaxis]*Ay[ix1:ix2,:,iz1:iz2]))/rr[:,newaxis,newaxis]

        elif self.l_2dxz:
            Ex = (esolver.phi[ix1:ix2,:,iz1:iz2] - esolver.phi[ix1+1:ix2+1,:,iz1:iz2])/esolver.dx
            Ey = zeros_like(Ex)
            Ez = (esolver.phi[ix1:ix2,:,iz1:iz2] - esolver.phi[ix1:ix2,:,iz1+1:iz2+1])/esolver.dz/zfact/zfact

            Bx =      -((Ay[ix1:ix2,:,iz1+1:iz2+1] - Ay[ix1:ix2,:,iz1:iz2]))/Asolver_dz
            By = 0.5*(-((Az[ix1+1:ix2+1,:,iz1  :iz2  ] - Az[ix1:ix2,:,iz1  :iz2  ]) +
                        (Az[ix1+1:ix2+1,:,iz1+1:iz2+1] - Az[ix1:ix2,:,iz1+1:iz2+1]))/Asolver_dx
                      +((Ax[ix1  :ix2  ,:,iz1+1:iz2+1] - Ax[ix1  :ix2  ,:,iz1:iz2]) +
                        (Ax[ix1+1:ix2+1,:,iz1+1:iz2+1] - Ax[ix1+1:ix2+1,:,iz1:iz2]))/Asolver_dz)
            Bz =      +((Ay[ix1+1:ix2+1,:,iz1:iz2] - Ay[ix1:ix2,:,iz1:iz2]))/Asolver_dx

        self.fields.Bx[...] = Bx
        self.fields.By[...] = By
        self.fields.Bz[...] = Bz
        self.fields.Ex[...] = Ex
        self.fields.Ey[...] = Ey
        self.fields.Ez[...] = Ez

        # --- This is copied from push_b_part_2
        # --- novercycle and icycle need to be defined for setebp.
        if self.ntsub<1.:
            self.novercycle = nint(1./self.ntsub)
            self.icycle = (top.it-1)%self.novercycle
        else:
            self.novercycle = 1
            self.icycle = 0

        # --- This is copied from the end of dosolve.
        self.setebp()
        if not self.l_nodalgrid and top.efetch[0] != 4:self.yee2node3d()
        if self.l_nodalgrid and top.efetch[0] == 4:
            self.fields.l_nodecentered=True
            self.node2yee3d()
        if self.l_smooth_particle_fields and any(self.npass_smooth>0):
            self.smoothfields()

    def set_num_Cherenkov_cor_coefs(self):
        if top.efetch[0]==1:self.gather_method="Momentum"
        if top.efetch[0]>1 and self.l_lower_order_in_v:self.gather_method="Galerkin"
        if top.efetch[0]>1 and not self.l_lower_order_in_v:self.gather_method="Uniform"

        self.num_Cherenkov_cor_coefs = {}
        # *************************** Galerkin gather ***************************
        self.num_Cherenkov_cor_coefs["Galerkin"]=AppendableArray(unitshape=[2,4])
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47536,2.04288,-0.598163,0.0314711],[-2.80862,2.80104,-1.14615,0.154077]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47536,2.04288,-0.598163,0.0314711],[-2.80862,2.80104,-1.14615,0.154077]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47545,2.04309,-0.598307,0.0315029],[-2.80851,2.80078,-1.14595,0.154027]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.4756,2.04342,-0.598549,0.0315558],[-2.80832,2.80034,-1.14561,0.153945]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47581,2.0439,-0.598886,0.0316298],[-2.80807,2.79973,-1.14514,0.153829]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47608,2.0445,-0.59932,0.031725],[-2.80774,2.79894,-1.14454,0.15368]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47641,2.04525,-0.59985,0.0318412],[-2.80733,2.79798,-1.1438,0.153498]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.4768,2.04612,-0.600477,0.0319785],[-2.80685,2.79685,-1.14292,0.153284]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47725,2.04714,-0.6012,0.0321367],[-2.8063,2.79554,-1.14192,0.153036]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47776,2.04829,-0.602019,0.0323158],[-2.80568,2.79405,-1.14077,0.152756]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47833,2.04957,-0.602934,0.0325158],[-2.80498,2.79239,-1.1395,0.152443]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47896,2.05099,-0.603944,0.0327364],[-2.80421,2.79056,-1.13809,0.152098]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.47965,2.05254,-0.605051,0.0329777],[-2.80337,2.78856,-1.13656,0.151721]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.4804,2.05423,-0.606253,0.0332396],[-2.80246,2.78638,-1.13488,0.151312]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48121,2.05606,-0.60755,0.0335218],[-2.80147,2.78404,-1.13308,0.150871]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48208,2.05802,-0.608942,0.0338243],[-2.80041,2.78152,-1.13115,0.150397]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48301,2.06012,-0.610429,0.0341469],[-2.79927,2.77882,-1.12908,0.149893]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48401,2.06235,-0.61201,0.0344895],[-2.79807,2.77596,-1.12689,0.149356]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48506,2.06471,-0.613685,0.0348519],[-2.79679,2.77292,-1.12456,0.148789]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48618,2.06721,-0.615453,0.0352339],[-2.79543,2.76972,-1.12211,0.14819]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48735,2.06984,-0.617314,0.0356353],[-2.79401,2.76634,-1.11953,0.14756]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48859,2.07261,-0.619268,0.0360559],[-2.79251,2.76279,-1.11681,0.1469]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.48988,2.0755,-0.621312,0.0364954],[-2.79094,2.75907,-1.11397,0.146208]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.49123,2.07853,-0.623447,0.0369536],[-2.78929,2.75517,-1.111,0.145486]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.49265,2.08169,-0.625672,0.0374302],[-2.78757,2.7511,-1.10789,0.144733]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.49412,2.08498,-0.627986,0.0379248],[-2.78578,2.74686,-1.10466,0.14395]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.49565,2.0884,-0.630386,0.0384372],[-2.78391,2.74245,-1.1013,0.143137]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.49724,2.09194,-0.632873,0.0389669],[-2.78196,2.73786,-1.09781,0.142293]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.49888,2.09561,-0.635443,0.0395135],[-2.77994,2.73309,-1.09419,0.141419]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.50058,2.09939,-0.638096,0.0400766],[-2.77784,2.72814,-1.09043,0.140514]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.50234,2.1033,-0.640829,0.0406557],[-2.77566,2.72301,-1.08654,0.139578]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.50415,2.10732,-0.64364,0.0412502],[-2.7734,2.7177,-1.08252,0.138612]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.50601,2.11145,-0.646526,0.0418594],[-2.77106,2.7122,-1.07836,0.137614]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.50791,2.1157,-0.649485,0.0424828],[-2.76864,2.70651,-1.07406,0.136586]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.50987,2.12004,-0.652512,0.0431196],[-2.76613,2.70062,-1.06962,0.135525]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.51187,2.12448,-0.655604,0.0437688],[-2.76353,2.69453,-1.06503,0.134432]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.51392,2.12901,-0.658756,0.0444297],[-2.76084,2.68824,-1.0603,0.133307]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.516,2.13363,-0.661964,0.0451011],[-2.75806,2.68173,-1.05541,0.132148]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.51812,2.13832,-0.665221,0.0457818],[-2.75518,2.675,-1.05037,0.130954]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.52027,2.14308,-0.668521,0.0464705],[-2.75219,2.66804,-1.04516,0.129725]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.52244,2.14789,-0.671856,0.0471658],[-2.7491,2.66084,-1.03978,0.12846]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.52464,2.15274,-0.675218,0.0478658],[-2.7459,2.65339,-1.03423,0.127156]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.52684,2.15762,-0.678596,0.0485687],[-2.74257,2.64566,-1.02848,0.125813]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.52906,2.16251,-0.68198,0.0492723],[-2.73912,2.63765,-1.02254,0.124428]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.53126,2.16738,-0.685355,0.049974],[-2.73552,2.62934,-1.01638,0.122999]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.53345,2.17222,-0.688706,0.0506708],[-2.73178,2.62069,-1.01,0.121523]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.53561,2.177,-0.692015,0.0513594],[-2.72787,2.61169,-1.00337,0.119996]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.53773,2.18168,-0.69526,0.0520359],[-2.72379,2.6023,-0.996479,0.118417]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.53978,2.18623,-0.698416,0.0526955],[-2.71951,2.59248,-0.989294,0.116778]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54175,2.19059,-0.701452,0.053333],[-2.71501,2.58218,-0.981786,0.115076]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.5436,2.19471,-0.704331,0.0539417],[-2.71026,2.57135,-0.97392,0.113303]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54531,2.19852,-0.70701,0.0545141],[-2.70524,2.55991,-0.965651,0.111453]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54683,2.20193,-0.709433,0.0550409],[-2.69989,2.54778,-0.956922,0.109514]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.5481,2.20483,-0.711533,0.0555106],[-2.69416,2.53484,-0.947666,0.107476]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54906,2.20709,-0.713224,0.0559094],[-2.68799,2.52096,-0.937795,0.105324]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54963,2.20852,-0.714397,0.0562198],[-2.68129,2.50596,-0.927197,0.103039]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54968,2.20888,-0.714907,0.0564196],[-2.67394,2.48959,-0.915724,0.100597]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54905,2.20785,-0.714562,0.0564797],[-2.66578,2.47153,-0.903179,0.097968]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54751,2.20496,-0.713094,0.0563618],[-2.65657,2.4513,-0.889283,0.0951084]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54472,2.19955,-0.710118,0.0560124],[-2.64598,2.42824,-0.873638,0.0919592]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.54014,2.19058,-0.705048,0.0553544],[-2.63347,2.40127,-0.855632,0.0884325]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.53286,2.1763,-0.69693,0.0542684],[-2.61813,2.36864,-0.834261,0.0843898]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.52115,2.15344,-0.684027,0.05255],[-2.59821,2.32701,-0.807691,0.0795876]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.50098,2.11466,-0.66255,0.0497817],[-2.56971,2.26887,-0.77188,0.0735132]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.45797,2.03459,-0.620099,0.0446889],[-2.51823,2.16823,-0.713448,0.0645399]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.28371,1.72254,-0.465905,0.0283268],[-2.33537,1.8294,-0.533852,0.0409941]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.4885,2.04899,-0.599292,0.0390466],[-2.53143,2.14818,-0.670502,0.053982]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.1433,1.36735,-0.220924,-0.00215633],[-2.17737,1.43641,-0.259095,0.00101255]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.4943,2.07019,-0.610552,0.035166],[-2.51929,2.12931,-0.654743,0.0452381]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.84529,2.77303,-1.00018,0.0724884],[-2.86122,2.82221,-1.05039,0.0894636]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.72242,2.51888,-0.847226,0.0509964],[-2.72908,2.54506,-0.87834,0.0626188]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.65633,2.3744,-0.750392,0.0326366],[-2.6536,2.37954,-0.7665,0.0409117]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.59601,2.23412,-0.646421,0.00868027],[-2.58374,2.21923,-0.649738,0.0146791]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.51477,2.0369,-0.491066,-0.0306397],[-2.49284,2.00346,-0.48457,-0.0255348]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.35935,1.65155,-0.178971,-0.112713],[-2.32762,1.60337,-0.1698,-0.105287]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-1.84315,0.361693,0.876104,-0.393844],[-1.80149,0.316787,0.855414,-0.369652]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.65422,2.39262,-0.789663,0.0516265],[-2.60242,2.28418,-0.721378,0.040091]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-3.46529,4.42354,-2.45543,0.497097],[-3.40335,4.25157,-2.29817,0.449834]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-3.15747,3.65311,-1.824,0.328432],[-3.0852,3.47341,-1.67791,0.28982]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-3.04694,3.37613,-1.59668,0.267631],[-2.9642,3.17856,-1.44399,0.229852]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.99205,3.23814,-1.48302,0.237103],[-2.89872,3.01966,-1.31861,0.197945]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.96075,3.15894,-1.41733,0.219317],[-2.85668,2.91811,-1.23894,0.17783]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.94172,3.11028,-1.37649,0.20811],[-2.82679,2.84621,-1.18287,0.163785]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.92994,3.07962,-1.35025,0.200755],[-2.80401,2.79167,-1.14058,0.153278]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.92283,3.06054,-1.33338,0.195859],[-2.78577,2.74819,-1.10706,0.145015]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.91894,3.04938,-1.3229,0.192637],[-2.77061,2.7122,-1.07946,0.138267]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.91736,3.04394,-1.31702,0.190612],[-2.75764,2.68152,-1.05606,0.132589]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.91753,3.04278,-1.31456,0.189477],[-2.74627,2.65475,-1.03575,0.127695]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.91905,3.04494,-1.31475,0.189026],[-2.73612,2.63093,-1.01777,0.123395]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.92165,3.04973,-1.31705,0.189117],[-2.72692,2.6094,-1.00159,0.119553]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.92512,3.05667,-1.32105,0.189646],[-2.71846,2.58968,-0.986841,0.116074]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.92933,3.06539,-1.32646,0.190538],[-2.71061,2.57142,-0.973239,0.112887]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.93416,3.07562,-1.33308,0.191735],[-2.70323,2.55434,-0.960573,0.109937]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.93952,3.08715,-1.34072,0.193194],[-2.69626,2.53824,-0.948678,0.107185]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.94535,3.09982,-1.34925,0.194881],[-2.68962,2.52294,-0.937429,0.104598]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.95159,3.11349,-1.35858,0.196769],[-2.68327,2.50833,-0.926722,0.102151]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.9582,3.12805,-1.36861,0.198838],[-2.67714,2.4943,-0.916477,0.0998223]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.96514,3.14342,-1.37929,0.201068],[-2.67122,2.48076,-0.906627,0.0975966]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.97239,3.15953,-1.39055,0.203448],[-2.66546,2.46764,-0.897118,0.0954599]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.97991,3.17632,-1.40234,0.205964],[-2.65985,2.45489,-0.887903,0.0934011]]))
        self.num_Cherenkov_cor_coefs["Galerkin"].append(array([[-2.98769,3.19374,-1.41463,0.208607],[-2.65437,2.44244,-0.878945,0.0914107]]))
        # *************************** Momentum gather ***************************
        self.num_Cherenkov_cor_coefs["Momentum"]=AppendableArray(unitshape=[2,4])
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98767,3.19368,-1.41458,0.208594],[-2.65428,2.44224,-0.878796,0.0913764]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98767,3.19368,-1.41458,0.208594],[-2.65428,2.44224,-0.878796,0.0913764]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.9876,3.19351,-1.41444,0.208555],[-2.65401,2.44163,-0.878347,0.0912737]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98749,3.19323,-1.41421,0.208491],[-2.65357,2.44061,-0.877602,0.0911031]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98734,3.19285,-1.41388,0.208402],[-2.65296,2.43919,-0.876563,0.0908654]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98716,3.19237,-1.41348,0.20829],[-2.65217,2.43738,-0.875236,0.0905616]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98693,3.19179,-1.41299,0.208157],[-2.65121,2.43518,-0.873624,0.0901933]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98668,3.19114,-1.41244,0.208006],[-2.65009,2.43261,-0.871736,0.089762]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98639,3.19041,-1.41183,0.207837],[-2.6488,2.42966,-0.869579,0.0892697]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98608,3.18962,-1.41116,0.207655],[-2.64736,2.42635,-0.86716,0.0887186]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98576,3.18878,-1.41046,0.207463],[-2.64577,2.4227,-0.86449,0.088111]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98542,3.18791,-1.40972,0.207262],[-2.64403,2.41871,-0.861579,0.0874494]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98507,3.18701,-1.40898,0.207058],[-2.64215,2.4144,-0.858436,0.0867364]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98472,3.18612,-1.40822,0.206852],[-2.64013,2.40979,-0.855074,0.0859749]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98437,3.18523,-1.40748,0.20665],[-2.63798,2.40488,-0.851503,0.0851676]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98403,3.18437,-1.40676,0.206454],[-2.63571,2.39969,-0.847734,0.0843175]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98371,3.18355,-1.40608,0.206269],[-2.63333,2.39425,-0.84378,0.0834274]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98341,3.18279,-1.40544,0.206097],[-2.63083,2.38855,-0.839652,0.0825004]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98314,3.1821,-1.40487,0.205943],[-2.62822,2.38262,-0.835361,0.0815392]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98291,3.1815,-1.40437,0.20581],[-2.62552,2.37647,-0.830918,0.0805468]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98271,3.181,-1.40397,0.205702],[-2.62272,2.37011,-0.826336,0.0795259]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98256,3.18062,-1.40366,0.205622],[-2.61984,2.36357,-0.821624,0.0784793]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98246,3.18037,-1.40346,0.205572],[-2.61687,2.35684,-0.816793,0.0774095]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98241,3.18027,-1.40339,0.205557],[-2.61383,2.34995,-0.811854,0.076319]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98243,3.18033,-1.40345,0.20558],[-2.61071,2.34291,-0.806815,0.0752104]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98251,3.18056,-1.40366,0.205642],[-2.60753,2.33573,-0.801685,0.0740857]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98267,3.18097,-1.40402,0.205747],[-2.60428,2.32842,-0.796474,0.0729473]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.9829,3.18157,-1.40455,0.205896],[-2.60098,2.32098,-0.791189,0.0717971]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98321,3.18238,-1.40524,0.206093],[-2.59762,2.31344,-0.785838,0.070637]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.9836,3.1834,-1.40612,0.206339],[-2.59421,2.3058,-0.780428,0.0694688]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98407,3.18465,-1.40718,0.206637],[-2.59076,2.29806,-0.774966,0.0682941]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98464,3.18612,-1.40844,0.206987],[-2.58726,2.29024,-0.769458,0.0671146]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.9853,3.18783,-1.40989,0.207392],[-2.58372,2.28234,-0.76391,0.0659314]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98606,3.18979,-1.41156,0.207853],[-2.58014,2.27437,-0.758326,0.064746]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98691,3.19199,-1.41343,0.208372],[-2.57653,2.26633,-0.752711,0.0635596]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98787,3.19446,-1.41551,0.208948],[-2.57288,2.25824,-0.74707,0.0623731]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.98892,3.19718,-1.41782,0.209585],[-2.5692,2.25009,-0.741407,0.0611876]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.99008,3.20017,-1.42034,0.210281],[-2.5655,2.24189,-0.735726,0.060004]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.99135,3.20343,-1.42309,0.211039],[-2.56177,2.23365,-0.730029,0.058823]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.99273,3.20697,-1.42608,0.211859],[-2.55801,2.22537,-0.724319,0.0576454]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.99421,3.21078,-1.42929,0.212741],[-2.55422,2.21704,-0.7186,0.0564718]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.9958,3.21487,-1.43273,0.213686],[-2.55042,2.20869,-0.712873,0.0553027]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.99751,3.21925,-1.43641,0.214695],[-2.54659,2.2003,-0.707141,0.0541387]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-2.99933,3.22391,-1.44033,0.215767],[-2.54274,2.19188,-0.701405,0.0529802]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.00126,3.22886,-1.44449,0.216904],[-2.53887,2.18343,-0.695668,0.0518276]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.0033,3.2341,-1.44888,0.218105],[-2.53498,2.17496,-0.68993,0.0506813]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.00546,3.23964,-1.45352,0.219371],[-2.53108,2.16647,-0.684193,0.0495415]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.00774,3.24546,-1.4584,0.220701],[-2.52715,2.15795,-0.678458,0.0484085]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.01013,3.25158,-1.46352,0.222096],[-2.52321,2.14941,-0.672726,0.0472826]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.01264,3.25799,-1.46888,0.223557],[-2.51925,2.14085,-0.666998,0.0461638]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.01526,3.26469,-1.47449,0.225082],[-2.51528,2.13227,-0.661275,0.0450525]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.018,3.2717,-1.48034,0.226672],[-2.51129,2.12367,-0.655557,0.0439488]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.02086,3.279,-1.48644,0.228327],[-2.50728,2.11506,-0.649845,0.0428526]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.02384,3.28659,-1.49278,0.230048],[-2.50326,2.10642,-0.644139,0.0417643]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.02694,3.29449,-1.49936,0.231833],[-2.49922,2.09778,-0.63844,0.0406838]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.03015,3.30268,-1.50619,0.233684],[-2.49516,2.08911,-0.632749,0.0396112]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.03348,3.31117,-1.51327,0.2356],[-2.4911,2.08043,-0.627065,0.0385466]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.03693,3.31996,-1.52059,0.237581],[-2.48701,2.07174,-0.621388,0.03749]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.0405,3.32904,-1.52815,0.239627],[-2.48291,2.06303,-0.61572,0.0364415]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.04418,3.33843,-1.53596,0.241738],[-2.4788,2.05431,-0.610059,0.0354011]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.04799,3.34812,-1.54402,0.243915],[-2.47467,2.04557,-0.604408,0.0343687]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.05191,3.3581,-1.55232,0.246156],[-2.47053,2.03682,-0.598764,0.0333446]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.05596,3.36838,-1.56087,0.248463],[-2.46637,2.02805,-0.593129,0.0323286]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.06012,3.37896,-1.56966,0.250835],[-2.4622,2.01927,-0.587503,0.0313207]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.0644,3.38984,-1.5787,0.253272],[-2.45802,2.01047,-0.581886,0.0303211]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.0688,3.40102,-1.58798,0.255774],[-2.45381,2.00167,-0.576277,0.0293296]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.07332,3.4125,-1.59751,0.258342],[-2.4496,1.99285,-0.570678,0.0283464]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.07795,3.42428,-1.60729,0.260976],[-2.44537,1.98401,-0.565088,0.0273714]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.08271,3.43636,-1.61731,0.263675],[-2.44112,1.97516,-0.559506,0.0264046]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.08758,3.44874,-1.62758,0.26644],[-2.43686,1.9663,-0.553934,0.0254461]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.09257,3.46142,-1.6381,0.269271],[-2.43259,1.95742,-0.548372,0.0244959]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.09769,3.47439,-1.64886,0.272168],[-2.4283,1.94853,-0.542818,0.0235539]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.10292,3.48767,-1.65987,0.275131],[-2.424,1.93963,-0.537274,0.0226202]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.10826,3.50125,-1.67113,0.278161],[-2.41968,1.93072,-0.53174,0.0216949]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.11373,3.51513,-1.68264,0.281258],[-2.41534,1.92179,-0.526216,0.0207778]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.11932,3.52931,-1.6944,0.284421],[-2.411,1.91284,-0.520701,0.0198692]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.12502,3.54379,-1.7064,0.287652],[-2.40663,1.90389,-0.515196,0.0189689]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.13084,3.55857,-1.71866,0.29095],[-2.40225,1.89492,-0.5097,0.0180771]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.13678,3.57365,-1.73117,0.294316],[-2.39786,1.88594,-0.504215,0.0171938]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.14284,3.58903,-1.74393,0.297751],[-2.39345,1.87694,-0.49874,0.0163189]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.14902,3.60472,-1.75694,0.301254],[-2.38903,1.86793,-0.493276,0.0154526]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.15531,3.6207,-1.7702,0.304825],[-2.38459,1.85891,-0.487822,0.0145948]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.16173,3.63699,-1.78372,0.308466],[-2.38014,1.84988,-0.482378,0.0137456]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.16826,3.65359,-1.79749,0.312177],[-2.37567,1.84083,-0.476946,0.0129052]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.17491,3.67048,-1.81152,0.315958],[-2.37119,1.83177,-0.471524,0.0120734]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.18168,3.68768,-1.82581,0.31981],[-2.36669,1.8227,-0.466113,0.0112504]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.18856,3.70519,-1.84035,0.323732],[-2.36217,1.81361,-0.460713,0.0104362]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.19557,3.72299,-1.85515,0.327726],[-2.35764,1.80451,-0.455325,0.00963088]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.20269,3.74111,-1.8702,0.331792],[-2.3531,1.7954,-0.449949,0.0088345]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.20993,3.75953,-1.88552,0.33593],[-2.34854,1.78628,-0.444584,0.00804712]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.21729,3.77825,-1.9011,0.340142],[-2.34396,1.77714,-0.439231,0.00726881]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.22476,3.79729,-1.91694,0.344427],[-2.33937,1.768,-0.43389,0.00649963]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.23236,3.81662,-1.93305,0.348786],[-2.33476,1.75884,-0.428562,0.00573964]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.24007,3.83627,-1.94942,0.353221],[-2.33014,1.74966,-0.423246,0.00498893]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.2479,3.85623,-1.96606,0.35773],[-2.3255,1.74048,-0.417943,0.00424756]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.25584,3.87649,-1.98296,0.362316],[-2.32085,1.73128,-0.412654,0.00351561]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.26391,3.89706,-2.00013,0.366978],[-2.31618,1.72208,-0.407377,0.00279317]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.27209,3.91795,-2.01758,0.371718],[-2.31149,1.71286,-0.402114,0.0020803]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.28039,3.93914,-2.03529,0.376536],[-2.30679,1.70363,-0.396865,0.00137709]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.2888,3.96065,-2.05328,0.381432],[-2.30207,1.69438,-0.391629,0.000683629]]))
        self.num_Cherenkov_cor_coefs["Momentum"].append(array([[-3.2888,3.96065,-2.05328,0.381432],[-2.30207,1.69438,-0.391629,0.000683629]]))
        # *************************** Uniform gather ***************************
        self.num_Cherenkov_cor_coefs["Uniform"]=AppendableArray(unitshape=[2,4])
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54365,2.19481,-0.704398,0.0539562],[-2.71024,2.57129,-0.973879,0.113294]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54365,2.19481,-0.704398,0.0539562],[-2.71024,2.57129,-0.973879,0.113294]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54377,2.19509,-0.704601,0.0539996],[-2.71017,2.57112,-0.973754,0.113264]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54399,2.19557,-0.704937,0.054072],[-2.71005,2.57084,-0.973545,0.113215]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54429,2.19624,-0.705409,0.0541734],[-2.70988,2.57045,-0.973254,0.113147]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54467,2.1971,-0.706015,0.0543038],[-2.70966,2.56995,-0.972879,0.113059]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54514,2.19815,-0.706757,0.0544631],[-2.7094,2.56934,-0.972422,0.112952]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54569,2.19939,-0.707634,0.0546515],[-2.70909,2.56861,-0.971882,0.112825]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54634,2.20083,-0.708646,0.054869],[-2.70873,2.56778,-0.97126,0.112679]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54706,2.20245,-0.709793,0.0551155],[-2.70832,2.56683,-0.970556,0.112514]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54788,2.20427,-0.711076,0.0553912],[-2.70787,2.56578,-0.96977,0.112329]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54878,2.20628,-0.712495,0.055696],[-2.70737,2.56462,-0.968902,0.112126]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.54976,2.20849,-0.714051,0.05603],[-2.70682,2.56334,-0.967954,0.111904]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.55083,2.21089,-0.715742,0.0563933],[-2.70623,2.56196,-0.966926,0.111663]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.55199,2.21348,-0.717571,0.0567859],[-2.70559,2.56047,-0.965817,0.111403]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.55324,2.21627,-0.719536,0.0572079],[-2.7049,2.55888,-0.964629,0.111125]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.55457,2.21925,-0.721639,0.0576593],[-2.70417,2.55717,-0.963362,0.110828]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.55599,2.22243,-0.72388,0.0581403],[-2.70339,2.55536,-0.962017,0.110514]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.5575,2.2258,-0.726259,0.0586508],[-2.70256,2.55345,-0.960594,0.110181]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.5591,2.22938,-0.728777,0.0591909],[-2.70169,2.55143,-0.959095,0.109831]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.56078,2.23315,-0.731433,0.0597608],[-2.70078,2.54931,-0.957519,0.109463]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.56256,2.23711,-0.73423,0.0603606],[-2.69982,2.54709,-0.955867,0.109077]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.56442,2.24128,-0.737166,0.0609903],[-2.69882,2.54476,-0.954141,0.108675]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.56638,2.24565,-0.740243,0.0616501],[-2.69777,2.54234,-0.95234,0.108255]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.56842,2.25022,-0.743461,0.06234],[-2.69668,2.53981,-0.950466,0.107819]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.57055,2.25499,-0.74682,0.0630602],[-2.69555,2.53718,-0.94852,0.107366]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.57278,2.25996,-0.750322,0.0638108],[-2.69437,2.53446,-0.946501,0.106897]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.57509,2.26513,-0.753967,0.0645919],[-2.69315,2.53164,-0.944412,0.106411]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.5775,2.27051,-0.757755,0.0654037],[-2.69189,2.52872,-0.942253,0.10591]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.57999,2.27609,-0.761687,0.0662464],[-2.69059,2.52571,-0.940025,0.105393]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.58258,2.28188,-0.765764,0.06712],[-2.68925,2.52261,-0.937729,0.104861]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.58527,2.28788,-0.769987,0.0680247],[-2.68786,2.51941,-0.935366,0.104314]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.58804,2.29408,-0.774356,0.0689608],[-2.68644,2.51612,-0.932936,0.103752]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.59091,2.30049,-0.778872,0.0699283],[-2.68497,2.51274,-0.930441,0.103175]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.59387,2.30712,-0.783535,0.0709276],[-2.68347,2.50926,-0.927881,0.102584]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.59693,2.31395,-0.788347,0.0719586],[-2.68193,2.5057,-0.925257,0.101979]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.60008,2.32099,-0.793308,0.0730218],[-2.68034,2.50206,-0.922572,0.101361]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.60333,2.32825,-0.79842,0.0741173],[-2.67872,2.49832,-0.919824,0.100729]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.60667,2.33572,-0.803682,0.0752453],[-2.67707,2.49451,-0.917016,0.100083]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.61011,2.3434,-0.809096,0.076406],[-2.67537,2.4906,-0.914148,0.0994249]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.61364,2.35131,-0.814663,0.0775997],[-2.67364,2.48662,-0.911222,0.098754]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.61727,2.35942,-0.820384,0.0788266],[-2.67187,2.48255,-0.908238,0.0980708]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.621,2.36776,-0.826259,0.080087],[-2.67006,2.4784,-0.905197,0.0973754]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.62482,2.37632,-0.83229,0.0813812],[-2.66822,2.47417,-0.9021,0.0966682]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.62875,2.38509,-0.838477,0.0827094],[-2.66635,2.46987,-0.898949,0.0959495]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.63277,2.39409,-0.844822,0.084072],[-2.66443,2.46548,-0.895744,0.0952195]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.63689,2.40331,-0.851325,0.0854692],[-2.66249,2.46102,-0.892486,0.0944786]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.64111,2.41275,-0.857988,0.0869014],[-2.66051,2.45649,-0.889176,0.093727]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.64543,2.42242,-0.864812,0.0883688],[-2.65849,2.45188,-0.885815,0.0929649]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.64985,2.43232,-0.871797,0.0898718],[-2.65645,2.44719,-0.882405,0.0921928]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.65442,2.44255,-0.879027,0.0914289],[-2.65435,2.4424,-0.878921,0.0914059]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.65899,2.45279,-0.886257,0.0929859],[-2.65226,2.43761,-0.875437,0.0906191]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.66371,2.46337,-0.893734,0.0945978],[-2.65011,2.43272,-0.871883,0.0898181]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.66853,2.47418,-0.901378,0.0962468],[-2.64794,2.42776,-0.868282,0.089008]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.67346,2.48522,-0.909189,0.0979331],[-2.64573,2.42272,-0.864635,0.0881892]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.67849,2.4965,-0.917168,0.0996572],[-2.64349,2.41763,-0.860945,0.0873619]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.68362,2.50801,-0.925317,0.10142],[-2.64122,2.41246,-0.85721,0.0865262]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.68885,2.51975,-0.933638,0.103221],[-2.63892,2.40723,-0.853434,0.0856826]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.69419,2.53173,-0.942131,0.105061],[-2.63659,2.40194,-0.849615,0.0848312]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.69963,2.54395,-0.950797,0.10694],[-2.63423,2.39658,-0.845756,0.0839724]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.70518,2.55641,-0.959639,0.10886],[-2.63185,2.39116,-0.841856,0.0831062]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.71083,2.56911,-0.968657,0.11082],[-2.62943,2.38569,-0.837917,0.0822331]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.71658,2.58205,-0.977852,0.112821],[-2.62698,2.38015,-0.83394,0.0813533]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.72244,2.59523,-0.987227,0.114863],[-2.62451,2.37455,-0.829926,0.0804669]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.7284,2.60866,-0.996782,0.116947],[-2.62201,2.36889,-0.825874,0.0795743]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.73447,2.62233,-1.00652,0.119074],[-2.61948,2.36318,-0.821787,0.0786756]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.74065,2.63625,-1.01644,0.121244],[-2.61692,2.35741,-0.817664,0.0777711]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.74693,2.65041,-1.02654,0.123457],[-2.61434,2.35158,-0.813507,0.0768611]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.75332,2.66482,-1.03684,0.125714],[-2.61173,2.3457,-0.809317,0.0759458]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.75982,2.67949,-1.04731,0.128015],[-2.60909,2.33976,-0.805093,0.0750253]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.76642,2.6944,-1.05798,0.130362],[-2.60643,2.33377,-0.800837,0.0740999]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.77313,2.70957,-1.06884,0.132755],[-2.60374,2.32773,-0.79655,0.0731699]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.77995,2.72499,-1.07989,0.135194],[-2.60102,2.32163,-0.792231,0.0722354]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.78687,2.74066,-1.09114,0.137679],[-2.59828,2.31549,-0.787883,0.0712967]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.79391,2.75659,-1.10258,0.140213],[-2.59551,2.30929,-0.783505,0.0703539]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.80105,2.77278,-1.11422,0.142794],[-2.59272,2.30304,-0.779099,0.0694073]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.8083,2.78923,-1.12605,0.145425],[-2.5899,2.29675,-0.774664,0.0684571]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.81565,2.80593,-1.13809,0.148105],[-2.58706,2.2904,-0.770201,0.0675035]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.82312,2.8229,-1.15033,0.150835],[-2.5842,2.28401,-0.765712,0.0665466]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.8307,2.84013,-1.16278,0.153615],[-2.58131,2.27757,-0.761196,0.0655867]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.83838,2.85762,-1.17543,0.156448],[-2.57839,2.27108,-0.756655,0.064624]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.84618,2.87537,-1.18828,0.159332],[-2.57546,2.26454,-0.752089,0.0636587]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.85408,2.8934,-1.20135,0.16227],[-2.57249,2.25796,-0.747498,0.0626909]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.8621,2.91168,-1.21463,0.165261],[-2.56951,2.25134,-0.742883,0.0617208]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.87022,2.93024,-1.22812,0.168306],[-2.5665,2.24467,-0.738245,0.0607487]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.87846,2.94907,-1.24182,0.171407],[-2.56347,2.23796,-0.733584,0.0597746]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.8868,2.96817,-1.25574,0.174564],[-2.56041,2.2312,-0.728901,0.0587988]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.89526,2.98754,-1.26988,0.177777],[-2.55733,2.2244,-0.724196,0.0578215]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.90382,3.00718,-1.28424,0.181048],[-2.55423,2.21755,-0.719469,0.0568428]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.9125,3.0271,-1.29882,0.184377],[-2.55111,2.21067,-0.714722,0.0558629]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.92128,3.04729,-1.31363,0.187765],[-2.54796,2.20374,-0.709955,0.054882]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.93018,3.06776,-1.32866,0.191213],[-2.54479,2.19677,-0.705168,0.0539002]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.93919,3.08851,-1.34393,0.194722],[-2.5416,2.18976,-0.700362,0.0529177]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.94831,3.10954,-1.35942,0.198292],[-2.53839,2.18271,-0.695537,0.0519347]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.95754,3.13085,-1.37514,0.201925],[-2.53515,2.17562,-0.690693,0.0509514]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.96688,3.15245,-1.3911,0.205621],[-2.5319,2.1685,-0.685832,0.0499678]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.97633,3.17432,-1.4073,0.209381],[-2.52862,2.16133,-0.680954,0.0489842]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.9859,3.19648,-1.42373,0.213206],[-2.52531,2.15412,-0.676058,0.0480007]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-2.99557,3.21893,-1.44041,0.217097],[-2.52199,2.14687,-0.671146,0.0470175]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-3.00536,3.24167,-1.45733,0.221056],[-2.51864,2.13959,-0.666218,0.0460347]]))
        self.num_Cherenkov_cor_coefs["Uniform"].append(array([[-3.01526,3.26469,-1.47449,0.225082],[-2.51528,2.13227,-0.661275,0.0450525]]))

    def get_num_Cherenkov_cor_coefs(self):
        dtodz=clight*top.dt/self.dz
        allcoefs=self.num_Cherenkov_cor_coefs[self.gather_method]
        ndt=shape(allcoefs[...])[0]
        dtodz_unit=1./(ndt-1)
        dtodz_norm=dtodz/dtodz_unit
        i=max(0,min(ndt-2,int(dtodz_norm)))
        w=dtodz_norm-i
        excoef=(1.-w)*allcoefs[i,0,:]+w*allcoefs[i+1,0,:]
        bycoef=(1.-w)*allcoefs[i,1,:]+w*allcoefs[i+1,1,:]
        return excoef,bycoef

    def simpledump(self,filename='fields',dir='./',l_appendtimestep=True):
        if me==0:
            fname = dir+filename
            if l_appendtimestep:fname+='_%g'%top.it
            fname+='.pck'
            f=PW.PW(fname)
        ex=self.gatherex()
        if me==0:f.ex=ex
        ey=self.gatherey()
        if me==0:f.ey=ey
        ez=self.gatherez()
        if me==0:f.ez=ez
        bx=self.gatherbx()
        if me==0:f.bx=bx
        by=self.gatherby()
        if me==0:f.by=by
        bz=self.gatherbz()
        if me==0:f.bz=bz
        if self.l_pushf:
            ff=self.gatherf()
            if me==0:f.f=ff
        if self.l_pushg:
            g=self.gatherg()
            if me==0:f.g=g
        jx=self.gatherjx()
        if me==0:f.jx=jx
        jy=self.gatherjy()
        if me==0:f.jy=jy
        jz=self.gatherjz()
        if me==0:f.jz=jz
        if self.l_getrho:
            rho=self.gatherrho()
            if me==0:f.rho=rho
            rhoold=self.gatherrhoold()
            if me==0:f.rhoold=rhoold
        if me==0:f.close()
        del ex,ey,ez,bx,by,bz,jx,jy,jz
        if self.l_pushf:del f
        if self.l_pushg:del g
        if self.l_getrho:del rho,rhoold

    def simpleread(self,filename='fields',dir='./',l_appendtimestep=True):
        if me==0:
            fname = dir+filename
            if l_appendtimestep:fname+='_%g'%top.it
            fname+='.pck'
            f=PR.PR(fname)
        toreturn = {}
        toreturn['ex']=f.ex
        toreturn['ey']=f.ey
        toreturn['ez']=f.ez
        toreturn['jx']=f.jx
        toreturn['jy']=f.jy
        toreturn['jz']=f.jz
        toreturn['bx']=f.bx
        toreturn['by']=f.by
        toreturn['bz']=f.bz

        if self.l_pushf:toreturn['f']=f.f
        if self.l_pushg:toreturn['g']=f.g
        if self.l_getrho:toreturn['rho']=f.rho
        if self.l_getrho:toreturn['rhoold']=f.rhoold
        return toreturn

    def simplecompare(self,filename='fields',dir='./',l_appendtimestep=True,win=0):
        fields = self.simpleread(filename,dir,l_appendtimestep)
        myfields = {}
        myfields['ex']=self.gatherex()
        myfields['ey']=self.gatherey()
        myfields['ez']=self.gatherez()
        myfields['bx']=self.gatherbx()
        myfields['by']=self.gatherby()
        myfields['bz']=self.gatherbz()
        myfields['jx']=self.gatherjx()
        myfields['jy']=self.gatherjy()
        myfields['jz']=self.gatherjz()
        if self.l_pushf: myfields['f']=self.gatherf()
        if self.l_pushg: myfields['g']=self.gatherg()
        if self.l_getrho:
            myfields['rho']=self.gatherrho()
            myfields['rhoold']=self.gatherrhoold()
        flist = ['ex','ey','ez','bx','by','bz','jx','jy','jz']
        if self.l_pushf:flist+=['f']
        if self.l_pushg:flist+=['g']
        if self.l_getrho:
            flist+=['rho']
            flist+=['rhoold']
        if me==0:
            window(win)
            for ff in flist:
                print 'plotting '+ff
                ppgeneric(fields[ff],view=3)
                ppgeneric(myfields[ff],view=4)
                ppgeneric(fields[ff]-myfields[ff],view=5)
                plsys(1);ptitles(ff)
                raw_input('press return for next field...')
                fma(0)

def allocatesf(f,stencil):
    f.syf = EM3D_SPLITYEEFIELDtype()
    f.fieldtype = f.syf.fieldtype
    f.proc=me
    f.syf.stencil=stencil

def allocatef(f,stencil):
    f.yf = EM3D_YEEFIELDtype()
    f.fieldtype = f.yf.fieldtype
    f.proc=me
    f.yf.stencil=stencil

def pyinit_3dem_block(nx, ny, nz,
                      nbndx, nbndy, nbndz,
                      nxguard, nyguard, nzguard,
                      norderx,nordery,norderz,
                      dt, dx, dy, dz,
                      clight, mu0,
                      xmin, ymin, zmin,
                      xlb, ylb, zlb, xrb, yrb, zrb,
                      deposit_energy_density,
                      l_nodalgrid,
                      refinement,l_pushf,l_pushg,
                      l_getrho,
                      stencil,
                      npass_smooth,
                      l_smooth_particle_fields,
                      ntsub,
                      l_1dz,
                      l_2dxz,
                      l_2drz,
                      theta_damp,
                      sigmae,
                      sigmab,
                      excoef,
                      bycoef,
                      pml_method,
                      circ_m):

    if xlb == em3d.otherproc:
        procxl = top.procneighbors[0,0]
    else:
        procxl = me
    if xrb == em3d.otherproc:
        procxr = top.procneighbors[1,0]
    else:
        procxr = me
    if ylb == em3d.otherproc:
        procyl = top.procneighbors[0,1]
    else:
        procyl = me
    if yrb == em3d.otherproc:
        procyr = top.procneighbors[1,1]
    else:
        procyr = me
    if zlb == em3d.otherproc:
        proczl = top.procneighbors[0,2]
    else:
        proczl = me
    if zrb == em3d.otherproc:
        proczr = top.procneighbors[1,2]
    else:
        proczr = me

    b = EM3D_BLOCKtype()
    b.core = EM3D_FIELDtype()
    b.core.yf = EM3D_YEEFIELDtype()
    b.core.fieldtype = b.core.yf.fieldtype
    b.core.proc = me
    b.core.yf.stencil=stencil

    b.nx = nx
    b.ny = ny
    b.nz = nz
    b.nbndx = nbndx
    b.nbndy = nbndy
    b.nbndz = nbndz
    b.nxguard = nxguard
    b.nyguard = nyguard
    b.nzguard = nzguard
    b.xmin = xmin
    b.ymin = ymin
    b.zmin = zmin
    b.dx = dx
    b.dy = dy
    b.dz = dz
    b.xmax = xmin+dx*nx
    b.ymax = ymin+dy*ny
    b.zmax = zmin+dz*nz
    b.dxi = 1./dx
    b.dyi = 1./dy
    b.dzi = 1./dz
    b.xlbnd = xlb
    b.xrbnd = xrb
    b.ylbnd = ylb
    b.yrbnd = yrb
    b.zlbnd = zlb
    b.zrbnd = zrb

    f=b.core.yf
    f.l_macroscopic = False
    f.l_nodalgrid = l_nodalgrid
    f.nx = nx
    f.ny = ny
    f.nz = nz
    f.circ_m = circ_m
    if norderx is not inf:f.norderx = norderx
    if norderx is not inf:f.nordery = nordery
    if norderx is not inf:f.norderz = norderz
    f.theta_damp=theta_damp
    f.sigmae=sigmae
    f.sigmab=sigmab
    if l_getrho:
        f.nxr = f.nx
        f.nyr = f.ny
        f.nzr = f.nz
    if 0:#refinement is None and (all(npass_smooth==0) or not l_smooth_particle_fields):
        f.nxp = 0
        f.nyp = 0
        f.nzp = 0
    else:
        f.nxp = f.nx
        f.nyp = f.ny
        f.nzp = f.nz
    if ntsub<0.75:
        f.nxpnext = f.nx
        f.nypnext = f.ny
        f.nzpnext = f.nz
    if not l_pushf:
        f.nxf = 0
        f.nyf = 0
        f.nzf = 0
    else:
        f.nxf = f.nx
        f.nyf = f.ny
        f.nzf = f.nz
    if not l_pushg:
        f.nxg = 0
        f.nyg = 0
        f.nzg = 0
    else:
        f.nxg = f.nx
        f.nyg = f.ny
        f.nzg = f.nz
    if f.theta_damp != 0.:
        f.nxdamp=f.nx
        f.nydamp=f.ny
        f.nzdamp=f.nz
    if f.sigmae != 0. or f.sigmab != 0.:
        f.nxext=f.nx
        f.nyext=f.ny
        f.nzext=f.nz
    f.nxguard = nxguard
    f.nyguard = nyguard
    f.nzguard = nzguard
    # set min/max of cells positions with FORTRAN indexing
    f.ixmin = 0
    f.iymin = 0
    f.izmin = 0
    f.ixmax = f.nx
    f.iymax = f.ny
    f.izmax = f.nz
    f.ixming = -f.nxguard
    f.iyming = -f.nyguard
    f.izming = -f.nzguard
    f.ixmaxg = f.ixmax+f.nxguard
    f.iymaxg = f.iymax+f.nyguard
    f.izmaxg = f.izmax+f.nzguard
    # set min/max of cells positions with Python indexing
    f.jxmin = f.ixmin-f.ixming
    f.jymin = f.iymin-f.iyming
    f.jzmin = f.izmin-f.izming
    f.jxmax = f.ixmax-f.ixming
    f.jymax = f.iymax-f.iyming
    f.jzmax = f.izmax-f.izming
    f.jxming = 0
    f.jyming = 0
    f.jzming = 0
    f.jxmaxg = f.ixmaxg-f.ixming
    f.jymaxg = f.iymaxg-f.iyming
    f.jzmaxg = f.izmaxg-f.izming
    f.xmin = xmin
    f.ymin = ymin
    f.zmin = zmin
    f.dx = dx
    f.dy = dy
    f.dz = dz
    f.xmax = xmin+dx*nx
    f.ymax = ymin+dy*ny
    f.zmax = zmin+dz*nz
    f.dxi = 1./dx
    f.dyi = 1./dy
    f.dzi = 1./dz
    f.clight = clight
    f.mu0    = mu0
    f.l_1dz = l_1dz
    f.l_2dxz = l_2dxz
    f.l_2drz = l_2drz

    if deposit_energy_density:
        f.nxmp=f.nx
        f.nymp=f.ny
        f.nzmp=f.nz

    f.gchange()

    if f.norderx is not inf:f.xcoefs[:] = FD_weights_hvincenti(f.norderx,l_staggered=not l_nodalgrid)
    if f.nordery is not inf:f.ycoefs[:] = FD_weights_hvincenti(f.nordery,l_staggered=not l_nodalgrid)
    if f.norderx is not inf:f.zcoefs[:] = FD_weights_hvincenti(f.norderz,l_staggered=not l_nodalgrid)

    if excoef is not None:
        # --- sets field deposition transformation stencil
        f.ex_stencil[0] =  (256+128*excoef[0]+96*excoef[1]+80*excoef[2]+70*excoef[3])/256
        f.ex_stencil[1] = -(     64*excoef[0]+64*excoef[1]+60*excoef[2]+56*excoef[3])/256
        f.ex_stencil[2] =  (                  16*excoef[1]+24*excoef[2]+28*excoef[3])/256
        f.ex_stencil[3] = -(                                4*excoef[2]+ 8*excoef[3])/256
        f.ex_stencil[4] =  (                                             1*excoef[3])/256
    if bycoef is not None:
        f.by_stencil[0] =  (256+128*bycoef[0]+96*bycoef[1]+80*bycoef[2]+70*bycoef[3])/256
        f.by_stencil[1] = -(     64*bycoef[0]+64*bycoef[1]+60*bycoef[2]+56*bycoef[3])/256
        f.by_stencil[2] =  (                  16*bycoef[1]+24*bycoef[2]+28*bycoef[3])/256
        f.by_stencil[3] = -(                                4*bycoef[2]+ 8*bycoef[3])/256
        f.by_stencil[4] =  (                                             1*bycoef[3])/256

    if 0:#refinement is None and  (all(npass_smooth==0) or not l_smooth_particle_fields):
        f.nxp = f.nx
        f.nyp = f.ny
        f.nzp = f.nz
        f.Exp = f.Ex
        f.Eyp = f.Ey
        f.Ezp = f.Ez
        f.Bxp = f.Bx
        f.Byp = f.By
        f.Bzp = f.Bz

    f.Ex[...] = 0.
    f.Ey[...] = 0.
    f.Ez[...] = 0.
    f.Bx[...] = 0.
    f.By[...] = 0.
    f.Bz[...] = 0.
    f.Exp[...] = 0.
    f.Eyp[...] = 0.
    f.Ezp[...] = 0.
    f.Bxp[...] = 0.
    f.Byp[...] = 0.
    f.Bzp[...] = 0.
    f.Jxarray[...] = 0.
    f.Jyarray[...] = 0.
    f.Jzarray[...] = 0.
    if l_getrho :
        f.Rhoarray[...] = 0.
    if f.circ_m>0:
        f.Ex_circ[...] = 0.
        f.Ey_circ[...] = 0.
        f.Ez_circ[...] = 0.
        f.Bx_circ[...] = 0.
        f.By_circ[...] = 0.
        f.Bz_circ[...] = 0.
        f.Exp_circ[...] = 0.
        f.Eyp_circ[...] = 0.
        f.Ezp_circ[...] = 0.
        f.Bxp_circ[...] = 0.
        f.Byp_circ[...] = 0.
        f.Bzp_circ[...] = 0.
        f.Jx_circ[...] = 0.
        f.Jy_circ[...] = 0.
        f.Jz_circ[...] = 0.
        f.Jarray_circ[...] = 0.
        f.Rho_circ[...] = 0.

    nnx=em3d.nn
    nny=em3d.nn
    nnz=em3d.nn
    smaxx=em3d.s_max_init
    smaxy=em3d.s_max_init
    smaxz=em3d.s_max_init
    sdeltax=em3d.s_delta
    sdeltay=em3d.s_delta
    sdeltaz=em3d.s_delta

    xminm = b.xmin-nbndx*dx
    xmin0 = b.xmin
    xminp = b.xmax
    yminm = b.ymin-nbndy*dy
    ymin0 = b.ymin
    yminp = b.ymax
    zminm = b.zmin-nbndz*dz
    zmin0 = b.zmin
    zminp = b.zmax

# --- sides
# x
    b.sidexl = EM3D_FIELDtype()
    if xlb==openbc:
        allocatesf(b.sidexl,stencil)
        init_splitfield(b.sidexl.syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, ymin0, zmin0,
                        clight,-1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    b.sidexl.proc=me
    if xlb==periodic:
        b.sidexl=b.core
    if xlb==em3d.otherproc:
        allocatef(b.sidexl,stencil)
        b.sidexl.proc=procxl
    b.sidexr = EM3D_FIELDtype()
    if xrb==openbc:
        allocatesf(b.sidexr,stencil)
        init_splitfield(b.sidexr.syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, ymin0, zmin0,
                        clight, 1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    b.sidexr.proc=me
    if xrb==periodic:
        b.sidexr=b.core
    if xrb==em3d.otherproc:
        allocatef(b.sidexr,stencil)
        b.sidexr.proc=procxr
# y
    b.sideyl = EM3D_FIELDtype()
    if ylb==openbc:
        allocatesf(b.sideyl,stencil)
        init_splitfield(b.sideyl.syf,nx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xmin0, yminm, zmin0,
                        clight, 0,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    b.sideyl.proc=me
    if ylb==periodic:
        b.sideyl=b.core
    if ylb==em3d.otherproc:
        allocatef(b.sideyl,stencil)
        b.sideyl.proc=procyl
    b.sideyr = EM3D_FIELDtype()
    if yrb==openbc:
        allocatesf(b.sideyr,stencil)
        init_splitfield(b.sideyr.syf,nx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xmin0, yminp, zmin0,
                        clight, 0, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    b.sideyr.proc=me
    if yrb==periodic:
        b.sideyr=b.core
    if yrb==em3d.otherproc:
        allocatef(b.sideyr,stencil)
        b.sideyr.proc=procyr
# z
    b.sidezl = EM3D_FIELDtype()
    if zlb==openbc:
        allocatesf(b.sidezl,stencil)
        init_splitfield(b.sidezl.syf,nx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xmin0, ymin0, zminm,
                        clight, 0, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    b.sidezl.proc=me
    if zlb==periodic:
        b.sidezl=b.core
    if zlb==em3d.otherproc:
        allocatef(b.sidezl,stencil)
        b.sidezl.proc=proczl
    b.sidezr = EM3D_FIELDtype()
    if zrb==openbc:
        allocatesf(b.sidezr,stencil)
        init_splitfield(b.sidezr.syf,nx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xmin0, ymin0, zminp,
                        clight, 0, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    b.sidezr.proc=me
    if zrb==periodic:
        b.sidezr=b.core
    if zrb==em3d.otherproc:
        allocatef(b.sidezr,stencil)
        b.sidezr.proc=proczr

# --- edges
# xy
    b.edgexlyl = EM3D_FIELDtype()
    if xlb==openbc and ylb==openbc:
        allocatesf(b.edgexlyl,stencil)
        init_splitfield(b.edgexlyl.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, yminm, zmin0,
                        clight,-1,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xlb==openbc:
            if ylb==periodic:
                b.edgexlyl=b.sidexl
            if ylb==em3d.otherproc:
                allocatesf(b.edgexlyl,stencil)
                b.edgexlyl.proc=procyl
        if ylb==openbc:
            if xlb==periodic:
                b.edgexlyl=b.sideyl
            if xlb==em3d.otherproc:
                allocatesf(b.edgexlyl,stencil)
                b.edgexlyl.proc=procxl

    b.edgexryl = EM3D_FIELDtype()
    if xrb==openbc and ylb==openbc:
        allocatesf(b.edgexryl,stencil)
        init_splitfield(b.edgexryl.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, yminm, zmin0,
                        clight, 1,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xrb==openbc:
            if ylb==periodic:
                b.edgexryl=b.sidexr
            if ylb==em3d.otherproc:
                allocatesf(b.edgexryl,stencil)
                b.edgexryl.proc=procyl
        if ylb==openbc:
            if xrb==periodic:
                b.edgexryl=b.sideyl
            if xrb==em3d.otherproc:
                allocatesf(b.edgexryl,stencil)
                b.edgexryl.proc=procxr

    b.edgexlyr = EM3D_FIELDtype()
    if xlb==openbc and yrb==openbc:
        allocatesf(b.edgexlyr,stencil)
        init_splitfield(b.edgexlyr.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, yminp, zmin0,
                        clight,-1, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xlb==openbc:
            if yrb==periodic:
                b.edgexlyr=b.sidexl
            if yrb==em3d.otherproc:
                allocatesf(b.edgexlyr,stencil)
                b.edgexlyr.proc=procyr
        if yrb==openbc:
            if xlb==periodic:
                b.edgexlyr=b.sideyr
            if xlb==em3d.otherproc:
                allocatesf(b.edgexlyr,stencil)
                b.edgexlyr.proc=procxl

    b.edgexryr = EM3D_FIELDtype()
    if xrb==openbc and yrb==openbc:
        allocatesf(b.edgexryr,stencil)
        init_splitfield(b.edgexryr.syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, yminp, zmin0,
                        clight, 1, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xrb==openbc:
            if yrb==periodic:
                b.edgexryr=b.sidexr
            if yrb==em3d.otherproc:
                allocatesf(b.edgexryr,stencil)
                b.edgexryr.proc=procyr
        if yrb==openbc:
            if xrb==periodic:
                b.edgexryr=b.sideyr
            if xrb==em3d.otherproc:
                allocatesf(b.edgexryr,stencil)
                b.edgexryr.proc=procxr
# xz
    b.edgexlzl = EM3D_FIELDtype()
    if xlb==openbc and zlb==openbc:
        allocatesf(b.edgexlzl,stencil)
        init_splitfield(b.edgexlzl.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, ymin0, zminm,
                        clight,-1, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xlb==openbc:
            if zlb==periodic:
                b.edgexlzl=b.sidexl
            if zlb==em3d.otherproc:
                allocatesf(b.edgexlzl,stencil)
                b.edgexlzl.proc=proczl
        if zlb==openbc:
            if xlb==periodic:
                b.edgexlzl=b.sidezl
            if xlb==em3d.otherproc:
                allocatesf(b.edgexlzl,stencil)
                b.edgexlzl.proc=procxl

    b.edgexrzl = EM3D_FIELDtype()
    if xrb==openbc and zlb==openbc:
        allocatesf(b.edgexrzl,stencil)
        init_splitfield(b.edgexrzl.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, ymin0, zminm,
                        clight, 1, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xrb==openbc:
            if zlb==periodic:
                b.edgexrzl=b.sidexr
            if zlb==em3d.otherproc:
                allocatesf(b.edgexrzl,stencil)
                b.edgexrzl.proc=proczl
        if zlb==openbc:
            if xrb==periodic:
                b.edgexrzl=b.sidezl
            if xrb==em3d.otherproc:
                allocatesf(b.edgexrzl,stencil)
                b.edgexrzl.proc=procxr

    b.edgexlzr = EM3D_FIELDtype()
    if xlb==openbc and zrb==openbc:
        allocatesf(b.edgexlzr,stencil)
        init_splitfield(b.edgexlzr.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, ymin0, zminp,
                        clight,-1, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xlb==openbc:
            if zrb==periodic:
                b.edgexlzr=b.sidexl
            if zrb==em3d.otherproc:
                allocatesf(b.edgexlzr,stencil)
                b.edgexlzr.proc=proczr
        if zrb==openbc:
            if xlb==periodic:
                b.edgexlzr=b.sidezr
            if xlb==em3d.otherproc:
                allocatesf(b.edgexlzr,stencil)
                b.edgexlzr.proc=procxl

    b.edgexrzr = EM3D_FIELDtype()
    if xrb==openbc and zrb==openbc:
        allocatesf(b.edgexrzr,stencil)
        init_splitfield(b.edgexrzr.syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, ymin0, zminp,
                        clight, 1, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xrb==openbc:
            if zrb==periodic:
                b.edgexrzr=b.sidexr
            if zrb==em3d.otherproc:
                allocatesf(b.edgexrzr,stencil)
                b.edgexrzr.proc=proczr
        if zrb==openbc:
            if xrb==periodic:
                b.edgexrzr=b.sidezr
            if xrb==em3d.otherproc:
                allocatesf(b.edgexrzr,stencil)
                b.edgexrzr.proc=procxr

# yz
    b.edgeylzl = EM3D_FIELDtype()
    if ylb==openbc and zlb==openbc:
        allocatesf(b.edgeylzl,stencil)
        init_splitfield(b.edgeylzl.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xmin0, yminm, zminm,
                        clight, 0,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if ylb==openbc:
            if zlb==periodic:
                b.edgeylzl=b.sideyl
            if zlb==em3d.otherproc:
                allocatesf(b.edgeylzl,stencil)
                b.edgeylzl.proc=proczl
        if zlb==openbc:
            if ylb==periodic:
                b.edgeylzl=b.sidezl
            if ylb==em3d.otherproc:
                allocatesf(b.edgeylzl,stencil)
                b.edgeylzl.proc=procyl
    b.edgeyrzl = EM3D_FIELDtype()
    if yrb==openbc and zlb==openbc:
        allocatesf(b.edgeyrzl,stencil)
        init_splitfield(b.edgeyrzl.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xmin0, yminp, zminm,
                        clight, 0, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if yrb==openbc:
            if zlb==periodic:
                b.edgeyrzl=b.sideyr
            if zlb==em3d.otherproc:
                allocatesf(b.edgeyrzl,stencil)
                b.edgeyrzl.proc=proczl
        if zlb==openbc:
            if yrb==periodic:
                b.edgeyrzl=b.sidezl
            if yrb==em3d.otherproc:
                allocatesf(b.edgeyrzl,stencil)
                b.edgeyrzl.proc=procyr
    b.edgeylzr = EM3D_FIELDtype()
    if ylb==openbc and zrb==openbc:
        allocatesf(b.edgeylzr,stencil)
        init_splitfield(b.edgeylzr.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xmin0, yminm, zminp,
                        clight, 0,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if ylb==openbc:
            if zrb==periodic:
                b.edgeylzr=b.sideyl
            if zrb==em3d.otherproc:
                allocatesf(b.edgeylzr,stencil)
                b.edgeylzr.proc=proczr
        if zrb==openbc:
            if ylb==periodic:
                b.edgeylzr=b.sidezr
            if ylb==em3d.otherproc:
                allocatesf(b.edgeylzr,stencil)
                b.edgeylzr.proc=procyl
    b.edgeyrzr = EM3D_FIELDtype()
    if yrb==openbc and zrb==openbc:
        allocatesf(b.edgeyrzr,stencil)
        init_splitfield(b.edgeyrzr.syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xmin0, yminp, zminp,
                        clight, 0, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if yrb==openbc:
            if zrb==periodic:
                b.edgeyrzr=b.sideyr
            if zrb==em3d.otherproc:
                allocatesf(b.edgeyrzr,stencil)
                b.edgeyrzr.proc=proczr
        if zrb==openbc:
            if yrb==periodic:
                b.edgeyrzr=b.sidezr
            if yrb==em3d.otherproc:
                allocatesf(b.edgeyrzr,stencil)
                b.edgeyrzr.proc=procyr

# *** corners
    b.cornerxlylzl = EM3D_FIELDtype()
    if xlb==openbc and ylb==openbc and zlb==openbc:
        allocatesf(b.cornerxlylzl,stencil)
        init_splitfield(b.cornerxlylzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, yminm, zminm,
                        clight,-1,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xlb==openbc and ylb==openbc:
            if zlb==periodic:
                b.cornerxlylzl=b.edgexlyl
            if zlb==em3d.otherproc:
                allocatesf(b.cornerxlylzl,stencil)
                b.cornerxlylzl.proc=proczl
        if xlb==openbc and zlb==openbc:
            if ylb==periodic:
                b.cornerxlylzl=b.edgexlzl
            if ylb==em3d.otherproc:
                allocatesf(b.cornerxlylzl,stencil)
                b.cornerxlylzl.proc=procyl
        if ylb==openbc and zlb==openbc:
            if xlb==periodic:
                b.cornerxlylzl=b.edgeylzl
            if xlb==em3d.otherproc:
                allocatesf(b.cornerxlylzl,stencil)
                b.cornerxlylzl.proc=procxl
    b.cornerxrylzl = EM3D_FIELDtype()
    if xrb==openbc and ylb==openbc and zlb==openbc:
        allocatesf(b.cornerxrylzl,stencil)
        init_splitfield(b.cornerxrylzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, yminm, zminm,
                        clight, 1,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xrb==openbc and ylb==openbc:
            if zlb==periodic:
                b.cornerxrylzl=b.edgexryl
            if zlb==em3d.otherproc:
                allocatesf(b.cornerxrylzl,stencil)
                b.cornerxrylzl.proc=proczl
        if xrb==openbc and zlb==openbc:
            if ylb==periodic:
                b.cornerxrylzl=b.edgexrzl
            if ylb==em3d.otherproc:
                allocatesf(b.cornerxrylzl,stencil)
                b.cornerxrylzl.proc=procyl
        if ylb==openbc and zlb==openbc:
            if xrb==periodic:
                b.cornerxrylzl=b.edgeylzl
            if xrb==em3d.otherproc:
                allocatesf(b.cornerxrylzl,stencil)
                b.cornerxrylzl.proc=procxr
    b.cornerxlyrzl = EM3D_FIELDtype()
    if xlb==openbc and yrb==openbc and zlb==openbc:
        allocatesf(b.cornerxlyrzl,stencil)
        init_splitfield(b.cornerxlyrzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, yminp, zminm,
                        clight,-1, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xlb==openbc and yrb==openbc:
            if zlb==periodic:
                b.cornerxlyrzl=b.edgexlyr
            if zlb==em3d.otherproc:
                allocatesf(b.cornerxlyrzl,stencil)
                b.cornerxlyrzl.proc=proczl
        if xlb==openbc and zlb==openbc:
            if yrb==periodic:
                b.cornerxlyrzl=b.edgexlzl
            if yrb==em3d.otherproc:
                allocatesf(b.cornerxlyrzl,stencil)
                b.cornerxlyrzl.proc=procyr
        if yrb==openbc and zlb==openbc:
            if xlb==periodic:
                b.cornerxlyrzl=b.edgeyrzl
            if xlb==em3d.otherproc:
                allocatesf(b.cornerxlyrzl,stencil)
                b.cornerxlyrzl.proc=procxl
    b.cornerxryrzl = EM3D_FIELDtype()
    if xrb==openbc and yrb==openbc and zlb==openbc:
        allocatesf(b.cornerxryrzl,stencil)
        init_splitfield(b.cornerxryrzl.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, yminp, zminm,
                        clight, 1, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xrb==openbc and yrb==openbc:
            if zlb==periodic:
                b.cornerxryrzl=b.edgexryr
            if zlb==em3d.otherproc:
                allocatesf(b.cornerxryrzl,stencil)
                b.cornerxryrzl.proc=proczl
        if xrb==openbc and zlb==openbc:
            if yrb==periodic:
                b.cornerxryrzl=b.edgexrzl
            if yrb==em3d.otherproc:
                allocatesf(b.cornerxryrzl,stencil)
                b.cornerxryrzl.proc=procyr
        if yrb==openbc and zlb==openbc:
            if xrb==periodic:
                b.cornerxryrzl=b.edgeyrzl
            if xrb==em3d.otherproc:
                allocatesf(b.cornerxryrzl,stencil)
                b.cornerxryrzl.proc=procxr
    b.cornerxlylzr = EM3D_FIELDtype()
    if xlb==openbc and ylb==openbc and zrb==openbc:
        allocatesf(b.cornerxlylzr,stencil)
        init_splitfield(b.cornerxlylzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, yminm, zminp,
                        clight,-1,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xlb==openbc and ylb==openbc:
            if zrb==periodic:
                b.cornerxlylzr=b.edgexlyl
            if zrb==em3d.otherproc:
                allocatesf(b.cornerxlylzr,stencil)
                b.cornerxlylzr.proc=proczr
        if xlb==openbc and zrb==openbc:
            if ylb==periodic:
                b.cornerxlylzr=b.edgexlzr
            if ylb==em3d.otherproc:
                allocatesf(b.cornerxlylzr,stencil)
                b.cornerxlylzr.proc=procyl
        if ylb==openbc and zrb==openbc:
            if xlb==periodic:
                b.cornerxlylzr=b.edgeylzr
            if xlb==em3d.otherproc:
                allocatesf(b.cornerxlylzr,stencil)
                b.cornerxlylzr.proc=procxl
    b.cornerxrylzr = EM3D_FIELDtype()
    if xrb==openbc and ylb==openbc and zrb==openbc:
        allocatesf(b.cornerxrylzr,stencil)
        init_splitfield(b.cornerxrylzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, yminm, zminp,
                        clight, 1,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xrb==openbc and ylb==openbc:
            if zrb==periodic:
                b.cornerxrylzr=b.edgexryl
            if zrb==em3d.otherproc:
                allocatesf(b.cornerxrylzr,stencil)
                b.cornerxrylzr.proc=proczr
        if xrb==openbc and zrb==openbc:
            if ylb==periodic:
                b.cornerxrylzr=b.edgexrzr
            if ylb==em3d.otherproc:
                allocatesf(b.cornerxrylzr,stencil)
                b.cornerxrylzr.proc=procyl
        if ylb==openbc and zrb==openbc:
            if xrb==periodic:
                b.cornerxrylzr=b.edgeylzr
            if xrb==em3d.otherproc:
                allocatesf(b.cornerxrylzr,stencil)
                b.cornerxrylzr.proc=procxr
    b.cornerxlyrzr = EM3D_FIELDtype()
    if xlb==openbc and yrb==openbc and zrb==openbc:
        allocatesf(b.cornerxlyrzr,stencil)
        init_splitfield(b.cornerxlyrzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminm, yminp, zminp,
                        clight,-1, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xlb==openbc and yrb==openbc:
            if zrb==periodic:
                b.cornerxlyrzr=b.edgexlyr
            if zrb==em3d.otherproc:
                allocatesf(b.cornerxlyrzr,stencil)
                b.cornerxlyrzr.proc=proczr
        if xlb==openbc and zrb==openbc:
            if yrb==periodic:
                b.cornerxlyrzr=b.edgexlzr
            if yrb==em3d.otherproc:
                allocatesf(b.cornerxlyrzr,stencil)
                b.cornerxlyrzr.proc=procyr
        if yrb==openbc and zrb==openbc:
            if xlb==periodic:
                b.cornerxlyrzr=b.edgeyrzr
            if xlb==em3d.otherproc:
                allocatesf(b.cornerxlyrzr,stencil)
                b.cornerxlyrzr.proc=procxl
    b.cornerxryrzr = EM3D_FIELDtype()
    if xrb==openbc and yrb==openbc and zrb==openbc:
        allocatesf(b.cornerxryrzr,stencil)
        init_splitfield(b.cornerxryrzr.syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, xminp, yminp, zminp,
                        clight, 1, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz,
                        f.norderx,f.nordery,f.norderz,f.xcoefs,f.ycoefs,f.zcoefs,
                        l_nodalgrid,pml_method)
    else:
        if xrb==openbc and yrb==openbc:
            if zrb==periodic:
                b.cornerxryrzr=b.edgexryr
            if zrb==em3d.otherproc:
                allocatesf(b.cornerxryrzr,stencil)
                b.cornerxryrzr.proc=proczr
        if xrb==openbc and zrb==openbc:
            if yrb==periodic:
                b.cornerxryrzr=b.edgexrzr
            if yrb==em3d.otherproc:
                allocatesf(b.cornerxryrzr,stencil)
                b.cornerxryrzr.proc=procyr
        if yrb==openbc and zrb==openbc:
            if xrb==periodic:
                b.cornerxryrzr=b.edgeyrzr
            if xrb==em3d.otherproc:
                allocatesf(b.cornerxryrzr,stencil)
                b.cornerxryrzr.proc=procxr

    return b

def FD_weights(z,n,m):
    """
 adapted from Matlab code from Fornberg (1998)
 Calculates FD weights. The parameters are:
  z   location where approximations are to be accurate,
  n   number of grid points,
  m   highest derivative that we want to find weights for
  c   array size m+1,lentgh(x) containing (as output) in
      successive rows the weights for derivatives 0,1,...,m.
    """

    x = arange(-n/2+1,n/2+1)*1.

    c=zeros([m+1,n]); c1=1.; c4=x[0]-z; c[0,0]=1.;
    for i in range(1,n):
        mn=min(i+1,m+1); c2=1.; c5=c4; c4=x[i]-z;
        for j in range(0,i-0):
            c3=x[i]-x[j];  c2=c2*c3;
            if j==i-1:
                c[1:mn,i]=c1*(arange(1,mn)*c[0:mn-1,i-1]-c5*c[1:mn,i-1])/c2;
                c[0,i]=-c1*c5*c[0,i-1]/c2;
            c[1:mn,j]=(c4*c[1:mn,j]-arange(1,mn)*c[0:mn-1,j])/c3;
            c[0,j]=c4*c[0,j]/c3;
        c1=c2;

    return c

def FD_weights_hvincenti(p,l_staggered=False):
    # --- from Henri Vincenti's formulas
    factorial = math.factorial

    c = zeros(int(p/2))
    for i in range(int(p/2)):
        l=i+1
        if l_staggered:
            lognumer = math.log(16.)*(1.-p/2.)+math.log(factorial(p-1.))*2
            logdenom = math.log(2.*l-1.)*2.+math.log(factorial(p/2.+l-1.))+math.log(factorial(p/2.-l))+2*math.log(factorial(p/2.-1.))
        else:
            lognumer = math.log(factorial(p/2.))*2
            logdenom = math.log(factorial(p/2.+l))+math.log(factorial(p/2.-l))+math.log(l)
        c[i] = (-1.)**(l+1)*exp(lognumer-logdenom)
    return c



##############################################################################
# --- This can only be done after the class is defined.
#try:
#  psyco.bind(EM3D)
#except NameError:
#  pass
