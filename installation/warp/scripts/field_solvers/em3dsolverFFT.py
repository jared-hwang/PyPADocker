"""Class for 2D & 3D FFT-based electromagnetic solver """
from GPSTD import *
from em3dsolver import *

class EM3DFFT(EM3D):

    __em3dfftinputs__ = []
    __flaginputs__ = {'spectral':1, # on/off spectral solver
                      'spectral_current':0, # flag for deposition that generalizes Esirkepov in k-space for any order
                      'current_cor':1,  # flag for correction current to verify Gauss' Law (not needed if spectral_current=1; more local than Boris correction on fields)
                      'boris_cor':0, # flag for correction of fields verify Gauss' Law (i.e. Boris correction; not needed if spectral_current=1 or spectral_current=1)
                      'spectral_mix':0.,'Jmult':False,
                      'l_spectral_staggered':False,
                      'l_theta_with_modif_k':True, # Whether to use the modified (finite-order) k in the Galilean Theta coefficients
                      'sigmab_x':0.,'sigmab_y':0.,'sigmab_z':0.}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(EM3DFFT.__flaginputs__,kw)
        EM3D.__init__(self,kwdict=kw)
        self.l_spectral_staggered = not self.l_nodalgrid

        self.l_getrho=True
        self.pml_method=2

    def finalize(self,lforce=False):
        if self.finalized and not lforce: return
        EM3D.finalize(self)
        self.allocatefieldarraysFFT()

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
                     'V_galilean':s.V_galilean,\
                     'V_pseudogalilean':s.V_pseudogalilean,\
                     'clight':clight}

            if s.ntsub is np.inf:
                if not self.l_getrho:
                    self.l_getrho = True
                    f.nxr = f.nx
                    f.nyr = f.ny
                    f.nzr = f.nz
                    f.gchange()

                self.GPSTDMaxwell = PSATD_Maxwell(yf=self.fields,
                            eps0=eps0, bc_periodic=bc_periodic,
                            l_theta_with_modif_k=self.l_theta_with_modif_k,
                            **kwGPSTD)
            else:
                if self.l_pushf and not self.l_getrho:
                    self.l_getrho = True
                    f.nxr = f.nx
                    f.nyr = f.ny
                    f.nzr = f.nz
                    f.gchange()
                self.GPSTDMaxwell = GPSTD_Maxwell(yf=self.fields,
                                                  eps0=eps0,
                                                  bc_periodic=bc_periodic,
                                                  l_theta_with_modif_k=self.l_theta_with_modif_k,
                                                  **kwGPSTD)

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

        if self.spectral:
            kwPML = kwGPSTD

            if s.ntsub==inf:
                GPSTD_PML = PSATD_Maxwell_PML
            else:
                GPSTD_PML = GPSTD_Maxwell_PML

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


################################################################################
#                                   LASER
################################################################################

    def depose_rho_old_laser(self,f, laser_antenna_i):
        for q in [1.,-1.]:# q represents the sign of the charged macroparticles
            # The antenna is made of two types of fictious particles : positive
            # and negative
            # --- deposit dRho/dt for correction of current to verify
            # --- the discrete continuity equation dRho/dt+div J=0.

            if top.wpid==0:
                wfact = ones((1,),'d')
                l_particles_weight = false
            else:
                l_particles_weight = true
            # computes position shift of grid in one time step with
            # Galilean frame
            dxgal = self.V_galilean[0]*top.dt
            dygal = self.V_galilean[1]*top.dt
            dzgal = self.V_galilean[2]*top.dt

            if self.l_2dxz:
                depose_rhoold_n_2dxz(
                        self.fields.Rhoold_local,laser_antenna_i.nn,
                        laser_antenna_i.xx + q*laser_antenna_i.xdx,
                        laser_antenna_i.zz + q*laser_antenna_i.zdz,
                        laser_antenna_i.vx + q*laser_antenna_i.ux,
                        laser_antenna_i.vy + q*laser_antenna_i.uy,
                        laser_antenna_i.vz + q*laser_antenna_i.uz,
                        laser_antenna_i.gi,
                        laser_antenna_i.weights, q,
                        f.xmin-dxgal, f.zmin+self.zgrid-dzgal,
                        top.dt,
                        f.dx,f.dz,
                        f.nx,f.nz,
                        f.nxguard,f.nzguard,
                        self.laser_depos_order_x,
                        self.laser_depos_order_z,
                        l_particles_weight, w3d.l4symtry)
            else:
                depose_rhoold_n_3d(
                        self.fields.Rhoold_local,laser_antenna_i.nn,
                        laser_antenna_i.xx + q*laser_antenna_i.xdx,
                        laser_antenna_i.yy + q*laser_antenna_i.ydy,
                        laser_antenna_i.zz + q*laser_antenna_i.zdz,
                        laser_antenna_i.vx + q*laser_antenna_i.ux,
                        laser_antenna_i.vy + q*laser_antenna_i.uy,
                        laser_antenna_i.vz + q*laser_antenna_i.uz,
                        laser_antenna_i.gi,
                        laser_antenna_i.weights, q,
                        f.xmin-dxgal,f.ymin-dygal,f.zmin+self.zgrid-dzgal,
                        top.dt,
                        f.dx,f.dy,f.dz,
                        f.nx,f.ny,f.nz,
                        f.nxguard,f.nyguard,f.nzguard,
                        self.laser_depos_order_x,
                        self.laser_depos_order_y,
                        self.laser_depos_order_z,
                        l_particles_weight, w3d.l4symtry)

################################################################################
# CHARGE/CURRENT DEPOSITION
################################################################################

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

        if self.spectral_current:
            self.depose_current_density_spectral(n,f,x,y,z,ux,uy,uz,gaminv,dt,
                                        wfact,zgrid,
                                        q,w,nox,noy,noz,l_particles_weight)
        else:
            EM3D.depose_current_density(self,n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,
                                        q,w,nox,noy,noz,l_particles_weight)

        if self.l_getrho :
            EM3D.depose_charge_density(self,n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight)

        if self.current_cor:
            # --- deposit dRho/dt for correction of current to verify
            # --- the discrete continuity equation dRho/dt+div J=0.
            if top.wpid==0:
                wfact = ones((1,),'d')
                l_particles_weight = false
            else:
                l_particles_weight = true
            # --- computes position shift of grid in one time step with Galilean frame
            dxgal = self.V_galilean[0]*top.dt
            dygal = self.V_galilean[1]*top.dt
            dzgal = self.V_galilean[2]*top.dt
            if self.l_2dxz:
                depose_rhoold_n_2dxz(self.fields.Rhoold_local,n,
                                                    x,z,ux,uy,uz,
                                                    gaminv,wfact,q*w,
                                                    f.xmin-dxgal,f.zmin+self.zgrid-dzgal,
                                                    top.dt*top.pgroup.ndts[js],
                                                    f.dx,f.dz,
                                                    f.nx,f.nz,
                                                    f.nxguard,f.nzguard,
                                                    top.depos_order[0,js],
                                                    top.depos_order[2,js],
                                                    l_particles_weight,w3d.l4symtry)
            else:
                depose_rhoold_n_3d(self.fields.Rhoold_local,n,
                                                    x,y,z,ux,uy,uz,
                                                    gaminv,wfact,q*w,
                                                    f.xmin-dxgal,f.ymin-dygal,f.zmin+self.zgrid-dzgal,
                                                    top.dt*top.pgroup.ndts[js],
                                                    f.dx,f.dy,f.dz,
                                                    f.nx,f.ny,f.nz,
                                                    f.nxguard,f.nyguard,f.nzguard,
                                                    top.depos_order[0,js],
                                                    top.depos_order[1,js],
                                                    top.depos_order[2,js],
                                                    l_particles_weight,w3d.l4symtry)

    def depose_current_density(self,n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight):
        if not self.spectral_current:
            EM3D.depose_current_density(self,n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight)
        else:
            self.depose_current_density_spectral(n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight)

    def depose_current_density_spectral(self,n,f,x,y,z,ux,uy,uz,gaminv,dt,wfact,zgrid,q,w,nox,noy,noz,l_particles_weight):
        if self.l_1dz:
            raise Exception('Need to add spectral current deposition in 1-D')
        elif self.l_2dxz:
            jx = self.fields.Jx[:,self.fields.nyguard,:]
            jy = self.fields.Jy[:,self.fields.nyguard,:]
            jz = self.fields.Jz[:,self.fields.nyguard,:]
            depose_j_n_2dxz_spectral(jx,jy,jz,n,
                                                            x,z,ux,uy,uz,
                                                            gaminv,wfact,q*w,
                                                            f.xmin-0.5*top.dt*self.V_galilean[0],
                                                            f.zmin+self.zgrid-0.5*top.dt*self.V_galilean[2],
                                                            dt,
                                                            f.dx,f.dz,
                                                            f.nx,f.nz,
                                                            f.nxguard,f.nzguard,
                                                            nox,
                                                            noz,
                                                            l_particles_weight,w3d.l4symtry)
        else:
            raise Exception('Need to add spectral current deposition in 3-D')

    def zerosourcep(self):
        EM3D.zerosourcep(self)

        # --- zero proper portion of Jarray
        for indts in range(top.nsndts):
            if top.ldts[indts]:
                if self.current_cor:
                    self.fields.Rhoold_local[...] = 0.

    def finalizesourcep(self):
        if self.sourcepfinalized: return
        self.sourcepfinalized = 1
        if self.l_verbose:print 'finalizesourcep'
        # --- add slices
        self.add_source_ndts_slices()
        self.aftersetsourcep()
        # -- add laser, loop on the different antennas
        for i in range(len(self.laser_antenna)):
            self.add_laser(self.block.core.yf, self.laser_antenna[i])
        if not self.spectral_current:self.Jyee2node3d()
        if self.spectral_current:self.getcurrent_spectral()
        # --- smooth current density
        if any(self.npass_smooth>0):self.smoothdensity()
        # --- apply current correction (before exchange between domains for locality)
        if self.current_cor:
            self.current_cor_spectral()
        # --- apply boundary conditions
        self.applysourceboundaryconditions()
        # --- exchange guard cells data across domains
        self.exchange_bc(self.block)
        if self.l_verbose:print 'finalizesourcep done'

    def boris_correction(self):
# --- apply Boris correction
        l_plotdive=1
        emK = self.FSpace
        em = self
        f = self.fields
        ixl,ixu,iyl,iyu,izl,izu = emK.get_ius()
        fields_shape = [ixu-ixl,iyu-iyl,izu-izl]

        j = 1j

        if emK.nx>1:ExF = emK.rfftn(squeeze(f.Ex[ixl:ixu,iyl:iyu,izl:izu]))
        if emK.ny>1:EyF = emK.rfftn(squeeze(f.Ey[ixl:ixu,iyl:iyu,izl:izu]))
        if emK.nz>1:EzF = emK.rfftn(squeeze(f.Ez[ixl:ixu,iyl:iyu,izl:izu]))

        RhoF = emK.rfftn(squeeze(f.Rho[ixl:ixu,iyl:iyu,izl:izu]))

        divemrho = -j*RhoF/eps0
        if emK.nx>1: divemrho -= emK.kxm*ExF
        if emK.ny>1: divemrho -= emK.kym*EyF
        if emK.nz>1: divemrho -= emK.kzm*EzF

        if l_plotdive:
                window(4);fma();ppg(abs(divemrho)*clight*top.dt,view=3);
                ppg(abs(emK.irfftn(divemrho)),view=4)
#                ppg(abs(j*RhoF*clight*top.dt/eps0),view=4)
                refresh()

        if emK.nx>1: ExF+=divemrho*emK.kxpn/emK.kmag
        if emK.ny>1: EyF+=divemrho*emK.kypn/emK.kmag
        if emK.nz>1: EzF+=divemrho*emK.kzpn/emK.kmag

        if l_plotdive:
                divemrho = -(emK.kxm*ExF+emK.kzm*EzF)-j*RhoF/eps0
                divemrho *= clight*top.dt
                ppg(abs(divemrho),view=5);ppg(abs(j*RhoF*clight*top.dt/eps0),view=6)

        if emK.nx>1:
            Ex = emK.irfftn(ExF)
            Ex.resize(fields_shape)
            f.Ex[ixl:ixu,iyl:iyu,izl:izu] = Ex.real
        if emK.ny>1:
            Ey = emK.irfftn(EyF)
            Ey.resize(fields_shape)
            f.Ey[ixl:ixu,iyl:iyu,izl:izu] = Ey.real
        if emK.nx>1:
            Ez = emK.irfftn(EzF)
            Ez.resize(fields_shape)
            f.Ez[ixl:ixu,iyl:iyu,izl:izu] = Ez.real


    def boris_correctionold(self):
# --- apply Boris correction
        l_plotdive=1
        em = self = getregisteredsolver()
        emK = self.GPSTDMaxwell
        j = 1j
        if self.boris_cor:
            ExF = emK.Ffields['ex']
            EzF = emK.Ffields['ez']
            RhoF = emK.Ffields['rho']
            divemrho = -(emK.kxm*ExF+emK.kzm*EzF)-j*RhoF/eps0
            if l_plotdive:
                window(4);fma();ppg(abs(divemrho),view=3);ppg(abs(j*RhoF),view=4)
            ExF+=divemrho*emK.kxpn/emK.kmag
            EzF+=divemrho*emK.kzpn/emK.kmag
            if l_plotdive:
                divemrho = -(emK.kxm*ExF+emK.kzm*EzF)-j*RhoF/eps0
                ppg(abs(divemrho),view=5);ppg(abs(j*RhoF),view=6)

    def wrap_periodic_BC(self,flist=[]):
        emK = self.FSpace

        if self.bounds[0]==periodic:
            ngx = self.nxguard
        else:
            ngx = 0
        if self.bounds[2]==periodic:
            ngy = self.nyguard
        else:
            ngy = 0
        if self.bounds[4]==periodic:
            ngz = self.nzguard
        else:
            ngz = 0

        if emK.nx>1 and self.bounds[0]==periodic:
            for J in flist:
                J[ngx:2*ngx+1,...]+=J[-ngx-1:,...]
                J[-ngx-1:,...]=0.
                J[-2*ngx-1:-ngx-1,...]+=J[:ngx,...]
                J[:ngx,...]=0.
        if emK.ny>1 and self.bounds[2]==periodic:
            for J in flist:
                J[:,ngy:2*ngy+1,:]+=J[:,-ngy-1:,:]
                J[:,-ngy-1:,:]=0.
                J[:,-2*ngy-1:-ngy-1,:]+=J[:,:ngy,:]
                J[:,:ngy,:]=0.
        if emK.nz>1 and self.bounds[4]==periodic:
            for J in flist:
                J[...,ngz:2*ngz+1]+=J[...,-ngz-1:]
                J[...,-ngz-1:]=0.
                J[...,-2*ngz-1:-ngz-1]+=J[...,:ngz]
                J[...,:ngz]=0.

    def current_cor_spectral(self):
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

        if all(self.V_galilean==0.) and all(self.V_pseudogalilean==0.):
            em.dRhoodtF = emK.rfftn(squeeze((f.Rho-f.Rhoold_local)[ixl:ixu,iyl:iyu,izl:izu]/top.dt),plan=emK.planj_rfftn)
        else:
            em.RhoF = emK.rfftn(squeeze(f.Rho[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)
            em.RhooldF = emK.rfftn(squeeze(f.Rhoold_local[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)

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

        if all(self.V_galilean==0.) and all(self.V_pseudogalilean==0.):
            if emK.nx>1:
                Jxl = j*em.dRhoodtF*emK.kxpn/emK.kmag
            if emK.ny>1:
                Jyl = j*em.dRhoodtF*emK.kypn/emK.kmag
            if emK.nz>1:
                Jzl = j*em.dRhoodtF*emK.kzpn/emK.kmag
        else:
            if emK.nx>1:
                Jxl = em.RhoF*emK.JxCorRhomult-em.RhooldF*emK.JxCorRhooldmult
            if emK.ny>1:
                Jyl = em.RhoF*emK.JyCorRhomult-em.RhooldF*emK.JyCorRhooldmult
            if emK.nz>1:
                Jzl = em.RhoF*emK.JzCorRhomult-em.RhooldF*emK.JzCorRhooldmult

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

    def getcurrent_spectral(self):
        if not self.spectral_current:return
        j=1j      # imaginary number
        emK = self.FSpace
        f = self.fields
        ixl,ixu,iyl,iyu,izl,izu = emK.get_ius()
        fields_shape = [ixu-ixl,iyu-iyl,izu-izl]

        self.wrap_periodic_BC([f.Jx,f.Jy,f.Jz])

        if emK.planj_rfftn is None:
            emK.planj_rfftn= emK.create_plan_rfftn(np.asarray(fields_shape))
        if emK.planj_irfftn is None:
            emK.planj_irfftn= emK.create_plan_irfftn(np.asarray(fields_shape))

        if emK.nx>1:JxF = emK.rfftn(squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)
#        if emK.ny>1:JyF = emK.fft(squeeze(f.Jy[ixl:ixu,iyl:iyu,izl:izu]))
        if emK.nz>1:JzF = emK.rfftn(squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)

#         if emK.nx>1:JxF = emK.fft(squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]),axis=0)
# #        if emK.ny>1:JyF = emK.fft(squeeze(f.Jy[ixl:ixu,iyl:iyu,izl:izu]))
#         if emK.nz>1:JzF = emK.fft(squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]),axis=1)


        if self.l_spectral_staggered:
            JxF = 1j*JxF/where(emK.kx==0.,1j,emK.kx*exp(-j*emK.kx_unmod*self.dx/2))
            JzF = 1j*JzF/where(emK.kz==0.,1j,emK.kz*exp(-j*emK.kz_unmod*self.dz/2))
        else:
            JxF = 1j*JxF/where(emK.kx==0.,1j,emK.kx)
            JzF = 1j*JzF/where(emK.kz==0.,1j,emK.kz)

        if any(self.V_pseudogalilean<>0.):
            JxF*=emK.CDcoef
#            JyF*=emK.CDcoef
            JzF*=emK.CDcoef

        # --- sets currents to zero at the Nyquist wavelength
        if emK.nx>1:
            dkx = 2*pi/(self.nxlocal+2*self.nxguard)
            index_Nyquist = compress(abs(emK.kxunit*self.dx)>=pi-0.1*dkx,arange(self.nxlocal+2*self.nxguard))
            JxF[index_Nyquist,...]=0.

        if emK.ny>1:
            dky = 2*pi/(self.nylocal+2*self.nyguard)
            index_Nyquist = compress(abs(emK.kyunit*self.dy)>=pi-0.1*dky,arange(self.nylocal+2*self.nyguard))
            JyF[:,index_Nyquist,:]=0.

        if emK.nz>1:
            dkz = 2*pi/(self.nzlocal+2*self.nzguard)
            index_Nyquist = compress(abs(emK.kzunit*self.dz)>=pi-0.1*dkz,arange(self.nzlocal+2*self.nzguard))
            JzF[...,index_Nyquist]=0.

        # --- adjust average current values
        if emK.nx>1:
            Jxcs=ave(cumsum(f.Jx[ixl:ixu,iyl:iyu,izl:izu],0),0)
#             Jx = emK.ifft(JxF,axis=0)
            Jx = emK.irfftn(JxF, np.asarray(np.shape(squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]))), plan=emK.planj_irfftn)#, field_out=squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]))
            Jx.resize(fields_shape)
            Jx=Jx.real
            for i in range(shape(Jx)[1]):
                Jx[:,i]-=Jxcs[i]*self.dx

        if emK.nz>1:
            Jzcs=ave(cumsum(f.Jz[ixl:ixu,iyl:iyu,izl:izu],2),2)
#             Jz = emK.ifft(JzF,axis=1)
            Jz = emK.irfftn(JzF, np.asarray(np.shape(squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]))), plan=emK.planj_irfftn)#, field_out=squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]))
            Jz.resize(fields_shape)
            Jz=Jz.real
            for i in range(shape(Jz)[0]):
                Jz[i,:]-=Jzcs[i]*self.dz

        if emK.nx>1:
            f.Jx[ixl:ixu,iyl:iyu,izl:izu]=Jx
        if emK.nz>1:
            f.Jz[ixl:ixu,iyl:iyu,izl:izu]=Jz

#        self.JxF=JxF
#        self.JzF=JzF

        del Jx, Jz, Jxcs, Jzcs, JxF, JzF

    def smoothdensity(self):
        if all(self.npass_smooth==0):return
        EM3D.smoothdensity(self)
        nx,ny,nz = shape(self.fields.Jx)
        nsm = shape(self.npass_smooth)[1]
        l_mask_method=1
        if self.mask_smooth is None:
            self.mask_smooth=[]
            for js in range(nsm):
                self.mask_smooth.append(None)

        if l_mask_method==2:
            if self.current_cor:
                DRhooldcopy = self.fields.Rhoold_local.copy()
        for js in range(nsm):
            if self.mask_smooth[js] is None or l_mask_method==2:
                if self.current_cor:
                    smooth3d_121_stride(self.fields.Rhoold_local[...],nx-1,ny-1,nz-1,
                                        self.npass_smooth[:,js].copy(),
                                        self.alpha_smooth[:,js].copy(),
                                        self.stride_smooth[:,js].copy())
            else:
                if self.current_cor:
                    smooth3d_121_stride_mask(self.fields.Rhoold_local[...],self.mask_smooth[js],nx-1,ny-1,nz-1,
                                           self.npass_smooth[:,js].copy(),
                                           self.alpha_smooth[:,js].copy(),
                                           self.stride_smooth[:,js].copy())
        if l_mask_method==2:
            if self.current_cor:
                self.fields.Rhoold_local *= self.mask_smooth
                self.fields.Rhoold_local += DRhooldcopy*(1.-self.mask_smooth)

    def Jyee2node3d(self):
        if self.l_nodalgrid and not self.l_deposit_nodal:
            Jyee2node3d(self.block.core.yf)
            if self.refinement is not None:
                self.__class__.__bases__[1].Jyee2node3d(self.field_coarse)

        if self.spectral and not self.l_nodalgrid and self.l_deposit_nodal:
            self.Jnode2yeeK()

    def apply_KFilter(self,F,KFilter):
        emK = self.FSpace
        ixl,ixu,iyl,iyu,izl,izu = emK.get_ius()
        KF = emK.rfftn(np.squeeze(F[ixl:ixu,iyl:iyu,izl:izu]))
        f = emK.irfftn(KF*KFilter,np.asarray(np.shape(np.squeeze(F[ixl:ixu,iyl:iyu,izl:izu]))))
        f.resize([ixu-ixl,iyu-iyl,izu-izl])
        F[ixl:ixu,iyl:iyu,izl:izu] = f.real
        del KF

    def smoothfields_poly(self):
        if not self.spectral:
            EM3D.smoothfields_poly(self)
            return
        else:
            self.apply_KFilter(self.fields.Exp,self.Exmultiplier)
            self.apply_KFilter(self.fields.Eyp,self.Exmultiplier)
            if self.ebcor==2:
                self.apply_KFilter(self.fields.Bxp,self.Bymultiplier)
                self.apply_KFilter(self.fields.Byp,self.Bymultiplier)

    def apply_KShift(self,F,KShift):
        emK = self.FSpace
        KF = emK.rfftn(F)
        F[...] = emK.irfftn(KF*KShift, np.asarray(np.shape(F))).real
        del KF
        return

    def fieldcenterK(self):
        emK = self.FSpace

        nx,ny,nz = shape(self.fields.Exp)
        kx_unmod = np.ones(shape(self.fields.Exp))
        kz_unmod = np.ones(shape(self.fields.Exp))
        kxunit = 2.*np.pi*(emK.rfftfreq(nx))/self.dx
        kzunit = 2.*np.pi*(emK.fftfreq(nz))/self.dz
        for i in range(nz):
            kx_unmod[:,0,i] *= kxunit
        for i in range(nx):
            kz_unmod[i,0,:] *= kzunit

        self.kx_unmod = kx_unmod
        self.kz_unmod = kz_unmod

        self.apply_KShift(self.fields.Exp,exp(-1j*kx_unmod*self.dx/2))
        self.apply_KShift(self.fields.Ezp,exp(-1j*kz_unmod*self.dz/2))
        self.apply_KShift(self.fields.Bxp,exp(-1j*(kz_unmod*self.dz/2)))
        self.apply_KShift(self.fields.Byp,exp(-1j*(kx_unmod*self.dx/2+kz_unmod*self.dz/2)))
        self.apply_KShift(self.fields.Bzp,exp(-1j*(kx_unmod*self.dx/2)))

        self.fields.l_nodecentered=True
        return
        self.apply_KFilter(self.fields.Exp,exp(-1j*kx_unmod*self.dx/2))
        self.apply_KFilter(self.fields.Eyp,exp(-1j*ky_unmod*self.dy/2))
        self.apply_KFilter(self.fields.Ezp,exp(-1j*kz_unmod*self.dz/2))
        self.apply_KFilter(self.fields.Bxp,exp(-1j*(ky_unmod*self.dy/2+kz_unmod*self.dz/2)))
        self.apply_KFilter(self.fields.Byp,exp(-1j*(kx_unmod*self.dx/2+kz_unmod*self.dz/2)))
        self.apply_KFilter(self.fields.Bzp,exp(-1j*(kx_unmod*self.dx/2+ky_unmod*self.dy/2)))

    def fieldcenterK(self):
        # --- centers fields at nodes by shifting in K-space
        emK = self.FSpace

        Cx = exp(-1j*emK.kx_unmod*self.dx/2)
        Cy = exp(-1j*emK.ky_unmod*self.dy/2)
        Cz = exp(-1j*emK.kz_unmod*self.dz/2)

        self.apply_KFilter(self.fields.Exp,Cx)
        self.apply_KFilter(self.fields.Eyp,Cy)
        self.apply_KFilter(self.fields.Ezp,Cz)
        self.apply_KFilter(self.fields.Bxp,Cy*Cz)
        self.apply_KFilter(self.fields.Byp,Cx*Cz)
        self.apply_KFilter(self.fields.Bzp,Cx*Cy)

        self.fields.l_nodecentered=True
        return

    def Jnode2yeeK(self):
        emK = self.FSpace

        nx,ny,nz = shape(self.fields.Exp)
        kx_unmod = np.ones(shape(self.fields.Exp))
        kz_unmod = np.ones(shape(self.fields.Exp))
        kxunit = 2.*np.pi*(emK.rfftfreq(nx))/self.dx
        kzunit = 2.*np.pi*(emK.fftfreq(nz))/self.dz
        for i in range(nz):
            kx_unmod[:,0,i] *= kxunit
        for i in range(nx):
            kz_unmod[i,0,:] *= kzunit

#        self.kx_unmod = kx_unmod
#        self.kz_unmod = kz_unmod

#        self.apply_KShift(self.fields.J[...,0],exp(1j*kx_unmod*self.dx/2))
#        self.apply_KShift(self.fields.J[...,1],exp(1j*ky_unmod*self.dz/2))
#        self.apply_KShift(self.fields.J[...,2],exp(1j*kz_unmod*self.dz/2))

        self.apply_KShift(self.fields.J[...,0],cos(kx_unmod*self.dx/2))
#        self.apply_KShift(self.fields.J[...,1],exp(1j*ky_unmod*self.dz/2))
        self.apply_KShift(self.fields.J[...,2],cos(kz_unmod*self.dz/2))

#        self.apply_KShift(self.fields.J[...,0],0.5*(1.+cos(kx_unmod*self.dx)))
#        self.apply_KShift(self.fields.J[...,1],exp(1j*ky_unmod*self.dz/2))
#        self.apply_KShift(self.fields.J[...,2],0.5*(1.+cos(kz_unmod*self.dz)))


#        self.apply_KShift(self.fields.J[...,0],0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))
#        self.apply_KShift(self.fields.J[...,1],0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))
#        self.apply_KShift(self.fields.J[...,2],0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))
#        self.apply_KShift(self.fields.Rho,0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))
#        self.apply_KShift(self.fields.DRhoodt,0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))

#        self.apply_KShift(self.fields.J[...,0],cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))
#        self.apply_KShift(self.fields.J[...,1],cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))
#        self.apply_KShift(self.fields.J[...,2],cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))
#        self.apply_KShift(self.fields.Rho,cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))
#        self.apply_KShift(self.fields.DRhoodt,cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))

        return


    def Jnode2yeeK(self):
        emK = self.FSpace

#        self.kx_unmod = kx_unmod
#        self.kz_unmod = kz_unmod

#        self.apply_KShift(self.fields.J[...,0],exp(1j*kx_unmod*self.dx/2))
#        self.apply_KShift(self.fields.J[...,1],exp(1j*ky_unmod*self.dz/2))
#        self.apply_KShift(self.fields.J[...,2],exp(1j*kz_unmod*self.dz/2))

        self.apply_KFilter(self.fields.J[...,0],cos(emK.kx_unmod*self.dx/2))
#        self.apply_KShift(self.fields.J[...,1],exp(1j*ky_unmod*self.dz/2))
        self.apply_KFilter(self.fields.J[...,2],cos(emK.kz_unmod*self.dz/2))

#        self.apply_KShift(self.fields.J[...,0],0.5*(1.+cos(kx_unmod*self.dx)))
#        self.apply_KShift(self.fields.J[...,1],exp(1j*ky_unmod*self.dz/2))
#        self.apply_KShift(self.fields.J[...,2],0.5*(1.+cos(kz_unmod*self.dz)))


#        self.apply_KShift(self.fields.J[...,0],0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))
#        self.apply_KShift(self.fields.J[...,1],0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))
#        self.apply_KShift(self.fields.J[...,2],0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))
#        self.apply_KShift(self.fields.Rho,0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))
#        self.apply_KShift(self.fields.DRhoodt,0.5*(1.+cos(kx_unmod*self.dx))*0.5*(1.+cos(kz_unmod*self.dz)))

#        self.apply_KShift(self.fields.J[...,0],cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))
#        self.apply_KShift(self.fields.J[...,1],cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))
#        self.apply_KShift(self.fields.J[...,2],cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))
#        self.apply_KShift(self.fields.Rho,cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))
#        self.apply_KShift(self.fields.DRhoodt,cos(kx_unmod*self.dx/2)*cos(kz_unmod*self.dz/2))

        return

    def push_spectral_psaotd(self):
    ###################################
    #              PSAOTD             #
    ###################################

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
