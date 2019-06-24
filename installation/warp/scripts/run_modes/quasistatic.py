"""
This script contains the class for performing runs in quasistatic approximation.
In parallel, there is a caveat that each processor carries particles at different time steps. 
"""

# TODO: Mesh Refinement; get Ez from beam using dA/dt=vz*dA/dz (for now, Ez=0).

from warp import *
from ..diagnostics.getzmom import *
#from ..field_solvers.AMR import *
from ..utils.appendablearray import *
from ..particles.ionization import *
import __main__
try:
  from pos import *
except:
  pass
try:
  import psyco
except ImportError:
  print 'Warning:  psyco not found'

class Quasistatic(SubcycledPoissonSolver):
  # this class differs from the previous one that it does not use the 3-D solver for the ions, but the 2-D one.
  def __init__(self,ions,MRroot=None,l_verbose=0,l_findmgparam=0,
               nparpgrp=top.nparpgrp,Ninit=1000,l_mode=1,l_selfe=1,l_selfi=1,maps=None,pboundxy=None,
               conductors=[],scraper=None,l_weakstrong=0,nelecperiod=1,
               l_elecuniform=0,xminelec=None,xmaxelec=None,yminelec=None,ymaxelec=None,nxelec=None,nyelec=None,
               backgroundtype=Electron,l_push_z=true,l_inject_elec_MR=false,l_warpzmmnt=false,
               npushzperiod=1,lattice=None,l_freeze_xelec=false,
               dispx=None,dispy=None,disppx=None,disppy=None,lcollapsemaps=1,
               xinitelec=None,yinitelec=None,winitelec=None,
               vxinitelec=None,vyinitelec=None,vzinitelec=None,strideinitelec=1,
               l_elec4sym=0,l_randinit=False,l_reuseelec=False,l_randomflip=False,
               l_beam_1d=False,sigmax=0.,sigmay=0.,l_planarYZ=0,
               l_elec_greensum=False,elec_aellipse=0.,elec_bellipse=0.,
               Ekick=0.,Lkick=0.,Akick=None,Bkick=1,
               nbuckets=1,tparallel=1,
               l_independentbunches=0,nsaverhophi=-1,
               l_periodicz=0,
               sec=None,l_posinst_track_electrons=False,
               l_gas_ionization = 0,          # turn on/off ionization of background gas
               gas_species      = Dihydrogen,  # background gas species
               gas_pressure     = 1.e-5*10,    # background gas pressure
               gas_temperature  = 300.,        # background gas temperature
               cross_section    = 2.e6,        # background gas ionization cross-section (by beam)
               emitted_energy0  = 100.,        # ionized electrons mean energy
               emitted_energy_sigma = 50.,     # ionized electrons energy sigma
               freeze_beam_list = [],
               nde=None,demax=0.1,l_zdediag=0,ke0=0.,
               **kw):
    assert top.grid_overlap != 2,"Error: the quasistatic class needs top.grid_overlap==1."
    assert (float(npes)/float(nbuckets))==float(npes/nbuckets),"Error in Quasistatic:the number of processors must be proportional to the number of buckets."
    self.nbuckets=nbuckets
    ntgroups=top.nxprocs*top.nyprocs
    self.ntgroups=ntgroups
    bucketmpiid=int(me*nbuckets/npes)
    self.bucketid=nbuckets-bucketmpiid
    if not lparallel:# or nbuckets==1:
      if not lparallel:
        self.mympigroup=None
      else:
        self.mympigroup=comm_world
      self.gme=me/self.ntgroups
      self.gnpes=npes/self.ntgroups
      print 'gme = ',self.gme
      mygroupinit = arange(npes)
      mygroup = mygroupinit[me%self.ntgroups::self.ntgroups]
    else:
     if lparallel:
      procs = arange(npes)
      mygroupinit = compress(int(procs*nbuckets/npes)==bucketmpiid,procs)
      mygroup = mygroupinit[me%self.ntgroups::self.ntgroups]
      try:
        self.mympigroup=comm_world.comm_create(mygroup)
      except:
        ranks_group_list = list(mygroup)
        mpi_group = comm_world.Get_group()
        newgroup = mpi_group.Incl(ranks_group_list)
        self.mympigroup=comm_world.Create(newgroup)

      self.gme = self.mympigroup.rank
      self.gnpes = self.mympigroup.size
    if ntgroups>1:
      mytgroup = mygroupinit[ntgroups*((me-mygroupinit[0])/ntgroups)]+arange(ntgroups)
      print mygroupinit
      print "me,mygroup,mytgroup = ",me,mygroup,mytgroup
      try:
        self.mympitgroup=comm_world.comm_create(mytgroup)
      except:
        self.mympitgroup=comm_world.Create(mytgroup)
      self.gtme = self.mympitgroup.rank
      self.gtnpes = self.mympitgroup.size
    else:
      self.gtme = 0
      self.gtnpes=1
    length=(w3d.zmmax-w3d.zmmin)/self.nbuckets
    self.zminbucket=-length/2-(self.bucketid-1)*length
    self.zmaxbucket=self.zminbucket+length
    self.nzbucket=w3d.nz/self.nbuckets
    self.l_independentbunches=l_independentbunches
    self.grid_overlap = 1
    FieldSolver.__init__(self,kwdict=kw)
    if l_planarYZ:
      w3d.solvergeom=w3d.Ygeom
      frz.l_get_fields_on_grid=False
    else:
      w3d.solvergeom=w3d.XYgeom
    self.gridelecs=[]
    self.gridions=[]
    self.gridionscp=[]
    self.grids=[]
    self.tparallel = tparallel
    if l_selfi:
      frz.l_bgrid=true
    for gxy in [self.gridions,self.gridelecs,self.gridionscp]:
     for i in range(2):
      gxy.append(GRIDtype())
      self.grids.append(gxy[-1])
      frz.basegrid=gxy[-1]
      if l_planarYZ:
        frz.init_base(0,w3d.ny,w3d.xmmax-w3d.xmmin,w3d.dy,w3d.xmmin,w3d.ymmin,false)
      else:
        if ntgroups==1 or self.tparallel in [0,1]:
          if self.tparallel==0:
            assert top.nxprocs==top.nyprocs,"Error: nxprocs must be equal to nyprocs."
            frz.init_base(w3d.nx/top.nxprocs,w3d.ny/top.nyprocs,w3d.dx*top.nxprocs,w3d.dy*top.nyprocs,w3d.xmmin,w3d.ymmin,false)
            g=frz.basegrid
            ntransist=2*2
            for i in range(nint(log2(top.nxprocs))):
              add_patch(g.gid,
                        w3d.xmminlocal,
                        w3d.xmmaxlocal,
                        w3d.ymminlocal,
                        w3d.ymmaxlocal,
                        2,2,
                        ntransist,ntransist,ntransist,ntransist)
              g=g.down
          else:
            frz.init_base(w3d.nx,w3d.ny,w3d.dx,w3d.dy,w3d.xmmin,w3d.ymmin,false)
        else:
          frz.init_base(w3d.nxlocal,w3d.nylocal,w3d.dx,w3d.dy,w3d.xmminlocal,w3d.ymminlocal,true)
      frz.mgridrz_accuracy = f3d.mgtol
    self.conductors=conductors
    self.l_condinstalled=0
    self.MRroot = MRroot
    self.dt = w3d.dz/top.vbeam # time step for pushing electrons
    self.nparpgrp = nparpgrp
    self.l_verbose = l_verbose
    self.l_parallelverbose = 0
    self.l_findmgparam = l_findmgparam
    self.l_periodicz = l_periodicz
    self.Ninit = Ninit
    self.l_mode = l_mode
    if npes>1:
      self.l_mode=1
    else:
      self.l_mode = l_mode
    self.l_selfe = l_selfe
    self.l_selfi = l_selfi
    self.l_elecuniform=l_elecuniform
    if xminelec is None:
      self.xminelec=w3d.xmmin
    else:  
      self.xminelec=xminelec
    if xmaxelec is None:
      self.xmaxelec=w3d.xmmax
    else:  
      self.xmaxelec=xmaxelec
    if yminelec is None:
      self.yminelec=w3d.ymmin
    else:  
      self.yminelec=yminelec
    if ymaxelec is None:
      self.ymaxelec=w3d.ymmax
    else:  
      self.ymaxelec=ymaxelec
    self.nxelec=nxelec
    self.nyelec=nyelec
    self.l_elec4sym=l_elec4sym
    self.l_randinit=l_randinit
    self.l_reuseelec=l_reuseelec
    self.l_randomflip=l_randomflip
    self.l_weakstrong=l_weakstrong
    self.l_beam_1d=l_beam_1d
    self.sigmax=sigmax
    self.sigmay=sigmay
    self.l_elec_greensum=l_elec_greensum
    self.elec_aellipse=elec_aellipse
    self.elec_bellipse=elec_bellipse
    self.nelecperiod=nelecperiod
    self.l_inject_elec_MR=l_inject_elec_MR
    if self.l_inject_elec_MR:
      top.wpid=nextpid()
    self.scraper=scraper
    self.sec=sec
    self.l_posinst_track_electrons=l_posinst_track_electrons
    if maps is not None:
      self.l_maps=true
      self.maps=maps
    else:
      self.l_maps=false
    self.l_push_z=l_push_z and not l_beam_1d
    if pboundxy is None:
      self.pboundxy = top.pboundxy
    else:
      self.pboundxy = pboundxy
    self.l_planarYZ=l_planarYZ
    # --- sets z range
    if self.l_mode==2:
      self.izmin = w3d.nz/4
      self.izmax = w3d.nz-self.izmin
    else:
      self.izmin = 0
      self.izmax = w3d.nzp
    self.znodes = top.zpminlocal+arange(self.izmax+1)*w3d.dz
    self.wz0=self.wz1=None
    self.pgions = top.pgroup
    self.ions = [ions]
    top.lspeciesmoments=False
    for iz in range(globalmax(w3d.nzp)-1):
      self.ions.append(Species(type=ions.type,fselfb=top.pgroup.fselfb[0]))
      top.pgroup.sq[-1] = top.pgroup.sq[0]
    self.pgelec = ParticleGroup()
    self.pgelec.npid=top.npid
    self.pgelec.gchange()
    top.pgroup = self.pgelec
    self.electrons = Species(type=backgroundtype)
    self.backgroundtype=backgroundtype
    if sec is not None:
      for cond in conductors:
        sec.add(incident_species=self.electrons,
                emitted_species=self.electrons,
                conductor=cond,
                interaction_type = 0) 
    self.l_gas_ionization=l_gas_ionization
    if l_gas_ionization:
      self.ioniz=[]
      for iz in range(self.izmin,self.izmax):
        self.ioniz.append(Ionization(stride=1))#stride=10,nx=30,ny=30,xmax=0.06,ymax=0.06)
      torr_to_MKS=133.3224  # 1 Torr=133.3224  N/m**2
      # ideal gas law: rho=p/(k*T)
      gas_prefactor = (torr_to_MKS/boltzmann)/294.
      self.gas_density   = gas_prefactor*gas_pressure*(294./gas_temperature)
      self.cross_section = cross_section
      self.emitted_energy_sigma = emitted_energy_sigma
      self.emitted_energy0 = emitted_energy0
      self.ioniz_installed=0
    top.pgroup = self.pgions
    self.ionstoprev = [0]
    self.fullsortdone = 0
    self.pnum   = AppendableArray(typecode='d')
    self.xbar   = AppendableArray(typecode='d')
    self.ybar   = AppendableArray(typecode='d')
    self.zbar   = AppendableArray(typecode='d')
    self.xpbar  = AppendableArray(typecode='d')
    self.ypbar  = AppendableArray(typecode='d')
    self.xpnbar = AppendableArray(typecode='d')
    self.ypnbar = AppendableArray(typecode='d')
    self.x2     = AppendableArray(typecode='d')
    self.y2     = AppendableArray(typecode='d')
    self.z2     = AppendableArray(typecode='d')
    self.xp2    = AppendableArray(typecode='d')
    self.yp2    = AppendableArray(typecode='d')
    self.xxpbar = AppendableArray(typecode='d')
    self.yypbar = AppendableArray(typecode='d')
    self.xpn2   = AppendableArray(typecode='d')
    self.ypn2   = AppendableArray(typecode='d')
    self.xxpnbar   = AppendableArray(typecode='d')
    self.yypnbar   = AppendableArray(typecode='d')
    self.pnumztmp   = zeros(self.izmax+1,'d')
    self.xbarztmp   = zeros(self.izmax+1,'d')
    self.ybarztmp   = zeros(self.izmax+1,'d')
    self.zbarztmp   = zeros(self.izmax+1,'d')
    self.xpbarztmp  = zeros(self.izmax+1,'d')
    self.ypbarztmp  = zeros(self.izmax+1,'d')
    self.xpnbarztmp = zeros(self.izmax+1,'d')
    self.ypnbarztmp = zeros(self.izmax+1,'d')
    self.x2ztmp     = zeros(self.izmax+1,'d')
    self.y2ztmp     = zeros(self.izmax+1,'d')
    self.z2ztmp     = zeros(self.izmax+1,'d')
    self.xp2ztmp    = zeros(self.izmax+1,'d')
    self.yp2ztmp    = zeros(self.izmax+1,'d')
    self.xxpbarztmp = zeros(self.izmax+1,'d')
    self.yypbarztmp = zeros(self.izmax+1,'d')
    self.xpn2ztmp   = zeros(self.izmax+1,'d')
    self.ypn2ztmp   = zeros(self.izmax+1,'d')
    self.xxpnbarztmp   = zeros(self.izmax+1,'d')
    self.yypnbarztmp   = zeros(self.izmax+1,'d')
    self.pnumz   = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.xbarz   = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.ybarz   = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.zbarz   = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.xpbarz  = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.ypbarz  = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.xpnbarz = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.ypnbarz = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.x2z     = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.y2z     = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.z2z     = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.xp2z    = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.yp2z    = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.xxpbarz = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.yypbarz = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.xpn2z   = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.ypn2z   = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.xxpnbarz   = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.yypnbarz   = AppendableArray(typecode='d',unitshape=(self.izmax+1,))
    self.timemmnts = AppendableArray(typecode='d')
    self.l_zdediag=l_zdediag
    if self.l_zdediag:
      if nde is None:
        self.nde=self.nzbucket
      else:
        self.nde=nde
      self.ke0 = ke0
      self.demax=demax
      self.dketmp0 = zeros(self.nde+1,'d')
      self.denstmp0 = zeros(self.nde+1,'d')
      self.dketmp1 = zeros(self.nde+1,'d')
      self.denstmp1 = zeros(self.nde+1,'d')
      self.dketmp  = zeros((self.izmax+1,self.nde+1,),'d')
      self.denstmp = zeros((self.izmax+1,self.nde+1,),'d')
      self.zdedipole = AppendableArray(typecode='d',unitshape=(self.izmax+1,self.nde+1,))
      self.zdedens   = AppendableArray(typecode='d',unitshape=(self.izmax+1,self.nde+1,))
    self.iz=-1
    self.nelechist=AppendableArray()
    self.l_timing=false
    self.l_warpzmmnt=l_warpzmmnt
    self.reset_timers()
    self.l_freeze_xelec=l_freeze_xelec
    self.freeze_beam_list=freeze_beam_list
    self.npushzperiod=npushzperiod
    if lattice is not None and lcollapsemaps:
      collapsedlattice=[]
      trmap = None
      lmap=0.
      nmap=""
      for elmnt in lattice:
        if trmap is None:
          trmap=elmnt.Map.copy()
          lmap=elmnt.L
          nmap=elmnt.name
        else:
          trmap=dot(elmnt.Map,trmap)
          lmap+=elmnt.L
          nmap+='+'+elmnt.name
        if elmnt.ecflag or elmnt is lattice[-1]:
          elmnt.Map=trmap.copy()
          elmnt.L=lmap
          elmnt.name=nmap
          collapsedlattice.append(elmnt)
          trmap=None
          lmap=0.
          nmap=''
      lattice = collapsedlattice
    self.lattice=lattice
    if self.lattice is not None:
      self.ist=0
      self.nst=len(lattice)
      top.dt=lattice[0].L/top.vbeam
    self.ilcount=zeros(w3d.nzp+1)
    self.ionstonext = [0]
    self.xinitelec=xinitelec
    self.yinitelec=yinitelec
    self.vxinitelec=vxinitelec
    self.vyinitelec=vyinitelec
    self.vzinitelec=vzinitelec
    self.winitelec=winitelec
    self.strideinitelec=strideinitelec
    self.ielecstride=0
    self.fbgainy=0.
    self.Ekick=Ekick
    self.Lkick=Lkick
    self.Akick=Akick
    self.Bkick=Bkick
    self.mystation=zeros(self.izmax,'d')
    self.nsaverhophi=nsaverhophi
    self.timehist=AppendableArray()
    self.smoother_ions=None
    self.smoother_elec=None
    self.store_rhotoprev()
    self.store_rhotonext()
    if self.l_posinst_track_electrons:
         top.pgroup=self.pgelec
         pp_pos2warp(0)
         top.pgroup=self.pgions
    
  def reset_timers(self):
    self.time_loop=0.
    self.time_sort=0.
    self.time_getmmnts=0.
    self.time_compute_sw=0.
    self.time_sendrecv_storeions=0.
    self.time_set_gamma=0.
    self.time_deposit_ions=0.
    self.time_solve=0.
    self.time_deposit_electrons=0.
    self.time_add_ei_fields=0.
    self.time_stack_rhophi=0.
    self.time_push_electrons=0.
    self.time_push_ions=0.
    self.time_reset_rho1=0.
    self.time_store_ionstoprev=0.

  def install_conductors(self):
   for g in self.gridions+self.gridelecs+self.gridionscp:
      for conductor in self.conductors:
          try:
            cond = conductor.cond
          except AttributeError:
            cond = conductor
          for ig in range(frz.ngrids):
#          installconductors(cond,conductors=ConductorType(),gridrz=g,nz=0)
            installconductors(cond,conductors=ConductorType(),
                              nx=g.nr,ny=g.nz,nz=0,
                              xmmin=g.rmin,xmmax=g.rmax,
                              ymmin=g.zmin,ymmax=g.zmax,
                              gridrz=g)
            if ig<frz.ngrids-1:g=g.down
   self.l_condinstalled=1
  
  def install_ionization(self):
      for iz in range(self.izmin,self.izmax):
        self.ioniz[iz].add(incident_species = self.ions[iz],
                     ndens            = self.gas_density,
                     emitted_species  = [self.electrons],
#                    emitted_species  = [gas_ions,electrons],
                     incident_pgroup  = self.pgions,
                     emitted_pgroup   = self.pgelec,
                     cross_section    = self.cross_section*1.e-28*self.pgions.sw[0]/(self.pgelec.sw[0]*w3d.dz),
                     emitted_energy_sigma = self.emitted_energy_sigma,
                     emitted_energy0 = self.emitted_energy0)
      self.ioniz_installed=1
      
  def step(self,nt=1,l_return_dist=false,l_plotelec=0,l_savedist=0,l_printsteps=1):
   # --- if running in parallel, enforcing l_mode=1
   if npes>1:self.l_mode=1
   if self.nsaverhophi>0:self.setuprhophi()
   if not self.l_condinstalled:self.install_conductors()
   if self.l_gas_ionization and not self.ioniz_installed:self.install_ionization()
   for it in range(nt):
     ptimeloop = wtime()
     
     # --- sets flag for pushing electrons
     if self.lattice is not None:
       l_push_elec = self.lattice[self.ist].ecflag
       l_noelec = not l_push_elec
       if self.l_verbose:print 'l_push_elec',l_push_elec
     else:
       l_push_elec = (not self.l_weakstrong) or ((top.it-(npes-me))%self.nelecperiod==0)

     # --- sets flag for pushing beam particles
     if self.bucketid in self.freeze_beam_list:
       self.l_push_ions = False
     else:
       self.l_push_ions = True

     # --- call beforestep functions
     callbeforestepfuncs.callfuncsinlist()

     # --- sort ions along z
     if self.l_timing:ptime = wtime()
     self.sort_ions_along_z()
     if self.l_timing: self.time_sort += wtime()-ptime

     # --- gather moments
#     if top.it==0 or (me>=(npes-top.it) and (top.it-(npes-me))%top.nhist==0):
#     if top.it==0 or (me>=(npes-top.it) and (top.it-(self.gnpes-self.gme))%top.nhist==0):
#     self.itbucket = top.it-(self.bucketid-1)*npes/self.nbuckets
     self.itbucket = top.it-(self.bucketid-1)*npes/self.nbuckets
     if top.it==0 or (me>=(npes-top.it*self.ntgroups) and (self.itbucket-(self.gnpes-self.gme))%top.nhist==0):
       if self.l_timing:ptime = wtime()
       if self.l_timing: self.time_getmmnts += wtime()-ptime

     # --- compute scatter/gather weights for beam ions (linear weighting)
     if self.l_timing:ptime = wtime()
     self.compute_sw()
     if self.l_timing: self.time_compute_sw += wtime()-ptime

     # --- save distribution in file
     if l_savedist:
       if self.nbuckets==1:
         f=PW.PW('dist%g.pdb'%top.it)
       else:
         f=PW.PW('dist_%g_%g.pdb'%(self.bucketid,top.it))
       for iz in range(w3d.nzp):
         exec('f.x%i=getx(js=iz)'%iz)
         exec('f.y%i=gety(js=iz)'%iz)
         exec('f.z%i=getz(js=iz)'%iz)
         exec('f.vx%i=getvx(js=iz)'%iz)
         exec('f.vy%i=getvy(js=iz)'%iz)
         exec('f.vz%i=getvz(js=iz)'%iz)
         exec('f.ex%i=getex(js=iz)'%iz)
         exec('f.ey%i=getey(js=iz)'%iz)
         exec('f.ez%i=getez(js=iz)'%iz)
       f.close()

     if l_push_elec:
       # --- generate electrons (on last processor only)
       if self.l_independentbunches:
          if self.gme>=(self.gnpes-self.ntgroups) and not(l_plotelec and top.it>0):
            self.pgelec.nps[0]=0
            self.create_electrons()
       else:
          if me>=(npes-self.ntgroups) and not(l_plotelec and top.it>0):self.create_electrons()
#       if me==max(0,npes-1):self.create_electrons()

       # --- (diagnostic) stacking of electron distribution
       if l_return_dist:
         self.dist=[]
         sp=self.electrons
         self.dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])

     # --- must be called right before entering the main loop
     # --- stores ions to be passed to previous processor in temporary arrays
#     iz = 0
#     if self.l_timing:ptime = wtime()
#     self.store_ionstoprev()
#     if self.l_timing: self.time_store_ionstoprev += wtime()-ptime
#     self.apply_ions_bndconditions(iz)
     self.elec=[]
     # --- loop over "3-D" grid mesh in z
     for iz in range(self.izmax-1,self.izmin-1,-1):
       self.iz=iz
         
       if self.iz==(self.izmax-3):
         if self.l_timing:ptime = wtime()
         self.sendrecv_storedions_toprev()       
         if self.l_timing: self.time_sendrecv_storeions += wtime()-ptime
       if self.l_timing:ptime = wtime()

       # --- push ions velocity (2nd half)
       if self.l_push_ions and me>=(npes-(top.it+1)*self.ntgroups):
         if iz<w3d.nzp-1:
           self.set_gamma(iz+1)
#           self.push_ions_velocity_second_half(iz+1)
         if iz==0:
           self.set_gamma(iz)
#           self.push_ions_velocity_second_half(iz)
       if self.l_timing: self.time_set_gamma += wtime()-ptime
 
       if not self.l_beam_1d:
         if self.l_timing:ptime = wtime()
         # --- deposit ions charge
         self.deposit_ions()
         if self.l_timing: self.time_deposit_ions += wtime()-ptime

       # --- call 2-D field solver for ions
       if (self.l_selfi or l_push_elec) and not self.l_beam_1d:
         if self.l_timing:ptime = wtime()
         frz.basegrid=self.gridions[1];mk_grids_ptr()
         # --- optimize mgparam
#         if  self.l_findmgparam and top.it-1==npes-me-1 and self.iz in [w3d.nzlocal/2,w3d.nzlocal/2+1]:
         if  self.l_findmgparam and top.it==0 and self.iz in [w3d.nzlocal/2,w3d.nzlocal/2+1]:
          if me>=npes-self.ntgroups:find_mgparam_rz(true)
         self.gridions[1].phi=0.
         if self.ntgroups>1:self.sum_rho_in_tgroups(self.gridions[1])
         if self.smoother_ions is not None:self.smooth_rho(self.gridions[1],self.smoother_ions)
         solve_mgridrz(self.gridions[1],frz.mgridrz_accuracy,true)
         if self.iz==0:
           frz.basegrid=self.gridions[0];mk_grids_ptr()
           self.gridions[0].phi=0.
           if self.ntgroups>1:self.sum_rho_in_tgroups(self.gridions[0])
           if self.smoother_ions is not None:self.smooth_rho(self.gridions[0],self.smoother_ions)
           solve_mgridrz(self.gridions[0],frz.mgridrz_accuracy,true)
         if self.l_timing: self.time_solve += wtime()-ptime

       if l_push_elec and not self.l_elec_greensum:
         if self.l_timing:ptime = wtime()
         # --- deposit electrons charge
         self.deposit_electrons(1)
         if self.l_timing: self.time_deposit_electrons += wtime()-ptime

       if not self.l_elec_greensum:
         if self.l_timing:ptime = wtime()
         # --- call 2-D field solver for electrons
         frz.basegrid=self.gridelecs[1];mk_grids_ptr()
           # --- optimize mgparam
#           if  self.l_findmgparam and top.it-1==npes-me-1 and self.iz in [w3d.nzlocal/2,w3d.nzlocal/2+1]:
         if l_push_elec:
           if  self.l_findmgparam and top.it==0 and self.iz in [w3d.nzlocal/2,w3d.nzlocal/2+1]:
             if me==npes-1:find_mgparam_rz(true)
           self.gridelecs[1].phi=0.
           if self.ntgroups>1:self.sum_rho_in_tgroups(self.gridelecs[1])
           if self.smoother_elec is not None:self.smooth_rho(self.gridelecs[1],self.smoother_elec)
           solve_mgridrz(self.gridelecs[1],frz.mgridrz_accuracy,true)
         else:
           if l_noelec:
             ge = self.gridelecs[1]
             for ig in range(frz.ngrids):
               if ig>0:
                 ge = ge.down
               ge.phi[...]=0.
           else:
             ge = self.gridelecs[1]
             for ig in range(frz.ngrids):
               if ig>0:
                 ge = ge.down
               ge.phi[...]=self.phie[ig][:,:,self.iz+1]
         if self.l_timing: self.time_solve += wtime()-ptime

       if not self.l_elec_greensum:
         if self.l_timing:ptime = wtime()
         # --- add electron and ion fields
         self.add_ei_fields(1)
         if self.l_timing: self.time_add_ei_fields += wtime()-ptime

       if self.l_timing:ptime = wtime()

       # --- gather rho and phi copies0
       if self.nsaverhophi>0:
#         if (me>=(npes-(top.it+1)*self.ntgroups) and (self.itbucket-(self.gnpes-self.gme))%self.nsaverhophi==0):
         if (me>=(npes-(top.it+1)*self.ntgroups) and (self.itbucket-(self.gnpes-self.gme))%self.nsaverhophi==0):
           if l_push_elec and not self.l_elec_greensum:
             ge = self.gridelecs[1]
             for ig in range(frz.ngrids):
               if ig>0:
                 ge = ge.down
               self.rhoe[ig][:,:,self.iz+1] = ge.rho.copy()
               self.phie[ig][:,:,self.iz+1] = ge.phi.copy()
           if not self.l_beam_1d:
             gi = self.gridions[1]
             for ig in range(frz.ngrids):
               if ig>0:
                 gi = gi.down
               self.rhoi[ig][:,:,self.iz+1] = gi.rho.copy()
               self.phii[ig][:,:,self.iz+1] = gi.phi.copy()
       if self.l_timing: self.time_stack_rhophi += wtime()-ptime

       if l_push_elec:
         if self.l_timing:ptime = wtime()
         # --- push electrons 
         self.push_electrons()
         if iz==0 and not self.l_elec_greensum:
           self.deposit_electrons(0)
           frz.basegrid=self.gridelecs[0];mk_grids_ptr()
           self.gridelecs[0].phi=0.
           if self.ntgroups>1:self.sum_rho_in_tgroups(self.gridelecs[0])
           if self.smoother_elec is not None:self.smooth_rho(self.gridelecs[0],self.smoother_elec)
           solve_mgridrz(self.gridelecs[0],frz.mgridrz_accuracy,true)
           self.add_ei_fields(0)
         if self.l_timing: self.time_push_electrons += wtime()-ptime

         # --- plot electrons
#         if l_plotelec :self.plot_electrons()
         if l_plotelec and (npes==1 or iz<w3d.nz/max(1,npes)-1):self.plot_electrons()
#         if l_plotelec and iz<w3d.nz/max(1,npes):self.plot_electrons()

       # --- gather moments
#       if iz<w3d.nzp-1 and (top.it==0 or (me>=(npes-top.it) and (top.it-(npes-me))%top.nhist==0)):
#       if iz<w3d.nzp-1 and (top.it==0 or (me>=(npes-top.it) and (top.it-(self.gnpes-self.gme))%top.nhist==0)):
       if iz<w3d.nzp-1 and (top.it==0 or (me>=(npes-top.it*self.ntgroups) and (self.itbucket-(self.gnpes-self.gme-1))%top.nhist==0)):

         if self.l_timing:ptime = wtime()
         self.getmmnts(iz+1)
         if iz==0:self.getmmnts(iz)
         if self.l_timing: self.time_getmmnts += wtime()-ptime

       if self.l_timing:ptime = wtime()
       # --- push ions
       # WARNING: under current configuration, velocity push MUST be BEFORE positions push
       if self.l_push_ions:
         if me>=(npes-(top.it+1)*self.ntgroups):
           self.gather_ions_fields()
           if iz<self.izmax-1:
             self.push_ions_velocity_full(iz+1)
             self.push_ions_positions(iz+1)
             if iz==self.izmax-2:self.store_ionstonext(iz+1)
         if iz<self.izmax-1:
           self.apply_ions_bndconditions(iz+1)
         if me>=(npes-(top.it+1)*self.ntgroups):
           if iz==0:
             self.push_ions_velocity_full(iz)
             self.push_ions_positions(iz)
             if self.l_timing:ptime = wtime()
             if self.l_timing: self.time_store_ionstoprev += wtime()-ptime
             self.store_ionstoprev()
         if iz==0:
             self.apply_ions_bndconditions(iz)
         if me>=(npes-(top.it+1)*self.ntgroups):
           if iz==0:
             if self.lattice is not None:
               self.ist = (self.ist+1)%self.nst

       if self.l_timing: self.time_push_ions += wtime()-ptime

       if l_push_elec:
         # --- (diagnostic) stacking of electron distribution
         if l_return_dist:
           sp=self.slist[0]
           self.dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])

       # --- switch 2-D solver grid
       for g in [self.gridelecs,self.gridions]:
         gtmp=g.pop(0)
         g.append(gtmp)
#       for iz in range(globalmax(w3d.nzp)-w3d.nzp):
#         if l_plotelec :self.plot_electrons()
         
    # --- END of FOR loop
     if self.l_timing:ptime = wtime()

    # --- gather latest rho and phi copies and save into file
     if self.nsaverhophi>0:
       if me==0 and top.it==0:
           f=PW.PW('rhophidata.pdb')
           f.ngrids=frz.ngrids
           f.npes=npes
           f.nz=w3d.nz
           f.zmin=w3d.zmmin
           f.zmax=w3d.zmmax
           nx=[]
           ny=[]
           xmin=[]
           xmax=[]
           ymin=[]
           ymax=[]
           ge = self.gridelecs[1]
           for ig in range(frz.ngrids):
             if ig>0:
               ge = ge.down
             nx.append(ge.nr)
             ny.append(ge.nz)
             xmin.append(ge.rmin)
             xmax.append(ge.rmax)
             ymin.append(ge.zmin)
             ymax.append(ge.zmax)
           f.nx=nx
           f.ny=ny
           f.xmin=xmin
           f.xmax=xmax
           f.ymin=ymin
           f.ymax=ymax
           f.close()
       if (me>=(npes-(top.it+1)*self.ntgroups) and (self.itbucket-(self.gnpes-self.gme))%self.nsaverhophi==0):
         if l_push_elec and not self.l_elec_greensum:
           ge = self.gridelecs[1]
           for ig in range(frz.ngrids):
             if ig>0:
               ge = ge.down
             self.rhoe[ig][:,:,0] = ge.rho.copy()
             self.phie[ig][:,:,0] = ge.phi.copy()
         if not self.l_beam_1d:
           gi = self.gridions[1]
           for ig in range(frz.ngrids):
             if ig>0:
               gi = gi.down
             self.rhoi[ig][:,:,0] = gi.rho.copy()
             self.phii[ig][:,:,0] = gi.phi.copy()
         f=PW.PW('rhophi%g_%g.pdb'%(me,self.itbucket-(self.gnpes-1-self.gme)))
         f.rhoe=self.rhoe
         f.phie=self.phie
         f.rhoi=self.rhoi
         f.phii=self.phii
         f.close()
     if self.l_timing: self.time_stack_rhophi += wtime()-ptime

     if self.l_timing:ptime = wtime()
     self.reset_rho1()
     if self.l_timing: self.time_reset_rho1 += wtime()-ptime
#       if lparallel:comm_world.barrier()
       
     # --- clear electrons on processor 0, shift them to the previous ones for me>0
     self.clear_electrons()
     top.pgroup = self.pgions
     
     # --- broadcast of mgparams from last processor
     if top.it==0:self.sendrecv_mgparams()

     # --- store moments data
#     if top.it==0 or (me>=(npes-top.it) and (top.it-(npes-me))%top.nhist==0):
#     if top.it==0 or (me>=(npes-top.it) and (top.it-(self.gnpes-self.gme))%top.nhist==0):
     if top.it==0 or (me>=(npes-top.it*self.ntgroups) and (self.itbucket-(self.gnpes-self.gme-1))%top.nhist==0):

       if self.l_timing:ptime = wtime()
       self.getmmnts_store()
       if self.l_timing: self.time_getmmnts += wtime()-ptime

     # --- update time, time counter
     if me>=(npes-(top.it+1)*self.ntgroups):top.time+=top.dt
     if self.lattice is not None:
       top.dt=self.lattice[self.ist].L/top.vbeam
     if self.l_verbose: print me,top.it,self.iz,'me = ',me,';compute zmmnt'
     if self.l_warpzmmnt and top.it%top.nhist==0:
       zmmnt()
       minidiag(top.it,top.time,top.lspecial)
     top.it+=1
#     self.itbucket+=1

     if self.l_verbose: print me,top.it,self.iz,'me = ',me,';call afterstep functions'
     # --- call afterstep functions
     callafterstepfuncs.callfuncsinlist()
     if self.l_verbose:print me,self.iz,'it = %i, time = %gs.'%(top.it,top.time)

     if l_printsteps and me==0:print me,self.iz,'it = %i, time = %gs.'%(top.it,top.time)

     self.time_loop = wtime()-ptimeloop
#     print 'time loop = ',self.time_loop
#     self.timehist.append(self.time_loop)
#     fma();pla(self.timehist[:]);refresh()
     
  def print_timers(self):
    print 'loop               ',self.time_loop
    print 'sort               ',self.time_sort#,' = ',self.time_sort/time_loop,'\%'
    print 'getmmnts           ',self.time_getmmnts
    print 'compute_sw         ',self.time_compute_sw
    print 'sendrecv_storeions ',self.time_sendrecv_storeions
    print 'set_gamma          ',self.time_set_gamma
    print 'deposit_ions       ',self.time_deposit_ions
    print 'solve              ',self.time_solve
    print 'deposit_electrons  ',self.time_deposit_electrons
    print 'add_ei_fields      ',self.time_add_ei_fields
    print 'stack_rhophi       ',self.time_stack_rhophi
    print 'push_electrons     ',self.time_push_electrons
    print 'push_ions          ',self.time_push_ions
    print 'reset_rho1         ',self.time_reset_rho1
    print 'store_ionstoprev   ',self.time_store_ionstoprev
  
  def setuprhophi(self):
    try:self.rhoe
    except:
      self.rhoe = []
      self.phie = []
      self.rhoi = []
      self.phii = []
      ge = self.gridelecs[0]
      for ig in range(frz.ngrids):
        if ig>0:
          ge = ge.down
        self.rhoe.append(zeros([ge.nr+1,ge.nz+1,w3d.nzp+1],'d'))
        self.phie.append(zeros([ge.nr+2*ge.nguardx+1,ge.nz+2*ge.nguardz+1,w3d.nzp+1],'d'))
        self.rhoi.append(zeros([ge.nr+1,ge.nz+1,w3d.nzp+1],'d'))
        self.phii.append(zeros([ge.nr+2*ge.nguardx+1,ge.nz+2*ge.nguardz+1,w3d.nzp+1],'d'))
  
  def reset_rho1(self):
    if self.l_verbose:print me,top.it,self.iz,'enter reset_rho1'
    # --- reset rho1 on proc 0
    if me == npes-1:
     for i in range(2):
      bg=frz.basegrid=self.gridions[i]
      mk_grids_ptr()
      reset_rzmgrid_rho()
    if self.l_verbose:print me,top.it,self.iz,'exit reset_rho1'

  def store_rhotonext(self):
    gi = self.gridionscp[0]
    self.rhotonext=[]
    for ig in range(frz.ngrids):
      if ig>0:
        gi = gi.down
      self.rhotonext.append(gi.rho.copy())
  
  def store_rhotoprev(self):
    gi = self.gridions[0]
    self.rhotoprev=[]
    for ig in range(frz.ngrids):
      if ig>0:
        gi = gi.down
      self.rhotoprev.append(gi.rho.copy())
  
  def send_rhotonext(self):
      if me<npes-self.ntgroups:
        mpisend(self.rhotonext, dest = me+self.ntgroups)
      if me>self.ntgroups-1:
        recved = mpirecv(source = me-self.ntgroups)
        g = self.gridions[0]
        for ig in range(frz.ngrids):
          if ig>0:
            g = g.down
          g.rho[...] += recved[ig]
  
  def send_rhotoprev(self):
      if me>self.ntgroups-1:
        mpisend(self.rhotoprev, dest = me-self.ntgroups)
      if me<npes-self.ntgroups:
        recved = mpirecv(source = me+self.ntgroups)
        g = self.gridions[1]
        for ig in range(frz.ngrids):
          if ig>0:
            g = g.down
          g.rho[...] += recved[ig]
    
  def store_ionstoprev(self):
    if self.l_verbose:print me,top.it,self.iz,'enter store_ionstoprev'
    js = 0
    pg = self.pgions
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    if pg.nps[js]==0:
      self.ionstoprev = [0]
      return    
    ii = il+compress(pg.zp[il:iu]<w3d.zmminp,arange(pg.nps[js]))   
    if len(ii)==0:
      self.ionstoprev = [0]
    else:
#      print 'me',me,len(ii),(take(pg.zp,ii,0)-w3d.zmminp)/w3d.dz
      self.ionstoprev = [len(ii)]
      self.ionstoprev.append(take(pg.xp,ii,0).copy())
      self.ionstoprev.append(take(pg.yp,ii,0).copy())
      self.ionstoprev.append(take(pg.zp,ii,0).copy())
      self.ionstoprev.append(take(pg.uxp,ii,0).copy())
      self.ionstoprev.append(take(pg.uyp,ii,0).copy())
      self.ionstoprev.append(take(pg.uzp,ii,0).copy())
      self.ionstoprev.append(take(pg.gaminv,ii,0).copy())
      self.ionstoprev.append(take(pg.ex,ii,0).copy())
      self.ionstoprev.append(take(pg.ey,ii,0).copy())
      self.ionstoprev.append(take(pg.ez,ii,0).copy())
      self.ionstoprev.append(take(pg.bx,ii,0).copy())
      self.ionstoprev.append(take(pg.by,ii,0).copy())
      self.ionstoprev.append(take(pg.bz,ii,0).copy())
      if pg.npid>0:
        self.ionstoprev.append(take(pg.pid,ii,0).copy())
    put(pg.gaminv,ii,0.)
    processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
    self.set_sw(js)
    if self.l_verbose:print me,top.it,self.iz,'exit store_ionstoprev'
    
  def store_ionstonext(self,js):
    if self.l_verbose:print me,top.it,js,'enter store_ionstonext'
    pg = self.pgions
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    if pg.nps[js]==0:
      self.ionstonext = [0]
      return    
    ii = il+compress(pg.zp[il:iu]>=w3d.zmmaxp,arange(pg.nps[js]))    
    if len(ii)==0:
      self.ionstonext = [0]
    else:
#      print 'me',me,len(ii),(take(pg.zp,ii,0)-w3d.zmminp)/w3d.dz
      self.ionstonext = [len(ii)]
      self.ionstonext.append(take(pg.xp,ii,0).copy())
      self.ionstonext.append(take(pg.yp,ii,0).copy())
      self.ionstonext.append(take(pg.zp,ii,0).copy())
      self.ionstonext.append(take(pg.uxp,ii,0).copy())
      self.ionstonext.append(take(pg.uyp,ii,0).copy())
      self.ionstonext.append(take(pg.uzp,ii,0).copy())
      self.ionstonext.append(take(pg.gaminv,ii,0).copy())
      self.ionstonext.append(take(pg.ex,ii,0).copy())
      self.ionstonext.append(take(pg.ey,ii,0).copy())
      self.ionstonext.append(take(pg.ez,ii,0).copy())
      self.ionstonext.append(take(pg.bx,ii,0).copy())
      self.ionstonext.append(take(pg.by,ii,0).copy())
      self.ionstonext.append(take(pg.bz,ii,0).copy())
      if pg.npid>0:
        self.ionstonext.append(take(pg.pid,ii,0).copy())
    put(pg.gaminv,ii,0.)
    processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
    self.set_sw(js)
    if self.l_verbose:print me,top.it,self.iz,'exit store_ionstonext'
    
  def sendrecv_mgparams(self):
    if self.l_verbose:print me,'enter sendrecv_mgparams'
    if not lparallel:return
    tosend = []
    if me==npes-1:
      for glist in [self.gridions,self.gridelecs]:
        for ig in range(len(glist)):
          g = glist[ig]
          for i in range(frz.ngrids):
            if i>0:g=g.down
            tosend.append(g.npre)
            tosend.append(g.npost)
            tosend.append(g.npmin)
            tosend.append(g.mgparam)
#      for ip in range(npes-1):
#        comm_world.send(tosend,ip)
      if self.l_verbose:print me,tosend
    recved = warp_parallel.broadcast(tosend,npes-1)
    if self.l_verbose:print me,recved
    if me<npes-1:
#     recved,status = comm_world.recv(npes-1)
      i = 0
      for glist in [self.gridions,self.gridelecs]:
        for ig in range(len(glist)):
          g = glist[ig]
          for ig in range(frz.ngrids):
            if ig>0:g=g.down
            g.npre = recved[i];i+=1
            g.npost = recved[i];i+=1
            g.npmin = recved[i];i+=1
            g.mgparam = recved[i];i+=1
    if self.l_verbose:print me,'exit sendrecv_mgparams'

  def sendrecv_storedions_toprev(self):
    if not lparallel:return
    if self.l_verbose:print me,top.it,self.iz,'enter sendrecv_storedions_toprev'
    # --- sends stored ions
    if me>self.ntgroups-1:
      mpisend(self.ionstoprev, dest = me-self.ntgroups)
      if self.ionstoprev[0]>0:
        if self.l_parallelverbose:print me, 'sends ',self.ionstoprev[0],' ions to ',me-1
#        print 'send itn',self.ionstoprev[0],self.ionstoprev[3]
    # --- receives stored ions
    if me<npes-self.ntgroups:
      js = w3d.nzp-1
      pg = self.pgions
      recved = mpirecv(source = me+self.ntgroups)
      np = recved[0]
      if np>0:
        if self.l_parallelverbose:print me, 'recvs ',np,' ions from ',me+1
#        print 'recved itn',recved[0],recved[3]
        if pg.npid>0:
          pid=recved[-1]
        else:
          pid=0.
        addparticles(js = js,
                     x  = recved[1],
                     y  = recved[2],
                     z  = recved[3],
                     vx = recved[4],
                     vy = recved[5],
                     vz = recved[6],
                     gi = recved[7],
                     ex = recved[8],
                     ey = recved[9],
                     ez = recved[10],
                     bx = recved[11],
                     by = recved[12],
                     bz = recved[13],                     
                     pid=pid,
                     lmomentum=1,
                     lallindomain=1,
                     lfields=1,
                     pgroup=self.pgions)
        self.set_sw(js)
    if self.l_verbose:print me,top.it,self.iz,'exit sendrecv_storedions_toprev'
             
  def sendrecv_storedions_tonext(self):
    if not lparallel:return
    if self.l_verbose:print me,top.it,self.iz,'enter sendrecv_storedions_tonext'
    # --- sends stored ions
    if me<npes-self.ntgroups:
      mpisend(self.ionstonext, dest = me+self.ntgroups)
      if self.ionstonext[0]>0:
        if self.l_parallelverbose:print me, 'sends ',self.ionstonext[0],' ions to ',me+1
#        print 'send',self.ionstonext[0],self.ionstonext[3]
    # --- receives stored ions
    if me>self.ntgroups-1:
      js = 0
      pg = self.pgions
      recved = mpirecv(source = me-self.ntgroups)
      np = recved[0]
      if np>0:
        if self.l_parallelverbose:print me, 'recvs ',np,' ions from ',me-1
#        print 'recved',recved[0],recved[3]
        if pg.npid>0:
          pid=recved[-1]
        else:
          pid=0.
        addparticles(js = js,
                     x  = recved[1],
                     y  = recved[2],
                     z  = recved[3],
                     vx = recved[4],
                     vy = recved[5],
                     vz = recved[6],
                     gi = recved[7],
                     ex = recved[8],
                     ey = recved[9],
                     ez = recved[10],
                     bx = recved[11],
                     by = recved[12],
                     bz = recved[13],                     
                     pid=pid,
                     lmomentum=1,
                     lallindomain=1,
                     lfields=1,
                     pgroup=self.pgions)
        self.set_sw(js)
    if self.l_verbose:print me,top.it,self.iz,'exit sendrecv_storedions_tonext'
             
  def sendparticlestonextold(self,js,ii):
    if self.l_verbose:print me,top.it,self.iz,'enter sendparticlestonext'
    pg = self.pgions
    top.pgroup = pg
    tosend = []
    tosend.append(len(ii))
    if len(ii)>0:
      tosend.append(take(pg.xp,ii))
      tosend.append(take(pg.yp,ii))
      tosend.append(take(pg.zp,ii))
      tosend.append(take(pg.uxp,ii))
      tosend.append(take(pg.uyp,ii))
      tosend.append(take(pg.uzp,ii))
      tosend.append(take(pg.gaminv,ii))
      tosend.append(take(pg.ex,ii))
      tosend.append(take(pg.ey,ii))
      tosend.append(take(pg.ez,ii))
      tosend.append(take(pg.bx,ii))
      tosend.append(take(pg.by,ii))
      tosend.append(take(pg.bz,ii))
      if top.npid>0:tosend.append(take(pg.pid,ii,0))
    mpisend(tosend, dest = me+1)
    if tosend[0]>0:
      if self.l_parallelverbose:print me, 'sends ',tosend[0],' ions to ',me+1
    if len(ii)>0:
      put(pg.gaminv,ii,0.)
    processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
    self.set_sw(js)
    if self.l_verbose:print me,top.it,self.iz,'exit sendparticlestonext'
    
  def recvparticlesfrompreviousold(self):
    if self.l_verbose:print me,top.it,self.iz,'enter recvparticlesfromprevious'
    pg = self.pgions
    top.pgroup = pg
    recved = mpirecv(source = me-1)
    np = recved[0]
    if np>0:
      if self.l_parallelverbose:print me, 'recvs ',np,' ions from ',me-1
      if top.npid>0:
        pid=self.recved[-1]
      else:
        pid=0.
      addparticles(js = 0,
                   x  = recved[1],
                   y  = recved[2],
                   z  = recved[3],
                   vx = recved[4],
                   vy = recved[5],
                   vz = recved[6],
                   gi = recved[7],
                   ex = recved[8],
                   ey = recved[9],
                   ez = recved[10],
                   bx = recved[11],
                   by = recved[12],
                   bz = recved[13],                     
                   pid=pid,
                   lmomentum=1,
                   lallindomain=1,
                   lfields=1,
                   pgroup=self.pgions)
      self.set_sw(0)
    if self.l_verbose:print me,top.it,self.iz,'exit recvparticlesfromprevious'
             
  def deposit_ions(self):
    if self.l_verbose:print me,top.it,self.iz,'enter deposit_ions'
    if self.iz==self.izmax-3:self.sendrecv_storedions_tonext()
    if self.iz==0:self.deposit_ions_last_step()
    if self.iz==0:self.store_rhotonext()
    if self.iz==self.izmax-1:
      bg=frz.basegrid=self.gridions[1]
      mk_grids_ptr()
      reset_rzmgrid_rho()
    pg = self.pgions
    top.pgroup = pg
    js = self.iz
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    bg=frz.basegrid=self.gridions[0]
    mk_grids_ptr()
    # --- reset 2-D rho arrays
    reset_rzmgrid_rho()
    if self.iz==0:self.send_rhotonext()
    if me>=(npes-(top.it+1)*self.ntgroups):
    # --- send rho to next proc on last step
     if iu-il>0:
      xp = pg.xp[il:iu]
      yp = pg.yp[il:iu]
      q = pg.sq[js]*pg.sw[js]/w3d.dz
      # --- deposit ion charge in first grid
      if self.l_planarYZ:
        rhoweightz_weight(yp,self.wz0[js],iu-il,q,bg.nz,bg.dz,0.)
      else:
        rhoweightrz_weights(xp,yp,yp,self.wz0[js],iu-il,q,bg.nr,bg.nz,bg.dr,bg.dz,bg.rmin,0.)
     if self.iz<self.izmax-1 or me==npes-1:
      bg=frz.basegrid=self.gridions[1]
      mk_grids_ptr()
      if iu-il>0:
        # --- deposit ion charge in second grid
        if self.l_planarYZ:
          rhoweightz_weight(yp,self.wz1[js],iu-il,q,bg.nz,bg.dz,0.)
        else:
          rhoweightrz_weights(xp,yp,yp,self.wz1[js],iu-il,q,bg.nr,bg.nz,bg.dr,bg.dz,bg.rmin,0.)
    if self.iz==0:           self.store_rhotoprev()
    if self.iz==self.izmax-1:self.send_rhotoprev()
    # --- distribute rho among patches (WARNING: must be called only once after everything has been deposited.)
    bg=frz.basegrid=self.gridions[1]
    mk_grids_ptr()
    # --- exchange data between procs
    if self.ntgroups>1:get_rho_from_rhop(frz.basegrid)
    distribute_rho_rz()
    if self.iz==0:
      bg=frz.basegrid=self.gridions[0]
      mk_grids_ptr()
      # --- exchange data between procs
      if self.ntgroups>1:get_rho_from_rhop(frz.basegrid)
      distribute_rho_rz()
    if self.l_verbose:print me,top.it,self.iz,'exit deposit_ions'
                      
  def deposit_ions_last_step(self):
  # deposits ions on last cell at last step, using ions newly advanced, from last cell.
    if self.l_verbose:print me,top.it,self.iz,'enter deposit_ions_last_step'
    pg = self.pgions
    top.pgroup = pg
    iz = w3d.nzp-1
    js = iz
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    np=iu-il
    # --- reset 2-D rho arrays
    bg=frz.basegrid=self.gridionscp[0]
    mk_grids_ptr()
    reset_rzmgrid_rho()
    if np>0:
      q = pg.sq[js]*pg.sw[js]/w3d.dz
      zmin=w3d.zmminp+js*w3d.dz
      wz1 = (pg.zp[il:iu]-zmin)/w3d.dz
      # --- select particles in slice js
      ii = compress(floor(wz1)==0,arange(np))
      xp = take(pg.xp[il:iu],ii)
      yp = take(pg.yp[il:iu],ii)
      wz1 = take(wz1,ii)
      np = len(ii)
      if np>0:
        # --- deposit ion charge in first grid
#        print 'dep q,wz1,js',q,wz1,js
        if self.l_planarYZ:
          rhoweightz_weight(yp,wz1,np,q,bg.nz,bg.dz,0.)
        else:
          rhoweightrz_weights(xp,yp,yp,wz1,np,q,bg.nr,bg.nz,bg.dr,bg.dz,bg.rmin,0.)
    il = pg.ins[js-1]-1
    iu = il+pg.nps[js-1]
    np=iu-il
    # --- reset 2-D rho arrays
    if np>0:
      q = pg.sq[js-1]*pg.sw[js-1]/w3d.dz
      zmin=w3d.zmminp+js*w3d.dz
      wz1 = (pg.zp[il:iu]-zmin)/w3d.dz
      # --- select particles in slice js
      ii = compress(floor(wz1)==0,arange(np))
      xp = take(pg.xp[il:iu],ii)
      yp = take(pg.yp[il:iu],ii)
      wz1 = take(wz1,ii)
      np = len(ii)
      if np>0:
        # --- deposit ion charge in first grid
#        print 'dep q,wz1,js',q,wz1,js
        if self.l_planarYZ:
          rhoweightz_weight(yp,wz1,np,q,bg.nz,bg.dz,0.)
        else:
          rhoweightrz_weights(xp,yp,yp,wz1,np,q,bg.nr,bg.nz,bg.dr,bg.dz,bg.rmin,0.)
    if self.l_verbose:print me,top.it,self.iz,'exit deposit_ions_last_step'
#    if lparallel:comm_world.barrier()
                      
  def sort_ions_along_z(self):
    if self.l_verbose:print me,top.it,self.iz,'enter sort_ions_along_z'
    if not self.fullsortdone:
      self.full_sort_ions_along_z()
      return
    pg = self.pgions
    if sum(pg.nps)==0:return
    for js in range(w3d.nzp):
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      zmin = js*w3d.dz+w3d.zmminp
      iz = floor((pg.zp[il:iu]-zmin)/w3d.dz) 
      if len(iz)>0:
        if min(iz)<-1 or max(iz)>1: 
          print 'error in sort_ions_along_z:iz',js,min(iz),max(iz),pg.zp[il:iu]
          raise Exception('error in sort_ions_along_z')
      izleft = il+compress(iz<0,arange(pg.nps[js]))
      izright = il+compress(iz>0,arange(pg.nps[js]))
#      if js==0 and len(izleft)>0:raise Exception('Error in sort_ions_along_z:js==0 and len(izleft)>0')
#      if js==w3d.nzp-1 and len(izright)>0:raise Exception('Error in sort_ions_along_z:js==w3d.nzp-1 and len(izright)>0')
      # --- Due to roundoff errors, it is possible that particles might be flagged as 
      if js==0:izleft=[]
      if js==w3d.nzp-1:izright=[]
      if len(izleft)>0:
        toaddleft = [take(pg.xp,izleft).copy(),
                     take(pg.yp,izleft).copy(),
                     take(pg.zp,izleft).copy(),
                     take(pg.uxp,izleft).copy(),
                     take(pg.uyp,izleft).copy(),
                     take(pg.uzp,izleft).copy(),
                     take(pg.ex,izleft).copy(),
                     take(pg.ey,izleft).copy(),
                     take(pg.ez,izleft).copy(),
                     take(pg.bx,izleft).copy(),
                     take(pg.by,izleft).copy(),
                     take(pg.bz,izleft).copy(),
                     take(pg.gaminv,izleft).copy()]
        if pg.npid>0:toaddleft.append(take(pg.pid,izleft,0).copy())
        put(pg.gaminv,izleft,0.)
      if len(izright)>0:
        toaddright = [take(pg.xp,izright).copy(),
                     take(pg.yp,izright).copy(),
                     take(pg.zp,izright).copy(),
                     take(pg.uxp,izright).copy(),
                     take(pg.uyp,izright).copy(),
                     take(pg.uzp,izright).copy(),
                     take(pg.ex,izright).copy(),
                     take(pg.ey,izright).copy(),
                     take(pg.ez,izright).copy(),
                     take(pg.bx,izright).copy(),
                     take(pg.by,izright).copy(),
                     take(pg.bz,izright).copy(),
                     take(pg.gaminv,izright).copy()]
        if pg.npid>0:toaddright.append(take(pg.pid,izright,0).copy())
        put(pg.gaminv,izright,0.)
      processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
      if len(izleft)>0:
        if pg.npid==0:
          pid = 0.
        else:
          pid = toaddleft[-1]
#        print 'add ',len(izleft),'from up' 
        addparticles(js=js-1,
                     x = toaddleft[0],
                     y = toaddleft[1],
                     z = toaddleft[2],
                     vx = toaddleft[3],
                     vy = toaddleft[4],
                     vz = toaddleft[5],
                     ex = toaddleft[6],
                     ey = toaddleft[7],
                     ez = toaddleft[8],
                     bx = toaddleft[9],
                     by = toaddleft[10],
                     bz = toaddleft[11],
                     gi = toaddleft[12],
                     pid = pid,
                     lmomentum=1,
                     lallindomain=1,
                     lfields=1)
      if len(izright)>0:
        if pg.npid==0:
          pid = 0.
        else:
          pid = toaddright[-1]
#        print 'add ',len(izright),'from down' 
        addparticles(js=js+1,
                     x = toaddright[0],
                     y = toaddright[1],
                     z = toaddright[2],
                     vx = toaddright[3],
                     vy = toaddright[4],
                     vz = toaddright[5],
                     ex = toaddright[6],
                     ey = toaddright[7],
                     ez = toaddright[8],
                     bx = toaddright[9],
                     by = toaddright[10],
                     bz = toaddright[11],
                     gi = toaddright[12],
                     pid = pid,
                     lmomentum=1,
                     lallindomain=1,
                     lfields=1)
    if self.l_verbose:print me,top.it,self.iz,'exit sort_ions_along_z'

  def full_sort_ions_along_z(self):
    if self.l_verbose:print me,top.it,self.iz,'enter full_sort_ions_along_z'
    pg = self.pgions
    top.pgroup = self.pgions
    if sum(pg.nps)==0:return
    js = 0
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    self.set_gamma(0)
    x=pg.xp[il:iu].copy()
    y=pg.yp[il:iu].copy()
    z=pg.zp[il:iu].copy()
    ux=pg.uxp[il:iu].copy()
    uy=pg.uyp[il:iu].copy()
    uz=pg.uzp[il:iu].copy()
    gi=pg.gaminv[il:iu].copy()
    if pg.npid>0:
      pid = pg.pid[il:iu,:].copy()
    pg.nps[0]=0
    for iz in range(w3d.nzp):
      zmin = iz*w3d.dz+w3d.zmminp+top.zgrid
      ii = compress((z>=zmin) & (z<zmin+w3d.dz),arange(shape(x)[0]))
      if len(ii)>0:
        if pg.npid>0:
          pida = pid.take(ii,0)
        else:
          pida = 0.
        addparticles(x=take(x,ii),
                     y=take(y,ii),
                     z=take(z,ii),
                     vx=take(ux,ii),
                     vy=take(uy,ii),
                     vz=take(uz,ii),
                     gi=take(gi,ii),
                     js=iz,
                     pid = pida,
                     lmomentum=1,
                     lallindomain=1)
      self.set_gamma(iz)
    self.fullsortdone = 1
    if self.l_verbose:print me,top.it,self.iz,'enter full_sort_ions_along_z'

  def compute_sw(self):
    if self.l_verbose:print me,top.it,self.iz,'enter compute_sw'
    del self.wz0,self.wz1
    self.wz0=[]
    self.wz1=[]
    for iz in range(w3d.nzp):
      wz0,wz1 = self.get_sw(iz)
      self.wz0.append(wz0)
      self.wz1.append(wz1)
    if self.l_verbose:print me,top.it,self.iz,'exit compute_sw'

  def get_sw(self,js):
      if self.l_verbose:print me,top.it,self.iz,'enter get_sw'
      pg = self.pgions
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      if pg.nps[js]==0:
        return None,None
      else:
        zmin=w3d.zmminp+js*w3d.dz
        wz1 = (pg.zp[il:iu]-zmin)/w3d.dz
        return 1.-wz1,wz1
      if self.l_verbose:print me,top.it,self.iz,'exit get_sw'

  def set_sw(self,js):
      if self.l_verbose:print me,top.it,self.iz,'enter set_sw'
      wz0,wz1 = self.get_sw(js)
      self.wz1[js] = wz1
      self.wz0[js] = wz0
      if self.l_verbose:print me,top.it,self.iz,'exit set_sw'
        
  def gather_ions_fields_old(self):
    if self.l_verbose:print me,top.it,self.iz,'enter push_ions'
    pg = self.pgions
    top.pgroup = pg
    if me>=(npes-(top.it+1)*self.ntgroups):
#    if me==(npes-(top.it+1)*self.ntgroups):
      js = self.iz
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      np = pg.nps[js]
      if np>0:
        pg.ex[il:iu]=0.
        pg.ey[il:iu]=0.
        pg.ez[il:iu]=0.
        pg.bx[il:iu]=0.
        pg.by[il:iu]=0.
        pg.bz[il:iu]=0.
        frz.basegrid = self.gridelecs[1]; g = frz.basegrid
        mk_grids_ptr()
        if self.l_planarYZ:
          fieldweightz(pg.yp[il:iu],pg.ey[il:iu],np,top.zgrid)
        else:
          fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,top.zgrid,top.efetch[0])
        pg.ex[il:iu]*=self.wz1[js]
        pg.ey[il:iu]*=self.wz1[js]
        if self.l_selfi:
          frz.basegrid = self.gridions[1]; g = frz.basegrid
          mk_grids_ptr()
          if self.l_planarYZ:
            fieldweightzb(pg.yp[il:iu],pg.bx[il:iu],np,top.zgrid)
          else:
            fieldweightxzb(pg.xp[il:iu],pg.yp[il:iu],pg.bx[il:iu],pg.by[il:iu],np,top.zgrid,top.efetch[0])
          pg.bx[il:iu]*=self.wz1[js]
          pg.by[il:iu]*=self.wz1[js]
        if self.iz==0:
          ex = zeros(np,'d')
          ey = zeros(np,'d')
          frz.basegrid = self.gridelecs[0]; g = frz.basegrid
          mk_grids_ptr()
          if self.l_planarYZ:
            fieldweightz(pg.yp[il:iu],ey,np,top.zgrid)
          else:
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],ex,ey,np,top.zgrid,top.efetch[0])
          ex*=self.wz0[js]
          ey*=self.wz0[js]
          pg.ex[il:iu]+=ex
          pg.ey[il:iu]+=ey
          if self.l_selfi:
            bx = zeros(np,'d')
            by = zeros(np,'d')
            frz.basegrid = self.gridions[0]; g = frz.basegrid
            mk_grids_ptr()
            if self.l_planarYZ:
              fieldweightzb(pg.yp[il:iu],bx,np,top.zgrid)
            else:
              fieldweightxzb(pg.xp[il:iu],pg.yp[il:iu],bx,by,np,top.zgrid,top.efetch[0])
            bx*=self.wz0[js]
            by*=self.wz0[js]
            pg.bx[il:iu]+=bx
            pg.by[il:iu]+=by
          self.add_other_fields(js)
      if self.iz<w3d.nzp-1:        
        js = self.iz+1
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        np = pg.nps[js]
        if np>0:
          ex = zeros(np,'d')
          ey = zeros(np,'d')
          frz.basegrid = self.gridelecs[1]; g = frz.basegrid
          mk_grids_ptr()
          if self.l_planarYZ:
            fieldweightz(pg.yp[il:iu],ey,np,top.zgrid)
          else:
            fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],ex,ey,np,top.zgrid,top.efetch[0])
          ex*=self.wz0[js]
          ey*=self.wz0[js]
          pg.ex[il:iu]+=ex
          pg.ey[il:iu]+=ey
          if self.l_selfi:
            bx = zeros(np,'d')
            by = zeros(np,'d')
            frz.basegrid = self.gridions[1]; g = frz.basegrid
            mk_grids_ptr()
            if self.l_planarYZ:
              fieldweightzb(pg.yp[il:iu],bx,np,top.zgrid)
            else:
              fieldweightxzb(pg.xp[il:iu],pg.yp[il:iu],bx,by,np,top.zgrid,top.efetch[0])
            bx*=self.wz0[js]
            by*=self.wz0[js]
            pg.bx[il:iu]+=bx
            pg.by[il:iu]+=by
          self.add_other_fields(js)
    if self.l_verbose:print me,top.it,self.iz,'exit push_ions'
             
  def gather_ions_fields(self):
    if self.l_verbose:print me,top.it,self.iz,'enter gather_ions_fields'
    pg = self.pgions
    top.pgroup = pg
    if me>=(npes-(top.it+1)*self.ntgroups):
#    if me==(npes-(top.it+1)*self.ntgroups):
      js = self.iz
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      np = pg.nps[js]
      if np>0:
        pg.ex[il:iu]=0.
        pg.ey[il:iu]=0.
        pg.ez[il:iu]=0.
        pg.bx[il:iu]=0.
        pg.by[il:iu]=0.
        pg.bz[il:iu]=0.
        if self.l_elec_greensum:
          self.ions_addefield_greensum(js)
        else:
          self.ions_fieldweight(js,self.gridelecs[1],self.wz1[js],l_add=0,l_bfield=0)
        if self.l_selfi:
          self.ions_fieldweight(js,self.gridions[1],self.wz1[js],l_add=0,l_bfield=1)
        if self.iz==0:
          if not self.l_elec_greensum:
            self.ions_fieldweight(js,self.gridelecs[0],self.wz0[js],l_add=1,l_bfield=0)
          if self.l_selfi:
            self.ions_fieldweight(js,self.gridions[0],self.wz0[js],l_add=1,l_bfield=1)
          self.add_other_fields(js)
      if self.iz<w3d.nzp-1:        
        js = self.iz+1
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        np = pg.nps[js]
        if np>0:
          if not self.l_elec_greensum:
            self.ions_fieldweight(js,self.gridelecs[1],self.wz0[js],l_add=1,l_bfield=0)
          if self.l_selfi:
            self.ions_fieldweight(js,self.gridions[1],self.wz0[js],l_add=1,l_bfield=1)
          self.add_other_fields(js)
    if self.l_verbose:print me,top.it,self.iz,'exit gather_ions_fields'
             
  def ions_fieldweight(self,js,g,w,l_bfield=0,l_add=0):
    pg = self.pgions
    top.pgroup = pg
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    np = pg.nps[js]
    if np==0:return
#    print 'fwb: iz,np,ex,ey,w',js,np,mean(pg.ex[il:iu]),std(pg.ex[il:iu]),mean(pg.ey[il:iu]),std(pg.ey[il:iu]),ave(w),std(w)
    frz.basegrid = g; mk_grids_ptr()
    if self.l_planarYZ:
      if l_bfield:
        fwxz = fieldweightzb
      else:
        fwxz = fieldweightz
    else:
      if l_bfield:
        fwxz = fieldweightxzb
      else:
        fwxz = fieldweightxz
    if l_add:
     ex=zeros(np,'d')
     ey=zeros(np,'d')
    else:
      if l_bfield:
        ex=pg.bx[il:iu]
        ey=pg.by[il:iu]
      else:
        ex=pg.ex[il:iu]
        ey=pg.ey[il:iu]
    if self.l_planarYZ:
      if l_bfield:
        fwxz(pg.yp[il:iu],ex,np,0.) # in this case, ex stands for bx
      else:
        fwxz(pg.yp[il:iu],ey,np,0.)
    else:
      fwxz(pg.xp[il:iu],pg.yp[il:iu],ex,ey,np,top.zgrid,top.efetch[0])
#    xe=getx(bcast=0,gather=0,pgroup=self.pgelec);ye=gety(bcast=0,gather=0,pgroup=self.pgelec);
#    print 'fw: iz,phi,rho,xe,ye',js,mean(g.phi),std(g.phi),mean(g.rho),std(g.rho),mean(xe),std(xe),mean(ye),std(ye)
    if l_add:
      if l_bfield:
        pg.bx[il:iu]+=ex*w
        pg.by[il:iu]+=ey*w
      else:
        pg.ex[il:iu]+=ex*w
        pg.ey[il:iu]+=ey*w
    else: 
      ex*=w
      ey*=w
#    print 'fwa: iz,ex,ey',js,mean(pg.ex[il:iu]),std(pg.ex[il:iu]),mean(pg.ey[il:iu]),std(pg.ey[il:iu])

  def ions_addefield_greensum(self,js):
    pgi = self.pgions
    npi = pgi.nps[js]
    pge = self.pgelec
    npe = pge.nps[0]
    if npi==0 or npe==0:return
    ili = pgi.ins[js]-1
    iui = ili+npi
    ile = pge.ins[0]-1
    iue = ile+npe
    w = pge.pid[ile:iue,top.wpid-1]*pge.sw[0]
    Lambda = w*pge.sq[0]  # note that the weight of the electrons is a line density
    for i in range(npi):
      x = pgi.xp[ili+i]-pge.xp[ile:iue]
      y = pgi.yp[ili+i]-pge.yp[ile:iue]
      aellipse = self.elec_aellipse
      bellipse = self.elec_bellipse
#      Ex, Ey = field_uniform_ellipse(Lambda,x,y,0.,0.,aellipse,bellipse)
      Ex, Ey = Bassetti_Erskine(Lambda,x,y,0.,0.,aellipse,bellipse)
      pgi.ex[ili+i]=sum(Ex)
      pgi.ey[ili+i]=sum(Ey)
      pgi.ez[ili+i]=0.


  def add_other_fields(self,js):
    if self.l_verbose:print me,top.it,self.iz,'add_other_fields'
    pg = self.pgions
    np = pg.nps[js]
    if np==0:return
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    if not self.l_maps:
      # --- Add in ears and uniform focusing E field pieces
      othere3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                  top.zbeam,top.zimax,top.zimin,top.straight,top.ifeears,top.eears,
                  top.eearsofz,top.dzzi,top.nzzarr,top.zzmin,
                  top.dedr,top.dexdx,top.deydy,top.dbdr,top.dbxdy,top.dbydx,
                  pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                  pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu])
      # --- Apply cancellation of E+VxB correction
#           --- Set quad, dipole E and B; All: Bz
#            call exteb3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
#     &                   gaminv(ipmin),-halfdt_s,halfdt_s,
#     &                   bx(ipmin),by(ipmin),bz(ipmin),
#     &                   ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),
#     &                   bendres,bendradi,fulldt_s)
#c           --- Correction to z on entry/exit to accelerator gap
#            call zgapcorr(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
#     &                    -halfdt_s, halfdt_s, fulldt_s, sm(1), sq(1), time)

  def push_ions_velocity_full(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter push_ions_velocity_first_half'
    pg = self.pgions
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
      self.set_gamma(js)
      # --- push velocity from magnetic field
      bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu], 
                  pg.sq[js],pg.sm[js],top.dt, top.ibpush)
      # --- push velocity from electric field (half step)
      epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                 pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu], 
                 pg.sq[js],pg.sm[js],0.5*top.dt)
      # --- update gamma
      self.set_gamma(js)
    if self.l_verbose:print me,top.it,self.iz,'exit push_ions_velocity_first_half'
    
  def push_ions_velocity_first_half(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter push_ions_velocity_first_half'
    pg = self.pgions
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
      self.set_gamma(js)
      # --- push velocity from magnetic field
      bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                  pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu], 
                  pg.sq[js],pg.sm[js],0.5*top.dt, top.ibpush)

    if self.l_verbose:print me,top.it,self.iz,'exit push_ions_velocity_first_half'
    
  def push_ions_velocity_second_half(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter push_ions_velocity_second_half'
    pg = self.pgions
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
      self.set_gamma(js)

    if self.l_verbose:print me,top.it,self.iz,'exit push_ions_velocity_second_half'
    
  def set_gamma(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter push_ions_velocity_second_half'
    pg = self.pgions
    np = pg.nps[js]
    if np==0:return
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    # --- update gamma
    gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
             top.gamadv,top.lrelativ)

    if self.l_verbose:print me,top.it,self.iz,'exit push_ions_velocity_second_half'
    
  def push_ions_positions(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter push_ions_positions'
    if me>=(npes-(top.it+1)*self.ntgroups):
      pg = self.pgions
      np = pg.nps[js]
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      if self.fbgainy>0.:
        ybar=ave(pg.yp[il:iu])
        pg.yp[il:iu]-=self.fbgainy*ybar
      if self.lattice is not None:
        if np>0:
          zp=pg.zp[il:iu].copy()
          self.lattice[self.ist].apply_transfer_map(pg,il,iu)
        if self.iz==w3d.nzp-2 and self.l_verbose:print 'push beam in ',self.lattice[self.ist].name
        self.ilcount[js]+=1
        istadd=1
        while self.lattice[(self.ist+istadd)%self.nst].L==0.:
          if np>0:
            self.lattice[(self.ist+istadd)%self.nst].apply_transfer_map(pg,il,iu)
          if self.iz==w3d.nzp-2:print 'push beam in ',self.lattice[(self.ist+istadd)%self.nst].name
          self.ilcount[js]+=1
          istadd+=1
        if js==0:self.ist=(self.ist+istadd-1)%self.nst
        if np>0:
          dz = max(abs(pg.zp[il:iu]-zp))
          if dz>w3d.dz:print 'Error in push_ions_positions: dz>w3d.dz',dz/w3d.dz
      else:
       if np==0:return
       if self.l_maps:
        if 0:
          apply_simple_map(np,
                           pg.xp[il:iu],
                           pg.yp[il:iu],
                           pg.uxp[il:iu],
                           pg.uyp[il:iu],
                           pg.uzp[il:iu],
                           self.maps.Mtx,
                           self.maps.Mty)
          if self.l_push_z and top.it%self.npushzperiod==0:
            self.maps.apply_synchrotron_motion(pg,il,iu,top.dt*self.npushzperiod)
#            pg.zp[il:iu]+=(pg.uzp[il:iu]*pg.gaminv[il:iu]-top.vbeam)*top.dt
        else:
          if self.bucketid==1:
            zbeam=0.
          else:
            zbeam=-(self.bucketid-1)*(w3d.zmmax-w3d.zmmin)/self.nbuckets
          self.maps.apply_transfer_map(pg,il,iu,top.dt,self.l_push_z,zbeam=zbeam)
        self.set_gamma(js)
       else:
        # --- push positions (average longitudinal velocity is removed)
        if np==0:return
        if self.l_push_z:
           xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                   pg.uxp[il:iu],pg.uyp[il:iu],
                   (pg.uzp[il:iu]-top.vbeam/pg.gaminv[il:iu]),
                   pg.gaminv[il:iu],top.dt)      
        else:
           xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                   pg.uxp[il:iu],pg.uyp[il:iu],
                   zeros(np,'d'),
                   pg.gaminv[il:iu],top.dt)      

      # --- apply transverse electric field kick
      if self.Ekick != 0.:
        if self.Bkick==self.bucketid and (top.it-(npes-1-me))%self.maps.nstations==0:
          E = self.Ekick
          if self.Akick is not None:
            try:
              E *= self.Akick[int(self.mystation[js]/self.maps.nstations)]
            except:
              pass
          dt = self.Lkick/(pg.gaminv[il:iu]*pg.uzp[il:iu])
          pg.yp[il:iu] += pg.uyp[il:iu]*pg.gaminv[il:iu]*dt \
                        + 0.5*pg.sq[js]*E*dt**2*pg.gaminv[il:iu]/pg.sm[js]
          pg.uyp[il:iu] += pg.sq[js]*E*dt/pg.sm[js]
          self.set_gamma(js)
        
      self.mystation[js]+=1

    if self.l_verbose:print me,top.it,self.iz,'exit push_ions_positions'

  def apply_ions_bndconditions(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter apply_ions_bndconditions'
    # --- apply boundary conditions
    pg = self.pgions
    if self.tparallel==1 and pg.nps[js]==0:return
    self.apply_bnd_conditions(js)
    if self.l_verbose:print me,top.it,self.iz,'exit apply_ions_bndconditions'
    
  def apply_bnd_conditions(self,js):
    if self.l_verbose:print me,top.it,self.iz,'enter apply_bnd_conditions'
    pg=top.pgroup=self.pgions
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    stckxy3d(top.pgroup,js,top.zbeam,true)
    if self.tparallel in [0,2]:
      xparticleboundaries(pg,js,js,w3d.xmmaxlocal,w3d.xmminlocal,True,w3d.l4symtry,False)
      yparticleboundaries(pg,js,js,w3d.ymmaxlocal,w3d.ymminlocal,True,w3d.l2symtry or w3d.l4symtry,False)
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
    else:
      partbndwithdata(pg.nps[js],pg.xp[il:iu],pg.uxp[il:iu],pg.gaminv[il:iu],
                      w3d.xmmax,w3d.xmmin,0.,
                      top.pboundxy,top.pboundxy)
      partbndwithdata(pg.nps[js],pg.yp[il:iu],pg.uyp[il:iu],pg.gaminv[il:iu],
                      w3d.ymmax,w3d.ymmin,0.,
                      top.pboundxy,top.pboundxy)
    if (js==0 or js==w3d.nzp-1) and pg.nps[js]>0:
      if js==0:top.pboundnz=-1
      if js==w3d.nzp-1:top.pbound0=-1
      partbndwithdata(pg.nps[js],pg.zp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                      w3d.zmmax,w3d.zmmin,top.zgrid,
                      top.pbound0,top.pboundnz)
      if js==0:top.pboundnz=0
      if js==w3d.nzp-1:top.pbound0=0
    if self.scraper is not None:
      self.scraper.updategrid()
      self.scraper.scrape(js)
    processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
    if self.l_verbose:print me,top.it,self.iz,'enter apply_bnd_conditions'

  def clear_electrons(self):
    # --- clear electrons on processor 0, shift them to the previous ones for me>0
    if self.l_verbose:print me,top.it,self.iz,'enter clear_electrons'
    pg = self.pgelec
    top.pgroup = self.pgelec
    js = 0
    # --- send electrons to previous processor
    if lparallel and me>self.ntgroups-1:
      tosend=[pg.nps[js]]
      if pg.nps[js]>0:
        il=top.pgroup.ins[js]-1
        iu=il+top.pgroup.nps[js]
        tosend.append(pg.xp[il:iu])
        tosend.append(pg.yp[il:iu])
        tosend.append(pg.zp[il:iu])
        tosend.append(pg.uxp[il:iu])
        tosend.append(pg.uyp[il:iu])
        tosend.append(pg.uzp[il:iu])
        tosend.append(pg.gaminv[il:iu])
        if top.npid>0:tosend.append(pg.pid[il:iu,:])
      mpisend(tosend, dest = me-self.ntgroups)
    # --- store electrons for diagnostics on all processors
    self.npelec_laststep=pg.nps[js]
    if 0:#pg.nps[js]>0:
      il=top.pgroup.ins[js]-1
      iu=il+top.pgroup.nps[js]
      self.xelec_laststep=pg.xp[il:iu].copy()
      self.yelec_laststep=pg.yp[il:iu].copy()
      self.uxelec_laststep=pg.uxp[il:iu].copy()
      self.uyelec_laststep=pg.uyp[il:iu].copy()
      self.uzelec_laststep=pg.uyp[il:iu].copy()
    # --- clear electrons on all processors
    if not self.l_periodicz:pg.nps[0]=0      
    # --- receive electrons from next processor
    if lparallel and me<npes-self.ntgroups:
      if self.l_posinst_track_electrons:
        pos.nlast=0
        pp_pos2warp(0)
        pos.time=0.
        pos.tbirth[:]=0.
      self.recved = mpirecv(source = me+self.ntgroups)
      np=self.recved[0]
      if np>0:
        if top.npid>0:
          pid=self.recved[8]
        else:
          pid=0.
        self.electrons.addpart(self.recved[1],
                              self.recved[2],
                              w3d.zmmaxp-1.e-10*w3d.dz,#
                              #self.recved[3],
                              self.recved[4],
                              self.recved[5],
                              self.recved[6],
                              self.recved[7],
                              pid=pid,
                              lmomentum=1,
                              lallindomain=1)
    # --- if periodic in z, first processor sends electrons to last
    if self.l_periodicz:      
      # --- send electrons to previous processor
      if lparallel and me==0:
        tosend=[pg.nps[js]]
        if pg.nps[js]>0:
          il=top.pgroup.ins[js]-1
          iu=il+top.pgroup.nps[js]
          tosend.append(pg.xp[il:iu])
          tosend.append(pg.yp[il:iu])
          tosend.append(pg.zp[il:iu])
          tosend.append(pg.uxp[il:iu])
          tosend.append(pg.uyp[il:iu])
          tosend.append(pg.uzp[il:iu])
          tosend.append(pg.gaminv[il:iu])
          if top.npid>0:tosend.append(pg.pid[il:iu,:])
        mpisend(tosend, dest = npes-1)
        pg.nps[0]=0      
      if self.l_verbose:print me,top.it,self.iz,'exit clear_electrons'
      # --- receive electrons from next processor
      if lparallel and me==npes-1:
        if self.l_posinst_track_electrons:
          pos.nlast=0
          pp_pos2warp(0)
          pos.time=0.
          pos.tbirth[:]=0.
        self.recved = mpirecv(source = 0)
        np=self.recved[0]
        if np>0:
          if top.npid>0:
            pid=self.recved[8]
          else:
            pid=0.
          self.electrons.addpart(self.recved[1],
                                self.recved[2],
                                w3d.zmmaxp-1.e-10*w3d.dz,#
                                #self.recved[3],
                                self.recved[4],
                                self.recved[5],
                                self.recved[6],
                                self.recved[7],
                                pid=pid,
                                lmomentum=1,
                                lallindomain=1)
#    if lparallel:comm_world.barrier()
    
  def plot_electrons(self):
        if self.l_verbose:print me,top.it,self.iz,'enter plot_electrons'
#        ppgeneric(getvx(js=1),getx(js=1))
        window(0);fma();
#        frz.basegrid=self.gridions[1]
#        mk_grids_ptr()
#        plphirz();refresh()
#        window(0);fma();
#        n=getn(js=1)
#        if n>0:
#          x = self.slist[0].getx()
#          y = self.slist[0].gety()
#          ex = self.slist[0].getex()
#          if me==0:ppco(x,y,ex,msize=3,ncolor=100);refresh()
#        if iz==0:
#        self.electrons.ppxy(msize=2);refresh()
        ppxy(pgroup=self.pgelec,msize=4)
        print 'plot_elec',self.iz,'nb electrons = ',self.pgelec.nps[0],self.electrons.getn()
#        self.electrons.ppxy(msize=5);refresh()#msize=1,color='density',ncolor=100,bcast=0,gather=0);refresh()
#        else:
#          self.slist[0].ppxex(msize=2);
#        limits(w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
#        ppxex(js=1)
#        pfxy(iz=iz);
#        fma()
        if self.l_verbose:print me,top.it,self.iz,'enter plot_electrons'

  def create_electrons(self):
     if self.l_verbose:print me,top.it,self.iz,'enter create_electrons'
     pg = self.pgelec
     top.pgroup = self.pgelec
     if self.l_posinst_track_electrons:
       pp_pos2warp(0)
       pos.time=0.
       pos.tbirth[:]=0.
     # --- generate electrons (on last processor only)
     if not self.l_periodicz:pg.nps[0]=0
     if self.l_mode==2:
       zmax = w3d.zmmax/2+top.zgrid
     else:
       zmax = self.zmmaxlocal-1.e-10*w3d.dz+top.zgrid
     if self.xinitelec is not None: 
       if self.winitelec is None:
         weight=1.
       else:
         weight=self.winitelec[self.ielecstride::self.strideinitelec]
       xinit = self.xinitelec[self.ielecstride::self.strideinitelec]
       yinit = self.yinitelec[self.ielecstride::self.strideinitelec]
       if self.l_randomflip or self.l_randinit:
         xinit=xinit.copy()       
         yinit=yinit.copy()       
       if self.l_randomflip:
        if 0:
         if ranf()>0.5:
           xinit*=-1.
         if ranf()>0.5:
           yinit*=-1.
        else:
         xinit = where(ranf(xinit)>0.5,xinit,-xinit)
         yinit = where(ranf(yinit)>0.5,yinit,-yinit)
       if self.l_randinit:
#         xinit+=(ranf(xinit)-0.5)*w3d.dx*0.25
         yinit+=(ranf(yinit)-0.5)*w3d.dy*0.25
#         xinit = where(xinit>=w3d.xmmax,w3d.xmmax-1.e-10*w3d.dx,xinit)
#         xinit = where(xinit<=w3d.xmmin,w3d.xmmin+1.e-10*w3d.dx,xinit)
         yinit = where(yinit>=w3d.ymmax,w3d.ymmax-1.e-10*w3d.dy,yinit)
         yinit = where(yinit<=w3d.ymmin,w3d.ymmin+1.e-10*w3d.dy,yinit)
       self.electrons.addpart(xinit,yinit,xinit*0.+zmax,
                         self.vxinitelec[self.ielecstride::self.strideinitelec],
                         self.vyinitelec[self.ielecstride::self.strideinitelec],
                         self.vzinitelec[self.ielecstride::self.strideinitelec],
                         w=weight,
                         lmomentum=false)
       if self.l_elec4sym:
        for isym in range(3):
         if isym==0:
           sx=-1.;sy=1.
         if isym==1:
           sx=1.;sy=-1.
         if isym==2:
           sx=-1.;sy=-1.
         self.electrons.addpart(sx*xinit,sy*yinit,xinit*0.+zmax,
                                sx*self.vxinitelec[self.ielecstride::self.strideinitelec],
                                sy*self.vyinitelec[self.ielecstride::self.strideinitelec],
                                self.vzinitelec[self.ielecstride::self.strideinitelec],
                                w=weight,
                                lmomentum=false)
       self.ielecstride = (self.ielecstride+1)%self.strideinitelec
     else: 
       if self.l_inject_elec_MR:
         js = self.electrons.jslist[0]
         g = frz.basegrid; mk_grids_ptr()
         xmin = g.xmin+g.transit_min_r*g.dr
         xmax = g.xmax-g.transit_max_r*g.dr
         ymin = g.zmin+g.transit_min_z*g.dz
         ymax = g.zmax-g.transit_max_z*g.dz
         for i in range(frz.ngrids):
           x,y = getmesh2d(xmin+g.dr/2,g.dr,g.nr-g.transit_min_r-g.transit_max_r-1, 
                           ymin+g.dz/2,g.dz,g.nz-g.transit_min_z-g.transit_max_z-1)
           np = product(shape(x))
           x.shape = (np,)
           y.shape = (np,)
           if i<frz.ngrids-1:
             g=g.down
             xmin = g.xmin+g.transit_min_r*g.dr
             xmax = g.xmax-g.transit_max_r*g.dr
             ymin = g.zmin+g.transit_min_z*g.dz
             ymax = g.zmax-g.transit_max_z*g.dz
             ii = compress( (x>xmax) | (x<xmin) | (y>ymax) | (y<ymin), arange(np))
             x = take(x,ii)
             y = take(y,ii)
           self.electrons.addpart(x,y,zmax,1.e-10,1.e-10,1.e-10,w=1./(4.**i))
       else:
         if top.prwall<sqrt(w3d.xmmax**2+w3d.ymmax**2):
           self.electrons.add_uniform_cylinder(self.Ninit,min(top.prwall,w3d.xmmax),zmax,zmax,1.e-10,1.e-10,1.e-10,lallindomain=true)
         else:
           if self.l_elecuniform:
             spacing='uniform'
           else:
             spacing='random'
           self.electrons.add_uniform_box(self.Ninit,self.xminelec,self.xmaxelec,
                                          self.yminelec,self.ymaxelec,zmax,zmax,
                                          1.e-10,1.e-10,1.e-10,
                                          nx=self.nxelec,ny=self.nyelec,spacing=spacing,lallindomain=false)
           
       if self.l_reuseelec:
         js = self.electrons.jslist[0]
         il = pg.ins[js]-1
         iu = il+pg.nps[js]
         self.xinitelec=pg.xp[il:iu].copy()
         self.yinitelec=pg.yp[il:iu].copy()
         self.vxinitelec=pg.uxp[il:iu]*pg.gaminv[il:iu]
         self.vyinitelec=pg.uyp[il:iu]*pg.gaminv[il:iu]
         self.vzinitelec=pg.uzp[il:iu]*pg.gaminv[il:iu]
     # --- scrape electrons 
     if self.scraper is not None:
       self.scraper.scrapeall(local=1,clear=1)
#     shrinkpart(pg)
     if self.l_verbose:print me,top.it,self.iz,'exit create_electrons'

  def deposit_electrons(self,i=1):
      if self.l_verbose:print me,top.it,self.iz,'enter deposit_electrons'
      pg = self.pgelec
      top.pgroup = self.pgelec
      bg = frz.basegrid = self.gridelecs[i]
      mk_grids_ptr()
      # --- reset 2-D rho arrays
      reset_rzmgrid_rho()
      bg.rho=0.
      # --- loop over species list
      for sp in [self.electrons]:
        # --- loop over species index
        for js in sp.jslist:
          np = pg.nps[js]
          if np==0:continue
          il = pg.ins[js]-1
          iu = il+pg.nps[js]
          if self.l_verbose:print me,top.it,self.iz,'me = ',me,';deposit rho' # MR OK
          # --- deposit rho
          if top.wpid==0:
            if self.l_planarYZ:
              rhoweightz(pg.yp[il:iu],np,pg.sq[js]*pg.sw[js],bg.nz,bg.dz,0.)
            else:
              rhoweightrz(pg.xp[il:iu],pg.yp[il:iu],pg.yp[il:iu],np,
                        pg.sq[js]*pg.sw[js],bg.nr,bg.nz,
                        bg.dr,bg.dz,bg.rmin,0.)
          else:
            if self.l_planarYZ:
              rhoweightz_weight(pg.yp[il:iu],pg.pid[:,top.wpid-1],np,pg.sq[js]*pg.sw[js],bg.nz,bg.dz,0.)
            else:
              rhoweightrz_weights(pg.xp[il:iu],pg.yp[il:iu],pg.yp[il:iu],
                        pg.pid[:,top.wpid-1],np,
                        pg.sq[js]*pg.sw[js],bg.nr,bg.nz,
                        bg.dr,bg.dz,bg.rmin,0.)

      # --- exchange data between procs
      if self.ntgroups>1:get_rho_from_rhop(frz.basegrid)
      # --- distribute rho among patches
      distribute_rho_rz()
      if self.l_verbose:print me,top.it,self.iz,'exit deposit_electrons'

  def push_electrons(self):
       if self.l_verbose:print me,top.it,self.iz,'enter push_electrons'
       self.nelechist.append(self.pgelec.nps[0])
       pg = self.pgelec
       top.pgroup = self.pgelec
       frz.basegrid = self.gridions[1]
       mk_grids_ptr()
       if top.bx0 != 0. or top.by0 != 0.:
         bx0 = ones(self.nparpgrp,'d')*top.bx0
         by0 = ones(self.nparpgrp,'d')*top.by0
         bz0 = ones(self.nparpgrp,'d')*top.bz0
         self.l_bfield=1
       else:
         self.l_bfield=0
       # --- loop over species list
       for sp in [self.electrons]:
        # --- loop over species index
        for js in sp.jslist:
          ng = 1+pg.nps[js]/self.nparpgrp
          if pg.nps[js]>0:
           for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if np==0:continue
            # --- gather self-forces
            if self.l_bfield:
              if self.l_verbose:print me,top.it,self.iz,'me = ',me,';bpush'
              bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                          bx0[:np], by0[:np], bz0[:np], pg.sq[js],pg.sm[js],0.5*self.dt, top.ibpush)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';fieldweight' # MR OK
            if self.l_beam_1d:
#              Lambda = self.pgions.sw[self.iz]*self.pgions.sq[self.iz]/w3d.dz
              Lambda = self.pgions.sw[self.iz]*self.pgions.sq[self.iz]*self.pgions.nps[self.iz]/w3d.dz
              x0 = self.pgions.xp[self.iz]
              y0 = self.pgions.yp[self.iz]
              pos.sigx = self.sigmax[self.iz]
              pos.sigy = self.sigmay[self.iz]
              Exbe, Eybe = allefield(Lambda,pg.xp[il:iu],pg.yp[il:iu],x0,y0)
#              Exbe, Eybe = Bassetti_Erskine(Lambda,pg.xp[il:iu],pg.yp[il:iu],x0,y0,sigmax,sigmay)
              pg.ex[il:iu]=Exbe
              pg.ey[il:iu]=Eybe
              pg.ez[il:iu]=0.
            else:
              pg.ex[il:iu]=0.
              pg.ey[il:iu]=0.
              pg.ez[il:iu]=0.
              if self.l_planarYZ:
                fieldweightz(pg.yp[il:iu],pg.ey[il:iu],np,0.)
              else:
                fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,0.,top.efetch[js])
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';xpush'
#            xe=getx(bcast=0,gather=0,pgroup=self.pgelec);ye=gety(bcast=0,gather=0,pgroup=self.pgelec)
#            print 'peb:iz,xe,ye',self.iz,ave(xe),std(xe),ave(ye),std(ye)
            if not self.l_posinst_track_electrons:
              if self.l_freeze_xelec:
                xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                        pg.uxp[il:iu]*0.,pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],self.dt)
              else:
                xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                        pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],self.dt)
#            xe=getx(bcast=0,gather=0,pgroup=self.pgelec);ye=gety(bcast=0,gather=0,pgroup=self.pgelec)
#            print 'pea:iz,xe,ye',self.iz,ave(xe),std(xe),ave(ye),std(ye)
            pg.zp[il:iu]-=w3d.dz
#            window(2)
#            ppg(frz.basegrid.phi);refresh()
#            fma();ppxy(pgroup=pg,js=js,color=red,msize=2);refresh()
#            print pg.nps
#            window(0)
#            window(3);ppzvx(js=1,msize=2);window(0)

          if not self.l_posinst_track_electrons:
            if self.scraper is not None:
#              top.ns=self.pgelec.ns
#              gchange('LostParticles')
              self.scraper.scrapeall(local=1,clear=self.sec is None)
              if self.sec is not None:
                self.sec.generate(local=1,l_accumulate_hist=0)
#                processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+self.dt*top.pgroup.ndts[js],top.zbeam)
            if self.tparallel in [0,2]:
#              particleboundariesxy(top.pgroup,js,lcallcontrollers=false)
              stckxy3d(top.pgroup,js,top.zbeam,true)
              xparticleboundaries(top.pgroup,js,js,w3d.xmmaxlocal,w3d.xmminlocal,True,w3d.l4symtry,False)
              yparticleboundaries(top.pgroup,js,js,w3d.ymmaxlocal,w3d.ymminlocal,True,w3d.l2symtry or w3d.l4symtry,False)
              processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+self.dt*top.pgroup.ndts[js],top.zbeam)
            else:
              if self.l_verbose:print me,top.it,self.iz,'me = ',me,';stckxy3d'
              if pg.nps[js]>0:
                stckxy3d(top.pgroup,js,top.zbeam,true)
                il = pg.ins[js]-1
                iu = pg.ins[js]-1+pg.nps[js]
                partbndwithdata(pg.nps[js],pg.xp[il:iu],pg.uxp[il:iu],pg.gaminv[il:iu],
                                w3d.xmmax,w3d.xmmin,0.,
                                top.pboundxy,top.pboundxy)
                partbndwithdata(pg.nps[js],pg.yp[il:iu],pg.uyp[il:iu],pg.gaminv[il:iu],
                                w3d.ymmax,w3d.ymmin,0.,
                                top.pboundxy,top.pboundxy)
                processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+self.dt*top.pgroup.ndts[js],top.zbeam)
 #             top.ns=self.pgions.ns
 #             gchange('LostParticles')
          else: # track electrons with Posinst
            if pg.nps[js]==0:continue
            pp_warp2pos()
            pos.dtim=self.dt
            reset_tl()
            pos.Bfunifx=top.bx0
            pos.Bfunify=top.by0
            pos.Bfunifz=top.bz0
            if pos.clight==0.:pos.define_physical_constants()
#            pos.om0=(pos.clight**2/pos.emass)*pos.
            pytrack_electrons()
            pos.time+=self.dt
            pp_pos2warp(0)
          top.npslost=0
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if np==0:continue
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';fieldweight' # MR OK
            if self.l_beam_1d:
#              Lambda = self.pgions.sw[self.iz]*self.pgions.sq[self.iz]/w3d.dz
              Lambda = self.pgions.sw[self.iz]*self.pgions.sq[self.iz]*self.pgions.nps[self.iz]/w3d.dz
              x0 = self.pgions.xp[self.iz]
              y0 = self.pgions.yp[self.iz]
              pos.sigx = self.sigmax[self.iz]
              pos.sigy = self.sigmay[self.iz]
              Exbe, Eybe = allefield(Lambda,pg.xp[il:iu],pg.yp[il:iu],x0,y0)
#              Exbe, Eybe = Bassetti_Erskine(Lambda,pg.xp[il:iu],pg.yp[il:iu],x0,y0,sigmax,sigmay)
              pg.ex[il:iu]=Exbe
              pg.ey[il:iu]=Eybe
              pg.ez[il:iu]=0.
            else:
              pg.ex[il:iu]=0.
              pg.ey[il:iu]=0.
              pg.ez[il:iu]=0.
              if self.l_planarYZ:
                fieldweightz(pg.yp[il:iu],pg.ey[il:iu],np,0.)
              else:
                fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,top.zgrid,top.efetch[js])
#            fetche3d(pg,il+1,np,js+1)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
            if self.l_bfield:
              if self.l_verbose:print me,top.it,self.iz,'me = ',me,';bpush'
              bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                          bx0[:np], by0[:np], bz0[:np], pg.sq[js],pg.sm[js],0.5*self.dt, top.ibpush)

       if self.l_gas_ionization:
            pgbf=self.pgelec.nps.copy()
            top.pgroup = self.pgions
            self.ioniz[self.iz].generate(dt=self.dt)
#            print 'Ioniz',pgbf,self.pgelec.nps

       if self.l_verbose:print me,top.it,self.iz,'exit push_electrons'
#       js = 0
#       il = pg.ins[js]-1
#       iu = pg.ins[js]-1+pg.nps[js]
#       if pg.nps[js]>0:
#         self.elec.append(pg.yp[il:iu].copy())

  def smooth_rho(self,grid,smoother):
    g=grid
    for ig in range(frz.ngrids):
      if ig>0:
        g = g.down
      smoother.apply(g.rho)
      
  def sum_rho_in_tgroups(self,g):
    if self.tparallel==2:return
    if self.tparallel==0:
      # --- exchange data between procs
      g.rho[...]=parallelsum(g.rho,comm=self.mympitgroup)
      for ig in range(frz.ngrids):
        if ig>0:
          g = g.down
          self.exchange_rho(g)
    else:
      for ig in range(frz.ngrids):
        if ig>0:
          g = g.down
        g.rho[...]=parallelsum(g.rho,comm=self.mympitgroup)
      
  def exchange_rho(self,g):
    # send data up in X
    if top.procneighbors[1,0]>me:
      nin = 2*(g.transit_max_r)+1
      comm_world.isend(g.rho[-2*nin+1:-nin+1,:],top.procneighbors[1,0])
    # send data down in X
    if top.procneighbors[0,0]<me:
      nin = 2*(g.transit_min_r)+1
      comm_world.isend(g.rho[nin-1:2*nin-1,:],top.procneighbors[0,0])
    # recv data from down in X
    if top.procneighbors[0,0]<me:
      nin = 2*(g.transit_min_r)+1
      g.rho[:nin,:]+=mpirecv(source = top.procneighbors[0,0])
    # recv data from up in X
    if top.procneighbors[1,0]>me:
      nin = 2*(g.transit_max_r)+1
      g.rho[-nin:,:]+=mpirecv(source = top.procneighbors[1,0])
    # send data up in Y
    if top.procneighbors[1,1]>me:
      nin = 2*(g.transit_max_z)+1
      comm_world.isend(g.rho[:,-2*nin+1:-nin+1],top.procneighbors[1,1])
    # send data down in Y
    if top.procneighbors[0,1]<me:
      nin = 2*(g.transit_min_z)+1
      comm_world.isend(g.rho[:,nin-1:2*nin-1],top.procneighbors[0,1])
    # recv data from down in Y
    if top.procneighbors[0,1]<me:
      nin = 2*(g.transit_min_z)+1
      g.rho[:,:nin]+=mpirecv(source = top.procneighbors[0,1])
    # recv data from up in Y
    if top.procneighbors[1,1]>me:
      nin = 2*(g.transit_max_z)+1
      g.rho[:,-nin:]+=mpirecv(source = top.procneighbors[1,1])

  def add_ei_fields(self,i=1):
    if self.l_verbose:print me,top.it,self.iz,'enter add_ei_fields'
    if self.l_selfi:self.getselfb(i)
    
    if self.l_selfe:
        ge = self.gridelecs[i]
        gi = self.gridions[i]
        for ig in range(frz.ngrids):
          if ig>0:
            ge = ge.down
            gi = gi.down
          gi.phi[...] += ge.phi[...] 
    if self.l_selfi:
        ge = self.gridelecs[i]
        gi = self.gridions[i]
        for ig in range(frz.ngrids):
          if ig>0:
            ge = ge.down
            gi = gi.down
          if self.l_selfe:
            ge.phi[...] = gi.phi[...]
          else:
            ge.phi[...] += gi.phi[...] 
    if frz.l_get_fields_on_grid:
        ge = self.gridelecs[i]
        gi = self.gridions[i]
        for g in [ge,gi]:
          frz.basegrid = g
          mk_grids_ptr()
          if self.lparallel:
            getphiforparticlesrz()
          else:
            getallfieldsfromphip()

    if self.l_verbose:print me,top.it,self.iz,'exit add_ei_fields'

  def getselfb(self,i=1):
    gi = self.gridions[i]
    for ig in range(frz.ngrids):
      if ig>0:gi=gi.down
      fact = (self.pgions.fselfb[0]/(clight*clight))
      if self.l_planarYZ:
        gi.brp[:,:] = fact*(gi.phi[0,2:] - gi.phi[0,:-2])/(2.*gi.dz)
        gi.bzp[...] = 0.
      else:
        gi.brp[:,:] =  fact*(gi.phi[1:-1,2:] - gi.phi[1:-1,:-2])/(2.*gi.dz)
        gi.bzp[:,:] = -fact*(gi.phi[2:,1:-1] - gi.phi[:-2,1:-1])/(2.*gi.dr)
      
#  def adddadttoe(self):
#    """Ez = -dA/dt = -beta**2 dphi/dz"""
#    # --- This assumes that nzguard is always 1
#    self.Ez[...] = (self.pgions.fselfb[0]/clight)**2*(potentialp[ix,iy,2:]-potentialp[ix,iy,:-2])/(2.*self.dz)
#    fieldp[2,:,:,:,0] += Ez

  def zerommnts(self):
        if self.l_verbose:print me,top.it,self.iz,'enter zerommnts'
        pg = self.pgions
        self.pnumztmp[:]=0.
        self.xbarztmp[:]=0.     
        self.ybarztmp[:]=0.      
        self.zbarztmp[:]=0.      
        self.xpbarztmp[:]=0.     
        self.ypbarztmp[:]=0.    
        self.xpnbarztmp[:]=0.      
        self.ypnbarztmp[:]=0.      
        self.x2ztmp[:]=0.     
        self.y2ztmp[:]=0.     
        self.z2ztmp[:]=0.      
        self.xp2ztmp[:]=0.      
        self.yp2ztmp[:]=0.     
        self.xpn2ztmp[:]=0.    
        self.ypn2ztmp[:]=0.     
        self.xxpbarztmp[:]=0.      
        self.yypbarztmp[:]=0.  
        self.xxpnbarztmp[:]=0.   
        self.yypnbarztmp[:]=0. 
        if self.l_zdediag:
          self.denstmp[...]=0.
          self.dketmp[...]=0.

  def getmmnts(self,js):
      if self.l_verbose:print me,top.it,self.iz,'enter getmmnts'
      pg = self.pgions
      if pg.nps[js]==0:return
      il = pg.ins[js]-1
      iu = il+pg.nps[js]
      x = pg.xp[il:iu]
      y = pg.yp[il:iu]
      z = pg.zp[il:iu]
      gaminv = pg.gaminv[il:iu]
      ux = pg.uxp[il:iu]
      uy = pg.uyp[il:iu]
      uz = pg.uzp[il:iu]
      xp = ux/uz
      yp = uy/uz
      beta = sqrt((1.-gaminv)*(1.+gaminv))
      gamma = 1./gaminv
      if self.l_zdediag:
        dke = gamma*pg.sm[js]*clight**2/self.ke0-1.
      if self.lattice is not None:
       ist = (self.ist-1)%self.nst
       if self.lattice[ist].dispx is not None:
        dpp = getuz(js=js,gather=0,pgroup=pg)/(top.vbeam*top.gammabar)-1.
        x = x-self.lattice[ist].dispx*dpp
        y = y-self.lattice[ist].dispy*dpp
        xp = xp-self.lattice[ist].disppx*dpp
        yp = yp-self.lattice[ist].disppy*dpp
      xpn = xp*beta*gamma
      ypn = yp*beta*gamma
      x2 = x*x
      y2 = y*y
      z2 = z*z
      xp2 = xp*xp
      yp2 = yp*yp
      xpn2 = xpn*xpn
      ypn2 = ypn*ypn
      xxp = x*xp
      yyp = y*yp
      xxpn = x*xpn
      yypn = y*ypn
      w1 = (z-self.znodes[js])/w3d.dz
      w0 = 1.-w1
      if self.l_zdediag:
        self.denstmp0[...]=0.
        self.dketmp0[...]=0.
        self.denstmp1[...]=0.
        self.dketmp1[...]=0.
        deposgrid1d(1,pg.nps[js],dke,y*w0,self.nde,self.dketmp0,self.denstmp0,-self.demax,self.demax)
        deposgrid1d(1,pg.nps[js],dke,y*w1,self.nde,self.dketmp1,self.denstmp1,-self.demax,self.demax)
      self.pnumztmp[js]=sum(w0)
      self.xbarztmp[js]=sum(x*w0)      
      self.ybarztmp[js]=sum(y*w0)      
      self.zbarztmp[js]=sum(z*w0)      
      self.xpbarztmp[js]=sum(xp*w0)      
      self.ypbarztmp[js]=sum(yp*w0)      
      self.xpnbarztmp[js]=sum(xpn*w0)      
      self.ypnbarztmp[js]=sum(ypn*w0)      
      self.x2ztmp[js]=sum(x2*w0)      
      self.y2ztmp[js]=sum(y2*w0)      
      self.z2ztmp[js]=sum(z2*w0)      
      self.xp2ztmp[js]=sum(xp2*w0)      
      self.yp2ztmp[js]=sum(yp2*w0)      
      self.xpn2ztmp[js]=sum(xpn2*w0)      
      self.ypn2ztmp[js]=sum(ypn2*w0)      
      self.xxpbarztmp[js]=sum(xxp*w0)      
      self.yypbarztmp[js]=sum(yyp*w0)      
      self.xxpnbarztmp[js]=sum(xxpn*w0)      
      self.yypnbarztmp[js]=sum(yypn*w0)      
      if self.l_zdediag:
        self.dketmp[js,:]=self.dketmp0
        self.denstmp[js,:]=self.denstmp0
      self.pnumztmp[js+1]+=sum(w1)
      self.xbarztmp[js+1]+=sum(x*w1)      
      self.ybarztmp[js+1]+=sum(y*w1)      
      self.zbarztmp[js+1]+=sum(z*w1)      
      self.xpbarztmp[js+1]+=sum(xp*w1)      
      self.ypbarztmp[js+1]+=sum(yp*w1)      
      self.xpnbarztmp[js+1]+=sum(xpn*w1)      
      self.ypnbarztmp[js+1]+=sum(ypn*w1)      
      self.x2ztmp[js+1]+=sum(x2*w1)      
      self.y2ztmp[js+1]+=sum(y2*w1)      
      self.z2ztmp[js+1]+=sum(z2*w1)      
      self.xp2ztmp[js+1]+=sum(xp2*w1)      
      self.yp2ztmp[js+1]+=sum(yp2*w1)      
      self.xpn2ztmp[js+1]+=sum(xpn2*w1)      
      self.ypn2ztmp[js+1]+=sum(ypn2*w1)      
      self.xxpbarztmp[js+1]+=sum(xxp*w1)      
      self.yypbarztmp[js+1]+=sum(yyp*w1)      
      self.xxpnbarztmp[js+1]+=sum(xxpn*w1)      
      self.yypnbarztmp[js+1]+=sum(yypn*w1)      
      if self.l_zdediag:
        self.dketmp[js+1,:]+=self.dketmp1
        self.denstmp[js+1,:]+=self.denstmp1
      if self.l_verbose:print me,top.it,self.iz,'exit getmmnts'

  def getmmnts_store(self):
    if self.l_verbose:print me,top.it,self.iz,'enter getmmnts_store'
    pg = self.pgions
    self.pnum.append(sum(self.pnumztmp))
    self.xbar.append(sum(self.xbarztmp))
    self.ybar.append(sum(self.ybarztmp))
    self.zbar.append(sum(self.zbarztmp))
    self.xpbar.append(sum(self.xpbarztmp))
    self.ypbar.append(sum(self.ypbarztmp))
    self.x2.append(sum(self.x2ztmp))
    self.y2.append(sum(self.y2ztmp))
    self.z2.append(sum(self.z2ztmp))
    self.xp2.append(sum(self.xp2ztmp))
    self.yp2.append(sum(self.yp2ztmp))
    self.xxpbar.append(sum(self.xxpbarztmp))
    self.yypbar.append(sum(self.yypbarztmp))
    self.xpnbar.append(sum(self.xpnbarztmp))
    self.ypnbar.append(sum(self.ypnbarztmp))
    self.xpn2.append(sum(self.xpn2ztmp))
    self.ypn2.append(sum(self.ypn2ztmp))
    self.xxpnbar.append(sum(self.xxpnbarztmp))
    self.yypnbar.append(sum(self.yypnbarztmp))
    self.pnumz.append(self.pnumztmp)
    self.xbarz.append(self.xbarztmp)
    self.ybarz.append(self.ybarztmp)
    self.zbarz.append(self.zbarztmp)
    self.xpbarz.append(self.xpbarztmp)
    self.ypbarz.append(self.ypbarztmp)
    self.x2z.append(self.x2ztmp)
    self.y2z.append(self.y2ztmp)
    self.z2z.append(self.z2ztmp)
    self.xp2z.append(self.xp2ztmp)
    self.yp2z.append(self.yp2ztmp)
    self.xxpbarz.append(self.xxpbarztmp)
    self.yypbarz.append(self.yypbarztmp)
    self.xpnbarz.append(self.xpnbarztmp)
    self.ypnbarz.append(self.ypnbarztmp)
    self.xpn2z.append(self.xpn2ztmp)
    self.ypn2z.append(self.ypn2ztmp)
    self.xxpnbarz.append(self.xxpnbarztmp)
    self.yypnbarz.append(self.yypnbarztmp)
    if self.l_zdediag:
      self.zdedipole.append(self.dketmp)
      self.zdedens.append(self.denstmp)
    self.timemmnts.append(top.time)
#    if top.it==0 or not lparallel:
#      self.timemmnts.append(top.time)
#    else:
#      self.timemmnts.append(top.time+(me-npes+1)*top.dt)
    self.zerommnts()
    if self.l_verbose:print me,top.it,self.iz,'exit getmmnts_store'
      
  def parallelsum(self,a):
      if self.ntgroups>1:a = parallelsum(a,comm=self.mympitgroup)
      if self.mympigroup is None:
         return a
      else:
         return parallelsum(a,comm=self.mympigroup)
  
  def parallelmin(self,a):
      if self.ntgroups>1:a = parallelmin(a,comm=self.mympitgroup)
      if self.mympigroup is None:
         return a
      else:
         return parallelmin(a,comm=self.mympigroup)
  
  def getpnum(self):
      n = self.parallelmin(len(self.pnum.data()))
      return self.parallelsum(self.pnum.data()[:n])

  def getxrms(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.x2.data()/pnum)
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.parallelsum(self.x2.data()[:n])/pnum)   

  def getyrms(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.y2.data()/pnum)
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.parallelsum(self.y2.data()[:n])/pnum)   

  def getzrms(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.z2.data()/pnum)
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.parallelsum(self.z2.data()[:n])/pnum)   

  def getxbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.xbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.xbar.data()[:n])/pnum)   

  def getybar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.ybar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.ybar.data()[:n])/pnum)   

  def getzbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.zbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.zbar.data()[:n])/pnum)   

  def getxpbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.xpbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.xpbar.data()[:n])/pnum)   

  def getypbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.ypbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.ypbar.data()[:n])/pnum)   

  def getxprms(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.xp2.data()/pnum)
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.parallelsum(self.xp2.data()[:n])/pnum)   

  def getyprms(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.yp2.data()/pnum)
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.parallelsum(self.yp2.data()[:n])/pnum)   

  def getxxpbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.xxpbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.xxpbar.data()[:n])/pnum)   

  def getyypbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.yypbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.yypbar.data()[:n])/pnum)   

  def getxpnbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.xpnbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.xpnbar.data()[:n])/pnum)   

  def getypnbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.ypnbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.ypnbar.data()[:n])/pnum)   

  def getxpnrms(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.xpn2.data()/pnum)
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.parallelsum(self.xpn2.data()[:n])/pnum)   

  def getypnrms(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.ypn2.data()/pnum)
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return sqrt(self.parallelsum(self.ypn2.data()[:n])/pnum)   

  def getxxpnbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.xxpnbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.xxpnbar.data()[:n])/pnum)   

  def getyypnbar(self):
    if not lparallel:
      pnum = self.pnum.data()
      pnum = where(pnum==0.,1.,pnum)
      return self.yypnbar.data()/pnum
    else:
      pnum = self.pnum.data()
      n = self.parallelmin(len(pnum))
      pnum = self.parallelsum(pnum[:n])
      pnum = where(pnum==0.,1.,pnum)
      return (self.parallelsum(self.yypnbar.data()[:n])/pnum)   

  def getemitxrms(self):
    xbar = self.getxbar()
    xpbar = self.getxpbar()
    if not lparallel:
      pnum = where(self.pnum.data()==0.,1.,self.pnum.data())
      x2  = self.x2.data()/pnum
      xp2 = self.xp2.data()/pnum
      xxp = self.xxpbar.data()/pnum
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnum = self.parallelsum(self.pnum.data()[:n])
      pnum = where(pnum==0.,1.,pnum)
      x2  = self.parallelsum(self.x2.data()[:n])/pnum
      xp2 = self.parallelsum(self.xp2.data()[:n])/pnum
      xxp = self.parallelsum(self.xxpbar.data()[:n])/pnum
    return sqrt((x2-xbar*xbar)*(xp2-xpbar*xpbar)-(xxp-xbar*xpbar)**2)
        
  def getemityrms(self):
    ybar = self.getybar()
    ypbar = self.getypbar()
    if not lparallel:
      pnum = where(self.pnum.data()==0.,1.,self.pnum.data())
      y2  = self.y2.data()/pnum
      yp2 = self.yp2.data()/pnum
      yyp = self.yypbar.data()/pnum
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnum = self.parallelsum(self.pnum.data()[:n])
      pnum = where(pnum==0.,1.,pnum)
      y2  = self.parallelsum(self.y2.data()[:n])/pnum
      yp2 = self.parallelsum(self.yp2.data()[:n])/pnum
      yyp = self.parallelsum(self.yypbar.data()[:n])/pnum
    return sqrt((y2-ybar*ybar)*(yp2-ypbar*ypbar)-(yyp-ybar*ypbar)**2)
        
  def getemitxnrms(self):
    xbar = self.getxbar()
    xpnbar = self.getxpnbar()
    if not lparallel:
      pnum = where(self.pnum.data()==0.,1.,self.pnum.data())
      x2  = self.x2.data()/pnum
      xpn2 = self.xpn2.data()/pnum
      xxpn = self.xxpnbar.data()/pnum
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnum = self.parallelsum(self.pnum.data()[:n])
      pnum = where(pnum==0.,1.,pnum)
      x2  = self.parallelsum(self.x2.data()[:n])/pnum
      xpn2 = self.parallelsum(self.xpn2.data()[:n])/pnum
      xxpn = self.parallelsum(self.xxpnbar.data()[:n])/pnum
    return sqrt((x2-xbar*xbar)*(xpn2-xpnbar*xpnbar)-(xxpn-xbar*xpnbar)**2)
        
  def getemitynrms(self):
    ybar = self.getybar()
    ypnbar = self.getypnbar()
    if not lparallel:
      pnum = where(self.pnum.data()==0.,1.,self.pnum.data())
      y2  = self.y2.data()/pnum
      ypn2 = self.ypn2.data()/pnum
      yypn = self.yypnbar.data()/pnum
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnum = self.parallelsum(self.pnum.data()[:n])
      pnum = where(pnum==0.,1.,pnum)
      y2  = self.parallelsum(self.y2.data()[:n])/pnum
      ypn2 = self.parallelsum(self.ypn2.data()[:n])/pnum
      yypn = self.parallelsum(self.yypnbar.data()[:n])/pnum
    return sqrt((y2-ybar*ybar)*(ypn2-ypnbar*ypnbar)-(yypn-ybar*ypnbar)**2)
        
  def gatherarray(self,a,zerotoone=0):
      # --- sum tparallel arrays from procs in one bucket
      if self.ntgroups>1:
        b = parallelsum(a,comm=self.mympitgroup)
      else:
        b = a.copy()
      if self.mympigroup is not None and self.gnpes>1:
        # --- update boundary
        if self.gme>0:
          mpisend(data = b[0,:], dest = self.gme-1, comm = self.mympigroup)
          b = b[1:,:]
        if self.gme<self.gnpes-1:
          b[-1,:] += mpirecv(source = self.gme+1, comm = self.mympigroup)
        # --- gather arrays from procs in one bucket
        b = gatherarray(b,bcast=0,comm=self.mympigroup)
      if zerotoone:
        b = where(b==0.,1.,b)
      return b
      
  def getpnumzhist(self):
      n = self.parallelmin(len(self.pnum.data()))
      return self.gatherarray(transpose(self.pnumz.data()[:n,:]))  

  def getxbarzhist(self):
    if not lparallel:
      pnumz = transpose(self.pnumz.data())
      pnumz = where(pnumz==0.,1.,pnumz)
      return transpose(self.xbarz.data())/pnumz
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      return (self.gatherarray(transpose(self.xbarz.data()[:n,:]))/self.gatherarray(pnumz,1))   

  def getybarzhist(self):
    if not lparallel:
      pnumz = transpose(self.pnumz.data())
      pnumz = where(pnumz==0.,1.,pnumz)
      return transpose(self.ybarz.data())/pnumz
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      return (self.gatherarray(transpose(self.ybarz.data()[:n,:]))/self.gatherarray(pnumz,1))   

  def getzbarzhist(self):
    if not lparallel:
      pnumz = transpose(self.pnumz.data())
      pnumz = where(pnumz==0.,1.,pnumz)
      return transpose(self.zbarz.data())/pnumz
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      return (self.gatherarray(transpose(self.zbarz.data()[:n,:]))/self.gatherarray(pnumz,1))   

  def getxpbarzhist(self):
    if not lparallel:
      pnumz = transpose(self.pnumz.data())
      pnumz = where(pnumz==0.,1.,pnumz)
      return transpose(self.xpbarz.data())/pnumz
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      return (self.gatherarray(transpose(self.xpbarz.data()[:n,:]))/self.gatherarray(pnumz,1))   

  def getypbarzhist(self):
    if not lparallel:
      pnumz = transpose(self.pnumz.data())
      pnumz = where(pnumz==0.,1.,pnumz)
      return transpose(self.ypbarz.data())/pnumz
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      return (self.gatherarray(transpose(self.ypbarz.data()[:n,:]))/self.gatherarray(pnumz,1))   

  def getxrmszhist(self):
    if not lparallel:
      pnumz = self.pnumz.data()
      pnumz = where(pnumz==0.,1.,pnumz)
      return transpose(sqrt(self.x2z.data()/pnumz))
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      return sqrt(self.gatherarray(transpose(self.x2z.data()[:n,:]))/self.gatherarray(pnumz,1))   

  def getyrmszhist(self):
    if not lparallel:
      pnumz = self.pnumz.data()
      pnumz = where(pnumz==0.,1.,pnumz)
      return transpose(sqrt(self.y2z.data()/pnumz))
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      return sqrt(self.gatherarray(transpose(self.y2z.data()[:n,:]))/self.gatherarray(pnumz,1))   

  def getzrmszhist(self):
    if not lparallel:
      pnumz = transpose(self.pnumz.data())
      pnumz = where(pnumz==0.,1.,pnumz)
      return transpose(sqrt(self.z2z.data()/pnumz))
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      return sqrt(self.gatherarray(transpose(self.z2z.data()[:n,:]))/self.gatherarray(pnumz,1))   

  def getemitxrmszhist(self):
    xbarz = self.getxbarzhist()
    xpbarz = self.getxpbarzhist()
    if not lparallel:
      pnumz = transpose(self.pnumz.data())
      pnumz = where(pnumz==0.,1.,pnumz)
      x2z  = transpose(self.x2z.data())/pnumz
      xp2z = transpose(self.xp2z.data())/pnumz
      xxpz = transpose(self.xxpbarz.data())/pnumz
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      pnumz = self.gatherarray(pnumz,1)
      x2z  = self.gatherarray(transpose(self.x2z.data()[:n,:]))/pnumz
      xp2z = self.gatherarray(transpose(self.xp2z.data()[:n,:]))/pnumz
      xxpz = self.gatherarray(transpose(self.xxpbarz.data()[:n,:]))/pnumz
    return sqrt(abs((x2z-xbarz*xbarz)*(xp2z-xpbarz*xpbarz)-(xxpz-xbarz*xpbarz)**2))

  def getemityrmszhist(self):
    ybarz = self.getybarzhist()
    ypbarz = self.getypbarzhist()
    if not lparallel:
      pnumz = transpose(self.pnumz.data())
      pnumz = where(pnumz==0.,1.,pnumz)
      y2z  = transpose(self.y2z.data())/pnumz
      yp2z = transpose(self.yp2z.data())/pnumz
      yypz = transpose(self.yypbarz.data())/pnumz
    else:
      n = self.parallelmin(len(self.pnum.data()))
      pnumz = transpose(self.pnumz.data()[:n,:])
      pnumz = self.gatherarray(pnumz,1)
      y2z  = self.gatherarray(transpose(self.y2z.data()[:n,:]))/pnumz
      yp2z = self.gatherarray(transpose(self.yp2z.data()[:n,:]))/pnumz
      yypz = self.gatherarray(transpose(self.yypbarz.data()[:n,:]))/pnumz
    return sqrt(abs((y2z-ybarz*ybarz)*(yp2z-ypbarz*ypbarz)-(yypz-ybarz*ypbarz)**2))

  def getzdedens(self):
    if not lparallel:
      return transpose(self.zdedens.data(),[1,2,0])
    else:
      n = self.parallelmin(len(self.pnum.data()))
      return self.gatherarray(transpose(self.zdedens.data()[:n,...],[1,2,0]))

  def getzdedipole(self):
    if not lparallel:
      zdedens = transpose(self.zdedens.data(),[1,2,0])
      zdedens = where(zdedens==0.,1.,zdedens)
      return transpose(self.zdedipole.data())/zdedens
    else:
      n = self.parallelmin(len(self.pnum.data()))
      zdedens = transpose(self.zdedens.data()[:n,...],[1,2,0])
      return (self.gatherarray(transpose(self.zdedipole.data()[:n,:],[1,2,0]))/self.gatherarray(zdedens,1))   

  def getrhoe(self,i=0):
    if me==0:
      rhoe=transpose(self.rhoe[i])
    else:
      rhoe=transpose(self.rhoe[i][:,:,1:])
    return transpose(gatherarray(rhoe))
      
  def getrhoi(self):
    if me==0:
      rhoi=transpose(self.rhoi[i])
    else:
      rhoi=transpose(self.rhoi[i][:,:,1:])
    return transpose(gatherarray(rhoi))
      
  def getphie(self):
    if me==0:
      phie=transpose(self.phie[i])
    else:
      phie=transpose(self.phie[i][:,:,1:])
    return transpose(gatherarray(phie))
      
  def getphii(self):
    if me==0:
      phii=transpose(self.phii[i])
    else:
      phii=transpose(self.phii[i][:,:,1:])
    return transpose(gatherarray(phii))
      
  def plfrac(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qspnum = self.getpnum()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qspnum/qspnum[0],xscale*qst,color=color,width=width)

  def plxrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsxrms = self.getxrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsxrms,xscale*qst,color=color,width=width)

  def plyrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsyrms = self.getyrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsyrms,xscale*qst,color=color,width=width)

  def plzrms(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qszrms = self.getzrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qszrms,xscale*qst,color=color,width=width)

  def plxbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsxbar = self.getxbar()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsxbar,xscale*qst,color=color,width=width)

  def plybar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsybar = self.getybar()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsybar,xscale*qst,color=color,width=width)

  def plzbar(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qszbar = self.getzbar()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qszbar,xscale*qst,color=color,width=width)

  def plemitx(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsexrms = self.getemitxrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsexrms,xscale*qst,color=color,width=width)

  def plemity(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qseyrms = self.getemityrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qseyrms,xscale*qst,color=color,width=width)

  def plemitxn(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qsexrms = self.getemitxnrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qsexrms,xscale*qst,color=color,width=width)

  def plemityn(self,color=black,width=1,type='solid',xscale=1.,yscale=1.):  
    qseyrms = self.getemitynrms()
    if me==0:
      qst = self.timemmnts.data()
      pla(yscale*qseyrms,xscale*qst,color=color,width=width)

class Quasistaticold:
  def __init__(self,slist=None,MRroot=None,beam=None,l_verbose=0,l_findmgparam=0,
               nparpgrp=top.nparpgrp,Ninit=1000,l_mode=1,l_selfe=1,maps=None,pboundxy=None,
               conductors=[],l_elecuniform=0,scraper=None,l_weakstrong=0,nelecperiod=1):
    assert beam is not None, 'beam must be passed as argument of Quasistatic'
    self.beam=beam
    if w3d.solvergeom==w3d.RZgeom:
      self.l_rz=1
      self.gridr=GRIDtype()
      gridcp=frz.basegrid
      frz.basegrid=self.gridr
      w3d.solvergeom=w3d.Rgeom
      nguardz=frz.nguardz+0
      init_base(w3d.nx,0,w3d.dx,0.,0.,0.,false)
      frz.basegrid=gridcp
      w3d.solvergeom=w3d.RZgeom
      frz.nguardz=nguardz
      frz.init_base(w3d.nx,w3d.nz,w3d.dx,w3d.dz,w3d.xmmin,w3d.zmminlocal,false)
      self.l_MR=0
    else:
      self.l_rz=0
      solvergeomcp=w3d.solvergeom
      w3d.solvergeom=w3d.XYgeom
      frz.init_base(w3d.nx,w3d.ny,w3d.dx,w3d.dy,w3d.xmmin,w3d.ymmin,false)
      if MRroot is None:
        self.l_MR = 0
        frz.mgridrz_accuracy = f3d.mgtol
      else:
        self.l_MR = 1
        frz.mgridrz_accuracy = MRroot.mgtol
        self.MRroot = MRroot
        b=MRroot
        brz=frz.basegrid
        l_addblock=true
        while l_addblock:
          try:
            b=b.children[0]
            add_subgrid(brz.gid[0],
                        b.nx,#-2*b.nguard*b.refinement[0],
                        b.ny,#-2*b.nguard*b.refinement[1],
                        b.dx,b.dy,
                        b.xmmin,#+b.nguard*b.dx*b.refinement[0],
                        b.ymmin,#+b.nguard*b.dy*b.refinement[1],
                        b.nguard*b.refinement[0],
                        b.nguard*b.refinement[0],
                        b.nguard*b.refinement[1],
                        b.nguard*b.refinement[1])
            brz=brz.down
          except:
            l_addblock=false
        g=frz.basegrid
        for conductor in conductors:
          if l_addblock:raise Exception('conductors need to be installed in quasistatic MR')
          try:
            cond = conductor.cond
          except AttributeError:
            cond = conductor
          self.conductors=ConductorType()
          installconductors(cond,conductors=self.conductors,gridrz=g,nz=0)
      w3d.solvergeom=solvergeomcp
    self.MRroot = MRroot
    self.dt = w3d.dz/top.vbeam # time step for pushing electrons
    self.nparpgrp = nparpgrp
    self.l_verbose = l_verbose
    self.slist = slist
    self.l_findmgparam = l_findmgparam and not self.l_rz
    self.Ninit = Ninit
    self.l_mode = l_mode
    if npes>1:
      self.l_mode=1
    else:
      self.l_mode = l_mode
    self.l_selfe = l_selfe
    self.l_elecuniform=l_elecuniform
    self.l_weakstrong=l_weakstrong
    self.nelecperiod=nelecperiod
    self.scraper=scraper
    if maps is not None:
      self.l_maps=true
      self.maps=maps
    if pboundxy is None:
      self.pboundxy = top.pboundxy
    else:
      self.pboundxy = pboundxy
      
  def step(self,nt=1,l_plotelec=0):
   # --- if running in parallel, enforcing l_mode=1
   if npes>1:self.l_mode=1
   for it in range(nt):
     # --- call beforestep functions
     callbeforestepfuncs.callfuncsinlist()
     if (not self.l_weakstrong or top.it==0) and (top.it%self.nelecperiod==0):
       loadrho()
       fieldsol()
       # --- switch to 2-D solver
       self.solver_3d_to_2d()
       # --- generate electrons (on last processor only)
       if me==max(0,npes-1) and (top.it==0  or not l_plotelec):self.create_electrons()
       # --- push 2-D electrons slice along bunch
       self.push_electrons(l_plotelec=l_plotelec)
       # --- switch back to 3-D solver
       self.solver_2d_to_3d()
     # --- push ions
     self.push_ions()
     # --- update time
     top.time+=top.dt
     # --- compute moments
     if self.l_verbose: print me,top.it,self.iz,'me = ',me,';compute zmmnt'
     if top.it%top.nhist==0:
       zmmnt()
       minidiag(top.it,top.time,top.lspecial)
     # --- call afterstep functions
     callafterstepfuncs.callfuncsinlist()
     # --- update time counter
     top.it+=1

  def push_ions(self):
    if me>=(npes-(top.it+1)*self.ntgroups):
#    if me==(npes-(top.it+1)*self.ntgroups):
     if self.l_maps:
      if self.l_verbose: print me,top.it,self.iz,'me = ',me,';apply space_charge kick'
      self.maps.apply_space_charge_kick(self.beam)
      if self.l_verbose: print me,top.it,self.iz,'me = ',me,';apply transfer map'
      self.maps.apply_transfer_map(self.beam)
      if self.l_verbose: print me,top.it,self.iz,'me = ',me,';apply bnd conditions'
      self.maps.apply_bnd_conditions(self.beam)

  def save_ions(self,zmax):
    # -- save data for ions with z<zmax
    # --- set shortcuts
    pg = top.pgroup
    js = 0
    il = pg.ins[js]-1
    iu = il+pg.nps[js]
    self.isaved = compress(pg.zp[il:iu]<zmax,il+arange(pg.nps[js]))
    self.xsaved = take(pg.xp[il:iu],self.isaved)
    self.ysaved = take(pg.yp[il:iu],self.isaved)
    self.zsaved = take(pg.zp[il:iu],self.isaved)
    self.uxsaved = take(pg.uxp[il:iu],self.isaved)
    self.uysaved = take(pg.uyp[il:iu],self.isaved)
    self.uzsaved = take(pg.uzp[il:iu],self.isaved)
    self.gisaved = take(pg.gaminv[il:iu],self.isaved)
    if pg.npid>0:self.pidsaved = take(pg.pid[il:iu,:],self.isaved,0)
    
  def push_electrons(self,l_return_dist=false,l_plotelec=0):
    # --- set shortcuts
    pg = top.pgroup
    bg = frz.basegrid

    # --- (diagnostic) stacking of electron distribution
    if l_return_dist:
      self.dist=[]
      sp=self.slist[0]
      self.dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])

    # --- sets z range
    if self.l_mode==2:
      izmin = w3d.nz/4
      izmax = w3d.nz-izmin
    else:
      izmin = 0
      izmax = w3d.nzp

    # --- loop over 3-D grid mesh in z
    for iz in range(izmax,izmin-1,-1):
      self.iz=iz
#      if self solve 2d field
      self.solve2dfield()
      # --- switch 2d and 3d slice of potentials
      self.switch2d3d(iz)          
      if l_plotelec and iz<w3d.nz/max(1,npes):self.plot_electrons()
      # --- push particles but on last step
      if(iz>0):self.push_electrons_2d()
      # --- (diagnostic) stacking of electron distribution
      if l_return_dist:
        sp=self.slist[0]
        self.dist.append([sp.getx(),sp.gety(),sp.getvx(),sp.getvy()])

    # --- clear electrons on processor 0, shift them to the previous ones for me>0
    self.clear_electrons()

    # --- fills end of 3-D phi arrays
    if self.l_mode==1:
      self.fillphizends(1)
    else:
      self.fillphizends(izmin)

#    if self.MRroot is not None:
#      self.setidosolve(self.MRroot,0)
##      self.MRroot.getallpotentialpforparticles()    
##      fieldsol()
#      self.MRroot.solve()
#      self.setidosolve(self.MRroot,1)      
      
  def clear_electrons(self):
    pg=top.pgroup
    # --- clear electrons on processor 0, shift them to the previous ones for me>0
    if me==0:
      pg.nps[1]=0      
    if npes>1:
      js = 1
      if pg.nps[js]>0:
        il=top.pgroup.ins[js]-1
        iu=il+top.pgroup.nps[js]
        top.pgroup.zp[il:iu]=w3d.zmminlocal-1.e-10*w3d.dz
      zparticleboundaries(top.pgroup,js,js,w3d.zmmax,w3d.zmmin,true)

  def plot_electrons(self):
#        ppgeneric(getvx(js=1),getx(js=1))
        window(0);fma();
#        ppco(getx(js=1),gety(js=1),getex(js=1),msize=2,ncolor=100);refresh()
#        if iz==0:
        self.slist[0].ppxy(msize=2);refresh()
#        else:
#          self.slist[0].ppxex(msize=2);
#        limits(w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
#        ppxex(js=1)
#        pfxy(iz=iz);
#        fma()

  def create_electrons(self):
     # --- define shortcuts
     pg = top.pgroup
     bg = frz.basegrid
     # --- generate electrons (on last processor only)
     pg.nps[1]=0
     if self.l_mode==2:
      zmax = w3d.zmmax/2+top.zgrid
     else:
      zmax = w3d.zmmax-1.e-10*w3d.dz+top.zgrid
     if top.prwall<sqrt(w3d.xmmax**2+w3d.ymmax**2) or self.l_rz:
      self.slist[0].add_uniform_cylinder(self.Ninit,min(top.prwall,w3d.xmmax),zmax,zmax,1.e-10,1.e-10,1.e-10,lallindomain=true)
     else:
      if self.l_elecuniform:
        spacing='uniform'
      else:
        spacing='random'
      self.slist[0].add_uniform_box(self.Ninit,w3d.xmmin,w3d.xmmax,
                                    w3d.ymmin,w3d.ymmax,zmax,zmax,
                                    1.e-10,1.e-10,1.e-10,spacing=spacing,lallindomain=true)
     # --- scrape electrons 
     if self.scraper is not None:self.scraper.scrapeall(local=1,clear=1)
#     shrinkpart(pg)

  def solve2dfield(self,l_findmgparam=false):
      if self.l_verbose:print me,top.it,self.iz,'me = ',me,'; enter fieldsolve'   # MR OK
      # --- define shortcuts
      pg = top.pgroup
      bg = frz.basegrid
      # --- reset 2-D rho arrays
      if self.l_rz:
        self.gridr.rho[...]=0.
      else:
        reset_rzmgrid_rho()
      # --- loop over species list
      for sp in self.slist:
        # --- loop over species index
        for js in sp.jslist:
          np = pg.nps[js]
          if np==0:continue
          il = pg.ins[js]-1
          iu = il+pg.nps[js]
          if self.l_verbose:print me,top.it,self.iz,'me = ',me,';deposit rho' # MR OK
          # --- deposit rho
          if self.l_rz:
            rhoweightr(pg.xp[il:iu],pg.yp[il:iu],np,pg.sq[js]*pg.sw[js],
                       bg.nr,bg.dr,bg.rmin)
          else:
            rhoweightrz(pg.xp[il:iu],pg.yp[il:iu],pg.yp[il:iu],np,
                        pg.sq[js]*pg.sw[js],bg.nr,bg.nz,
                        bg.dr,bg.dz,bg.rmin,0.)
      # --- distribute rho among patches
      distribute_rho_rz()
      
      # --- optimize mgparam
      if self.l_findmgparam and top.it==0 and self.iz==w3d.nz/2:find_mgparam_rz(false)

      # --- call 2-D field solver
      solve_mgridrz(bg,frz.mgridrz_accuracy,true)
      if self.l_verbose:print me,top.it,self.iz,'me = ',me,'; exit fieldsolve'   # MR OK

  def solver_3d_to_2d(self):
    # --- switch to 2-D solver
    if self.l_rz:
      self.gridcp=frz.basegrid
      frz.basegrid=self.gridr
      bg=self.gridr
      self.solvergeomcp = w3d.solvergeom
      self.fstypecp = top.fstype
      w3d.solvergeom = w3d.Rgeom
    else:
      self.solvergeomcp = w3d.solvergeom
      self.fstypecp = top.fstype
      w3d.solvergeom = w3d.XYgeom

  def solver_2d_to_3d(self):
    # --- switch to 3-D solver
    w3d.solvergeom = self.solvergeomcp
    top.fstype = self.fstypecp
    if self.l_rz:
      frz.basegrid=self.gridcp
      bg=frz.basegrid

  def switch2d3d(self,iz):
      # --- define shortcuts
      pg = top.pgroup
      bg = frz.basegrid
      # --- stack 2-D phi from electrons into 3-D array
      # --- and set 2-D phi from precalculated 3-D potential
      if self.l_MR:
        b3d = self.MRroot
        brz = bg
        for ig in range(frz.ngrids):
          if ig>0:
            b3d = b3d.children[0]
            brz = brz.down
            enzl=b3d.extradimslower[2]
            enzu=b3d.extradimsupper[2]
          else:
            enzl=0
            enzu=0
          if self.l_maps:
            phitmp = brz.phi[:,:].copy()
            if self.l_selfe:
#              brz.phi[1:-1,1:-1] += b3d.phi[...,iz+1+enzl] 
              brz.phi[1:-1,1:-1] += b3d.potentialarray[:,:,iz+1+enzl,0,0]
            else:
              brz.phi[:,:] = b3d.phi[...,iz+1+enzl].copy()
#              brz.phi[1:-1,1:-1] = b3d.phi[...,iz+1+enzl].copy()
            b3d.potentialparray[:,:,iz+1+enzl,0,0] = phitmp.copy()
#            b3d.phi[:,:,iz+1+enzl] = phitmp.copy()
            del phitmp
          else:
            if self.l_selfe:
#              brz.phi[1:-1,1:-1] += b3d.phi[...,iz+1+enzl] 
#              b3d.phi[:,:,iz+1+enzl] = brz.phi[1:-1,1:-1]
              b3d.potentialparray[:,:,iz+1+enzl,0,1] = brz.phi[1:-1,1:-1]
              brz.phi[1:-1,1:-1] += b3d.potentialparray[:,:,iz+1+enzl,0,0]
              b3d.potentialparray[:,:,iz+1+enzl,0,0]=0. # would be added again by getallpotentialpforparticles() otherwise
            else:
              phitmp = brz.phi[1:-1,1:-1].copy()
              brz.phi[1:-1,1:-1] = b3d.phi[...,iz+1+enzl].copy()
              b3d.phi[:,:,iz+1+enzl] += phitmp.copy()
      else:
       if self.l_rz:
        phitmp = bg.phi.copy()
        if self.l_selfe:
          bg.phi[:,0] += gridcp.phi[...,iz+1] 
        else:
          bg.phi[:,0] = gridcp.phi[...,iz+1].copy() 
        gridcp.phi[...,iz+1] = phitmp[:,0].copy()
       else:
        if self.l_maps:
          phitmp = bg.phi[1:-1,1:-1].copy()
          if self.l_selfe:
            bg.phi[1:-1,1:-1] += w3d.phi[:,:,iz+1] 
          else:
            bg.phi[1:-1,1:-1] = w3d.phi[:,:,iz+1].copy() 
          w3d.phi[:,:,iz+1] = phitmp.copy()
        else:
          raise Exception('this part needs to be programmed: Quasistatic+Boris-l_MR')
      # --- update guard cells of 2-D solver
      updateguardcells2d()
      if frz.l_get_fields_on_grid:getallfieldsfromphip()

  def push_electrons_2d(self):
       # --- define shortcuts
       pg = top.pgroup
       bg = frz.basegrid
       # --- loop over species list
       for sp in self.slist:
        # --- loop over species index
        for js in sp.jslist:
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if np==0:continue
            # --- gather self-forces
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';fieldweight' # MR OK
            pg.ex[il:iu]=0.
            pg.ey[il:iu]=0.
            pg.ez[il:iu]=0.
            if self.l_rz:
              fieldweightr(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np)
            else:
              fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,0.,top.efetch[js])
#            fetche3d(pg,il+1,np,js+1)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';xpush'
            xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                    pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],self.dt)
            pg.zp[il:iu]-=w3d.dz
#            pli(frz.basegrid.phi);refresh()
#            ppzx(js=1,color=red,msize=2);refresh()
#            window(3);ppzvx(js=1,msize=2);window(0)
          if self.l_verbose:print me,top.it,self.iz,'me = ',me,';stckxy3d'
          stckxy3d(top.pgroup,js,top.zbeam,true)
          # --- This may need to call partbndwithdata since stckxy3d does not
          # --- check against the grid anymore!
          if self.scraper is not None:
            self.scraper.scrapeall(local=1,clear=1)
          else:
            processlostpart(top.pgroup,js+1,top.clearlostpart,top.time+self.dt*top.pgroup.ndts[js],top.zbeam)
          top.npslost=0
          ng = 1+pg.nps[js]/self.nparpgrp
          for ig in range(ng):
            il = pg.ins[js]-1+self.nparpgrp*ig
            iu = min(il+self.nparpgrp,pg.ins[js]-1+pg.nps[js])
            np = iu-il
            if np==0:continue
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';fieldweight' # MR OK
            pg.ex[il:iu]=0.
            pg.ey[il:iu]=0.
            pg.ez[il:iu]=0.
            if self.l_rz:
              fieldweightr(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np)
            else:
              fieldweightxz(pg.xp[il:iu],pg.yp[il:iu],pg.ex[il:iu],pg.ey[il:iu],np,top.zgrid,top.efetch[js])
#            fetche3d(pg,il+1,np,js+1)
            if self.l_verbose:print me,top.it,self.iz,'me = ',me,';epush'
            epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                    pg.ex[il:iu],pg.ey[il:iu],pg.ez[il:iu],pg.sq[js],pg.sm[js],0.5*self.dt)
            gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     top.gamadv,top.lrelativ)

  def setefetch(self,fsolver,efetch):
        fsolver.efetch=efetch
        try:
          self.setefetch(fsolver.children[0],efetch)
        except:
          return
 
  def setidosolve(self,fsolver,idosolve):
        fsolver.l_internal_dosolve=idosolve
        try:
          self.setidosolve(fsolver.children[0],idosolve)
        except:
          return
 
  def fillphizends(self,nz):
      if self.l_MR:
        b3d = self.MRroot
        for ig in range(frz.ngrids):
          if ig>0:
            b3d = b3d.children[0]
            enzl=b3d.extradimslower[2]
            enzu=b3d.extradimsupper[2]+w3d.nzlocal-w3d.nzp
          else:
            enzl=0
            enzu=w3d.nzlocal-w3d.nzp
#            enzu=1
          for iz in range(nz+enzl):
#            b3d.phi[:,:,iz] = b3d.phi[:,:,nz+enzl]
            b3d.potentialparray[:,:,iz,0,0] = b3d.potentialparray[:,:,nz+enzl,0,0]
          for iz in range(nz+enzu):
#            b3d.phi[:,:,-iz-1] = b3d.phi[:,:,-nz-1-enzu]
            b3d.potentialparray[:,:,-iz-1,0,0] = b3d.potentialparray[:,:,-nz-1-enzu,0,0]
      else:
        if self.l_rz:
          for iz in range(nz):
            bg.phi[:,iz] = bg.phi[:,nz]
            bg.phi[:,-iz-1] = bg.phi[:,-nz-1]
        else:
          for iz in range(nz):
            w3d.phi[:,:,iz] = w3d.phi[:,:,nz]
            w3d.phi[:,:,-iz-1] = w3d.phi[:,:,-nz-1]

# --- This can only be done after the class is defined.
try:
  psyco.bind(Quasistatic)
except NameError:
  print 'Warning:psyco binding of Quasistatic class failed.'


