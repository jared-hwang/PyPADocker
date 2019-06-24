"""
ParticleScraper: class for creating particle scraping
"""
from ..warp import *
#import decorators


def particlescraperdoc():
    from ..particles import particlescraper
    print particlescraper.__doc__


class ParticleScraper(object):
    """
  Class for creating particle scraper for conductors
   - conductors: a conductor or list of conductors which act as particle scrapers
                 Note that each conductor MUST have a unique id that is > 0.
   - lsavecondid: when true, the id of the conductor where the particle is
                  lost is save. The id is saved in the array top.pidlost[:,-1].
   - lsaveintercept: when true, the location and surface normal where the
                     particle intercepted the conductor surface is calculated.
                     The location is overwritten onto the xplost, yplost, and
                     zplost arrays. The angles describing the surface normal are
                     put into pidlost[:,-3] and pidlost[:,-2]. The spherical
                     coordinate angle theta is in -3, and phi is in -2.
                     The time at which the particles are lost is put into
                     pidlost[:,-4].
                     Note that the condid where the particle is lost is also
                     saved in pidlost[:,-1].
   - lrefineallintercept: when true, with lsaveintercept, when determining if
                          particles are lost, particles near conductors are
                          advanced, starting from the old positions, using a
                          refined time step small compared to the cyclotron
                          gyroperiod.
                          This option is useful when using the Drift-Lorentz
                          mover, i.e. when interpdk is turned on, and when the
                          time step is larger than the gyroperiod. In that case,
                          the calculation of the intercept can be inaccurate
                          since the gyromotion is not resolved. The time-step
                          refinement resolves the gyromotion so the incident
                          angle of the particle onto the conductor can be
                          correctly obtained.
   - lrefineintercept: when true, with lsaveintercept, only lost particles are
                       advanced from the old positions using a time step small
                       compared to the cyclotron gyroperiod to calculate a
                       better value for the intercept. This option is not
                       recommended, being superceded by lrefineallintercept.
   - nstepsperorbit=8: number of refined time steps to take when using
                       lrefineallintercept or lrefineintercept.
   - lcollectlpdata: When true, the lost particles statistics will be collected
                     for each conductor in the list lostparticles_data (Assembly
                     class).
   - mglevel=0: Coarsening level for index grid which is used to determine
                which conductors particles are near. This grid is a full size,
                3d (or 2d) array and can require a not insignificant amount of
                time to compute. If it is expected that few particles will be
                lost, using a coarser grid can substantially reduce the memory
                and time requirements for this grid. However, a coarser grid
                will mean that more particles will be flagged as being near the
                conductors, and the detailed check required for these particles
                will be more expensive.  A trade off is required that can really
                only be optimized empirically. A value of 0, 1, or 2 are
                probably optimal.
   - install=1: flag whether or not to install the scraper so that the scraping
                automatically happens every time step.
   - lbeforescraper=0: when true, the scraper is installed before grid scraping,
                      otherwise it is installed after. This is not recommended
                      since the scraping algorithms implemented here will not
                      work for particles out of the grid bounds anyway.
                      This option is retain only for legacy purposes.
   - lfastscraper=0: A faster, but approximate, version of the scraper. Note
                     that with this option, which conductor the particle is
                     scraped on is not calculated (and so cannot be used with
                     lsavecondid nor lsaveintercept)
   - grid=None: A instance of the Grid class can be supplied, allowing control
                over the region where the scraping is done and the resolution
                of the scraping data.
   - gridmode=0: Flag which sets whether the grid is fixed or automatically
                 updated as the beam frame moves. When 0, the grid is
                 automatically updated, when 1, the grid is fixed.
                 Note that when gridmode=0 and the beam frame is moving,
                 there will be some overhead since the internal data related
                 to the conductor location will be updated every time step.
   - nxscale=1,nyscale=1,nzscale=1: Scale factor on the number of grid cells
                                    on the workig mesh. These should be integer
                                    values.
   - interceptvelocitymethod='finitedifference':
          This sets how the velocity is calculated when finding the intercept
          where lost particles entered a conductor. When 'finitedifference', use
          a finite difference of the new minus the old particle positions. This
          gives more robust results in the presence of magnetic fields.
          When 'actualvelocity', use the actual particle velocity. This option
          should only ever be used in special cases and for testing.
   - species: List of species that should be scraped. Defaults to all species.

  After an instance is created, additional conductors can be added by calling
  the method registerconductors which takes either a conductor or a list of
  conductors are an argument. Otherwise, nothing further needs to be done - the
  scraping will happen automatically.

  Various methods are accessible if additional fine control is needed.
    """
    #__metaclass__ = decorators.TimedClass
    def __init__(self,conductors=None,lsavecondid=0,lsaveintercept=0,
                      lrefineintercept=0,lrefineallintercept=0,nstepsperorbit=8,
                      lcollectlpdata=0,mglevel=0,aura=0.,
                      install=1,lbeforescraper=0,lfastscraper=0,
                      grid=None,gridmode=0,nxscale=1,nyscale=1,nzscale=1,
                      interceptvelocitymethod='finitedifference',
                      species=None):
        self.mglevel = mglevel
        self.aura = aura
        self.lbeforescraper = lbeforescraper
        self.lfastscraper = lfastscraper
        self.interceptvelocitymethod = interceptvelocitymethod
        self.species = species
        # --- First set so install is false. Reset later with input value.
        # --- This is needed since in some cases registerconductors may want
        # --- to do the install. This just skips it in that case.
        self.install = 0
        # --- Remember if the user specified the grid.
        self.usergrid = (grid is not None)
        # --- Don't create the grid until it is needed.
        self.grid = grid
        self.gridmode = gridmode
        self.nxscale = nxscale
        self.nyscale = nyscale
        self.nzscale = nzscale
        # --- By default, don't save the old positions or velocities.
        self.lsaveoldpositions = false
        self.lsaveoldvelocities = false
        # --- register any initial conductors
        self.conductors = []
        self.newconductors = false
        self.reflectiveconductors = false
        self.registerconductors(conductors)
        # --- Allocate arrays for lost particles
        gchange("LostParticles")
        # --- If the conductor id where particles are lost is being saved,
        # --- need to turn on saving of lost particles.
        self.lsaveintercept = lsaveintercept or lrefineintercept or lrefineallintercept
        self.lrefineintercept = lrefineintercept
        self.lrefineallintercept = lrefineallintercept
        self.nstepsperorbit = nstepsperorbit
        self.lsavecondid = (lsavecondid or lsaveintercept or
                            lrefineintercept or lrefineallintercept or
                            lcollectlpdata)
        assert not (self.lsavecondid and self.lfastscraper),"With the fast scraper, the conductor information where the particles are lost is not calculated and so lsavecondid will not work"
        self.lcollectlpdata = lcollectlpdata
        if self.lsavecondid:
            top.lsavelostpart = true
        if self.lsaveintercept:
            self.lsaveoldpositions = true
            if self.lrefineintercept or self.lrefineallintercept:
                self.lsaveoldvelocities = true
        self.l_print_timing=0
        # --- If the user specified the grid, then add the conductors
        if self.usergrid:
            # --- Make sure that the grid has some finite extent. If not, then just
            # --- return so that the scraper will not be operational. This can
            # --- in parallel since the user supplied grid may not extend over the
            # --- full system and so some processors will not overlap the grid.
            if ((grid.nxlocal == 0 and grid.nylocal == 0 and grid.nzlocal == 0) or
                (grid.nxlocal < 0 or grid.nylocal < 0 or grid.nzlocal < 0)):
                if not lparallel:
                    print "Grid: warning: the user supplied grid has zero extent."
                return
            self.updateconductors()
        # --- Install the call to scrape particles if requested
        self.install = install
        self.installscraper()
        # --- This is needed but is not necessarily the correct code.
        # --- xoldpid etc aren't defined until saveolddata is called, and it
        # --- isn't called until after scraping happens. But the xoldpid etc
        # --- are needed the first time scraping happens - this produces an
        # --- error. This is a conceptual error in general when particles
        # --- are being added. On their first time being scraped, the xold
        # --- etc have not yet been saved, but those quantities are needed.
        # --- Adding a call to saveolddata here partially fixes the first
        # --- problem. Fixing the other will be more complicated - perhaps
        # --- requiring a flag that says whether the xold has been saved
        # --- yet.
        self.saveolddata()

    def installscraper(self):
        """
    Install the scraper so that it is called during at the appropriate place
    in a time step. This is normally done automically."""
        if not self.install: return
        # --- Install the call to scrape particles
        if self.lbeforescraper:
            if not isinstalledbeforescraper(self.scrapeall):
                installbeforescraper(self.scrapeall)
        else:
            if not isinstalledparticlescraper(self.scrapeall):
                installparticlescraper(self.scrapeall)

    def disable(self):
        """Uninstall the scraper so that it will not be called.
        """
        if self.lbeforescraper:
            if isinstalledbeforescraper(self.scrapeall):
                uninstallbeforescraper(self.scrapeall)
        else:
            if isinstalledparticlescraper(self.scrapeall):
                uninstallparticlescraper(self.scrapeall)

    def __setstate__(self,dict):
        """This is called when the instance is unpickled."""
        self.__dict__.update(dict)
        self.installscraper()

        # --- This is needed so that new variables get values when restoring
        # --- from an older version.
        if 'reducedisinside' not in self.__dict__:
            #self.reducedisinside = self.grid.isinside.copy()
            try:
                self.reducedisinside = self.grid.isinside
            except AttributeError:
                pass
        if 'lrefineintercept' not in self.__dict__:
            self.lrefineintercept = 0
        if 'lrefineallintercept' not in self.__dict__:
            self.lrefineallintercept = 0

        if 'lfastscraper' not in self.__dict__:
            self.lfastscraper = 0
        if 'reflectiveconductors' not in self.__dict__:
            self.reflectiveconductors = false
            for c in self.conductors:
                if c.material == 'reflector':
                    # --- Set the flag signifying the presence of reflective conductors
                    self.reflectiveconductors = true

    def registerconductors(self,newconductors):
        """Adds conductors to the list of conductors that particles are scraped on.
        """
        if newconductors is None: return
        if not isinstance(newconductors,list): newconductors = [newconductors]
        for c in newconductors:
            assert c.condid > 0,"The conductor id must be greater than zero in order for the particle scraping to work."
            self.conductors.append(c)
            if c.material == 'reflector':
                # --- For reflector materials, the old position is saved so that when
                # --- particles reflect, their current position will be replaced
                # --- with the old.
                self.lsaveoldpositions = true
                #self.lsaveoldvelocities = true
                 # --- Set the flag signifying the presence of reflective conductors
                self.reflectiveconductors = true
        self.newconductors = true

    def unregisterconductors(self,conductor,nooverlap=0):
        """Remove the conductor from the list of conductors that particles are
    scraped on.
        """
        self.conductors.remove(conductor)
        if not nooverlap or self.lfastscraper:
            # --- This is horribly inefficient!!!
            #self.grid.resetisinside()
            #self.updateconductors()
            # --- This causes the conductor data to be regenerated for all conductors
            self.newconductors = true
        else:
            self.grid.removeisinside(conductor)

    def updategrid(self,lforce=0):
        """
    Update the grid to match any changes to the underlying grid, for example
    after load balancing. This also does the initialization of the grid.
    The code is structured this way so that the grid is only created when needed.
    This allows the particles scraper instance to be created before all of the
    data needed by the grid is set up. This will normally be called automatically.
        """
        if self.l_print_timing:tstart=wtime()
        if self.grid is None: lforce = 1

        gridchanged = 0
        if lforce or not self.usergrid:

            # --- Check if self.grid.decomp is defined. If not, then force
            # --- an update.
            try:
                self.grid.decomp
            except AttributeError:
                lforce = 1

            if not lforce:
                # --- Check if the solver's grid has changed. If not, then
                # --- return since nothing needs to be done.
                gdc = self.grid.decomp
                tdc = top.ppdecomp
                gridchanged = not (
                    self.grid.nxlocal == w3d.nxp*self.nxscale and
                    self.grid.nylocal == w3d.nyp*self.nyscale and
                    self.grid.nzlocal == w3d.nzp*self.nzscale and
                    self.grid.xmmin == w3d.xmmin and
                    self.grid.xmmax == w3d.xmmax and
                    self.grid.ymmin == w3d.ymmin and
                    self.grid.ymmax == w3d.ymmax and
                    self.grid.zmmin == w3d.zmmin and
                    self.grid.zmmax == w3d.zmmax and
                    gdc.ix[gdc.ixproc] == tdc.ix[tdc.ixproc]*self.nxscale and
                    gdc.iy[gdc.iyproc] == tdc.iy[tdc.iyproc]*self.nyscale and
                    gdc.iz[gdc.izproc] == tdc.iz[tdc.izproc]*self.nzscale and
                    gdc.nx[gdc.ixproc] == tdc.nx[tdc.ixproc]*self.nxscale and
                    gdc.ny[gdc.iyproc] == tdc.ny[tdc.iyproc]*self.nyscale and
                    gdc.nz[gdc.izproc] == tdc.nz[tdc.izproc]*self.nzscale and
                    (self.zbeamprev == top.zbeam or self.gridmode != 0))
            else:
                gridchanged = 1

            if gridchanged:
                # --- Note that a copy of the decomposition is passed in.
                # --- The decomposition in top may be changed the next time
                # --- loadbalancing is done, but the decomposition in self.grid
                # --- should not be changed. Instead, a whole new grid is created.
                self.grid = Grid(decomp=copy.deepcopy(top.ppdecomp),
                                 nxscale=self.nxscale,nyscale=self.nyscale,
                                 nzscale=self.nzscale)
                self.zbeamprev = self.grid.zbeam

        if self.newconductors or gridchanged:
            self.newconductors = false
            self.updateconductors()

            if top.chdtspid>0:
                if (w3d.nxc != self.grid.nxlocal or
                    w3d.nyc != self.grid.nylocal or
                    w3d.nzc != self.grid.nzlocal):
                    w3d.nxc=self.grid.nxlocal
                    w3d.nyc=self.grid.nylocal
                    w3d.nzc=self.grid.nzlocal
                    gchange('Fields3dParticles')
                    sum_neighbors3d(nint(self.grid.isinside),w3d.isnearbycond,
                                    w3d.nxc,w3d.nyc,w3d.nzc)

        if self.l_print_timing:
            tend=wtime()
            print 'updategrid',tend-tstart

    def updateconductors(self):
        """Generate the data on the grid from the scraping conductors. This will normally be called automatically."""
        if not self.lfastscraper:
            self.grid.resetisinside()
            for c in self.conductors:
                self.grid.getisinside(c,mglevel=self.mglevel,aura=self.aura)
            # --- reducedisinside is a copy of isinside but will be modified to remove
            # --- redundant information. This provides an optimization of the routines
            # --- which find intersections with conductors. Normally, a particle is
            # --- compared against the conductors that the grid point surrounding it
            # --- are in. If more than one of those grid points are in the same
            # --- conductor, the particle will be checked against that conductor
            # --- multiple times. This is a waste of CPU time. The reducing routine
            # --- checks if a grid point is between two grid points that are in the
            # --- same conductor as itself. If so, then the fact that the grid point
            # --- is inside that conductor can be ignored, since particles nearby
            # --- will get a reference to the conductor from the neighboring grid
            # --- points. Note that the routine never ignores grid points that have
            # --- nx,ny,nz all even.
            #self.reducedisinside = fzeros(self.grid.isinside.shape,'d')
            #self.reducedisinside[...] = self.grid.isinside
            #reduceisinsidegrid(self.grid.isinside,self.reducedisinside,
            #                   self.grid.nx,self.grid.ny,self.grid.nz)
            # --- There is a problem with the above so don't use for now
            # --- Just make a reference. Similarly in setstate.
            self.reducedisinside = self.grid.isinside
        else:
            for c in self.conductors:
                self.grid.getdistances(c)

    def jslist(self,jslist=None):
        "Return a list of species indices that are scraped"
        if jslist is not None:
            # --- If a value if given, just return it
            result = jslist
        elif self.species is None:
            # --- The default is all species
            result = range(top.pgroup.ns)
        else:
            # --- Gather the jslists from each specified species
            result = []
            for s in self.species:
                result += s.jslist
        return result

    def saveolddata(self):
        """Saves old particle data. In some cases, the location of the particle
    before it was lost is needed."""
        # --- If no data is to be saved, then do nothing.
        if not (self.lsaveoldpositions or self.lsaveoldvelocities): return

        # --- Check if the pid indices have been setup.
        if self.lsaveoldpositions and 'xoldpid' not in self.__dict__:
            # --- Note that nextpid returns numbers based on 1 based indexing
            self.xoldpid = nextpid() - 1
            self.yoldpid = nextpid() - 1
            self.zoldpid = nextpid() - 1
            self.oldisOK = nextpid() - 1
            setuppgroup(top.pgroup)
        if self.lsaveoldvelocities and 'uxoldpid' not in self.__dict__:
            self.uxoldpid = nextpid() - 1
            self.uyoldpid = nextpid() - 1
            self.uzoldpid = nextpid() - 1
            setuppgroup(top.pgroup)

        # --- Do the saving.
        for js in self.jslist():
            if top.pgroup.ldts[js] and top.pgroup.nps[js] > 0:
                i1 = top.pgroup.ins[js] - 1
                i2 = i1 + top.pgroup.nps[js]
                if self.lsaveoldpositions:
                    # --- The code could be written this way now since the get routines
                    # --- can now return a direct reference to the data.
                    #getpid(id=self.xoldpid,js=js,gather=0)[:] = getx(js=js,gather=0)
                    #getpid(id=self.yoldpid,js=js,gather=0)[:] = gety(js=js,gather=0)
                    #getpid(id=self.zoldpid,js=js,gather=0)[:] = getz(js=js,gather=0)
                    # --- But this is ~3 times faster for small numbers of particles
                    # --- due to overhead of the get functions.
                    top.pgroup.pid[i1:i2,self.xoldpid] = top.pgroup.xp[i1:i2]
                    top.pgroup.pid[i1:i2,self.yoldpid] = top.pgroup.yp[i1:i2]
                    top.pgroup.pid[i1:i2,self.zoldpid] = top.pgroup.zp[i1:i2]
                    top.pgroup.pid[i1:i2,self.oldisOK] = 1.
                if self.lsaveoldvelocities:
                    #getpid(id=self.uxoldpid,js=js,gather=0)[:] = getux(js=js,gather=0)
                    #getpid(id=self.uyoldpid,js=js,gather=0)[:] = getuy(js=js,gather=0)
                    #getpid(id=self.uzoldpid,js=js,gather=0)[:] = getuz(js=js,gather=0)
                    top.pgroup.pid[i1:i2,self.uxoldpid] = top.pgroup.uxp[i1:i2]
                    top.pgroup.pid[i1:i2,self.uyoldpid] = top.pgroup.uyp[i1:i2]
                    top.pgroup.pid[i1:i2,self.uzoldpid] = top.pgroup.uzp[i1:i2]

    def applysymmetry(self,xc,yc):
        """Apply symmetry conditions to the positions so that the data passed
    into isinside is consistent with that obtained from the grid.
        """
        if self.grid.l4symtry:
            xcsym = abs(xc)
            ycsym = abs(yc)
        elif self.grid.l2symtry:
            xcsym = xc
            ycsym = abs(yc)
        else:
            xcsym = xc
            ycsym = yc
        return xcsym,ycsym

    def scrapeall(self,clear=0,local=0,jslist=None):
        """
    Apply scraping to all of the species. This will normally be called automatically.
      - clear=0: when true, lost particles are removed from the particle arrays.
                 Note that this routine is normally called during the course of
                 a time step, where the removal of lost particles is handled
                 elsewhere.
      - local=0: This only affects the lcollectlpdata option.
                 When true, this collection of data is turned off
                 (avoiding a parallel operation).
      - jslist=None: Optional list of species to scrape. It defaults to all species.
        """
        if len(self.conductors)==0: return
        self.updategrid()
        for js in self.jslist(jslist):
            if top.pgroup.ldts[js]:
                if self.l_print_timing:tstart=wtime()
                if self.lfastscraper:
                    self.fastscrape(js);
                else:
                    self.scrape(js);
                if self.l_print_timing:tend=wtime()
                if self.l_print_timing:print js,'scrape',tend-tstart
                if self.l_print_timing:tstart=wtime()
                if clear or self.lsavecondid:
                    processlostpart(top.pgroup,js+1,top.clearlostpart,top.time,top.zbeam)
                if self.l_print_timing:tend=wtime()
                if self.l_print_timing:print js,'processlosspart',tend-tstart
                if self.l_print_timing:tstart=wtime()
                if self.lsavecondid:
                    self.savecondid(js,local=local)
                if self.l_print_timing:tend=wtime()
                if self.l_print_timing:print js,'savecondid',tend-tstart
        if self.reflectiveconductors:
            # --- If there are any reflecting conductors, then redo the parallel
            # --- boundary conditions, since any reflected particles will have
            # --- their position replaced by their old position, which may be in
            # --- a different domain.
            # --- Note that this is a global operation, so all processors must
            # --- make this call.
            particlegridboundaries3d(top.pgroup,-1)
        self.saveolddata()
    #scrapeall = decorators.timedmethod(scrapeall)

    def scrape(self,js):
        """Apply scraping to species js. It is better to call scrapeall. This will normally be called automatically."""
        # --- If there are no particles in this species, that nothing needs to be done
        if top.pgroup.nps[js] == 0: return

        # --- Get mesh information into local variables
        dx,dy,dz,nx,ny,nz,ix,iy,iz = self.grid.getmeshsize(self.mglevel)
        xmin = self.grid.xmmin + ix*dx
        xmax = self.grid.xmmin + (ix+nx)*dx
        ymin = self.grid.ymmin + iy*dy
        ymax = self.grid.ymmin + (iy+ny)*dy
        zmin = self.grid.zmmin + iz*dz + top.zbeam
        zmax = self.grid.zmmin + (iz+nz)*dz + top.zbeam

        # --- The ixa etc are the location of the x=0 plane. This is needed
        # --- since in certain cases, the size of the collapsed dimension
        # --- can be greater than one and the 0 plane needs to be found, it
        # --- is the plane where the data is stored. This is true for example
        # --- with the RZ EM solver.
        # --- The check of nx > 0 etc is done since in certain cases the xmin
        # --- can be nonzero when nx is zero, giving an erroneous value for ixa.
        # --- This is the case when the quasistatic solver is being used.
        ixa = iya = iza = 0
        if nx > 0: ixa = nint(-xmin/dx)
        if ny > 0: iya = nint(-ymin/dy)
        if nz > 0: iza = nint(-zmin/dz)
        isinside = self.grid.isinside

        # --- Get handy references to the particles in the species
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        xx = top.pgroup.xp[i1:i2]
        yy = top.pgroup.yp[i1:i2]
        zz = top.pgroup.zp[i1:i2]
        pp = zeros(top.pgroup.nps[js],'d')
        #if js==1:print js,i1,i2,top.pgroup.zp[i1:i2],top.zbeam

        # --- Find which particles are close to a conductor. This
        # --- interpolates from the isinside grid. The results are
        # --- put into the array pp.
        if w3d.solvergeom in [w3d.XYZgeom]:
            getgrid3d(top.pgroup.nps[js],xx,yy,zz,pp,
                      nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                      w3d.l2symtry,w3d.l4symtry)
        elif w3d.solvergeom == w3d.RZgeom:
            # --- Note that for RZ, the radius is calculated for this, but
            # --- the original particle position is used below.
            rr = sqrt(xx**2 + yy**2)
            getgrid2d(top.pgroup.nps[js],rr,zz,pp,nx,nz,isinside[:,iya,:],
                      xmin,xmax,zmin,zmax)
        elif w3d.solvergeom == w3d.XZgeom:
            xsym,ysym = self.applysymmetry(xx,0)
            getgrid2d(top.pgroup.nps[js],xsym,zz,pp,nx,nz,isinside[:,iya,:],
                      xmin,xmax,zmin,zmax)
        elif w3d.solvergeom == w3d.XYgeom:
            xsym,ysym = self.applysymmetry(xx,yy)
            getgrid2d(top.pgroup.nps[js],xsym,ysym,pp,nx,ny,isinside[:,:,iza],
                      xmin,xmax,ymin,ymax)
        elif w3d.solvergeom == w3d.Rgeom:
            # --- Note that for R, the radius is calculated for this, but
            # --- the original particle position is used below.
            rr = sqrt(xx**2 + yy**2)
            getgrid1d(top.pgroup.nps[js],rr,pp,nx,isinside[:,iya,iza],
                      xmin,xmax)
        elif w3d.solvergeom == w3d.Ygeom:
            xsym,ysym = self.applysymmetry(0,yy)
            getgrid1d(top.pgroup.nps[js],ysym,pp,ny,isinside[ixa,:,iza],
                      ymin,ymax)
        elif w3d.solvergeom == w3d.Zgeom:
            getgrid1d(top.pgroup.nps[js],zz,pp,nz,isinside[ixa,iya,:],
                      ymin,ymax)
        else:
            raise Exception("The particle scraping only works for XYZ, XY, RZ, R, Y and Z geometry")

        # --- Get indices for all of the particles which are close to a
        # --- conductor. If there are none, then immediately return.
        # --- Note, of course, that close may mean inside.
        iclose = compress(pp>0.,arange(i1,i2))
        if len(iclose) == 0: return

        # --- Get the positions of particles which are close to a conductor.
        xx = take(xx,iclose-i1)
        yy = take(yy,iclose-i1)
        zz = take(zz,iclose-i1)

        # --- The 'g' lists give the locations of the corners of the grid cell
        # --- relative to the grid location of the particles close to a
        # --- conductor. Also, get those grid locations.
        if w3d.solvergeom in [w3d.XYZgeom]:
            nd = 3
            gdx = [0.,dx,0.,dx,0.,dx,0.,dx]
            gdy = [0.,0.,dy,dy,0.,0.,dy,dy]
            gdz = [0.,0.,0.,0.,dz,dz,dz,dz]
            xg = xmin+aint(abs(xx-xmin)/dx)*dx
            yg = ymin+aint(abs(yy-ymin)/dy)*dy
            zg = zmin+aint(abs(zz-zmin)/dz)*dz
        elif w3d.solvergeom in [w3d.RZgeom]:
            nd = 2
            gdx = [0.,dx,0.,dx]
            gdz = [0.,0.,dz,dz]
            # --- Like above, the radius is calculated in the temporary, but the
            # --- original particle position is used below.
            # --- These two lines calculating rr give the same result, but the second
            # --- is probably faster
            #rr = sqrt(xx**2 + yy**2)
            rr = take(rr,iclose-i1)
            xg = xmin+aint(abs(rr-xmin)/dx)*dx
            zg = zmin+aint(abs(zz-zmin)/dz)*dz
        elif w3d.solvergeom in [w3d.XZgeom]:
            nd = 2
            gdx = [0.,dx,0.,dx]
            gdz = [0.,0.,dz,dz]
            xg = xmin+aint(abs(xx-xmin)/dx)*dx
            zg = zmin+aint(abs(zz-zmin)/dz)*dz
        elif w3d.solvergeom == w3d.XYgeom:
            nd = 2
            gdx = [0.,dx,0.,dx]
            gdy = [0.,0.,dy,dy]
            xg = xmin+aint(abs(xx-xmin)/dx)*dx
            yg = ymin+aint(abs(yy-ymin)/dy)*dy
        elif w3d.solvergeom in [w3d.Rgeom]:
            nd = 1
            gdx = [0.,dx]
            # --- Like above, the radius is calculated in the temporary, but the
            # --- original particle position is used below.
            # --- These two lines calculating rr give the same result, but the second
            # --- is probably faster
            #rr = sqrt(xx**2 + yy**2)
            rr = take(rr,iclose-i1)
            xg = xmin+aint(abs(rr-xmin)/dx)*dx
        elif w3d.solvergeom == w3d.Ygeom:
            nd = 1
            gdy = [0.,dy]
            yg = ymin+aint(abs(yy-ymin)/dy)*dy
        elif w3d.solvergeom == w3d.Zgeom:
            nd = 1
            gdz = [0.,dz]
            zg = zmin+aint(abs(zz-zmin)/dz)*dz

        nn = len(iclose)
        pp = zeros(nn,'d')

        # --- Loop over the corners of the grid cell
        for i in range(2**nd):

            # --- Get id of the conductor that the particles are near
            # --- See comments in updateconductors regarding reducedisinside
            # --- An optimization trick is to shift the grid rather than
            # --- the particles, avoiding adding scalars to arrays.
            if w3d.solvergeom in [w3d.XYZgeom]:
                getgridngp3d(nn,xg,yg,zg,pp,nx,ny,nz,self.reducedisinside,
                             xmin-gdx[i],xmax-gdx[i],ymin-gdy[i],ymax-gdy[i],
                             zmin-gdz[i],zmax-gdz[i],0.,w3d.l2symtry,w3d.l4symtry)
            elif w3d.solvergeom in [w3d.XZgeom,w3d.RZgeom]:
                xgsym,ygsym = self.applysymmetry(xg,0)
                getgridngp2d(nn,xgsym,zg,pp,nx,nz,self.reducedisinside[:,iya,:],
                             xmin-gdx[i],xmax-gdx[i],zmin-gdz[i],zmax-gdz[i])
            elif w3d.solvergeom == w3d.XYgeom:
                xgsym,ygsym = self.applysymmetry(xg,yg)
                getgridngp2d(nn,xgsym,ygsym,pp,nx,ny,self.reducedisinside[:,:,iza],
                             xmin-gdx[i],xmax-gdx[i],ymin-gdy[i],ymax-gdy[i])
            elif w3d.solvergeom == w3d.Rgeom:
                xgsym,ygsym = self.applysymmetry(xg,yg)
                getgridngp1d(nn,xgsym,pp,nx,self.reducedisinside[:,iya,iza],
                             xmin-gdx[i],xmax-gdx[i])
            elif w3d.solvergeom == w3d.Ygeom:
                xgsym,ygsym = self.applysymmetry(0,yg)
                getgridngp1d(nn,ygsym,pp,ny,self.reducedisinside[ixa,:,iza],
                             ymin-gdy[i],ymax-gdy[i])
            elif w3d.solvergeom == w3d.zgeom:
                getgridngp1d(nn,zgsym,pp,nz,self.reducedisinside[ixa,iya,:],
                             zmin-gdz[i],zmax-gdz[i])

            # --- Loop over the conductors, removing particles that are found inside
            # --- of each.
            for c in self.conductors:

                # --- Get indices relative to the temporary arrays.
                # --- Note that iclose is relative to the full particle arrays.
                itempclose=arange(nn)

                # --- Get indices of particles that are close to the conductor
                ii = compress(pp == c.condid,itempclose)

                # --- If there are no particles close, then skip to the next conductor
                if len(ii) == 0: continue

                # --- Get positions of the particles that are close
                xc = take(xx,ii)
                yc = take(yy,ii)
                zc = take(zz,ii)

                # --- Find the particles that are currently inside and down-select
                # --- the indices. The nint is needed since the quantities is used in
                # --- logical expressions below which require ints.
                xcsym,ycsym = self.applysymmetry(xc,yc)
                currentisinside = nint(c.isinside(xcsym,ycsym,zc).isinside)
                iic = compress(currentisinside,ii)
                ic = take(iclose,iic)

                if self.lrefineallintercept:
                    # --- Refine whether or not particles are lost by taking small time
                    # --- steps, starting from the old position. Note that it is possible
                    # --- that particles that were lost may not be lost upon refinement,
                    # --- and similarly, particles that were not lost, may be lost upon
                    # --- refinement.
                    # --- Get the old coordinates of particles that are close.
                    iclose1 = take(iclose,ii)
                    xo = take(top.pgroup.pid[:,self.xoldpid],iclose1)
                    yo = take(top.pgroup.pid[:,self.yoldpid],iclose1)
                    zo = take(top.pgroup.pid[:,self.zoldpid],iclose1)
                    uxo = take(top.pgroup.pid[:,self.uxoldpid],iclose1)
                    uyo = take(top.pgroup.pid[:,self.uyoldpid],iclose1)
                    uzo = take(top.pgroup.pid[:,self.uzoldpid],iclose1)
                    oldisOK = nint(take(top.pgroup.pid[:,self.oldisOK],iclose1))

                    # --- Get the current fields
                    ex = take(top.pgroup.ex,iclose1)
                    ey = take(top.pgroup.ey,iclose1)
                    ez = take(top.pgroup.ez,iclose1)
                    bx = take(top.pgroup.bx,iclose1)
                    by = take(top.pgroup.by,iclose1)
                    bz = take(top.pgroup.bz,iclose1)

                    # --- This is a possible optmization, which skips particles to
                    # --- far from the conductor to possibly hit it.

                    # --- Get the largest distance that the particles could travel
                    # --- in one time step.
                    qom = top.pgroup.sq[js]/top.pgroup.sm[js]
                    xchange = abs(uxo*top.dt) + abs(0.5*qom*ex*top.dt**2)
                    ychange = abs(uyo*top.dt) + abs(0.5*qom*ey*top.dt**2)
                    zchange = abs(uzo*top.dt) + abs(0.5*qom*ez*top.dt**2)
                    maxchange = sqrt(xchange**2 + ychange**2 + zchange**2)

                    # --- Compare the largest travel distance to the distance from the
                    # --- conductor, and skip particles that are far enough away that
                    # --- they would not hit the conductor. Do a logical_or with
                    # --- currentisinside just to be sure that scraped particles are not
                    # --- missed.
                    distance = c.distance(xo,yo,zo)
                    closeenough = logical_or((maxchange > distance.distance),
                                             currentisinside)

                    # --- Downselect the particles which are close enough to the
                    # --- coductor that they could be lost.
                    ii = compress(closeenough,ii)
                    if len(ii) == 0: continue
                    iclose1 = take(iclose,ii)
                    xc = compress(closeenough,xc)
                    yc = compress(closeenough,yc)
                    zc = compress(closeenough,zc)
                    xo = compress(closeenough,xo)
                    yo = compress(closeenough,yo)
                    zo = compress(closeenough,zo)
                    uxo = compress(closeenough,uxo)
                    uyo = compress(closeenough,uyo)
                    uzo = compress(closeenough,uzo)
                    ex = compress(closeenough,ex)
                    ey = compress(closeenough,ey)
                    ez = compress(closeenough,ez)
                    bx = compress(closeenough,bx)
                    by = compress(closeenough,by)
                    bz = compress(closeenough,bz)
                    currentisinside = compress(closeenough,currentisinside)
                    oldisOK = compress(closeenough,oldisOK)

                    # --- Create some temporaries
                    itime = None
                    dt = top.dt*top.pgroup.ndts[js]*top.pgroup.dtscale[js]*ones(len(ii))
                    q = top.pgroup.sq[js]
                    m = top.pgroup.sm[js]

                    # --- Do the refinement calculation. The currentisinside argument sets
                    # --- when the current position is replaced by the refined position.
                    # --- If the particle is currently lost but in the refined
                    # --- calculation is not lost, then the replace the current position
                    # --- with that refined position that is not lost. Similarly, if the
                    # --- particle is currently not lost, but in the refined calculation
                    # --- is lost, then replace the current position with the refined
                    # --- position.
                    refinedisinside = zeros(len(xc),'l')
                    self.refineintercept(c,xc,yc,zc,xo,yo,zo,uxo,uyo,uzo,
                                         ex,ey,ez,bx,by,bz,itime,dt,q,m,currentisinside,
                                         refinedisinside)

                    # --- For some newly created particles, there was no old saved data.
                    # --- In those cases, ignore the result of the refined calculation,
                    # --- and use the result obtained originally in currentisinside.
                    refinedisinside = where(oldisOK,refinedisinside,currentisinside)

                    # --- iic lists the particles that are lost in the refined
                    # --- calculation. These will be scraped. Particles which were
                    # --- considered lost but where not lost based on the refined
                    # --- calculation still need to have their refined positions checked
                    # --- against other conductors. There is a possible problem here.
                    # --- The refined trajectory could put the particle in a different
                    # --- grid cell than the original, and it could be inside a conductor
                    # --- that the original wasn't considered close too. This would leave
                    # --- that particle unscraped at that refined position but inside
                    # --- a conductor. This case would be messy to deal with, requiring
                    # --- a second loop over conductors.
                    iic = compress(refinedisinside,ii)
                    ic = take(iclose,iic)

                    # --- Do the replacements as described above. Note that for lost
                    # --- particles, xc,yc,zc hold the positions of the particles one
                    # --- small time step into the conductor. Don't do the replacement
                    # --- for new particles since the old data is no good.
                    iio         = (currentisinside | refinedisinside) & oldisOK
                    iiu         = compress(iio,arange(shape(xc)[0]))
                    iuserefined = compress(iio,iclose1)
                    put(top.pgroup.xp, iuserefined,take(xc, iiu))
                    put(top.pgroup.yp, iuserefined,take(yc, iiu))
                    put(top.pgroup.zp, iuserefined,take(zc, iiu))
                    put(top.pgroup.uxp,iuserefined,take(uxo,iiu))
                    put(top.pgroup.uyp,iuserefined,take(uyo,iiu))
                    put(top.pgroup.uzp,iuserefined,take(uzo,iiu))

                    # --- Note that the old values of the positions are changed
                    # --- only for particles for which the refined calculation
                    # --- shows they are lost. This is needed for the interception
                    # --- calculation done in savecondid.
                    iclose2 = compress(refinedisinside,iclose1)
                    if len(iclose2) > 0:
                        put(top.pgroup.pid[:,self.xoldpid],iclose2,compress(refinedisinside,xo))
                        put(top.pgroup.pid[:,self.yoldpid],iclose2,compress(refinedisinside,yo))
                        put(top.pgroup.pid[:,self.zoldpid],iclose2,compress(refinedisinside,zo))

                # --- If no particle are inside the conductor, then skip to the next one
                if len(iic) == 0: continue

                if c.material == 'reflector':
                    # --- For lost new particles, which have no old data, not much can be
                    # --- done, so they are set to be removed.
                    oldisOK = nint(take(top.pgroup.pid[:,self.oldisOK],ic))
                    icnew = compress(logical_not(oldisOK),ic)
                    if len(icnew) > 0:
                        put(top.pgroup.gaminv,icnew,0.)
                        # --- Only old particles will be reflected.
                        ic = compress(oldisOK,ic)

                    # --- For particles which are inside, replace the position with
                    # --- the old position and reverse the velocity.
                    put(top.pgroup.xp,ic,take(top.pgroup.pid[:,self.xoldpid],ic))
                    put(top.pgroup.yp,ic,take(top.pgroup.pid[:,self.yoldpid],ic))
                    put(top.pgroup.zp,ic,take(top.pgroup.pid[:,self.zoldpid],ic))
                    if self.lsaveoldvelocities:
                        # --- If its available, use the old velocity.
                        # --- Should this be the default?
                        put(top.pgroup.uxp,ic,-take(top.pgroup.pid[:,self.uxoldpid],ic))
                        put(top.pgroup.uyp,ic,-take(top.pgroup.pid[:,self.uyoldpid],ic))
                        put(top.pgroup.uzp,ic,-take(top.pgroup.pid[:,self.uzoldpid],ic))
                    else:
                        # --- Otherwise use the new velocity. Can this lead to errors?
                        put(top.pgroup.uxp,ic,-take(top.pgroup.uxp,ic))
                        put(top.pgroup.uyp,ic,-take(top.pgroup.uyp,ic))
                        put(top.pgroup.uzp,ic,-take(top.pgroup.uzp,ic))

                else:

                    # --- For particles which are inside, set gaminv to 0, the lost
                    # --- particle flag
                    put(top.pgroup.gaminv,ic,0.)

                # --- Remove the already handled particles, returning if there
                # --- are no more.
                put(iclose,iic,-1)
                iclose = compress(iclose>=0,iclose)
                nn = len(iclose)
                if nn == 0: return
                put(itempclose,iic,-1)
                itempclose = compress(itempclose>=0,itempclose)
                xx = take(xx,itempclose)
                yy = take(yy,itempclose)
                zz = take(zz,itempclose)
                pp = take(pp,itempclose)
                if w3d.solvergeom in [w3d.XYZgeom,w3d.XZgeom,w3d.XYgeom,w3d.RZgeom,w3d.Rgeom]:
                    xg = take(xg,itempclose)
                if w3d.solvergeom in [w3d.XYZgeom,w3d.XYgeom,w3d.Ygeom]:
                    yg = take(yg,itempclose)
                if w3d.solvergeom in [w3d.XYZgeom,w3d.XZgeom,w3d.RZgeom,w3d.Zgeom]:
                    zg = take(zg,itempclose)


    def savecondid(self,js,local=0):
        """Saves information about lost particles, including the conductor id
    where the particles are lost, the intercept point and angle of incidence where
    the particle struck the conductor, and integrated data for the conductors,
    counting the current lost on the conductor.
        """
        jsid = top.pgroup.sid[js]

        # --- Just return if there are no lost particles.
        if top.npslost[jsid] == 0:
            if self.lcollectlpdata and not local:
                # --- If data is being collected, the 0 from this processor must still
                # --- be added to the sum.
                for c in self.conductors:
                    # --- This parallelsum coordinates with the ones below.
                    w=parallelsum(0.)
                    if w != 0.:
                        c.lostparticles_data.append(array([top.time,
                                                           w*top.pgroup.sq[js]*top.pgroup.sw[js],
                                                           top.dt,
                                                           jsid]))
            return

        # --- First make sure there is extra space in the pidlost array.
        pidspace = 1
        if self.lsaveintercept: pidspace = 4
        if top.npidlost < top.npid+pidspace:
            top.npidlost = top.npid + pidspace
            gchange("LostParticles")

        # --- Much of this code is duplicated from scrape above so if it changes,
        # --- this should change as well.
        dx,dy,dz,nx,ny,nz,ix,iy,iz = self.grid.getmeshsize(self.mglevel)
        xmin = self.grid.xmmin + ix*dx
        xmax = self.grid.xmmin + (ix+nx)*dx
        ymin = self.grid.ymmin + iy*dy
        ymax = self.grid.ymmin + (iy+ny)*dy
        zmin = self.grid.zmmin + iz*dz + top.zbeam
        zmax = self.grid.zmmin + (iz+nz)*dz + top.zbeam
        isinside = self.grid.isinside

        i1 = top.inslost[jsid] - 1
        i2 = top.inslost[jsid] + top.npslost[jsid] - 1
        xx = top.xplost[i1:i2]
        yy = top.yplost[i1:i2]
        zz = top.zplost[i1:i2]

        if w3d.solvergeom in [w3d.RZgeom,w3d.Rgeom]:
            xx = sqrt(xx**2 + yy**2)
            yy = zeros(len(xx),'d')

        # --- Get the indices of all lost particles that havn't been localized
        # --- to a conductor.
        iscrape = compress(top.pidlost[i1:i2,-1]==0,arange(i1,i2))
        if self.lcollectlpdata:iscrape1=iscrape.copy()

        # --- Duplicate the particle list eight times, once for each corner.
        iscrape = numpy.repeat(iscrape,8)
        nn = len(iscrape)
        x8 = take(xx,iscrape-i1)
        y8 = take(yy,iscrape-i1)
        z8 = take(zz,iscrape-i1)
        xg = xmin+aint(abs(x8-xmin)/dx)*dx + array(int(nn/8)*[0.,dx,0.,dx,0.,dx,0.,dx])
        yg = ymin+aint(abs(y8-ymin)/dy)*dy + array(int(nn/8)*[0.,0.,dy,dy,0.,0.,dy,dy])
        zg = zmin+aint(abs(z8-zmin)/dz)*dz + array(int(nn/8)*[0.,0.,0.,0.,dz,dz,dz,dz])
        pp = zeros(nn,'d')

        # --- The ixa etc are the location of the x=0 plane. This is needed
        # --- since in certain cases, the size of the collapsed dimension
        # --- can be greater than one and the 0 plane needs to be found, it
        # --- is the plane where the data is stored. This is true for example
        # --- with the RZ EM solver.
        # --- The check of nx > 0 etc is done since in certain cases the xmin
        # --- can be nonzero when nx is zero, giving an erroneous value for ixa.
        # --- This is the case when the quasistatic solver is being used.
        ixa = iya = iza = 0
        if nx > 0: ixa = nint(-xmin/dx)
        if ny > 0: iya = nint(-ymin/dy)
        if nz > 0: iza = nint(-zmin/dz)

        # --- Get conductor id that particles are near
        # --- See comments in updateconductors regarding reducedisinside
        if w3d.solvergeom in [w3d.XYZgeom]:
            getgridngp3d(nn,xg,yg,zg,pp,
                         nx,ny,nz,self.reducedisinside,xmin,xmax,ymin,ymax,zmin,zmax,0.,
                         w3d.l2symtry,w3d.l4symtry)
        elif w3d.solvergeom == w3d.RZgeom or w3d.solvergeom == w3d.XZgeom:
            xgsym,ygsym = self.applysymmetry(xg,0)
            getgridngp2d(nn,xgsym,zg,pp,nx,nz,self.reducedisinside[:,iya,:],xmin,xmax,zmin,zmax)
        elif w3d.solvergeom == w3d.XYgeom:
            xgsym,ygsym = self.applysymmetry(xg,yg)
            getgridngp2d(nn,xgsym,ygsym,pp,nx,ny,self.reducedisinside[:,:,iza],xmin,xmax,ymin,ymax)
        elif w3d.solvergeom == w3d.Rgeom:
            xgsym,ygsym = self.applysymmetry(xg,0)
            getgridngp1d(nn,xgsym,pp,nx,self.reducedisinside[:,iya,iza],xmin,xmax)
        elif w3d.solvergeom == w3d.Ygeom:
            xgsym,ygsym = self.applysymmetry(0,yg)
            getgridngp1d(nn,ygsym,pp,ny,self.reducedisinside[ixa,:,iza],ymin,ymax)
        elif w3d.solvergeom == w3d.Zgeom:
            getgridngp1d(nn,zg,pp,nz,self.reducedisinside[ixa,iya,:],zmin,zmax)
        else:
            raise Exception("The particle scraping only works for XYZ, XZ, XY and RZ geometry")

        if w3d.solvergeom in [w3d.RZgeom,w3d.Rgeom]:
            xx = top.xplost[i1:i2]
            yy = top.yplost[i1:i2]
            x8 = take(xx,iscrape-i1)
            y8 = take(yy,iscrape-i1)

        # --- Loop over the conductors, removing particles inside of each.
        for c in self.conductors:
            ii = compress(pp == c.condid,arange(nn))
            if len(ii) == 0:
                if self.lcollectlpdata and not local:
                    # --- This parallelsum coordinates with the other processors
                    w=parallelsum(0.)
                    if w != 0.:
                        c.lostparticles_data.append(array([top.time,
                                                           w*top.pgroup.sq[js]*top.pgroup.sw[js],
                                                           top.dt,
                                                           jsid]))
                continue
            xc = take(x8,ii)
            yc = take(y8,ii)
            zc = take(z8,ii)

            xcsym,ycsym = self.applysymmetry(xc,yc)
            ic = take(iscrape,ii)
            ic = compress(c.isinside(xcsym,ycsym,zc).isinside,ic)
            if len(ic) == 0:
                if self.lcollectlpdata and not local:
                    # --- This parallelsum coordinates with the other processors
                    w=parallelsum(0.)
                    if w != 0.:
                        c.lostparticles_data.append(array([top.time,
                                                           w*top.pgroup.sq[js]*top.pgroup.sw[js],
                                                           top.dt,
                                                           jsid]))
                continue
            # --- For particles which are inside, set pid to the id of the conductor
            # --- where the particle is lost.
            put(top.pidlost[:,-1],ic,c.condid)
            # --- Save location and surface normal where particle intercepted the
            # --- conductor.
            if self.lsaveintercept:

                # --- Don't calculate the intercept for new particles, for which there
                # --- is no old data saved.
                oldisOK = nint(take(top.pidlost[:,self.oldisOK],ic))
                icnew = compress(logical_not(oldisOK),ic)
                if len(icnew) > 0:
                    # --- Set the conductor ID to zero, so that these particles are
                    # --- ignored.
                    put(top.pidlost[:,-1],icnew,0)
                    # --- Downselect to only include the older particles.
                    ic = compress(oldisOK,ic)

                xc = take(xx,ic-i1)
                yc = take(yy,ic-i1)
                zc = take(zz,ic-i1)
                xo = take(top.pidlost[:,self.xoldpid],ic)
                yo = take(top.pidlost[:,self.yoldpid],ic)
                zo = take(top.pidlost[:,self.zoldpid],ic)

                dt = top.dt*top.pgroup.ndts[js]*top.pgroup.dtscale[js]
                if self.lrefineintercept:
                    uxo = take(top.pidlost[:,self.uxoldpid],ic)
                    uyo = take(top.pidlost[:,self.uyoldpid],ic)
                    uzo = take(top.pidlost[:,self.uzoldpid],ic)
                    ex = take(top.exlost,ic)
                    ey = take(top.eylost,ic)
                    ez = take(top.ezlost,ic)
                    bx = take(top.bxlost,ic)
                    by = take(top.bylost,ic)
                    bz = take(top.bzlost,ic)
                    itime = zeros(len(ic),'d')
                    dt *= ones(len(ic))
                    q = top.pgroup.sq[js]
                    m = top.pgroup.sm[js]
                    self.refineintercept(c,xc,yc,zc,xo,yo,zo,uxo,uyo,uzo,
                                         ex,ey,ez,bx,by,bz,itime,dt,q,m,0,
                                         zeros(len(xc),'l'))
                else:
                    itime = 0.

                if self.lrefineallintercept:
                    # --- In this case, the old and new positions are the points
                    # --- just outside and inside of the conductor, differing by the
                    # --- refined time step size. That refined step size is needed
                    # --- to get the correct approximation to the velocity.
                    bx = take(top.bxlost,ic)
                    by = take(top.bylost,ic)
                    bz = take(top.bzlost,ic)
                    q = top.pgroup.sq[js]
                    m = top.pgroup.sm[js]
                    dt = dt/self.getrefinedtimestepnumber(dt,bx,by,bz,q,m)

                # --- use an approximate calculation.
                if self.interceptvelocitymethod == 'finitedifference':
                    vx = (xc-xo)/dt
                    vy = (yc-yo)/dt
                    vz = (zc-zo)/dt
                elif self.interceptvelocitymethod == 'actualvelocity':
                    ux = take(top.uxplost,ic)
                    uy = take(top.uyplost,ic)
                    uz = take(top.uzplost,ic)
                    gi = 1./sqrt(1.+(ux**2+uy**2+uz**2)/clight**2)
                    vx = ux*gi
                    vy = uy*gi
                    vz = uz*gi
                else:
                    raise ValueError("interceptvelocitymethod has an incorrect value, %s"%self.interceptvelocitymethod)

                # --- get v in lab frame
                if top.boost_gamma>1.:
                    boost_beta  = -sqrt(1.-1./top.boost_gamma**2)
                    fact = 1./(1.-boost_beta*vz/clight)
                    vx = vx*fact/top.boost_gamma
                    vy = vy*fact/top.boost_gamma
                    vz = (vz-boost_beta*clight)*fact

                intercept = c.intercept(xc,yc,zc,vx,vy,vz)
                dtintercept = (sqrt((xc - intercept.xi)**2 +
                                    (yc - intercept.yi)**2 +
                                    (zc - intercept.zi)**2)/
                      dvnz(sqrt(vx**2 + vy**2 + vz**2))) + itime

                put(top.xplost,ic,intercept.xi)
                put(top.yplost,ic,intercept.yi)
                put(top.zplost,ic,intercept.zi)

                # --- Also, reset the velocities
                if top.lrelativ:
                    beta = sqrt(vx**2 + vy**2 + vz**2)/clight
                    # --- If beta is too large, then reset the velocities. Note that
                    # --- there may be some other but that is making beta too large.
                    # --- It looks like in some cases, the x and xold etc positions are
                    # --- inconsistent and very far from each each, giving an errorneous
                    # --- value of vx etc.
                    # --- betacorrection is written this way to avoid dividing by zero
                    # --- when beta == 0.
                    betacorrection = 1./where(beta >= 0.99999,beta/0.99999,1.)
                    beta = minimum(beta,0.99999)
                    gamma = 1./sqrt((1.-beta)*(1.+beta))
                    ux = vx*(gamma*betacorrection)
                    uy = vy*(gamma*betacorrection)
                    uz = vz*(gamma*betacorrection)
                else:
                    gamma = 1.
                    ux = vx
                    uy = vy
                    uz = vz
                put(top.uxplost,ic,ux)
                put(top.uyplost,ic,uy)
                put(top.uzplost,ic,uz)
                put(top.gaminvlost,ic,1./gamma)

                # --- Set the angle of incidence and time of interception
                put(top.pidlost[:,-3],ic,intercept.itheta)
                put(top.pidlost[:,-2],ic,intercept.iphi)
                put(top.pidlost[:,-4],ic,top.time - dtintercept)

            if self.lcollectlpdata:
                pidlostcondid = take(top.pidlost[:,-1],iscrape1)
                pidtoconsider = compress(pidlostcondid==c.condid,iscrape1)
                if top.wpid==0:
                    w = len(pidtoconsider)
                else:
                    w = sum(take(top.pidlost[:,top.wpid-1],pidtoconsider))
                # --- This parallelsum coordinates with the ones above
                if not local:w=parallelsum(w)
                c.lostparticles_data.append(array([top.time,
                                                   w*top.pgroup.sq[js]*top.pgroup.sw[js],
                                                   top.dt,
                                                   jsid]))


    def getrefinedtimestepnumber(self,dt,bx,by,bz,q,m):
        """Calculates a refined time step size for each particle that is a
    fraction of the cyclotron period (calculated from the B field of each
    particle).
        """
        # --- The cyclotron frequency for each particle
        magB = sqrt(bx**2 + by**2 + bz**2)
        omegac = q/m*magB

        # --- Get the number of steps for each particle.
        # --- This is set by self.nstepsperorbit which is the number of steps
        # --- per cyclotron orbit.
        isteps = nint(dt*omegac/(2.*pi)*self.nstepsperorbit)

        # --- So that the orbit is always refined, at least a little, set the
        # --- minimum number of steps to be self.nstepsperorbit. This is helpful
        # --- for example if the B field is zero.
        isteps = maximum(isteps,self.nstepsperorbit)

        # --- Now, return the refined step number
        return isteps

    def refineintercept(self,c,xc,yc,zc,xo,yo,zo,uxo,uyo,uzo,ex,ey,ez,bx,by,bz,
                        itime,dt,q,m,luserefinedifnotlost,isinside):
        """Refine the location of the intercept, advancing the particle using a
    time step that is small compared to the cyclotron period of each particle,
    starting from the old position before the particle was lost.
    c: the conductor
    xc,yc,zc: input holding the current particle position
              output holding the point just inside the conductor
    xo,yo,zo: input holding the old particle position
              output holding the point just outside the conductor
    uxo,uyo,uzo: input holding old velocity, time synchronized with xo,yo,zo
                 output holding velocity near time when particle entered the
                 conductor
    ex,ey,ez,bx,by,bz: fixed E and B fields
    itime: output holding time that the particle reached xo,yo,zo
           If input value is None, then the time is not saved
    dt: output holding particle time step sizes
    q,m: charge and mass of the particles
    luserefinedifnotlost: when true, if the refined particle orbit is not lost,
                          then replace replace the current position with the
                          refined position
        """
        # --- Note that this routine is written in a non-pythonic way, where the
        # --- changes are made directly in the input arrays.

        # --- Get the number of particles
        nn = len(xo)

        # --- Get the number of refined time steps for each particle
        isteps = self.getrefinedtimestepnumber(dt,bx,by,bz,q,m)

        # --- Get the maximum number of steps needed
        nsteps  = max(isteps)

        # --- Calculate the overcycled step size for each particle
        # --- Note that dtover will be modified in the loop below,
        # --- and dt is used to return the time step size.
        dtover = dt/isteps
        dt[:] = dtover

        # --- Recalculate gaminv. The code could also save the old gaminv, but
        # --- this should be equivalent.
        if top.lrelativ:
            usq = (uxo**2 + uyo**2 + uzo**2)
            gamma = sqrt(1. + usq/clight**2)
            gaminv = 1./gamma
        else:
            gaminv = ones(nn,'d')

        # --- Save the positions. This is needed so that the data can be
        # --- restored if no intercept is found below.
        xcsave = xc.copy()
        ycsave = yc.copy()
        zcsave = zc.copy()
        xosave = xo.copy()
        yosave = yo.copy()
        zosave = zo.copy()
        uxosave = uxo.copy()
        uyosave = uyo.copy()
        uzosave = uzo.copy()

        # --- Get the starting positions of the advance. These should all be
        # --- outside of the conductor. The loop below advances the particles
        # --- using these 'c' arrays.
        xc[:] = xo
        yc[:] = yo
        zc[:] = zo
        if itime is not None:
            itime[:] = 0.

        # --- Do the over cycling loop. Note that all particles are advanced by
        # --- the maximum number of steps, though once a particle goes inside
        # --- the conductor, its time step is set to zero so the coordinates
        # --- don't change anymore.
        # --- One possible optimization is to have the fortran advancing
        # --- routines skip particles that have a zero time step size.
        for it in range(nsteps):

            # --- Do a full split leap-frog advance (with constant E and B fields)
            # --- Note that this does the advance in place, directly changing the
            # --- input arrays.
            bpusht3d(nn,uxo,uyo,uzo,gaminv,bx,by,bz,q,m,dtover,0.5,top.ibpush)
            epusht3d(nn,uxo,uyo,uzo,ex,ey,ez,q,m,dtover,0.5)
            gammaadv(nn,gaminv,uxo,uyo,uzo,top.gamadv,top.lrelativ)
            xpusht3d(nn,xc,yc,zc,uxo,uyo,uzo,gaminv,dtover)
            epusht3d(nn,uxo,uyo,uzo,ex,ey,ez,q,m,dtover,0.5)
            gammaadv(nn,gaminv,uxo,uyo,uzo,top.gamadv,top.lrelativ)
            bpusht3d(nn,uxo,uyo,uzo,gaminv,bx,by,bz,q,m,dtover,0.5,top.ibpush)

            # --- This provides a nice diagnostic for testing
            #plp(yc,xc,marker=circle,color=green)
            #pldj(xo,yo,xc,yc,color=green)

            # --- Check whether the new positions are inside of the conductor.
            xcsym,ycsym = self.applysymmetry(xc,yc)
            isinside[:] = c.isinside(xcsym,ycsym,zc).isinside

            # --- Kludgy code to handle some boundary conditions
            # --- This code is OK in serial, but in parallel, it will break with
            # --- periodic b.c.s and if a particle crosses a parallel domain
            # --- boundary.
            zmmin = w3d.zmmin + top.zbeam
            zmmax = w3d.zmmax + top.zbeam
            if top.pboundxy == periodic:
                xc[:] = where(xc > w3d.xmmax,xc-(w3d.xmmax-w3d.xmmin),xc)
                xc[:] = where(xc < w3d.xmmin,xc+(w3d.xmmax-w3d.xmmin),xc)
                yc[:] = where(yc > w3d.ymmax,yc-(w3d.ymmax-w3d.ymmin),yc)
                yc[:] = where(yc < w3d.ymmin,yc+(w3d.ymmax-w3d.ymmin),yc)
            elif top.pboundxy == reflect:
                uxo[:] = where(xc > w3d.xmmax,-uxo,uxo)
                uxo[:] = where(xc < w3d.xmmin,-uxo,uxo)
                uyo[:] = where(yc > w3d.ymmax,-uyo,uyo)
                uyo[:] = where(yc < w3d.ymmin,-uyo,uyo)
                xc[:] = where(xc > w3d.xmmax,2.*w3d.xmmax-xc,xc)
                xc[:] = where(xc < w3d.xmmin,2.*w3d.xmmin-xc,xc)
                yc[:] = where(yc > w3d.ymmax,2.*w3d.ymmax-yc,yc)
                yc[:] = where(yc < w3d.ymmin,2.*w3d.ymmin-yc,yc)
            elif top.pboundxy == absorb:
                isinside[:] = where(xc > w3d.xmmax,1,isinside)
                isinside[:] = where(xc < w3d.xmmin,1,isinside)
                isinside[:] = where(yc > w3d.ymmax,1,isinside)
                isinside[:] = where(yc < w3d.ymmin,1,isinside)
            if top.pbound0 == periodic or top.pboundnz == periodic:
                zc[:] = where(zc > zmax,zc-(zmax-zmmin),zc)
                zc[:] = where(zc < zmmin,zc+(zmax-zmmin),zc)
            if top.pboundnz == reflect:
                uzo[:] = where(zc > zmmax,-uzo,uzo)
                zc[:] = where(zc > zmmax,2.*zmmax-zc,zc)
            if top.pbound0 == reflect:
                uzo[:] = where(zc < zmmin,-uzo,uzo)
                zc[:] = where(zc < zmmin,2.*zmmin-zc,zc)
            if top.pboundnz == absorb:
                isinside[:] = where(zc > zmmax,1,isinside)
            if top.pbound0 == absorb:
                isinside[:] = where(zc < zmmin,1,isinside)

            # --- For the particles that are still outside, set the old positions
            # --- to be the updated positions
            xo[:] = where(isinside,xo,xc)
            yo[:] = where(isinside,yo,yc)
            zo[:] = where(isinside,zo,zc)

            # --- If a particle is inside the conductor, then stop advancing it,
            # --- setting its time step size to zero.
            dtover = where(isinside,0.,dtover)

            # --- Now, advance itime. Note that for particles that are now inside,
            # --- itime is not advanced since it is the time just before the particle
            # --- enters the conductor.
            if itime is not None:
                itime[:] = itime + dtover

            # --- Quit the loop if all intercepts have been found.
            if alltrue(dtover==0.): break

        # --- Check for cases where no interception was found. In those cases,
        # --- restore the original data since that will at least not cause
        # --- a code problem. This checks if dtover hasn't been zeroed out
        # --- which means that the particle was never flagged as being inside
        # --- in the loop above.
        if sometrue(dtover > 0.):
            userefined = ((dtover == 0.) | luserefinedifnotlost)
            xc[:] = where(userefined,xc,xcsave)
            yc[:] = where(userefined,yc,ycsave)
            zc[:] = where(userefined,zc,zcsave)
            xo[:] = where(userefined,xo,xosave)
            yo[:] = where(userefined,yo,yosave)
            zo[:] = where(userefined,zo,zosave)
            uxo[:] = where(userefined,uxo,uxosave)
            uyo[:] = where(userefined,uyo,uyosave)
            uzo[:] = where(userefined,uzo,uzosave)
            if itime is not None:
                itime[:] = where(userefined,itime,0.)

    def fastscrape(self,js):
        """A fast but not precise method of scraping particles. In this method,
    the grid is setup so that each grid point stores the distance of that grid
    point to the nearest conductor surface, with negative values indicating
    that the grid point is inside of the conductor. This data is interpolated to
    the particles and particles that get a negative value are considered to be
    inside the conductor and are scraped. It is only approximate due to
    interpolating errors from the grid. This will normally be called automatically.
        """
        # --- If there are no particles in this species, that nothing needs
        # --- to be done
        if top.pgroup.nps[js] == 0: return

        # --- Get mesh information into local variables
        dx,dy,dz,nx,ny,nz,ix,iy,iz = self.grid.getmeshsize(self.mglevel)
        xmin = self.grid.xmmin + ix*dx
        xmax = self.grid.xmmin + (ix+nx)*dx
        ymin = self.grid.ymmin + iy*dy
        ymax = self.grid.ymmin + (iy+ny)*dy
        zmin = self.grid.zmmin + iz*dz + top.zbeam
        zmax = self.grid.zmmin + (iz+nz)*dz + top.zbeam
        distances = self.grid.distances

        # --- Get handy references to the particles in the species
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        xx = top.pgroup.xp[i1:i2]
        yy = top.pgroup.yp[i1:i2]
        zz = top.pgroup.zp[i1:i2]
        pp = zeros(top.pgroup.nps[js],'d')

        # --- Find the distances of each particle to the conductors.
        # --- The interpolates from the distances grid. The results are
        # --- put into the array pp.
        if w3d.solvergeom in [w3d.XYZgeom]:
            getgrid3d(top.pgroup.nps[js],xx,yy,zz,pp,
                      nx,ny,nz,distances,xmin,xmax,ymin,ymax,zmin,zmax,
                      w3d.l2symtry,w3d.l4symtry)
        elif w3d.solvergeom == w3d.RZgeom:
            # --- Note that for RZ, the radius is calculated for this, but
            # --- the original particle position is used below.
            rr = sqrt(xx**2 + yy**2)
            getgrid2d(top.pgroup.nps[js],rr,zz,pp,nx,nz,distances[:,0,:],
                      xmin,xmax,zmin,zmax)
        elif w3d.solvergeom == w3d.XZgeom:
            xsym,ysym = self.applysymmetry(xx,0)
            getgrid2d(top.pgroup.nps[js],xsym,zz,pp,nx,nz,distances[:,0,:],
                      xmin,xmax,zmin,zmax)
        elif w3d.solvergeom == w3d.XYgeom:
            xsym,ysym = self.applysymmetry(xx,yy)
            getgrid2d(top.pgroup.nps[js],xsym,ysym,pp,nx,ny,distances[:,:,0],
                      xmin,xmax,ymin,ymax)
        elif w3d.solvergeom == w3d.Rgeom:
            # --- Note that for R, the radius is calculated for this, but
            # --- the original particle position is used below.
            rr = sqrt(xx**2 + yy**2)
            getgrid1d(top.pgroup.nps[js],rr,pp,nx,distances[:,0,0],
                      xmin,xmax)
        elif w3d.solvergeom == w3d.Ygeom:
            xsym,ysym = self.applysymmetry(0,yy)
            getgrid2d(top.pgroup.nps[js],ysym,pp,ny,distances[0,:,0],
                      ymin,ymax)
        elif w3d.solvergeom == w3d.Zgeom:
            getgrid2d(top.pgroup.nps[js],zsym,pp,nz,distances[0,0,:],
                      zmin,zmax)
        else:
            raise Exception("The particle scraping only works for XYZ, XY and RZ geometry")

        # --- Any particles which have a negative distance are approximately
        # --- inside of the conductors. Those are considered lost or reflected.
        ilost = compress(pp<0.,arange(i1,i2))
        if len(ilost) == 0: return

        if self.reflectiveconductors:
            # --- Note that since conductors are not distinguished, all conductors
            # --- are considered reflective if any are.
            # --- For lost new particles, which have no old data, not much can be
            # --- done, so they are set to be removed.
            oldisOK = nint(take(top.pgroup.pid[:,self.oldisOK],ilost))
            icnew = compress(logical_not(oldisOK),ilost)
            if len(icnew) > 0:
                put(top.pgroup.gaminv,icnew,0.)
                # --- Only old particles will be reflected.
                ilost = compress(oldisOK,ilost)

            # --- For particles which are inside, replace the position with
            # --- the old position and reverse the velocity.
            put(top.pgroup.xp,ilost,take(top.pgroup.pid[:,self.xoldpid],ilost))
            put(top.pgroup.yp,ilost,take(top.pgroup.pid[:,self.yoldpid],ilost))
            put(top.pgroup.zp,ilost,take(top.pgroup.pid[:,self.zoldpid],ilost))
            if self.lsaveoldvelocities:
                # --- If its available, use the old velocity.
                # --- Should this be the default?
                put(top.pgroup.uxp,ilost,-take(top.pgroup.pid[:,self.uxoldpid],ilost))
                put(top.pgroup.uyp,ilost,-take(top.pgroup.pid[:,self.uyoldpid],ilost))
                put(top.pgroup.uzp,ilost,-take(top.pgroup.pid[:,self.uzoldpid],ilost))
            else:
                # --- Otherwise use the new velocity. Can this lead to errors?
                put(top.pgroup.uxp,ilost,-take(top.pgroup.uxp,ilost))
                put(top.pgroup.uyp,ilost,-take(top.pgroup.uyp,ilost))
                put(top.pgroup.uzp,ilost,-take(top.pgroup.uzp,ilost))

        else:
            # --- For particles which are inside, set gaminv to 0, the lost
            # --- particle flag
            put(top.pgroup.gaminv,ilost,0.)

    def pdxy(self,iz=0,fullplane=0,xyantisymmetric=0,**kw):
        """Makes a plot of the internal data used to keep track of the location
        of conductors. Extra keyword arguments are passed to ppgeneric."""
        self.updategrid()
        if self.lfastscraper:
            data = self.grid.distances
        else:
            data = self.grid.isinside
        self.grid.fsdecomp = self.grid.decomp
        self.grid.solvergeom = w3d.solvergeom
        data = getdecomposedarray(data,iz=iz,bcast=0,local=0,
                                  fullplane=fullplane,
                                  xyantisymmetric=xyantisymmetric,
                                  solver=self.grid)
        del self.grid.fsdecomp
        del self.grid.solvergeom
        ppgeneric(grid=data,
                  xmin=self.grid.xmmin,xmax=self.grid.xmmax,
                  ymin=self.grid.ymmin,ymax=self.grid.ymmax,**kw)

    def pdzx(self,iy=0,fullplane=0,**kw):
        """Makes a plot of the internal data used to keep track of the location
        of conductors. Extra keyword arguments are passed to ppgeneric."""
        self.updategrid()
        if self.lfastscraper:
            data = self.grid.distances
        else:
            data = self.grid.isinside
        self.grid.fsdecomp = self.grid.decomp
        self.grid.solvergeom = w3d.solvergeom
        data = getdecomposedarray(data,iy=iy,bcast=0,local=0,
                                  fullplane=fullplane,
                                  solver=self.grid)
        del self.grid.fsdecomp
        del self.grid.solvergeom
        ppgeneric(gridt=data,
                  xmin=self.grid.zmminlocal,xmax=self.grid.zmmaxlocal,
                  ymin=self.grid.xmminlocal,ymax=self.grid.xmmaxlocal,**kw)

    def pdzy(self,ix=0,fullplane=0,**kw):
        """Makes a plot of the internal data used to keep track of the location
        of conductors. Extra keyword arguments are passed to ppgeneric."""
        self.updategrid()
        if self.lfastscraper:
            data = self.grid.distances
        else:
            data = self.grid.isinside
        self.grid.fsdecomp = self.grid.decomp
        self.grid.solvergeom = w3d.solvergeom
        data = getdecomposedarray(data,ix=ix,bcast=0,local=0,
                                  fullplane=fullplane,
                                  solver=self.grid)
        del self.grid.fsdecomp
        del self.grid.solvergeom
        ppgeneric(gridt=data,
                  xmin=self.grid.zmminlocal,xmax=self.grid.zmmaxlocal,
                  ymin=self.grid.ymminlocal,ymax=self.grid.ymmaxlocal,**kw)
