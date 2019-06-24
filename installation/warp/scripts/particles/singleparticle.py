from ..warp import *
from ..utils.appendablearray import *


class TraceParticle(object):
    """
  Class for adding a trace particle
  This must be called by all processes when running in parallel and the inputs
  (including the coordinates) must be the same across all processors.
  Creator arguments...
   - x,y,z,vx,vy,vz: 6 coordinates. They can either be
       scalars (for one particle) or sequences (for multiple particles). The
       sequences must all be of the same length. If some are sequences and some
       scalars, the scalar values will be used for all particles. If vz
       is not specified, all the particles will use top.vbeam.
   - maxsteps=1000: an estimate of the number of steps taken. OK if value too
                    small
   - savedata=1: frequency (in time steps) for saving of particle trajectories
   - js=0: species of particles
   - enforceinitboundaries=true: when true, the particle boundaries are
                                 enforced for the initial particle positions.
                                 It is a good idea to turn this off if the
                                 trace particles are setup before a generate
                                 is done. However, care is needed when running
                                 in parallel.
   - lsavefields=0: when true, save the time histories of the field data

  Available methods...
   - gett(i=0):  returns history of time for i'th particle
   - getx(i=0):  returns history of x for i'th particle
   - gety(i=0):  returns history of y for i'th particle
   - getz(i=0):  returns history of z for i'th particle
   - getvx(i=0): returns history of vx for i'th particle
   - getvy(i=0): returns history of vy for i'th particle
   - getvz(i=0): returns history of vz for i'th particle
   - getgi(i=0): returns history of gamma inverse for i'th particle
   - getr(i=0):  returns history of r for i'th particle

   if lsavefields is true
   - getex(i=0):  returns history of Ex for i'th particle
   - getey(i=0):  returns history of Ey for i'th particle
   - getez(i=0):  returns history of Ez for i'th particle
   - getbx(i=0):  returns history of Bx for i'th particle
   - getby(i=0):  returns history of By for i'th particle
   - getbz(i=0):  returns history of Bz for i'th particle

   - pxt(i=0), pyt(i=0), prt(i=0), pzt(i=0)
     pvxt(i=0), pvyt(i=0), pvzt(i=0), pgit(i=0)
     pxy(i=0), pzx(i=0), pzy(i=0), pzr(i=0), pzvx(i=0), pzvy(i=0), pzvz(i=0)
     plots pairs of quantities for the i'th particle. All of these also
     take the same optional arguments as the plg command. If the keyword
     'trajectory' is set to false, then only the most recent particle position is
     plotted.

   - disable(): Clears out the particles. The history data is maintained, but
                the particles are no longer advanced.
   - enable(): Re-enables particles which have been disabled. They will
               continue from the place they were disabled.
    """
    # --- Class attribute to keep track of species instances
    # --- Needed so that ins and nps can be set properly the first
    # --- time a species is used.
    _instance_dict = {}

    # ----------------------------------------------------------------------
    def __init__(self, x=0., y=0., z=0., vx=0., vy=0., vz=None,
                 maxsteps=1000, savedata=1, js=0,
                 enforceinitboundaries=true, lsavefields=false):
        assert js < top.ns, "species must be already existing"

        # --- store value of allspecl
        self.allspecl = top.allspecl

        # --- Save some if the input
        self.js = js
        self.savedata = savedata
        self.enabled = 0
        self.enforceinitboundaries = enforceinitboundaries
        self.lsavefields = lsavefields

        # --- Check if species needs to be setup
        if js not in TraceParticle._instance_dict:
            TraceParticle._instance_dict[js] = 1
            if top.pgroup.sq[js] == 0.:
                top.pgroup.sq[js] = top.zion*top.echarge
            if top.pgroup.sm[js] == 0.:
                top.pgroup.sm[js] = top.aion*top.amu
            if top.pgroup.sid[js] == -1:
                top.pgroup.sid[js] = js

        # --- Use the particle's ssn to keep track of them
        if top.ssnpid == 0:
            top.ssnpid = nextpid()
        setuppgroup(top.pgroup)

        # --- Setup particles
        self.spinit(x, y, z, vx, vy, vz)
        # --- Setup history arrays
        self.setuphistory(maxsteps)

    # ----------------------------------------------------------------------
    def __del__(self):
        try:
            # --- If this is happening when python is quitting, WARP packages
            # --- may not exist anymore and errors would happen.
            self.disable()
        except:
            pass

    # ----------------------------------------------------------------------
    # --- Initialize the single particle.
    def spinit(self, x=0., y=0., z=0., vx=0., vy=0., vz=None,
               maxsteps=1000):
        "Initializes one of more particles to run independently"
        # --- Set default value of vz and make sure it is not zero
        if vz is None:
            vz = top.vbeam
        # --- Make sure that the coordinates are not of list or tuple type
        if isinstance(x, (list, tuple)):
            x = array(x)
        if isinstance(y, (list, tuple)):
            y = array(y)
        if isinstance(z, (list, tuple)):
            z = array(z)
        if isinstance(vx, (list, tuple)):
            vx = array(vx)
        if isinstance(vy, (list, tuple)):
            vy = array(vy)
        if isinstance(vz, (list, tuple)):
            vz = array(vz)
        # --- Find number of particles
        self.nn = 1
        if isinstance(x, ndarray):
            self.nn = len(x)
        if isinstance(y, ndarray):
            self.nn = len(y)
        if isinstance(z, ndarray):
            self.nn = len(z)
        if isinstance(vx, ndarray):
            self.nn = len(vx)
        if isinstance(vy, ndarray):
            self.nn = len(vy)
        if isinstance(vz, ndarray):
            self.nn = len(vz)
        # --- Store the starting values
        self.xinit = x
        self.yinit = y
        self.zinit = z
        self.vxinit = vx
        self.vyinit = vy
        self.vzinit = vz
        # --- Set gamma inverse
        if top.lrelativ:
            vsqinit = self.vxinit**2+self.vyinit**2+self.vzinit**2
            self.giinit = sqrt(1. - vsqinit/clight**2)
        else:
            self.giinit = 1.
        self.uxinit = self.vxinit/self.giinit
        self.uyinit = self.vyinit/self.giinit
        self.uzinit = self.vzinit/self.giinit
        self.pidinit = zeros((self.nn, top.npid), 'd')
        # --- Create new ssn's, making sure that the values are the same across
        # --- all processors.
        if me == 0:
            self.ssn = top.ssn + iota(self.nn)
            top.ssn = top.ssn + self.nn
        else:
            self.ssn = None
        self.ssn = warp_parallel.broadcast(self.ssn)
        self.pidinit[:,top.ssnpid-1] = self.ssn
        # --- Set current values
        self.x = self.xinit*ones(self.nn)
        self.y = self.yinit*ones(self.nn)
        self.z = self.zinit*ones(self.nn)
        self.ux = self.uxinit*ones(self.nn)
        self.uy = self.uyinit*ones(self.nn)
        self.uz = self.uzinit*ones(self.nn)
        self.gi = self.giinit*ones(self.nn)
        self.pid = self.pidinit

        # --- Initialize particle coordinates
        self.enable()
        # self.checklive()

    # ----------------------------------------------------------------------
    def enable(self):
        """Load data into fortran arrays"""
        if self.enabled:
            return
        self.enabled = 1
        self.startit = top.it
        # --- load the data
        addparticles(x=self.x, y=self.y, z=self.z, vx=self.ux, vy=self.uy, vz=self.uz,
                     gi=self.gi, pid=self.pid, js=self.js, lmomentum=true,
                     lallindomain=not self.enforceinitboundaries)
        # --- Enforce the transverse particle boundary conditions
        if self.enforceinitboundaries:
            particleboundaries3d(top.pgroup, self.js, true)
        # self.checklive()
        # --- Add routine after step to save data
        installafterstep(self.spsavedata)

    # ----------------------------------------------------------------------
    def disable(self):
        """Disables the particles"""
        if not self.enabled:
            return
        self.enabled = 0
        # --- Save the last value in case the particles are re-enabled
        for i in range(self.nn):
            ii = selectparticles(js=self.js, ssn=self.ssn[i])
            if len(ii) == 0:
                self.live[i] = 0
                self.uz[i] = 0.
                continue
            self.x[i] = getx(js=self.js, ii=ii)
            self.y[i] = gety(js=self.js, ii=ii)
            self.z[i] = getz(js=self.js, ii=ii)
            self.ux[i] = getux(js=self.js, ii=ii)
            self.uy[i] = getuy(js=self.js, ii=ii)
            self.uz[i] = getuz(js=self.js, ii=ii)
            self.gi[i] = getgaminv(js=self.js, ii=ii)
            self.pid[i, :] = getpid(js=self.js, ii=ii, id=-1)[0,:]
            # --- Set gaminv to zero, signal of dead particles
            top.pgroup.gaminv[ii] = 0.
        # --- Clear out the tracer particles
        clearpart(top.pgroup, self.js+1, 1)
        # --- remove routine from after step
        uninstallafterstep(self.spsavedata)

    # ----------------------------------------------------------------------
    def checklive(self):
        """Check which particles are still alive. Note that this refers directly
    to the WARP particle database since the particle's data is not saved if it
    is not alive."""
        self.live = zeros(self.nn, 'l')
        for i in range(self.nn):
            ii = selectparticles(js=self.js, ssn=self.ssn[i])
            # --- If the particle is not live, then there is no particle with
            # --- that ssn.
            if len(ii) > 0:
                self.live[i] = 1

    # ----------------------------------------------------------------------
    def reset(self, clearhistory=0):
        """Reset back to the starting conditions.
      - clearhistory=0: when true, the history data is cleared.
        """
        # --- The disable and enable removes the particles from the fortran
        # --- arrays and reinstalls the initial values.
        self.disable()
        # --- Set current values to the initial values
        self.x = self.xinit*ones(self.nn)
        self.y = self.yinit*ones(self.nn)
        self.z = self.zinit*ones(self.nn)
        self.ux = self.uxinit*ones(self.nn)
        self.uy = self.uyinit*ones(self.nn)
        self.uz = self.uzinit*ones(self.nn)
        self.gi = self.giinit*ones(self.nn)
        self.pid = self.pidinit
        self.enable()
        if self.savedata and clearhistory:
            self.setuphistory(maxsteps=len(self.spt[0]))
        else:
            self.spsavedata()

    # ----------------------------------------------------------------------
    def setuphistory(self, maxsteps=1000):
        # --- Create arrays for saving trajectory
        if self.savedata:
            self.spt = []
            self.spx = []
            self.spy = []
            self.spz = []
            self.spvx = []
            self.spvy = []
            self.spvz = []
            self.spgi = []
            if package()[0] == 'wxy':
                self.spdt = []
            if self.lsavefields:
                self.spex = []
                self.spey = []
                self.spez = []
                self.spbx = []
                self.spby = []
                self.spbz = []
            for i in range(self.nn):
                self.spt.append(AppendableArray(maxsteps, typecode='d'))
                self.spx.append(AppendableArray(maxsteps, typecode='d'))
                self.spy.append(AppendableArray(maxsteps, typecode='d'))
                self.spz.append(AppendableArray(maxsteps, typecode='d'))
                self.spvx.append(AppendableArray(maxsteps, typecode='d'))
                self.spvy.append(AppendableArray(maxsteps, typecode='d'))
                self.spvz.append(AppendableArray(maxsteps, typecode='d'))
                self.spgi.append(AppendableArray(maxsteps, typecode='d'))
                if package()[0] == 'wxy':
                    self.spdt.append(AppendableArray(maxsteps, typecode='d'))
                if self.lsavefields:
                    self.spex.append(AppendableArray(maxsteps, typecode='d'))
                    self.spey.append(AppendableArray(maxsteps, typecode='d'))
                    self.spez.append(AppendableArray(maxsteps, typecode='d'))
                    self.spbx.append(AppendableArray(maxsteps, typecode='d'))
                    self.spby.append(AppendableArray(maxsteps, typecode='d'))
                    self.spbz.append(AppendableArray(maxsteps, typecode='d'))
            self.spsavedata()

    # ----------------------------------------------------------------------
    def spsavedata(self):
        """Saves data"""
        self.checklive()
        if not self.savedata:
            return
        # --- set top.allspecl to true to ensure that v and x are synchronized at next time step
        if (top.it+1 - self.startit) % self.savedata == 0:
            self.allspecl = top.allspecl
            top.allspecl = 1

        if (top.it - self.startit) % self.savedata != 0:
            return
        top.allspecl = self.allspecl
        for i in range(self.nn):
            ii = selectparticles(js=self.js, ssn=self.ssn[i])
            np = globalsum(len(ii))
            if np > 0:
                self.spt[i].append(top.time)
                self.spx[i].append(getx(js=self.js, ii=ii))
                self.spy[i].append(gety(js=self.js, ii=ii))
                self.spz[i].append(getz(js=self.js, ii=ii))
                self.spvx[i].append(getux(js=self.js, ii=ii))
                self.spvy[i].append(getuy(js=self.js, ii=ii))
                self.spvz[i].append(getuz(js=self.js, ii=ii))
                self.spgi[i].append(getgaminv(js=self.js, ii=ii))
                if package()[0] == 'wxy':
                    self.spdt[i].append(getpid(js=self.js, ii=ii, id=wxy.dtpid-1))
                if self.lsavefields:
                    self.spex[i].append(getex(js=self.js, ii=ii))
                    self.spey[i].append(getey(js=self.js, ii=ii))
                    self.spez[i].append(getez(js=self.js, ii=ii))
                    self.spbx[i].append(getbx(js=self.js, ii=ii))
                    self.spby[i].append(getby(js=self.js, ii=ii))
                    self.spbz[i].append(getbz(js=self.js, ii=ii))

    # ----------------------------------------------------------------------
    def getsavedata(self):
        "Retrieves the save particle trajectories on a single array"
        if self.nn == 1:
            return array((self.spx[0].data(), self.spy[0].data(), self.spz[0].data(),
                         self.spvx[0].data(), self.spvy[0].data(), self.spvz[0].data(),
                         self.spgi[0].data()))
        else:
            r = zeros((self.nn, 7, len(self.spx[0])), 'd')
            for i in range(self.nn):
                r[i,0,:] = self.spx[i].data()
                r[i,1,:] = self.spy[i].data()
                r[i,2,:] = self.spz[i].data()
                r[i,3,:] = self.spvx[i].data()
                r[i,4,:] = self.spvy[i].data()
                r[i,5,:] = self.spvz[i].data()
                r[i,6,:] = self.spgi[i].data()
            return r

    # ----------------------------------------------------------------------
    def gett(self, i=0):
        return self.spt[i].data()

    def getx(self, i=0):
        return self.spx[i].data()

    def gety(self, i=0):
        return self.spy[i].data()

    def getz(self, i=0):
        return self.spz[i].data()

    def getux(self, i=0):
        return self.spvx[i].data()

    def getuy(self, i=0):
        return self.spvy[i].data()

    def getuz(self, i=0):
        return self.spvz[i].data()

    def getvx(self, i=0):
        return self.spvx[i].data()*self.spgi[i].data()

    def getvy(self, i=0):
        return self.spvy[i].data()*self.spgi[i].data()

    def getvz(self, i=0):
        return self.spvz[i].data()*self.spgi[i].data()

    def getgi(self, i=0):
        return self.spgi[i].data()

    def getdt(self, i=0):
        return self.spdt[i].data()

    def getr(self, i=0):
        return sqrt(self.getx(i)**2 + self.gety(i)**2)

    def getex(self, i=0):
        if not self.lsavefields:
            raise NameError('lsavefields must be set for the fields to be saved')
        return self.spex[i].data()

    def getey(self, i=0):
        if not self.lsavefields:
            raise NameError('lsavefields must be set for the fields to be saved')
        return self.spey[i].data()

    def getez(self, i=0):
        if not self.lsavefields:
            raise NameError('lsavefields must be set for the fields to be saved')
        return self.spez[i].data()

    def getbx(self, i=0):
        if not self.lsavefields:
            raise NameError('lsavefields must be set for the fields to be saved')
        return self.spbx[i].data()

    def getby(self, i=0):
        if not self.lsavefields:
            raise NameError('lsavefields must be set for the fields to be saved')
        return self.spby[i].data()

    def getbz(self, i=0):
        if not self.lsavefields:
            raise NameError('lsavefields must be set for the fields to be saved')
        return self.spbz[i].data()

    # ----------------------------------------------------------------------
    def plotparticle(self, y, x, kw):
        if kw.pop('trajectory', True):
            plg(y, x, **kw)
        else:
            if len(x) > 0:
                plp(y[-1], x[-1], **kw)

    # ----------------------------------------------------------------------
    def pxt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "time (s)", "x (m)")
        self.plotparticle(self.getx(i), self.gett(i), kw)

    def pyt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "time (s)", "y (m)")
        self.plotparticle(self.gety(i), self.gett(i), kw)

    def prt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "time (s)", "r (m)")
        self.plotparticle(self.getr(i), self.gett(i), kw)

    def pzt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "time (s)", "z (m)")
        self.plotparticle(self.getz(i), self.gett(i), kw)

    def pvxt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "time (s)", "Vx (m/s)")
        self.plotparticle(self.getvx(i), self.gett(i), kw)

    def pvyt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "time (s)", "Vy (m/s)")
        self.plotparticle(self.getvy(i), self.gett(i), kw)

    def pvzt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "time (s)", "Vz (m/s)")
        self.plotparticle(self.getvz(i), self.gett(i), kw)

    def pgit(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "time (s)", "gamma inverse")
        self.plotparticle(self.getgi(i), self.gett(i), kw)

    def pxy(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "x (m)", "y (m)")
        self.plotparticle(self.gety(i), self.getx(i), kw)

    def pzx(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "z (m)", "x (m)")
        self.plotparticle(self.getx(i), self.getz(i), kw)

    def pzy(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "z (m)", "y (m)")
        self.plotparticle(self.gety(i), self.getz(i), kw)

    def pzr(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "z (m)", "r (m)")
        self.plotparticle(self.getr(i), self.getz(i), kw)

    def pzvx(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "z (m)", "Vx (m/s)")
        self.plotparticle(self.getvx(i), self.getz(i), kw)

    def pzvy(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "z (m)", "Vy (m/s)")
        self.plotparticle(self.getvy(i), self.getz(i), kw)

    def pzvz(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Trace particle", "z (m)", "Vz (m/s)")
        self.plotparticle(self.getvz(i), self.getz(i), kw)


class NoninteractingParticles(TraceParticle):
    """
  Class for running WARP with non-interacting particles - no self-fields are
  done and no diagnostic moments are calculated. Creator arguments...
   - x,y,z,vx,vy,vz: 6 coordinates. They can either be
       scalars (for one particle) or sequences (for multiple particles). The
       sequences must all be of the same length. If some are sequences and some
       scalars, the scalar values will be used for all particles. If vz
       is not specified, all the particles will use top.vbeam.
   - maxsteps=1000: an estimate of the number of steps taken. OK if value too
                    small
   - savedata=1: frequency (in time steps) for saving of particle trajectories
   - zerophi=0: when true, w3d.phi is zero out
   - resettime=0: when true, time and beam frame location reset to initial values
   - js=0: species of particles
   - enforceinitboundaries=true: when true, the particle boundaries are
                                 enforced for the initial particle positions.
                                 It is a good idea to turn this off if the
                                 trace particles are setup before a generate
                                 is done.
   - lsavefields=0: when true, save the time histories of the field data

  Available methods...
   - gett(i=0):  returns history of time for i'th particle
   - getx(i=0):  returns history of x for i'th particle
   - gety(i=0):  returns history of y for i'th particle
   - getz(i=0):  returns history of z for i'th particle
   - getvx(i=0): returns history of vx for i'th particle
   - getvy(i=0): returns history of vy for i'th particle
   - getvz(i=0): returns history of vz for i'th particle
   - getgi(i=0): returns history of gamma inverse for i'th particle

   if lsavefields is true
   - getex(i=0):  returns history of Ex for i'th particle
   - getey(i=0):  returns history of Ey for i'th particle
   - getez(i=0):  returns history of Ez for i'th particle
   - getbx(i=0):  returns history of Bx for i'th particle
   - getby(i=0):  returns history of By for i'th particle
   - getbz(i=0):  returns history of Bz for i'th particle

   - pxt(i=0), pyt(i=0), pzt(i=0), pvxt(i=0), pvyt(i=0), pvzt(i=0), pgit(i=0)
     pxy(i=0), pzx(i=0), pzy(i=0), pzvx(i=0), pzvy(i=0), pzvz(i=0)
     plots pairs of quantities for the i'th particle. All of these also
     take the same optional arguments as the plg command.

   - disable(): Clears out the particles. The history data is maintained, but
                the particles are no longer advanced.
   - enable(): Re-enables particles which have been disabled. They will
               continue from the place they were disabled.
    """

    # ----------------------------------------------------------------------
    def __init__(self, x=0., y=0., z=0., vx=0., vy=0., vz=None,
                 maxsteps=1000, savedata=1, zerophi=0, resettime=0, js=0,
                 enforceinitboundaries=true, lsavefields=false):
        if js not in TraceParticle._instance_dict:
            TraceParticle._instance_dict[js] = 1
            if js == 0:
                top.pgroup.ins[js] = 1
            else:
                top.pgroup.ins[js] = top.pgroup.ins[js-1] + top.pgroup.nps[js-1]
            top.pgroup.nps[js] = 0
            top.pgroup.sq[js] = top.zion*top.echarge
            top.pgroup.sm[js] = top.aion*top.amu
            top.pgroup.sw[js] = 0.
        TraceParticle.__init__(self, x, y, z, vx, vy, vz, maxsteps, savedata, js,
                               enforceinitboundaries, lsavefields)
        # --- Do some initialization
        self.spsetup(zerophi)
        # --- Setup the lattice
        resetlat()
        setlatt()
        # --- Reset time if requested
        if resettime:
            self.resettime()

    # ----------------------------------------------------------------------
    def spsetup(self, zerophi=0):
        """Sets up for a single particle run: turns diagnostics off and saves some
    initial data.
      - zerophi=0: when true, zeros out the w3d.phi array"""
        # --- Turn off all of the diagnostics
        wxy.ldiag = 0
        top.ifzmmnt = 0
        top.itmomnts = 0
        top.itplps = 0
        top.itplfreq = 0
        top.zzmomnts = 0
        top.zzplps = 0
        top.zzplfreq = 0
        top.nhist = top.nt
        top.iflabwn = 0
        w3d.lrhodia3d = false
        w3d.lgetese3d = false
        w3d.lgtlchg3d = false
        # --- Turn off the field solvers
        top.fstype = -1
        top.bfstype = -1
        # --- Zero out phi if requested.
        if zerophi:
            w3d.phi = 0.
        # --- Turn off charge deposition. The laccumulate_rho is set to true
        # --- to turn off the zeroing of rho.
        top.depos = 'none'
        top.laccumulate_rho = true
        # --- Save initial time-step data
        self.itsave = top.it
        self.timesave = top.time
        self.zbeamsave = top.zbeam
        self.zgridsave = top.zgrid
        self.zgridprvsave = top.zgridprv
        self.vbeamsave = top.vbeam
        self.vbeamfrmsave = top.vbeamfrm
        self.dtsave = top.dt
        self.dssave = wxy.ds

    # ----------------------------------------------------------------------
    def resettime(self):
        # --- Initialize time stepping
        top.it = self.itsave
        top.time = self.timesave
        top.zbeam = self.zbeamsave
        top.zgrid = self.zgridsave
        top.zgridprv = self.zgridprvsave
        top.vbeam = self.vbeamsave
        top.vbeamfrm = self.vbeamfrmsave
        top.dt = self.dtsave
        wxy.ds = self.dssave
        setlatt()

    # ----------------------------------------------------------------------
    def reset(self, clearhistory=0):
        """Reset back to the starting conditions.
      - clearhistory=0: when true, the history data is cleared.
        """
        self.__class__.reset(clearhistory)
        self.resettime()


##############################################################################
# Keep original until new classes have been tested.
class SingleParticle:
    """
  Class for running WARP in a single particle mode - no self-fields are
  done and no diagnostic moments are calculated. Creator arguments...
   - x,y,z,vx,vy,vz: 6 coordinates. They can either be
       scalars (for one particle) or sequences (for multiple particles). The
       sequences must all be of the same length. If some are sequences and some
       scalars, the scalar values will be used for all particles. If vz
       is not specified, all the particles will use top.vbeam.
   - maxsteps=1000: an estimate of the number of steps taken. OK if value too
                    small
   - savedata=1: frequency (in time steps) for saving of particle trajectories
   - zerophi=0: when true, w3d.phi is zero out
   - resettime=0: when true, time and beam frame location reset to initial values
   - js=0: species of particles

  Available methods...
   - gett(i=0):  returns history of time for i'th particle
   - getx(i=0):  returns history of x for i'th particle
   - gety(i=0):  returns history of y for i'th particle
   - getz(i=0):  returns history of z for i'th particle
   - getvx(i=0): returns history of vx for i'th particle
   - getvy(i=0): returns history of vy for i'th particle
   - getvz(i=0): returns history of vz for i'th particle
   - getgi(i=0): returns history of gamma inverse for i'th particle

   - pxt(i=0), pyt(i=0), pzt(i=0), pvxt(i=0), pvyt(i=0), pvzt(i=0), pgit(i=0)
     pxy(i=0), pzx(i=0), pzy(i=0), pzvx(i=0), pzvy(i=0), pzvz(i=0)
     plots pairs of quantities for the i'th particle. All of these also
     take the same optional arguments as the plg command.

   - disable(): Clears out the particles. The history data is maintained, but
                the particles are no longer advanced.
   - enable(): Re-enables particles which have been disabled. They will
               continue from the place they were disabled.
    """
    # --- Class attribute to keep track of species instances
    # --- Needed so that ins and nps can be set properly the first
    # --- time a species is used.
    _instance_dict = {}

    # ----------------------------------------------------------------------
    def __init__(self, x=0., y=0., z=0., vx=0., vy=0., vz=None,
                 maxsteps=1000, savedata=1, zerophi=0, resettime=0, js=0):
        # --- Do some global initialization

        # --- store value of allspecl
        self.allspecl = top.allspecl

        self.js = js
        if js not in SingleParticle._instance_dict:
            SingleParticle._instance_dict[js] = 1
            if js == 0:
                top.pgroup.ins[js] = 1
            else:
                top.pgroup.ins[js] = top.pgroup.ins[js-1] + top.pgroup.nps[js-1]
            top.pgroup.nps[js] = 0
            top.pgroup.sq[js] = top.zion*top.echarge
            top.pgroup.sm[js] = top.aion*top.amu
            top.pgroup.sw[js] = 0.
            if top.pgroup.sid[js] == -1:
                top.pgroup.sid[js] = js
        self.savedata = savedata
        self.enabled = 0
        # --- Do some initialization
        self.spsetup(zerophi)
        # --- Setup particles
        self.spinit(x, y, z, vx, vy, vz)
        # --- Setup the lattice
        resetlat()
        setlatt()
        # --- Reset time if requested
        if resettime:
            self.resettime()
        # --- Setup history arrays
        self.setuphistory(maxsteps)

    # ----------------------------------------------------------------------
    def __del__(self):
        try:
            # --- If this is happening when python is quitting, WARP packages
            # --- may not exist anymore and errors would happen.
            self.disable()
        except:
            pass

    # ----------------------------------------------------------------------
    def spsetup(self, zerophi=0):
        """Sets up for a single particle run: turns diagnostics off and saves some
    initial data.
      - zerophi=0: when true, zeros out the w3d.phi array"""
        # --- Turn off all of the diagnostics
        wxy.ldiag = 0
        top.ifzmmnt = 0
        top.itmomnts = 0
        top.itplps = 0
        top.itplfreq = 0
        top.zzmomnts = 0
        top.zzplps = 0
        top.zzplfreq = 0
        top.nhist = top.nt
        top.iflabwn = 0
        w3d.lrhodia3d = false
        w3d.lgetese3d = false
        w3d.lgtlchg3d = false
        # --- Turn off the field solvers
        top.fstype = -1
        top.bfstype = -1
        # --- Zero out phi if requested.
        if zerophi:
            w3d.phi = 0.
        # --- Turn off charge deposition. The laccumulate_rho is set to true
        # --- to turn off the zeroing of rho.
        top.depos = 'none'
        top.laccumulate_rho = true
        # --- Save initial time-step data
        self.itsave = top.it
        self.timesave = top.time
        self.zbeamsave = top.zbeam
        self.zgridsave = top.zgrid
        self.zgridprvsave = top.zgridprv
        self.vbeamsave = top.vbeam
        self.vbeamfrmsave = top.vbeamfrm
        self.dtsave = top.dt
        self.dssave = wxy.ds

    # ----------------------------------------------------------------------
    # --- Initialize the single particle.
    def spinit(self, x=0., y=0., z=0., vx=0., vy=0., vz=None,
               maxsteps=1000):
        "Initializes one of more particles to run independently"
        # --- Set default value of vz and make sure it is not zero
        if vz is None:
            vz = top.vbeam
        # --- Make sure that the coordinates are not of list or tuple type
        if isinstance(x, (list, tuple)):
            x = array(x)
        if isinstance(y, (list, tuple)):
            y = array(y)
        if isinstance(z, (list, tuple)):
            z = array(z)
        if isinstance(vx, (list, tuple)):
            vx = array(vx)
        if isinstance(vy, (list, tuple)):
            vy = array(vy)
        if isinstance(vz, (list, tuple)):
            vz = array(vz)
        # --- Find number of particles
        self.nn = 1
        if isinstance(x, ndarray):
            self.nn = len(x)
        if isinstance(y, ndarray):
            self.nn = len(y)
        if isinstance(z, ndarray):
            self.nn = len(z)
        if isinstance(vx, ndarray):
            self.nn = len(vx)
        if isinstance(vy, ndarray):
            self.nn = len(vy)
        if isinstance(vz, ndarray):
            self.nn = len(vz)
        # --- Store the starting values
        self.xinit = x
        self.yinit = y
        self.zinit = z
        self.vxinit = vx
        self.vyinit = vy
        self.vzinit = vz
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        # --- Set gamma inverse
        if top.lrelativ:
            self.gi = sqrt(1. - (self.vx**2+self.vy**2+self.vz**2)/clight**2)
            gamma = 1./self.gi
            self.vx = self.vx*gamma
            self.vy = self.vy*gamma
            self.vz = self.vz*gamma
        else:
            self.gi = self.nn*[1.]
        # --- Initialize particle coordinates
        self.enable()
        self.checklive()

    # ----------------------------------------------------------------------
    def enable(self):
        """Load data into fortran arrays"""
        if self.enabled:
            return
        self.enabled = 1
        self.startit = top.it
        # --- make sure there is space
        chckpart(top.pgroup, self.js+1, 0, self.nn)
        # --- load the data
        ip1 = top.pgroup.ins[self.js] - 1 + top.pgroup.nps[self.js]
        top.pgroup.nps[self.js] = top.pgroup.nps[self.js] + self.nn
        ip2 = top.pgroup.ins[self.js] + top.pgroup.nps[self.js] - 1
        self.ip1 = ip1
        self.ip2 = ip2
        top.pgroup.xp[ip1:ip2] = self.x
        top.pgroup.yp[ip1:ip2] = self.y
        top.pgroup.zp[ip1:ip2] = self.z
        top.pgroup.uxp[ip1:ip2] = self.vx
        top.pgroup.uyp[ip1:ip2] = self.vy
        top.pgroup.uzp[ip1:ip2] = self.vz
        top.pgroup.gaminv[ip1:ip2] = self.gi
        # --- Enforce the particle boundary conditions
        particleboundaries3d(top.pgroup, self.js, true)
        # --- Add routine after step to save data
        installafterstep(self.spsavedata)

    # ----------------------------------------------------------------------
    def disable(self):
        """Disables the particles"""
        if not self.enabled:
            return
        self.enabled = 0
        # --- Save the last value in case the particles are re-enabled
        self.x = top.pgroup.xp[self.ip1:self.ip2] + 0.
        self.y = top.pgroup.yp[self.ip1:self.ip2] + 0.
        self.z = top.pgroup.zp[self.ip1:self.ip2] + 0.
        self.vx = top.pgroup.uxp[self.ip1:self.ip2] + 0.
        self.vy = top.pgroup.uyp[self.ip1:self.ip2] + 0.
        self.vz = top.pgroup.uzp[self.ip1:self.ip2] + 0.
        self.gi = top.pgroup.gaminv[self.ip1:self.ip2] + 0.
        # --- Set gaminv to zero, signal of dead particles
        top.pgroup.gaminv[self.ip1:self.ip2] = 0.
        # --- Remove particles from looping if they are at either end of the
        # --- particle data in the arrays
        if self.ip1 == top.pgroup.ins[self.js] - 1:
            top.pgroup.ins[self.js] = top.pgroup.ins[self.js] + self.nn
            top.pgroup.nps[self.js] = top.pgroup.nps[self.js] - self.nn
        elif self.ip2 == top.pgroup.ins[self.js] + top.pgroup.nps[self.js] - 1:
            top.pgroup.nps[self.js] = top.pgroup.nps[self.js] - self.nn
        # --- remove routine from after step
        uninstallafterstep(self.spsavedata)

    # ----------------------------------------------------------------------
    def checklive(self):
        """Check which particles are still alive. Note that this refers directly
    to the WARP particle database since the particle's data is not saved if it
    is not alive."""
        self.live = ones(self.ip2 - self.ip1+1)

    # ----------------------------------------------------------------------
    def resettime(self):
        # --- Initialize time stepping
        top.it = self.itsave
        top.time = self.timesave
        top.zbeam = self.zbeamsave
        top.zgrid = self.zgridsave
        top.zgridprv = self.zgridprvsave
        top.vbeam = self.vbeamsave
        top.vbeamfrm = self.vbeamfrmsave
        top.dt = self.dtsave
        wxy.ds = self.dssave
        setlatt()

    # ----------------------------------------------------------------------
    def reset(self, clearhistory=0):
        """Reset back to the starting conditions.
      - clearhistory=0: when true, the history data is cleared.
        """
        # --- The disable and enable removes the particles from the fortran
        # --- arrays and reinstalls the initial values.
        self.disable()
        self.x = self.xinit
        self.y = self.yinit
        self.z = self.zinit
        self.vx = self.vxinit
        self.vy = self.vyinit
        self.vz = self.vzinit
        self.resettime()
        self.enable()
        if self.savedata and clearhistory:
            self.setuphistory(maxsteps=len(self.spt[0]))
        else:
            self.spsavedata()

    # ----------------------------------------------------------------------
    def setuphistory(self, maxsteps=1000):
        # --- Create arrays for saving trajectory
        if self.savedata:
            self.spt = []
            self.spx = []
            self.spy = []
            self.spz = []
            self.spvx = []
            self.spvy = []
            self.spvz = []
            self.spgi = []
            if package()[0] == 'wxy':
                self.spdt = []
            for i in range(self.nn):
                self.spt.append(AppendableArray(maxsteps, typecode='d'))
                self.spx.append(AppendableArray(maxsteps, typecode='d'))
                self.spy.append(AppendableArray(maxsteps, typecode='d'))
                self.spz.append(AppendableArray(maxsteps, typecode='d'))
                self.spvx.append(AppendableArray(maxsteps, typecode='d'))
                self.spvy.append(AppendableArray(maxsteps, typecode='d'))
                self.spvz.append(AppendableArray(maxsteps, typecode='d'))
                self.spgi.append(AppendableArray(maxsteps, typecode='d'))
                if package()[0] == 'wxy':
                    self.spdt.append(AppendableArray(maxsteps, typecode='d'))
            self.spsavedata()

    # ----------------------------------------------------------------------
    def spsavedata(self):
        """Saves data"""
        self.checklive()
        if not self.savedata:
            return
        # --- set top.allspecl to true to ensure that v and x are synchronized at next time step
        if (top.it+1 - self.startit) % self.savedata == 0:
            self.allspecl = top.allspecl
            top.allspecl = 1

        if (top.it - self.startit) % self.savedata != 0:
            return
        top.allspecl = self.allspecl
        for i in range(self.nn):
            if self.live[i]:
                self.spt[i].append(top.time)
                self.spx[i].append(top.pgroup.xp[self.ip1 + i])
                self.spy[i].append(top.pgroup.yp[self.ip1 + i])
                self.spz[i].append(top.pgroup.zp[self.ip1 + i])
                self.spvx[i].append(top.pgroup.uxp[self.ip1 + i])
                self.spvy[i].append(top.pgroup.uyp[self.ip1 + i])
                self.spvz[i].append(top.pgroup.uzp[self.ip1 + i])
                self.spgi[i].append(top.pgroup.gaminv[self.ip1 + i])
                if package()[0] == 'wxy':
                    self.spdt[i].append(top.pgroup.pid[self.ip1 + i, wxy.dtpid-1])

    # ----------------------------------------------------------------------
    def getsavedata(self):
        "Retrieves the save particle trajectories on a single array"
        if self.nn == 1:
            return array((self.spx[0].data(), self.spy[0].data(), self.spz[0].data(),
                          self.spvx[0].data(), self.spvy[0].data(), self.spvz[0].data(),
                          self.spgi[0].data()))
        else:
            r = zeros((self.nn, 7, len(self.spx[0])), 'd')
            for i in range(self.nn):
                r[i,0,:] = self.spx[i].data()
                r[i,1,:] = self.spy[i].data()
                r[i,2,:] = self.spz[i].data()
                r[i,3,:] = self.spvx[i].data()
                r[i,4,:] = self.spvy[i].data()
                r[i,5,:] = self.spvz[i].data()
                r[i,6,:] = self.spgi[i].data()
            return r
 
    # ----------------------------------------------------------------------
    def gett(self, i=0):
        return self.spt[i].data()

    def getx(self, i=0):
        return self.spx[i].data()

    def gety(self, i=0):
        return self.spy[i].data()

    def getz(self, i=0):
        return self.spz[i].data()

    def getvx(self, i=0):
        return self.spvx[i].data()*self.spgi[i].data()

    def getvy(self, i=0):
        return self.spvy[i].data()*self.spgi[i].data()

    def getvz(self, i=0):
        return self.spvz[i].data()*self.spgi[i].data()

    def getgi(self, i=0):
        return self.spgi[i].data()

    def getdt(self, i=0):
        return self.spdt[i].data()

    def getr(self, i=0):
        return sqrt(self.getx(i)**2 + self.gety(i)**2)

    # ----------------------------------------------------------------------
    def pxt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "time (s)", "x (m)")
        plg(self.getx(i), self.gett(i), **kw)

    def pyt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "time (s)", "y (m)")
        plg(self.gety(i), self.gett(i), **kw)

    def prt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "time (s)", "r (m)")
        plg(self.getr(i), self.gett(i), **kw)

    def pzt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "time (s)", "z (m)")
        plg(self.getz(i), self.gett(i), **kw)

    def pvxt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "time (s)", "Vx (m/s)")
        plg(self.getvx(i), self.gett(i), **kw)

    def pvyt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "time (s)", "Vy (m/s)")
        plg(self.getvy(i), self.gett(i), **kw)

    def pvzt(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "time (s)", "Vz (m/s)")
        plg(self.getvz(i), self.gett(i), **kw)

    def pgit(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "time (s)", "gamma inverse")
        plg(self.getgi(i), self.gett(i), **kw)

    def pxy(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "x (m)", "y (m)")
        plg(self.gety(i), self.getx(i), **kw)

    def pzx(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "z (m)", "x (m)")
        plg(self.getx(i), self.getz(i), **kw)

    def pzy(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "z (m)", "y (m)")
        plg(self.gety(i), self.getz(i), **kw)

    def pzr(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "z (m)", "r (m)")
        plg(self.getr(i), self.getz(i), **kw)

    def pzvx(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "z (m)", "Vx (m/s)")
        plg(self.getvx(i), self.getz(i), **kw)

    def pzvy(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "z (m)", "Vy (m/s)")
        plg(self.getvy(i), self.getz(i), **kw)

    def pzvz(self, i=0, **kw):
        if kw.pop('titles', True):
            ptitles("Single particle", "z (m)", "Vz (m/s)")
        plg(self.getvz(i), self.getz(i), **kw)
