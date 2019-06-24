"""
Defines the PICAPI and TimeStepper classes.
"""

class PICAPI:
    """
    PIC code API
    This specifies all of the pieces that are needed to carray out time steps
    in a typical PIC code. Code that implements this API can be controlled
    by the generic time stepping class, which allows easy switching between
    various way of organizing the time step. It also allows insertion of
    user specified controllers to perform additional actions during a time step.
    """
    def get_time(self): return 0.
    def set_time(self, time): self.time = time
    def get_step_size(self): return 1.
    def get_step_number(self): return 0
    def set_step_number(self, it): self.it = it

    def push_positions(self, dt): pass
    def push_velocities_withE(self, dt): pass
    def push_velocities_withB(self, dt): pass
    def get_self_fields(self): pass
    def get_applied_fields(self, dt=None, dtl=None, dtr=None): pass
    def calculate_source(self): pass
    def push_Efields(self, dt): pass
    def push_Bfields(self, dt): pass
    def apply_particle_boundary_conditions(self): pass


class TimeStepper(object):
    """Generic time stepping class
    This relies on instances that use the PICAPI.
    
    splitting='velocity': Species the organization of a time step.
                          'velocity' means that the step will be split so
                              that the position and velocity will be synchronized
                              after a full position advance.
                          'position' means that the step will be split so
                              that the position and velocity will be synchronized
                              after a full velocity advance.
    alwayssplit=True: When True, split steps are always done, otherwise full steps
                      are done for efficiency.
    """

    def __init__(self, PICcode, splitting='velocity', alwayssplit=True):
        self.PICcode = PICcode
        self.splitting = splitting
        self.alwayssplit = alwayssplit

    def step(self, nsteps=1):
        for i in range(nsteps):
            firststep = ((i == 0) or self.alwayssplit)
            laststep = ((i+1 == nsteps) or self.alwayssplit)
            if self.splitting == 'position':
                self.splitpositionstep(firststep, laststep)
            elif self.splitting == 'velocity':
                self.splitvelocitystep(firststep, laststep)
            else:
                raise Exception('Unsupported splitting option')

    def splitvelocitystep(self, firststep, laststep):
        """Split leap-frog with velocity advance split.
        Step ends with position and velocity synchronized at the end of a full position advance.
        """
        PICcode = self.PICcode

        time = PICcode.get_time()
        istep = PICcode.get_step_number()
        dt = PICcode.get_step_size()

        if firststep:
            PICcode.get_self_fields()
            PICcode.get_applied_fields(dt, 0., 1./2.)
            PICcode.push_velocities_withB(dt/2.)
            PICcode.push_velocities_withE(dt/2.)
        else:
            PICcode.push_velocities_withE(dt/2.)
            PICcode.push_velocities_withB(dt)
            PICcode.push_velocities_withE(dt/2.)

        PICcode.push_positions(dt)
        PICcode.apply_particle_boundary_conditions()
        PICcode.calculate_source()
        PICcode.push_Bfields(dt/2.)
        PICcode.push_Efields(dt)
        PICcode.push_Bfields(dt/2.)

        istep += 1
        time += dt

        PICcode.set_step_number(istep)
        PICcode.set_time(time)

        PICcode.get_self_fields()

        if laststep:
            PICcode.get_applied_fields(dt, -1./2., 0.)
            PICcode.push_velocities_withE(dt/2.)
            PICcode.push_velocities_withB(dt/2.)
        else:
            PICcode.get_applied_fields(dt, -1./2., 1./2.)

    def splitpositionstep(self, firststep, laststep):
        """Split leap-frog with position advance split.
        Step ends with position and velocity synchronized at the end of a full velocity advance.
        """
        PICcode = self.PICcode

        time = PICcode.get_time()
        istep = PICcode.get_step_number()
        dt = PICcode.get_step_size()

        if firststep:
            PICcode.push_positions(dt/2.)
        else:
            PICcode.push_positions(dt)

        PICcode.apply_particle_boundary_conditions()
        PICcode.calculate_source()
        PICcode.push_Bfields(dt/2.)
        PICcode.push_Efields(dt)
        PICcode.push_Bfields(dt/2.)

        istep += 1
        time += dt

        PICcode.set_step_number(istep)
        PICcode.set_time(time)

        PICcode.get_self_fields()
        PICcode.get_applied_fields(dt, -1./2., 1./2.)
        PICcode.push_velocities_withE(dt/2.)
        PICcode.push_velocities_withB(dt)
        PICcode.push_velocities_withE(dt/2.)

        if laststep:
            PICcode.push_positions(dt/2.)
            PICcode.apply_particle_boundary_conditions()
