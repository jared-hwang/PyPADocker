"""WarpPIC API for time stepping"""
from .. import warp
from ..run_modes.timestepper import PICAPI

class WarpPIC(PICAPI):

    def __init__(self, species=None):
        if species is not None:
            self.species = [species]
        else:
            self.species = warp.listofallspecies

    def get_time(self):
        return warp.top.time

    def set_time(self, time):
        warp.top.time = time

    def get_step_size(self):
        return warp.top.dt
        
    def get_step_number(self):
        return warp.top.it

    def set_step_number(self, it):
        warp.top.it = it

    def push_positions(self, dt):
        for specie in self.species:
            specie.xpush3d(dt)

    def push_velocities_withE(self, dt):
        for specie in self.species:
            specie.epush3d(dt)

    def push_velocities_withB(self, dt):
        for specie in self.species:
            specie.bpush3d(dt)

    def get_self_fields(self):
        for specie in self.species:
            specie.getselffields()

    def get_applied_fields(self, dt, dtl, dtr):
        for specie in self.species:
            specie.getappliedfields(dt, dtl, dtr)

    def calculate_source(self):
        for specie in self.species:
            specie.loadrho(lfinalize_rho=False)
        warp.finalizerho()

    def push_Efields(self, dt):
        warp.fieldsolve()

    def push_Bfields(self, dt):
        pass

    def apply_particle_boundary_conditions(self):
        for specie in self.species:
            specie.particleboundaries()

