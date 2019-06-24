"""Classes following the PICMI standard
"""
import re
import picmistandard
import numpy as np
from .warp import *
from .init_tools.boost_tools import BoostConverter
from .init_tools.plasma_initialization import PlasmaInjector
from .init_tools.beam_tools import initialize_beam_fields
from .field_solvers.em3dsolverFFT import *
from .data_dumping import openpmd_diag
import warp

codename = 'warp'
picmistandard.register_codename(codename)

c = warp.clight
ep0 = warp.eps0
mu0 = warp.mu0
q_e = warp.echarge
m_e = warp.Electron.mass
m_p = warp.Proton.mass

# --- Always assume that relativity should be turned on.
warp.top.lrelativ = warp.true

# --- Species name dictionary.
# --- This maps the PICMI type (a string) to the Warp species types
# --- This will go in a separate file, since it will be big.
species_type_dict = {'electron':warp.Electron,
                     'positron':warp.Positron,
                     'proton':warp.Proton}
# --- Add all of the elements to the dictionary
for name,d in warp.periodic_table.items():
    symbol = d['Symbol']
    m = re.search('[0-9]', symbol)
    if m is not None:
        symbol = '#' + m.group() + symbol
    species_type_dict[symbol] = getattr(warp, name)


class Species(picmistandard.PICMI_Species):
    def init(self, kw):

        # --- When using the PICMI standard, always make the particles
        # --- have variable weights. This will greatly simplify the
        # --- handling of particle initialization. The species weight
        # --- will be set to 1, and the variable weight will hold
        # --- the actual weight.

        wtype = species_type_dict.get(self.particle_type, None)
        self.wspecies = warp.Species(type = wtype,
                                     name = self.name,
                                     charge = self.charge,
                                     charge_state = self.charge_state,
                                     mass = self.mass,
                                     lvariableweights = True,
                                     weight = 1.)

    def initialize_species_inputs(self, layout):

        if self.particle_shape is not None:
            if isinstance(self.particle_shape, str):
                interpolation_order = {'NGP':0, 'linear':1, 'quadratic':2, 'cubic':3}[self.particle_shape]
            else:
                interpolation_order = self.particle_shape
            self.wspecies.depos_order = interpolation_order

        if self.initial_distribution is not None:
            self.layout = layout
            installparticleloader(self.particle_loader)

    def particle_loader(self):
        self.initial_distribution.loaddistribution(self.wspecies, self.layout)


picmistandard.PICMI_MultiSpecies.Species_class = Species
class MultiSpecies(picmistandard.PICMI_MultiSpecies):
    pass


class GaussianBunchDistribution(picmistandard.PICMI_GaussianBunchDistribution):
    def loaddistribution(self, species, layout):
        assert isinstance(layout, PseudoRandomLayout), Exception('Warp only supports PseudoRandomLayout with GaussianBunchDistribution')
        assert layout.n_macroparticles is not None, Exception('Warp only support n_macroparticles with PseudoRandomLayout with GaussianBunchDistribution')
        np = layout.n_macroparticles
        deltax = self.rms_bunch_size[0]
        deltay = self.rms_bunch_size[1]
        deltaz = self.rms_bunch_size[2]
        vthx = self.rms_velocity[0]
        vthy = self.rms_velocity[1]
        vthz = self.rms_velocity[2]
        xmean = self.centroid_position[0]
        ymean = self.centroid_position[1]
        zmean = self.centroid_position[2]
        vxmean = self.centroid_velocity[0]
        vymean = self.centroid_velocity[1]
        vzmean = self.centroid_velocity[2]
        vxdiv = self.velocity_divergence[0]
        vydiv = self.velocity_divergence[1]
        vzdiv = self.velocity_divergence[2]
        w = self.n_physical_particles/np
        species.add_gaussian_dist(np, deltax, deltay, deltaz, vthx, vthy, vthz,
                                  xmean, ymean, zmean, vxmean, vymean, vzmean, vxdiv, vydiv, vzdiv,
                                  zdist='random', rdist='linear', fourfold=False, lmomentum=True, w=w)

        if top.boost_gamma > 1.:
            # --- Code copied from Boosted_Frame class
            # --- This maps all particles from this species, which will not be
            # --- correct if some had be created elsewhere.
            boost_beta = sqrt(1. - 1./top.boost_gamma**2)
            for jspr,js in enumerate(species.jslist):
                if getn(pgroup=top.pgroup, js=js, bcast=0, gather=0)>0:
                    il = top.pgroup.ins[js] - 1
                    iu = il + top.pgroup.nps[js]
                    uzfrm = top.boost_gamma*boost_beta*warp.clight
                    tpr = top.boost_gamma*top.time - uzfrm*top.pgroup.zp[il:iu]/warp.clight**2
                    top.pgroup.zp[il:iu] = top.boost_gamma*top.pgroup.zp[il:iu] - uzfrm*top.time
                    vx = getvx(pgroup=top.pgroup, js=js, bcast=0, gather=0)
                    vy = getvy(pgroup=top.pgroup, js=js, bcast=0, gather=0)
                    vz = getvz(pgroup=top.pgroup, js=js, bcast=0, gather=0)
                    fact = 1./(1. - boost_beta*vz/warp.clight)
                    vxpr = vx*fact/top.boost_gamma
                    vypr = vy*fact/top.boost_gamma
                    vzpr = (vz - boost_beta*warp.clight)*fact
                    top.pgroup.xp[il:iu] = top.pgroup.xp[il:iu] - tpr*vxpr
                    top.pgroup.yp[il:iu] = top.pgroup.yp[il:iu] - tpr*vypr
                    top.pgroup.zp[il:iu] = top.pgroup.zp[il:iu] - tpr*vzpr
                    gammapr = 1./sqrt(1. - (vxpr*vxpr + vypr*vypr + vzpr*vzpr)/warp.clight**2)
                    top.pgroup.uxp[il:iu] = vxpr*gammapr
                    top.pgroup.uyp[il:iu] = vypr*gammapr
                    top.pgroup.uzp[il:iu] = vzpr*gammapr
                    top.pgroup.gaminv[il:iu] = 1./gammapr
                    if top.uxoldpid > 0:
                        top.pgroup.pid[il:iu,top.uxoldpid-1] = top.pgroup.uxp[il:iu]
                    if top.uyoldpid > 0:
                        top.pgroup.pid[il:iu,top.uyoldpid-1] = top.pgroup.uyp[il:iu]
                    if top.uzoldpid > 0:
                        top.pgroup.pid[il:iu,top.uzoldpid-1] = top.pgroup.uzp[il:iu]
            particleboundaries3d(top.pgroup, -1, False)


class UniformDistribution(picmistandard.PICMI_UniformDistribution):
    def loaddistribution(self, species, layout):
        xmin = self.lower_bound[0]
        if xmin is None:
            xmin = -largepos
        ymin = self.lower_bound[1]
        if ymin is None:
            ymin = -largepos
        zmin = self.lower_bound[2]
        if zmin is None:
            zmin = -largepos
        xmax = self.upper_bound[0]
        if xmax is None:
            xmax = +largepos
        ymax = self.upper_bound[1]
        if ymax is None:
            ymax = +largepos
        zmax = self.upper_bound[2]
        if zmax is None:
            zmax = +largepos
        ux_th = self.rms_velocity[0]
        uy_th = self.rms_velocity[1]
        uz_th = self.rms_velocity[2]
        ux_m = self.directed_velocity[0]
        uy_m = self.directed_velocity[1]
        uz_m = self.directed_velocity[2]

        density = self.density

        if top.boost_gamma > 1.:
            boost = BoostConverter(top.boost_gamma)
            beta = uz_m/sqrt(1. + uz_m**2)
            zmin, zmax = boost.copropag_length([zmin, zmax], beta_object=beta)
            density, = boost.copropag_density([density], beta_object=beta)
            uz_m, = boost.longitudinal_momentum([uz_m])

        if isinstance(layout, GriddedLayout):
            # --- Note that layout.grid is ignored
            p_nx = layout.n_macroparticle_per_cell[0]
            if len(layout.n_macroparticle_per_cell) > 2:
                p_ny = layout.n_macroparticle_per_cell[1]
            else:
                p_ny = 1
            p_nz = layout.n_macroparticle_per_cell[-1]

            npreal_per_cell = density*w3d.dx*w3d.dy*w3d.dz
            w = npreal_per_cell/(p_nx*p_ny*p_nz)
            def dens_func(x, y, z, w=w):
                return w

            if top.vbeamfrm > 0:
                injection_direction = +1
            else:
                injection_direction = -1
            plasmainjector = PlasmaInjector(elec=species, ions=None, w3d=w3d, top=top,
                                            dim='%dd'%layout.grid.number_of_dimensions,
                                            p_nx=p_nx, p_ny=p_ny, p_nz=p_nz,
                                            p_xmin=xmin, p_ymin=ymin, p_zmin=zmin,
                                            p_xmax=xmax, p_ymax=ymax, p_zmax=zmax,
                                            dens_func=dens_func,
                                            ux_m=ux_m, uy_m=uy_m, uz_m=uz_m,
                                            ux_th=ux_th, uy_th=uy_th, uz_th=uz_th,
                                            injection_direction=injection_direction)
            if self.fill_in:
                installuserinjection(plasmainjector.continuous_injection)

        elif isinstance(layout, PseudoRandomLayout):
            assert not self.fill_in, Exception('Warp does not support fill_in with PseudoRandomLayout')
            # For this option, the mins and maxs must be bounded by the domain
            xmin = max(xmin, w3d.xmmin)
            xmax = min(xmax, w3d.xmmax)
            ymin = max(ymin, w3d.ymmin)
            ymax = min(ymax, w3d.ymmax)
            zmin = max(zmin, w3d.zmmin)
            zmax = min(zmax, w3d.zmmax)
            # Determine the number of particles to load
            if layout.n_macroparticles_per_cell is not None:
                np = layout.n_macroparticles_per_cell*w3d.nx*w3d.ny*w3d.nz
            elif layout.n_macroparticles is not None:
                np = layout.n_macroparticles
            # The particle weight
            npreal = density*(xmax - xmin)*(ymax - ymin)*(zmax - zmin)
            w = np.full(np, npreal/np)
            species.add_uniform_box(np=np, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                                    vthx=ux_th, vthy=uy_th, vthz=uz_th,
                                    vxmean=ux_m, vymean=uy_m, vzmean=uz_m,
                                    lmomentum=1, spacing='random',
                                    lallindomain=warp.true,
                                    w=w)


class AnalyticDistribution(picmistandard.PICMI_AnalyticDistribution):
    def loaddistribution(self, species, layout):
        xmin = self.lower_bound[0]
        if xmin is None:
            xmin = -largepos
        ymin = self.lower_bound[1]
        if ymin is None:
            ymin = -largepos
        zmin = self.lower_bound[2]
        if zmin is None:
            zmin = -largepos
        xmax = self.upper_bound[0]
        if xmax is None:
            xmax = +largepos
        ymax = self.upper_bound[1]
        if ymax is None:
            ymax = +largepos
        zmax = self.upper_bound[2]
        if zmax is None:
            zmax = +largepos
        ux_th = self.rms_velocity[0]
        uy_th = self.rms_velocity[1]
        uz_th = self.rms_velocity[2]
        ux_m = self.directed_velocity[0]
        uy_m = self.directed_velocity[1]
        uz_m = self.directed_velocity[2]

        if top.boost_gamma > 1.:
            boost = BoostConverter(top.boost_gamma)
            beta = uz_m/sqrt(1. + uz_m**2)
            zmin, zmax = boost.copropag_length([zmin, zmax], beta_object=beta)
            density_boost_converter, = boost.copropag_density([1.], beta_object=beta)
            z_boost_converter, = boost.copropag_length([1.], beta_object=beta)
            uz_m, = boost.longitudinal_momentum([uz_m])
        else:
            density_boost_converter = 1.
            z_boost_converter = 1.
            
        if isinstance(layout, GriddedLayout):
            # --- Note that layout.grid is ignored
            p_nx = layout.n_macroparticle_per_cell[0]
            if len(layout.n_macroparticle_per_cell) > 2:
                p_ny = layout.n_macroparticle_per_cell[1]
            else:
                p_ny = 1
            p_nz = layout.n_macroparticle_per_cell[-1]

            if isinstance(self.density_expression, str):
                cell_volume_per_particle = w3d.dx*w3d.dy*w3d.dz/(p_nx*p_ny*p_nz)*density_boost_converter
                def dens_func(x, y, z, density=self.density_expression, cell_volume_per_particle=cell_volume_per_particle, z_boost_converter=z_boost_converter):
                    # --- Include globals so that numpy is available
                    if top.boost_gamma > 1.:
                        z = z/z_boost_converter
                    dct = locals()
                    dct.update(self.user_defined_kw)
                    d = eval(density, globals(), dct)
                    return d*cell_volume_per_particle
            else:
                npreal_per_cell = self.density_expression*w3d.dx*w3d.dy*w3d.dz*density_boost_converter
                w = npreal_per_cell/(p_nx*p_ny*p_nz)
                def dens_func(x, y, z, w=w):
                    return w

            if top.vbeamfrm > 0:
                injection_direction = +1
            else:
                injection_direction = -1
            plasmainjector = PlasmaInjector(elec=species, ions=None, w3d=w3d, top=top, dim='3d',
                                            p_nx=p_nx, p_ny=p_ny, p_nz=p_nz,
                                            p_xmin=xmin, p_ymin=ymin, p_zmin=zmin,
                                            p_xmax=xmax, p_ymax=ymax, p_zmax=zmax,
                                            dens_func=dens_func,
                                            ux_m=ux_m, uy_m=uy_m, uz_m=uz_m,
                                            ux_th=ux_th, uy_th=uy_th, uz_th=uz_th,
                                            injection_direction=injection_direction)
            if self.fill_in:
                installuserinjection(plasmainjector.continuous_injection)

        elif isinstance(layout, PseudoRandomLayout):
            assert not self.fill_in, Exception('Warp does not support fill_in with PseudoRandomLayout')
            assert not isinstance(self.density_expression, str), Exception('Warp does not support PseudoRandomLayout with nonuniform distribution')
            # For this option, the mins and maxs must be bounded by the domain
            xmin = max(xmin, w3d.xmmin)
            xmax = min(xmax, w3d.xmmax)
            ymin = max(ymin, w3d.ymmin)
            ymax = min(ymax, w3d.ymmax)
            zmin = max(zmin, w3d.zmmin)
            zmax = min(zmax, w3d.zmmax)
            # Determine the number of particles to load
            if layout.n_macroparticles_per_cell is not None:
                np = layout.n_macroparticles_per_cell*w3d.nx*w3d.ny*w3d.nz
            elif layout.n_macroparticles is not None:
                np = layout.n_macroparticles
            # The particle weight
            npreal = self.density_expression*(xmax - xmin)*(ymax - ymin)*(zmax - zmin)*density_boost_converter
            w = np.full(np, npreal/np)
            species.add_uniform_box(np=np, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                                    vthx=ux_th, vthy=uy_th, vthz=uz_th,
                                    vxmean=ux_m, vymean=uy_m, vzmean=uz_m,
                                    lmomentum=1, spacing='random',
                                    lallindomain=warp.true, w=w)


class ParticleListDistribution(picmistandard.PICMI_ParticleListDistribution):
    def loaddistribution(self, species, layout):
        species.addparticles(x=self.x, y=self.y, z=self.z, ux=self.ux, uy=self.uy, uz=self.uz,
                             vx=None, vy=None, vz=None)


class ParticleDistributionPlanarInjector(picmistandard.PICMI_ParticleDistributionPlanarInjector):
    pass


class GriddedLayout(picmistandard.PICMI_GriddedLayout):
    pass


class PseudoRandomLayout(picmistandard.PICMI_PseudoRandomLayout):
    pass


class BinomialSmoother(picmistandard.PICMI_BinomialSmoother):
    pass


class CylindricalGrid(picmistandard.PICMI_CylindricalGrid):
    def init(self, kw):
        w3d.nx = self.nr
        w3d.ny = 0
        w3d.nz = self.nz
        w3d.xmmin = self.rmin
        w3d.xmmax = self.rmax
        w3d.ymmin = 0.
        w3d.ymmax = 1.
        w3d.zmmin = self.zmin
        w3d.zmmax = self.zmax

        bc_dict = {'dirichlet':warp.dirichlet,
                     'neumann':warp.neumann,
                    'periodic':warp.periodic,
                        'open':warp.openbc}
        self.bounds = [bc_dict[self.lower_boundary_conditions[0]], bc_dict[self.upper_boundary_conditions[0]],
                       bc_dict[self.lower_boundary_conditions[-1]], bc_dict[self.upper_boundary_conditions[-1]]]
        w3d.boundxy = self.bounds[1]
        w3d.bound0 = self.bounds[2]
        w3d.boundnz = self.bounds[3]
        top.pboundxy = self.bounds[1]
        top.pbound0 = self.bounds[2]
        top.pboundnz = self.bounds[3]
        if top.pboundxy == warp.openbc: top.pboundxy = warp.absorb
        if top.pbound0 == warp.openbc: top.pbound0 = warp.absorb
        if top.pboundnz == warp.openbc: top.pboundnz = warp.absorb

        if self.moving_window_velocity is not None:
            top.vbeam = top.vbeamfrm = self.moving_window_velocity[-1]
            top.lgridqnt = true
    

class Cartesian2DGrid(picmistandard.PICMI_Cartesian2DGrid):
    def init(self, kw):
        w3d.nx = self.nx
        w3d.ny = 2
        w3d.nz = self.ny
        w3d.xmmin = self.xmin
        w3d.xmmax = self.xmax
        w3d.ymmin = -float(w3d.ny)/2.
        w3d.ymmax = float(w3d.ny)/2.
        w3d.zmmin = self.ymin
        w3d.zmmax = self.ymax

        bc_dict = {'dirichlet':warp.dirichlet,
                        'neumann':warp.neumann,
                        'periodic':warp.periodic,
                        'open':warp.openbc}
        self.bounds = [bc_dict[self.lower_boundary_conditions[0]], bc_dict[self.upper_boundary_conditions[0]],
                       bc_dict[self.lower_boundary_conditions[1]], bc_dict[self.upper_boundary_conditions[1]]]
        w3d.boundxy = self.bounds[1]
        w3d.bound0 = self.bounds[2]
        w3d.boundnz = self.bounds[3]
        top.pboundxy = self.bounds[1]
        top.pbound0 = self.bounds[2]
        top.pboundnz = self.bounds[3]
        if top.pboundxy == warp.openbc: top.pboundxy = warp.absorb
        if top.pbound0 == warp.openbc: top.pbound0 = warp.absorb
        if top.pboundnz == warp.openbc: top.pboundnz = warp.absorb

        if self.moving_window_velocity is not None:
            top.vbeam = top.vbeamfrm = self.moving_window_velocity[-1]
            top.lgridqnt = true


class Cartesian3DGrid(picmistandard.PICMI_Cartesian3DGrid):
    def init(self, kw):
        w3d.nx = self.nx
        w3d.ny = self.ny
        w3d.nz = self.nz
        w3d.xmmin = self.xmin
        w3d.xmmax = self.xmax
        w3d.ymmin = self.ymin
        w3d.ymmax = self.ymax
        w3d.zmmin = self.zmin
        w3d.zmmax = self.zmax

        bc_dict = {'dirichlet':warp.dirichlet,
                        'neumann':warp.neumann,
                        'periodic':warp.periodic,
                        'open':warp.openbc}
        self.bounds = [bc_dict[self.lower_boundary_conditions[0]], bc_dict[self.upper_boundary_conditions[0]],
                       bc_dict[self.lower_boundary_conditions[1]], bc_dict[self.upper_boundary_conditions[1]],
                       bc_dict[self.lower_boundary_conditions[2]], bc_dict[self.upper_boundary_conditions[2]]]
        w3d.boundxy = self.bounds[1]
        w3d.bound0 = self.bounds[4]
        w3d.boundnz = self.bounds[5]
        top.pboundxy = self.bounds[1]
        top.pbound0 = self.bounds[4]
        top.pboundnz = self.bounds[5]
        if top.pboundxy == warp.openbc: top.pboundxy = warp.absorb
        if top.pbound0 == warp.openbc: top.pbound0 = warp.absorb
        if top.pboundnz == warp.openbc: top.pboundnz = warp.absorb

        if self.moving_window_velocity is not None:
            top.vbeam = top.vbeamfrm = self.moving_window_velocity[2]
            top.lgridqnt = true


class ElectromagneticSolver(picmistandard.PICMI_ElectromagneticSolver):
    def init(self, kw):
        self.l_correct_num_Cherenkov = kw.pop('warp_l_correct_num_Cherenkov', EM3D.__flaginputs__['l_correct_num_Cherenkov'])
        self.type_rz_depose = kw.pop('warp_type_rz_depose', EM3D.__flaginputs__['type_rz_depose'])
        self.l_setcowancoefs = kw.pop('warp_l_setcowancoefs', EM3D.__flaginputs__['l_setcowancoefs'])
        self.l_getrho = kw.pop('warp_l_getrho', EM3D.__flaginputs__['l_getrho'])

    def initialize_solver_inputs(self):
        if self.method is not None:
            # Stencil controls the courant condition for each solver in WARP
            # Yee solver: stencil = 0 (cdt=1./sqrt(d) dx, where d is dimensionality) 
            # CKC, PSATD solvers: stencil=1 (cdt=dx)
            # N.B: For PSTD solvers and GPSTD solvers the courant condition depends on 
            # the solver order and is currently not automatically implemented in WARP. 
            # In this case and by default, we currently set stencil=0 for these solvers 
            # as well as other solvers. 
            stencil = {'Yee':0, 'CKC':1, 'PSATD':1, 'PSTD':0, 'GPSTD':0}[self.method]
        else:
            stencil = 0

        if self.method in ['PSATD','GPSTD','PSTD']: 
            spectral = 1
            if self.stencil_order is None: 
                # If stencil order is not defined
                # By default, use infinite order stencil 
                self.stencil_order = [-1, -1, -1]
            if self.method == 'PSATD': 
                ntsub = np.inf 
            elif self.method == 'PSTD':
                ntsub = 1
            elif self.method == 'GPSTD':
                ntsub = 2
        else: 
            ntsub = 1
            spectral = 0
            # If stencil_order not defined, 
            # use stencil_order = [2, 2, 2] by default
            if self.stencil_order is None: 
                self.stencil_order = [2, 2, 2]
                
        if isinstance(self.source_smoother, BinomialSmoother): 
            npass_smooth = self.source_smoother.n_pass
            alpha_smooth = self.source_smoother.alpha
            stride_smooth = self.source_smoother.stride
            if (self.grid.number_of_dimensions == 2): 
                for i in range(len(npass_smooth[0])):
                    npass_smooth[1][i] = 0
        else: 
            npass_smooth = [[ 0 ], [ 0 ], [ 0 ]]
            alpha_smooth = [[ 1.], [ 1.], [ 1.]]
            stride_smooth = [[ 1 ], [ 1 ], [ 1 ]]      
            
        self.solver = EM3DFFT(stencil=stencil, 
                              norderx = self.stencil_order[0], 
                              nordery = self.stencil_order[1], 
                              norderz = self.stencil_order[2], 
                              ntsub = ntsub, 
                              l_2dxz = self.grid.number_of_dimensions == 2, 
                              l_2drz = isinstance(self.grid, CylindricalGrid),
                              l_1dz = self.grid.number_of_dimensions == 1, 
                              spectral = spectral, 
                              npass_smooth = npass_smooth,
                              alpha_smooth = alpha_smooth, 
                              stride_smooth = stride_smooth,
                              l_correct_num_Cherenkov = self.l_correct_num_Cherenkov,
                              type_rz_depose = self.type_rz_depose,
                              l_setcowancoefs = self.l_setcowancoefs,
                              current_cor = spectral,
                              l_getrho = self.l_getrho)


class ElectrostaticSolver(picmistandard.PICMI_ElectrostaticSolver):
    def initialize_solver_inputs(self):
        self.solver = MultiGrid3d()


class GaussianLaser(picmistandard.PICMI_GaussianLaser):
    def initialize_laser_inputs(self, solver, antenna, gamma_boost):
        from .init_tools import add_laser
        dim = '%dd'%solver.grid.number_of_dimensions
        if self.zeta is None: self.zeta = 0.
        if self.beta is None: self.beta = 0.
        if self.phi2 is None: self.phi2 = 0.
        antenna_z0 = antenna.position[2]
        add_laser(solver.solver, dim, self.a0, self.waist, self.duration*warp.clight,
                  self.centroid_position[2], self.focal_position[2],
                  lambda0=self.wavelength, theta_pol=self.polarization_angle, source_z=antenna_z0,
                  zeta=self.zeta, beta=self.beta, phi2=self.phi2, 
                  gamma_boost=gamma_boost, laser_file=None, laser_file_energy=None)


class LaserAntenna(picmistandard.PICMI_LaserAntenna):
    pass


class Simulation(picmistandard.PICMI_Simulation):
    def init(self, kw):
        if self.verbose is not None:
            top.verbosity = self.verbose + 1

        if not hasattr(setup, 'pname'):
            # --- Setup the graphics if needed.
            setup()

        self.lebcancel_pusher = kw.pop('warp_lebcancel_pusher', None)
        self.initialize_solver_after_generate = kw.pop('warp_initialize_solver_after_generate', False)

        self.inputs_initialized = False

    def initialize_sim_inputs(self):
        if self.inputs_initialized:
            return

        self.inputs_initialized = True

        top.boost_gamma = self.gamma_boost

        if self.gamma_boost > 1:
            boost = BoostConverter(top.boost_gamma)
            w3d.zmmin, w3d.zmmax = boost.copropag_length([w3d.zmmin, w3d.zmmax],
                                                         beta_object = self.solver.grid.moving_window_velocity[-1]/c)

        for i in range(len(self.species)):
            if self.species[i].particle_shape is None:
                self.species[i].particle_shape = self.particle_shape
            self.species[i].initialize_species_inputs(self.layouts[i])

        if self.lebcancel_pusher is not None:
            top.pgroup.lebcancel_pusher = lebcancel_pusher

        if self.initialize_solver_after_generate:
            # --- Turn off the field solver during the generate
            top.fstype = -1
        else:
            self.solver.initialize_solver_inputs()
            registersolver(self.solver.solver)
            installparticleloader(self.initialize_beam_fields)

        package('w3d')
        generate()

        if self.initialize_solver_after_generate:
            self.solver.initialize_solver_inputs()
            registersolver(self.solver.solver)
            self.initialize_beam_fields()

        for i in range(len(self.lasers)):
            self.lasers[i].initialize_laser_inputs(self.solver, self.laser_injection_methods[i], self.gamma_boost)

        for diag in self.diagnostics:
            diag.initialize_diag_inputs(sim = self)

    def initialize_beam_fields(self):
        """This generates the self-fields of any species with initialize_self_field True
        It'll be called after all of the species are loaded, and before the field solver starts.
        """
        for i in range(len(self.species)):
            if self.initialize_self_fields[i]:
                dim = '%dd'%self.solver.grid.number_of_dimensions
                initialize_beam_fields(self.solver.solver, dim, self.species[i].wspecies, w3d, top)

    def step(self, nsteps=None):
        self.initialize_sim_inputs()
        if nsteps is None:
            nsteps = self.max_steps
        step(nsteps)

    def write_input_file(self, file_name='inputs'):
        self.initialize_sim_inputs()
        pass


class FieldDiagnostic(picmistandard.PICMI_FieldDiagnostic): 
    def initialize_diag_inputs(self, sim):
        if any(self.lower_bound != self.grid.lower_bound) or any(self.upper_bound != self.grid.upper_bound):
            print('Warning: Warp cannot return a subdomain. Bounds set to grid bounds')
        sub_sampling = (np.array(self.grid.number_of_cells)/np.array(self.number_of_cells)).astype('l')
        if self.grid.number_of_dimensions == 2:
            # A list of length 3 is expected.
            sub_sampling = [sub_sampling[0], 1, sub_sampling[1]]
        diag_field = openpmd_diag.FieldDiagnostic(period = self.period,
                                                  top = top,
                                                  w3d = w3d,
                                                  em = sim.solver.solver,
                                                  comm_world = comm_world,
                                                  fieldtypes = self.data_list,
                                                  iteration_min = self.step_min,
                                                  iteration_max = self.step_max,
                                                  sub_sampling = sub_sampling,
                                                  lparallel_output = self.parallelio,
                                                  write_dir = self.write_dir)
        # Install after step 
        installafterstep(diag_field.write)


class ParticleDiagnostic(picmistandard.PICMI_ParticleDiagnostic): 
    def initialize_diag_inputs(self, sim):
        species_dict = dict()
        # Check if self.species is a Species object or an iterable of Specie
        if np.iterable(self.species): 
            for sp in self.species: 
                if isinstance(sp, Species):
                    species_dict[sp.name] = sp.wspecies
        else: 
            if isinstance(self.species, Species):
                species_dict[self.species.name] = self.species.wspecies
        self.species_dict = species_dict 
        
        # Init Warp diag
        diag_part = openpmd_diag.ParticleDiagnostic(period = self.period,
                                                    top = top,
                                                    w3d = w3d,
                                                    species = species_dict,
                                                    comm_world = comm_world,
                                                    particle_data = self.data_list,
                                                    iteration_min = self.step_min,
                                                    iteration_max = self.step_max,
                                                    lparallel_output = self.parallelio,
                                                    write_dir = self.write_dir)
        # Install after step 
        installafterstep(diag_part.write)


class LabFrameFieldDiagnostic(picmistandard.PICMI_LabFrameFieldDiagnostic):
    def init(self, kw):
        self.diag_period = kw.pop('warp_diag_period', 20)

    def initialize_diag_inputs(self, sim):
        if self.grid.moving_window_velocity is not None:
            v_lab = self.grid.moving_window_velocity[-1]
        else:
            v_lab = 0.
        diag_field = openpmd_diag.BoostedFieldDiagnostic(zmin_lab = self.grid.lower_bound[-1],
                                                         zmax_lab = self.grid.upper_bound[-1],
                                                         v_lab = v_lab,
                                                         dt_snapshots_lab = self.dt_snapshots,
                                                         Ntot_snapshots_lab = self.num_snapshots,
                                                         gamma_boost = sim.gamma_boost,
                                                         period = self.diag_period,
                                                         em = sim.solver.solver,
                                                         top = top,
                                                         w3d = w3d,
                                                         comm_world = comm_world,
                                                         fieldtypes = self.data_list,
                                                         z_subsampling = self.z_subsampling,
                                                         write_dir = self.write_dir,
                                                         boost_dir = 1,
                                                         lparallel_output = self.parallelio,
                                                         t_min_lab = self.time_start,
                                                         xmin_lab = None,
                                                         xmax_lab = None,
                                                         ymin_lab = None,
                                                         ymax_lab = None)
        # Install after step 
        installafterstep(diag_field.write)


class LabFrameParticleDiagnostic(picmistandard.PICMI_LabFrameParticleDiagnostic):
    def init(self, kw):
        self.diag_period = kw.pop('warp_diag_period', 20)

    def initialize_diag_inputs(self, sim):
        species_dict = dict()
        # Check if self.species is a Species object or an iterable of Specie
        if np.iterable(self.species): 
            for sp in self.species: 
                if isinstance(sp, Species):
                    species_dict[sp.name] = sp.wspecies
        else: 
            if isinstance(self.species, Species):
                species_dict[self.species.name] = self.species.wspecies
        self.species_dict = species_dict 

        if self.grid.moving_window_velocity is not None:
            v_lab = self.grid.moving_window_velocity[-1]
        else:
            v_lab = 0.

        diag_part = openpmd_diag.BoostedParticleDiagnostic(zmin_lab = self.grid.lower_bound[-1],
                                                           zmax_lab = self.grid.upper_bound[-1],
                                                           v_lab = v_lab,
                                                           dt_snapshots_lab = self.dt_snapshots,
                                                           Ntot_snapshots_lab = self.num_snapshots,
                                                           gamma_boost = sim.gamma_boost,
                                                           period = self.diag_period,
                                                           em = sim.solver.solver,
                                                           top = top,
                                                           w3d = w3d,
                                                           comm_world = comm_world,
                                                           particle_data = self.data_list,
                                                           select = None,
                                                           write_dir = self.write_dir,
                                                           species = species_dict,
                                                           boost_dir = 1,
                                                           lparallel_output = self.parallelio,
                                                           t_min_lab = self.time_start)

