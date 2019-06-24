import numpy as np
#from pywarpx import picmi
from warp import picmi

##########################
# physics parameters
##########################

gamma_boost = 10.

# --- laser

laser_a0              = 1.        # Normalized potential vector
laser_wavelength      = 8e-07     # Wavelength of the laser (in meters)
laser_waist           = 8e-06     # Waist of the laser (in meters)
laser_length          = 3.e-6     # Length of the laser (in meters)
laser_polarization    = 0.        # Polarization angle (in rad)
laser_z0              = -2*laser_length  # Initial position of the centroid (meters)
laser_focal_position  = 0.        # Laser focal position (meters)
laser_antenna_z       = -0.1e-6   # Position of the antenna (meters)

# --- plasma

plasma_density = 2.5e25
plasma_min     = [-25.e-6, -25.e-6,  5.0e-6]
plasma_max     = [ 25.e-6,  25.e-6,  3000.e-6]

# --- Beam
beam_density = 1.e26
beam_uz = 1000.
beam_xmin = -3.e-6
beam_xmax = +3.e-6
beam_ymin = -3.e-6
beam_ymax = +3.e-6
beam_zmin = -12.e-6
beam_zmax = -10.e-6


##########################
# numerics parameters
##########################

# --- Nb time steps

max_steps = 500

# --- grid

nx = 100
ny = 100
nz = 800

xmin = -30.e-6
xmax = +30.e-6
ymin = -30.e-6
ymax = +30.e-6
zmin = -15.e-6
zmax = +5.e-6

moving_window_velocity = [0., 0., picmi.c]

plasma_number_per_cell_each_dim = [2, 4, 2]
beam_number_per_cell_each_dim = [4, 8, 4]

##########################
# physics components
##########################

# --- laser

laser = picmi.GaussianLaser(wavelength            = laser_wavelength,
                            waist                 = laser_waist,
                            duration              = laser_length/picmi.c,
                            focal_position        = [0., 0., laser_focal_position],
                            centroid_position     = [0., 0., laser_z0],
                            polarization_angle    = laser_polarization,
                            propagation_direction = [0,0,1],
                            E0 = laser_a0*2.*np.pi*picmi.m_e*picmi.c**2/(picmi.q_e*laser_wavelength)) # Maximum amplitude of the laser field (in V/m)

laser_antenna = picmi.LaserAntenna(position = [0., 0., laser_antenna_z],  # This point is on the laser plane
                                   normal_vector = [0., 0., 1.])  # The plane normal direction

# --- plasma

uniform_plasma = picmi.UniformDistribution(density     = plasma_density,
                                           lower_bound = plasma_min,
                                           upper_bound = plasma_max,
                                           fill_in     = True)

electrons = picmi.Species(particle_type = 'electron',
                          particle_shape = 'linear',
                          name = 'electrons',
                          initial_distribution = uniform_plasma)

ions = picmi.Species(particle_type = 'H',
                     charge_state = +1,
                     particle_shape = 'linear',
                     name = 'ions',
                     initial_distribution = uniform_plasma)

# --- beam

parabolic_beam = picmi.AnalyticDistribution(density_expression = "beam_density*np.maximum(0., (z - beam_zmin)*(beam_zmax - z)*4/(beam_zmax - beam_zmin)**2*(1. - (sqrt(x**2 + y**2)/beam_rmax)**2))",
                                            beam_density = beam_density,
                                            beam_rmax = beam_xmax,
                                            beam_zmin = beam_zmin,
                                            beam_zmax = beam_zmax,
                                            lower_bound = [beam_xmin, beam_ymin, beam_zmin],
                                            upper_bound = [beam_xmax, beam_ymax, beam_zmax],
                                            directed_velocity = [0., 0., beam_uz])

beam = picmi.Species(particle_type = 'electron',
                     particle_shape = 'linear',
                     name = 'beam',
                     initial_distribution = parabolic_beam)


##########################
# numerics components
##########################

grid = picmi.Cartesian3DGrid(number_of_cells = [nx, ny, nz],
                             lower_bound = [xmin, ymin, zmin],
                             upper_bound = [xmax, ymax, zmax],
                             lower_boundary_conditions = ['open', 'open', 'open'],
                             upper_boundary_conditions = ['open', 'open', 'open'],
                             moving_window_velocity = moving_window_velocity,
                             warpx_max_grid_size=32)

smoother = picmi.BinomialSmoother(n_pass = [[1], [1], [1]],
                                  compensation = [[False], [False], [False]],
                                  stride = [[1], [1], [1]],
                                  alpha = [[0.5], [0.5], [0.5]])

solver = picmi.ElectromagneticSolver(grid = grid,
                                     method = 'Yee',
                                     cfl = 1.,
                                     source_smoother = smoother,
                                     field_smoother = smoother,
                                     warp_l_correct_num_Cherenkov = True,
                                     warp_type_rz_depose = 1,
                                     warp_l_setcowancoefs = True,
                                     warp_l_getrho = True)


##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(grid = grid,
                                   period = 20,
                                   write_dir = 'diags')

part_diag = picmi.ParticleDiagnostic(period = 100,
                                     species = [electrons, ions, beam],
                                     write_dir = 'diags')

field_diag_lab = picmi.LabFrameFieldDiagnostic(grid = grid,
                                               num_snapshots = 20,
                                               dt_snapshots = 0.5*(zmax - zmin)/picmi.c,
                                               data_list = ["rho", "E", "B", "J"],
                                               write_dir = 'lab_diags')

part_diag_lab = picmi.LabFrameParticleDiagnostic(grid = grid,
                                                 num_snapshots = 20,
                                                 dt_snapshots = 0.5*(zmax - zmin)/picmi.c,
                                                 species = [electrons, ions, beam],
                                                 write_dir = 'lab_diags')


##########################
# simulation setup
##########################

sim = picmi.Simulation(solver = solver,
                       max_steps = max_steps,
                       gamma_boost = gamma_boost,
                       verbose = 1,
                       cfl = 1.0,
                       warp_initialize_solver_after_generate = True)

sim.add_species(electrons, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=plasma_number_per_cell_each_dim))
sim.add_species(ions, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=plasma_number_per_cell_each_dim))
sim.add_species(beam, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=beam_number_per_cell_each_dim),
                                                 initialize_self_field = True)

sim.add_laser(laser, injection_method=laser_antenna)

sim.add_diagnostic(field_diag)
#sim.add_diagnostic(part_diag)
sim.add_diagnostic(field_diag_lab)
sim.add_diagnostic(part_diag_lab)

##########################
# simulation run
##########################

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
#sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step(max_steps)
