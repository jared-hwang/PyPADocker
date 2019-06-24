"""
Example Pierce diode calculation.
Hot plate source emitting singly ionized potassium
"""
import warp as wp
from warp.run_modes.egun_like import gun

# --- Set four-character run id, comment lines, user's name.
wp.top.pline2 = "Pierce diode example"
wp.top.pline1 = "Injected beam. Semi-Gaus."
wp.top.runmaker = "DPG"

# --- Invoke setup routine for the plotting
wp.setup()

# --- Set the dimensionality
wp.w3d.solvergeom = wp.w3d.RZgeom

# --- Sets method of running
# ---   Steady state gun mode
# ---   Time dependent simulation (when False)
steady_state_gun = True

# --- Basic parameters
channel_radius = 15.*wp.cm

diode_voltage = 93.*wp.kV

# --- Setup source plate
source_radius = 5.5*wp.cm
source_temperature = 0.1  # in eV
source_curvature_radius = 30.*wp.cm  # --- radius of curvature of emitting surface
pierce_angle = 67.

# --- Setup diode aperture plate
zplate = 8.*wp.cm  # --- plate location
rplate = 5.5*wp.cm  # --- aperture radius
plate_width = 2.5*wp.cm  # --- thickness of aperture plate

# --- Setup simulation species
beam = wp.Species(type=wp.Potassium, charge_state=+1, name='beam')

# --- Child-Langmuir current between parallel plates
j = 4./9.*wp.eps0*wp.sqrt(2.*wp.echarge*beam.charge_state/beam.mass)*diode_voltage**1.5/zplate**2
diode_current = wp.pi*source_radius**2*j

print("Child-Langmuir current density = ", j)
print("Child-Langmuir current = ", diode_current)

# --- Set basic beam parameters
beam.a0 = source_radius
beam.b0 = source_radius
beam.ap0 = .0e0
beam.bp0 = .0e0
beam.ibeam = diode_current
beam.vthz = wp.sqrt(source_temperature*wp.jperev/beam.mass)
beam.vthperp = wp.sqrt(source_temperature*wp.jperev/beam.mass)
wp.derivqty()

# --- Length of simulation box
runlen = zplate + 5.*wp.cm

# --- Set boundary conditions
# ---   for field solve
wp.w3d.bound0 = wp.dirichlet
wp.w3d.boundnz = wp.neumann
wp.w3d.boundxy = wp.neumann
# ---   for particles
wp.top.pbound0 = wp.absorb
wp.top.pboundnz = wp.absorb
wp.top.prwall = channel_radius

# --- Set field grid size
wp.w3d.xmmin = -channel_radius
wp.w3d.xmmax = +channel_radius
wp.w3d.ymmin = -channel_radius
wp.w3d.ymmax = +channel_radius
wp.w3d.zmmin = 0.
wp.w3d.zmmax = runlen

# --- Field grid dimensions - note that nx and ny must be even.
wp.w3d.nx = wp.w3d.ny = 32
wp.w3d.nz = 32

# --- Set the time step size. This needs to be small enough to satisfy the Courant limit.
dz = (wp.w3d.zmmax - wp.w3d.zmmin)/wp.w3d.nz
vzfinal = wp.sqrt(2.*diode_voltage*wp.jperev/beam.mass)
wp.top.dt = 0.4*(dz/vzfinal)

# --- Specify injection of the particles
wp.top.inject = 2  # 2 means space-charge limited injection
wp.top.rinject = source_curvature_radius  # Source radius of curvature
wp.top.npinject = 150  # Approximate number of particles injected each step
wp.top.vinject = diode_voltage
wp.w3d.l_inj_exact = True

# --- If using the RZ geometry, set so injection uses the same geometry
wp.w3d.l_inj_rz = (wp.w3d.solvergeom == wp.w3d.RZgeom)

# --- Set up fieldsolver
wp.f3d.mgtol = 1.e-1  # Multigrid solver convergence tolerance, in volts

solver = wp.MultiGrid2D()
wp.registersolver(solver)

piercezlen = (channel_radius - source_radius)*wp.tan((90.-pierce_angle)*wp.pi/180.)
piercezlen = 0.04
rround = plate_width/2.

# --- Create source conductors

# --- Outer radius of Pierce cone
rpierce = source_radius + piercezlen*wp.tan(pierce_angle*wp.pi/180.)

# --- Depth of curved emitting surface
sourcezlen = (source_radius**2/(source_curvature_radius +
                                wp.sqrt(source_curvature_radius**2 - source_radius**2)))

# --- the rsrf and zsrf specify the line in RZ describing the shape of the source and Pierce cone.
# --- The first segment is an arc, the curved emitting surface.
source = wp.ZSrfrv(rsrf=[0., source_radius, rpierce, channel_radius, channel_radius],
                   zsrf=[0., sourcezlen, sourcezlen + piercezlen, sourcezlen + piercezlen, 0.],
                   zc=[source_curvature_radius, None, None, None, None],
                   rc=[0., None, None, None, None],
                   voltage=diode_voltage)

wp.installconductor(source, dfill=wp.largepos)

# --- Create aperture plate
plate = wp.ZRoundedCylinderOut(radius=rplate, length=plate_width, radius2=rround, voltage=0., zcent=zplate)

wp.installconductor(plate, dfill=wp.largepos)

# --- Setup the particle scraper
scraper = wp.ParticleScraper([source, plate])

# --- Set pline1 to include appropriate parameters
if wp.w3d.solvergeom == wp.w3d.RZgeom:
    wp.top.pline1 = ("Injected beam. Semi-Gaus. %dx%d. npinject=%d, dt=%d" %
                     (wp.w3d.nx, wp.w3d.nz, wp.top.npinject, wp.top.dt))
else:
    wp.top.pline1 = ("Injected beam. Semi-Gaus. %dx%dx%d. npinject=%d, dt=%d" %
                     (wp.w3d.nx, wp.w3d.ny, wp.w3d.nz, wp.top.npinject, wp.top.dt))

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.)
wp.package("w3d")
wp.generate()

# --- Open up plotting windows
wp.winon()
wp.winon(1, suffix='current')


def beamplots():
    wp.window(0)
    wp.fma()
    wp.pfzr(plotsg=0, cond=0, titles=False)
    source.draw(filled=150, fullplane=False)
    plate.draw(filled=100, fullplane=False)
    wp.ppzr(titles=False)
    wp.limits(wp.w3d.zmminglobal, wp.w3d.zmmaxglobal, 0., channel_radius)
    wp.ptitles('Hot plate source', 'Z (m)', 'R (m)')
    wp.refresh()

    wp.window(1)
    wp.fma()
    wp.pzcurr()
    wp.limits(wp.w3d.zmminglobal, wp.w3d.zmmaxglobal, 0., diode_current*1.5)
    wp.refresh()


if steady_state_gun:
    # --- Steady-state operation
    # --- This does steady-state gun iterations, plotting the z versus r
    # --- after each iteration.
    wp.top.inj_param = 0.2
    for iter in range(10):
        gun(1, ipstep=1, lvariabletimestep=1)
        beamplots()

else:

    # --- Call beamplots after every 20 steps
    @cwp.allfromafterstep
    def makebeamplots():
        if top.it % 20 == 0:
            beamplots()

    wp.step(700)

# --- Make sure that last plot frames get sent to the cgm file
wp.window(0)
wp.hcp()
wp.window(1)
wp.hcp()
