import os
import shutil
from warp import *
from warp.data_dumping.openpmd_diag import FieldDiagnostic, \
    ParticleDiagnostic, BoostedFieldDiagnostic, \
    BoostedParticleDiagnostic, ProbeParticleDiagnostic

def set_numerics( depos_order, efetch, particle_pusher, dim ):
    """
    Set up some of the parameters of the numerical scheme
    """
    # Generate the grid and simulation structures
    saved_vbeamfrm = top.vbeamfrm
    package('w3d')
    generate()
    top.vbeamfrm = saved_vbeamfrm

    # Particle pusher
    top.lrelativ = true
    top.pgroup.lebcancel_pusher = particle_pusher

    # Gathering and current deposition
    top.depos_order[...] = depos_order
    if dim in ["1d", "2d"]:
        top.depos_order[1,:] = 1
    if dim in  ["1d"]:
        top.depos_order[0,:] = 1
    top.efetch[...] = efetch
    if dim == "circ":
        # Only efetch=1 is implemented in circ
        top.efetch[...] = 1

    # Specific boundary conditions for 1d
    if dim == "1d":
        w3d.boundxy = periodic
        top.pboundxy = periodic

def set_boundary_conditions( f_boundz, f_boundxy, p_boundz, p_boundxy ):
    """
    Set the boundary conditions for the fields and particles
    """
    # Fields
    w3d.bound0 = f_boundz
    w3d.boundnz = f_boundz
    w3d.boundxy = f_boundxy

    # Particle
    top.pbound0 = p_boundz
    top.pboundnz = p_boundz
    top.pboundxy = p_boundxy

def set_simulation_box( Nz, Nx, Ny, zmin, zmax, xmax, ymax, dim,
                        xmin=None, ymin=None):
    """
    Calculate the bounds of simulation box and the cell sizes.
    """
    # Deactivates electrostatic solver
    top.fstype = -1

    # Longitudinal direction
    w3d.zmmin = zmin
    w3d.zmmax = zmax
    w3d.nz = Nz
    w3d.dz = (w3d.zmmax-w3d.zmmin)*1./w3d.nz

    # Transverse direction (correspond to the r direction for "circ")
    w3d.xmmax = xmax
    if dim == "1d":
        w3d.xmmin = -1.
        w3d.xmmax = 1.
        w3d.ymmin = -1.
        w3d.ymmax = 1.
    elif dim in ["2d", "3d"]:
        if xmin is None:
            w3d.xmmin = -xmax
        else:
            w3d.xmmin = xmin
    elif dim == "circ":
        w3d.xmmin = 0
    w3d.nx = Nx
    w3d.dx = (w3d.xmmax-w3d.xmmin)*1./w3d.nx

    # 3rd direction
    if dim=="1d":
        w3d.nx = w3d.ny = 2
    elif dim == "3d":
        w3d.ymmax = ymax
        if ymin is None:
            w3d.ymmin = -ymax
        else:
            w3d.ymmin = ymin
        w3d.ny = Ny
    elif dim in ["2d", "circ"] :
        w3d.ymmax = 1
        w3d.ymmin = -1
        w3d.ny = 2
    w3d.dy = (w3d.ymmax-w3d.ymmin)*1./w3d.ny

    # Specific case of cylindrical geometry
    if dim == "circ":
        w3d.solvergeom=w3d.RZgeom

def set_moving_window( l_moving_window, v_moving_window ):
    """
    Set the attributes of the moving window
    """
    if l_moving_window:
        top.vbeamfrm = v_moving_window

def set_diagnostics( interactive ):
    """
    Set up some general parameters for diagnostics
    """
    # Ensures that no extra saving is done
    top.nhist = top.nt

    # The next switches save computational time by
    # disactivating some diagnostics
    w3d.lrhodia3d = false
    w3d.lgetese3d = false
    w3d.lgtlchg3d = false
    top.ifzmmnt = 0
    top.itmomnts = 0
    top.itplps = 0
    top.itplfreq = 0
    top.iflabwn = 0
    top.zzmomnts = 0
    top.zzplps = 0
    top.zzplfreq = 0
    top.lprntpara = false
    top.lpsplots = false

    # Prepare the gist output in interactive mode
    if interactive==1:
        setup()
        winon(0,dpi=100)


def set_smoothing_parameters(l_smooth, dim, npass_smooth,
                             alpha_smooth, stride_smooth):
    """
    Correct the user-specified smoothing parameters
    """
    if l_smooth:
        # In the 2D case, suppress smoothing in the transverse direction
        if dim=='1d':
            for i in range(len(npass_smooth[0])):
                npass_smooth[0][i]=0
        if dim in ["1d", "2d", "circ"]:
            for i in range( len(npass_smooth[0]) ):
                npass_smooth[1][i]=0
    else:
        # If l_smooth is 0, then set the elements of the list npass_smooth to 0
        for i in range(len(npass_smooth)):
            for j in range(len(npass_smooth[i])):
                npass_smooth[i][j] = 0

def prepare_weights( n_e, nppcellx, nppcelly, nppcellz,
                            dim, circ_m ):
    """
    Return the maximum weight of the macroparticles
    Do a few checks and prepare the weight array of the macroparticles

    Parameters:
    -----------
    n_e : float
        Electron density (in number of physical particles per m^3)

    nppcellx, nppcelly, nppcellz : floats
        Number of macroparticles per cell, along each direction
        When dim="2d", nppcelly is ignored
        When dim="circ", nppcellx corresponds to the radial number of particles
        and nppcelly should be a multiple of 4*circ_m

    dim : str
        Either "2d", "circ" or "3d"

    circ_m : int
        The number of azimuthal modes above m=0

    Returns:
    --------
    A single real number, corresponding to the maximum weight of a species
    """
    # Do a few checks of the input parameters
    if dim=="1d":
        nppcellx = nppcelly =1
    elif dim=="2d":
        nppcelly = 1
    elif dim=="circ" and circ_m > 0:
        if nppcelly % (4*circ_m) != 0:
            raise ValueError('nppcelly should be a multiple of 4*circ_m')

    # Allocate space for variable weights
    if top.wpid == 0:
        top.wpid = nextpid()

    # Determine the maximum volume of a cell
    if dim == "circ":
        Vcell = 2*pi * w3d.xmmax * w3d.dx * w3d.dz
    else:
        Vcell = w3d.dx * w3d.dy * w3d.dz

    # Determine the corresponding maximal weight
    weight = n_e * Vcell / (nppcellx * nppcelly * nppcellz )

    return(weight)

def remove_existing_directory( directory_list ):
    """
    Remove the directories in `directory_list`, if they exist

    Parameters:
    -----------
    directory_list: string or list of string
        Each string is a local path to a directory
    """
    # Only the first proc will delete the directories
    if me==0:
        # If the user passes a string, transform it to a list
        if type(directory_list) is str:
            directory_list = [directory_list]
        # Go through the list and remove directories
        for directory in directory_list:
            if os.path.exists( directory ):
                print 'Removing previous directory %s' %directory
                shutil.rmtree( directory )
