import numpy as np
from scipy.constants import c, m_e, e
from .boost_tools import BoostConverter
# Import laser antenna and laser profiles
from ..field_solvers.laser.laser_profiles import *
from ..field_solvers.laser.laser_antenna import LaserAntenna
from warp import openbc, w3d

def add_laser( em, dim, a0, w0, ctau, z0, zf=None, lambda0=0.8e-6,
               theta_pol=0., source_z=0., zeta=0, beta=0, phi2=0,
               gamma_boost=None, cep=0.,
               laser_file=None, laser_file_energy=None):
    """
    Add a linearly-polarized, Gaussian laser pulse in the em object,
    by setting the correct laser_func, laser_emax, laser_source_z
    and laser_polangle

    NB: When using this interface, the antenna is necessarily
    motionless in the lab-frame.

    Parameters
    ----------
    em : an EM3D object
       The structure that contains the fields of the simulation

    dim: str
       Either "2d", "3d" or "circ"

    a0 : float (unitless)
       *Used only if no laser_file is provided, i.e. for a Gaussian pulse*
       The a0 of a Gaussian pulse at focus

    w0 : float (in meters)
       *Used only if no laser_file is provided, i.e. for a Gaussian pulse*
       The waist of the Gaussian pulse at focus

    ctau : float (in meters)
       *Used only if no laser_file is provided, i.e. for a Gaussian pulse*
       The "longitudinal waist" (or pulse length) of the Gaussian pulse

    z0 : float (in meters)
       *Used only if no laser_file is provided, i.e. for a Gaussian pulse*
       The position of the laser centroid relative to z=0.

    zf : float (in meters), optional
       *Used only if no laser_file is provided, i.e. for a Gaussian pulse*
       The position of the laser focus relative to z=0.
       If not provided, then the laser focus is at z0

    lambda0 : float (in meters), optional
       The central wavelength of the laser
       Default : 0.8 microns (Ti:Sapph laser)

    theta_pol : float (in radians), optional
       The angle of polarization with respect to the x axis
       Default : 0 rad

    source_z : float (in meters), optional
       The position of the antenna that launches the laser

    zeta: float (in m.s), optional
       *Used only if no laser_file is provided, i.e. for a Gaussian pulse*
       Spatial chirp, at focus,
       as defined in Akturk et al., Opt Express, vol 12, no 19 (2014)

    beta: float (in s), optional
       *Used only if no laser_file is provided, i.e. for a Gaussian pulse*
       Angular dispersion, at focus,
       as defined in Akturk et al., Opt Express, vol 12, no 19 (2014)

    phi2: float (in s^2), optional
       *Used only if no laser_file is provided, i.e. for a Gaussian pulse*
       Temporal chirp, at focus,
       as defined in Akturk et al., Opt Express, vol 12, no 19 (2014)

    gamma_boost : float, optional
        When initializing the laser in a boosted frame, set the value of
        `gamma_boost` to the corresponding Lorentz factor. All the other
        quantities (ctau, zf, source_z, etc.) are to be given in the lab frame.

    cep: float (in rad), optional
        Carrier-Envelope Phase

    laser_file: string, optional
        name of the hdf5 file containing the data. The file objects names should be:
        - x <vector> for 2d and 3d
        - y <vector> for 3d
        - r <vector> for circ
        - t <vector> time vector for 2d, 3d and circ
        - Ereal <2d or 3d matrix> real part of the envelope of the laser field:
            2d matrix (t, x) for 2d
            2d matrix (t, r) for circles
            3d matrix (t, x, y) for 3d
        - Eimag: same as Ereal with the imaginary part of E

        Note that the longitudinal coordinate used here is time, so the front of the
        pulse is in the first elements of Ereal and Eimag along the first dimension.

        To calculate the laser field from this envelope, we use the '+' convention:
        E(t) = ( Ereal + i*Eimag ) * exp( +i*omega0*t )

    laser_file_energy: float (in J), optional
        pulse energy (in Joules). The laser field is rescaled using this factor
    """
    # Wavevector and speed of the antenna
    k0 = 2*np.pi/lambda0
    source_v = 0.
    inv_c = 1./c
    tau = ctau * inv_c
    t_peak = (source_z - z0) * inv_c
    if zf is None:
        focal_length = source_z - z0
    else:
        focal_length = source_z - zf

    # Create a laser_profile object
    # Note that the laser_profile needs to be a callable instance of a class,
    # i.e. an instance of a class with the __call__ method. This avoids the
    # problem of the EM solver not being picklable if laser_func were an
    # instance method, which is not picklable.

    # When running a simulation in boosted frame, convert these parameters
    boost = None
    if (gamma_boost is not None):
        boost = BoostConverter( gamma_boost )
        source_z, = boost.copropag_length([ source_z ],
                                          beta_object=source_v/c)
        source_v, = boost.velocity([ source_v ])

    # - Case of a Gaussian pulse
    if laser_file is None:
        # Create a laser profile object to store these parameters
        if (beta == 0) and (zeta == 0) and (phi2 == 0):
            # Without spatio-temporal correlations
            laser_profile = GaussianProfile( k0, w0, tau, t_peak, a0, dim,
                focal_length=focal_length, boost=boost, source_v=source_v, cep=cep )
        else:
            # With spatio-temporal correlations
            laser_profile = GaussianSTCProfile( k0, w0, tau, t_peak, a0, zeta,
                                   beta, phi2, dim, focal_length=focal_length,
                                   boost=boost, source_v=source_v )

    # - Case of an experimental profile
    else:
        # Create a laser profile object
        laser_profile = ExperimentalProfile( k0, laser_file, laser_file_energy, dim,
                                             boost=boost, source_v=source_v
                                           )

    # Additional parameters for the order of deposition of the laser
    em.laser_depos_order_x=1
    em.laser_depos_order_y=1
    em.laser_depos_order_z=1

    # Add laser thanks to add_extra_antenna
    add_extra_antenna(em, w3d, dim, laser_profile,
                      laser_emax=laser_profile.E0, laser_source_z=source_z,
                      laser_source_v=source_v, laser_polangle=theta_pol)


#===============================================================================
def add_extra_antenna(em, w3d, dim, laser_func, laser_polvector=None,
                      laser_vector=np.array([0.,0.,1.]), laser_spot=None,
                      laser_emax=None, laser_source_z=None,
                      laser_source_v=np.array([0., 0., 0.]),
                      laser_polangle=None):
    """
    Add an extra antenna to the EM3D class.
    The user is free to use one of the two different paradigms to introduce a
    laser :
         - laser_func, laser_vector, laser_polvector, laser_spot, laser_emax
         - laser_func, laser_source_z, laser_polangle, laser_emax
    within or without an antenna velocity.

    In a case where both options are used, the second overwrites the first one.

    Parameters
    ----------
    em : an EM3D object
       The structure that contains the fields of the simulation

    dim: str
       Either "2d", "3d" or "circ"

    laser_func : a laser profile object
        Fonction of (x, y, t) which defines the laser amplitude
        Typically one of the profiles defined in laser_profiles.py

    laser_polvector : 1darray of floats (unitless), of shape 3, optional
        Laser polarization vector

    laser_vector : 1darray of floats (unitless), of shape 3, optional
        Laser directional vector
        Default: [0.,0.,1.]

    laser_spot : 1darray of floats (in meter), of shape 3, optional
        Initial position of the laser centroid

    laser_emax : floats (in V/m), optional
        Maximal expected amplitude generated by the antenna

    laser_source_z : float (in meters), optional
       The position of the antenna that launches the laser

    laser_source_v : 1darray of floats (m/s), of shape 3, optional
       The velocity of the antenna that launches the laser
       Default: [0.,0.,0.]

    laser_polangle : float (in rad), optional
       The angle of polarization with respect to the x axis
    """

    laser_antenna = LaserAntenna(laser_func, laser_vector,
                                 laser_polvector, laser_spot,
                                 laser_emax, laser_source_z, laser_source_v,
                                 laser_polangle, w3d, dim, em.circ_m)
    em.laser_antenna.append( laser_antenna )


#===============================================================================
def retropropagation(em, w3d, negative_propagation=False):
    """
    This routine is used to retropropagate a laser.

        When the function is called, the B field sign is changed to inverse the
        direction of propagation.
        Then, after considering the plane of the antenna, which separates the box
        into 2 half-spaces, all the fields in one half-space are set to zero.
        Only one pulse generated by the antenna is thus kept, depending on the
        value of the flag negative_propagation.
        
    Note that it is meant to work with a single antenna. If em.laser_antenna contains 
    several antennas, only the first one will be retropropagated.

        Parameter:
    -----------

        negative_propagation: boolean
                Indicate the half-space set to 0. If False, it suppresses the pulse
                propagating along laser_vector. If None, none of the spaces are set to
                0.
    """
    f = em.fields

    # Change the sign of B
    f.Bx = - f.Bx
    f.By = - f.By
    f.Bz = - f.Bz

    # Put zero values in the half space where the propagation was positive by
    # default.
    nbpoints = f.Ex.shape
    xmin = w3d.xmminlocal - em.nxguard*em.dx
    xmax = w3d.xmmaxlocal + em.nxguard*em.dx
    ymin = w3d.ymminlocal - em.nyguard*em.dy
    ymax = w3d.ymmaxlocal + em.nyguard*em.dy
    zmin = w3d.zmminlocal - em.nzguard*em.dz
    zmax = w3d.zmmaxlocal + em.nzguard*em.dz

    x = np.linspace(xmin, xmax, nbpoints[0])
    y = np.linspace(ymin, ymax, nbpoints[1])
    z = np.linspace(zmin, zmax, nbpoints[2])
    x,y,z = np.meshgrid(x,y,z,indexing='ij')

    vect = em.laser_antenna[0].vector
    spot = em.laser_antenna[0].spot

    mesh_points_antenna_frame = (x-spot[0]) * vect[0] + (y-spot[1]) * vect[1] \
                                + (z-spot[2]) * vect[2]

    # Set the fields to 0 if negative_propagation is defined
    if negative_propagation is not None :
        # Condition to find the corresponding halfspace depending on the value
        # of negative_propagation.
        if negative_propagation:
            zero_condition = (mesh_points_antenna_frame < 0 )
        else:
            zero_condition = (mesh_points_antenna_frame > 0 )

        # All the field arrays are put to 0 in this halfspace.
        # The Rho array is not reset assuming there were no particles before
        # and then no charges.
        f.Ex[zero_condition] = 0
        f.Ey[zero_condition] = 0
        f.Ez[zero_condition] = 0
        f.Bx[zero_condition] = 0
        f.By[zero_condition] = 0
        f.Bz[zero_condition] = 0
        f.Jx[zero_condition] = 0
        f.Jy[zero_condition] = 0
        f.Jz[zero_condition] = 0

        # Set the fields in the PML to 0 if existing
        b = em.bounds
        list_boundaries = []

        # side
        if b[0] == openbc:
            list_boundaries.append(em.block.sidexl.syf)
        if b[1] == openbc:
            list_boundaries.append(em.block.sidexr.syf)
        if b[2] == openbc:
            list_boundaries.append(em.block.sideyl.syf)
        if b[3] == openbc:
            list_boundaries.append(em.block.sideyr.syf)
        if b[4] == openbc:
            list_boundaries.append(em.block.sidezl.syf)
        if b[5] == openbc:
            list_boundaries.append(em.block.sidezr.syf)

        # edge
        if b[0] == openbc and b[2] == openbc:
            list_boundaries.append(em.block.edgexlyl.syf)
        if b[0] == openbc and b[3] == openbc:
            list_boundaries.append(em.block.edgexlyr.syf)
        if b[0] == openbc and b[4] == openbc:
            list_boundaries.append(em.block.edgexlzl.syf)
        if b[0] == openbc and b[5] == openbc:
            list_boundaries.append(em.block.edgexlzr.syf)
        if b[1] == openbc and b[2] == openbc:
            list_boundaries.append(em.block.edgexryl.syf)
        if b[1] == openbc and b[3] == openbc:
            list_boundaries.append(em.block.edgexryr.syf)
        if b[1] == openbc and b[4] == openbc:
            list_boundaries.append(em.block.edgexrzl.syf)
        if b[1] == openbc and b[5] == openbc:
            list_boundaries.append(em.block.edgexrzr.syf)
        if b[2] == openbc and b[4] == openbc:
            list_boundaries.append(em.block.edgeylzl.syf)
        if b[2] == openbc and b[5] == openbc:
            list_boundaries.append(em.block.edgeylzr.syf)
        if b[3] == openbc and b[4] == openbc:
            list_boundaries.append(em.block.edgeyrzl.syf)
        if b[3] == openbc and b[5] == openbc:
            list_boundaries.append(em.block.edgeyrzr.syf)

        # corner
        if b[0] == openbc and b[2] == openbc and b[4] == openbc:
            list_boundaries.append(em.block.cornerxlylzl.syf)
        if b[0] == openbc and b[2] == openbc and b[5] == openbc:
            list_boundaries.append(em.block.cornerxlylzr.syf)
        if b[0] == openbc and b[3] == openbc and b[4] == openbc:
            list_boundaries.append(em.block.cornerxlyrzl.syf)
        if b[0] == openbc and b[3] == openbc and b[5] == openbc:
            list_boundaries.append(em.block.cornerxlyrzr.syf)
        if b[1] == openbc and b[2] == openbc and b[4] == openbc:
            list_boundaries.append(em.block.cornerxrylzl.syf)
        if b[1] == openbc and b[2] == openbc and b[5] == openbc:
            list_boundaries.append(em.block.cornerxrylzr.syf)
        if b[1] == openbc and b[3] == openbc and b[4] == openbc:
            list_boundaries.append(em.block.cornerxryrzl.syf)
        if b[1] == openbc and b[3] == openbc and b[5] == openbc:
            list_boundaries.append(em.block.cornerxryrzr.syf)

        for syf in list_boundaries :
            syf.exx[...] = 0.;           syf.bxx[...] = 0.
            syf.exy[...] = 0.;           syf.bxy[...] = 0.
            syf.exz[...] = 0.;           syf.bxz[...] = 0.
            syf.eyx[...] = 0.;           syf.byx[...] = 0.
            syf.eyy[...] = 0.;           syf.byy[...] = 0.
            syf.eyz[...] = 0.;           syf.byz[...] = 0.
            syf.ezx[...] = 0.;           syf.bzx[...] = 0.
            syf.ezy[...] = 0.;           syf.bzy[...] = 0.
            syf.ezz[...] = 0.;           syf.bzz[...] = 0.

    print "================================================"
    print " Retropropagation completed."
    print "================================================"
