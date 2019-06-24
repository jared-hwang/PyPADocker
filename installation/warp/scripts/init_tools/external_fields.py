"""
This file defines a set of external fields, i.e. fields that are
calculated analytically and applied to the particles at each timestep.

These fields are not actually evolved on the simulation grid.
"""
from warp import *
import numpy as np
from scipy.constants import m_e, c, e

def add_external_sstf_laser( a0, w0, ctau, zf, tf, beta=0,
                lambda0=0.8e-6, theta_pol=0., forward_propagating=True ):
    """
    Add a linearly-polarized laser with simultaneous spatial
    and temporal focusing (SSTF).

    See Durst et al., Opt Commun. 2008, 281(7):1796-1805 for
    the analytical formula.

    Parameters
    ----------
    a0: float (dimensionless)
        The a0 at focus

    w0: float (meter)
        The waist in the y direction
        This is also the waist in the x direction when beta=0

    ctau: float (meter)
        The length of the pulse when beta=0
        With the notations of Durst et al., ctau = 2*c/\Omega

    zf: float (meter)
        The position of the focal plane

    tf: float (seconds)
        The time at which the on-axis intensity in the focal plane
        is maximal.

    beta: float (dimensionless)
        The dimensionless rate of chirp
        beta = (2*alpha)/(tau*win) where win is the waist before focusing
        and alpha is the rate of spatial chirp (meters.second)

    lambda0: float (dimensionless)
        The wavelength of the laser

    theta_pol: float (radian)
        The angle of polarization of the laser, with respect to the
        x axis.

    forward_propagating: bool (optional)
        Wether the pulse is propagating forward (towards positive z)
        or backward (towards negative z)
    """
    # Create an SstfLaser object to store the parameters
    sstf_laser = SstfLaser( a0, w0, ctau, zf, tf, beta, lambda0,
                theta_pol, forward_propagating=forward_propagating )

    # Install the object in Warp, so that it calls the function
    # add_field_on_particles at each timestep
    installothereuser( sstf_laser.add_field_on_particles )


class SstfLaser( object ):
    """Class that calculates analytically the field of an SSTF laser"""

    def __init__( self, a0, w0, ctau, zf, tf, beta=0,
                lambda0=0.8e-6, theta_pol=0., forward_propagating=True ) :
        """
        Register the parameters

        See the docstring of add_external_sstf_laser for the parameters
        """
        # Intermediate variables
        k0 = 2*np.pi/lambda0

        # Amplitude of the transverse electric and magnetic
        # field, along each axis
        # (These numbers scale the spatial profile which is calculated
        # by self.compute_profile)
        E0 = a0*m_e*c**2*k0/e
        self.E0x = E0 * np.cos(theta_pol)
        self.E0y = E0 * np.sin(theta_pol)
        self.B0x = -self.E0y / c
        self.B0y = self.E0x / c

        # Precalculate and register useful scalars
        self.tf = tf
        self.zf = zf
        self.beta_ctau_over_w0 = beta * ctau / w0
        self.inv_ZR = 2./(k0*w0**2)
        self.i_beta2 = 1.j*beta**2
        self.i_k0 = 1.j*k0
        self.inv_ctau2 = 1./ctau**2
        self.inv_w02 = 1./w0**2

        # Propagation speed
        if forward_propagating is True:
            self.vg = c
        else:
            self.vg = -c


    def add_field_on_particles( self ):
        """
        Function to be called at each timestep, through `installothereuser`

        This function adds the external, analytical SSTF field to
        the fields felt by the macroparticles.
        """
        # Extract the time
        t = top.time

        # Get the indices of the current group of particles
        jmin = w3d.jmin
        jmax = w3d.jmax
        if jmax <= jmin:
            return()

        # Extract the arrays of particle positions
        x = top.pgroup.xp[ jmin:jmax ]
        y = top.pgroup.yp[ jmin:jmax ]
        z = top.pgroup.zp[ jmin:jmax ]

        # Calculate the laser profile (this normalized to 1 at focus)
        field_profile = self.compute_profile( x, y, z, t )

        # Add the fields (with proper dimension and amplitude)
        # to the particle field array
        top.pgroup.ex[ jmin:jmax ] += self.E0x * field_profile
        top.pgroup.ey[ jmin:jmax ] += self.E0y * field_profile
        top.pgroup.bx[ jmin:jmax ] += self.B0x * field_profile
        top.pgroup.by[ jmin:jmax ] += self.B0y * field_profile


    def compute_profile( self, x, y, z, t ):
        """
        Return the normalized laser profile (normalized to 1 at focus)
        at the position of the particles

        Parameters
        ----------
        x, y, z: 1darray of floats (meter)
            Arrays of shape (n_particles,) which contain the particle positions

        t: float (second)
            Time at which the field is calculated

        Returns
        -------
        A 1darray of reals, of shape (n_particles,)
        """
        # Single scalar, that calculates the position of the centroid
        z_centroid = self.zf + self.vg*(t- self.tf)

        # Precalculate different 1d arrays
        # - Normalized position with respect to focus
        z_prime = self.inv_ZR * (z - self.zf)
        # - Diffraction factor
        inv_spatial = 1./(1 + 1.j * z_prime )
        # - Terms for the temporal envelope
        inv_temporal = 1./( 1. + self.i_beta2 * (z_prime*inv_spatial) )
        numerator = z - z_centroid - self.beta_ctau_over_w0 * (x*inv_spatial)
        # - Prefactor
        prefactor = inv_spatial * inv_temporal**0.5

        # Get the complex profile
        profile = prefactor
        # - Introduce propagation phase
        profile *= np.exp( self.i_k0*( z - z_centroid ) )
        # - Introduce envelope
        profile *= np.exp( - self.inv_w02 * ( (x**2+y**2)*inv_spatial ) \
                           - self.inv_ctau2 * ( numerator**2 * inv_temporal ) )

        # Return the real field
        return( profile.real )
