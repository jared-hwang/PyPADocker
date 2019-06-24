"""
This file defines a set of standard laser profiles, that can be passed
as `laser_func` to the EM3D class or to the `LaserAntenna` class
"""
import numpy as np
from scipy.constants import c, m_e, e, epsilon_0
from scipy.interpolate import RegularGridInterpolator
from scipy.special import genlaguerre, j1
from math import factorial
import h5py
# Try importing parallel functions, in order to broadcast
# the experimental laser file, if required
try:
    from warp_parallel import mpibcast, me
except ImportError:
    # Single-proc simulation
    mpibcast = lambda x:x
    me = 0

class ExperimentalProfile( object ):
    """Class that calculates the laser from a data file."""

    def __init__( self, k0, laser_file, laser_file_energy, dim, boost=None, source_v=0 ):

        """
        k0: float (in rad/m)
            wavenumber of the laser
           
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

        laser_file_energy: float (in J)
            pulse energy (in Joules). The laser field is rescaled using this factor
            
        dim: string '2d', '3d' or 'circ'
            geometry used for the simulation

        boost: a BoostConverter object
            If not None, the laser is emitted in the corresponding boosted-frame
            (even though all parameters are passed in the lab frame)

        source_v: float (in meters/second)
            The speed of the antenna in the direction normal to its plane    
        """

        # The first processor loads the file and sends it to the others
        # (This prevents all the processors from accessing the same file,
        # which can substantially slow down the simulation)
        
        self.v_antenna = source_v
        self.boost = boost
        self.dim = dim
                
        if me==0:
            with h5py.File(laser_file, 'r') as f:
                if self.dim == '2d':
                    t = f['t'][:]
                    x = f['x'][:]
                    Ereal = f['Ereal'][:,:]
                    Eimag = f['Eimag'][:,:]
                elif self.dim == 'circ':
                    t = f['t'][:]
                    r = f['r'][:]
                    Ereal = f['Ereal'][:,:]
                    Eimag = f['Eimag'][:,:]
                elif self.dim == '3d':
                    t = f['t'][:]
                    x = f['x'][:]
                    y = f['y'][:]
                    Ereal = f['Ereal'][:,:,:]
                    Eimag = f['Eimag'][:,:,:]
        else:
            x = None
            y = None
            t = None
            Ereal = None
            Eimag = None
        # Broadcast the data to all procs
        if self.dim == '2d':
            x = mpibcast( x )
        elif self.dim == 'circ':
            r = mpibcast( r )
        elif self.dim == '3d':
            x = mpibcast( x )
            y = mpibcast( y )
        t = mpibcast( t )
        Ereal = mpibcast( Ereal )
        Eimag = mpibcast( Eimag )

        # Recover the complex field
        E_data = Ereal + 1.j*Eimag
        
        # Register the wavevector
        self.k0 = k0

        # Rescale field to energy = laser_file_energy and define interpolation function
        # We do the Slowly-Varying Envelope Approximation. That explains the factor 1/2 
        # We also assume that B = E/c           
        if self.dim == '2d':
            # Calculate the pulse energy in laser_file
            dt = t[1] - t[0]
            dx = x[1] - x[0]
            # In 2d, the energy must be given in J/m
            # This is equivalent to assuming that the pulse is 1m-long in y
            energy_data = epsilon_0*dx*c*dt*np.sum(Ereal**2+Eimag**2) / 2
            # Rescale E
            E_data = E_data * np.sqrt(laser_file_energy/energy_data)
            # Interpolation object
            self.interp_func = RegularGridInterpolator( (t, x), E_data,
                                bounds_error=False, fill_value=0. )
        elif self.dim == 'circ':
            # Calculate the pulse energy in laser_file
            dt = t[1] - t[0]
            dr = r[1] - r[0]
            r_matrix = r.reshape(1, r.shape[0])
            energy_data = np.pi*epsilon_0*dr*c*dt*np.sum(r_matrix*(Ereal**2+Eimag**2))
            # Rescale E
            E_data = E_data * np.sqrt(laser_file_energy/energy_data)
            # Interpolation object
            self.interp_func = RegularGridInterpolator( (t, r), E_data,
                                bounds_error=False, fill_value=0. )
        if self.dim == '3d':
            # Calculate the pulse energy in laser_file
            dt = t[1] - t[0]
            dx = x[1] - x[0]
            dy = y[1] - y[0]
            energy_data = epsilon_0*dx*dy*c*dt*np.sum(Ereal**2+Eimag**2) / 2
            # Rescale E
            E_data = E_data * np.sqrt(laser_file_energy/energy_data)
            # Interpolation object
            self.interp_func = RegularGridInterpolator( (t, x, y), E_data,
                                bounds_error=False, fill_value=0. )

        self.E0 = np.abs( E_data ).max()

    def __call__( self, x, y, t_modified ):
        """
        Return the transverse profile of the laser at the position
        of the antenna

        Parameters:
        -----------
        x: float or ndarray
            First transverse direction in meters

        y: float or ndarray
            Second transverse direction in meters

        t_modified: float
            Time in seconds, multiplied by (1-v_antenna/c)
            This multiplication is done in em3dsolver.py, when
            calling the present function.
        """
        # Get the true time
        # (The top.time has been multiplied by (1-v_antenna/c)
        # in em3dsolver.py, before calling the present function)        
        t = t_modified/(1.-self.v_antenna/c)
        # Get the position of the antenna at this time
        z_source = self.v_antenna * t
        # When running in the boosted frame, convert these position to
        # the lab frame, so as to use the lab-frame formula of the laser
        if self.boost is not None:
            zlab_source = self.boost.gamma0*( z_source + self.boost.beta0*c*t )
            tlab_source = self.boost.gamma0*( t + self.boost.beta0*z_source/c )
            # Overwrite boosted frame values, within the scope of this function
            z_source = zlab_source
            t = tlab_source
            
        # Boosted-frame: convert the laser amplitude
        # These formula assume that the antenna is motionless in the lab frame
        conversion_factor = 1.
        if self.boost is not None:
            conversion_factor = 1./self.boost.gamma0
            # The line below is to compensate the fact that the laser
            # amplitude is multiplied by (1-v_antenna/c) in em3dsolver.py
            conversion_factor *= 1./(1. - self.v_antenna/c)
            
        # Interpolate to find the complex amplitude
        if self.dim == '2d':
            Ecomplex = self.interp_func( (t, x) )
        elif self.dim == 'circ':
            # x is used as the radial coordinate r
            Ecomplex = self.interp_func( (t, x) )
        elif self.dim == '3d':
            Ecomplex = self.interp_func( (t, x, y) )
        # Add laser oscillations and temporal chirp

        Eosc = ( Ecomplex * np.exp( 1.j*self.k0*c*t ) ).real

        return( Eosc * conversion_factor )

class GaussianProfile( object ):
    """Class that calculates a Gaussian laser pulse."""

    def __init__( self, k0, waist, tau, t_peak, a0, dim,
        focal_length=0, temporal_order=2, boost=None, source_v=0, cep=0. ):
        """
        Define a Gaussian laser profile.
        (Gaussian transversally, hypergaussian longitudinally)

        This object can then be passed to the `EM3D` class, as the argument
        `laser_func`, in order to have a Gaussian laser emitted by the antenna.

        Parameters:
        -----------
        k0: float (in meters^-1)
            Laser wavevector (in the lab frame)

        waist: float (in meters)
            Laser waist in the focal plane (in the lab frame)

        tau: float (in seconds)
            Laser temporal waist (in the lab frame)

        t_peak: float (in seconds)
            The time at which the peak of the laser pulse is emitted
            by the antenna (in the lab frame)

        a0: float (dimensionless)
            The peak normalized vector potential, in the focal plane

        dim: string
            The dimension of the simulation. Either "1d", "2d", "circ", or "3d"

        focal_length: float (in meters)
            Distance from the laser antenna to the focal plane
            (along the direction normal to the plane of the antenna)
            Use a positive number for a laser that is focusing as it is
            being emitted by the antenna ; use a negative number for a
            laser that is defocusing as it is being emitted.

        temporal order: int
            The order of the hypergaussian temporal profile
            (Use 2 for a Gaussian temporal profile)

        boost: a BoostConverter object
            If not None, the laser is emitted in the corresponding boosted-frame
            (even though all parameters are passed in the lab frame)

        source_v: float (meters/second)
            The speed of the antenna in the direction normal to its plane
            
        cep: float (rad)
            Carrier-Envelope Phase
        """
        # Set a number of parameters for the laser
        E0 = a0*m_e*c**2*k0/e
        zr = 0.5*k0*waist**2

        # Store the parameters
        self.k0 = k0
        self.inv_waist2 = 1./waist**2
        self.inv_zr = 1./zr
        self.inv_tau = 1./tau
        self.t_peak = t_peak
        self.E0 = E0
        self.v_antenna = source_v
        self.focal_length = focal_length
        self.boost = boost
        self.temporal_order = temporal_order
        self.cep = cep

        # Geometric coefficient (for the evolution of the amplitude)
        self.geom_coeff = get_geometric_coeff( dim )

    def __call__( self, x, y, t_modified ):
        """
        Return the transverse profile of the laser at the position
        of the antenna

        Parameters:
        -----------
        x: float or ndarray
            First transverse direction in meters

        y: float or ndarray
            Second transverse direction in meters

        t_modified: float
            Time in seconds, multiplied by (1-v_antenna/c)
            This multiplication is done in em3dsolver.py, when
            calling the present function.
        """
        # Get the true time
        # (The top.time has been multiplied by (1-v_antenna/c)
        # in em3dsolver.py, before calling the present function)
        t = t_modified/(1.-self.v_antenna/c)
        # Get the position of the antenna at this time
        z_source = self.v_antenna * t
        focal_length = self.focal_length
        # Change focal length in time, when using a moving antenna
        if self.boost is not None:
            focal_length = self.focal_length - self.v_antenna * t

        # When running in the boosted frame, convert these position to
        # the lab frame, so as to use the lab-frame formula of the laser
        if self.boost is not None:
            zlab_source = self.boost.gamma0*( z_source + self.boost.beta0*c*t )
            tlab_source = self.boost.gamma0*( t + self.boost.beta0*z_source/c )
            # Overwrite boosted frame values, within the scope of this function
            z_source = zlab_source
            t = tlab_source

        # Lab-frame formula for the laser:
        # - Waist and curvature and the position of the source
        diffract_factor = 1 - 1j*focal_length*self.inv_zr

        # Calculate the argument of the complex exponential
        exp_argument = 1j * self.k0*c*( t - self.t_peak - z_source/c ) + 1j * self.cep \
          - (x**2 + y**2) * self.inv_waist2 / diffract_factor \
          - ((t - self.t_peak - z_source/c ) * self.inv_tau)**self.temporal_order

        # - Combine profiles
        profile =  np.exp(exp_argument) / ( diffract_factor**self.geom_coeff )

        # Boosted-frame: convert the laser amplitude
        # These formula assume that the antenna is motionless in the lab frame
        if self.boost is not None:
            conversion_factor = 1./self.boost.gamma0
            # The line below is to compensate the fact that the laser
            # amplitude is multiplied by (1-v_antenna/c) in em3dsolver.py
            conversion_factor *= 1./(1. - self.v_antenna/c)
            E0 = conversion_factor * self.E0
        else:
            E0 = self.E0

        # Return the combination of profile and field amplitude
        return( E0 * profile.real )

class JincGaussianAngleProfile( object ):
    """Class that calculates a laser pulse with transverse Jinc profile, longitudinal
       Gaussian profile and propagating at an arbitrary angle with respect to z."""

    def __init__( self, k0, waist, tau, t_peak, a0, dim,
        temporal_order=2, boost=None, source_v=0, theta_zx=0., x_center=0. ):
        """
        Define a laser profile with is a Jinc function transversely and a 
        supergaussian function longitudinally. The antenna is orthogonal to z, but 
        this profile supports a non-zero angle in the (z, x) plane, and is compatible
        with the boosted frame and a moving window.
        
        Note 1: The focal plane is the plane of the antenna. No general analytical formula 
        exist otherwise.
        Note2: This profile only works for small angles theta_zx.

        This object can then be passed to the `EM3D` class, as the argument
        `laser_func`, in order to have a Gaussian laser emitted by the antenna.

        Parameters:
        -----------
        k0: float (in meters^-1)
            Laser wavevector (in the lab frame)

        waist: float (in meters)
            Laser waist in the focal plane (in the lab frame)

        tau: float (in seconds)
            Laser temporal waist (in the lab frame)

        t_peak: float (in seconds)
            The time at which the peak of the laser pulse is emitted
            by the antenna (in the lab frame)

        a0: float (dimensionless)
            The peak normalized vector potential, in the focal plane

        dim: string
            The dimension of the simulation. Either "1d", "2d", "circ", or "3d"

        temporal order: int
            The order of the hypergaussian temporal profile
            (Use 2 for a Gaussian temporal profile)

        boost: a BoostConverter object
            If not None, the laser is emitted in the corresponding boosted-frame
            (even though all parameters are passed in the lab frame)

        source_v: float (meters/second)
            The speed of the antenna in the direction normal to its plane
            
        theta_zx: fload (rad)
            Angle of the k vector in the (z,x) plane with respect to the z axis.
            For instance, 
              theta_zx = 0: along +z
              theta_zx = pi/2: along +x
              
        x_center: float (m)
            Laser centroid position in the transverse direction x.
            Used to inject a laser pulse off-axis.
        """
        # Set a number of parameters for the laser
        E0 = a0*m_e*c**2*k0/e

        # Store the parameters
        self.k0 = k0
        self.waist = waist
        self.inv_tau = 1./tau
        self.t_peak = t_peak
        self.E0 = E0
        self.v_antenna = source_v
        self.boost = boost
        self.temporal_order = temporal_order
        self.theta_zx = theta_zx
        self.x_center = x_center

        # Geometric coefficient (for the evolution of the amplitude)
        self.geom_coeff = get_geometric_coeff( dim )

    def __call__( self, x, y, t_modified ):
        """
        Return the transverse profile of the laser at the position
        of the antenna

        Parameters:
        -----------
        x: float or ndarray
            First transverse direction in meters

        y: float or ndarray
            Second transverse direction in meters

        t_modified: float
            Time in seconds, multiplied by (1-v_antenna/c)
            This multiplication is done in em3dsolver.py, when
            calling the present function.
        """
        # Get the true time
        # (The top.time has been multiplied by (1-v_antenna/c)
        # in em3dsolver.py, before calling the present function)
        t = t_modified/(1.-self.v_antenna/c)
        # Get the position of the antenna at this time
        z_source = self.v_antenna * t

        # When running in the boosted frame, convert these position to
        # the lab frame, so as to use the lab-frame formula of the laser
        if self.boost is not None:
            zlab_source = self.boost.gamma0*( z_source + self.boost.beta0*c*t )
            tlab_source = self.boost.gamma0*( t + self.boost.beta0*z_source/c )
            # Overwrite boosted frame values, within the scope of this function
            z_source = zlab_source
            t = tlab_source
#         Rotated coordinate, to allow for propagation angle
#         x_rotated = (x-self.x_center)*np.cos(self.theta_zx) + c*t*np.sin(self.theta_zx)
#         t_rotated = t*np.cos(self.theta_zx) - (x-self.x_center)/c * np.sin(self.theta_zx)
#         Define spatio-temporal profile
#         r = np.maximum(np.sqrt(x_rotated**2+y**2), self.waist*1.e-8)
#         space_profile = 2*j1(r/self.waist) / (r/self.waist)
#         time_profile  = np.exp(-((t_rotated - self.t_peak - z_source/c ) * self.inv_tau)
#                                       **self.temporal_order)
#         phase         = self.k0*c*(t_rotated-self.t_peak-z_source/c )

        # Rotated coordinate, to allow for propagation angle
        x_rotated = (x-self.x_center)*np.cos(self.theta_zx) + \
                     c*(t-self.t_peak-z_source/c)*np.sin(self.theta_zx)
        t_rotated = (t-self.t_peak-z_source/c)*np.cos(self.theta_zx) -\
                    (x-self.x_center)/c * np.sin(self.theta_zx)
        # Define spatio-temporal profile
        r = np.maximum(np.sqrt(x_rotated**2+y**2), self.waist*1.e-8)
        space_profile = 2*j1(r/self.waist) / (r/self.waist)
        time_profile  = np.exp(-(t_rotated * self.inv_tau)
                                      **self.temporal_order)
        phase         = self.k0*c*t_rotated
        
        # Boosted-frame: convert the laser amplitude
        # These formula assume that the antenna is motionless in the lab frame
        if self.boost is not None:
            conversion_factor = 1./self.boost.gamma0
            # The line below is to compensate the fact that the laser
            # amplitude is multiplied by (1-v_antenna/c) in em3dsolver.py
            conversion_factor *= 1./(1. - self.v_antenna/c)
            E0 = conversion_factor * self.E0            
        else:
            E0 = self.E0
        
        # Return the combination of profile and field amplitude
        return( E0 * space_profile * time_profile * np.cos(phase) )

class GaussianSTCProfile( object ):
    """Class that calculates a Gaussian laser pulse with spatio-temporal
    correlations (i.e. spatial chirp, angular dispersion and temporal chirp)"""

    def __init__( self, k0, waist, tau, t_peak, a0, zeta, beta, phi2, dim,
              focal_length=0, boost=None, source_v=0 ):
        """
        Define a Gaussian laser profile with spatio-temporal correlations
        (i.e. with spatial chirp, angular dispersion and temporal chirp)

        This object can then be passed to the `EM3D` class, as the argument
        `laser_func`, in order to have a Gaussian with spatio-temporal
        correlations laser emitted by the antenna.

        See Akturk et al., Pulse-front tilt caused by spatial and temporal
        chirp, Optics Express, Vol. 12 No. 19 (2014) for more details.

        Parameters:
        -----------
        k0: float (in meters^-1)
            Laser wavevector (in the lab frame)

        waist: float (in meters)
            Laser waist in the focal plane (in the lab frame)

        tau: float (in seconds)
            Laser temporal waist (in the lab frame)

        t_peak: float (in seconds)
            The time at which the peak of the laser pulse is emitted
            by the antenna (in the lab frame)

        a0: float (dimensionless)
            The peak normalized vector potential, in the focal plane

        dim: string
            The dimension of the simulation. Either "1d", "2d", "circ", or "3d"

        zeta: float (in meter.second)
            The amount of spatial chirp, at focus (in the lab frame)
            Namely, a wave packet centered on the frequency w0 + dw
            (where w0 is the central frequency corresponding to k0)
            will reach its peak intensity at an off-axis position
            x(dw) = zeta * dw

        beta: float (in radian.second)
            The amount of angular dispersion, at focus (in the lab frame)
            Namely, a wave packet centered on the frequency (w0 + dw) has
            its wavefronts tilted by an angle theta(dw) = beta*dw
            with respect to the normal of the plane of the antenna.

        phi2: float (in second^2)
            The amount of temporal chirp, at focuse (in the lab frame)
            Namely, a wave packet centered on the frequency (w0 + dw) will
            reach its peak intensity at a time t(dw) = t_peak + phi2*dw.
            Thus, a positive phi2 corresponds to positive chirp, i.e. red part
            of the spectrum in the front of the pulse and blue part of the
            spectrum in the back.

        focal_length: float (in meters)
            Distance from the laser antenna to the focal plane
            (along the direction normal to the plane of the antenna)
            Use a positive number for a laser that is focusing as it is
            being emitted by the antenna ; use a negative number for a
            laser that is defocusing as it is being emitted.

        boost: a BoostConverter object
            If not None, the laser is emitted in the corresponding boosted-frame
            (even though all parameters are passed in the lab frame)

        source_v: float (meters/second)
            The speed of the antenna in the direction normal to its plane
        """

        # Set a number of parameters for the laser
        E0 = a0*m_e*c**2*k0/e
        zr = 0.5*k0*waist**2

        # Store the parameters
        self.k0 = k0
        self.inv_zr = 1./zr
        self.inv_waist2 = 1./waist**2
        self.inv_tau2 = 1/tau**2
        self.focal_length = focal_length
        self.t_peak = t_peak
        self.v_antenna = source_v
        self.E0 = E0
        self.beta = beta
        self.zeta = zeta
        self.phi2 = phi2
        self.boost = boost

        # Geometric coefficient (for the evolution of the amplitude)
        self.geom_coeff = get_geometric_coeff( dim )

    def __call__( self, x, y, t_modified ):
        """
        Return the transverse profile of the laser at the position
        of the antenna

        Parameters:
        -----------
        x: float or ndarray
            First transverse direction in meters

        y: float or ndarray
            Second transverse direction in meters

        t_modified: float
            Time in seconds, multiplied by (1-v_antenna/c)
            This multiplication is done in em3dsolver.py, when
            calling the present function.
        """
        # Get the true time
        # (The top.time has been multiplied by (1-v_antenna/c)
        # in em3dsolver.py, before calling the present function)
        t = t_modified/(1.-self.v_antenna/c)
        # Get the position of the antenna at this time
        z_source = self.v_antenna * t
        focal_length = self.focal_length
        # Change focal length in time, when using a moving antenna
        if self.boost is not None:
            focal_length = self.focal_length - self.v_antenna * t

        # When running in the boosted frame, convert these position to
        # the lab frame, so as to use the lab-frame formula of the laser
        if self.boost is not None:
            zlab_source = self.boost.gamma0*( z_source + self.boost.beta0*c*t )
            tlab_source = self.boost.gamma0*( t + self.boost.beta0*z_source/c )
            # Overwrite boosted frame values, within the scope of this function
            z_source = zlab_source
            t = tlab_source

        # Diffraction and stretching factor
        diffract_factor = 1 - 1j*focal_length*self.inv_zr
        stretch_factor = 1 + \
          4*(self.zeta + self.beta*focal_length)**2 * \
            (self.inv_tau2*self.inv_waist2) / diffract_factor \
        + 2j*(self.phi2 - self.beta**2*self.k0*focal_length) * self.inv_tau2

        # Calculate the argument of the complex exponential
        exp_argument = 1j * self.k0*c*( t - self.t_peak ) \
          - (x**2 + y**2) * self.inv_waist2 / diffract_factor \
          - 1./stretch_factor * self.inv_tau2 * \
            ( t  - self.t_peak - z_source/c - self.beta*self.k0*x \
            - 2j*x*(self.zeta + self.beta*focal_length) \
                *self.inv_waist2/diffract_factor )**2

        # Get the profile
        profile = np.exp(exp_argument) / \
          ( diffract_factor**self.geom_coeff * stretch_factor**0.5 )

        # Boosted-frame: convert the laser amplitude
        # These formula assume that the antenna is motionless in the lab frame
        if self.boost is not None:
            conversion_factor = 1./self.boost.gamma0
            # The line below is to compensate the fact that the laser
            # amplitude is multiplied by (1-v_antenna/c) in em3dsolver.py
            conversion_factor *= 1./(1. - self.v_antenna/c)
            E0 = conversion_factor * self.E0
        else:
            E0 = self.E0

        return( E0 * profile.real )

class LaguerreGaussianProfile(object):
    """
    Class that calculates a Laguerre-Gaussian laser pulse.
    A typical LG pulse is defined as :

    E(x,y,z) = \left(\frac{r \sqrt{2}}{w} \right)^n L_{mn} \left[
                 \frac{2 r^2}{w^2} \right] e^{- i n \varphi} \; GaussianProfile

    where  r = (x^2 + y^2)^(1/2) and \phi = arctan(y/x).
    n and m are specific parameters to calculate the Laguerre
    polynomial :

    L_{mn} \left[  x \right] = \frac{e^x}{m! \: x^n}  \frac{d^m}{dx^m} \left[
                               e^{-x} x^{n+m} \right]

    Be careful, a new Gouy phase is defined such as :
                \psi_{LG}  = (2m + n+ 1) \psi_{G}

    In order to keep a normalized energy, it is necessary to divide the LG
    pulse energy expression by a coefficient alpha (depending on m and n). This
    coefficient can be found analytically and is equal to (m+n)!/m!.
    To take it into account, we divide the electric field by sqrt(alpha).


    Note than when n and m are both equal to 0, this function returns the same
    results as GaussianProfile.
    """

    def __init__( self, m, n, k0, waist, tau, t_peak, a0, dim,
        focal_length=0, temporal_order=2, boost=None, source_v=0 ):
        """
        Define a Laguerre-Gaussian laser profile.
        (Laguerre-Gaussian transversally, hypergaussian longitudinally)

        This object can then be passed to the `EM3D` class, as the argument
        `laser_func`, in order to have a LG laser emitted by the antenna.

        Parameters:
        -----------

        m, n: integer (dimensionless)
            Laguerre polynomial coefficients, cf definition in the docstring
            of this class

        k0: float (in meters^-1)
            Laser wavevector (in the lab frame)

        waist: float (in meters)
            Laser waist in the focal plane (in the lab frame)

        tau: float (in seconds)
            Laser temporal waist (in the lab frame)

        t_peak: float (in seconds)
            The time at which the peak of the laser pulse is emitted
            by the antenna (in the lab frame)

        a0: float (dimensionless)
            The peak normalized vector potential, in the focal plane

        dim: string
            The dimension of the simulation. Either "1d", "2d", "circ", or "3d"

        focal_length: float (in meters)
            Distance from the laser antenna to the focal plane
            (along the direction normal to the plane of the antenna)
            Use a positive number for a laser that is focusing as it is
            being emitted by the antenna ; use a negative number for a
            laser that is defocusing as it is being emitted.

        temporal order: int
            The order of the hypergaussian temporal profile
            (Use 2 for a Gaussian temporal profile)

        boost: a BoostConverter object
            If not None, the laser is emitted in the corresponding boosted-frame
            (even though all parameters are passed in the lab frame)

        source_v: float (meters/second)
            The speed of the antenna in the direction normal to its plane
        """
        # Set a number of parameters for the laser
        E0 = a0*m_e*c**2*k0/e
        zr = 0.5*k0*waist**2

        # Store the parameters
        self.m = m
        self.n = n
        self.k0 = k0
        self.waist = waist
        self.zr = zr
        self.inv_tau = 1./tau
        self.t_peak = t_peak
        self.E0 = E0
        self.v_antenna = source_v
        self.focal_length = focal_length
        self.boost = boost
        self.temporal_order = temporal_order

        # Geometric coefficient (for the evolution of the amplitude)
        self.geom_coeff = get_geometric_coeff( dim )

    def __call__( self, x, y, t_modified ):
        """
        Return the transverse profile of the laser at the position
        of the antenna

        Parameters:
        -----------
        x: float or ndarray
            First transverse direction in meters

        y: float or ndarray
            Second transverse direction in meters

        t_modified: float
            Time in seconds, multiplied by (1-v_antenna/c)
            This multiplication is done in em3dsolver.py, when
            calling the present function.
        """
        # Calculate the array of radius
        r2 = x**2 + y**2
        r = np.sqrt(r2)

        # phi is the argument of x+i*y
        phi = np.angle(x+1j*y)

        # Get the true time
        # (The top.time has been multiplied by (1-v_antenna/c)
        # in em3dsolver.py, before calling the present function)
        t = t_modified/(1.-self.v_antenna/c)
        # Get the position of the antenna at this time
        z_source = self.v_antenna * t
        z = self.focal_length
        # Change focal length in time, when using a moving antenna
        if self.boost is not None:
            z = self.focal_length - self.v_antenna * t

        # When running in the boosted frame, convert these position to
        # the lab frame, so as to use the lab-frame formula of the laser
        if self.boost is not None:
            zlab_source = self.boost.gamma0*( z_source + self.boost.beta0*c*t )
            tlab_source = self.boost.gamma0*( t + self.boost.beta0*z_source/c )
            # Overwrite boosted frame values, within the scope of this function
            z_source = zlab_source
            t = tlab_source

        # Lab-frame formula for the laser:
        # - Waist and curvature and the position of the source
        w = self.waist * np.sqrt( 1 + ( z/self.zr )**2 )
        R = z *( 1 + ( self.zr/z )**2 )

        # Generate the Laguerre function via the scipy function
        L_mn = genlaguerre(self.m, self.n)

        # Calculate alpha, the normalisation coefficient (cf docstring)
        alpha = factorial( self.n + self.m ) / factorial( self.m )

        # - Propagation phase at the position of the source
        propag_phase = self.k0*c*( t - self.t_peak ) \
             + self.k0 * r2 / (2*R) \
             - self.geom_coeff*(self.n + 2*self.m + 1) *np.arctan( z/self.zr )\
             - phi * self.n

        # - Longitudinal and transverse profile
        trans_profile = np.exp( - r2 / w**2 ) \
                        * (r*np.sqrt(2)/w)**self.n * L_mn(2*(r/w)**2)

        long_profile = np.exp(
        - ((t - self.t_peak - z_source/c ) * self.inv_tau)**self.temporal_order)

        # -Curvature oscillations
        curvature_oscillations = np.cos( propag_phase )

        # - Prefactor
        prefactor = (self.waist/w)**self.geom_coeff / np.sqrt(alpha)


        # - Combine profiles
        profile =  prefactor * long_profile * trans_profile \
                    * curvature_oscillations

        # Boosted-frame: convert the laser amplitude
        # These formula assume that the antenna is motionless in the lab frame
        if self.boost is not None:
            conversion_factor = 1./self.boost.gamma0
            # The line below is to compensate the fact that the laser
            # amplitude is multiplied by (1-v_antenna/c) in em3dsolver.py
            conversion_factor *= 1./(1. - self.v_antenna/c)
            E0 = conversion_factor * self.E0
        else:
            E0 = self.E0

        return( E0*profile )


def get_geometric_coeff( dim ):
    """
    Calculate geometric coefficient for the Gouy phase and the evolution
    of the amplitude.

    Parameters:
    -----------
    dim: string
        Either "1d", "2d", "circ" or "3d"
    """
    if  dim=="1d":
        geom_coeff = 0.
    elif dim=="2d":
        geom_coeff = 0.5
    elif dim in ["circ", "3d"]:
        geom_coeff = 1.
    else:
        raise ValueError(
            "`dim` is not recognized: use '1d', '2d', 'circ' or '3d'")
    return( geom_coeff )
