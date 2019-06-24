import numpy as np
from ...warp import eps0 as epsilon_0
from ...warp import PicklableFunction

class LaserAntenna(object):

    def __init__(self, laser_func, vector, polvector, spot, emax,
                 source_z, source_v, polangle, w3d, dim, circ_m ):

        # Initialize the variable self.spot
        if spot is None:
            if source_z is None:
                source_z = w3d.zmmin
            source_z = max(min(source_z,w3d.zmmax),w3d.zmmin)
            self.spot = np.array([ 0, 0, source_z ])
        else:
            self.spot = spot

        # Initialize the variable polvector
        if polvector is None:
            if polangle is None:
                polangle = 0.
            polvector = np.array([ np.cos(polangle), np.sin(polangle), 0.] )

        # Error if the 2 main vectors are not orthognal
        assert np.isclose(np.dot(vector, polvector), 0.), \
            "Error : laser_vector and laser_polvector must be orthogonal. "

        # Normalisation of vector and polvector
        vect_norm = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
        polvect_norm = np.sqrt(polvector[0]**2+polvector[1]**2+polvector[2]**2)
        self.vector = vector/vect_norm
        self.polvector = polvector/polvect_norm

        # Creation of a third vector orthogonal to both vector and polvector
        self.polvector_2 = np.cross(self.vector, self.polvector)

        # Maximal amplitude
        self.emax = emax

        # Function that describes the laser evolution
        self.laser_func = PicklableFunction( laser_func )

        # Dimension
        self.circ_m = circ_m
        self.dim = dim

        # Particle initialization
        self.initialize_virtual_particles( w3d )

        # Antenna velocity
        # If dim(source_v)=1, source_v is consider as a norm and the velocity
        # is normal to the plane of the antenna. If not, source_v is a vector
        # and the antenna moves along this vector.
        if np.size(source_v) == 1:
            self.vx = source_v * self.vector[0]
            self.vy = source_v * self.vector[1]
            self.vz = source_v * self.vector[2]
            self.v  = source_v

        elif np.size(source_v) == 3:
            self.vx = source_v[0]
            self.vy = source_v[1]
            self.vz = source_v[2]
            self.v  = np.sqrt(self.vx**2 + self.vy**2 + self.vz**2)

        else:
            raise TypeError("source_v must be a scalar or a 3-components \
            vector")

    def initialize_virtual_particles( self, w3d ):
        """
        Initialization of the antenna particles depending on the dimension and
        the laser propagation vector.
        """
        # Shortcut definition
        x0 = self.spot[0]
        y0 = self.spot[1]
        z0 = self.spot[2]
        xmin = w3d.xmmin
        xmax = w3d.xmmax
        ymin = w3d.ymmin
        ymax = w3d.ymmax
        zmin = w3d.zmmin
        zmax = w3d.zmmax

        if self.dim == "1d":
            # Ux is chosen orthogonal to self.vector in the plane (x,z)
            Uy = np.array([0.,1.,0.])
            Ux = np.cross(Uy,self.vector)
            self.Ux = Ux
            self.Uy = Uy

            # 1D injection along x
            self.nn_global = 1
            self.xx_global = x0 + np.zeros(self.nn_global)
            self.yy_global = y0 + np.zeros(self.nn_global)
            self.zz_global = z0 + np.zeros(self.nn_global)

        elif self.dim == "circ":
            # 2D circ

            # Check that the normal vector is along z and that the
            # position of the antenna is on the axis
            # (Otherwise the simulation cannot be performed in cylindrical
            # coordinate)
            assert np.allclose( self.vector, np.array([0,0,1]) )
            assert np.allclose( self.spot[:1], np.array([0,0]) )
            # Get the vectors that give the coordinate system of the antenna
            Ux = self.polvector
            Uy = self.polvector_2
            self.Ux = Ux
            self.Uy = Uy

            # The points of the antenna are along a star-pattern
            imin = np.floor( xmin/w3d.dx )
            imax = np.floor( xmax/w3d.dx )
            rr = w3d.dx * np.arange( imin, imax+1 )
            self.weights_circ = 2 * np.pi * rr / w3d.dx
            self.weights_circ /= 4 * self.circ_m
            w0 = self.weights_circ.copy()
            self.xx_global = rr.copy()
            self.yy_global = np.zeros_like( self.xx_global  )
            for i in range( 1, 4*self.circ_m ):
                phase = 0.5*np.pi*float(i)/self.circ_m
                self.xx_global = np.concatenate( (self.xx_global,
                                                    rr*np.cos(phase)) )
                self.yy_global = np.concatenate( (self.yy_global,
                                                    rr*np.sin(phase)) )
                self.weights_circ = np.concatenate((self.weights_circ,w0))
            self.nn_global = np.shape(self.xx_global)[0]
            self.zz_global = z0 + np.zeros(self.nn_global)

        elif self.dim == "2d":
            # 2D plane

            # Ux is chosen orthogonal to self.vector in the plane (x,z)
            Uy = np.array([0.,1.,0.])
            Ux = np.cross(Uy,self.vector)
            self.Ux = Ux
            self.Uy = Uy

            # Spacing between virtual particles to ensure at least one
            # particle per cell
            # select only the Ux components different from 0
            list_Ux = []
            if not Ux[0] == 0.: list_Ux.append(w3d.dx/np.abs(Ux[0]))
            if not Ux[1] == 0.: list_Ux.append(w3d.dy/np.abs(Ux[1]))
            if not Ux[2] == 0.: list_Ux.append(w3d.dz/np.abs(Ux[2]))
            self.Sx = min(list_Ux)

            # Boundaries of the box, depending on sign of the components of Ux
            xmin_i = switch_min_max(xmin, xmax, Ux[0])
            ymin_i = switch_min_max(ymin, ymax, Ux[1])
            zmin_i = switch_min_max(zmin, zmax, Ux[2])
            xmax_i = switch_min_max(xmax, xmin, Ux[0])
            ymax_i = switch_min_max(ymax, ymin, Ux[1])
            zmax_i = switch_min_max(zmax, zmin, Ux[2])

            # Find the range of integer with which the particles will be
            # initialized
            imin = Ux[0]*(xmin_i-x0) + Ux[1]*(ymin_i-y0) + Ux[2]*(zmin_i-z0)
            imax = Ux[0]*(xmax_i-x0) + Ux[1]*(ymax_i-y0) + Ux[2]*(zmax_i-z0)
            imin = np.floor(imin/self.Sx)
            imax = np.floor(imax/self.Sx)+1
            antenna_i = np.arange(imin, imax)

            # Initialize the particle positions
            self.xx_global = x0 + self.Sx*Ux[0]*antenna_i
            self.zz_global = z0 + self.Sx*Ux[2]*antenna_i

            # Keep only the particles that are inside the global box
            is_in_global_box = (self.xx_global >= xmin) \
                                & (self.xx_global < xmax) \
                                & (self.zz_global >= zmin) \
                                & (self.zz_global < zmax)
            self.zz_global = self.zz_global[is_in_global_box]
            self.xx_global = self.xx_global[is_in_global_box]
            self.yy_global = np.zeros(len(self.xx_global))
            # Number of virtual particles
            self.nn_global = np.shape(self.xx_global)[0]

        else:
            # 3D case, Ux = polvector and Uy = polvector_2
            Ux = self.polvector
            Uy = self.polvector_2
            self.Ux = Ux
            self.Uy = Uy

            # Spacing between virtual particles to ensure at least
            # one particle per cell
            # select only the components of Ux and Uy different from 0
            list_Ux = []; list_Uy = []
            if not Ux[0] == 0.: list_Ux.append( w3d.dx/np.abs(Ux[0]) )
            if not Ux[1] == 0.: list_Ux.append( w3d.dy/np.abs(Ux[1]) )
            if not Ux[2] == 0.: list_Ux.append( w3d.dz/np.abs(Ux[2]) )
            if not Uy[0] == 0.: list_Uy.append( w3d.dx/np.abs(Uy[0]) )
            if not Uy[1] == 0.: list_Uy.append( w3d.dy/np.abs(Uy[1]) )
            if not Uy[2] == 0.: list_Uy.append( w3d.dz/np.abs(Uy[2]) )
            self.Sx = min(list_Ux)
            self.Sy = min(list_Uy)

            # Boundaries of the box, depending on sign of the components of Ux
            xmin_i = switch_min_max(xmin, xmax, Ux[0])
            ymin_i = switch_min_max(ymin, ymax, Ux[1])
            zmin_i = switch_min_max(zmin, zmax, Ux[2])
            xmax_i = switch_min_max(xmax, xmin, Ux[0])
            ymax_i = switch_min_max(ymax, ymin, Ux[1])
            zmax_i = switch_min_max(zmax, zmin, Ux[2])
            # Boundaries of the box, depending on sign of the components of Uy
            xmin_j = switch_min_max(xmin, xmax, Uy[0])
            ymin_j = switch_min_max(ymin, ymax, Uy[1])
            zmin_j = switch_min_max(zmin, zmax, Uy[2])
            xmax_j = switch_min_max(xmax, xmin, Uy[0])
            ymax_j = switch_min_max(ymax, ymin, Uy[1])
            zmax_j = switch_min_max(zmax, zmin, Uy[2])

            # Find the range of integer with which the particles will be
            # initialized
            imin = Ux[0]*(xmin_i-x0) + Ux[1]*(ymin_i-y0) + Ux[2]*(zmin_i-z0)
            imax = Ux[0]*(xmax_i-x0) + Ux[1]*(ymax_i-y0) + Ux[2]*(zmax_i-z0)
            jmin = Uy[0]*(xmin_j-x0) + Uy[1]*(ymin_j-y0) + Uy[2]*(zmin_j-z0)
            jmax = Uy[0]*(xmax_j-x0) + Uy[1]*(ymax_j-y0) + Uy[2]*(zmax_j-z0)
            imin = np.floor(imin/self.Sx)
            imax = np.floor(imax/self.Sx)+1
            jmin = np.floor(jmin/self.Sy)
            jmax = np.floor(jmax/self.Sy)+1
            array_i = np.arange(imin, imax)
            array_j = np.arange(jmin, jmax)
            antenna_i, antenna_j = np.meshgrid(array_i,array_j)

            # Initialize the particle positions
            self.xx_global = x0 + self.Sx*Ux[0]*antenna_i + self.Sy*Uy[0]*antenna_j
            self.yy_global = y0 + self.Sx*Ux[1]*antenna_i + self.Sy*Uy[1]*antenna_j
            self.zz_global = z0 + self.Sx*Ux[2]*antenna_i + self.Sy*Uy[2]*antenna_j
            self.xx_global = self.xx_global.flatten()
            self.yy_global = self.yy_global.flatten()
            self.zz_global = self.zz_global.flatten()

            # Keep only the particles that are inside the global box
            is_in_global_box = (self.xx_global >= xmin) \
                                & (self.xx_global < xmax) \
                                & (self.yy_global >= ymin) \
                                & (self.yy_global < ymax) \
                                & (self.zz_global >= zmin) \
                                & (self.zz_global < zmax)
            self.zz_global = self.zz_global[is_in_global_box]
            self.yy_global = self.yy_global[is_in_global_box]
            self.xx_global = self.xx_global[is_in_global_box]
            # Number of virtual particles
            self.nn_global = np.shape(self.xx_global)[0]

        # Set the deplacement around the initial position and normalized momenta
        # variation of each macroparticles to 0
        self.xdx_global = np.zeros(self.nn_global)
        self.ydy_global = np.zeros(self.nn_global)
        self.zdz_global = np.zeros(self.nn_global)
        self.ux_global = np.zeros(self.nn_global)
        self.uy_global = np.zeros(self.nn_global)
        self.uz_global = np.zeros(self.nn_global)
        self.gi_global = np.ones(self.nn_global)

        # Calculate the weights
        self.weights_global = np.ones(self.nn_global) * epsilon_0*self.emax/0.01

        if self.dim == "2d":
            self.weights_global *= self.Sx
        elif self.dim == "3d" :
            self.weights_global *= self.Sy*self.Sx
        elif self.circ_m > 0 : # Circ
            # Laser initialized with particles in a star-pattern
            self.weights_global *= w3d.dx**2 * self.weights_circ

    def push_virtual_particles(self, top, f, clight ):
        """
        Calculate the motion parameters of the laser antenna at a given
        timestep
        """
        dt = top.dt

        self.spot[0] += self.vx * dt
        self.spot[1] += self.vy * dt
        self.spot[2] += self.vz * dt

        x0 = self.spot[0]
        y0 = self.spot[1]
        z0 = self.spot[2]

        Ux = self.Ux
        Uy = self.Uy

        self.xx_global += self.vx * dt
        self.yy_global += self.vy * dt
        self.zz_global += self.vz * dt

        # Coordinate of the antenna in the plane (Ux,Uy)
        x = (self.xx_global-x0)*Ux[0] + (self.yy_global-y0)*Ux[1] + (self.zz_global-z0)*Ux[2]
        y = (self.xx_global-x0)*Uy[0] + (self.yy_global-y0)*Uy[1] + (self.zz_global-z0)*Uy[2]
        t = top.time*(1.-self.v/clight)
        amp = self.laser_func(x,y,t)

        # --- displaces fixed weight particles on "continuous" trajectories
        dispmax = 0.01*clight
        coef_ampli = dispmax * (1.-self.v/clight) / self.emax

        if isinstance(amp,list): #elliptic polarization
            amp_x = amp[0]*self.polvector[0] + amp[1]*self.polvector_2[0]
            amp_y = amp[0]*self.polvector[1] + amp[1]*self.polvector_2[1]
            amp_z = amp[0]*self.polvector[2] + amp[1]*self.polvector_2[2]

            amplitude_x = coef_ampli * amp_x
            amplitude_y = coef_ampli * amp_y
            amplitude_z = coef_ampli * amp_z

        else: #linear polarization
            amplitude_x = coef_ampli * amp * self.polvector[0]
            amplitude_y = coef_ampli * amp * self.polvector[1]
            amplitude_z = coef_ampli * amp * self.polvector[2]

        # Set the amplitude of the normalized momenta of the fictious
        # macroparticles
        self.ux_global[...] = amplitude_x
        self.uy_global[...] = amplitude_y
        self.uz_global[...] = amplitude_z

        # Set the corresponding displacement of the fictious macroparticles
        self.xdx_global[...] += self.ux_global * dt
        self.ydy_global[...] += self.uy_global * dt
        self.zdz_global[...] += self.uz_global * dt

    def select_particles_in_local_box( self, w3d, zgrid ):
        """
        Among the global arrays of particles, select only the particles
        that are currently in the local box.

        These arrays are then used for current deposition.
        """
        # Extract local boundaries
        xmin = w3d.xmminlocal
        xmax = w3d.xmmaxlocal
        ymin = w3d.ymminlocal
        ymax = w3d.ymmaxlocal
        zmin = w3d.zmminlocal + zgrid
        zmax = w3d.zmmaxlocal + zgrid

        # Select the particles that are in the local box
        if self.dim == "1d":
            in_local_box = (self.zz_global >= zmin) & (self.zz_global < zmax)
        elif self.dim == "circ":
            r2 = self.xx_global**2 + self.yy_global**2
            in_local_box = (self.zz_global >= zmin) & (self.zz_global < zmax) \
                         & (r2 >= xmin**2) & (r2 < xmax**2)
        elif self.dim == "2d":
            in_local_box = (self.zz_global >= zmin) & (self.zz_global < zmax) \
                         & (self.xx_global >= xmin) & (self.xx_global < xmax)
        elif self.dim == "3d":
            in_local_box = (self.zz_global >= zmin) & (self.zz_global < zmax) \
                         & (self.xx_global >= xmin) & (self.xx_global < xmax) \
                         & (self.yy_global >= ymin) & (self.yy_global < ymax)

        # Build the arrays of selected particles
        self.nn = in_local_box.sum()
        self.xx = self.xx_global[ in_local_box ]
        self.yy = self.yy_global[ in_local_box ]
        self.zz = self.zz_global[ in_local_box ]
        self.xdx = self.xdx_global[ in_local_box ]
        self.ydy = self.ydy_global[ in_local_box ]
        self.zdz = self.zdz_global[ in_local_box ]
        self.ux = self.ux_global[ in_local_box ]
        self.uy = self.uy_global[ in_local_box ]
        self.uz = self.uz_global[ in_local_box ]
        self.gi = self.gi_global[ in_local_box ]
        self.weights = self.weights_global[ in_local_box ]

# Additional routines
def switch_min_max( x1, x2, u ):
    """
    Return x1 or x2 depending on the sign of u
    """
    if u >= 0 :
        return x1
    else:
        return x2
