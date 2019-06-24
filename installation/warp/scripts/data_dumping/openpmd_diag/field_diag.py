"""
This file defines the class FieldDiagnostic
"""
import os
import numpy as np
from generic_diag import OpenPMDDiagnostic
from field_extraction import get_circ_dataset, \
     get_cart3d_dataset, get_cart2d_dataset, get_cart1d_dataset, \
     get_global_indices

# Import a number of useful dictionaries
from data_dict import field_boundary_dict, particle_boundary_dict, \
     field_solver_dict, x_offset_dict, y_offset_dict, z_offset_dict

class FieldDiagnostic(OpenPMDDiagnostic):
    """
    Class that defines the field diagnostics to be done.

    Usage
    -----
    After initialization, the diagnostic is called by using the
    `write` method.
    """

    def __init__(self, period, em, top, w3d, comm_world=None,
                 iteration_min=None, iteration_max=None,
                 fieldtypes=["rho", "E", "B", "J"],
                 sub_sampling=[1,1,1], write_dir=None,
                 lparallel_output=False, write_metadata_parallel=False ):
        """
        Initialize the field diagnostic.

        Parameters
        ----------
        period: int
            The period of the diagnostics, in number of timesteps.
            (i.e. the diagnostics are written whenever the number
            of iterations is divisible by `period`)

        em: an EM3D object (as defined in em3dsolver)
            Contains the fields data and the different methods to extract it

        top: the object representing the `top` package in Warp
            Contains information on the time.

        w3d: the object representing the `w3d` package in Warp
            Contains the dimensions of the grid.

        comm_world: a communicator object
            Either an mpi4py or a pyMPI object, or None (single-proc)

        iteration_min: int, optional
            iteration at which the diagnostic starts to be active

        iteration_max: int, optional
            iteration at which the diagnostic stops being active

        fieldtypes: a list of strings, optional
            The strings are either "rho", "E", "B" or "J"
            and indicate which field should be written.
            Default: all fields are written

        sub_sampling: a list of integers, optional
            contains subsampling periods of grid quantities along
            each direction x, y, z. Default is 1,1,1 (i.e all points)
            N.B: THIS DOES NOT WORK WITH BOOSTED FRAME

        write_dir: string, optional
            The POSIX path to the directory where the results are
            to be written. If none is provided, this will be the path
            of the current working directory.

        lparallel_output: boolean, optional
            Switch to set output mode (parallel or gathering)
            If "True": Parallel output

        write_metadata_parallel : boolean
            If "True" : file metadata are written in parallel
        """
        # General setup
        OpenPMDDiagnostic.__init__(self, period, top, w3d, comm_world,
            iteration_min=iteration_min, iteration_max=iteration_max,
            lparallel_output=lparallel_output, write_dir=write_dir,
            write_metadata_parallel=write_metadata_parallel )

        # Register the arguments
        self.em = em
        self.fieldtypes = fieldtypes

        # Determine the dimensions (Cartesian or cylindrical)
        if self.em.l_2drz is True:
            self.dim = "circ"
        elif self.em.l_1dz is True:
            self.dim = "1d"
        elif self.em.l_2dxz is True:
            self.dim = "2d"
        else:
            self.dim = "3d"

        # Determine the coordinates
        if self.dim == "circ":
            self.coords = ['r', 't', 'z']
        else:
            self.coords = ['x', 'y', 'z']

        # Compute global indices of the domain
        # With subsampling periods sub_sampling[0,1,2] in x, y, z
        istartx, ixsub, nxsub = get_global_indices(top.fsdecomp.ix,
                                    top.fsdecomp.nx, sub_sampling[0])
        istarty, iysub, nysub = get_global_indices(top.fsdecomp.iy,
                                    top.fsdecomp.ny, sub_sampling[1])
        istartz, izsub, nzsub = get_global_indices(top.fsdecomp.iz,
                                    top.fsdecomp.nz, sub_sampling[2])

        # Determine the global indices of the local domain
        if (self.dim == "1d"):
            global_indices = np.zeros([2,1],dtype=np.int)
        elif (self.dim == "2d") or (self.dim=="circ"):
            global_indices = np.zeros([2,2], dtype = np.int)
        elif self.dim == "3d":
            global_indices = np.zeros([2,3], dtype = np.int)

        # In the z direction
        # Indices within [global_indices[0,-1],global_indices[1,-1][ are dumped
        global_indices[0,-1] = istartz[top.iprocgrid[2]]
        global_indices[1,-1] = global_indices[0,-1] + \
          nzsub[top.iprocgrid[2]]+1

        # In the x direction
        # Indices within [global_indices[0,0],global_indices[1,0][ are dumped
        if self.dim in ["2d","circ","3d"]:
            global_indices[0,0] = istartx[top.iprocgrid[0]]
            global_indices[1,0] = global_indices[0,0] + nxsub[top.iprocgrid[0]] + 1

        # In the y direction
        # Indices within [global_indices[0,1],global_indices[1,1][ are dumped
        if self.dim == "3d":
            global_indices[0,1] = istarty[top.iprocgrid[1]]
            global_indices[1,1] = global_indices[0,1] + \
              nysub[top.iprocgrid[1]] + 1

        # Set FieldDiagnostic attributes for dumping
        self.global_indices = global_indices

        # Dumped data array dimensions in x,y,z
        self.nx = ((em.nx)-(em.nx)%sub_sampling[0])//sub_sampling[0]
        self.ny = ((em.ny)-(em.ny)%sub_sampling[1])//sub_sampling[1]
        self.nz = ((em.nz)-(em.nz)%sub_sampling[2])//sub_sampling[2]

        # Mesh cell sizes of dumped grid arrays
        self.dx = em.dx*sub_sampling[0]
        self.dy = em.dy*sub_sampling[1]
        self.dz = em.dz*sub_sampling[2]

        # Start indices of data to be dumped in non-sampled grid array in x,y,z
        self.start=[istartx[top.iprocgrid[0]]-top.fsdecomp.ix[top.iprocgrid[0]],
                    istarty[top.iprocgrid[1]]-top.fsdecomp.iy[top.iprocgrid[1]],
                    istartz[top.iprocgrid[2]]-top.fsdecomp.iz[top.iprocgrid[2]]]
        # Subsampling period
        self.sub_sampling=sub_sampling


    def write_hdf5( self, iteration ):
        """
        Write an HDF5 file that complies with the OpenPMD standard

        Parameter
        ---------
        iteration: int
             The current iteration number of the simulation.
        """
        # Find the file name
        filename = "data%08d.h5" %iteration
        fullpath = os.path.join( self.write_dir, "hdf5", filename )

        # Create the file and setup its attributes
        zmin = self.top.zgrid + self.w3d.zmmin
        self.create_file_empty_meshes( fullpath, self.top.it,
                self.top.time, self.nz, zmin, self.dz, self.top.dt )

        # Open the file again (possibly in parallel), and get the field path
        f = self.open_file( fullpath, parallel_open=self.lparallel_output )
        # (f is None if this processor does not participate in writing data)
        if f is not None:
            field_path = "/data/%d/fields/" %iteration
            field_grp = f[field_path]
        else:
            field_grp = None

        # Loop over the different quantities that should be written
        for fieldtype in self.fieldtypes:
            # Scalar field
            if fieldtype == "rho":
                self.write_dataset( field_grp, "rho", "rho" )
            # Vector field
            elif fieldtype in ["E", "B", "J"]:
                for coord in self.coords:
                    quantity = "%s%s" %(fieldtype, coord)
                    path = "%s/%s" %(fieldtype, coord)
                    self.write_dataset( field_grp, path, quantity )

        # Close the file
        if f is not None:
            f.close()

    # Writing methods
    # ---------------
    def write_dataset( self, field_grp, path, quantity ):
        """
        Write a given dataset

        Parameters
        ----------
        field_grp: an h5py.Group object
            The group that corresponds to the path indicated in meshesPath

        path: string
            The relative path where to write the dataset, in field_grp

        quantity: string
            Describes which field is being written.
            (Either rho, Er, Et, Ez, Br, Bz, Bt, Jr, Jt or Jz)
        """
        # Extract the correct dataset
        if field_grp is not None:
            dset = field_grp[path]
        else:
            dset = None

        # Circ case
        if self.dim == "circ":
            self.write_circ_dataset( dset, quantity )
        # 1D cartesian
        elif self.dim == "1d":
            self.write_cart1d_dataset( dset, quantity )
        # 2D Cartesian case
        elif self.dim == "2d":
            self.write_cart2d_dataset( dset, quantity )
        # 3D Cartesian case
        elif self.dim == "3d":
            self.write_cart3d_dataset( dset, quantity )

    def write_circ_dataset( self, dset, quantity):
        """
        Write a dataset in Circ coordinates
        """
        # Fill the dataset with these quantities
        # Gathering mode
        if not self.lparallel_output:
            F = get_circ_dataset( self.em, quantity, lgather=True, sub_sampling=self.sub_sampling, start=self.start )
            if self.rank == 0:
                dset[:,:,:] = F
        # Parallel mode
        else:
            F = get_circ_dataset( self.em, quantity, lgather=False, sub_sampling=self.sub_sampling, start=self.start )
            indices = self.global_indices
            with dset.collective:
                dset[ :, indices[0,0]:indices[1,0],
                        indices[0,1]:indices[1,1] ] = F




    def write_cart1d_dataset( self, dset, quantity ):
        """
        Write a 1D dataset in Cartesian coordinates
        """
        # Fill the dataset with these quantities
        # Gathering mode
        if not self.lparallel_output:
            F = get_cart1d_dataset( self.em, quantity, lgather=True,
                    sub_sampling=self.sub_sampling, start=self.start )
            if self.rank == 0:
                dset[:] = F
        # Parallel mode
        else:
            F = get_cart1d_dataset( self.em, quantity, lgather=False,
                    sub_sampling=self.sub_sampling, start=self.start )
            indices = self.global_indices
            with dset.collective:
                dset[ indices[0,0]:indices[1,0]] = F

    def write_cart2d_dataset( self, dset, quantity ):
        """
        Write a dataset in Cartesian coordinates
        """
        # Fill the dataset with these quantities
        # Gathering mode
        if not self.lparallel_output:
            F = get_cart2d_dataset( self.em, quantity, lgather=True,
                    sub_sampling=self.sub_sampling, start=self.start )
            if self.rank == 0:
                dset[:,:] = F
        # Parallel mode
        else:
            F = get_cart2d_dataset( self.em, quantity, lgather=False,
                    sub_sampling=self.sub_sampling, start=self.start )
            indices = self.global_indices
            with dset.collective:
                dset[ indices[0,0]:indices[1,0],
                        indices[0,1]:indices[1,1] ] = F

    def write_cart3d_dataset( self, dset, quantity ):
        """
        Write a dataset in Cartesian coordinates
        """
        # Fill the dataset with these quantities
        # Gathering mode
        if not self.lparallel_output:
            F = get_cart3d_dataset( self.em, quantity, lgather=True,
                    sub_sampling=self.sub_sampling, start=self.start )
            if self.rank == 0:
                dset[:,:,:] = F
        # Parallel mode
        else:
            F = get_cart3d_dataset( self.em, quantity, lgather=False,
                    sub_sampling=self.sub_sampling, start=self.start )
            indices = self.global_indices
            with dset.collective:
                dset[ indices[0,0]:indices[1,0],
                    indices[0,1]:indices[1,1],
                    indices[0,2]:indices[1,2] ] = F

    # OpenPMD setup methods
    # ---------------------

    def create_file_empty_meshes( self, fullpath, iteration,
                                   time, Nz, zmin, dz, dt, Nx = None, xmin = None, Ny= None, ymin = None ):
        """
        Create an openPMD file with empty meshes and setup all its attributes

        Parameters
        ----------
        fullpath: string
            The absolute path to the file to be created

        iteration: int
            The iteration number of this diagnostic

        time: float (seconds)
            The physical time at this iteration

        Nz: int
            The number of gridpoints along z in this diagnostics

        zmin: float (meters)
            The position of the left end of the box

        dz: float (meters)
            The resolution in z of this diagnostic

        dt: float (seconds)
            The timestep of the simulation

        Nx, Ny: int
            The number of gridpoints along x, y in this diagnostics 
            If None then act as if Nx = nx_globa, Ny = ny_global

        xmin: float (meters)
           The position of the lower boundary of the box along x 
           If None, then act as if xmin = w3d.xmmin

        ymin: float (meters)
           The position of the lower boundary of the box along y
           If None, then act as if ymin = w3d.ymmin

  
        """
        # Determine the shape of the datasets that will be written
        # Circ case
        if(Nx is None): 
            nx = self.nx+1
        else: 
            nx = Nx
        # 3d case
        if(self.dim == "3d"):
            if(Ny is None): 
               ny = self.ny + 1
            else:
               ny = Ny 
        if self.dim == "circ":
            data_shape = ( 2*self.em.circ_m+1, nx, Nz+1 )
        # 1D case
        elif self.dim == "1d":
            data_shape = ( Nz+1, )
        # 2D case
        elif self.dim == "2d":
            data_shape = ( nx, Nz+1 )
        # 3D case
        elif self.dim == "3d":
            data_shape = ( nx, ny, Nz+1 )

        # Create the file (can be done by one proc or in parallel)
        f = self.open_file( fullpath,
            parallel_open=self.write_metadata_parallel)

        # Setup the different layers of the openPMD file
        # (f is None if this processor does not participate is writing data)
        if f is not None:

            # Setup the attributes of the top level of the file
            self.setup_openpmd_file( f, iteration, time, dt )

            # Setup the meshes group (contains all the fields)
            f.attrs["meshesPath"] = np.string_("fields/")
            field_path = "/data/%d/fields/" %iteration
            field_grp = f.require_group(field_path)
            self.setup_openpmd_meshes_group(field_grp)

            # Loop over the different quantities that should be written
            # and setup the corresponding datasets
            for fieldtype in self.fieldtypes:

                # Scalar field
                if fieldtype == "rho":
                    # Setup the dataset
                    dset = field_grp.require_dataset(
                        "rho", data_shape, dtype='f8')
                    self.setup_openpmd_mesh_component( dset, "rho" )
                    # Setup the record to which it belongs
                    self.setup_openpmd_mesh_record( dset, "rho", dz, zmin, xmin, ymin )

                # Vector field
                elif fieldtype in ["E", "B", "J"]:
                    # Setup the datasets
                    for coord in self.coords:
                        quantity = "%s%s" %(fieldtype, coord)
                        path = "%s/%s" %(fieldtype, coord)
                        dset = field_grp.require_dataset(
                            path, data_shape, dtype='f8')
                        self.setup_openpmd_mesh_component( dset, quantity )
                    # Setup the record to which they belong
                    self.setup_openpmd_mesh_record(
                        field_grp[fieldtype], fieldtype, dz, zmin, xmin, ymin )

                # Unknown field
                else:
                    raise ValueError(
                        "Invalid string in fieldtypes: %s" %fieldtype)

            # Close the file
            f.close()

    def setup_openpmd_meshes_group( self, dset ):
        """
        Set the attributes that are specific to the mesh path

        Parameter
        ---------
        dset: an h5py.Group object that contains all the mesh quantities
        """
        # Field Solver
        dset.attrs["fieldSolver"] = field_solver_dict[ self.em.stencil ]
        # Field and particle boundary
        # - 1D and Circ
        if self.em.l_1dz:
            dset.attrs["fieldBoundary"] = np.array([
                field_boundary_dict[ self.w3d.bound0 ],
                field_boundary_dict[ self.w3d.boundnz ] ])
            dset.attrs["particleBoundary"] = np.array([
                particle_boundary_dict[ self.top.pbound0 ],
                particle_boundary_dict[ self.top.pboundnz ] ])
        # - 2D and Circ
        elif self.em.l_2dxz:
            dset.attrs["fieldBoundary"] = np.array([
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.bound0 ],
                field_boundary_dict[ self.w3d.boundnz ] ])
            dset.attrs["particleBoundary"] = np.array([
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pbound0 ],
                particle_boundary_dict[ self.top.pboundnz ] ])
        # - 3D
        else:
            dset.attrs["fieldBoundary"] = np.array([
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.bound0 ],
                field_boundary_dict[ self.w3d.boundnz ] ])
            dset.attrs["particleBoundary"] = np.array([
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pbound0 ],
                particle_boundary_dict[ self.top.pboundnz ] ])

        # Current Smoothing
        if np.all( self.em.npass_smooth == 0 ):
            dset.attrs["currentSmoothing"] = np.string_("none")
        else:
            dset.attrs["currentSmoothing"] = np.string_("Binomial")
            dset.attrs["currentSmoothingParameters"] = \
              np.string_(str(self.em.npass_smooth))
        # Charge correction
        dset.attrs["chargeCorrection"] = np.string_("none")


    def setup_openpmd_mesh_record( self, dset, quantity, dz, zmin, xmin = None, ymin = None ):
        """
        Sets the attributes that are specific to a mesh record

        Parameter
        ---------
        dset: an h5py.Dataset or h5py.Group object

        quantity: string
            The name of the record (e.g. "rho", "J", "E" or "B")

        dz: float (meters)
            The resolution in z of this diagnostic

        xmin, ymin, zmin: float (meters)
            The position of the left end of the grid
        """
        # Generic record attributes
        self.setup_openpmd_record( dset, quantity )

        if(xmin is None): 
            xmin_d = self.w3d.xmmin
        else: 
            xmin_d = xmin
        if(ymin is None):
            ymin_d = self.w3d.ymmin
        else:
            ymin_d = ymin 

        # Geometry parameters
        # - thetaMode
        if self.dim == "circ":
            dset.attrs['geometry']  = np.string_("thetaMode")
            dset.attrs['geometryParameters'] = \
              np.string_("m=%d;imag=+" %(self.em.circ_m + 1))
            dset.attrs['gridSpacing'] = np.array([ self.dx, dz ])
            dset.attrs['axisLabels'] = np.array([ b'r', b'z' ])
            dset.attrs["gridGlobalOffset"] = np.array([xmin_d, zmin])
        # - 1D Cartesian
        elif self.dim == "1d":
            dset.attrs['geometry'] = np.string_("cartesian")
            dset.attrs['gridSpacing'] = np.array([ dz ])
            dset.attrs['axisLabels'] = np.array([ b'z' ])
            dset.attrs["gridGlobalOffset"] = np.array([zmin])
        # - 2D Cartesian
        elif self.dim == "2d":
            dset.attrs['geometry'] = np.string_("cartesian")
            dset.attrs['gridSpacing'] = np.array([ self.dx, dz ])
            dset.attrs['axisLabels'] = np.array([ b'x', b'z' ])
            dset.attrs["gridGlobalOffset"] = np.array([xmin_d, zmin])
        # - 3D Cartesian
        elif self.dim == "3d":
            dset.attrs['geometry'] = np.string_("cartesian")
            dset.attrs['gridSpacing'] = np.array([ self.dx, self.dy, dz ])
            dset.attrs['axisLabels'] = np.array([ b'x', b'y', b'z' ])
            dset.attrs["gridGlobalOffset"] = np.array([ xmin_d,
                                                    ymin_d, zmin])
        # Generic attributes
        dset.attrs["dataOrder"] = np.string_("C")
        dset.attrs["gridUnitSI"] = 1.
        dset.attrs["fieldSmoothing"] = np.string_("none")

    def setup_openpmd_mesh_component( self, dset, quantity ):
        """
        Set up the attributes of a mesh component

        Parameter
        ---------
        dset: an h5py.Dataset

        quantity: string
            The field that is being written
        """
        # Generic setup of the component
        self.setup_openpmd_component( dset )

        # Field positions
        if (self.em.l_1dz==True):
            positions = np.array([0.])
        elif (self.em.l_2dxz==True) :
            positions = np.array([0., 0.])
        else:
            positions = np.array([0.,0.,0.])
        # Along z
        positions[-1] = z_offset_dict[ quantity ]
        # Along x (2D xz and 3D xyz cartesian only)
        if (self.em.l_1dz==False) :
            positions[0] = x_offset_dict[ quantity ]
        # Along y (3D xyz Cartesian only)
        if (self.em.l_2dxz==False) :
            positions[1] = y_offset_dict[quantity]

        dset.attrs['position'] = positions
