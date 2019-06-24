"""
This file defines the generic class OpenPMDDiagnostic.

This class is a common class from which both ParticleDiagnostic
and FieldDiagnostic inherit
"""
import os
import h5py
import datetime
import shutil
from dateutil.tz import tzlocal
import numpy as np

# Dictionaries of correspondance for openPMD
from data_dict import unit_dimension_dict

class OpenPMDDiagnostic(object) :
    """
    Generic class that contains methods which are common
    to both FieldDiagnostic and ParticleDiagnostic
    """

    def __init__(self, period, top, w3d, comm_world,
        iteration_min=None, iteration_max=None, lparallel_output=False,
        write_metadata_parallel=False, write_dir=None ) :
        """
        General setup of the diagnostic

        Parameters
        ----------
        period : int
            The period of the diagnostics, in number of timesteps.
            (i.e. the diagnostics are written whenever the number
            of iterations is divisible by `period`)

        top : the object representing the `top` package in Warp
            Contains information on the time.

        w3d : the object representing the `w3d` package in Warp
            Contains the dimensions of the grid.

        comm_world : a communicator object
            Either an mpi4py or a pyMPI object, or None (single-proc)

        iteration_min: int, optional
            iteration at which the diagnostic starts to be active

        iteration_max: int, optional
            iteration at which the diagnostic stops being active

        lparallel_output : boolean, optional
            Switch to set output mode (parallel or gathering)
            If "True" : Parallel output

        write_metadata_parallel : boolean
            If "True" : file metadata are written in parallel

        write_dir : string, optional
            The POSIX path to the directory where the results are
            to be written. If none is provided, this will be the path
            of the current working directory
        """
        # Get the rank of this processor
        if comm_world is not None :
            self.rank = comm_world.rank
        else :
            self.rank = 0

        # Register the arguments
        if (iteration_min is None):
            self.iteration_min = 0
        else:
            self.iteration_min = iteration_min
        if (iteration_max is None):
            self.iteration_max = np.inf
        else:
            self.iteration_max = iteration_max
        self.top = top
        self.w3d = w3d
        self.period = period
        self.comm_world = comm_world
        self.lparallel_output = lparallel_output
        self.write_metadata_parallel = write_metadata_parallel
        if (self.comm_world is None) or (self.comm_world.size==1):
            self.lparallel_output = False
            self.write_metadata_parallel = False

        # Get the directory in which to write the data
        if write_dir is None :
            self.write_dir = os.path.join( os.getcwd(), 'diags' )
        else :
            self.write_dir = os.path.abspath(write_dir)

        # Create a few addiditional directories within self.write_dir
        self.create_dir("")
        # If the directory hdf5 exists, remove it. (Otherwise the code
        # may crash if the directory hdf5 contains preexisting files.)
        self.create_dir("hdf5")

    def open_file( self, fullpath, parallel_open,comm=None ):
        """
        Open a file in parallel on several processors, depending
        on the flag parallel_open and self.rank

        If a processor does not participate in the opening of
        the file, this returns None, for that processor

        Parameter
        ---------
        fullpath: string
            The absolute path to the openPMD file

        parallel_open:
            Whether the file is opened in parallel, or whether
            only processor 0 opens the file

        Returns
        -------
        An h5py.File object, or None
        """
        if(parallel_open):
          if (comm is None): 
            comm = self.comm_world
        # In serial mode, only the first proc opens/creates the file.
        if parallel_open == False and self.rank == 0 :
            # Create the filename and open hdf5 file
            f = h5py.File( fullpath, mode="a" )
        # In parallel mode, all proc open the file
        elif parallel_open == True :
            # Create the filename and open hdf5 file
            f = h5py.File( fullpath, mode="a", driver='mpio',
                           comm=comm)
        else:
            f = None

        return(f)

    def write( self ) :
        """
        Check if the data should be written at this iteration
        (based on self.period) and if yes, write it.

        The variable top.it should be defined in the Python
        environment, and should represent the total number of
        timesteps in the simulation.
        """
        # Check if the fields should be written at this iteration
        if ((self.top.it % self.period == 0) and \
        (self.top.it>=self.iteration_min)    and \
        (self.top.it<=self.iteration_max)):
            # Write the hdf5 file if needed
            self.write_hdf5( self.top.it )

    def create_dir( self, dir_path, remove_existing=False ) :
        """
        Check whether the directory exists, and if not create it.

        Parameter
        ---------
        dir_path : string
           Relative path from the directory where the diagnostics
           are written
        """
        # The following operations are done only by the first processor.
        if self.rank == 0 :

            # Get the full path
            full_path = os.path.join( self.write_dir, dir_path )

            # Check wether it exists, and create it if needed
            if os.path.exists(full_path) and remove_existing==True:
                shutil.rmtree(full_path)
            if os.path.exists(full_path) == False :
                try:
                    os.makedirs(full_path)
                except OSError :
                    pass

    def setup_openpmd_file( self, f, iteration, time, dt ) :
        """
        Sets the attributes of the hdf5 file, that comply with OpenPMD

        Parameter
        ---------
        f : an h5py.File object

        iteration: int
            The iteration number of this diagnostic

        time: float (seconds)
            The physical time at this iteration

        dt: float (seconds)
            The timestep of the simulation
        """
        # Set the attributes of the HDF5 file

        # General attributes
        f.attrs["openPMD"] = np.string_("1.1.0")
        f.attrs["openPMDextension"] = np.uint32(1)
        f.attrs["software"] = np.string_("warp")
        f.attrs["softwareVersion"] = np.string_("4")
        f.attrs["date"] = np.string_(
            datetime.datetime.now(tzlocal()).strftime('%Y-%m-%d %H:%M:%S %z'))
        f.attrs["iterationEncoding"] = np.string_("fileBased")
        f.attrs["iterationFormat"] =  np.string_("data%T.h5")

        # Setup the basePath
        f.attrs["basePath"] = np.string_("/data/%T/")
        base_path = "/data/%d/" %iteration
        bp = f.require_group( base_path )
        bp.attrs["time"] = time
        bp.attrs["dt"] = dt
        bp.attrs["timeUnitSI"] = 1.

    def setup_openpmd_record( self, dset, quantity ) :
        """
        Sets the attributes of a record, that comply with OpenPMD

        Parameter
        ---------
        dset : an h5py.Dataset or h5py.Group object

        quantity : string
           The name of the record considered
        """
        dset.attrs["unitDimension"] = unit_dimension_dict[quantity]
        # timeOffset is set to zero (approximation)
        dset.attrs["timeOffset"] = 0.

    def setup_openpmd_component( self, dset ) :
        """
        Sets the attributes of a component, that comply with OpenPMD

        Parameter
        ---------
        dset : an h5py.Dataset or h5py.Group object
        """
        dset.attrs["unitSI"] = 1.
