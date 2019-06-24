"""
This file defines the class BoostedFieldDiagnostic

Major features:
- The class reuses the existing methods of FieldDiagnostic
  as much as possible, through class inheritance
- The class implements memory buffering of the slices, so as
  not to write to disk at every timestep
"""
import os
import numpy as np
import time
from scipy.constants import c
from field_diag import FieldDiagnostic
from field_extraction import get_dataset
from data_dict import z_offset_dict
from warp_parallel import gather, me, mpiallgather
try:
    from mpi4py import MPI
except ImportError:
    MPI = None
    pass

class BoostedFieldDiagnostic(FieldDiagnostic):
    """
    Class that writes the fields *in the lab frame*, from
    a simulation in the boosted frame

    Usage
    -----
    After initialization, the diagnostic is called by using the
    `write` method.
    """
    def __init__(self, zmin_lab, zmax_lab, v_lab, dt_snapshots_lab,
                 Ntot_snapshots_lab, gamma_boost, period, em, top, w3d,
                 comm_world=None, fieldtypes=["rho", "E", "B", "J"],
                 z_subsampling=1, write_dir=None, boost_dir=1,
                 lparallel_output=False,t_min_lab=0., xmin_lab = None,
                 xmax_lab=None, ymin_lab=None, ymax_lab=None ):
        """
        Initialize diagnostics that retrieve the data in the lab frame,
        as a series of snapshot (one file per snapshot),
        within a virtual moving window defined by zmin_lab, zmax_lab, v_lab.

        Note: In the current implementation, these diagnostics do not
        use parallel HDF5 output. Rank 0 creates and writes all the files.

        Parameters
        ----------
        zmin_lab, zmax_lab: floats (meters)
            Positions of the minimum and maximum of the virtual moving window,
            *in the lab frame*, at t=0


        v_lab: float (m.s^-1)
            Speed of the moving window *in the lab frame*

        dt_snapshots_lab: float (seconds)
            Time interval *in the lab frame* between two successive snapshots

        Ntot_snapshots_lab: int
            Total number of snapshots that this diagnostic will produce
     
        t_min_lab: real (seconds)
            Time for the first snapshot in the lab frame.
            Snapshots are given at t = t_min_lab + i * dt_snapshot_lab -- with i = 0:Ntot_snapshots_lab-1

        period: int
            Number of iterations for which the data is accumulated in memory,
            before finally writing it to the disk.

        z_subsampling: int
            A factor which is applied on the resolution of the lab frame
            reconstruction.

        boost_dir: int (1 or -1)
            The direction of the Lorentz transformation from the lab frame
            to the boosted frame (along the z axis)
                         
        lparallel_output: boolean 
            Enable/disable parallel IO 

        xmin_lab, xmax_lab, ymin_lab, ymax_lab: floats (meters)
            Positions of the minimum and maximum of the virtual moving window,
            *in the lab frame*, at t=0. If None then suppose all the sim box


        See the documentation of FieldDiagnostic for the other parameters
        """
        # Do not leave write_dir as None, as this may conflict with
        # the default directory ('./diags') in which diagnostics in the
        # boosted frame are written.
        if write_dir is None:
            write_dir='lab_diags'

        # Initialize the normal attributes of a FieldDiagnostic
        FieldDiagnostic.__init__(self, period, em, top, w3d,
                comm_world, fieldtypes=fieldtypes, write_dir=write_dir,
                lparallel_output=lparallel_output)


        # Check user input
        boost_dir = int(boost_dir)
        assert boost_dir in [1,-1]

        # Compute global and local indices for transverse directions (x and y)
        self.get_indices_transverse_directions(xmin_lab, xmax_lab, ymin_lab,\
                                               ymax_lab)

        # Register the boost quantities
        self.gamma_boost = gamma_boost
        self.boost_dir = boost_dir
        self.inv_gamma_boost = 1./gamma_boost
        self.beta_boost = np.sqrt( 1. - self.inv_gamma_boost**2 ) * boost_dir
        self.inv_beta_boost = 1./self.beta_boost
        self.lparallel_output = lparallel_output
        self.t_min_lab = t_min_lab
   
        #if parallel output , needs to store the mpi group of comm_world
        if(lparallel_output):  self.mpi_group = self.comm_world.Get_group()


        # Find the z resolution and size of the diagnostic *in the lab frame*
        # (Needed to initialize metadata in the openPMD file)
        dz_lab = np.abs(c*self.top.dt * self.inv_beta_boost*self.inv_gamma_boost)
        Nz = int(round( (zmax_lab - zmin_lab)/dz_lab ))

        # In case of subsampling along z, increase dz and reduce Nz
        if z_subsampling > 1:
            dz_lab = dz_lab * z_subsampling
            Nz = Nz // z_subsampling
        self.inv_dz_lab = 1./dz_lab

        # Create the list of LabSnapshot objects
        self.snapshots = []
        # Record the time it takes
        if self.rank == 0:
            measured_start = time.clock()
            print('\nInitializing the lab-frame diagnostics: %d files...' %(
                Ntot_snapshots_lab) )
        self.Ntot_snapshots_lab = Ntot_snapshots_lab
        # Loop through the lab snapshots and create the corresponding files
        for i in range( Ntot_snapshots_lab ):
            t_lab = i * dt_snapshots_lab + self.t_min_lab
            snapshot = LabSnapshot( t_lab,
                                    zmin_lab + v_lab*t_lab,
                                    zmax_lab + v_lab*t_lab,
                                    self.write_dir, i, self.rank,
                                    boost_dir,lparallel_output)
            self.snapshots.append( snapshot )
            # Initialize a corresponding empty file
            if self.rank == 0:
                self.create_file_empty_meshes( snapshot.filename, i,
                snapshot.t_lab, Nz, snapshot.zmin_lab, dz_lab, self.top.dt, 
                self.Nx_total , xmin_lab, self.Ny_total, ymin_lab   )

        # Print a message that records the time for initialization
        if self.rank == 0:
            measured_end = time.clock()
            print('Time taken for initialization of the files: %.5f s' %(
                measured_end-measured_start) )

        # Create a slice handler, which will do all the extraction, Lorentz
        # transformation, etc for each slice to be registered in a
        # LabSnapshot, and abstracts the dimension

        start = np.array([self.shift_x_min,self.shift_y_min,0])
        self.slice_handler = SliceHandler(self.gamma_boost, self.beta_boost, 
                                          self.dim,start, self.nx_dump,
                                          self.ny_dump )

    def get_indices_transverse_directions(self, xmin_lab, xmax_lab, 
                                          ymin_lab, ymax_lab):
        
        """ This routine compute global and local indices that needed when 
        using xmin_lab ... ymax_lab parameters to dump only a portion of the 
        fields with this diagnostics
        Since the Lorentz transform is done along the z direction, 
        indices along transverse directions remain the same for each dump during
        the simulation """

        # Get local indices of each mpi task
        self.indices = np.copy(self.global_indices)
        self.shift_x_min = 0 
        self.shift_y_min = 0
        self.shift_x_max = 0
        self.shift_y_max = 0

        if(xmin_lab is not None):
            # Compute the shift introduced by xmin_lab for each proc
            self.shift_x_min = int(max(0,(xmin_lab-self.em.xmminlocal)/self.dx))
            self.indices[0,0] += self.shift_x_min
            # Compute the shift introduced by xmin_lab for the whole domain
            self.ix_start_g   = max(0,int((xmin_lab-self.w3d.xmmin)/self.dx))
        else:
            self.ix_start_g = 0
        if self.dim == "3d":
            if(ymin_lab is not None):
                # Compute the shift introduced by ymin_lab for each proc
                self.shift_y_min = int(max(0,(ymin_lab-self.em.ymminlocal)/self.dy))
                self.indices[0,1] += self.shift_y_min   
                # Compute the shift introduced by ymin_lab for the whole domain     
                self.iy_start_g = max(0,int((ymin_lab-self.w3d.ymmin)/self.dy))
            else:
                self.iy_start_g = 0
        
        if(xmax_lab is not None):
            # Compute the shift introduced by xmax_lab for each proc
            self.shift_x_max = int(max(0,(self.em.xmmaxlocal-xmax_lab)/self.dx))
            self.indices[1,0] -= self.shift_x_max
        if self.dim == "3d":
            if(ymax_lab is not None):
                # Compute the shift introduced by ymax_lab for each proc
                shift_y_min = int(max(0,(self.em.ymmaxlocal-ymax_lab)/self.dy))
                self.indices[1,1] -= self.shift_y_max

        # If current mpi subdomain does not intersect with the diag window
        # then set all to 0
        if(xmax_lab is not None and xmax_lab < self.em.xmminlocal): 
            self.indices[0,0] = 0
            self.indices[1,0] = 0
        if(xmin_lab is not None and xmin_lab > self.em.xmmaxlocal):
                self.indices[0,0] = 0
                self.indices[1,0] = 0


        if(self.dim=="3d"):
            if(ymax_lab is not None and ymax_lab < self.em.ymminlocal):
                self.indices[0,0] = 0
                self.indices[1,0] = 0
            if(ymin_lab is not None and  ymin_lab > self.em.ymmaxlocal):
                self.indices[0,0] = 0
                self.indices[1,0] = 0

        # Compute number of data points to dump along x  by current mpi task
        self.nx_dump = max(0,self.indices[1,0] - self.indices[0,0])

        # Compute number of data points to dump along y  by current mpi task
        if self.dim == "3d" :   
            self.ny_dump = max(0,self.indices[1,1] - self.indices[0,1])
        else: 
             self.ny_dump = 0

        self.Nx_total = None
        self.Ny_total = None

        # Compute the total number of points to dump along x and y globally
        if (self.comm_world is not None) and (self.comm_world.size > 1):
            self.global_indices_list = gather( self.indices,
                                               comm=self.comm_world )
            self.Nx_total = gather(self.nx_dump,comm=self.comm_world)
            self.Nx_total = sum(self.Nx_total)//(self.top.fsdecomp.nyprocs*self.top.fsdecomp.nzprocs)

            
            if self.dim == "3d" : 
                self.Ny_total = gather(self.ny_dump, comm=self.comm_world)
                self.Ny_total = sum(self.Ny_total)//(self.top.fsdecomp.nxprocs*self.top.fsdecomp.nzprocs)
        
       
        self.indices[:,0] -= self.ix_start_g
        if(self.dim == "3d"): 
            self.indices[:,1] -= self.iy_start_g

    def write( self ):
        """
        Redefines the method write of the parent class FieldDiagnostic

        Should be registered with installafterstep in Warp
        """
        # At each timestep, store a slices of the fields in memory buffers
        self.store_snapshot_slices()

        # Every self.period, write the buffered slices to disk
        if self.top.it % self.period == 0:
            if not self.lparallel_output: 
              self.flush_to_disk()
            else: 
              self.flush_to_disk_parallel()

    def store_snapshot_slices( self ):
        """
        Store slices of the fields in the memory buffers of the
        corresponding lab snapshots
        """
        # Find the limits of the local subdomain at this iteration
        zmin_boost = self.top.zgrid + self.em.zmminlocal
        zmax_boost = self.top.zgrid + self.em.zmmaxlocal

        # Loop through the labsnapshots
        for snapshot in self.snapshots:

            # Update the positions of the output slice of this snapshot
            # in the lab and boosted frame (current_z_lab and current_z_boost)
            snapshot.update_current_output_positions( self.top.time,
                            self.inv_gamma_boost, self.inv_beta_boost )

            # For this snapshot:
            # - check if the output position *in the boosted frame*
            #   is in the current local domain
            # - check if the output position *in the lab frame*
            #   is within the lab-frame boundaries of the current snapshot
            if ( (snapshot.current_z_boost > zmin_boost) and \
                 (snapshot.current_z_boost < zmax_boost) and \
                 (snapshot.current_z_lab > snapshot.zmin_lab) and \
                 (snapshot.current_z_lab < snapshot.zmax_lab) ):

                # In this case, extract the proper slice from the field array,
                # perform a Lorentz transform to the lab frame, and store
                # the results in a properly-formed array
                slice_array = self.slice_handler.extract_slice(
                    self.em, snapshot.current_z_boost, zmin_boost )
                # Register this in the buffers of this snapshot
                snapshot.register_slice( slice_array, self.inv_dz_lab )


    def flush_to_disk_parallel(self): 
        """ 
        Write the buffered slices of fields to the disk using parallel h5py
        Erease the buffered slices of LabSnapshot objects
        Dada is NOT gathered to proc 0 before being save to disk
        Instead, at each data dump, Ntot_snapshot communicators are created,
        and each hdf5 file is opened calling the corresponding communicator.
        Then each mpi dumps fields on the relevent files.
        """
      
        f2i = self.slice_handler.field_to_index
 
        # allocates Ntot_snapshot arrays of different kinds

        field_array = [None]*self.Ntot_snapshots_lab
          
        iz_min = [None]*self.Ntot_snapshots_lab
        iz_max = [None]*self.Ntot_snapshots_lab
        f = [None]*self.Ntot_snapshots_lab
        
        # this flag is turned true if current mpi is dumping on the i_th snapshot during this flush  
        write_on = [None]*self.Ntot_snapshots_lab
        newgroup = [None]*self.Ntot_snapshots_lab
        dump_comm = [None]*self.Ntot_snapshots_lab
        ranks_group_list = [None]*self.Ntot_snapshots_lab
        field_grp =[None]*self.Ntot_snapshots_lab

        #limits of data dumping along x and y dirctions
        indices = self.indices
        
        # Loop over boosted frame snapshots in order to build an mpi sub comm for each snapshot
        # each communicator encodes informations about which processors need to dump data for this snapshot 
        # Then each mpi open relevent h5files in parallel, all h5files stay open until the flush is completed
        
        
        for i,snapshot in enumerate(self.snapshots):
            #Compact succesive slices that have been buffered 
            #over time into a single array 
            # This returns None, None, None for proc which has no slices
            field_array[i], iz_min[i], iz_max[i] = snapshot.compact_slices()
 
            # if field_array is not None , then this proc will have to dump data
            write_on[i] = False
            #if write_on[i] == True then the current mpi needs to dump data for the i_th snap
            if (field_array[i] is not None): write_on[i] =True
            if(self.nx_dump <= 0) : write_on[i] = False
            if(self.dim == "3d") :  
               if(self.ny_dump <= 0) : write_on[i] = False 
            # Erase the memory buffers
            snapshot.buffered_slices = []
            snapshot.buffer_z_indices = []
            # Creates the subcommunicator that will open the h5 file and dump data
            in_list = -1 
            if ( write_on[i] ):
                in_list = me
            ranks_group_list[i] = mpiallgather( in_list )
            
            ranks_group_list[i] = list(set(ranks_group_list[i]))

            # deletes -1 from the list of ranks
            ranks_group_list[i] = [x for x in ranks_group_list[i] if x >= 0 ]

            newgroup[i] = self.mpi_group.Incl(ranks_group_list[i])
            dump_comm[i] = self.comm_world.Create(newgroup[i])
            #Each mpi opens relevent snapshot files for himself  
            if(dump_comm[i]!= MPI.COMM_NULL):
                f[i] = self.open_file( snapshot.filename, parallel_open=True, comm=dump_comm[i] )
            else: 
                f[i] = None
             
            newgroup[i].Free()
        # Cleans unused data 
        ranks_group_list = []
        newgroup = []
        
        # Dumps data on each snapshot using previously initiized mpi communicators. 
        
        for i,snapshot in enumerate(self.snapshots):  
            if(write_on[i]):
                if f[i] is not None:
                    field_path = "/data/%d/fields/" %snapshot.iteration
                    field_grp[i] = f[i][field_path]
                else:
                    field_grp[i] = None
                # Loop over the different quantities that should be written
                for fieldtype in self.fieldtypes:
                    if fieldtype == "rho":
                        quantity = "rho"
                        path = "rho"
                        if field_grp[i] is not None:
                            dset = field_grp[i][path]
                        else: 
                            dset = None
                        if self.dim == "2d":
                            data = field_array[i][ f2i[ quantity ] ]
                        elif self.dim == "3d":
                            data = field_array[i][ f2i[ quantity ] ]
                        if self.dim == "2d":
                            with dset.collective:
                                dset[ indices[0,0]:indices[1,0],iz_min[i]:iz_max[i] ] = \
                                data[:self.nx_dump,:iz_max[i] -iz_min[i]]

                        elif self.dim == "3d":
                            with dset.collective:
                                dset[indices[0,0]:indices[1,0],indices[0,1]:indices[1,1],iz_min[i]:iz_max[i] ] = \
                                data[:self.nx_dump,:self.ny_dump,:iz_max[i] -iz_min[i]]


                    else:
                        if fieldtype in ["E", "B", "J"]:
                            for coord in self.coords:
                                quantity = "%s%s" %(fieldtype, coord)
                                path = "%s/%s" %(fieldtype, coord)
                                t_i = time.clock()
                                if field_grp is not None:
                                    dset = field_grp[i][path]
                                else: 
                                    dset = None
                                if self.dim == "2d": 
                                    data = field_array[i][ f2i[ quantity ] ]
                                elif self.dim == "3d":
                                    data = field_array[i][ f2i[ quantity ] ]
                                if self.dim == "2d":
                                    with dset.collective:
                                        dset[ indices[0,0]:indices[1,0],iz_min[i]:iz_max[i] ] = \
                                        data[:self.nx_dump,:iz_max[i] -iz_min[i]]
                                elif self.dim == "3d":
                                    with dset.collective:
                                        dset[indices[0,0]:indices[1,0],indices[0,1]:indices[1,1],iz_min[i]:iz_max[i] ] = \
                                        data[:self.nx_dump,:self.ny_dump,:iz_max[i] -iz_min[i]]
        #closes current snapshot file
        for i in range(self.Ntot_snapshots_lab):
            if f[i] is not None:
               f[i].close()
            if dump_comm[i] != MPI.COMM_NULL:
                dump_comm[i].Free()
        #Further cleaning
        data = []
        field_array = []
        f = [] 
        iz_min = []
        iz_max = []
        write_on = []
        dump_comm = []
        field_grp = []
        indices = []
        

        




    def flush_to_disk( self ):
        """
        Writes the buffered slices of fields to the disk

        Erase the buffered slices of the LabSnapshot objects

        Notice: In parallel version, data are gathered to proc 0
        before being saved to disk
        """
        # Loop through the labsnapshots and flush the data
        for snapshot in self.snapshots:
            
            # Compact the successive slices that have been buffered
            # over time into a single array
            # This returns None, None, None for proc which has no slices
            field_array, iz_min, iz_max = snapshot.compact_slices()

            # Erase the memory buffers
            snapshot.buffered_slices = []
            snapshot.buffer_z_indices = []

            # Gather the compacted slices from several proc
            if (self.comm_world is None) or (self.comm_world.size == 1):
                # Serial simulation
                global_field_array = field_array
                global_iz_min = iz_min
                global_iz_max = iz_max
            else:
                # Create new communicator with procs that have non-empty data to send
                in_list = 0
                if (field_array is not None) or (me == 0):
                    in_list = me
                ranks_group_list = mpiallgather( in_list )
                ranks_group_list = list(set(ranks_group_list))
                mpi_group = self.comm_world.Get_group()
                self.ranks_group_list = ranks_group_list
                newgroup = mpi_group.Incl(ranks_group_list)
                dump_comm = self.comm_world.Create(newgroup)
                # Gather data on proc 0 into this communicator
                if dump_comm != MPI.COMM_NULL:
                    field_array_list_comm = dump_comm.gather( field_array )
                    iz_min_list_comm = dump_comm.gather( iz_min )
                    iz_max_list_comm = dump_comm.gather( iz_max )
                
                # First proc: merge the field arrays from each proc
                if self.rank == 0:
                    # Check whether any processor had some slices
                    no_slices = True
                    for i_proc in xrange(dump_comm.Get_size()):
                        if field_array_list_comm[i_proc] is not None:
                            no_slices = False
                    # If there are no slices, set global quantities to None
                    if no_slices:
                        global_field_array = None
                        global_iz_min = None
                        global_iz_max = None
                    # If there are some slices, gather them
                    else:
                        global_field_array, global_iz_min, global_iz_max = \
                          self.gather_slices(field_array_list_comm, 
                              iz_min_list_comm, iz_max_list_comm, dump_comm.Get_size())

                # Free the dump communicator
                mpi_group.Free()
                newgroup.Free()
                if dump_comm != MPI.COMM_NULL:
                    dump_comm.Free()

            # Write the gathered slices to disk
            if (self.rank == 0) and (global_field_array is not None):
                self.write_slices( global_field_array, global_iz_min,
                global_iz_max, snapshot, self.slice_handler.field_to_index )

            # Free gathered arrays
            field_array_list_comm = []
            iz_min_list_comm = []
            iz_max_list_comm = []

    def gather_slices( self, field_array_list, iz_min_list, iz_max_list, size_list ):
        """
        Merge the arrays in field_array_list (one array per proc) into
        a single array

        Parameters
        ----------
        field_array_list: list of arrays
           One element per proc (the element is None if the proc had no data)
           The shape of the array is the one returned by `compact_slices`

        iz_min_list, iz_max_list: lists of integers
           One element per proc (the element is None if the proc had no data)
           The index along z between which each field_array should be included
           into the global array (iz_min is inclusive, iz_max exclusive)
        """
        # Find the global iz_min and global iz_max
        global_iz_min = min([n for n in iz_min_list if n is not None])
        global_iz_max = max([n for n in iz_max_list if n is not None])

        # Allocate the global field array, with the proper size
        nslice = global_iz_max - global_iz_min
        # Circ case
        if self.dim == "circ":
            data_shape = ( 10, 2*self.em.circ_m+1, self.nx+1, nslice )
        # 2D case
        elif self.dim == "2d":
            data_shape = ( 10, self.nx+1, nslice )
        # 3D case
        elif self.dim == "3d":
            data_shape = ( 10, self.nx+1, self.ny+1, nslice )
        global_array = np.zeros( data_shape )
        # Loop through all the processors
        # Fit the field arrays one by one into the global_array
        for i_proc in xrange(size_list):

            i_proc_commworld = self.ranks_group_list[ i_proc ]

            # If this proc has no data, skip it
            if field_array_list[ i_proc ] is None:
                continue

            # Find the indices where the array will be fitted
            ix_min = self.global_indices_list[ i_proc_commworld ][0,0]
            ix_max = self.global_indices_list[ i_proc_commworld ][1,0]
            iy_min = self.global_indices_list[ i_proc_commworld ][0,1]
            iy_max = self.global_indices_list[ i_proc_commworld ][1,1]
            # Longitudinal indices within the array global_array
            s_min = iz_min_list[ i_proc ] - global_iz_min
            s_max = iz_max_list[ i_proc ] - global_iz_min

            # Copy the arrays to the proper position
            if self.dim == "2d":
                global_array[ :, ix_min:ix_max, s_min:s_max ] \
                  = field_array_list[i_proc][ :, :ix_max-ix_min]
            elif self.dim == "3d":
                global_array[ :, ix_min:ix_max, iy_min:iy_max, s_min:s_max ] \
                  = field_array_list[i_proc][ :, :ix_max-ix_min, :iy_max-iy_min]
            elif self.dim == "circ":
                # The second index corresponds to the azimuthal mode
                global_array[ :, :, ix_min:ix_max, s_min:s_max ] \
                  = field_array_list[i_proc][ :, :, :ix_max-ix_min]

        return( global_array, global_iz_min, global_iz_max )


    def write_slices( self, field_array, iz_min, iz_max, snapshot, f2i ):
        """
        For one given snapshot, write the slices of the
        different fields to an openPMD file

        Parameters
        ----------
        field_array: array of reals
            Array of shape
            - (10, em.nxlocal+1, nslices) if dim="2d"
            - (10, em.nxlocal+1, em.nylocal+1, nslices) if dim="3d"
            - (10, 2*em.circ_m+1, em.nxlocal+1, nslices) if dim="circ"

        iz_min, iz_max: integers
            The indices between which the slices will be written
            iz_min is inclusice and iz_max is exclusive

        snapshot: a LabSnaphot object

        f2i: dict
            Dictionary of correspondance between the field names
            and the integer index in the field_array
        """
        # Open the file without parallel I/O in this implementation
        f = self.open_file( snapshot.filename, parallel_open=False )
        field_path = "/data/%d/fields/" %snapshot.iteration
        field_grp = f[field_path]

        # Loop over the different quantities that should be written
        for fieldtype in self.fieldtypes:
            # Scalar field
            if fieldtype == "rho":
                data = field_array[ f2i[ "rho" ] ]
                self.write_field_slices( field_grp, data, "rho",
                                         "rho", iz_min, iz_max )
            # Vector field
            elif fieldtype in ["E", "B", "J"]:
                for coord in self.coords:
                    quantity = "%s%s" %(fieldtype, coord)
                    path = "%s/%s" %(fieldtype, coord)
                    data = field_array[ f2i[ quantity ] ]
                    self.write_field_slices( field_grp, data, path,
                                        quantity, iz_min, iz_max )

        # Close the file
        f.close()

    def write_field_slices( self, field_grp, data, path,
                            quantity, iz_min, iz_max ):
        """
        Writes the slices of a given field into the openPMD file

        Parameters
        ----------
        field_grp: an hdf5.Group
            The h5py group that contains all the meshes

        data: array of reals
            An array containing the slices for one given field

        path: string
            The path of the dataset to write within field_grp

        quantity: string
            A string that indicates which field is being written
            (e.g. 'Ex', 'Br', or 'rho')

        iz_min, iz_max: integers
            The indices between which the slices will be written
            iz_min is inclusice and iz_max is exclusive
        """
        dset = field_grp[ path ]
        indices = self.indices

        # Write the fields depending on the geometry
        if self.dim == "2d":
            dset[ :, iz_min:iz_max ] = data
        elif self.dim == "3d":
            dset[ :, :, iz_min:iz_max ] = data
        elif self.dim == "circ":
            # The first index corresponds to the azimuthal mode
            dset[ :, :, iz_min:iz_max ] = data

class LabSnapshot:
    """
    Class that stores data relative to one given snapshot
    in the lab frame (i.e. one given *time* in the lab frame)
    """

    def __init__(self, t_lab, zmin_lab, zmax_lab,
                 write_dir, i, rank, boost_dir,lparallel=False):
        """
        Initialize a LabSnapshot

        Parameters
        ----------
        t_lab: float (seconds)
            Time of this snapshot *in the lab frame*

        zmin_lab, zmax_lab: floats
            Longitudinal limits of this snapshot

        write_dir: string
            Absolute path to the directory where the data for
            this snapshot is to be written

        i: int
            Number of the file where this snapshot is to be written

        rank: int
            Index number of the processor

        boost_dir: int (1 or -1)
            The direction of the Lorentz transformation from the lab frame
            to the boosted frame (along the z axis)
        """
        # Deduce the name of the filename where this snapshot writes
        if(lparallel == False):
            if rank == 0:
              self.filename = os.path.join( write_dir, 'hdf5/data%08d.h5' %i)
        else: 
            self.filename = os.path.join( write_dir, 'hdf5/data%08d.h5' %i)
        self.iteration = i

        # Time and boundaries in the lab frame (constants quantities)
        self.zmin_lab = zmin_lab
        self.zmax_lab = zmax_lab
        self.t_lab = t_lab
        self.boost_dir = boost_dir

        # Positions where the fields are to be registered
        # (Change at every iteration)
        self.current_z_lab = 0
        self.current_z_boost = 0

        # Buffered field slice and corresponding array index in z
        self.buffered_slices = []
        self.buffer_z_indices = []

        self.boost_dir = boost_dir

    def update_current_output_positions( self, t_boost, inv_gamma, inv_beta ):
        """
        Update the positions of output for this snapshot, so that
        if corresponds to the time t_boost in the boosted frame

        Parameters
        ----------
        t_boost: float (seconds)
            Time of the current iteration, in the boosted frame

        inv_gamma, inv_beta: floats
            Inverse of the Lorentz factor of the boost, and inverse
            of the corresponding beta.
        """
        t_lab = self.t_lab

        # This implements the Lorentz transformation formulas,
        # for a snapshot having a fixed t_lab
        self.current_z_boost = ( t_lab*inv_gamma - t_boost )*c*inv_beta
        self.current_z_lab = ( t_lab - t_boost*inv_gamma )*c*inv_beta

    def register_slice( self, slice_array, inv_dz_lab ):
        """
        Store the slice of fields represented by slice_array
        and also store the z index at which this slice should be
        written in the final lab frame array

        Parameters
        ----------
        slice_array: array of reals
            An array of packed fields that corresponds to one slice,
            as given by the SliceHandler object

        inv_dz_lab: float
            Inverse of the grid spacing in z, *in the lab frame*
        """
        # Find the index of the slice in the lab frame
        iz_lab = int(round( (self.current_z_lab - self.zmin_lab)*inv_dz_lab ))

        # Store the slice, if it was not already previously stored
        # (when dt is small and dz is large, this can happen)
        if (iz_lab in self.buffer_z_indices) == False:
            self.buffered_slices.append( slice_array )
            self.buffer_z_indices.append( iz_lab )

    def compact_slices(self):
        """
        Compact the successive slices that have been buffered
        over time into a single array, and return the indices
        at which this array should be written.

        Returns
        -------
        field_array: an array of reals of shape
        - (10, em.nxlocal+1, nslices) if dim is "2d"
        - (10, em.nxlocal+1, em.nylocal+1, nslices) if dim is "3d"
        - (10, 2*em.circ_m+1, em.nxlocal+1, nslices) if dim is "circ"
        In the above nslices is the number of buffered slices

        iz_min, iz_max: integers
        The indices between which the slices should be written
        (iz_min is inclusive, iz_max is exclusive)

        Returns None if the slices are empty
        """
        # Return None if the slices are empty
        if len(self.buffer_z_indices) == 0:
            return( None, None, None )

        # Check that the indices of the slices are contiguous
        # (This should be a consequence of the transformation implemented
        # in update_current_output_positions, and of the calculation
        # of inv_dz_lab.)
        iz_old = self.buffer_z_indices[0]
        for iz in self.buffer_z_indices[1:]:
            if iz != iz_old - self.boost_dir:
                raise UserWarning('In the boosted frame diagnostic, '
                        'the buffered slices are not contiguous in z.\n'
                        'The boosted frame diagnostics may be inaccurate.')
                break
            iz_old = iz

        # Pack the different slices together
        if self.boost_dir == 1:
            # Reverse the order of the slices when stacking the array,
            # since the slices where registered for right to left
            try:
                field_array = np.stack( self.buffered_slices[::-1], axis=-1 )
            except AttributeError:
                # If the version of numpy is older than 1.10, stack does
                # not exist. In this case, do it by hand:
                index  = np.array(np.shape( self.buffered_slices[::-1] ))
                rolled_index = np.roll(index, -1)
                field_array = np.dstack( self.buffered_slices[::-1] )
                field_array = field_array.reshape( (rolled_index), order="F" )
        elif self.boost_dir == -1:
            try:
                field_array = np.stack( self.buffered_slices, axis=-1 )
            except AttributeError:
                # If the version of numpy is older than 1.10, stack does
                # not exist. In this case, do it by hand:
                index  = np.array(np.shape( self.buffered_slices ))
                rolled_index = np.roll(index, -1)
                field_array = np.dstack( self.buffered_slices )
                field_array = field_array.reshape( (rolled_index), order="F" )

        # Get the first and last index in z
        # (Following Python conventions, iz_min is inclusive,
        # iz_max is exclusive)
        if self.boost_dir == 1:
            iz_min = self.buffer_z_indices[-1]
            iz_max = self.buffer_z_indices[0] + 1
        elif self.boost_dir == -1:
            iz_min = self.buffer_z_indices[0]
            iz_max = self.buffer_z_indices[-1] + 1

        return( field_array, iz_min, iz_max )

class SliceHandler:
    """
    Class that extracts, Lorentz-transforms and writes slices of the fields
    """
    def __init__( self, gamma_boost, beta_boost, dim, start, nx_dump, ny_dump ):
        """
        Initialize the SliceHandler object

        Parameters
        ----------
        gamma_boost, beta_boost: float
            The Lorentz factor of the boost and the corresponding beta

        dim: string
            Either "2d", "3d", or "circ"
            Indicates the geometry of the fields
        """
        # Store the arguments
        self.dim = dim
        self.gamma_boost = gamma_boost
        self.beta_boost = beta_boost
        self.start = start
        self.nx_dump = nx_dump
        self.ny_dump = ny_dump

        # Create a dictionary that contains the correspondance
        # between the field names and array index
        if (dim=="2d") or (dim=="3d") :
            self.field_to_index = {'Ex':0, 'Ey':1, 'Ez':2, 'Bx':3,
                'By':4, 'Bz':5, 'Jx':6, 'Jy':7, 'Jz':8, 'rho':9}
        elif dim=="circ":
            self.field_to_index = {'Er':0, 'Et':1, 'Ez':2, 'Br':3,
                'Bt':4, 'Bz':5, 'Jr':6, 'Jt':7, 'Jz':8, 'rho':9}

    def extract_slice( self, em, z_boost, zmin_boost ):
        """
        Returns an array that contains the slice of the fields at
        z_boost (the fields returned are already transformed to the lab frame)

        Parameters
        ----------
        em: an EM3DSolver object
            The object from which to extract the fields

        z_boost: float (meters)
            Position of the slice in the boosted frame

        zmin_boost: float (meters)
            Position of the left end of physical part of the local subdomain
            (i.e. excludes guard cells)

        Returns
        -------
        An array of reals that packs together the slices of the
        different fields.

        The first index of this array corresponds to the field type
        (10 different field types), and the correspondance
        between the field type and integer index is given self.field_to_index

        The shape of this arrays is:
        - (10, em.nxlocal+1,) for dim="2d"
        - (10, em.nxlocal+1, em.nylocal+1) for dim="3d"
        - (10, 2*em.circ_m+1, em.nxlocal+1) for dim="circ"
        """
        # Extract a slice of the fields *in the boosted frame*
        # at z_boost, using interpolation, and store them in an array
        # (See the docstring of the extract_slice_boosted_frame for
        # the shape of this array.)
        slice_array = self.extract_slice_boosted_frame(
            em, z_boost, zmin_boost )

        # Perform the Lorentz transformation of the fields *from
        # the boosted frame to the lab frame*
        self.transform_fields_to_lab_frame( slice_array,em )

        return( slice_array )

    def extract_slice_boosted_frame( self, em, z_boost, zmin_boost ):
        """
        Extract a slice of the fields at z_boost, using interpolation in z

        See the docstring of extract_slice for the parameters.

        Returns
        -------
        An array that packs together the slices of the different fields.
            The shape of this arrays is:
            - (10, em.nxlocal+1,) for dim="2d"
            - (10, em.nxlocal+1, em.nylocal+1) for dim="3d"
            - (10, 2*em.circ_m+1, em.nxlocal+1) for dim="circ"
        """
        # Allocate an array of the proper shape
        if(self.nx_dump is None): nx_ = em.nxlocal+1
        else : nx_ = self.nx_dump
        if(self.ny_dump is None): ny_ = em.nylocal+1
        else : ny_ = self.ny_dump
        if self.dim=="2d":
            slice_array = np.empty( (10, nx_,) ,order="F")
        elif self.dim=="3d":
            slice_array = np.empty( (10, nx_, ny_), order = "F" )
        elif self.dim=="circ":
            slice_array = np.empty( (10, 2*em.circ_m+1, em.nxlocal+1) ,order="F")

        # Find the index of the slice in the boosted frame
        # and the corresponding interpolation shape factor
        dz = em.dz
        # Centered
        z_centered_gridunits = ( z_boost - zmin_boost )/dz
        iz_centered = int( z_centered_gridunits )
        Sz_centered = iz_centered + 1 - z_centered_gridunits
        # Staggered
        z_staggered_gridunits = ( z_boost - zmin_boost - 0.5*dz )/dz
        iz_staggered = int( z_staggered_gridunits )
        Sz_staggered = iz_staggered + 1 - z_staggered_gridunits

        # Shortcut for the correspondance between field and integer index
        f2i = self.field_to_index

        # Loop through the fields, and extract the proper slice for each field
        for quantity in self.field_to_index.keys():
            # Here typical values for `quantity` are e.g. 'Er', 'Bx', 'rho'

            # Choose the index and interpolating factor, depending
            # on whether the field is centered in z or staggered
            # - Centered field in z
            if z_offset_dict[quantity] == 0:
                iz = iz_centered
                Sz = Sz_centered
            # - Staggered field in z
            elif z_offset_dict[quantity] == 0.5:
                iz = iz_staggered
                Sz = Sz_staggered
            else:
                raise ValueError( 'Unknown staggered offset for %s: %f' %(
                    quantity, z_offset_dict[quantity] ))

            # Interpolate the centered field in z
            # (Transversally-staggered fields are also interpolated
            # to the nodes of the grid, thanks to the flag transverse_centered)
            slice_array[ f2i[quantity], ... ] = Sz * get_dataset(
                self.dim, em, quantity, lgather=False, iz_slice=iz,
                transverse_centered=True, start = self.start, nx_d = self.nx_dump,  ny_d = self.ny_dump )
            slice_array[ f2i[quantity], ... ] += (1.-Sz) * get_dataset(
                self.dim, em, quantity, lgather=False, iz_slice=iz+1,
                transverse_centered=True, start = self.start, nx_d = self.nx_dump,  ny_d = self.ny_dump )

        return( slice_array )

    def transform_fields_to_lab_frame( self, fields ,em):
        """
        Modifies the array `fields` in place, to transform the field values
        from the boosted frame to the lab frame.

        The transformation is a transformation with -beta_boost, thus
        the corresponding formulas are:
        - for the transverse part of E and B:
        $\vec{E}_{lab} = \gamma(\vec{E} - c\vec{\beta} \times\vec{B})$
        $\vec{B}_{lab} = \gamma(\vec{B} + \vec{\beta}/c \times\vec{E})$
        - for rho and Jz:
        $\rho_{lab} = \gamma(\rho + \beta J_{z}/c)$
        $J_{z,lab} = \gamma(J_z + c\beta \rho)$

        Parameter
        ---------
        fields: array of floats
            An array that packs together the slices of the different fields.
            The shape of this arrays is:
            - (10, em.nxlocal+1,) for dim="2d"
            - (10, em.nxlocal+1, em.nylocal+1) for dim="3d"
            - (10, 2*em.circ_m+1, em.nxlocal+1) for dim="circ"
        """
        # Some shortcuts
        gamma = self.gamma_boost
        cbeta = c*self.beta_boost
        beta_c = self.beta_boost/c
        # Shortcut to give the correspondance between field name
        # (e.g. 'Ex', 'rho') and integer index in the array
        f2i = self.field_to_index

        # Lorentz transformations
        # For E and B
        # (NB: Ez and Bz are unchanged by the Lorentz transform)
        if(self.nx_dump is None) : n2 = em.nxlocal+1
        else : n2 = self.nx_dump
        if(self.dim  == "3d"):
            if(self.ny_dump is None) : n3 = em.nylocal+1
            else: n3 = self.ny_dump

        if self.dim in ["2d", "3d"]:

            if (hasattr(em,"l_pxr")) :
                if(em.l_pxr == True):
                    n1 = 10
                    if(self.dim == "2d"):
                        em.lorentz_transform2d(n1,n2,fields,gamma,cbeta,beta_c)
                    else: 
                        em.lorentz_transform3d(n1,n2,n3,fields,gamma,cbeta,beta_c)         
            else :
                # Use temporary arrays when changing Ex and By in place
                ex_lab = gamma*( fields[f2i['Ex']] + cbeta * fields[f2i['By']] )
                by_lab = gamma*( fields[f2i['By']] + beta_c * fields[f2i['Ex']] )
                fields[ f2i['Ex'], ... ] = ex_lab
                fields[ f2i['By'], ... ] = by_lab
                # Use temporary arrays when changing Ey and Bx in place
                ey_lab = gamma*( fields[f2i['Ey']] - cbeta * fields[f2i['Bx']] )
                bx_lab = gamma*( fields[f2i['Bx']] - beta_c * fields[f2i['Ey']] )
                fields[ f2i['Ey'], ... ] = ey_lab
                fields[ f2i['Bx'], ... ] = bx_lab
                # For rho and J
                # (NB: the transverse components of J are unchanged)
                # Use temporary arrays when changing rho and Jz in place

                rho_lab = gamma*( fields[f2i['rho']] + beta_c * fields[f2i['Jz']] )
                Jz_lab =  gamma*( fields[f2i['Jz']] + cbeta * fields[f2i['rho']] )
                fields[ f2i['rho'], ... ] = rho_lab
                fields[ f2i['Jz'], ... ] = Jz_lab
    
        elif self.dim=="circ":
            # Use temporary arrays when changing Er and Bt in place
            er_lab = gamma*( fields[f2i['Er']] + cbeta * fields[f2i['Bt']] )
            bt_lab = gamma*( fields[f2i['Bt']] + beta_c * fields[f2i['Er']] )
            fields[ f2i['Er'], ... ] = er_lab
            fields[ f2i['Bt'], ... ] = bt_lab
            # Use temporary arrays when changing Et and Br in place
            et_lab = gamma*( fields[f2i['Et']] - cbeta * fields[f2i['Br']] )
            br_lab = gamma*( fields[f2i['Br']] - beta_c * fields[f2i['Et']] )
            fields[ f2i['Et'], ... ] = et_lab
            fields[ f2i['Br'], ... ] = br_lab

            # For rho and J
            # (NB: the transverse components of J are unchanged)
            # Use temporary arrays when changing rho and Jz in place

            rho_lab = gamma*( fields[f2i['rho']] + beta_c * fields[f2i['Jz']] )
            Jz_lab =  gamma*( fields[f2i['Jz']] + cbeta * fields[f2i['rho']] )
            fields[ f2i['rho'], ... ] = rho_lab
            fields[ f2i['Jz'], ... ] = Jz_lab

