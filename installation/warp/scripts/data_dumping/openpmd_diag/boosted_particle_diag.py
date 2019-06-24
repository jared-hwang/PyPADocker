"""
This file defines the class BoostedParticleDiagnostic

Major features:
- The class reuses the existing methods of ParticleDiagnostic
  as much as possible, through class inheritance
- The class implements memory buffering of the slices, so as
  not to write to disk at every timestep
"""
import os
import numpy as np
import time
from scipy.constants import c
from particle_diag import ParticleDiagnostic
from warp_parallel import me, mpiallgather
try:
    from mpi4py import MPI
except ImportError:
    MPI = None
    pass

class BoostedParticleDiagnostic(ParticleDiagnostic):
    """
    Class that writes the particles *in the lab frame*,
    from a simulation in the boosted frame

    Usage
    -----
    After initialization, the diagnostic is called by using
    the 'write' method.
    """
    def __init__(self, zmin_lab, zmax_lab, v_lab, dt_snapshots_lab,
                 Ntot_snapshots_lab, gamma_boost, period,
                 em, top, w3d, comm_world=None,
                 particle_data=["position", "momentum", "weighting"],
                 select=None, write_dir=None,
                 species={"electrons": None}, boost_dir=1,lparallel_output=False,t_min_lab=0. ):
        """
        Initialize diagnostics that retrieve the data in the lab frame,
        as a series of snapshot (one file per snapshot),
        within a virtual moving window defined by zmin_lab, zmax_lab, v_lab.


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

        boost_dir: int (1 or -1)
            The direction of the Lorentz transformation from the lab frame
            to the boosted frame (along the z axis)

        See the documentation of ParticleDiagnostic for the other parameters
        """
        # Do not leave write_dir as None, as this may conflict with
        # the default directory ('./diags') in which diagnostics in the
        # boosted frame are written
        if write_dir is None:
            write_dir = 'lab_diags'
        self.lparallel_output = lparallel_output
        # Initialize Particle diagnostic normal attributes
        ParticleDiagnostic.__init__(self, period, top, w3d, comm_world,
            species=species, particle_data=particle_data, select=select,
            write_dir=write_dir,write_metadata_parallel=lparallel_output, lparallel_output = lparallel_output)
        # Note: The boosted frame diagnostics cannot use parallel HDF5 output

        # Check user input
        boost_dir = int(boost_dir)
        assert boost_dir in [1,-1]
        
        # Register the boost quantities
        self.em = em
        self.gamma_boost = gamma_boost
        self.inv_gamma_boost = 1./gamma_boost
        self.beta_boost = np.sqrt(1. - self.inv_gamma_boost**2) * boost_dir
        self.inv_beta_boost = 1./self.beta_boost
        self.Ntot_snapshots_lab = Ntot_snapshots_lab
        self.t_min_lab = t_min_lab
        self.dump_p_fields = False
        if("E" in self.particle_data or "B" in self.particle_data): self.dump_p_fields = True

        # Create the list of LabSnapshot objects
        self.snapshots = []

        # Record the time it takes
        if self.rank == 0:
            measured_start = time.clock()
            print('\nInitializing the lab-frame diagnostics: %d files...' %(
                Ntot_snapshots_lab) )

        # Loop through the lab snapshots and create the corresponding files
        self.particle_catcher = ParticleCatcher(
            self.gamma_boost, self.beta_boost, top,self.dump_p_fields,em )
        self.particle_catcher.allocate_previous_instant()

        for i in range( Ntot_snapshots_lab ):
            t_lab = i*dt_snapshots_lab + self.t_min_lab
            snapshot = LabSnapshot( t_lab, zmin_lab + v_lab*t_lab,
                            top.dt, zmax_lab + v_lab*t_lab,
                            self.write_dir, i, self.species_dict, self.rank,self.dump_p_fields,self.em )
            self.snapshots.append( snapshot )
            # Initialize a corresponding empty file
            self.create_file_empty_particles(
                snapshot.filename, i, snapshot.t_lab, self.top.dt)
 

        if(self.lparallel_output) :
            self.mpi_group = self.comm_world.Get_group() 

        # Print a message that records the time for initialization
        if self.rank == 0:
            measured_end = time.clock()
            print('Time taken for initialization of the files: %.5f s' %(
                measured_end - measured_start) )


    def write( self ):
        """
        Redefines the method write of the parent class ParticleDiagnostic

        Should be registered with installafterstep in Warp
        """
        # At each timestep, store a slice of the particles in memory buffers
        self.store_snapshot_slices()

        # Every self.period, write the buffered slices to disk
        if self.top.it % self.period == 0:
            if not self.lparallel_output:
                self.flush_to_disk()
            else:
                self.flush_to_disk_parallel() 

    def store_snapshot_slices( self ):
        """
        Store slices of the particles in the memory buffers of the
        corresponding lab snapshots

        The particles (in a slice) are also selected in accordance with
        the selection rules provided as argument to BoostedParticleDiagnostic
        """
        # Loop through the labsnapshots
        for snapshot in self.snapshots:

            # Update the positions of the output slice of this snapshot
            # in the lab and boosted frame (current_z_lab and current_z_boost)
            snapshot.update_current_output_positions( self.top.time,
                            self.inv_gamma_boost, self.inv_beta_boost)

            # For this snapshot:
            # - check if the output position *in the lab frame*
            #   is within the lab-frame boundaries of the current snapshot
            if ( (snapshot.current_z_lab > snapshot.zmin_lab) and \
                 (snapshot.current_z_lab < snapshot.zmax_lab) ):

                # Loop through the particle species and register the
                # particle arrays in the snapshot objects (buffering)
                for species_name, species in self.species_dict.iteritems():

                    slice_array = self.particle_catcher.extract_slice(
                        species, self.select, snapshot.prev_z_boost,
                        snapshot.current_z_boost, snapshot.t_lab )
                    snapshot.register_slice( slice_array, species_name )
  
    def flush_to_disk_parallel(self):
        """
        Writes the buffered slices of particles to the disk using parallel IO. Erase the
        buffered slices of the LabSnapshot objects
        """


        #Init empty ntot_snapshot arrays of h5_file dictionnaries 
        #The dictionary is incremented with species_names
        
        f = [None]*self.Ntot_snapshots_lab


        #Init empty ntot_snapshot arrays of particle_arrays dictionnaries 
        #The dictionary is incremented with species_names
        particle_array = [None]*self.Ntot_snapshots_lab
 
        #Init empty ntot_snapshot arrays of boolean dictionnaries 
        #The dictionary is incremented with species_names
        write_on =[None]*self.Ntot_snapshots_lab
  
        #Init empty ntot_snapshot arrays of integer_array dictionnaries 
        #The dictionary is incremented with species_names
        n_rank =[None]*self.Ntot_snapshots_lab


        #Init empty ntot_snapshot arrays of integer dictionnaries 
        #The dictionary is incremented with species_names
        #local_number number of particles to be dumped 
        nlocals_dict = [None]*self.Ntot_snapshots_lab


        #Init empty ntot_snapshot arrays of ineteger dictionnaries 
        #The dictionary is incremented with species_names
        #Global number of particles to be dumped 
        nglobal_dict = [None]*self.Ntot_snapshots_lab


        #Init empty ntot_snapshot array of MPI.COMM dictionnaries 
        #The dictionary is incremented with species_names
        #For each snapshot and for each species, dump_comm is passed to h5py to dump data
        dump_comm = [None]*self.Ntot_snapshots_lab

        for i in range(self.Ntot_snapshots_lab):
                f[i] = dict()
                particle_array[i] = dict()
                write_on[i] = dict()
                n_rank[i] = dict()
                dump_comm[i] = dict()       
                nlocals_dict[i] = dict()
                nglobal_dict[i] = dict()
                for species_name in self.species_dict:
                    f[i][species_name] = 0
                    particle_array[i][species_name] = None 
                    write_on[i][species_name] = False
                    n_rank[i][species_name] = None 
                    dump_comm[i][species_name] =  MPI.COMM_NULL
                    nlocals_dict[i][species_name] = None
                    nglobal_dict[i][species_name] = None  





        # Loop over snapshots and species to construct dump_comms and open h5_files
        for snapshot in self.snapshots:
            i = snapshot.iteration 
            for species_name in self.species_dict:

                particle_array[i][species_name] = snapshot.compact_slices(species_name)
                n = np.size(particle_array[i][species_name][0])
                nlocals_dict[i][species_name]= mpiallgather( n )
                nglobal_dict[i][species_name]=np.sum(nlocals_dict[i][species_name])

                n_rank[i][species_name] = nlocals_dict[i][species_name]

                # Create a communicator containing ranks that have
                # particles to dump
                in_list = -1
                if (np.shape(particle_array[i][species_name])[1] != 0) :
                    in_list = me
                    write_on[i][species_name] = True
                ranks_group_list = mpiallgather( in_list )
                # deletes -1 from the list of ranks
                ranks_group_list = [x for x in ranks_group_list if x >= 0 ]
                ranks_group_list = list(set(ranks_group_list))
                newgroup = self.mpi_group.Incl(ranks_group_list)
                # Create communicator
                dump_comm[i][species_name] = self.comm_world.Create(newgroup)
                newgroup.Free()
                ranks_group_list = []
                if(write_on[i][species_name]): 
                    if(dump_comm[i][species_name] is not None and dump_comm[i][species_name] != MPI.COMM_NULL  ):
                        #each MPI opens relevent h5 files with adequate dump_com 
                        f[i][species_name]=self.open_file(snapshot.filename, parallel_open= \
                                                          self.lparallel_output,comm=dump_comm[i][species_name] )
                 
         
   
        # loop over snapshots and species to flush data
        for snapshot in self.snapshots:
            i= snapshot.iteration 
            for species_name in self.species_dict:
                if(write_on[i][species_name]):
                      n_part_to_dump = np.shape(particle_array[i][species_name])[1]
                      self.write_slices(particle_array[i][species_name], species_name, snapshot,
                                        self.particle_catcher.particle_to_index,
                                        comm=dump_comm[i][species_name],n_rank=n_rank[i][species_name],
                                        n_global=nglobal_dict[i],h5_file=f[i][species_name])
                      snapshot.buffered_slices[species_name] = []

        # loop over snapshots and species to  close files and free dump_comm
        for i in range(self.Ntot_snapshots_lab):
            for species_name in self.species_dict: 
                if(write_on[i][species_name]):
                    f[i][species_name].close()
                if dump_comm[i][species_name] != MPI.COMM_NULL:
                    dump_comm[i][species_name].Free()


       # cleaning 
        f = []
        dump_comm = []        
        particle_array = []
        write_on = []
        n_rank = [] 
        nlocals_dict = []
        nglobal_dict = []
      




    def flush_to_disk(self):
        """
        Writes the buffered slices of particles to the disk. Erase the
        buffered slices of the LabSnapshot objects

        Notice: In parallel version, data are gathered to proc 0
        before being saved to disk
        """
        # Loop through the labsnapshots and flush the data

        for snapshot in self.snapshots:

            # Compact the successive slices that have been buffered
            # over time into a single array
            for species_name in self.species_dict:
                particle_array = snapshot.compact_slices(species_name)
                if self.comm_world is not None:
                    # Create a communicator containing ranks that have
                    # particles to dump
                    in_list = 0
                    if (np.shape(particle_array)[1] != 0) or (me == 0):
                        in_list = me
                    ranks_group_list = mpiallgather( in_list )
                    ranks_group_list = list(set(ranks_group_list))
                    mpi_group = self.comm_world.Get_group()
                    self.ranks_group_list = ranks_group_list
                    # Create group
                    newgroup = mpi_group.Incl(ranks_group_list)
                    # Create communicator
                    dump_comm = self.comm_world.Create(newgroup)
                    # Gather data on proc 0 into this communicator
                    if dump_comm != MPI.COMM_NULL:
                        list_part_array = dump_comm.gather( particle_array )
                    # Free the dump communicator
                    mpi_group.Free()
                    newgroup.Free()
                    if dump_comm != MPI.COMM_NULL:
                        dump_comm.Free()

                    # Rank 0 concatenates all lists into an array with all particles
                    if self.rank == 0:
                        p_array = np.concatenate(list_part_array, axis=1)
                else:
                    p_array = particle_array
                    
                # Write this array to disk (if this snapshot has new slices)
                if self.rank == 0 and p_array.size:
                    self.write_slices(p_array, species_name, snapshot,
                        self.particle_catcher.particle_to_index)
                # Erase the buffers
                snapshot.buffered_slices[species_name] = []

    def write_boosted_dataset(self, species_grp, path, data, quantity,n_rank=None,n_global=None,comm=None):
        """
        Writes each quantity of the buffered dataset to the disk, the
        final step of the writing
        """
        if not self.lparallel_output: 
            dset = species_grp[path]
            index = dset.shape[0]

            # Resize the h5py dataset
            dset.resize(index+len(data), axis=0)

            # Write the data to the dataset at correct indices
            dset[index:] = data
        else:
            dset = species_grp[path] 
            index = dset.shape[0]
            dset.resize(index+n_global, axis=0)            

            if n_rank is not None:
                iold = index+sum(n_rank[0:self.rank])
                # Calculate the last index occupied by the current rank
                inew = iold+n_rank[self.rank]
                # Write the local data to the global array
                dset[iold:inew] = data
            #One proc writes the data (serial and lparallel_output=False)
            else:
                if (self.rank==0):
                    # Write the data to the dataset at correct indices
                    dset[index:] = data

              

    def write_slices( self, particle_array, species_name, snapshot, p2i,n_rank=None,n_global=None,comm=None,h5_file=None ):
        """
        For one given snapshot, write the slices of the
        different species to an openPMD file

        Parameters
        ----------
        particle_array: array of reals
            Array of shape (8, num_part)

        species_name: String
            A String that acts as the key for the buffered_slices dictionary

        snapshot: a LabSnaphot object

        p2i: dict
            Dictionary of correspondance between the particle quantities
            and the integer index in the particle_array
 
        n_rank:  array of integer
            array of local number of particles that need to be dumped for each mpi task

        n_global: integer
            global number of particles  for this flush.

        comm : MPI_COMMUNICATOR
             mpi communicator used to dump the data when parallel IO

        h5_file : An h5py.File object
            if parallel IO then  this routine does not open h5 files, 
            instead h5 files are opened in flush_to_disk_parallel, and h5_file the returned object from open_file 
            for current snapshot 

        """
        
        # Open the file without parallel I/O in this implementation
        # If using parallel IO then files have already been opened
        if not self.lparallel_output:
           f = self.open_file( snapshot.filename, parallel_open=False)
        else: 
           f = h5_file
        particle_path = "/data/%d/particles/%s" %(snapshot.iteration,
                                                    species_name)
        species_grp = f[particle_path]

        if(self.lparallel_output): 
            ng=n_global[species_name] 
        else: 
            ng = None   
        # Loop over the different quantities that should be written
        for particle_var in self.particle_data:

            if particle_var == "position":
                for coord in ["x","y","z"]:
                    quantity= coord
                    path = "%s/%s" %(particle_var, quantity)
                    data = particle_array[ p2i[ quantity ] ]
                    self.write_boosted_dataset(species_grp, path, data, quantity,n_rank=n_rank,\
                                               n_global=ng,comm=comm)


            elif particle_var == "momentum":
                for coord in ["x","y","z"]:
                    quantity= "u%s" %coord
                    path = "%s/%s" %(particle_var,coord)
                    data = particle_array[ p2i[ quantity ] ]
                    self.write_boosted_dataset(species_grp, path, data, quantity,n_rank=n_rank,\
                                               n_global=ng,comm=comm)

            elif particle_var == "B":
                for coord in ["x","y","z"]:
                    quantity= "b%s" %coord
                    path = "%s/%s" %(particle_var,coord)
                    data = particle_array[ p2i[ quantity ] ]
                    self.write_boosted_dataset(species_grp, path, data, quantity,n_rank=n_rank,\
                                               n_global=ng,comm=comm)

            elif particle_var == "E":
                for coord in ["x","y","z"]:
                    quantity= "e%s" %coord
                    path = "%s/%s" %(particle_var,coord)
                    data = particle_array[ p2i[ quantity ] ]
                    self.write_boosted_dataset(species_grp, path, data, quantity,n_rank=n_rank,\
                                               n_global=ng,comm=comm)


            elif particle_var == "weighting":
               quantity= "w"
               path = 'weighting'
               data = particle_array[ p2i[ quantity ] ]
               self.write_boosted_dataset(species_grp, path, data, quantity,n_rank=n_rank,\
                                          n_global=ng,comm=comm)


            elif particle_var == "id":
               quantity= "id"
               path = 'id'
               data = particle_array[ p2i[ quantity ] ]
               self.write_boosted_dataset(species_grp, path, data, quantity,n_rank=n_rank,\
                                          n_global=ng,comm=comm)

        data = []

        #If serial IO then close files here
        
        if not self.lparallel_output:
            f.close()


class LabSnapshot:
    """
    Class that stores data relative to one given snapshot
    in the lab frame (i.e. one given *time* in the lab frame)
    """
    def __init__(self, t_lab, zmin_lab, dt, zmax_lab, write_dir, i,
        species_dict, rank, dump_p_fields, em):
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

        species_dict: dict
            Contains all the species name of the species object
            (inherited from Warp)

        dump_p_fields: boolean
            Flag to dump particle fields
        """
        # Deduce the name of the filename where this snapshot writes
        self.filename = os.path.join( write_dir, 'hdf5/data%08d.h5' %i)
        self.iteration = i
        self.dt = dt
        self.dump_p_fields = dump_p_fields
        self.n_quantities = 0

        if(dump_p_fields) : 
            self.n_quantities = 15
        else: 
            self.n_quantities = 9
 
        self.em = em
        # Time and boundaries in the lab frame (constants quantities)
        self.zmin_lab = zmin_lab
        self.zmax_lab = zmax_lab
        self.t_lab = t_lab

        # Positions where the fields are to be registered
        # (Change at every iteration)
        self.current_z_lab = 0
        self.current_z_boost = 0

        # Prepare buffered slices
        self.buffered_slices = {}
        for species_name in species_dict:
            self.buffered_slices[species_name] = []

    def update_current_output_positions( self, t_boost, inv_gamma, inv_beta ):
        """
        Update the current and previous positions of output for this snapshot,
        so that it corresponds to the time t_boost in the boosted frame

        Parameters
        ----------
        t_boost: float (seconds)
            Time of the current iteration, in the boosted frame

        inv_gamma, inv_beta: floats
            Inverse of the Lorentz factor of the boost, and inverse
            of the corresponding beta
        """
        # Some shorcuts for further calculation's purposes
        t_lab = self.t_lab
        t_boost_prev = t_boost - self.dt

        # This implements the Lorentz transformation formulas,
        # for a snapshot having a fixed t_lab
        self.current_z_boost = (t_lab*inv_gamma - t_boost)*c*inv_beta
        self.prev_z_boost = (t_lab*inv_gamma - t_boost_prev)*c*inv_beta
        self.current_z_lab = (t_lab - t_boost*inv_gamma)*c*inv_beta
        self.prev_z_lab = (t_lab - t_boost_prev*inv_gamma)*c*inv_beta

    def register_slice(self, slice_array, species):
        """
        Store the slice of particles represented by slice_array

        Parameters
        ----------
        slice_array: array of reals
            An array of packed fields that corresponds to one slice,
            as given by the ParticleCatcher object

        species: String, key of the species_dict
            Act as the key for the buffered_slices dictionary
        """
        # Store the values
        self.buffered_slices[species].append(slice_array)

    def compact_slices(self, species):
        """
        Compact the successive slices that have been buffered
        over time into a single array.

        Parameters
        ----------
        species: String, key of the species_dict
            Act as the key for the buffered_slices dictionary

        Returns
        -------
        paticle_array: an array of reals of shape (9, numPart)
        regardless of the dimension

        Returns None if the slices are empty
        """
        if self.buffered_slices[species] != []:
            particle_array = np.concatenate(
                self.buffered_slices[species], axis=1)
        else:
            particle_array = np.empty((self.n_quantities,0))

        return particle_array

class ParticleCatcher:
    """
    Class that extracts, Lorentz-transforms and gathers particles
    """
    def __init__(self, gamma_boost, beta_boost, top,dump_f=False, em=None):
        """
        Initialize the ParticleCatcher object

        Parameters
        ----------
        gamma_boost, beta_boost: float
            The Lorentz factor of the boost and the corresponding beta

        top: WARP object

        dump_f : boolean
            Flag for field dumping
        
        em : EM Object
  
       
        """
        # Some attributes neccessary for particle selections
        self.gamma_boost = gamma_boost
        self.beta_boost = beta_boost
        self.top = top
        self.dump_p_fields = dump_f
        self.em = em

        # Create a dictionary that contains the correspondance
        # between the particles quantity and array index
        if(dump_f == False):
            self.particle_to_index = {'x':0, 'y':1, 'z':2, 'ux':3,
                    'uy':4, 'uz':5, 'w':6, 'gamma':7, 't':8}
            if self.top.ssnpid > 0:
                self.particle_to_index['id'] = 9

        else: 
            self.particle_to_index = {'x':0, 'y':1, 'z':2, 'ux':3,
                    'uy':4, 'uz':5, 'w':6, 'gamma':7, 't':8,'ex':9, 'ey':10, 'ez':11,'bx':12, 'by':13, 'bz':14} 
            if self.top.ssnpid > 0:
                self.particle_to_index['id'] = 15
    def get_particle_slice( self, species, prev_z_boost, current_z_boost ):
        """
        Select the particles for the current slice, and extract their
        positions and momenta at the current and previous timestep

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which one slice will be extracted

        current_z_boost, prev_z_boost: floats
            The position, in the boosted frame, of the plane that
            corresponds to the lab snapshot, at the current iteration and
            at the previous iteration.
            Only particles that have crossed this plane, between the previous
            iteration and the current iteration, are extracted

        Returns
        -------
        num_part: int
            Number of selected particles
        """
        # Quantities at current time step
        current_x = self.get_quantity( species, "x" )
        current_y = self.get_quantity( species, "y" )
        current_z = self.get_quantity( species, "z" )
        current_ux = self.get_quantity( species, "ux" )
        current_uy = self.get_quantity( species, "uy" )
        current_uz = self.get_quantity( species, "uz" )
        current_weights = self.get_quantity( species, "w" )
        if(self.dump_p_fields): 
            current_ex = self.get_quantity( species, "ex" )
            current_ey = self.get_quantity( species, "ey" )
            current_ez = self.get_quantity( species, "ez" )
            current_bx = self.get_quantity( species, "bx" )
            current_by = self.get_quantity( species, "by" )
            current_bz = self.get_quantity( species, "bz" )
        
        if self.top.ssnpid > 0:
            current_id = self.get_quantity( species, "id" )

        # Quantities at previous time step
        previous_x = self.get_quantity( species, "x", l_prev=True )
        previous_y = self.get_quantity( species, "y", l_prev=True )
        previous_z = self.get_quantity( species, "z", l_prev=True )
        previous_ux = self.get_quantity( species, "ux", l_prev=True )
        previous_uy = self.get_quantity( species, "uy", l_prev=True )
        previous_uz = self.get_quantity( species, "uz", l_prev=True )
        if(self.dump_p_fields):
            previous_ex = self.get_quantity( species, "ex",l_prev=True )
            previous_ey = self.get_quantity( species, "ey",l_prev=True )
            previous_ez = self.get_quantity( species, "ez",l_prev=True )
            previous_bx = self.get_quantity( species, "bx",l_prev=True )
            previous_by = self.get_quantity( species, "by",l_prev=True )
            previous_bz = self.get_quantity( species, "bz",l_prev=True )


        # A particle array for mapping purposes
        particle_indices = np.arange( len(current_z) )

        # For this snapshot:
        # - check if the output position *in the boosted frame*
        #   crosses the zboost in a forward motion
        # - check if the output position *in the boosted frame*
        #   crosses the zboost_prev in a backward motion
        selected_indices = np.compress((
            ((current_z >= current_z_boost ) & (previous_z <= prev_z_boost )) |
            ((current_z <= current_z_boost ) & (previous_z >= prev_z_boost ))
            ), particle_indices)

        num_part = np.shape(selected_indices)[0]

        ## Particle quantities that satisfy the aforementioned condition
        self.mass = species.mass

        self.x_captured = np.take(current_x, selected_indices)
        self.y_captured = np.take(current_y, selected_indices)
        self.z_captured = np.take(current_z, selected_indices)
        self.ux_captured = np.take(current_ux, selected_indices)
        self.uy_captured = np.take(current_uy, selected_indices)
        self.uz_captured = np.take(current_uz, selected_indices)
        self.w_captured = np.take(current_weights, selected_indices)
        self.gamma_captured = np.sqrt(1. + (self.ux_captured**2+\
            self.uy_captured**2 + self.uz_captured**2)/c**2)
        if self.top.ssnpid > 0:
            self.id_captured = np.take(current_id, selected_indices)
        if(self.dump_p_fields):
            self.ex_captured = np.take(current_ex, selected_indices)
            self.ey_captured = np.take(current_ey, selected_indices)
            self.ez_captured = np.take(current_ez, selected_indices)
            self.bx_captured = np.take(current_bx, selected_indices)
            self.by_captured = np.take(current_by, selected_indices)
            self.bz_captured = np.take(current_bz, selected_indices) 

        self.x_prev_captured = np.take(previous_x, selected_indices)
        self.y_prev_captured = np.take(previous_y, selected_indices)
        self.z_prev_captured = np.take(previous_z, selected_indices)
        self.ux_prev_captured = np.take(previous_ux, selected_indices)
        self.uy_prev_captured = np.take(previous_uy, selected_indices)
        self.uz_prev_captured = np.take(previous_uz, selected_indices)
        self.gamma_prev_captured = np.sqrt(1. + (self.ux_prev_captured**2+\
            self.uy_prev_captured**2 + self.uz_prev_captured**2)/c**2)
        if(self.dump_p_fields):
            self.ex_prev_captured = np.take(previous_ex, selected_indices)
            self.ey_prev_captured = np.take(previous_ey, selected_indices)
            self.ez_prev_captured = np.take(previous_ez, selected_indices)
            self.bx_prev_captured = np.take(previous_bx, selected_indices)
            self.by_prev_captured = np.take(previous_by, selected_indices)
            self.bz_prev_captured = np.take(previous_bz, selected_indices)


        return( num_part )

    def transform_particles_to_lab_frame(self):
        """
        Transform the particle quantities from the boosted frame to the
        lab frame. These are classical Lorentz transformation equations
        """
        uzfrm = -self.beta_boost*self.gamma_boost*c
        ic = 1./c
        ic2 = 1./c**2

        # Time in lab frame
        self.t = self.gamma_boost*self.top.time - uzfrm*self.z_captured*ic2
        self.t_prev = self.gamma_boost*(self.top.time - self.top.dt) \
          - uzfrm*self.z_prev_captured*ic2

        # Position in lab frame
        self.z_captured = self.gamma_boost*(self.z_captured + \
            self.beta_boost*c*self.top.time)
        self.z_prev_captured = self.gamma_boost*(self.z_prev_captured \
            + self.beta_boost*c*(self.top.time-self.top.dt))

        # Momentum in lab frame
        self.uz_captured = self.gamma_boost*self.uz_captured \
        - self.gamma_captured*uzfrm
        self.uz_prev_captured = self.gamma_boost*self.uz_prev_captured \
        - self.gamma_prev_captured*uzfrm

        # Field in lab frame
        if(self.dump_p_fields): 
            cbeta = self.beta_boost*c
            beta_ov_c = self.beta_boost*ic

            temp = np.copy(self.ey_captured)
            self.ey_captured = self.gamma_boost*(self.ey_captured - cbeta*self.bx_captured)
            self.bx_captured = self.gamma_boost*(self.bx_captured - beta_ov_c*temp)

            temp = np.copy(self.ex_captured)
            self.ex_captured = self.gamma_boost*(self.ex_captured + cbeta*self.by_captured)
            self.by_captured = self.gamma_boost*(self.by_captured + beta_ov_c*temp)


            temp = np.copy(self.ey_prev_captured)
            self.ey_prev_captured = self.gamma_boost*(self.ey_prev_captured - cbeta*self.bx_prev_captured)
            self.bx_prev_captured = self.gamma_boost*(self.bx_prev_captured - beta_ov_c*temp)

            temp = np.copy(self.ex_prev_captured)
            self.ex_prev_captured = self.gamma_boost*(self.ex_prev_captured + cbeta*self.by_prev_captured)
            self.by_prev_captured = self.gamma_boost*(self.by_prev_captured + beta_ov_c*temp)
            temp = []



    def interpolate_to_time(self, t_output):
        """
        Interpolate the particle quantities in time to the instant t

        Parameter
        ---------
        t_output: float
            Time **in the lab frame** at which the quantities should be
            interpolated. (Typically, the time of a given LabSnapshot)
        """
        # Calculation interpolation weights, for the time interpolation
        weight_prev = (self.t - t_output)/(self.t - self.t_prev)
        weight_next = (t_output - self.t_prev)/(self.t - self.t_prev)

        # Perform the interpolation
        self.x_captured = \
          self.x_prev_captured * weight_prev + self.x_captured * weight_next
        self.y_captured = \
          self.y_prev_captured * weight_prev + self.y_captured * weight_next
        self.z_captured = \
          self.z_prev_captured * weight_prev + self.z_captured * weight_next
        self.ux_captured = \
          self.ux_prev_captured * weight_prev + self.ux_captured * weight_next
        self.uy_captured = \
          self.uy_prev_captured * weight_prev + self.uy_captured * weight_next
        self.uz_captured = \
          self.uz_prev_captured * weight_prev + self.uz_captured * weight_next
        self.gamma_captured = self.gamma_prev_captured * weight_prev + \
          self.gamma_captured * weight_next
        if(self.dump_p_fields):
            self.ex_captured = \
              self.ex_prev_captured * weight_prev + self.ex_captured * weight_next
            self.ey_captured = \
              self.ey_prev_captured * weight_prev + self.ey_captured * weight_next
            self.ez_captured = \
              self.ez_prev_captured * weight_prev + self.ez_captured * weight_next
            self.bx_captured = \
              self.bx_prev_captured * weight_prev + self.bx_captured * weight_next
            self.by_captured = \
              self.by_prev_captured * weight_prev + self.by_captured * weight_next
            self.bz_captured = \
              self.bz_prev_captured * weight_prev + self.bz_captured * weight_next

    

    def gather_array(self, quantity):
        """
        Gather the quantity arrays and normalize the momenta
        Parameters
        ----------
        quantity: String
            Quantity of the particles that is wished to be gathered

        Returns
        -------
        ar: array of reals
            An array of gathered particle's quantity
        """
        ar = np.zeros(np.shape(self.x_captured)[0])

        if quantity == "x":
            ar = np.array(self.x_captured)
        elif quantity == "y":
            ar = np.array(self.y_captured)
        elif quantity == "z":
            ar = np.array(self.z_captured)
        elif quantity == "ux":
            ar = np.array(self.ux_captured)
        elif quantity == "uy":
            ar = np.array(self.uy_captured)
        elif quantity == "uz":
            ar = np.array(self.uz_captured)
        elif quantity == "w":
            ar = np.array(self.w_captured)
        elif quantity == "gamma":
            ar = np.array(self.gamma_captured)
        elif quantity == "id":
            ar = np.array(self.id_captured)
        elif quantity == "ex":
            ar = np.array(self.ex_captured)
        elif quantity == "ey":
            ar = np.array(self.ey_captured)
        elif quantity == "ez":
            ar = np.array(self.ez_captured)
        elif quantity == "bx":
            ar = np.array(self.bx_captured)
        elif quantity == "by":
            ar = np.array(self.by_captured)
        elif quantity == "bz":
            ar = np.array(self.bz_captured)

        return ar

    def extract_slice(self, species, select, prev_z_boost,
                      current_z_boost, t_output ):
        """
        Extract a slice of the particles that corresponds to a given
        lab snapshot, and transform them to the lab frame

        If select is present, extract only particles that satisfy the criteria

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which one slice will be extracted

        select: dict
            A set of rules defined by the users in selecting the particles
            Ex: {"uz": [50, 100]} for particles which have normalized
            values between 50 and 100

        current_z_boost, prev_z_boost: floats
            The position, in the boosted frame, of the plane that
            corresponds to the lab snapshot, at the current iteration and
            at the previous iteration.
            Only particles that have crossed this plane, between the previous
            iteration and the current iteration, are extracted

        t_output: float (s)
            Time **in the lab frame** at which the quantities should be
            obtained. (Typically, the time of a given LabSnapshot.) This
            is used when interpolating between two simulation timesteps.

        Returns
        -------
        slice_array: An array of reals of shape (8, num_part)
            An array that packs together the different particle quantities
            (x, y, z, ux, uy, uz, weight)
        """
        # Declare an attribute for convenience
        p2i = self.particle_to_index

        # Get the particles
        num_part = self.get_particle_slice( species, prev_z_boost,
                        current_z_boost )

        if (hasattr(self.em,"l_pxr")):
            if(self.em.l_pxr): 
                if(self.dump_p_fields):
                    self.em.lorentz_transform_parts_with_fields(num_part, self.gamma_boost, self.beta_boost, self.top.time, self.top.dt,t_output,\
                                     self.x_captured, self.x_prev_captured,\
                                     self.y_captured, self.y_prev_captured,\
                                     self.z_captured, self.z_prev_captured,\
                                     self.ux_captured, self.ux_prev_captured,\
                                     self.uy_captured, self.uy_prev_captured,\
                                     self.uz_captured, self.uz_prev_captured,\
                                     self.gamma_captured, self.gamma_prev_captured,\
                                     self.ex_captured, self.ex_prev_captured,\
                                     self.ey_captured, self.ey_prev_captured,\
                                     self.ez_captured, self.ez_prev_captured,\
                                     self.bx_captured, self.bx_prev_captured,\
                                     self.by_captured, self.by_prev_captured,\
                                     self.bz_captured, self.bz_prev_captured)

                else: 
                    self.em.lorentz_transform_parts_without_fields(num_part, self.gamma_boost, self.beta_boost, self.top.time, self.top.dt,t_output,\
                                     self.x_captured, self.x_prev_captured,\
                                     self.y_captured, self.y_prev_captured,\
                                     self.z_captured, self.z_prev_captured,\
                                     self.ux_captured, self.ux_prev_captured,\
                                     self.uy_captured, self.uy_prev_captured,\
                                     self.uz_captured, self.uz_prev_captured,\
                                     self.gamma_captured, self.gamma_prev_captured
                                     )
            else:
                # Transform the particles from boosted frame back to lab frame
                self.transform_particles_to_lab_frame()
                # Interpolate the particle quantities in time, to t_output
                self.interpolate_to_time( t_output )
        else:
            # Transform the particles from boosted frame back to lab frame
            self.transform_particles_to_lab_frame()

            # Interpolate the particle quantities in time, to t_output
            self.interpolate_to_time( t_output )
        slice_array = np.empty((np.shape(p2i.keys())[0], num_part,))

        for quantity in self.particle_to_index.keys():
            # Here typical values for 'quantity' are e.g. 'z', 'ux', 'gamma'
            # you should just gather array locally
            slice_array[ p2i[quantity], ... ] = self.gather_array(quantity)

        # Choose the particles based on the select criteria defined by the
        # users. Notice: this implementation still comes with a cost,
        # one way to optimize it would be to do the selection before Lorentz
        # transformation back to the lab frame
        if (select is not None) and slice_array.size:
            select_array = self.apply_selection(select, slice_array)
            row, column =  np.where(select_array==True)
            temp_slice_array = slice_array[row,column]

            # Temp_slice_array is a 1D numpy array, we reshape it so that it
            # has the same size as slice_array
            slice_array = np.reshape(
                temp_slice_array,(np.shape(p2i.keys())[0],-1))

        # Multiplying momenta by the species mass to make them unitless
        for quantity in ["ux", "uy", "uz"]:
            slice_array[p2i[quantity]] *= species.mass

        return slice_array

    def get_quantity(self, species, quantity, l_prev=False):
        """
        Get a given particle quantity

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which the quantity is extracted

        quantity: string
            Describes which quantity is queried
            Either "x", "y", "z", "ux", "uy", "uz", "w", "ex", "ey", "ez","bx","by","bz" or "id"

        l_prev: boolean
            If True, then return the quantities of the previous timestep;
            else return quantities of the current timestep
        """
        # Extract the chosen quantities

        # At current timestep
        if not(l_prev):
            if quantity == "x":
                quantity_array = species.getx( gather=False )
            elif quantity == "y":
                quantity_array = species.gety( gather=False )
            elif quantity == "z":
                quantity_array = species.getz( gather=False )
            elif quantity == "ux":
                quantity_array = species.getux( gather=False )
            elif quantity == "uy":
                quantity_array = species.getuy( gather=False )
            elif quantity == "uz":
                quantity_array = species.getuz( gather=False )
            elif quantity == "w":
                quantity_array = species.getweights( gather=False )
            elif quantity == "id":
                quantity_array = species.getssn( gather=False )
            elif(quantity == "ex"):
                quantity_array = species.getex(gather = False)
            elif(quantity == "ey"):
                quantity_array = species.getey(gather = False)
            elif(quantity == "ez"):
                quantity_array = species.getez(gather = False)
            elif(quantity == "bx"):
                quantity_array = species.getbx(gather = False)
            elif(quantity == "by"):
                quantity_array = species.getby(gather = False)
            elif(quantity == "bz"):
                quantity_array = species.getbz(gather = False)


        # Or at previous timestep
        else:
            if quantity == "x":
                quantity_array = species.getxold( gather=False )
            elif quantity == "y":
                quantity_array = species.getyold( gather=False )
            elif quantity == "z":
                quantity_array = species.getzold( gather=False )
            elif quantity == "ux":
                quantity_array = species.getuxold( gather=False )
            elif quantity == "uy":
                quantity_array = species.getuyold( gather=False )
            elif quantity == "uz":
                quantity_array = species.getuzold( gather=False )
            elif(quantity == "ex"): 
                quantity_array = species.getexold(gather = False)
            elif(quantity == "ey"):
                quantity_array = species.geteyold(gather = False)
            elif(quantity == "ez"):
                quantity_array = species.getezold(gather = False)
            elif(quantity == "bx"):
                quantity_array = species.getbxold(gather = False)
            elif(quantity == "by"):
                quantity_array = species.getbyold(gather = False)
            elif(quantity == "bz"):
                quantity_array = species.getbzold(gather = False)          

        return( quantity_array )

    def allocate_previous_instant(self):
        """
        Allocate the top.'quantity'oldpid arrays. This is used to store
        the previous values of the quantities.
        """
        if not self.top.xoldpid:
            self.top.xoldpid = self.top.nextpid()
        if not self.top.yoldpid:
            self.top.yoldpid = self.top.nextpid()
        if not self.top.zoldpid:
            self.top.zoldpid = self.top.nextpid()
        if not self.top.uxoldpid:
            self.top.uxoldpid = self.top.nextpid()
        if not self.top.uyoldpid:
            self.top.uyoldpid = self.top.nextpid()
        if not self.top.uzoldpid:
            self.top.uzoldpid = self.top.nextpid()
        if(self.dump_p_fields):
            if not self.top.exoldpid:
                self.top.exoldpid = self.top.nextpid() 
            if not self.top.eyoldpid:
                self.top.eyoldpid = self.top.nextpid()
            if not self.top.ezoldpid:
                self.top.ezoldpid = self.top.nextpid()
            if not self.top.bxoldpid:
                self.top.bxoldpid = self.top.nextpid()
            if not self.top.byoldpid:
                self.top.byoldpid = self.top.nextpid()
            if not self.top.bzoldpid:
                self.top.bzoldpid = self.top.nextpid()


    def apply_selection(self, select, slice_array):
        """
        Apply the rules of self.select to determine which
        particles should be written

        Parameters
        ----------
        select: a dictionary that defines all selection rules based
        on the quantities

        Returns
        -------
        A 1d array of the same shape as that particle array
        containing True for the particles that satify all
        the rules of self.select
        """
        p2i = self.particle_to_index

        # Initialize an array filled with True
        select_array = np.ones( np.shape(slice_array), dtype='bool' )

        # Apply the rules successively
        # Go through the quantities on which a rule applies
        for quantity in select.keys():
            # Lower bound
            if select[quantity][0] is not None:
                select_array = np.logical_and(
                    slice_array[p2i[quantity]] >\
                     select[quantity][0], select_array )
            # Upper bound
            if select[quantity][1] is not None:
                select_array = np.logical_and(
                    slice_array[p2i[quantity]] <\
                    select[quantity][1], select_array )

        return select_array
