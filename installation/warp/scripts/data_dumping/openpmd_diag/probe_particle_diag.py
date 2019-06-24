"""
This file defines the class ProbeParticleDiagnostic

This diagnostic records the particles that go through a given plane
in time and saves them in a single openPMD file (while doing caching, in
order not to access the disk at every single timestep)

Note: this is a hack the openPMD standard in that it stores particles from
different time in the same openPMD file. (The particles have an addiditional
attribute `t` which is stored in the openPMD file.)
"""
import os
import numpy as np
import time
from scipy.constants import c
from particle_diag import ParticleDiagnostic
from warp_parallel import gatherarray, mpiallgather
from data_dict import particle_quantity_dict

class ParticleAccumulator(ParticleDiagnostic):
    """
    Class that allows buffering of particle quantities

    Usage
    -----
    After initialization, the diagnostic is called by using
    the 'write' method.
    """
    def __init__(self,period_flush,period_diag, top, w3d, comm_world=None,
                 particle_data=["position", "momentum", "weighting", "t"],
                 select=None, write_dir=None, lparallel_output=False,
                 write_metadata_parallel=False,
                 species={"electrons": None},iteration_min=None,iteration_max=None,
                 onefile_per_flush=False):
        """
        Initialization

        Parameters
        ----------

        period_flush: int
            Number of iterations for which the data is accumulated in memory,
            before finally writing it to the disk.

        period_diag: int
            Period at which the diagnostic looks for new particle
            to be accumulated in memory.

        onefile_per_flush: boolean
            if False (default), produces one file for the entire run.
            if True, produces one file per flush (useful for very large dumps
            -e.g in 3D-, where resizing large datasets can be really costly)

        See the documentation of ParticleDiagnostic for the other parameters
        """
        # Do not leave write_dir as None, as this may conflict with
        # the default directory ('./diags')
        if write_dir is None:
            write_dir = 'probe_diags'

        # Initialize Particle diagnostic normal attributes
        ParticleDiagnostic.__init__(self, period_flush, top, w3d,
            comm_world=comm_world,
            species=species, particle_data=particle_data, select=select,
            write_dir=write_dir, lparallel_output=lparallel_output,
            write_metadata_parallel=write_metadata_parallel,
            iteration_min=iteration_min,iteration_max=iteration_max)
        self.period_diag = period_diag
        self.onefile_per_flush=onefile_per_flush

        # Initialize proper helper objects
        self.particle_storer = ParticleStorer( top.dt, self.write_dir,
            self.species_dict, self.lparallel_output, self.rank )
        # Sanity check
        if ("t" not in self.particle_data):
            particle_data.append("t")
        # Init particle catcher object
        self.init_catcher_object()

        # Initialize a corresponding empty file
        if (not self.onefile_per_flush) and \
            (self.write_metadata_parallel or self.rank == 0):
            self.create_file_empty_particles(
                self.particle_storer.filename, 0, 0, self.top.dt )

    def init_catcher_object (self):
        self.particle_catcher = ParticleCatcher( self.top, self.particle_data)


    def write( self ):
        """
        Redefines the method write of the parent class ParticleDiagnostic

        Should be registered with installafterstep in Warp
        """
        if ((self.top.it>=self.iteration_min) and \
            (self.top.it<=self.iteration_max)):
            # At each period_diag, store new particles in memory buffers
            if self.top.it % self.period_diag == 0:
                self.store_new_particles()
            # Every self.period, write the buffered slices to disk
            if (self.top.it % self.period == 0 or \
                self.top.it==self.iteration_max):
                self.flush_to_disk()

    def store_new_particles( self ):
        """
        Store new of the particles in the memory buffers of the
        particle storer

        The stored particles are also selected in accordance with
        the selection rules provided as argument to ProbeParticleDiagnostic
        """
        # Loop through the particle species and register the
        # particle arrays in the particle storer object (buffering)
        for species_name, species in self.species_dict.iteritems():

            slice_array = self.particle_catcher.extract_slice(
                        species, self.select )
            self.particle_storer.register_slice( slice_array, species_name )

    def flush_to_disk(self):
        """
        Writes the buffered slices of particles to the disk. Erase the
        buffered slices of the ParticleStorer object

        """
        # Prepare dictionary that contain, for each species, the list of the
        # local number of macroparticles to be dumped on each proc,
        # the total number of particles to be dumped across all procs,
        # and the compact 2d arrays of particle quantities (with shape
        # (n_quantity, n_particle))
        nlocals_dict = dict()
        nglobal_dict = dict()
        parray_dict  = dict()

        # Compact the successive slices that have been buffered
        # over time into a single array
        for species_name in self.species_dict:
            particle_array = self.particle_storer.compact_slices(species_name)

            if self.comm_world is not None:
                if (self.lparallel_output):
                    # Prepare parallel HDF5 output
                    parray_dict[species_name]=particle_array
                    n = np.size(particle_array[0])
                    nlocals_dict[species_name]= mpiallgather( n )
                    nglobal_dict[species_name]=np.sum(nlocals_dict[species_name])

                else:
                    # Prepare HDF5 output by the first proc, using MPI gathering
                    nlocals_dict[species_name]= None
                    n_rank = self.comm_world.allgather(np.shape(particle_array)[1])

                    # Note that gatherarray routine in parallel.py only works
                    # with 1D array. Here we flatten the 2D particle arrays
                    # before gathering.
                    g_curr = gatherarray(particle_array.flatten(),
                        root=0, comm=self.comm_world )

                    if self.rank == 0:
                        # Get the number of quantities
                        nquant = np.shape(self.particle_catcher.particle_to_index.keys())[0]

                        # Prepare an empty array for reshaping purposes. The
                        # final shape of the array is (8, total_num_particles)
                        parray_dict[species_name]= np.empty((nquant, 0))

                        # Index needed in reshaping process
                        n_ind = 0

                        # Loop over all the processors, if the processor
                        # contains particles, we reshape the gathered_array
                        # and reconstruct by concatenation
                        for i in xrange(self.top.nprocs):

                            if n_rank[i] != 0:
                                parray_dict[species_name] = \
                                np.concatenate((parray_dict[species_name], np.reshape( \
                                    g_curr[n_ind:n_ind+nquant*n_rank[i]], \
                                    (nquant,n_rank[i]))),axis=1)

                                # Update the index
                                n_ind += nquant*n_rank[i]
                    else:
                        parray_dict[species_name] = particle_array
                    # Get global size on all procs
                    n = np.sum(n_rank)
                    nglobal_dict[species_name]= n

            else:
                # Prepare single-proc output (for single-proc simulation)
                parray_dict[species_name] = particle_array
                n = np.size(parray_dict[species_name][0])
                nlocals_dict[species_name]= None
                nglobal_dict[species_name]= n

        if self.onefile_per_flush:
            # Create the file for current flush
            iteration = self.top.it
            file_suffix = "data%08d.h5" %iteration
            curr_filename = os.path.join( self.write_dir, "hdf5", file_suffix  )
            self.create_file_empty_particles( curr_filename, iteration, \
                     self.top.time, self.top.dt, select_nglobal_dict=nglobal_dict )
        else:
            iteration = self.particle_storer.iteration
            # File already created (same file for all flushes)
            curr_filename = self.particle_storer.filename

        # Open the file with or without parallel I/O depending on self.lparallel_output
        f = self.open_file( curr_filename, parallel_open=self.lparallel_output)

        for species_name in self.species_dict:
            species_path = "/data/%d/particles/%s" %(iteration,species_name)
            if f is not None:
                species_grp = f[species_path]
            else:
                species_grp = None
            # Write this array to disk (if this self.particle_storer has new slices)
            self.write_slices(species_grp, parray_dict[species_name], \
            self.particle_catcher.particle_to_index, nlocals_dict[species_name],nglobal_dict[species_name])
            # Erase the buffers
            self.particle_storer.buffered_slices[species_name] = []

        # Close the file
        if f is not None:
            f.close()

    def write_probe_dataset(self, species_grp, path, data, quantity, n_rank, nglobal):
        """
        Writes each quantity of the buffered dataset to the disk, the
        final step of the writing
        """
        if (species_grp is not None) and (nglobal>0) :
            dset = species_grp[path]
            # Resize the h5py dataset if one file for entire run
            if not self.onefile_per_flush:
                index = dset.shape[0]
                dset.resize(index+nglobal, axis=0)
            else:
                index=0
            # All procs write the data
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

    def write_slices( self, species_grp, particle_array, p2i, n_locals, nglobal ):
        """
        Write the slices of the different species to an openPMD file

        Parameters
        ----------
        species_grp: an h5py.Group
            Represent the group in which to write the current species

        particle_array: array of reals
            Array of shape (8, num_part)

        p2i: dict
            Dictionary of correspondance between the particle quantities
            and the integer index in the particle_array

        n_locals: list or None:
            A list with one element per MPI rank, containing the number of
            particles to be dumped from each rank.
            Necessary for parallel HDF5 output.
            If None: indicates that the output is serial (either single-proc
            simulation, or output by first proc using MPI gather)

        nglobal: int
            The total number of particles to be dumped, across all procs.
        """

        # Loop over the different quantities that should be written
        for particle_var in self.particle_data:

            if particle_var in  ["position","momentum","E", "B"]:
                for coord in ["x","y","z"]:
                    quantity= "%s%s" %(particle_quantity_dict[particle_var],coord)
                    path = "%s/%s" %(particle_var, coord)
                    data = particle_array[ p2i[ quantity ] ]
                    self.write_probe_dataset(
                            species_grp, path, data, quantity, n_locals, nglobal)

            elif particle_var == "t":
               quantity= "t"
               path = "t"
               data = particle_array[ p2i[ quantity ] ]
               self.write_probe_dataset(species_grp, path, data, quantity, n_locals, \
               nglobal)

            elif particle_var == "weighting":
               quantity= "w"
               path = "weighting"
               data = particle_array[ p2i[ quantity ] ]
               self.write_probe_dataset(species_grp, path, data, quantity, n_locals, \
               nglobal )

            elif particle_var == "id":
               quantity= "id"
               path = "id"
               data = (np.rint(particle_array[ p2i[ quantity ] ])).astype('uint64')
               self.write_probe_dataset(species_grp, path, data, quantity, n_locals, \
               nglobal )


class ProbeParticleDiagnostic(ParticleAccumulator):
    """
    Class that writes the particles that go across a given plane, in
    the direction given by `plane_normal_vector`
    (The particles that cross the plane in the other direction are not saved.)

    Usage
    -----
    After initialization, the diagnostic is called by using
    the 'write' method.
    """
    def __init__(self, plane_position, plane_normal_vector,
                 period, top, w3d, comm_world=None,
                 particle_data=["position", "momentum", "weighting", "t"],
                 select=None, write_dir=None, lparallel_output=False,
                 write_metadata_parallel=False, onefile_per_flush=False,
                 species={"electrons": None},iteration_min=None,
                 iteration_max=None, plane_velocity=None):
        """
        Initialize diagnostics that retrieve the particles crossing a given
        plane.

        Parameters
        ----------
        plane_position: a 1darray containing 3 floats (in meters)
            The position (in x, y, z) of one of the points of the plane

        plane_normal_vector: a 1darray containing 3 floats
            The coordinates (in x, y, z) of one of the vectors of the plane

        period: int
            Number of iterations for which the data is accumulated in memory,
            before finally writing it to the disk.

        plane_velocity : a 1darray containing 3 floats
            Speed of the plane (in x ,y, z) in m/s

        See the documentation of ParticleDiagnostic and ParticleAccumulator
        for the other parameters
        """
        # Do not leave write_dir as None, as this may conflict with
        # the default directory ('./diags')
        if write_dir is None:
            write_dir = 'probe_diags'
        self.plane_position = plane_position
        self.plane_normal_vector = plane_normal_vector
        self.plane_velocity = plane_velocity

        # Initialize Particle Accumulator normal attributes
        ParticleAccumulator.__init__(self, period, 1, top, w3d, comm_world,
            species=species, particle_data=particle_data, select=select,
            write_dir=write_dir, lparallel_output=lparallel_output,
            write_metadata_parallel=write_metadata_parallel,
            onefile_per_flush=onefile_per_flush,
            iteration_min=iteration_min,iteration_max=iteration_max)


        # Initialize proper helper objects
        self.particle_storer = ParticleStorer( top.dt, self.write_dir,
            self.species_dict, self.lparallel_output, self.rank )
        self.init_catcher_object()
        self.particle_catcher.allocate_previous_instant()


    def init_catcher_object (self):
        self.particle_catcher = ParticleProbeCatcher( self.top, self.plane_position, \
                                self.plane_normal_vector,  self.plane_velocity,
                                self.particle_data )

    def write( self ):
        """
        Redefines the method write of the parent class ParticleDiagnostic

        Should be registered with installafterstep in Warp
        """
        if self.plane_velocity is not None:
            increment = self.plane_velocity * self.top.dt
            self.plane_position += increment
            self.particle_catcher.plane_position +=  increment

        ParticleAccumulator.write( self )

class ParticleStorer:
    """
    Class that stores data relative to the particles that are crossing the plane
    """
    def __init__(self, dt, write_dir, species_dict, lparallel_output, rank):
        """
        Initialize a ParticleStorer object

        Parameters
        ----------
        write_dir: string
            Absolute path to the directory where the data for
            this snapshot is to be written

        species_dict: dict
            Contains all the species name of the species object
            (inherited from Warp)
        """
        # Deduce the name of the filename where this snapshot writes
        self.filename = os.path.join( write_dir, 'hdf5/data%08d.h5' %0)
        self.iteration = 0
        self.dt = dt

        # Prepare buffered slices
        self.buffered_slices = {}
        for species_name in species_dict:
            self.buffered_slices[species_name] = []

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
            particle_array = np.empty((8,0))

        return particle_array

class ParticleCatcher:
    """
    Class that extracts, and gathers particle quantities
    Provides tools to select and extract particle quantities according
    to selection rules. Selection rules can be customized by defining
    a derived class of ParticleCatcher and defining a new method
    get_particle_slice()
    """
    def __init__(self, top, particle_data):
        """
        Initialize the ParticleCatcher object

        Parameters
        ----------
        particle_data: list of particle data to "catch"

        top: WARP object
        """

        # Get list of particle quantities to catch
        # for current species
        list_of_quantities=[]
        for particle_var in particle_data:
            if particle_var in ["position","momentum","E","B"]:
                for coord in ["x", "y", "z"]:
                    quantity = "%s%s" %(particle_quantity_dict[particle_var],coord)
                    list_of_quantities.append(quantity)
            elif particle_var=="weighting":
                list_of_quantities.append("w")
            elif particle_var=="id":
                list_of_quantities.append("id")
            elif particle_var=="t":
                list_of_quantities.append("t")

        # Some attributes neccessary for particle selections
        self.top = top

        # Create a dictionary that contains the correspondance
        # between the particles quantity and array index
        self.list_of_quantities=list_of_quantities
        nquants=len(self.list_of_quantities)
        self.nquants=nquants
        particle_to_index=dict()
        for i in range(self.nquants):
            particle_to_index[self.list_of_quantities[i]]=i
            if self.list_of_quantities[i]=="t":
                self.t_index=i
        self.particle_to_index = particle_to_index
        self.captured_quantities = dict()

    def get_particle_slice( self, species ):
        """
        Select the particles for the current slice, and extract their
        particle quantities. By default, all particles are taken in the
        generic class. User have to redefine this function in a derived
        class to implement custom selection rule adapted to its diag

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which one slice will be extracted

        Returns
        -------
        num_part: int
            Number of selected particles
        """

        # By default all particles are chosen
        for i in range(self.nquants):
            # Quantities at current time step
            if (i is not self.t_index):
                self.captured_quantities[self.list_of_quantities[i]]= \
                    self.get_quantity( species,self.list_of_quantities[i] )

        i_not_t = np.delete(np.arange(self.nquants),self.t_index)[0]
        num_part= \
        np.size(self.captured_quantities[self.list_of_quantities[i_not_t]])
        self.captured_quantities[self.list_of_quantities[self.t_index]]= \
            np.ones(num_part)*self.top.time

        return( num_part )

    def gather_array(self, quantity):
        """
        Get quantity arrays to be gathered
        ----------
        quantity: String
            Quantity of the particles that is wished to be gathered

        Returns
        -------
        ar: array of reals
            An array of gathered particle's quantity
        """

        return self.captured_quantities[quantity]

    def extract_slice(self, species, select ):
        """
        Extract a slice of the particles

        If select is present, extract only particles that satisfy the criteria

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which one slice will be extracted

        select: dict
            A set of rules defined by the users in selecting the particles
            Ex: {"uz": [50, 100]} for particles which have normalized
            values between 50 and 100

        Returns
        -------
        slice_array: An array of reals of shape (8, num_part)
            An array that packs together the different particle quantities
            (x, y, z, ux, uy, uz, weight, t)
        """
        # Declare an attribute for convenience
        p2i = self.particle_to_index

        # Get the particles
        num_part = self.get_particle_slice( species )
        slice_array = np.empty((np.shape(p2i.keys())[0], num_part,))

        # Get the particle quantities
        for quantity in self.particle_to_index.keys():
            # Here typical values for 'quantity' are e.g. 'z', 'ux', 'gamma'
            # you should just gather array locally
            slice_array[ p2i[quantity], ... ] = self.gather_array(quantity)

        # Choose the particles based on the select criteria defined by the
        # users.
        if (select is not None) and slice_array.size:
            select_array = self.apply_selection(select, slice_array)
            row, column =  np.where(select_array==True)
            temp_slice_array = slice_array[row,column]

            # Temp_slice_array is a 1D numpy array, we reshape it so that it
            # has the same size as slice_array
            slice_array = np.reshape(
                temp_slice_array,(np.shape(p2i.keys())[0],-1))

        # Multiplying momenta by the species mass to make them unitless
        for quantity in self.particle_to_index.keys():
             if quantity in ["ux", "uy", "uz"]:
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
            Either "x", "y", "z", "ux", "uy", "uz", "w"

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

         # Quantities that do not depend on time step
        if quantity == "w":
             quantity_array = species.getweights( gather=False )
        elif quantity == "id":
             quantity_array = \
             (np.rint( species.getssn(gather=False) )).astype('uint64')
        elif quantity == "ex":
            quantity_array = species.getex( gather=False )
        elif quantity == "ey":
            quantity_array = species.getey( gather=False )
        elif quantity == "ez":
            quantity_array = species.getez( gather=False )
        elif quantity == "bx":
            quantity_array = species.getbx( gather=False )
        elif quantity == "by":
            quantity_array = species.getby( gather=False )
        elif quantity == "bz":
            quantity_array = species.getbz( gather=False )

        return( quantity_array )

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

class ParticleProbeCatcher(ParticleCatcher):
    """
    Class that extracts, interpolates and gathers particles
    """
    def __init__(self, top, plane_position, plane_normal_vector, plane_velocity,
                  particle_data ):
        """
        Initialize the ParticleProbeCatcher object

        Parameters
        ----------
        plane_position: a list of 3 floats (in meters)
            The position (in x, y, z) of one of the points of the plane

        plane_normal_vector: a list of 3 floats
            The coordinates (in x, y, z) of one of the vectors of the plane

        top: WARP object
        """

        # Init ParticleCatcher Normal attributes
        ParticleCatcher.__init__(self, top, particle_data)

        # Some attributes neccessary for particle selections
        self.plane_position = plane_position
        self.plane_normal_vector = plane_normal_vector
        self.plane_velocity  = plane_velocity

    def get_particle_slice( self, species ):
        """
        Select the particles for the current slice, and extract their
        positions and momenta at the current and previous timestep

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which one slice will be extracted

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

        # Quantities at previous time step
        previous_x = self.get_quantity( species, "x", l_prev=True )
        previous_y = self.get_quantity( species, "y", l_prev=True )
        previous_z = self.get_quantity( species, "z", l_prev=True )
        previous_ux = self.get_quantity( species, "ux", l_prev=True )
        previous_uy = self.get_quantity( species, "uy", l_prev=True )
        previous_uz = self.get_quantity( species, "uz", l_prev=True )

        # Quantities non related to the time
        if self.top.wpid:
            weights = self.get_quantity( species, "w" )
        if self.top.ssnpid:
            pid = self.get_quantity( species, "id" )

        if 'bx' in self.list_of_quantities:
            bx = self.get_quantity(species, 'bx')
        if 'by' in self.list_of_quantities:
            by = self.get_quantity(species, 'by')
        if 'bz' in self.list_of_quantities:
            bz = self.get_quantity(species, 'bz')
        if 'ex' in self.list_of_quantities:
            ex = self.get_quantity(species, 'ex')
        if 'ey' in self.list_of_quantities:
            ey = self.get_quantity(species, 'ey')
        if 'ez' in self.list_of_quantities:
            ez = self.get_quantity(species, 'ez')

        # This part is then common for WARP or WARP + PICSAR
        # A particle array for mapping purposes
        particle_indices = np.arange( len(current_z) )

        # For this snapshot:
        # - check if the particles where before the plane at the previous timestep
        # - check if the particle are beyond the plane at the current timestep
        if self.plane_velocity is not None:
            r_old = self.plane_position - self.plane_velocity * self.top.dt
        else:
            r_old = self.plane_position
        r = self.plane_position
        n = self.plane_normal_vector
        previous_position_relative_to_plane = \
              n[0]*(previous_x - r_old[0]) \
            + n[1]*(previous_y - r_old[1]) \
            + n[2]*(previous_z - r_old[2])
        current_position_relative_to_plane = \
              n[0]*(current_x - r[0]) \
            + n[1]*(current_y - r[1]) \
            + n[2]*(current_z - r[2])
        selected_indices = np.compress(
            (previous_position_relative_to_plane <= 0 ) &
            (current_position_relative_to_plane > 0 )  , particle_indices)

        num_part = np.shape(selected_indices)[0]

        self.mass = species.mass
        if self.top.wpid:
            self.captured_quantities['w'] = np.take(weights, selected_indices)
        if self.top.ssnpid:
            self.captured_quantities['id'] = np.take( pid, selected_indices)
        if 'bx' in self.list_of_quantities:
            self.captured_quantities['bx'] = np.take( bx, selected_indices)
        if 'by' in self.list_of_quantities:
            self.captured_quantities['by'] = np.take( by, selected_indices)
        if 'bz' in self.list_of_quantities:
            self.captured_quantities['bz'] = np.take( bz, selected_indices)
        if 'ex' in self.list_of_quantities:
            self.captured_quantities['ex'] = np.take( ex, selected_indices)
        if 'ey' in self.list_of_quantities:
            self.captured_quantities['ey'] = np.take( ey, selected_indices)
        if 'ez' in self.list_of_quantities:
            self.captured_quantities['ez'] = np.take( ez, selected_indices)

        ## Select the particle quantities that satisfy the
        ## aforementioned condition
        current_x = np.take(current_x, selected_indices)
        current_y = np.take(current_y, selected_indices)
        current_z = np.take(current_z, selected_indices)
        current_ux = np.take(current_ux, selected_indices)
        current_uy = np.take(current_uy, selected_indices)
        current_uz = np.take(current_uz, selected_indices)
        current_position_relative_to_plane = np.take(
            current_position_relative_to_plane, selected_indices )

        previous_x = np.take(previous_x, selected_indices)
        previous_y = np.take(previous_y, selected_indices)
        previous_z = np.take(previous_z, selected_indices)
        previous_ux = np.take(previous_ux, selected_indices)
        previous_uy = np.take(previous_uy, selected_indices)
        previous_uz = np.take(previous_uz, selected_indices)
        previous_position_relative_to_plane = np.take(
            previous_position_relative_to_plane, selected_indices )

        # Interpolate particle quantity to the time when they cross the plane
        norm_factor = 1 / ( np.abs(previous_position_relative_to_plane) \
                + current_position_relative_to_plane )
        interp_current = np.abs(previous_position_relative_to_plane) * norm_factor
        interp_previous = current_position_relative_to_plane * norm_factor

        self.captured_quantities['t']= interp_current * self.top.time + \
                            interp_previous * (self.top.time - self.top.dt)
        self.captured_quantities['x'] = interp_current * current_x + \
                            interp_previous * previous_x
        self.captured_quantities['y'] = interp_current * current_y + \
                            interp_previous * previous_y
        self.captured_quantities['z'] = interp_current * current_z + \
                            interp_previous * previous_z
        self.captured_quantities['ux'] = interp_current * current_ux + \
                            interp_previous * previous_ux
        self.captured_quantities['uy'] = interp_current * current_uy + \
                            interp_previous * previous_uy
        self.captured_quantities['uz'] = interp_current * current_uz + \
                            interp_previous * previous_uz

        return( num_part )
