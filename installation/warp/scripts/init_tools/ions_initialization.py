from warp import Species, Proton, Electron
from warp.particles import tunnel_ionization

def initialize_ion_dict( ion_states, weight, group_elec_by_element=False ):
    """
    Initialize the species for ions (and electron from ionization), for
    different elements and different charge states
    Return two dictionaries: one for ions and one for the electrons

    The tunnel ionization calculation is also automatically put into place
    after initializing the ions.

    Parameters:
    -----------
    ion_states: dict
        A dictionary of the form:
        { 'Hydrogen': {'relative_density':0.75, 'q_start':1, 'q_max':1 },
            'Helium': {'relative_density':0.25, 'q_start':0, 'q_max':2 } }
        that describes how to initialize the ions.
        - 'relative_density' sets the relative weight of the ions compared
           to the reference weight provided
        - 'q_start' is the charge state of the ions when they are initialized
        - 'q_max' is the maximum charge state at which ionization will get
           (for faster calculation, this may be lower than the element's Z)

    weight: double
        A reference weight that correspond to a relative_density of 1.

    group_elec_by_element: bool, optional
        Whether to return all the electrons from a single element together
        or to separate them according to the ion charge state from which
        they come. 

    Returns:
    --------
    Two dictionaries containing a list of species
    (one list per element in the dictionary;
    one species per charge state in each list)
    e.g. for the `ion_states` given above, this returns
    {   'Hydrogen': [ <ion species Hydrogen 1> ],
          'Helium': [ <ion species Helium 0>, <ion species Helium 1>,
                    <ion species Helium 2> ]    },
    {   'Hydrogen': [],
          'Helium': [ <elec species from Helium 0>,
                      <elec species from Helium 1> ]  }
    """
    # Initialize the dictionary that will be returned
    ions = {}
    elec_from_ions = {}

    # Case without an ionization file
    if tunnel_ionization is None:

        # Loop over the elements (here, element is a string)
        for element in ion_states.keys():

            ion_weight = ion_states[ element ][ 'relative_density' ] * weight
            q = ion_states[ element ][ 'q_start' ]
            ion = Species( type=Proton, weight=ion_weight,
                           charge_state=q, name=element+str(q)+'+' )

            # Add the ion list for the current element to the ions dictionary
            ions[ element ] = [ ion ]
            elec_from_ions[ element ] = []

    # Case with an ionization file
    if tunnel_ionization is not None:
                    
        # Initialize the ionizer object, which keeps track of which
        # species ionizes into which other species
        tunnel_ionizer = tunnel_ionization.TunnelIonization(stride=1)

        # Loop over the elements (here, element is a string)
        for element in ion_states.keys():

            # Create the ion and electron list
            ion_list = []
            elec_list = []
            # Weight is that of the electrons, multiplied by relative density
            ion_weight = ion_states[ element ][ 'relative_density' ] * weight
            # Extract the element object from tunnel_ionization
            ion_element = getattr( tunnel_ionization, element )
            # Extract the min and max of the charge states
            q_max = min( ion_states[element]['q_max'], ion_element.Z )
            q_start = ion_states[element]['q_start']

            # Loop over the charge states and create the ions
            for q in range(q_start, q_max+1):
                ion = Species( type=ion_element, weight=ion_weight,
                               charge_state=q, name=element+str(q)+'+' )
                ion_list.append( ion )
            # Loop over the charge states and create the electrons
            if group_elec_by_element:
                elec_list = [ Species( type=Electron, weight=ion_weight,
                                name='electron from '+element ) ]
            else:
                for q in range(q_start, q_max):
                    elec = Species( type=Electron, weight=ion_weight,
                                name='electron from '+element+str(q)+'+' )
                    elec_list.append( elec )
            # Introduce tunnel ionization between the newly created species
            for q in range(q_start, q_max):
                if group_elec_by_element:
                    q_elec = 0
                else:
                    q_elec = q-q_start
                tunnel_ionizer.add( incident_species=ion_list[q-q_start],
                    emitted_species=[ion_list[q+1-q_start], elec_list[q_elec]] )

            # Add the ion list for the current element to the ions dictionary
            ions[ element ] = ion_list
            elec_from_ions[ element ] = elec_list

    # Return the dictionary object
    return( ions, elec_from_ions )
