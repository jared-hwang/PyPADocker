from warp import installafterstep
from warp.field_solvers.em3dsolver import EM3D
from warp.field_solvers.em3dsolverFFT import EM3DFFT
import numpy as np

def initialize_em_solver( stencil, dim, 
    npass_smooth, alpha_smooth, stride_smooth,
    circ_m=0, norder=None, ntsub=None, nguard_spectral=50,
    centered_spectral=False, **kwargs ):
    """
    Return an instance of an electromagnetic solver, with a
    configuration that corresponds to the provided parameters

    This is a simplified wrapper around EM3D and EM3DFFT. To have access to
    the full list of options, please initialize EM3D and EM3DFFT by hand.

    Parameters
    ----------
    stencil: int
       Describes the type of stencil that is used for the field solver
       0: Yee
       1: Cole-Karkkainen
       -1: Spectral
       For stencil 0 and 1, the default current deposition is Esirkepov
       For stencil -1, the default is direct deposition (rho v) + correction

    npass_smooth: 2darray of integers, of shape (3,2)
        An array which indicates the number of smoother and
        compensator along each direction
        e.g. np.array([ [2,1], [2,1], [2,1] ]) for 2 smoother and 1 compensator
        along each direction

    alpha_smooth: 2darray of floats, of shape (3,2)
        Coefficients along each direction, for the smoother and compensator

    stride_smooth: 2darray of integers, of shape (3,2)
        The stride of the smoother and compensator along each direction

    circ_m: int, optional
        The number of azimuthal modes beyond m=0

    norder: int, optional
        Order of the stencil, which is used for the calculation of
        the spatial derivatives (only for the Yee solver and spectral solver)
        The default is 2 for the Yee solver, and infinity for spectral solver

    ntsub: int, optional
        Number of subcycles per PIC cycle.
        (On each subcycle, the Maxwell equations are advanced in time,
        but the particles are not moved and the current is not redeposited)
        The default is 1 for the finite-difference solvers, and infinity
        (i.e. analytical integration) for the spectral solver.

    nguard_spectral: int, optional
        Number of guard cells around the domain, when stencil==-1
        (When stencil is different than -1, than the number of guard cells
        is determined automatically)

    centered_spectral: bool, optional
        Whether to use a centered scheme for the spectral solver
        
    kwargs: additional optional arguments to be passed to EM3D
        (For instance, dtcoef=0.7)

    Returns
    -------
    An instance of EM3D or EM3DFFT
    """
    # Run a few checks
    if (norder is not None) and ((stencil in [0, -1]) is False):
        raise ValueError("The parameter norder can only be set \n"
            "for the Yee solver (stencil 0) or spectral solver (stencil -1)")

    # Automatically determine a set of geometric parameters
    l_2dxz = (dim in ["2d", "circ"])
    l_2drz = (dim in ["circ"])
    l_1dz = (dim=="1d")
    circ_m = (dim =="circ")*circ_m

    # Configure the solver
    if stencil == -1:

        # Spectral: automatically determine the order and subcycling
        # The default is infinite order
        if norder is None:
            norder = np.inf
        # The default is infinite subcycling
        if ntsub is None:
            ntsub = np.inf

        # Initialize a spectral solver
        em = EM3DFFT(
            # Here the stencil is set only for the purpose of
            # automatically determining dt_courant.
            # Stencil 3 sets dt = w3d.dz/c automatically
            stencil=3, 
            npass_smooth=npass_smooth,
            alpha_smooth=alpha_smooth,
            stride_smooth=stride_smooth,
            l_2dxz=l_2dxz,
            l_1dz=l_1dz,
            l_getrho=True,
            pml_method=2,
            nxguard=nguard_spectral,
            nzguard=nguard_spectral,
            spectral_current=0,
            current_cor=1,
            ntsub=ntsub,
            norderx=norder,
            norderz=norder,
            l_esirkepov=0,
            l_fieldcenterK=True,
            l_nodalgrid=int(centered_spectral),
            l_deposit_nodal=int(centered_spectral),
            **kwargs )

        # Damp the field in the guard cells, in order to prevent
        # them from wrapping around
        def zerofieldsinguardz():
            fields = [em.fields.Ex,em.fields.Ey,em.fields.Ez,
                em.fields.Bx,em.fields.By,em.fields.Bz]
            for field in fields:
                field[:,:,-em.nzguard/2:] = 0
                field[:,:,:em.nzguard/2] = 0
        installafterstep( zerofieldsinguardz )

    else:
        # Initialize a finite-difference solver        
        em = EM3D(
            stencil=stencil,
            npass_smooth=npass_smooth,
            alpha_smooth=alpha_smooth,
            stride_smooth=stride_smooth,
            l_2dxz=l_2dxz,
            l_2drz=l_2drz,
            l_1dz=l_1dz,
            l_getrho=True,
            circ_m=circ_m,
            type_rz_depose=1,
            l_correct_num_Cherenkov=True,
            **kwargs )
    
    # Return the solver
    return( em )
    

