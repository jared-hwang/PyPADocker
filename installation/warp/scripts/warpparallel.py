"""Functions used in parallel version.
Most important ones are the paralleldump and parallelrestore functions.
"""
from warp import *
import __main__
import copy


def warpparalleldoc():
    import warpparallel
    print warpparallel.__doc__

try:
    #pyMPI version
    import mpi
    top.comm_world = comm_world.comm_fortran() #pyMPI
    top.lcomm_world_initted = true

except:
    try:
        #mpi4py version
        from mpi4py import MPI as mpi
        top.comm_world = comm_world.py2f() #mpi4py
        top.lcomm_world_initted = true
    except:
        pass

top.my_index = me
top.nprocs = npes
top.nslaves = top.nprocs

# ---------------------------------------------------------------------------
def gatherallzarray(a,zaxis=0):
    """Gathers and broadcasts the data in a z-array which is decomposed
  the same way as the particle domains. Each processor contributes the
  data from within the particle decomposition region it owns. This works
  with any array from the groups Z_arrays and Z_Moments.
   - first argument is the z-array
   - zaxis: axis which is decomposed in z
    """
    if not lparallel: return a
    # --- Get start and end of particle decomposition region
    iz1 = 0
    if me < npes-1: iz2 = top.izpslave[me+1] - 1 - top.izpslave[me]
    else:           iz2 = w3d.nz - top.izpslave[me]
    # --- Rearrange array to put the decomposed axis first
    if zaxis != 0: a = swapaxes(a,0,zaxis)
    # --- Gather and broadcast it
    result = gatherarray(a[iz1:iz2+1,...],bcast=1)
    # --- Rearrange array to put the decomposed axis back where it started
    if zaxis != 0: result = swapaxes(result,0,zaxis)
    return result

# ---------------------------------------------------------------------------
def scatterallzarray(a,zaxis=0):
    """Scatters the data in a z-array which is decomposed the same way as
  the particle domains. Each processor contributes the data from within
  the particle decomposition region it owns. This works with any array
  from the groups Z_arrays and Z_Moments.
   - first argument is the z-array
   - zaxis: axis which is decomposed in z
    """
    if not lparallel: return a
    # --- Rearrange array to put the decomposed axis first
    if zaxis != 0: a = swapaxes(a,0,zaxis)
    # --- Get the appropriate subsection
    result = a[top.izpslave[me]:top.izpslave[me]+top.nzpslave[me] + 1,...]
    # --- Rearrange array to put the decomposed axis back where it started
    if zaxis != 0: result = swapaxes(result,0,zaxis)
    return result

# ---------------------------------------------------------------------------
def gatherallzfsarray(a,zaxis=0):
    """Gathers and broadcasts the data in a z-array decomposed in the same
  way as the field grid. Each processor contributes the data from within
  the field-solve decomposition region it owns.
   - first argument is the z-array
   - zaxis: axis which is decomposed in z
    """
    if not lparallel: return a
    # --- Get start and end of field-solve decomposition region
    iz1 = 0
    if me < npes-1: iz2 = top.izfsslave[me+1] - 1 - top.izfsslave[me]
    else:           iz2 = w3d.nz - top.izfsslave[me]
    # --- Rearrange array to put the decomposed axis first
    if zaxis != 0: a = swapaxes(a,0,zaxis)
    # --- Gather and broadcast it
    result = gatherarray(a[iz1:iz2+1,...],bcast=1)
    # --- Rearrange array to put the decomposed axis back where it started
    if zaxis != 0: result = swapaxes(result,0,zaxis)
    return result

# ---------------------------------------------------------------------------
def scatterallzfsarray(a,zaxis=0):
    """Scatters the data in a z-array decomposed in the same way as the
  field grid. Each processor contributes the data from within the
  field-solve decomposition region it owns.
   - first argument is the z-array
   - zaxis: axis which is decomposed in z
    """
    if not lparallel: return a
    # --- Rearrange array to put the decomposed axis first
    if zaxis != 0: a = swapaxes(a,0,zaxis)
    # --- Get the appropriate subsection
    result = a[top.izfsslave[me]:top.izfsslave[me]+top.nzfsslave[me] + 1,...]
    # --- Rearrange array to put the decomposed axis back where it started
    if zaxis != 0: result = swapaxes(result,0,zaxis)
    return result

#-------------------------------------------------------------------------
def convertiztope(iz):
    """Given an iz value, returns the processor number whose particle region
  contains that value."""
    if 0 <= iz <= w3d.nz:
        # --- This finds all of the processors for which have iz within their
        # --- domain. The last one is selected since in the regions which
        # --- overlap, the standard is that the processor to the right has
        # --- priority for that region, i.e. the processor which has the
        # --- overlapping on it left hand edge.
        pe = compress(logical_and(less_equal(top.izpslave,iz),
                        less_equal(iz,top.izpslave+top.nzpslave)),arange(npes))[-1]
    else:
        pe = None
    return pe
convertizptope = convertiztope

def convertizfstope(iz):
    """Given an iz value, returns the processor number whose field solve region
  contains that value."""
    if 0 <= iz <= w3d.nz:
        # --- This finds all of the processors for which have iz within their
        # --- domain. The last one is selected since in the regions which
        # --- overlap, the standard is that the processor to the right has
        # --- priority for that region, i.e. the processor which has the
        # --- overlapping on it left hand edge.
        pe = compress(logical_and(less_equal(top.izfsslave,iz),
                      less_equal(iz,top.izfsslave+top.nzfsslave)),arange(npes))[-1]
    else:
        pe = None
    return pe

#-------------------------------------------------------------------------
def broadcastgroupHist():
    '''Broadcasts the history data, group Hist, from processor 0 to all of the
  other processors. This is needed since normally, the history data is only
  saved on processor 0.'''
    varlist = top.varlist("Hist")
    for vname in varlist:
        x = top.getpyobject(vname)
        x = parallel.broadcast(x)
        if x is not None:
            setattr(top,vname,x)

