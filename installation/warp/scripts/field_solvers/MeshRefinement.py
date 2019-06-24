"""Implements adaptive mesh refinement in 3d
"""
from __future__ import generators
__all__ = ['MeshRefinement',
           'MRBlock3D','MRBlock','MRBlock2D','MRBlockRZ','MRBlock2DDielectric',
           'MRBlockImplicit2D','EMMRBlock']
from ..warp import *
from find_mgparam import find_mgparam
import types
import collections
try:
    import Opyndx
    VisualizableClass = Opyndx.Visualizable
except ImportError:
    # --- If Opyndx is not available, then use object as the base class,
    # --- disabling any visualization.
    VisualizableClass = object

#import threading

try:
    import psyco
except ImportError:
    pass

#########################################################################
# Note that MRBlock3D is psyco.bind at the end of the file
class MeshRefinement(VisualizableClass):
    """
  Implements adaptive mesh refinement in 3d
   - refinement=None: amount of refinement along each axis
   - lower,upper: extent of domain in relative to parent, in its own grid
                  cell size, and not including any guard cells.
                  Only applies to refinement block.
   - dims: dimensions of the grid, only used for root block, the one with
           no parents.
   - mins,maxs: locations of the grid lower and upper bounds in the beam frame
   - children: list of tuples, each containing three elements,
               (lower,upper,refinement). It is recommended that refinement
               blocks be added later using addchild.
   - lreducedpickle=true: when true, a small pickle is made by removing all of
                          the big arrays. The information can be regenerated
                          upon restart. (Expert use only.)
   - nguard=[1,1,1]: number of coarse guard cells where rho is deposited but
                     the fields not gather. This is needed to reduce the
                     self-field effect of particles near refinement boundaries.
                     The default value should be sufficient.
   - nguarddepos=[0,0,0]: number of coarse guard cells where rho (or current)
                          is not deposited. This is needed for the refinement
                          of the EM solver.
    """
    def __init__(self,parent=None,refinement=None,
                      lower=None,upper=None,
                      fulllower=None,fullupper=None,
                      dims=None,mins=None,maxs=None,
                      nguard=[1,1,1],
                      nguarddepos=[0,0,0],
                      l_EM=0,
                      children=None,lchild=False,**kw):

        self.lchild = lchild

        # --- Save the input dictionary to pass to children
        self.kw = kw

        # --- By default, use electrostatic solver
        self.l_EM=l_EM

        # --- The value of solvergeom needs to be fetched here (normally, is
        # --- would be set in the call to the __init__ of the field solver.
        self.solvergeom = kw.get('solvergeom',w3d.solvergeom)
        if 'solvergeom' in kw: del kw['solvergeom']

        # --- Check the dimensionality. If less than 3D, set the appropriate
        # --- dimension to zero. (The code only works with XYZ and RZ now.)
        if self.solvergeom == w3d.RZgeom and 0:
            if refinement is not None:
                if isinstance(refinement,(ndarray,collections.Sequence)):
                    refinement = [refinement[0],1,refinement[-1]]
                else:
                    refinement = [refinement,1,refinement]
            if lower is not None:
                lower = [lower[0],0,lower[-1]]
            if upper is not None:
                upper = [upper[0],0,upper[-1]]
            if dims is not None:
                dims = [dims[0],0,dims[-1]]
            if mins is not None:
                mins = [mins[0],0.,mins[-1]]
            if maxs is not None:
                maxs = [maxs[0],0.,maxs[-1]]

        if parent is None:

            # --- No parents, so just create empty lists
            self.parents = []
            self.root = self

            # --- It is assumed that the root block will be the first one created.
            # --- So clear out the global block list and count.
            self.totalnumberofblocks = 0
            self.listofblocks = []
            self.finalized = 0
            self.my_index = me
            self.my_indexbase = me

            # --- For the root, the dimensions and extent of the grid should
            # --- be specified. If not, they will be taken from w3d.
            self.dims = dims
            self.mins = mins
            self.maxs = maxs
            self.totalrefinement = ones(3,'l')
            self.forcesymmetries = 1
            self.nguard = nguard
            self.nguarddepos = nguarddepos
            self.deltas = None

            self.refinement = None

            # --- The root block is always active
            self.isactive = True

            # --- top.nparpgrp is now set in the init of the base class.

        else:

            # --- Save the parent and the index number. These are saved in lists
            # --- since a block can have multiple parents.
            self.parents = [parent.blocknumber]
            self.root = parent.root
            self.my_index = 0
            self.my_indexbase = me
            self.nprocs = 1
            self.nxprocs = 1
            self.nyprocs = 1
            self.nzprocs = 1
            self.ixproc = 0
            self.iyproc = 0
            self.izproc = 0

            # --- Make sure that refinement is an array of length three. If a scalar
            # --- is input, it is broadcast to all three axis.
            if len(shape(refinement)) == 0:
                refinement = 3*[refinement]
            self.refinement = array(refinement)
            if len(shape(nguard)) == 0:
                nguard = 3*[nguard]
            self.nguard = array(nguard)
            self.nguarddepos = array(nguarddepos)

            self.totalrefinement = parent.totalrefinement*self.refinement
            self.deltas = parent.deltas/self.refinement
            self.rootdims = self.root.dims*self.totalrefinement
            self.forcesymmetries = 0

            if lower is None and upper is None:
                # --- The grid mins and maxs are input. Convert them to arrays.
                self.mins = array(mins)
                self.maxs = array(maxs)

                # --- Make sure that the input makes sense.
                # --- Make sure that the block has a finite extent in all dimensions
                assert alltrue(self.maxs>=self.mins),\
                     "The child must have a finite extent in all dimensions"

                # --- Save the input values. the actual values used will be modified
                # --- by various constraints.
                self.minsinput = self.mins.copy()
                self.maxsinput = self.maxs.copy()
                self.lowerinput = None
                self.upperinput = None

                # --- The lower and upper are calculated to be an integer number of
                # --- parent grid cells times the refinement factor. The lower is
                # --- rounded down and upper rounded up to ensure that the block
                # --- includes the entire extent specified by mins and maxs.
                self.mins = maximum(self.mins,self.root.mins)
                self.maxs = minimum(self.maxs,self.root.maxs)
                self.lower = (nint(floor((self.mins - self.root.minsglobal)/parent.deltas))*
                             self.refinement)
                self.upper = (nint(ceil((self.maxs - self.root.minsglobal)/parent.deltas))*
                             self.refinement)
                self.lower = maximum(self.root.lower*self.totalrefinement,self.lower)
                self.upper = minimum(self.root.upper*self.totalrefinement,self.upper)

            else:
                # --- The grid lower and upper bounds are input. The bounds are
                # --- relative to the root grid, but scaled by the total refinement.
                self.lower = nint(array(lower))
                self.upper = nint(array(upper))

                # --- Make sure that the input makes sense.
                # --- Make sure that the block has a finite extent in all dimensions
                assert alltrue(self.upper>=self.lower),\
                     "The child must have a finite extent in all dimensions"

                # --- Save the input values. the actual values used will be modified
                # --- by various constraints.
                self.minsinput = None
                self.maxsinput = None
                self.lowerinput = self.lower.copy()
                self.upperinput = self.upper.copy()

            # --- In parallel, an extra grid cell in z can be added since the
            # --- information in the z guard planes of potential is correct.
            self.extradimslower = zeros(3,'l')
            self.extradimsupper = zeros(3,'l')
            if npes > 1 and not self.l_EM:
                if self.root.ixproc > 0:                   self.extradimslower[0] = 1
                if self.root.ixproc < self.root.nxprocs-1: self.extradimsupper[0] = 1
                if self.root.iyproc > 0:                   self.extradimslower[1] = 1
                if self.root.iyproc < self.root.nyprocs-1: self.extradimsupper[1] = 1
                if self.root.izproc > 0:                   self.extradimslower[2] = 1
                if self.root.izproc < self.root.nzprocs-1: self.extradimsupper[2] = 1

            rootfulllower = self.root.fulllower*self.totalrefinement - self.extradimslower*self.totalrefinement
            rootfullupper = self.root.fullupper*self.totalrefinement + self.extradimsupper*self.totalrefinement

            # --- Now, extend the domain by the given number of guard cells. Checks
            # --- are made so that the domain doesn't extend beyond the root grid.
            if fulllower is None:
                self.fulllower = maximum(rootfulllower,self.lower - self.nguard*self.refinement)
            else:
                self.fulllower = array(fulllower)
            if fullupper is None:
                self.fullupper = minimum(rootfullupper,self.upper + self.nguard*self.refinement)
            else:
                self.fullupper = array(fullupper)

            # --- Get the number of grid points along each dimension
            self.dims = self.fullupper - self.fulllower

            # --- Make sure that the number of grid points is even.
            # --- If it is odd, then enough cells are added to extend to the next
            # --- grid cell of the parent. It is then cutoff at the root grid.
            self.fulllower = where(self.dims%2==1,self.fulllower-self.refinement,
                                                  self.fulllower)
            self.fulllower = maximum(rootfulllower,self.fulllower)
            self.dims = self.fullupper - self.fulllower

            # --- If it is still odd (which means that the cells added above
            # --- where cutoff at zero) then add some at the top.
            self.fullupper = where(self.dims%2==1,self.fullupper+self.refinement,
                                                  self.fullupper)
            self.fullupper = minimum(rootfullupper,self.fullupper)
            self.dims = self.fullupper - self.fulllower

            # --- If it is still odd, then there is some serious problem. The number
            # --- in the base grid may be odd.
            assert alltrue(self.dims[:2]%2 == 0),\
                   """The number of grid cells in one of the transverse dimensions is odd - they all must be even. Check that the number of cells in the base grid is even."""

            # --- Now calculate the extent of the grid
            self.mins = self.root.minsglobal + self.fulllower*self.deltas
            self.maxs = self.root.minsglobal + self.fullupper*self.deltas

            # --- Recalculate the deltas so that it is done is a manner that is
            # --- consistent with the calculation elsewhere. This is needed to
            # --- avoid problems with round off. This ensures that the values
            # --- of the deltas calculated here and elsewhere will be the same
            # --- out to machine precision.
            self.deltas = where(self.dims>0.,
                                (self.maxs - self.mins)/
                                    where(self.dims>0.,self.dims,1),
                                self.deltas)

            # --- Note that it is not needed to check if the child overlaps the
            # --- parent. This allows a child to be added to any parent and avoids
            # --- problems in the parallel version where a child may not intersect
            # --- the root block on some processors.

            # --- Check if the block is active, i.e. has a finite extent.
            # --- This should only be false in parallel, where blocks will not
            # --- necessarily intersect to domain of some processors.
            # --- This flag can also provide a way for a user to turn off blocks
            # --- for testing purposes or otherwise.
            self.isactive = alltrue(self.upper>=self.lower)

            # --- First, just use same boundary conditions as root.
            self.bounds = self.root.bounds.copy()
            self.pbounds = self.root.pbounds.copy()

            # --- Check if the mesh doesn't reach the edge of the root grid.
            # --- If not, switch to Dirichlet/Open boundary for ES/EM solver. Also check to make sure
            # --- that one side isn't periodic when the other is.
            if self.l_EM:
                interiorbc = openbc
            else:
                interiorbc = dirichlet
            self.bounds[::2] = where(self.fulllower > 0,interiorbc,self.bounds[::2])
            self.bounds[1::2] = where(self.fullupper < self.root.dimsglobal*self.totalrefinement,
                                      interiorbc,self.bounds[1::2])
            self.bounds[::2] = where((self.bounds[::2] == 2)&(self.bounds[1::2] != 2),
                                     interiorbc,self.bounds[::2])
            self.bounds[1::2] = where((self.bounds[::2] == 2)&(self.bounds[1::2] != 2),
                                      interiorbc,self.bounds[1::2])

            self.l2symtry = self.root.l2symtry
            self.l4symtry = self.root.l4symtry

            self.pbounds[::2] = where(self.fulllower > 0,interiorbc,self.pbounds[::2])
            self.pbounds[1::2] = where(self.fullupper < self.root.dimsglobal*self.totalrefinement,
                                      interiorbc,self.pbounds[1::2])
            self.pbounds[::2] = where((self.pbounds[::2] == 2)&(self.pbounds[1::2] != 2),
                                      interiorbc,self.pbounds[::2])
            self.pbounds[1::2] = where((self.pbounds[::2] == 2)&(self.pbounds[1::2] != 2),
                                       interiorbc,self.pbounds[1::2])

            # --- Create some temporaries for optimization
            self.fullloweroverrefinement = aint(self.fulllower/self.refinement)
            self.fullupperoverrefinement = aint(self.fullupper/self.refinement)

        # --- Set individual quantities based on the values in the arrays,
        # --- if they have been set.
        if self.deltas is not None:
            self.dx = self.deltas[0]
            self.dy = self.deltas[1]
            self.dz = self.deltas[2]
        if self.dims is not None:
            self.nx = self.dims[0]
            self.ny = self.dims[1]
            self.nz = self.dims[2]
        if self.mins is not None:
            self.xmmin = self.mins[0]
            self.ymmin = self.mins[1]
            self.zmmin = self.mins[2]
        if self.maxs is not None:
            self.xmmax = self.maxs[0]
            self.ymmax = self.maxs[1]
            self.zmmax = self.maxs[2]

        # --- Do some further initialization. This assumes that the first
        # --- class inherited from is this class, so the second is the class
        # --- for the solver.
        self.__class__.__bases__[1].__init__(self,**kw)

        if parent is None:
            # --- This is only needed by the root grid in cases when the grid
            # --- parameters are obtained from w3d instead of the argument list.
            self.dims = array([self.nxlocal,self.nylocal,self.nzlocal])
            self.dimsglobal = array([self.nx,self.ny,self.nz])
            self.deltas = array([self.dx,self.dy,self.dz])
            self.mins = array([self.xmminlocal,self.ymminlocal,self.zmminlocal])
            self.maxs = array([self.xmmaxlocal,self.ymmaxlocal,self.zmmaxlocal])
            self.minsglobal = array([self.xmmin,self.ymmin,self.zmmin])
            self.maxsglobal = array([self.xmmax,self.ymmax,self.zmmax])
            self.lower = nint(((self.mins - self.minsglobal)/self.deltas))
            self.upper = nint(((self.maxs - self.minsglobal)/self.deltas))
            self.fulllower = self.lower.copy()
            self.fullupper = self.upper.copy()
            self.rootdims = self.dims

        # --- childdomains is the node centered grid which keeps track of which
        # --- cells are owned by which children. If there are no children,
        # --- then it is not needed.
        self.childdomains = None

        # --- Set up variable time steps
        if top.chdtspid > 0:
            if top.dxpid == 0: top.dxpid = nextpid()
            if top.dypid == 0: top.dypid = nextpid()
            if top.dzpid == 0: top.dzpid = nextpid()

        # --- Get the current global block number and increment the counter.
        # --- Also, add self to the global list of blocks.
        self.blocknumber = self.root.totalnumberofblocks
        self.root.totalnumberofblocks += 1
        self.root.listofblocks.append(self)

        # --- Note that a dictionary is used for the overlaps so that lookups
        # --- are faster, also, the values in the dictionary are list containing
        # --- the domain of the overlap. The blocks with lower and higher block
        # --- number are treated differently, so create separate lists for each.
        # --- Two separate lists is likely only a small optimization.
        self.overlapslower = {}
        self.overlapshigher = {}
        self.overlapsparallelleft = {}
        self.overlapsparallelright = {}

        # --- Make sure that self.lcylindrical is set
        try:
            self.lcylindrical
        except AttributeError:
            self.lcylindrical = (self.solvergeom == w3d.RZgeom)

        # --- Now add any specified children
        self.children = []
        if children is not None:
            for l,u,r in children:
                self.addchild(l,u,refinement=r)

    def __getstate__(self):
        """
    Check whether this instance is the registered solver so that upon unpickling
    it knows whether to re-register itself.
        """
        # --- Make sure the the lreducedpickle option gets propagated to all
        # --- of the blocks.
        self.lreducedpickle = self.root.lreducedpickle
        self.lnorestoreonpickle = self.root.lnorestoreonpickle
        dict = self.__class__.__bases__[1].__getstate__(self)
        dict['_unpicklingcount'] = len(self.root.listofblocks)
        if self.lreducedpickle:
            # --- Remove the big objects from the dictionary. This can be
            # --- regenerated upon the restore.
            dict['childdomains'] = None
        return dict

    def __setstate__(self,dict):
        _unpicklingcount = dict['_unpicklingcount']
        del dict['_unpicklingcount']
        self.__class__.__bases__[1].__setstate__(self,dict)
        try:
            self.root._unpicklingcount -= 1
        except AttributeError:
            self.root._unpicklingcount = _unpicklingcount - 1
        if (self.root._unpicklingcount == 0 and
            self.lreducedpickle and not self.lnorestoreonpickle):
            del self.root._unpicklingcount
            # --- At this point, all of the children have been restored.
            # --- Regenerate childdomains
            self.root.initializechilddomains()
            # --- If source and potential weren't saved, make sure that they are setup.
            # --- Though, this may not always be the right thing to do.
            # --- These can only be done at the end of the restart since only then
            # --- is it gauranteed that the particles are read in.
            if 'sourceparray' not in dict:
                # --- Note that sourcep is only deleted when subsycling is not being
                # --- done and so doesn't then need to be restored.
                installafterrestart(self.root.loadsource)
            installafterrestart(self.root.solve)
        # --- Set the isactive attribute in case this is restored from an old
        # --- dump file.
        if 'isactive' not in self.__dict__:
            self.isactive = 1
        if 'l_EM' not in self.__dict__:
            self.l_EM = 0

    def addchild(self,lower=None,upper=None,fulllower=None,fullupper=None,
                      mins=None,maxs=None,
                      refinement=[2,2,2],nguard=None,nguarddepos=None):
        """
    Add a mesh refined block to this block.
    It is recommended to use mins and maxs to specify the extent of the block.
     - mins,maxs: locations of the grid lower and upper bounds in the beam frame
    Alternatively, lower and upper can be used, but are less flexible and more
    error prone.
     - lower,upper: extent of domain in relative to parent, in its own grid
                    cell size, and not including any guard cells

     - refinement=[2,2,2]: amount of refinement along each axis
     - nguard=[1,1,1] (or value set in base block): number of coarse guard cells
              where rho is deposited but the fields not gather. This is needed to
              reduce the self-field effect of particles near refinement
              boundaries. The default value should be sufficient.
     - nguarddepos=[0,0,0] (or value set in base block): number of coarse guard
              cells where rho (or current) is not deposited. This is needed
              for the refinement of the EM solver.
        """
        # --- A possible thing to do here is set self.root.finalized = false so
        # --- that if a new child is added later, the finalization will be redone
        # --- so that the new child is properly incorporated. This could lead to
        # --- potential problems though, like improperly setup conductors and
        # --- charge density in the child. So it would be better if this
        # --- was forced to be explicit so it is clear what is happening.
        # --- For now, raise an error if a child is being added after finalization.
        # --- This can be circumvented as needed though by expliciting setting
        # --- self.root.finalized = false.
        if self.root.finalized:
            raise RuntimeError('addchild must be called before any operations, such as a field solve, are done')

        if nguard is None: nguard = self.nguard
        if nguarddepos is None: nguarddepos = self.nguarddepos
        # --- Note that all exceptions should be reported.
        child = self.__class__(parent=self,lower=lower,upper=upper,
                               fulllower=fulllower,fullupper=fullupper,
                               mins=mins,maxs=maxs,
                               refinement=refinement,nguard=nguard,
                               nguarddepos=nguarddepos,lchild=True,
                               **self.kw)
        self.children.append(child)
        return child

    def resetroot(self):
        # --- No parents, so just create empty lists
        self.parents = []
        self.root = self
        self.totalnumberofblocks = 1
        self.listofblocks = [self]
        self.childdomains = None
        self.children = []
        self.finalized = false

    #--------------------------------------------------------------------------
    # --- The next several methods handle conductors
    #--------------------------------------------------------------------------

    def installconductor(self,conductor,
                              xmin=None,xmax=None,
                              ymin=None,ymax=None,
                              zmin=None,zmax=None,
                              dfill=top.largepos):

        # --- This check is needed since sometimes during a restore from a pickle,
        # --- this routine may be called by a parent, when it is being restored,
        # --- before this instance is restored, in which case no attributes,
        # --- including 'parents', has been set yet.
        if 'parents' not in self.__dict__: return False

        # --- Call the installconductor from the inherited field solver class
        self.__class__.__bases__[1].installconductor(self,conductor,
                                                          xmin,xmax,ymin,ymax,zmin,zmax,dfill)

        # --- Call installconductor for all of the children.
        for child in self.children:
            child.installconductor(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill)

        return True

    def clearconductors(self,fselfblist=None):
        if not self.isfirstcall(): return
        self.__class__.__bases__[1].clearconductors(self,fselfblist)
        for child in self.children:
            child.clearconductors(fselfblist)

    def getconductors(self,alllevels=1,result=None):
        if result is None: result = []
        result.append(self.getconductorobject())
        if alllevels:
            for child in self.children:
                child.getconductors(alllevels,result)
        return result

    def setconductorvoltage(self,voltage,condid=0,discrete=false,
                            setvinject=false):
        'Recursively calls setconductorvoltage for base and all children'
        if not self.isfirstcall(): return
        self.__class__.__bases__[1].setconductorvoltage(self,voltage,condid,discrete,setvinject)
        for child in self.children:
            child.setconductorvoltage(voltage,condid,discrete,setvinject=false)

    #--------------------------------------------------------------------------
    # --- The next several methods handle initialization that is done after
    # --- all blocks have been added.
    #--------------------------------------------------------------------------

    def finalize(self,lforce=0):
        # --- This should only be called at the top level.
        if self != self.root: return
        if self.finalized and not lforce: return
        blocklists = self.generateblocklevellists()
        self.blocklists = blocklists
        neighborblocklists = self.swapblocklistswithprocessneighbors()
        self.clearparentsandchildren()
        self.findallchildren(blocklists)
        self.initializechilddomains()
        self.findoverlappingsiblings(blocklists[1:],
                                     neighborblocklists[0][1:],
                                     neighborblocklists[1][1:],
                                     neighborblocklists[2][1:],
                                     neighborblocklists[3][1:],
                                     neighborblocklists[4][1:],
                                     neighborblocklists[5][1:])
        self.finalized = 1

    def generateblocklevellists(self,blocklists=None):
        NMAXLEVELS = 32
        if blocklists is None:
            # --- This will only happen at the top level.
            # --- Create a list of empty lists. Each empty list will get the blocks
            # --- at the appropriate levels appended to it. Note that 32 is
            # --- assumed to be a large enough number - there almost certainly
            # --- will never be 32 levels of refinement.
            blocklists = [[] for i in range(NMAXLEVELS)]
        # --- Add this instance to the top level of the list and pass the rest
        # --- of it to the children
        if self not in blocklists[0]:
            blocklists[0].append(self)
            for child in self.children:
                b = child.generateblocklevellists(blocklists[1:])
        return blocklists

    def generatedummyblocklevellists(self,dummyblocklists=None):
        NMAXLEVELS = 32
        if dummyblocklists is None:
            # --- This will only happen at the top level.
            # --- Create a list of empty lists. Each empty list will get the blocks
            # --- at the appropriate levels appended to it. Note that 32 is
            # --- assumed to be a large enough number - there almost certainly
            # --- will never be 32 levels of refinement.
            dummyblocklists = [[] for i in range(NMAXLEVELS)]
        # --- Add this instance to the top level of the list and pass the rest
        # --- of it to the children
        if self not in dummyblocklists[0]:
            dummyself = {}
            dummyself['isactive'] = self.isactive
            dummyself['fulllower'] = self.fulllower
            dummyself['fullupper'] = self.fullupper
            dummyself['blocknumber'] = self.blocknumber
            dummyblocklists[0].append(dummyself)
            for child in self.children:
                b = child.generatedummyblocklevellists(dummyblocklists[1:])
        return dummyblocklists

    def swapblocklistswithprocessneighbors(self):
        NMAXLEVELS = 32
        # --- Note that this will only be called when self = self.root
        if not lparallel:
            # --- Return a bunch of empty lists
            #return [[[] for i in range(NMAXLEVELS)] for i in range(6)]
            return 6*[NMAXLEVELS*[[]]]
        else:
            # --- Create a dummyblocklist
            dummyblocklists = self.generatedummyblocklevellists()

            neighborpes = self.neighborpes

            # --- Exchange lists with processes neighboring along X
            if self.ixproc > 0 and neighborpes[0] >= 0:
                mpisend(dummyblocklists, dest = neighborpes[0])
            if self.ixproc < self.nxprocs-1 and neighborpes[1] >= 0:
                blocklistsxp = mpirecv(source = neighborpes[1])
            else:
                blocklistsxp = NMAXLEVELS*[[]]

            if self.ixproc < self.nxprocs-1 and neighborpes[1] >= 0:
                mpisend(dummyblocklists, dest = neighborpes[1])
            if self.ixproc > 0 and neighborpes[0] >= 0:
                blocklistsxm = mpirecv(source = neighborpes[0])
            else:
                blocklistsxm = NMAXLEVELS*[[]]

            # --- Exchange lists with processes neighboring along Y
            if self.iyproc > 0 and neighborpes[2] >= 0:
                mpisend(dummyblocklists, dest = neighborpes[2])
            if self.iyproc < self.nyprocs-1 and neighborpes[3] >= 0:
                blocklistsyp = mpirecv(source = neighborpes[3])
            else:
                blocklistsyp = NMAXLEVELS*[[]]

            if self.iyproc < self.nyprocs-1 and neighborpes[3] >= 0:
                mpisend(dummyblocklists, dest = neighborpes[3])
            if self.iyproc > 0 and neighborpes[2] >= 0:
                blocklistsym = mpirecv(source = neighborpes[2])
            else:
                blocklistsym = NMAXLEVELS*[[]]

            # --- Exchange lists with processes neighboring along Z
            if self.izproc > 0 and neighborpes[4] >= 0:
                mpisend(dummyblocklists, dest = neighborpes[4])
            if self.izproc < self.nzprocs-1 and neighborpes[5] >= 0:
                blocklistszp = mpirecv(source = neighborpes[5])
            else:
                blocklistszp = NMAXLEVELS*[[]]

            if self.izproc < self.nzprocs-1 and neighborpes[5] >= 0:
                mpisend(dummyblocklists, dest = neighborpes[5])
            if self.izproc > 0 and neighborpes[4] >= 0:
                blocklistszm = mpirecv(source = neighborpes[4])
            else:
                blocklistszm = NMAXLEVELS*[[]]

            return (blocklistsxm,blocklistsxp,
                    blocklistsym,blocklistsyp,
                    blocklistszm,blocklistszp)

    def clearparentsandchildren(self):
        for block in self.listofblocks:
            block.parents = []
            block.children = []
            block.overlapslower = {}
            block.overlapshigher = {}
            block.overlapsparallelleft = {}
            block.overlapsparallelright = {}

    def findallchildren(self,blocklists):
        for block in blocklists[1]:
            if not block.isactive: continue
            # --- Get extent of possible overlapping domain
            l = maximum(block.fullloweroverrefinement,self.fulllower)
            u = minimum(block.fullupperoverrefinement,self.fullupper)
            #if alltrue(u >= l):
            if (u[0] >= l[0] and u[1] >= l[1] and u[2] >= l[2]):
                self.children.append(block)
                block.parents.append(self.blocknumber)

        # --- Only the first block in the list makes the call for the next level.
        # --- This guarantees that this method is called only once for each block.
        if blocklists[0][0] == self:
            for block in blocklists[1]:
                block.findallchildren(blocklists[1:])

    def initializechilddomains(self):
        """
    Sets the regions that are covered by the children.
        """
        if not self.isfirstcall(): return
        self.childdomains = None
        # --- Loop over the children, first calling each, then setting
        # --- childdomain appropriately.
        for child in self.children:
            child.initializechilddomains()

            # --- Set full domain to negative of child number first.
            if self.l_EM:
                l = maximum(self.fulllower,child.fullloweroverrefinement+self.nguarddepos)
                u = child.fullupperoverrefinement-self.nguarddepos
                # --- If the child extends to the edge of the parent mesh, it claims the
                # --- grid points on the upper edges.
                u += where(child.fullupperoverrefinement == self.fullupper,1+array(self.nguarddepos),0)
                l -= where(child.fullloweroverrefinement == self.fulllower,array(self.nguarddepos),0)
            else:
                l = maximum(self.fulllower,child.fullloweroverrefinement)
                u = child.fullupperoverrefinement
                # --- If the child extends to the edge of the parent mesh, it claims the
                # --- grid points on the upper edges.
                #u = u + where(u == self.fullupper,1,0)
                u[u == self.fullupper] += 1
            # --- The child claims all unclaimed areas.
            ii = self.getchilddomains(l,u)
            iit = transpose(ii)
            iit[...] = where(iit==self.blocknumber,-child.blocknumber,iit)

            # --- Set interior to positive child number.
            l = maximum(self.fulllower,(child.lower/child.refinement).astype('l'))
            u = (child.upper/child.refinement).astype('l')
            # --- Check against the case where only the guard cells of the child
            # --- overlap the parent.
            if (u[0] < l[0] or u[1] < l[1] or u[2] < l[2]): continue
            # --- If the child extends to the edge of the parent mesh, it claims the
            # --- grid points on the upper edges.
            #u = u + where(u == self.fullupper,1,0)
            u[u == self.fullupper] += 1
            # --- The child claims its full interior area
            ii = self.getchilddomains(l,u)
            iit = transpose(ii)
            iit[...] = +child.blocknumber

    def findoverlappingsiblings(self,blocklists,blocklistsxm,blocklistsxp,
                                                blocklistsym,blocklistsyp,
                                                blocklistszm,blocklistszp):
        """
    Recursive routine to find, at each level of refinement, all overlapping
    siblings.
    This is faster than the original routine since each pair of blocks is checked
    only once, rather than twice for each parent as in the original.
        """
        # --- When the list is empty, there are no blocks, so just return.
        if (len(blocklists[0]) == 0 and
            len(blocklistsxm[0]) == 0 and len(blocklistsxp[0]) == 0 and
            len(blocklistsym[0]) == 0 and len(blocklistsyp[0]) == 0 and
            len(blocklistszm[0]) == 0 and len(blocklistszp[0]) == 0): return

        # --- Make the call for the next level.
        self.findoverlappingsiblings(blocklists[1:],
                                     blocklistsxm[1:],blocklistsxp[1:],
                                     blocklistsym[1:],blocklistsyp[1:],
                                     blocklistszm[1:],blocklistszp[1:])

        # --- Get a copy of the list (which will be mangled below).
        blocklistscopy = copy.copy(blocklists[0])

        # --- Loop over all blocks.
        for block in blocklists[0]:

            # --- Loop only over blocks that havn't been checked yet.
            del blocklistscopy[0]
            if not block.isactive: continue

            for sibling in blocklistscopy:
                if not sibling.isactive: continue

                # --- Get the area common to the block and its sibling.
                # --- Don't do anything if there is no overlap.
                sl = maximum(sibling.fulllower,block.fulllower)
                su = minimum(sibling.fullupper,block.fullupper)
                if sl[0] > su[0] or sl[1] > su[1] or sl[2] > su[2]: continue

                # --- The ordering ofthe blocks in the list matter since lower blocks
                # --- get precedence, and the source is accumualated there.
                if block.blocknumber < sibling.blocknumber:
                    block.overlapshigher[sibling.blocknumber] = [sl,su]
                    sibling.overlapslower[block.blocknumber] = [sl,su]
                else:
                    block.overlapslower[sibling.blocknumber] = [sl,su]
                    sibling.overlapshigher[block.blocknumber] = [sl,su]

        # --- Now find overlapping blocks in neighboring processors.
        # --- For the precedence, even numbered processors get precdence
        # --- over odd numbered processors. Otherwise, the logic is the same
        # --- as above.
        neighborblocklists = [blocklistsxm[0],blocklistsxp[0],
                              blocklistsym[0],blocklistsyp[0],
                              blocklistszm[0],blocklistszp[0]]
        for block in blocklists[0]:
            if not block.isactive: continue
            for neighborlists,pe in zip(neighborblocklists,self.neighborpes):
                if pe < 0 or pe == npes or pe == me: continue
                for neighborblock in neighborlists:
                    if not neighborblock['isactive']: continue
                    sl = maximum(neighborblock['fulllower'],block.fulllower)
                    su = minimum(neighborblock['fullupper'],block.fullupper)
                    # --- It probably doesn't hurt anything in 3D but will be a small
                    # --- waste of time to include blocks that only overlap on an edge,
                    # --- but those overlaps must be included in 1D and 2D.
                    if sl[0] > su[0] or sl[1] > su[1] or sl[2] > su[2]: continue
                    if pe < me:
                        block.overlapsparallelleft[pe,neighborblock['blocknumber']] = [sl,su]
                    else:
                        block.overlapsparallelright[pe,neighborblock['blocknumber']] = [sl,su]

    def clearinactiveregions(self,nbcells,parent=None,level=1):
        """
    For regions which should not be refined but are included because of the
    coalescing of blocks, the childdomains is set so that the E-fields will
    not be fetched from there (it is set negative).
     - nbcells: array containing the refinement level of the grid cells.
                The shape will be the same as the intersection of self and the
                calling parent.
     - parent: the calling parent
     - level: the total amount refinement
        """
        if parent is None:
            # --- If there is no parent, it is assumed that this is the root block
            l = self.fulllower
            u = self.fullupper
        else:
            # --- Find intersection of parent and self.
            l = maximum(parent.fulllower*self.refinement,self.fulllower)
            u = minimum(parent.fullupper*self.refinement,self.fullupper)

        for child in self.children:
            # --- Find intersection of parent, self, and child
            cl = maximum(child.fullloweroverrefinement,l)
            cu = minimum(child.fullupperoverrefinement,u)
            #if sometrue(cl > cu): continue
            if cl[0] > cu[0] or cl[1] > cu[1] or cl[2] > cu[2]: continue

            # --- Get childdomains in the intersection region, Wherever the
            # --- the refinement level is lower than the childs, force childdomains
            # --- to be negative.
            ii = self.getchilddomains(cl,cu,1)
            nbc = self.getlocalarray(nbcells,cl,cu,fulllower=l)
            # --- Note that if the upper edge of the child does not extend to the
            # --- upper edge of self, then the childdomains at cu will have the
            # --- blocknumber of self, and so should  not be reset. The check
            # --- for (ii == child.blocknumber) prevents this. The check also
            # --- prevents one child from clearing out the domain of another, though
            # --- that wouldn't break anything since that other child would be
            # --- clearing out the same region itself. Also, with the check,
            # --- the second argument of the where can just be -ii instead
            # --- of -abs(ii) since only positive values of ii will have the sign
            # --- changed.
            r = child.refinement
            ii[...] = where((nbc<max(level*r)) & (ii == child.blocknumber),-ii,ii)

            # --- Stretch out the array so it has the refined cell size of the child
            nbcstretched = zeros(1+(cu-cl)*r,'l')
            for k in range(r[2]):
                if k == 0: ksl = slice(None)
                else:      ksl = slice(-1)
                for j in range(r[1]):
                    if j == 0: jsl = slice(None)
                    else:      jsl = slice(-1)
                    for i in range(r[0]):
                        if i == 0: isl = slice(None)
                        else:      isl = slice(-1)
                        nbcstretched[i::r[0],j::r[1],k::r[2]] = nbc[isl,jsl,ksl]

            child.clearinactiveregions(nbcstretched,self,level*r)

    #--------------------------------------------------------------------------
    # --- The next several methods handle the charge density calculation.
    #--------------------------------------------------------------------------

    def setsourcepforparticles(self,*args):
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].setsourcepforparticles(block,*args)

    def setsourceforfieldsolve(self,*args):
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].setsourceforfieldsolve(block,*args)

    def allocatedataarrays(self):
        # --- If not root, than only allocate the arrays of this block
        if self != self.root:
            if self.isactive:
                self.__class__.__bases__[1].allocatedataarrays(self)
            return
        # --- Otherwise, do all blocks.
        # --- Make sure that the final setup was done. This is put here
        # --- since this routine is called by loadsource and solve, which
        # --- call finalize anyway.
        self.finalize()
        # --- Now loop over blocks, calling allocatedataarrays
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].allocatedataarrays(block)

    def zerosource(self):
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].zerosource(block)

    def zerosourcep(self):
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].zerosourcep(block)

    def applysourceboundaryconditions(self):
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].applysourceboundaryconditions(block)

    def averagesourcepwithsubcycling(self):
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].averagesourcepwithsubcycling(block)

    def saveprevioussource(self):
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].saveprevioussource(block)

    def averagewithprevioussource(self,param):
        for block in self.listofblocks:
            if not block.isactive: continue
            self.__class__.__bases__[1].averagewithprevioussource(block,param)

    def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wfact,zgrid,*args,**kw):
        """
    Given the list of particles, a charge and a weight, deposits the charge
    density of the mesh structure.
    This first gets the blocknumber of the block where each of the particles are
    to be deposited. This is then sorted once. The loop is then over the list
    of blocks, rather than walking through the tree structure.
    The sort still takes up about 40% of the time. Note that this depends on
    having ichilddomains filled with the blocknumber rather than the child number
    relative to the parent.
        """
        lrootonly = kw.get('lrootonly',0)
        if len(self.children) > 0 and not lrootonly:
            ichild = zeros(len(x),'l')
            self.getichild(x,y,z,ichild,zgrid)
            x,y,z,ux,uy,uz,gaminv,wfact,nperchild = self.sortbyichild(ichild,x,y,z,ux,uy,uz,gaminv,wfact)

        else:
            nperchild = [len(x)]

        # --- For each block, pass to it the particles in it's domain.
        i = 0
        for block,n in zip(self.root.listofblocks,nperchild):
            if n == 0: continue # does this cause problems? XXX
            self.__class__.__bases__[1].setsourcepatposition(block,x[i:i+n],y[i:i+n],
                                                             z[i:i+n],
                                                             ux[i:i+n],uy[i:i+n],
                                                             uz[i:i+n],gaminv[i:i+n],
                                                             wfact[i:i+n],zgrid,
                                                             *args)
            i = i + n

    def aftersetsourcep(self):
        # --- distribute charge density among blocks
        self.getsourcepfromoverlaps()
        self.gathersourcepfromchildren()
        self.exchangesourcepwithneighbors()
        self.restoresourcepinoverlaps()

    def exchangesourcepwithneighbors(self):
        """
    Exchange sourcep in blocks overlapping blocks on neighboring processors.
        """
        assert self is self.root,"This should only be called by the root block"
        if npes <= 1 or self.l_EM: return
        #if len(self.children) == 0: return

        # --- All blocks first send the overlapping sourcep to the neighbors.
        # --- Note that the precedences for sourcep is handled locally so
        # --- here, all exchanges are made. Note that all data to be
        # --- sent is gathered into a dictionary and sent together.

        # --- Gather all of the data to be sent from all of the refinement blocks
        # --- The separation of the left and right dicts is needed for the end
        # --- case of nxprocs==2 and periodic boundaries, where the processor to
        # --- the left and right is the same.
        senddictsleft = {}
        senddictsright = {}
        for pe in self.neighborpeslist: senddictsleft[pe] = {}
        for pe in self.neighborpeslist: senddictsright[pe] = {}
        for block in self.listofblocks:
            for (pe,othernumber),data in block.overlapsparallelleft.iteritems():
                l,u = data
                sourcep = block.getsourcepslice(l,u)
                senddictsleft[pe].setdefault(othernumber,[]).append((l,u,sourcep))
            for (pe,othernumber),data in block.overlapsparallelright.iteritems():
                l,u = data
                sourcep = block.getsourcepslice(l,u)
                senddictsright[pe].setdefault(othernumber,[]).append((l,u,sourcep))

        # --- Loop over the three axis, in the order x, y, z. Transferring the data
        # --- this way will implicitly handle diagonal neighbors.
        for nel,ner in [[0,1],[2,3],[4,5]]:
            pel = self.neighborpes[nel]
            per = self.neighborpes[ner]

            # --- Loop over the parity - evens send to odds first, then vice versa
            # --- The send/recv pattern is done to avoid lock up since the
            # --- operations are blocking.
            for parity in [0,1]:
                if (self.ixproc+self.iyproc+self.izproc)%2 == parity:
                    if pel in senddictsleft:  mpisend(senddictsleft[pel], dest = pel)
                    if per in senddictsright: mpisend(senddictsright[per], dest = per)
                else:
                    if per in senddictsright: dictfromright = mpirecv(source = per)
                    else:                     dictfromright = {}
                    if pel in senddictsleft:  dictfromleft = mpirecv(source = pel)
                    else:                     dictfromleft = {}

            # --- Create a list of the blocks that receive data.
            blocksreceivingdata = []

            # --- Add in the data from the right
            for blocknumber,data in dictfromright.iteritems():
                block = self.getblockfromnumber(blocknumber)
                for l,u,osourcep in data:
                    ssourcep = block.getsourcepslice(l,u)
                    add(ssourcep,osourcep,ssourcep)
                blocksreceivingdata.append(block)

            # --- The from the left
            for blocknumber,data in dictfromleft.iteritems():
                block = self.getblockfromnumber(blocknumber)
                for l,u,osourcep in data:
                    ssourcep = block.getsourcepslice(l,u)
                    add(ssourcep,osourcep,ssourcep)
                blocksreceivingdata.append(block)

            # --- For the blocks that have received data, go back and zero out
            # --- the regions which overlap blocks with lower number, since those
            # --- lower blocks "own" that data. This is needed to avoid the
            # --- duplication of data that would occur if multiple blocks would
            # --- send the same data onto the next processors. This is only needed
            # --- after the data transfers in x and y, and so is skipped when
            # --- nel == 4.
            if nel < 4:
                for block in blocksreceivingdata:
                    block.zerosourcepinoverlap()

    def getsourcepfromoverlaps(self):
        """
    Add in the sourcep from overlaping areas. The sourcep is gathered into the
    block with the lowerest number. Later on, the sourcep will be copied back to
    the higher numbered blocks. Note that overlaps from neighboring processors
    has already been taken care of. This should only ever be called by the root
    block.
        """
        if self.l_EM: return
        assert self is self.root,"This should only be called by the root block"

        # --- This loops over the blocks in ascending order to ensure that in any
        # --- area with overlap, the block with the lowest number is the one that
        # --- gets the sourcep. This avoids problems of double counting sourcep.
        # --- This could also be done be zeroing out osourcep, but that is extra
        # --- (unecessary) computational work, since it already will be done
        # --- at the end of gathersourcepfromchildren.
        for block in self.listofblocks:
            for othernumber,overlapdomain in block.overlapshigher.iteritems():
                other = block.getblockfromnumber(othernumber)
                l,u = overlapdomain
                ssourcep = block.getsourcepslice(l,u)
                osourcep = other.getsourcepslice(l,u)
                add(ssourcep,osourcep,ssourcep)
                osourcep[...] = 0.

    def sortbyichild(self,ichild,x,y,z,ux,uy,uz,gaminv,wfact):
        nw = len(wfact)
        wfactout = zeros(nw,'d')
        if len(ux) == 0:
            xout,yout,zout,uzout = empty((4,len(x)),'d')
            nperchild = zeros(self.root.totalnumberofblocks,'l')
            sortparticlesbyindex1(len(x),ichild,x,y,z,uz,nw,wfact,
                                  self.root.totalnumberofblocks,
                                  xout,yout,zout,uzout,wfactout,nperchild)
            return xout,yout,zout,ux,uy,uzout,gaminv,wfactout,nperchild
        else:
            xout,yout,zout,uxout,uyout,uzout,gaminvout = zeros((7,len(x)),'d')
            nperchild = zeros(self.root.totalnumberofblocks,'l')
            sortparticlesbyindex2(len(x),ichild,x,y,z,ux,uy,uz,gaminv,nw,wfact,
                                  self.root.totalnumberofblocks,
                                  xout,yout,zout,
                                  uxout,uyout,uzout,gaminvout,wfactout,nperchild)
            return xout,yout,zout,uxout,uyout,uzout,gaminvout,wfactout,nperchild

    def getichild(self,x,y,z,ichild,zgrid):
        """
    Gathers the ichild for the setsourcep.
        """
        # --- This must wait until all of the parents have have set ichild
        # --- so that the value in the children takes precedence.
        if not self.islastcall(): return
        if len(x) == 0: return
        if len(self.children) > 0:
            # --- Convert particles to axisymmetry if needed
            if self is self.root and self.solvergeom == w3d.RZgeom:
                x = sqrt(x**2 + y**2)
                #y = 0.*y
                y = zeros(y.shape,'d')
            # --- Find out whether the particles are in the local domain or one of
            # --- the children's.
            getichild(self.blocknumber,len(x),x,y,z,ichild,
                      self.nxlocal,self.nylocal,self.nzlocal,self.childdomains,
                      self.xmminlocal,self.xmmaxlocal,
                      self.ymminlocal,self.ymmaxlocal,
                      self.zmminlocal,self.zmmaxlocal,zgrid,
                      self.l2symtry,self.l4symtry)
            for child in self.children:
                child.getichild(x,y,z,ichild,zgrid)

    def zerosourcepinoverlap(self):
        """
    This zeros out sourcep in overlapping regions for higher numbered blocks.
    When sourcep is passed from child to parent, in any overlapping regions only
    one child needs to pass the data to the parent.  The choice as to which does
    the passing is determined by the blocknumber - the lower gets to do the
    passing.  For the others, the sourcep in the overlapping region is cleared
    out. That sourcep will be restored later by a call to
    restoresourcepinoverlaps.  Note that this is not recursive, since it is
    called separately by each block from gathersourcepfromchildren.
        """
        if self.l_EM: return
        for othernumber,overlapdomain in self.overlapslower.iteritems():
            l,u = overlapdomain
            ssourcep = self.getsourcepslice(l,u)
            ssourcep[...] = 0.

    def gathersourcepfromchildren(self):
        """
    Fortran version
        """
        # --- Do this only the first time this is called. This should only be
        # --- done once and since each parent requires that this be done
        # --- before it can get its sourcep from here, it must be done on the
        # --- first call.
        if not self.isfirstcall(): return
        # --- Loop over the children
        for child in self.children:

            # --- Make sure that the child has gathered sourcep from its children.
            child.gathersourcepfromchildren()

            # --- Get coordinates of child relative to this domain
            if self.lparallel and self is self.root:
                # --- On the root processor, the extent of the sourcep array can be
                # --- different from the extent of the source array, which is what
                # --- the fulllower and fullupper describe. This difference must be
                # --- taken into account.
                ppdecomp = self.ppdecomp
                fulllower = [ppdecomp.ix[self.ixproc],
                             ppdecomp.iy[self.iyproc],
                             ppdecomp.iz[self.izproc]]
                fullupper = [ppdecomp.ix[self.ixproc]+ppdecomp.nx[self.ixproc],
                             ppdecomp.iy[self.iyproc]+ppdecomp.ny[self.iyproc],
                             ppdecomp.iz[self.izproc]+ppdecomp.nz[self.izproc]]
            else:
                fulllower = self.fulllower
                fullupper = self.fullupper
            l = maximum(child.fullloweroverrefinement,fulllower)
            u = minimum(child.fullupperoverrefinement,fullupper)

            w = child.getwarrayforsourcep()

            pdims = array([self.nxp,self.nyp,self.nzp])
            childpdims = array([child.nxp,child.nyp,child.nzp])
            gathersourcefromchild(self.sourcep,self.ncomponents,
                                  [self.nxguardrho,self.nyguardrho,self.nzguardrho],
                                  pdims,
                                  child.sourcep,childpdims,
                                  l,u,fulllower,child.fulllower,child.fullupper,
                                  child.refinement,w,
                                  self.xmeshlocal,child.xmeshlocal,
                                  self.lcylindrical)

        # --- zerosourcepinoverlap is call here so that any contribution from
        # --- the children in the overlap regions will get zeroed as necessary.
        self.zerosourcepinoverlap()

    def restoresourcepinoverlaps(self):
        """
    Restore sourcep in overlapping areas for blocks which had the sourcep zeroed
    out, the higher numbered blocks. This should only ever be called by the root
    block.
        """
        if self.l_EM: return
        assert self is self.root,"This should only be called by the root block"
        # --- The loop does not need to be in ascending order, but this just
        # --- matches the getsourcepfromoverlaps routine.
        for block in self.listofblocks:
            for othernumber,overlapdomain in block.overlapslower.iteritems():
                other = block.getblockfromnumber(othernumber)
                l,u = overlapdomain
                ssourcep = block.getsourcepslice(l,u)
                osourcep = other.getsourcepslice(l,u)
                ssourcep[...] = osourcep

    #--------------------------------------------------------------------------
    # --- Methods to carry out the field solve
    #--------------------------------------------------------------------------

    def dosolve(self,iwhich=0,zfact=None,isourcepndtscopies=None,indts=None,iselfb=None):
        # --- Make sure that the final setup was done.
        self.finalize()

        # --- Wait until all of the parents have called here until actually
        # --- doing the solve. This ensures that the potential in all of the parents
        # --- which is needed on the boundaries will be up to date.
        if not self.islastcall(): return

        # --- solve on potential, first getting potential from the parents - both the
        # --- boundary conditions and the interior values as the initial
        # --- value.
        if self.l_internal_dosolve:self.setpotentialfromparents()
        self.__class__.__bases__[1].dosolve(self,iwhich,zfact,isourcepndtscopies,indts,iselfb)

        # --- solve for children, using the routine which does the correct
        # --- referencing for subcycling and self-B correction
        for child in self.children:
            child.dosolveonpotential(iwhich,zfact,isourcepndtscopies,indts,iselfb)

        """
        if self == self.root:
          self.__class__.__bases__[1].dosolve(self,iwhich,zfact,isourcepndtscopies,indts,iselfb)
          for blocklist in self.blocklists[1:]:
            if len(blocklist) == 0: break
            tlist = []
            for block in blocklist:
              t = threading.Thread(target=block.dosolveonpotential,args=tuple([iwhich,zfact,isourcepndtscopies,indts,iselfb]))
              t.start()
              tlist.append(t)
              print self.blocknumber,len(tlist)
            for t in tlist:
              t.join()
        else:
          self.setpotentialfromparents()
          self.__class__.__bases__[1].dosolve(self,iwhich,zfact,isourcepndtscopies,indts,iselfb)
        """

    def setpotentialfromparents(self):
        """
    Sets potential, using the values from the parent grid. Setting the full
    potential array gives a better initial guess for the field solver.
        """
        if self.l_EM:return
        for parentnumber in self.parents:
            parent = self.getblockfromnumber(parentnumber)
            # --- Coordinates of mesh relative to parent's mesh location
            # --- and refinement. The minimum and maximum are needed in case
            # --- this mesh extends beyond the parent's.
            plower = (parent.fulllower - self.extradimslower)*self.refinement
            pupper = (parent.fullupper + self.extradimsupper)*self.refinement
            l = maximum(plower,self.fulllower)
            u = minimum(pupper,self.fullupper)
            # --- The full potential arrays are passed in to avoid copying the subsets
            # --- since the fortran needs contiguous arrays.
            gatherpotentialfromparents(self.potential,self.ncomponents,
                                       [self.nxguardphi,self.nyguardphi,self.nzguardphi],
                                       self.dims,l,u,
                                       self.fulllower,
                                       parent.potential,parent.dims,parent.fulllower,
                                       self.refinement)

    #--------------------------------------------------------------------------
    # --- Methods to fetch fields and potential
    #--------------------------------------------------------------------------

    def setpotentialpforparticles(self,*args):
        if self is not self.root:
            self.__class__.__bases__[1].setpotentialpforparticles(self,*args)
        else:
            for block in self.listofblocks:
                if not block.isactive: continue
                self.__class__.__bases__[1].setpotentialpforparticles(block,*args)

    def setfieldpforparticles(self,*args):
        if self is not self.root:
            self.__class__.__bases__[1].setfieldpforparticles(self,*args)
        else:
            for block in self.listofblocks:
                if not block.isactive: continue
                self.__class__.__bases__[1].setfieldpforparticles(block,*args)

    def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
        if not self.finalized: return
        # --- The fetchfield without sorting everything is faster, so use it.
        # --- It is faster because the extra sorting takes a not insignificant
        # --- amount of time, more than unsorting the E arrays.
        self.fetchfieldfrompositionswithoutsort(x,y,z,ex,ey,ez,bx,by,bz,js,pgroup)
        #self.fetchfieldfrompositionswithpsort(x,y,z,ex,ey,ez,bx,by,bz,js,pgroup)

    def fetchfieldfrompositionswithsort(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,
                                        pgroup=None):
        """
    Given the list of particles, fetch the E fields.
    This first gets the blocknumber of the block where each of the particles are
    to be deposited. This is then sorted once. The loop is then over the list
    of blocks, rather than walking through the tree structure.
    The sort takes up about 40% of the time. It is significantly faster
    using the fortran sort.
    Note that this depends on having ichilddomains filled with the
    blocknumber rather than the child number relative to the parent.
    Also, this ends up with the input data remaining sorted.
        """
        if len(self.children) > 0:

            ichild = zeros(len(x),'l')
            # --- This assumes that the root block has blocknumber zero.
            self.getichild_positiveonly(x,y,z,ichild)

            # --- This sorts the particle data in place, including
            # --- the velocities, gaminv, and pid.
            nn = self.root.totalnumberofblocks
            nperchild = zeros(nn,'l')
            particlesortbyindex(pgroup,ichild,0,w3d.ipminfsapi,w3d.npfsapi,
                                nn,nperchild)

        else:
            nperchild = [len(x)]

        # --- For each block, pass to it the particles in it's domain.
        i = 0
        for block,n in zip(self.root.listofblocks,nperchild):
            if n == 0: continue
            self.__class__.__bases__[1].fetchfieldfrompositions(block,
                                               x[i:i+n],y[i:i+n],z[i:i+n],
                                               ex[i:i+n],ey[i:i+n],ez[i:i+n],
                                               bx[i:i+n],by[i:i+n],bz[i:i+n],js,
                                               pgroup)
            if pgroup is not None and top.chdtspid > 0:
                ipmin = w3d.ipminfsapi - 1 + i
                pgroup.pid[ipmin:ipmin+n,top.dxpid-1] = block.dx
                pgroup.pid[ipmin:ipmin+n,top.dypid-1] = block.dy
                pgroup.pid[ipmin:ipmin+n,top.dzpid-1] = block.dz
            i = i + n

    def fetchfieldfrompositionswithoutsort(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,
                                           pgroup=None):
        """
    This is the old version of fetchfieldfrompositions that doesn't rely on having
    access to the particle group and does not sort the input data.
        """

        if pgroup is not None and top.chdtspid > 0:
            ipmin = w3d.ipminfsapi - 1
            n = len(x)
            dx = pgroup.pid[ipmin:ipmin+n,top.dxpid-1]
            dy = pgroup.pid[ipmin:ipmin+n,top.dypid-1]
            dz = pgroup.pid[ipmin:ipmin+n,top.dzpid-1]

        if len(self.children) > 0:

            ichild = zeros(len(x),'l')
            # --- This assumes that the root block has blocknumber zero.
            self.getichild_positiveonly(x,y,z,ichild)

            x,y,z,isort,nperchild = self.sortbyichildgetisort(ichild,x,y,z)

            # --- Create temporary arrays to hold the E field
            tex,tey,tez = zeros((3,len(ex)),'d')
            tbx,tby,tbz = zeros((3,len(bx)),'d')
            if pgroup is not None and top.chdtspid > 0:
                tdx,tdy,tdz = zeros((3,len(bx)),'d')

        else:
            isort = None
            nperchild = [len(x)]
            tex,tey,tez = ex,ey,ez
            tbx,tby,tbz = bx,by,bz
            if pgroup is not None and top.chdtspid > 0:
                tdx,tdy,tdz = dx,dy,dz

        # --- For each block, pass to it the particles in it's domain.
        i = 0
        for block,n in zip(self.root.listofblocks,nperchild):
            if n == 0: continue
            self.__class__.__bases__[1].fetchfieldfrompositions(block,
                                            x[i:i+n],y[i:i+n],z[i:i+n],
                                            tex[i:i+n],tey[i:i+n],tez[i:i+n],
                                            tbx[i:i+n],tby[i:i+n],tbz[i:i+n],js,
                                            pgroup)
            if pgroup is not None and top.chdtspid > 0:
                tdx[i:i+n] = block.dx
                tdy[i:i+n] = block.dy
                tdz[i:i+n] = block.dz
            i = i + n

        # --- Now, put the E fields back into the original arrays, unsorting
        # --- the data
        if isort is not None:
            if len(tex) > 0: addsortedefield(len(tex),isort,tex,tey,tez,ex,ey,ez)
            if len(tbx) > 0: addsortedefield(len(tbx),isort,tbx,tby,tbz,bx,by,bz)
            if pgroup is not None and top.chdtspid > 0:
                # --- Note that these need to be unsorted just like the fields.
                if len(tdx) > 0: addsortedefield(len(tdx),isort,tdx,tdy,tdz,dx,dy,dz)

    def fetchpotentialfrompositions(self,x,y,z,potential):
        if not self.finalized: return
        self.fetchpotentialfrompositionswithoutsort(x,y,z,potential)

    def fetchpotentialfrompositions_old(self,x,y,z,potential):
        """
    Fetches the potential, given a list of positions
        """
        if not self.finalized: return
        if len(x) == 0: return
        if len(self.children) > 0:

            # --- Convert particles to axisymmetry if needed
            if self is self.root and self.solvergeom == w3d.RZgeom:
                x = sqrt(x**2 + y**2)
                y = 0.*y

            # --- Find out whether the particles are in the local domain or one of
            # --- the children's. It is assumed at first to be from the local
            # --- domain, and is only set to one of the childs domains where
            # --- childdomains is positive (which does not include any guard cells).
            ichild = zeros(len(x),'l')
            add(ichild,self.blocknumber,ichild)
            getichildpositiveonly(self.blocknumber,len(x),x,y,z,ichild,
                                  self.nxlocal,self.nylocal,self.nzlocal,
                                  self.childdomains,
                                  self.xmminlocal,self.xmmaxlocal,
                                  self.ymminlocal,self.ymmaxlocal,
                                  self.zmminlocal,self.zmmaxlocal,
                                  self.root.getzgridprv(),
                                  self.l2symtry,self.l4symtry)

            for block in [self]+self.children:
                # --- Get list of particles within the ith domain
                # --- Note that when block==self, the particles selected are the
                # --- ones from this instances domain.
                ii = oldnonzero(ichild==block.blocknumber)
                if len(ii) == 0: continue
                # --- Get positions of those particles.
                xc = take(x,ii)
                yc = take(y,ii)
                zc = take(z,ii)
                # --- Create temporary arrays to hold the potential
                tpotential = zeros(shape(xc),'d')
                # --- Now get the field
                if block == self:
                    self.__class__.__bases__[1].fetchpotentialfrompositions(self,xc,yc,zc,tpotential)
                else:
                    block.fetchpotentialfrompositions(xc,yc,zc,tpotential)
                # --- Put the potential into the passed in arrays
                put(potential,ii,tpotential) # won't work for magnetostatic

        else:

            # --- Get potential from this domain
            self.__class__.__bases__[1].fetchpotentialfrompositions(self,x,y,z,potential)

    def fetchpotentialfrompositionswithoutsort(self,x,y,z,potential):
        """
    This is the version of fetchpotentialfrompositions that doesn't rely on having
    access to the particle group and does not sort the input data.
        """

        if len(self.children) > 0:

            ichild = zeros(len(x),'l')
            # --- This assumes that the root block has blocknumber zero.
            self.getichild_positiveonly(x,y,z,ichild)

            x,y,z,isort,nperchild = self.sortbyichildgetisort(ichild,x,y,z)

            # --- Create temporary arrays to hold the potential
            tpotential = zeros(len(potential),'d')

        else:
            isort = None
            nperchild = [len(x)]
            tpotential = potential

        # --- For each block, pass to it the particles in it's domain.
        i = 0
        for block,n in zip(self.root.listofblocks,nperchild):
            if n == 0: continue
            self.__class__.__bases__[1].fetchpotentialfrompositions(block,
                                            x[i:i+n],y[i:i+n],z[i:i+n],
                                            tpotential[i:i+n])
            i = i + n

        # --- Now, put the potentials back into the original arrays, unsorting
        # --- the data
        if isort is not None:
            if len(tpotential) > 0: addsortedpotential(len(tpotential),isort,tpotential,potential)

    def sortbyichildgetisort(self,ichild,x,y,z):
        xout,yout,zout = empty((3,len(x)),'d')
        isort = empty(len(x),'l')
        nperchild = empty(self.root.totalnumberofblocks,'l')
        sortparticlesbyindexgetisort(len(x),ichild,x,y,z,
                                     self.root.totalnumberofblocks,
                                     xout,yout,zout,isort,nperchild)
        return xout,yout,zout,isort,nperchild

    def getichild_positiveonly(self,x,y,z,ichild):
        """
    Gathers the ichild for the fetchfield_allsort.
        """
        # --- This must wait until all of the parents have have set ichild
        # --- so that the value in the children takes precedence.
        if not self.islastcall(): return
        if len(x) == 0: return
        if len(self.children) > 0:

            # --- Convert particles to axisymmetry if needed
            if self is self.root and self.solvergeom == w3d.RZgeom:
                x = sqrt(x**2 + y**2)
                y = 0.*y

            # --- Find out whether the particles are in the local domain or one of
            # --- the children's.
            getichildpositiveonly(self.blocknumber,len(x),x,y,z,ichild,
                                  self.nxlocal,self.nylocal,self.nzlocal,
                                  self.childdomains,
                                  self.xmminlocal,self.xmmaxlocal,
                                  self.ymminlocal,self.ymmaxlocal,
                                  self.zmminlocal,self.zmmaxlocal,
                                  self.root.getzgridprv(),
                                  self.l2symtry,self.l4symtry)
            for child in self.children:
                child.getichild_positiveonly(x,y,z,ichild)

    #--------------------------------------------------------------------------
    # --- Extra user interface methods
    #--------------------------------------------------------------------------

    def setadvectionvelocity(self,v):
        """Sets the velocity at which the refined grids move relative to this grid,
    Note also that position relative to which the advection is done is reset
    to zero."""
        self.advectionvelocity = v
        self.advectionpostion = 0.
        if self.advectionvelocity == 0:
            if isinstalledbeforeloadrho(self.advectioncontrol):
                uninstallbeforeloadrho(self.advectioncontrol)
        else:
            if not isinstalledbeforeloadrho(self.advectioncontrol):
                installbeforeloadrho(self.advectioncontrol)

    def advectioncontrol(self):
        i1 = nint((self.advectionpostion)/self.dz)
        i2 = nint((self.advectionpostion + self.advectionvelocity*top.dt)/self.dz)
        zzadvect = self.dz*(i2-i1)
        self.advectchildren(zzadvect,i2-i1,root=self)
        self.advectionpostion = self.advectionpostion + self.advectionvelocity*top.dt

    def advectchildren(self,zzadvect,nzadvect,root):
        if nzadvect == 0: return
        if self != root:
            # --- Advect all of the variables of the children holding z positions.
            if not self.isfirstcall(): return
            self.fullloweroverrefinement[-1] += nzadvect
            self.fullupperoverrefinement[-1] += nzadvect
            # --- The nzadvect for the parent was passed in so increase it
            # --- appropriately.
            nzadvect *= self.refinement[-1]
            self.zmmin += zzadvect
            self.zmmax += zzadvect
            self.zmminlocal += zzadvect
            self.zmmaxlocal += zzadvect
            self.zmminp += zzadvect
            self.zmmaxp += zzadvect
            self.fsdecomp.zmin += zzadvect
            self.fsdecomp.zmax += zzadvect
            self.ppdecomp.zmin += zzadvect
            self.ppdecomp.zmax += zzadvect
            self.mins[-1] += zzadvect
            self.maxs[-1] += zzadvect
            self.zmesh[:] += zzadvect
            self.lower[-1] += nzadvect
            self.upper[-1] += nzadvect
            self.fulllower[-1] += nzadvect
            self.fullupper[-1] += nzadvect

        else:
            # --- Advect the data in the childdomains on the root grid.
            childdomains = self.getchilddomains(self.fulllower,self.fullupper,
                                                upperedge=1)
            if nzadvect > 0:
                for iz in range(self.nz,nzadvect-1,-1):
                    childdomains[...,iz] = childdomains[...,iz-nzadvect]
                childdomains[...,:nzadvect] = self.blocknumber
            elif nzadvect < 0:
                for iz in range(0,self.nz-nzadvect+1):
                    childdomains[...,iz] = childdomains[...,iz+nzadvect]
                childdomains[...,self.nz-nzadvect+1:] = self.blocknumber

        # --- Now do the advection in the children
        for child in self.children:
            child.advectchildren(zzadvect,nzadvect,root)

    #--------------------------------------------------------------------------
    # --- Utility methods
    #--------------------------------------------------------------------------

    def getblockfromnumber(self,number):
        return self.root.listofblocks[number]

    def islastcall(self):
        "Returns true when last parent has called"
        try:                   self.ncallsfromparents
        except AttributeError: self.ncallsfromparents = 0
        self.ncallsfromparents += 1
        if self.ncallsfromparents < len(self.parents): return 0
        self.ncallsfromparents = 0
        return 1

    def isfirstcall(self):
        "Returns true when first parent has called"
        try:                   self.ncallsfromparents
        except AttributeError: self.ncallsfromparents = 0
        self.ncallsfromparents += 1
        if self.ncallsfromparents > 1:
            if self.ncallsfromparents == len(self.parents):
                self.ncallsfromparents = 0
            return 0
        # --- This extra check is needed in case there is one or no parent.
        if self.ncallsfromparents >= len(self.parents):
            self.ncallsfromparents = 0
        return 1

    def setname(self,name='c',ichild=None):
        if not self.isfirstcall(): return
        import __main__
        if ichild is not None:
            name = name + '%d'%ichild
            __main__.__dict__[name] = self
        self.mainname = name
        for child,ichild in zip(self.children,range(1,1+len(self.children))):
            child.setname(name,ichild)

    def getmem(self):
        if not self.isfirstcall(): return
        memtot = product(self.dims + 1)
        for child in self.children:
            memtot = memtot + child.getmem()
        return memtot

    def getwarrayforsourcep(self):
        try:
            # --- Return the array if it has already been calculated
            return self.warrayforsourcep
        except AttributeError:
            pass

        # --- Create weight array needed for sourcep deposition.
        # --- Is linear falloff in the weights correct for r > 2?
        r = self.refinement
        wi = [0,0,0]
        for i in range(3):
            if self.dims[i] > 0:
                wi[i] = r[i] - abs(1.*iota(-r[i]+1,+r[i]-1))
                wi[i] = wi[i]/sum(wi[i])
            else:
                wi[i] = ones(2*r[i]-1,'d')
        # --- Expand into 3 dimensions
        w = outer(wi[0],outer(wi[1],wi[2]))
        w.shape = (2*r[0]-1,2*r[1]-1,2*r[2]-1)
        result = fzeros(2*r-1,'d')
        result[...] = w
        self.warrayforsourcep = result
        return result

    def setzgrid(self,zgrid):
        if self == self.root:
            # --- This calls the setzgrid from the field solver base class
            self.__class__.__bases__[1].setzgrid(self,zgrid)
        else:
            # --- This makes sure that for children, the value of zgrid is
            # --- determined by the root grid.
            self.root.setzgrid(zgrid)

    def getzgrid(self):
        if self == self.root:
            # --- This calls the getzgrid from the field solver base class
            return self.__class__.__bases__[1].getzgrid(self)
        else:
            return self.root.getzgrid()

    def getzgridprv(self):
        if self == self.root:
            # --- This calls the getzgridprv from the field solver base class
            return self.__class__.__bases__[1].getzgridprv(self)
        else:
            return self.root.getzgridprv()

    def getzgridndts(self):
        if self == self.root:
            # --- This calls the getzgridndts from the field solver base class
            return self.__class__.__bases__[1].getzgridndts(self)
        else:
            return self.root.getzgridndts()

    def getpotentialpslice(self,lower=None,upper=None,**kw):
#   if lower is None: lower = self.fulllower - array([0,0,1])
#   if upper is None: upper = self.fullupper + array([0,0,1])
        # --- Note that this takes into account the guard cells in z.
        ix1,iy1,iz1 = lower - self.fulllower
        ix2,iy2,iz2 = upper - self.fulllower + 1
        iz1 = iz1 + 1
        iz2 = iz2 + 1
        return self.potentialp[...,ix1:ix2,iy1:iy2,iz1:iz2]
    def getsourcepslice(self,lower=None,upper=None,r=[1,1,1]):
#   if lower is None: lower = self.fulllower
#   if upper is None: upper = self.fullupper
        ix1,iy1,iz1 = lower - self.fulllower
        ix2,iy2,iz2 = upper - self.fulllower + 1
        return self.sourcep[...,ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
    def getpotentialslice(self,lower,upper):
        # --- Note that this takes into account the guard cells in z.
        ix1,iy1,iz1 = (lower - self.fulllower +
                       array([self.nxguardphi,self.nyguardphi,self.nzguardphi]))
        ix2,iy2,iz2 = (upper - self.fulllower + 1 +
                       array([self.nxguardphi,self.nyguardphi,self.nzguardphi]))
        return self.potential[...,ix1:ix2,iy1:iy2,iz1:iz2]
    def getsourceslice(self,lower,upper,r=[1,1,1]):
        ix1,iy1,iz1 = lower - self.fulllower
        ix2,iy2,iz2 = upper - self.fulllower + 1
        return self.source[...,ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
    def getfieldslice(self,lower=None,upper=None,comp=slice(None),r=[1,1,1]):
        if lower is None: lower = self.lower
        if upper is None: upper = self.upper
        if isinstance(comp, basestring):
            ic = ['x','y','z'].index(comp)
        else:
            ic = comp
        ix1,iy1,iz1 = lower - self.fulllower
        ix2,iy2,iz2 = upper - self.fulllower + 1
        self.calcselfe(recalculate=0)
        return self.field[ic,ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
    def getchilddomains(self,lower,upper,upperedge=0):
        if self.childdomains is None:
            #self.childdomains = fzeros(1+self.dims)  + self.blocknumber
            self.childdomains = fzeros(1+self.dims,'l')
            add(self.childdomains,self.blocknumber,self.childdomains)
        ix1,iy1,iz1 = lower - self.fulllower
        ix2,iy2,iz2 = upper - self.fulllower + upperedge
        return self.childdomains[ix1:ix2,iy1:iy2,iz1:iz2]
    def getlocalarray(self,array,lower,upper,r=[1,1,1],fulllower=None,
                      upperedge=1):
        if fulllower is None: fulllower = self.fulllower
        ix1,iy1,iz1 = lower - fulllower
        ix2,iy2,iz2 = upper - fulllower + upperedge
        return array[ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]

    def setmgtol(self,mgtol=None):
        """
    Sets the convergence tolerance for all blocks. If mgtol is not given, it uses
    f3d.mgtol.
        """
        if mgtol is None: mgtol = f3d.mgtol
        self.mgtol = mgtol*ones(self.ncomponents,'l')
        for child in self.children:
            child.setmgtol(mgtol)

    def setmgmaxiters(self,mgmaxiters=None):
        """
    Sets the maximum number of iterations for all blocks. If mgmaxiters is
    not given, it uses f3d.mgmaxiters.
        """
        if mgmaxiters is None: mgmaxiters = f3d.mgmaxiters
        self.mgmaxiters = mgmaxiters
        for child in self.children:
            child.setmgmaxiters(mgmaxiters)

    def find_mgparam(self,lsavephi=false,resetpasses=0):
        for block in self.listofblocks:
            if not block.isactive: continue
            print "Finding mgparam for block number ",block.blocknumber,me
            # --- Temporarily remove the children so the solve is only done
            # --- on this block
            childrensave = block.children
            block.children = []
            self.__class__.__bases__[1].find_mgparam(block,lsavephi=lsavephi,resetpasses=resetpasses)
            # --- Restore the children
            block.children = childrensave

    def walkblocks(self):
        yield self
        for child in self.children:
            for b in child.walkblocks():
                yield b

    def setattrrecursive(self,name,value):
        for block in self.root.listofblocks:
            setattr(block,name,value)

    def arraysliceoperation(self,ip,idim,getdataname,op,opnd,null,comp=None):
        """
    Applies the operator to the array at the specified plane. The blocks only
    contribute within their domains of ownership.
        """
        # --- Each block only needs to check once
        # --- XXX This call breaks something
        #if not self.islastcall(): return null
        # --- Don't do anything if the ip is outside the block or if the child
        # --- is not active
        if (ip < self.fulllower[idim] or ip > self.fullupper[idim] or
            not self.isactive):
            return null
        # --- Get the appropriate slice of potential and the childdomains array
        ii = [slice(None),slice(None),slice(None)]
        ii[idim] = ip - self.fulllower[idim]
        ix,iy,iz = ii
        getdata = getattr(self,getdataname)
        array = getdata(self.fulllower,self.fullupper)

        if not isinstance(comp,types.IntType):
            # --- 'E','B','J','A' will give the field magnitude
            try:
                ic = ['x','y','z','E','B','J','A'].index(comp)
            except ValueError:
                pass
            if self.lcylindrical:
                try:
                    ic = ['r','theta','z','E','B','J','A'].index(comp)
                except ValueError:
                    pass
        else:
            ic = comp
        assert isinstance(ic,types.IntType),"Unrecognized component was input"

        if ic is not None and len(shape(array)) == 4:
            if ic > 2:
                arrayx = array[0,...]
                arrayy = array[1,...]
                arrayz = array[2,...]
                array = sqrt(arrayx**2 + arrayy**2 + arrayz**2)
            else:
                array = array[ic,...]

        if len(self.children) > 0:
            # --- Skip points that self doesn't own
            c = self.getchilddomains(self.fulllower,self.fullupper,1)
            array = where(c[ix,iy,iz]==self.blocknumber,array[ix,iy,iz],null)
        else:
            # --- The transpose is taken so that the array is in C ordering so
            # --- the opnd will be faster.
            array = transpose(array[ix,iy,iz])
        # --- Find the max of self's and the children's potential
        result = opnd(array)
        for child in self.children:
            ipc = ip*child.refinement[idim]
            cresult = child.arraysliceoperation(ipc,idim,getdataname,op,opnd,null,ic)
            result = op(result,cresult)
        return result

    #--------------------------------------------------------------------------
    # --- The following are used for plotting.
    #--------------------------------------------------------------------------

    def genericpf(self,getdataname,kw,idim,pffunc,ip=None):
        """
    Generic plotting routine. This plots only the local domain. Domains of the
    children are also skipped, but the same call is made for them so they will
    be plotted.
        """
        # --- The the block is not active, then don't do anything.
        # --- Note that the root block must execute this function so that
        # --- in parallel, plotlistofthings is called, which is a global
        # --- operations (and must be called by all processors).
        if self is not self.root and not self.isactive: return

        # --- Wait until all parents have called so that the child's domain
        # --- if not overlapped by a parent. This only affects the cellaray plots.
        # --- This also avoids the child plotting multiple times.
        if not self.islastcall(): return

        # --- Get the plane to be plotted
        if ip is None:
            ip = kw.get(('ix','iy','iz')[idim],None)
            if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
        else:
            ip = ip*self.refinement[idim]

        kw[('ix','iy','iz')[idim]] = ip - self.fulllower[idim]

        # --- Set the values of cmin and cmax for all levels. This must be
        # --- done by the root level.
        if self is self.root:
            comp = kw.get('comp','z')
            cmin = kw.get('cmin',None)
            cmax = kw.get('cmax',None)
            cmin_in = cmin
            cmax_in = cmax
            cmin = self.arraysliceoperation(ip,idim,getdataname,min,minnd,+largepos,
                                            comp)
            cmax = self.arraysliceoperation(ip,idim,getdataname,max,maxnd,-largepos,
                                            comp)
            cmin = globalmin(cmin)
            cmax = globalmax(cmax)

            # --- If cmin and/or cmax were not given, but are calculated from the
            # --- data, then apply the same transform to then that will be done to
            # --- the data. Do the transforms after the global min and maxes so
            # --- that possible zero local values don't mess things up.
            try:
                gridscale = kw['gridscale']
            except KeyError:
                gridscale = 1.
            try:
                uselog = kw['uselog']
                if uselog == 'e' or uselog == 1.:
                    logscale = 1.
                else:
                    assert uselog>0,'uselog must be greater than zero'
                    logscale = log(uselog)
            except KeyError:
                logscale = None

            if gridscale < 0.:
                cmin,cmax = cmax,cmin

            if cmin_in is None:
                cmin *= gridscale
                if logscale is not None:
                    if cmin > 0.:
                        cmin = log(cmin)/logscale
                    else:
                        cmin = 0.
            else:
                cmin = cmin_in

            if cmax_in is None:
                cmax *= gridscale
                if logscale is not None:
                    if cmax > 0.:
                        cmax = log(cmax)/logscale
                    else:
                        cmax = 0.
            else:
                cmax = cmax_in

            kw['cmin'] = cmin
            kw['cmax'] = cmax
            kw['local'] = 1

            accumulateplotlists()

        try:
            # --- Only make the plot if the plane is included in the domain.
            # --- Even if there is no overlap, the children must be called since
            # --- they may overlap (in the domain of a different parent).
            # --- Note that the full extent is used so that the childdomains slice
            # --- will have the same shape as ireg.
            if self.fulllower[idim] <= ip and ip <= self.fullupper[idim]:
                if self.childdomains is not None:
                    # --- Create the ireg array, which will be set to zero in the domain
                    # --- of the children.
                    ss = list(shape(self.childdomains))
                    del ss[idim]
                    ireg = zeros(ss,'l')
                    ii = [slice(-1),slice(-1),slice(-1)]
                    ii[idim] = ip - self.fulllower[idim]
                    ix,iy,iz = ii
                    ireg[1:,1:]=equal(self.childdomains[ix,iy,iz],self.blocknumber)
                    if idim != 2: ireg = transpose(ireg)
                    kw['ireg'] = ireg
                else:
                    kw['ireg'] = None
                self.__class__.__bases__[1].genericpf(self,kw,pffunc)
                kw['titles'] = 0
                kw['lcolorbar'] = 0

            for child in self.children:
                child.genericpf(getdataname,kw,idim,pffunc,ip)

        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def drawbox(self,ip=None,idim=2,withguards=1,color=[],selfonly=0):
        if len(color)==0: color=['red', 'green', 'blue', 'cyan', 'magenta','yellow']
        if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
        if self is self.root: accumulateplotlists()
        try:
            if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return
            ii = [0,1,2]
            del ii[idim]
            if withguards:
                i01 = self.mins[ii[0]]
                i02 = self.maxs[ii[0]]
                i11 = self.mins[ii[1]]
                i12 = self.maxs[ii[1]]
            else:
                i01 = self.root.mins[ii[0]] + self.lower[ii[0]]*self.deltas[ii[0]]
                i02 = self.root.mins[ii[0]] + self.upper[ii[0]]*self.deltas[ii[0]]
                i11 = self.root.mins[ii[1]] + self.lower[ii[1]]*self.deltas[ii[1]]
                i12 = self.root.mins[ii[1]] + self.upper[ii[1]]*self.deltas[ii[1]]
            yy = [i01,i01,i02,i02,i01]
            xx = [i11,i12,i12,i11,i11]
            if idim==2:
                yy,xx = xx,yy
            else:
                xx = array(xx) + self.root.getzgrid()
            plg(yy,xx,color=color[0])
            if not selfonly:
                for child in self.children:
                    child.drawbox(ip=ip*child.refinement[idim],idim=idim,
                                  withguards=withguards,color=color[1:])
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def drawboxzy(self,ix=None,withguards=1,color=[],selfonly=0):
        self.drawbox(ip=ix,idim=0,withguards=withguards,color=color,
                     selfonly=selfonly)
    def drawboxzx(self,iy=None,withguards=1,color=[],selfonly=0):
        self.drawbox(ip=iy,idim=1,withguards=withguards,color=color,
                     selfonly=selfonly)
    def drawboxxy(self,iz=None,withguards=1,color=[],selfonly=0):
        self.drawbox(ip=iz,idim=2,withguards=withguards,color=color,
                     selfonly=selfonly)

    def drawfilledbox(self,ip=None,idim=2,withguards=1,ibox=None,selfonly=0):
        if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
        if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return
        ii = [0,1,2]
        del ii[idim]
        if withguards:
            i01 = self.mins[ii[0]]
            i02 = self.maxs[ii[0]]
            i11 = self.mins[ii[1]]
            i12 = self.maxs[ii[1]]
        else:
            i01 = self.root.mins[ii[0]] + self.lower[ii[0]]*self.deltas[ii[0]]
            i02 = self.root.mins[ii[0]] + self.upper[ii[0]]*self.deltas[ii[0]]
            i11 = self.root.mins[ii[1]] + self.lower[ii[1]]*self.deltas[ii[1]]
            i12 = self.root.mins[ii[1]] + self.upper[ii[1]]*self.deltas[ii[1]]
        xx = [i01,i01,i02,i02,i01]
        yy = [i11,i12,i12,i11,i11]
        if idim==2: xx,yy = yy,xx
        if ibox is None: ibox = ones(1,ubyte)
        else:            ibox = (ibox+1).astype(ubyte)
        if self is self.root: accumulateplotlists()
        try:
            plfp(ibox,yy,xx,[5])
            if not selfonly:
                for child in self.children:
                    child.drawfilledbox(ip=ip*child.refinement[idim],idim=idim,
                                        withguards=withguards,ibox=ibox)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def createdxobject(self,kwdict={},**kw):
        """
    Create DX object drawing the object.
      - withguards=1: when true, the guard cells are included in the bounding box
        """
        kw.update(kwdict)
        withguards = kw.get('withguards',1)
        xmin,xmax = self.xmminlocal,self.xmmaxlocal
        ymin,ymax = self.ymminlocal,self.ymmaxlocal
        zmin,zmax = self.zmminlocal,self.zmmaxlocal
        if not withguards:
            ng = self.nguard*self.refinement
            xmin,xmax = xmin+ng[0]*self.dx, xmax-ng[0]*self.dx
            ymin,ymax = ymin+ng[1]*self.dy, ymax-ng[1]*self.dy
            zmin,zmax = zmin+ng[2]*self.dz, zmax-ng[2]*self.dz
        dxlist = [Opyndx.viewboundingbox(xmin,xmax,ymin,ymax,zmin,zmax)]
        for child in self.children:
            dxlist.append(child.getdxobject(kwdict=kw))
        self.dxobject = Opyndx.DXCollection(*dxlist)

    def plpotentialz(self,comp=None,ix=None,iy=None,colors=None,selfonly=0,
                     scale=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if self is self.root: accumulateplotlists()
        try:
            if ix < self.fulllower[0]: return
            if iy < self.fulllower[1]: return
            if ix > self.fullupper[0]: return
            if iy > self.fullupper[1]: return
            ix1 = ix - self.fulllower[0] + self.nxguardphi
            iy1 = iy - self.fulllower[1] + self.nyguardphi
            if scale:
                mesh = self.zmeshlocal
            else:
                iz1 = self.fulllower[2]
                iz2 = self.fullupper[2]
                mesh = arange(iz1,iz2+1)/self.totalrefinement[2]
            if comp is None: ppp = self.potential
            else:            ppp = self.potential[comp,...]
            plg(ppp[ix1,iy1,1:-1],mesh,color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plpotentialz(comp,ix*child.refinement[0],iy*child.refinement[1],
                                       colors=colors,scale=scale)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plpotentialx(self,comp=None,iy=None,iz=None,colors=None,selfonly=0,
                     scale=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if self is self.root: accumulateplotlists()
        try:
            if iy < self.fulllower[1]: return
            if iz < self.fulllower[2]: return
            if iy > self.fullupper[1]: return
            if iz > self.fullupper[2]: return
            iy1 = iy - self.fulllower[1] + self.nyguard
            iz1 = iz - self.fulllower[2] + self.nzguard
            if scale:
                mesh = self.xmeshlocal
            else:
                ix1 = self.fulllower[0]
                ix2 = self.fullupper[0]
                mesh = arange(ix1,ix2+1)/self.totalrefinement[0]
            if comp is None: ppp = self.potential
            else:            ppp = self.potential[comp,...]
            plg(ppp[1:-1,iy1,iz1],mesh,color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plpotentialx(comp,iy*child.refinement[1],iz*child.refinement[2],
                                       colors=colors,scale=scale)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plpotentialy(self,comp=None,ix=None,iz=None,colors=None,selfonly=0,
                     scale=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if self is self.root: accumulateplotlists()
        try:
            if ix < self.fulllower[0]: return
            if iz < self.fulllower[2]: return
            if ix > self.fullupper[0]: return
            if iz > self.fullupper[2]: return
            ix1 = ix - self.fulllower[0] + self.nxguardphi
            iz1 = iz - self.fulllower[2] + self.nzguardphi
            if scale:
                mesh = self.ymeshlocal
            else:
                iy1 = self.fulllower[1]
                iy2 = self.fullupper[1]
                mesh = arange(iy1,iy2+1)/self.totalrefinement[1]
            if comp is None: ppp = self.potential
            else:            ppp = self.potential[comp,...]
            plg(ppp[ix1,1:-1,iz1],mesh,color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plpotentialy(comp,ix*child.refinement[0],iz*child.refinement[2],
                                       colors=colors,scale=scale)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plsourcez(self,comp=None,ix=None,iy=None,colors=None,selfonly=0,scale=1,
                  withboundary=0):
        # --- Note that rho at the child boundaries is incorrect and not used, so
        # --- don't plot it.
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if self is self.root: accumulateplotlists()
        try:
            if ix < self.fulllower[0]: return
            if iy < self.fulllower[1]: return
            if ix > self.fullupper[0]: return
            if iy > self.fullupper[1]: return
            if self != self.root and not withboundary:
                if ix == self.fulllower[0] and self.nx > 0: return
                if iy == self.fulllower[1] and self.ny > 0: return
                if ix == self.fullupper[0] and self.nx > 0: return
                if iy == self.fullupper[1] and self.ny > 0: return
                zslice = slice(1,-1)
            else:
                zslice = slice(None)
            ix1 = ix - self.fulllower[0]
            iy1 = iy - self.fulllower[1]
            if scale:
                mesh = self.zmeshlocal[zslice]
            else:
                iz1 = self.fulllower[2]
                iz2 = self.fullupper[2]
                if self != self.root and not withboundary:
                    iz1 += 1
                    iz2 -= 1
                mesh = arange(iz1,iz2+1)/self.totalrefinement[2]
            if comp is None: rrr = self.source
            else:            rrr = self.source[comp,...]
            plg(rrr[ix1,iy1,zslice],mesh,color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plsourcez(comp,ix*child.refinement[0],iy*child.refinement[1],
                                    colors=colors,scale=scale,withboundary=withboundary)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plsourcex(self,comp=None,iy=None,iz=None,colors=None,selfonly=0,scale=1,
                  withboundary=0):
        # --- Note that source at the child boundaries is incorrect and not used, so
        # --- don't plot it.
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if self is self.root: accumulateplotlists()
        try:
            if iy < self.fulllower[1]: return
            if iz < self.fulllower[2]: return
            if iy > self.fullupper[1]: return
            if iz > self.fullupper[2]: return
            if self != self.root and not withboundary:
                if iy == self.fulllower[1] and self.ny > 0: return
                if iz == self.fulllower[2] and self.nz > 0: return
                if iy == self.fullupper[1] and self.ny > 0: return
                if iz == self.fullupper[2] and self.nz > 0: return
                xslice = slice(1,-1)
            else:
                xslice = slice(None)
            iy1 = iy - self.fulllower[1]
            iz1 = iz - self.fulllower[2]
            if scale:
                mesh = self.xmeshlocal[xslice]
            else:
                ix1 = self.fulllower[0]
                ix2 = self.fullupper[0]
                if self != self.root and not withboundary:
                    ix1 += 1
                    ix2 -= 1
                mesh = arange(ix1,ix2+1)/self.totalrefinement[0]
            if comp is None: rrr = self.source
            else:            rrr = self.source[comp,...]
            plg(rrr[xslice,iy1,iz1],mesh,color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plsourcex(comp,iy*child.refinement[1],iz*child.refinement[2],
                                    colors=colors,scale=scale,withboundary=withboundary)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plsourcey(self,comp=None,ix=None,iz=None,colors=None,selfonly=0,scale=1,
                  withboundary=0):
        # --- Note that source at the child boundaries is incorrect and not used, so
        # --- don't plot it.
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if self is self.root: accumulateplotlists()
        try:
            if ix < self.fulllower[0]: return
            if iz < self.fulllower[2]: return
            if ix > self.fullupper[0]: return
            if iz > self.fullupper[2]: return
            if self != self.root and not withboundary:
                if ix == self.fulllower[0] and self.nx > 0: return
                if iz == self.fulllower[2] and self.nz > 0: return
                if ix == self.fullupper[0] and self.nx > 0: return
                if iz == self.fullupper[2] and self.nz > 0: return
                yslice = slice(1,-1)
            else:
                yslice = slice(None)
            ix1 = ix - self.fulllower[0]
            iz1 = iz - self.fulllower[2]
            if scale:
                mesh = self.ymeshlocal[yslice]
            else:
                iy1 = self.fulllower[1]
                iy2 = self.fullupper[1]
                if self != self.root and not withboundary:
                    iy1 += 1
                    iy2 -= 1
                mesh = arange(iy1,iy2+1)/self.totalrefinement[1]
            if comp is None: rrr = self.source
            else:            rrr = self.source[comp,...]
            plg(rrr[ix1,yslice,iz1],mesh,color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plsourcey(comp,ix*child.refinement[0],iz*child.refinement[2],
                                    colors=colors,scale=scale,withboundary=withboundary)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plfieldz(self,comp=2,ix=None,iy=None,colors=None,selfonly=0,scale=1,
                 withguard=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if withguard:
            lower,upper = self.fulllower,self.fullupper
            iz = slice(None)
        else:
            lower,upper = self.lower,self.upper
            iz = slice(self.lower[2] - self.fulllower[2],
                       self.upper[2] - self.fulllower[2] + 1)
        if self is self.root: accumulateplotlists()
        try:
            if ix < lower[0]: return
            if iy < lower[1]: return
            if ix > upper[0]: return
            if iy > upper[1]: return
            ix1 = ix - self.fulllower[0]
            iy1 = iy - self.fulllower[1]
            if scale:
                mesh = self.zmeshlocal
            else:
                iz1 = self.fulllower[2]
                iz2 = self.fullupper[2]
                mesh = arange(iz1,iz2+1)/self.totalrefinement[2]
            plg(self.field[comp,ix1,iy1,iz],mesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plfieldz(comp,ix*child.refinement[0],iy*child.refinement[1],
                                   colors=colors,scale=scale,withguard=withguard)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plfieldx(self,comp=2,iy=None,iz=None,colors=None,selfonly=0,scale=1,
                 withguard=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if withguard:
            lower,upper = self.fulllower,self.fullupper
            ix = slice(None)
        else:
            lower,upper = self.lower,self.upper
            ix = slice(self.lower[0] - self.fulllower[0],
                       self.upper[0] - self.fulllower[0] + 1)
        if self is self.root: accumulateplotlists()
        try:
            if iy < lower[1]: return
            if iz < lower[2]: return
            if iy > upper[1]: return
            if iz > upper[2]: return
            iy1 = iy - self.fulllower[1]
            iz1 = iz - self.fulllower[2]
            if scale:
                mesh = self.xmeshlocal
            else:
                ix1 = self.fulllower[0]
                ix2 = self.fullupper[0]
                mesh = arange(ix1,ix2+1)/self.totalrefinement[0]
            plg(self.field[comp,ix,iy1,iz1],mesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plfieldx(comp,iy*child.refinement[1],iz*child.refinement[2],
                                   colors=colors,scale=scale,withguard=withguard)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plfieldy(self,comp=2,ix=None,iz=None,colors=None,selfonly=0,scale=1,
                 withguard=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if withguard:
            lower,upper = self.fulllower,self.fullupper
            iy = slice(None)
        else:
            lower,upper = self.lower,self.upper
            iy = slice(self.lower[1] - self.fulllower[1],
                       self.upper[1] - self.fulllower[1] + 1)
        if self is self.root: accumulateplotlists()
        try:
            if ix < lower[0]: return
            if iz < lower[2]: return
            if ix > upper[0]: return
            if iz > upper[2]: return
            ix1 = ix - self.fulllower[0]
            iz1 = iz - self.fulllower[2]
            if scale:
                mesh = self.ymeshlocal
            else:
                iy1 = self.fulllower[1]
                iy2 = self.fullupper[1]
                mesh = arange(iy1,iy2+1)/self.totalrefinement[1]
            plg(self.field[comp,ix1,iy,iz1],mesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plfieldy(comp,ix*child.refinement[0],iz*child.refinement[2],
                                   colors=colors,scale=scale,withguard=withguard)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)














##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MRBlock3D(MeshRefinement,MultiGrid3D):
    """
  Implements adaptive mesh refinement in 3d for the electrostatic field solver
   - parent:
   - refinement=None: amount of refinement along each axis
   - lower,upper: extent of domain in relative to parent, in its own grid
                  cell size, and not including any guard cells
   - dims: dimensions of the grid, only used for root block, the one with
           no parents
   - mins,maxs: locations of the grid lower and upper bounds in the beam frame
   - root: coarsest level grid
   - children: list of tuples, each containing three elements,
               (lower,upper,refinement). Children can also be added later
               using addchild.
   - lreducedpickle=true: when true, a small pickle is made by removing all of
                          the big arrays. The information can be regenerated
                          upon restart.
    """
    def __init__(self,parent=None,refinement=None,
                      lower=None,upper=None,
                      fulllower=None,fullupper=None,
                      dims=None,mins=None,maxs=None,
                      nguard=[1,1,1],
                      children=None,**kw):

        # --- Note that this calls the MultiGrid3D __init__ as well.
        self.__class__.__bases__[0].__init__(self,
                        parent=parent,refinement=refinement,
                        lower=lower,upper=upper,
                        fulllower=fulllower,fullupper=fullupper,
                        dims=dims,mins=mins,maxs=maxs,
                        nguard=nguard,
                        children=children,**kw)

    def getgetdataname(self,kw):
        if kw.get('plotselfe',0):
            return 'getfieldslice'
        elif kw.get('plotrho',0):
            return 'getsourceslice'
        else:
            return 'getpotentialslice'

    def pfxy(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,2,pfxy)
    def pfzx(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzx)
    def pfzy(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,0,pfzy)
    def pfxyg(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,2,pfxyg)
    def pfzxg(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzxg)
    def pfzyg(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,0,pfzyg)

    def plphiz(self,ix=None,iy=None,colors=None,selfonly=0,scale=1):
        self.plpotentialz(ix=ix,iy=iy,colors=colors,selfonly=selfonly,scale=scale)

    def plphix(self,iy=None,iz=None,colors=None,selfonly=0,scale=1):
        self.plpotentialx(iy=iy,iz=iz,colors=colors,selfonly=selfonly,scale=scale)

    def plphiy(self,ix=None,iz=None,colors=None,selfonly=0,scale=1):
        self.plpotentialy(ix=ix,iz=iz,colors=colors,selfonly=selfonly,scale=scale)

    def plrhoz(self,ix=None,iy=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcez(ix=ix,iy=iy,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plrhox(self,iy=None,iz=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcex(iy=iy,iz=iz,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plrhoy(self,ix=None,iz=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcey(ix=ix,iz=iz,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plselfez(self,comp=2,ix=None,iy=None,colors=None,selfonly=0,scale=1,
                 withguard=1):
        self.plfieldz(ix=ix,iy=iy,colors=colors,selfonly=selfonly,scale=scale,
                      withguard=withguard)

    def plselfex(self,comp=2,iy=None,iz=None,colors=None,selfonly=0,scale=1,
                 withguard=1):
        self.plfieldx(iy=iy,iz=iz,colors=colors,selfonly=selfonly,scale=scale,
                      withguard=withguard)

    def plselfey(self,comp=2,ix=None,iz=None,colors=None,selfonly=0,scale=1,
                 withguard=1):
        self.plfieldy(ix=ix,iz=iz,colors=colors,selfonly=selfonly,scale=scale,
                      withguard=withguard)

MRBlock = MRBlock3D


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MRBlock2D(MeshRefinement,MultiGrid2D):
    """
  Implements adaptive mesh refinement in 2d for the electrostatic field solver
   - parent:
   - refinement=None: amount of refinement along each axis
   - lower,upper: extent of domain in relative to parent, in its own grid
                  cell size, and not including any guard cells
   - dims: dimensions of the grid, only used for root block, the one with
           no parents
   - mins,maxs: locations of the grid lower and upper bounds in the beam frame
   - root: coarsest level grid
   - children: list of tuples, each containing three elements,
               (lower,upper,refinement). Children can also be added later
               using addchild.
   - lreducedpickle=true: when true, a small pickle is made by removing all of
                          the big arrays. The information can be regenerated
                          upon restart.
    """
    def __init__(self,parent=None,refinement=None,
                      lower=None,upper=None,
                      fulllower=None,fullupper=None,
                      dims=None,mins=None,maxs=None,
                      nguard=[1,1,1],
                      children=None,**kw):

        # --- Note that this calls the MultiGrid2D __init__ as well.
        self.__class__.__bases__[0].__init__(self,
                        parent=parent,refinement=refinement,
                        lower=lower,upper=upper,
                        fulllower=fulllower,fullupper=fullupper,
                        dims=dims,mins=mins,maxs=maxs,
                        nguard=nguard,
                        children=children,**kw)

    def getgetdataname(self,kw):
        if kw.get('plotselfe',0):
            return 'getfieldslice'
        elif kw.get('plotrho',0):
            return 'getsourceslice'
        else:
            return 'getpotentialslice'

    def drawboxzr(self,iy=None,withguards=1,color=[],selfonly=0):
        self.drawbox(ip=iy,idim=1,withguards=withguards,color=color,
                     selfonly=selfonly)

    def pfzx(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzx)
    def pfzxg(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzxg)

    def plphiz(self,ix=None,colors=None,selfonly=0,scale=1):
        self.plpotentialz(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale)

    def plphix(self,iz=None,colors=None,selfonly=0,scale=1):
        self.plpotentialx(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale)

    def plrhoz(self,ix=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcez(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plrhox(self,iz=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcex(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plselfez(self,comp=2,ix=None,colors=None,selfonly=0,scale=1,withguard=1):
        self.plfieldz(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                      withguard=withguard)

    def plselfex(self,comp=2,iz=None,colors=None,selfonly=0,scale=1,withguard=1):
        self.plfieldx(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                      withguard=withguard)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MRBlockRZ(MeshRefinement,MultiGridRZ):
    """
  Implements adaptive mesh refinement in RZ for the electrostatic field solver
   - parent:
   - refinement=None: amount of refinement along each axis
   - lower,upper: extent of domain in relative to parent, in its own grid
                  cell size, and not including any guard cells
   - dims: dimensions of the grid, only used for root block, the one with
           no parents
   - mins,maxs: locations of the grid lower and upper bounds in the beam frame
   - root: coarsest level grid
   - children: list of tuples, each containing three elements,
               (lower,upper,refinement). Children can also be added later
               using addchild.
   - lreducedpickle=true: when true, a small pickle is made by removing all of
                          the big arrays. The information can be regenerated
                          upon restart.
    """
    def __init__(self,parent=None,refinement=None,
                      lower=None,upper=None,
                      fulllower=None,fullupper=None,
                      dims=None,mins=None,maxs=None,
                      nguard=[1,1,1],
                      children=None,**kw):

        # --- Note that this calls the MultiGridRZ __init__ as well.
        self.__class__.__bases__[0].__init__(self,
                        parent=parent,refinement=refinement,
                        lower=lower,upper=upper,
                        fulllower=fulllower,fullupper=fullupper,
                        dims=dims,mins=mins,maxs=maxs,
                        nguard=nguard,
                        children=children,**kw)

    def getgetdataname(self,kw):
        if kw.get('plotselfe',0):
            return 'getfieldslice'
        elif kw.get('plotrho',0):
            return 'getsourceslice'
        else:
            return 'getpotentialslice'

    def drawboxzr(self,iy=None,withguards=1,color=[],selfonly=0):
        self.drawbox(ip=iy,idim=1,withguards=withguards,color=color,
                     selfonly=selfonly)

    def pfzx(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzx)
    def pfzxg(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzxg)

    def plphiz(self,ix=None,colors=None,selfonly=0,scale=1):
        self.plpotentialz(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale)

    def plphix(self,iz=None,colors=None,selfonly=0,scale=1):
        self.plpotentialx(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale)

    def plrhoz(self,ix=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcez(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plrhox(self,iz=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcex(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plselfez(self,comp=2,ix=None,colors=None,selfonly=0,scale=1,withguard=1):
        self.plfieldz(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                      withguard=withguard)

    def plselfex(self,comp=2,iz=None,colors=None,selfonly=0,scale=1,withguard=1):
        self.plfieldx(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                      withguard=withguard)




##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MRBlock2DDielectric(MeshRefinement,MultiGrid2DDielectric):
    """
  Implements adaptive mesh refinement in 2d for the electrostatic field solver
  with variable dielectric
   - parent:
   - refinement=None: amount of refinement along each axis
   - lower,upper: extent of domain in relative to parent, in its own grid
                  cell size, and not including any guard cells
   - dims: dimensions of the grid, only used for root block, the one with
           no parents
   - mins,maxs: locations of the grid lower and upper bounds in the beam frame
   - root: coarsest level grid
   - children: list of tuples, each containing three elements,
               (lower,upper,refinement). Children can also be added later
               using addchild.
   - lreducedpickle=true: when true, a small pickle is made by removing all of
                          the big arrays. The information can be regenerated
                          upon restart.
    """
    def __init__(self,parent=None,refinement=None,
                      lower=None,upper=None,
                      fulllower=None,fullupper=None,
                      dims=None,mins=None,maxs=None,
                      nguard=[1,1,1],
                      children=None,**kw):

        # --- Note that this calls the MultiGrid2DDielectric __init__ as well.
        self.__class__.__bases__[0].__init__(self,
                        parent=parent,refinement=refinement,
                        lower=lower,upper=upper,
                        fulllower=fulllower,fullupper=fullupper,
                        dims=dims,mins=mins,maxs=maxs,
                        nguard=nguard,
                        children=children,**kw)

    def getgetdataname(self,kw):
        if kw.get('plotselfe',0):
            return 'getfieldslice'
        elif kw.get('plotrho',0):
            return 'getsourceslice'
        else:
            return 'getpotentialslice'

    def drawboxzr(self,iy=None,withguards=1,color=[],selfonly=0):
        self.drawbox(ip=iy,idim=1,withguards=withguards,color=color,
                     selfonly=selfonly)

    def pfzx(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzx)
    def pfzxg(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzxg)

    def plphiz(self,ix=None,colors=None,selfonly=0,scale=1):
        self.plpotentialz(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale)

    def plphix(self,iz=None,colors=None,selfonly=0,scale=1):
        self.plpotentialx(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale)

    def plrhoz(self,ix=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcez(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plrhox(self,iz=None,colors=None,selfonly=0,scale=1,
               withboundary=0):
        self.plsourcex(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale,
                       withboundary=withboundary)

    def plselfez(self,comp=2,ix=None,colors=None,selfonly=0,scale=1,withguard=1):
        self.plfieldz(ix=ix,iy=0,colors=colors,selfonly=selfonly,scale=scale)

    def plselfex(self,comp=2,iz=None,colors=None,selfonly=0,scale=1,withguard=1):
        self.plfieldx(iz=iz,iy=0,colors=colors,selfonly=selfonly,scale=scale)

    def plepsilonz(self,ix=None,colors=None,selfonly=0,scale=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if self is self.root: accumulateplotlists()
        try:
            if ix < self.fulllower[0]: return
            if ix > self.fullupper[0]: return
            ix1 = ix - self.fulllower[0] + self.nxguardphi
            iz1 = self.fulllower[2]
            iz2 = self.fullupper[2] - 1
            if scale:
                mesh = self.zmeshlocal[:-1] + self.dz/2.
            else:
                mesh = arange(iz1,iz2+1)/self.totalrefinement[2] + 0.5
            iz1 += -self.fulllower[2] + self.nzguard
            iz2 += -self.fulllower[2] + self.nzguard
            print self.epsilon[ix1,iz1:iz2+1].shape,mesh.shape
            plg(self.epsilon[ix1,iz1:iz2+1],mesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plepsilonz(ix*child.refinement[0],colors=colors,scale=scale)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plepsilonx(self,iz=None,colors=None,selfonly=0,scale=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if self is self.root: accumulateplotlists()
        try:
            if iz < self.fulllower[2]: return
            if iz > self.fullupper[2]: return
            iz1 = iz - self.fulllower[2] + self.nzguard
            ix1 = self.fulllower[0]
            ix2 = self.fullupper[0] - 1
            if scale:
                mesh = self.xmeshlocal[:-1] + self.dx/2.
            else:
                mesh = arange(ix1,ix2+1)/self.totalrefinement[0] + 0.5
            ix1 += -self.fulllower[0] + self.nxguardphi
            ix2 += -self.fulllower[0] + self.nxguardphi
            plg(self.epsilon[ix1:ix2+1,iz1],mesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plepsilonx(iz*child.refinement[2],colors=colors,scale=scale)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
class MRBlockImplicit2D(MeshRefinement,MultiGridImplicit2D):
    """
  Implements adaptive mesh refinement in 2d for the implicit electrostatic
  field solver
   - parent:
   - refinement=None: amount of refinement along each axis
   - lower,upper: extent of domain in relative to parent, in its own grid
                  cell size, and not including any guard cells
   - dims: dimensions of the grid, only used for root block, the one with
           no parents
   - mins,maxs: locations of the grid lower and upper bounds in the beam frame
   - root: coarsest level grid
   - children: list of tuples, each containing three elements,
               (lower,upper,refinement). Children can also be added later
               using addchild.
   - lreducedpickle=true: when true, a small pickle is made by removing all of
                          the big arrays. The information can be regenerated
                          upon restart.
    """
    def __init__(self,parent=None,refinement=None,
                      lower=None,upper=None,
                      fulllower=None,fullupper=None,
                      dims=None,mins=None,maxs=None,
                      nguard=[1,1,1],
                      children=None,**kw):

        # --- Note that this calls the MultiGridImplicit2D __init__ as well.
        self.__class__.__bases__[0].__init__(self,
                        parent=parent,refinement=refinement,
                        lower=lower,upper=upper,
                        fulllower=fulllower,fullupper=fullupper,
                        dims=dims,mins=mins,maxs=maxs,
                        nguard=nguard,
                        children=children,**kw)

    def getgetdataname(self,kw):
        if kw.get('plotselfe',0):
            return 'getfieldslice'
        elif kw.get('plotrho',0):
            return 'getsourceslice'
        else:
            return 'getpotentialslice'

    def drawboxzr(self,iy=None,withguards=1,color=[],selfonly=0):
        self.drawbox(ip=iy,idim=1,withguards=withguards,color=color,
                     selfonly=selfonly)

    def pfzx(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzx)
    def pfzxg(self,kwdict=None,**kw):
        if kwdict is None: kwdict = {}
        kwdict.update(kw)
        self.genericpf(self.getgetdataname(kw),kwdict,1,pfzxg)

    def plphiz(self,ix=None,colors=None,selfonly=0):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if ix < self.fulllower[0]: return
        if ix > self.fullupper[0]: return
        if self is self.root: accumulateplotlists()
        try:
            plg(self.potential[ix-self.fulllower[0]+1,0,1:-1],self.zmesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plphiz(ix*child.refinement[0],colors=colors)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plphix(self,iz=None,colors=None,selfonly=0):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if iz < self.fulllower[2]: return
        if iz > self.fullupper[2]: return
        if self is self.root: accumulateplotlists()
        try:
            plg(self.potential[1:-1,0,iz-self.fulllower[2]+1],self.xmesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plphix(iz*child.refinement[2],colors=colors)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plrhoz(self,ix=None,colors=None,selfonly=0):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if ix < self.fulllower[0]: return
        if ix > self.fullupper[0]: return
        if self is self.root: accumulateplotlists()
        try:
            plg(self.source[ix-self.fulllower[0],0,:],self.zmesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plrhoz(ix*child.refinement[0],colors=colors)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plrhox(self,iz=None,colors=None,selfonly=0):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if iz < self.fulllower[2]: return
        if iz > self.fullupper[2]: return
        if self is self.root: accumulateplotlists()
        try:
            plg(self.source[:,0,iz-self.fulllower[2]],self.xmesh,
                color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plrhox(iz*child.refinement[2],colors=colors)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plselfez(self,comp=2,ix=None,colors=None,selfonly=0,withguard=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if withguard:
            lower,upper = self.fulllower,self.fullupper
            iz = slice(None)
        else:
            lower,upper = self.lower,self.upper
            iz = slice(self.lower[2] - self.fulllower[2],
                       self.upper[2] - self.fulllower[2] + 1)
        if ix < lower[0]: return
        if ix > upper[0]: return
        if self is self.root: accumulateplotlists()
        try:
            plg(self.field[comp,ix-self.fulllower[0],0,iz],
                self.zmesh[iz],color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plselfez(comp,ix*child.refinement[0],
                                   colors=colors,withguard=withguard)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)

    def plselfex(self,comp=2,iz=None,colors=None,selfonly=0,withguard=1):
        if colors is None: colors = color
        elif not isinstance(colors,collections.Sequence): colors = list([colors])
        if withguard:
            lower,upper = self.fulllower,self.fullupper
            ix = slice(None)
        else:
            lower,upper = self.lower,self.upper
            ix = slice(self.lower[0] - self.fulllower[0],
                       self.upper[0] - self.fulllower[0] + 1)
        if iz < lower[2]: return
        if iz > upper[2]: return
        if self is self.root: accumulateplotlists()
        try:
            plg(self.field[comp,ix,0,iz-self.fulllower[2]],
                           self.xmesh[ix],color=colors[self.blocknumber%len(colors)])
            if not selfonly:
                for child in self.children:
                    child.plselfex(comp,iz*child.refinement[2],
                                   colors=colors,withguard=withguard)
        finally:
            if self is self.root: plotlistofthings(lturnofflist=1)


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
class EMMRBlock(MeshRefinement,EM3D):
    """
  Implements adaptive mesh refinement in 3d for the electromagnetic field solver
   - parent:
   - refinement=None: amount of refinement along each axis
   - lower,upper: extent of domain in relative to parent, in its own grid
                  cell size, and not including any guard cells
   - dims: dimensions of the grid, only used for root block, the one with
           no parents
   - mins,maxs: locations of the grid lower and upper bounds in the beam frame
   - root: coarsest level grid
   - children: list of tuples, each containing three elements,
               (lower,upper,refinement). Children can also be added later
               using addchild.
   - lreducedpickle=true: when true, a small pickle is made by removing all of
                          the big arrays. The information can be regenerated
                          upon restart.
    """
    def __init__(self,parent=None,refinement=None,
                      lower=None,upper=None,
                      fulllower=None,fullupper=None,
                      dims=None,mins=None,maxs=None,
                      nguard=[2,2,2],
                      nguarddepos=[1,1,1],
                      children=None,**kw):
        # --- Note that this calls the EM3D __init__ as well.
        self.__class__.__bases__[0].__init__(self,
                        parent=parent,refinement=refinement,
                        lower=lower,upper=upper,
                        fulllower=fulllower,fullupper=fullupper,
                        dims=dims,mins=mins,maxs=maxs,
                        nguard=nguard,l_EM=1,
                        nguarddepos=nguarddepos,
                        children=children,**kw)

    def solve2ndhalf(self):
        self.allocatedataarrays()
        if self.solveroff:return
        if self.fields.spectral:
            self.move_window_fields()
            return
        if self.mode==2:
            self.solve2ndhalfmode2()
            return
        if self.l_verbose:print 'solve 2nd half',self
        if top.dt != self.dtinit:raise Exception('Time step has been changed since initialization of EM3D.')
        self.push_b_part_2()
        if self.l_pushf:self.exchange_f()
        self.exchange_b()
        self.move_window_fields()
        if self.l_verbose:print 'solve 2nd half done'

    def dosolve(self,iwhich=0,*args):
        self.getconductorobject()
        if self.solveroff:return
        if self.ntsub<1.:
            self.novercycle = nint(1./self.ntsub)
            self.icycle = (top.it-1)%self.novercycle
        else:
            self.novercycle = 1
            self.icycle = 0
        if self.mode==2:
            self.dosolvemode2()
            return
        if any(top.fselfb != 0.):raise Exception('Error:EM solver does not work if fselfb != 0.')
        if self.l_verbose:print 'solve 1st half'
        if top.dt != self.dtinit:raise Exception('Time step has been changed since initialization of EM3D.')
        if self.fields.spectral:
#            self.move_window_fields()
            self.push_spectral_psaotd()
        else:
            self.push_e()
            self.exchange_e()
            for i in range(int(self.ntsub)-1):
                self.push_b_full()
                if self.l_pushf:self.exchange_f()
                self.exchange_b()
                self.push_e_full(i)
                self.exchange_e()
            self.push_b_part_1()
        if self.pml_method==2:
            scale_em3d_bnd_fields(self.block,top.dt,self.l_pushf,self.l_pushg)
        if self.fields.spectral:
            self.exchange_e()
        if self.l_pushf:self.exchange_f()
        self.exchange_b()
        self.setebp()
        if not self.l_nodalgrid and top.efetch[0] != 4:
            if self.l_fieldcenterK:
                self.fieldcenterK()
            else:
                self.yee2node3d()
        if self.l_nodalgrid and top.efetch[0] == 4:
            self.fields.l_nodecentered=True
            self.node2yee3d()
        if self.l_smooth_particle_fields and any(self.npass_smooth>0):
            self.smoothfields()
        self.addsubstractfieldfromparent()
        if self.l_correct_num_Cherenkov:self.smoothfields_poly()
        # --- for fields that are overcycled, they need to be pushed backward every ntsub
        self.push_e(dir=-1)
        self.exchange_e(dir=-1)
        if self.l_verbose:print 'solve 1st half done'


    def step(self,n=1,freq_print=10,lallspecl=0,l_alawarpx=False):
        for i in range(n):
            if top.it%freq_print==0:print 'it = %g time = %g'%(top.it,top.time)
            if lallspecl:
                l_first=l_last=1
            else:
                if i==0:
                    l_first=1
                else:
                    l_first=0
                if i==n-1:
                    l_last=1
                else:
                    l_last=0
            if l_alawarpx:
                self.onestep_alawarpx(l_first,l_last)
            else:
                self.onestep(l_first,l_last)

    def onestep_alawarpx(self,l_first,l_last):
        if self.ntsub<1.:
            self.novercycle = nint(1./self.ntsub)
            self.icycle = (top.it-1)%self.novercycle
        else:
            self.novercycle = 1
            self.icycle = 0
        self.__class__.__bases__[1].novercycle=self.novercycle
        self.__class__.__bases__[1].icycle=self.icycle
        
        # --- call beforestep functions
        callbeforestepfuncs.callfuncsinlist()

        # --- push by half step if first step (where all are synchronized)
        if l_first:
            if not self.solveroff:
                self.allocatedataarrays()
                self.push_e(l_half=True)
                self.exchange_e()
            for js in range(top.pgroup.ns):
                self.push_positions(js,dtmult=0.5)

        if not self.solveroff:
            for i in range(int(self.ntsub)-1):
                self.push_b_full()
                if self.l_pushf:self.exchange_f()
                self.exchange_b()
                self.push_e_full(i)
                self.exchange_e()
            self.push_b_part_1()

            if self.pml_method==2:
                scale_em3d_bnd_fields(self.block,top.dt,self.l_pushf,self.l_pushg)

            if self.l_pushf:self.exchange_f()
            self.exchange_b()
            self.setebp()

            if not self.l_nodalgrid and top.efetch[0] != 4:
                if self.l_fieldcenterK:
                    self.fieldcenterK()
                else:
                    self.yee2node3d()
            if self.l_nodalgrid and top.efetch[0] == 4:
                self.fields.l_nodecentered=True
                self.node2yee3d()
            if self.l_smooth_particle_fields and any(self.npass_smooth>0):
                self.smoothfields()
            self.addsubstractfieldfromparent()
            if self.l_correct_num_Cherenkov:self.smoothfields_poly()

        top.zgrid+=top.vbeamfrm*top.dt
        top.zbeam=top.zgrid

        w3d.pgroupfsapi = top.pgroup
        for js in range(top.pgroup.ns):
            self.fetcheb(js)
            self.push_velocity_full(js)
            self.record_old_positions(js)
            self.push_positions(js)

        inject3d(1, top.pgroup)

        # --- call user-defined injection routines
        userinjection.callfuncsinlist()

        particleboundaries3d(top.pgroup,-1,False)

        self.solve2ndhalf()

        # --- call beforeloadrho functions
        if isinstalledbeforeloadrho(self.solve2ndhalf):
            uninstallbeforeloadrho(self.solve2ndhalf)
        beforeloadrho.callfuncsinlist()

#        self.loadsource()
        self.loadrho()
        self.loadj()

        self.getconductorobject()

        if l_last:
            for js in range(top.pgroup.ns):
                self.record_old_positions(js)
                self.push_positions(js,dtmult=-0.5)
            if not self.solveroff:self.push_e(l_half=True)
        else:
            if not self.solveroff:
                self.push_e()
                self.exchange_e()

        # --- update time, time counter
        top.time+=top.dt
        if top.it%top.nhist==0:
#           zmmnt()
           minidiag(top.it,top.time,top.lspecial)
        top.it+=1

        # --- call afterstep functions
        callafterstepfuncs.callfuncsinlist()
        
    def push_e(self,dir=1,l_half=False):
        for child in self.children:
            child.push_e(dir,l_half)
        self.__class__.__bases__[1].push_e(self,dir,l_half)

    def exchange_e(self,dir=1.):
        for child in self.children:
            child.exchange_e(dir)
        self.__class__.__bases__[1].exchange_e(self,dir)

    def push_b_part_1(self,dir=1):
        for child in self.children:
            child.push_b_part_1(dir)
        self.__class__.__bases__[1].push_b_part_1(self,dir)

    def push_b_part_2(self):
        for child in self.children:
            child.push_b_part_2()
        self.__class__.__bases__[1].push_b_part_2(self)

    def exchange_b(self,dir=1):
        for child in self.children:
            child.exchange_b(dir)
        self.__class__.__bases__[1].exchange_b(self,dir)

    def exchange_f(self,dir=1):
        for child in self.children:
            if child.l_pushf:
                child.exchange_f(dir)
        if self.l_pushf:
            self.__class__.__bases__[1].exchange_f(self,dir)

    def setebp(self):
        for child in self.children:
            child.setebp()
        self.__class__.__bases__[1].setebp(self)

    def yee2node3d(self):
        for child in self.children:
            child.yee2node3d()
        self.__class__.__bases__[1].yee2node3d(self)

    def node2yee3d(self):
        for child in self.children:
            child.node2yee3d()
        self.__class__.__bases__[1].node2yee3d(self)

    def push_eb_subcycle(self):
        if self is self.root:
            self.push_ebsubcycle_children([self])
        if len(self.children)>0:
            self.push_ebsubcycle_children(self.children)
            for child in self.children:
                child.push_eb_subcycle()

    def push_ebsubcycle_children(self,children):
        c0 = children[0]
        ntsub = c0.ntsub
        for i in range(int(ntsub)-1):
            for c in children:
                self.__class__.__bases__[1].push_b_full(c)
                if c.l_pushf:
                    self.__class__.__bases__[1].exchange_f(c)
                self.__class__.__bases__[1].exchange_b(c)
                self.__class__.__bases__[1].push_e_full(c,i)
                self.__class__.__bases__[1].exchange_e(c)
        if c0.refinement is not None:
            ntsub = c0.field_coarse.ntsub
            for i in range(int(ntsub)-1):
                for c in children:
                    self.__class__.__bases__[1].push_b_full(c.field_coarse)
                    self.__class__.__bases__[1].exchange_f(c.field_coarse)
                    self.__class__.__bases__[1].exchange_b(c.field_coarse)
                    self.__class__.__bases__[1].push_e_full(c.field_coarse,i)
                    self.__class__.__bases__[1].exchange_e(c.field_coarse)

    def installconductor(self,conductor,
                              xmin=None,xmax=None,
                              ymin=None,ymax=None,
                              zmin=None,zmax=None,
                              dfill=top.largepos):

        # --- Call the installconductor from the inherited MeshRefinement class.
        # --- This will return True only if it is Ok calling it at this time.
        result = self.__class__.__bases__[0].installconductor(self,conductor,
                                                          xmin,xmax,ymin,ymax,zmin,zmax,dfill)
        if result:
            # --- Call installconductor for the field_coarse
            if self.refinement is not None:
                self.field_coarse.installconductor(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill)
        return result

    def getconductors(self,alllevels=1,result=None):
        result = self.__class__.__bases__[0].getconductors(self,alllevels,result)
        if self.refinement is not None:
            result.append(self.field_coarse.getconductorobject())
        return result

    def addsubstractfieldfromparent(self):
        """
    Add own field and field from parent, substracting field from field_coarse.block, and
    putting the result in Exp, Eyp, Ezp, Bxp, Byp and Bzp.
        """
        for parentnumber in self.parents:
            parent = self.getblockfromnumber(parentnumber)
            # --- Coordinates of mesh relative to parent's mesh location
            # --- and refinement. The minimum and maximum are needed in case
            # --- this mesh extends beyond the parent's.
            plower = (parent.fulllower - self.extradimslower)*self.refinement
            pupper = (parent.fullupper + self.extradimsupper)*self.refinement
            l = maximum(plower,self.fulllower)
            u = minimum(pupper,self.fullupper)
            lp = (l-plower)/self.refinement
            if top.efetch[0] != 4:
                addsubstractfields_nodal(self.block,self.field_coarse.block,parent.block,lp,self.refinement,self.l_2dxz)
            else:
                addsubstractfields(self.block,self.field_coarse.block,parent.block,lp,self.refinement,self.l_2dxz)
        for child in self.children:
            child.addsubstractfieldfromparent()

    def move_cells_x(self,n):
        self.__class__.__bases__[1].move_cells_x(self,n)
        if self.refinement is not None:
            self.__class__.__bases__[1].move_cells_x(self.field_coarse,n/self.refinement[0])
        dx = self.dx*n
        self.mins[0]+=dx
        self.maxs[0]+=dx
        for child in self.children:
            child.move_cells_x(n*child.refinement[0])

    def move_cells_y(self,n):
        self.__class__.__bases__[1].move_cells_y(self,n)
        if self.refinement is not None:
            self.__class__.__bases__[1].move_cells_y(self.field_coarse,n/self.refinement[1])
        dy = self.dy*n
        self.mins[1]+=dy
        self.maxs[1]+=dy
        for child in self.children:
            child.move_cells_y(n*child.refinement[1])

    def move_cells_z(self,n):
        self.__class__.__bases__[1].move_cells_z(self,n)
        if self.refinement is not None:
            self.__class__.__bases__[1].move_cells_z(self.field_coarse,n/self.refinement[2])
        dz = self.dz*n
        self.mins[2]+=dz
        self.maxs[2]+=dz
        for child in self.children:
            child.move_cells_z(n*child.refinement[2])

    def shift_cells_z(self,n):
        self.__class__.__bases__[1].shift_cells_z(self,n)
        if self.refinement is not None:
            self.__class__.__bases__[1].shift_cells_z(self.field_coarse,n/self.refinement[2])
        dz = self.dz*n
        self.mins[2]+=dz
        self.maxs[2]+=dz
        for child in self.children:
            child.shift_cells_z(n*child.refinement[2])

    def finalize(self,lforce=False):
        if self != self.root: return
        if self.finalized and not lforce: return
        self.__class__.__bases__[0].finalize(self) # MeshRefinement.finalize
        for block in self.listofblocks:
            self.__class__.__bases__[1].finalize(block,lforce=True) # EM3D.finalize
            if block != self.root:
                self.__class__.__bases__[1].finalize(block.field_coarse,lforce=True) # EM3D.finalize
        self.setbcoverlaps()
        self.checkconnections()
        self.fillchilddomains()
        self.setdtinit()
        self.finalized = True # Redundant, but here for clarity

    def aftersetsourcep(self):
        # --- distribute charge density among blocks
        self.gathersourcepfromchildren()

    def finalizesourcep(self):
        if self.l_verbose:print 'finalizesourcep',self.sourcepfinalized
        if self.sourcepfinalized: return
        self.sourcepfinalized = True
        self.add_source_ndts_slices()
        self.aftersetsourcep()
        # --- smooth current density
        self.smoothdensity()
        # -- add laser 
        self.add_laser(self.block.core.yf)
        self.applysourceboundaryconditions()
        if self.l_verbose:print 'finalizesourcep done'

    def gathersourcepfromchildren(self):
        """
    Fortran version
        """
        # --- Do this only the first time this is called. This should only be
        # --- done once and since each parent requires that this be done
        # --- before it can get its sourcep from here, it must be done on the
        # --- first call.
        if not self.isfirstcall(): return
        # --- Loop over the children
        for child in self.children:

            # --- Make sure that the child has gathered sourcep from its children.
            child.gathersourcepfromchildren()

            # --- Get coordinates of child relative to this domain
            fulllower = self.fulllower
            fullupper = self.fullupper
            l = maximum(child.fullloweroverrefinement,fulllower)
            u = minimum(child.fullupperoverrefinement,fullupper)

            w = child.getwarrayforsourcep()

            # --- project charge and current density from fine patch to coarse twin + parent
            cb = child.block
            cbc = child.field_coarse.block
            lp = l-fulllower
            project_jxjyjz(cb.core.yf.Jx,cb.core.yf.Jy,cb.core.yf.Jz,
                           cbc.core.yf.Jxarray[...,0],cbc.core.yf.Jyarray[...,0],
                           cbc.core.yf.Jzarray[...,0],
                           self.block.core.yf.Jxarray[...,0],
                           self.block.core.yf.Jyarray[...,0],
                           self.block.core.yf.Jzarray[...,0],
                           cb.nx,cb.ny,cb.nz,
                           self.block.nx,self.block.ny,self.block.nz,
                           cb.nxguard,cb.nyguard,cb.nzguard,
                           child.refinement[0],
                           child.refinement[1],
                           child.refinement[2],
                           lp[0],lp[1],lp[2],self.l_2dxz,
                           self.icycle,self.novercycle)
            cbc.core.yf.Jx = cbc.core.yf.Jxarray[...,0]
            cbc.core.yf.Jy = cbc.core.yf.Jyarray[...,0]
            cbc.core.yf.Jz = cbc.core.yf.Jzarray[...,0]
            self.block.core.yf.Jx = self.block.core.yf.Jxarray[...,0]
            self.block.core.yf.Jy = self.block.core.yf.Jyarray[...,0]
            self.block.core.yf.Jz = self.block.core.yf.Jzarray[...,0]
            if child.l_pushf:
                project_rho(cb.core.yf.Rho,
                             cbc.core.yf.Rhoarray[...,0],
                             self.block.core.yf.Rho,
                             cb.nx,cb.ny,cb.nz,
                             self.block.nx,self.block.ny,self.block.nz,
                             cb.nxguard,cb.nyguard,cb.nzguard,
                             child.refinement[0],
                             child.refinement[1],
                             child.refinement[2],
                             lp[0],lp[1],lp[2],self.l_2dxz)
                cbc.core.yf.Rho = cbc.core.yf.Rhoarray[...,0]

            # --- apply density mask to fine patch and twin charge and current densities
            if 0:
                ntrans = 4
                child.Jxbf = cb.core.yf.Jx.copy()
                child.Jybf = cb.core.yf.Jy.copy()
                child.Jzbf = cb.core.yf.Jz.copy()
                child.Jxbfc = cbc.core.yf.Jxarray[...,0].copy()
                child.Jybfc = cbc.core.yf.Jyarray[...,0].copy()
                child.Jzbfc = cbc.core.yf.Jzarray[...,0].copy()
                apply_dmask(cb.core.yf.Rho,
                            cb.core.yf.Jx,
                            cb.core.yf.Jy,
                            cb.core.yf.Jz,
                            cb.core.yf.dmaskx,
                            cb.core.yf.dmasky,
                            cb.core.yf.dmaskz,
                            child.bounds,child.nguarddepos*child.refinement,child.refinement*ntrans,
                            cb.nx,cb.ny,cb.nz,
                            cb.nxguardphi,cb.nyguardphi,cb.nzguardphi,
                            self.l_pushf,self.l_2dxz)
                apply_dmask(cbc.core.yf.Rhoarray[...,0],
                            cbc.core.yf.Jxarray[...,0],
                            cbc.core.yf.Jyarray[...,0],
                            cbc.core.yf.Jzarray[...,0],
                            cbc.core.yf.dmaskx,
                            cbc.core.yf.dmasky,
                            cbc.core.yf.dmaskz,
                            child.bounds,child.nguarddepos,array([1,1,1])*ntrans,
                            cbc.nx,cbc.ny,cbc.nz,
                            cbc.nxguardphi,cbc.nyguardphi,cbc.nzguardphi,
                            self.l_pushf,self.l_2dxz)
                child.Jxaf = cb.core.yf.Jx.copy()
                child.Jyaf = cb.core.yf.Jy.copy()
                child.Jzaf = cb.core.yf.Jz.copy()
                child.Jxafc = cbc.core.yf.Jxarray[...,0].copy()
                child.Jyafc = cbc.core.yf.Jyarray[...,0].copy()
                child.Jzafc = cbc.core.yf.Jzarray[...,0].copy()

    def add_source_ndts_slices(self):
        for child in self.children:
            child.add_source_ndts_slices()
        self.__class__.__bases__[1].add_source_ndts_slices(self)

    def smoothdensity(self):
        for child in self.children:
            child.smoothdensity()
        self.__class__.__bases__[1].smoothdensity(self)
        if self.refinement is not None:
            self.__class__.__bases__[1].smoothdensity(self.field_coarse)

    def smoothfields(self):
        for child in self.children:
            if child.l_smooth_particle_fields and any(child.npass_smooth>0):
                child.smoothfields()
        if self.l_smooth_particle_fields and any(self.npass_smooth>0):
            self.__class__.__bases__[1].smoothfields(self)

    def setbcoverlaps(self):
        self.__class__.__bases__[1].setbcoverlaps(self)
        for child in self.children:
            child.setbcoverlaps()

    def checkconnections(self):
        self.__class__.__bases__[1].checkconnections(self)
        for child in self.children:
            child.checkconnections()

    def fillchilddomains(self):
        if self is not self.root:self.__class__.__bases__[1].fillchilddomains(self)
        for child in self.children:
            child.fillchilddomains()

    def setdtinit(self):
        self.__class__.__bases__[1].setdtinit(self)
        for child in self.children:
            child.setdtinit()

    ##########################################################################
    # Define the basic plot commands
    def fetchcmincmax(self,dataname,l_children,guards,kw):
        overlap = True
        direction = kw.get('direction',None)
        slice = kw.get('slice',None)
        l_abs = kw.get('l_abs',False)
        cmin_in = kw.get('cmin',None)
        cmax_in = kw.get('cmax',None)
        if cmin_in is not None: cmin = cmin_in
        if cmax_in is not None: cmax = cmax_in
        if cmin_in is None or cmax_in is None:
            data = getattr(self,dataname)(guards,overlap)
            if self.l_2dxz and guards:data=data[:,0,:]
            slice,dataslice = self.getdataatslice(data,direction,slice,l_abs)
            if cmin_in is None: cmin = minnd(dataslice)
            if cmax_in is None: cmax = maxnd(dataslice)
            if l_children:
                for i in range(1,len(self.blocklists)):
                    for c in self.blocklists[i]:
                        if c.isactive:
                            data = getattr(c,dataname)(guards,overlap)
                            if self.l_2dxz and guards:data=data[:,0,:]
                            slice,dataslice = c.getdataatslice(data,direction,slice,l_abs)
                            if cmin_in is None: cmin = min(cmin,minnd(dataslice))
                            if cmax_in is None: cmax = max(cmax,maxnd(dataslice))
            if cmin_in is None: cmin = globalmin(cmin)
            if cmax_in is None: cmax = globalmax(cmax)
            gridscale = kw.get('gridscale',None)
            if cmin is not None and gridscale is not None: cmin *= gridscale
            if cmax is not None and gridscale is not None: cmax *= gridscale
        return slice,cmin,cmax

    def genericpfem3dMR(self,title='',dataname='',l_children=1,guards=0,guardMR=0,**kw):
        overlap=True
        slice,cmin,cmax = self.fetchcmincmax(dataname,l_children,guards,kw)
        kw['cmin'] = cmin
        kw['cmax'] = cmax
        kw['slice'] = slice

        data = getattr(self,dataname)(guards,overlap)
        self.genericpfem3d(data,title,**kw)

        if l_children:
            for i in range(1,len(self.blocklists)):
                for c in self.blocklists[i]:
                    if c.isactive:
                        if guardMR:
                            guardMR=array([0,0,0])
                        else:
                            guardMR=c.nguard*c.refinement
                        xmin = c.block.xmin+guardMR[0]*c.block.dx
                        xmax = c.block.xmax-guardMR[0]*c.block.dx
                        ymin = c.block.ymin+guardMR[1]*c.block.dy
                        ymax = c.block.ymax-guardMR[1]*c.block.dy
                        zmin = c.block.zmin+guardMR[2]*c.block.dz
                        zmax = c.block.zmax-guardMR[2]*c.block.dz
                        data = getattr(c,dataname)(guards,overlap,guardMR)
                        c.genericpfem3d(data, \
                        title,xmmin=xmin,xmmax=xmax,\
                        ymmin=ymin,ymmax=ymax,\
                        zmmin=zmin,zmmax=zmax,**kw)
                    else:
                        c.genericpfem3d(None,title,**kw)


    def getarray(self,g,guards=0,overlap=0,guardMR=[0,0,0]):
        if guards:
            return g
        else:
            f=self.fields
            ox=oy=oz=0
            if not overlap:
                if self.block.xrbnd==em3d.otherproc:ox=1
                if self.block.yrbnd==em3d.otherproc:oy=1
                if self.block.zrbnd==em3d.otherproc:oz=1
            if self.l_1dz:
                return g[0,0,f.nzguard+guardMR[2]:-f.nzguard-oz-guardMR[2]]
            elif self.l_2dxz:
                return g[f.nxguard+guardMR[0]:-f.nxguard-ox-guardMR[0],0, \
                         f.nzguard+guardMR[2]:-f.nzguard-oz-guardMR[2]]
            else:
                return g[f.nxguard+guardMR[0]:-f.nxguard-ox-guardMR[0], \
                         f.nyguard+guardMR[1]:-f.nyguard-oy-guardMR[1], \
                         f.nzguard+guardMR[2]:-f.nzguard-oz-guardMR[2]]

    def getex(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Exp,guards,overlap,guardMR)

    def getey(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Eyp,guards,overlap,guardMR)

    def getez(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Ezp,guards,overlap,guardMR)

    def getbx(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Bxp,guards,overlap,guardMR)

    def getby(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Byp,guards,overlap,guardMR)

    def getbz(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Bzp,guards,overlap,guardMR)

    def getexg(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Ex,guards,overlap,guardMR)

    def geteyg(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Ey,guards,overlap,guardMR)

    def getezg(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Ez,guards,overlap,guardMR)

    def getbxg(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Bx,guards,overlap,guardMR)

    def getbyg(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.By,guards,overlap,guardMR)

    def getbzg(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Bz,guards,overlap,guardMR)

    def getjx(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Jx,guards,overlap,guardMR)

    def getjy(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Jy,guards,overlap,guardMR)

    def getjz(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Jz,guards,overlap,guardMR)

    def getrho(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.Rho,guards,overlap,guardMR)

    def getf(self,guards=0,overlap=0,guardMR=[0,0,0]):
        return self.getarray(self.fields.F,guards,overlap,guardMR)

    def getdive(self,guards=0,overlap=0,guardMR=[0,0,0]):
        dive = zeros(shape(self.fields.Ex),'d')
        f = self.fields
        if top.efetch[0] != 4:node2yee3d(f)
        getdive(f.Ex,f.Ey,f.Ez,dive,f.dx,f.dy,f.dz,
                f.nx,f.ny,f.nz,f.nxguard,f.nyguard,f.nzguard,
                f.xmin,
                self.l_2dxz,self.l_2drz,self.l_nodalgrid)
        if top.efetch[0] != 4:yee2node3d(f)
        return self.getarray(dive,guards,overlap,guardMR)

    def pfex(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('E_x','getex',l_children,guards,guardMR,**kw)

    def pfey(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('E_y','getey',l_children,guards,guardMR,**kw)

    def pfez(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('E_z','getez',l_children,guards,guardMR,**kw)

    def pfbx(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('B_x','getbx',l_children,guards,guardMR,**kw)

    def pfby(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('B_y','getby',l_children,guards,guardMR,**kw)

    def pfbz(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('B_z','getbz',l_children,guards,guardMR,**kw)

    def pfexg(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('E_x','getexg',l_children,guards,guardMR,**kw)

    def pfeyg(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('E_y','geteyg',l_children,guards,guardMR,**kw)

    def pfezg(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('E_z','getezg',l_children,guards,guardMR,**kw)

    def pfbxg(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('B_x','getbxg',l_children,guards,guardMR,**kw)

    def pfbyg(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('B_y','getbyg',l_children,guards,guardMR,**kw)

    def pfbzg(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('B_z','getbzg',l_children,guards,guardMR,**kw)

    def pfjx(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('J_x','getjx',l_children,guards,guardMR,**kw)

    def pfjy(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('J_y','getjy',l_children,guards,guardMR,**kw)

    def pfjz(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('J_z','getjz',l_children,guards,guardMR,**kw)

    def pfrho(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('Rho','getrho',l_children,guards,guardMR,**kw)

    def pff(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('F','getf',l_children,guards,guardMR,**kw)

    def pfdive(self,l_children=1,guards=0,guardMR=0,**kw):
        self.genericpfem3dMR('Div E','getdive',l_children,guards,guardMR,**kw)

# --- This can only be done after MRBlock3D is defined.
try:
    psyco.bind(MRBlock3D)
except NameError:
    pass
