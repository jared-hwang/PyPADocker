"""
The module supplies functions for dealing with particles.

These commands returns particle info based on selection criteria.
selectparticles(): return list of indices of particles selected
getn(): get number of particles selected
getx(), gety(), getz(), getr(), gettheta(): get particle position
getvx(), getvy(), getvz(): get particle velocity
getux(), getuy(), getuz(): get particle momentum/mass
getxp(), getyp(), getrp(): get tranverse normalized velocities
getgaminv(): get gamma inverse
getpid(): get particle identification number

getxxpslope(), getyypslope(): get x-x' and y-y' slopes

add_particles(): add particles to the simulation

setup_subsets(): Create subsets for particle plots (negative window numbers)
clear_subsets(): Clears the subsets for particle plots (negative window
numbers)
"""
from ..warp import *

# -------------------------------------------------------------------------
def particlesdoc():
    print __doc__

# -------------------------------------------------------------------------
def setpid(pidname, pkg=None, force=False):
    """Sets the specified pid.
    If it has already been set (i.e. it is nonzero), then do nothing.
     - pidname: name of the pid
     - pkg=top: package where pid is defined
     - force=False: If True, set it even if it has already be set
    """
    if pkg is None:
        pkg = top
    if getattr(pkg, pidname) == 0 or force:
        setattr(pkg, pidname, nextpid())

##########################################################################
# Setup the random subsets. This should only be used for special purposes.
# This routine must be explicitly called by the user for this mechanism to
# be operational.
# These routines can be called at any time by the user to add more subsets,
# for example subsets with more or fewer particles, or subsets based on
# the number of particles in a species other than 0.
psubset={}

def setup_subsets(js=0, pgroup=None):
    """
  Use of this function is deprecated and is recommended only for special purposes.
  Adds plotting subset to the list. The size of the subsets is controlled by
  top.npplot.
    - js=0 is the species to create a subset for
    - pgroup=top.pgroup: particle group to base subsets on
    """
    global psubset
    if pgroup is None:
        pgroup = top.pgroup
    if lparallel:
        totalnp = parallelsum(pgroup.nps[js])
        if totalnp == 0:
            totalnp = 1
        fracnp = float(pgroup.nps[js])/float(totalnp)
    else:
        fracnp = 1.
    psubset['n'] = pgroup.nps[js]
    for i, np in enumerate(top.npplot):
        if pgroup.nps[js] == 0:
            psubset[i] = zeros(0, 'i')
            continue
        ntopick = min(1., (np*fracnp)/pgroup.nps[js])
        ii = arange(pgroup.nps[js])
        rr = random.random(pgroup.nps[js])
        # rr = ranf(zeros(pgroup.nps[js]))
        psubset[i] = ii[rr<ntopick].astype('i')

# ----------------------------------------------------------------------------
def clear_subsets():
    "Clears the particle subsets so that they can be updated."
    global psubset
    psubset = {}

# ----------------------------------------------------------------------------
# Modified version of the sample routine from the random module of Python2.3
def populationsample(population, k):
    """Chooses k unique random elements from a population sequence.

    Returns a new list containing elements from the population while
    leaving the original population unchanged.  The resulting list is
    in selection order so that all sub-slices will also be valid random
    samples.  This allows raffle winners (the sample) to be partitioned
    into grand prize and second place winners (the subslices).

    Members of the population need not be hashable or unique.  If the
    population contains repeats, then each occurrence is a possible
    selection in the sample.

    To choose a sample in a range of integers, use range as an argument.
    This is especially fast and space efficient for sampling from a
    large population:   sample(range(10000000), 60)
    """

    # Sampling without replacement entails tracking either potential
    # selections (the pool) in a list or previous selections in a
    # dictionary.

    # When the number of selections is small compared to the population,
    # then tracking selections is efficient, requiring only a small
    # dictionary and an occasional reselection.  For a larger number of
    # selections, the pool tracking method is preferred since the list takes
    # less space than the dictionary and it doesn't suffer from frequent
    # reselections.

    n = population.size
    if k > n:
        return population
    if k < 0:
        raise ValueError("sample size must be non-negative")
    rand = random.random
    _int = int
    # --- Result is an array rather than a list ---
    result = zeros(k, gettypecode(array([population[0]])))
    if n < 6 * k:     # if n len list takes less space than a k len dict
        pool = list(population)
        for i in range(k):         # invariant:  non-selected at [0,n-i)
            j = _int(rand() * (n - i))
            result[i] = pool[j]
            pool[j] = pool[n - i - 1]   # move non-selected item into vacancy
    else:
        try:
            n > 0 and (population[0], population[n//2], population[n-1])
        except (TypeError, KeyError):   # handle sets and dictionaries
            population = tuple(population)
        selected = {}
        for i in range(k):
            j = _int(rand() * n)
            while j in selected:
                j = _int(rand() * n)
            result[i] = selected[j] = population[j]
    return result

##########################################################################
# -------------------------------------------------------------------------
def getattrwithsuffix(object=top, name=None, suffix='', checkmain=1, pkg='top',
                      default=None):
    """
  Helper function modeled on getattr. Like getattr, given an object and a name,
  it will return object.name. Optionally, a suffix can be given, in which case,
  it will return the equivalent of getattr(object,name + suffix).  The object can
  be None, in which case the main dictionary will be searched, returning
  __main__.__dict__[name + suffix]. The object can also be a dictionary or an
  opened PDB file. Finally, if the object is specified but name + suffix is not
  found and if checkmain is true, then the main dictionary will again be
  searched. If nothing is found and default is specified, it is returned,
  otherwise an exception is raised.
    """
    assert name is not None, "A name must be supplied"
    if object is not None:
        # --- This will work if object is a Warp package or a generic instance,
        # --- including a pdb file.
        try:
            result = getattr(object, name + suffix)
            return result
        except (AttributeError, KeyError):
            pass
    if type(object) is DictionaryType:
        try:
            result = object[name + suffix]
            return result
        except KeyError:
            pass
    if type(object) is InstanceType and object.__class__ is PR.PR:
        # --- Try to get data from a pdb file that was written as part of the
        # --- specified package.
        try:
            result = object.read(name + suffix + '@' + pkg)
            return result
        except:
            pass
    if object is None or checkmain:
        import __main__
        try:
            result = __main__.__dict__[name + suffix]
            return result
        except KeyError:
            pass

    # --- If all else fails, try finding the name without the suffix.
    if suffix != '':
        return getattrwithsuffix(object=object, name=name, suffix='',
                                 checkmain=checkmain, pkg=pkg, default=default)

    if default is not None:
        return default
    else:
        raise AttributeError

# -------------------------------------------------------------------------
def _getobjectpgroup(kw):
    suffix = ((kw.get('lost', 0) and 'lost') or '') + kw.get('suffix', '')
    object = kw.get('object', top)
    if suffix == '':
        try:
            pgroup = kw.get('pgroup', object.pgroup)
            if pgroup is None:
                pgroup = object.pgroup
        except (AttributeError, KeyError):
            pgroup = object
    else:
        pgroup = object
    return suffix, object, pgroup

# -------------------------------------------------------------------------
def _setindices(data, ii, i1, i2):
    if len(data) > i2 - i1:
        # --- This does the down select to get the data for the specific
        # --- species, the extent of which is set by i1 and i2.
        # --- The check assumes that if the length of the data array is longer
        # --- than the extent of the species, then data contains all of the
        # --- species and the down select is needed.
        data = data[slice(i1, i2)]
    if ii is None:
        # --- If ii was not already set, set it to have the indices covering
        # --- the full range of the species. No further down selection is needed.
        ii = arange(i1, i2)
    else:
        # --- If ii is a slice, convert it to indices.
        if isinstance(ii, slice):
            ii = arange(ii.start, ii.stop)
        # --- Down select data based on the existing ii. Note that the indexing
        # --- uses ii-i1 since data has already bean sliced using i1 and i2.
        data = data[ii-i1]
    return data, ii

# -------------------------------------------------------------------------
# --- Complete dictionary of possible keywords and their default values
_selectparticles_kwdefaults = {"js":0, "jslist":None, "win":None,
                               "x":None, "y":None, "z":None,
                               "ix":None, "wx":1., "iy":None, "wy":1., "iz":None, "wz":1.,
                               "xl":None, "xu":None, "yl":None, "yu":None, "zl":None, "zu":None,
                               "zc":None, "xc":None, "yc":None,
                               "ssn":None, "ssnid":None, "ii":None,
                               "lost":false, "suffix":'', "object":top, "pgroup":None,
                               "w3dobject":None,
                               'checkargs':0, 'allowbadargs':0}

def selectparticles(iw=0, kwdict={}, **kw):
    """
  Selects particles based on either subsets or windows. By default all particles
  are selected.
  Multiple selection criteria are supported.
    - iw=0: Window to chose from
            For positive values, the z range is determined from the win argument.
            For negative values, a random subset is chosen, with the number
            of particles set by top.npplot[-iw-1].
    - js=0: Species to select particles from
    - jslist=None: List of Species to choose from, e.g. [0,3,4].
                   Alternatively, -1 will include all species.
    - win=top.zwindows + top.zbeam: Windows to use (in lab frame)
    - x=xp: Coordinate to use for x range selection
    - y=yp: Coordinate to use for y range selection
    - z=zp: Coordinate to use for z range selection
    - ix=None: Picks particles within xmesh[ix] +- wx*dx
    - iy=None: Picks particles within ymesh[iy] +- wy*dy
    - iz=None: Picks particles within zbeam + zmesh[iz] +- wz*dz
    - xc=None: Picks particles within xc +- wx*dx
    - yc=None: Picks particles within yc +- wy*dy
    - zc=None: Picks particles within zc +- wz*dz
    - wx=1.: Width of x selection window in grid cells
    - wy=1.: Width of y selection window in grid cells
    - wz=1.: Width of z selection window in grid cells
    - xl=None: Lower range in x of selection region
    - xu=None: Upper range in x of selection region
    - yl=None: Lower range in y of selection region
    - yu=None: Upper range in y of selection region
    - zl=None: Lower range in z of selection region
    - zu=None: Upper range in z of selection region
    - ssn=None: Returns the particle or particles with the specified ssn
                Raises and error if ssn's are not saved - top.ssnpid must be setup
                or the ssnid must be specified.
    - ssnid=top.ssnpid: ID to use when selecting particles with a specified ssn.
    - ii=None: Particle index list
               If supplied, use it for the list of particles instead
               of choosing particles from the given species.
    - lost=false: Flags whether lost particles are returned
                  When true, returns indices to the lost particles rather than
                  the live particles. Note that the lost particle data is
                  saved in different arrays than the live particles.
    - suffix=None: Variable name suffix
                   When specified, python variables with the specified suffix
                   will be used rather than the arrays from top.pgroup.
    - object=top: Object to get particle data from
                  Besides top, this can be an open PDB file, or a dictionary.
    - pgroup=top.pgroup: Particle group to get particles from
    """

    """
  Note, that the following kind of code doesn't work, and I'm not sure of an
  easy way to fix it.
  r = getr(jslist=-1)
  ppxy(jslist=-1, x=r, xu=.1)
  The problem is that the array r is only as long as there are many particles,
  but the ii index will be expecting that r have length npmax.
    """

    # --- Create dictionary of local values and copy it into local dictionary,
    # --- ignoring keywords not listed in _selectparticles_kwdefaults.
    kwvalues = _selectparticles_kwdefaults.copy()
    kwvalues.update(kw)
    kwvalues.update(kwdict)
    # for arg in _selectparticles_kwdefaults: exec arg + " = kwvalues['" + arg + "']"
    # --- This is MUCH faster than the loop, about 100x, but not as nice looking.
    js = kwvalues['js']
    jslist = kwvalues['jslist']
    win = kwvalues['win']
    x = kwvalues['x']
    y = kwvalues['y']
    z = kwvalues['z']
    ix = kwvalues['ix']
    wx = kwvalues['wx']
    iy = kwvalues['iy']
    wy = kwvalues['wy']
    iz = kwvalues['iz']
    wz = kwvalues['wz']
    xl = kwvalues['xl']
    xu = kwvalues['xu']
    yl = kwvalues['yl']
    yu = kwvalues['yu']
    zl = kwvalues['zl']
    zu = kwvalues['zu']
    zc = kwvalues['zc']
    xc = kwvalues['xc']
    yc = kwvalues['yc']
    ssn = kwvalues['ssn']
    ssnid = kwvalues['ssnid']
    ii = kwvalues['ii']
    lost = kwvalues['lost']
    suffix = kwvalues['suffix']
    object = kwvalues['object']
    w3dobject = kwvalues['w3dobject']
    checkargs = kwvalues['checkargs']
    allowbadargs = kwvalues['allowbadargs']

    # --- Check the argument list for bad arguments.
    # --- 'checkargs' allows this routine to be called only to check the
    # --- input for bad arguments.
    # --- 'allowbadargs' allows this routine to be called with bad arguments.
    # --- These are intentionally undocumented features.
    badargs = {}
    for k, v in kwvalues.iteritems():
        if k not in _selectparticles_kwdefaults:
            badargs[k] = v
    if checkargs:
        return badargs
    if badargs and not allowbadargs:
        raise TypeError, "bad argument " + ' '.join(badargs.keys())

    suffix, object, pgroup = _getobjectpgroup(kwvalues)

    # --- If lost is true, then strip off the 'lost' part of the suffix
    # --- for nonparticle quantities
    suffixparticle = suffix
    if lost:
        suffix = suffix[4:]

    ins = getattrwithsuffix(pgroup, 'ins', suffixparticle)
    nps = getattrwithsuffix(pgroup, 'nps', suffixparticle)
    ns = ins.size

    # --- If jslist defined, call selectparticles repeatedly for each species
    # --- on the list.
    # --- Note that if ii is not None, then jslist is ignored since the list
    # --- of particles is already selected and it presumably already contains
    # --- particles from the species in jslist.
    # --- Remove list so selectparticles is subsequently called with one species at a time
    del kwvalues['jslist']
    if jslist is not None and ii is None:
        if jslist == -1:
            jslist = range(0, ns)
        if len(jslist) == 1:
            # --- If jslist has only one species listed, then it is faster
            # --- to act as if js had been set and continue through the rest of
            # --- the function normally.
            js = jslist[0]
        else:
            # --- Create the ii array, which will contain the indices of all of the
            # --- particles in the species in jslist.
            ii = []
            for js in jslist:
                ii += range(ins[js] - 1, ins[js] + nps[js] - 1)
            ii = array(ii, 'l')

    # --- If the w3dobject was not passed in, use w3d, or if the object was
    # --- set to a PDB file, use it.
    if w3dobject is None:
        if object is not top:
            w3dobject = object
        else:
            w3dobject = w3d

    if ii is None:
        i1 = ins[js] - 1
        i2 = ins[js] + nps[js] - 1
    else:
        # --- If ii is given, it already includes the particles from the desired
        # --- species, so no down select to species js is needed.
        i1 = 0
        i2 = getattrwithsuffix(pgroup, 'npmax', suffixparticle)

    if i2 <= i1:
        return array([], 'l')

    if ssn is not None:
        if ssnid is None:
            ssnid = top.ssnpid
        assert ssnid > 0, "ssn's are not used"
        id = getattrwithsuffix(pgroup, 'pid', suffixparticle)[:,ssnid-1]
        id, ii = _setindices(id, ii, i1, i2)
        ii = compress(nint(id)==ssn, ii)

    if xl is not None or xu is not None:
        if x is None:
            x = getattrwithsuffix(pgroup, 'xp', suffixparticle)
        if xl is None:
            xl = -largepos
        if xu is None:
            xu = +largepos
        if xl > xu:
            print "Warning: xl > xu"
        x, ii = _setindices(x, ii, i1, i2)
        ii = compress(logical_and(less(xl, x), less(x, xu)), ii)
    if yl is not None or yu is not None:
        if y is None:
            y = getattrwithsuffix(pgroup, 'yp', suffixparticle)
        if yl is None:
            yl = -largepos
        if yu is None:
            yu = +largepos
        if yl > yu:
            print "Warning: yl > yu"
        y, ii = _setindices(y, ii, i1, i2)
        ii = compress(logical_and(less(yl, y), less(y, yu)), ii)
    if zl is not None or zu is not None:
        if z is None:
            z = getattrwithsuffix(pgroup, 'zp', suffixparticle)
        if zl is None:
            zl = -largepos
        if zu is None:
            zu = +largepos
        if zl > zu:
            print "Warning: zl > zu"
        z, ii = _setindices(z, ii, i1, i2)
        ii = compress(logical_and(less(zl, z), less(z, zu)), ii)

    if ix is not None:
        xmmin = getattrwithsuffix(w3dobject, 'xmmin', suffix, pkg='w3d')
        dx = getattrwithsuffix(w3dobject, 'dx', suffix, pkg='w3d')
        xl = xmmin + ix*dx - wx*dx
        xu = xmmin + ix*dx + wx*dx
        if x is None:
            x = getattrwithsuffix(pgroup, 'xp', suffixparticle)
        x, ii = _setindices(x, ii, i1, i2)
        ii = compress(logical_and(less(xl, x), less(x, xu)), ii)
    if iy is not None:
        ymmin = getattrwithsuffix(w3dobject, 'ymmin', suffix, pkg='w3d')
        dy = getattrwithsuffix(w3dobject, 'dy', suffix, pkg='w3d')
        yl = ymmin + iy*dy - wy*dy
        yu = ymmin + iy*dy + wy*dy
        if y is None:
            y = getattrwithsuffix(pgroup, 'yp', suffixparticle)
        y, ii = _setindices(y, ii, i1, i2)
        ii = compress(logical_and(less(yl, y), less(y, yu)), ii)
    if iz is not None:
        zbeam = getattrwithsuffix(object, 'zbeam', suffix)
        zmmin = getattrwithsuffix(w3dobject, 'zmmin', suffix, pkg='w3d')
        dz = getattrwithsuffix(w3dobject, 'dz', suffix, pkg='w3d')
        zl = zmmin + iz*dz - wz*dz + zbeam
        zu = zmmin + iz*dz + wz*dz + zbeam
        if z is None:
            z = getattrwithsuffix(pgroup, 'zp', suffixparticle)
        z, ii = _setindices(z, ii, i1, i2)
        ii = compress(logical_and(less(zl, z), less(z, zu)), ii)

    if xc is not None:
        if x is None:
            x = getattrwithsuffix(pgroup, 'xp', suffixparticle)
        dx = getattrwithsuffix(w3dobject, 'dx', suffix, pkg='w3d')
        xl = xc - wx*dx
        xu = xc + wx*dx
        x, ii = _setindices(x, ii, i1, i2)
        ii = compress(logical_and(less(xl, x), less(x, xu)), ii)
    if yc is not None:
        if y is None:
            y = getattrwithsuffix(pgroup, 'yp', suffixparticle)
        dy = getattrwithsuffix(w3dobject, 'dy', suffix, pkg='w3d')
        yl = yc - wy*dy
        yu = yc + wy*dy
        y, ii = _setindices(y, ii, i1, i2)
        ii = compress(logical_and(less(yl, y), less(y, yu)), ii)
    if zc is not None:
        if z is None:
            z = getattrwithsuffix(pgroup, 'zp', suffixparticle)
        dz = getattrwithsuffix(w3dobject, 'dz', suffix, pkg='w3d')
        zl = zc - wz*dz
        zu = zc + wz*dz
        z, ii = _setindices(z, ii, i1, i2)
        ii = compress(logical_and(less(zl, z), less(z, zu)), ii)

    if iw < 0:
        if len(psubset) > 0 and psubset['n'] == nps[js]:
            # --- The use of the precalculated subset is not recommended. It is now
            # --- only used if the number of particles is the same as the number
            # --- used when generating the subsets. Also, the routine setup_subsets
            # --- must have been explicitly called by the user for psubset to be
            # --- setup.
            try:
                subset = psubset[-iw-1]
            except KeyError:
                raise IndexError, "Bad window number"
            if ii is None:
                nn = nps[js]
            else:
                nn = ii.size
            subset = compress(less(subset, nn), subset)
            if ii is None:
                ii = i1 + subset
            else:
                ii = ii[subset]
        else:
            # --- This chooses a random sampling of the particles. It will be
            # --- different each time.
            npplot = getattrwithsuffix(object, 'npplot', suffix)
            if ii is None:
                ii = arange(i1, i2)
            ii = populationsample(ii, npplot[-iw-1])
    elif iw == 0:
        if ii is None:
            # --- If all particles are being selected, then a slice object is returned
            # --- so that a reference to the particle array is returned, avoiding
            # --- memory copying.
            ii = slice(i1, i2)
    else:
        zbeam = getattrwithsuffix(object, 'zbeam', suffix)
        zwindows = getattrwithsuffix(object, 'zwindows', suffix)
        if win is None:
            win = zwindows[:,iw] + zbeam
        if win.ndim == 2:
            win = win[:,iw]
        if z is None:
            z = getattrwithsuffix(pgroup, 'zp', suffixparticle)
        z, ii = _setindices(z, ii, i1, i2)
        zslice = z
        ii = ii[logical_and(less(win[0], zslice), less(zslice, win[1]))]

    return ii

# -------------------------------------------------------------------------
# The following return a specific coordinate of the selected particles
# More documetation added after they are declared.
# --- The checks of the len of ii are there in case the particle arrays
# --- are unallocated.
# -------------------------------------------------------------------------
_particlebcastdefault = [1]

def setgetparticlebcastdefault(defval=1):
    """Change the default value of the bcast argument in the particle selection
  routines.
   - defval=1: The default value.
    """
    _particlebcastdefault[0] = defval

# This function handles list of pgroups
# It returns a concatenated list resultlist
# Which contains a collection of attributes of all pgroups
# that are obtained using getfunc
def _handlepgroups(getfunc, iw, gather, bcast, **kw):
    if bcast is None:
        bcast = _particlebcastdefault[0]
    pgroups = kw.pop('pgroups')
    resultlist = []
    for pgroup in pgroups:
        kw['pgroup'] = pgroup
        resultlist.append(getfunc(iw=iw, gather=False, bcast=False, **kw))
    result = concatenate(resultlist)
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
def getn(iw=0, gather=1, bcast=None, **kw):
    "Returns number of particles in selection. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        pgroups = kw.pop('pgroups')
        n = 0
        for pgroup in pgroups:
            kw['pgroup'] = pgroup
            n += getn(iw, gather=False, bcast=False, **kw)
    else:
        if bcast is None:
            bcast = _particlebcastdefault[0]
        ii = selectparticles(iw=iw, kwdict=kw)
        if isinstance(ii, slice):
            n = ii.stop - ii.start
        else:
            n = ii.size
    if lparallel and gather:
        return globalsum(n)
    else:
        return n

# ----------------------------------------------------
# -------------------------------------------------------------------------
def getx(iw=0, gather=1, bcast=None, **kw):
    "Returns the X positions. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getx, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = x[ii]
    elif ii.size > 0:
        result = x[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def gety(iw=0, gather=1, bcast=None, **kw):
    "Returns the Y positions. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(gety, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        y = getattrwithsuffix(pgroup, 'yp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = y[ii]
    elif ii.size > 0:
        result = y[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getz(iw=0, gather=1, bcast=None, **kw):
    "Returns the Z positions. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getz, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        z = getattrwithsuffix(pgroup, 'zp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = z[ii]
    elif ii.size > 0:
        result = z[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getr(iw=0, gather=1, bcast=None, **kw):
    "Returns the R postions. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getr, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = sqrt(x[ii]**2 + y[ii]**2)
    elif ii.size > 0:
        result = sqrt(x[ii]**2 + y[ii]**2)
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def gettheta(iw=0, gather=1, bcast=None, **kw):
    "Returns the theta postions. For particle selection options, see :py:func:`~particles.selectparticles`."
    if bcast is None:
        bcast = _particlebcastdefault[0]
    if 'pgroups' in kw:
        return _handlepgroups(gettheta, iw, gather, bcast, **kw)
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = arctan2(y[ii], x[ii])
    elif ii.size > 0:
        result = arctan2(y[ii], x[ii])
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getvx(iw=0, gather=1, bcast=None, **kw):
    "Returns the X velocity. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getvx, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        ux = getattrwithsuffix(pgroup, 'uxp', suffix)
        gaminv = getattrwithsuffix(pgroup, 'gaminv', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        if top.lrelativ:
            result = ux[ii]*gaminv[ii]
        else:
            result = ux[ii]
    elif ii.size > 0:
        result = ux[ii]*gaminv[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getvy(iw=0, gather=1, bcast=None, **kw):
    "Returns the Y velocity. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getvy, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        uy = getattrwithsuffix(pgroup, 'uyp', suffix)
        gaminv = getattrwithsuffix(pgroup, 'gaminv', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        if top.lrelativ:
            result = uy[ii]*gaminv[ii]
        else:
            result = uy[ii]
    elif ii.size > 0:
        result = uy[ii]*gaminv[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getvz(iw=0, gather=1, bcast=None, **kw):
    "Returns the Z velocity. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getvz, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        uz = getattrwithsuffix(pgroup, 'uzp', suffix)
        gaminv = getattrwithsuffix(pgroup, 'gaminv', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        if top.lrelativ:
            result = uz[ii]*gaminv[ii]
        else:
            result = uz[ii]
    elif ii.size > 0:
        result = uz[ii]*gaminv[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getvr(iw=0, gather=1, bcast=None, **kw):
    "Returns the radial velocity. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getvr, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
        ux = getattrwithsuffix(pgroup, 'uxp', suffix)
        uy = getattrwithsuffix(pgroup, 'uyp', suffix)
        gaminv = getattrwithsuffix(pgroup, 'gaminv', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        tt = arctan2(y[ii], x[ii])
        result = (ux[ii]*cos(tt) + uy[ii]*sin(tt))*gaminv[ii]
    elif ii.size > 0:
        tt = arctan2(y[ii], x[ii])
        result = (ux[ii]*cos(tt) + uy[ii]*sin(tt))*gaminv[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getvtheta(iw=0, gather=1, bcast=None, **kw):
    "Returns the azimuthal velocity. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getvtheta, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
        ux = getattrwithsuffix(pgroup, 'uxp', suffix)
        uy = getattrwithsuffix(pgroup, 'uyp', suffix)
        gaminv = getattrwithsuffix(pgroup, 'gaminv', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        tt = arctan2(y[ii], x[ii])
        result = (-ux[ii]*sin(tt) + uy[ii]*cos(tt))*gaminv[ii]
    elif ii.size > 0:
        tt = arctan2(y[ii], x[ii])
        result = (-ux[ii]*sin(tt) + uy[ii]*cos(tt))*gaminv[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getux(iw=0, gather=1, bcast=None, **kw):
    "Returns the X momentum over mass. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getux, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        ux = getattrwithsuffix(pgroup, 'uxp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = ux[ii]
    elif ii.size > 0:
        result = ux[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getuy(iw=0, gather=1, bcast=None, **kw):
    "Returns the Y momentum over mass. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getuy, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        uy = getattrwithsuffix(pgroup, 'uyp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = uy[ii]
    elif ii.size > 0:
        result = uy[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getuz(iw=0, gather=1, bcast=None, **kw):
    "Returns the Z momentum over mass. For particle selection options, see :py:func:`~particles.selectparticles`. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getuz, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        uz = getattrwithsuffix(pgroup, 'uzp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = uz[ii]
    elif ii.size > 0:
        result = uz[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getxp(iw=0, gather=1, bcast=None, **kw):
    "Returns the X velocity over the Z velocity (X'). For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getxp, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        ux = getattrwithsuffix(pgroup, 'uxp', suffix)
        uz = getattrwithsuffix(pgroup, 'uzp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = ux[ii]/uz[ii]
    elif ii.size > 0:
        result = ux[ii]/uz[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getyp(iw=0, gather=1, bcast=None, **kw):
    "Returns the Y velocity over the Z velocity (Y'). For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getyp, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        uy = getattrwithsuffix(pgroup, 'uyp', suffix)
        uz = getattrwithsuffix(pgroup, 'uzp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = uy[ii]/uz[ii]
    elif ii.size > 0:
        result = uy[ii]/uz[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getrp(iw=0, gather=1, bcast=None, **kw):
    "Returns the radial velocity over the Z velocity (R'). For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getrp, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
        ux = getattrwithsuffix(pgroup, 'uxp', suffix)
        uy = getattrwithsuffix(pgroup, 'uyp', suffix)
        uz = getattrwithsuffix(pgroup, 'uzp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        tt = arctan2(y[ii], x[ii])
        result = ((ux[ii]*cos(tt) + uy[ii]*sin(tt))/uz[ii])
    elif ii.size > 0:
        tt = arctan2(y[ii], x[ii])
        result = ((ux[ii]*cos(tt) + uy[ii]*sin(tt))/uz[ii])
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def gettp(iw=0, gather=1, bcast=None, **kw):
    "Returns the azimuthal velocity over the Z velocity (R'). For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(gettp, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
        ux = getattrwithsuffix(pgroup, 'uxp', suffix)
        uy = getattrwithsuffix(pgroup, 'uyp', suffix)
        uz = getattrwithsuffix(pgroup, 'uzp', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        tt = arctan2(y[ii], x[ii])
        result = ((-ux[ii]*sin(tt) + uy[ii]*cos(tt))/uz[ii])
    elif ii.size > 0:
        tt = arctan2(y[ii], x[ii])
        result = ((-ux[ii]*sin(tt) + uy[ii]*cos(tt))/uz[ii])
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getgaminv(iw=0, gather=1, bcast=None, **kw):
    "Returns the gamma inverse. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getgaminv, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        gaminv = getattrwithsuffix(pgroup, 'gaminv', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = gaminv[ii]
    elif ii.size > 0:
        result = gaminv[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getex(iw=0, gather=1, bcast=None, **kw):
    "Returns the Ex field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getex, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        ex = getattrwithsuffix(pgroup, 'ex', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = ex[ii]
    elif ii.size > 0:
        result = ex[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getey(iw=0, gather=1, bcast=None, **kw):
    "Returns the Ey field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getey, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        ey = getattrwithsuffix(pgroup, 'ey', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = ey[ii]
    elif ii.size > 0:
        result = ey[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getez(iw=0, gather=1, bcast=None, **kw):
    "Returns the Ez field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getez, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        ez = getattrwithsuffix(pgroup, 'ez', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = ez[ii]
    elif ii.size > 0:
        result = ez[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def geter(iw=0, gather=1, bcast=None, **kw):
    "Returns the Er field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(geter, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
        ex = getattrwithsuffix(pgroup, 'ex', suffix)
        ey = getattrwithsuffix(pgroup, 'ey', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        theta = arctan2(y[ii], x[ii])
        result = ex[ii]*cos(theta) + ey[ii]*sin(theta)
    elif ii.size > 0:
        theta = arctan2(y[ii], x[ii])
        result = ex[ii]*cos(theta) + ey[ii]*sin(theta)
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getetheta(iw=0, gather=1, bcast=None, **kw):
    "Returns the Etheta field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getetheta, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
        ex = getattrwithsuffix(pgroup, 'ex', suffix)
        ey = getattrwithsuffix(pgroup, 'ey', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        theta = arctan2(y[ii], x[ii])
        result = -ex[ii]*sin(theta) + ey[ii]*cos(theta)
    elif ii.size > 0:
        theta = arctan2(y[ii], x[ii])
        result = -ex[ii]*sin(theta) + ey[ii]*cos(theta)
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getbx(iw=0, gather=1, bcast=None, **kw):
    "Returns the Bx field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getbx, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        bx = getattrwithsuffix(pgroup, 'bx', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = bx[ii]
    elif ii.size > 0:
        result = bx[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getby(iw=0, gather=1, bcast=None, **kw):
    "Returns the By field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getby, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        by = getattrwithsuffix(pgroup, 'by', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = by[ii]
    elif ii.size > 0:
        result = by[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getbz(iw=0, gather=1, bcast=None, **kw):
    "Returns the Bz field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getbz, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        bz = getattrwithsuffix(pgroup, 'bz', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        result = bz[ii]
    elif ii.size > 0:
        result = bz[ii]
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getbr(iw=0, gather=1, bcast=None, **kw):
    "Returns the Br field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getbr, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
        bx = getattrwithsuffix(pgroup, 'bx', suffix)
        by = getattrwithsuffix(pgroup, 'by', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        theta = arctan2(y[ii], x[ii])
        result = bx[ii]*cos(theta) + by[ii]*sin(theta)
    elif ii.size > 0:
        theta = arctan2(y[ii], x[ii])
        result = bx[ii]*cos(theta) + by[ii]*sin(theta)
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getbtheta(iw=0, gather=1, bcast=None, **kw):
    "Returns the Btheta field applied to the particles. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getbtheta, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    ii = selectparticles(iw=iw, kwdict=kw)
    suffix, object, pgroup = _getobjectpgroup(kw)
    if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
        x = getattrwithsuffix(pgroup, 'xp', suffix)
        y = getattrwithsuffix(pgroup, 'yp', suffix)
        bx = getattrwithsuffix(pgroup, 'bx', suffix)
        by = getattrwithsuffix(pgroup, 'by', suffix)
    if isinstance(ii, slice) and ii.stop > ii.start:
        theta = arctan2(y[ii], x[ii])
        result = -bx[ii]*sin(theta) + by[ii]*cos(theta)
    elif ii.size > 0:
        theta = arctan2(y[ii], x[ii])
        result = -bx[ii]*sin(theta) + by[ii]*cos(theta)
    else:
        result = array([], 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getpid(id=0, iw=0, gather=1, bcast=None, **kw):
    """Returns particle id information.
    - id=0: which pid value to return
            Note that when top.wpid or other id's from the top package are used,
            1 must be subtracted when passed in, i.e. id=top.wpid-1.
            If id=-1, returns all pids.
  For particle selection options, see :py:func:`~particles.selectparticles`.
    """
    if 'pgroups' in kw:
        return _handlepgroups(getpid, iw, gather, bcast, idpgroups=id, **kw)
    if 'idpgroups' in kw:
        id=kw.pop('idpgroups')
    if bcast is None:
        bcast = _particlebcastdefault[0]
    suffix, object, pgroup = _getobjectpgroup(kw)
    lost = kw.get('lost', 0)
    if lost:
        npidlost = getattrwithsuffix(object, 'npidlost')
        dopid = (npidlost > 0)
    else:
        npid = getattrwithsuffix(pgroup, 'npid', suffix)
        dopid = (npid > 0)
    if dopid:
        ii = selectparticles(iw=iw, kwdict=kw)
        if (isinstance(ii, slice) and ii.stop > ii.start) or ii.size > 0:
            pid = getattrwithsuffix(pgroup, 'pid', suffix)
        if isinstance(ii, slice) and ii.stop > ii.start:
            if id >= 0:
                result = pid[ii,id]
            else:
                result = pid[ii,:]
        elif ii.size > 0:
            if id >= 0:
                result = take(pid[:,id], ii)
            else:
                result = take(pid[:,:], ii, 0)
        else:
            if id >= 0:
                result = zeros(0, 'd')
            else:
                result = zeros((0, top.npid), 'd')
    else:
        if id >= 0:
            result = zeros(0, 'd')
        else:
            result = zeros((0, top.npid), 'd')
    if lparallel and gather:
        return gatherarray(result, bcast=bcast)
    else:
        return result

# -------------------------------------------------------------------------
def getvdrifts(iw=0, js=0, jslist=None, gather=1, bcast=None, edrift=1, bdrift=1, **kw):
    "Returns the velocity drifts used in the drift-Lorentz mover. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getvdrifts, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    if jslist is None:
        jslist=[js]
    nptot=0
    vxl=[]
    vyl=[]
    vzl=[]
    suffix, object, pgroup = _getobjectpgroup(kw)
    for js in jslist:
        ii = selectparticles(iw=iw, js=js, jslist=None, kwdict=kw)
        if isinstance(ii, slice) and ii.stop <= ii.start:
            continue
        if not isinstance(ii, slice) and ii.size == 0:
            continue
        if isinstance(ii, slice):
            ii = arange(ii.start, ii.stop)
        np=ii.size
        nptot += np
        x = take(getattrwithsuffix(object, 'xp', suffix), ii)
        y = take(getattrwithsuffix(object, 'yp', suffix), ii)
        z = take(getattrwithsuffix(object, 'zp', suffix), ii)
        ux = take(getattrwithsuffix(object, 'uxp', suffix), ii)
        uy = take(getattrwithsuffix(object, 'uyp', suffix), ii)
        uz = take(getattrwithsuffix(object, 'uzp', suffix), ii)
        npinttmp=w3d.npint + 0
        npfieldtmp=w3d.npfield + 0
        w3d.npint=np
        w3d.npfield=np
        gchange('DKInterptmp')
        setvdrifts(pgroup, np, js + 1, x, y, z, ux, uy, uz, "corrector")
        if edrift:
            if bdrift:
                vxl.append(w3d.vdx[:np])
                vyl.append(w3d.vdy[:np])
                vzl.append(w3d.vdz[:np])
            else:
                vxl.append(w3d.vex[:np])
                vyl.append(w3d.vey[:np])
                vzl.append(w3d.vez[:np])
        else:
            if bdrift:
                vxl.append(w3d.vbx[:np])
                vyl.append(w3d.vby[:np])
                vzl.append(w3d.vbz[:np])
        w3d.npint=npinttmp
        w3d.npfield=npfieldtmp
        gchange('DKInterptmp')
    if nptot == 0:
        return array([], 'd'), array([], 'd'), array([], 'd')
    else:
        vx = zeros(nptot, 'd')
        vy = zeros(nptot, 'd')
        vz = zeros(nptot, 'd')
        ip=0
        for il in range(len(vxl)):
            np=shape(vxl[il])
            vx[ip:ip+np]=vxl[il]
            vy[ip:ip+np]=vyl[il]
            vz[ip:ip+np]=vzl[il]
            ip += np
    if lparallel and gather:
        return (gatherarray(vx, bcast=bcast),
                gatherarray(vy, bcast=bcast),
                gatherarray(vz, bcast=bcast))
    else:
        return vx, vy, vz

# -------------------------------------------------------------------------
# def getw(iw=0, gather=1, bcast=None, **kw):
#     """Returns particle weights.
#   This is equivalent to getpid(id=top.wpid-1). For particle selection options, see :py:func:`~particles.selectparticles`.
#     """
#     return getpid(id=top.wpid-1, gather=gather, bcast=bcast, pgroups=pgroups, **kw)
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
def getw(iw=0, gather=1, bcast=None, **kw):
    """Returns particle weights.
  This is equivalent to getpid(id=top.wpid-1). For particle selection options, see :py:func:`~particles.selectparticles`.
    """
    return getpid(id=top.wpid-1, gather=gather, bcast=bcast, **kw)

def getke(iw=0, js=0, jslist=None, gather=1, bcast=None, **kw):
    "Returns the particles kinetic energy. For particle selection options, see :py:func:`~particles.selectparticles`."
    if 'pgroups' in kw:
        return _handlepgroups(getke, iw, gather, bcast, **kw)
    if bcast is None:
        bcast = _particlebcastdefault[0]
    if jslist is None:
        jslist=[js]
    nptot=0
    kel=[]
    suffix, object, pgroup = _getobjectpgroup(kw)
    masses = getattrwithsuffix(pgroup, 'sm', suffix)
    for js in jslist:
        ii = selectparticles(iw=iw, js=js, jslist=None, kwdict=kw)
        if isinstance(ii, slice) and ii.stop <= ii.start:
            continue
        if not isinstance(ii, slice) and ii.size == 0:
            continue
        if isinstance(ii, slice):
            ii = arange(ii.start, ii.stop)
        np=ii.size
        nptot += np
        mass = masses[js]
        if top.lrelativ:
            gamma = 1./take(getattrwithsuffix(pgroup, 'gaminv', suffix), ii)
            kes = (gamma - 1.)*mass*clight**2/echarge
        else:
            vx = take(getattrwithsuffix(pgroup, 'uxp', suffix), ii)
            vy = take(getattrwithsuffix(pgroup, 'uyp', suffix), ii)
            vz = take(getattrwithsuffix(pgroup, 'uzp', suffix), ii)
            kes = 0.5*mass*(vx**2 + vy**2 + vz**2)/echarge
        kel.append(kes[:np])
    if nptot == 0:
        ke = array([], 'd')
    else:
        ke = zeros(nptot, 'd')
        ip=0
        for il in range(len(kel)):
            np=shape(kel[il])[0]
            ke[ip:ip+np]=kel[il]
            ip += np
    if lparallel and gather:
        return gatherarray(ke, bcast=bcast)
    else:
        return ke

# -------------------------------------------------------------------------
# Add the some parallel documentation to each of the routines.
_paralleldoc = """

  - gather=true: Flags whether results are gathered to processor 0
              If true, all processors must make the call.
              What is returned by the other processors is determined by the
              bcast argument.
              When false, each processor does the calculation locally, only
              considers particles within its domain, and returns the
              local value.
              The argument is ignored in serial.
  - bcast=1: Flag whether results are broadcast to all processors
             Only used when gather is true. When true, the result is broadcast
             from processor 0 to all other processors. The return value is the
             same on all processors. Otherwise, only processor 0 will have the
             correct result - the return values on other processors should not
             be used.
             Note that the setgetparticlebcastdefault can be used to change
             the default value to 0, which is more efficient if the bcast is
             not needed.
             The argument is ignored in serial.
"""
getn.__doc__ += _paralleldoc
getx.__doc__ += _paralleldoc
gety.__doc__ += _paralleldoc
getz.__doc__ += _paralleldoc
getr.__doc__ += _paralleldoc
gettheta.__doc__ += _paralleldoc
getvx.__doc__ += _paralleldoc
getvy.__doc__ += _paralleldoc
getvz.__doc__ += _paralleldoc
getvr.__doc__ += _paralleldoc
getvtheta.__doc__ += _paralleldoc
getux.__doc__ += _paralleldoc
getuy.__doc__ += _paralleldoc
getuz.__doc__ += _paralleldoc
getxp.__doc__ += _paralleldoc
getyp.__doc__ += _paralleldoc
getrp.__doc__ += _paralleldoc
gettp.__doc__ += _paralleldoc
getgaminv.__doc__ += _paralleldoc
getex.__doc__ += _paralleldoc
getey.__doc__ += _paralleldoc
getez.__doc__ += _paralleldoc
geter.__doc__ += _paralleldoc
getetheta.__doc__ += _paralleldoc
getbx.__doc__ += _paralleldoc
getby.__doc__ += _paralleldoc
getbz.__doc__ += _paralleldoc
getbr.__doc__ += _paralleldoc
getbtheta.__doc__ += _paralleldoc
getpid.__doc__ += _paralleldoc
getvdrifts.__doc__ += _paralleldoc
getw.__doc__ += _paralleldoc
del _paralleldoc
# -------------------------------------------------------------------------

##########################################################################
def getxxpslope(iw=0, iz=-1, kwdict={}, checkargs=0):
    """
  Calculates the x-x' slope from particle moments.
  This returns a tuple containing (slope, xoffset, xpoffset, vz).
  The product slope*vz gives the slope for x-vx.
    - iw=0: Window number
    - iz=-1: Z grid location
    - zl=None: When specified, lower range of selection region
    - zu=None: When specified, upper range of selection region
    - zc=None: When specified, center of range of selection region
    - slopejs=-1: Species whose moments are used to calculate the slope
                  -1 means use data combined from all species.
    - js, jslist: can also specify a species or list of species.
                  slopejs takes precedence, then jslist, then js.
    """
    if checkargs:
        # --- This is only needed since no other routines take the slopejs
        # --- argument.
        kwdictcopy = kwdict.copy()
        if 'slopejs' in kwdictcopy:
            del kwdictcopy['slopejs']
        return kwdictcopy
    zl = kwdict.get('zl')
    zu = kwdict.get('zu')
    zc = kwdict.get('zc')
    js = kwdict.get('js', -1)
    jslist = kwdict.get('jslist', [js])
    slopejs = kwdict.get('slopejs', None)
    if slopejs is not None:
        jslist = [slopejs]
    suffix = kwdict.get('suffix', '')
    object = kwdict.get('object', top)
    zwindows = getattrwithsuffix(object, 'zwindows', suffix)
    if zl is not None and zu is not None:
        zc = 0.5*(zl + zu)
    if iw > 0:
        zc = 0.5*(zwindows[0,iw] + zwindows[1,iw])
    if zc is not None:
        zmmntmin = getattrwithsuffix(object, 'zmmntmin', suffix)
        dzm = getattrwithsuffix(object, 'dzm', suffix)
        iz = nint(floor((zc - zmmntmin)/dzm))
        wz =            (zc - zmmntmin)/dzm - iz
    else:
        wz = 0.
    nzmmnt = getattrwithsuffix(object, 'nzmmnt', suffix)
    nsum = 0.
    xxpbarsum = 0.
    xbarsum = 0.
    xpbarsum = 0.
    xsqbarsum = 0.
    vzbarsum = 0.
    for js in jslist:
        if iz is not None and 0 <= iz <= nzmmnt:
            izp1 = iz + 1
            if iz == nzmmnt:
                izp1 = iz
            nz = getattrwithsuffix(object, 'pnumz', suffix)
            xxpbarz = getattrwithsuffix(object, 'xxpbarz', suffix)
            xbarz = getattrwithsuffix(object, 'xbarz', suffix)
            xpbarz = getattrwithsuffix(object, 'xpbarz', suffix)
            xsqbarz = getattrwithsuffix(object, 'xsqbarz', suffix)
            vzbarz = getattrwithsuffix(object, 'vzbarz', suffix)
            n = nz[iz,js]*(1. - wz) + nz[izp1,js]*wz
            xxpbar = xxpbarz[iz,js]*(1. - wz) + xxpbarz[izp1,js]*wz
            xbar = xbarz[iz,js]*(1. - wz) + xbarz[izp1,js]*wz
            xpbar = xpbarz[iz,js]*(1. - wz) + xpbarz[izp1,js]*wz
            xsqbar = xsqbarz[iz,js]*(1. - wz) + xsqbarz[izp1,js]*wz
            vzbar = vzbarz[iz,js]*(1. - wz) + vzbarz[izp1,js]*wz
        elif iw <= 0:
            n = getattrwithsuffix(object, 'pnum', suffix)[0,js]
            xxpbar = getattrwithsuffix(object, 'xxpbar', suffix)[0,js]
            xbar = getattrwithsuffix(object, 'xbar', suffix)[0,js]
            xpbar = getattrwithsuffix(object, 'xpbar', suffix)[0,js]
            xsqbar = getattrwithsuffix(object, 'xsqbar', suffix)[0,js]
            vzbar = getattrwithsuffix(object, 'vzbar', suffix)[0,js]
        else:
            n = 0.
        nsum += n
        xxpbarsum += xxpbar*n
        xbarsum += xbar*n
        xpbarsum += xpbar*n
        xsqbarsum += xsqbar*n
        vzbarsum += vzbar*n

    if nsum > 0.:
        xxpbar = xxpbarsum/nsum
        xbar = xbarsum/nsum
        xpbar = xpbarsum/nsum
        xsqbar = xsqbarsum/nsum
        vzbar = vzbarsum/nsum
        xrms = sqrt(xsqbar - xbar**2)
        slope = (xxpbar - xbar*xpbar)/dvnz(xrms)**2
        xoffset = xbar
        xpoffset = xpbar
        vz = vzbar
    else:
        slope = 0.
        xoffset = 0.
        xpoffset = 0.
        vz = 0.
    return (slope, xoffset, xpoffset, vz)

# -------------------------------------------------------------------------
def getyypslope(iw=0, iz=-1, kwdict={}, checkargs=0):
    """
  Calculates the y-y' slope from particle moments.
  This returns a tuple containing (slope, yoffset, ypoffset, vz).
  The product slope*vz gives the slope for y-vy.
    - iw=0: Window number
    - iz=-1: Z grid location
    - zl=None: When specified, lower range of selection region
    - zu=None: When specified, upper range of selection region
    - zc=None: When specified, center of range of selection region
    - slopejs=-1: Species whose moments are used to calculate the slope
                  -1 means use data combined from all species.
    - js, jslist: can also specify a species or list of species.
                  slopejs takes precedence, then jslist, then js.
    """
    if checkargs:
        # --- This is only needed since no other routines take the slopejs
        # --- argument.
        kwdictcopy = kwdict.copy()
        if 'slopejs' in kwdictcopy:
            del kwdictcopy['slopejs']
        return kwdictcopy
    zl = kwdict.get('zl')
    zu = kwdict.get('zu')
    zc = kwdict.get('zc')
    js = kwdict.get('js', -1)
    jslist = kwdict.get('jslist', [js])
    slopejs = kwdict.get('slopejs', None)
    if slopejs is not None:
        jslist = [slopejs]
    suffix = kwdict.get('suffix', '')
    object = kwdict.get('object', top)
    zwindows = getattrwithsuffix(object, 'zwindows', suffix)
    if zl is not None and zu is not None:
        zc = 0.5*(zl + zu)
    if iw > 0:
        zc = 0.5*(zwindows[0,iw] + zwindows[1,iw])
    if zc is not None:
        zmmntmin = getattrwithsuffix(object, 'zmmntmin', suffix)
        dzm = getattrwithsuffix(object, 'dzm', suffix)
        iz = nint(floor((zc - zmmntmin)/dzm))
        wz =            (zc - zmmntmin)/dzm - iz
    else:
        wz = 0.
    nzmmnt = getattrwithsuffix(object, 'nzmmnt', suffix)
    nsum = 0.
    yypbarsum = 0.
    ybarsum = 0.
    ypbarsum = 0.
    ysqbarsum = 0.
    vzbarsum = 0.
    for js in jslist:
        if iz is not None and 0 <= iz <= nzmmnt:
            izp1 = iz + 1
            if iz == nzmmnt:
                izp1 = iz
            nz = getattrwithsuffix(object, 'pnumz', suffix)
            yypbarz = getattrwithsuffix(object, 'yypbarz', suffix)
            ybarz = getattrwithsuffix(object, 'ybarz', suffix)
            ypbarz = getattrwithsuffix(object, 'ypbarz', suffix)
            ysqbarz = getattrwithsuffix(object, 'ysqbarz', suffix)
            vzbarz = getattrwithsuffix(object, 'vzbarz', suffix)
            n = nz[iz,js]*(1. - wz) + nz[izp1,js]*wz
            yypbar = yypbarz[iz,js]*(1. - wz) + yypbarz[izp1,js]*wz
            ybar = ybarz[iz,js]*(1. - wz) + ybarz[izp1,js]*wz
            ypbar = ypbarz[iz,js]*(1. - wz) + ypbarz[izp1,js]*wz
            ysqbar = ysqbarz[iz,js]*(1. - wz) + ysqbarz[izp1,js]*wz
            vzbar = vzbarz[iz,js]*(1. - wz) + vzbarz[izp1,js]*wz
        elif iw <= 0:
            n = getattrwithsuffix(object, 'pnum', suffix)[0,js]
            yypbar = getattrwithsuffix(object, 'yypbar', suffix)[0,js]
            ybar = getattrwithsuffix(object, 'ybar', suffix)[0,js]
            ypbar = getattrwithsuffix(object, 'ypbar', suffix)[0,js]
            ysqbar = getattrwithsuffix(object, 'ysqbar', suffix)[0,js]
            vzbar = getattrwithsuffix(object, 'vzbar', suffix)[0,js]
        else:
            n = 0.
        nsum += n
        yypbarsum += yypbar*n
        ybarsum += ybar*n
        ypbarsum += ypbar*n
        ysqbarsum += ysqbar*n
        vzbarsum += vzbar*n

    if nsum > 0.:
        yypbar = yypbarsum/nsum
        ybar = ybarsum/nsum
        ypbar = ypbarsum/nsum
        ysqbar = ysqbarsum/nsum
        vzbar = vzbarsum/nsum
        yrms = sqrt(ysqbar - ybar**2)
        slope = (yypbar - ybar*ypbar)/dvnz(yrms)**2
        yoffset = ybar
        ypoffset = ypbar
        vz = vzbar
    else:
        slope = 0.
        yoffset = 0.
        ypoffset = 0.
        vz = 0.
    return (slope, yoffset, ypoffset, vz)

# -------------------------------------------------------------------------
def getvzrange(kwdict={}):
    "Returns a tuple containing the Vz range for plots"
    suffix = kwdict.get('suffix', '')
    object = kwdict.get('object', top)
    vzrng = getattrwithsuffix(object, 'vzrng', suffix)
    if (vzrng != 0.):
        vtilt = getattrwithsuffix(object, 'vtilt', suffix)
        vbeam = getattrwithsuffix(object, 'vbeam', suffix)
        vzshift = getattrwithsuffix(object, 'vzshift', suffix)
        vzmax = (1. + vtilt)*vbeam*(1. + vzrng) - vzshift
        vzmin = (1. - vtilt)*vbeam*(1. - vzrng) - vzshift
    else:
        js = kwdict.get('js', -1)
        slopejs = kwdict.get('slopejs', js)
        vzminp = getattrwithsuffix(object, 'vzminp', suffix)
        vzmaxp = getattrwithsuffix(object, 'vzmaxp', suffix)
        vzmax = vzmaxp[slopejs] + 0.1*(vzmaxp[slopejs] - vzminp[slopejs])
        vzmin = vzminp[slopejs] - 0.1*(vzmaxp[slopejs] - vzminp[slopejs])
    return (vzmin, vzmax)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
def add_particles(x=0., y=0., z=0., vx=0., vy=0., vz=0., gi=1.,
                  pid=None, w=None, ux=None, uy=None, uz=None,
                  js=0, species_number=None, sid=None,
                  lallindomain=None, unique_particles=None,
                  xmmin=None, xmmax=None,
                  ymmin=None, ymmax=None,
                  zmmin=None, zmmax=None,
                  lmomentum=false,
                  l2symtry=None, l4symtry=None, lrz=None,
                  resetrho=false, dofieldsol=false, resetmoments=false,
                  pgroup=None,
                  ex=0., ey=0., ez=0., bx=0., by=0., bz=0.,
                  pidpairs=None,
                  lfields=false, lnewparticles=true, lusespaceabove=true,
                  lreturndata=false):
    """
  Adds particles to the simulation
    - x,y,z,vx,vy,vz,gi: particle coordinates and velocities.
                         Can be arrays or scalars. Scalars are broadcast to all
                         particles. Any that are unsupplied default to zero,
                         except gi, which defaults to 1. (gi is 1/gamma, the
                         relatistic paramter)
    - ux, uy, uz: Instead of vx, vy, vz, these can be used, passing in momentums.
                  The flag lmomentum will be set to true.
    - pid: additional particle information, such as an ID number or weight.
    - pidpairs=None: Allows setting specific pid columns. Argument must be a
                     list of lists, each having the format [id,pidvalue].
                     id is the one-based index returned by nextpid.
                     The assigment pid[:,id-1] = pidvalue is done.
    - w=1.: particle weight
            this is only used if top.wpid > 0 and if lnewparticles is true.
    - js=0, or species_number=0: species to which new particles will be added
    - lallindomain=false: Flags whether particles are within the parallel domains
                          When true, all particles are assumed to be with in the
                          extent of the domain so particle scraping is not done.
                          This is automatically set to true when the code is in
                          slice mode, i.e. after package('wxy'). Except if the
                          option is explicitly set. If false, the particles that
                          are outside the parallel domain are not added.
    - unique_particles: Specifies whether all of the particles passed in are unique.
                        If this is set, it overides lallindomain, and sets
                        lallindomain=unique_particles. If this is False, then only
                        particles within the parallel subdomains are added, otherwise
                        all particles are added.
    - xmmin=top.xpminlocal,xmmax=top.xpmaxlocal
                       x extent of the domain - should only be set in unusual
                       circumstances.
    - ymmin=top.ypminlocal,ymmax=top.ypmaxlocal
                       y extent of the domain - should only be set in unusual
                       circumstances.
    - zmmin=top.zpminlocal + top.zgrid,zmmax=top.zpmaxlocal + top.zgrid:
                       z extent of the domain - should only be set in unusual
                       circumstances.
    - l2symtry,l4symtry,lrz: System symmetries, default to w3d values
    - lmomentum=false: Flags whether momentum or velocities are being input.
                       Set to false when velocities are input as velocities, true
                       when input as massless momentum (as WARP stores them).
                       Only used when top.lrelativ is true.
    - lnewparticles=true: Flag whether the particles are treated as new.
                          If so, the ssn will be set if needed, and the
                          position saved as the birth location.
                          Set this to false if using add_particles to restore
                          particles.
    - lusespaceabove=true: Flag where to put the particles in memory.
                           When true, the new particles are preferentially
                           placed in memory in the space above the existing
                           particles, otherwise below. This is only needed for
                           special cases.
    - lreturndata=false: When true, the data from the newly created particles
                         is returned. This will include only particles that
                         were actually added (i.e. does not include any
                         particles outside the domain when lallindomain is
                         false). Also, if lmomentum if false and top.lrelativ
                         is true, the velocites will have been converted into
                         massless momentum. If lfields is true, the fields are
                         also returned. The data is returned in the order
                         x, y, z, vx, vy, vz, gi, pid, ex, ey, ez, bx, by, bz.
    - pgroup=top.pgroup: Particle group to add particles too
    """

    if pgroup is None:
        pgroup = top.pgroup

    # --- Keep track if pid was input. If not, set it to the default value.
    # --- This is needed to fix a possible conflict with the w input argument
    pid_unset = (pid is None)
    if pid is None:
        pid = 0.

    # --- Check if this is a new species
    if species_number is not None:
        js = species_number
    if js + 1 > top.ns:
        setnspecies(js + 1, pgroup)

    # --- Set the sid if it hasn't already be done
    if sid is None:
        sid = js
    if pgroup.sid[js] == -1:
        pgroup.sid[js] = sid

    # --- Use momentum quantities if specified. These take
    # --- precedence over vx, vy, and vz.
    if ux is not None or uy is not None or uz is not None:
        lmomentum = true
        if ux is not None:
            vx = ux
        if uy is not None:
            vy = uy
        if uz is not None:
            vz = uz

    # --- Get length of arrays, set to one for scalars
    lenx = size(x)
    leny = size(y)
    lenz = size(z)
    lenvx = size(vx)
    lenvy = size(vy)
    lenvz = size(vz)
    lengi = size(gi)
    lenex = size(ex)
    leney = size(ey)
    lenez = size(ez)
    lenbx = size(bx)
    lenby = size(by)
    lenbz = size(bz)
    try:
        lenpid = len(pid)
    except TypeError:
        lenpid = 1

    # --- If any of the inputs are arrays that are zero length, then return
    if (lenx == 0 or leny == 0 or lenz == 0 or lenvx == 0 or lenvy == 0 or lenvz == 0 or
        lengi == 0 or lenpid == 0 or
        lenex == 0 or leney == 0 or lenez == 0 or lenbx == 0 or lenby == 0 or lenbz == 0):
        return

    # --- Max length of input arrays
    maxlen = max(lenx, leny, lenz, lenvx, lenvy, lenvz, lengi, lenpid, lenex, leney, lenez, lenbx, lenby, lenbz)
    assert lenx==maxlen or lenx==1, "Length of x doesn't match len of others"
    assert leny==maxlen or leny==1, "Length of y doesn't match len of others"
    assert lenz==maxlen or lenz==1, "Length of z doesn't match len of others"
    assert lenvx==maxlen or lenvx==1, "Length of vx doesn't match len of others"
    assert lenvy==maxlen or lenvy==1, "Length of vy doesn't match len of others"
    assert lenvz==maxlen or lenvz==1, "Length of vz doesn't match len of others"
    assert lengi==maxlen or lengi==1, "Length of gi doesn't match len of others"
    assert lenex==maxlen or lenex==1, "Length of ex doesn't match len of others"
    assert leney==maxlen or leney==1, "Length of ey doesn't match len of others"
    assert lenez==maxlen or lenez==1, "Length of ez doesn't match len of others"
    assert lenbx==maxlen or lenbx==1, "Length of bx doesn't match len of others"
    assert lenby==maxlen or lenby==1, "Length of by doesn't match len of others"
    assert lenbz==maxlen or lenbz==1, "Length of bz doesn't match len of others"
    assert lenpid==maxlen or lenpid==1, "Length of pid doesn't match len of others"

    # --- Convert all to arrays of length maxlen, broadcasting scalars
    x = array(x)*ones(maxlen, 'd')
    y = array(y)*ones(maxlen, 'd')
    z = array(z)*ones(maxlen, 'd')
    vx = array(vx)*ones(maxlen, 'd')
    vy = array(vy)*ones(maxlen, 'd')
    vz = array(vz)*ones(maxlen, 'd')
    gi = array(gi)*ones(maxlen, 'd')
    ex = array(ex)*ones(maxlen, 'd')
    ey = array(ey)*ones(maxlen, 'd')
    ez = array(ez)*ones(maxlen, 'd')
    bx = array(bx)*ones(maxlen, 'd')
    by = array(by)*ones(maxlen, 'd')
    bz = array(bz)*ones(maxlen, 'd')
    pid = array(pid)*ones([maxlen, top.npid], 'd')

    if lnewparticles:
        # --- Set time of creation
        if top.tbirthpid > 0:
            pid[:,top.tbirthpid-1] = top.time

        # --- Set weights
        # --- The weights can be set either by passing in a pid array or by the w argument.
        # --- The w argument takes precedence. If neither are supplied, a default value
        # --- of 1 is used.
        if pid_unset and w is None:
            w = 1.
        if w is not None and top.wpid > 0:
            pid[:,top.wpid-1] = array(w)*ones(maxlen, 'd')
        # --- Note that ssn is set in addpart

        # --- Set xyz old
        if top.xoldpid > 0:
            pid[:,top.xoldpid-1] = x
        if top.yoldpid > 0:
            pid[:,top.yoldpid-1] = y
        if top.zoldpid > 0:
            pid[:,top.zoldpid-1] = z

    # --- Load in any pid data passed in
    if pidpairs is not None:
        for id, pp in pidpairs:
            pid[:,id-1] = pp

    # --- Set extent of domain
    if xmmin is None:
        xmmin = top.xpminlocal
    if xmmax is None:
        xmmax = top.xpmaxlocal
    if w3d.solvergeom == w3d.XZgeom:
        if ymmin is None:
            ymmin = -largepos
        if ymmax is None:
            ymmax = +largepos
    else:
        if ymmin is None:
            ymmin = top.ypminlocal
        if ymmax is None:
            ymmax = top.ypmaxlocal
    if zmmin is None:
        zmmin = top.zpminlocal + top.zgrid
    if zmmax is None:
        zmmax = top.zpmaxlocal + top.zgrid

    if l2symtry is None:
        l2symtry = w3d.l2symtry
    if l4symtry is None:
        l4symtry = w3d.l4symtry
    if lrz is None:
        lrz = (w3d.solvergeom in [w3d.RZgeom, w3d.Rgeom])

    if unique_particles is not None:
        lallindomain = unique_particles

    # --- When running in slice mode, automatically set lallindomain to true.
    # --- It is assumed that all particles will be within the specified domain,
    # --- since in the slice mode, the z of the particles is ignored.
    # --- The user can still set lallindomain to false to override this.
    if lallindomain is None:
        if getcurrpkg() == 'wxy':
            lallindomain = true
        else:
            lallindomain = false

    # --- Do some error checking
    if not lallindomain and zmmin == zmmax:
        print "=================================================================="
        print "Add_particles: warning - no particles will be loaded - you should"
        print "either set lallindomain=true or set zmmin and zmmax so they are"
        print "different from each other."
        print "=================================================================="

    # --- Get the number of particles before adding the new ones.
    if lreturndata:
        nbefore = pgroup.nps[js]

    # --- Now data can be passed into the fortran add_particles routine.
    addpart(pgroup, maxlen, top.npid, x, y, z, vx, vy, vz, gi, ex, ey, ez, bx, by, bz, pid, js+1,
            lallindomain, xmmin, xmmax, ymmin, ymmax, zmmin, zmmax,
            l2symtry, l4symtry, lrz,
            lmomentum, lfields, lnewparticles, lusespaceabove)

    # --- Get the number of particles after adding the new ones.
    if lreturndata:
        nafter = pgroup.nps[js]

    # --- If the slice code is active, then call initdtp
    if package()[0] == 'wxy':
        initdtp(pgroup)

    # --- Do followup work if requested
    if resetrho:
        loadrho()
    if dofieldsol:
        fieldsol(-1)
    if resetmoments:
        from ..diagnostics import getzmom
        getzmom.zmmnt()
        if top.it%top.nhist == 0:
            top.jhist = top.jhist - 1
            savehist(top.time)

    if lreturndata:
        if lallindomain:
            if lfields:
                return (x, y, z, vx, vy, vz, gi, pid, ex, ey, ez, bx, by, bz)
            else:
                return (x, y, z, vx, vy, vz, gi, pid)
        else:
            # --- In this case, the particles outside the domain have been removed
            # --- and new data compressed to the bottom of the arrays.
            nn = nafter - nbefore
            if lfields:
                return (x[:nn], y[:nn], z[:nn], vx[:nn], vy[:nn], vz[:nn], gi[:nn], pid[:nn,:], ex[:nn], ey[:nn], ez[:nn], bx[:nn], by[:nn], bz[:nn])
            else:
                return (x[:nn], y[:nn], z[:nn], vx[:nn], vy[:nn], vz[:nn], gi[:nn], pid[:nn,:])


addparticles = add_particles
