# File MPLOT.PY --- standard post-processing for module-impedance runs

from ..warp import *


### MPLOT - setup plots
def mplot(dumpfile):
    restore(dumpfile)
    # --- to set cgm file name, default marks=0, run hcpon
    setup(prefix=dumpfile[0:4]+"h")
    winon(dpi=75)

### MOUNTAINPLOT1 - Mountain-range plots of quantities saved vs. z at
### every timestep
# A much faster version of the original mountainplot1. The speed is
# gained in the averaging by two ways.  First, the loop over time (2nd
# dimension) is vectorized, creating an array with the qty at all times
# instead of one at a time. The second speed up is minimizing the work
# needed for the sum. Each ave is the same as the neighboring average
# except for the addition and subtraction of one point. Taking advantage
# of that minimizes the number of additions.
def mountainplot1(qtyname,qty,kwdict={},**kw):
    """
  Mountain-range plots of quantities saved vs. z at every timestep
    - qtyname: String used in axis labels and title
    - qty: Data to be plotted
    - ifdelta=1: Turns on subtraction from initial value
    - ifordt=0: Makes ordinate in time (otherwise units of qty)
    - ifvst=1: Makes abscissa in time (otherwise z)
    - ifneg=0: Plots negative of data (or delta)
    - nlines=100: Number of lines to overlay
    - istart=0: Time index at which to start the plots
    - iend=shape(qty)[1]-1: Time index at which to end the plots
    - istep=jhist/nlines: Can be specified instead of nlines
    - navg=0: Turns on averaging (over 2*navg+1 points)
    - offset=0: Abscissa offset of overlaid lines
    - ordoffset=0: Ordinate offset of overlaid lines. Can be an array such as
                   top.hzbeam.
    - color='fg': Line color
    - nz=shape(qty)[0]-1: Number of z points
    - dz=w3d.dz: Size of z points
    - zmmin=w3d.zmmin: Start of z points
    - hvbeam=top.hvbeam: Velocity used for horizontal offsets
    - titles=1: When 0, no titles are plotted
    - titlet=None: Top title
    - titleb=None: Bottom title
    - titlel=None: Left title
    - titler=None: Right title
    """
    kwdefaults = {'ord':None,'ifdelta':0,'ifordt':0,'ifvst':1,'ifneg':0,
                  'nlines':100,'istart':0,'iend':None,'istep':None,
                  'navg':0,'offset':0,'ordoffset':0,
                  'color':'fg','nz':None,
                  'dz':w3d.dz,'zmmin':w3d.zmmin,'hvbeam':top.hvbeam,'titles':1,
                  'titlet':None,'titleb':None,'titlel':None,'titler':None}
    kwvalues = kwdefaults.copy()
    kwvalues.update(kw)
    kwvalues.update(kwdict)
    for arg in kwdefaults: exec(arg+" = kwvalues['"+arg+"']")
    badargs = checkarguments(kwvalues,kwdefaults)
    if badargs: raise Exception("bad argument ",' '.join(badargs.keys()))

    # --- Special arguments
    if iend is None: iend = shape(qty)[1] - 1
    if nz is None: nz = shape(qty)[0] - 1

    if ord is None:
        if ifvst and hvbeam[0] != 0.:
            ord = dz/hvbeam[0]*(iota(nz+2,2,-1)-2)
        else:
            ord = dz*iota(0,nz)+zmmin
    if istep is None: istep = max(1,(iend-istart)/nlines)
    sign = 1
    shift = offset
    titlet = qtyname
    titlel = qtyname + " (offset for later curves)"
    if ifdelta:
        titlet = "delta " + titlet
        titlel = "delta " + titlel
    if ifneg:
        sign = -1
        titlet = "- " + titlet
        titlel = "- " + titlel
    titleb = "z (beam frame)"
    if ifvst and hvbeam[0] != 0.:
        titleb = "time (relative to beam arrival at gap)"
    if ifordt:
        titlel = "time       "
        shift = top.dt * top.nhist
    if offset != 0.:
        abscissascale = shift/offset
    else:
        abscissascale = 1.
    titler="nlines = %d  navg = %d  offset = %6.2e" % (nlines,navg,offset)
    ptitles(titlet,titleb,titlel,titler)
    if navg:
        hl = qty[:,istart:iend+1:istep] + 0.
        hl[navg,:] = ave(qty[navg-navg:navg+navg+1,::istep])
        for j in range(navg+1,nz-navg-1):
            hl[j,:] = hl[j-1,:] + (qty[j+navg,::istep] -
                                   qty[j-navg-1,::istep])/(2*navg+1)
    else:
        hl = qty[:,istart:iend+1:istep]
    if ifdelta:
        hl0 = hl[:,0]
    else:
        hl0 = zeros(shape(hl[:,0]),'d')
    if len(shape(ordoffset)) == 0:
        ordoffset = iota(0,iend)*ordoffset
    else:
        assert len(ordoffset) == shape(qty)[1],\
               "ordoffset must have same length is data's second dimension"
    j = -1
    for i in range(istart,iend+1,istep):
        j = j + 1
        plg(i*shift+abscissascale*sign*(hl[:,j]-hl0),ord+ordoffset[i],color=color)


### PMLCHG - plot line charge
def pmlchg(ifdelta=1,ifordt=0,ifneg=0,nlines=100,ifvst=1,navg=0,offset=5.e-9,
          ordoffset=0.,istep=None,title="Line charge"):
    mountainplot1(title,top.hlinechg,ifdelta=ifdelta,ifordt=ifordt,
                  ifneg=ifneg,ifvst=ifvst,nlines=nlines,navg=navg,
                  offset=offset,ordoffset=ordoffset,istep=istep)

### PMCURR - plot current
def pmcurr(ifdelta=1,ifordt=0,ifneg=0,nlines=100,ifvst=1,navg=0,offset=1.e-3,
          ordoffset=0.,istep=None,title="Current"):
    mountainplot1(title,top.hcurrz[...,-1],ifdelta=ifdelta,ifordt=ifordt,
                  ifneg=ifneg,ifvst=ifvst,nlines=nlines,navg=navg,
                  offset=offset,ordoffset=ordoffset,istep=istep)

### PMVZOFZ - plot axial velocity
def pmvzofz(ifdelta=1,ifordt=0,ifneg=0,nlines=100,ifvst=1,navg=0,offset=5.e-9,
            ordoffset=0.,istep=None,title="Axial velocity"):
    mountainplot1(title,top.hvzofz,ifdelta=ifdelta,ifordt=ifordt,
                  ifneg=ifneg,ifvst=ifvst,nlines=nlines,navg=navg,
                  offset=offset,ordoffset=ordoffset,istep=istep)

### PVGAP - Mountain-range plot of vgap, which is saved at every 10th accel gap
def pvgap(_hvgap=None,ifzt=0,ifneg=0,nincr=1,offset=5000,color='fg'):
    if not _hvgap: _hvgap = hvgap
    assert top.hvbeam[0] != 0.,'top.hvbeam[0] must be nonzero'
    loi = (0.5+(top.zlatstrt+.5*(top.acclzs[::10]+top.acclze[::10])-top.zzmax)/
           (top.hvbeam[0]*top.dt)).astype("i")
    hii = (0.5+(top.zlatstrt+.5*(top.acclzs[::10]+top.acclze[::10])-top.zzmin)/
           (top.hvbeam[0]*top.dt)).astype("i")
    abscissa=iota(0,hii[0]-loi[0]-1)*top.dt
    labscissa = len(abscissa)
    start = 0
    sign = 1
    shift = offset
    titlet =  "Gap Voltage"
    titlel = "Gap Voltage (offset for later curves)"
    if ifneg:
        sign = -1
        titlet = "- " + titlet
        titlel = "- " + titlel
    titleb = "time (relative to beam arrival at gap)"
    if ifzt:
        titlel = "time (of beam arrival at gap)"
        titleb = "z (relative to beam head)"
        abscissa = top.hvbeam[0]*abscissa
        shift = top.dt * steps_p_perd / 2 * 10
        start = (top.zlatstrt+top.acclzs[0]-top.zzmax) / top.hvbeam[0]
    titler="nincr = %d  offset = %6.2e" % (nincr,offset)
    ptitles(titlet,titleb,titlel,titler)
    imax = 0
    for i in range(0,shape(_hvgap)[0]):
        if max(abs(_hvgap[i,:])):
            imax = i
    for i in range(0,imax+1,nincr):
        plg(start+i*shift+shift/offset*sign*_hvgap[i,loi[i]:loi[i]+labscissa],
            abscissa,color=color)
