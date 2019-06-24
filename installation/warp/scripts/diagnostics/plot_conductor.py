from ..warp import *
import __main__
import copy


def plot_conductordoc():
    print """
  The following functions plot contours of the potential in various planes
  along with the conductors in that plane. The first three plot with the axis
  having units of meters. The suffix 'g' means that it plots with the axis
  having units of number of grid cells. The 'box' suffix means
  that it plots a box around each grid point inside of a conductor.

  pfxy, pfzx, pfzy, pfzr
  pfxyg, pfzxg, pfzyg, pfzrg
  pfxybox, pfzxbox, pfzybox, pfzxboxi, pfzyboxi

  plotgrid: plots the x-z mesh in the lab frame (including any bends)
  pfzxlab: makes the pfzx plot in the lab frame (including any bends)
  plotsrfrv: handy command to plot r versus z for a suface of revolution, giving
             the function describing it

  plotquadoutline: plots outline of quadrupole structure for quad elements
  plotheleoutline: plots outline of quadrupole structure for hele elements
  plotemltoutline: plots outline of quadrupole structure for emlt elements
  plotmmltoutline: plots outline of pipe structure for mmlt elements
  plotpgrdoutline: plots outline of quadrupole structure for pgrd elements
  plotaccloutline: plots outline of pipe structure for accl elements
  plotdrftoutline: plots outline of pipe structure for drft elements

  visualizeconductors: create 3-D visualization of conductors based on subgrid
                       data

  cleanconductors: not a plot routine, buts removes conductor points not
                   within the the range of the field solve
    """

######################################################################
# functions to plot the conductor points and subgrid data            #
# in MKS units                                                       #
######################################################################

# --- Convenience function to plot the sub-grid data
def plotcond(iy,ix,iz,izp,numb,ymin,xmin,dy,dx,color,mglevel,yscale,xscale,
             conductors,local):
    if conductors is None: return
    interior = conductors.interior
    nn = interior.n
    if nn > 0:
        sx = ['x','y','z'][ix]
        sy = ['x','y','z'][iy]
        sz = ['x','y','z'][iz]
        lx = getattr(conductors,'levell'+sx)[mglevel]
        ly = getattr(conductors,'levell'+sy)[mglevel]
        lz = getattr(conductors,'levell'+sz)[mglevel]
        ixc = interior.indx[ix,:]*lx
        iyc = interior.indx[iy,:]*ly
        izc = interior.indx[iz,:]*lz
        if ix == 0: ixc = ixc + conductors.levelix[mglevel]*lx
        if iy == 0: iyc = iyc + conductors.levelix[mglevel]*ly
        if iz == 0: izc = izc + conductors.levelix[mglevel]*lz
        if ix == 1: ixc = ixc + conductors.leveliy[mglevel]*lx
        if iy == 1: iyc = iyc + conductors.leveliy[mglevel]*ly
        if iz == 1: izc = izc + conductors.leveliy[mglevel]*lz
        if ix == 2: ixc = ixc + conductors.leveliz[mglevel]*lx
        if iy == 2: iyc = iyc + conductors.leveliz[mglevel]*ly
        if iz == 2: izc = izc + conductors.leveliz[mglevel]*lz
        cnumb = interior.numb
        level = equal(mglevel,interior.ilevel)
    else:
        ixc = array([])
        iyc = array([])
        izc = array([])
        cnumb = array([])
        level = array([])
    if w3d.solvergeom == w3d.XYgeom and iz == 2:
        izequal = 1
    else:
        izequal = equal(izc[:nn],izp)
    ii = compress(logical_and(izequal,level[:nn]),arange(nn))
    if numb is not None:
        cnumb = take(cnumb,ii)
        ii = compress(equal(cnumb,numb),ii)
    xx = (take(ixc,ii)*dx+xmin)*xscale
    yy = (take(iyc,ii)*dy+ymin)*yscale
    plp(yy,xx,color=color,local=local)

def plotsubgrid(iy,ix,iz,pp,izp,numb,ymin,xmin,dy,dx,color,subgridlen,mglevel,
                yscale,xscale,inverted,conductors,local):
    if conductors is None: return
    evensubgrid = conductors.evensubgrid
    oddsubgrid = conductors.oddsubgrid
    pp = ((pp + conductors.leveliz[mglevel]) % 2)
    subgrid = [evensubgrid,oddsubgrid][pp]
    nn = subgrid.n
    if nn > 0:
        sx = ['x','y','z'][ix]
        sy = ['x','y','z'][iy]
        sz = ['x','y','z'][iz]
        lx = getattr(conductors,'levell'+sx)[mglevel]
        ly = getattr(conductors,'levell'+sy)[mglevel]
        lz = getattr(conductors,'levell'+sz)[mglevel]
        ixc = subgrid.indx[ix,:]
        iyc = subgrid.indx[iy,:]
        izc = subgrid.indx[iz,:]*lz
        if ix == 0: ixc = ixc + conductors.levelix[mglevel]
        if iy == 0: iyc = iyc + conductors.levelix[mglevel]
        if iz == 0: izc = izc + conductors.levelix[mglevel]
        if ix == 1: ixc = ixc + conductors.leveliy[mglevel]
        if iy == 1: iyc = iyc + conductors.leveliy[mglevel]
        if iz == 1: izc = izc + conductors.leveliy[mglevel]
        if ix == 2: ixc = ixc + conductors.leveliz[mglevel]
        if iy == 2: iyc = iyc + conductors.leveliz[mglevel]
        if iz == 2: izc = izc + conductors.leveliz[mglevel]
        delmx = abs(subgrid.dels[2*ix  ])
        delpx = abs(subgrid.dels[2*ix+1])
        delmy = abs(subgrid.dels[2*iy  ])
        delpy = abs(subgrid.dels[2*iy+1])
        numbmx = subgrid.numb[2*ix  ]
        numbpx = subgrid.numb[2*ix+1]
        numbmy = subgrid.numb[2*iy  ]
        numbpy = subgrid.numb[2*iy+1]
        level = equal(mglevel,subgrid.ilevel)
    else:
        lx = array([1])
        ly = array([1])
        lz = array([1])
        ixc = array([])
        iyc = array([])
        izc = array([])
        delmx = array([])
        delpx = array([])
        delmy = array([])
        delpy = array([])
        numbmx = array([])
        numbpx = array([])
        numbmy = array([])
        numbpy = array([])
        level = array([])
    if w3d.solvergeom == w3d.XYgeom and iz == 2:
        izequal = 1
    else:
        izequal = equal(izc[:nn],izp)
    ii = compress(logical_and(izequal,equal(level[:nn],1)),arange(nn))
    dx = dx*lx
    dy = dy*ly
    xx = (take(ixc,ii)*dx+xmin)*xscale
    yy = (take(iyc,ii)*dy+ymin)*yscale
    delmx = take(delmx,ii)*dx*xscale
    delpx = take(delpx,ii)*dx*xscale
    delmy = take(delmy,ii)*dy*yscale
    delpy = take(delpy,ii)*dy*yscale
    if numb is not None:
        numbmx = take(numbmx,ii)
        numbpx = take(numbpx,ii)
        numbmy = take(numbmy,ii)
        numbpy = take(numbpy,ii)
        xxmx = compress(numbmx == numb,xx)
        xxpx = compress(numbpx == numb,xx)
        xxmy = compress(numbmy == numb,xx)
        xxpy = compress(numbpy == numb,xx)
        yymx = compress(numbmx == numb,yy)
        yypx = compress(numbpx == numb,yy)
        yymy = compress(numbmy == numb,yy)
        yypy = compress(numbpy == numb,yy)
        delmx = compress(numbmx == numb,delmx)
        delpx = compress(numbmx == numb,delpx)
        delmy = compress(numbmy == numb,delmy)
        delpy = compress(numbmy == numb,delpy)
    else:
        xxmx = xx
        xxpx = xx
        xxmy = xx
        xxpy = xx
        yymx = yy
        yypx = yy
        yymy = yy
        yypy = yy
    # --- This code combines all of the individual lines into the list pp.
    # --- This vectorized code avoids slower explicit loops.
    pp = []
    dx = dx*xscale
    dy = dy*yscale
    if inverted:
        sgmy = lambda x,y,d,dg:[y,y-d,x,x]
        sgmx = lambda x,y,d,dg:[y,y,x,x-d]
        sgpy = lambda x,y,d,dg:[y,y+d,x,x]
        sgpx = lambda x,y,d,dg:[y,y,x,x+d]
    else:
        sgmy = lambda x,y,d,dg:[y-dg*ceil(abs(d)),y-d,x,x]
        sgmx = lambda x,y,d,dg:[y,y,x-dg*ceil(abs(d)),x-d]
        sgpy = lambda x,y,d,dg:[y+dg*ceil(abs(d)),y+d,x,x]
        sgpx = lambda x,y,d,dg:[y,y,x+dg*ceil(abs(d)),x+d]
    for mapfunc,xx,yy,dl,dd in [(sgmy,xxmy,yymy,delmy,dy),
                                (sgmx,xxmx,yymx,delmx,dx),
                                (sgpy,xxpy,yypy,delpy,dy),
                                (sgpx,xxpx,yypx,delpx,dx)]:
        ii = compress(less(abs(dl),abs(dd)*subgridlen),arange(len(xx)))
        xx = take(xx,ii)
        yy = take(yy,ii)
        dl = take(dl,ii)
        pp = pp + map(mapfunc,xx,yy,dl,dd*ones(len(xx)))
    # --- Convert the list to an array and plot it
    # --- If the length is zero, create an empty array that conforms
    if len(pp) == 0: pp = transpose(array([[],[],[],[]]))
    pp = array(pp)
    pldj(pp[:,2],pp[:,0],pp[:,3],pp[:,1],color=color,local=local)

def plotcondfill(iy,ix,iz,izp,ymin,xmin,dy,dx,mglevel,yscale,xscale,
                 conductors,local):
    """
  Plots conductors, filling them in with a solid color. The color is given
  by the conductor number.
    """
    if conductors is None: return
    interior = conductors.interior
    evensubgrid = conductors.evensubgrid
    oddsubgrid = conductors.oddsubgrid
    sx = ['x','y','z'][ix]
    sy = ['x','y','z'][iy]
    sz = ['x','y','z'][iz]
    # --- Get the numbers of conductor points
    nc = interior.n
    ne = evensubgrid.n
    no = oddsubgrid.n
    ns = ne + no
    lx = getattr(conductors,'levell'+sx)[mglevel]
    ly = getattr(conductors,'levell'+sy)[mglevel]
    lz = getattr(conductors,'levell'+sz)[mglevel]
    if nc > 0:
        ixc = interior.indx[ix,:nc]
        iyc = interior.indx[iy,:nc]
        izc = interior.indx[iz,:nc]
        numb = interior.numb[:nc]
        # --- Add z offset of data. This only applies for the parallel version
        if ix == 0: ixc = ixc + conductors.levelix[mglevel]
        if iy == 0: iyc = iyc + conductors.levelix[mglevel]
        if iz == 0: izc = izc + conductors.levelix[mglevel]
        if ix == 1: ixc = ixc + conductors.leveliy[mglevel]
        if iy == 1: iyc = iyc + conductors.leveliy[mglevel]
        if iz == 1: izc = izc + conductors.leveliy[mglevel]
        if ix == 2: ixc = ixc + conductors.leveliz[mglevel]
        if iy == 2: iyc = iyc + conductors.leveliz[mglevel]
        if iz == 2: izc = izc + conductors.leveliz[mglevel]
        levelc = equal(mglevel,interior.ilevel[:nc])
    else:
        ixc = array([])
        iyc = array([])
        izc = array([])
        numb = array([])
        levelc = array([])
    if ne > 0:
        iexs = evensubgrid.indx[ix,:ne]
        ieys = evensubgrid.indx[iy,:ne]
        iezs = evensubgrid.indx[iz,:ne]
        ecdelmx = abs(evensubgrid.dels[2*ix  ,:ne])
        ecdelpx = abs(evensubgrid.dels[2*ix+1,:ne])
        ecdelmy = abs(evensubgrid.dels[2*iy  ,:ne])
        ecdelpy = abs(evensubgrid.dels[2*iy+1,:ne])
        enumbmx = evensubgrid.numb[2*ix  ,:ne]
        enumbpx = evensubgrid.numb[2*ix+1,:ne]
        enumbmy = evensubgrid.numb[2*iy  ,:ne]
        enumbpy = evensubgrid.numb[2*iy+1,:ne]
        elevels = equal(mglevel,evensubgrid.ilevel[:ne])
    else:
        iexs = array([])
        ieys = array([])
        iezs = array([])
        ecdelmx = array([])
        ecdelpx = array([])
        ecdelmy = array([])
        ecdelpy = array([])
        enumbmx = array([])
        enumbpx = array([])
        enumbmy = array([])
        enumbpy = array([])
        elevels = array([])
    if no > 0:
        ioxs = oddsubgrid.indx[ix,:no]
        ioys = oddsubgrid.indx[iy,:no]
        iozs = oddsubgrid.indx[iz,:no]
        ocdelmx = abs(oddsubgrid.dels[2*ix  ,:no])
        ocdelpx = abs(oddsubgrid.dels[2*ix+1,:no])
        ocdelmy = abs(oddsubgrid.dels[2*iy  ,:no])
        ocdelpy = abs(oddsubgrid.dels[2*iy+1,:no])
        onumbmx = oddsubgrid.numb[2*ix  ,:no]
        onumbpx = oddsubgrid.numb[2*ix+1,:no]
        onumbmy = oddsubgrid.numb[2*iy  ,:no]
        onumbpy = oddsubgrid.numb[2*iy+1,:no]
        olevels = equal(mglevel,oddsubgrid.ilevel[:no])
    else:
        ioxs = array([])
        ioys = array([])
        iozs = array([])
        ocdelmx = array([])
        ocdelpx = array([])
        ocdelmy = array([])
        ocdelpy = array([])
        onumbmx = array([])
        onumbpx = array([])
        onumbmy = array([])
        onumbpy = array([])
        olevels = array([])

    def arrayjoin(a,b):
        r = zeros(len(a)+len(b),gettypecode(a))
        r[:len(a)] = a
        r[len(a):] = b
        return r

    # --- The even and odd data are merged into the same list.
    ixs = arrayjoin(iexs,ioxs)
    iys = arrayjoin(ieys,ioys)
    izs = arrayjoin(iezs,iozs)*lz
    # --- Add z offset of data. This only applies for the parallel version
    if ix == 0: ixs = ixs + conductors.levelix[mglevel]
    if iy == 0: iys = iys + conductors.levelix[mglevel]
    if iz == 0: izs = izs + conductors.levelix[mglevel]
    if ix == 1: ixs = ixs + conductors.leveliy[mglevel]
    if iy == 1: iys = iys + conductors.leveliy[mglevel]
    if iz == 1: izs = izs + conductors.leveliy[mglevel]
    if ix == 2: ixs = ixs + conductors.leveliz[mglevel]
    if iy == 2: iys = iys + conductors.leveliz[mglevel]
    if iz == 2: izs = izs + conductors.leveliz[mglevel]

    delmx = arrayjoin(ecdelmx,ocdelmx)
    delpx = arrayjoin(ecdelpx,ocdelpx)
    delmy = arrayjoin(ecdelmy,ocdelmy)
    delpy = arrayjoin(ecdelpy,ocdelpy)
    numbmx = arrayjoin(enumbmx,onumbmx)
    numbpx = arrayjoin(enumbpx,onumbpx)
    numbmy = arrayjoin(enumbmy,onumbmy)
    numbpy = arrayjoin(enumbpy,onumbpy)
    levels = arrayjoin(elevels,olevels)

    # --- Select out the conductor points in the appropriate slice and in
    # --- the appropriate refinement level.
    if w3d.solvergeom == w3d.XYgeom and iz == 2:
        izcequal = 1
        izsequal = 1
    else:
        izcequal = equal(izc,izp)
        izsequal = equal(izs,izp)
    iic = compress(logical_and(izcequal,equal(levelc,1)),arange(nc))
    iis = compress(logical_and(izsequal,equal(levels,1)),arange(ns))
    dx = dx*lx*xscale
    dy = dy*ly*yscale
    ixc = take(ixc,iic)
    iyc = take(iyc,iic)
    xxc = ixc*dx+xmin*xscale
    yyc = iyc*dy+ymin*yscale
    numb = take(numb,iic)
    ixs = take(ixs,iis)
    iys = take(iys,iis)
    xxs = ixs*dx+xmin
    yys = iys*dy+ymin
    delmx = take(delmx,iis)*dx
    delpx = take(delpx,iis)*dx
    delmy = take(delmy,iis)*dy
    delpy = take(delpy,iis)*dy
    numbmx = take(numbmx,iis)
    numbpx = take(numbpx,iis)
    numbmy = take(numbmy,iis)
    numbpy = take(numbpy,iis)
# # --- Now, after the global gather, if there is no data then quit.
# if len(ixc) + len(ixs) == 0: return
    # --- Get max grid point so that an array can be created which covers all
    # --- of the data.
    if len(ixc) == 0:
        maxixc = 0
        maxiyc = 0
    else:
        maxixc = max(ixc)
        maxiyc = max(iyc)
    if len(ixs) == 0:
        maxixs = 0
        maxiys = 0
    else:
        maxixs = max(ixs)
        maxiys = max(iys)
    nx = max(maxixc,maxixs) + 1
    ny = max(maxiyc,maxiys) + 1
    iii = zeros((5,1+nx,1+ny),'l') - 1
    mx,px,my,py = 1,2,3,4
    # --- Flag grid points where the conductors are.
    for i in range(len(ixc)): iii[0,ixc[i],iyc[i]] = i
    for i in range(len(ixs)):
        if abs(delmx[i]) < abs(dx): iii[mx,ixs[i],iys[i]] = i
        if abs(delpx[i]) < abs(dx): iii[px,ixs[i],iys[i]] = i
        if abs(delmy[i]) < abs(dy): iii[my,ixs[i],iys[i]] = i
        if abs(delpy[i]) < abs(dy): iii[py,ixs[i],iys[i]] = i
    x = []
    y = []
    z = []
    n = []

    # --- Loop over all conductor points, drawing a fill polygon for each.
    ixall = list(ixc) + list(ixs)
    iyall = list(iyc) + list(iys)
    for ix,iy in zip(ixall,iyall):
        # --- For each point, deal with the grid cell up and to the right.
        # --- This is done so that no two grid points draw over the same spot.
        i1 = _corner1(ix  ,iy  ,px,py,+delpx,+delpy,numb,numbpx,numbpy,0,
                      iii,xxc,yyc,xxs,yys,x,y,z,n)
        i2 = _corner1(ix  ,iy+1,px,my,+delpx,-delmy,numb,numbpx,numbmy,1,
                      iii,xxc,yyc,xxs,yys,x,y,z,n)
        i3 = _corner1(ix+1,iy+1,mx,my,-delmx,-delmy,numb,numbmx,numbmy,0,
                      iii,xxc,yyc,xxs,yys,x,y,z,n)
        i4 = _corner1(ix+1,iy  ,mx,py,-delmx,+delpy,numb,numbmx,numbpy,1,
                      iii,xxc,yyc,xxs,yys,x,y,z,n)
        if i1[0] + i2[0] + i3[0] + i4[0] > 0:
            n.append(i1[0] + i2[0] + i3[0] + i4[0])
            if   i1[1] is not None: z.append(i1[1])
            elif i2[1] is not None: z.append(i2[1])
            elif i3[1] is not None: z.append(i3[1])
            elif i4[1] is not None: z.append(i4[1])
        # --- One special case needs to be dealt with. If there is no grid
        # --- point in the lower left corner of a grid cell, it would otherwise
        # --- be skipped.
        if (ix > 0 and iii[py,ix,iy] >= 0 and max(iii[:,ix-1,iy]) == -1 and
                       iii[px,ix-1,iy+1] >= 0 and iii[0,ix,iy+1] >= 0):
            i2 = _corner1(ix-1,iy+1,px,my,+delpx,-delmy,numb,numbpx,numbmy,1,
                          iii,xxc,yyc,xxs,yys,x,y,z,n)
            i3 = _corner1(ix  ,iy+1,mx,my,-delmx,-delmy,numb,numbmx,numbmy,0,
                          iii,xxc,yyc,xxs,yys,x,y,z,n)
            i4 = _corner1(ix  ,iy  ,mx,py,-delmx,+delpy,numb,numbmx,numbpy,1,
                          iii,xxc,yyc,xxs,yys,x,y,z,n)
            if i2[0] + i3[0] + i4[0] > 0:
                n.append(i2[0] + i3[0] + i4[0])
                if   i2[1] is not None: z.append(i2[1])
                elif i3[1] is not None: z.append(i3[1])
                elif i4[1] is not None: z.append(i4[1])

    # --- Now that the data is gathered, make the plot.
    z = array(z).astype(ubyte)
    n = array(n).astype('l')
    plfp(z,y,x,n,local=local)

def _corner1(ix,iy,sx,sy,delsx,delsy,nm,nmsx,nmsy,iparity,
             iii,xxc,yyc,xxs,yys,x,y,z,n):
    """
  Special routine which gathers the data to draw a polygon in the local grid
  cell. It must be called for each of the four corners of a grid cell.
    """
    nn = 0
    numb = None
    if iii[0,ix,iy] >= 0:
        nn = nn + 1
        i0 = iii[0,ix,iy]
        x.append(xxc[i0])
        y.append(yyc[i0])
        numb = nm[i0]
    for i in range(2):
        if i == iparity:
            if iii[sx,ix,iy] >= 0:
                nn = nn + 1
                i0 = iii[sx,ix,iy]
                x.append(xxs[i0]+delsx[i0])
                y.append(yys[i0])
                numb = nmsx[i0]
        else:
            if iii[sy,ix,iy] >= 0:
                nn = nn + 1
                i0 = iii[sy,ix,iy]
                x.append(xxs[i0])
                y.append(yys[i0]+delsy[i0])
                numb = nmsy[i0]
    return nn,numb

def plotcondfillnew(yy,xx,zz,iz,ymin,xmin,dy,dx,mglevel,yscale,xscale,
                    f3dcond,f3dmg,local):
    """
  Plots conductors, filling them in with a solid color. The color is given
  by the conductor number.
    """
    if f3dcond is None: return
    # --- Get the numbers of conductor points
    nc = f3dcond.ncond
    ne = f3dcond.necndbdy
    no = f3dcond.nocndbdy
    ns = ne + no
    lx = getattr(f3dmg,'mglevelsl'+xx)[mglevel]
    ly = getattr(f3dmg,'mglevelsl'+yy)[mglevel]
    lz = getattr(f3dmg,'mglevelsl'+zz)[mglevel]
    if nc > 0:
        ixc = getattr(f3dcond,'i'+xx+'cond')[:nc]
        iyc = getattr(f3dcond,'i'+yy+'cond')[:nc]
        izc = getattr(f3dcond,'i'+zz+'cond')[:nc]
        numb = f3dcond.condnumb[:nc]
        # --- Add z offset of data. This only applies for the parallel version
        if xx == 0: ixc = ixc + f3dmg.mglevelsix[mglevel]
        if yy == 0: iyc = iyc + f3dmg.mglevelsix[mglevel]
        if zz == 0: izc = izc + f3dmg.mglevelsix[mglevel]
        if xx == 1: ixc = ixc + f3dmg.mglevelsiy[mglevel]
        if yy == 1: iyc = iyc + f3dmg.mglevelsiy[mglevel]
        if zz == 1: izc = izc + f3dmg.mglevelsiy[mglevel]
        if xx == 2: ixc = ixc + f3dmg.mglevelsiz[mglevel]
        if yy == 2: iyc = iyc + f3dmg.mglevelsiz[mglevel]
        if zz == 2: izc = izc + f3dmg.mglevelsiz[mglevel]
        try:
            levc = f3dcond.icondlevel[:nc]
            levelc = equal(mglevel,levc)
        except AttributeError:
            levelc = ones(nc,'l')
    else:
        ixc = array([])
        iyc = array([])
        izc = array([])
        numb = array([])
        levelc = array([])
    if ne > 0:
        iexs = getattr(f3dcond,'iecnd'+xx)[:ne]
        ieys = getattr(f3dcond,'iecnd'+yy)[:ne]
        iezs = getattr(f3dcond,'iecnd'+zz)[:ne]
        ecdelmx = getattr(f3dcond,'ecdelm'+xx)[:ne]
        ecdelpx = getattr(f3dcond,'ecdelp'+xx)[:ne]
        ecdelmy = getattr(f3dcond,'ecdelm'+yy)[:ne]
        ecdelpy = getattr(f3dcond,'ecdelp'+yy)[:ne]
        enumbmx = getattr(f3dcond,'ecnumbm'+xx)[:ne]
        enumbpx = getattr(f3dcond,'ecnumbp'+xx)[:ne]
        enumbmy = getattr(f3dcond,'ecnumbm'+yy)[:ne]
        enumbpy = getattr(f3dcond,'ecnumbp'+yy)[:ne]
        try:
            elevs = f3dcond.iecndlevel[:ne]
            elevels = equal(mglevel,elevs)
        except AttributeError:
            elevels = ones(ne,'l')
    else:
        iexs = array([])
        ieys = array([])
        iezs = array([])
        ecdelmx = array([])
        ecdelpx = array([])
        ecdelmy = array([])
        ecdelpy = array([])
        enumbmx = array([])
        enumbpx = array([])
        enumbmy = array([])
        enumbpy = array([])
        elevels = array([])
    if no > 0:
        ioxs = getattr(f3dcond,'iocnd'+xx)[:no]
        ioys = getattr(f3dcond,'iocnd'+yy)[:no]
        iozs = getattr(f3dcond,'iocnd'+zz)[:no]
        ocdelmx = getattr(f3dcond,'ocdelm'+xx)[:no]
        ocdelpx = getattr(f3dcond,'ocdelp'+xx)[:no]
        ocdelmy = getattr(f3dcond,'ocdelm'+yy)[:no]
        ocdelpy = getattr(f3dcond,'ocdelp'+yy)[:no]
        onumbmx = getattr(f3dcond,'ocnumbm'+xx)[:no]
        onumbpx = getattr(f3dcond,'ocnumbp'+xx)[:no]
        onumbmy = getattr(f3dcond,'ocnumbm'+yy)[:no]
        onumbpy = getattr(f3dcond,'ocnumbp'+yy)[:no]
        try:
            olevs = f3dcond.iocndlevel[:no]
            olevels = equal(mglevel,olevs)
        except AttributeError:
            olevels = ones(no,'l')
    else:
        ioxs = array([])
        ioys = array([])
        iozs = array([])
        ocdelmx = array([])
        ocdelpx = array([])
        ocdelmy = array([])
        ocdelpy = array([])
        onumbmx = array([])
        onumbpx = array([])
        onumbmy = array([])
        onumbpy = array([])
        olevels = array([])

    # --- The even and odd data are merged into the same list.
    ixs = array(list(iexs) + list(ioxs))
    iys = array(list(ieys) + list(ioys))
    izs = array(list(iezs) + list(iozs))*lz
    # --- Add z offset of data. This only applies for the parallel version
    if xx == 0: ixs = ixs + f3dmg.mglevelsix[mglevel]
    if yy == 0: iys = iys + f3dmg.mglevelsix[mglevel]
    if zz == 0: izs = izs + f3dmg.mglevelsix[mglevel]
    if xx == 1: ixs = ixs + f3dmg.mglevelsiy[mglevel]
    if yy == 1: iys = iys + f3dmg.mglevelsiy[mglevel]
    if zz == 1: izs = izs + f3dmg.mglevelsiy[mglevel]
    if xx == 2: ixs = ixs + f3dmg.mglevelsiz[mglevel]
    if yy == 2: iys = iys + f3dmg.mglevelsiz[mglevel]
    if zz == 2: izs = izs + f3dmg.mglevelsiz[mglevel]

    delmx = array(list(ecdelmx) + list(ocdelmx))
    delpx = array(list(ecdelpx) + list(ocdelpx))
    delmy = array(list(ecdelmy) + list(ocdelmy))
    delpy = array(list(ecdelpy) + list(ocdelpy))
    numbmx = array(list(enumbmx) + list(onumbmx))
    numbpx = array(list(enumbpx) + list(onumbpx))
    numbmy = array(list(enumbmy) + list(onumbmy))
    numbpy = array(list(enumbpy) + list(onumbpy))
    levels = array(list(elevels) + list(olevels))

    # --- Select out the conductor points in the appropriate slice and in
    # --- the appropriate refinement level.
    if w3d.solvergeom == w3d.XYgeom and iz == 2:
        izcequal = 1
        izsequal = 1
    else:
        izcequal = equal(izc,izp)
        izsequal = equal(izs,izp)
    iic = compress(logical_and(izcequal,equal(levelc,1)),arange(nc))
    iis = compress(logical_and(izsequal,equal(levels,1)),arange(ns))
    dx = dx*lx*xscale
    dy = dy*ly*yscale
    ixc = take(ixc,iic)
    iyc = take(iyc,iic)
    xxc = ixc*dx+xmin*xscale
    yyc = iyc*dy+ymin*yscale
    numb = take(numb,iic)
    ixs = take(ixs,iis)
    iys = take(iys,iis)
    xxs = ixs*dx+xmin
    yys = iys*dy+ymin
    delmx = take(delmx,iis)*dx
    delpx = take(delpx,iis)*dx
    delmy = take(delmy,iis)*dy
    delpy = take(delpy,iis)*dy
    numbmx = take(numbmx,iis)
    numbpx = take(numbpx,iis)
    numbmy = take(numbmy,iis)
    numbpy = take(numbpy,iis)
# # --- Now, after the global gather, if there is no data then quit.
# if len(ixc) + len(ixs) == 0: return
    # --- Get max grid point so that an array can be created which covers all
    # --- of the data.
    if len(ixc) == 0:
        maxixc = 0
        maxiyc = 0
    else:
        maxixc = max(ixc)
        maxiyc = max(iyc)
    if len(ixs) == 0:
        maxixs = 0
        maxiys = 0
    else:
        maxixs = max(ixs)
        maxiys = max(iys)
    nx = max(maxixc,maxixs) + 1
    ny = max(maxiyc,maxiys) + 1
    iii = zeros((5,1+nx,1+ny),'l')
    mx,px,my,py = 1,2,3,4
    # --- Flag grid points where the conductors are.
    for i in range(len(ixc)): iii[0,ixc[i],iyc[i]] = i+1
    for i in range(len(ixs)):
        if abs(delmx[i]) < abs(dx): iii[mx,ixs[i],iys[i]] = i+1
        if abs(delpx[i]) < abs(dx): iii[px,ixs[i],iys[i]] = i+1
        if abs(delmy[i]) < abs(dy): iii[my,ixs[i],iys[i]] = i+1
        if abs(delpy[i]) < abs(dy): iii[py,ixs[i],iys[i]] = i+1
    # --- Zero out data for all points internal to a conductor
    iiisum = zeros((3+nx,3+ny),'l')
    iiisum[1:-1,1:-1] = sum(iii,axis=0)
    iiisums = zeros((3+nx,3+ny),'l')
    iiisums[1:-1,1:-1] = sum(iii[1:,:,:],axis=0)
    iiisurf = where(((iiisum[:-2,1:-1]>0)&(iiisum[1:-1,:-2]>0)&
                     (iiisum[2:,1:-1]>0)&(iiisum[1:-1,2:]>0)),
                    iiisums[1:-1,1:-1],iiisum[1:-1,1:-1])
    iiis = zeros((5,1+nx,1+ny),'l')
    for i in range(5):
        iiis[i,:,:] = where(iiisurf>0,iii[i,:,:],0)

    x = []
    y = []
    z = []
    n = []

    # --- Loop over all conductor points, drawing a fill polygon for each.
    ixall = list(ixc) + list(ixs)
    iyall = list(iyc) + list(iys)
    for ix,iy in zip(ixall,iyall):
        if iiisurf[ix,iy] == 0 or iiis[0,ix,iy] > 0: continue
        iix = []
        iiy = []
        z.append(0)
        n.append(0)
        while 1:
            if iiis[0,ix,iy] > 0:
                ixprev = ix
                iyprev = iy
                i0 = iiis[0,ix,iy] - 1
                iiis[0,ix,iy] = 0
                iiisurf[ix,iy] = 0
                x.append(xxc[i0])
                y.append(yyc[i0])
                z[-1] = numb[i0]
                n[-1] = n[-1] + 1
                if ix>0 and iiis[0,ix-1,iy] > 0: ix = ix - 1
                elif ix>0 and iiis[px,ix-1,iy] > 0: ix = ix - 1
                elif ix>0 and iiis[my,ix-1,iy] > 0: ix = ix - 1
                elif ix>0 and iy>0 and iiis[py,ix-1,iy-1] > 0: ix,iy = ix - 1,iy - 1
                elif ix>0 and iy>0 and iiis[0,ix-1,iy-1] > 0: ix,iy = ix - 1,iy - 1
                elif ix>0 and iy>0 and iiis[px,ix-1,iy-1] > 0: ix,iy = ix - 1,iy - 1
                elif iy>0 and iiis[mx,ix,iy-1] > 0: iy = iy - 1
                elif iy>0 and iiis[py,ix,iy-1] > 0: iy = iy - 1
                elif iy>0 and iiis[0,ix,iy-1] > 0: iy = iy - 1
                elif iy>0 and iiis[px,ix,iy-1] > 0: iy = iy - 1
                elif ix<nx and iy>0 and iiis[mx,ix+1,iy-1] > 0: ix,iy = ix + 1,iy - 1
                elif ix<nx and iy>0 and iiis[0,ix+1,iy-1] > 0: ix,iy = ix + 1,iy - 1
                elif ix<nx and iy>0 and iiis[py,ix+1,iy-1] > 0: ix,iy = ix + 1,iy - 1
                elif ix<nx and iiis[my,ix+1,iy] > 0: ix = ix + 1
                elif ix<nx and iiis[mx,ix+1,iy] > 0: ix = ix + 1
                elif ix<nx and iiis[0,ix+1,iy] > 0: ix = ix + 1
                elif ix<nx and iiis[py,ix+1,iy] > 0: ix = ix + 1
                elif ix<nx and iy<ny and iiis[my,ix+1,iy+1] > 0: ix,iy = ix + 1,iy + 1
                elif ix<nx and iy<ny and iiis[0,ix+1,iy+1] > 0: ix,iy = ix + 1,iy + 1
                elif ix<nx and iy<ny and iiis[mx,ix+1,iy+1] > 0: ix,iy = ix + 1,iy + 1
                elif iy<ny and iiis[px,ix,iy+1] > 0: iy = iy + 1
                elif iy<ny and iiis[my,ix,iy+1] > 0: iy = iy + 1
                elif iy<ny and iiis[0,ix,iy+1] > 0: iy = iy + 1
                elif iy<ny and iiis[mx,ix,iy+1] > 0: iy = iy + 1
                elif ix>0 and iy<ny and iiis[px,ix-1,iy+1] > 0: ix,iy = ix - 1,iy + 1
                elif ix>0 and iy<ny and iiis[0,ix-1,iy+1] > 0: ix,iy = ix - 1,iy + 1
                elif ix>0 and iy<ny and iiis[my,ix-1,iy+1] > 0: ix,iy = ix - 1,iy + 1
                elif ix>0 and iiis[py,ix-1,iy] > 0: ix = ix - 1
                continue
            if iiis[mx,ix,iy] > 0 and not (ixprev == ix+1 or
               (ixprev == ix and iiis[my,ix,iy] > 0)):
                ixprev = ix
                iyprev = iy
                i0 = iiis[mx,ix,iy] - 1
                iiis[mx,ix,iy] = 0
                iiisurf[ix,iy] = 0
                x.append(xxs[i0]-delmx[i0])
                y.append(yys[i0])
                z[-1] = numbmx[i0]
                n[-1] = n[-1] + 1
                if iiis[0,ix,iy] > 0: pass
                elif iiis[py,ix,iy] > 0: pass
                elif iy<ny and iiis[0,ix,iy+1] > 0: iy = iy + 1
                elif iy<ny and iiis[mx,ix,iy+1] > 0: iy = iy + 1
                elif iy<ny and ix>0 and iiis[0,ix-1,iy+1] > 0: ix,iy = ix - 1,iy + 1
                elif iy<ny and ix>0 and iiis[my,ix-1,iy+1] > 0: ix,iy = ix - 1,iy + 1
                elif ix>0 and iiis[0,ix-1,iy] > 0: ix = ix - 1
                elif ix>0 and iiis[px,ix-1,iy] > 0: ix = ix - 1
                continue
            if iiis[py,ix,iy] > 0:
                ixprev = ix
                iyprev = iy
                i0 = iiis[py,ix,iy] - 1
                iiis[py,ix,iy] = 0
                iiisurf[ix,iy] = 0
                x.append(xxs[i0])
                y.append(yys[i0]+delpy[i0])
                z[-1] = numbpy[i0]
                n[-1] = n[-1] + 1
                if iiis[0,ix,iy] > 0: pass
                elif iiis[px,ix,iy] > 0: pass
                elif ix<nx and iiis[0,ix+1,iy] > 0: ix = ix + 1
                elif ix<nx and iiis[py,ix+1,iy] > 0: ix = ix + 1
                elif ix<nx and iy<ny and iiis[0,ix+1,iy+1] > 0: ix,iy = ix + 1,iy + 1
                elif ix<nx and iy<ny and iiis[mx,ix+1,iy+1] > 0: ix,iy = ix + 1,iy + 1
                elif iy<ny and iiis[0,ix,iy+1] > 0: iy = iy + 1
                elif iy<ny and iiis[my,ix,iy+1] > 0: iy = iy + 1
                continue
            if iiis[px,ix,iy] > 0:
                ixprev = ix
                iyprev = iy
                i0 = iiis[px,ix,iy] - 1
                iiis[px,ix,iy] = 0
                iiisurf[ix,iy] = 0
                x.append(xxs[i0]+delpx[i0])
                y.append(yys[i0])
                z[-1] = numbpx[i0]
                n[-1] = n[-1] + 1
                if iiis[0,ix,iy] > 0: pass
                elif iiis[my,ix,iy] > 0: pass
                elif iy>0 and iiis[0,ix,iy-1] > 0: iy = iy - 1
                elif iy>0 and iiis[px,ix,iy-1] > 0: iy = iy - 1
                elif iy>0 and ix<nx and iiis[0,ix+1,iy-1] > 0: ix,iy = ix + 1,iy - 1
                elif iy>0 and ix<nx and iiis[py,ix+1,iy-1] > 0: ix,iy = ix + 1,iy - 1
                elif ix<nx and iiis[0,ix+1,iy] > 0: ix = ix + 1
                elif ix<nx and iiis[mx,ix+1,iy] > 0: ix = ix + 1
                continue
            if iiis[my,ix,iy] > 0:
                ixprev = ix
                iyprev = iy
                i0 = iiis[my,ix,iy] - 1
                iiis[my,ix,iy] = 0
                iiisurf[ix,iy] = 0
                x.append(xxs[i0])
                y.append(yys[i0]-delmy[i0])
                z[-1] = numbmy[i0]
                n[-1] = n[-1] + 1
                if iiis[0,ix,iy] > 0: pass
                elif iiis[mx,ix,iy] > 0: pass
                elif ix>0 and iiis[0,ix-1,iy] > 0: ix = ix - 1
                elif ix>0 and iiis[my,ix-1,iy] > 0: ix = ix - 1
                elif ix>0 and iy>0 and iiis[0,ix-1,iy-1] > 0: ix,iy = ix - 1,iy - 1
                elif ix>0 and iy>0 and iiis[px,ix-1,iy-1] > 0: ix,iy = ix - 1,iy - 1
                elif iy>0 and iiis[0,ix,iy-1] > 0: iy = iy - 1
                elif iy>0 and iiis[py,ix,iy-1] > 0: iy = iy - 1
                continue
            break
        if n[-1] == 0:
            del n[-1]
            del z[-1]

    # --- Note that this may not work in parallel!!
    print n
    # --- Now that the data is gathered, make the plot.
    z = array(z).astype(ubyte)
    n = array(n).astype('l')
    plfp(z,y,x,n,local=local)
    plp([y[0]],[x[0]],marker=circle,color=green,local=local)
    plg(y,x,color=red,local=local)

######################################################################
######################################################################
# x-y plane
def pfxy(iz=None,fullplane=1,cond=1,plotsg=1,fill=0,scale=1,
         plotphi=1,plotrho=0,plotselfe=0,comp='z',
         subgridlen=1.,phicolor=blue,rhocolor=red,selfecolor=green,
         condcolor='fg',oddcolor=red,evencolor=green,numb=None,mglevel=0,
         inverted=1,conductors=None,solver=None,kwdict=None,**kw):
    """
  Plots conductors and contours of electrostatic potential in X-Y plane
    - iz=nint(-zmmin/dz): z index of plane
    - fullplane=1: when true, plot all quandrants regardless of symmetries
    - cond=1: when true, plot grid point inside of conductors
    - plotsg=1: when true, plots subgrid data
    - fill=0: when true, fills in conductor
    - scale=1 when true, plots data in lab frame, otherwise grid frame
    - xscale,yscale=1: scaling factor applied to x and y axis
    - plotphi=1: when true, plot contours of potential
    - plotrho=0: when true, plot contours of the charge density
    - plotselfe=0: when true, plot contours or vectors of the electric field
    - comp='z': the component of the electric field to plot, 'x', 'y', 'z' or 'E',
                use 'E' for the magnitude.
    - phicolor=blue: color of phi contours
    - rhocolor=red: color of rho contours
    - selfecolor=green: color of selfe contours or vectors
    - condcolor='fg' color of conductor points inside conductors
    - oddcolor=red color of odd subgrid points
    - evencolor=green color of even subgrid points
    - subgridlen=1 maximum length of subgrid line which are plotted
    - numb: specify which conductors to plot based on the conductor number
    - mglevel=0: level of multigrid to plot data for
    - inverted=1: when false, draws subgrid lines starting inside the conductor
    - See :py:func:`~warpplots.pcphixy` and :py:func:`~warpplots.pcrhoxy` for more plotting options.
    """
    if kwdict is None: kwdict = {}
    kwdict.update(kw)

    # --- These input arguments are also accepted by ppgeneric so are treated
    # --- differently.
    xscale = kwdict.get('xscale',1)
    yscale = kwdict.get('yscale',1)

    if solver is None:
        # --- Check is a solver is registered, and if so, call the appropriate
        # --- method of that instance.
        solver = getregisteredsolver()
        if solver is not None:
            if hasattr(solver,'pfxy'):
                # --- If the solver has pfxy defined, call that instead.
                # --- Note that kw is not a valid keyword and so must be removed. Its
                # --- contents have been put into kwdict.
                kw = copy.copy(locals())
                kw.update(kwdict)
                del kw['kw']
                del kw['kwdict']
                solver.pfxy(**kw)
                return
        else:
            solver = w3d
    if conductors is None:
        if solver is w3d:
            conductors = f3d.conductors
        elif hasattr(solver,'conductors'):
            conductors = solver.conductors
        elif hasattr(solver,'getconductorobject'):
            conductors = solver.getconductorobject()

    # --- This routine by default operates in parallel
    local = kwdict.setdefault('local',0)

    if iz is None:
        if solver.dz != 0.:
            iz = nint(-solver.zmmin/solver.dz)
        else:
            iz = 0
    if w3d.solvergeom != w3d.XYgeom:
        if iz < 0 or solver.nz < iz: return
    if scale:
        dx = solver.dx
        dy = solver.dy
        xmmin = solver.xmmin
        ymmin = solver.ymmin
    else:
        dx = 1.
        dy = 1.
        xmmin = 0.
        ymmin = 0.
    if plotphi or plotrho or plotselfe:
        if not scale:
            kwdict['xmin'] = 0
            kwdict['xmax'] = solver.nx
            kwdict['ymin'] = 0
            kwdict['ymax'] = solver.ny
        if plotphi:
            kwdict['ccolor'] = phicolor
            pcphixy(*(iz,fullplane,solver),**kwdict)
        if plotrho:
            kwdict['ccolor'] = rhocolor
            pcrhoxy(*(iz,fullplane,solver),**kwdict)
        if plotselfe:
            kwdict['ccolor'] = selfecolor
            pcselfexy(*(comp,iz,fullplane,solver),**kwdict)
    if fill:
        plotcondfill(1,0,2,iz,ymmin,xmmin,dy,dx,mglevel,yscale,xscale,conductors,
                     local)
        if fullplane and (solver.l2symtry or solver.l4symtry):
            plotcondfill(1,0,2,iz,ymmin,xmmin,dy,dx,mglevel,yscale,-xscale,conductors,
                         local)
        if fullplane and solver.l4symtry:
            plotcondfill(1,0,2,iz,ymmin,xmmin,dy,dx,mglevel,-yscale,xscale,conductors,
                         local)
            plotcondfill(1,0,2,iz,ymmin,xmmin,dy,dx,mglevel,-yscale,-xscale,
                         conductors,local)
    if cond:
        plotcond(1,0,2,iz,numb,ymmin,xmmin,dy,dx,condcolor,mglevel,yscale,xscale,
                 conductors,local)
        if fullplane and (solver.l2symtry or solver.l4symtry):
            plotcond(1,0,2,iz,numb,ymmin,xmmin,dy,dx,condcolor,mglevel,-yscale,xscale,
                     conductors,local)
        if fullplane and solver.l4symtry:
            plotcond(1,0,2,iz,numb,ymmin,xmmin,dy,dx,condcolor,mglevel,yscale,-xscale,
                     conductors,local)
            plotcond(1,0,2,iz,numb,ymmin,xmmin,dy,dx,condcolor,mglevel,
                     -yscale,-xscale,conductors,local)
    if plotsg:
        plotsubgrid(1,0,2,0,iz,numb,ymmin,xmmin,dy,dx,evencolor,
                    subgridlen,mglevel,yscale,xscale,inverted,conductors,local)
        plotsubgrid(1,0,2,1,iz,numb,ymmin,xmmin,dy,dx,oddcolor,
                    subgridlen,mglevel,yscale,xscale,inverted,conductors,local)
        if fullplane and (solver.l2symtry or solver.l4symtry):
            plotsubgrid(1,0,2,0,iz,numb,ymmin,xmmin,dy,dx,evencolor,
                        subgridlen,mglevel,yscale,-xscale,inverted,conductors,local)
            plotsubgrid(1,0,2,1,iz,numb,ymmin,xmmin,dy,dx,oddcolor,
                        subgridlen,mglevel,yscale,-xscale,inverted,conductors,local)
        if fullplane and solver.l4symtry:
            plotsubgrid(1,0,2,0,iz,numb,ymmin,xmmin,dy,dx,evencolor,
                        subgridlen,mglevel,-yscale,xscale,inverted,conductors,local)
            plotsubgrid(1,0,2,1,iz,numb,ymmin,xmmin,dy,dx,oddcolor,
                        subgridlen,mglevel,-yscale,xscale,inverted,conductors,local)
            plotsubgrid(1,0,2,0,iz,numb,ymmin,xmmin,dy,dx,evencolor,
                        subgridlen,mglevel,-yscale,-xscale,inverted,conductors,local)
            plotsubgrid(1,0,2,1,iz,numb,ymmin,xmmin,dy,dx,oddcolor,
                        subgridlen,mglevel,-yscale,-xscale,inverted,conductors,local)

# z-x plane
def pfzx(iy=None,fullplane=1,lbeamframe=0,
         cond=1,plotsg=1,fill=0,scale=1,
         plotphi=1,plotrho=0,plotselfe=0,comp='z',
         subgridlen=1.,phicolor=blue,rhocolor=red,selfecolor=green,
         condcolor='fg',oddcolor=red,evencolor=green,numb=None,mglevel=0,
         inverted=1,conductors=None,solver=None,kwdict=None,**kw):
    """
  Plots conductors and contours of electrostatic potential in Z-X plane
    - iy=nint(-ymmin/dy): y index of plane
    - fullplane=1: when true, plot all quadrants regardless of symmetries
    - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
    - cond=1: when true, plot grid point inside of conductors
    - plotsg=1: when true, plots subgrid data
    - fill=0: when true, fills in conductor
    - scale=1 when true, plots data in lab frame, otherwise grid frame
    - xscale,yscale=1: scaling factor applied to x and y axis
    - plotphi=1: when true, plot contours of potential
    - plotrho=0: when true, plot contours of the charge density
    - plotselfe=0: when true, plot contours or vectors of the electric field
    - comp='z': the component of the electric field to plot, 'x', 'y', 'z' or 'E',
                use 'E' for the magnitude.
    - phicolor=blue: color of phi contours
    - rhocolor=red: color of rho contours
    - selfecolor=green: color of selfe contours or vectors
    - condcolor='fg' color of conductor points inside conductors
    - oddcolor=red color of odd subgrid points
    - evencolor=green color of even subgrid points
    - subgridlen=1 maximum length of subgrid line which are plotted
    - numb: specify which conductors to plot based on the conductor number
    - mglevel=0: level of multigrid to plot data for
    - inverted=1: when false, draws subgrid lines starting inside the conductor
    - See :py:func:`~warpplots.pcphizx` and :py:func:`~warpplots.pcrhozx` for more plotting options.
    """
    if kwdict is None: kwdict = {}
    kwdict.update(kw)

    # --- These input arguments are also accepted by ppgeneric so are treated
    # --- differently.
    xscale = kwdict.get('xscale',1)
    yscale = kwdict.get('yscale',1)

    if solver is None:
        # --- Check is a solver is registered, and if so, call the appropriate
        # --- method of that instance.
        solver = getregisteredsolver()
        if solver is not None:
            if hasattr(solver,'pfzx'):
                # --- If the solver has pfzx defined, call that instead.
                # --- Note that kw is not a valid keyword and so must be removed. Its
                # --- contents have been put into kwdict.
                kw = copy.copy(locals())
                kw.update(kwdict)
                del kw['kw']
                del kw['kwdict']
                solver.pfzx(**kw)
                return
        else:
            solver = w3d
    if conductors is None:
        if solver is w3d:
            conductors = f3d.conductors
        elif hasattr(solver,'conductors'):
            conductors = solver.conductors
        elif hasattr(solver,'getconductorobject'):
            conductors = solver.getconductorobject()

    # --- This routine by default operates in parallel
    local = kwdict.setdefault('local',0)

    if iy is None:
        if solver.dy != 0.:
            iy = nint(-solver.ymmin/solver.dy)
        else:
            iy = 0
    if iy < 0 or solver.ny < iy: return
    if lbeamframe: zbeam = 0.
    else:          zbeam = top.zbeam
    if scale:
        dx = solver.dx
        dz = solver.dz
        xmmin = solver.xmmin
        zmmin = solver.zmmin + zbeam
    else:
        dx = 1.
        dz = 1.
        xmmin = 0.
        zmmin = 0.
    if plotphi or plotrho or plotselfe:
        if not scale:
            kwdict['xmin'] = 0
            kwdict['xmax'] = solver.nz
            kwdict['ymin'] = 0.
            kwdict['ymax'] = solver.nx
        if plotphi:
            kwdict['ccolor'] = phicolor
            pcphizx(*(iy,fullplane,lbeamframe,solver),**kwdict)
        if plotrho:
            kwdict['ccolor'] = rhocolor
            pcrhozx(*(iy,fullplane,lbeamframe,solver),**kwdict)
        if plotselfe:
            kwdict['ccolor'] = selfecolor
            pcselfezx(*(comp,iy,fullplane,solver),**kwdict)
    if fill:
        plotcondfill(0,2,1,iy,xmmin,zmmin,dx,dz,mglevel,yscale,xscale,conductors,
                     local)
        if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
            plotcondfill(0,2,1,iy,xmmin,zmmin,dx,dz,mglevel,-yscale,xscale,conductors,
                         local)
    if plotsg:
        plotsubgrid(0,2,1,0,iy,numb,xmmin,zmmin,dx,dz,evencolor,
                    subgridlen,mglevel,yscale,xscale,inverted,conductors,local)
        plotsubgrid(0,2,1,1,iy,numb,xmmin,zmmin,dx,dz,oddcolor,
                    subgridlen,mglevel,yscale,xscale,inverted,conductors,local)
        if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
            plotsubgrid(0,2,1,0,iy,numb,xmmin,zmmin,dx,dz,evencolor,
                        subgridlen,mglevel,-yscale,xscale,inverted,conductors,local)
            plotsubgrid(0,2,1,1,iy,numb,xmmin,zmmin,dx,dz,oddcolor,
                        subgridlen,mglevel,-yscale,xscale,inverted,conductors,local)
    if cond:
        plotcond(0,2,1,iy,numb,xmmin,zmmin,dx,dz,condcolor,mglevel,yscale,xscale,
                 conductors,local)
        if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
            plotcond(0,2,1,iy,numb,xmmin,zmmin,dx,dz,condcolor,mglevel,-yscale,xscale,
                     conductors,local)

# z-r plane
def pfzr(**kw):
    if 'fullplane' not in kw: kw['fullplane'] = 0
    pfzx(**kw)

# z-y plane
def pfzy(ix=None,fullplane=1,lbeamframe=0,
         cond=1,plotsg=1,fill=0,scale=1,
         plotphi=1,plotrho=0,plotselfe=0,comp='z',
         subgridlen=1.,phicolor=blue,rhocolor=red,selfecolor=green,
         condcolor='fg',oddcolor=red,evencolor=green,numb=None,mglevel=0,
         inverted=1,conductors=None,solver=None,kwdict=None,**kw):
    """
  Plots conductors and contours of electrostatic potential in Z-Y plane
    - ix=nint(-xmmin/dx): x index of plane
    - fullplane=1: when true, plot all quadrants regardless of symmetries
    - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
    - cond=1: when true, plot grid point inside of conductors
    - plotsg=1: when true, plots subgrid data
    - fill=0: when true, fills in conductor
    - scale=1 when true, plots data in lab frame, otherwise grid frame
    - xscale,yscale=1: scaling factor applied to x and y axis
    - plotphi=1: when true, plot contours of potential
    - plotrho=0: when true, plot contours of the charge density
    - plotselfe=0: when true, plot contours or vectors of the electric field
    - comp='z': the component of the electric field to plot, 'x', 'y', 'z' or 'E',
                use 'E' for the magnitude.
    - phicolor=blue: color of phi contours
    - rhocolor=red: color of rho contours
    - selfecolor=green: color of selfe contours or vectors
    - condcolor='fg' color of conductor points inside conductors
    - oddcolor=red color of odd subgrid points
    - evencolor=green color of even subgrid points
    - subgridlen=1 maximum length of subgrid line which are plotted
    - numb: specify which conductors to plot based on the conductor number
    - mglevel=0: level of multigrid to plot data for
    - inverted=1: when false, draws subgrid lines starting inside the conductor
    - See :py:func:`~warpplots.pcphizy` and :py:func:`~warpplots.pcrhozy` for more plotting options.
    """
    if kwdict is None: kwdict = {}
    kwdict.update(kw)

    # --- These input arguments are also accepted by ppgeneric so are treated
    # --- differently.
    xscale = kwdict.get('xscale',1)
    yscale = kwdict.get('yscale',1)

    if solver is None:
        # --- Check is a solver is registered, and if so, call the appropriate
        # --- method of that instance.
        solver = getregisteredsolver()
        if solver is not None:
            if hasattr(solver,'pfzy'):
                # --- If the solver has pfzy defined, call that instead.
                # --- Note that kw is not a valid keyword and so must be removed. Its
                # --- contents have been put into kwdict.
                kw = copy.copy(locals())
                kw.update(kwdict)
                del kw['kw']
                del kw['kwdict']
                solver.pfzy(**kw)
                return
        else:
            solver = w3d
    if conductors is None:
        if solver is w3d:
            conductors = f3d.conductors
        elif hasattr(solver,'conductors'):
            conductors = solver.conductors
        elif hasattr(solver,'getconductorobject'):
            conductors = solver.getconductorobject()

    # --- This routine by default operates in parallel
    local = kwdict.setdefault('local',0)

    if ix is None:
        if solver.dx != 0.:
            ix = nint(-solver.xmmin/solver.dx)
        else:
            ix = 0
    if ix < 0 or solver.nx < ix: return
    if lbeamframe: zbeam = 0.
    else:          zbeam = top.zbeam
    if scale:
        dy = solver.dy
        dz = solver.dz
        ymmin = solver.ymmin
        zmmin = solver.zmmin + zbeam
    else:
        dy = 1.
        dz = 1.
        ymmin = 0.
        zmmin = 0.
    if plotphi or plotrho or plotselfe:
        if not scale:
            kwdict['xmin'] = 0
            kwdict['xmax'] = solver.nz
            kwdict['ymin'] = 0
            kwdict['ymax'] = solver.ny
        if plotphi:
            kwdict['ccolor'] = phicolor
            pcphizy(*(ix,fullplane,lbeamframe,solver),**kwdict)
        if plotrho:
            kwdict['ccolor'] = rhocolor
            pcrhozy(*(ix,fullplane,lbeamframe,solver),**kwdict)
        if plotselfe:
            kwdict['ccolor'] = selfecolor
            pcselfezy(*(comp,ix,fullplane,solver),**kwdict)
    if fill:
        plotcondfill(1,2,0,ix,ymmin,zmmin,dy,dz,mglevel,yscale,xscale,conductors,
                     local)
        if fullplane and (solver.l2symtry or solver.l4symtry):
            plotcondfill(1,2,0,ix,ymmin,zmmin,dy,dz,mglevel,-yscale,xscale,conductors,
                         local)
    if cond:
        plotcond(1,2,0,ix,numb,ymmin,zmmin,dy,dz,condcolor,mglevel,yscale,xscale,
                 conductors,local)
        if fullplane and (solver.l2symtry or solver.l4symtry):
            plotcond(1,2,0,ix,numb,ymmin,zmmin,dy,dz,condcolor,mglevel,-yscale,xscale,
                     conductors,local)
    if plotsg:
        plotsubgrid(1,2,0,0,ix,numb,ymmin,zmmin,dy,dz,evencolor,
                    subgridlen,mglevel,yscale,xscale,inverted,conductors,local)
        plotsubgrid(1,2,0,1,ix,numb,ymmin,zmmin,dy,dz,oddcolor,
                    subgridlen,mglevel,yscale,xscale,inverted,conductors,local)
        if fullplane and (solver.l2symtry or solver.l4symtry):
            plotsubgrid(1,2,0,0,ix,numb,ymmin,zmmin,dy,dz,evencolor,
                        subgridlen,mglevel,-yscale,xscale,inverted,conductors,local)
            plotsubgrid(1,2,0,1,ix,numb,ymmin,zmmin,dy,dz,oddcolor,
                        subgridlen,mglevel,-yscale,xscale,inverted,conductors,local)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in grid units                                                      #
######################################################################

# x-y plane
def pfxyg(iz=None,fullplane=1,
          cond=1,plotsg=1,fill=0,plotphi=1,plotrho=0,plotselfe=0,comp='z',
          phicolor=blue,rhocolor=red,selfecolor=green,
          subgridlen=1.,condcolor='fg',oddcolor=red,evencolor=green,
          numb=None,mglevel=0,inverted=1,conductors=None,solver=None,**kw):
    """
  Plots conductors and contours of electrostatic potential in X-Y plane in grid
  frame
  Same arguments as py:func:`pfxy`
    """
    pfxy(iz=iz,fullplane=fullplane,scale=0,cond=cond,plotsg=plotsg,fill=fill,
         plotphi=plotphi,plotrho=plotrho,plotselfe=plotselfe,comp=comp,
         subgridlen=subgridlen,
         phicolor=phicolor,rhocolor=rhocolor,selfecolor=selfecolor,
         condcolor=condcolor,oddcolor=oddcolor,evencolor=evencolor,
         numb=numb,mglevel=mglevel,inverted=inverted,conductors=conductors,
         solver=solver,kwdict=kw)

# z-x plane
def pfzxg(iy=None,fullplane=1,lbeamframe=0,
          cond=1,plotsg=1,fill=0,plotphi=1,plotrho=0,plotselfe=0,comp='z',
          subgridlen=1.,phicolor=blue,rhocolor=red,selfecolor=green,
          condcolor='fg',oddcolor=red,evencolor=green,numb=None,mglevel=0,
          inverted=1,conductors=None,solver=None,**kw):
    """
  Plots conductors and contours of electrostatic potential in Z-X plane in grid
  frame
  Same arguments as py:func:`pfzx`
    """
    pfzx(iy=iy,fullplane=fullplane,lbeamframe=lbeamframe,scale=0,
         cond=cond,plotsg=plotsg,fill=fill,
         plotphi=plotphi,plotrho=plotrho,plotselfe=plotselfe,comp=comp,
         subgridlen=subgridlen,
         phicolor=phicolor,rhocolor=rhocolor,selfecolor=selfecolor,
         condcolor=condcolor,oddcolor=oddcolor,evencolor=evencolor,
         numb=numb,mglevel=mglevel,inverted=inverted,conductors=conductors,
         solver=solver,kwdict=kw)

# z-r plane
pfzrg = pfzxg

# z-y plane
def pfzyg(ix=None,fullplane=1,lbeamframe=0,
          cond=1,plotsg=1,fill=0,plotphi=1,plotrho=0,plotselfe=0,comp='z',
          subgridlen=1.,phicolor=blue,rhocolor=red,selfecolor=green,
          condcolor='fg',oddcolor=red,evencolor=green,numb=None,mglevel=0,
          inverted=1,conductors=None,solver=None,**kw):
    """
  Plots conductors and contours of electrostatic potential in Z-Y plane in grid
  frame
  Same arguments as py:func:`pfzy`
    """
    pfzy(ix=ix,fullplane=fullplane,lbeamframe=lbeamframe,scale=0,
         cond=cond,plotsg=plotsg,fill=fill,
         plotphi=plotphi,plotrho=plotrho,plotselfe=plotselfe,comp=comp,
         subgridlen=subgridlen,
         phicolor=phicolor,rhocolor=rhocolor,selfecolor=selfecolor,
         condcolor=condcolor,oddcolor=oddcolor,evencolor=evencolor,
         numb=numb,mglevel=mglevel,inverted=inverted,conductors=conductors,
         solver=solver,kwdict=kw)

######################################################################

############################################################################
# These plot a box at each conductor point
############################################################################

# x-y plane
def pfxybox(iz=None,contours=8,plotsg=1,scale=1,signx=1,signy=1,
            plotphi=1,plotrho=0,plotselfe=0,comp='z',filled=0,
            phicolor=blue,rhocolor=red,selfecolor=green,
            condcolor='fg',conductors=None,solver=None,
            kwdict=None,**kw):
    """
  Plots square at conductor points and contours of electrostatic potential
  in X-Y plane
    - iz=nint(-zmmin/dz): z index of plane
    - contours=8 optional number of or list of contours
    - plotsg=1 when true, plots subgrid data
    - scale=1 when true, plots data in lab frame, otherwise grid frame
    - signx=1 sign of x, used for plotting symmetry planes
    - signy=1 sign of y, used for plotting symmetry planes
    - plotphi=1: when true, plot contours of potential
    - plotrho=0: when true, plot contours of the charge density
    - plotselfe=0: when true, plot contours or vectors of the electric field
    - comp='z': the component of the electric field to plot, 'x', 'y', 'z' or 'E',
                use 'E' for the magnitude.
    - filled=0 when true, plots filled contours
    - phicolor=blue: color of phi contours
    - rhocolor=red: color of rho contours
    - selfecolor=green: color of selfe contours or vectors
    - condcolor='fg' color of conductor points inside conductors
    - subgridlen=1 maximum length of subgrid line which are plotted
    """
    if kwdict is None: kwdict = {}
    kwdict.update(kw)

    if solver is None:
        # --- Check is a solver is registered, and if so, call the appropriate
        # --- method of that instance.
        solver = getregisteredsolver()
        if solver is not None:
            if hasattr(solver,'pfxybox'):
                # --- If the solver has pfxybox defined, call that instead.
                # --- Note that kw is not a valid keyword and so must be removed. Its
                # --- contents have been put into kwdict.
                kw = copy.copy(locals())
                kw.update(kwdict)
                del kw['kw']
                del kw['kwdict']
                solver.pfxybox(**kw)
                return
        else:
            solver = w3d
    if conductors is None:
        if solver is w3d:
            conductors = f3d.conductors
        elif hasattr(solver,'conductors'):
            conductors = solver.conductors
        elif hasattr(solver,'getconductorobject'):
            conductors = solver.getconductorobject()

    # --- This routine by default operates in parallel
    local = kwdict.setdefault('local',0)

    if iz is None: iz = solver.iz_axis
    if iz < 0 or solver.nz < iz: return
    if scale:
        dy = solver.dy*signy
        dx = solver.dx*signx
        ymmin = solver.ymmin
        xmmin = solver.xmmin
    else:
        dy = 1.*signy
        dx = 1.*signx
        ymmin = 0.
        xmmin = 0.
    if plotphi or plotrho or plotselfe:
        if not scale:
            kwdict['xmin'] = 0
            kwdict['xmax'] = solver.nx
            kwdict['ymin'] = 0
            kwdict['ymax'] = solver.ny
        if plotphi:
            kwdict['ccolor'] = phicolor
            pcphixy(*(iz,fullplane,solver),**kwdict)
        if plotrho:
            kwdict['ccolor'] = rhocolor
            pcrhoxy(*(iz,fullplane,solver),**kwdict)
        if plotselfe:
            kwdict['ccolor'] = selfecolor
            pcselfexy(*(comp,iz,fullplane,solver),**kwdict)
    if conductors.interior.n > 0:
        n = conductors.interior.n
        izc = conductors.indx[2,0:n]+conductors.leveliz[0]
        ii = compress(equal(izc,iz),arange(n))
        x = take(conductors.indx[0,0:n],ii)*dx+xmmin
        y = take(conductors.indx[1,0:n],ii)*dy+ymmin
    else:
        x = array([])
        y = array([])
    if len(x) > 0:
        pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
            array([x-dx/2,x+dx/2,x+dx/2,x-dx/2,x-dx/2]),
            color=condcolor,local=local)

# z-x plane
def pfzxbox(iy=None,contours=8,plotsg=1,scale=1,signz=1,signx=1,
            plotphi=1,plotrho=0,plotselfe=0,comp='z',filled=0,
            phicolor=blue,rhocolor=red,selfecolor=green,
            condcolor='fg',conductors=None,solver=None,
            kwdict=None,**kw):
    """
  Plots square at conductor points and contours of electrostatic potential
  in Z-X plane
    - iy=nint(-ymmin/dy): y index of plane
    - contours=8 optional number of or list of contours
    - plotsg=1 when true, plots subgrid data
    - scale=1 when true, plots data in lab frame, otherwise grid frame
    - signz=1 sign of z, used for plotting symmetry planes
    - signx=1 sign of x, used for plotting symmetry planes
    - plotphi=1: when true, plot contours of potential
    - plotrho=0: when true, plot contours of the charge density
    - plotselfe=0: when true, plot contours or vectors of the electric field
    - comp='z': the component of the electric field to plot, 'x', 'y', 'z' or 'E',
                use 'E' for the magnitude.
    - filled=0 when true, plots filled contours
    - phicolor=blue: color of phi contours
    - rhocolor=red: color of rho contours
    - selfecolor=green: color of selfe contours or vectors
    - condcolor='fg' color of conductor points inside conductors
    """
    if kwdict is None: kwdict = {}
    kwdict.update(kw)

    if solver is None:
        # --- Check is a solver is registered, and if so, call the appropriate
        # --- method of that instance.
        solver = getregisteredsolver()
        if solver is not None:
            if hasattr(solver,'pfzxbox'):
                # --- If the solver has pfzxbox defined, call that instead.
                # --- Note that kw is not a valid keyword and so must be removed. Its
                # --- contents have been put into kwdict.
                kw = copy.copy(locals())
                kw.update(kwdict)
                del kw['kw']
                del kw['kwdict']
                solver.pfzxbox(**kw)
                return
        else:
            solver = w3d
    if conductors is None:
        if solver is w3d:
            conductors = f3d.conductors
        elif hasattr(solver,'conductors'):
            conductors = solver.conductors
        elif hasattr(solver,'getconductorobject'):
            conductors = solver.getconductorobject()

    # --- This routine by default operates in parallel
    local = kwdict.setdefault('local',0)

    if not iy:
        if solver.dy != 0.:
            iy = nint(-solver.ymmin/solver.dy)
        else:
            iy = 0
    if iy < 0 or solver.ny < iy: return
    if scale:
        dx = solver.dx*signx
        dz = solver.dz*signz
        xmmin = solver.xmmin
        zmmin = solver.zmmin
    else:
        dx = 1.*signx
        dz = 1.*signz
        xmmin = 0.
        zmmin = 0.
    if plotphi or plotrho or plotselfe:
        if not scale:
            kwdict['xmin'] = 0
            kwdict['xmax'] = solver.nx
            kwdict['ymin'] = 0
            kwdict['ymax'] = solver.ny
        if plotphi:
            kwdict['ccolor'] = phicolor
            pcphizx(*(iy,fullplane,solver),**kwdict)
        if plotrho:
            kwdict['ccolor'] = rhocolor
            pcrhozx(*(iy,fullplane,solver),**kwdict)
        if plotselfe:
            kwdict['ccolor'] = selfecolor
            pcselfezx(*(comp,iy,fullplane,solver),**kwdict)
    n = conductors.interior.n
    if (n > 0):
        ii = compress(equal(conductors.interior.indx[1,0:n],iy),arange(n))
        x = take(conductors.interior.indx[0,0:n],ii)*dx+xmmin
        z = take(conductors.interior.indx[2,0:n],ii)*dz+zmmin
    else:
        x = array([])
        z = array([])
    if len(x) > 0:
        pla(array([x-dx/2,x-dx/2,x+dx/2,x+dx/2,x-dx/2]),
            array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
            color=condcolor,local=local)

# z-y plane
def pfzybox(ix=None,contours=8,plotsg=1,scale=1,signz=1,signy=1,
            plotphi=1,plotrho=0,plotselfe=0,comp='z',filled=0,
            phicolor=blue,rhocolor=red,selfecolor=green,
            condcolor='fg',conductors=None,solver=None,
            kwdict=None,**kw):
    """
  Plots square at conductor points and contours of electrostatic potential
  in Z-Y plane
    - ix=nint(-xmmin/dx): x index of plane
    - contours=8 optional number of or list of contours
    - plotsg=1 when true, plots subgrid data
    - scale=1 when true, plots data in lab frame, otherwise grid frame
    - signz=1 sign of z, used for plotting symmetry planes
    - signy=1 sign of y, used for plotting symmetry planes
    - plotphi=1: when true, plot contours of potential
    - plotrho=0: when true, plot contours of the charge density
    - plotselfe=0: when true, plot contours or vectors of the electric field
    - comp='z': the component of the electric field to plot, 'x', 'y', 'z' or 'E',
                use 'E' for the magnitude.
    - filled=0 when true, plots filled contours
    - phicolor=blue: color of phi contours
    - rhocolor=red: color of rho contours
    - selfecolor=green: color of selfe contours or vectors
    - condcolor='fg' color of conductor points inside conductors
    """
    if kwdict is None: kwdict = {}
    kwdict.update(kw)

    if solver is None:
        # --- Check is a solver is registered, and if so, call the appropriate
        # --- method of that instance.
        solver = getregisteredsolver()
        if solver is not None:
            if hasattr(solver,'pfzybox'):
                # --- If the solver has pfzybox defined, call that instead.
                # --- Note that kw is not a valid keyword and so must be removed. Its
                # --- contents have been put into kwdict.
                kw = copy.copy(locals())
                kw.update(kwdict)
                del kw['kw']
                del kw['kwdict']
                solver.pfzybox(**kw)
                return
        else:
            solver = w3d
    if conductors is None:
        if solver is w3d:
            conductors = f3d.conductors
        elif hasattr(solver,'conductors'):
            conductors = solver.conductors
        elif hasattr(solver,'getconductorobject'):
            conductors = solver.getconductorobject()

    # --- This routine by default operates in parallel
    local = kwdict.setdefault('local',0)

    if not ix:
        if solver.dx != 0.:
            ix = nint(-solver.xmmin/solver.dx)
        else:
            ix = 0
    if ix < 0 or solver.nx < ix: return
    if scale:
        dy = solver.dy*signy
        dz = solver.dz*signz
        ymmin = solver.ymmin
        zmmin = solver.zmmin
    else:
        dy = 1.*signy
        dz = 1.*signz
        ymmin = 0.
        zmmin = 0.
    if plotphi or plotrho or plotselfe:
        if not scale:
            kwdict['xmin'] = 0
            kwdict['xmax'] = solver.nx
            kwdict['ymin'] = 0
            kwdict['ymax'] = solver.ny
        if plotphi:
            kwdict['ccolor'] = phicolor
            pcphizy(*(ix,fullplane,solver),**kwdict)
        if plotrho:
            kwdict['ccolor'] = rhocolor
            pcrhozy(*(ix,fullplane,solver),**kwdict)
        if plotselfe:
            kwdict['ccolor'] = selfecolor
            pcselfezy(*(comp,ix,fullplane,solver),**kwdict)
    n = conductors.interior.n
    if (n > 0):
        ii = compress(equal(conductors.interior.indx[0,0:n],ix),arange(n))
        y = take(conductors.interior.indx[1,0:n],ii)*dy+ymmin
        z = take(conductors.interior.indx[2,0:n],ii)*dz+zmmin
    else:
        y = array([])
        z = array([])
    if len(y) > 0:
        pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
            array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
            color=condcolor,local=local)

# z-x plane
def pfzxboxi(iy=None,contours=8,plotsg=1,scale=1,signz=1,
             plotphi=1,plotrho=0,plotselfe=0,comp='z',
             filled=0,phicolor=blue,rhocolor=red,selfecolor=green,
             condcolor='fg',conductors=None,solver=None,**kw):
    """
  Plots square at conductor points and contours of electrostatic potential
  in Z-(-X) plane
    - iy=nint(-ymmin/dy): y index of plane
    - contours=8 optional number of or list of contours
    - plotsg=1 when true, plots subgrid data
    - scale=1 when true, plots data in lab frame, otherwise grid frame
    - signz=1 sign of z, used for plotting symmetry planes
    - plotphi=1: when true, plot contours of potential
    - plotrho=0: when true, plot contours of the charge density
    - plotselfe=0: when true, plot contours or vectors of the electric field
    - comp='z': the component of the electric field to plot, 'x', 'y', 'z' or 'E',
                use 'E' for the magnitude.
    - filled=0 when true, plots filled contours
    - phicolor=blue: color of phi contours
    - rhocolor=red: color of rho contours
    - selfecolor=green: color of selfe contours or vectors
    - condcolor='fg' color of conductor points inside conductors
    """
    pfzxbox(iy=iy,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
            signx=-1,
            plotphi=plotphi,plotrho=plotrho,plotselfe=plotselfe,comp=comp,
            filled=filled,
            phicolor=phicolor,rhocolor=rhocolor,selfecolor=selfecolor,
            condcolor=condcolor,conductors=conductors,solver=w3d,kwdict=kw)

# z-y plane
def pfzyboxi(ix=None,contours=8,plotsg=1,scale=1,signz=1,signy=-1,
             plotphi=1,plotrho=0,plotselfe=0,comp='z',
             filled=0,phicolor=blue,rhocolor=red,selfecolor=green,
             condcolor='fg',conductors=None,solver=None,**kw):
    """
  Plots square at conductor points and contours of electrostatic potential
  in Z-(-Y) plane
    - ix=nint(-xmmin/dx): x index of plane
    - contours=8 optional number of or list of contours
    - plotsg=1 when true, plots subgrid data
    - scale=1 when true, plots data in lab frame, otherwise grid frame
    - signz=1 sign of z, used for plotting symmetry planes
    - plotphi=1: when true, plot contours of potential
    - plotrho=0: when true, plot contours of the charge density
    - plotselfe=0: when true, plot contours or vectors of the electric field
    - comp='z': the component of the electric field to plot, 'x', 'y', 'z' or 'E',
                use 'E' for the magnitude.
    - filled=0 when true, plots filled contours
    - phicolor=blue: color of phi contours
    - rhocolor=red: color of rho contours
    - selfecolor=green: color of selfe contours or vectors
    - condcolor='fg' color of conductor points inside conductors
    """
    pfzybox(ix=ix,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
            signy=-1,
            plotphi=plotphi,plotrho=plotrho,plotselfe=plotselfe,comp=comp,
            filled=filled,
            phicolor=phicolor,rhocolor=rhocolor,selfecolor=selfecolor,
            condcolor=condcolor,conductors=conductors,solver=w3d,kwdict=kw)




############################################################################
# These plots plot the conductor points colored based on the conductor
# number. The list of colors is input by the user.

# --- convenience function
def findunique(i):
    ii = sort(i)
    result = list(compress(ii[:-1]!=ii[1:],ii[:-1])) + [ii[-1]]
    return result

def plotcondn(yy,xx,zz,iz,ymmin,xmmin,dy,dx,mglevel,signy,signx,conductors,
              local):
    if conductors is None: return
    ncolor = len(color)
    if conductors.interior.n > 0: nn = conductors.interior.numb
    else:             nn = array([])
    nlist = gatherarray(nn)
    nlist = findunique(nlist)
    nlist = warp_parallel.broadcast(nlist)
    for i in nlist:
        plotcond(yy,xx,zz,iz,i,ymmin,xmmin,dy,dx,color[i%ncolor],
                 mglevel,signy,signx,conductors,local)

def pfzxn(iy=None,numbs=None,colors=None,cmarker=point,smarker=circle,
          scale=1,signz=1,signx=1,subgridlen=1.,fullplane=1,mglevel=0,
          inverted=1,conductors=f3d.conductors,solver=w3d,local=0):
    if iy is None:
        if solver.dy != 0.:
            iy = nint(-solver.ymmin/solver.dy)
        else:
            iy = 0
    if iy < 0 or solver.ny < iy: return
    if colors is None: colors = color
    if scale:
        dx = solver.dx*signx
        dz = solver.dz*signz
        xmmin = solver.xmmin
        zmmin = solver.zmmin
    else:
        dx = 1.*signx
        dz = 1.*signz
        xmmin = 0.
        zmmin = 0.
    plotcondn(0,2,1,iy,xmmin,zmmin,dx,dz,mglevel,1,1,conductors,local)
    if fullplane and solver.l4symtry:
        plotcondn(0,2,1,iy,xmmin,zmmin,dx,dz,mglevel,-1,1,conductors,local)
    ncolor = len(colors)
    nlist = gatherarray(conductors.evensubgrid.numb[0,:conductors.evensubgrid.n])
    nlist = findunique(nlist)
    #nlist.remove(0)
    nlist = warp_parallel.broadcast(nlist)
    for i in nlist:
        plotsubgrid(0,2,1,0,iy,i,xmmin,zmmin,dx,dz,
                    colors[i%ncolor],subgridlen,mglevel,1,1,inverted,conductors,local)
        if fullplane and solver.l4symtry:
            plotsubgrid(0,2,1,0,iy,i,xmmin,zmmin,dx,dz,
                        colors[i%ncolor],subgridlen,mglevel,-1,1,inverted,conductors,local)
    nlist = gatherarray(conductors.oddsubgrid.numb[0,:conductors.oddsubgrid.n])
    nlist = findunique(nlist)
    nlist = warp_parallel.broadcast(nlist)
    for i in nlist:
        plotsubgrid(0,2,1,1,iy,i,xmmin,zmmin,dx,dz,
                    colors[i%ncolor],subgridlen,mglevel,1,1,inverted,conductors,local)
        if fullplane and solver.l4symtry:
            plotsubgrid(0,2,1,1,iy,i,xmmin,zmmin,dx,dz,
                        colors[i%ncolor],subgridlen,mglevel,-1,1,inverted,conductors,local)


############################################################################
# These plot the conductors in laboratory frame, using the tolabfrm routine
# to convert from code frame to lab frame.  There is also a routine to plot
# the computational x-z grid in lab frame.
############################################################################

# --- plot grid in lab frame (including bends)
def plotgrid(zbeam=None,ii=2,plotcond=1,solver=w3d,zcent=None):
    """
  Plots Z-X grid in the lab frame (including bends)
    - zbeam=top.zbeam is the center position
    - ii=2 is the step size in the grid points plotted
      2 means that every other grid line is plotted
    - plotcond=1 when true, plots conductors
    """
    if zbeam is None: zbeam=top.zbeam
    if zcent is None: zcent=top.zbeam
    # --- declare temporary data space, 2 2-D arrays to hold grid coordinates
    nx = max(1,nint(1.*solver.nx/ii))
    dx = 1.*solver.nx/nx*solver.dx
    nz = max(1,nint(1.*solver.nz/ii))
    dz = 1.*solver.nz/nz*solver.dz
    xg,zg = getmesh2d(solver.xmmin,dx,nx,solver.zmmin+zbeam,dz,nz)

    # --- If in a bend, convert the grid data to the lab frame
    if top.bends:
        # --- reshape arrays to make 1-D arrays to pass to tolabfrm
        nn = (nx+1)*(nz+1)
        xg.shape = (nn,)
        zg.shape = (nn,)

        # --- Convert data to lab frame
        tolabfrm(zcent,nn,xg,zg)

        # --- Reshape back into 2-D arrays
        xg.shape = (nx+1,nz+1)
        zg.shape = (nx+1,nz+1)

    # --- Make plots
    pla(xg,zg,marks=0)
    pla(transpose(xg),transpose(zg),marks=0)

    if plotcond: pfzxlab(zbeam)

# --- Make pfzx plot in lab frame
def pfzxlab(zz=None,iy=None,condcolor='fg',conductors=f3d.conductors,
            solver=w3d,mglevel=0):
    """Plots conductors in Z-X lab frame (including bends)
    - zz=top.zbeam is the center position
    - condcolor='fg' color of conductor points inside conductors
    """
    if not zz: zz=top.zbeam
    if iy is None:
        if solver.dy != 0.:
            iy = nint(-solver.ymmin/solver.dy)
        else:
            iy = 0
    # --- if zz is not equal to zbeam, then calculate conductors for new location
    if (zz != top.zbeam):
        z = top.zbeam
        g = top.zgrid
        top.zbeam = zz
        top.zgrid = zz
        setlatt()
        fieldsol(1)
    # --- gather conductor data
    if conductors.interior.n > 0:
        interior = conductors.interior
        nn = interior.n
        if nn > 0:
            lx = conductors.levellx[0]
            ly = conductors.levelly[0]
            lz = conductors.levellz[0]
            ixc = interior.indx[0,:] + conductors.levelix[0]
            iyc = interior.indx[1,:] + conductors.leveliy[0]
            izc = interior.indx[2,:] + conductors.leveliz[0]
            level = equal(mglevel,interior.ilevel)
        ii = compress(logical_and(equal(iyc[:nn],iy),level[:nn]),arange(nn))

        xl=take(conductors.interior.indx[0,:],ii)*solver.dx+solver.xmmin
        zl=take(conductors.interior.indx[2,:],ii)*solver.dz+solver.zmmin+zz
        # --- convert to lab frame
        tolabfrm(zz,len(xl),xl,zl)
        # --- make plot
        plp(xl,zl,color=condcolor)
    # --- restore original conductor data at zbeam
    if (zz != top.zbeam):
        top.zbeam = z
        top.zgrid = g
        setlatt()
        fieldsol(1)


#####################################################################
def plotsrfrvinout(srfrvin,srfrvout,zmin,zmax,n=1000,color='fg',gridframe=0,
                   rscale=1,zscale=1,
                   roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0,
                   outline=1,fillcolor=None):
    """Handy function for plotting the r versus z for a surface of revolution
   - srfrvin: surface of revolution function to plot
   - srfrvou: surface of revolution function to plot
   - zmin,zmax: z range to plot
   - n=1000: number of points to plot
   - color='fg': color of line
   - gridframe=0: when true, plots in grid frame
   - rscale=1: scaling for radius
   - zscale=1: scaling for z
   - roff=0: offset for radius
   - zoff=0: offset for z
   - rmin=0: minimum value of r plotted (before applying rscale and roff)
   - rmax=0: maximum value of r plotted (before applying rscale and roff)
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    zz = iota(0,n)*(zmax - zmin)/n + zmin
    rrin = ones(n+1,'d')
    rrout = ones(n+1,'d')
    for i in range(n+1):
        f3d.srfrv_z = zz[i]
        srfrvin()
        rrin[i] = f3d.srfrv_r
        srfrvout()
        rrout[i] = f3d.srfrv_r
    if gridframe:
        zz = (zz - w3d.zmmin)/w3d.dz
        rrin = (rrin)/w3d.dx + ir_axis
        rrout = (rrout)/w3d.dx + ir_axis
    rr = array(list(rrin) + list(rrout[::-1]))
    zz = array(list(zz) + list(zz[::-1]))
    rr = where(less(rr,rmin),rmin,rr)
    rr = where(greater(rr,rmax),rmax,rr)
    if outline:
        plg(rscale*rr+roff,zscale*zz+zoff,color=color)
    if fillcolor is not None:
        cc = array([fillcolor]).astype(ubyte)
        plfp(cc,rscale*rr+roff,zscale*zz+zoff,[2*n+2])

#####################################################################
def plotsrfrvin(srfrvin,zmin,zmax,n=1000,color='fg',gridframe=0,
                rscale=1,zscale=1,
                roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0,
                outline=1,fillcolor=None):
    """Handy function for plotting the r versus z for a surface of revolution
   - srfrvin: surface of revolution function to plot
   - zmin,zmax: z range to plot
   - n=1000: number of points to plot
   - color='fg': color of line
   - gridframe=0: when true, plots in grid frame
   - rscale=1: scaling for radius
   - zscale=1: scaling for z
   - roff=0: offset for radius
   - zoff=0: offset for z
   - rmin=0: minimum value of r plotted (before applying rscale and roff)
   - rmax=0: maximum value of r plotted (before applying rscale and roff)
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    zz = iota(0,n)*(zmax - zmin)/n + zmin
    rr = ones(n+1,'d')
    for i in range(n+1):
        f3d.srfrv_z = zz[i]
        srfrvin()
        rr[i] = f3d.srfrv_r
    if gridframe:
        zz = (zz - w3d.zmmin)/w3d.dz
        rr = (rr)/w3d.dx + ir_axis
    rr = where(less(rr,rmin),rmin,rr)
    rr = where(greater(rr,rmax),rmax,rr)
    rr = rscale*rr+roff
    zz = zscale*zz+zoff
    if outline:
        plg(rr,zz,color=color)
    if fillcolor is not None:
        cc = array([fillcolor]).astype(ubyte)
        rr = array(list(rr) + [rmin,rmin])
        zz = array(list(zz) + [zmax,zmin])
        plfp(cc,rr,zz,[n+3])

#####################################################################
def plotsrfrvout(srfrvin,zmin,zmax,n=1000,color='fg',gridframe=0,
                 rscale=1,zscale=1,
                 roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0,
                 outline=1,fillcolor=None):
    """Handy function for plotting the r versus z for a surface of revolution
   - srfrvout: surface of revolution function to plot
   - zmin,zmax: z range to plot
   - n=1000: number of points to plot
   - color='fg': color of line
   - gridframe=0: when true, plots in grid frame
   - rscale=1: scaling for radius
   - zscale=1: scaling for z
   - roff=0: offset for radius
   - zoff=0: offset for z
   - rmin=0: minimum value of r plotted (before applying rscale and roff)
   - rmax=0: maximum value of r plotted (before applying rscale and roff)
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    zz = iota(0,n)*(zmax - zmin)/n + zmin
    rr = ones(n+1,'d')
    for i in range(n+1):
        f3d.srfrv_z = zz[i]
        srfrvout()
        rr[i] = f3d.srfrv_r
    if gridframe:
        zz = (zz - w3d.zmmin)/w3d.dz
        rr = (rr)/w3d.dx + ir_axis
    rr = where(less(rr,rmin),rmin,rr)
    rr = where(greater(rr,rmax),rmax,rr)
    rr = rscale*rr+roff
    zz = zscale*zz+zoff
    if outline:
        plg(rr,zz,color=color)
    if fillcolor is not None:
        cc = array([fillcolor]).astype(ubyte)
        rmax = min(rmax,w3d.xmmax)
        rmax = max(rmax,max(rr))
        rr = array(list(rr) + [rmax,rmax])
        zz = array(list(zz) + [zmax,zmin])
        plfp(cc,rr,zz,[n+3])

#####################################################################
def plotsrfrv(srfrv,zmin,zmax,n=1000,color='fg',gridframe=0,rscale=1,zscale=1,
              roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0,
              outline=1,fillcolor=None):
    """Handy function for plotting the r versus z for a surface of revolution
   - srfrv: surface of revolution function to plot
   - zmin,zmax: z range to plot
   - n=1000: number of points to plot
   - color='fg': color of line
   - gridframe=0: when true, plots in grid frame
   - rscale=1: scaling for radius
   - zscale=1: scaling for z
   - roff=0: offset for radius
   - zoff=0: offset for z
   - rmin=0: minimum value of r plotted (before applying rscale and roff)
   - rmax=0: maximum value of r plotted (before applying rscale and roff)
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    zz = iota(0,n)*(zmax - zmin)/n + zmin
    rr = ones(n+1,'d')
    for i in range(n+1):
        f3d.srfrv_z = zz[i]
        srfrv()
        rr[i] = f3d.srfrv_r
    if gridframe:
        zz = (zz - w3d.zmmin)/w3d.dz
        rr = (rr)/w3d.dx + ir_axis
    rr = where(less(rr,rmin),rmin,rr)
    rr = where(greater(rr,rmax),rmax,rr)
    if outline:
        plg(rscale*rr+roff,zscale*zz+zoff,color=color)
    if fillcolor is not None:
        cc = array([fillcolor]).astype(ubyte)
        plfp(cc,rscale*rr+roff,zscale*zz+zoff,[n+1])


#####################################################################
#####################################################################
def plotelementoutline(color,gridframe,axis,zl,zu,ie,ne,outline,fillcolor,
                       ezs,eze,eap,eax,eay,eox,eoy,
                       err=None,erl=None,egl=None,egp=None,
                       epa=None,epr=None,epw=None,dpal=None,dpar=None,zoffset=None):
    """Plots the outline of electrostatic elements
    - color: line color
    - gridframe: when true, make plot in grid coordinates
    - axis: selects axis to plot, either 'x' or 'y'
    """
    if zu is None and zl is None:
        zl = w3d.zmmin + top.zbeam
        zu = w3d.zmmax + top.zbeam
    if zu is not None and zoffset is None and top.zlatperi > 0.:
        # --- This allows multiple periodic repeats of the lattice to be
        # --- plotted.
        if zl is None: zl = top.zlatstrt
        zoffset = floor((zl-top.zlatstrt)/top.zlatperi)*top.zlatperi
        while zoffset < zu:
            z1 = max(zl,zoffset) - zoffset
            z2 = min(zu,zoffset+top.zlatperi) - zoffset
            plotelementoutline(color,gridframe,axis,z1,z2,ie,ne,outline,fillcolor,
                               ezs,eze,eap,eax,eay,eox,eoy,err,erl,egl,egp,
                               epa,epr,epw,dpal,dpar,zoffset=zoffset)
            zoffset = zoffset + top.zlatperi
        return
    # --- Set the aperture radius.
    if max(eax) == 0.: eax = eap
    if max(eay) == 0.: eay = eap
    if axis == 'x': eap = eax
    else:           eap = eay
    # --- Set defaults for z-range
    if zl is None: zl = -largepos
    if zu is None: zu = +largepos
    #if zoffset is None: zoffset = top.zlatstrt
    if zoffset is None: zoffset = 0.
    if axis == 'x': gpsign = 1
    else:           gpsign = -1
    for i in range(ie,ie+ne):
        # --- Skip if outside the z-range
        if (eze[i]+top.zlatstrt+zoffset < zl or
            ezs[i]+top.zlatstrt+zoffset > zu): continue
        # --- plot rods
        # --- If aperture is zero, then this quad is skipped
        rodap = eap[i]
        if erl is not None and erl[i] > 0.:
            rodlen = erl[i]
            gp = egp[i]*gpsign
            gaplen = egl[i]
        else:
            rodlen = (eze[i] - ezs[i])
            gp = 1*gpsign
            gaplen = 0.
        if err is not None and err[i] > 0.: rodrr = err[i]
        else:                               rodrr = 8./7.*eap[i]
        if axis == 'x': offset = eox[i]
        else:           offset = eoy[i]
        if rodap > 0. and rodlen > 0.:
            rr = rodap + rodrr + rodrr*array([1.,1.,-1.,-1.,1.])
            zz = gp*(-0.5*(rodlen+gaplen) + rodlen*array([0.,1.,1.,0.,0.]))
            rr1 = offset + rr
            rr2 = offset - rr
            zz = 0.5*(eze[i] + ezs[i]) + top.zlatstrt + zoffset + zz
            if gridframe:
                rr1 = rr1/w3d.dx
                rr2 = rr2/w3d.dx
                zz = (zz - w3d.zmmin)/w3d.dz
            if outline:
                plg(rr1,zz,color=color)
                plg(rr2,zz,color=color)
            if fillcolor is not None:
                cc = array([fillcolor]).astype(ubyte)
                plfp(cc,rr1,zz,[5])
                plfp(cc,rr2,zz,[5])
        # --- Plot end plates
        if epw is not None and epw[i] > 0.:
            pw = epw[i]
            if epa[i] > 0.: pa = epa[i]
            else:           pa = eap[i]
            if epr[i] > 0.: pr = epr[i]
            else:           pr = rodap + 2.*rodrr
            pal = pa + dpal[i]
            par = pa + dpar[i]
            rrl = array([pr,pr,pal,pal,pr])
            rrr = array([pr,pr,par,par,pr])
            zz = pw*array([0.,1.,1.,0.,0.])
            rrl1 = offset + rrl
            rrl2 = offset - rrl
            rrr1 = offset + rrr
            rrr2 = offset - rrr
            zzl = (0.5*(eze[i] + ezs[i]) - 0.5*(rodlen+gaplen) - zz +
                  top.zlatstrt + zoffset)
            zzr = (0.5*(eze[i] + ezs[i]) + 0.5*(rodlen+gaplen) + zz +
                  top.zlatstrt + zoffset)
            if gridframe:
                rrl1 = rrl1/w3d.dx
                rrl2 = rrl2/w3d.dx
                rrr1 = rrr1/w3d.dx
                rrr2 = rrr2/w3d.dx
                zzl = (zzl - w3d.zmmin)/w3d.dz
                zzr = (zzr - w3d.zmmin)/w3d.dz
            if outline:
                plg(rrl1,zzl,color=color)
                plg(rrl2,zzl,color=color)
                plg(rrr1,zzr,color=color)
                plg(rrr2,zzr,color=color)
            if fillcolor is not None:
                cc = array([fillcolor]).astype(ubyte)
                plfp(cc,rrl1,zzl,[5])
                plfp(cc,rrl2,zzl,[5])
                plfp(cc,rrr1,zzr,[5])
                plfp(cc,rrr2,zzr,[5])


#---------------------------------------------------------------------------
def plotquadoutline(zl=None,zu=None,iq=0,nq=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
    """Plots the outline of quadrupole elements
   - zl,zu: range in lab frame to plot, defaults to extent of phi field.
   - iq=0: starting quad to plot
   - nq=top.nquad+1: number of quads to plot
   - color='fg': line color
   - gridframe=0: when true, make plot in grid coordinates
   - axis='x': selects axis to plot, either 'x' or 'y'
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    if top.nquad < 0: return
    if nq is None: nq = top.nquad + 1
    plotelementoutline(color,gridframe,axis,zl,zu,iq,nq,outline,fillcolor,
                       top.quadzs,top.quadze,top.quadap,top.quadax,top.quaday,
                       top.qoffx,top.qoffy,
                       top.quadrr,top.quadrl,top.quadgl,top.quadgp,
                       top.quadpa,top.quadpr,top.quadpw,
                       top.qdelpal,top.qdelpar)

#---------------------------------------------------------------------------
def plotheleoutline(zl=None,zu=None,ih=0,nh=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
    """Plots the outline of hele elements
   - zl,zu: range in lab frame to plot, defaults to extent of phi field.
   - ih=0: starting hele to plot
   - nh=top.nhele+1: number of heles to plot
   - color='fg': line color
   - gridframe=0: when true, make plot in grid coordinates
   - axis='x': selects axis to plot, either 'x' or 'y'
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    if top.nhele < 0: return
    if nh is None: nh = top.nhele + 1
    plotelementoutline(color,gridframe,axis,zl,zu,ih,nh,outline,fillcolor,
                       top.helezs,top.heleze,top.heleap,top.heleax,top.heleay,
                       top.heleox,top.heleoy,
                       top.helerr,top.helerl,top.helegl,top.helegp,
                       top.helepa,zeros(top.nhele+1,'d'),top.helepw,
                       zeros(top.nhele+1,'d'),zeros(top.nhele+1,'d'))

#---------------------------------------------------------------------------
def plotemltoutline(zl=None,zu=None,ie=0,ne=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
    """Plots the outline of emlt elements
   - zl,zu: range in lab frame to plot, defaults to extent of phi field.
   - ie=0: starting emlt to plot
   - ne=top.nemlt+1: number of emlts to plot
   - color='fg': line color
   - gridframe=0: when true, make plot in grid coordinates
   - axis='x': selects axis to plot, either 'x' or 'y'
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    if top.nemlt < 0: return
    if ne is None: ne = top.nemlt + 1
    plotelementoutline(color,gridframe,axis,zl,zu,ie,ne,outline,fillcolor,
                       top.emltzs,top.emltze,top.emltap,top.emltax,top.emltay,
                       top.emltox,top.emltoy,
                       top.emltrr,top.emltrl,top.emltgl,top.emltgp,
                       top.emltpa,zeros(top.nemlt+1,'d'),top.emltpw,
                       zeros(top.nemlt+1,'d'),zeros(top.nemlt+1,'d'))

#---------------------------------------------------------------------------
def plotmmltoutline(zl=None,zu=None,ie=0,nm=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
    """Plots the outline of mmlt elements
   - zl,zu: range in lab frame to plot, defaults to extent of phi field.
   - ie=0: starting mmlt to plot
   - nm=top.nmmlt+1: number of mmlts to plot
   - color='fg': line color
   - gridframe=0: when true, make plot in grid coordinates
   - axis='x': selects axis to plot, either 'x' or 'y'
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    if top.nmmlt < 0: return
    if nm is None: nm = top.nmmlt + 1
    plotelementoutline(color,gridframe,axis,zl,zu,ie,nm,outline,fillcolor,
                       top.mmltzs,top.mmltze,top.mmltap,top.mmltax,top.mmltay,
                       top.mmltox,top.mmltoy)

#---------------------------------------------------------------------------
def plotpgrdoutline(zl=None,zu=None,ip=0,np=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
    """Plots the outline of pgrd elements
   - zl,zu: range in lab frame to plot, defaults to extent of phi field.
   - ip=0: starting pgrd to plot
   - np=top.npgrd+1: number of pgrds to plot
   - color='fg': line color
   - gridframe=0: when true, make plot in grid coordinates
   - axis='x': selects axis to plot, either 'x' or 'y'
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    if top.npgrd < 0: return
    if np is None: np = top.npgrd + 1
    plotelementoutline(color,gridframe,axis,zl,zu,ip,np,outline,fillcolor,
                       top.pgrdzs,top.pgrdze,top.pgrdap,top.pgrdax,top.pgrday,
                       top.pgrdox,top.pgrdoy,
                       top.pgrdrr,top.pgrdrl,top.pgrdgl,top.pgrdgp,
                       top.pgrdpa,zeros(top.npgrd+1,'d'),top.pgrdpw,
                       zeros(top.npgrd+1,'d'),zeros(top.npgrd+1,'d'))

#---------------------------------------------------------------------------
def plotaccloutline(zl=None,zu=None,ia=0,na=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
    """Plots the outline of accl elements
   - zl,zu: range in lab frame to plot, defaults to extent of phi field.
   - ia=0: starting accl to plot
   - na=top.naccl+1: number of accls to plot
   - color='fg': line color
   - gridframe=0: when true, make plot in grid coordinates
   - axis='x': selects axis to plot, either 'x' or 'y'
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    if top.naccl < 0: return
    if na is None: na = top.naccl + 1
    plotelementoutline(color,gridframe,axis,zl,zu,ia,na,outline,fillcolor,
                       top.acclzs,top.acclze,top.acclap,top.acclax,top.acclay,
                       top.acclox,top.accloy)

#---------------------------------------------------------------------------
def plotdrftoutline(zl=None,zu=None,id=0,nd=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
    """Plots the outline of drft elements
   - zl,zu: range in lab frame to plot, defaults to extent of phi field.
   - id=0: starting drft to plot
   - nd=top.ndrft+1: number of drfts to plot
   - color='fg': line color
   - gridframe=0: when true, make plot in grid coordinates
   - axis='x': selects axis to plot, either 'x' or 'y'
   - outline=1: when true, draw outline
   - fillcolor=None: optionally sets fill color
    """
    if top.ndrft < 0: return
    if nd is None: nd = top.ndrft + 1
    plotelementoutline(color,gridframe,axis,zl,zu,id,nd,outline,fillcolor,
                       top.drftzs,top.drftze,top.drftap,top.drftax,top.drftay,
                       top.drftox,top.drftoy)

#########################################################################
#########################################################################
#########################################################################
def updatemgconductors(dumpfilename):
    """
  This routine updates conductors from old dump files that were written
  before the conductor data was put into derived types.
   - dumpfilename: name of old dump file
    """
    ff = PR.PR(dumpfilename)

    f3d.conductors.interior.nmax = ff.read('ncondmax@f3d')
    f3d.conductors.evensubgrid.nmax = ff.read('ncndmax@f3d')
    f3d.conductors.oddsubgrid.nmax = ff.read('ncndmax@f3d')
    gchange('Conductor3d')

    f3d.conductors.interior.n = ff.read('ncond@f3d')
    if f3d.conductors.interior.n > 0:
        f3d.conductors.interior.indx[0,:] = ff.read('ixcond@f3d')
        f3d.conductors.interior.indx[1,:] = ff.read('iycond@f3d')
        f3d.conductors.interior.indx[2,:] = ff.read('izcond@f3d')
        f3d.conductors.interior.volt[:] = ff.read('condvolt@f3d')
        f3d.conductors.interior.numb[:] = ff.read('condnumb@f3d')
        f3d.conductors.interior.ilevel[:] = ff.read('icondlevel@f3d')

    f3d.conductors.evensubgrid.n = ff.read('necndbdy@f3d')
    if f3d.conductors.evensubgrid.n > 0:
        f3d.conductors.evensubgrid.indx[0,:] = ff.read('iecndx@f3d')
        f3d.conductors.evensubgrid.indx[1,:] = ff.read('iecndy@f3d')
        f3d.conductors.evensubgrid.indx[2,:] = ff.read('iecndz@f3d')
        f3d.conductors.evensubgrid.dels[0,:] = ff.read('ecdelmx@f3d')
        f3d.conductors.evensubgrid.dels[1,:] = ff.read('ecdelpx@f3d')
        f3d.conductors.evensubgrid.dels[2,:] = ff.read('ecdelmy@f3d')
        f3d.conductors.evensubgrid.dels[3,:] = ff.read('ecdelpy@f3d')
        f3d.conductors.evensubgrid.dels[4,:] = ff.read('ecdelmz@f3d')
        f3d.conductors.evensubgrid.dels[5,:] = ff.read('ecdelpz@f3d')
        f3d.conductors.evensubgrid.volt[0,:] = ff.read('ecvoltmx@f3d')
        f3d.conductors.evensubgrid.volt[1,:] = ff.read('ecvoltpx@f3d')
        f3d.conductors.evensubgrid.volt[2,:] = ff.read('ecvoltmy@f3d')
        f3d.conductors.evensubgrid.volt[3,:] = ff.read('ecvoltpy@f3d')
        f3d.conductors.evensubgrid.volt[4,:] = ff.read('ecvoltmz@f3d')
        f3d.conductors.evensubgrid.volt[5,:] = ff.read('ecvoltpz@f3d')
        f3d.conductors.evensubgrid.numb[0,:] = ff.read('ecnumbmx@f3d')
        f3d.conductors.evensubgrid.numb[1,:] = ff.read('ecnumbpx@f3d')
        f3d.conductors.evensubgrid.numb[2,:] = ff.read('ecnumbmy@f3d')
        f3d.conductors.evensubgrid.numb[3,:] = ff.read('ecnumbpy@f3d')
        f3d.conductors.evensubgrid.numb[4,:] = ff.read('ecnumbmz@f3d')
        f3d.conductors.evensubgrid.numb[5,:] = ff.read('ecnumbpz@f3d')
        f3d.conductors.evensubgrid.ilevel[:] = ff.read('iecndlevel@f3d')

    f3d.conductors.oddsubgrid.n = ff.read('nocndbdy@f3d')
    if f3d.conductors.oddsubgrid.n > 0:
        f3d.conductors.oddsubgrid.indx[0,:] = ff.read('iocndx@f3d')
        f3d.conductors.oddsubgrid.indx[1,:] = ff.read('iocndy@f3d')
        f3d.conductors.oddsubgrid.indx[2,:] = ff.read('iocndz@f3d')
        f3d.conductors.oddsubgrid.dels[0,:] = ff.read('ocdelmx@f3d')
        f3d.conductors.oddsubgrid.dels[1,:] = ff.read('ocdelpx@f3d')
        f3d.conductors.oddsubgrid.dels[2,:] = ff.read('ocdelmy@f3d')
        f3d.conductors.oddsubgrid.dels[3,:] = ff.read('ocdelpy@f3d')
        f3d.conductors.oddsubgrid.dels[4,:] = ff.read('ocdelmz@f3d')
        f3d.conductors.oddsubgrid.dels[5,:] = ff.read('ocdelpz@f3d')
        f3d.conductors.oddsubgrid.volt[0,:] = ff.read('ocvoltmx@f3d')
        f3d.conductors.oddsubgrid.volt[1,:] = ff.read('ocvoltpx@f3d')
        f3d.conductors.oddsubgrid.volt[2,:] = ff.read('ocvoltmy@f3d')
        f3d.conductors.oddsubgrid.volt[3,:] = ff.read('ocvoltpy@f3d')
        f3d.conductors.oddsubgrid.volt[4,:] = ff.read('ocvoltmz@f3d')
        f3d.conductors.oddsubgrid.volt[5,:] = ff.read('ocvoltpz@f3d')
        f3d.conductors.oddsubgrid.numb[0,:] = ff.read('ocnumbmx@f3d')
        f3d.conductors.oddsubgrid.numb[1,:] = ff.read('ocnumbpx@f3d')
        f3d.conductors.oddsubgrid.numb[2,:] = ff.read('ocnumbmy@f3d')
        f3d.conductors.oddsubgrid.numb[3,:] = ff.read('ocnumbpy@f3d')
        f3d.conductors.oddsubgrid.numb[4,:] = ff.read('ocnumbmz@f3d')
        f3d.conductors.oddsubgrid.numb[5,:] = ff.read('ocnumbpz@f3d')
        f3d.conductors.oddsubgrid.ilevel[:] = ff.read('iocndlevel@f3d')

    ff.close()

def updatemgconductorsold():
    """
  This routine updates conductors which from old dump files.
  of the code, the conductor coordinates and deltas were stored relative to
  the finest grid. In the current version, the data is stored relative the to
  grid that data is to be used for.
    """
    # --- Check if one of the MG conductor arrays is allocated. If not, then
    # --- return since the conductors have not been generated for the MG solver.
    # --- This check is used instead of fstype==7 since when the field solve
    # --- is done periodically, fstype may be -1 in the restart dump.
    test = f3d.getpyobject("ecdelmx")
    if test is None: return

    # --- First, copy all of the data
    ixcond = f3d.ixcond[0:f3d.ncond] + 0.
    iycond = f3d.iycond[0:f3d.ncond] + 0.
    izcond = f3d.izcond[0:f3d.ncond] + 0.
    condvolt = f3d.condvolt[0:f3d.ncond] + 0.
    condnumb = f3d.condnumb[0:f3d.ncond] + 0.
    icondlxy = f3d.icondlxy[0:f3d.ncond] + 0
    icondlz  = f3d.icondlz[0:f3d.ncond] + 0
    iecndx = f3d.iecndx[0:f3d.necndbdy] + 0.
    iecndy = f3d.iecndy[0:f3d.necndbdy] + 0.
    iecndz = f3d.iecndz[0:f3d.necndbdy] + 0.
    ecdelmx = f3d.ecdelmx[0:f3d.necndbdy] + 0.
    ecdelpx = f3d.ecdelpx[0:f3d.necndbdy] + 0.
    ecdelmy = f3d.ecdelmy[0:f3d.necndbdy] + 0.
    ecdelpy = f3d.ecdelpy[0:f3d.necndbdy] + 0.
    ecdelmz = f3d.ecdelmz[0:f3d.necndbdy] + 0.
    ecdelpz = f3d.ecdelpz[0:f3d.necndbdy] + 0.
    ecvolt = f3d.ecvolt[0:f3d.necndbdy] + 0.
    ecnumb = f3d.ecnumb[0:f3d.necndbdy] + 0.
    ecvoltmx = f3d.ecvoltmx[0:f3d.necndbdy] + 0.
    ecvoltpx = f3d.ecvoltpx[0:f3d.necndbdy] + 0.
    ecvoltmy = f3d.ecvoltmy[0:f3d.necndbdy] + 0.
    ecvoltpy = f3d.ecvoltpy[0:f3d.necndbdy] + 0.
    ecvoltmz = f3d.ecvoltmz[0:f3d.necndbdy] + 0.
    ecvoltpz = f3d.ecvoltpz[0:f3d.necndbdy] + 0.
    ecnumbmx = f3d.ecnumbmx[0:f3d.necndbdy] + 0.
    ecnumbpx = f3d.ecnumbpx[0:f3d.necndbdy] + 0.
    ecnumbmy = f3d.ecnumbmy[0:f3d.necndbdy] + 0.
    ecnumbpy = f3d.ecnumbpy[0:f3d.necndbdy] + 0.
    ecnumbmz = f3d.ecnumbmz[0:f3d.necndbdy] + 0.
    ecnumbpz = f3d.ecnumbpz[0:f3d.necndbdy] + 0.
    iecndlxy = f3d.iecndlxy[0:f3d.necndbdy] + 0.
    iecndlz = f3d.iecndlz[0:f3d.necndbdy] + 0.
    iocndx = f3d.iocndx[0:f3d.nocndbdy] + 0.
    iocndy = f3d.iocndy[0:f3d.nocndbdy] + 0.
    iocndz = f3d.iocndz[0:f3d.nocndbdy] + 0.
    ocdelmx = f3d.ocdelmx[0:f3d.nocndbdy] + 0.
    ocdelpx = f3d.ocdelpx[0:f3d.nocndbdy] + 0.
    ocdelmy = f3d.ocdelmy[0:f3d.nocndbdy] + 0.
    ocdelpy = f3d.ocdelpy[0:f3d.nocndbdy] + 0.
    ocdelmz = f3d.ocdelmz[0:f3d.nocndbdy] + 0.
    ocdelpz = f3d.ocdelpz[0:f3d.nocndbdy] + 0.
    ocvolt = f3d.ocvolt[0:f3d.nocndbdy] + 0.
    ocnumb = f3d.ocnumb[0:f3d.nocndbdy] + 0.
    ocvoltmx = f3d.ocvoltmx[0:f3d.nocndbdy] + 0.
    ocvoltpx = f3d.ocvoltpx[0:f3d.nocndbdy] + 0.
    ocvoltmy = f3d.ocvoltmy[0:f3d.nocndbdy] + 0.
    ocvoltpy = f3d.ocvoltpy[0:f3d.nocndbdy] + 0.
    ocvoltmz = f3d.ocvoltmz[0:f3d.nocndbdy] + 0.
    ocvoltpz = f3d.ocvoltpz[0:f3d.nocndbdy] + 0.
    ocnumbmx = f3d.ocnumbmx[0:f3d.nocndbdy] + 0.
    ocnumbpx = f3d.ocnumbpx[0:f3d.nocndbdy] + 0.
    ocnumbmy = f3d.ocnumbmy[0:f3d.nocndbdy] + 0.
    ocnumbpy = f3d.ocnumbpy[0:f3d.nocndbdy] + 0.
    ocnumbmz = f3d.ocnumbmz[0:f3d.nocndbdy] + 0.
    ocnumbpz = f3d.ocnumbpz[0:f3d.nocndbdy] + 0.
    iocndlxy = f3d.iocndlxy[0:f3d.nocndbdy] + 0.
    iocndlz = f3d.iocndlz[0:f3d.nocndbdy] + 0.

    # --- Get the list of coarsening levels
    llxy = findunique(icondlxy)
    llz  = findunique(icondlz)

    # --- Create lists to hold the converted data.
    ixcondnew = []
    iycondnew = []
    izcondnew = []
    condvoltnew = []
    icondlxynew = []
    icondlznew = []
    icndxnew = []
    icndynew = []
    icndznew = []
    cdelmxnew = []
    cdelmynew = []
    cdelmznew = []
    cdelpxnew = []
    cdelpynew = []
    cdelpznew = []
    cvoltnew = []
    cvoltmxnew = []
    cvoltpxnew = []
    cvoltmynew = []
    cvoltpynew = []
    cvoltmznew = []
    cvoltpznew = []
    icndlxynew = []
    icndlznew = []

    # --- Get number of old data points
    ncond = len(ixcond)
    necndbdy = len(iecndx)
    nocndbdy = len(iocndx)

    # --- Loop over coarsening levels, collecting and converting the data
    # --- for each level.
    for j in range(len(llxy)):
        ii = compress(logical_and(icondlxy >= llxy[j],icondlz >= llz[j]), \
                      arange(ncond))
        ixcondnew = ixcondnew + list(take(ixcond/llxy[j],ii))
        iycondnew = iycondnew + list(take(iycond/llxy[j],ii))
        izcondnew = izcondnew + list(take(izcond/llz[j],ii))
        condvoltnew = condvoltnew + list(take(condvolt,ii))
        icondlxynew = icondlxynew + list(llxy[j]*ones(len(ii),'l'))
        icondlznew = icondlznew + list(llz[j]*ones(len(ii),'l'))

        ii = compress(logical_and(iecndlxy >= llxy[j],iecndlz >= llz[j]), \
                      arange(necndbdy))
        icndxnew = icndxnew + list(take(iecndx/llxy[j],ii))
        icndynew = icndynew + list(take(iecndy/llxy[j],ii))
        icndznew = icndznew + list(take(iecndz/llz[j],ii))
        cdelmxnew = cdelmxnew + list(take(ecdelmx/llxy[j],ii))
        cdelmynew = cdelmynew + list(take(ecdelmy/llxy[j],ii))
        cdelmznew = cdelmznew + list(take(ecdelmz/llz[j],ii))
        cdelpxnew = cdelpxnew + list(take(ecdelpx/llxy[j],ii))
        cdelpynew = cdelpynew + list(take(ecdelpy/llxy[j],ii))
        cdelpznew = cdelpznew + list(take(ecdelpz/llz[j],ii))
        cvoltnew = cvoltnew + list(take(ecvolt,ii))
        cvoltmxnew = cvoltmxnew + list(take(ecvoltmx,ii))
        cvoltpxnew = cvoltpxnew + list(take(ecvoltpx,ii))
        cvoltmynew = cvoltmynew + list(take(ecvoltmy,ii))
        cvoltpynew = cvoltpynew + list(take(ecvoltpy,ii))
        cvoltmznew = cvoltmznew + list(take(ecvoltmz,ii))
        cvoltpznew = cvoltpznew + list(take(ecvoltpz,ii))
        icndlxynew = icndlxynew + list(llxy[j]*ones(len(ii),'l'))
        icndlznew = icndlznew + list(llz[j]*ones(len(ii),'l'))

        ii = compress(logical_and(iocndlxy >= llxy[j],iocndlz >= llz[j]), \
                      arange(nocndbdy))
        icndxnew = icndxnew + list(take(iocndx/llxy[j],ii))
        icndynew = icndynew + list(take(iocndy/llxy[j],ii))
        icndznew = icndznew + list(take(iocndz/llz[j],ii))
        cdelmxnew = cdelmxnew + list(take(ocdelmx/llxy[j],ii))
        cdelmynew = cdelmynew + list(take(ocdelmy/llxy[j],ii))
        cdelmznew = cdelmznew + list(take(ocdelmz/llz[j],ii))
        cdelpxnew = cdelpxnew + list(take(ocdelpx/llxy[j],ii))
        cdelpynew = cdelpynew + list(take(ocdelpy/llxy[j],ii))
        cdelpznew = cdelpznew + list(take(ocdelpz/llz[j],ii))
        cvoltnew = cvoltnew + list(take(ocvolt,ii))
        cvoltmxnew = cvoltmxnew + list(take(ocvoltmx,ii))
        cvoltpxnew = cvoltpxnew + list(take(ocvoltpx,ii))
        cvoltmynew = cvoltmynew + list(take(ocvoltmy,ii))
        cvoltpynew = cvoltpynew + list(take(ocvoltpy,ii))
        cvoltmznew = cvoltmznew + list(take(ocvoltmz,ii))
        cvoltpznew = cvoltpznew + list(take(ocvoltpz,ii))
        icndlxynew = icndlxynew + list(llxy[j]*ones(len(ii),'l'))
        icndlznew = icndlznew + list(llz[j]*ones(len(ii),'l'))

    # --- Find odd and even points
    icndxnew = array(icndxnew,dtype='i')
    icndynew = array(icndynew,dtype='i')
    icndznew = array(icndznew,dtype='i')
    ee = compress((icndxnew+icndynew+icndznew)%2==0,arange(len(icndxnew)))
    oo = compress((icndxnew+icndynew+icndznew)%2==1,arange(len(icndxnew)))

    # --- Reset counters and reallocate arrays for the converted data.
    f3d.ncondmax = len(ixcondnew)
    f3d.ncndmax = max(len(ee),len(oo))
    gchange("Conductor3d")
    f3d.ncond = len(ixcondnew)
    f3d.necndbdy = len(ee)
    f3d.nocndbdy = len(oo)

    # --- Copy converted data into WARP arrays.
    f3d.ixcond[:f3d.ncond] = ixcondnew
    f3d.iycond[:f3d.ncond] = iycondnew
    f3d.izcond[:f3d.ncond] = izcondnew
    f3d.condvolt[:f3d.ncond] = condvoltnew
    f3d.icondlevel[:f3d.ncond] = nint(log(icondlxynew)/log(2.)+0.5)
    #f3d.icondlxy[:f3d.ncond] = icondlxynew
    #f3d.icondlz[:f3d.ncond] = icondlznew

    f3d.iecndx[:f3d.necndbdy] = take(icndxnew,ee)
    f3d.iecndy[:f3d.necndbdy] = take(icndynew,ee)
    f3d.iecndz[:f3d.necndbdy] = take(icndznew,ee)
    f3d.ecdelmx[:f3d.necndbdy] = take(cdelmxnew,ee)
    f3d.ecdelmy[:f3d.necndbdy] = take(cdelmynew,ee)
    f3d.ecdelmz[:f3d.necndbdy] = take(cdelmznew,ee)
    f3d.ecdelpx[:f3d.necndbdy] = take(cdelpxnew,ee)
    f3d.ecdelpy[:f3d.necndbdy] = take(cdelpynew,ee)
    f3d.ecdelpz[:f3d.necndbdy] = take(cdelpznew,ee)
    f3d.ecvolt[:f3d.necndbdy] = take(cvoltnew,ee)
    f3d.ecvoltmx[:f3d.necndbdy] = take(cvoltmxnew,ee)
    f3d.ecvoltpx[:f3d.necndbdy] = take(cvoltpxnew,ee)
    f3d.ecvoltmy[:f3d.necndbdy] = take(cvoltmynew,ee)
    f3d.ecvoltpy[:f3d.necndbdy] = take(cvoltpynew,ee)
    f3d.ecvoltmz[:f3d.necndbdy] = take(cvoltmznew,ee)
    f3d.ecvoltpz[:f3d.necndbdy] = take(cvoltpznew,ee)
    f3d.iecndlevel[:f3d.ncond] = nint(log(iecndlxynew)/log(2.)+0.5)
    #f3d.iecndlxy[:f3d.necndbdy] = take(icndlxynew,ee)
    #f3d.iecndlz[:f3d.necndbdy] = take(icndlznew,ee)

    f3d.iocndx[:f3d.nocndbdy] = take(icndxnew,oo)
    f3d.iocndy[:f3d.nocndbdy] = take(icndynew,oo)
    f3d.iocndz[:f3d.nocndbdy] = take(icndznew,oo)
    f3d.ocdelmx[:f3d.nocndbdy] = take(cdelmxnew,oo)
    f3d.ocdelmy[:f3d.nocndbdy] = take(cdelmynew,oo)
    f3d.ocdelmz[:f3d.nocndbdy] = take(cdelmznew,oo)
    f3d.ocdelpx[:f3d.nocndbdy] = take(cdelpxnew,oo)
    f3d.ocdelpy[:f3d.nocndbdy] = take(cdelpynew,oo)
    f3d.ocdelpz[:f3d.nocndbdy] = take(cdelpznew,oo)
    f3d.ocvolt[:f3d.nocndbdy] = take(cvoltnew,oo)
    f3d.ocvoltmx[:f3d.nocndbdy] = take(cvoltmxnew,oo)
    f3d.ocvoltpx[:f3d.nocndbdy] = take(cvoltpxnew,oo)
    f3d.ocvoltmy[:f3d.nocndbdy] = take(cvoltmynew,oo)
    f3d.ocvoltpy[:f3d.nocndbdy] = take(cvoltpynew,oo)
    f3d.ocvoltmz[:f3d.nocndbdy] = take(cvoltmznew,oo)
    f3d.ocvoltpz[:f3d.nocndbdy] = take(cvoltpznew,oo)
    f3d.iocndlevel[:f3d.ncond] = nint(log(iocndlxynew)/log(2.)+0.5)
    #f3d.iocndlxy[:f3d.nocndbdy] = take(icndlxynew,oo)
    #f3d.iocndlz[:f3d.nocndbdy] = take(icndlznew,oo)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
# Convenience routines for calling the surface of revolution routines.
# The presence of these allows changing the argument lists of the fortran
# versions without breaking code.
#########################################################################
def srfrvout(rofzfunc=" ",volt=0.,zmin=None,zmax=None,xcent=0.,ycent=0.,
             rmax=top.largepos,lfill=false,
             xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
             zmminlocal=None,zmmaxlocal=None,zmmin=None,zbeam=None,
             dx=None,dy=None,dz=None,
             nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
             xmesh=None,ymesh=None,l2symtry=None,l4symtry=None,
             srfrv_pernz=0,condid=0):
    """
  Sets up a conductor represented by the outside of a surface of revolution.
  The routine rofzfunc should be of the form
     def rofz():
       f3d.srfrv_r = f(f3d.srfrv_z)
  where f() is the radius as a function of z and srfrv_z and srfrv_r are
  compiled variables which pass data into and out of rofz.
  The begining of the calling sequence would then be
  srfrvout(rofz,...)

  Input:
    rofzfunc=" ": routine which calculates the radius as function of z.
    volt=0.: voltage on the conductor.
    zmin=w3d.zmminlocal: minimum z of the conductor.
    zmax=w3d.zmmaxlocal: maximum z of the conductor.
    xcent=0.: x center of the conductor.
    ycent=0.: y center of the conductor.
    rmax=LargePos: maximum radius of the conductor.
    lfill=false: logical requesting that the whole conductor be filled
                 with points.
    lshell=true: logical requesting that the shell be subgrid resolved
    xmin,xmax,ymin,ymax: min and max transverse extent of conductor.
                         default to w3d.xmmin,xmmax,ymmin,ymmax
    zmminlocal,zmmaxlocal,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
    xmesh,ymesh,l2symtry,l4symtry:
             are all variables describing the grid. Default to variables in w3d
             and top with the same name.
    srfrv_pernz=0: when non-zero, a piece-wise linear approximation is made of
                   the surface (for optimization)
    condid=0: Id number to identify this conductor
  Output is put directly into the conductor arrays of Conductor3d.
    """
    print """
  Warning: srfrvout is obsolete and should no longer be used.
  use ZSrfrvOut instead.
    """
    if xmin is None: xmin = w3d.xmminlocal
    if xmax is None: xmax = w3d.xmmaxlocal
    if ymin is None: ymin = w3d.ymminlocal
    if ymax is None: ymax = w3d.ymmaxlocal
    if zmin is None: zmin = w3d.zmminlocal
    if zmax is None: zmax = w3d.zmmaxlocal
    if xmminlocal is None: xmminlocal = w3d.xmminlocal
    if xmmaxlocal is None: xmmaxlocal = w3d.xmmaxlocal
    if ymminlocal is None: ymminlocal = w3d.ymminlocal
    if ymmaxlocal is None: ymmaxlocal = w3d.ymmaxlocal
    if zmminlocal is None: zmminlocal = w3d.zmminlocal
    if zmmaxlocal is None: zmmaxlocal = w3d.zmmaxlocal
    if zmmin is None: zmmin = w3d.zmmin
    if zbeam is None: zbeam = top.zbeam
    if dx is None: dx = w3d.dx
    if dy is None: dy = w3d.dy
    if dz is None: dz = w3d.dz
    if nx is None: nx = w3d.nx
    if ny is None: ny = w3d.ny
    if nz is None: nz = w3d.nz
    if ix_axis is None: ix_axis = nint(-w3d.xmmin/w3d.dx)
    if iy_axis is None: iy_axis = nint(-w3d.ymmin/w3d.dy)
    if xmesh is None: xmesh = w3d.xmesh
    if ymesh is None: ymesh = w3d.ymesh
    if l2symtry is None: l2symtry = w3d.l2symtry
    if l4symtry is None: l4symtry = w3d.l4symtry
    if srfrv_pernz > 0:
        save_srfrv_pernz = f3d.srfrv_pernz
        f3d.srfrv_pernz = srfrv_pernz

    # --- Make sure the rofzfunc is in main.
    # --- Note that this can only really work if a reference to the function
    # --- is passed in (instead of the name).
    if not f3d.lsrlinr and callable(rofzfunc):
        __main__.__dict__[rofzfunc.__name__] = rofzfunc

    # --- Get the name of the input function if a reference to the function
    # --- was passed in.
    if callable(rofzfunc): rofzfunc = rofzfunc.__name__

    # --- Get srfrv function - in older versions of the code the final 'f'
    # --- was not there
    try:
        srfrvfunc = f3d.srfrvoutf
    except AttributeError:
        srfrvfunc = f3d.srfrvout

    # --- Now call the fortran version
    srfrvfunc(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
              xmin,xmax,ymin,ymax,lshell,zmminlocal,zmmaxlocal,zmmin,zbeam,dx,dy,dz,
              nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)

    # --- Reset srfrv_pernz if needed
    if srfrv_pernz > 0: f3d.srfrv_pernz = save_srfrv_pernz

#---------------------------------------------------------------------------
def srfrvin(rofzfunc=" ",volt=0.,zmin=None,zmax=None,xcent=0.,ycent=0.,
            rmin=0.,lfill=false,
            xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
            zmminlocal=None,zmmaxlocal=None,zmmin=None,zbeam=None,
            dx=None,dy=None,dz=None,
            nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
            xmesh=None,ymesh=None,l2symtry=None,l4symtry=None,
            srfrv_pernz=0,condid=0):
    """
  Sets up a conductor represented by the inside of a surface of revolution.
  The routine rofzfunc should be of the form
     def rofz():
       f3d.srfrv_r = f(f3d.srfrv_z)
  where f() is the radius as a function of z and srfrv_z and srfrv_r are
  compiled variables which pass data into and out of rofz.
  The begining of the calling sequence would then be
  srfrvout(rofz,...)

  Input:
    rofzfunc=" ": routine which calculates the radius as function of z.
    volt=0.: voltage on the conductor.
    zmin=w3d.zmminlocal: minimum z of the conductor.
    zmax=w3d.zmmaxlocal: maximum z of the conductor.
    xcent=0.: x center of the conductor.
    ycent=0.: y center of the conductor.
    rmin=0.: minimum radius of the conductor.
    lfill=false: logical requesting that the whole conductor be filled
                 with points.
    lshell=true: logical requesting that the shell be subgrid resolved
    xmin,xmax,ymin,ymax: min and max transverse extent of conductor.
                         default to w3d.xmmin,xmmax,ymmin,ymmax
    zmminlocal,zmmaxlocal,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
    xmesh,ymesh,l2symtry,l4symtry:
             are all variables describing the grid. Default to variables in w3d
             and top with the same name.
    srfrv_pernz=0: when non-zero, a piece-wise linear approximation is made of
                   the surface (for optimization)
    condid=0: Id number to identify this conductor
  Output is put directly into the conductor arrays of Conductor3d.
    """
    print """
  Warning: srfrvin is obsolete and should no longer be used.
  use ZSrfrvIn instead.
    """
    if xmin is None: xmin = w3d.xmminlocal
    if xmax is None: xmax = w3d.xmmaxlocal
    if ymin is None: ymin = w3d.ymminlocal
    if ymax is None: ymax = w3d.ymmaxlocal
    if zmin is None: zmin = w3d.zmminlocal
    if zmax is None: zmax = w3d.zmmaxlocal
    if xmminlocal is None: xmminlocal = w3d.xmminlocal
    if xmmaxlocal is None: xmmaxlocal = w3d.xmmaxlocal
    if ymminlocal is None: ymminlocal = w3d.ymminlocal
    if ymmaxlocal is None: ymmaxlocal = w3d.ymmaxlocal
    if zmminlocal is None: zmminlocal = w3d.zmminlocal
    if zmmaxlocal is None: zmmaxlocal = w3d.zmmaxlocal
    if zmmin is None: zmmin = w3d.zmmin
    if zbeam is None: zbeam = top.zbeam
    if dx is None: dx = w3d.dx
    if dy is None: dy = w3d.dy
    if dz is None: dz = w3d.dz
    if nx is None: nx = w3d.nx
    if ny is None: ny = w3d.ny
    if nz is None: nz = w3d.nz
    if ix_axis is None: ix_axis = nint(-w3d.xmmin/w3d.dx)
    if iy_axis is None: iy_axis = nint(-w3d.ymmin/w3d.dy)
    if xmesh is None: xmesh = w3d.xmesh
    if ymesh is None: ymesh = w3d.ymesh
    if l2symtry is None: l2symtry = w3d.l2symtry
    if l4symtry is None: l4symtry = w3d.l4symtry
    if srfrv_pernz > 0:
        save_srfrv_pernz = f3d.srfrv_pernz
        f3d.srfrv_pernz = srfrv_pernz

    # --- Make sure the rofzfunc is in main.
    # --- Note that this can only really work if a reference to the function
    # --- is passed in (instead of the name).
    if not f3d.lsrlinr and callable(rofzfunc):
        __main__.__dict__[rofzfunc.__name__] = rofzfunc

    # --- Get the name of the input function if a reference to the function
    # --- was passed in.
    if callable(rofzfunc): rofzfunc = rofzfunc.__name__

    # --- Get srfrv function - in older versions of the code the final 'f'
    # --- was not there
    try:
        srfrvfunc = f3d.srfrvinf
    except AttributeError:
        srfrvfunc = f3d.srfrvin

    # --- Now call the fortran version
    srfrvfunc(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
              xmin,xmax,ymin,ymax,lshell,zmminlocal,zmmaxlocal,zmmin,zbeam,dx,dy,dz,
              nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)

    # --- Reset srfrv_pernz if needed
    if srfrv_pernz > 0: f3d.srfrv_pernz = save_srfrv_pernz

#---------------------------------------------------------------------------
def srfrvinout(rminofz=" ",rmaxofz=" ",volt=0.,zmin=None,zmax=None,
               xcent=0.,ycent=0.,lzend=true,
               xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
               zmminlocal=None,zmmaxlocal=None,zmmin=None,zbeam=None,
               dx=None,dy=None,dz=None,
               nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
               xmesh=None,ymesh=None,l2symtry=None,l4symtry=None,
               srfrv_pernz=0,condid=0):
    """
  Sets up a conductor between two surfaces of revolution.
  The routines rminofz and rmaxofz should be of the form
     def rofz():
       f3d.srfrv_r = f(f3d.srfrv_z)
  where f() is the radius as a function of z and srfrv_z and srfrv_r are
  compiled variables which pass data into and out of rofz.

  Input:
    rminofz=" ": routine which calculates the inner radius as function of z.
    rmaxofz=" ": routine which calculates the outer radius as function of z.
    volt=0.: voltage on the conductor.
    zmin=w3d.zmminlocal: minimum z of the conductor.
    zmax=w3d.zmmaxlocal: maximum z of the conductor.
    xcent=0.: x center of the conductor.
    ycent=0.: y center of the conductor.
    rmin=0.: minimum radius of the conductor.
    lzend=true: logical requesting that the end of the conductor be included
    lshell=true: logical requesting that the shell be subgrid resolved
    xmin,xmax,ymin,ymax: min and max transverse extent of conductor.
                         default to w3d.xmmin,xmmax,ymmin,ymmax
    zmminlocal,zmmaxlocal,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
    xmesh,ymesh,l2symtry,l4symtry:
             are all variables describing the grid. Default to variables in w3d
             and top with the same name.
    srfrv_pernz=0: when non-zero, a piece-wise linear approximation is made of
                   the surface (for optimization)
    condid=0: Id number to identify this conductor
  Output is put directly into the conductor arrays of Conductor3d.
    """
    print """
  Warning: srfrvinout is obsolete and should no longer be used.
  use ZSrfrvInOut instead.
    """
    if xmin is None: xmin = w3d.xmminlocal
    if xmax is None: xmax = w3d.xmmaxlocal
    if ymin is None: ymin = w3d.ymminlocal
    if ymax is None: ymax = w3d.ymmaxlocal
    if zmin is None: zmin = w3d.zmminlocal
    if zmax is None: zmax = w3d.zmmaxlocal
    if xmminlocal is None: xmminlocal = w3d.xmminlocal
    if xmmaxlocal is None: xmmaxlocal = w3d.xmmaxlocal
    if ymminlocal is None: ymminlocal = w3d.ymminlocal
    if ymmaxlocal is None: ymmaxlocal = w3d.ymmaxlocal
    if zmminlocal is None: zmminlocal = w3d.zmminlocal
    if zmmaxlocal is None: zmmaxlocal = w3d.zmmaxlocal
    if zmmin is None: zmmin = w3d.zmmin
    if zbeam is None: zbeam = top.zbeam
    if dx is None: dx = w3d.dx
    if dy is None: dy = w3d.dy
    if dz is None: dz = w3d.dz
    if nx is None: nx = w3d.nx
    if ny is None: ny = w3d.ny
    if nz is None: nz = w3d.nz
    if ix_axis is None: ix_axis = nint(-w3d.xmmin/w3d.dx)
    if iy_axis is None: iy_axis = nint(-w3d.ymmin/w3d.dy)
    if xmesh is None: xmesh = w3d.xmesh
    if ymesh is None: ymesh = w3d.ymesh
    if l2symtry is None: l2symtry = w3d.l2symtry
    if l4symtry is None: l4symtry = w3d.l4symtry
    if srfrv_pernz > 0:
        save_srfrv_pernz = f3d.srfrv_pernz
        f3d.srfrv_pernz = srfrv_pernz

    # --- Make sure the rofzfunc is in main.
    # --- Note that this can only really work if a reference to the function
    # --- is passed in (instead of the name).
    if not f3d.lsrlinr and callable(rminofz):
        __main__.__dict__[rminofz.__name__] = rminofz
    if not f3d.lsrlinr and callable(rmaxofz):
        __main__.__dict__[rmaxofz.__name__] = rmaxofz

    # --- Get the name of the input function if a reference to the function
    # --- was passed in.
    if callable(rminofz): rminofz = rminofz.__name__
    if callable(rmaxofz): rmaxofz = rmaxofz.__name__

    # --- Get srfrv function - in older versions of the code the final 'f'
    # --- was not there
    try:
        srfrvfunc = f3d.srfrvinoutf
    except AttributeError:
        srfrvfunc = f3d.srfrvinout

    # --- Now call the fortran version
    srfrvfunc(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,lzend,
              xmin,xmax,ymin,ymax,lshell,zmminlocal,zmmaxlocal,zmmin,zbeam,dx,dy,dz,
              nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)

    # --- Reset srfrv_pernz if needed
    if srfrv_pernz > 0: f3d.srfrv_pernz = save_srfrv_pernz

#---------------------------------------------------------------------------
def  platepnt(ixmin=None,ixmax=None,iymin=None,iymax=None,
              ix_axis=None,iy_axis=None,dx=None,dy=None,
              aper=None,rmax=None,vvv=None,xoff=None,yoff=None,
              delz_in=None,iz=None,lz_in_plate=None,fuzz=None,condid=0):
    """
  Python interface for the platepnt routine. This now just calls
  the srfrvout routine. Note that the option lz_in_plate is now ignored.
    ixmin=0: minimum value of ix
    ixmax=w3d.nx: maximum value of ix
    iymin=0: minimum value of iy
    iymax=w3d.ny: maximum value of iy
    ix_axis=nint(-xmmin/dx): x grid location of beam center
    iy_axis=nint(-ymmin/dy): y grid location of beam center
    dx=w3d.dx: grid cell size in x
    dy=w3d.dy: grid cell size in y
    aper=0.: inner aperture of the plate
    rmax=LargePos: maximum radius of the plate
    vvv=0.: voltage on the plate
    xoff=0.: x offset of aperture
    yoff=0.: y offset of aperture
    delz_in=0.: fraction of cell, when outside of plate, to edge of plate
    iz=0: axial grid location of the plate
    fuzz=1.e-5*w3d.dz: number smalled compared to grid cell size, used
                       to prevent precision problems
    condid=0: Id number to identify this conductor
    """
    if ixmin is None: ixmin = 0
    if ixmax is None: ixmax = w3d.nx
    if iymin is None: iymin = 0
    if iymax is None: iymax = w3d.ny
    if ix_axis is None: ix_axis = nint(-w3d.xmmin/w3d.dx)
    if iy_axis is None: iy_axis = nint(-w3d.ymmin/w3d.dy)
    if dx is None: dx = w3d.dx
    if dy is None: dy = w3d.dy
    if aper is None: aper = 0.
    if rmax is None: rmax = top.largepos
    if vvv is None: vvv = 0.
    if xoff is None: xoff = 0.
    if yoff is None: yoff = 0.
    if delz_in is None: delz_in = 0.
    if iz is None: iz = 0
    if lz_in_plate is None: lz_in_plate = false
    if fuzz is None: fuzz = 1.e-5*w3d.dz

    f3d.lsrlinr = true
    f3d.npnts_sr = 2
    gchange("Surface_of_Rev")
    f3d.z_sr[0] = w3d.zmminlocal + iz*w3d.dz - 1.e-11*w3d.dz
    f3d.z_sr[1] = w3d.zmminlocal + iz*w3d.dz + 1.e-11*w3d.dz
    f3d.r_sr[:] = aper
    srfrvout(" ",vvv,f3d.z_sr[0],f3d.z_sr[1],xoff,yoff,rmax,true,
             w3d.xmmin+ixmin*w3d.dx,w3d.xmmin+ixmax*w3d.dx,
             w3d.ymmin+iymin*w3d.dy,w3d.ymmin+iymax*w3d.dy,true,
             w3d.zmminlocal,w3d.zmmaxlocal,top.zbeam,dx,dy,w3d.dz,
             w3d.nx,w3d.ny,w3d.nz,ix_axis,iy_axis,condid=condid)

    f3d.lsrlinr = false

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
def setconductorvoltage(voltage,condid=0,discrete=false,setvinject=false,
                        conductors=None):
    """
  Sets the voltage on a conductor, given an id.
   - voltage: voltage on conductor. Can be one of the following...
       - list of voltages at z grid cell locations, the list must have the same
         number of points along z as the grid, i.e. w3d.nz.
       - function which takes three arguments, x, y, z positions relative to lab
         frame and returns the voltage (only works in 3d)
       - a scalar voltage value
   - condid=0: conductor object or id number
   - discrete=false: when true, z locations for plus/minus z subgrid
                     points are round up/down.
   - setvinject=false: when true, sets top.vinject
   - conductors=f3d.conductors: allows alternate conductor to be set other
                                than the default ones
    """
    if isinstance(condid,Assembly): condid = condid.condid
    solver = getregisteredsolver()
    if conductors is None:
        if solver is not None:
            try:
                solver.setconductorvoltage(voltage,condid,discrete,setvinject)
            except AttributeError:
                print "Warning: setconductorvoltage not implemented for the resgisterd solver"
            return
        else:
            conductors = f3d.conductors

    interior = conductors.interior
    evensubgrid = conductors.evensubgrid
    oddsubgrid = conductors.oddsubgrid
    n = interior.n
    ne = evensubgrid.n
    no = oddsubgrid.n

    # --- Set vinject first is requested.
    if setvinject:
        if type(voltage) in [list,tuple,ndarray]:
            # --- Set it to the voltage on the left edge
            top.vinject = voltage[0]
        elif callable(voltage):
            # --- Set it to the voltage at the source center
            top.vinject = voltage(top.xinject,top.yinject,top.zinject)
        else:
            top.vinject = voltage

    if w3d.solvergeom in [w3d.RZgeom, w3d.XZgeom, w3d.XYgeom] and solver is None:
        if type(voltage) in [list,tuple,ndarray]:
        # --- Voltage is assumed to be the voltages are the z grid cell locations
        # --- (in the global beam frame).
            setconductorvoltagerz(voltage,w3d.nz,w3d.zmmin,w3d.dz,discrete,
                                  condid)
        else:
            setconductorvoltagerz_id(condid,voltage)
        return

    if type(voltage) in [list,tuple,ndarray]:
        # --- Voltage is assumed to be the voltages are the z grid cell locations
        # --- (in the global beam frame).

        if interior.n > 0:
            # --- Get z location of the conductor points (taking into account
            # --- the differing grid cell sizes for coarse levels). Use that to
            # --- gather the value of voltage at those locations.
            lz = take(conductors.levellz,interior.ilevel)
            iz = take(conductors.leveliz,interior.ilevel)
            icz = nint(interior.indx[2,:]*lz + iz)
            icz = minimum(icz,len(voltage)-2,icz)
            cv = take(voltage,icz)

        if evensubgrid.n > 0:
            # --- Get z location of even points. For conductors in the transverse
            # --- plane, use the voltage at the grid cell locations.
            iecl = take(conductors.levellz,evensubgrid.ilevel)
            iecliz = take(conductors.leveliz,evensubgrid.ilevel)
            iecz = nint(evensubgrid.indx[2,:]*iecl + iecliz)
            iecz = minimum(iecz,len(voltage)-2,iecz)
            ecv = take(voltage,iecz)
            ecvmx = ecv
            ecvpx = ecv
            ecvmy = ecv
            ecvpy = ecv

            # --- For conductors to the left, find the location of the conductor
            # --- and linear interpolate from the voltage data. If the discrete flag
            # --- is set, then round down to the nearest grid point.
            ecmz = iecz + where(logical_and(0 < evensubgrid.dels[4,:],evensubgrid.dels[4,:] < 1.),
                                -evensubgrid.dels[4,:],0)*iecl
            iecmz = ecmz.astype(long)
            if discrete: wecmz = 0.
            else:        wecmz = ecmz - iecmz
            ecvmz = take(voltage,iecmz)*(1.-wecmz) + take(voltage,iecmz+1)*wecmz

            # --- Same for conductors to the right. If discrete is set, round up.
            ecpz = iecz + where(logical_and(0 < evensubgrid.dels[5,:],evensubgrid.dels[5,:] < 1.),
                                -evensubgrid.dels[5,:],0)*iecl
            iecpz = ecpz.astype(long)
            if discrete: wecpz = 1.
            else:        wecpz = ecpz - iecpz
            ecvpz = take(voltage,iecpz)*(1.-wecpz) + take(voltage,iecpz+1)*wecpz

        if oddsubgrid.n > 0:
            # --- Repeat for odd conductor points.
            iocl = take(conductors.levellz,oddsubgrid.ilevel)
            iocliz = take(conductors.leveliz,oddsubgrid.ilevel)
            iocz = nint(oddsubgrid.indx[2,:]*iocl + iocliz)
            iocz = minimum(iocz,len(voltage)-2,iocz)
            ocv = take(voltage,iocz)
            ocvmx = ocv
            ocvpx = ocv
            ocvmy = ocv
            ocvpy = ocv

            ocmz = iocz + where(logical_and(0 < oddsubgrid.dels[4,:],oddsubgrid.dels[4,:] < 1.),
                                -oddsubgrid.dels[4,:],0)*iocl
            iocmz = ocmz.astype(long)
            if discrete: wocmz = 0.
            else:        wocmz = ocmz - iocmz
            ocvmz = take(voltage,iocmz)*(1.-wocmz) + take(voltage,iocmz+1)*wocmz

            ocpz = iocz + where(logical_and(0 < oddsubgrid.dels[5,:],oddsubgrid.dels[5,:] < 1.),
                                -oddsubgrid.dels[5,:],0)*iocl
            iocpz = ocpz.astype(long)
            if discrete: wocpz = 1.
            else:        wocpz = ocpz - iocpz
            ocvpz = take(voltage,iocpz)*(1.-wocpz) + take(voltage,iocpz+1)*wocpz

    elif callable(voltage):
        # --- Assumes that voltage is a function which takes 3 arguments, the
        # --- coordinates x, y, z, in meters relative to the beam frame.
        solver = getregisteredsolver()
        if solver is None: solver = w3d

        if interior.n > 0:
            icx = (interior.indx[0,:]*take(nint(conductors.levellx),interior.ilevel)+
                             take(conductors.levelix,interior.ilevel))
            icy = (interior.indx[1,:]*take(nint(conductors.levelly),interior.ilevel)+
                             take(conductors.leveliy,interior.ilevel))
            icz = (interior.indx[2,:]*take(nint(conductors.levellz),interior.ilevel)+
                             take(conductors.leveliz,interior.ilevel))
            cx = solver.xmmin + icx*solver.dx
            cy = solver.ymmin + icy*solver.dy
            cz = solver.zmmin + icz*solver.dz
            cv = voltage(cx,cy,cz)

        def between0and1(x):
            return logical_and(0.<x,x<1.)

        if evensubgrid.n > 0:
            ieclx = take(nint(conductors.levellx),evensubgrid.ilevel)
            iecly = take(nint(conductors.levelly),evensubgrid.ilevel)
            ieclz = take(nint(conductors.levellz),evensubgrid.ilevel)
            ecx = solver.xmmin + solver.dx*(evensubgrid.indx[0,:]*ieclx +
                                      take(conductors.levelix,evensubgrid.ilevel))
            ecy = solver.ymmin + solver.dy*(evensubgrid.indx[1,:]*iecly +
                                      take(conductors.leveliy,evensubgrid.ilevel))
            ecz = solver.zmmin + solver.dz*(evensubgrid.indx[2,:]*ieclz +
                                      take(conductors.leveliz,evensubgrid.ilevel))

            edels = evensubgrid.dels
            ecxmx = ecx - where(between0and1(edels[0,:]),edels[0,:],0)*ieclx*solver.dx
            ecxpx = ecx + where(between0and1(edels[1,:]),edels[1,:],0)*ieclx*solver.dx
            ecymy = ecy - where(between0and1(edels[2,:]),edels[2,:],0)*iecly*solver.dy
            ecypy = ecy + where(between0and1(edels[3,:]),edels[3,:],0)*iecly*solver.dy
            eczmz = ecz - where(between0and1(edels[4,:]),edels[4,:],0)*ieclz*solver.dz
            eczpz = ecz + where(between0and1(edels[5,:]),edels[5,:],0)*ieclz*solver.dz
            ecvmx = voltage(ecxmx,ecy  ,ecz  )
            ecvpx = voltage(ecxpx,ecy  ,ecz  )
            ecvmy = voltage(ecx  ,ecymy,ecz  )
            ecvpy = voltage(ecx  ,ecypy,ecz  )
            ecvmz = voltage(ecx  ,ecy  ,eczmz)
            ecvpz = voltage(ecx  ,ecy  ,eczpz)

        if oddsubgrid.n > 0:
            ioclx = take(nint(conductors.levellx),oddsubgrid.ilevel)
            iocly = take(nint(conductors.levelly),oddsubgrid.ilevel)
            ioclz = take(nint(conductors.levellz),oddsubgrid.ilevel)
            ocx = solver.xmmin + solver.dx*(oddsubgrid.indx[0,:]*ioclx +
                                      take(conductors.levelix,oddsubgrid.ilevel))
            ocy = solver.ymmin + solver.dy*(oddsubgrid.indx[1,:]*iocly +
                                      take(conductors.leveliy,oddsubgrid.ilevel))
            ocz = solver.zmmin + solver.dz*(oddsubgrid.indx[2,:]*ioclz +
                                      take(conductors.leveliz,oddsubgrid.ilevel))

            odels = oddsubgrid.dels
            ocxmx = ocx - where(between0and1(odels[0,:]),odels[0,:],0)*ioclx*solver.dx
            ocxpx = ocx + where(between0and1(odels[1,:]),odels[1,:],0)*ioclx*solver.dx
            ocymy = ocy - where(between0and1(odels[2,:]),odels[2,:],0)*iocly*solver.dy
            ocypy = ocy + where(between0and1(odels[3,:]),odels[3,:],0)*iocly*solver.dy
            oczmz = ocz - where(between0and1(odels[4,:]),odels[4,:],0)*ioclz*solver.dz
            oczpz = ocz + where(between0and1(odels[5,:]),odels[5,:],0)*ioclz*solver.dz
            ocvmx = voltage(ocxmx,ocy  ,ocz  )
            ocvpx = voltage(ocxpx,ocy  ,ocz  )
            ocvmy = voltage(ocx  ,ocymy,ocz  )
            ocvpy = voltage(ocx  ,ocypy,ocz  )
            ocvmz = voltage(ocx  ,ocy  ,oczmz)
            ocvpz = voltage(ocx  ,ocy  ,oczpz)

    else:
        cv = voltage
        ecvmx = ecvpx = ecvmy = ecvpy = ecvmz = ecvpz = voltage
        ocvmx = ocvpx = ocvmy = ocvpy = ocvmz = ocvpz = voltage

    # --- Now, put the voltage data into the fortran arrays.
    if conductors.interior.n > 0:
        f = conductors.interior
        f.volt[:] = where(equal(f.numb,condid),cv,f.volt)
    if conductors.evensubgrid.n > 0:
        f = conductors.evensubgrid
        f.volt[0,:] = where(equal(f.numb[0,:],condid),ecvmx,f.volt[0,:])
        f.volt[1,:] = where(equal(f.numb[1,:],condid),ecvpx,f.volt[1,:])
        f.volt[2,:] = where(equal(f.numb[2,:],condid),ecvmy,f.volt[2,:])
        f.volt[3,:] = where(equal(f.numb[3,:],condid),ecvpy,f.volt[3,:])
        f.volt[4,:] = where(equal(f.numb[4,:],condid),ecvmz,f.volt[4,:])
        f.volt[5,:] = where(equal(f.numb[5,:],condid),ecvpz,f.volt[5,:])
    if conductors.oddsubgrid.n > 0:
        f = conductors.oddsubgrid
        f.volt[0,:] = where(equal(f.numb[0,:],condid),ocvmx,f.volt[0,:])
        f.volt[1,:] = where(equal(f.numb[1,:],condid),ocvpx,f.volt[1,:])
        f.volt[2,:] = where(equal(f.numb[2,:],condid),ocvmy,f.volt[2,:])
        f.volt[3,:] = where(equal(f.numb[3,:],condid),ocvpy,f.volt[3,:])
        f.volt[4,:] = where(equal(f.numb[4,:],condid),ocvmz,f.volt[4,:])
        f.volt[5,:] = where(equal(f.numb[5,:],condid),ocvpz,f.volt[5,:])

#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
def visualizeconductors(condid=None,color=None,mglevel=0,
                        scene=None,title="Conductor geometry",vrange=None,
                        smooth=0,showgrid=0,showaxes=0,
                        conductors=f3d.conductors):
    """
  Creates 3-D visualization of the conductors based off of the subgrid data.
   - condid=None: optional conductor ID of object to draw
   - color=None: color, in the form of a list, [r,g,b], e.g. [1,0,0] to get red
   - mglevel=0: multigrid level to draw
   - scene=None: scene to use - by default creates a new one.
              the scene is returned by the function
   - title="Conductor geometry": window title
   - vrange=None: range of each dimension - used to scale size of image, in form
               [x,y,z]. e.g. to decrease z by 10, use [1,1,10]
   - smooth=0: not yet supported
   - showgrid=0: When true, displays field grid box
   - showaxes=0: When true, displays the axes
   - conductors=f3d.conductors: allows alternate conductor to be visualized other
                                than the default ones

  Returns the scene use to draw the image
    """
    if conductors is None: return
    if lparallel: return
    try:
        import Opyndx
    except ImportError:
        return

    solver = getregisteredsolver()
    if solver is None:
        solver = w3d
        solvertop = top
    else:
        solvertop = solver
    # --- Make sure that the conductor data is properly installed.
    checkconductors(solver.nx,solver.ny,solver.nz,
                    solver.nxlocal,solver.nylocal,solver.nzlocal,
                    solver.dx,solver.dy,solver.dz,conductors,solvertop.fsdecomp)

    # --- Save grid size
    nx = solver.nxlocal
    ny = solver.nylocal
    nz = solver.nzlocal

    # --- Get conductors
    ie = conductors.evensubgrid.istart[mglevel  ] - 1
    io = conductors.oddsubgrid.istart[mglevel  ] - 1
    ne = conductors.evensubgrid.istart[mglevel+1] - 1
    no = conductors.oddsubgrid.istart[mglevel+1] - 1
    if condid is None:
        # --- Get all conductors
        f = conductors.evensubgrid
        iecndx = f.indx[0,ie:ne]
        iecndy = f.indx[1,ie:ne]
        iecndz = f.indx[2,ie:ne]
        ecdelmx = abs(f.dels[0,ie:ne])
        ecdelpx = abs(f.dels[1,ie:ne])
        ecdelmy = abs(f.dels[2,ie:ne])
        ecdelpy = abs(f.dels[3,ie:ne])
        ecdelmz = abs(f.dels[4,ie:ne])
        ecdelpz = abs(f.dels[5,ie:ne])
        f = conductors.oddsubgrid
        iocndx = f.indx[0,io:no]
        iocndy = f.indx[1,io:no]
        iocndz = f.indx[2,io:no]
        ocdelmx = abs(f.dels[0,io:no])
        ocdelpx = abs(f.dels[1,io:no])
        ocdelmy = abs(f.dels[2,io:no])
        ocdelpy = abs(f.dels[3,io:no])
        ocdelmz = abs(f.dels[4,io:no])
        ocdelpz = abs(f.dels[5,io:no])
    else:
        # --- Get only points matching the specified condid
        f = conductors.evensubgrid
        iecndx = compress(f.numb[0,ie:ne]==condid,f.indx[0,ie:ne])
        iecndy = compress(f.numb[0,ie:ne]==condid,f.indx[1,ie:ne])
        iecndz = compress(f.numb[0,ie:ne]==condid,f.indx[2,ie:ne])
        ecdelmx = compress(f.numb[0,ie:ne]==condid,f.abs(dels[0,ie:ne]))
        ecdelpx = compress(f.numb[1,ie:ne]==condid,f.abs(dels[1,ie:ne]))
        ecdelmy = compress(f.numb[2,ie:ne]==condid,f.abs(dels[2,ie:ne]))
        ecdelpy = compress(f.numb[3,ie:ne]==condid,f.abs(dels[3,ie:ne]))
        ecdelmz = compress(f.numb[4,ie:ne]==condid,f.abs(dels[4,ie:ne]))
        ecdelpz = compress(f.numb[5,ie:ne]==condid,f.abs(dels[5,ie:ne]))
        f = conductors.oddsubgrid
        iocndx = compress(f.numb[0,io:no]==condid,f.indx[0,io:no])
        iocndy = compress(f.numb[0,io:no]==condid,f.indx[1,io:no])
        iocndz = compress(f.numb[0,io:no]==condid,f.indx[2,io:no])
        ocdelmx = compress(f.numb[0,io:no]==condid,f.abs(dels[0,io:no]))
        ocdelpx = compress(f.numb[1,io:no]==condid,f.abs(dels[1,io:no]))
        ocdelmy = compress(f.numb[2,io:no]==condid,f.abs(dels[2,io:no]))
        ocdelpy = compress(f.numb[3,io:no]==condid,f.abs(dels[3,io:no]))
        ocdelmz = compress(f.numb[4,io:no]==condid,f.abs(dels[4,io:no]))
        ocdelpz = compress(f.numb[5,io:no]==condid,f.abs(dels[5,io:no]))

    nn = len(iecndx) + len(iocndx)
    if nn == 0: return
    icndx = concatenate((iecndx,iocndx))
    icndy = concatenate((iecndy,iocndy))
    icndz = concatenate((iecndz,iocndz))
    delmx = concatenate((ecdelmx,ocdelmx))
    delpx = concatenate((ecdelpx,ocdelpx))
    delmy = concatenate((ecdelmy,ocdelmy))
    delpy = concatenate((ecdelpy,ocdelpy))
    delmz = concatenate((ecdelmz,ocdelmz))
    delpz = concatenate((ecdelpz,ocdelpz))
    icnd = array([icndx,icndy,icndz])
    dels = array([delmx,delpx,delmy,delpy,delmz,delpz])

    model = Opyndx.VisualModel(twoSided=true,scene=scene,title=title,
                               vrange=vrange)

    gridmin = array([solver.xmmin,solver.ymmin,solver.zmmin])
    griddd = array([solver.dx*conductors.levellx[mglevel],
                    solver.dy*conductors.levelly[mglevel],
                    solver.dz*conductors.levellz[mglevel]])
    gridnn = array([nint(solver.nx/conductors.levellx[mglevel]),
                    nint(solver.ny/conductors.levelly[mglevel]),
                    nint(solver.nz/conductors.levellz[mglevel])])

    # --- This fortran routine generates the triangulated surface. It was
    # --- converted to fortran for speed.
    f3d.ntriangles = 0
    getconductorfacets(nn,icnd,dels,gridnn,griddd,gridmin)
    if smooth:
        f3d.maxtriangles = f3d.ntrianges
        gchange("ConductorGeometryVisualization")
        tt = f3d.triangles - gridmin[:,newaxis,newaxis]
        tt = tt[0,:,:]**2 + tt[1,:,:]**2 + tt[2,:,:]**2
        tt = (tt/max(tt)*100000000).astype(long)
        tt.shape = (3*f3d.ntriangles,)
        ii = argsort(tt)
        conductorsmoothshading(tt,ii)

    model.triangles = reshape(transpose(f3d.triangles),(3*f3d.ntriangles,3))
    model.normals = reshape(transpose(f3d.normals),(3*f3d.ntriangles,3))
    if color is not None:
        model.colors = (3*f3d.ntriangles)*[color]

    model.Display(showgrid=showgrid,showaxes=showaxes)
    return model
