from __future__ import generators
from warp import *
from MeshRefinement import *
try:
    import Opyndx
    VisualizableClass = Opyndx.Visualizable
except ImportError:
    # --- If Opyndx is not available, then use object as the base class,
    # --- disabling any visualization.
    VisualizableClass = object
import time
import gc # Garbage collection

try:
    import psyco
except ImportError:
    pass


class AMRTree(VisualizableClass):
    """
  Adaptive Mesh Refinement class.
    """
    def __init__(self,lremoveblockswithoutparticles=0):
        self.lremoveblockswithoutparticles = lremoveblockswithoutparticles
        self.nblocks      = 0
        self.f            = None
        self.f_user       = None
        self.nbcells_user = None
        if getcurrpkg()=='w3d' and w3d.solvergeom==w3d.XYZgeom:
            self.solvergeom = w3d.XYZgeom
            self.blocks=MRBlock()
        else:
            self.solvergeom = w3d.solvergeom
            try:
                self.blocks={}
                self.blocks[frz.basegrid.gid[0]]={'grid':frz.basegrid,'installed_conductors':[]}
            except:
                self.blocks=MRBlock()
        self.colors       = ['red','blue','yellow','green','cyan','magenta','white']
        self.conductors   = []
        self.conductorsdfill   = []
        self.beforefs = ControllerFunction()
        self.dfill = 2
        self.enable()

    def enable(self):
        if self.solvergeom == w3d.XYZgeom:
            registersolver(self)
        else:
            try:
                bg=frz.basegrid
                installbeforefs(self.generate)
            except:
                registersolver(self)
    def __getstate__(self):
        if type(self.blocks) != DictType:
            # --- In this case, there is no trouble with pickle
            return self.__dict__
        else:
            # --- In this case, the complexly nested frz.basegrid structure
            # --- causes pickle to have fits and it often reachs the
            # --- recursion depth limit. So, remove all references to grid here,
            # --- replacing them with the gids so they can be restored later.
            # --- This also removes the problem of having extra instances of
            # --- the grids saved in the dump file.
            dict = self.__dict__.copy()
            cleanblocks = self.blocks.copy()
            del dict['blocks']
            dict['cleanblocks'] = cleanblocks
            for gid in cleanblocks:
                cleanblocks[gid]['grid'] = gid
            return dict

    def __setstate__(self,dict):
        # --- This is called when the instance is unpickled.
        # --- Note that the restorefrzgrid method below is not called now since
        # --- this instance may be restored before frz.basegrid is restored.
        self.__dict__.update(dict)
        self.enable()
        installafterrestart(self.restorefrzgrid)

    def restorefrzgrid(self):
        # --- This should be called whenever 'grid' is used to ensure that it
        # --- has been properly restored after being unpickled.
        # --- Note that it can't be called right away since the frz.basegrid
        # --- may not have been restored yet.
        try:
            cleanblocks = self.__dict__['cleanblocks']
            for g in walkfrzgrid():
                cleanblocks[g.gid[0]]['grid'] = g
            self.__dict__['blocks'] = cleanblocks
            del self.__dict__['cleanblocks']
        except KeyError:
            pass

    # --------------------------------------------------------------------
    # --- Field solver API routines.
    # --- Currently, these only support the 3d version.
    def loadrho(self,lzero=true,lfinalize_rho=true):
        # --- If the mesh refinement is going to be recalculated this step,
        # --- then rho is only needed on the root block. Skip loading rho
        # --- on the mesh refined blocks to save time.
        ifcond = w3d.AMRgenerate_periodicity==0
        if not ifcond:
            ifcond = top.it%w3d.AMRgenerate_periodicity != 0
        if (ifcond or not lzero): lrootonly = 0
        else:                     lrootonly = 1
        if 1:#self.solvergeom == w3d.XYZgeom:
            self.blocks.loadrho(lzero=lzero,lfinalize_rho=lfinalize_rho,lrootonly=lrootonly)
#        else:
#            raise Exception("Only 3d supported as registered solver.")

    def loadj(self,lzero=true,lfinalize_rho=true):
        # --- If the mesh refinement is going to be recalculated this step,
        # --- then rho is only needed on the root block. Skip loading rho
        # --- on the mesh refined blocks to save time.
        ifcond = w3d.AMRgenerate_periodicity==0
        if not ifcond:
            ifcond = top.it%w3d.AMRgenerate_periodicity != 0
        if (ifcond or not lzero): lrootonly = 0
        else:                     lrootonly = 1
        if 1:#self.solvergeom == w3d.XYZgeom:
            self.blocks.loadj(lzero=lzero,lfinalize_rho=lfinalize_rho,lrootonly=lrootonly)

    # --- Diagnostic routines
    def rhodia(self):
        pass
    def gtlchg(self):
        pass
    def srhoax(self):
        pass
    def getese(self):
        pass
    def sphiax(self):
        pass
    def sezax(self):
        pass

    def solve(self):
#        if self.solvergeom == w3d.XYZgeom:
            self.generate()
            self.blocks.solve()
#        else:
#            raise Exception("Only 3d supported as registered solver.")

    def fetche(self):
#        if self.solvergeom == w3d.XYZgeom:
            self.blocks.fetche()
#        else:
#            raise Exception("Only 3d supported as registered solver.")

    def fetchb(self):
#        if self.solvergeom == w3d.XYZgeom:
            self.blocks.fetchb()
#        else:
#            raise Exception("Only 3d supported as registered solver.")

    def fetchphi(self):
#        if self.solvergeom == w3d.XYZgeom:
            self.blocks.fetchphi()
#        else:
#            raise Exception("Only 3d supported as registered solver.")

    def installconductor(self,conductor,dfill=2):
        self.addconductor(conductor,dfill)
        if self.solvergeom==w3d.XYZgeom:
            self.blocks.installconductor(conductor,dfill=self.dfill)
        else:
            try:
                bg=frz.basegrid
                self.restorefrzgrid()
                for block in self.blocks.itervalues():
                    g = block['grid']
                    if conductor not in block['installed_conductors']:
                        try:
                            cond = conductor.cond
                        except AttributeError:
                            cond = conductor
                        installconductors(cond,nx=g.nr,ny=0,nzlocal=g.nz,nz=g.nz,
                                            xmmin=g.xmin,xmmax=g.xmax,
                                            zmmin=g.zmin,zmmax=g.zmax,dfill=dfill,
                                            gridrz=g)
                        block['installed_conductors'].append(cond)
                get_cond_rz(1)
            except:
                self.blocks.installconductor(conductor,dfill=self.dfill)
            
    def hasconductors(self):
        if self.solvergeom == w3d.XYZgeom:
            return self.blocks.hasconductors()
        else:
            try:
                bg=frz.basegrid
                return f3d.conductors.interior.n > 0
            except:
                return self.blocks.hasconductors()
            
    def getconductors(self):
        if self.solvergeom == w3d.XYZgeom:
            return self.blocks.getconductors(alllevels=1)
        else:
            try:
                bg=frz.basegrid
                return f3d.conductors
            except:
                return self.blocks.getconductors(alllevels=1)

    def setconductorvoltage(self,voltage,condid=0,discrete=false,setvinject=false):
        if self.solvergeom == w3d.XYZgeom:
            self.blocks.setconductorvoltage(voltage,condid,discrete,setvinject)
        else:
            try:
                bg=frz.basegrid
                setconductorvoltage(voltage,condid,discrete,setvinject)
            except:
                self.blocks.setconductorvoltage(voltage,condid,discrete,setvinject)

    def installbeforefs(self,f):
        self.beforefs.installfuncinlist(f)

    def uninstallbeforefs(self,f):
        self.beforefs.uninstallfuncinlist(f)

    def pfxy(self,**kw): self.blocks.pfxy(kwdict=kw)
    def pfzx(self,**kw): self.blocks.pfzx(kwdict=kw)
    def pfzy(self,**kw): self.blocks.pfzy(kwdict=kw)
    def pfxyg(self,**kw): self.blocks.pfxyg(kwdict=kw)
    def pfzxg(self,**kw): self.blocks.pfzxg(kwdict=kw)
    def pfzyg(self,**kw): self.blocks.pfzyg(kwdict=kw)

    # --------------------------------------------------------------------


    def addconductor(self,conductor,dfill=2):
        self.conductors.append(conductor)
        self.conductorsdfill.append(dfill)

    def getabsgrad(self,f,dx,dy,dz):
        """
      get average of absolute value of grad(f).
      f is supposed to be in Fortran ordering.
        """
        if rank(f)==2:
            s=shape(f)
            nx = s[0]
            ny = s[1]
            g = fzeros([nx+2,ny+2],float64)
            # fill interior
            g[1:-1,1:-1] = f
            # set boundaries
            g[0,   1:-1] = 2.*f[0, :]-f[1, :]
            g[-1,  1:-1] = 2.*f[-1,:]-f[-2,:]
            g[1:-1,0   ] = 2.*f[:, 0]-f[:, 1]
            g[1:-1,-1  ] = 2.*f[:,-1]-f[:,-2]
            # computes average of gradients
            gr = 0.5*abs((g[1:-1,1:-1]-g[ :-2,1:-1])/dy) \
               + 0.5*abs((g[2:,  1:-1]-g[1:-1,1:-1])/dy) \
               + 0.5*abs((g[1:-1,1:-1]-g[1:-1, :-2])/dx) \
               + 0.5*abs((g[1:-1,2:  ]-g[1:-1,1:-1])/dx)
        else:
            #s=shape(f)
            #g = fzeros(array(s)+2,float64)
            ## fill interior
            #g[1:-1,1:-1,1:-1] = f
            ## set boundaries
            #g[0   ,1:-1,1:-1] = f[1 ,: ,: ]
            #g[-1  ,1:-1,1:-1] = f[-2,: ,: ]
            #g[1:-1,0   ,1:-1] = f[: ,1 ,: ]
            #g[1:-1,-1  ,1:-1] = f[: ,-2,: ]
            #g[1:-1,1:-1,0   ] = f[: ,: ,1 ]
            #g[1:-1,1:-1,-1  ] = f[: ,: ,-2]
            ## computes average of gradients
            ## --- The factor of 0.5 does affect anything since the gradient
            ## --- is scaled by its maximum.
            #gr = abs((g[1:-1,1:-1,1:-1] - g[ :-2,1:-1,1:-1])/dx) \
            #   + abs((g[2:,  1:-1,1:-1] - g[1:-1,1:-1,1:-1])/dx) \
            #   + abs((g[1:-1,1:-1,1:-1] - g[1:-1, :-2,1:-1])/dy) \
            #   + abs((g[1:-1,2:  ,1:-1] - g[1:-1,1:-1,1:-1])/dy) \
            #   + abs((g[1:-1,1:-1,1:-1] - g[1:-1,1:-1, :-2])/dz) \
            #   + abs((g[1:-1,1:-1,2:  ] - g[1:-1,1:-1,1:-1])/dz)
            #gx = abs((g[1:-1,1:-1,1:-1] - g[ :-2,1:-1,1:-1]))
            #add(gx,abs((g[2:,  1:-1,1:-1] - g[1:-1,1:-1,1:-1])),gx)
            #divide(gx,dx,gx)
            #gy = abs((g[1:-1,1:-1,1:-1] - g[1:-1, :-2,1:-1]))
            #add(gy,abs((g[1:-1,2:  ,1:-1] - g[1:-1,1:-1,1:-1])),gy)
            #divide(gy,dy,gy)
            #gz = abs((g[1:-1,1:-1,1:-1] - g[1:-1,1:-1, :-2]))
            #add(gz,abs((g[1:-1,1:-1,2:  ] - g[1:-1,1:-1,1:-1])),gz)
            #divide(gz,dz,gz)
            #add(gx,gy,gx)
            #add(gx,gz,gx)
            gr = zeros(shape(f),'d')
            nx,ny,nz = array(shape(f)) - 1
            getabsgrad3d(nx,ny,nz,f,gr,dx,dy,dz)
        return gr

    def getedges_byslice(self,f,dx,dy,dz,threshold):
        """
      Returns array with non-zero values only at edges.
      For each line (horizontals and verticals), the values which
      are above threshold*max(values of f in line) are selected as edges.
        """
        # get average of absolute value of gradient of f
        fg = self.getabsgrad(f,dx,dy,dz)

        dim = rank(f)

        # get edges using vertical lines
        maxfg = fg.max(0)
        g1 = where(fg>threshold*maxfg[newaxis,...],fg,0.)

        # get edges using horizontal lines
        maxfg = fg.max(1)
        if dim==2:
            g2 = where(fg>threshold*maxfg[:,newaxis],fg,0.)
        elif dim==3:
            g2 = where(fg>threshold*maxfg[:,newaxis,:],fg,0.)

        # take max of g1 and g2
        g = where(g1>g2,g1,g2)

        if dim==3:
            maxfg = fg.max(2)
            g3 = where(fg>threshold*maxfg[...,newaxis],fg,0.)
            # returns max of g and g3
            g = where(g>g3,g,g3)

        return g

    def getnbcell_edges(self,f,dx,dy,dz,threshold,Rgrad):
        """
      Returns array with non-zero value RMR at edges, zero elsewhere.
      For each line (horizontals and verticals), the values which
      are above threshold*max(values of f in line) are selected as edges.
        """
        if Rgrad<=1: return 0*f
        fg = self.getedges_byslice(f,dx,dy,dz,threshold)
        m = maxnd(fg)
        return where(fg>1.e-10*m,aint(Rgrad),0)

    def getnbcell_rho(self,f,Rdens,MRfact):
        """
      returns nb cells proportional to density f, with mesh refinement factor of b
      (default=2) and maximum number of refined cells per coarse cell nmax (along
      one dimension).
        """
        if Rdens<=1: return 0*f
        # get dimension (2-D or 3-D)
        dim = rank(f)
        # get number of refinement levels
        n = nint(log(Rdens)/log(MRfact))
        # get nb cells proportional to f
#      fg=MRfact**(dim*(n+1))*f/maxnd(f)
        if 1:
            fg=MRfact**(dim*n)*f/maxnd(f)
            fg=where(fg>1,fg,1)
#        return MRfact**aint(log(fg)/log(dim**MRfact))
            return MRfact**aint(log(fg)/log(MRfact**dim))
#        return b**aint(log(fg)/log(b**dim))
        else:
            nbpcell=4
            if dim==2:
                fg = aint(f*w3d.dx*w3d.dy/(top.pgroup.sw[0]*echarge))
            else:
                fg = aint(f*w3d.dx*w3d.dy*w3d.dz/(top.pgroup.sw[0]*echarge))
            fg = where(fg>1,fg,1)
            fg = MRfact**aint(log(fg)/log(nbpcell*MRfact**dim))
            return where(fg>n,n,fg)

    def getnbcells(self,f,dx,dy,dz,Rdens,threshold,Rgrad,MRfact=2,lmax=4):
        if f is None: return None
        if maxnd(abs(self.f))==0: return None
        fg1 = self.getnbcell_edges(f,dx,dy,dz,threshold,Rgrad)
        fg2 = self.getnbcell_rho(f,Rdens,MRfact)
        f = aint(where(fg1>fg2,fg1,fg2))
        # next loop removes isolated blocks of lmax cells or less
        # this needs improvements
        if lmax>0:
            if rank(f)==3:
                nx=shape(f)[0]
                ny=shape(f)[1]
                nz=shape(f)[2]
                t = fzeros([nx,ny,nz],'l')
                r = maxnd(f)
                while r>=1:
                    sum_neighbors3d(where(f==r,0,1),t,nx-1,ny-1,nz-1)
                    fl = where(f==r and t<=lmax,1,0)
                    f = where(fl==1,r,f)
                    r=r/MRfact
            if rank(f)==2:
                nr=shape(f)[0]
                nz=shape(f)[1]
                t = fzeros([nr,nz],'l')
                r = maxnd(f)
                while r>=1:
                    sum_neighbors(where(f==r,0,1),t,nr-1,nz-1)
                    fl = where(f==r and t<=lmax,1,0)
                    f = where(fl==1,r,f)
                    r=r/MRfact
        return f

    def getconds2d(self,f0,fno,ib,lx,ux,ly,uy,progressive,nooverlap):
        f0t = f0[lx:ux,ly:uy]
        if(progressive):
            cond = sum(f0t>=ib)>0
        else:
            cond = sum(f0t==ib)>0
        if(nooverlap):
            cond2 = sum(fno[lx:ux,ly:uy])==0
        else:
            cond2 = true
        return cond,cond2

    def getconds3d(self,f0,fno,ib,lx,ux,ly,uy,lz,uz,progressive,nooverlap):
        f0t = f0[lx:ux,ly:uy,lz:uz]
        if(progressive):
            cond = sum(sum(f0t>=ib))>0
        else:
            cond = sum(sum(f0t==ib))>0
        if(nooverlap):
            cond2 = sum(sum(fno[lx:ux,ly:uy,lz:uz]))==0
        else:
            cond2 = true
        return cond,cond2

    def get_area_fraction(self,f,j,k,l,ix,iy,iz,ixm,iym,izm,ib,progressive=true):
        if rank(f)==2:
            if(progressive):
                return float(sum(sum(f[j-ixm:j+ix,k-iym:k+iy]>=ib)))/((ix+ixm)*(iy+iym))
            else:
                return float(sum(sum(f[j-ixm:j+ix,k-iym:k+iy]==ib)))/((ix+ixm)*(iy+iym))
        else:
            if(progressive):
                return float(sum(sum(sum(f[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz]>=ib))))/((ix+ixm)*(iy+iym)*(iz+izm))
            else:
                return float(sum(sum(sum(f[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz]==ib))))/((ix+ixm)*(iy+iym)*(iz+izm))

    def setlist(self,f,rl,b,progressive=true,nooverlap=true):
        self.listblocks = [0]
        if maxnd(abs(f)) == 0.: return
        p = progressive
        nlevels = nint(log(maxnd(f))/log(b))+1
        dim = rank(f)
        nx = shape(f)[0]
        ny = shape(f)[1]
        if dim==3:
            nz = shape(f)[2]
        else:
            nz=0; iz=0; izm=0; l=0

        try:
            self.f0
        except:
            self.f0=[]
        # loop all refinement levels
        for i in range(nlevels-1,0,-1):
            r = rl[min(i-1,len(rl)-1)]
            ib = b**i
            if(progressive and i<nlevels-1):
                f0 = where(self.sumpatch(listpatches,nx,ny,nz,dim)>0,ib,f)
            else:
                f0 = f.copy()
            f1 = ravel(f0)
            if(progressive):
                listnodes = oldnonzero(f1>=ib)
            else:
                listnodes = oldnonzero(f1==ib)
            listpatches = []
            if(nooverlap):
                if dim==2:
                    fno = zeros([nx,ny],'l')
                else:
                    fno = zeros([nx,ny,nz],'l')
            # loop all nodes where refinement is needed
            for n in listnodes:
                ix  = 1
                iy  = 1
                ixm = 1#0
                iym = 1#0
                if w3d.l2symtry:
                    iym=0
                if w3d.l4symtry:
                    ixm=0
                    iym=0
                if dim==2:
                    j = n/ny
                    k = n%ny
                    f0t = f0[j,k]
                else: #dim=3
                    j = n/(ny*nz)
                    k = n%(ny*nz)/nz
                    l = n%nz
                    iz  = 1
                    izm = 0
                    f0t = f0[j,k,l]
                # check if a patch is present
                if(progressive):
                    cond = f0t>=ib
                else:
                    cond = f0t==ib
                if(cond):
                    tryit = true
                    # try to expand the patch
                    while tryit:
                        ix0  = ix+0
                        iy0  = iy+0
                        ixm0 = ixm+0
                        iym0 = iym+0
                        if dim==2:
                        # x up
                            if (j+ix<nx):
                                cond,cond2=self.getconds2d(f0,fno,ib,j+ix,j+ix+1,k-iym,k+iy,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix+1,iy,iz,ixm,iym,izm,ib,p)>=r):ix+=1
                            # y up
                            if(k+iy<ny):
                                cond,cond2=self.getconds2d(f0,fno,ib,j-ixm,j+ix,k+iy,k+iy+1,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy+1,iz,ixm,iym,izm,ib,p)>=r):iy+=1
                            # x down
                            if(j-ixm-1>-1):
                                cond,cond2=self.getconds2d(f0,fno,ib,j-ixm-1,j-ixm,k-iym,k+iy,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm+1,iym,izm,ib,p)>=r):ixm+=1
                            # y down
                            if(k-iym-1>-1):
                                cond,cond2=self.getconds2d(f0,fno,ib,j-ixm,j+ix,k-iym-1,k-iym,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm,iym+1,izm,ib,p)>=r):iym+=1
                            if(ix==ix0 and iy==iy0 and ixm==ixm0 and iym==iym0):tryit=false
                        # end if dim==2

                        if dim==3:
                            iz0  = iz+0
                            izm0 = izm+0
                            # x up
                            if(j+ix<nx):
                                cond,cond2=self.getconds3d(f0,fno,ib,j+ix,j+ix+1,k-iym,k+iy,l-izm,l+iz,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix+1,iy,iz,ixm,iym,izm,ib,p)>=r):ix+=1
                            # y up
                            if(k+iy<ny):
                                cond,cond2=self.getconds3d(f0,fno,ib,j-ixm,j+ix,k+iy,k+iy+1,l-izm,l+iz,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy+1,iz,ixm,iym,izm,ib,p)>=r):iy+=1
                            # z up
                            if(l+iz<nz):
                                cond,cond2=self.getconds3d(f0,fno,ib,j-ixm,j+ix,k-iym,k+iy,l+iz,l+iz+1,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy,iz+1,ixm,iym,izm,ib,p)>=r):iz+=1
                            # x down
                            if(j-ixm-1>-1):
                                cond,cond2=self.getconds3d(f0,fno,ib,j-ixm-1,j-ixm,k-iym,k+iy,l-izm,l+iz,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm+1,iym,izm,ib,p)>=r):ixm+=1
                            # y down
                            if(k-iym-1>-1):
                                cond,cond2=self.getconds3d(f0,fno,ib,j-ixm,j+ix,k-iym-1,k-iym,l-izm,l+iz,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm,iym+1,izm,ib,p)>=r):iym+=1
                            # z down
                            if(l-izm-1>-1):
                                cond,cond2=self.getconds3d(f0,fno,ib,j-ixm,j+ix,k-iym,k+iy,l-izm-1,l-izm,progressive,nooverlap)
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm,iym,izm+1,ib,p)>=r):izm+=1

                            if(ix==ix0 and iy==iy0 and iz==iz0 and
                               ixm==ixm0 and iym==iym0 and izm==izm0):tryit=false

                    if dim==2:
                        f0[j-ixm:j+ix,k-iym:k+iy] = 0
                        if(nooverlap):fno[j-ixm:j+ix,k-iym:k+iy] = 1
                        listpatches.append([j-ixm,k-iym,ixm+ix,iym+iy])
                    else:
                        f0[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz] = 0
                        if(nooverlap):fno[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz] = 1
                        listpatches.append([j-ixm,k-iym,l-izm,ixm+ix,iym+iy,izm+iz])
            self.listblocks.insert(1,listpatches)

    def setlistold(self,f,rl,b,progressive=true,nooverlap=true):
        p = progressive
        nlevels = nint(log(maxnd(f))/log(b))+1
        dim = rank(f)
        nx = shape(f)[0]
        ny = shape(f)[1]
        if dim==3:
            nz = shape(f)[2]
        else:
            nz=0; iz=0; izm=0; l=0
        self.listblocks = [0]

        # loop all refinement levels
        for i in range(nlevels-1,0,-1):
            r = rl[min(i-1,len(rl)-1)]
            ib = b**i
            if(progressive and i<nlevels-1):
                f0 = where(self.sumpatch(listpatches,nx,ny,nz,dim)>0,ib,f)
            else:
                f0 = f.copy()
            f1 = ravel(f0)
            if(progressive):
                listnodes = oldnonzero(f1>=ib)
            else:
                listnodes = oldnonzero(f1==ib)
            listpatches = []
            if(nooverlap):
                if dim==2:
                    fno = zeros([nx,ny],'l')
                else:
                    fno = zeros([nx,ny,nz],'l')
            # loop all nodes where refinement is needed
            for n in listnodes:
                ix  = 1
                iy  = 1
                ixm = 0
                iym = 0
                if dim==2:
                    j = n/ny
                    k = n%ny
                    f0t = f0[j,k]
                else: #dim=3
                    j = n/(ny*nz)
                    k = n%(ny*nz)/nz
                    l = n%nz
                    iz  = 1
                    izm = 0
                    f0t = f0[j,k,l]
                # check if a patch is present
                if(progressive):
                    cond = f0t>=ib
                else:
                    cond = f0t==ib
                if(cond):
                    tryit = true
                    # try to expand the patch
                    while tryit:
                        ix0  = ix+0
                        iy0  = iy+0
                        ixm0 = ixm+0
                        iym0 = iym+0
                        if dim==3:
                            iz0  = iz+0
                            izm0 = izm+0
                        # x up
                        if(j+ix<nx):
                            if dim==2:
                                f0t = f0[j+ix,k-iym:k+iy]
                            else:
                                f0t = f0[j+ix,k-iym:k+iy,l-izm:l+iz]
                            if(progressive):
                                cond = sum(f0t>=ib)>0
                            else:
                                cond = sum(f0t==ib)>0
                            if(nooverlap):
                                if dim==2:
                                    cond2 = sum(fno[j+ix,k-iym:k+iy])==0
                                else:
                                    cond2 = sum(sum(fno[j+ix,k-iym:k+iy,l-izm:l+iz]))==0
                            else:
                                cond2 = true
                            if(cond and cond2):
                                if(self.get_area_fraction(f0,j,k,l,ix+1,iy,iz,ixm,iym,izm,ib,p)>=r):ix+=1
                        # y up
                        if(k+iy<ny):
                            if dim==2:
                                f0t = f0[j-ixm:j+ix,k+iy]
                            else:
                                f0t = f0[j-ixm:j+ix,k+iy,l-izm:l+iz]
                            if(progressive):
                                cond = sum(f0t>=ib)>0
                            else:
                                cond = sum(f0t==ib)>0
                            if(nooverlap):
                                if dim==2:
                                    cond2 = sum(fno[j-ixm:j+ix,k+iy])==0
                                else:
                                    cond2 = sum(sum(fno[j-ixm:j+ix,k+iy,l-izm:l+iz]))==0
                            else:
                                cond2 = true
                            if(cond and cond2):
                                if(self.get_area_fraction(f0,j,k,l,ix,iy+1,iz,ixm,iym,izm,ib,p)>=r):iy+=1
                        # z up
                        if dim==3:
                            if(l+iz<nz):
                                f0t = f0[j-ixm:j+ix,k-iym:k+iy,l+iz]
                                if(progressive):
                                    cond = sum(f0t>=ib)>0
                                else:
                                    cond = sum(f0t==ib)>0
                                if(nooverlap):
                                    cond2 = sum(sum(fno[j-ixm:j+ix,k-iym:k+iy,l+iz]))==0
                                else:
                                    cond2 = true
                                if(cond and cond2):
                                    if(self.get_area_fraction(f0,j,k,l,ix,iy,iz+1,ixm,iym,izm,ib,p)>=r):iz+=1
                        # x down
                        if(j-ixm-1>-1):
                            if dim==2:
                                f0t = f0[j-ixm-1,k-iym:k+iy]
                            else:
                                f0t = f0[j-ixm-1,k-iym:k+iy,l-izm:l+iz]
                            if(progressive):
                                cond = sum(f0t>=ib)>0
                            else:
                                cond = sum(f0t==ib)>0
                            if(nooverlap):
                                if dim==2:
                                    cond2 = sum(fno[j-ixm-1,k-iym:k+iy])==0
                                else:
                                    cond2 = sum(sum(fno[j-ixm-1,k-iym:k+iy,l-izm:l+iz]))==0
                            else:
                                cond2 = true
                            if(cond and cond2):
                                if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm+1,iym,izm,ib,p)>=r):ixm+=1
                        # y down
                        if(k-iym-1>-1):
                            if dim==2:
                                f0t = f0[j-ixm:j+ix,k-iym-1]
                            else:
                                f0t = f0[j-ixm:j+ix,k-iym-1,l-izm:l+iz]
                            if(progressive):
                                cond = sum(f0t>=ib)>0
                            else:
                                cond = sum(f0t==ib)>0
                            if(nooverlap):
                                if dim==2:
                                    cond2 = sum(fno[j-ixm:j+ix,k-iym-1])==0
                                else:
                                    cond2 = sum(sum(fno[j-ixm:j+ix,k-iym-1,l-izm:l+iz]))==0
                            else:
                                cond2 = true
                            if(cond and cond2):
                                if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm,iym+1,izm,ib,p)>=r):iym+=1
                        # z down
                        if(l-izm-1>-1):
                            f0t = f0[j-ixm:j+ix,k-iym:k+iy,l-izm-1]
                            if(progressive):
                                cond = sum(f0t>=ib)>0
                            else:
                                cond = sum(f0t==ib)>0
                            if(nooverlap):
                                cond2 = sum(sum(fno[j-ixm:j+ix,k-iym:k+iy,l-izm-1]))==0
                            else:
                                cond2 = true
                            if(cond and cond2):
                                if(self.get_area_fraction(f0,j,k,l,ix,iy,iz,ixm,iym,izm+1,ib,p)>=r):izm+=1
                        if dim==2:
                            if(ix==ix0 and iy==iy0 and ixm==ixm0 and iym==iym0):tryit=false
                        else:
                            if(ix==ix0 and iy==iy0 and iz==iz0 and
                               ixm==ixm0 and iym==iym0 and izm==izm0):tryit=false
                    if dim==2:
                        f0[j-ixm:j+ix,k-iym:k+iy] = 0
                        if(nooverlap):fno[j-ixm:j+ix,k-iym:k+iy] = 1
                        listpatches.append([j-ixm,k-iym,ixm+ix,iym+iy])
                    else:
                        f0[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz] = 0
                        if(nooverlap):fno[j-ixm:j+ix,k-iym:k+iy,l-izm:l+iz] = 1
                        listpatches.append([j-ixm,k-iym,l-izm,ixm+ix,iym+iy,izm+iz])
            self.listblocks.insert(1,listpatches)

    def sumpatch(self,listpatches,nx,ny,nz,dim):
        if dim==2:
            f = fzeros([nx,ny],'l')
            for patch in listpatches:
                j=patch[0]
                k=patch[1]
                ix = patch[2]
                iy = patch[3]
                f[j:j+ix,k:k+iy] += 1
        else:
            f = zeros([nx,ny,nz],'l')
            for patch in listpatches:
                j=patch[0]
                k=patch[1]
                l=patch[2]
                ix = patch[3]
                iy = patch[4]
                iz = patch[5]
                f[j:j+ix,k:k+iy,l:l+iz] += 1
        return f

    def setblocks(self):
        self.nblocks=0
        if self.solvergeom == w3d.XYZgeom:
            self.blocks.resetroot()
            mothergrid = self.blocks
        else:
            try:
                self.restorefrzgrid()
                mothergrid = frz.basegrid
                self.del_blocks2d()
            except:
                self.blocks.resetroot()
                mothergrid = self.blocks
        # --- Enforce garbage collection since the MR structures
        # --- have complex linkages and may not be immediately removed.
        gc.collect()
        xmin0 = w3d.xmmin; xmax0 = w3d.xmmax; dx = w3d.dx
        if self.solvergeom == w3d.XYZgeom:
            ymin0 = w3d.ymmin; ymax0 = w3d.ymmax; dy = w3d.dy
            zmin0 = w3d.zmmin; zmax0 = w3d.zmmax; dz = w3d.dz
        if self.solvergeom == w3d.RZgeom or self.solvergeom == w3d.XZgeom:
            ymin0 = w3d.zmmin; ymax0 = w3d.zmmax; dy = w3d.dz
        if self.solvergeom == w3d.XYgeom:
            ymin0 = w3d.ymmin; ymax0 = w3d.ymmax; dy = w3d.dy
        for ii,blocks in enumerate(self.listblocks[1:]):
            i = ii+1
            if i>1:
                if self.solvergeom == w3d.XYZgeom:
                    mothergrid = mothergrid.children[0]
                else:
                    try:
                        bg=frz.basegrid
                        mothergrid = mothergrid.down
                    except:
                        mothergrid = mothergrid.children[0]
            r=self.MRfact**i
            for patch in blocks:
                if (self.lremoveblockswithoutparticles and
                   not self.patchhasparticles(patch)): continue
                self.nblocks+=1
                if self.solvergeom == w3d.XYZgeom:
                    lower = nint(array(patch[:3])*r)
                    upper = lower + nint(array(patch[3:])*r)
                    print lower,upper,upper-lower
                    mothergrid.addchild(lower,upper,refinement=self.MRfact)
                else:
                    try:
                        bg=frz.basegrid
                        nx = nint(patch[2]*r)
                        ny = nint(patch[3]*r)
                        dxnew = dx/r
                        dynew = dy/r
                        dxmother = dxnew*self.MRfact
                        dymother = dynew*self.MRfact
                        xmin = xmin0 + patch[0]*dx
                        xmax = xmin  + nx*dxnew
                        ymin = ymin0 + patch[1]*dy
                        ymax = ymin  + ny*dynew
                        nx, xmin, xmax, t_xmin, t_xmax = self.add_transit(nx, xmin, xmax, dxmother, xmin0, xmax0)
                        ny, ymin, ymax, t_ymin, t_ymax = self.add_transit(ny, ymin, ymax, dymother, ymin0, ymax0)
                        add_subgrid(mothergrid.gid[0],nx,ny,dxnew,dynew,xmin,ymin,
                                    t_xmin*self.MRfact,t_xmax*self.MRfact,
                                    t_ymin*self.MRfact,t_ymax*self.MRfact)
                    except:
                        lower = nint(array(patch[:3])*r)
                        upper = lower + nint(array(patch[3:])*r)
                        mothergrid.addchild(lower,upper,refinement=self.MRfact)

        if self.solvergeom == w3d.XYZgeom:
            self.blocks.finalize()
        else:
            try:
                g = frz.basegrid
                for i in range(1,frz.ngrids):
                    try:
                        g = g.next
                    except:
                        g=g.down
                    self.blocks[g.gid[0]] = {'grid':g,'installed_conductors':[]}
            except:
                self.blocks.finalize()
            
    def add_transit(self, nx, xmin, xmax, dxmother, xmin0, xmax0):
        nt = self.ntransit
        n = self.MRfact

        xmin_try = xmin-nt*dxmother
        xmax_try = xmax+nt*dxmother
        t_xmin   = min(nt, max(0,nt-int((xmin0-xmin_try)/dxmother)))
        t_xmax   = min(nt, max(0,nt-int((xmax_try-xmax0)/dxmother)))
        xmin    -= t_xmin*dxmother
        xmax    += t_xmax*dxmother
        nx      += n*(t_xmin+t_xmax)

        return nx, xmin, xmax, t_xmin, t_xmax

    def del_blocks2d(self,g=None):
        if g==None: g=frz.basegrid
        try:
            self.del_blocks2d(g.next)
        except:
            try:
                self.del_blocks2d(g.down)
            except:
                pass
        if g is not frz.basegrid:
            id = g.gid[0]
            del_subgrid(id)
            self.blocks.__delitem__(id)
        else:
            frz.ngrids=1
            self.nblocks=0
            g.loc_part=g.gid[0]
            g.loc_part_fd=g.gid[0]

    def patchhasparticles(self,patch):
        if self.solvergeom==w3d.XYZgeom:
            i1,j1,k1 = patch[:3]
            i2,j2,k2 = array([i1,j1,k1]) + patch[3:]
            rho = self.blocks.rho[i1:i2+1,j1:j2+1,k1:k2+1]
        else:
            i1,j1 = patch[:2]
            i2,j2 = array([i1,j1]) + patch[2:]
            try:
                rho = frz.basegrid.rho[i1:i2+1,j1:j2+1]
            except:
                rho = self.blocks.rho[i1:i2+1,k1:k2+1]
        if maxnd(abs(rho)) == 0.: return 0
        else:                     return 1

    def generate(self,l_timing=0,l_allocate_blocks=1):
        """
      Generate AMR blocks based on values and gradients of self.f. If self.f is None, its default is
      the charge density. If self.f is null, an error message is raised.
        """
        # return if not time to generate a new set of blocks
        # Note that beforefs is still called.

        ifcond = w3d.AMRgenerate_periodicity==0
        if not ifcond:
            ifcond = top.it%w3d.AMRgenerate_periodicity != 0
        if ifcond:
            self.beforefs()
            return
        print 'generate grids at it = ',top.it

        # check if w3d.AMRlevels set properly and set defaut variables
        if w3d.AMRlevels<=0: raise Exception('Error: AMRTree.generate called with w3d.AMRlevels<=0/')
        self.MRfact   = w3d.AMRrefinement
        self.ntransit = w3d.AMRtransit

        # generate nbcells from self.f or use self.nbcells_user if provided
        if self.nbcells_user is None:
            if l_timing: starttime = time.clock()
            l_nbcellsnone=1
            # set self.f to the charge density array by default
            if self.f_user is not None:
                if callable(self.f_user):
                    self.f = self.f_user()
                else:
                    self.f = self.f_user
            elif self.f is None:
                if self.solvergeom==w3d.XYZgeom:
                    self.f = self.blocks.rho
                else:
                    try:
                        self.f = frz.basegrid.rho
                    except:
                        self.f = self.blocks.rho
                    
            # set cell dimensions of mother grid according to geometry
            dx = w3d.dx
            dz = w3d.dz
            if self.solvergeom==w3d.XYZgeom or self.solvergeom==w3d.XYgeom:
                dy = w3d.dy
            if self.solvergeom==w3d.XZgeom or self.solvergeom==w3d.RZgeom:
                dy = w3d.dz

            # set parameters controlling the automatic generation of blocks
            if w3d.AMRmaxlevel_density ==-1: w3d.AMRmaxlevel_density =w3d.AMRlevels
            if w3d.AMRmaxlevel_gradient==-1: w3d.AMRmaxlevel_gradient=w3d.AMRlevels
            self.Rdens = w3d.AMRrefinement**w3d.AMRmaxlevel_density
            self.Rgrad = w3d.AMRrefinement**w3d.AMRmaxlevel_gradient
            self.threshold = w3d.AMRthreshold_gradient
            if self.solvergeom==w3d.XYZgeom:
                self.maxcells_isolated_blocks = w3d.AMRmaxsize_isolated_blocks**3
            else:
                self.maxcells_isolated_blocks = w3d.AMRmaxsize_isolated_blocks**2
            self.nbcells=self.getnbcells(self.f,dx,dy,dz,self.Rdens,self.threshold,self.Rgrad,
                                         MRfact=self.MRfact,lmax=self.maxcells_isolated_blocks)
            self.f = None
            if l_timing:
                endtime = time.clock()
                print 'created nbcells in ',endtime-starttime,' seconds.'
        else:
            if callable(self.nbcells_user):
                self.nbcells = self.nbcells_user()
            else:
                self.nbcells = self.nbcells_user

        if self.nbcells is not None:
            # generate list of blocks from array nbcells
            if l_timing: starttime = time.clock()
            if self.solvergeom==w3d.XYZgeom:
                self.setlist(self.nbcells[:-1,:-1,:-1],w3d.AMRcoalescing,self.MRfact,true)
            else:
                self.setlist(self.nbcells[:-1,:-1],w3d.AMRcoalescing,self.MRfact,true)
            if l_timing:
                endtime = time.clock()
                print 'generated list in ',endtime-starttime,' seconds.'

            if not l_allocate_blocks:return

            # allocate blocks from list self.listblocks
            if l_timing: starttime = time.clock()
            self.setblocks()
            if l_timing:
                endtime = time.clock()
                print 'generated blocks in ',endtime-starttime,' seconds.'

            # clear inactive regions in each blocks
            if not w3d.AMRuse_inactive_regions:
                if l_timing: starttime = time.clock()
                if self.solvergeom==w3d.XYZgeom:
                    self.blocks.clearinactiveregions(self.nbcells)
                else:
                    try:
                        g = frz.basegrid
                        adjust_lpfd(self.nbcells,g.nr,g.nz,g.rmin,g.rmax,g.zmin,g.zmax)
                    except:
                        self.blocks.clearinactiveregions(self.nbcells)
                if l_timing:
                    endtime = time.clock()
                    print 'Cleared inactive regions in ',endtime-starttime,' seconds.'

        # set conductor data
        if l_timing: starttime = time.clock()
        if self.solvergeom==w3d.XYZgeom:
            for cond,dfill in zip(self.conductors,self.conductorsdfill):
                self.blocks.installconductor(cond,dfill=dfill)
        else:
            try:
                bg = frz.basegrid
                # this loop is needed so that the grids are correctly registered when installing conductors
                if self.nblocks>0:
                    g = frz.basegrid
                    for idummy in range(frz.ngrids-1):
                        rdummy=g.nr
                        try:
                            g=g.next
                        except:
                            try:
                                g=g.down
                            except:
                                pass
                # install conductors
                for cond,dfill in zip(self.conductors,self.conductorsdfill):
                    for block in self.blocks.itervalues():
                        g=block['grid']
                        if g is not frz.basegrid:
                            if cond not in block['installed_conductors']:
                                installconductors(cond,nx=g.nr,ny=0,nzlocal=g.nz,nz=g.nz,
                                                  xmmin=g.xmin,xmmax=g.xmax,
                                                  zmmin=g.zmin,zmmax=g.zmax,dfill=dfill,
                                                  gridrz=g)
                                block['installed_conductors'].append(cond)
                get_cond_rz(1)
            except:
                for cond,dfill in zip(self.conductors,self.conductorsdfill):
                    self.blocks.installconductor(cond,dfill=dfill)
        if l_timing:
            endtime = time.clock()
            print 'generated conductors in ',endtime-starttime,' seconds.'

        if l_timing: starttime = time.clock()
        # load charge density on new set of blocks
        if self.solvergeom == w3d.XYZgeom:
            # --- Call the solvers loadrho routine directly since AMR instance
            # --- is already the registered solver.
            self.blocks.loadrho(lzero=true,lfinalize_rho=true,lrootonly=0)
        else:
            try:
                bg=frz.basegrid
                loadrho()
            except:
                self.blocks.loadrho(lzero=true,lfinalize_rho=true,lrootonly=0)
        if l_timing:
            endtime = time.clock()
            print 'loaded rho in ',endtime-starttime,' seconds.'

        self.beforefs()
        print 'Generated ',self.nblocks,' blocks.'

    def draw_blocks2d(self,level=None,color='black',width=1.,allmesh=0,f=1):
        for ii,blocks in enumerate(self.listblocks[1:]):
            i=ii+1
            if level is None or level==i:
                for patch in blocks:
                    j=patch[0]
                    k=patch[1]
                    l = patch[2]
                    h = patch[3]
                    r=(float(self.MRfact)**i)/f
                    nx = nint(l*r)
                    ny = nint(h*r)
                    if self.solvergeom != w3d.RZgeom:
                        xmin=w3d.xmmin
                        ymin=w3d.ymmin
                        dx=w3d.dx
                        dy=w3d.dy
                        if(allmesh):
                            self.draw_mesh(nx,ny,xmin+j*dx,ymin+k*dy,dx/r,dy/r,color=self.colors[i%len(self.colors)],width=width)
                        else:
                            self.draw_box(ymin+k*dy, ymin+k*dy+h*dy, xmin+j*dx, xmin+j*dx+l*dx, color=self.colors[i%len(self.colors)],width=width)
                    else:
                        xmin=w3d.xmmin
                        ymin=w3d.zmmin
                        dx=w3d.dx
                        dy=w3d.dz
                        if(allmesh):
                            self.draw_mesh(ny,nx,ymin+k*dy,xmin+j*dx,dy/r,dx/r,color=self.colors[i%len(self.colors)],width=width)
                        else:
                            self.draw_box(xmin+j*dx, xmin+j*dx+l*dx, ymin+k*dy, ymin+k*dy+h*dy, color=self.colors[i%len(self.colors)],width=width)

    def draw_mesh(self,nx,ny,xmin,ymin,dx,dy,color='black',width=1):
        x = xmin+arange(nx+1)*dx
        y = ymin+arange(ny+1)*dy
        xxmin = xmin*ones(ny+1)
        yymin = ymin*ones(nx+1)
        pldj(x,yymin,x,yymin+ny*dy,color=color,width=width)
        pldj(xxmin,y,xxmin+nx*dx,y,color=color,width=width)


    def draw_box(self,rmin, rmax, zmin, zmax, color='blue',width=1):
        pldj([zmin,zmin,zmin,zmax],
             [rmin,rmax,rmin,rmin],
             [zmax,zmax,zmin,zmax],
             [rmin,rmax,rmax,rmax],color=color,width=width)

    def createdxobject(self,kwdict={},**kw):
        """
      Create DX object drawing the object.
      - withguards=1: when true, the guard cells are included in the bounding box

        """
        kw.update(kwdict)
        withguards = kw.get('withguards',1)
        level = kw.get('level',None)
        xmin = kw.get('xmin')
        xmax = kw.get('xmax')
        ymin = kw.get('ymin')
        ymax = kw.get('ymax')
        zmin = kw.get('zmin')
        zmax = kw.get('zmax')
        if xmin is None:xmin=w3d.xmmin
        if xmax is None:xmax=w3d.xmmax
        if ymin is None:ymin=w3d.ymmin
        if ymax is None:ymax=w3d.ymmax
        if zmin is None:zmin=w3d.zmmin
        if zmax is None:zmax=w3d.zmmax
        dxlist = []
        for i,blocks in enumerate(self.listblocks[1:]):
            if level is None or level==i:
                for patch in blocks:
                    xminp = w3d.xmmin + patch[0]*w3d.dx
                    xmaxp =     xminp + patch[3]*w3d.dx
                    yminp = w3d.ymmin + patch[1]*w3d.dy
                    ymaxp =     yminp + patch[4]*w3d.dy
                    zminp = w3d.zmmin + patch[2]*w3d.dz
                    zmaxp =     zminp + patch[5]*w3d.dz
                    if( not( (xminp>xmax) or (xmaxp<xmin) or \
                             (yminp>ymax) or (ymaxp<ymin) or \
                             (zminp>zmax) or (zmaxp<zmin))):
                        dxlist.append(Opyndx.viewboundingbox(
                                                      max(xmin,xminp),min(xmax,xmaxp),
                                                      max(ymin,yminp),min(ymax,ymaxp),
                                                      max(zmin,zminp),min(zmax,zmaxp),
                                                      self.colors[i]))
        self.dxobject = Opyndx.DXCollection(*dxlist)
    def draw(self,level=None,xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None, \
             color='black',width=1.,allmesh=0):
        if self.solvergeom==w3d.XYZgeom:#MR:
            self.createdxobject(level=level,
                                xmin=xmin,xmax=xmax,
                                ymin=ymin,ymax=ymax,
                                zmin=zmin,zmax=zmax)
            Opyndx.DXImage(self)
        else:
            self.draw_blocks2d(level=level,color=color,width=width,allmesh=allmesh)


def draw_mesh(nx,ny,xmin,ymin,dx,dy,color='black',width=1):
    x = xmin+arange(nx+1)*dx
    y = ymin+arange(ny+1)*dy
    xxmin = xmin*ones(ny+1)
    yymin = ymin*ones(nx+1)
    pldj(x,yymin,x,yymin+ny*dy,color=color,width=width)
    pldj(xxmin,y,xxmin+nx*dx,y,color=color,width=width)

def plphirz(grid=None,which='phi',cmin=None,cmax=None,
                 border=1,bordercolor=None,borderwidth=1,
                 mesh=0,meshcolor=None,meshwidth=1,meshr=1,
                 siblings=1,children=1,firstcall=1,level=1,maxlevel=0,
                 delay=0,transit=0,l_transpose=0):
    if grid is None:
        g = frz.basegrid
    else:
        g = grid

    colors = ['red','blue','yellow','green','cyan','magenta','white']
    i=g.levelref
    if bordercolor is None:
        bordercoloru=colors[i%len(colors)]
    else:
        bordercoloru=bordercolor
    if meshcolor is None:
        meshcoloru=colors[i%len(colors)]
    else:
        meshcoloru=meshcolor

    zmin = g.zmin-0.5*g.dz
    rmin = g.rmin-0.5*g.dr
    zmax = zmin+(g.nz+1)*g.dz
    rmax = rmin+(g.nr+1)*g.dr
    dr=g.dr
    dz=g.dz
    nr=g.nr
    nz=g.nz

    if not transit:
        zmin = zmin+g.transit_min_z*g.dz
        zmax = zmax-g.transit_max_z*g.dz
        rmin = rmin+g.transit_min_r*g.dr
        rmax = rmax-g.transit_max_r*g.dr
        nz -= 2*g.transit_min_z
        nr -= 2*g.transit_min_r
    if(which=='phi'):
        f = g.phi[g.nguardx:-g.nguardx,g.nguardz:-g.nguardz]
    if(which=='rho'):
        f = g.rho
    if(which=='lp'):
        f = g.loc_part
    if(which=='lpfd'):
        f = g.loc_part_fd
    if(firstcall):
        if(which=='phi'):
            if cmin is None:cmin = minnd(frz.basegrid.phi)
            if cmax is None:cmax = maxnd(frz.basegrid.phi)
        elif(which=='rho'):
            if cmin is None:cmin = minnd(frz.basegrid.rho)
            if cmax is None:cmax = maxnd(frz.basegrid.rho)
        else:
            cmin=0
            cmax=frz.ngrids
    if not transit:
        f = f[g.transit_min_r:g.nr+1-g.transit_max_r,g.transit_min_z:g.nz+1-g.transit_max_z]
    if not l_transpose:
        pli(f,zmin,rmin,zmax,rmax,cmin=cmin,cmax=cmax)
    else:
        pli(transpose(f),rmin,zmin,rmax,zmax,cmin=cmin,cmax=cmax)
    if(mesh):
#        nr = nint(float(g.nr)/meshr)
#        nz = nint(float(g.nz)/meshr)
#        dr = g.dr*meshr
#        dz = g.dz*meshr
        nr/=meshr
        nz/=meshr
        dr*=meshr
        dz*=meshr
        if not l_transpose:
#          draw_mesh(nz,nr,g.zmin,g.rmin,dz,dr,color=meshcoloru,width=meshwidth)
            draw_mesh(nz+1,nr+1,zmin,rmin,dz,dr,color=meshcoloru,width=meshwidth)
        else:
#          draw_mesh(nr,nz,g.rmin,g.zmin,dr,dz,color=meshcoloru,width=meshwidth)
            draw_mesh(nr+1,nz+1,rmin,zmin,dr,dz,color=meshcoloru,width=meshwidth)
    if(border):
        if not l_transpose:
            draw_box(rmin, rmax, zmin, zmax, color=bordercoloru,width=borderwidth)
        else:
            draw_box(zmin, zmax, rmin, rmax, color=bordercoloru,width=borderwidth)
    time.sleep(delay)
    pyg_pending()
    pyg_idler()
    if(siblings):
        try:
            plphirz(g.next,which,cmin,cmax,border,bordercolor,borderwidth,mesh,meshcolor,meshwidth,meshr,
                       siblings,children=0,firstcall=0,level=level,maxlevel=maxlevel,delay=delay,transit=transit,
                       l_transpose=l_transpose)
        except:
            pass
    if(children):
        if maxlevel==0 or level<maxlevel:
            try:
                plphirz(g.down,which,cmin,cmax,border,bordercolor,borderwidth,mesh,meshcolor,meshwidth,meshr,
                           siblings,children,firstcall=0,level=level+1,maxlevel=maxlevel,delay=delay,transit=transit,
                          l_transpose=l_transpose)
            except:
                pass
    if(firstcall):
        colorbar(cmin,cmax,view=plsys())

def plrhorz(**args):
    plphirz(which='rho',**args)

def plcondrz(grid=None,border=1,bordercolor='yellow',mesh=0,meshcolor='white',meshr=1,
                 siblings=1,children=1,firstcall=1,level=1,maxlevel=0,delay=0):
    if grid is None:
        g = frz.basegrid
    else:
        g = grid

    zmin = g.zmin-g.dz
    rmin = g.rmin-g.dr
    zmax = g.zmax+g.dz
    rmax = g.rmax+g.dr
    b = g.bndfirst
    for i in range(b.nb_conductors):
        if(i==0):
            c = b.cndfirst
        else:
            c = c.next
        color=red
        for ic in range(c.nbbnd):
            if ic>=c.nbbndred:color=green
            z=zmin+c.kk[ic]*g.dz
            x=rmin+c.jj[ic]*g.dr
            if(c.dxm[ic]<g.dr):pldj([z],[x],[z],[x-c.dxm[ic]],color=color)
            if(c.dxp[ic]<g.dr):pldj([z],[x],[z],[x+c.dxp[ic]],color=color)
            if(c.dzm[ic]<g.dz):pldj([z],[x],[z-c.dzm[ic]],[x],color=color)
            if(c.dzp[ic]<g.dz):pldj([z],[x],[z+c.dzp[ic]],[x],color=color)

    if(mesh):
        nr = nint(float(g.nr)/meshr)
        nz = nint(float(g.nz)/meshr)
        dr = g.dr*meshr
        dz = g.dz*meshr
        draw_mesh(nz,nr,g.zmin,g.rmin,dz,dr,color=meshcolor)
    if(border):
        draw_box(rmin, rmax, zmin, zmax, color=bordercolor)
    time.sleep(delay)
    pyg_pending()
    pyg_idler()
    if(siblings):
        try:
            plcondrz(g.next,border,bordercolor,mesh,meshcolor,meshr,
                       siblings,children=0,firstcall=0,level=level,maxlevel=maxlevel,delay=delay)
        except:
            pass
    if(children):
        if maxlevel==0 or level<maxlevel:
            try:
                plcondrz(g.down,border,bordercolor,mesh,meshcolor,meshr,
                           siblings,children,firstcall=0,level=level+1,maxlevel=maxlevel,delay=delay)
            except:
                pass

def draw_box(rmin, rmax, zmin, zmax, color='blue',width=1):
    pldj([zmin,zmin,zmin,zmax],
          [rmin,rmax,rmin,rmin],
          [zmax,zmax,zmin,zmax],
          [rmin,rmax,rmax,rmax],color=color,width=width)



def walkfrzgrid(bnd=0,cnd=0,base=None):
    """Walk through the frz grid structure, returning each element once.
   - bnd=0: when true, the objects returned will be the BNDtypes
   - cnd=0: when true, the objects returned will be the CONDtypes
  Note that the base argument should never be specified.
    """
    # --- Create list which hold the objects already returned. This is used
    # --- to make sure that an object is not returned multiple times.
    try: walkfrzgrid.glist
    except AttributeError: walkfrzgrid.glist = []

    # --- If this is the top call, set the base.
    # --- Also make sure that bnd is set if cnd is set.
    if base is None:
        base = frz.basegrid
        if cnd: bnd = 1

    if not bnd and not cnd:
        # --- If not return a BND or CONDtype, return the base
        # --- Note that base may still be a BND or COND type
        # --- passed in from above.
        if base not in walkfrzgrid.glist:
        # --- Only return it if it has not been returned already.
            walkfrzgrid.glist.append(base)
            yield base
    elif bnd:
        # --- If BND, walk through the BND types, returning each one.
        # --- cnd is passed in since the BND must be walked through
        # --- to return the CONDs.
        for b in walkfrzgrid(bnd=0,cnd=cnd,base=base.bndfirst):
            if b is not None: yield b
            else: break
    else:
        # --- Walk through the CONDs
        if base.nb_conductors > 0:
            for c in walkfrzgrid(bnd=0,cnd=0,base=base.cndfirst):
                if c is not None: yield c
                else: break

    try:
        # --- Check for a level down. Note that only grids have a down
        # --- attribute, hence the try/except.
        # --- If down is unallocated, getpyobject returns None.
        if base.getpyobject('down') is not None:
            for g in walkfrzgrid(bnd=bnd,cnd=cnd,base=base.down):
                if g is not None: yield g
                else: break
    except AttributeError:
        pass

    # --- Check for a next object. This applies to grids, BNDs, and CONDs.
    if base.getpyobject('next') is not None:
        for g in walkfrzgrid(bnd=bnd,cnd=cnd,base=base.next):
            if g is not None: yield g
            else: break

    # --- If this is the end of the top call, then reset glist.
    # --- This is way this routine will only work (multiple times)
    # --- if the starting point is frz.basegrid.
    if base == frz.basegrid: walkfrzgrid.glist = []




# --- This can only be done after AMRTree is defined.
try:
    psyco.bind(AMRTree)
except NameError:
    pass
