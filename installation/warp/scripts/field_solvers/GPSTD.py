""" Class for FFT-based Pseudo-Spectral solver """
import numpy as np
import math
import copy
import os

try:
    # Try to import fortran wrapper of FFTW
    import fastfftforpy as fftpy
    import fastfftpy as fstpy
    fst=fstpy.fastfft
    fft=fftpy
    l_fftw_fort=True
except:
    fft = np.fft
    l_fftw_fort=False

print('l_fftw_fort',l_fftw_fort)

def iszero(f):
    """
    Returns True if f==0.

    Parameters
    ----------
    f : any
        Input argument, can be of any type and shape.

    Returns
    -------
    val : bool
        True if f==0., False if it is not.

    """
    if type(f) is type(0.):
        if f==0.:
            return True
    return False

def getmatcompress(mat):
    matcompress = copy.deepcopy(mat)
    matlist = []
    matlistindices = []
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if mat[i][j] in matlist:
                l = matlist.index(mat[i][j])
                matcompress[i][j]=matlistindices[l]
            else:
                matlist.append(mat[i][j])
#              matlistabs.append(abs(mat[i][j]))
                matlistindices.append([i,j])
                matcompress[i][j]=None

    return None#matcompress

def multmat(matrix1,matrix2,matcompress=None):
    """
    Returns dot product of two 2-D lists.
    Each element of the list is a scalar or a real 1-D, 2-D or 3-D array.

    Parameters
    ----------
    matrix1 : 2-D matrix list
        First argument.
    matrix2 : 2-D matrix list
        Second argument.
    matcompress : 2-D matrix list for compression
        Optional third argument

    Returns
    -------
    output : 2-D matrix list
        Returns the dot product of matrix1 and matrix 2.

    """
        # Matrix multiplication
    if len(matrix1[0]) != len(matrix2):
    # Check matrix dimensions
        print 'Matrices must be m*n and n*p to multiply!'
    else:
        # Multiply if correct dimensions
        new_matrix = copy.deepcopy(matrix2)
        for i in range(len(matrix2)):
            for j in range(len(matrix2[0])):
                new_matrix[i][j]=0.
        for i in range(len(matrix2)):
            for j in range(len(matrix1[0])):
                if matcompress is not None:
                    if matcompress[i][j] is not None:
                        ii = matcompress[i][j][0]
                        jj = matcompress[i][j][1]
                        new_matrix[i][j]=new_matrix[ii][jj]
                        doadd=0
                    else:
                        doadd=1
                else:
                    doadd=1
                if doadd:
                    for k in range(len(matrix1)):
                        if not iszero(matrix2[i][k]) and not iszero(matrix1[k][j]):
                            new_matrix[i][j] += matrix2[i][k]*matrix1[k][j]
        return new_matrix

def exp_by_squaring_matrixlist(x, n, matcompress=None):
    if n < 1:
        raise Exception('Error in exp_by_squaring_matrixlist: n<1.')

    if n == 1:return x
    y = np.identity(len(x)).tolist()
    while n > 1:
      if float(n/2)==float(n)/2: # if n is even then
        x = multmat(x,x,matcompress=matcompress)
        n /= 2
      else:
        y = multmat(x,y,matcompress=matcompress)
        x = multmat(x,x,matcompress=matcompress)
        n = (n-1)/2
    return multmat(x,y,matcompress=matcompress)

def FD_weights(z,n,m):
    """
 adapted from Matlab code from Fornberg (1998)
 Calculates FD weights. The parameters are:
  z   location np.where approximations are to be accurate,
  n   number of grid points,
  m   highest derivative that we want to find weights for
  c   array size m+1,lentgh(x) containing (as output) in
      successive rows the weights for derivatives 0,1,...,m.
    """

    x = np.arange(-n/2+1,n/2+1)*1.

    c=np.zeros([m+1,n]); c1=1.; c4=x[0]-z; c[0,0]=1.;
    for i in range(1,n):
        mn=min(i+1,m+1); c2=1.; c5=c4; c4=x[i]-z;
        for j in range(0,i-0):
            c3=x[i]-x[j];  c2=c2*c3;
            if j==i-1:
                c[1:mn,i]=c1*(np.arange(1,mn)*c[0:mn-1,i-1]-c5*c[1:mn,i-1])/c2;
                c[0,i]=-c1*c5*c[0,i-1]/c2;
            c[1:mn,j]=(c4*c[1:mn,j]-np.arange(1,mn)*c[0:mn-1,j])/c3;
            c[0,j]=c4*c[0,j]/c3;
        c1=c2;

    return c

def FD_weights_hvincenti(p,l_staggered=False):
    # --- from Henri Vincenti's formulas
    factorial = math.factorial

    c = np.zeros(p//2)
    for i in range(p//2):
        l=i+1
        if l_staggered:
            lognumer = math.log(16.)*(1.-p/2.)+math.log(factorial(p-1.))*2
            logdenom = math.log(2.*l-1.)*2.+math.log(factorial(p/2.+l-1.))+math.log(factorial(p/2.-l))+2*math.log(factorial(p/2.-1.))
        else:
            lognumer = math.log(factorial(p/2.))*2
            logdenom = math.log(factorial(p/2.+l))+math.log(factorial(p/2.-l))+math.log(l)
        c[i] = (-1.)**(l+1)*np.exp(lognumer-logdenom)
    return c

class Fourier_Space():

    __flaginputs__ = {'nx':1,'ny':1,'nz':1,
                      'nxguard':0,'nyguard':0,'nzguard':0,
                      'norderx':np.inf,'nordery':np.inf,'norderz':np.inf,
                      'dt':1.,'dx':1.,'dy':1.,'dz':1.,
                      'l_staggered':False,
                      'bc_periodic':[0,0,0],
                      'l_fftw': l_fftw_fort,
                      'nthreads': None}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(Fourier_Space.__flaginputs__,kw)



        # Dimensions of real space arrays
        nx = self.nx
        ny = self.ny
        nz = self.nz
        if not self.bc_periodic[0]:
            nx+=2*self.nxguard
        if not self.bc_periodic[1]:
            ny+=2*self.nyguard
        if not self.bc_periodic[2]:
            nz+=2*self.nzguard
        self.dim_r2c=dims_r2c=[]
        if nx>1:
            dims_r2c+=[nx]
        if ny>1:
            dims_r2c+=[ny]
        if nz>1:
            dims_r2c+=[nz]

        # Dimensions of Fourier space arrays
        if (len(dims_r2c)==3):
            nxf = nx//2+1
            nyf = ny
            nzf = nz
        elif (len(dims_r2c) ==2):
            if (nx>1):
                nxf=nx//2+1
                nyf=ny
                nzf=nz
            else:
                nxf=nx
                nyf=ny//2
                nzf=nz
        self.dims=dims=[]
        if nxf>1:
            dims+=[nxf]
        if nyf>1:
            dims+=[nyf]
        if nzf>1:
            dims+=[nzf]

        # --- sets kx, ky, kz, k
        if nx>1:
            self.kxn = np.ones(dims)
            self.kx_unmod = np.ones(dims)
            kxunit = 2.*np.pi*(self.rfftfreq(nx))/self.dx
            kxunit_mod = kxunit.copy()
            self.kxunit = kxunit
            if self.norderx is not np.inf:
                xcoefs=0.
                if self.l_staggered:
#                    xc = 0.5
                    xi = 2
                else:
#                    xc = 0.
                    xi = 1
#                w = FD_weights(xc, self.norderx+1,1)[-1,self.norderx/2+1:]
                w = FD_weights_hvincenti(self.norderx,l_staggered=self.l_staggered)

                for i in range(self.norderx//2):
                    xcoefs+=w[i]*2*np.sin(kxunit*(xi*i+1)*self.dx/xi)

                kxunit_mod*=xcoefs/np.where(kxunit==0.,1.,kxunit*self.dx)

        if ny>1:
            self.kyn = np.ones(dims)
            self.ky_unmod = np.ones(dims)
            if (nx == 1): # Case 2D YZ
                kyunit = 2.*np.pi*(self.rfftfreq(ny))/self.dy
            else:         # CASE 3D XYZ or 2D XY
                kyunit = 2.*np.pi*(self.fftfreq(ny))/self.dy
            kyunit_mod = kyunit.copy()
            self.kyunit = kyunit
            if self.nordery is not np.inf:
                ycoefs=0.
                if self.l_staggered:
                    yc = 0.5
                    yi = 2
                else:
                    yc = 0.
                    yi = 1
#                w = FD_weights(yc, self.nordery+1,1)[-1,self.nordery/2+1:]
                w = FD_weights_hvincenti(self.nordery,l_staggered=self.l_staggered)

                for i in range(self.nordery//2):
                    ycoefs+=w[i]*2*np.sin(kyunit*(yi*i+1)*self.dy/yi)

                kyunit_mod*=ycoefs/np.where(kyunit==0.,1.,kyunit*self.dy)

        if nz>1:
            self.kzn = np.ones(dims)
            self.kz_unmod = np.ones(dims)
            kzunit = 2.*np.pi*(self.fftfreq(nz))/self.dz
            kzunit_mod = kzunit.copy()
            self.kzunit = kzunit
            if self.norderz is not np.inf:
                zcoefs=0.
                if self.l_staggered:
                    zc = 0.5
                    zi = 2
                else:
                    zc = 0.
                    zi = 1
#                w = FD_weights(zc, self.norderz+1,1)[-1,self.norderz/2+1:]
                w = FD_weights_hvincenti(self.norderz,l_staggered=self.l_staggered)

                for i in range(self.norderz//2):
                    zcoefs+=w[i]*2*np.sin(kzunit*(zi*i+1)*self.dz/zi)

                kzunit_mod*=zcoefs/np.where(kzunit==0.,1.,kzunit*self.dz)

        if len(dims)==3:
            for k in range(nzf):
                for j in range(nyf):
                    self.kxn[:,j,k] *= kxunit_mod
                    self.kx_unmod[:,j,k] *= kxunit

            for i in range(nxf):
                for k in range(nzf):
                    self.kyn[i,:,k] *= kyunit_mod
                    self.ky_unmod[i,:,k] *= kyunit

            for i in range(nxf):
                for j in range(nyf):
                    self.kzn[i,j,:] *= kzunit_mod
                    self.kz_unmod[i,j,:] *= kzunit

            self.kx=self.kxn.copy()
            self.ky=self.kyn.copy()
            self.kz=self.kzn.copy()
            self.k = np.sqrt(self.kxn*self.kxn+self.kyn*self.kyn+self.kzn*self.kzn)
            self.kmag = abs(np.where(self.k==0.,1.,self.k))
            self.kxn = self.kxn/self.kmag
            self.kyn = self.kyn/self.kmag
            self.kzn = self.kzn/self.kmag


            j = 1j
            if self.l_staggered:
                self.kxmn = self.kxn*np.exp(-j*self.kx_unmod*self.dx/2)
                self.kxpn = self.kxn*np.exp( j*self.kx_unmod*self.dx/2)
                self.kymn = self.kyn*np.exp(-j*self.ky_unmod*self.dy/2)
                self.kypn = self.kyn*np.exp( j*self.ky_unmod*self.dy/2)
                self.kzmn = self.kzn*np.exp(-j*self.kz_unmod*self.dz/2)
                self.kzpn = self.kzn*np.exp( j*self.kz_unmod*self.dz/2)
                self.kxm = self.kx*np.exp(-j*self.kx_unmod*self.dx/2)
                self.kxp = self.kx*np.exp( j*self.kx_unmod*self.dx/2)
                self.kym = self.ky*np.exp(-j*self.ky_unmod*self.dy/2)
                self.kyp = self.ky*np.exp( j*self.ky_unmod*self.dy/2)
                self.kzm = self.kz*np.exp(-j*self.kz_unmod*self.dz/2)
                self.kzp = self.kz*np.exp( j*self.kz_unmod*self.dz/2)
            else:
                self.kxmn = self.kxn
                self.kxpn = self.kxn
                self.kymn = self.kyn
                self.kypn = self.kyn
                self.kzmn = self.kzn
                self.kzpn = self.kzn
                self.kxm = self.kx
                self.kxp = self.kx
                self.kym = self.ky
                self.kyp = self.ky
                self.kzm = self.kz
                self.kzp = self.kz

        if len(dims)==2:
            j = 1j
            # --- 2D YZ
            if nx==1:
                for i in range(nzf):
                    self.kyn[:,i] *= kyunit_mod
                    self.ky_unmod[:,i] *= kyunit
                for i in range(nyf):
                    self.kzn[i,:] *= kzunit_mod
                    self.kz_unmod[i,:] *= kzunit
                self.kx_unmod = 0.
                self.ky=self.kyn.copy()
                self.kz=self.kzn.copy()
                self.k = np.sqrt(self.kyn*self.kyn+self.kzn*self.kzn)
                self.kmag = abs(np.where(self.k==0.,1.,self.k))
                self.kyn = self.kyn/self.kmag
                self.kzn = self.kzn/self.kmag
                self.kxp = 0.
                self.kxm = 0.
                self.kxpn = 0.
                self.kxmn = 0.
                if self.l_staggered:
                    self.kymn = self.kyn*np.exp(-j*self.ky_unmod*self.dy/2)
                    self.kypn = self.kyn*np.exp( j*self.ky_unmod*self.dy/2)
                    self.kzmn = self.kzn*np.exp(-j*self.kz_unmod*self.dz/2)
                    self.kzpn = self.kzn*np.exp( j*self.kz_unmod*self.dz/2)
                    self.kym = self.ky*np.exp(-j*self.ky_unmod*self.dy/2)
                    self.kyp = self.ky*np.exp( j*self.ky_unmod*self.dy/2)
                    self.kzm = self.kz*np.exp(-j*self.kz_unmod*self.dz/2)
                    self.kzp = self.kz*np.exp( j*self.kz_unmod*self.dz/2)
                else:
                    self.kymn = self.kyn
                    self.kypn = self.kyn
                    self.kzmn = self.kzn
                    self.kzpn = self.kzn
                    self.kym = self.ky
                    self.kyp = self.ky
                    self.kzm = self.kz
                    self.kzp = self.kz

            if ny==1:
            # --- 2D XZ
                for i in range(nzf):
                    self.kxn[:,i] *= kxunit_mod
                    self.kx_unmod[:,i] *= kxunit
                for i in range(nxf):
                    self.kzn[i,:] *= kzunit_mod
                    self.kz_unmod[i,:] *= kzunit
                self.ky_unmod = 0.
                self.kx=self.kxn.copy()
                self.kz=self.kzn.copy()
                self.k = np.sqrt(self.kxn*self.kxn+self.kzn*self.kzn)
                self.kmag = abs(np.where(self.k==0.,1.,self.k))
                self.kxn = self.kxn/self.kmag
                self.kzn = self.kzn/self.kmag
                self.kyp = 0.
                self.kym = 0.
                self.kypn = 0.
                self.kymn = 0.
                if self.l_staggered:
                    self.kxmn = self.kxn*np.exp(-j*self.kx_unmod*self.dx/2)
                    self.kxpn = self.kxn*np.exp( j*self.kx_unmod*self.dx/2)
                    self.kzmn = self.kzn*np.exp(-j*self.kz_unmod*self.dz/2)
                    self.kzpn = self.kzn*np.exp( j*self.kz_unmod*self.dz/2)
                    self.kxm = self.kx*np.exp(-j*self.kx_unmod*self.dx/2)
                    self.kxp = self.kx*np.exp( j*self.kx_unmod*self.dx/2)
                    self.kzm = self.kz*np.exp(-j*self.kz_unmod*self.dz/2)
                    self.kzp = self.kz*np.exp( j*self.kz_unmod*self.dz/2)
                else:
                    self.kxmn = self.kxn
                    self.kxpn = self.kxn
                    self.kzmn = self.kzn
                    self.kzpn = self.kzn
                    self.kxm = self.kx
                    self.kxp = self.kx
                    self.kzm = self.kz
                    self.kzp = self.kz

            if nz==1:
            # --- 2D XY
                for i in range(nyf):
                    self.kxn[:,i] *= kxunit_mod
                    self.kx_unmod[:,i] *= kxunit
                for i in range(nxf):
                    self.kyn[i,:] *= kyunit_mod
                    self.ky_unmod[i,:] *= kyunit
                self.kz_unmod = 0.
                self.kx=self.kxn.copy()
                self.ky=self.kyn.copy()
                self.k = np.sqrt(self.kxn*self.kxn+self.kyn*self.kyn)
                self.kmag = abs(np.where(self.k==0.,1.,self.k))
                self.kxn = self.kxn/self.kmag
                self.kyn = self.kyn/self.kmag
                self.kzp = 0.
                self.kzm = 0.
                self.kzpn = 0.
                self.kzmn = 0.
                if self.l_staggered:
                    self.kxmn = self.kxn*np.exp(-j*self.kx_unmod*self.dx/2)
                    self.kxpn = self.kxn*np.exp( j*self.kx_unmod*self.dx/2)
                    self.kymn = self.kyn*np.exp(-j*self.ky_unmod*self.dy/2)
                    self.kypn = self.kyn*np.exp( j*self.ky_unmod*self.dy/2)
                    self.kxm = self.kx*np.exp(-j*self.kx_unmod*self.dx/2)
                    self.kxp = self.kx*np.exp( j*self.kx_unmod*self.dx/2)
                    self.kym = self.ky*np.exp(-j*self.ky_unmod*self.dy/2)
                    self.kyp = self.ky*np.exp( j*self.ky_unmod*self.dy/2)
                else:
                    self.kxmn = self.kxn
                    self.kxpn = self.kxn
                    self.kymn = self.kyn
                    self.kypn = self.kyn
                    self.kxm = self.kx
                    self.kxp = self.kx
                    self.kym = self.ky
                    self.kyp = self.ky
        # compute FFTW plans if FFTW loaded
        if (self.l_fftw):
            if self.nthreads is None:
                self.nthreads=int(os.getenv('OMP_NUM_THREADS',1))
        # Init plans
        self.plan_rfftn={}
        self.plan_irfftn={}
        self.planj_rfftn=None
        self.planj_irfftn=None

    def create_plan_rfftn(self, dims):
        if self.l_fftw:
            return fftpy.compute_plan_rfftn(dims,nthreads=self.nthreads, plan_opt=fst.fftw_measure)
        else:
            return None

    def create_plan_irfftn(self, dims):
        if self.l_fftw:
            return fftpy.compute_plan_irfftn(dims,nthreads=self.nthreads, plan_opt=fst.fftw_measure)
        else:
            return None

    def fftn(self, a, field_out=None, plan=None):
        if (self.l_fftw):
            return fftpy.fftn(a,plan=plan,nthreads=self.nthreads, field_out=field_out)
        else:
            return np.fft.fftn(a)

    def rfftn(self, a, field_out=None, plan =None):
        if (self.l_fftw):
            return fftpy.rfftn(a,plan=plan,nthreads=self.nthreads, field_out=field_out)
        else:
            axes = np.arange(0,a.ndim)
            last_axe = axes[-1]
            axes[-1] = axes[0]
            axes[0]  = last_axe
            return np.fft.rfftn(a, axes=axes)

    def ifftn(self, a, field_out=None, plan=None):
        if (self.l_fftw):
            return fftpy.ifftn(a,plan=plan,nthreads=self.nthreads, field_out=field_out)
        else:
            return np.fft.ifftn(a)
    def irfftn(self, a, dims, field_out=None, plan=None):
        if (self.l_fftw):
            return fftpy.irfftn(a,dims, plan=plan,nthreads=self.nthreads, field_out=field_out)
        else:
            axes = np.arange(0,a.ndim)
            last_axe = axes[-1]
            axes[-1] = axes[0]
            axes[0]  = last_axe
            rdims = dims.copy()
            rdims[0] = dims[-1]
            rdims[-1]=dims[0]
            return np.fft.irfftn(a, axes=axes, s=rdims)

    def fft(self, a, axis=0):
        if (self.l_fftw):
            return fftpy.fft(a, axis=axis,nthreads=self.nthreads)
        else:
            return np.fft.fft(a, axis=axis)

    def ifft(self, a, axis=0):
        if (self.l_fftw):
            return fftpy.ifft(a, axis=axis,nthreads=self.nthreads)
        else:
            return np.fft.ifft(a, axis=axis)

    def fftfreq(self, a):
        if (self.l_fftw):
            return fftpy.fftfreq(a)
        else:
            return np.fft.fftfreq(a)

    def rfftfreq(self, a):
        if (self.l_fftw):
            return fftpy.rfftfreq(a)
        else:
            return np.fft.rfftfreq(a)

    def processdefaultsfromdict(self,dict,kw):
        for name,defvalue in dict.iteritems():
            if name not in self.__dict__:
                self.__dict__[name] = kw.get(name,defvalue)
            if name in kw: del kw[name]

    def get_ius(self):
        if self.bc_periodic[0]:
            ixl = self.nxguard
            ixu = max(1,self.nx+self.nxguard)
        else:
            ixl = 0
            ixu = max(1,self.nx+2*self.nxguard)
        if self.bc_periodic[1]:
           iyl = self.nyguard
           iyu = max(1,self.ny+self.nyguard)
        else:
            iyl = 0
            iyu = max(1,self.ny+2*self.nyguard)
        if self.bc_periodic[2]:
           izl = self.nzguard
           izu = max(1,self.nz+self.nzguard)
        else:
            izl = 0
            izu = max(1,self.nz+2*self.nzguard)
        return ixl,ixu,iyl,iyu,izl,izu

    def divsetorig(self,num,denom,orig):
        """
        Divide num by denom when denom is 0 at origin. Value at origin is provided by orig.
        """
        denommag = denom.copy()
        if len(self.dims)==1:denommag[0]=1.
        if len(self.dims)==2:denommag[0,0]=1.
        if len(self.dims)==3:denommag[0,0,0]=1.
        r = num/denommag
        if len(self.dims)==1:r[0]=orig
        if len(self.dims)==2:r[0,0]=orig
        if len(self.dims)==3:r[0,0,0]=orig
        del denommag
        return r

class GPSTD_Matrix():

    def __init__(self,fields={}):
        self.fields=fields
        n = len(self.fields.keys())

        self.mat = np.zeros([n,n]).tolist()
        for i in range(n):
            self.mat[i][i]=1.

        self.fields_name = {}
        self.fields_order = {}
        for i,fname in enumerate(self.fields.keys()):
            self.fields_name[i]=fname
            self.fields_order[fname]=i

    def add_op(self,fname,ops):
        assert fname in self.fields.keys(), "Error on GPSTD/add_op: field not found in dictionary."
        ifield = self.fields_order[fname]
        for op in ops.keys():
            assert op in self.fields.keys(), "Error on GPSTD/add_op: field not found in dictionary."
            iop = self.fields_order[op]
            self.mat[ifield][iop] = ops[op]

class GPSTD(Fourier_Space):

    __flaginputs__ = {'dt':1.,'ntsub':1}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(GPSTD.__flaginputs__,kw)
        Fourier_Space.__init__(self,kwdict=kw)

        self.fields = {}        # dict. of fields in real space
        self.Ffields = {}       # dict. of fields in Fourier space
        self.LSource = {}       # dict. of logical for sources
        self.Sfilters = {}      # dict. of filters to apply to sources before push
        self.Ffilters = {}      # dict. of filters to apply to fields after push

    def add_fields(self,f,l_source=False):
        self.fields.update(f)

        for k in f.keys():
            self.LSource[k] = l_source

        self.fields_name = {}
        self.fields_order = {}
        for i,fname in enumerate(self.fields.keys()):
            self.fields_name[i]=fname
            self.fields_order[fname]=i

    def add_Ffilter(self,Fname,Ffilter):
        self.Ffilters[Fname]=Ffilter

    def add_Sfilter(self,Sname,Sfilter):
        self.Sfilters[Sname]=Sfilter

    def get_Ffields(self):
        ixl,ixu,iyl,iyu,izl,izu = self.get_ius()
        if self.Ffields=={}:
            self.fields_shape = [ixu-ixl,iyu-iyl,izu-izl]
            for k in self.fields.keys():
                self.plan_rfftn[k] = self.create_plan_rfftn(np.asarray(self.fields_shape))
                self.Ffields[k]    =self.rfftn(np.squeeze(self.fields[k][ixl:ixu,iyl:iyu,izl:izu]),plan=self.plan_rfftn[k])
        else:
            for k in self.fields.keys():
                self.Ffields[k]=self.rfftn(np.squeeze(self.fields[k][ixl:ixu,iyl:iyu,izl:izu]),field_out=self.Ffields[k],plan=self.plan_rfftn[k])

    def get_fields(self):
        ixl,ixu,iyl,iyu,izl,izu = self.get_ius()
        if (self.plan_irfftn=={}):
            for k in self.fields.keys():
                    self.plan_irfftn[k] = self.create_plan_irfftn(np.asarray(self.fields_shape))
        for k in self.fields.keys():
            if not self.LSource[k]:
                shapek = np.asarray(np.shape(np.squeeze(self.fields[k][ixl:ixu,iyl:iyu,izl:izu])))
                f = self.irfftn(self.Ffields[k], shapek, field_out=np.squeeze(self.fields[k][ixl:ixu,iyl:iyu,izl:izu]), plan=self.plan_irfftn[k])
                f.resize(self.fields_shape)
                self.fields[k][ixl:ixu,iyl:iyu,izl:izu] = f.real

    def push_fields(self):

        self.get_Ffields()

        # --- filter sources before push
        for k in self.Sfilters.keys():
            self.Ffields[k]*=self.Sfilters[k]

        mymat = self.mymat
        n = len(mymat)

        # --- set dictionary of field values before time step
        oldfields = {}
        for k in self.Ffields.keys():
            oldfields[k] = self.Ffields[k].copy()

        # --- set dictionary of field flags for update
        updated_fields = {}
        for k in self.Ffields.keys():
            updated_fields[k] = False

        # --- fields update
        # --- multiply vector 'fields' by matrix 'mymat', returning result in 'fields'
        for i in range(n):
            # --- get key of field of rank i
            ki = self.fields_name[i]
            # --- cycle to next item if current field is a source
            if self.LSource[ki]:continue

            # --- update diagonal first
            # --- if matrix value is 1., do nothing
            if not mymat[i][i] is 1.:
                updated_fields[ki] = True
                if mymat[i][i] is 0.:
                # --- if matrix value is 0., zero out array
                    self.Ffields[i][...] = 0.
                else:
                # --- otherwise, multiply by matrix value
                    self.Ffields[ki] *= mymat[i][i]

            # --- update field for non-diagonal matrix elements.
            for j in range(n):
                if i!=j and not mymat[i][j] is 0.:
                    # --- update only if matrix element is non-zero
                    updated_fields[ki] = True
                    kj = self.fields_name[j]
                    self.Ffields[ki] += mymat[i][j]*oldfields[kj]

        del oldfields

        # --- filter fields after push
        for k in self.Ffilters.keys():
            self.Ffields[k]*=self.Ffilters[k]

        self.get_fields()

        # --- set periodic BC
        if self.bc_periodic[0]:
            ngx = self.nxguard
        else:
            ngx = 0
        if self.bc_periodic[1]:
            ngy = self.nyguard
        else:
            ngy = 0
        if self.bc_periodic[2]:
            ngz = self.nzguard
        else:
            ngz = 0

        if self.bc_periodic[0]:
            for k in self.fields.keys():
                if updated_fields[k]:
                    f = self.fields[k]
                    f[-ngx-1:,...]=f[ngx:2*ngx+1,...]
                    f[:ngx,...]=f[-2*ngx-1:-ngx-1,...]
        if self.bc_periodic[1]:
            for k in self.fields.keys():
                if updated_fields[k]:
                    f = self.fields[k]
                    f[:,-ngy-1:,:]=f[:,ngy:2*ngy+1,:]
                    f[:,:ngy,:]=f[:,-2*ngy-1:-ngy-1,:]
        if self.bc_periodic[2]:
            for k in self.fields.keys():
                if updated_fields[k]:
                    f = self.fields[k]
                    f[...,-ngz-1:]=f[...,ngz:2*ngz+1]
                    f[...,:ngz]=f[...,-2*ngz-1:-ngz-1]

        del updated_fields

    def getm(self,a='ex',b='ex'):
        o=self.fields_order
        m=self.mymat
        return m[o[a]][o[b]]

class GPSTD_Maxwell_PML(GPSTD):

    __flaginputs__ = {'syf':None,'l_pushf':False,'l_pushg':False,'clight':299792458.0}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(GPSTD_Maxwell_PML.__flaginputs__,kw)

        syf=self.syf
        nx = np.max([1,syf.nx])
        ny = np.max([1,syf.ny])
        nz = np.max([1,syf.nz])
        kw['nx']=nx
        kw['ny']=ny
        kw['nz']=nz

        GPSTD.__init__(self,kwdict=kw)

        j = 1j

        if self.l_pushf:
            self.add_fields({"exx":syf.exx, \
                             "exy":syf.exy, \
                             "exz":syf.exz, \
                             "eyx":syf.eyx, \
                             "eyy":syf.eyy, \
                             "eyz":syf.eyz, \
                             "ezx":syf.ezx, \
                             "ezy":syf.ezy, \
                             "ezz":syf.ezz})
        else:
            self.add_fields({"exy":syf.exy, \
                             "exz":syf.exz, \
                             "eyx":syf.eyx, \
                             "eyz":syf.eyz, \
                             "ezx":syf.ezx, \
                             "ezy":syf.ezy})
        if self.l_pushg:
            self.add_fields({"bxx":syf.bxx, \
                             "bxy":syf.bxy, \
                             "bxz":syf.bxz, \
                             "byx":syf.byx, \
                             "byy":syf.byy, \
                             "byz":syf.byz, \
                             "bzx":syf.bzx, \
                             "bzy":syf.bzy, \
                             "bzz":syf.bzz})
        else:
            self.add_fields({"bxy":syf.bxy, \
                             "bxz":syf.bxz, \
                             "byx":syf.byx, \
                             "byz":syf.byz, \
                             "bzx":syf.bzx, \
                             "bzy":syf.bzy})
        if self.l_pushf:
            self.add_fields({"fx":syf.fx, \
                             "fy":syf.fy, \
                             "fz":syf.fz})

        if self.l_pushg:
            self.add_fields({"gx":syf.gx, \
                             "gy":syf.gy, \
                             "gz":syf.gz})

        self.get_Ffields()

        m0 = 0.
        m1 = 1.
        dt=self.dt/self.ntsub
        cdt=dt*self.clight

        if self.nx>1:
            axm = j*dt*self.clight*self.kxm
            axp = j*dt*self.clight*self.kxp
        else:
            axm = axp = 0.

        if self.ny>1:
            aym = j*dt*self.clight*self.kym
            ayp = j*dt*self.clight*self.kyp
        else:
            aym = ayp = 0.

        if self.nz>1:
            azm = j*dt*self.clight*self.kzm
            azp = j*dt*self.clight*self.kzp
        else:
            azm = azp = 0.

        if self.nx>1:
            axp0 = 0.5/self.ntsub
            axm0 = 0.65/self.ntsub
        else:
            axm0 = axp0 = 0.

        if self.ny>1:
            ayp0 = 0.55/self.ntsub
            aym0 = 0.45/self.ntsub
        else:
            aym0 = ayp0 = 0.

        if self.nz>1:
            azp0 = 0.35/self.ntsub
            azm0 = 0.25/self.ntsub
        else:
            azm0 = azp0 = 0.

        self.mymatref = self.getmaxwellmat_pml(axp0,ayp0,azp0,axm0,aym0,azm0, \
                            0.1/self.ntsub,0.11/self.ntsub,m0,m1, \
                            0.5*self.dx,0.5*self.dy,0.5*self.dz,l_matref=1)

        matcompress = getmatcompress(self.mymatref)

        self.mymatref = exp_by_squaring_matrixlist(self.mymatref, self.ntsub, matcompress=matcompress)

        self.mymat = self.getmaxwellmat_pml(axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,\
                     self.kx_unmod,self.ky_unmod,self.kz_unmod,l_matref=0,matcompress=matcompress)

        self.mymat = exp_by_squaring_matrixlist(self.mymat, self.ntsub, matcompress=matcompress)

    def getmaxwellmat_pml(self,axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,
                      kx_unmod,ky_unmod,kz_unmod,l_matref=0,
                      matcompress=None):
        c=self.clight

        matpushb = GPSTD_Matrix(self.fields)
        if self.l_pushf:
            # --- bx
            if self.l_pushg:matpushb.add_op('bxx',{'bxx':1.,'gx':axm/2,'gy':axm/2,'gz':axm/2})
            matpushb.add_op('bxy',{'bxy':1.,'ezx':-ayp/2,'ezy':-ayp/2,'ezz':-ayp/2})
            matpushb.add_op('bxz',{'bxz':1.,'eyx': azp/2,'eyy': azp/2,'eyz': azp/2})
            # --- by
            matpushb.add_op('byx',{'byx':1.,'ezx': axp/2,'ezy': axp/2,'ezz': axp/2})
            if self.l_pushg:matpushb.add_op('byy',{'byy':1.,'gx':aym/2,'gy':aym/2,'gz':aym/2})
            matpushb.add_op('byz',{'byz':1.,'exx':-azp/2,'exy':-azp/2,'exz':-azp/2})
            # --- bz
            matpushb.add_op('bzx',{'bzx':1.,'eyx':-axp/2,'eyy':-axp/2,'eyz':-axp/2})
            matpushb.add_op('bzy',{'bzy':1.,'exx': ayp/2,'exy': ayp/2,'exz': ayp/2})
            if self.l_pushg:matpushb.add_op('bzz',{'bzz':1.,'gx':azm/2,'gy':azm/2,'gz':azm/2})
        else:
            # --- bx
            if self.l_pushg:matpushb.add_op('bxx',{'bxx':1.,'gx':axm/2,'gy':axm/2,'gz':axm/2})
            matpushb.add_op('bxy',{'bxy':1.,'ezx':-ayp/2,'ezy':-ayp/2})
            matpushb.add_op('bxz',{'bxz':1.,'eyx': azp/2,'eyz': azp/2})
            # --- by
            matpushb.add_op('byx',{'byx':1.,'ezx': axp/2,'ezy': axp/2})
            if self.l_pushg:matpushb.add_op('byy',{'byy':1.,'gx':aym/2,'gy':aym/2,'gz':aym/2})
            matpushb.add_op('byz',{'byz':1.,'exy':-azp/2,'exz':-azp/2})
            # --- bz
            matpushb.add_op('bzx',{'bzx':1.,'eyx':-axp/2,'eyz':-axp/2})
            matpushb.add_op('bzy',{'bzy':1.,'exy': ayp/2,'exz': ayp/2})
            if self.l_pushg:matpushb.add_op('bzz',{'bzz':1.,'gx':azm/2,'gy':azm/2,'gz':azm/2})

        matpushe = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            # --- ex
            if self.l_pushf:matpushe.add_op('exx',{'exx':1.,'fx':axp,'fy':axp,'fz':axp})
            matpushe.add_op('exy',{'exy':1.,'bzx': aym,'bzy': aym,'bzz': aym})
            matpushe.add_op('exz',{'exz':1.,'byx':-azm,'byy':-azm,'byz':-azm})
            # --- ey
            matpushe.add_op('eyx',{'eyx':1.,'bzx':-axm,'bzy':-axm,'bzz':-axm})
            if self.l_pushf:matpushe.add_op('eyy',{'eyy':1.,'fx':ayp,'fy':ayp,'fz':ayp})
            matpushe.add_op('eyz',{'eyz':1.,'bxx': azm,'bxy': azm,'bxz': azm})
            # --- ez
            matpushe.add_op('ezx',{'ezx':1.,'byx': axm,'byy': axm,'byz': axm})
            matpushe.add_op('ezy',{'ezy':1.,'bxx':-aym,'bxy':-aym,'bxz':-aym})
            if self.l_pushf:matpushe.add_op('ezz',{'ezz':1.,'fx':azp,'fy':azp,'fz':azp})
        else:
            # --- ex
            if self.l_pushf:matpushe.add_op('exx',{'exx':1.,'fx':axp,'fy':axp,'fz':axp})
            matpushe.add_op('exy',{'exy':1.,'bzx': aym,'bzy': aym})
            matpushe.add_op('exz',{'exz':1.,'byx':-azm,'byz':-azm})
            # --- ey
            matpushe.add_op('eyx',{'eyx':1.,'bzx':-axm,'bzy':-axm})
            if self.l_pushf:matpushe.add_op('eyy',{'eyy':1.,'fx':ayp,'fy':ayp,'fz':ayp})
            matpushe.add_op('eyz',{'eyz':1.,'bxy': azm,'bxz': azm})
            # --- ez
            matpushe.add_op('ezx',{'ezx':1.,'byx': axm,'byz': axm})
            matpushe.add_op('ezy',{'ezy':1.,'bxy':-aym,'bxz':-aym})
            if self.l_pushf:matpushe.add_op('ezz',{'ezz':1.,'fx':azp,'fy':azp,'fz':azp})

        if self.l_pushf:
            matpushf = GPSTD_Matrix(self.fields)
            matpushf.add_op('fx',{'fx':1.,'exx':axm/2,'exy':axm/2,'exz':axm/2})
            matpushf.add_op('fy',{'fy':1.,'eyx':aym/2,'eyy':aym/2,'eyz':aym/2})
            matpushf.add_op('fz',{'fz':1.,'ezx':azm/2,'ezy':azm/2,'ezz':azm/2})

        if self.l_pushg:
            matpushg = GPSTD_Matrix(self.fields)
            matpushg.add_op('gx',{'gx':1.,'bxx':axp,'bxy':axp,'bxz':axp})
            matpushg.add_op('gy',{'gy':1.,'byx':ayp,'byy':ayp,'byz':ayp})
            matpushg.add_op('gz',{'gz':1.,'bzx':azp,'bzy':azp,'bzz':azp})

        if self.l_pushf:
            mymat = multmat(matpushf.mat,matpushb.mat,matcompress=matcompress)
            if self.l_pushg:mymat = multmat(mymat,matpushg.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushe.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushf.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushb.mat,matcompress=matcompress)
        else:
            mymat = multmat(matpushb.mat,matpushe.mat,matcompress=matcompress)
            if self.l_pushg:mymat = multmat(mymat,matpushg.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushb.mat,matcompress=matcompress)

        return mymat

    def push(self):
        syf = self.syf

        self.push_fields()

        for f in self.fields.values():
            if self.nx>1:
                f[:self.nxguard/2,...]=0.
                f[-self.nxguard/2:,...]=0.
            if self.ny>1:
                f[:,:self.nyguard/2,:]=0.
                f[:,-self.nyguard/2:,:]=0.
            if self.nz>1:
                f[...,:self.nzguard/2]=0.
                f[...,-self.nzguard/2:]=0.

#      scale_em3d_split_fields(syf,top.dt,self.l_pushf)

        return

class PSATD_Maxwell_PML(GPSTD):

    __flaginputs__ = {'syf':None,'l_pushf':False,'l_pushg':False,'clight':299792458.0}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(GPSTD_Maxwell_PML.__flaginputs__,kw)

        syf=self.syf
        nx = np.max([1,syf.nx])
        ny = np.max([1,syf.ny])
        nz = np.max([1,syf.nz])
        kw['nx']=nx
        kw['ny']=ny
        kw['nz']=nz

        GPSTD.__init__(self,kwdict=kw)

        dt=self.dt
        cdt=dt*self.clight
        self.wdt = self.k*cdt
        self.coswdt=np.cos(self.wdt)
        self.sinwdt=np.sin(self.wdt)

        j = 1j

        if self.l_pushf:
            self.add_fields({"exx":syf.exx, \
                             "exy":syf.exy, \
                             "exz":syf.exz, \
                             "eyx":syf.eyx, \
                             "eyy":syf.eyy, \
                             "eyz":syf.eyz, \
                             "ezx":syf.ezx, \
                             "ezy":syf.ezy, \
                             "ezz":syf.ezz})
        else:
            self.add_fields({"exy":syf.exy, \
                             "exz":syf.exz, \
                             "eyx":syf.eyx, \
                             "eyz":syf.eyz, \
                             "ezx":syf.ezx, \
                             "ezy":syf.ezy})
        if self.l_pushg:
            self.add_fields({"bxx":syf.bxx, \
                             "bxy":syf.bxy, \
                             "bxz":syf.bxz, \
                             "byx":syf.byx, \
                             "byy":syf.byy, \
                             "byz":syf.byz, \
                             "bzx":syf.bzx, \
                             "bzy":syf.bzy, \
                             "bzz":syf.bzz})
        else:
            self.add_fields({"bxy":syf.bxy, \
                             "bxz":syf.bxz, \
                             "byx":syf.byx, \
                             "byz":syf.byz, \
                             "bzx":syf.bzx, \
                             "bzy":syf.bzy})
        if self.l_pushf:
            self.add_fields({"fx":syf.fx, \
                             "fy":syf.fy, \
                             "fz":syf.fz})

        if self.l_pushg:
            self.add_fields({"gx":syf.gx, \
                             "gy":syf.gy, \
                             "gz":syf.gz})

        self.get_Ffields()

        m0 = 0.
        m1 = 1.
        dt=self.dt
        cdt=dt*self.clight
        C=self.coswdt
        S=self.sinwdt

        if self.nx>1:
            axm = j*S*self.kxmn
            axp = j*S*self.kxpn
        else:
            axm = axp = 0.

        if self.ny>1:
            aym = j*S*self.kymn
            ayp = j*S*self.kypn
        else:
            aym = ayp = 0.

        if self.nz>1:
            azm = j*S*self.kzmn
            azp = j*S*self.kzpn
        else:
            azm = azp = 0.

        self.mymat = self.getmaxwellmat_pml(C,S,axp,ayp,azp,axm,aym,azm)

    def getmaxwellmat_pml(self,C,S,axp,ayp,azp,axm,aym,azm):
        mymat = GPSTD_Matrix(self.fields)
        if self.l_pushf:
            # --- bx
            if self.l_pushg:mymat.add_op('bxx',{'bxx':C,'gx':axm,'gy':axm,'gz':axm})
            mymat.add_op('bxy',{'bxy':C,'ezx':-ayp,'ezy':-ayp,'ezz':-ayp})
            mymat.add_op('bxz',{'bxz':C,'eyx': azp,'eyy': azp,'eyz': azp})
            # --- by
            mymat.add_op('byx',{'byx':C,'ezx': axp,'ezy': axp,'ezz': axp})
            if self.l_pushg:mymat.add_op('byy',{'byy':C,'gx':aym,'gy':aym,'gz':aym})
            mymat.add_op('byz',{'byz':C,'exx':-azp,'exy':-azp,'exz':-azp})
            # --- bz
            mymat.add_op('bzx',{'bzx':C,'eyx':-axp,'eyy':-axp,'eyz':-axp})
            mymat.add_op('bzy',{'bzy':C,'exx': ayp,'exy': ayp,'exz': ayp})
            if self.l_pushg:mymat.add_op('bzz',{'bzz':C,'gx':azm,'gy':azm,'gz':azm})
        else:
            # --- bx
            if self.l_pushg:mymat.add_op('bxx',{'bxx':C,'gx':axm,'gy':axm,'gz':axm})
            mymat.add_op('bxy',{'bxy':C,'ezx':-ayp,'ezy':-ayp})
            mymat.add_op('bxz',{'bxz':C,'eyx': azp,'eyz': azp})
            # --- by
            mymat.add_op('byx',{'byx':C,'ezx': axp,'ezy': axp})
            if self.l_pushg:mymat.add_op('byy',{'byy':C,'gx':aym,'gy':aym,'gz':aym})
            mymat.add_op('byz',{'byz':C,'exy':-azp,'exz':-azp})
            # --- bz
            mymat.add_op('bzx',{'bzx':C,'eyx':-axp,'eyz':-axp})
            mymat.add_op('bzy',{'bzy':C,'exy': ayp,'exz': ayp})
            if self.l_pushg:mymat.add_op('bzz',{'bzz':C,'gx':azm,'gy':azm,'gz':azm})

        if self.l_pushg:
            # --- ex
            if self.l_pushf:mymat.add_op('exx',{'exx':C,'fx':axp,'fy':axp,'fz':axp})
            mymat.add_op('exy',{'exy':C,'bzx': aym,'bzy': aym,'bzz': aym})
            mymat.add_op('exz',{'exz':C,'byx':-azm,'byy':-azm,'byz':-azm})
            # --- ey
            mymat.add_op('eyx',{'eyx':C,'bzx':-axm,'bzy':-axm,'bzz':-axm})
            if self.l_pushf:mymat.add_op('eyy',{'eyy':C,'fx':ayp,'fy':ayp,'fz':ayp})
            mymat.add_op('eyz',{'eyz':C,'bxx': azm,'bxy': azm,'bxz': azm})
            # --- ez
            mymat.add_op('ezx',{'ezx':C,'byx': axm,'byy': axm,'byz': axm})
            mymat.add_op('ezy',{'ezy':C,'bxx':-aym,'bxy':-aym,'bxz':-aym})
            if self.l_pushf:mymat.add_op('ezz',{'ezz':C,'fx':azp,'fy':azp,'fz':azp})
        else:
            # --- ex
            if self.l_pushf:mymat.add_op('exx',{'exx':C,'fx':axp,'fy':axp,'fz':axp})
            mymat.add_op('exy',{'exy':C,'bzx': aym,'bzy': aym})
            mymat.add_op('exz',{'exz':C,'byx':-azm,'byz':-azm})
            # --- ey
            mymat.add_op('eyx',{'eyx':C,'bzx':-axm,'bzy':-axm})
            if self.l_pushf:mymat.add_op('eyy',{'eyy':C,'fx':ayp,'fy':ayp,'fz':ayp})
            mymat.add_op('eyz',{'eyz':C,'bxy': azm,'bxz': azm})
            # --- ez
            mymat.add_op('ezx',{'ezx':C,'byx': axm,'byz': axm})
            mymat.add_op('ezy',{'ezy':C,'bxy':-aym,'bxz':-aym})
            if self.l_pushf:mymat.add_op('ezz',{'ezz':C,'fx':azp,'fy':azp,'fz':azp})

        if self.l_pushf:
            mymat.add_op('fx',{'fx':C,'exx':axm,'exy':axm,'exz':axm})
            mymat.add_op('fy',{'fy':C,'eyx':aym,'eyy':aym,'eyz':aym})
            mymat.add_op('fz',{'fz':C,'ezx':azm,'ezy':azm,'ezz':azm})

        if self.l_pushg:
            mymat.add_op('gx',{'gx':C,'bxx':axp,'bxy':axp,'bxz':axp})
            mymat.add_op('gy',{'gy':C,'byx':ayp,'byy':ayp,'byz':ayp})
            mymat.add_op('gz',{'gz':C,'bzx':azp,'bzy':azp,'bzz':azp})

        return mymat.mat

    def push(self):
        syf = self.syf

        self.push_fields()

        for f in self.fields.values():
            if self.nx>1:
                f[:self.nxguard,...]=0.
                f[-self.nxguard/2:,...]=0.
            if self.ny>1:
                f[:,:self.nyguard/2,:]=0.
                f[:,-self.nyguard/2:,:]=0.
            if self.nz>1:
                f[...,:self.nzguard/2]=0.
                f[...,-self.nzguard/2:]=0.

#      scale_em3d_split_fields(syf,top.dt,self.l_pushf)

        return

class GPSTD_Maxwell(GPSTD):

    __flaginputs__ = {'yf':None,
                      'l_pushf':False,
                      'l_pushg':False,
                      'clight':299792458.0,
                      'eps0':8.854187817620389e-12,
                      'V_galilean':np.array([0.,0.,0.]),
                      'V_pseudogalilean':np.array([0.,0.,0.]),
                      'l_theta_with_modif_k':True}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(GPSTD_Maxwell.__flaginputs__,kw)

        yf=self.yf
        nx = np.max([1,yf.nx])
        ny = np.max([1,yf.ny])
        nz = np.max([1,yf.nz])
        kw['nx']=nx
        kw['ny']=ny
        kw['nz']=nz

        GPSTD.__init__(self,kwdict=kw)

        j = 1j

        self.add_fields({"bx":yf.Bx,"by":yf.By,"bz":yf.Bz, \
                         "ex":yf.Ex,"ey":yf.Ey,"ez":yf.Ez})
        self.add_fields({"rho":yf.Rho},True)

        if self.l_pushf:
            self.add_fields({"f":yf.F})
            self.add_fields({"rhoold":yf.Rhoold},True)
            self.add_fields({"rhonew":yf.Rho},True)
            self.add_fields({"drho":yf.Rho},True)
        if self.l_pushg:
            self.add_fields({"g":yf.G})
#        if self.l_pushf or self.l_getrho:
#            self.add_fields({"rho":yf.Rho},True)
        self.add_fields({"rho":yf.Rho},True)
        self.add_fields({"drho":yf.Rhoold},True)
        self.add_fields({"jx":yf.Jx,"jy":yf.Jy,"jz":yf.Jz},True)

        self.get_Ffields()

        m0 = 0.
        m1 = 1.
        dt=self.dt/self.ntsub
        cdt=dt*self.clight

        if self.nx>1:
            axm = j*dt*self.clight*self.kxm
            axp = j*dt*self.clight*self.kxp
        else:
            axm = axp = 0.

        if self.ny>1:
            aym = j*dt*self.clight*self.kym
            ayp = j*dt*self.clight*self.kyp
        else:
            aym = ayp = 0.

        if self.nz>1:
            azm = j*dt*self.clight*self.kzm
            azp = j*dt*self.clight*self.kzp
        else:
            azm = azp = 0.

        if self.nx>1:
            axp0 = 0.5/self.ntsub
            axm0 = 0.65/self.ntsub
        else:
            axm0 = axp0 = 0.

        if self.ny>1:
            ayp0 = 0.55/self.ntsub
            aym0 = 0.45/self.ntsub
        else:
            aym0 = ayp0 = 0.

        if self.nz>1:
            azp0 = 0.35/self.ntsub
            azm0 = 0.25/self.ntsub
        else:
            azm0 = azp0 = 0.

        self.mymatref = self.getmaxwellmat(axp0,ayp0,azp0,axm0,aym0,azm0, \
                            0.1/self.ntsub,0.11/self.ntsub,m0,m1, \
                            0.5*self.dx,0.5*self.dy,0.5*self.dz,l_matref=1)

        matcompress = getmatcompress(self.mymatref)

        self.mymatref = exp_by_squaring_matrixlist(self.mymatref, self.ntsub, matcompress=matcompress)

        if np.all(self.V_galilean==0.):
            self.mymat = self.getmaxwellmat(axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,\
                         self.kx_unmod,self.ky_unmod,self.kz_unmod,l_matref=0,matcompress=matcompress)
        else:
            self.mymat = self.getmaxwellmat_galilean(axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,\
                         self.kx_unmod,self.ky_unmod,self.kz_unmod,l_matref=0,matcompress=matcompress,
                         V_galilean=self.V_galilean)

    def getmaxwellmat(self,axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,
                      kx_unmod,ky_unmod,kz_unmod,l_matref=0,
                      matcompress=None):
        c=self.clight

        if self.l_pushf:
            matpushrho = GPSTD_Matrix(self.fields)
#            matpushrho.add_op('rho',{'rho':1.,'jx':-axm/c,'jy':-aym/c,'jz':-azm/c})
            matpushrho.add_op('rho',{'rho':1.,'drho':1./self.ntsub})

        matpushb = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            matpushb.add_op('bx',{'bx':1.,'ey': azp/(2*c),'ez':-ayp/(2*c),'g':axm/(2*c)})
            matpushb.add_op('by',{'by':1.,'ex':-azp/(2*c),'ez': axp/(2*c),'g':aym/(2*c)})
            matpushb.add_op('bz',{'bz':1.,'ex': ayp/(2*c),'ey':-axp/(2*c),'g':azm/(2*c)})
        else:
            matpushb.add_op('bx',{'bx':1.,'ey': azp/(2*c),'ez':-ayp/(2*c)})
            matpushb.add_op('by',{'by':1.,'ex':-azp/(2*c),'ez': axp/(2*c)})
            matpushb.add_op('bz',{'bz':1.,'ex': ayp/(2*c),'ey':-axp/(2*c)})

        matpushe = GPSTD_Matrix(self.fields)
        if self.l_pushf:
            matpushe.add_op('ex',{'ex':1.,'by':-azm*c,'bz': aym*c,'f':axp,'jx':-dt/self.eps0})
            matpushe.add_op('ey',{'ey':1.,'bx': azm*c,'bz':-axm*c,'f':ayp,'jy':-dt/self.eps0})
            matpushe.add_op('ez',{'ez':1.,'bx':-aym*c,'by': axm*c,'f':azp,'jz':-dt/self.eps0})
        else:
            matpushe.add_op('ex',{'ex':1.,'by':-azm*c,'bz': aym*c,'jx':-dt/self.eps0})
            matpushe.add_op('ey',{'ey':1.,'bx': azm*c,'bz':-axm*c,'jy':-dt/self.eps0})
            matpushe.add_op('ez',{'ez':1.,'bx':-aym*c,'by': axm*c,'jz':-dt/self.eps0})

        if self.l_pushf:
            matpushf = GPSTD_Matrix(self.fields)
            matpushf.add_op('f',{'f':1.,'ex':axm/2,'ey':aym/2,'ez':azm/2,'rho':-0.5*cdt/self.eps0})

        if self.l_pushg:
            matpushg = GPSTD_Matrix(self.fields)
            matpushg.add_op('g',{'g':1.,'bx':axp*c,'by':ayp*c,'bz':azp*c})

        if self.l_pushf:
            mymat_init = GPSTD_Matrix(self.fields)
            mymat_init.add_op('rho',{'rho':0.,'rhoold':1.})
            mymat_init.add_op('drho',{'drho':0.,'rhonew':1.,'rhoold':-1.})

        if self.l_pushf:
            mymat = multmat(matpushf.mat,matpushb.mat,matcompress=matcompress)
            if self.l_pushg:mymat = multmat(mymat,matpushg.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushe.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushrho.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushf.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushb.mat,matcompress=matcompress)
        else:
            mymat = multmat(matpushb.mat,matpushe.mat,matcompress=matcompress)
            if self.l_pushg:mymat = multmat(mymat,matpushg.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushb.mat,matcompress=matcompress)

        self.mymat = exp_by_squaring_matrixlist(mymat, self.ntsub, matcompress=matcompress)
        if self.l_pushf:
            self.mymat = multmat(mymat_init.mat,self.mymat)

        if l_matref:
            self.matpushb=matpushb
            self.matpushe=matpushe
            if self.l_pushf:self.matpushf=matpushf
            if self.l_pushg:self.matpushg=matpushg
            if self.l_pushf:self.matpushrho=matpushrho

        return self.mymat

    def getmaxwellmat_galilean(self,axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,
                      kx_unmod,ky_unmod,kz_unmod,l_matref=0,
                      matcompress=None,V_galilean=[0.,0.,0.]):

        # --- equivalent to J constant in lab frame grid when 1
        # --- equivalent to J constant in moving frame grid when 0
        l_matpushj = 0

        j=1j
        c=self.clight
        V0  = np.sqrt(V_galilean[0]*V_galilean[0]+V_galilean[1]*V_galilean[1]+V_galilean[2]*V_galilean[2])
        w=self.k*c
        kV=self.kx_unmod*V_galilean[0]+self.ky_unmod*V_galilean[1]+self.kz_unmod*V_galilean[2]
        T=np.exp(j*kV*dt)
        T2=np.exp(j*kV*dt/2)

        Theta=np.exp(j*kV*self.dt)
        coef = -j*self.divsetorig(j*kV,1.-Theta,-1./self.dt)

        is_singular=(w!=0) & (kV==0)
        coef[is_singular] = j/self.dt
        self.cr=coef.copy()
        self.cro=coef.copy()*Theta
        coef /= self.kmag
        self.JxCorRhomult = coef*self.kxpn
        self.JyCorRhomult = coef*self.kypn
        self.JzCorRhomult = coef*self.kzpn

        self.JxCorRhooldmult = coef*self.kxpn*Theta
        self.JyCorRhooldmult = coef*self.kypn*Theta
        self.JzCorRhooldmult = coef*self.kzpn*Theta

        axp = axp
        ayp = ayp
        azp = azp
        axm = axm
        aym = aym
        azm = azm

        if self.l_pushf:
            matpushrho = GPSTD_Matrix(self.fields)
#            matpushrho.add_op('rho',{'rho':T,'jx':-axm/c,'jy':-aym/c,'jz':-azm/c})
#            matpushrho.add_op('rho',{'rho':T,'drho':1./self.ntsub})
            alpha = 0.5*j*kV*dt
            matpushrho.add_op('rho',{'rho':(1.+alpha)/(1.-alpha),'drho':1./(self.ntsub*(1.-alpha))})
            matpushdrho = GPSTD_Matrix(self.fields)
            matpushdrho.add_op('drho',{'drho':T})

        matpushb = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            matpushb.add_op('bx',{'bx':T2,'ey': azp*T2/(2*c),'ez':-ayp*T2/(2*c),'g':axm*T2/(2*c)})
            matpushb.add_op('by',{'by':T2,'ex':-azp*T2/(2*c),'ez': axp*T2/(2*c),'g':aym*T2/(2*c)})
            matpushb.add_op('bz',{'bz':T2,'ex': ayp*T2/(2*c),'ey':-axp*T2/(2*c),'g':azm*T2/(2*c)})
        else:
            matpushb.add_op('bx',{'bx':T2,'ey': azp*T2/(2*c),'ez':-ayp*T2/(2*c)})
            matpushb.add_op('by',{'by':T2,'ex':-azp*T2/(2*c),'ez': axp*T2/(2*c)})
            matpushb.add_op('bz',{'bz':T2,'ex': ayp*T2/(2*c),'ey':-axp*T2/(2*c)})

        matpushe = GPSTD_Matrix(self.fields)
        if self.l_pushf:
            matpushe.add_op('ex',{'ex':T,'by':-azm*T2*c,'bz': aym*T2*c,'f':axp*T2,'jx':-dt*T/self.eps0})
            matpushe.add_op('ey',{'ey':T,'bx': azm*T2*c,'bz':-axm*T2*c,'f':ayp*T2,'jy':-dt*T/self.eps0})
            matpushe.add_op('ez',{'ez':T,'bx':-aym*T2*c,'by': axm*T2*c,'f':azp*T2,'jz':-dt*T/self.eps0})
        else:
            matpushe.add_op('ex',{'ex':T,'by':-azm*T2*c,'bz': aym*T2*c,'jx':-dt*T/self.eps0})
            matpushe.add_op('ey',{'ey':T,'bx': azm*T2*c,'bz':-axm*T2*c,'jy':-dt*T/self.eps0})
            matpushe.add_op('ez',{'ez':T,'bx':-aym*T2*c,'by': axm*T2*c,'jz':-dt*T/self.eps0})

        if self.l_pushf:
            matpushf = GPSTD_Matrix(self.fields)
            matpushf.add_op('f',{'f':T2,'ex':axm*T2/2,'ey':aym*T2/2,'ez':azm*T2/2,'rho':-0.5*cdt/self.eps0})

        if self.l_pushg:
            matpushg = GPSTD_Matrix(self.fields)
            matpushg.add_op('g',{'g':T,'bx':axp*T2*c,'by':ayp*T2*c,'bz':azp*T2*c})

        mymat_init = GPSTD_Matrix(self.fields)
        if l_matpushj:
            shift_init = np.exp((-0.5+0*0.5/self.ntsub)*j*kV*self.dt)
        else:
            shift_init = 1.
        mymat_init.add_op('jx',{'jx':shift_init})
        mymat_init.add_op('jy',{'jy':shift_init})
        mymat_init.add_op('jz',{'jz':shift_init})
        if self.l_pushf:
            mymat_init.add_op('rho',{'rhoold':1.})
            if l_matpushj:
                mymat_init.add_op('drho',{'rhonew':np.exp(-j*kV*self.dt),'rhoold':-1.})
            else:
                mymat_init.add_op('drho',{'rhonew':np.exp(-j*kV*self.dt),'rhoold':-1.})

        matpushj = GPSTD_Matrix(self.fields)
        matpushj.add_op('jx',{'jx':T})
        matpushj.add_op('jy',{'jy':T})
        matpushj.add_op('jz',{'jz':T})

        if self.l_pushf:
            mymat = multmat(matpushf.mat,matpushb.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushe.mat,matcompress=matcompress)
        else:
            mymat = multmat(matpushb.mat,matpushe.mat,matcompress=matcompress)
        if self.l_pushg:mymat = multmat(mymat,matpushg.mat,matcompress=matcompress)
        if self.l_pushf:mymat = multmat(mymat,matpushdrho.mat,matcompress=matcompress)
        if self.l_pushf:mymat = multmat(mymat,matpushrho.mat,matcompress=matcompress)
        if self.l_pushf:mymat = multmat(mymat,matpushf.mat,matcompress=matcompress)
        mymat = multmat(mymat,matpushb.mat,matcompress=matcompress)

        if l_matpushj:mymat = multmat(mymat,matpushj.mat,matcompress=matcompress)

        mymat = exp_by_squaring_matrixlist(mymat, self.ntsub, matcompress=matcompress)

        self.mymat = multmat(mymat_init.mat,mymat)

        if l_matref:
            self.matpushb=matpushb
            self.matpushe=matpushe
            if self.l_pushf:self.matpushf=matpushf
            if self.l_pushg:self.matpushg=matpushg
            if self.l_pushf:self.matpushrho=matpushrho

        return self.mymat

    def push(self):

        self.push_fields()

        return

class PSATD_Maxwell(GPSTD):

    __flaginputs__ = {'yf':None,
                      'l_pushf':False,
                      'l_pushg':False,
                      'clight':299792458.0,
                      'eps0':8.854187817620389e-12,
                      'V_galilean':np.array([0.,0.,0.]),
                      'V_pseudogalilean':np.array([0.,0.,0.]),
                      'l_theta_with_modif_k':True}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(PSATD_Maxwell.__flaginputs__,kw)

        yf=self.yf

        nx = np.max([1,yf.nx])
        ny = np.max([1,yf.ny])
        nz = np.max([1,yf.nz])
        kw['nx']=nx
        kw['ny']=ny
        kw['nz']=nz

        GPSTD.__init__(self,kwdict=kw)

        j = 1j

        self.add_fields({"bx":yf.Bx,"by":yf.By,"bz":yf.Bz, \
                         "ex":yf.Ex,"ey":yf.Ey,"ez":yf.Ez})
        if self.l_pushf:
            self.add_fields({"f":yf.F})
        if self.l_pushg:
            self.add_fields({"g":yf.G})
        self.add_fields({"rhoold":yf.Rhoold},True)
        self.add_fields({"rhonew":yf.Rho},True)
        self.add_fields({"jx":yf.Jx,"jy":yf.Jy,"jz":yf.Jz},True)

        self.get_Ffields()

        m0 = 0.
        m1 = 1.
        dt=self.dt
        cdt=dt*self.clight
        self.wdt = self.k*cdt
        self.coswdt=np.cos(self.wdt)
        self.sinwdt=np.sin(self.wdt)
        C=self.coswdt
        S=self.sinwdt

        if np.all(self.V_galilean==0.) and np.all(self.V_pseudogalilean==0.):
            self.mymat = self.getmaxwellmat(self.kxpn,self.kypn,self.kzpn,\
                         self.kxmn,self.kymn,self.kzmn,dt,cdt)
        else:
            if np.any(self.V_galilean != 0.):
                self.mymat = self.getmaxwellmat_galilean(self.kxpn,self.kypn,self.kzpn,\
                             self.kxmn,self.kymn,self.kzmn,dt,cdt,self.V_galilean,
                             l_theta_with_modif_k=self.l_theta_with_modif_k)
            if np.any(self.V_pseudogalilean != 0.):
                self.mymat = self.getmaxwellmat_pseudogalilean(self.kxpn,self.kypn,self.kzpn,\
                             self.kxmn,self.kymn,self.kzmn,dt,cdt,self.V_pseudogalilean)

    def getmaxwellmat(self,kxpn,kypn,kzpn,kxmn,kymn,kzmn,dt,cdt):

        j = 1j
        c=self.clight
        C=self.coswdt
        S=self.sinwdt

        Soverk = self.divsetorig(S,self.kmag,self.dt*self.clight)
        Jmult = 1./(self.kmag*self.clight*self.eps0)

        EJmult = -self.divsetorig(S,self.kmag*self.clight*self.eps0,self.dt/self.eps0)

        ERhomult = j*(-EJmult/dt-1./self.eps0)/self.kmag
        ERhooldmult = j*(C/self.eps0+EJmult/dt) /self.kmag

        BJmult = j*(C-1.)*Jmult/self.clight

        FJmult = j*(C-1.)*Jmult
        FRhomult = self.divsetorig(C-1.,dt*self.kmag**2*self.clight*self.eps0,-0.5*self.dt*self.clight/self.eps0)

        if self.nx>1:
            axm = j*S*self.kxmn
            axp = j*S*self.kxpn
            kxpn = self.kxpn
            kxmn = self.kxmn
        else:
            axm = axp = 0.
            bxm = bxp = 0.
            kxpn = kxmn = 0.

        if self.ny>1:
            aym = j*S*self.kymn
            ayp = j*S*self.kypn
            kypn = self.kypn
            kymn = self.kymn
        else:
            aym = ayp = 0.
            bym = byp = 0.
            kypn = kymn = 0.

        if self.nz>1:
            azm = j*S*self.kzmn
            azp = j*S*self.kzpn
            kzpn = self.kzpn
            kzmn = self.kzmn
        else:
            azm = azp = 0.
            bzm = bzp = 0.
            kzpn = kzmn = 0.

        self.BJmult = BJmult
        self.EJmult = EJmult
        self.ERhomult = ERhomult
        self.ERhooldmult = ERhooldmult
        self.Jmult = Jmult
        self.Soverk = Soverk

        mymat = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'g':axm/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'g':aym/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'g':azm/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})
        else:
            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})

        if self.l_pushf:
            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'f':axp,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'f':ayp,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'f':azp,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})
        else:
            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})

        if self.l_pushf:
            mymat.add_op('f',{'f':C,'ex':axm,'ey':aym,'ez':azm, \
                                    'jx': kxmn*FJmult,'jy': kymn*FJmult,'jz': kzmn*FJmult, \
                                    'rhonew':FRhomult,\
                                    'rhoold':-FRhomult - Soverk/self.eps0})

        if self.l_pushg:
            mymat.add_op('g',{'g':C,'bx':axp*c,'by':ayp*c,'bz':azp*c})

        return mymat.mat

    def getmaxwellmat_galilean(self,kxpn,kypn,kzpn,kxmn,kymn,kzmn,dt,cdt,
                V_galilean=np.array([0.,0.,0.]), l_theta_with_modif_k=True):

        j = 1j
        V0 = np.linalg.norm(V_galilean)
        c=self.clight
        C=self.coswdt
        S=self.sinwdt
        if l_theta_with_modif_k:
            kV = self.kxp*V_galilean[0] + self.kyp*V_galilean[1] + self.kzp*V_galilean[2]
        else:
            kV=self.kx_unmod*V_galilean[0]+self.ky_unmod*V_galilean[1]+self.kz_unmod*V_galilean[2]
        Theta=T=np.exp(j*kV*dt)
        CT = C*Theta
        ST = S*Theta
        w = self.k*c
        kVow = self.divsetorig(kV,self.kmag*c,V0/c)
        So1mT = self.divsetorig(S,1.-T,j*c/V0)
        onemCo1mT = self.divsetorig(1.-C,1.-T,0.)

        denom = (w*w-kV*kV)
        self.denom=denom

        X1 = self.divsetorig(1.-CT+j*kVow*ST, denom, dt**2*0.5)
        X2 = self.divsetorig(1.+j*kVow*T*So1mT+kVow**2*T*onemCo1mT, denom, dt**2/6)
        X3 = T*self.divsetorig(C+j*kVow*T*So1mT+kVow**2*onemCo1mT, denom, -dt**2/3)

        #Apply a special limit when kV=0 but w!=0
        is_singular=(w!=0) & (kV==0)
        X2[is_singular]=(1.-S[is_singular]/(w[is_singular]*dt))/w[is_singular]**2
        X3[is_singular] = T[is_singular]*(C[is_singular]-S[is_singular]/(w[is_singular]*dt))/w[is_singular]**2

        Soverk = self.divsetorig(S,self.kmag,self.dt*self.clight)
        Jmult = 1./(self.kmag*self.clight*self.eps0)

        EJmult = -self.divsetorig(ST,self.kmag*self.clight*self.eps0,dt/self.eps0)+j*X1*kV/self.eps0

        ERhomult = -j*c**2*X2*self.k/self.eps0
        ERhooldmult = j*c**2*X3*self.k/self.eps0

        BJmult = -self.k*j*X1/self.eps0

        FJmult = j*(C-1.)*Jmult
        FRhomult = (C-1.)/(dt*self.kmag**2*self.clight*self.eps0)

        coef = -j*self.divsetorig(j*kV,1.-T,-1./dt)
        coef[is_singular]=j/dt
        self.CDcoef = coef.copy()*np.exp(0.5*j*kV*dt)

        self.cr=coef.copy()
        self.cro=coef.copy()*T
        coef /= self.kmag
        self.JxCorRhomult = coef*self.kxpn
        self.JyCorRhomult = coef*self.kypn
        self.JzCorRhomult = coef*self.kzpn

        self.JxCorRhooldmult = coef*self.kxpn*T
        self.JyCorRhooldmult = coef*self.kypn*T
        self.JzCorRhooldmult = coef*self.kzpn*T

        if len(self.dims)==1:
            FRhomult[0] = -0.5*self.dt*self.clight/self.eps0
        if len(self.dims)==2:
            FRhomult[0,0] = -0.5*self.dt*self.clight/self.eps0
        if len(self.dims)==3:
            FRhomult[0,0,0] = -0.5*self.dt*self.clight/self.eps0

        if self.nx>1:
            axm = j*ST*self.kxmn
            axp = j*ST*self.kxpn
            kxpn = self.kxpn
            kxmn = self.kxmn
        else:
            axm = axp = 0.
            bxm = bxp = 0.
            kxpn = kxmn = 0.

        if self.ny>1:
            aym = j*ST*self.kymn
            ayp = j*ST*self.kypn
            kypn = self.kypn
            kymn = self.kymn
        else:
            aym = ayp = 0.
            bym = byp = 0.
            kypn = kymn = 0.

        if self.nz>1:
            azm = j*ST*self.kzmn
            azp = j*ST*self.kzpn
            kzpn = self.kzpn
            kzmn = self.kzmn
        else:
            azm = azp = 0.
            bzm = bzp = 0.
            kzpn = kzmn = 0.

        self.BJmult = BJmult
        self.EJmult = EJmult
        self.ERhomult = ERhomult
        self.ERhooldmult = ERhooldmult
        self.Jmult = Jmult
        self.Soverk = Soverk
        self.X1=X1
        self.X2=X2
        self.X3=X3
        self.kVow = kVow
        self.So1mT = So1mT
        self.onemCo1mT = onemCo1mT
        self.kV=kV
        self.T=T
        self.CT = CT
        self.ST = ST

        mymat = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            mymat.add_op('bx',{'bx':CT,'ey': azp/c,'ez':-ayp/c,'g':axm/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':CT,'ex':-azp/c,'ez': axp/c,'g':aym/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':CT,'ex': ayp/c,'ey':-axp/c,'g':azm/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})
        else:
            mymat.add_op('bx',{'bx':CT,'ey': azp/c,'ez':-ayp/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':CT,'ex':-azp/c,'ez': axp/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':CT,'ex': ayp/c,'ey':-axp/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})

        if self.l_pushf:
            mymat.add_op('ex',{'ex':CT,'by':-azm*c,'bz': aym*c,'jx':EJmult,'f':axp,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':CT,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'f':ayp,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':CT,'bx':-aym*c,'by': axm*c,'jz':EJmult,'f':azp,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})
        else:
            mymat.add_op('ex',{'ex':CT,'by':-azm*c,'bz': aym*c,'jx':EJmult,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':CT,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':CT,'bx':-aym*c,'by': axm*c,'jz':EJmult,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})

        if self.l_pushf:
            print 'l_pushf not yet implemented in PSATD Galilean'
            raise
            mymat.add_op('f',{'f':CT,'ex':axm,'ey':aym,'ez':azm, \
                                    'jx': kxmn*FJmult,'jy': kymn*FJmult,'jz': kzmn*FJmult, \
                                    'rhonew':FRhomult,\
                                    'rhoold':-FRhomult - Soverk/self.eps0})

        if self.l_pushg:
            print 'l_pushg not yet implemented in PSATD Galilean'
            raise
            mymat.add_op('g',{'g':CT,'bx':axp*c,'by':ayp*c,'bz':azp*c})

        return mymat.mat

    def getmaxwellmat_pseudogalilean(self,kxpn,kypn,kzpn,kxmn,kymn,kzmn,dt,cdt,V_galilean=np.array([0.,0.,0.])):

        j = 1j
        V0 = np.sqrt(V_galilean[0]*V_galilean[0]+V_galilean[1]*V_galilean[1]+V_galilean[2]*V_galilean[2])
        c=self.clight
        C=self.coswdt
        S=self.sinwdt
        kV=self.kx_unmod*V_galilean[0]+self.ky_unmod*V_galilean[1]+self.kz_unmod*V_galilean[2]
        Theta=T=np.exp(j*kV*dt)
        invT = np.exp(-j*kV*dt)
        CT = C*Theta
        ST = S*Theta
        w = self.k*c
        kVow = self.divsetorig(kV,self.kmag*c,V0/c)
        So1mT = self.divsetorig(S,1.-T,j*c/V0)
        onemCo1mT = self.divsetorig(1.-C,1.-T,0.)

        denom = (w*w-kV*kV)
        self.denom=denom

        X1 = self.divsetorig(1.-CT+j*kVow*ST, denom, dt**2*0.5)
        X2 = self.divsetorig(1.+j*kVow*T*So1mT+kVow**2*T*onemCo1mT, denom, dt**2/6)
        X3 = self.divsetorig(C+j*kVow*T*So1mT+kVow**2*onemCo1mT, denom, -dt**2/3)

        #Apply a special limit when kV=0 but w!=0
        is_singular=(w!=0.) & (kV==0.)
        X2[is_singular]=(1.-S[is_singular]/(w[is_singular]*dt))/w[is_singular]**2
        X3[is_singular] = (C[is_singular]-S[is_singular]/(w[is_singular]*dt))/w[is_singular]**2

        X1 *= np.exp(-0.5*j*kV*dt)

        Soverk = self.divsetorig(S,self.kmag,self.dt*self.clight)
        Jmult = 1./(self.kmag*self.clight*self.eps0)

        EJmult = -np.exp(-0.5*j*kV*dt)*self.divsetorig(ST,self.kmag*self.clight*self.eps0,dt/self.eps0)+j*X1*kV/self.eps0

        ERhomult = -j*c**2*X2*self.k/self.eps0
        ERhooldmult = j*c**2*X3*self.k/self.eps0

        BJmult = -self.k*j*X1/self.eps0

        FJmult = j*(C-1.)*Jmult
        FRhomult = (C-1.)/(dt*self.kmag**2*self.clight*self.eps0)

        coef = self.divsetorig(j*kV*dt,T-1.,1.)
        coef[is_singular]=1.

        self.CDcoef = coef.copy()*np.exp(0.5*j*kV*dt)
        coef *= j/(dt*self.kmag)

        self.JxCorRhomult = coef*self.kxpn*np.exp(0.5*j*kV*dt)
        self.JyCorRhomult = coef*self.kypn*np.exp(0.5*j*kV*dt)
        self.JzCorRhomult = coef*self.kzpn*np.exp(0.5*j*kV*dt)

        self.JxCorRhooldmult = coef*self.kxpn*np.exp(0.5*j*kV*dt)
        self.JyCorRhooldmult = coef*self.kypn*np.exp(0.5*j*kV*dt)
        self.JzCorRhooldmult = coef*self.kzpn*np.exp(0.5*j*kV*dt)

        if len(self.dims)==1:
            FRhomult[0] = -0.5*self.dt*self.clight/self.eps0
        if len(self.dims)==2:
            FRhomult[0,0] = -0.5*self.dt*self.clight/self.eps0
        if len(self.dims)==3:
            FRhomult[0,0,0] = -0.5*self.dt*self.clight/self.eps0

        if self.nx>1:
            axm = j*S*self.kxmn
            axp = j*S*self.kxpn
            kxpn = self.kxpn
            kxmn = self.kxmn
        else:
            axm = axp = 0.
            bxm = bxp = 0.
            kxpn = kxmn = 0.

        if self.ny>1:
            aym = j*S*self.kymn
            ayp = j*S*self.kypn
            kypn = self.kypn
            kymn = self.kymn
        else:
            aym = ayp = 0.
            bym = byp = 0.
            kypn = kymn = 0.

        if self.nz>1:
            azm = j*S*self.kzmn
            azp = j*S*self.kzpn
            kzpn = self.kzpn
            kzmn = self.kzmn
        else:
            azm = azp = 0.
            bzm = bzp = 0.
            kzpn = kzmn = 0.

        self.BJmult = BJmult
        self.EJmult = EJmult
        self.ERhomult = ERhomult
        self.ERhooldmult = ERhooldmult
        self.Jmult = Jmult
        self.Soverk = Soverk
        self.X1=X1
        self.X2=X2
        self.X3=X3
        self.kVow = kVow
        self.So1mT = So1mT
        self.onemCo1mT = onemCo1mT
        self.kV=kV
        self.T=T
        self.CT = CT
        self.ST = ST

        mymat = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'g':axm/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'g':aym/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'g':azm/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})
        else:
            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})

        if self.l_pushf:
            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'f':axp,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'f':ayp,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'f':azp,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})
        else:
            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})

        if self.l_pushf:
            print 'l_pushf not yet implemented in PSATD PseudoGalilean'
            raise
            mymat.add_op('f',{'f':C,'ex':axm,'ey':aym,'ez':azm, \
                                    'jx': kxmn*FJmult,'jy': kymn*FJmult,'jz': kzmn*FJmult, \
                                    'rhonew':FRhomult,\
                                    'rhoold':-FRhomult - Soverk/self.eps0})

        if self.l_pushg:
            mymat.add_op('g',{'g':C,'bx':axp*c,'by':ayp*c,'bz':azp*c})

        return mymat.mat

    def push(self):

        self.push_fields()

        return
