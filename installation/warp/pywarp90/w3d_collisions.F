#include "top.h"
c=============================================================================
c# Copyright (c) 2007-2007, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This has routines handling particle collisions
c  Alex Friedman, LLNL, (925)422-0827
c  David P. Grote, LLNL, (925)423-7194
c  Bruce Cohen, LLNL, (925)422-9823
c=============================================================================
      subroutine langevincollisions3d(lself,np,ux1,uy1,uz1,
     &                                vxmi,vymi,vzmi,density2,vthsqi2,
     &                                q1,q2,m1,m2,vthsqinit1,vthsqinit2,
     &                                dt,loglambda,epvth)
      use Constant
      use Beam_acc, Only: lrelativ
      logical(ISZ):: lself
      integer(ISZ):: np
      real(kind=8):: ux1(np),uy1(np),uz1(np)
      real(kind=8):: vxmi(np),vymi(np),vzmi(np),density2(np),vthsqi2(np)
      real(kind=8):: q1,q2,m1,m2,vthsqinit1,vthsqinit2
      real(kind=8):: dt,loglambda(np),epvth

c CIC collision operator, 2 space, 3 velocity dimensions
c Implements Manheimer, Lampe, Joyce, JCP 138, 563 (1997)
c test particle species i1 scatters of field particle species i2

c nudt input to this routine corresponds to:
c     nudt = dt * (2^1.5 /3pi^.5) nu_NRL (i1,i1) for species i1
c where  nu_NRL = (4pi Z_1^4 e^4 n_1 ln Lambda)/(m_1^2 v_th_1^3)

c For collisions of i1 on i2, this routine internally converts
c nudt to nudt = dt * (2^1.5 /3pi^.5) nu_NRL (i1,i2) for species i1 on i2
c where  nu_NRL(1,2) = (4pi Z_1^2 Z_2^2 e^4 n_2 ln Lambda)/(m_1^2 v_th_2^3)

c This routine is to be called every ncint timesteps.
c ncint = (1/nudt)*(1/ncoll) where nudt=nu_i *dt and ncoll=no. of
c of times we wish to apply collision operator over the period 1/nudt
c nudt and ncoll are input by user.

c Note:
c v_thermal near edges of plasma is prone to big errors that
c can cause a runaway of particle energy.--DON'T NEED TO WORRY THIS.

      integer(ISZ):: i
      real(kind=8):: sqrt2,twopi,twoorootpi
      real(kind=8):: nudtcoeff,nudt,eta,vthsq2,vmx2sq,v2,gamma2,vx,vy,vz,udotv
      real(kind=8):: uxf,uyf,uzf,uxf0
      real(kind=8):: uu,uui,uu1,uu2,vv,fdv,d11,d33,gamma,gamsq,G,P0
      real(kind=8):: xrand1
      real(kind=8):: aa,norm,const,ww,erf
      real(kind=8):: wranf

c     --- First, calculate the part of the collision frequency that
c     --- depends on constants or global quantities.
      nudtcoeff = dt*2**1.5/(3.*sqrt(pi))*
     &            q1**2*q2**2/m1**2/(4.*pi*eps0**2)

c     --- Weight factor to the thermal velocity incorporated to deal with
c     --- a cooling instability.
c     --- See Cohen et al PoP 13, 22705 (2006) after eq.(4).
      if(epvth == 0.) epvth=0.95

c     --- Useful constants
      sqrt2 = sqrt(2.)
      twopi = 2.*pi
      twoorootpi = 2./sqrt(pi)

c     --- Mass scaling factor. Note that for self collisions, eta=1
      eta = 0.5*(1.+m1/m2)

c     --- Needed for Sherlock form of approximation
      aa = (8./(3.*pi))*(pi-3.)/(4.-pi)
      norm = (3.*sqrt(pi))/4.
      const = 2./(3.*sqrt(pi))

c     --- Now loop over the test particles.
      do i=1,np

        if (vthsqi2(i) < 1.e-6) cycle

        vthsq2 = epvth*vthsqi2(i) + (1. - epvth)*vthsqinit2
        nudt = nudtcoeff*density2(i)*loglambda(i)/vthsq2**1.5
        nudt = abs(nudt)

        if(nudt >= 0.99) then
c         --- This should print some kind of warning perhaps?
          write(*,*) "nudt=",nudt
        endif

c Compute particle velocities in local drift frame

        if (lrelativ) then

c         --- Lorentz
          vmx2sq = vxmi(i)**2 + vymi(i)**2 + vzmi(i)**2
          v2 = sqrt(vmx2sq)
          if (v2 .ne. 0.) then
            gamma2 = sqrt(1. + vmx2sq/clight**2)
            vx = vxmi(i)/v2
            vy = vymi(i)/v2
            vz = vzmi(i)/v2
            P0 = sqrt(1. + (ux1(i)**2 + uy1(i)**2 + uz1(i)**2)/clight**2)
            udotv = ux1(i)*vx + uy1(i)*vy + uz1(i)*vz
            uxf = ux1(i) + (gamma2 - 1.)*udotv*vx - vxmi(i)*P0
            uyf = uy1(i) + (gamma2 - 1.)*udotv*vy - vymi(i)*P0
            uzf = uz1(i) + (gamma2 - 1.)*udotv*vz - vzmi(i)*P0
          else
            uxf = ux1(i)
            uyf = uy1(i)
            uzf = uz1(i)
          endif

        else

c         --- Gallilean
          uxf = ux1(i) - vxmi(i)
          uyf = uy1(i) - vymi(i)
          uzf = uz1(i) - vzmi(i)

        endif

c Compute 3D speed and some other speed normalizations for relative
c angles in coordinate frame oriented parallel to local velocity

        uu  =  sqrt(uxf**2+uyf**2+uzf**2) + 1.e-06*sqrt(vthsqinit1)
        if (uu == 0.) cycle

        if (lrelativ) then
          gamsq = 1. + (uxf*uxf + uyf*uyf + uzf*uzf)/clight**2
          gamma = sqrt(gamsq)
        else
          gamsq = 1.
          gamma = 1.
        endif
        vv = uu/gamma
        ww = uu/sqrt(2.*vthsqi2(i))

c-----------------------------------------------------------------------------
c       --- Calculate diffusion and drag coefficients for isotropic F(v)
c       --- and resulting kicks.

c       --- Note that the relativistic corrections derived by Braams and
c       --- Karney are not important assuming that the thermal velocity of the
c       --- field species is small compared to c. The corrections (ratios
c       --- of K Bessel functions) are all close to 1.

c       --- The quantity used in the calculation of F and D is the velocity, v,
c       --- instead of the momentum, u. This matches the use of v for the
c       --- test particles in the form of F and D from Braams and Karney.
c       --- Note that this gives the curious result that F and D have finite
c       --- values as the test particle energy grows large (since v -> c).
c       --- Also, following the form of Braams and Karney, since they write the
c       --- distribution function f in terms of u, the corresponding Langevin
c       --- equations are also written in terms of u, i.e. du/dt = F*dt + Q.

c       --- Approximation to the error function due to S. Winitzki.
c       --- This expression is very good, given an accuracy of
c       --- approximately 1 part in 10**4 or better over all ww.
        erf = sqrt(1. - exp(-ww**2*((4./pi) + aa*ww**2)/(1. + aa*ww**2)))

c       --- This approximate form is due to Sherlock,
c       --- JCP 227 (2008) 2286-2292.
c       --- Note that the exact form below is preferred.
c       fdv = abs(eta*nudt/(1. + (2.*const)*ww**3))

c       d33 = sqrt(abs(2.*vthsqi2(i)*(fdv/eta)*abs(-log(1. - wranf()))))
c    &        *cos(twopi*wranf())

c       d11 = sqrt(nudt*2.*(vthsqi2(i)*norm)
c    &             *(erf/ww - (1.*const)/(1. + (2.*const)*ww**3))
c    &             *abs(-log(1. - wranf())))

c       --- The erf is calculated anyway for the Sherlock form, so there is
c       --- no reason not to use it and use the exact expressions.
c       --- This seems to give a noticably better result for the two species,
c       --- unequal temperature test case, giving results closer to the analytic
c       --- expressions for equilibration over long times, nu*t ~ 1. This is
c       --- unexpected since the difference in the expressions is small, ~5%.

c       --- In Braams and Karney, the expression for d11 contains the term
c       --- (1/u**2 + 1/(gamma**2*cvac**2)). For small v << c, this gives
c       --- approximately 1/v**2. For large v ~ c, this is << 1. So, for this
c       --- term, the difference between using v and u is insignificant, so the
c       --- expression here only uses v.

        G = (erf - twoorootpi*ww*exp(-ww**2))/(2*ww**2)

        fdv = abs(eta*nudt*G/ww/const)

        d33 = sqrt(abs(2.*vthsqi2(i)*(fdv/eta)*abs(-log(1. - wranf()))))
     &        *cos(twopi*wranf())

        d11 = sqrt(nudt*2.*vthsqi2(i)*norm*((erf - G)/ww)*abs(-log(1.-wranf())))

c       --- This is a Pade approximation from B. Cohen.
c       --- The Sherlock and exact forms are better. This does not have the
c       --- correct form at large v.
c       fdv = abs( eta*nudt/(1.+1.06483*(uu/(sqrt2*sqrt(vthsqi2(i))))**2.488) )
c       d33 = sqrt(abs(2.*vthsqi2(i)*(fdv/eta)*abs(-log(1.-wranf()))))
c    &      *cos(twopi*wranf())
c       d11 = sqrt(nudt*2.*vthsqi2(i)/(1.+0.1856*
c    &      (uu/(sqrt2*sqrt(vthsqi2(i))))**1.848)*abs(-log(1.-wranf())))

c need to multiply d11 by cos and sin of ranf*twopi for kicks
c-----------------------------------------------------------------------------

c       --- Includes energy cons. correction ~dt**2
c       --- Note that the extra gamma factor comes about when the Langevin
c       --- equation is written in terms of u, and so the correction term
c       --- becomes dF = -dt/(2*u)*F**2. Also note that the fdv calculated
c       --- above is actually Fd/v. */
        fdv = -vv*fdv*(1.+0.5*fdv/gamma)

c       --- Do the Monte Carlo collisional kicks in the particle velocity
c       --- local frame and transform back to x-y-z velocity coordinates
c       --- in drift frame
        xrand1 = wranf()

c       --- Later on we need a duplicate copy of uxf (after uxf is
c       --- overwritten with a new value after the collision) for computing a
c       --- cosine = uxf/uu1.
        uxf0 = uxf

c       --- reduce number of divides by precalculating 1./uu
        uui = 1./uu
        uu1 =  sqrt(uxf**2+uyf**2)+ 1.e-06*sqrt(vthsqinit1)

        if (uu1 .ne. 0.) then
          uu2 =  1./(uu*uu1 +1.e-06*vthsqinit1)

          uxf = uxf + (fdv+d33)*uxf*uui + d11*cos(twopi*xrand1)*uyf/uu1 
     &        + d11*sin(twopi*xrand1)*uzf*uxf*uu2

          uyf = uyf + (fdv+d33)*uyf*uui + d11*cos(twopi*xrand1)*(-uxf0)/uu1 
     &        + d11*sin(twopi*xrand1)*uzf*uyf*uu2

          uzf = uzf + (fdv+d33)*uzf*uui + d11*sin(twopi*xrand1)*(-uu1**2)*uu2

        else

          uzf = uzf + (fdv+d33)*uzf*uui

        endif

c Galilean tranformation to lab velocity coordinates from drift frame

        if (lrelativ) then

c         --- Lorentz
          if (v2 .ne. 0.) then
            P0 = sqrt(1. + (uxf**2 + uyf**2 + uzf**2)/clight**2)
            udotv = uxf*vx + uyf*vy + uzf*vz
            ux1(i) = uxf + (gamma2 - 1.)*udotv*vx + vxmi(i)*P0
            uy1(i) = uyf + (gamma2 - 1.)*udotv*vy + vymi(i)*P0
            uz1(i) = uzf + (gamma2 - 1.)*udotv*vz + vzmi(i)*P0
          else
            ux1(i) = uxf
            uy1(i) = uyf
            uz1(i) = uzf
          endif

        else

c         --- Gallilean
          ux1(i) = uxf + vxmi(i)
          uy1(i) = uyf + vymi(i)
          uz1(i) = uzf + vzmi(i)

        endif

      enddo   !!loop over particle blocks

      return
      end  !!langevincollisions3d
c=============================================================================

c=============================================================================
c=============================================================================
c=============================================================================

c=============================================================================

      subroutine ijcoll2d(i1,i2,np1,np2,nxp,nyp,ux1,ux2,uy1,uy2,
     &   uz1,uz2,epvth,x1,y1,x2,y2,lx,ly,aion1,zion1,aion2,zion2,nion1,nion2,
     &   vthmisq1,vthmisq2,nudt,ncint,ixbc,izbc)
      use Constant
cc CIC collision operator, 2 space, 3 velocity dimensions
cc Implements Manheimer, Lampe, Joyce, JCP 138, 563 (1997)
cc test particle species i1 scatters of field particle species i2

cc nudt input to this routine corresponds to:
cc     nudt = dt * (2^1.5 /3pi^.5) nu_NRL (i1,i1) for species i1
cc where  nu_NRL = (4pi Z_1^4 e^4 n_1 ln Lambda)/(m_1^2 v_th_1^3)

cc For collisions of i1 on i2, this routine internally converts
cc nudt to nudt = dt * (2^1.5 /3pi^.5) nu_NRL (i1,i2) for species i1 on i2
cc where  nu_NRL(1,2) = (4pi Z_1^2 Z_2^2 e^4 n_2 ln Lambda)/(m_1^2 v_th_2^3)

cc This routine is to be called every ncint timesteps.

cc ixbc=0,1,2  izbc=0,1,2  bounded,bounded,periodic

cc Add whatever "Use" statements as needed here

      integer i,ip,ix,iy,ix1,iy1,ncint,nx1ny1,ixbc,izbc
      integer i1,i2,np1,np2,nxp,nyp

      real escale,eta,vthmisqa,vthmisqb,aiona,ziona,aionb,zionb
      real vthmisqhyd,vtha,vthb,na,nb,vth1,lx,ly
      real ke1tot,ke2tot,delke12,twopi,wxa,wxb,wya,wyb,s0i

      real ux1(np1),ux2(np2),uy1(np1),uy2(np2),uz1(np1),uz2(np2)
      real x1(np1),y1(np1),x2(np2),y2(np2),xf,yf

      real nuadt,s00
      real ddxp,ddyp,epvth	
      real etotf1,etotf2,vv,vv1,vv2,fdv,d33,d11,xrand1,uxf,uxf0,uyf,uzf
      real sqrt2
      real aion1,zion1,aion2,zion2,nion1,nion2,vthmisq1,vthmisq2
      real nudt,vthsqi
      real sddx2d,wranf
      real aa,norm,const,ww,erf



cc Declare Allocatable new variables

      real, allocatable, dimension(:) :: pdmy2          ! dummy array
      real, allocatable, dimension(:) :: pdmy1          ! dummy array
      real, allocatable, dimension(:) :: vxmi           ! local mean x velocity of background
      real, allocatable, dimension(:) :: vymi           ! local mean y velocity of background
      real, allocatable, dimension(:) :: vzmi           ! local mean z velocity of background


      real, allocatable, dimension(:,:) :: s0             ! density array
      real, allocatable, dimension(:,:) :: svx             ! x-momentum array
      real, allocatable, dimension(:,:) :: svy             ! y-momentum array
      real, allocatable, dimension(:,:) :: svz             ! z-momentum array
      real, allocatable, dimension(:,:) :: se             ! energy array
      integer:: alloc_status

      allocate(pdmy1(np1), pdmy2(np2), vxmi(np1), vymi(np1), vzmi(np1),
     &         s0(0:nxp,0:nyp), se(0:nxp,0:nyp),
     &         svx(0:nxp,0:nyp), svy(0:nxp,0:nyp), svz(0:nxp,0:nyp), stat=alloc_status)
      if (alloc_status /= 0) then
        print*,"Error:ijcoll2d: pdmy1 et al could not be allocated"
      endif

cc The strategy to use here is to pass the velocity and position
cc arrays for both species needed in the collisions

cc The test particle is species 1 and the field particle is species 2


c     --- Needed for Sherlock for of approximation
      aa = (8./(3.*pi))*(pi-3.)/(4.-pi)
      norm = (3.*sqrt(pi))/4.
      const = 2./(3.*sqrt(pi))

      ddxp=sddx2d(lx,nxp)
      ddyp=sddx2d(ly,nyp)
      sqrt2 = sqrt(2.)

      nx1ny1 = (nxp+1)*(nyp+1)

cncoll ; ncint; nudt ; step ; i1; i2
c call this routine every ncint timesteps

c ncint = (1/nudt)*(1/ncoll) where nudt=nu_i *dt and ncoll=no. of
c of times we wish to apply collision operator over the period 1/nudt
c nudt and ncoll are input by user. Then ncint is calculated in ib_deck
c before calling this routine.

      twopi=8.*atan(1.)

      ke1tot=sum(ux1**2+uy1**2+uz1**2)/(3.*np1*vthmisq1)  ! ke1tot
      ke2tot=sum(ux2**2+uy2**2+uz2**2)/(3.*np2*vthmisq2); ! ke2tot


c call pdiag to compute moments of 2D particle array on grid

      pdmy1=0.
      pdmy2=0.



      call pdiag2d(x2,y2,ux2,uy2,uz2,pdmy2,np2,lx,ly,nxp,nyp,s0,
     &        svx,svy,svz,se,ixbc,izbc)

cc      write(*,*) "sum(s0)/((nxp+1)*(nyp+1))=",sum(s0)/((nxp+1)*(nyp+1))

cc      write(*,*) "x2(1:np2)=",x2(1:np2)
cc      write(*,*) "y2(1:np2)=",y2(1:np2)
cc      write(*,*) "ux2(1:np2)=",ux2(1:np2)
cc      write(*,*) "uy2(1:np2)=",uy2(1:np2)
cc      write(*,*) "uz2(1:np2)=",x2(1:np2)
cc      write(*,*) "s0(0:nxp,0:nyp)=",s0(0:nxp,0:nyp)
cc      write(*,*) "se(0:nxp,0:nyp)=",se(0:nxp,0:nyp)


      vth1 = sqrt(vthmisq1)

cc Need to define charges and masses of species 1 and 2, depends
cc on the particular code, nion1,2 & vthmisq1,2 are species 1,2 (absolute index)
cc while i1 and i2 are relative indices for collisions

      vtha = sqrt(vthmisq1)
      na=nion1
      ziona=zion1
      aiona=aion1

      vthb = sqrt(vthmisq2)
      nb=nion2
      zionb=zion2
      aionb=aion2


      if(i1.eq.i2) then
         eta = 1.
      else
         eta = 0.5*(1.+aion1/aion2)
       endif

c pdiag returns <(v_x-<v_x>)^2+"y & z">/3.=v_th^2


c v_thermal near edges of plasma is prone to big errors that
c can cause a runaway of particle energy.--DON'T NEED TO WORRY THIS.

      vthmisqb=vthmisq2
      vthmisqa=vthmisq1

cCollision rate nudt is input = nu_eff (1st species, whatever it is)
c and frquency ncint is computed elsewhere


      etotf1 = 0.
      etotf2 = 0.

cFor i1 on i2, with nudt defined relative to first species:

      do i=1,np1

         nuadt = nudt*(nb/nion1)*float(ncint)*(ziona/zion1)**2*
     &         (zionb/zion1)**2*(vth1/vthb)**3*(aion1/aiona)**2


         xf = x1(i)
         yf = y1(i)

         If(i.eq.1) write(*,*) "nudt,nuadt=",nudt,nuadt
 
c Compute grid indices near particle from (x,y), indices span 0 to nx, 0 to ny

         ix = mod(int(xf*ddxp),nxp+1)
         ix1 = mod(int(xf*ddxp)+1,nxp+1)

         iy = mod(int(yf*ddyp),nyp+1)
         iy1 = mod(int(yf*ddyp)+1,nyp+1)

c Compute linear interpolation factors

         wxb = xf*ddxp-aint(xf*ddxp)
         wxa = 1.-wxb
         wyb = yf*ddyp-aint(yf*ddyp)
         wya = 1.-wyb



c Compute local average velocities and vx**2+vy**2 using linear interpolation

         vxmi(i) = wxa*wya*svx(ix,iy)+wxb*wya*svx(ix1,iy)
     &          + wxa*wyb*svx(ix,iy1) + wxb*wyb*svx(ix1,iy1)

         vymi(i) = wxa*wya*svy(ix,iy)+wxb*wya*svy(ix1,iy)
     &          + wxa*wyb*svy(ix,iy1) + wxb*wyb*svy(ix1,iy1)
         vzmi(i) = wxa*wya*svz(ix,iy)+wxb*wya*svz(ix1,iy)
     &          + wxa*wyb*svz(ix,iy1) + wxb*wyb*svz(ix1,iy1)

         vthsqi = wxa*wya*se(ix,iy)+wxb*wya*se(ix1,iy)+1.e-6*vthmisqb
     &          + wxa*wyb*se(ix,iy1) + wxb*wyb*se(ix1,iy1)
         vthsqi = abs(vthsqi)

         s0i = wxa*wya*s0(ix,iy)+wxb*wya*s0(ix1,iy)
     &          + wxa*wyb*s0(ix,iy1)+ wxb*wyb*s0(ix1,iy1)

 
         s00 = sum(s0)/nx1ny1 

cc         if(i.eq.10) write(*,*) "s0i,s0(ix,iy),s00,nx1ny1=",
cc     &          s0i,s0(ix,iy),s00,nx1ny1
cc         if(i.eq.10) write(*,*) "s0(1:nxp,1:nyp)=",s0(1:nxp,1:nyp)
cc         if(i.eq.10) write(*,*) "se(1:nxp,1:nyp)=",se(1:nxp,1:nyp)
   
cassumes that the average density is const and moats are negligible


         if(epvth.eq.0.) epvth=0.95		!!epvth=0.95 works fine
c Apply collisions with local correction to the temperature and density
c The nu_0 NRL = 2**1.5 * nu_alpha Jones et al.

         nuadt = nuadt*(s0i/s00)*(vthmisqb/(epvth*vthsqi
     &        +(1.-epvth)*vthmisqb))**1.5
         nuadt= abs(nuadt)

         if(i.eq.1) write(*,*) "s0i,s00,vthsqi,vthmisqb,nuadt=",
     &              s0i,s00,vthsqi,vthmisqb,nuadt

cc         if(i.eq.10) write(*,*) "vthmisqb,vthsqi,nuadt=",vthmisqb,vthsqi,nuadt
 
         if(nuadt.ge.0.99) then
            write(*,*) "nuadt=",nuadt
         endif

c Compute particle velocities in local drift frame

         uxf = ux1(i)-vxmi(i)
         uyf = uy1(i)-vymi(i)
         uzf = uz1(i)-vzmi(i)

c compute the temperature before the collision

         etotf1=(uxf**2+uyf**2+uzf**2)/(3.*np1*vthmisqa)+etotf1


c Compute 3D speed and some other speed normalizations for relative
c angles in coordinate frame oriented parallel to local velocity

         vv  =  sqrt(uxf**2+uyf**2+uzf**2) + 1.e-06*sqrt(vthmisqa)
         vv1 =  sqrt(uxf**2+uyf**2)+ 1.e-06*sqrt(vthmisqa)
         vv2 =  1./(vv*vv1 +1.e-06*vthmisqa)

c-----------------------------------------------------------------------------
c calculate diffusion and drag coefficients for isotropic F(v) and
c resulting kicks, using approximation from Sherlock, JCP 227 (2008) 2286-2292.

         ww = vv/sqrt(2.*vthsqi)

c        --- Approximation to the error function due to S. Winitzki.
         erf = sqrt(1.-exp(-ww**2 * ((4./pi)+aa*ww**2)/(1.+aa*ww**2)) )

         fdv = abs( eta*nuadt/(1.+(2.*const)*ww**3) )

         d33 = sqrt(abs(2.*vthsqi*(fdv/eta)*abs(-log(1.-wranf()))))
     &        *cos(twopi*wranf())

         d11 = sqrt(nuadt*2.*(vthsqi*norm)*(  erf/ww - (1.*const)
     &       /(1.+(2.*const)*ww**3 ) )*abs(-log(1.-wranf())))

c-----------------------------------------------------------------------------
c calculate diffusion and drag coefficients for isotropic F(v) and
c resulting kicks, using Pade approximation from B. Cohen.
c The Sherlock form is better and has correct form at large v.

cc        fdv = abs( eta*nuadt/(1.+1.08*(vv/(sqrt2*sqrt(vthsqi)))**2.76) )
c        fdv = abs( eta*nuadt/(1.+1.06483*(vv/(sqrt2*sqrt(vthsqi)))**2.488) )

c        d33 = sqrt(abs(2.*vthsqi*(fdv/eta)*abs(-log(1.-wranf()))))
c    &       *cos(twopi*wranf())

c        d11 = sqrt(nuadt*2.*vthsqi/(1.+0.1856*
c    &       (vv/(sqrt2*sqrt(vthsqi)))**1.848)*abs(-log(1.-wranf())))
c need to multiply d11 by cos and sin of ranf*twopi for kicks
c-----------------------------------------------------------------------------

         fdv = -vv*fdv*(1.+0.5*fdv)      !! Includes energy cons. correction ~dt**2

c Do the Monte Carlo collisional kicks in the particle velocity local frame
c and transform back to x-y-z velocity coordinates in drift frame

         xrand1=wranf()

c Later on we need a duplicate copy of uxf (after uxf is
c overwritten with a new value after the collision) for computing a
c cosine = uxf/vv1, so stash a copy of uxf in s0i



         uxf0=uxf

c reduce number of divides by replacing vv with 1./vv

         vv=1./vv

         uxf = uxf + (fdv+d33)*uxf*vv + d11*cos(twopi*xrand1)*uyf/vv1 
     &       + d11*sin(twopi*xrand1)*uzf*uxf*vv2

         uyf = uyf + (fdv+d33)*uyf*vv + d11*cos(twopi*xrand1)*(-uxf0)/vv1 
     &       + d11*sin(twopi*xrand1)*uzf*uyf*vv2

         uzf = uzf + (fdv+d33)*uzf*vv + d11*sin(twopi*xrand1)*(-vv1**2)*vv2


c Galilean tranformation to lab velocity coordinates from drift frame

         ux1(i) = uxf
         uy1(i) = uyf
         uz1(i) = uzf

c compute temperature after the collision
         etotf2=(uxf**2+uyf**2+uzf**2)/(3.*np1*vthmisqa)+etotf2

      enddo   !!loop over particle blocks

      escale=sqrt(etotf1/etotf2)

      write(*,*) "escale=",escale

      do i = 1,np1

         uxf = ux1(i)
         uyf = uy1(i)
         uzf = uz1(i)



cScale the particle velocities to conserve energy after the
ccollisions.  Manheimer et al. suggest doing this.
cThis collision scheme does NOT conserve energy like Jones
calgorithm does, except statistically for an ideal Maxwellian
cdistribution of test particles, and erf function drag-diff coeffs,
cnone of which we satisfy.


         if(i1.eq.i2) then
            uxf=uxf*escale
            uyf=uyf*escale
            uzf=uzf*escale
         endif

c Galilean tranformation to lab velocity coordinates from drift frame

         ux1(i) = uxf + vxmi(i)
         uy1(i) = uyf + vymi(i)
         uz1(i) = uzf + vzmi(i)

      enddo



!! optional writes of key collision quantities
      write(*,*)  "nuadt=",nuadt
      write(*,*)   "etotf1=",etotf1
      write(*,*)   "etotf2=",etotf2
      write(*,*)   "escale=",escale


      ke1tot=sum(ux1**2+uy1**2+uz1**2)/(3.*np1*vthmisq1)
      write(*,*) "ke1tot=",ke1tot

      ke2tot=sum(ux2**2+uy2**2+uz2**2)/(3.*np2*vthmisq2); 
      write(*,*) "ke2tot=",ke2tot
!!      delke12 = ke1tot-(mion2/mion1)*ke2tot*(vthmisq2/vthmisq1) !! delke12


ccpaws

cc Deallocate arrays

      deallocate(pdmy1, pdmy2, vxmi, vymi, vzmi)
      deallocate(s0, se, svx, svy, svz)



      end  !!ijcoll2d


c pdiag.m

c   nxp,nzp here can differ from field arrays nx,nz
c    need to be able to do only a part of the system?

      subroutine pdiag2d(x,z,ux,uy,uz,gammai,np,lxp,lzp,nxp,nzp,
     & s0,sx,sy,sz,se,ixbc,izbc)
      
      Implicit None

      integer np, nxp, nzp
      real x(np), z(np), ux(np), uy(np), uz(np), gammai(np)
c  only need for nx,nz

      real s0(0:nxp,0:nzp), sx(0:nxp,0:nzp), sy(0:nxp,0:nzp), sz(0:nxp,0:nzp), se(0:nxp,0:nzp)
      real lxp, lzp

      integer i, ix, iz,ixbc,izbc
      real ddxp, ddzp, rx, rz, wxa, wxb, wza, wzb
      real wxaza, wxazb, wxbza, wxbzb
      real qvc, tx, ty, tz, tem
      real sddx2d

      ddxp = sddx2d(lxp,nxp)
      ddzp = sddx2d(lzp,nzp)

cc      write(*,*) "ddxp=","ddzp=",ddxp,ddzp
cc      write(*,*) "lxp,lzp,nxp,nzp=",lxp,lzp,nxp,nzp
      qvc = ddxp*ddzp

      do ix = 0,nxp
      do iz = 0,nzp
       s0(ix,iz) = 0.
       sx(ix,iz) = 0.
       sy(ix,iz) = 0.
       sz(ix,iz) = 0.
       se(ix,iz) = 0.
      enddo
      enddo

      do i = 1,np
       rx = x(i)*ddxp
       rz = z(i)*ddzp
       ix = rx
       iz = rz
c      Bilinear (area) weighting.
       wxb = rx - ix
       wzb = rz - iz
       wxa = 1. - wxb
       wza = 1. - wzb
       wxaza = wxa*wza
       wxazb = wxa*wzb
       wxbza = wxb*wza
       wxbzb = wxb*wzb
       s0(ix  ,iz  ) = s0(ix  ,iz  ) + qvc*wxaza
       s0(ix  ,iz+1) = s0(ix  ,iz+1) + qvc*wxazb
       s0(ix+1,iz  ) = s0(ix+1,iz  ) + qvc*wxbza
       s0(ix+1,iz+1) = s0(ix+1,iz+1) + qvc*wxbzb
       tx = qvc*ux(i)
       sx(ix  ,iz  ) = sx(ix  ,iz  ) + tx*wxaza
       sx(ix  ,iz+1) = sx(ix  ,iz+1) + tx*wxazb
       sx(ix+1,iz  ) = sx(ix+1,iz  ) + tx*wxbza
       sx(ix+1,iz+1) = sx(ix+1,iz+1) + tx*wxbzb
       ty = qvc*uy(i)
       sy(ix  ,iz  ) = sy(ix  ,iz  ) + ty*wxaza
       sy(ix  ,iz+1) = sy(ix  ,iz+1) + ty*wxazb
       sy(ix+1,iz  ) = sy(ix+1,iz  ) + ty*wxbza
       sy(ix+1,iz+1) = sy(ix+1,iz+1) + ty*wxbzb
       tz = qvc*uz(i)
       sz(ix  ,iz  ) = sz(ix  ,iz  ) + tz*wxaza
       sz(ix  ,iz+1) = sz(ix  ,iz+1) + tz*wxazb
       sz(ix+1,iz  ) = sz(ix+1,iz  ) + tz*wxbza
       sz(ix+1,iz+1) = sz(ix+1,iz+1) + tz*wxbzb
      enddo

c  Establish periodicity in z.-- assume periodic in z
cc      if(nyp = ny) then

      if(izbc.eq.2) then

         do ix = 0,nxp
            s0(ix, 0) = s0(ix,0) + s0(ix,nzp)
            s0(ix,nzp) = s0(ix,0)
            sx(ix, 0) = sx(ix,0) + sx(ix,nzp)
            sx(ix,nzp) = sx(ix,0)
            sy(ix, 0) = sy(ix,0) + sy(ix,nzp)
            sy(ix,nzp) = sy(ix,0)
            sz(ix, 0) = sz(ix,0) + sz(ix,nzp)
            sz(ix,nzp) = sz(ix,0)
         enddo
      endif

c  Establish periodicity in x.
c***not if open sided bc
       if(ixbc.eq.2) then
          do iz = 0,nzp
             s0(0, iz) = s0(0,iz) + s0(nxp,iz)
             s0(nxp,iz) = s0(0,iz)
             sx(0, iz) = sx(0,iz) + sx(nxp,iz)
             sx(nxp,iz) = sx(0,iz)
             sy(0, iz) = sy(0,iz) + sy(nxp,iz)
             sy(nxp,iz) = sy(0,iz)
             sz(0, iz) = sz(0,iz) + sz(nxp,iz)
             sz(nxp,iz) = sz(0,iz)
          enddo
       endif

c  For now, just divide by s0 to get mean velocities.
      do ix = 0,nxp
      do iz = 0,nzp
       if(s0(ix,iz).ne. 0.) then
        sx(ix,iz) = sx(ix,iz)/s0(ix,iz)
        sy(ix,iz) = sy(ix,iz)/s0(ix,iz)
        sz(ix,iz) = sz(ix,iz)/s0(ix,iz)
       endif
      enddo
      enddo

      do i = 1,np
       rx = x(i)*ddxp
       rz = z(i)*ddzp
       ix = rx
       iz = rz
c      Bilinear (area) weighting.
       wxb = rx - ix
       wzb = rz - iz
       wxa = 1. - wxb
       wza = 1. - wzb
       wxaza = wxa*wza
       wxazb = wxa*wzb
       wxbza = wxb*wza
       wxbzb = wxb*wzb
       tx = wxaza*sx(ix,iz  ) + wxbza*sx(ix+1,iz) +
     &            wxazb*sx(ix,iz+1) + wxbzb*sx(ix+1,iz+1)
       ty = wxaza*sy(ix,iz  ) + wxbza*sy(ix+1,iz) +
     &            wxazb*sy(ix,iz+1) + wxbzb*sy(ix+1,iz+1)
       tz = wxaza*sz(ix,iz  ) + wxbza*sz(ix+1,iz) +
     &            wxazb*sz(ix,iz+1) + wxbzb*sz(ix+1,iz+1)
       tem = ((ux(i) - tx)**2 + (uy(i) - ty)**2 + (uz(i) - tz)**2)/3.
       se(ix  ,iz  ) = se(ix  ,iz  ) + tem*wxaza*qvc
       se(ix  ,iz+1) = se(ix  ,iz+1) + tem*wxazb*qvc
       se(ix+1,iz  ) = se(ix+1,iz  ) + tem*wxbza*qvc
       se(ix+1,iz+1) = se(ix+1,iz+1) + tem*wxbzb*qvc
cc       write(*,*) "i,x(i),z(i)=",i,x(i),z(i)
cc       write(*,*) "ix,iz=",ix,iz
cc       write(*,*) "i,rx,rz,wxb,wzb,se(ix,iz),se(ix,iz+1),se(ix+1,iz),
cc     &     se(ix+1,iz+1)="
cc       write(*,*) i,rx,rz,wxb,wzb,se(ix,iz),se(ix,iz+1),se(ix+1,iz),
cc     &     se(ix+1,iz+1)
  
 
      enddo

c  Establish periodicity in z.  assume periodic in z
cc      if(nyp = ny) then
      if(izbc.eq.2) then
         do ix = 0,nxp
            se(ix, 0) = se(ix,0) + se(ix,nzp)
            se(ix,nzp) = se(ix,0)
         enddo
      endif

c  Establish periodicity in x.
c***not if open sided bc
      if(ixbc.eq.2) then
         do iz = 0,nzp
            se(0, iz) = se(0,iz) + se(nxp,iz)
            se(nxp,iz) = se(0,iz)
         enddo
      endif

c  Now have thermal energy, except for missing mass factor

c  For now, just divide by s0 to get mean particle thermal energy.
      do ix = 0,nxp
      do iz = 0,nzp
       if(s0(ix,iz).ne. 0.) se(ix,iz) = se(ix,iz)/s0(ix,iz)
      enddo
      enddo

cc      do ix = 0,nxp
cc         do iz = 0,nzp
cc            write(*,*) "ix,iz,s0,se,svx,svy,svz=",
cc     &              ix,iz,s0(ix,iz),se(ix,iz),sx(ix,iz)
cc     &              ,sy(ix,iz),sz(ix,iz)
cc            write(*,*)

cc         enddo
cc      enddo

      return
      end

      real function sddx2d(lx,nx)

      Implicit None
      real lx
      integer nx
      real ddx, eps, loddx, hiddx
c  sddx = nx/lx, or is reduced to guarantee lx*ddx <= nx.

      sddx2d = nx/lx
cc      write(*,*) "sddx2d=",sddx2d
      eps = lx*sddx2d - nx
cc      write(*,*) "eps=",eps
      if(eps .lt. 0.) return

c  Sigh. Bisection method to find a suitable ddx.
      ddx = sddx2d
      loddx = ddx*(1-2.e-7)
      hiddx = ddx
      do
        ddx = 0.5*(loddx + hiddx)
        eps = lx*ddx - nx
        if(ddx .eq. loddx)goto 666
        if(ddx .eq. hiddx) then
          ddx = loddx
          eps = lx*ddx - nx
          goto 666
        endif
        if(eps < 0) then
          loddx = ddx
        else
          hiddx = ddx
          ddx = loddx
        endif
      enddo
 666  continue
      sddx2d = ddx
cc      write(*,*) "sddx2d=",sddx2d

      if(lx*ddx .ge. nx) then
        write(*,*) 'sddx2d error'
      endif
      return
      end

c=============================================================================
      subroutine compton_scatter_electrons(photon_energy, photon_density, photon_px, photon_py, photon_pz,
     &                                     nx, ny, nz, dx, dy, dz, xmmin, ymmin, zmmin, zgrid0, dt,
     &                                     boost_gamma, time, weightingmethod, anglethreshold,
     &                                     electrongroup, electrons, photongroup, photons, tbirthpid, scattercountpid,
     &                                     withgammaeffective, scatterelectron)
      use Constant
      use ParticleGroupModule
      use Particles,Only: wpid
      integer(ISZ):: nx, ny, nz
      real(kind=8):: photon_energy
      real(kind=8):: photon_density(0:nx,0:ny,0:nz)
      real(kind=8):: photon_px(0:nx,0:ny,0:nz)
      real(kind=8):: photon_py(0:nx,0:ny,0:nz)
      real(kind=8):: photon_pz(0:nx,0:ny,0:nz)
      real(kind=8):: dx, dy, dz, dt
      real(kind=8):: xmmin, ymmin, zmmin, zgrid0
      real(kind=8):: boost_gamma, time
      integer(ISZ):: weightingmethod
      real(kind=8):: anglethreshold
      type(ParticleGroup):: electrongroup
      type(ParticleGroup):: photongroup
      integer(ISZ):: photons
      integer(ISZ):: electrons
      integer(ISZ):: tbirthpid, scattercountpid
      logical(ISZ):: withgammaeffective
      logical(ISZ):: scatterelectron

c zgrid0 is the location of the grid at time = 0

      integer(ISZ):: ie0, ie, ig
      integer(ISZ):: ix, iy, iz
      integer(ISZ):: nnewspace
      integer(ISZ):: nnewphotons, inew

      real(kind=8):: x
      real(kind=8):: rn1, rn2, c0, b, a
      real(kind=8):: k, kprime
      real(kind=8):: ratio, arg
      real(kind=8):: theta, phi
      real(kind=8):: e_radius, e_area, sigma, vz, zswept
      real(kind=8):: numer, denom
      real(kind=8):: factor1, factor2
      real(kind=8):: prob
      real(kind=8):: f, h, factor

      real(kind=8):: boost_uz
      real(kind=8):: xlab, ylab, zlab, tlab
      real(kind=8):: uxlab, uylab, uzlab, gamma_lab
      real(kind=8):: ux, uy, uz
      real(kind=8):: gamma_frame, umag
      real(kind=8):: uxframe, uyframe, uzframe
      real(kind=8):: uxhat, uyhat, uzhat, udotv, pdotv
      real(kind=8):: zgrid
      real(kind=8):: pdens, px, py, pz, pe, p_tot, p_perp
      real(kind=8):: laser_a0sq
      real(kind=8):: laser_Aamp, laser_Ax, laser_Ay
      real(kind=8):: uxinlaser, uyinlaser, uzinlaser, gammainlaser
      !real(kind=8):: pi0, pz_effective, gamma_effective
      real(kind=8):: omega, kz, eta, gammae, gammabetae, gamma_effective
      real(kind=8):: px_old, py_old, pz_old
      real(kind=8):: p_old, e_old, p_new, e_new
      real(kind=8):: pxtemp, pytemp, pztemp, ptemp
      real(kind=8):: pxlab, pylab, pzlab, pelab
      real(kind=8):: delta_px, delta_py, delta_pz, delta_e, gamma_new
      real(kind=8):: psx, psy, psz
      real(kind=8):: wx, wy, wz
      real(kind=8):: cphi, sphi, ctheta, stheta
      real(kind=8):: rrr

      integer(ISZ):: nscattered
      integer(ISZ), save:: nnewphotonstotal = 10000

      if (electrongroup%nps(electrons) == 0) return

      boost_uz = clight*sqrt(boost_gamma*boost_gamma - 1.)

      e_radius = echarge**2/(4.0*pi*eps0*emass*clight**2)
      e_area = pi*e_radius**2

      ! --- Make some space to start with.
      nnewspace = 2*nnewphotonstotal
      call chckpart(photongroup,photons,0,nnewspace)

      ! --- Loop over electrons
      nscattered = 0
      nnewphotonstotal = 0
      do ie0=0, electrongroup%nps(electrons)-1

        ! --- Note that with handling the index this way, the correct
        ! --- electron index will be used even if the data is shifted
        ! --- around as new photons are added.
        ie = ie0 + electrongroup%ins(electrons)

        ! --- Boost particle to lab frame
        xlab = electrongroup%xp(ie)
        ylab = electrongroup%yp(ie)
        zlab = boost_gamma*electrongroup%zp(ie) + boost_uz*time
        tlab = boost_gamma*time + boost_uz*electrongroup%zp(ie)/clight**2
        uxlab = electrongroup%uxp(ie)
        uylab = electrongroup%uyp(ie)
        uzlab = boost_gamma*electrongroup%uzp(ie) + boost_uz/electrongroup%gaminv(ie)
        gamma_lab = sqrt(1. + (uxlab**2 + uylab**2 + uzlab**2)/clight**2)

        zgrid = zgrid0 - clight*tlab
        ix = int((xlab - xmmin)/dx)
        iy = int((ylab - ymmin)/dy)
        iz = int((zlab - zmmin - zgrid)/dz)
        wx =     (xlab - xmmin)/dx - ix
        wy =     (ylab - ymmin)/dy - iy
        wz =     (zlab - zmmin - zgrid)/dz - iz

        if (ix < 0 .or. ix >= nx .or.
     &      iy < 0 .or. iy >= ny .or.
     &      iz < 0 .or. iz >= nz) cycle

        ! --- Get the photon density at the electron location (in lab frame)
        pdens = photon_density(ix  ,iy  ,iz  )*(1. - wx)*(1. - wy)*(1. - wz) +
     &          photon_density(ix+1,iy  ,iz  )*      wx *(1. - wy)*(1. - wz) +
     &          photon_density(ix  ,iy+1,iz  )*(1. - wx)*      wy *(1. - wz) +
     &          photon_density(ix+1,iy+1,iz  )*      wx *      wy *(1. - wz) +
     &          photon_density(ix  ,iy  ,iz+1)*(1. - wx)*(1. - wy)*      wz  +
     &          photon_density(ix+1,iy  ,iz+1)*      wx *(1. - wy)*      wz  +
     &          photon_density(ix  ,iy+1,iz+1)*(1. - wx)*      wy *      wz  +
     &          photon_density(ix+1,iy+1,iz+1)*      wx *      wy *      wz 

        ! --- Get the momentum of the incoming photon (in lab frame)
        pxtemp = photon_px(ix  ,iy  ,iz  )*(1. - wx)*(1. - wy)*(1. - wz) +
     &           photon_px(ix+1,iy  ,iz  )*      wx *(1. - wy)*(1. - wz) +
     &           photon_px(ix  ,iy+1,iz  )*(1. - wx)*      wy *(1. - wz) +
     &           photon_px(ix+1,iy+1,iz  )*      wx *      wy *(1. - wz) +
     &           photon_px(ix  ,iy  ,iz+1)*(1. - wx)*(1. - wy)*      wz  +
     &           photon_px(ix+1,iy  ,iz+1)*      wx *(1. - wy)*      wz  +
     &           photon_px(ix  ,iy+1,iz+1)*(1. - wx)*      wy *      wz  +
     &           photon_px(ix+1,iy+1,iz+1)*      wx *      wy *      wz
        pytemp = photon_py(ix  ,iy  ,iz  )*(1. - wx)*(1. - wy)*(1. - wz) +
     &           photon_py(ix+1,iy  ,iz  )*      wx *(1. - wy)*(1. - wz) +
     &           photon_py(ix  ,iy+1,iz  )*(1. - wx)*      wy *(1. - wz) +
     &           photon_py(ix+1,iy+1,iz  )*      wx *      wy *(1. - wz) +
     &           photon_py(ix  ,iy  ,iz+1)*(1. - wx)*(1. - wy)*      wz  +
     &           photon_py(ix+1,iy  ,iz+1)*      wx *(1. - wy)*      wz  +
     &           photon_py(ix  ,iy+1,iz+1)*(1. - wx)*      wy *      wz  +
     &           photon_py(ix+1,iy+1,iz+1)*      wx *      wy *      wz
        pztemp = photon_pz(ix  ,iy  ,iz  )*(1. - wx)*(1. - wy)*(1. - wz) +
     &           photon_pz(ix+1,iy  ,iz  )*      wx *(1. - wy)*(1. - wz) +
     &           photon_pz(ix  ,iy+1,iz  )*(1. - wx)*      wy *(1. - wz) +
     &           photon_pz(ix+1,iy+1,iz  )*      wx *      wy *(1. - wz) +
     &           photon_pz(ix  ,iy  ,iz+1)*(1. - wx)*(1. - wy)*      wz  +
     &           photon_pz(ix+1,iy  ,iz+1)*      wx *(1. - wy)*      wz  +
     &           photon_pz(ix  ,iy+1,iz+1)*(1. - wx)*      wy *      wz  +
     &           photon_pz(ix+1,iy+1,iz+1)*      wx *      wy *      wz

        ! --- Photon energy in lab frame
        pelab = photon_energy

        ! --- The photon_px, y, and z are only used for directionality.
        ! --- While the momentum encapsulates the energy, the
        ! --- interpolated value may add an artificial numeric spread
        ! --- in the energy. Scale the momentum to be consistent with
        ! --- the input photon energy.
        ! --- Note that the px, y, z are momentum, with units kg m/s.
        ptemp = sqrt(pxtemp**2 + pytemp**2 + pztemp**2)
        pxlab = photon_energy/clight*(pxtemp/ptemp)
        pylab = photon_energy/clight*(pytemp/ptemp)
        pzlab = photon_energy/clight*(pztemp/ptemp)

        ! --- Quantities needed for frame transformations into and out
        ! --- the rest frame of the electron
        if (withgammaeffective) then
c         --- equations from Rykovanov
c         pi0 = gamma_lab - uzlab/clight
c         pz_effective = uzlab/clight + laser_a0sq/(2.*pi0)
c         gamma_effective = sqrt(1. + pz_effective**2)

c         --- Based on Sarachik and Schappert
c         laser_a0sq = echarge**2*planck**2*pdens/(2.*pi**2*eps0*emass**2*clight**2*photon_energy)
c         omega = 2.*pi*photon_energy/planck
c         kz = -omega/clight ! Note that this has the laser moving with vz = -clight
c         eta = omega*tlab - kz*zlab ! For now, use phase = 0
c         ! --- gammae and gammabetae are relative to the electron rest frame before the scatter
c         gammae = 1. + laser_a0sq*sin(eta)**2/2.
c         !gammabetae = -sqrt(gammae**2 - 1.) ! The sign assumes that the laser has vz = -clight
c         gammabetae = -laser_a0sq*sin(eta)**2/2.
c         ! --- Transform to lab frame, using gamma and u of the electron in the lab frame as the boost
c         gamma_effective = gamma_lab*gammae + uzlab/clight*gammabetae

c         --- Based on Bardsley, Penetrante, and Mittleman
c         --- Note that Aamp is the unnormalized vector potential
          laser_Aamp = sqrt(1.*planck**2*pdens/(4.*pi**2*eps0*photon_energy))
          omega = 2.*pi*photon_energy/planck
          kz = -omega/clight ! Note that this has the laser moving with vz = -clight
          eta = omega*tlab - kz*zlab ! For now, use phase = 0
          !call random_number(rn1)
          rn1 = 0.
          !laser_Ax = laser_Aamp*sin(eta + rn1*2.*pi) ! Linearly polarized in x
          !laser_Ay = 0.
          laser_Ax = sqrt(2.)*laser_Aamp*sin(eta + rn1*2.*pi) ! Circularly polarized
          laser_Ay = sqrt(2.)*laser_Aamp*cos(eta + rn1*2.*pi) !
          ! --- p = p0 + eA + zhat*(eA.eA + 2eA.p0)/(2(E0/c - p0||)) eq 2.5
          uxinlaser = uxlab + echarge*laser_Ax/emass
          uyinlaser = uylab + echarge*laser_Ay/emass
          uzinlaser = uzlab - (echarge**2*(laser_Ax**2 + laser_Ay**2)/emass + 2*echarge*(laser_Ax*uxlab + laser_Ay*uylab))
     &                     /(emass*2.*(gamma_lab*clight + uzlab))
          gammainlaser = sqrt(1. + (uxinlaser**2 + uyinlaser**2 + uzinlaser**2)/clight**2)
          gamma_effective = gammainlaser
          gamma_frame = gamma_effective
          uxframe = uxinlaser
          uyframe = uyinlaser
          uzframe = uzinlaser
        else
          gamma_frame = gamma_lab
          uxframe = uxlab
          uyframe = uylab
          uzframe = uzlab
        endif

        umag = sqrt(uxframe**2 + uyframe**2 + uzframe**2)
        if (umag > 0.) then
          uxhat = uxframe/umag
          uyhat = uyframe/umag
          uzhat = uzframe/umag
        endif

        ! --- Transform the photon into the electron rest frame
        if (umag > 0.) then
          pdotv = pxlab*uxhat + pylab*uyhat + pzlab*uzhat
          px = pxlab + (gamma_frame - 1.)*pdotv*uxhat - uxframe*pelab/clight**2
          py = pylab + (gamma_frame - 1.)*pdotv*uyhat - uyframe*pelab/clight**2
          pz = pzlab + (gamma_frame - 1.)*pdotv*uzhat - uzframe*pelab/clight**2
        else
          px = pxlab
          py = pylab
          pz = pzlab
        endif

        ! --- Photon energy in the electron rest frame
        pe = sqrt(px**2 + py**2 + pz**2)*clight
        k = pe/(emass*clight**2) ! scaled energy, relative to electron rest mass

        numer   = 2.0*(2.0 + k*(1.0 + k)*(8.0 + k))
        denom   = k**2*(1.0 + 2.0*k)**2
        factor1 = numer/denom

        numer   = (2.0 + k*(2.0 - k))*log(1.0 + 2.0*k)
        factor2 = numer/k**3

        ! --- sigma is the total cross section
        sigma = e_area*(factor1 - factor2)

        ! --- Probability of scattering
        vz = electrongroup%uzp(ie)*electrongroup%gaminv(ie)
        zswept = boost_gamma*dt*(vz + clight)
        arg = sigma*pdens*zswept
        ! --- For a single scatter event, the exp form is probably not
        ! --- correct. This code would be wrong anyway if arg > 1.
        !prob = 1.0 - exp(-arg)
        prob = arg
        call random_number(rn1)
        call random_number(rn2)

        if (weightingmethod == 1) then
            ! --- Every electron will produce a small number of photons
            nnewphotons = prob*electrongroup%sw(electrons)/photongroup%sw(photons) + (1. - rn1)
            ! --- Adjust the probability so that the appropriate
            ! --- fraction of electrons is scattered.
            !prob = photongroup%sw(photons)/electrongroup%sw(electrons)
        else
            ! --- Only a subset of electrons will produce all of the photons
            if (rn1 < prob) then
              nnewphotons = electrongroup%sw(electrons)/photongroup%sw(photons) + rn2
            else
              nnewphotons = 0
            endif
        endif

        if (nnewphotons == 0) cycle

        nnewphotonstotal = nnewphotonstotal + nnewphotons

        if (nnewphotonstotal > nnewspace) then
          nnewspace = 1.5*nnewspace + nnewphotons
          call chckpart(photongroup,photons,0,nnewspace)

          ! --- Update the electron index in case the data was moved during chckpart
          ie = ie0 + electrongroup%ins(electrons)

        endif

        ! --- rotate photon to align with z axis
        p_perp = sqrt(px**2 + py**2)

        p_tot = sqrt(px**2 + py**2 + pz**2) ! total momentum

        theta = acos(pz/p_tot)

        if (p_perp == 0.0) then
          phi = 0.0
        elseif (py >= 0.0) then
          phi = acos(px/p_perp)
        else
          phi = 2.*pi - acos(px/p_perp)
        endif

        cphi   = cos(phi) ! save rotation angles
        sphi   = sin(phi)
        ctheta = cos(theta)
        stheta = sin(theta)

        ! --- Save old values that are needed later
        ! --- Old means quantities in the electron rest frame before scattering
        px_old = px
        py_old = py
        pz_old = pz
        p_old  = p_tot
        e_old  = pe

        ! --- Get scattering angle from Klein-Nishina cross section
        ! --- The method used is described here:
        ! --- Ozmutl, E. N. "Sampling of Angular Distribution in Compton Scattering"
        ! --- Appl. Radiat. Isot. 43, 6, pp. 713-715 (1992)
        ! --- Int. J. Radiat. Appl. Instrum. Part A

        c0 = 2.0*(2.0*k*k + 2.0*k + 1.0)/(2.0*k + 1.0)**3
        b = (2.0 + c0)/(2.0 - c0)
        a = 2.0*b - 1.0

        do inew=0, nnewphotons-1

          kprime = 0.0

          do while (kprime == 0.0)

            call random_number(rn2)

            x = b - (b + 1.0)*(0.5*c0)**rn2 ! cos(theta)

            factor = 1.0 + k*(1.0 - x)
            denom  = factor*factor
            numer  = 1.0 + x*x + (k*k*(1.0 - x)**2)/factor
            f = numer/denom
            h = a/(b - x)

            ratio = f/h

            call random_number(rn2)

            if (ratio >= rn2) then
              kprime = 1.0/(1.0 - x + 1.0/k)
            endif

          enddo

          e_new = kprime*emass*clight**2 ! energy of scattered photons (in J)
          p_new = e_new/clight

          theta = acos(x) ! polar scattering angle

          call random_number(rn2)

          phi = 2.*pi*rn2 ! azimuthal scattering angle

          ! --- Rotation due to scattering
          psx = p_new*sin(theta)*cos(phi)
          psy = p_new*sin(theta)*sin(phi)
          psz = p_new*cos(theta)

          ! --- Rotate momentum vector back to old frame
          px =  ctheta*cphi*psx - sphi*psy + stheta*cphi*psz
          py =  ctheta*sphi*psx + cphi*psy + stheta*sphi*psz
          pz = -stheta*     psx +            ctheta*     psz

          ! --- Transform back into the lab frame (from the electron rest frame)
          if (umag > 0.) then
            pdotv = px*uxhat + py*uyhat + pz*uzhat
            pxlab = px + (gamma_frame - 1.)*pdotv*uxhat + uxframe*e_new/clight**2
            pylab = py + (gamma_frame - 1.)*pdotv*uyhat + uyframe*e_new/clight**2
            pzlab = pz + (gamma_frame - 1.)*pdotv*uzhat + uzframe*e_new/clight**2
          else
            pxlab = px
            pylab = py
            pzlab = pz
          endif
          pelab = sqrt(pxlab**2 + pylab**2 + pzlab**2)*clight

          ! --- Don't save the gammas at too large of an angle.
          theta = atan2(sqrt(pxlab**2 + pylab**2), pzlab)
          if (theta > anglethreshold) cycle

          ! --- Add the new photon to the end of the array.
          ig = photongroup%ins(photons) + photongroup%nps(photons)
          photongroup%nps(photons) = photongroup%nps(photons) + 1

          ! --- New photon is created at the electron position (in lab frame)
          photongroup%xp(ig) = xlab
          photongroup%yp(ig) = ylab
          photongroup%zp(ig) = zlab

          photongroup%uxp(ig) = pxlab
          photongroup%uyp(ig) = pylab
          photongroup%uzp(ig) = pzlab
          photongroup%gaminv(ig) = pelab/jperev ! A handy place to save the energy, in eV
          if (tbirthpid > 0) then
            photongroup%pid(ig,tbirthpid) = time
          endif

        enddo

        if (rn1 < prob .and. scatterelectron) then

          ! --- For weighting method 1, only a subset of electrons at
          ! --- this point are scattered.
          ! --- For weighting method 2, all electrons at
          ! --- this point are scattered (the if check is always true).

          ! --- Use the parameters from the last photon since it is convenient.
          ! --- Note that if photons%sw != electrons%sw, this algorithm does
          ! --- not exactly conserve energy and momentum.

          ! --- Modify electron, giving it the remainder of the momentum and energy
          delta_px = px_old - px
          delta_py = py_old - py
          delta_pz = pz_old - pz
          delta_e = e_old - e_new
          gamma_new = delta_e/(emass*clight**2) + 1.

          ux = delta_px/emass
          uy = delta_py/emass
          uz = delta_pz/emass

          ! --- Transform back into the lab frame (from the electron rest frame)
          if (umag > 0.) then
            udotv = ux*uxhat + uy*uyhat + uz*uzhat
            uxlab = ux + (gamma_frame - 1.)*udotv*uxhat + uxframe*gamma_new
            uylab = uy + (gamma_frame - 1.)*udotv*uyhat + uyframe*gamma_new
            uzlab = uz + (gamma_frame - 1.)*udotv*uzhat + uzframe*gamma_new
          else
            uxlab = ux
            uylab = uy
            uzlab = uz
          endif
          gamma_lab = sqrt(1. + (uxlab**2 + uylab**2 + uzlab**2)/clight**2)

          ! --- Boost electron back into boosted frame
          ux = uxlab
          uy = uylab
          uz = boost_gamma*uzlab - boost_uz*gamma_lab

          electrongroup%uxp(ie) = ux
          electrongroup%uyp(ie) = uy
          electrongroup%uzp(ie) = uz
          electrongroup%gaminv(ie) = 1./sqrt(1. + (ux**2 + uy**2 + uz**2)/clight**2)

          if (scattercountpid > 0) then
            electrongroup%pid(ie,scattercountpid) = electrongroup%pid(ie,scattercountpid) + 1
          endif

          nscattered = nscattered + 1

        endif

      enddo

      end subroutine compton_scatter_electrons
c=============================================================================
