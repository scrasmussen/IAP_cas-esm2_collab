      subroutine hadvppm(nn,dt,dx,con00,vel,area,areav,con11,
     &                   p1d,t1d,rhoa,mwgt,unitlab)
c  
c----CAMx v5.40 111010
c
c     HADVPPM performs advection using the one-dimensional implementation
c     of the piecewise parabolic method of Colella and Woodward (1984).
c     A piecewise continuous parabola is used as the intepolation polynomial.
c     The slope of the parabola at cell edges is computed from a cumulative
c     function of the advected quantity.  These slopes are further modified
c     so that the interpolation function is monotone.
c
c     This version based on CMAQ HPPM.F, v 1.1.1.1 9/14/98, written by
c     M.T. Odman (10/5/93), NCSC.  This version assumes constant grid cell
c     size.
c     
c     The following definitions are used:  
c  
c              |-----------> Positive direction  
c  
c     |Boundary|<-----------------Domain----------------->|Boundary| 
c  
c     | CON(1) | CON(2) |  ...  | CON(I) |  ...  |CON(N-1)| CON(N) |  
c  
c     VEL(1)-->|     VEL(I-1)-->|        |-->VEL(I)       |-->VEL(N-1)  
c  
c      FP(1)-->|      FP(I-1)-->|        |-->FP(I)        |-->FP(N-1)  
c  
c      FM(2)<--|        FM(I)<--|        |<--FM(I+1)      |<--FM(N)  
c             
c                            -->|   DX   |<-- 
c 
c     Copyright 1996 - 2010
c     ENVIRON International Corporation
c  
c     Modifications:  
c        5/17/00   small modification to flux1,2 to improve mass accounting
c       10/13/03   area weighting applied to winds rather than conc
c       11/04/03   only single FLXARR vector passed back for use in
c                  ZRATES and probling tools
c         9/1/09   Revised area weighting, removed from wind vector
c  
c     Input arguments:  
c        nn                  Number of cells  
c        dt                  Time step (s)  
c        dx                  Length of cell (m)
c        con                 Concentration vector (umol/m3) 
c        vel                 Wind speed vector (m/s)  
c        area                Cell area adjustment vector (1/m2) 
c        areav               Interfacial area adjustment vector (m2) 
c  
c     Output arguments:  
c        con                 Concentration vector (umol/m3)  
c        flxarr              Conc change from interfacial mass flux (umol/m3)
c  
c     Routines called:  
c        none  
c  
c     Called by:  
c        XYADVEC  
c        ZRATES
c 
c      include "camx.prm"
c

      character(len=*) :: unitlab      

      real dx,dt
      integer nn
      real con(nn),vel(nn),area(nn),areav(nn),flxarr(nn)

      real mwgt
      real p1d(nn),t1d(nn),rhoa(nn),trfac(nn)
      real con00(nn),con11(2:nn-1)
c
c-----Local parameters
c     STEEPEN is a flag for discontinuty capturing (steepening)
c     This is disabled in this version as marked by c***
c
c***      logical STEEPEN
c***      parameter (STEEPEN=.false.)
c***      parameter (ETA1=20.0, ETA2=0.05, EPS=0.01)
      parameter (TWO3RDS=2./3.)
c
c***      real fm(MX1D),fp(MX1D),cm(MX1D),cl(MX1D),cr(MX1D),dc(MX1D),
c***     &     c6(MX1D),d2c(MX1D),eta(MX1D),etabar(MX1D),cld(MX1D),crd(MX1D)
c
      real fm(nn)
      real fp(nn)
      real cm(nn)
      real cl(nn)
      real cr(nn)
      real dc(nn)
      real c6(nn)
c
c-----Entry point 
c
      ! conc is advected in umol/m3 or ug/m3

      if(trim(unitlab).eq.'ppbv') then
        do i=1,nn
          trfac(i)=mwgt/( 0.08206*t1d(i)/(p1d(i)/1000.0) )
          con(i)=con00(i)*trfac(i)
        enddo
      elseif(trim(unitlab).eq.'ugom3') then
        trfac=1.0
        con=con00
      else
        stop 'unitlab err in hadvppm.f'
      endif


c-----Set all fluxes to zero. Either positive or negative flux will
c     remain zero depending on the sign of the velocity
c
c***      zeta = dx*dx
      do i = 1,nn
        fm(i) = 0.
        fp(i) = 0.
      enddo
c
c-----Zero order polynomial at the boundary cells
c
      cm(2)  = con(2)
      cm(nn) = con(nn-1)
c
c-----First order polynomial at the next cells, no monotonicity constraint
c     needed
c
      cm(3)    = (con(3) + con(2))/2.
      cm(nn-1) = (con(nn-1) + con(nn-2))/2.
c
c-----Second order polynomial inside the domain
c
      do i = 3,nn-2
c
c-----Compute average slope in the i'th cell
c
        dc(i) = 0.5*(con(i+1) - con(i-1))
c      
c-----Guarantee that CM lies between CON(I) and CON(I+1)
c     monotonicity constraint
      
        if ((con(i+1) - con(i))*(con(i) - con(i-1)).gt.0.) then
          dc(i) = sign(1.,dc(i))*min(
     &                               abs(dc(i)),
     &                               2.*abs(con(i+1) - con(i)),
     &                               2.*abs(con(i) - con(i-1)))
        else
          dc(i) = 0.
        endif
      enddo
c
      do i = 3,nn-3
        cm(i+1) = con(i) + 
     &            0.5*(con(i+1) - con(i)) + (dc(i) - dc(i+1))/6.
      enddo
c
      do i = 2,nn-1
        cr(i) = cm(i+1)
        cl(i) = cm(i)
      enddo
c
c-----Optional discontinuty capturing
c     This is disbaled completely in this version 
c
c***      if (STEEPEN) then
c***        do i = 2,nn-1
c***          eta(i) = 0.
c***          cld(i) = con(i)
c***          crd(i) = con(i)
c***        enddo
c***c 
c***c-----Finite diff approximation to 2nd derivative
c***c
c***        do i = 3,nn-2
c***          d2c(i) = (con(i+1) - 2.*con(i) + con(i-1))/6.
c***        enddo
c***c
c***c-----No discontinuity detection near the boundary: cells 2, 3, NN-2, NN-1
c***c 
c***        do i = 4,nn-3  
c***c 
c***c-----Compute etabars
c***c 
c***          if ((-d2c(i+1)*d2c(i-1).gt.0.) .and.
c***     &        (abs(con(i+1) - con(i-1)) -
c***     &         EPS*min(abs(con(i+1)),abs(con(i-1))).gt.0.)) then
c***            etabar(i) = -zeta*(d2c(i+1) - d2c(i-1))/
c***     &                  (con(i+1) - con(i-1))
c***          else
c***            etabar(i) = 0.
c***          endif
c***          eta(i) = max(0.,min(ETA1*(etabar(i) - ETA2),1.)) 
c***          crd(i) = con(i+1) - 0.5*dc(i+1)
c***          cld(i) = con(i-1) + 0.5*dc(i-1)
c***        enddo
c***c
c***        do i = 2,nn-1
c***          cr(i) = cm(i+1) + eta(i)*(crd(i) - cm(i+1))
c***          cl(i) = cm(i) + eta(i)*(cld(i) - cm(i))
c***        enddo
c***      endif
c
c-----Generate piecewise parabolic distributions
c
      do i = 2,nn-1
c
c-----Monotonicity
c 
        if ((cr(i) - con(i))*(con(i) - cl(i)).gt.0.) then
          dc(i) = cr(i) - cl(i)
          c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))
c
c-----Overshoot cases
c
          if (dc(i)*c6(i) .gt. dc(i)*dc(i)) then
            cl(i) = 3.*con(i) - 2.*cr(i)
          elseif (-dc(i)*dc(i) .gt. dc(i)*c6(i)) then
            cr(i) = 3.*con(i) - 2.*cl(i)
          endif
        else
          cl(i) = con(i)
          cr(i) = con(i)
        endif
        dc(i) = cr(i) - cl(i)
        c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))
      enddo
c
c-----Compute fluxes from the parabolic distribution
c
      do i = 2,nn-1
        x = max(0., -vel(i-1)*dt/dx)
        fm(i) = x*(cl(i) + 0.5*x*(dc(i) + c6(i)*(1. - TWO3RDS*x)))
        x = max(0., vel(i)*dt/dx)
        fp(i) = x*(cr(i) - 0.5*x*(dc(i) - c6(i)*(1. - TWO3RDS*x)))
      enddo
c
c-----Compute fluxes from boundary cells assuming uniform distribution
c
      if (vel(1).gt.0.) then
        x = vel(1)*dt/dx
        fp(1) = x*con(1)
      endif
c
      if (vel(nn-1).lt.0.) then
        x = -vel(nn-1)*dt/dx
        fm(nn) = x*con(nn)
      endif
c
c-----Update concentrations
c
      flxarr(1) = (fp(1) - fm(2))
      do i = 2,nn-1
        flxarr(i) = (fp(i) - fm(i+1))
        con(i) = con(i) - 
     &           area(i)*(areav(i)*flxarr(i) - areav(i-1)*flxarr(i-1))
        con11(i)=con(i)/trfac(i)
      enddo     

c
      return
      end
