c
c     Copyright 1996, 2001 Systems Applications International
c       The REMSAD modeling system software is the property of SAI
c       and is protected by registered copyrights. Distribution of
c       the software outside of your organization is not permitted.
c       Preparation and distribution of any derivative software using
c       any portions of the REMSAD code, without the express permission
c       of SAI, is also not permitted. By accepting the REMSAD code
c       provided in this release, you are also accepting the
c       limitations, terms, and provisions of use as noted above.
      subroutine convmix(nz,cd,kbot,ktop,kfreez,deltat,dx,dy,
     &                pstar,xmd,rain,sigma,zi,ti,pi,alfa)
c
c-----ATDM v1.00
c
c     Calculates the vertical redistribution of a pollutant
c     within convective clouds
c
ctcm---05/06/14 - revised action on negative concentration to avoid
c                 artificially creating mass
c
cabh---00/06/20 - ABH changed comparison with 10-12 back to 10-6
cabh---00/06/21 - ABH changed comparison with 10-06 back to 10-12
c-----95/06/16:
c     C Emery SAI
c-----95/10/18:
c     -removed checks for lwet
c-----95/10/20:
c     -restructured sub-step loop for more efficiency
c-----95/10/30:
c     -test module with wetdep outside of sub-step loop for increased speed
c-----99/05/05: DFG
c     -removed the block for determining each sub-time step. This is because
c      the method is not stable and negative concentration occurs during updating.
c     -changed the explicit method to implicit iteration method. The error
c      tolerance is set to be 0.0005.  Maximum iteration number is 5 times.
c      The sub-time step is initial set to equal 15 minutes.  If it iterates
c      5 times and still not converge, the sub-time step will be cut half,
c      and so on.
c
c     Mixing in convective clouds follows the procedure outlined by
c     Lin et al. (JGR, 1994; p. 25615-25630). Concentration must be expressed
c     as a mixing ratio (g/g or ppb, etc), and not as mass per volume 
c     (g/m3 or mol/m3)
c
c     Input arguments:
c               lprnt:          diagnostic debug print flag
c               nz:             number of layers   
c               lind:           species index
c               kbot:           layer index of cloud bottom
c               ktop:           layer index of cloud top
c               kfreez:         layer index of freezing level 
c               deltat:         time step (s)
c               dx:             cell size in x-direction (km)
c               dy:             cell size in y-direction (km)
c               pstar:          column P* (mb)
c               xmd:            convective mixing coefficient (1/s)
c               rain:           rainfall rate (in/hr)
c               sigma:          vertical sigma structure 
c               zi:             profile of layer interface heights (m)
c               ti:             profile of temperature (K)
c               pi:             profile of pressure (mb)
c               wrat:           wet scavenging rate for sulfate (1/s)
c               alfa:           profile of concentrations (ppm)
c               idebug:         diagnostic debug flag
c               ierr:           error flag
c     Output arguments:
c               alfa:           profile of concentrations (ppm)
c               deltadp:        deposited mass (moles)
c     Files written:
c               6:              standard output file unit (*)         
c     Routines called:
c               wetdep
c     Called by:
c               rchem     rchemf
c
c      integer mxfz
c      parameter(mxfz=10)
      real   cd
      dimension ti(nz), pi(nz), zi(nz), sigma(nz+1)
      dimension alfa(nz),dalfdt(nz)
      dimension hold(nz)
c
c-----Bit 14 (8192's location) turns on debug writes in convmix
c
c      idiag = and(idebug,8192)
c      ldiag = .false.
c      if (idiag.ne.0 .and. lprnt .and. lind.eq.15) ldiag = .true.
c    &   (lind.eq.3 .or. lind.eq.13)) ldiag = .true.
c      if (idiag.ne.0 .and. lprnt )  ldiag = .true.
c      ldiag=.true.
c      if (ldiag) write(*,*)'Entering CONVMIX, deltat = ',deltat
c      ierr = 0
c
c-----First estimate some vertically integrated quantities
c
      xtol=0.02
      nitmax=15
      nstep = 1 
      timrem = deltat*60.
c
c set sub cloud mixing scale to 15 minutes
      dtmax=timrem
      dttot=0.
      dtmix=dtmax/float(nitmax)
      dt=dtmax
    5 continue
c
c   while setting up initial values for next step, also save max in
c   any layer and total mass in column.
c
      alfmax = 0.
      alfmass = 0.
      do k=1,nz
        hold(k)=alfa(k)
        alfmax = amax1(alfmax, alfa(k))
        alfmass = alfmass + alfa(k) * (sigma(k) - sigma(k+1))
      end do
c
      DO 1000 iter=1,nitmax
      alfasc = 0.
      do 10 k = 1,nz
        dalfdt(k) = 0.
        if (k.lt.kbot) then
          dsigma = sigma(k) - sigma(k+1)
          alfasc = alfasc + hold(k)*dsigma
        endif
   10 continue
      alfasc = alfasc/(1. - sigma(kbot))
c      if (ldiag) write(*,*)'alfasc ',alfasc
c
c-----Estimate rate of change within cloud layers
c
      dadtsum = 0.
      do 20 k = kbot,ktop
        dsigma = sigma(k) - sigma(k+1)
        dalfdt(k) = cd*xmd*(alfasc - hold(k))/4.
        dadtsum = dadtsum + dalfdt(k)*dsigma
c        if (ldiag) write(*,*)'dalfdt(k) ',k,dalfdt(k)
   20 continue
c
c-----Estimate rate of change within sub-cloud layers
c
      sumasc = 0.
      do 30 k = 1,kbot-1
        dsigma = sigma(k) - sigma(k+1)
        ak = hold(k) 
c-------------original-------------------------
c        if (dadtsum.lt.0.) ak = 1./hold(k)
c---------------original-----------------------
c++++++++++++++++lijie add+++++++++++++++++++++
         if(dadtsum.lt.0.) then
            dalfdt(k)=0.0
            goto 30
         endif
c++++++++++++++++finished++++++++++++++++++++++
        sumasc = sumasc + ak*dsigma
   30 continue
      dadtscs = 0.
      do 40 k = 1,kbot-1
        dsigma = sigma(k) - sigma(k+1)
        ak = hold(k) 
c--------------------original---------------------
c        if (dadtsum.lt.0.) ak = 1./hold(k)
c--------------------original---------------------
c+++++++++++++++++lijie add+++++++++++++++++++++++
           if (dadtsum.lt.0.) then
              goto 400
           endif  
c+++++++++++++++++++finished++++++++++++++++++++++
        dalfdt(k) = -ak*dadtsum/sumasc
  400   dadtscs = dadtscs + dalfdt(k)*dsigma
   40 continue
c+++++++++++++++++lijie add+++++++++++++++++++++++
      if(dadtsum.lt.0.) then
         
         dalfdt(ktop)=abs(dadtsum)/(sigma(ktop-1) - sigma(ktop))
c         print*,ktop,dalfdt(3),dalfdt(ktop),dadtsum
      endif 
c++++++++++++++++finished+++++++++++++++++++++++++
c
c-----Apply rates to each layer of the column
c
      tmpmax=0.
      do 60 k = 1,nz
        alfnew=hold(k)
cabh    if (hold(k) .gt. 1.0e-6) then
ctcm    if (hold(k) .gt. 1.0e-12) then
ctcm    if (hold(k) .gt. 1.0e-30) then   ! try an even smaller value
        if (hold(k) .gt. 1.0e-12) then   
           err=abs(dalfdt(k)/hold(k))
        else
           err=abs(dalfdt(k))
        endif
          
        hold(k)=alfnew + dtmix*dalfdt(k)
        if (err.gt.tmpmax) then
          tmpmax=err
        endif
c       if (hold(k).le.0.) hold(k)=alfnew/10.
        if (hold(k).le.-1.0e-12) tmpmax = 10.*xtol
   60 continue
                       
      if (tmpmax.le. xtol) then
        dttot=dttot+dtmix
      if (dttot.ge.timrem*0.99999) go to 100
      else 
c        dtmix=dtmix/2.0
c        nitmax=2*nitmax
        goto 200
      endif
  
 1000 CONTINUE           
 100  continue
c
c  don't let concentration go below lower bound or above the original
c  maximum in any layer
c
        alftest = 0.
        do k=1,nz
          alfa(k) = amin1(hold(k),alfmax)
          alftest = alftest + alfa(k) * (sigma(k) - sigma(k+1))
        end do

c
c  check mass and if it has changed, scale all concs.
c
        if ( abs(alftest - alfmass)/alfmass .gt. 0.02) then
c       if ( ix .eq. 79 .and. jx .eq.65) then
          alffac = alfmass/alftest
c         write(*,*) 'CONVMIX: Mass check at ', ix, jx, ' lind = ',lind
c         write(*,*) 'Initial mass: ',alfmass,' Final mass: ',alftest
c         write(*,*) 'Alfa before adjustment: '
          do k=1,nz
c           write(*,*) k, alfa(k)
            alfa(k) = alfa(k) * alffac
          enddo
        endif
      return
  200 continue
      end
