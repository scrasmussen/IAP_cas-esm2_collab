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
      subroutine strtmix(nz,csav,kbot,ktop,kfreez,deltat,dx,dy,
     &                   rain,pstar,sigma,zi,ti,pi,qi,alfa)
c
c-----ATDM v1.00
c
c     Calculates the vertical redistribution of a pollutant
c     within stratiform clouds
c
c-----95/06/16:
c     C Emery SAI
c-----95/06/23:
c     -added pstar to wetdep call
c-----95/08/17:
c     -added check for LWET
c     -added kbot to WETDEP argument list
c-----95/10/18:
c     -removed checks on lwet
c-----95/10/20:
c     -fixed coding errors
c-----95/10/30:
c     -test module with wetdep outside of sub-step loop for increased speed
c
c     Mixing in stratiform clouds is accomplished via an explicit K-theory
c     approach using a constant Kc of 50 m2/s. Concentration must be 
c     expressed as a mixing ratio (g/g or ppb, etc), and not as mass per 
c     volume (g/m3 or mol/m3)
c
c     Input arguments:
c               lprnt:          diagnostic debug print flag
c               nz:             number of layers
c               lind:           species index
c               kbot:           layer index of cloud bottom
c               ktop:           layer index of cloud top
c               kfreez:         layer index of freezing level
c               deltat:         time step (min)
c               dx:             cell size in x-direction (km)
c               dy:             cell size in x-direction (km)
c               rain:           rainfall rate (in/hr)
c               pstar:          column value of P* (mb)
c               sigma:          vertical sigma structure
c               zi:             profile of layer interface heights (m)
c               ti:             profile of temperature (K)
c               pi:             profile of pressure (mb)
c               qi:             profile of moisture (kg/kg)
c               wrat:           wet scavenging rate for sulfate (1/s)
c               alfa:           profile of concentrations (ppm)
c               idebug:         diagnostic debug flag
c     Output arguments:
c               alfa:           profile of concentrations (ppm)
c               deltadp:        deposited mass (moles)
c     Routines called:
c               wetdep
c     Called by:
c               rchem     rchemf
c
c
c      integer mxfz
c      parameter(mxfz=10)
      real sigmah(nz),dsigh(nz),dsigi(nz),wt(nz),csav
      dimension ti(nz), pi(nz), zi(nz), qi(nz), sigma(nz+1)
      dimension alfa(nz),wrat(nz)
      real kfact(nz), kc
      data eps /0.622/, rdry /287./, grav /9.8/
      data kc /50./
c
c-----Bit 14 (8192's location) turns on debug writes in strtmix
c
c      idiag = and(idebug,8192)
c      ldiag = .false.
c     if (idiag.ne.0 .and. lprnt .and. 
c    &   (lind.eq.3 .or. lind.eq.13)) ldiag = .true.
c      if (idiag.ne.0 .and. lprnt ) ldiag=.true.
c      if (ldiag) write(*,*)'Entering STRTMIX, deltat = ',deltat
c
c-----Calculate depth parameters
c
      do 10 k = kbot,ktop
        sigmah(k) = (sigma(k) + sigma(k+1))/2.
        dsigh(k) = sigma(k+1) - sigma(k)
   10 continue
      do 11 k = kbot,ktop-1
        dsigi(k) = 1./(sigmah(k+1) - sigmah(k))
   11 continue
c
c-----Calculate density and Kv-factor profile
c
c     write (*,*) 'before 120 in strtmix'
c
      do 120 k = kbot,ktop-1
        pbar = 100.*(pi(k) + pi(k+1))/2.
        tvk = ti(k)*(eps + qi(k))/(eps*(1. + qi(k)))
        tvkp1 = ti(k+1)*(eps + qi(k+1))/(eps*(1. + qi(k+1)))
        tvbar = (tvk + tvkp1)/2. 
        densf = pbar/(rdry*tvbar)
        kfact(k) = csav*kc*densf**2*60.*deltat*dsigi(k)*
     &             (grav/(100.*pstar))**2
c        if (ldiag) write(*,*)'k,densf,kfact',k,densf,kfact(k)
 120  continue
c
c-----Determine maximum time step and number of sub-time steps
c
      dtmax = 60.*deltat
      do 130 k = kbot,ktop-1
        dt = 60.*deltat*dsigi(k)*dsigh(k)**2/(2.*kfact(k))
        dtmax = amin1(dtmax,dt)
 130  continue
      if (dtmax.ge.60.*deltat) then
        ndt = 1
      else
        ndt = int(60.*deltat/dtmax) + 1
      endif
      dtmix = 60.*deltat/float(ndt)
c      if (ldiag) write(*,*)'ndt,dtmix,kbot,ktop',ndt,dtmix,kbot,ktop
c
c-----Start of vertical diffusion; loop over sub time-steps
c
      deltadp = 0.
      do 150 n = 1,ndt
c
        k = kbot
        wt(k) = alfa(k)
     &        + kfact(k)*(alfa(k+1) - alfa(k))/(dsigh(k)*float(ndt))
        wt(k) = amax1(wt(k),0.)
c
c       write (*,*) 'before 160 in strtmix'
c
        do 160 k = kbot+1,ktop-1
          wt(k) = alfa(k)
     &          + (kfact(k)*(alfa(k+1) - alfa(k))
     &          -  kfact(k-1)*(alfa(k) - alfa(k-1)))
     &             /(dsigh(k)*float(ndt))
          wt(k) = amax1(wt(k),0.)
 160    continue
c
        k = ktop
        wt(k) = alfa(k)
     &        - kfact(k-1)*(alfa(k) - alfa(k-1))/(dsigh(k)*float(ndt))
        wt(k) = amax1(wt(k),0.)
c
c-----Load concentration array and calculate wet deposition
c
        do 170 k = kbot,ktop
c          if (ldiag) write(*,*)'k,alfa,wt',k,alfa(k),wt(k)
          alfa(k) = amax1(wt(k),1.E-12)
 170    continue
c
c       call wetdep(ldiag,nz,lind,ktop,kbot,kfreez,dx,dy,dtmix,pstar,
c    &              rain,sigma,zi,ti,pi,wrat,alfa,deldp)
c       deltadp = deltadp + deldp
c
c-----End sub-time loop
c
 150  continue
c
c     dtdep = 60.*deltat
c     write (*,*) 'before wetdep in strtmix'
c     call wetdep(ldiag,nz,lind,ktop,kbot,kfreez,dx,dy,dtdep,pstar,
c    &            rain,sigma,zi,ti,pi,wrat,alfa,deldp,lqdh2o,ix,jx,lpro3)
c     deltadp = deltadp + deldp
c
      return
      end
