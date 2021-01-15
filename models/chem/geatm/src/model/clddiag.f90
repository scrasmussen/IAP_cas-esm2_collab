      subroutine clddiag(nz,lconpcp,pbl,rainc,dz_1d,zm_1d,pr,ta,qv,qc,ccov)
!
!-----CLDDIAG diagnoses un-resolved cloud fields for the purpose of calculating
!     cloud optical depth (for photolysis rates adjustment) and water contents
!     (for chemistry and wet deposition).
!
!     This code uses some aspects of cloud diagnostic routines in MCIP v3.3
!     and CMAQ v4.7.
!
!     Input:
!        nz       number of layers
!        lconpcp  Flag denoting sub-grid convective precip
!        pbl      PBL height (m)
!        rainc    convective rainfall rate (mm/hr)
!        zh       Layer interface heights (m)
!        pr       Layer pressure (Pa)
!        ta       Layer temperature (K)
!        qv       Layer humidity (kg/kg)
!        qc       Layer cloud water content (g/m3)
!     
!     Output:  
!        qc       Layer cloud water content (g/m3)
!        ccov     Layer cloud cover fraction
!
      implicit none
!
!-----Argument variables
!
      integer nz
      real pbl,rainc
      real zh(nz),pr(nz),ta(nz),qv(nz),qc(nz),ccov(nz)
      logical lconpcp

      real dz_1d(nz),zm_1d(nz)
!
!-----Local variables
!
      integer k,kmx,ksrc,istab,icount,kbot,ktop
      real rd,cp,lv,lf,mh2o,mair,evap0,rv,grav,tpert,qpert,rhc1,rhc2, &
           rdocp,lvocp,mvoma,gamma,evap
      real esat,qsat,rhc,sg1,tsrc,qsrc,psrc,theta,esrc,pp,tp,qp,etheta, &
           tdsrc,tdlps,delta,dzt,plcl,tlcl,tbar,dp,pbar,tbase,x1,dtdp, &
           frac,qxs,wl,pavg,tempa,tempb,tempc,qent,fent,tent,ti,dql,dqi, &
           fa,fb,htst,t1,qvc,prate,rhoair
      real zm(nz),dz(nz),rh(nz),tsat(nz),qad(nz),qwat(nz), &
           tcld(nz),twc(nz)

      data rd    /287./    !J/K/kg
      data cp    /1004./   !J/K/kg
      data lv    /2.5e6/   !J/kg
      data lf    /3.3e5/   !J/kg
      data mh2o  /18./     !g(h2o)/mol
      data mair  /28.8/    !g(air)/mol
      data evap0 /611./    !Pa
      data rv    /461./    !J/K/kg
      data grav  /9.8/     !m/s2
      data tpert /1.5/     !K
      data qpert /1.5e-3/  !kg/kg
      data rhc1  /0.7/
      data rhc2  /0.9/
!     
!-----Initialize cloud variables
!
      rdocp = rd/cp
      lvocp = lv/cp
      mvoma = mh2o/mair
      gamma = grav/cp
 
      kbot = 0
      ktop = 0
      do k = 1,nz
        ccov(k) = 0.
        twc(k) = 0.

! comment out by shun
!---------------
!        zm(k) = zh(k)/2.
!        dz(k) = zh(k)
!        if (k .gt. 1) then
!          zm(k) = (zh(k) + zh(k-1))/2.
!          dz(k) = zh(k) - zh(k-1)
!        endif
!----------------
        zm(k)=zm_1d(k)
        dz(k)=dz_1d(k)

      enddo
!
!-----Determine unresolved "stratus" cloud coverage based on RH profile
!
!     kmx = 1
      do k = 1,nz
!
!-----Calculate RH for this layer
!
        esat  = evap0*exp((lv/rv)*(1./273. - 1./ta(k)))
        qsat  = mvoma/(pr(k)/esat - 1.)
        rh(k) = min(qv(k)/qsat,1.)
!
!-----Set critical RH to 98% in PBL
!
!       if (zm(k) .lt. pbl) then
!         rhc  = 0.98
!         kmx = k
!         if (rh(k) .gt. rhc) ccov(k) = 0.34*(rh(k) - rhc)/(1. - rhc)
!       else
!
!-----Calculate critical RH profile above PBL
!
!         sg1 = pr(k)/pr(kmx)
!         rhc = 1. - (2.*sg1*(1. - sg1)*(1. + 1.732*(sg1 - 0.5)))
!         rhc = min(rhc,0.98)
!         if (rh(k) .gt. rhc) ccov(k) = ((rh(k) - rhc)/(1. - rhc))**2
!       endif
!       ccov(k) = max(min(ccov(k),1.),0.)
!
!-----If a resolved cloud exists here, use it.  Otherwise, set cloud water
!     for any un-resolved cloud coverage > 10%
!
        if (qc(k) .gt. 0.01) then           ! A resolved cloud exists here
          ccov(k) = 1.
        else
!         if (ccov(k) .gt. 0.1) then
!           qc(k) = 0.05e3*qv(k)*pr(k)/rd/ta(k)
!         else
            qc(k) = 0.
            ccov(k) = 0.
!         endif
        endif
      enddo
!     if (scmeth.eq.'GRELL') goto 999
!
!-----Start diagnosis of sub-grid convection
!
!-----Determine cloud source level by determining max equivalent potential
!     temperature between surface and 3000 m
!
      ksrc = 1
      tsrc = ta(1) + tpert
      qsrc = qv(1) + qpert
      psrc = pr(1)
      theta = tsrc*(1.e5/psrc)**rdocp
      esrc = theta*exp(lvocp*qsrc/tsrc)

      do k = 2,nz
        pp = pr(k)
        if (zm(k) .gt. 3000.) goto 10
        tp = ta(k) + tpert
        qp = qv(k) + qpert
        theta = tp*(1.e5/pp)**rdocp
        etheta = theta*exp(lvocp*qp/tp)
        if (etheta .gt. esrc) then
          ksrc = k
          tsrc = tp
          qsrc = qp
          psrc = pp
          esrc = etheta
        endif
      enddo
!
!-----If no covective precip and RH at source layer < 70%,
!     exit convective diagnosis
!
 10   if (.not.lconpcp .and. rh(ksrc) .le. rhc1) goto 999
!
!-----Compute lifting condensation level:
!     First, determine dew point lapse rate
!
      evap = qsrc*psrc/(mvoma + qsrc)
      tdsrc = 1./(1./273. - (rv/lv)*alog(evap/evap0))
      tdsrc = amin1(tdsrc,tsrc)
      tdlps = (grav*tdsrc*tdsrc)/(mvoma*lv*tsrc)
!
!-----Second, compute difference between dry adiabatic and dew point lapse
!     rate, height increment above source level to reach LCL, and pressure
!     at LCL
!
      delta = gamma - tdlps
      if (delta .le. 0.) then
        dzt = 0.
        plcl = psrc
        tlcl = tsrc
      else
        dzt = (tsrc - tdsrc)/delta
        tlcl = tsrc - gamma*dzt
        tbar = 0.5*(tsrc + tlcl)
        plcl = psrc*exp(-(grav/rd)*dzt/tbar)
      endif
!
!-----Determine layer containing cloud base
!
      do k = 2,nz
        if (pr(k) .le. plcl) then
          kbot = k
          goto 20
        endif
      enddo
!
!-----Cloud base layer never found, exit convective diagnosis
!
 20   if (kbot.eq.0) goto 999
!
!-----If no convective precip and base > PBL, exit convective diagnosis
!
      if (.not.lconpcp .and. zm(kbot) .gt. pbl) then
        kbot = 0
        goto 999
      endif
!
!-----Determine cloud top by following moist adiabat up from the base. 
!
      istab = 0
      do k = kbot,nz
        dp = pr(k-1) - pr(k)
        pbar = 0.5*(pr(k-1) + pr(k))
        if (k.eq.kbot) then
          dp = plcl - pr(k)
          pbar = 0.5*(plcl + pr(k))
          tbase = tlcl
        endif

        tbar = tbase - 6.5e-4*dp/2.
        esat = evap0*exp((lv/rv)*(1./273. - 1./tbar))
        qsat = mvoma/(pbar/esat - 1.)
        x1 = lv*qsat/(rd*tbar)
        dtdp = (rd*tbar)/(pbar*cp)*((1. + x1)     &
               /(1. + (mvoma*lvocp/tbar)*x1))
        tsat(k) = tbase - dp*dtdp
        esat = evap0*exp((lv/rv)*(1./273. - 1./tsat(k)))
        qad(k) = mvoma/(pr(k)/esat - 1.)
        tbase = tsat(k)

        if (istab.eq.0) then
          if (tsat(k).gt.ta(k)) then
            istab = 1
          endif
        else
          if (tsat(k).lt.ta(k)) then
            ktop = k - 1
            goto 30
          endif
        endif
      enddo
      ktop = nz - 1
!
!-----If ISTAB is still = 0, we have a "stable" cloud: exit convective diagnosis
!
 30   if (istab.eq.0) then
        kbot = 0
        ktop = 0
        goto 999
      endif
!
!-----Determine cloud coverage for non-precipitating convection
!
      if (.not.lconpcp) then
        frac = 0.5*(rh(ksrc) - rhc1)/(rhc2 - rhc1)
        if (frac.le.0.1) then
          kbot = 0
          ktop = 0
          goto 999
        endif
        do k = kbot,ktop
          ccov(k) = max(ccov(k),frac)
        enddo
      endif
!
!-----Determine entrainment factor using iterative method
!
      qxs = 0.
      do k = kbot,ktop
        wl = 0.7*exp((pr(k) - plcl)*0.000125) + 0.2
        if (k .eq. kbot) then
          pavg = (pr(k) + pr(k-1))/2.
          if (plcl .lt. pavg) then
            pavg = (pr(k+1) + pr(k))/2.
            pavg = (pavg + plcl)/2.
            wl = 0.7*exp((pavg - plcl)*0.000125) + 0.2
          endif
        endif

        qwat(k) = max(wl*(qsrc - qad(k)),1.e-20)
        tempa = tsat(k) - 20.
        tempb = tsat(k) + 10.

        qent = qv(k)
        x1   = max((qent - qsrc),1.e-10)
        esat = evap0*exp((lv/rv)*(1./273. - 1./tempa))
        qsat = mvoma/(pr(k)/esat - 1.)
        fent = (qsat + qwat(k) - qsrc)/x1
        fent = max(min(fent,1.),0.)

        tent = ta(k)
        ti = tsat(k)*(1. - fent) + tent*fent
        dql = (qsrc - qad(k))*(1. - fent - wl)
        dqi = 0.
        if (tempa .lt. 273.) then
          dqi = -qwat(k)*(tempa - 273.)/18.
          if (tempa .le. 255.) dqi = qwat(k)
        endif
        fa = cp*(tempa - ti) + lv*dql + lf*dqi
!
!-----Top of iterative loop: cut the temperature interval in half
!
        icount = 0
 599    continue
        htst = tempb - tempa
        if (htst .lt. 0.01) goto 595
        icount = icount + 1
        if (icount .gt. 100) then
          write(*,*) 'No convergence in entrainment solver!'
          stop
        endif

        tempc = (tempa + tempb)/2.
        qent  = qv(k)
        x1    = max((qent - qsrc),1.e-10)
        esat  = evap0*exp((lv/rv)*(1./273. - 1./tempc))
        qsat  = mvoma/(pr(k)/esat - 1.)
        fent  = (qsat + qwat(k) - qsrc)/x1
        fent  = max(min(fent,0.99),0.01)

        tent = ta(k)
        ti = tsat(k)*(1. - fent) + tent*fent
        dql = (qsrc - qad(k))*(1. - fent - wl)
        dqi = 0.
        if (tempc .lt. 273.) then
          dqi = -qwat(k)*(tempc - 273.)/18.
          if (tempc .le. 255.) dqi = qwat(k)
        endif
        fb = cp*(tempc - ti) + lv*dql + lf*dqi

        if (fa*fb) 590, 590, 591
590     tempb = tempc
        go to 599
591     tempa = tempc
        fa = fb
        go to 599
595     continue                  ! exit from iterator, convergence achieved
        tcld(k) = max(tempc,150.)
!
!-----Excess water for rainout (kg/m2)
!
        t1 = max(tcld(k),ta(k))
        esat  = evap0*exp((lv/rv)*(1./273. - 1./t1))
        qvc = mvoma/(pr(k)/esat - 1.)
        rhoair = pr(k)/(rd*tcld(k))
        qxs = qxs + (qwat(k) + qvc - qv(k))*rhoair*dz(k)
      enddo
!
!-----Determine cloud coverage for precipitating convection.
!     PRATE is storm rainout rate in mm/hour, noting that 1 kg
!     of water occupies a 1 mm thick layer of water in a square meter
!     of ground (accounts for density of water = 1000 kg/m3)
!
      if (lconpcp) then
        prate = 0.3*qxs
        frac = 1.
        if (prate .gt. 1.001*rainc) frac = rainc/prate
        if (frac .le. 0.1) then
          ktop = 0
          kbot = 0
          goto 999
        endif
        do k = kbot,ktop
          ccov(k) = max(ccov(k),frac)
        enddo
      endif
      do k = kbot,ktop
        rhoair = pr(k)/(rd*tcld(k))
        twc(k) = 1.e3*qwat(k)*rhoair
      enddo
!
 999  do k = 1,nz
        qc(k) = max(qc(k),twc(k))
      enddo
!
      return
      end
