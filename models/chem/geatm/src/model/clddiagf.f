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
      subroutine clddiagf(nx,ny,nz,lwet,delfx,delfy,uf,vf,htf,presf,
     &                    temprf,rainf,rh0,sigmaf,pstarf,psfcf,
     &                    lati,cloudf,xmdf,kcldbf,kcldtf,kfrzf,
     &                    fconvf,fstrtf,ffnonpr,dt,conc,iwb,ieb,jsb,jeb)
c
c-----ATDM v1.00
c
c     Diagnoses cloud parameters (top, bottom, freezing level, cloud cover)
c     and type (stratiform vs. convective).  Diagnosis based on input 
c     temperature, moisture, and pressure profiles.  Convection triggered if 
c     moisture convergence is above a threshold value and sounding through
c     a column is convectively unstable; following MM5 (Anthes et al, NCAR
c     Tech Note TN-282+STR, 1987).
c
c     Fine grid version
c
c-----95/05/19:
c     C Emery SAI
c-----95/10/04:
c     -put max limit on CS = 1.0
c     -put limit on CS such that if CS < 0.1, CS = 0.0
c-----95/10/18:
c     -restructured for more efficient handling of lwet and rain
c-----95/10/20:
c     -moved calculation of xmd to just after check of convective instability;
c      now checks for xmd > 0 for convection to occur
c     -checks for rain > 0.01 in/hr (instead of 0) to diagnose rain-clouds
c
c-----95/12/26:
c     -put min limit on kbot (>=3) for convective clouds at about line 199
c      (quick fix to try and prevent oscilation for strong convective
c       columns with low bottom)
c
c-----96/09/26:
c     -fixed argument type mismatch in scavrat call -- cd (real) passed 
c      as logical variable, fixed by defining lcd=.true. when cd .gt. 0.
c
c     Input arguments:
c               nx:             number of fine grid cells in x-direction     
c               ny:             number of fine grid cells in y-direction     
c               nz:             number of fine grid layers
c               lwet:           wet deposition flag
c               delfx:          fine grid cell size in x-direction (km)
c               delfy:          fine grid cell size in y-direction (km)
c               uf:             array of u-components winds (km/min)
c               vf:             array of v-components winds (km/min)
c               htf:            array of layer interface heights (m)      
c               presf:          array of pressure (mb)
c               temprf:         array of temperature (K)
c               rainf:          array of rain(inches/hr)
c               ch2of:          array of water vapor (kg/kg)
c               sigmaf:         fine grid sigma structure
c               pstarf:         array of psfc-ptop (mb)
c               psfcf:          array of surface pressure (mb)
c               tphioaf:        latitude convergence factor (1/km)
c               idebug:         diagnostic debug flag
c     Output arguments:
c               cloudf:         fraction of stratiform clouds
c               xmdf:           convective mixing coefficient (1/s)
c               kcldbf:         layer index of cloud base
c               kcldtf:         layer index of cloud top
c               kfrzf:          layer index of freezing level
c               fconvf:         fraction of precipitating convective coverage
c               fstrtf:         fraction of precipitating stratiform coverage
c               wratf:          sulfate scavenging rate (1/s)
c     Files written:
c               6:              standard output file unit (*)
c     Routines called:
c               scavrat
c     Called by:
c               update
c
c      integer :: mxfz
c      parameter (mxfz=10) !need to be modified to meet real layers
      integer :: iwb,ieb,jsb,jeb
      real :: dt 
      dimension uf(iwb:ieb,jsb:jeb,nz),vf(iwb:ieb,jsb:jeb,nz),
     &           htf(iwb:ieb,jsb:jeb,nz),presf(iwb:ieb,jsb:jeb,nz),
     &          temprf(iwb:ieb,jsb:jeb,nz),rainf(iwb:ieb,jsb:jeb),
     &          ch2of(iwb:ieb,jsb:jeb,nz),sigmaf(nz+1),
     &          pstarf(iwb:ieb,jsb:jeb),
     &          psfcf(iwb:ieb,jsb:jeb),tphioaf(iwb:ieb,jsb:jeb),
     &          rh0(iwb:ieb,jsb:jeb,nz) 
      dimension cloudf(iwb:ieb,jsb:jeb,nz),xmdf(iwb:ieb,jsb:jeb),
     &          kcldbf(iwb:ieb,jsb:jeb),
     &          kcldtf(iwb:ieb,jsb:jeb),kfrzf(iwb:ieb,jsb:jeb),
     &          fconvf(iwb:ieb,jsb:jeb),fstrtf(iwb:ieb,jsb:jeb)
      dimension delfx(iwb:ieb,jsb:jeb,nz),delfy(iwb:ieb,jsb:jeb,nz)
      real      lati(iwb:ieb,jsb:jeb),ffnonpr(iwb:ieb,jsb:jeb,nz)
      real lv
      dimension ti(nz), pi(nz), zi(nz),rhi(nz)
      dimension td(nz), theta(nz), qs(nz), rh(nz), rhc(nz)
      dimension qi(nz), thetae(nz), cs(nz), rate(nz)
      logical lwet, lprnt,lcd
      real conc(iwb:ieb,jsb:jeb,nz),alfa(nz)
      data p0 /1000./, gamma /0.286/, e0 /6.11/, eps /0.622/
      data lv /2.5e6/, rv /461./, cp /1004./, g /9.8/, rd /287./
      lprnt = .false.
      gamadry = g/cp

c     write(*,*)'Enter CLDDIAGF '
      do 100 j = jsb+1,jeb-1
        do 200 i = iwb+1,ieb-1
c
c-----Load variables for column
c
         tphioaf(i,j)=tan(lati(i,j)*0.01745329)/6367.
          lcd=.false.
          cd = 0.
          csav = 0.
          xmdf(i,j) = 0.
          ktop = 0
          kbot = 0
          kfreez = 0 
          do 10 k = 1,nz
            
            cs(k) = 0.
            zi(k) = htf(i,j,k)
            pi(k) = presf(i,j,k)
            ti(k) = temprf(i,j,k)
            rhi(k) =rh0(i,j,k)
            es=e0*exp((lv/rv)*(1./273. - 1./ti(k)))
            ch2of(i,j,k)=eps/(pi(k)/(rhi(k)*es/100.)-1)
            qi(k) = ch2of(i,j,k)
            alfa(k)=conc(i,j,k)
   10     continue
c         
c-----Loop over layers to find pressure measures at layer interface 
c     heights, and moisture parameters in layers
c
c         write(*,1001)
          k700 = 0
          rhav = 0.
          do 20 k = 1,nz
            if (pi(k).lt.500. .and. k700.eq.0) then
              k700 = k - 1
              if (k700.lt.1) k700 = 1
            endif
            theta(k) = ti(k)*(p0/pi(k))**gamma

            ev = qi(k)*pi(k)/(qi(k) + eps)
            tdk = 1./273 - (rv/lv)*alog(ev/e0)
            td(k) = 1./tdk
            es = e0*exp((lv/rv)*(1./273. - 1./ti(k)))
            qs(k) = (eps*es)/(pi(k) - es)
            rh(k) = 100.*ev/es
         

            rhav = rhav + rh(k)*(sigmaf(k) - sigmaf(k+1))
            prat = pi(k)/psfcf(i,j)
            term1 = 100.*2.*prat*(1. - prat)
            term2 = 1. + sqrt(3.)*(prat - 0.5)
            rhc(k) = 100. - term1*term2

c           write(*,1000) sigmaf(k),sigmaf(k+1),pi(k),ti(k),theta(k),
c    &                    qi(k),qs(k),td(k),rh(k),rhc(k)
   20     continue
 1000     format(5x,2f6.2,3f6.1,2f7.4,f6.1,2f5.0)
 1001     format('        d(sigma)    p     T    TH    qi     qs',
     &           '     Td      RH    RHc')
c         write(*,*)'Column mean RH: ',rhav
          if (rhav.lt.50.) rhav = 50.
          bfac = 2.*(1. - rhav/100.)

c
c-----Loop over levels and determine the extent of stratiform clouds
c
c         write(*,*)'Layer   Cloud Fraction'
          csmax = 0.
          if (ti(1).le.273.) kfreez = 1
          do 30 k = 1,nz
            if (rh(k).gt.rhc(k)) then
              cs(k) = ((rh(k) - rhc(k))/(100. - rhc(k)))**2 
              cs(k) = amin1(cs(k),1.)
              if (cs(k).lt.0.1) cs(k) = 0.
            endif
            csmax = amax1(cs(k),csmax)
            if (kfreez.eq.0 .and. ti(k).lt.270.) kfreez = k
c           write(*,1002) k,100.*cs(k)
   30     continue
 1002     format(i5,7x,f10.1)

c
c-----If it is raining in this column, determine whether cloud type is
c     convective or stratiform
c
c          if (.not.lwet .or. rainf(i,j).lt.0.01) goto 75
          if (.not.lwet) goto 75
c
c-----Loop over layers between surface and 700 mb
c     find layer having maximum equivalent potential temperature, and save
c     LCL (cloud base), T, and p at that level
c
          themax = 0.
          do 40 k = 1,k700
            gamaqs = (g*td(k)**2)/(eps*lv*ti(k))
            zs = (ti(k)-td(k))/(gamadry-gamaqs)
            tmps = ti(k) - gamadry*zs
            tbar = (tmps + ti(k))/2.
            dlnp = (g*zs)/(rd*tbar)
            prss = pi(k)*exp(-dlnp)
            thetae(k) = theta(k)*exp(lv*qi(k)/(cp*tmps))
            if (thetae(k).gt.themax) then
              if (k.eq.1) then
                zlcl = zi(k)/2. + zs
              else
                zlcl = zi(k-1) + (zi(k) - zi(k-1))/2. + zs
              endif
              tlcl = tmps
              plcl = prss
              themax = thetae(k)
            endif
   40     continue
c         write(*,*)'LCL at temperature: ',tlcl
c         write(*,*)'          pressure: ',plcl
c         write(*,*)'            height: ',zlcl
c         write(*,*)'           theta_e: ',themax
c
c-----Loop over layers from LCL to top; find saturated equivelent potential
c     temperature for each level, and convective cloud top, bottom, freezing 
c     level
c
          ipos = 0
          sum = 0.
          do 50 k = 1,nz
cgem            if (k.eq.1 .or. pi(k).gt.plcl) then
            if (k.le.2 .or. pi(k).gt.plcl) then
              kbot = k+1
              goto 50
            endif
            dsigma = sigmaf(k) - sigmaf(k+1)
            thetae(k) = theta(k)*exp(lv*qs(k)/(cp*ti(k)))
            if (ipos.eq.0 .and. themax.gt.thetae(k)) ipos = 1
            if (ipos.eq.1 .and. themax.le.thetae(k)) then
              ktop = k
              goto 51
            endif
            sum = sum + (themax - thetae(k))*dsigma
c           write(*,*)'themax,thetae,tdif,dsigma,sum',
c    &                 themax,thetae(k),themax-thetae(k),dsigma,sum
   50     continue
          if (ipos.eq.1) then
            ktop = nz
          else
c
c-----No convective clouds; no positive energy
c
c           write(*,*)'Sounding is convectively stable'   
c           write(*,*)'No positive area: energy = ',sum
            ktop = 0
            kbot = 0
            goto 61
          endif
   51     if (sum.lt.0) then
c           write(*,*)'Sounding is convectively stable'
c           write(*,*)'energy = ',sum
            ktop = 0
            kbot = 0
            goto 61
          else         
c
c-----Convective clouds possible; positive energy
c
            thsum = 0.
            thdfsum = 0.
            do 65 k = kbot,ktop
              dsigma = sigmaf(k) - sigmaf(k+1)
              thsum = thsum + thetae(k)*dsigma
              thdfsum = thdfsum + (thetae(kbot) - thetae(k))*dsigma
   65       continue
            thtebar = thsum/(sigmaf(kbot) - sigmaf(ktop+1))
            xmdf(i,j) = 1.e-3*(themax - thtebar)/(thdfsum+1.0e-10)
            if (xmdf(i,j) .gt. 1.0) xmdf(i,j) = 1.0
            if (xmdf(i,j).lt.0.) then
c             write(*,*)'Convection not strong enough'
c             write(*,*)'xmd < 0: ',i,j,xmd(i,j)
c             write(*,*)'themax,thtebar,thdfsum',
c    &                   themax,thtebar,thdfsum
              xmdf(i,j) = 0.
              ktop = 0
              kbot = 0
              goto 61
            endif
            delsig = sigmaf(kbot) - sigmaf(ktop+1)
            if (delsig.lt.0.3) then
c             write(*,*)'Convection not deep enough'
c             write(*,*)'energy = ',sum
              ktop = 0
              kbot = 0
              goto 61
            else
c             write(*,*)'Sounding is convective!'
c             write(*,*)'energy = ',sum
c             write(*,*)'Top at   : ',ktop,htop,sigmaf(ktop+1)
c             write(*,*)'Bottom at: ',kbot,hbot,sigmaf(kbot)
c             write(*,*)'Freezing : ',kfreez,htfreez
            endif
          endif
c
c-----If sounding is convectively unstable, calculate moisture convergence
c     through entire column and compare to critical value; then calculate 
c     fractional cover; convert pstarf from mb to Pa, and winds from km/min
c     to m/s
c
          istagr=.true. 
          sum = 0.
          do 60 k = 1,nz
            dsigma = sigmaf(k) - sigmaf(k+1)
            if (lstagr) then
              uip = uf(i,j,k)
              uim = uf(i-1,j,k)
              vjp = vf(i,j,k)
              vjm = vf(i,j-1,k)
              vmid =(vf(i,j,k) + vf(i,j-1,k))/2.
            else
              uip = (uf(i,j,k) + uf(i+1,j,k))/2.
              uim = (uf(i,j,k) + uf(i-1,j,k))/2.
              vjp = (vf(i,j,k) + vf(i,j+1,k))/2.
              vjm = (vf(i,j,k) + vf(i,j-1,k))/2.
              vmid = vf(i,j,k)
            endif
            pqip = 100.*pstarf(i,j)*ch2of(i,j,k)
            pqim = pqip
            pqjp = pqip
            pqjm = pqip
            if (uip.lt.0.) pqip = 100.*pstarf(i+1,j)*ch2of(i+1,j,k)
            if (uim.ge.0.) pqim = 100.*pstarf(i-1,j)*ch2of(i-1,j,k)
            if (vjp.lt.0.) pqjp = 100.*pstarf(i,j+1)*ch2of(i,j+1,k)
            if (vjm.ge.0.) pqjm = 100.*pstarf(i,j-1)*ch2of(i,j-1,k)
            qdivu = (pqip*uip - pqim*uim)/(1000.*delfx(i,j,k))
            qdivv = (pqjp*vjp - pqjm*vjm)/(1000.*delfy(i,j,k))
            qdivlat = -100.*pstarf(i,j)*ch2of(i,j,k)*vmid*
     &                 tphioaf(i,j)/1000.
            qdiv = qdivu + qdivv + qdivlat
            sum = sum - qdiv*dsigma
c           write(*,*)'qci,pui,qcim1,puim1'
c           write(*,*) qci,pui,qcim1,puim1
c           write(*,*)'qcj,pvj,qcjm1,pvjm1'
c           write(*,*) qcj,pvj,qcjm1,pvjm1
c           write(*,*)'qdivu,qdivv ',qdivu,qdivv,sum
   60     continue
c
c-----No convection; insufficient moisture convergence or cloud clover
c
          if (sum.lt.3.e-7) then
c           write(*,*)'Insufficient moisture convergence: ',sum
            kbot = 0
            ktop = 0
          else
            cd = (1. - bfac)*sum/4.3e-3
            if (cd.gt.1.) cd = 1.
            if (cd.lt.0.05) then
c             write(*,*)'Insufficient cloud cover: ',1.-bfac,cd
              cd = 0.
              kbot = 0
              ktop = 0
            else
c
c-----Convection is diagnosed!
c
c             write(*,*)'Convection is happening!'
c             write(*,*)'Moisture convergence: ',sum
c             write(*,*)'Convective cloud coverage: ',cd
              if (kfreez.gt.1) kfreez = kfreez + 1
              goto 75
            endif
          endif
c
c-----Find vertical distribution of precipitating stratus if convection not
c     occurring.  Ensure the following:
c       1) the lowest cloud deck found is responsible for rain
c       2) average cloud coverage for the deck is at least 10%, otherwise 
c          no rain
c       3) the depth of the deck is greater than 0.1 (in sigma coords, which 
c          translates to about 1000 m in the lowest model layers).
c       4) the bottom of the deck is no higher than 3000 m AGL
c
   61     continue
          kbot = 0
          ktop = 0
          csav = 0.
          if (csmax.lt.0.1) then
             rainf(i,j) = 0.
             goto 75
          endif
          do 70 k = 1,nz
            if (kbot.eq.0 .and. cs(k).ge.0.5*csmax) kbot = k
            if (kbot.ne.0 .and. cs(k).lt.0.5*csmax) then
              ktop = k - 1
              goto 71
            endif
   70     continue
          ktop = nz
   71     continue
          do 73 k = kbot,ktop
            dsigma = sigmaf(k) - sigmaf(k+1)
            csav = csav + cs(k)*dsigma
   73     continue
          sigtot = sigmaf(kbot) - sigmaf(ktop+1)
          csav = csav/sigtot
          if (sigtot.lt.0.1 .or. kbot.eq.ktop) then
             if (kbot.le.2) then
                cs(ktop+1) = cs(ktop)
             else
                cs(kbot-1) = cs(kbot)
             endif
             goto 61
          endif
          zbot = 0.
          if (kbot.gt.1) zbot = zi(kbot-1)
          if (csav.lt.0.1 .or. zbot.gt.3000.) then
            csav = 0.
            kbot = 0
            ktop = 0
            rainf(i,j) = 0.
          endif
c         write(*,*)'Top at   : ',ktop,htop,sigmaf(ktop+1)
c         write(*,*)'Bottom at: ',kbot,hbot,sigmaf(kbot)
c         write(*,*)'Height of freezing: ',kfreez,htfreez,sigmaf(kfreez)
c         write(*,*)'Average stratus coverage: ',csav

   75     continue
c
c-----Load common variables for use in cloud sub-model
c
          kcldbf(i,j) = kbot
          kcldtf(i,j) = ktop
          kfrzf(i,j) = kfreez
          fconvf(i,j) = cd
          fstrtf(i,j) = csav
          do 80 k = 1,nz
c           cldmax = 0.
c           do 85 kk = k,nz
c             cldmax = amax1(cldmax,cs(kk))
c  85       continue
c           cloudf(i,j,k) = cldmax
            if (k.ge.kbot .and. k.le.ktop) then
              if (cd.ne.0.) then
                cloudf(i,j,k) = amax1(cd,cs(k))
                if (cd .ge. cs(k)) ffnonpr(i,j,k) = 0.0
                if (cd .lt. cs(k)) ffnonpr(i,j,k) = 1.0-(cd/cs(k))
              else
                cloudf(i,j,k) = amax1(csav,cs(k))
                if (csav .ge. cs(k)) ffnonpr(i,j,k) = 0.0
                if (csav .lt. cs(k)) ffnonpr(i,j,k) = 1.0-(csav/cs(k))
              endif
            else
              cloudf(i,j,k) = cs(k)
              ffnonpr(i,j,k) = 1.0
            endif
c             write (*,*) 'clddiagf', i,j,k,ffnonpr(i,j,k)
   80     continue 
c
c-----Estimate convecivv results on the concentration (ppb) -----------------
       IF(fconvf(i,j).gt.0) then
!       call convmix(nz,cd,kbot,ktop,kfreez,dt,delfx(i,j,1),delfy(i,j,1),
!     &      pstarf(i,j),xmdf(i,j),rainf(i,j),sigmaf,zi,ti,pi,alfa)
       ENDIF


       IF(fstrtf(i,j).gt.0) then
!        call strtmix(nz,csav,kbot,ktop,kfreez,dt,delfx(i,j,1),
!     &        delfy(i,j,1),
!     &      rainf(i,j),pstarf(i,j),sigmaf,zi,ti,pi,qi,alfa)        
       ENDIF 

      DO 9001 k=1,nz
       conc(i,j,k)=alfa(k)
 9001 CONTINUE

 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  200   continue
  100 continue
      return
      end
