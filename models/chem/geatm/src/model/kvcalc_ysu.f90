      subroutine kv_ysu(nz,pbl,ustar,eli,wstar,zm,zz,thetav,uwind, &
            vwind,temp,qvap,kvmin,risfc,rkv)
!
!-----Determines Kv using the YSU (Yonsei University) methodology. 
!     YSU scheme developed by Song-You Hong, Yign Noh, and Jimy Dudhia
!     Monthly Weather Review Volume 134, Sept 2006.
!
!     PBLs are taken directly from WRF.
!
!     Arguments:
!       nz     number of layers
!       pbl    PBL height (m)
!       ustar  friction velocity (m/s)
!       eli    inverse Monin Obukhov length
!       wstar  convective velocity scale
!       zm     midlayer heights (m)
!       zz     layer interface height (m)
!       thetav virtual potential temperature (K)
!       uwind  U-wind component (m/s)
!       vwind  V-wind component (m/s)
!       temp   temperature (K)
!       qvap   water vapor (ppm)
!       kvmin  Minimum Kv value (m2/s)
!       risfc  bulk richardson number in surface layer
!       rkv    diffusivity (m2/s)

!     Local variables
!       rkloc  diffusivity above pbl (m2/s)
!       rkent  diffusivity in entrainment zone (m2/s)
!       phih,phim stability profile functions
!       prnum0 Prandtl number
!       prfac  
!       wscalek
!       prnumfac
!       zfac
!       wm2,wm3
!       delta
!       bfxpbl  buoyancy flux at pbl height
!       dthetav change in virtual potential temp
!       we      entrainment rate
!       entfac  entrainment factor

      implicit none

      integer nz
      real pbl,ustar,eli,wstar,kvmin,risfc
      real zm(nz),zz(nz),thetav(nz),uwind(nz),vwind(nz),temp(nz), &
           qvap(nz),rkv(nz)

      integer k,kpbl
      real vk,ss,du,dv,dtvdz,tvavg,qavg,tavg,xlv,alph, &
           cpair,chi,zk,rl2,rlam,ric,grav,rdry,rvap
      real rkloc,rkent,phim,phih,prnum,prnum0,prnumfac, &
           prfac,zfac,wscalek,wm2,wm3,xkzo,delta,bfxpbl,dthetav, &
           we,entfac,rig,rlamdz

      real, allocatable, dimension(:) :: qv(:),dzm(:)

      real phifac         !wind profile function evaluated at the
                          !top of the sfc layer
      real sfcfrac        !ratio of surface layer to PBL height
      real prmin,prmax    !min and max Prandtl number
      real zfmin          !minimum ratio of thickness to PBL
 
      data vk /0.4/, grav /9.8/, rdry /287./, &
           rlam /30./, ric /0.25/, rvap /461./,xlv /2.5e6/

      data phifac /8./, sfcfrac /0.1/, prmin /0.25/, prmax /4.0/, &
           zfmin /1.e-8/

      allocate (qv(nz),dzm(nz))

!
!-----Initialize
!
      do k=1,nz
        rkv(k) = 0.0
      enddo
!
!-----Get some PBL parameters
! 
      do k = 1,nz
         qv(k) = (qvap(k)/1.e6)*18./28.8     !mixing ratio(kg/kg)
      enddo
      do k = 2,nz
         dzm(k-1) = zm(k) - zm(k-1)
      enddo
      dzm(nz) = 0.
 
      do k = 2,nz
         if (zz(k).gt.pbl) then
           kpbl = k-1
           goto 10
         endif
      enddo
      kpbl = nz - 1
 10   continue

!
!-----Compute profile functions and the Prandtl number at the top of 
!     the surface layer.  Use the Monin Obukhov length from micromet.

      if (risfc .gt. 0.0) then
!       Stable
        phim = 1.0 + 5.0 *sfcfrac*pbl*eli
        phih = phim
        prnum0 = 1.0
      else
!       Unstable/neutral
        phim = (1.0 - 16.0 *sfcfrac*pbl*eli)**(-1./4.)
        phih = (1.0 - 16.0 *sfcfrac*pbl*eli)**(-1./2.)
        prfac = 6.8*vk*sfcfrac /phim /       &
          ((1. + 4.*vk*wstar**3./ustar**3.)**(1./3.))
        prnum0 = phih/phim + prfac
        prnum0 = min(prnum0,prmax)
        prnum0 = max(prnum0,prmin)
        prnum0 = prnum0/(1.0 + prfac)
      endif

! 
!-----Find kv profile within PBL.  Compute Km.  Then scale to a 
!     vertically-varying function of the Prandtl number to get Kh
!
      do k=1,kpbl
        zfac = min(max((1.0 - zz(k)/pbl),zfmin),1.)   
        if (risfc .gt. 0.0) then
!         Stable
          wscalek = ustar/(1.0 + (phim-1.0)*zz(k)/(sfcfrac*pbl))
          wscalek = max(wscalek, ustar/5.)
          prnum = 1.0
        else
!         Unstable
          wscalek =(ustar**3. +phim*vk*wstar**3.*(1.0-zfac))**(1./3.)   
          prnumfac = -3. * (max(zz(k)-sfcfrac*pbl,0.0))**2.0/   &
                   (pbl**2.0)
          prnum = 1.0 + (prnum0 - 1.0) * exp(prnumfac)
        endif

        rkv(k) = (0.001 * dzm(k)) + wscalek * vk * zz(k)* zfac**2.0 &
                 /prnum
        rkv(k) = max(rkv(k),kvmin)      
      enddo

!
!-----Prepare calculations inside entrainment zone 
!
      wm3   = (wstar**3. + 5.0 * ustar**3.)
      wm2   = wm3**(2./3.)
      delta = 0.02 * pbl + 0.05 * wm2 / (grav/thetav(1)*0.001*pbl)
      delta = min(delta,100.)

!     Assume buoyancy flux at PBL = 15% surface flux.
!     Calculate entrainment rate at PBL(we)
      bfxpbl = -0.15 * thetav(1)/grav * wm3/pbl
      if (kpbl .gt. 1) then
        dthetav = max(thetav(kpbl)-thetav(kpbl-1),0.01)
        we = max(bfxpbl/dthetav, -sqrt(wm2))
      else
        we = -sqrt(wm2)
      endif 
!
!-----Find kv profile above PBL -- 
!     Above entrainment: free atmosphere diffusion. Local K approach (Louis, 1979)
!     Within entrainment zone: geometric average of local Kv and entrainment Kv
! 
      if (kpbl .gt. nz-2) goto 100
      do k=kpbl+1,nz-1

!        Compute local kv 
         xkzo = 0.001 * dzm(k)
         du = uwind(k+1) - uwind(k)
         dv = vwind(k+1) - vwind(k)
         ss = (du**2. + dv**2.)/(dzm(k)**2.) + 1.e-9

         dtvdz = (thetav(k+1) - thetav(k))/dzm(k)
         tvavg = (thetav(k+1) + thetav(k))/2.
         rig = grav*dtvdz/(ss*tvavg)

         zk = vk*zz(k)
         rlamdz = min(max(0.1*dzm(k),rlam),300.)
         rl2 = (zk*rlamdz/(rlamdz + zk))**2
!
!------  Compute Kh
         if (rig .lt. 0.) then
!          Unstable 
           rkloc = xkzo + (rl2*sqrt(ss))*(1.-8.0*rig/(1.+1.286*  &
                    sqrt(-1.0*rig)))
         else
!          Stable
           rkloc = xkzo + (rl2*sqrt(ss))/(1.0+5.*rig)**2.
         endif
         rkloc = max(rkloc,kvmin)

!
!-----   Find kv profile above PBL, but within entrainment layer
!        Entrainment zone ends when the flux is reduced to 1%.  This
!        occurs when (z-pbl)**2/delta**2 >= 4.6  (e(-4.6) = 0.01)
     
         entfac = (zz(k) - pbl)**2./delta**2.
         if (entfac .lt. 4.6) then
           rkent = -we * dzm(kpbl) * exp(-entfac)
           rkv(k) = sqrt(rkent * rkloc)
           rkv(k) = max(rkv(k),kvmin)
         else
           rkv(k) = rkloc
           rkv(k) = max(rkv(k),kvmin)
         endif
       enddo
        
100    rkv(nz) = 0.

       deallocate (qv,dzm)

       return
       end
