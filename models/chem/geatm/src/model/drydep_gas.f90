

       subroutine drydep_gas(myid,cflag,fcloud, &
                       month,tsurf,cellat,cellon,pwc,cwc,height,&
                       press,windu,windv, &
                       solflux,landuse,water,&
                       lrdlai,wrflai, &
                       tempk,pbl_hgt,lrddrt,icddrt,lrdsno,icdsno,vdep,&
                       itzon,iyear,imonth,iday,ihour,iminute, &
                       sx,ex,sy,ey,nzz,ig )
!c
!c-----CAMx v4.42 070603
!c
!c
!c     DRYDEP is the driver for the calculation of gridded dry deposition
!c     velocities for a given grid. Deposition velocities are calculated for
!c     each gas species and for each aerosol size bin, weighted by the
!c     fractional land use specified for each cell.
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        4/4/00    Added aerosol deposition as f(size)
!c        4/4/01    Fixed a few bugs in the call for VD_AER
!c        1/9/02    Aerosol size and density now species-dependent
!c        3/26/03   Added scaling factor to surface resistance (provided
!c                  on chemparam file), and zero dry deposition for
!c                  insoluble gases
!c        4/9/03    Removed rain, added precip and cloud water contents:
!c                  surfaces can now be rain-wetted, fog-wetted, or dew-wetted
!c        6/6/03    Protect against divide by zero with totland
!c        6/11/03   Use optional surface roughness length, if available
!c        7/21/03   Use optional drought stress and snow cover, if available;
!c                  Introduced latitude-dependent specification of season;
!c                  Revised solar flux calculation to use RADM cloud adjustment
!c        8/18/03   Relate drought stress to Palmer Drought Index
!c        4/21/04   Incorporated sectional PM
!c        11/19/04  Incorporated season-dependent roughness length
!c
!c     Input arguments:
!c        ncol                number of columns
!c        nrow                number of rows
!c        nlay                number of layers
!c        tsurf               surface temperature field (K)
!c        cellat              cell centroid latitude (deg)
!c        cellon              cell centroid longitude (deg)
!c        pwc                 precipitation water content (g/m3)
!c        cwc                 cloud water content (g/m3)
!c        height              layer interface height field (m)
!c        press               layer pressure field (mb)
!c        windu               layer U-component wind field (m/s)
!c        windv               layer V-component wind field (m/s)
!c        solflux             solar flux(w/m2)
!c        landuse             landuse
!c        water               layer water vapor field (ppm)
!c        tempk               layer temperature field (K)
!c        lrddrt              flag that gridded drought stress is available
!c        icddrt              optional drought index
!c        lrdsno              flag that gridded snow cover is available
!c        icdsno              optional snow index
!c
!c     Output arguments:
!c        vdep                species-dependent deposition velocity field (m/s)
!c
!c     Routines called:
!c        CALDATE
!c        GETZNTH
!c        MICROMET
!c        VD_GAS
!c        VD_AER
!c        HENRYFNC
!c
!c     Called by:
!c        CAMx
!c
!c
      ! add 'z03' dry deposition scheme chenxsh@mail.iap.ac.cn,20151029

      use dep_parm

      implicit none


!      character(len=*),parameter :: cflag='W89'
      character(len=*) :: cflag ! W89

      integer :: sx,ex,sy,ey,nzz 
      real tsurf(sx-1:ex+1,sy-1:ey+1),cellat(sx-1:ex+1,sy-1:ey+1),&
         cellon(sx-1:ex+1,sy-1:ey+1)
      integer icddrt(sx-1:ex+1,sy-1:ey+1),&
         icdsno(sx-1:ex+1,sy-1:ey+1)
      real height(sx-1:ex+1,sy-1:ey+1,nzz),&
          press(sx-1:ex+1,sy-1:ey+1,nzz),&
          windu(sx-1:ex+1,sy-1:ey+1,nzz),&
          windv(sx-1:ex+1,sy-1:ey+1,nzz),&
          water(sx-1:ex+1,sy-1:ey+1,nzz),&
          tempk(sx-1:ex+1,sy-1:ey+1,nzz),pwc(sx-1:ex+1,sy-1:ey+1,nzz),&
          cwc(sx-1:ex+1,sy-1:ey+1,nzz)

      integer :: iyear,imonth,iday,ihour,iminute
      integer :: jdate
!      real zmol(sx-1:ex+1,sy-1:ey+1)
      real fcloud(sx-1:ex+1,sy-1:ey+1,nzz)


      real vdep(sx-1:ex+1,sy-1:ey+1),pbl_hgt(sx-1:ex+1,sy-1:ey+1)
      logical lstable
      logical ldark,lrdruf,lrddrt,lrdsnoi
      REAL :: solflux(sx-1:ex+1,sy-1:ey+1)
      INTEGER ::  landuse(sx-1:ex+1,sy-1:ey+1)
!      REAL :: henry0(74),tfact(74),f0(74),rscale(74)
!      REAL :: diffrat(74)
      real :: henry !,dstress(0:5)
!      REAL,DIMENSION(11,5) :: rj,rlu,rac,rgss,rgso,rlcs,rlco,z0lu
!c

      logical :: lrdlai
      real    :: wrflai(sx-1:ex+1,sy-1:ey+1)

      real :: el

      real :: zenith
      real :: lai_f
      real*8 :: srad_dbl,coszen_dbl


      integer :: myid
      integer :: month
      integer :: knh3,khno3,kso2
      integer :: i,j,ig,m,mbin,latbin,isesn
      integer :: itzon,itime,iwet,istress
      logical :: lrdsno
      real    :: deltaz,temp0,prss0,ucomp,vcomp,wind,qwatr,ev,es,rh,dew
      real    :: psih
      real    :: pi,cwmin
      real    :: totland,z0,coszen,pbl,ustar,wstar,zl
      real    :: vd


      real    :: henso2,henso20,tfactso2
      integer :: ko3,iflgso2,iflgo3
      



      integer, dimension(12) :: nday

      data nday/31,28,31,30,31,30,31,31,30,31,30,31/


      data pi/3.1415927/
      data cwmin /3.8E-5/



!c
!c-----Entry point
!c

      jdate=iyear*10000+imonth*100+iday
      call juldate(jdate) 

      itime=ihour*100+iminute

      month=imonth

!c
!c-----Loop over rows and columns
!c
      do 30 j = sy,ey
      do 20 i = sx,ex
!c
!c-----Determine season
!c
          
          if(diffrat(ig).eq.0) then
            vd=0.0
            vdep(i,j) = vd
            cycle
          endif


          mbin = month
          if (cellat(i,j).lt.0.) then
            mbin = mod(month+6,12)
            if (mbin.eq.0) mbin = 12
          endif
          latbin = 1
          if (abs(cellat(i,j)).gt.20.) then
            latbin = 2
          elseif (abs(cellat(i,j)).gt.35.) then
            latbin = 3
          elseif (abs(cellat(i,j)).gt.50.) then
            latbin = 4
          elseif (abs(cellat(i,j)).gt.75.) then
            latbin = 5
          endif
          if ((cellat(i,j).gt.50. .and. cellat(i,j).lt.75.) .and.&
             (cellon(i,j).gt.-15. .and. cellon(i,j).lt.15.)) latbin = 3
          isesn = iseason(latbin,mbin)
!c
!c-----Use input snow cover to set season, if specified
!c
          if (lrdsno .and. icdsno(i,j).eq.1) isesn = 4
!c
!c-----Calculate solar flux
!c
          solflux(i,j) = amax1(0.,solflux(i,j))
!c
!c-----Load local met variables
!c
          deltaz = height(i,j,1)/2.
          temp0 = tsurf(i,j) - 273.15
          prss0 = press(i,j,1) -&
                 2.*deltaz*(press(i,j,2) - press(i,j,1))/height(i,j,2)
          ucomp = (windu(i,j,1) + windu(i-1,j,1))/2.
          vcomp = (windv(i,j,1) + windv(i,j-1,1))/2.
          wind = sqrt(ucomp**2 + vcomp**2)
          wind = amax1(0.1,wind)

          if(deltaz.lt.0) then
            print*,'deltaz=',deltaz,i,j
            stop
          endif

!c
!c-----Determine surface wetness
!c
          qwatr = 1.e-6*water(i,j,1)*18./28.8
          ev = qwatr*prss0/(qwatr + eps)
          es = e0*exp((lv/rv)*(1./273. - 1./tsurf(i,j)))
          rh = amin1(1.,ev/es)

          iwet = 0

          if (pwc(i,j,1).gt.cwmin) then
            iwet = 2
          elseif (cwc(i,j,1).gt.cwmin) then
            iwet = 1
          else
            dew = (1. - rh)*(wind + 0.6)
            if (dew.lt.0.19) iwet = 1
          endif
!c
!c-----Use input drought stress, if specified
!c
          istress = 0
          if (lrddrt .and. icddrt(i,j).gt.0) istress = icddrt(i,j)
!c
!c-----Loop over land use; surface roughness for water is dependent on
!c     wind speed
         vdep(i,j) = 0.

          totland = 1.
!c
!            m=nint(landuse)
 
         if(trim(cflag).eq.'W89') then 

            m=landuse(i,j)
            z0 = z0lu(m,isesn)
            if (m.eq.7) z0 = amax1(z0,2.0e-6*wind**2.5)

          elseif(trim(cflag).eq.'Z03') then

            m=landuse(i,j)

!            print*,'m mbin=',m,mbin,i,j

!            print*,lai_ref(m,mbin)

            if(lrdlai) then   
             lai_f = wrflai(i,j)
             lai_f = amin1(lai_ref(m,15),lai_f)
             lai_f = amax1(lai_ref(m,14),lai_f)
            else
            lai_f = lai_ref(m,mbin) + &
     &                float(min(nday(mbin),iday))/float(nday(mbin))* &
     &                (lai_ref(m,mbin+1) - lai_ref(m,mbin))

            endif
     
!            print*,'lai_f=',lai_f            

            !if (lrdlai) then
            !    lai_f = lai_f*rlai
            !    lai_f = amin1(lai_ref(m,15),lai_f)
            !    lai_f = amax1(lai_ref(m,14),lai_f)
            !endif
            if (m.eq.1 .or. m.eq.3) then
               z0 = 2.0e-6*wind**2.5
            else
               if (z02(m).gt.z01(m)) then
                  z0 = z01(m) + (lai_f - lai_ref(m,14))/ &
     &                          (lai_ref(m,15) - lai_ref(m,14))* &
     &                          (z02(m) - z01(m))
               else
                  z0 = z01(m)
               endif
            endif

            call getznth(cellat(i,j),cellon(i,j),itime,jdate,itzon,zenith,ldark)
            coszen = cos(DBLE(zenith)*DBLE(pi)/180.)

          else
            stop 'naqpms-ddep-flag erro in drydep_gas.f90'
         
          endif
!c
!c-----Use input surface roughness, if specified
!c
!c
!c-----Get surface layer micrometeorological parameters for this cell and
!c     landuse type
!c
            if (prss0.lt.0) then
              write(*,'(//,a)') 'ERROR in DRYDEP:'
              write(*,*) 'Invalid pressure value'
              write(*,*) 'Cell   Height  Deltaz'
              write(*,*) i,j,height(i,j,1),deltaz,press(i,j,1),press(i,j,2),prss0
            endif
            
          pbl = pbl_hgt(i,j)
            
!            call micromet(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,&
!                         deltaz,wind,z0,pbl,ustar,psih,wstar,lstable)
             

            call micromet_el(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,&
                          deltaz,wind,z0,pbl,ustar,el,psih,wstar,lstable)

            zl = deltaz/el
            if (zl.lt.0.) then
              zl = min(-5.,zl)
            else
              zl = max(5.,zl)
            endif

!c
!c-----Loop over GAS species, and calculate deposition velocity for this cell,
!c     landuse, and current species.
!c     Use input drought stress code, if specified
!c
!           if(i==53.and.j==65.and.ig==11)print*,tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0
               knh3  =  4
               khno3 =  2
               kso2  = 18
               ko3   = 11
!               IF(ig==18) THEN   !SO2
                henso20  = 1.0e+5
                tfactso2 = -3156
              call henryfnc(0,henso20,tfactso2,tsurf(i,j),7.,knh3,khno3,&
                         kso2,henso2)
!               ENDIF
              if (henry0(ig).gt.1.e-6) then
                call henryfnc(ig,henry0(ig),tfact(ig),tsurf(i,j),7.,knh3,&
                             khno3,kso2,henry)
                iflgso2 = 0
                iflgo3 = 0
                if (ig.eq.kso2) then
                  iflgso2 = 1
                  henry = henso2
                endif
                if (ig.eq.ko3) iflgo3 = 1

               if(trim(cflag).eq.'W89') then
                call vd_gas(i,j,ig,m,istress,iwet,iflgso2,iflgo3,z0,&
                          deltaz,psih,ustar,diffrat(ig),henry,henso2,&
                          f0(ig),rscale(ig),temp0,dstress(0),solflux(i,j),&
                          rj(m,isesn),rlu(m,isesn),rac(m,isesn),&
                          rlcs(m,isesn),rlco(m,isesn),rgss(m,isesn),&
                          rgso(m,isesn),vd)
                elseif(trim(cflag).eq.'Z03') then
                  ! zl,pwc
                  !fcloud(i,j,1)=0.0
                  !print*,'henso2=',henso2

                  srad_dbl=solflux(i,j)
                  coszen_dbl=coszen

!                  print*,'lai_f to vd_gas_zhang=',lai_f               
   
                  call vd_gas_zhang(deltaz,zl,z0,ustar,tempk(i,j,1), &
     &                      tsurf(i,j),srad_dbl,rh,fcloud(i,j,1), &
     &                      pwc(i,j,1),coszen_dbl,m,isesn,henry, &
     &                      henso2,ig,lai_f,vd)

                  if(.not.(vd.ge.0.and.vd.le.1.0e3)) then
                    print*,'i&j&ig=',ig,i,j,vd
                    print*,deltaz,zl,z0,ustar
                    print*,tempk(i,j,1),tsurf(i,j)
                    print*,srad_dbl,rh,fcloud(i,j,1),pwc(i,j,1)
                    print*,coszen_dbl
                    print*,m,isesn
                    print*,henry,henso2,lai_f
                    stop 'vd-erro'
                  endif

                endif              
              else
                vd = 0.
              endif
              vdep(i,j) =  vd
!c
 20     continue
 30   continue
!c
      return
      end




     subroutine vd_gas(i,j,ig,ilu,istress,iwet,iso2,io3,z0,deltaz,psih,ustar,&
                       diffrat,henry,henso2,f0,rscale,ts,dstress,&
                       solflux,rj,rlu,rac,rlcs,rlco,rgss,rgso,vd)
!c
!c----CAMx v4.42 070603
!c
!c     VD_GAS calculates a deposition velocity for a specific gas species,
!c     grid cell, and land use category.  The parallel resistance approach of
!c     Wesely and Hicks (1977) is used with the improvements of Wesely (1989).
!c     Surface resistance (rs) to water (landuse 7) is determined following
!c     Sehmel (1980), as implemented in UAM-AERO (STI, 1996).
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        8/18/03       Relate drought stress to PDI index
!c        3/21/03       Modified Ra calculation to use layer 1 midpoint height
!c                      rather than default 10 m; changed drought stress effects
!c                      to stomatal resistance to be consistent with effects in
!c                      GLOBEIS
!c        3/26/03       Added scaling factor to surface resistance (provided
!c                      on chemparam file)
!c
!c     Input arguments:
!c        ilu                 land use index
!c        istress             vegetation drought stress index
!c        iwet                surface wetness index
!c                            0 = dry
!c                            1 = dew wetted
!c                            2 = rain wetted
!c        iso2                SO2 species flag (1=SO2,0=other)
!c        io3                 O3 species flag (1=O3,0=other)
!c        z0                  surface roughness length (m)
!c        deltaz              Layer 1 midpoint height (m)
!c        psih                similarity stability correction term
!c        ustar               friction velocity (m/s)
!c        diffrat             ratio of molecular diffusivity of water to species
!c        henry               Henry's Law constant (M/atm)
!c        henso2              Henry's Law constant for SO2 (M/atm)
!c        f0                  normalized reactivity parameter
!c        rscale              user-defined surface resistance scaling factor
!c        ts                  surface temperature (C)
!c        dstress             adjustment factors for drought stress
!c        solflux             Solar radiation flux (W/m2)
!c        rj                  baseline minimum stomatal resistance (s/m)
!c        rlu                 baseline upper canopy (cuticle) resistance (s/m)
!c        rac                 baseline canopy height/density resistance (s/m)
!c        rlcs                baseline SO2 lower canopy resistance (s/m)
!c        rlco                baseline O3 lower canopy resistance (s/m)
!c        rgss                baseline SO2 ground surface resistance (s/m)
!c        rgso                baseline O3 ground surface resistance (s/m)
!c
!c     Output arguments:
!c        vd                  deposition velocity (m/s)
!c
!c     Routines called:
!c        none
!c
!c     Called by:
!c        DRYDEP
     real    dstress(0:5)
!c
      data vk/0.4/, rmin/1.0/, rmax/1.e5/, d1/2./, d2/0.667/
      data vair/1.5e-5/, diffh2o/2.3e-5/
!c
!c-----Entry point
!c
!c-----Compute atmospheric resistance, RA
!c
      ra = (alog(deltaz/z0) - psih)/(vk*ustar)
      ra = amax1(ra,rmin)
!c
!c-----Compute the deposition layer resistance, RD
!c
      schmidt = vair*diffrat/diffh2o
      rd = d1*schmidt**d2/(vk*ustar)
      rd = amax1(rd,rmin)
!c
!c-----Compute the surface layer resistance over water, RS
!c
      if (ilu.eq.7) then
        rs = 1./(3.9e-5*henry*ustar*(ts + 273.15))
        rs = amax1(rs,rmin)
        goto 100
      endif
!c
!c
!c-----Compute stomatal resistance, RST
!c     Adjust for vegetation drought stress
!c
      rst = rmax
      if (ts.gt.0. .and. ts.lt.40.) then
        rst = diffrat*rj*(1. + (200./(solflux + 0.1))**2)*&
             (400./(ts*(40.-ts)))
        rst = rst * dstress(istress)
      endif
      if (iwet.gt.0) rst = 3.*rst
      rst = amin1(rmax,rst)
!c
!c-----Compute mesophyll resistance, RM
!c
      rm = 1./(henry/3000. + 100.*f0)
      rm = amax1(rmin,rm)
      rm = amin1(rmax,rm)
!c
!c-----Compute upper canopy resistance, RUC
!c     Adjust for surface wetness
!c
      if (iwet.eq.0) then
        ruc = rlu/(henry/henso2 + f0)
        ruc = amin1(rmax,ruc)
      else
        if (iwet.eq.1) then
          rlus = 100.
          if (ilu.eq.1) rlus = 50.
!c
!c  --- original equations from Wesely 89 ---
!c         rluo = 1./(1./3000. + 1./(3.*rlu))
!c
          rluo = 1000. + rlu
        else
!c
!c  --- original equations from Wesely 89 ---
!c         rlus = 1./(1./5000. + 1./(3.*rlu))
!c
          rlus = 2000. + rlu
          if (ilu.eq.1) rlus = 50.
          rluo = 1./(1./1000. + 1./(3.*rlu))
        endif
        if (iso2.eq.1) then
          ruc = rlus
        elseif (io3.eq.1) then
          ruc = rluo
        else
!c
!c  --- original equations from Wesely 89 ---
!         ruc = 1./(1./(3.*rlu) + 1.e-7*henry + f0/rluo)
!c
       ruc = 1./(henry/(henso2*rlus) + f0/rluo)

          ruc = amin1(rmax,ruc)
        endif
      endif
!         if(i==53.and.j==65.and.ig==11) print*,ruc,rmax,henso2
!c
!c-----Compute resistance to buoyant convection, RDC
!c     (assume effect of local terrain slope is non-resolvable; this
!c     factor is set to 1)
!c
      rdc = 100.*(1. + 1000./(solflux + 10.))
      rdc = amin1(rmax,rdc)
!c
!c-----Compute resistance of exposed surfaces in lower canopy, RCL
!c
      rlc = 1./(henry/(henso2*rlcs) + f0/rlco)
      rlc = amin1(rmax,rlc)
!c
!c-----Compute resistance of ground surface, RGS
!c
      rgs = 1./(henry/(henso2*rgss) + f0/rgso)
      rgs = amin1(rmax,rgs)
!c
!c-----Determine total surface resistance over land, RS
!c
      rs = 1./(rst + rm) + 1./ruc + 1./(rdc + rlc) + 1./(rac + rgs)
      rs = amax1(rmin,1./rs)
!c
!c-----Scale surface resistance for acids according to user-definition
!c
 100  continue
      rs = rs*rscale
!c
!c-----Final deposition velocity for this cell, land use, and species
      vd = 1./(ra + rd + rs)
!c
!      if(i==53.and.j==65.and.ig==11)print*,vd,ra,rd,rs,deltaz,psih,ustar
      return
      end
      subroutine micromet(temp,temp0,press,press0,deltaz,wind,z0,pbl,&
                         ustar,psih,wstar,lstable)
!c
!c----CAMx v4.42 070603
!c
!c     MICROMET calculates surface layer micro-meteorological flux-gradient
!c     relationships and variables based on Louis (1979)
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        8/15/02    Specified limits for denominators from 1e-10 to 1e-20
!c        3/21/03    Modified Z/L calculation to use input midpoint height
!c                   (instead of default 10 m) and to increase range from
!c                   (-1 to +1) to (-2.5 to +1.5)
!c        8/20/03    added calculation of w* and logical for stability
!c
!c     Input arguments:
!c        temp                Layer 1 temperature (K)
!c        temp0               Surface temperature (K)
!c        press               Layer 1 pressure (mb)
!c        press0              Surface pressure (mb)
!c        deltaz              Layer 1 midpoint height (m)
!c        wind                Layer 1 total wind speed (m/s)
!c        z0                  Surface roughness length (m)
!c        pbl                 PBL depth (m)
!c
!c     Output arguments:
!c        ustar               friction velocity (m/s)
!c        psih                Similarity stability correction term
!c        wstar               convective velocity scale (m/s)
!c        lstable             stability logical (T=stable)
!c
!c     Routines called:
!c        none
!c
!c     Called by:
!c        DRYDEP
!c        PIGGROW
!c
     logical lstable
      data vk/0.4/, g/9.8/, gamma/0.286/
!c
!c-----Entry point
!c
!c-----Calculate potential temperature and richardson number
!c
      theta = temp*(1000./press)**gamma
      theta0 = temp0*(1000./press0)**gamma
      dtheta = theta - theta0
      thetabar = (theta + theta0)/2.
      ri = (g/thetabar)*deltaz*dtheta/(wind**2 + 1.e-20)
!c
!c-----Determine stability functions
!c
      zscale = vk/alog(deltaz/z0)
      if (ri.le.0.) then
        lstable = .false.
        cm    = 69.56*sqrt(deltaz/z0)*zscale**2
        ch    = 49.82*sqrt(deltaz/z0)*zscale**2
        fm    = 1. - 9.4*ri/(1. + cm*sqrt(abs(ri)))
        fh    = 1. - 9.4*ri/(1. + ch*sqrt(abs(ri)))
      else
        lstable = .true.
        fm = 1./((1. + 4.7*ri)**2)
        fh = fm
      endif
!c
!c-----Calculate micromet variables
!c
      ustar2 = fm*(wind*zscale)**2
      ustar2 = amax1(1.e-20,ustar2)
      ustar = sqrt(ustar2)
      thstar = 1.35*zscale**2*wind*dtheta*fh/ustar
      el = ustar2*temp/(vk*g*thstar + 1.e-20)
      elabs = abs(el)
      wstar = 0.
      if (el.lt.0.) wstar = (pbl*ustar**3./(vk*elabs))**(1./3.)

      if (elabs.ge.1.e4) then
        psih = 0.0
      elseif (el.lt.0.) then
        zoverl = amin1(2.5,deltaz/elabs)
        zmel = alog(zoverl)
        psih = exp(0.598 + 0.39*zmel - 0.090*zmel*zmel)
      else
        zoverl = amin1(1.5,deltaz/elabs)
        psih = -5.*zoverl
      endif
!c
      return
      end
     subroutine henryfnc(ispc,hlaw0,tfact,temp,ph,knh3,khno3,kso2,hlaw)
!c
!c----CAMx v4.42 070603
!c
!c     HENRYFNC calculates temperature and dissociation adjustments to
!c     baseline Henry's Law constants.
!c
!c     Copyright 2006-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        None
!c
!c     Input arguments:
!c        ispc                Species index
!c        hlaw0               Baseline Henry's Law constant @298K (M/atm)
!c        tfact               temperature factor
!c        temp                ambient temperature (K)
!c        ph                  pH of liquid
!c        knh3                pointer to NH3
!c        khno3               pointer to HNO3
!c        kso2                pointer to SO2
!c
!c     Output arguments:
!c        hlaw                Adjusted Henry's Law constant (M/atm)
!c
!c     Routines called:
!c        None
!c
!c     Called by:
!c        WETDEP
!c        DRYDEP
!c
      implicit none
      integer ispc,knh3,khno3,kso2
      real hlaw0,tfact,temp,ph,hlaw

      real diss1,diss2
!c
      hlaw = hlaw0*exp(tfact*(1./298. - 1./temp))
      if (ispc.eq.knh3) then
        diss1 = 10.**(-189.1/temp - 4.117)
        diss2 = 10.**(-5839.5/temp - 9.7618*alog(temp) + 61.206)
        hlaw = hlaw*(1. + (diss1/diss2)*10.**(-ph))
      elseif (ispc.eq.khno3) then
       diss1 = 15.4
        hlaw = hlaw*(1. + diss1/(10.**(-ph)))
      elseif (ispc.eq.kso2) then
        diss1 = 10.**(853./temp)/54950.
        diss2 = 10.**(621.9/temp)/1.897e+9
        hlaw = hlaw*(1. + diss1/(10.**(-ph)) +&
                         diss1*diss2/(10.**(-2.*ph)))
      endif

      return
      end

