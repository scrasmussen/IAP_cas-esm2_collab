      subroutine drydep_aer_old(myid,cflag,month,tsurf,cellat,cellon,height,&
                            press,windu,windv,landuse,tempk,&
                            lrdlai,wrflai, &
                    pbl_hgt,diam2d,vdep_aer,sx,ex,sy,ey,nzz,ig)

! add 'z03' dry deposition scheme chenxsh@mail.iap.ac.cn,20151029

      use dep_parm

      character(len=*) :: cflag
!      character(len=*),parameter :: cflag='Z03' ! W89

      integer :: sx,ex,sy,ey
  
      real diam2d(sx-1:ex+1,sy-1:ey+1) 

      real tsurf(sx-1:ex+1,sy-1:ey+1),cellat(sx-1:ex+1,sy-1:ey+1),&
          cellon(sx-1:ex+1,sy-1:ey+1)
      real height(sx-1:ex+1,sy-1:ey+1,nzz),&     
          press(sx-1:ex+1,sy-1:ey+1,nzz),&
          windu(sx-1:ex+1,sy-1:ey+1,nzz),&
          windv(sx-1:ex+1,sy-1:ey+1,nzz),&
          tempk(sx-1:ex+1,sy-1:ey+1,nzz)
      real vdep_aer(sx-1:ex+1,sy-1:ey+1),pbl_hgt(sx-1:ex+1,sy-1:ey+1)

!       real zmol(sx-1:ex+1,sy-1:ey+1)

      logical lstable
!      REAL :: iseason(5,12)
      REAL :: diam  
! diam : log-mean sectional aerosol diameter (m)
      REAL :: rhop ! AEROSOL DENISITY in g/m3
      INTEGER ::  landuse(sx-1:ex+1,sy-1:ey+1)
!      REAL,DIMENSION(11,5) :: z0lu
!
      integer, dimension(12) :: nday

      real :: zl,el,lai_f

      logical :: lrdlai
      real    :: wrflai(sx-1:ex+1,sy-1:ey+1)


      data nday/31,28,31,30,31,30,31,31,30,31,30,31/


      data pi/3.1415927/
      data cwmin /3.8E-5/
!
!c-----Season indices by month and latitude band
!c     Season Indices            Latitude Bands
!c     1 = summer                1 = <20    Tropical
!c     2 = autumn                2 = 20-35  Sub-tropical
!c     3 = winter w/o snow       3 = 35-50  Temperate
!c     4 = winter w/ snow        4 = 50-75  Cool
!c     5 = spring                5 = >75    Polar
!c                    Latitude Band

      IF(ig == 75 ) rhop = 1.0E06 ! Primary PM25
      IF(ig == 76 ) rhop = 1.0E06 ! Primary PM10
      IF(ig == 77 ) rhop = 2.0E06 ! BC
      IF(ig == 78 ) rhop = 1.0E06 ! Primary OC
      IF(ig == 79 ) rhop = 0.0    ! H+
      IF(ig == 80 ) rhop = 2.0E06 ! Na+
      IF(ig == 81 ) rhop = 1.5E06 ! NH4+
      IF(ig == 82 ) rhop = 2.0E06 ! CL-
      IF(ig == 83 ) rhop = 1.5E06 ! SO42-
      IF(ig == 84 ) rhop = 1.5E06 ! HSO4-
      IF(ig == 85 ) rhop = 1.5E06 ! NO3-
      IF(ig == 86 ) rhop = 2.0E06 ! NACL
      IF(ig == 87 ) rhop = 1.5E06 ! NA2SO4
      IF(ig == 88 ) rhop = 1.5E06 ! NANO3
      IF(ig == 89 ) rhop = 1.5E06 ! NH42SO4
      IF(ig == 90 ) rhop = 1.5E06 ! NH4NO3
      IF(ig == 91 ) rhop = 1.5E06 ! NH4CL
      IF(ig == 92 ) rhop = 1.5E06 ! H2SO4
      IF(ig == 93 ) rhop = 1.5E06 ! NH4HSO4
      IF(ig == 94)  rhop = 1.5E06 ! NAHSO4
      IF(ig == 95)  rhop = 1.5E06 ! (NH4)4H(SO4)2(S)
      IF(ig == 96)  rhop = 1.0E06 ! SOA1
      IF(ig == 97)  rhop = 1.0E06 ! SOA2
      IF(ig == 98)  rhop = 1.0E06 ! SOA3
      IF(ig == 99)  rhop = 1.0E06 ! SOA4
      IF(ig == 100) rhop = 1.0E06 ! SOA5
      IF(ig == 101) rhop = 1.0E06 ! SOA6
      IF(ig == 102) rhop = 1.0E06 ! AH2O       

      do 30 j = sy,ey
        do 20 i = sx,ex

         diam=diam2d(i,j)
!c
!c-----Determine season
!c
!          zl=zmol(i,j)


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
!c
!c-----Load local met variables
          deltaz = height(i,j,1)/2.
          temp0 = tsurf(i,j) - 273.15
          prss0 = press(i,j,1) -&
                 2.*deltaz*(press(i,j,2) - press(i,j,1))/height(i,j,2)
          ucomp = (windu(i,j,1) + windu(i-1,j,1))/2.
          vcomp = (windv(i,j,1) + windv(i,j-1,1))/2.
          wind = sqrt(ucomp**2 + vcomp**2)
          wind = amax1(0.1,wind)
!c
!c
!c-----Loop over land use; surface roughness for water is dependent on        
!c     wind speed
         vdep_aer(i,j) = 0.
         totland = 1.

         if(trim(cflag).eq.'W89') then
                   
            m=landuse(i,j)
            z0 = z0lu(m,isesn)
            if (m.eq.7) z0 = amax1(z0,2.0e-6*wind**2.5)


          elseif(trim(cflag).eq.'Z03') then

            m=landuse(i,j)


            if(lrdlai) then
             lai_f = wrflai(i,j)
             lai_f = amin1(lai_ref(m,15),lai_f)
             lai_f = amax1(lai_ref(m,14),lai_f)
            else
            lai_f = lai_ref(m,mbin) + &
     &                float(min(nday(mbin),iday))/float(nday(mbin))* &
     &                (lai_ref(m,mbin+1) - lai_ref(m,mbin))
            endif
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

          else
            stop 'naqpms-ddep-flag erro in drydep_aer.f90'

          endif




!c

!c-----Use input surface roughness, if specified
!c
!c
!c-----Get surface layer micrometeorological parameters for this cell
!and
!c     landuse type
!c
         if (prss0.lt.0) then
           write(iout,'(//,a)') 'ERROR in DRYDEP:'
           write(iout,*) 'Invalid pressure value'
           write(iout,*) 'Cell   Height  Deltaz'
           write(iout,*) i,j,height(i,j,1),deltaz,press(i,j,1),press(i,j,2),prss0
         endif

         pbl = pbl_hgt(i,j)
!         call micromet(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,&
!                    deltaz,wind,z0,pbl,ustar,psih,wstar,lstable)


         call micromet_el(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,&
                    deltaz,wind,z0,pbl,ustar,el,psih,wstar,lstable)

         zl = deltaz/el
         if (zl.lt.0.) then
           zl = min(-5.,zl)
         else
           zl = max(5.,zl)
         endif
!c

         if(trim(cflag).eq.'W89') then
            call vd_aer(z0,deltaz,psih,ustar,diam,rhop,tsurf(i,j),vd)      

         elseif(trim(cflag).eq.'Z03') then
            call vd_aer_zhang(deltaz,zl,z0,ustar,diam,rhop &
                         ,tsurf(i,j),tempk(i,j,1),m,lai_f,vd)

            if(.not.(vd.ge.0.and.vd.le.1.0e3)) then
               print*,'i&j&ig=',ig,i,j,vd
               print*,deltaz,zl,z0,ustar
               print*,diam,rhop
               print*,tsurf(i,j),tempk(i,j,1)
               print*,m,lai_f
               stop
            endif

         endif

         vdep_aer(i,j) =  vd
        
   20     continue
   30   continue
    
                            
      return
      end      

      subroutine vd_aer_to_be_deleted(z0,deltaz,psih,ustar,diam,rhop,ts,vd)
!c
!c----CAMx v4.42 070603
!c
!c     VD_AER calculates a deposition velocity for a specific aerosol size
!c     bin, grid cell, and land use category.  The parallel resistance approach
!c     of Slinn and Slinn (1980) is used, as implemented in UAM-AERO
!c     (STI, 1996).
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        3/21/03       Modified Ra calculation to use layer 1 midpoint height
!c                      rather than default 10 m.
!c        9/18/03       Removed small bug in final Vd equation
!c
!c     Input arguments:
!c        z0                  surface roughness length (m)
!c        deltaz              Layer 1 midpoint height (m)
!c        psih                similarity stability correction term
!c        ustar               friction velocity (m/s)
!c        diam                log-mean sectional aerosol diameter (m)
!c        rhop                aerosol density (g/m3)
!c        ts                  surface temperature (K)
!c
!c     Output arguments:
!c        vd                  deposition velocity (m/s)
!c
!c     Routines called:
!c        none
!c
!c     Called by:
!c        DRYDEP

      data vk/0.4/, rmin/1.0/, xmfp/6.5e-8/, g/9.8/, vabs/1.81e-2/
      data boltz/1.38e-20/, pi/3.1415927/, vair/1.5e-5/
!c
!c-----Entry point
!c
!c-----Speed correction factor and sedimendation velocity
!c
      power = amin1(7.6,0.55*diam/xmfp)
      scf = 1. + (2.514 + 0.8*exp(-power))*xmfp/diam
      vsed = rhop*g*(diam*diam)*scf/(18.*vabs)
!c
!c-----Brownian diffusivity and Schmidt number
!c
      difbrwn = boltz*ts*scf/(3.*pi*vabs*diam)
      schmidt = vair/difbrwn
!c
!c-----Stokes number
!c
      stokes = vsed*(ustar*ustar)/(vair*g)
!c
!c-----Compute atmospheric resistance, RA
!c
      ra = (alog(deltaz/z0) - psih)/(vk*ustar)
      ra = amax1(ra,rmin)
!c
!c-----Compute the deposition layer resistance, RD
!c
      sc23 = schmidt**(-2./3.)
      power = -3./stokes
      if (power.lt.-37.) then
        xinert = 10.**(-37.)
      else
        xinert = 10.**(power)
      endif
      rd = 1./(ustar*(sc23 + xinert))
      rd = amax1(rd,rmin)
!c-----Final deposition velocity for this cell, land use, and aerosol size
!c
      vd = vsed + 1./(ra + rd + ra*rd*vsed)
!c
      return
      end
      
      
            

