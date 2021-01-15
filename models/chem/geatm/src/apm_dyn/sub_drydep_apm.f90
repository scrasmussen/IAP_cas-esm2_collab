      subroutine drydep_vel_apm(myid,month,tsurf,cellat,cellon,height,&
                            press,windu,windv,landuse,tempk,&
                    pbl_hgt,diam2d,rhop2d,vdep_aer,sx,ex,sy,ey,nzz)
      integer :: sx,ex,sy,ey
      real tsurf(sx-1:ex+1,sy-1:ey+1),cellat(sx-1:ex+1,sy-1:ey+1),&
          cellon(sx-1:ex+1,sy-1:ey+1)

      real diam2d(sx-1:ex+1,sy-1:ey+1),rhop2d(sx-1:ex+1,sy-1:ey+1)

      real height(sx-1:ex+1,sy-1:ey+1,nzz),&     
          press(sx-1:ex+1,sy-1:ey+1,nzz),&
          windu(sx-1:ex+1,sy-1:ey+1,nzz),&
          windv(sx-1:ex+1,sy-1:ey+1,nzz),&
          tempk(sx-1:ex+1,sy-1:ey+1,nzz)
      !real vdep_aer(sx-1:ex+1,sy-1:ey+1),pbl_hgt(sx-1:ex+1,sy-1:ey+1)
      real vdep_aer(sx:ex,sy:ey),pbl_hgt(sx-1:ex+1,sy-1:ey+1)


      logical lstable
      REAL :: iseason(5,12)
      REAL :: diam  
! diam : log-mean sectional aerosol diameter (m)
      REAL :: rhop ! AEROSOL DENISITY in g/m3
      INTEGER ::  landuse(sx-1:ex+1,sy-1:ey+1)
      REAL,DIMENSION(11,5) :: z0lu     
      INTEGER :: IA ! 2 DUST 1 SEASALT
!
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
      data iseason / 1, 3, 3, 3, 3,& ! Jan
                    1, 5, 3, 3, 3,& ! Feb
                    1, 5, 5, 3, 3,& ! Mar
                    1, 5, 5, 5, 3,& ! Apr
                    1, 1, 5, 5, 3,& ! May
                    1, 1, 1, 1, 5,& ! Jun
                    1, 1, 1, 1, 1,& ! Jul
                    1, 1, 1, 1, 2,& ! Aug
                    1, 1, 2, 2, 3,& ! Sep
                    1, 2, 2, 2, 3,& ! Oct
                    1, 2, 2, 3, 3,& ! Nov
                    1, 2, 3, 3, 3/ ! Dec
!c
!c-----Surface roughness (m) as a function of 11 landuse categories
!c     and 5 seasons; based on AERMET model (ref EPA SCRAM website)
!c
      data z0lu&
      /1.0,0.20,0.100,1.3,1.3,1.30,0.0001,0.002,0.20,0.150,0.30,&
       1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,&
       1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,&
       1.0,0.01,0.001,0.5,1.3,0.90,0.0001,0.002,0.05,0.006,0.15,&
       1.0,0.03,0.050,1.0,1.3,1.15,0.0001,0.002,0.20,0.040,0.30/

       !IF(IA==1) rhop = 2.2E06 ! SEASALT
       !IF(IA==2) rhop = 2.65E06  ! DUST

       !if(IA.eq.9999) rhop = 1.7E06 ! shun apm sulf

      do 30 j = sy,ey
        do 20 i = sx,ex
!c
!c-----Determine season
!c
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
                   
        m=landuse(i,j)
        z0 = z0lu(m,isesn)
        if (m.eq.7) z0 = amax1(z0,2.0e-6*wind**2.5)
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
         call micromet(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,&
                    deltaz,wind,z0,pbl,ustar,psih,wstar,lstable)
!cc

         diam=diam2d(i,j) ! shun add
         rhop=rhop2d(i,j) ! shun add

         call vd_dust(z0,deltaz,psih,ustar,diam,rhop,temp0,vd)      

         vdep_aer(i,j) =  vd
        
   20     continue
   30   continue
    
                            
      return
      end  subroutine drydep_vel_apm 

      
      
            

