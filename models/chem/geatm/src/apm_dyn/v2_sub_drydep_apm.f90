      subroutine drydep_vel_apm_v2(myid,cflag,month &
                            ,tsurf,cellat,cellon,height,press,windu,windv &
                            ,lrdlai,wrflai &
                            ,landuse,tempk,pbl_hgt,diam2d,rhop2d,vdep_aer &
                            ,ix,iy,nzz)
      use dep_parm

      character(len=*) :: cflag

      integer :: nzz
      integer :: ix,iy
      real tsurf,cellat,cellon

      real diam2d,rhop2d

      real height(nzz),press(nzz),windu(nzz),windv(nzz),tempk(nzz)
      real vdep_aer,pbl_hgt


      logical lstable
!      REAL :: iseason(5,12)
      REAL :: diam  
! diam : log-mean sectional aerosol diameter (m)
      REAL :: rhop ! AEROSOL DENISITY in g/m3
      INTEGER ::  landuse
!      REAL,DIMENSION(11,5) :: z0lu     
      INTEGER :: IA ! 2 DUST 1 SEASALT
!
      integer, dimension(12) :: nday
      real :: zl,el,lai_f

      logical :: lrdlai
      real    :: wrflai

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





       !if(IA.eq.9999) rhop = 1.7E06 ! shun apm sulf


        if(ix.eq.100.and.iy.eq.90.and..false.) then
           print*,'tsurf=',tsurf
           print*,'cellon&cellat=',cellon,cellat
           print*,'diam2d,rhop2d=',diam2d,rhop2d
           print*,'height=',height
           print*,'press=',press
           print*,'windu=',windu
           print*,'windv=',windv
           print*,'tempk=',tempk
        endif



!c
!c-----Determine season
!c


          mbin = month
          if (cellat.lt.0.) then
            mbin = mod(month+6,12)
            if (mbin.eq.0) mbin = 12
          endif
          latbin = 1
          if (abs(cellat).gt.20.) then
          latbin = 2
          elseif (abs(cellat).gt.35.) then
          latbin = 3
          elseif (abs(cellat).gt.50.) then
          latbin = 4
          elseif (abs(cellat).gt.75.) then
          latbin = 5
          endif
          if ((cellat.gt.50. .and. cellat.lt.75.) .and.&
             (cellon.gt.-15. .and. cellon.lt.15.)) latbin = 3
           isesn = iseason(latbin,mbin)
!c
!c-----Use input snow cover to set season, if specified
!c
!c
!c-----Load local met variables
          deltaz = height(1)/2.
          temp0 = tsurf - 273.15
          prss0 = press(1) -&
                 2.*deltaz*(press(2) - press(1))/height(2)
          ucomp = (windu(1) + windu(1))/2.
          vcomp = (windv(1) + windv(1))/2.
          wind = sqrt(ucomp**2 + vcomp**2)
          wind = amax1(0.1,wind)
!c
!c
!c-----Loop over land use; surface roughness for water is dependent on        
!c     wind speed
         vdep_aer = 0.
         totland = 1.


        if(trim(cflag).eq.'W89') then                   
           m=landuse
           z0 = z0lu(m,isesn)
           if (m.eq.7) z0 = amax1(z0,2.0e-6*wind**2.5)

        elseif(trim(cflag).eq.'Z03') then

           m=landuse

            if(lrdlai) then
             lai_f = wrflai
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
           stop 'apm-ddep-flag erro in v2_sub_drydep_apm.F90'
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
           write(iout,*) i,j,height(1),deltaz,press(1),press(2),prss0
         endif

         pbl = pbl_hgt
!         call micromet(tempk(1),tsurf,press(1),prss0,&
!                    deltaz,wind,z0,pbl,ustar,psih,wstar,lstable)
         call micromet_el(tempk(1),tsurf,press(1),prss0,&
                    deltaz,wind,z0,pbl,ustar,el,psih,wstar,lstable)
         zl = deltaz/el
         if (zl.lt.0.) then
           zl = min(-5.,zl)
         else
           zl = max(5.,zl)
         endif



         diam=diam2d ! shun add
         rhop=rhop2d ! shun add

         
         if(trim(cflag).eq.'W89') then

!            call vd_dust(z0,deltaz,psih,ustar,diam,rhop,temp0,vd)      
            call vd_aer(z0,deltaz,psih,ustar,diam,rhop,tsurf,vd)

         elseif(trim(cflag).eq.'Z03') then

            call vd_aer_zhang(deltaz,zl,z0,ustar,diam,rhop &
                         ,tsurf,tempk(1),m,lai_f,vd)

            if(.not.(vd.ge.0.and.vd.le.1.0e3)) then
               print*,'i&j&ig=',ig,i,j,vd
               print*,deltaz,zl,z0,ustar
               print*,diam,rhop
               print*,tsurf,tempk(1)
               print*,m,lai_f
               stop
            endif

         endif

         vdep_aer =  vd
        
    
                            
      return
      end  subroutine drydep_vel_apm_v2

      
      
            

