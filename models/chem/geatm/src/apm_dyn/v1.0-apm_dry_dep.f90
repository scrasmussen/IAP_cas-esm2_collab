subroutine apm_dry_dep_v1 &
 & ( myid &
 &  ,lapm &
 &  ,dt &
 &  ,imonth2 &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dz &
 &  ,Plev,temp,u_ws,v_ws,QVAPOR,hlayer,hgt1 &
 &  ,land_use,temp2m,LATITCRS,LONGICRS,SWDOWN,PBL_HGT &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

real :: dt

integer :: imonth2

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is,imode

integer :: mem2d,mem3d

integer :: ixy,i02,i03,iapm,iapmode

integer :: ip3mem(nzz,nest),ip2mem(nest)

real :: dryvel=0.0 ! need to be modified


real,dimension(mem3d) :: dz,Plev,temp,QVAPOR,u_ws,v_ws,hlayer
real,dimension(mem2d) :: HGT1,temp2m,LATITCRS,LONGICRS,SWDOWN,PBL_HGT

real,dimension(mem2d) :: land_use

integer :: iwb,ieb,jsb,jeb
real   ,allocatable,dimension(:,:)   :: tsurf,xlat,xlon,SWDOWN1,pblhgt !,vdep
integer,allocatable,dimension(:,:)   :: land
real   ,allocatable,dimension(:,:)   :: vdep,gas_vdep
real   ,allocatable,dimension(:,:,:) :: ppp,ttn,uu,vv,height1,QQ,water
real   ,allocatable,dimension(:,:)   :: diam2d,rhop2d

real,parameter :: rgfmin=1.0,rgfmax=100.0
real :: rgf_tmp
real :: diam

real,dimension(mem3d) :: wk,wks
real,dimension(mem2d) :: area0,area,ratio_area

real*8,allocatable,dimension(:) :: numcc,tmp,radius
real*8 :: totalbc,totaloc

! 1:FF 2:BB
integer,parameter :: ridx(1:8) = (/1,2,1,2,1,2,1,2/)
!integer,parameter :: bcidx(1:4) = (/1,2,5,6/)
!integer,parameter :: ocidx(1:4) = (/3,4,7,8/)

character*2,parameter :: flagbcoc(1:8) = (/'bc','bc','oc','oc' &
                                          ,'bc','bc','oc','oc'/)

real    :: vdepm
integer :: imlo,jmlo


!+++ sulferic acid vapor ++++!
logical :: lrddrt,lrdsno     !
integer :: icddrt,icdsno     !
!++++++++++++++++++++++++++++!



!IF(lapm) THEN ! apm flag
!print*,'dims in sub-ddep',sx(ne),ex(ne),sy(ne),ey(ne),nzz

iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1

allocate( tsurf(iwb:ieb,jsb:jeb) )
allocate( xlat(iwb:ieb,jsb:jeb) )
allocate( xlon(iwb:ieb,jsb:jeb) )
allocate( SWDOWN1(iwb:ieb,jsb:jeb) )
allocate( pblhgt(iwb:ieb,jsb:jeb) )

allocate( land(iwb:ieb,jsb:jeb) )

allocate( ppp(iwb:ieb,jsb:jeb,nzz) )
allocate( ttn(iwb:ieb,jsb:jeb,nzz) )
allocate( uu(iwb:ieb,jsb:jeb,nzz) )
allocate( vv(iwb:ieb,jsb:jeb,nzz) )
allocate( height1(iwb:ieb,jsb:jeb,nzz) )
allocate( QQ(iwb:ieb,jsb:jeb,nzz) )
allocate( water(iwb:ieb,jsb:jeb,nzz) )

allocate( diam2d(iwb:ieb,jsb:jeb) )
allocate( rhop2d(iwb:ieb,jsb:jeb) )

allocate( vdep(sx(ne):ex(ne),sy(ne):ey(ne)) )

allocate( gas_vdep(iwb:ieb,jsb:jeb) )


! get local variables

do j=sy(ne)-1,ey(ne)+1
do i=sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      ppp(i,j,k)=Plev(i03+ixy)
      ttn(i,j,k)=temp(i03+ixy)
      QQ(i,j,k)=QVAPOR(i03+ixy)
      uu(i,j,k)=u_ws(i03+ixy)
      vv(i,j,k)=v_ws(i03+ixy)
      water(i,j,k)=QVAPOR(i03+ixy)*1.E06*29./18.!ppm
      if(k==1) then
        height1(i,j,k) = 2*(hlayer(i03+ixy)*1000.-HGT1(i02+ixy))
      endif
      if(k/=1) then
        height1(i,j,k) = 2*(hlayer(i03+ixy)*1000.-HGT1(i02+ixy))- &
                       & height1(i,j,k-1)
      endif
   enddo    !k     

   CALL getland(LAND_USE(i02+ixy),land(i,j))

   tsurf(i,j)= temp2m(i02+ixy)
   xlat(i,j) = LATITCRS(i02+ixy)
   xlon(i,j) = LONGICRS(i02+ixy)
   SWDOWN1(i,j)=SWDOWN(i02+ixy)
   pblhgt(i,j) = PBL_HGT(i02+ixy)
enddo
enddo


!print*,'tsurf :',tsurf
!print*,'xlat :',xlat
!print*,'xlon :',xlon
!print*,size(land)
!print*,'land :',land

!print*,'pblhgt:',pblhgt

!print*,'height1',height1(:,:,1)
!print*,'ppp',ppp(:,:,1)
!print*,'uu',uu(:,:,1)
!print*,'vv',vv(:,:,1)
!print*,'ttn',ttn(:,:,1)


!return

!stop 'kkinputs in ddep'

!========================================================================
!> sulferic acid gas
 if(lfor_h2so4) then
   lrddrt=.false.;icddrt=0  !if droughtness index
   lrdsno=.false.;icdsno=0  !if snow index

   CALL drydep_gas(myid,imonth2,tsurf,xlat,xlon,QQ,QQ,height1,ppp,uu,vv &
               ,SWDOWN1,land,water,ttn,pblhgt,lrddrt,icddrt,lrdsno,icdsno &
               ,gas_vdep &
               ,sx(ne),ex(ne),sy(ne),ey(ne),nzz,1)

   i02=ip2mem(ne)
   do j=sy(ne)-1,ey(ne)+1
   do  i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      apmdryvel2d(i02+ixy) = gas_vdep(i,j)
   enddo
   enddo
  
   i03=ip3mem(1,ne) ! k=1
   call dry_dep_gas_zhu( myid, h2so4_gas(i03),apmdryvel2d(i02),dz(i03) &
                        ,sx(ne), ex(ne), sy(ne), ey(ne), dt )
 endif
!=======================================================================



!=======================================================================
!----------------------------- particles -------------------------------

  k=1

  i02=ip2mem(ne)
  i03=ip3mem(k,ne)

 !> shun : apm sulfate

  if(lfor_sulf) then
!print*,'sulf ddep'
   loop_so4 : do is=1,NSO4

!print*,'size sulf',is
     iapm=ip_sulf(k,is,ne)

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        if(lapm_wetsize) then
          diam2d(i,j)=rd_sulf(is)*2.0*rgf_sulf(i03+ixy)
          call check_rgf('rgf_sulf',i,j,k,rgfmin,rgfmax,rgf_sulf(i03+ixy))
        else
          diam2d(i,j)=rd_sulf(is)*2.0
        endif
!print*,'i&j&',i,j
!print*,'diam=',diam2d(i,j),rd_sulf(is)*2.0,rgf_sulf(i03+ixy)
        !diam2d(i,j)=rw_sulf(iapm+ixy)*2.0
        call check_value(i,j,k,0.0,1000.0,diam2d(i,j))
        rhop2d(i,j)=densulf*1.0e6 ! g/cm3 -> g/m3
     enddo
     enddo

     call drydep_vel_apm( myid,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv &
             &             ,land,ttn,pblhgt,diam2d,rhop2d,vdep &
             &             ,sx(ne),ex(ne),sy(ne),ey(ne),nzz )


     call max_ij_location(myid,sx(ne),ex(ne),sy(ne),ey(ne),imlo,jmlo,vdep,vdepm)
!print*,'max  vdep',vdepm

!print*,'sulf vdep',vdep


     if(vdepm.gt.10.0.and..false.) then
       print*,'nest=',ne
       print*,'imonth2=',imonth2
       do j=sy(ne),ey(ne)
       do i=sx(ne),ex(ne)
         if(vdep(i,j).gt.10) then
           print*,'checkvdel',i,j,vdep(i,j)
           print*,'temp=',tsurf(i,j)
           print*,'xlat=',xlat(i,j)
           print*,'xlon=',xlon(i,j)
           print*,'height1=',height1(i,j,:)
           print*,'uu=',uu(i,j,:)
           print*,'vv=',vv(i,j,:)
           print*,'land=',land(i,j)
           print*,'ttn=',ttn(i,j,:)
           print*,'pblhgt=',pblhgt(i,j)
           print*,'diam2d=',diam2d(i,j)
           print*,'rhop2d=',rhop2d(i,j)
           print*
         endif
       enddo
       enddo
       print*,'warning : vdep too large'
     endif



     !vdep=10.0

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_sulf(k,is,ne)
        apmdryvel2d(i02+ixy)=vdep(i,j)
        !print*,vdep(i,j),j,i
        !apmdryvel2d(i02+ixy)=0.1
        wk(i03+ixy)=apm_sulf(iapm+ixy)
     enddo
     enddo

     call dry_dep_gas_zhu( myid, wk(i03),apmdryvel2d(i02),dz(i03) &
             &            ,sx(ne), ex(ne), sy(ne), ey(ne), dt )


     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_sulf(k,is,ne)
        apm_sulf(iapm+ixy)=wk(i03+ixy)
     enddo
     enddo

   enddo loop_so4
  endif
 !< shun : end of apm sulfate 

!stop 'kksulf_ddep'


!deallocate(apmdryvel2d)


 !> shun : apm seasalt

if(lcoated_dyn) then

   if(.not.allocated(apm_msltsulf)) then
      allocate(apm_msltsulf(saltmem))
   else
      stop 'allocate(apm_msltsulf) erro'
   endif

   ! distribute bulk seasalt coated sulfate to bins
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius))) then
      allocate(numcc(NSEA))
      allocate(tmp(NSEA))
      allocate(radius(NSEA))
   else
      stop 'allocate(numcc) erro'
   endif
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NSEA
        iapm=ip_salt(k,is,ne)
        radius(is)=rd_salt(is)
        numcc(is)= apm_salt(iapm+ixy)/ &
                 & radius(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius(is)**3.0d0)
        tmp(is)=numcc(is)*radius(is)**2
      enddo
      do is=1,NSEA
        iapm=ip_salt(k,is,ne)
        if(sum(tmp(:)).gt.0) then
          apm_msltsulf(iapm+ixy)=msltsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
          apm_msltsulf(iapm+ixy)=0
        endif
      enddo
   enddo
   enddo
   !!!!!

endif

  if(lfor_salt) then
!print*,'salt ddep'
   loop_salt : do is=1,NSEA
!print*,'salt size',is
     iapm=ip_salt(k,is,ne)

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        if(lapm_wetsize) then
          diam2d(i,j)=rd_salt(is)*2.0*rgf_salt(i03+ixy)
          call check_rgf('rgf_salt',i,j,k,rgfmin,rgfmax,rgf_salt(i03+ixy))
        else
          diam2d(i,j)=rd_salt(is)*2.0
        endif
        !diam2d(i,j)=rw_salt(iapm+ixy)*2.0
        call check_value( i,j,k,0.0,1000.0,diam2d(i,j) )
        rhop2d(i,j)=densalt*1.0e6 ! g/cm3 -> g/m3
     enddo
     enddo

     call drydep_vel_apm( myid,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv &
             &             ,land,ttn,pblhgt,diam2d,rhop2d,vdep &
             &             ,sx(ne),ex(ne),sy(ne),ey(ne),nzz )


     call max_ij_location(myid,sx(ne),ex(ne),sy(ne),ey(ne),imlo,jmlo,vdep,vdepm)
!     print*,'max location',imlo,jmlo,vdepm

!print*,'salt vdep',vdep


     if(vdepm.gt.10.0.and..false.) then
       print*,'nest=',ne
       print*,'imonth2=',imonth2
       do j=sy(ne),ey(ne)
       do i=sx(ne),ex(ne)
         if(vdep(i,j).gt.10) then
           print*,'checkvdel',i,j,vdep(i,j)
           print*,'temp=',tsurf(i,j)
           print*,'xlat=',xlat(i,j)
           print*,'xlon=',xlon(i,j)
           print*,'height1=',height1(i,j,:)
           print*,'uu=',uu(i,j,:)
           print*,'vv=',vv(i,j,:)
           print*,'land=',land(i,j)
           print*,'ttn=',ttn(i,j,:)
           print*,'pblhgt=',pblhgt(i,j)
           print*,'diam2d=',diam2d(i,j)
           print*,'rhop2d=',rhop2d(i,j)
           print*
         endif
       enddo
       enddo
       print*,'warning : vdep too large'
     endif
!     print*,'saltvel=',vdep

!stop

!cycle

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_salt(k,is,ne)
        apmdryvel2d(i02+ixy)=vdep(i,j)
        !apmdryvel2d(i02+ixy)=0.1
        wk(i03+ixy)=apm_salt(iapm+ixy)
        if(lcoated_dyn) wks(i03+ixy)=apm_msltsulf(iapm+ixy)
     enddo
     enddo

     call dry_dep_gas_zhu( myid, wk(i03),apmdryvel2d(i02),dz(i03) &
             &            ,sx(ne), ex(ne), sy(ne), ey(ne), dt )


     if(lcoated_dyn) then
       call dry_dep_gas_zhu( myid, wks(i03),apmdryvel2d(i02),dz(i03) &
               &            ,sx(ne), ex(ne), sy(ne), ey(ne), dt )
     endif

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_salt(k,is,ne)
        apm_salt(iapm+ixy)=wk(i03+ixy)
        if(lcoated_dyn) apm_msltsulf(iapm+ixy)=wks(i03+ixy)
     enddo
     enddo
   enddo loop_salt
 endif

!stop 'kksaltvel'


if(lcoated_dyn) then
   ! bins seasalt coated sulfate to bulk
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      msltsulf(i03+ixy)=0.0d0
      do is=1,NSEA
        iapm=ip_salt(k,is,ne)
        msltsulf(i03+ixy)=msltsulf(i03+ixy)+apm_msltsulf(iapm+ixy)
      enddo
   enddo
   enddo
   !!!!!

   if(allocated(apm_msltsulf)) then
      deallocate(apm_msltsulf)
   else
      stop 'deallocate(apm_msltsulf) erro'
   endif

   if(allocated(numcc).and.allocated(tmp).and.allocated(radius)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius)
   else
      stop 'deallocate(numcc) erro'
   endif

endif

 !< shun : end of apm seasalt


 !> shun : apm dust

if(lcoated_dyn) then

   if(.not.allocated(apm_mdstsulf)) then
      allocate(apm_mdstsulf(dustmem))
   else
      stop 'allocate(apm_mdstsulf) erro'
   endif
   ! distribute bulk seasalt coated sulfate to bins
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius))) then
      allocate(numcc(NDSTB))
      allocate(tmp(NDSTB))
      allocate(radius(NDSTB))
   else
      stop 'allocate(numcc) erro'
   endif
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NDSTB
        iapm=ip_dust(k,is,ne)
        radius(is)=rd_dust(is)
        numcc(is)= apm_dust(iapm+ixy)/ &
                 & radius(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius(is)**3.0d0)
        tmp(is)=numcc(is)*radius(is)**2
      enddo
      do is=1,NDSTB
        iapm=ip_dust(k,is,ne)
        if(sum(tmp(:)).gt.0) then
          apm_mdstsulf(iapm+ixy)=mdstsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
          apm_mdstsulf(iapm+ixy)=0
        endif
      enddo
   enddo
   enddo
   !!!!!

endif

  if(lfor_dust) then
!   print*,'dust ddep'
   loop_dust : do is=1,NDSTB

     iapm=ip_dust(k,is,ne)
!print*,'dust size',is
     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        if(lapm_wetsize) then
          diam2d(i,j)=rd_dust(is)*2.0*rgf_dust(i03+ixy)
          call check_rgf('rgf_dust',i,j,k,rgfmin,rgfmax,rgf_dust(i03+ixy))
        else
          diam2d(i,j)=rd_dust(is)*2.0
        endif
        !diam2d(i,j)=rw_dust(iapm+ixy)*2.0
        call check_value(i,j,k,0.0,1000.0,diam2d(i,j))
        rhop2d(i,j)=dendust*1.0e6 ! g/cm3 -> g/m3
     enddo
     enddo

     call drydep_vel_apm( myid,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv &
             &             ,land,ttn,pblhgt,diam2d,rhop2d,vdep &
             &             ,sx(ne),ex(ne),sy(ne),ey(ne),nzz )


     call max_ij_location(myid,sx(ne),ex(ne),sy(ne),ey(ne),imlo,jmlo,vdep,vdepm)
!     print*,'max  vdep',imlo,jmlo,vdepm
!print*,'dust vdep',vdep


     if(vdepm.gt.10.0.and..false.) then
       print*,'nest=',ne
       print*,'imonth2=',imonth2
       do j=sy(ne),ey(ne)
       do i=sx(ne),ex(ne)
         if(vdep(i,j).gt.10) then
           print*,'checkvdel',i,j,vdep(i,j)
           print*,'temp=',tsurf(i,j)
           print*,'xlat=',xlat(i,j)
           print*,'xlon=',xlon(i,j)
           print*,'height1=',height1(i,j,:)
           print*,'uu=',uu(i,j,:)
           print*,'vv=',vv(i,j,:)
           print*,'land=',land(i,j)
           print*,'ttn=',ttn(i,j,:)
           print*,'pblhgt=',pblhgt(i,j)
           print*,'diam2d=',diam2d(i,j)
           print*,'rhop2d=',rhop2d(i,j)
           print*
         endif
       enddo
       enddo
       print*,'warning : vdep too large'
     endif

!print*,'dust vel=',vdep
!cycle


     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_dust(k,is,ne)
        apmdryvel2d(i02+ixy)=vdep(i,j)
        !apmdryvel2d(i02+ixy)=0.1
        wk(i03+ixy)=apm_dust(iapm+ixy)
        if(lcoated_dyn) wks(i03+ixy)=apm_mdstsulf(iapm+ixy)
     enddo
     enddo

     call dry_dep_gas_zhu( myid, wk(i03),apmdryvel2d(i02),dz(i03) &
             &            ,sx(ne), ex(ne), sy(ne), ey(ne), dt )

     if(lcoated_dyn) then
       call dry_dep_gas_zhu( myid, wks(i03),apmdryvel2d(i02),dz(i03) &
               &            ,sx(ne), ex(ne), sy(ne), ey(ne), dt )
     endif

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_dust(k,is,ne)
        apm_dust(iapm+ixy)=wk(i03+ixy)
        if(lcoated_dyn) apm_mdstsulf(iapm+ixy)=wks(i03+ixy)
     enddo
     enddo

   enddo loop_dust
  endif


!stop 'kkdustvel'

if(lcoated_dyn) then
   ! bins seasalt coated sulfate to bulk
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      mdstsulf(i03+ixy)=0.0d0
      do is=1,NDSTB
        iapm=ip_dust(k,is,ne)
        mdstsulf(i03+ixy)=mdstsulf(i03+ixy)+apm_mdstsulf(iapm+ixy)
      enddo
   enddo
   enddo
   !!!!!

   if(allocated(apm_mdstsulf)) then
      deallocate(apm_mdstsulf)
   else
      stop 'deallocate(apm_mdstsulf) erro'
   endif

   if(allocated(numcc).and.allocated(tmp).and.allocated(radius)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius)
   else
      stop 'deallocate(numcc) erro'
   endif

endif

 !< shun : end of apm dust


 !> shun : apm bcoc

if(lcoated_dyn) then
   if(.not.allocated(apm_mbcocsulf)) then
      allocate(apm_mbcocsulf(bcocmem))
   else
      stop 'allocate(apm_mbcocsulf) erro'
   endif
   ! distribute bulk bcoc coated sulfate to bins
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius))) then
      allocate(numcc(NBCOCT))
      allocate(tmp(NBCOCT))
      allocate(radius(NBCOCT))
   else
      stop 'allocate(numcc) erro'
   endif
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
        radius(is)=ref1d_bcoc(ridx(is))
        numcc(is)= apm_bcoc(iapm+ixy)/ &
                 & radius(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius(is)**3.0d0)
        tmp(is)=numcc(is)*radius(is)**2
      enddo
      totalbc=tmp(1)+tmp(5)+tmp(2)+tmp(6)
      totaloc=tmp(3)+tmp(7)+tmp(4)+tmp(8)
      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
        if(flagbcoc(is).eq.'bc') then
          if(totalbc.gt.0) then
            apm_mbcocsulf(iapm+ixy)=mbcsulf(i03+ixy)*tmp(is)/totalbc
          else
            apm_mbcocsulf(iapm+ixy)=0
          endif
        elseif(flagbcoc(is).eq.'oc') then
          if(totaloc.gt.0) then
            apm_mbcocsulf(iapm+ixy)=mocsulf(i03+ixy)*tmp(is)/totaloc
          else
            apm_mbcocsulf(iapm+ixy)=0
          endif
        else
          stop 'erro'
        endif
      enddo
   enddo
   enddo

endif



  if(lfor_bcoc) then
!   print*,'bcoc ddep'
   loop_bcoc : do is=1,NBCOCT

!print*,'bcoc categ',is
     iapm = ip_bcoc(k,is,ne)

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        !iapmode=ipmode_bcoc(k,ridx(is),ne)
        if(flagbcoc(is).eq.'bc') then
        !  diam=refw_bc(iapmode+ixy)*2.0
          rgf_tmp=rgf_bc(i03+ixy)
        elseif(flagbcoc(is).eq.'oc') then
        !  diam=refw_oc(iapmode+ixy)*2.0
          rgf_tmp=rgf_oc(i03+ixy)
        endif
        if(lapm_wetsize) then
          diam2d(i,j)=ref1d_bcoc(ridx(is))*2.0*rgf_tmp
          call check_rgf('rgf_bcoc',i,j,k,rgfmin,rgfmax,rgf_tmp)
        else
          diam2d(i,j)=ref1d_bcoc(ridx(is))*2.0
        endif
        call check_value(i,j,k,0.0,1000.0,diam2d(i,j))
        rhop2d(i,j)=denbcoc*1.0e6 ! g/cm3 -> g/m3
     enddo
     enddo

     call drydep_vel_apm( myid,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv &
             &             ,land,ttn,pblhgt,diam2d,rhop2d,vdep &
             &             ,sx(ne),ex(ne),sy(ne),ey(ne),nzz )


     call max_ij_location(myid,sx(ne),ex(ne),sy(ne),ey(ne),imlo,jmlo,vdep,vdepm)
!     print*,'max  vdep',vdepm

!print*,'bcoc vdep',vdep


     if(vdepm.gt.10.0.and..false.) then
       print*,'nest=',ne
       print*,'imonth2=',imonth2
       do j=sy(ne),ey(ne)
       do i=sx(ne),ex(ne)
         if(vdep(i,j).gt.10) then
           print*,'checkvdel',i,j,vdep(i,j)
           print*,'temp=',tsurf(i,j)
           print*,'xlat=',xlat(i,j)
           print*,'xlon=',xlon(i,j)
           print*,'height1=',height1(i,j,:)
           print*,'uu=',uu(i,j,:)
           print*,'vv=',vv(i,j,:)
           print*,'land=',land(i,j)
           print*,'ttn=',ttn(i,j,:)
           print*,'pblhgt=',pblhgt(i,j)
           print*,'diam2d=',diam2d(i,j)
           print*,'rhop2d=',rhop2d(i,j)
           print*
         endif
       enddo
       enddo
       print*,'warning : vdep too large'
     endif

!print*,'bcocvel=',vdep


     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_bcoc(k,is,ne)
        apmdryvel2d(i02+ixy)=vdep(i,j)
        !apmdryvel2d(i02+ixy)=0.1
        wk(i03+ixy)=apm_bcoc(iapm+ixy)
        if(lcoated_dyn) wks(i03+ixy)=apm_mbcocsulf(iapm+ixy)
     enddo
     enddo

     call dry_dep_gas_zhu( myid, wk(i03),apmdryvel2d(i02),dz(i03) &
             &            ,sx(ne), ex(ne), sy(ne), ey(ne), dt )


     if(lcoated_dyn) then
        call dry_dep_gas_zhu( myid, wks(i03),apmdryvel2d(i02),dz(i03) &
                &            ,sx(ne), ex(ne), sy(ne), ey(ne), dt )
     endif

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_bcoc(k,is,ne)
        apm_bcoc(iapm+ixy)=wk(i03+ixy)
        if(lcoated_dyn) apm_mbcocsulf(iapm+ixy)=wks(i03+ixy)
     enddo
     enddo

   enddo loop_bcoc
  endif

!return

!stop 'kkbcocvdep'


if(lcoated_dyn) then
   ! bins bcoc coated sulfate to bulk
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      mbcsulf(i03+ixy)=0.0d0
      mocsulf(i03+ixy)=0.0d0
      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
        if(flagbcoc(is).eq.'bc') then
          mbcsulf(i03+ixy)=mbcsulf(i03+ixy)+apm_mbcocsulf(iapm+ixy)
        elseif(flagbcoc(is).eq.'oc') then
          mocsulf(i03+ixy)=mocsulf(i03+ixy)+apm_mbcocsulf(iapm+ixy)
        else
          stop 'erro'
        endif
      enddo
   enddo
   enddo
   !!!!!
   if(allocated(apm_mbcocsulf)) then
      deallocate(apm_mbcocsulf)
   else
      stop 'deallocate(apm_mbcocsulf) erro'
   endif

   if(allocated(numcc).and.allocated(tmp).and.allocated(radius)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius)
   else
      stop 'deallocate(numcc) erro'
   endif
endif

 !< shun : end of apm bcoc


!ENDIF ! apm flag


!stop 'kkbcocvel'

end subroutine apm_dry_dep_v1



