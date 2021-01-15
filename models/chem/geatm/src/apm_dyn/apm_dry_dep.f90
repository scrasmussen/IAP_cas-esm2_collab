
subroutine apm_dry_dep &
 & ( myid &
 &  ,lapm &
 &  ,dt &
 &  ,imonth2 &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,iwb,ieb,jsb,jeb &
 &  ,xlat,xlon,land,tsurf,SWDOWN1,pblhgt &
 &  ,ppp,ttn,uu,vv,height1 &
 &  ,dz &
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

real,dimension(mem3d) :: dz

integer :: ixy,i02,i03,iapm,iapmode

integer :: ip3mem(nzz,nest),ip2mem(nest)

real :: dryvel=0.0 ! need to be modified

integer :: iwb,ieb,jsb,jeb
real,dimension(iwb:ieb,jsb:jeb) :: tsurf,xlat,xlon,SWDOWN1,pblhgt !,vdep
integer,dimension(iwb:ieb,jsb:jeb) :: land


real,dimension(sx(ne):ex(ne),sy(ne):ey(ne)) :: vdep

real,dimension(iwb:ieb,jsb:jeb,nzz) :: ppp,ttn,uu,vv,height1

real,dimension(iwb:ieb,jsb:jeb) :: diam2d,rhop2d

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



real  :: vdepm
integer :: imlo,jmlo


!IF(lapm) THEN ! apm flag
print*,'dims in sub',sx(ne),ex(ne),sy(ne),ey(ne),nzz
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

!stop 'kkradius in ddep'


  k=1

  i02=ip2mem(ne)
  i03=ip3mem(k,ne)

 !> shun : apm sulfate

  if(lfor_sulf) then
   print*,'sulf ddep'
   loop_so4 : do is=1,NSO4

!print*,'size sulf',is
     iapm=ip_sulf(k,is,ne)

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        if(lapm_wetsize) then
          diam2d(i,j)=rd_sulf(is)*2.0*rgf_sulf(i03+ixy)
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

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
!       call check_value(i,j,k,-1.0e-4,1000.0,vdep(i,j))
     enddo
     enddo

     call max_ij_location(iwb,ieb,jsb,jeb,imlo,jmlo,vdep,vdepm)
     print*,'max location',imlo,jmlo,vdepm

!     if(vdepm.gt.10) then
!       print*,'nest=',ne
!       print*,'vdep=',vdep
!     endif


!     vdep=10.0

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_sulf(k,is,ne)
        apmdryvel2d(i02+ixy)=vdep(i,j)
        !print*,vdep(i,j),j,i
        !apmdryvel2d(i02+ixy)=2.0
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
   print*,'salt ddep'
   loop_salt : do is=1,NSEA
     iapm=ip_salt(k,is,ne)

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        if(lapm_wetsize) then
          diam2d(i,j)=rd_salt(is)*2.0*rgf_salt(i03+ixy)
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


     call max_ij_location(iwb,ieb,jsb,jeb,imlo,jmlo,vdep,vdepm)
     print*,'max location',imlo,jmlo,vdepm

!     if(vdepm.gt.10) then
!       print*,'nest=',ne
!       print*,'vdep=',vdep
!     endif

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
!       call check_value(i,j,k,-1.0e-4,1000.0,vdep(i,j))
     enddo
     enddo

!     print*,'salt size',is
!     print*,'salt vel=',vdep

!cycle

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_salt(k,is,ne)
        apmdryvel2d(i02+ixy)=vdep(i,j)
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
   print*,'dust ddep'
   loop_dust : do is=1,NDSTB

     iapm=ip_dust(k,is,ne)
print*,'dust size',is
     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        if(lapm_wetsize) then
          diam2d(i,j)=rd_dust(is)*2.0*rgf_dust(i03+ixy)
        else
          diam2d(i,j)=rd_dust(is)*2.0
        endif
!print*,'i&j',i,j
!print*,diam2d(i,j),rd_dust(is)*2.0,rgf_dust(i03+ixy)
        !diam2d(i,j)=rw_dust(iapm+ixy)*2.0
        call check_value(i,j,k,0.0,1000.0,diam2d(i,j))
        rhop2d(i,j)=dendust*1.0e6 ! g/cm3 -> g/m3
     enddo
     enddo

     call drydep_vel_apm( myid,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv &
             &             ,land,ttn,pblhgt,diam2d,rhop2d,vdep &
             &             ,sx(ne),ex(ne),sy(ne),ey(ne),nzz )


     call max_ij_location(iwb,ieb,jsb,jeb,imlo,jmlo,vdep,vdepm)
     print*,'max location',imlo,jmlo,vdepm

!     if(vdepm.gt.10) then
!       print*,'nest=',ne
!       print*,'vdep=',vdep
!     endif

!print*,'dust vel=',vdep
!cycle

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
!       call check_value(i,j,k,-1.0e-4,1000.0,vdep(i,j))
     enddo
     enddo

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_dust(k,is,ne)
        apmdryvel2d(i02+ixy)=vdep(i,j)
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
   print*,'bcoc ddep'
   loop_bcoc : do is=1,NBCOCT

print*,'bcoc categ',is
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


     call max_ij_location(iwb,ieb,jsb,jeb,imlo,jmlo,vdep,vdepm)
     print*,'max location',imlo,jmlo,vdepm
     print*,vdep(35,64),is


     if(vdepm.gt.10.and..false.) then
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
           print*,'diam2d=',diam2d(i,j),ref1d_bcoc(ridx(is))*2.0
           print*,'rhop2d=',rhop2d(i,j)
           print*
         endif
       enddo
       enddo
     endif


!print*,'bcocvel=',vdep

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
!       call check_value(i,j,k,-1.0e-4,10.0,vdep(i,j))
     enddo
     enddo


     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_bcoc(k,is,ne)
        apmdryvel2d(i02+ixy)=vdep(i,j)
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

end subroutine apm_dry_dep



