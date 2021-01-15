
subroutine apm_gra_dep_v2 &
 & ( lapm &
 &  ,myid &
 &  ,dt,numTGRV &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dz &
 &  ,ip2mem &
 &  ,ip3mem,mem3d )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

real    :: dt
integer :: numTGRV

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is,IT


integer :: mem3d

real,dimension(mem3d) :: dz

integer :: ixy,i02,i03,iapm

integer :: ip2mem(nest)
integer :: ip3mem(nzz,nest)

real    :: gvel ! need to be modified


real*8,allocatable,dimension(:) :: numcc,tmp,radius
real*8 :: totalbc,totaloc

real :: diam,rhop

integer,parameter :: ridx(1:8) = (/1,2,1,2,1,2,1,2/)

character*2,parameter :: flagbcoc(1:8) = (/'bc','bc','oc','oc' &
                                          ,'bc','bc','oc','oc'/)


real :: subdt
real :: value1,value2

real,dimension(nzz) :: diam1d,rhop1d,gvel1d,con,deltz

real :: rgf_tmp

real,allocatable,dimension(:,:) :: consp


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 subdt=dt/numTGRV

LOOP_IT : DO IT=1,numTGRV

!=====================
!> shun : apm sulfate
  if(lfor_sulf.and..false.) then
!   print*,'sulf gdep'
   loop_so4 : do is=1,NSO4

!     print*,'sulfate_bin=',is

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

        do k=1,nzz
          i03=ip3mem(k,ne)
          iapm=ip_sulf(k,is,ne)
!apm_sulf(iapm+ixy)=10
          con(k)=apm_sulf(iapm+ixy)
          deltz(k)=dz(i03+ixy)
          if(lapm_wetsize) then
            diam1d(k)=rd_sulf(is)*2.0*rgf_sulf(i03+ixy)
          else
            diam1d(k)=rd_sulf(is)*2.0
          endif    
          rhop1d(k)=densulf*1.0e+6
          diam=diam1d(k)
          rhop=rhop1d(k)
          call vg_aer(diam,rhop,gvel)
!gvel=1.0
          gvel1d(k)=gvel
        enddo


        if(i.eq.33.and.j.eq.33.and.is.eq.40) then
!           print*,'gvel='
!           print*,gvel1d
!           print*,'sulf-con'
!           print*,con
        endif

        do k=1,nzz
          if(k==nzz) then
            value2=MAX(0.0,MIN(0.99,gvel1d(k)*subdt/deltz(k)))
            con(k)=con(k)-con(k)*value2
          else
            value1=MAX(0.0,MIN(0.99,gvel1d(k+1)*subdt/deltz(k+1))) ! shun
            value2=MAX(0.0,MIN(0.99,gvel1d(k)*subdt/deltz(k)))       ! shun
            if(k.ne.1) then
              con(k)=con(k)+con(k+1)*value1-con(k)*value2
            elseif(k.eq.1) then
              con(k)=con(k)+con(k+1)*value1 ! shun
            endif
          endif
        enddo

        if(i.eq.33.and.j.eq.33.and.is.eq.40) then
!           print*,'result'
!           print*,con
        endif


        do k=1,nzz
          iapm=ip_sulf(k,is,ne)
          apm_sulf(iapm+ixy)=con(k)
        enddo

     enddo
     enddo

    enddo loop_so4
   endif
!===================

!stop

!===================
!> shun : apm salt

  if(lfor_salt) then

!    print*,'salt gdep'

    if(lcoated_dyn) then
     allocate(consp(NSEA,nzz))
     allocate(numcc(NSEA))
     allocate(tmp(NSEA))
     allocate(radius(NSEA))
    endif

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
 
       if(lcoated_dyn) then
        do k=1,nzz
          i03=ip3mem(k,ne)
          do is=1,NSEA
            iapm=ip_salt(k,is,ne)
            radius(is)=rd_salt(is)
!apm_salt(iapm+ixy)=10
            numcc(is)= apm_salt(iapm+ixy)/ &
                 & radius(is)**3.0d0
            tmp(is)=numcc(is)*radius(is)**2
!tmp(is)=1
          enddo
          do is=1,NSEA
            iapm=ip_salt(k,is,ne)
!msltsulf(i03+ixy)=200
            if(sum(tmp(:)).gt.0) then
              consp(is,k)=msltsulf(i03+ixy)*tmp(is)/sum(tmp(:))
            else
              consp(is,k)=0
            endif
          enddo
        enddo
       endif

!consp=10
        loop_salt : do is=1,NSEA ! salt and its coating gdep
!        print*,'salt_bin=',is
        do k=1,nzz
          i03=ip3mem(k,ne)
          iapm=ip_salt(k,is,ne)
          con(k)=apm_salt(iapm+ixy)
!          con(k)=10.0
          deltz(k)=dz(i03+ixy)
          if(lapm_wetsize) then
            diam1d(k)=rd_salt(is)*2.0*rgf_salt(i03+ixy)
          else
            diam1d(k)=rd_salt(is)*2.0
          endif
          rhop1d(k)=densalt*1.0e+6
          diam=diam1d(k)
          rhop=rhop1d(k)
          call vg_aer(diam,rhop,gvel)
          gvel1d(k)=gvel
        enddo

        if(i.eq.33.and.j.eq.33.and.is.eq.20.and..false.) then
           print*,'gvel='
           print*,gvel1d
           print*,'sea-con'
           print*,con
        endif

        do k=1,nzz
          if(k==nzz) then
            value2=MAX(0.0,MIN(0.99,gvel1d(k)*subdt/deltz(k)))
            con(k)=con(k)-con(k)*value2
            if(lcoated_dyn) then
               consp(is,k)=consp(is,k)-consp(is,k)*value2
            endif
          else
            value1=MAX(0.0,MIN(0.99,gvel1d(k+1)*subdt/deltz(k+1))) ! shun
            value2=MAX(0.0,MIN(0.99,gvel1d(k)*subdt/deltz(k)))       ! shun
            if(k.ne.1) then
              con(k)=con(k)+con(k+1)*value1-con(k)*value2
              if(lcoated_dyn) then
                consp(is,k)=consp(is,k)+consp(is,k+1)*value1-consp(is,k)*value2 
              endif
            elseif(k.eq.1) then
              con(k)=con(k)+con(k+1)*value1 ! shun
              if(lcoated_dyn) then
                consp(is,k)=consp(is,k)+consp(is,k+1)*value1
              endif
            endif
          endif
        enddo

        if(i.eq.33.and.j.eq.33.and.is.eq.20.and..false.) then
           print*,'sea-dep'
           print*,con
           print*,'sea-coating'
           print*,consp(is,:)
        endif

        do k=1,nzz
          iapm=ip_salt(k,is,ne)
          apm_salt(iapm+ixy)=con(k)
        enddo

        enddo loop_salt ! salt and its coating gdep

        if(lcoated_dyn) then
        do k=1,nzz
          i03=ip3mem(k,ne) 
          msltsulf(i03+ixy)=sum(consp(:,k))
        enddo
        endif

        if(i.eq.33.and.j.eq.33.and..false.) then
          print*,'sea-coating'
          do k=1,20
          print*,sum(consp(:,k))
          enddo
          print*,'check',sum(consp(:,:)),sum(con(:))
        endif

     enddo
     enddo

     if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

  endif
!===================

!stop

!===================
!> shun : apm dust
if(lfor_dust) then

!     print*,'dust_gdep'

    if(lcoated_dyn) then
     allocate(consp(NDSTB,nzz))
     allocate(numcc(NDSTB))
     allocate(tmp(NDSTB))
     allocate(radius(NDSTB))
    endif

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

       if(lcoated_dyn) then
        do k=1,nzz
          i03=ip3mem(k,ne)
          do is=1,NDSTB
            iapm=ip_dust(k,is,ne)
            radius(is)=rd_dust(is)
            numcc(is)= apm_dust(iapm+ixy)/ &
                 & radius(is)**3.0d0
            tmp(is)=numcc(is)*radius(is)**2
          enddo
          do is=1,NDSTB
            iapm=ip_dust(k,is,ne)
            if(sum(tmp(:)).gt.0) then
              consp(is,k)=mdstsulf(i03+ixy)*tmp(is)/sum(tmp(:))
            else
              consp(is,k)=0
            endif
          enddo
        enddo
       endif

        loop_dst : do is=1,NDSTB

        do k=1,nzz
          i03=ip3mem(k,ne)
          iapm=ip_dust(k,is,ne)
          con(k)=apm_dust(iapm+ixy)
!          con(k)=10.0
          deltz(k)=dz(i03+ixy)
          if(lapm_wetsize) then
            diam1d(k)=rd_dust(is)*2.0*rgf_dust(i03+ixy)
          else
            diam1d(k)=rd_dust(is)*2.0
          endif   
          rhop1d(k)=dendust*1.0e+6
          diam=diam1d(k)
          rhop=rhop1d(k)
          call vg_aer(diam,rhop,gvel)
          gvel1d(k)=gvel
        enddo

        if(i.eq.100.and.j.eq.90) then
!           print*,con
        endif

        do k=1,nzz
          if(k==nzz) then
            value2=MAX(0.0,MIN(0.99,gvel1d(k)*subdt/deltz(k)))
            con(k)=con(k)-con(k)*value2
            if(lcoated_dyn) then
               consp(is,k)=consp(is,k)-consp(is,k)*value2
            endif
          else
            value1=MAX(0.0,MIN(0.99,gvel1d(k+1)*subdt/deltz(k+1))) ! shun
            value2=MAX(0.0,MIN(0.99,gvel1d(k)*subdt/deltz(k)))       ! shun
            if(k.ne.1) then
              con(k)=con(k)+con(k+1)*value1-con(k)*value2
              if(lcoated_dyn) then
                consp(is,k)=consp(is,k)+consp(is,k+1)*value1-consp(is,k)*value2
              endif
            elseif(k.eq.1) then
              con(k)=con(k)+con(k+1)*value1 ! shun
              if(lcoated_dyn) then
                consp(is,k)=consp(is,k)+consp(is,k+1)*value1
              endif
            endif
          endif
        enddo

        if(i.eq.100.and.j.eq.90) then
!           print*,con
        endif

        do k=1,nzz
          iapm=ip_dust(k,is,ne)
          apm_dust(iapm+ixy)=con(k)
        enddo

        enddo loop_dst

        if(lcoated_dyn) then
        do k=1,nzz
          i03=ip3mem(k,ne)
          mdstsulf(i03+ixy)=sum(consp(:,k))
        enddo
        endif

     enddo
     enddo

     if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

   endif
!==================


!==================
! shun : apm bcoc

if(lfor_bcoc.and..false.) then

!     print*,'bcoc_gdep'

    if(lcoated_dyn) then
     allocate(consp(NBCOCT,nzz))
     allocate(numcc(NBCOCT))
     allocate(tmp(NBCOCT))
     allocate(radius(NBCOCT))
    endif

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

       if(lcoated_dyn) then
        do k=1,nzz

          i03=ip3mem(k,ne)

          do is=1,NBCOCT
            iapm=ip_bcoc(k,is,ne)
            radius(is)=ref1d_bcoc(ridx(is))
            numcc(is)= apm_bcoc(iapm+ixy)/ &
                 & radius(is)**3.0d0
            tmp(is)=numcc(is)*radius(is)**2
          enddo
          totalbc=tmp(1)+tmp(5)+tmp(2)+tmp(6)
          totaloc=tmp(3)+tmp(7)+tmp(4)+tmp(8)
          do is=1,NBCOCT
            iapm=ip_bcoc(k,is,ne)
            if(flagbcoc(is).eq.'bc') then
              if(totalbc.gt.0) then
                consp(is,k)=mbcsulf(i03+ixy)*tmp(is)/totalbc
              else
                consp(is,k)=0
              endif
            elseif(flagbcoc(is).eq.'oc') then
              if(totaloc.gt.0) then
                 consp(is,k)=mocsulf(i03+ixy)*tmp(is)/totaloc
              else
                 consp(is,k)=0
              endif
            endif
          enddo

        enddo ! k
       endif

        loop_bcoc : do is=1,NBCOCT
  
        do k=1,nzz
          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)
          con(k)=apm_bcoc(iapm+ixy)
!          con(k)=10.0
          deltz(k)=dz(i03+ixy)
          if(lapm_wetsize) then
            if(flagbcoc(is).eq.'bc') then
               rgf_tmp=rgf_bc(i03+ixy)
            elseif(flagbcoc(is).eq.'oc') then
               rgf_tmp=rgf_oc(i03+ixy)
            endif
          else
            rgf_tmp=1.0
          endif
          diam1d(k)=ref1d_bcoc(ridx(is))*2.0*rgf_tmp
          rhop1d(k)=denbcoc*1.0e+6
          diam=diam1d(k)
          rhop=rhop1d(k)
          call vg_aer(diam,rhop,gvel)
          gvel1d(k)=gvel
        enddo

        if(i.eq.100.and.j.eq.90) then
!           print*,con
        endif

        do k=1,nzz
          if(k==nzz) then
            value2=MAX(0.0,MIN(0.99,gvel1d(k)*subdt/deltz(k)))
            con(k)=con(k)-con(k)*value2
            if(lcoated_dyn) then
               consp(is,k)=consp(is,k)-consp(is,k)*value2
            endif
          else
            value1=MAX(0.0,MIN(0.99,gvel1d(k+1)*subdt/deltz(k+1))) ! shun
            value2=MAX(0.0,MIN(0.99,gvel1d(k)*subdt/deltz(k)))       ! shun
            if(k.ne.1) then
              con(k)=con(k)+con(k+1)*value1-con(k)*value2
              if(lcoated_dyn) then
                consp(is,k)=consp(is,k)+consp(is,k+1)*value1-consp(is,k)*value2
              endif
            elseif(k.eq.1) then
              con(k)=con(k)+con(k+1)*value1 ! shun
              if(lcoated_dyn) then
                consp(is,k)=consp(is,k)+consp(is,k+1)*value1
              endif
            endif
          endif
        enddo

        if(i.eq.100.and.j.eq.90) then
!           print*,con
        endif


        do k=1,nzz
          iapm=ip_bcoc(k,is,ne)
          apm_bcoc(iapm+ixy)=con(k)
        enddo

        enddo loop_bcoc

       if(lcoated_dyn) then
        do k=1,nzz
           i03=ip3mem(k,ne)
           mbcsulf(i03+ixy)=consp(1,k)+consp(5,k)+consp(2,k)+consp(6,k)
           mocsulf(i03+ixy)=consp(3,k)+consp(7,k)+consp(4,k)+consp(8,k)
        enddo
       endif
!================

     enddo
     enddo

     if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

   endif
!====================

ENDDO LOOP_IT
 


end subroutine apm_gra_dep_v2





