
subroutine apm_v_dif &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,iwb,ieb,jsb,jeb &
 &  ,dzz &
 &  ,ip3mem &
 &  ,rkv,ttn,ppp,atm,kktop )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

real :: dt

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is

integer :: ixy,iapm,i03

integer :: ip3mem(nzz,nest)

integer :: iwb,ieb,jsb,jeb
real,dimension(iwb:ieb,jsb:jeb,nzz) :: ppp,ttn,conc,atm,rkv,dzz
real,dimension(iwb:ieb,jsb:jeb    ) :: kktop

integer :: ispflag


!IF(lapm) THEN ! apm flag

 !> shun : apm sulfate vdif
  if(lfor_sulf) then
   !print*,'sulf vdif'
   loop_so4_vdif : do is=1,NSO4

    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       iapm=ip_sulf(k,is,ne)
       conc(i,j,k)=apm_sulf(iapm+ixy)
     enddo
    enddo
    enddo
    ispflag=9999
    call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
               & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        iapm=ip_sulf(k,is,ne)
        apm_sulf(iapm+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo loop_so4_vdif
  endif
 !< shun : end of apm sulfate vdif


 !> shun : apm seasalt vdif
  if(lfor_salt) then
   !print*,'salt vdif'
   loop_salt_vdif : do is=1,NSEA
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       iapm=ip_salt(k,is,ne)
       conc(i,j,k)=apm_salt(iapm+ixy)
     enddo
    enddo
    enddo
    ispflag=9999
    call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
               & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        iapm=ip_salt(k,is,ne)
        apm_salt(iapm+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo loop_salt_vdif
  endif
 !< shun : end of apm seasalt vdif



 !> shun : apm dust vdif
  if(lfor_dust) then
   !print*,'dust vdif'
   loop_dust_vdif : do is=1,NDSTB
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       iapm=ip_dust(k,is,ne)
       conc(i,j,k)=apm_dust(iapm+ixy)
     enddo
    enddo
    enddo
    ispflag=9999
    call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
               & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        iapm=ip_dust(k,is,ne)
        apm_dust(iapm+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo loop_dust_vdif
  endif
 !< shun : end of apm dust hdif



 !> shun : apm bcoc vdif
  if(lfor_bcoc) then
   !print*,'bcoc vdif'
   loop_bcoc_vdif : do is=1,NBCOCT
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       iapm=ip_bcoc(k,is,ne)
       conc(i,j,k)=apm_bcoc(iapm+ixy)
     enddo
    enddo
    enddo
    ispflag=9999
    call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
               & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        iapm=ip_bcoc(k,is,ne)
        apm_bcoc(iapm+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo loop_bcoc_vdif
  endif
 !< shun : end of apm bcoc hdif

!return

!==============================================================
!==============================================================

! apm coated species
 if(lcoated_dyn) then

  !print*,'coated vdif'

   is=1 ! not used
!-> sulfate on seasalt
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       conc(i,j,k)=msltsulf(i03+ixy)
     enddo
   enddo
   enddo
   ispflag=9999
   call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
              & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       msltsulf(i03+ixy)=conc(i,j,k)
     enddo
   enddo
   enddo

!-> sulfate on dust
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       conc(i,j,k)=mdstsulf(i03+ixy)
     enddo
   enddo
   enddo
   ispflag=9999
   call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
              & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       mdstsulf(i03+ixy)=conc(i,j,k)
     enddo
   enddo
   enddo

!-> sulfate on BC
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       conc(i,j,k)=mbcsulf(i03+ixy)
     enddo
   enddo
   enddo
   ispflag=9999
   call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
              & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       mbcsulf(i03+ixy)=conc(i,j,k)
     enddo
   enddo
   enddo


!-> sulfate on POC
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       conc(i,j,k)=mocsulf(i03+ixy)
     enddo
   enddo
   enddo
   ispflag=9999
   call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
              & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       mocsulf(i03+ixy)=conc(i,j,k)
     enddo
   enddo
   enddo

 endif ! lcoated_dyn

!-> sulferic acid vapor
 if(lfor_h2so4) then
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       conc(i,j,k)=h2so4_gas(i03+ixy)
     enddo
   enddo
   enddo
   ispflag=1
   call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
              & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       h2so4_gas(i03+ixy)=conc(i,j,k)
     enddo
   enddo
   enddo
 endif


!ENDIF ! apm flag

!================================
   loop_binbc : do is=1,nbincb
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       iapm=ip_cbbin(k,is,ne)
       conc(i,j,k)=apm_binbc(iapm+ixy)
     enddo
    enddo
    enddo
    ispflag=9999
    call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
               & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        iapm=ip_cbbin(k,is,ne)
        apm_binbc(iapm+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo loop_binbc



   loop_binoc : do is=1,nbincb
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       iapm=ip_cbbin(k,is,ne)
       conc(i,j,k)=apm_binoc(iapm+ixy)
     enddo
    enddo
    enddo
    ispflag=9999
    call diffus( myid,ispflag,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dt &
               & ,rkv,dzz,ttn,ppp,conc,atm,kktop)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        iapm=ip_cbbin(k,is,ne)
        apm_binoc(iapm+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo loop_binoc



end subroutine apm_v_dif





