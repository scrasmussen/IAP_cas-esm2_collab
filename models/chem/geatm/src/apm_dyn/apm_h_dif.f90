
subroutine apm_h_dif &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dx,dy &
 &  ,kh &
 &  ,ip2mem,mem2d &
 &  ,ktop &
 &  ,ip3mem,mem3d )

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
!integer :: k,is

integer :: mem3d
real,dimension(mem3d) :: dx,dy,kh

integer :: mem2d
real,dimension(mem2d) :: ktop


integer :: ixy,i02,i03,iapm

integer :: ip2mem(nest)
integer :: ip3mem(nzz,nest)

real,dimension(mem3d) :: wk

!IF(lapm) THEN ! apm flag

 !> shun : apm sulfate hdif
if(lfor_sulf) then
   !print*,'sulf hdif'
   loop_so4_hdif : do is=1,NSO4

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_sulf(k,is,ne)
          wk(i03+ixy)=apm_sulf(iapm+ixy)
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_sulf(k,is,ne)
          apm_sulf(iapm+ixy)=wk(i03+ixy)
        enddo
     enddo
     enddo

   enddo loop_so4_hdif
  endif
 !< shun : end of apm sulfate hdif

!return


 !> shun : apm seasalt hdif
  if(lfor_salt) then
   !print*,'salt hdif'
   loop_salt_hdif : do is=1,NSEA

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_salt(k,is,ne)
          wk(i03+ixy)=apm_salt(iapm+ixy)
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_salt(k,is,ne)
          apm_salt(iapm+ixy)=wk(i03+ixy)
        enddo
     enddo
     enddo

   enddo loop_salt_hdif
  endif
 !< shun : end of apm seasalt hdif

 !> shun : apm dust hdif
  if(lfor_dust) then
   !print*,'dust hdif'
   loop_dust_hdif : do is=1,NDSTB

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_dust(k,is,ne)
          wk(i03+ixy)=apm_dust(iapm+ixy)
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_dust(k,is,ne)
          apm_dust(iapm+ixy)=wk(i03+ixy)
        enddo
     enddo
     enddo

   enddo loop_dust_hdif
  endif
 !< shun : end of apm dust hdif

 !> shun : apm bcoc hdif
  if(lfor_bcoc) then
   !print*,'bcoc hdif'
   loop_bcoc_hdif : do is=1,NBCOCT

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)
          wk(i03+ixy)=apm_bcoc(iapm+ixy)
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       iapm=ip_bcoc(k,is,ne)
       call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)
          apm_bcoc(iapm+ixy)=wk(i03+ixy)
        enddo
     enddo
     enddo

   enddo loop_bcoc_hdif
  endif
 !< shun : end of apm bcoc hdif

!return

!==========================================================================
!==========================================================================

! apm coated species

 if(lcoated_dyn) then
    !print*,'coated hdif'
!-> sulfate on seasalt
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=msltsulf(i03+ixy)
       enddo
    enddo
    enddo
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         msltsulf(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo

!-> sulfate on dust
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=mdstsulf(i03+ixy)
       enddo
    enddo
    enddo
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         mdstsulf(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo

!-> sulfate on BC
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=mbcsulf(i03+ixy)
       enddo
    enddo
    enddo
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         mbcsulf(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo

!-> sulfate on POC
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=mocsulf(i03+ixy)
       enddo
    enddo
    enddo
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         mocsulf(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo

  endif ! lcoated_dyn

!-> suferic acid vapor
  if(lfor_h2so4) then
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=h2so4_gas(i03+ixy)
       enddo
    enddo
    enddo
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, wk(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         h2so4_gas(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo
  endif
!!!!!!!!!



!ENDIF ! apm flag


end subroutine apm_h_dif





