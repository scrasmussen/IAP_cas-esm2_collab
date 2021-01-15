
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

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       iapm=ip_sulf(k,is,ne)
       call dif_hori( myid, apm_sulf(iapm), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo


   enddo loop_so4_hdif
  endif
 !< shun : end of apm sulfate hdif

!return


 !> shun : apm seasalt hdif
  if(lfor_salt) then
   !print*,'salt hdif'
   loop_salt_hdif : do is=1,NSEA

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       iapm=ip_salt(k,is,ne)
       call dif_hori( myid, apm_salt(iapm), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo

   enddo loop_salt_hdif
  endif
 !< shun : end of apm seasalt hdif

 !> shun : apm dust hdif
  if(lfor_dust) then
   !print*,'dust hdif'
   loop_dust_hdif : do is=1,NDSTB

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       iapm=ip_dust(k,is,ne)
       call dif_hori( myid, apm_dust(iapm), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo

   enddo loop_dust_hdif
  endif
 !< shun : end of apm dust hdif

 !> shun : apm bcoc hdif
  if(lfor_bcoc) then
   !print*,'bcoc hdif'
   loop_bcoc_hdif : do is=1,NBCOCT

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       iapm=ip_bcoc(k,is,ne)
       call dif_hori( myid, apm_bcoc(iapm), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
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
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, msltsulf(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo

!-> sulfate on dust
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, mdstsulf(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo

!-> sulfate on BC
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, mbcsulf(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo

!-> sulfate on POC
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, mocsulf(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo

  endif ! lcoated_dyn

!-> suferic acid vapor
  if(lfor_h2so4) then
    do k=1,nzz-1
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      call dif_hori( myid, h2so4_gas(i03), kh(i03), dx(i03), dy(i03), &
         & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
    enddo
  endif
!!!!!!!!!

   loop_binbc : do is=1,nbincb

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       iapm=ip_cbbin(k,is,ne)
       call dif_hori( myid, apm_binbc(iapm), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo

   enddo loop_binbc


   loop_binoc : do is=1,nbincb

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       i02=ip2mem(ne)
       iapm=ip_cbbin(k,is,ne)
       call dif_hori( myid, apm_binoc(iapm), kh(i03), dx(i03), dy(i03), &
          & sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,is)
     enddo

   enddo loop_binoc



!ENDIF ! apm flag


end subroutine apm_h_dif





