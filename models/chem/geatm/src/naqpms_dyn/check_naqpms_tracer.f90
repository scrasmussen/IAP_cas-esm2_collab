
subroutine check_naqpms_tracer &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,flag )

use naqpms_varlist
implicit none

integer :: myid

real :: dt

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

character(len=*) :: flag

integer :: i,j,k,ig,i04,i03,ixy

real,parameter :: cmax=1.0e6


do ig=1,iedgas
  do k=1,nzz
    i03=ip3mem(k,ne)
    i04=ip4mem(k,ig,ne)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       if(.not.(gas(i04+ixy).ge.0.and.gas(i04+ixy).le.cmax)) then
         print*,trim(flag),i,j,k,ig,gas(i04+ixy)
         stop
       endif
    enddo
    enddo
  enddo

enddo


end subroutine check_naqpms_tracer
