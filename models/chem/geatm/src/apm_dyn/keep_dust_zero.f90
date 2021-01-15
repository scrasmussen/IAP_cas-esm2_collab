
subroutine keep_dust_zero &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex )

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


integer :: ixy,i02,i03,iapm



 !> shun : apm dust hdif
   loop_dust_hdif : do is=1,NDSTB

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          iapm=ip_dust(k,is,ne)
          apm_dust(iapm+ixy)=0
        enddo
     enddo
     enddo

   enddo loop_dust_hdif
 !< shun : end of apm dust hdif

end subroutine keep_dust_zero





