
subroutine apm_clean_sp &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip3mem,mem3d &
 &  ,prgname )

use apm_varlist

use APM_INIT_MOD, only : IFNUCL

implicit none
include 'apm_parm.inc'

character*25 :: prgname

integer :: myid

real :: dt

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is
!integer :: k,is

integer :: mem3d

integer :: ixy,i02,i03,iapm

integer :: ip3mem(nzz,nest)

!real,parameter :: upper=1.0e4,lower=0

real :: tmp


!print*,sy(ne)-1,ey(ne)+1,sx(ne)-1,ex(ne)+1
!print*,ne

!return

!stop

if(IFNUCL.ne.0) then

return

else

   loop_so4 : do is=1,NSO4

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
          iapm=ip_sulf(k,is,ne)
          apm_sulf(iapm+ixy)=1.0d-30
        enddo
     enddo
     enddo

   enddo loop_so4

endif

end subroutine apm_clean_sp





