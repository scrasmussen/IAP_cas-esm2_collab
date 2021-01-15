subroutine get_apmsulf &
 & ( myid &
 &  ,lapm &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid


logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is,imode

integer :: mem2d,mem3d

integer :: ixy,i02,i03,iapm

integer :: ip3mem(nzz,nest),ip2mem(nest)


real :: tmp


loop_j : do j = sy(ne),ey(ne) 
loop_i : do i = sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

  loop_k : do k=1,nzz

    i03 = ip3mem(k,ne)

    tmp=0.0
    do is=1,NSO4
      iapm=ip_sulf(k,is,ne)
      tmp=tmp+apm_sulf(iapm+ixy)
    enddo

    apm_spsulfate(i03+ixy)=tmp

    apm_ppsulfate(i03+ixy)=msltsulf(i03+ixy)+mdstsulf(i03+ixy) &
                          +mbcsulf(i03+ixy)+mocsulf(i03+ixy)

    apm_sulfate(i03+ixy)=apm_spsulfate(i03+ixy)+apm_ppsulfate(i03+ixy)

  enddo loop_k

enddo loop_i
enddo loop_j


end subroutine get_apmsulf



