

subroutine apm_isor_tracer & 
                           & ( myid &
                           &  ,lapm &
                           &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
                           &  ,ip3mem,mem3d )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid
logical :: lapm
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: mem3d

integer :: ip3mem(nzz,nest)

integer :: i03,ixy,iapm

integer :: i,j,k,is


return


 loop_k : DO k=1,nzz-1

  i03 = ip3mem(k,ne)

  loop_j : DO j = sy(ne),ey(ne)
  loop_i : DO i = sx(ne),ex(ne)

    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

    do is=1,NSO4
      iapm=ip_sulf(k,is,ne)
      apm_sulf(iapm+ixy)=apm_sulf(iapm+ixy)*ch_ratio(i03+ixy) 
    enddo
 
    msltsulf(i03+ixy) = msltsulf(i03+ixy)*ch_ratio(i03+ixy)
    mdstsulf(i03+ixy) = mdstsulf(i03+ixy)*ch_ratio(i03+ixy)
    mbcsulf(i03+ixy)  = mbcsulf(i03+ixy)*ch_ratio(i03+ixy) 
    mocsulf(i03+ixy)  = mocsulf(i03+ixy)*ch_ratio(i03+ixy)

  enddo loop_i
  enddo loop_j

 enddo loop_k


end subroutine apm_isor_tracer 






