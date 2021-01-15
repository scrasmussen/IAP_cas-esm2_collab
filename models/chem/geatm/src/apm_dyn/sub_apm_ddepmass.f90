subroutine get_apmddepmass &
 & ( myid &
 &  ,lapm &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dx,dy,dz
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
integer :: mem2d,mem3d
integer :: ip3mem(nzz,nest),ip2mem(nest)
real,dimension(mem3d) :: dx,dy,dz

integer :: i,j,k,is,imode

integer :: ixy,i02,i03,iapm

real :: tmp


k=1
i03 = ip3mem(k,ne)
i02 = ip2mem(ne)
loop_j : do j = sy(ne),ey(ne) 
loop_i : do i = sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

  tmp1=0.0
  do is=1,NSO4
     iapm=ip_sulf(k,is,ne)
     tmp1=tmp1+apm_sulf(iapm+ixy)
  enddo

  tmp2 = msltsulf(i03+ixy)+mdstsulf(i03+ixy)+mbcsulf(i03+ixy)+mocsulf(i03+ixy)

  deltc1 = apm_spsulfate(i03+ixy) - tmp1
  deltc2 = apm_ppsulfate(i03+ixy) - tmp2
  deltct = deltc1+deltc2

  ddep_apm_spsulf(i02+ixy)=deltc1*dz(i03+ixy)+ddep_apm_spsulf(i02+ixy)
  ddep_apm_ppsulf(i02+ixy)=deltc2*dz(i03+ixy)+ddep_apm_ppsulf(i02+ixy)
  ddep_apm_sulf(i02+ixy)  =deltct*dz(i03+ixy)+ddep_apm_sulf(i02+ixy)

enddo loop_i
enddo loop_j


end subroutine get_apmddepmass



