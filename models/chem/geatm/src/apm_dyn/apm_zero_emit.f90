
subroutine apm_zero_emit ( myid,lapm,ne,nx,ny,nzz,nest,ip3mem,sx,ex,sy,ey )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid
logical :: lapm

integer :: i,j,k,is
integer :: ixy,i03,iapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: ip3mem(nzz,nest)


do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)
  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
  do k=1,nzz

    i03=ip3mem(k,ne)

    sulf_emit(i03+ixy) = 0.0d0 


    bc_emit(i03+ixy)   = 0.0d0

    oc_emit(i03+ixy)   = 0.0d0

    bc_emit_bb(i03+ixy)   = 0.0d0
 
    bc_emit_ff(i03+ixy)   = 0.0d0

    oc_emit_bb(i03+ixy)   = 0.0d0

    oc_emit_ff(i03+ixy)   = 0.0d0

    ppmfine_emit(i03+ixy) = 0.0d0

    do is=1,NSEA
      iapm=ip_salt(k,is,ne)
      salt_emit(iapm+ixy) = 0.0d0
    enddo

    dust_emit(i03+ixy)=0.0d0

  enddo
enddo
enddo


end subroutine apm_zero_emit





