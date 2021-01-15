
subroutine wr_gas( myid,gas,mem4d,ip4mem,nz,ngas,nest,ne &
                  ,sx,ex,sy,ey )

integer :: mem4d
real    :: gas(mem4d)

integer :: nz,ngas,nest

integer :: ip4mem(nz,ngas,nest)

integer :: sx(nest),ex(nest),sy(nest),ey(nest)

integer :: myid

integer :: i,j,k,ixy,ig,ip


real              :: msulf,msulf_pp,msulf_sp
integer,parameter :: nsulf = 8
integer,parameter :: sulf_o_gas(1:nsulf) = (/ 1, 1, 1, 1, 1, 1, 1, 1/)
integer,parameter :: sulf_index(1:nsulf) = (/83,84,87,89,92,93,94,95/)


integer,parameter :: nnx=69,nny=74


real :: test_so4(nnx,nny)


ig=1

do k=1,1 !nzz
ip=ip4mem(k,ig,ne)
do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)
 ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
 !print*,gas(ip+ixy)

enddo
enddo
enddo

!ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1




end
