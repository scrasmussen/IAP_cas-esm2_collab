
subroutine check_matrix(label,lower,upper, ist,ied,jst,jed,matrix)
implicit none
integer,parameter :: npoint=10000
character :: label*200
integer :: i,j,k
integer :: ist,jst,ied,jed
real    :: lower,upper
real    :: value
real    :: matrix(ist:ied,jst:jed)

logical :: lnormal

integer :: ico(npoint),jco(npoint)
real    :: vij(npoint)
integer :: itpoint,ipoint

lnormal=.true.
itpoint=0

do j=jst,jed
do i=ist,ied

 value=matrix(i,j)

 if(value.ge.lower.and.value.le.upper) then
  cycle
 else
  lnormal=.false.
  itpoint=itpoint+1
  ico(itpoint)=i
  jco(itpoint)=j
  vij=value
 endif

enddo
enddo

if(.not.lnormal) then
print*,'shun check: ',trim(label)
print*,'itpoint=',itpoint
do ipoint=1,itpoint
 print*,ico(ipoint),jco(ipoint),vij(ipoint)
enddo
print*,'in source code file sub_check_matrix.f90'
stop
endif


end subroutine check_matrix


