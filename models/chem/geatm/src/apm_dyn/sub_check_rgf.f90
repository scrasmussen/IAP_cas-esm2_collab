
subroutine check_rgf(label,i,j,k,lower,upper,value)
implicit none
character :: label*100
integer :: i,j,k
real    :: lower,upper
real    :: value
if(value.ge.lower.and.value.le.upper) then
 return
else
 print*,'check value :',trim(label)
 print*,i,j,k
 print*,lower,value,upper
 print*,'in source code file sub_check_rgf.f90'
 stop
endif
end subroutine check_rgf


