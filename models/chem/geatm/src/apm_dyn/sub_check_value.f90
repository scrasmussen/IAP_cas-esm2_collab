
subroutine check_value(i,j,k,lower,upper,value)
implicit none
integer :: i,j,k
real    :: lower,upper
real    :: value
if(value.gt.lower.and.value.le.upper) then
 return
else
 print*,'check diameter in deposition'
 print*,i,j,k
 print*,lower,value,upper
 print*,'in source code file sub_check_value.f90'
 stop
endif
end subroutine check_value


