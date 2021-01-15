
subroutine print_time(myid,lprint_time,timestr)
implicit none
!logical,save :: lprint_time

character*19 :: timestr
integer :: myid
logical :: lprint_time

if(lprint_time) then
 print*,'current time : ',timestr,myid
 !lprint_time=.false.
endif

return

end subroutine print_time
