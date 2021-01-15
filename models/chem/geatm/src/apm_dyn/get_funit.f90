subroutine get_funit(funit)
integer :: funit
logical :: is_used
do funit=10,1000
 inquire(unit=funit, opened=is_used)
 if (.not.is_used) exit
enddo
return
end
