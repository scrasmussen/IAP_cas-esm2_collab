subroutine get_funitnaqpms(funit)
integer :: funit
logical :: is_used
do funit=10,2000
 inquire(unit=funit, opened=is_used)
 if (.not.is_used) exit
enddo
return
end
