module isorropia_var_mod

integer,parameter :: isorropia_nvar = 22
integer,parameter :: isorropia_idx(1:isorropia_nvar) = (/   1,  2,  3,  4, 80 &
                                                       , 81, 82, 83, 84, 85 &
                                                       , 86, 87, 88, 89, 90 &
                                                       , 91, 92, 93, 94, 95 &
                                                       , 79,102  /)

real,dimension(isorropia_nvar) :: isorropia_var,iso_no3_ratio,iso_nh4_ratio



integer :: idx_iso

contains

subroutine  naqpms_aqueous_ratio
implicit none
integer :: i


end subroutine naqpms_aqueous_ratio



subroutine naqpms_aqueous_update( naer, caer, i03, ixy, mem4d, gas )
implicit none
integer :: naer, i03, ixy
real    :: caer(naer) ! naer=9
real    :: gas(mem4d)
integer :: i



end subroutine naqpms_aqueous_update



end module isoropia_var_mod
