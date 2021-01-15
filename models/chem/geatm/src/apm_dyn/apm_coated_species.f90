

! distribute mass of Sulfate, Nitrate, Ammonium, MSA, SOA
! according to surface area to each bin in each particle category
! 195426893

subroutine get_coated_tracers

 use apm_varlist
 include 'apm_parm.inc'
 


end subroutine get_coated_tracers


!
subroutine apm_sulf_1surface

end subroutine apm_sulf_1surface
!

!
subroutine apm_sulf_surface

end subroutine apm_sulf_surface
!



!
subroutine apm_coated_s(nx,ny,nbin,r,cts)
implicit none
integer :: nx,ny,nbin
real*8 :: r
real*8,dimension(nx,ny) :: cts

end subroutine apm_coated_s
