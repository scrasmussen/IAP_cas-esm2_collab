
module adv1d_comv

integer :: icpu

integer :: idx_dom

integer :: idx_spc

integer :: it_step

integer :: ixx,iyy

integer :: istx_dim,iedx_dim,isty_dim,iedy_dim

character :: fcheck*50

integer :: nt_step

real :: tmp1d_dat(20)

real,allocatable,dimension(:) :: jcbin

end module adv1d_comv


