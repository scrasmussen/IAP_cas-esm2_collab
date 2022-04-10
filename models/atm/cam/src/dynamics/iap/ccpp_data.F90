module ccpp_data
  !! \section arg_table_ccpp_data Argument Table
  !! \htmlinclude ccpp_data.html
  !!
  use shr_kind_mod,      only: r8 => SHR_KIND_R8
  use physics_types,     only: physics_state, physics_int_ephem, physics_int_pers, physics_global
  use ccpp_types,        only: ccpp_t
  use phys_grid,         only: ngcols
  use camsrfexch_types , only: cam_in_t  
    

  implicit none

  private
  
  public nchnks, &
         cdata_domain, &
         cdata_chunk, &
         ccpp_suite, &
         dt, &
         phys_int_ephem, &
         phys_int_pers, &
         phys_global, &
         phys_state, &
         cam_in
  
  integer                                :: nchnks
  type(ccpp_t),            save          :: cdata_domain
  type(ccpp_t),            save, pointer :: cdata_chunk(:)
  character(len=256)                     :: ccpp_suite='undefined'
  real(kind=r8)                          :: dt
  
  type(physics_int_ephem), save, pointer :: phys_int_ephem(:)
  type(physics_int_pers),  save, pointer :: phys_int_pers(:)
  type(physics_global),    save          :: phys_global
  type(physics_state),     save, pointer :: phys_state(:)

  type(cam_in_t),          save, pointer :: cam_in(:)
end module ccpp_data
