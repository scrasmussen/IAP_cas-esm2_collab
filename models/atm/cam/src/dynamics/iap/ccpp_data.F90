module ccpp_data
  !! \section arg_table_ccpp_data Argument Table
  !! \htmlinclude ccpp_data.html
  !!
  use shr_kind_mod,      only: r8 => SHR_KIND_R8
  use physics_types,     only: physics_state, physics_int_ephem, physics_int_pers, physics_global
  use phys_grid,         only: ngcols
  use camsrfexch_types , only: cam_in_t
    

  implicit none

  private
  
  public nchnks, &
         ngcols, &
         ccpp_suite, &
         dt, &
         phys_int_ephem, &
         phys_int_pers, &
         phys_global, &
         phys_state, &
         cam_in
  
  integer,            save               :: nchnks
  character(len=256), save               :: ccpp_suite='undefined'
  real(kind=r8),      save               :: dt
  
  type(physics_int_ephem), save, pointer :: phys_int_ephem(:) => null()
  type(physics_int_pers),  save, pointer :: phys_int_pers(:) => null()
  type(physics_global),    save          :: phys_global
  type(physics_state),     save, pointer :: phys_state(:) => null()

  type(cam_in_t),          save, pointer :: cam_in(:) => null()
end module ccpp_data
