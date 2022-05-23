module physics_register
  
  use phys_buffer,  only: pbuf_add
  use phys_control, only: phys_getopts
  use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
  
  implicit none
  
  private
  
  public :: convect_deep_register
  
  contains
  
  !=========================================================================================
  subroutine convect_deep_register()

  !----------------------------------------
  ! Purpose: register fields with the physics buffer
  !----------------------------------------

    implicit none
    
    integer :: idx
    character(len=16) :: deep_scheme    ! default set in phys_control.F90, use namelist to change
    
    ! get deep_scheme setting from phys_control
    call phys_getopts(deep_scheme_out = deep_scheme)

  ! zmh
    select case ( deep_scheme )
    case('ZM') !    Zhang-McFarlane (default)
       !call zm_conv_register
       call pbuf_add('DP_FLXPRC', 'global', 1, pverp, 1, idx) ! Flux of precipitation from deep convection (kg/m2/s)
       call pbuf_add('DP_FLXSNW', 'global', 1, pverp, 1, idx) ! Flux of snow from deep convection (kg/m2/s) 
       call pbuf_add('DP_CLDLIQ', 'global', 1, pver,  1, idx) ! deep gbm cloud liquid water (kg/kg)
       call pbuf_add('DP_CLDICE', 'global', 1, pver,  1, idx) ! deep gbm cloud liquid water (kg/kg)    
    case('ZYX1')
          ! if( trim(plume_model) == 'cam' .and. trim(shallow_scheme) =='off')then
          !     write(*,*)' When plume model is cam, shallow convection scheme needs to be set!'
          !     stop
          ! endif
          ! call zyx1_conv_intr_register

    case('ZYX2')
          ! if( trim(plume_model) == 'cam' .and. trim(shallow_scheme) =='off')then
          !     write(*,*)' When plume model is cam, shallow convection scheme needs to be set!'
          !     stop
          ! endif
          ! call zyx2_conv_intr_register

    end select

    call pbuf_add('ICWMRDP' , 'physpkg', 1,pver,      1,         idx)
    call pbuf_add('RPRDDP' , 'physpkg', 1,pver,       1,         idx)
    call pbuf_add('NEVAPR_DPCU' , 'physpkg', 1,pver,      1,     idx)

  !wxc zmh
    call pbuf_add('slflxdp' , 'physpkg', 1,pverp,      1,        idx)
    call pbuf_add('qtflxdp' , 'physpkg', 1,pverp,      1,        idx)


  end subroutine convect_deep_register

  !=========================================================================================
  
  
end module physics_register