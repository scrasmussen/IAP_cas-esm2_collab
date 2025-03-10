module camsrfexch_types
!-----------------------------------------------------------------------
!
! Module to handle data that is exchanged between the CAM atmosphere
! model and the surface models (land, sea-ice, and ocean).
!
!-----------------------------------------------------------------------
!
! USES:
!
!> \section arg_table_camsrfexch_types
!! \htmlinclude camsrfexch_types.html
!!
  use shr_kind_mod,  only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents,  only: pcnst
  use ppgrid,        only: pcols, begchunk, endchunk
  use phys_grid,     only: get_ncols_p, phys_grid_initialized
  use abortutils,    only: endrun
  use infnan,        only: inf
  use cam_logfile,   only: iulog

  implicit none

!----------------------------------------------------------------------- 
! PRIVATE: Make default data and interfaces private
!----------------------------------------------------------------------- 
  private     ! By default all data is private to this module
!
! Public interfaces
!
  public atm2hub_alloc              ! Atmosphere to surface data allocation method
  public hub2atm_alloc              ! Merged hub surface to atmosphere data allocation method
  public hub2atm_setopts            ! Set options to allocate optional parts of data type
  public atm2hub_deallocate
  public hub2atm_deallocate
!
! Public data types
!
  public cam_out_t                  ! Data from atmosphere
  public cam_in_t                   ! Merged surface data

!---------------------------------------------------------------------------
! This is the data that is sent from the atmosphere to the surface models
!---------------------------------------------------------------------------
!> \section arg_table_cam_out_t
!! \htmlinclude cam_out_t.html
!!
  type cam_out_t 
     integer  :: lchnk               ! chunk index
     integer  :: ncol                ! number of columns in chunk
     real(r8) :: tbot(pcols)         ! bot level temperature
     real(r8) :: zbot(pcols)         ! bot level height above surface
     real(r8) :: ubot(pcols)         ! bot level u wind
     real(r8) :: vbot(pcols)         ! bot level v wind
     real(r8) :: qbot(pcols,pcnst)   ! bot level specific humidity
     real(r8) :: pbot(pcols)         ! bot level pressure
     real(r8) :: rho(pcols)          ! bot level density	
     real(r8) :: netsw(pcols)        !	
     real(r8) :: flwds(pcols)        ! 
     real(r8) :: precsc(pcols)       !
     real(r8) :: precsl(pcols)       !
     real(r8) :: precc(pcols)        ! 
     real(r8) :: precl(pcols)        ! 
     real(r8) :: soll(pcols)         ! 
     real(r8) :: sols(pcols)         ! 
     real(r8) :: solld(pcols)        !
     real(r8) :: solsd(pcols)        !
     real(r8) :: srfrad(pcols)       !
     real(r8) :: thbot(pcols)        ! 
     real(r8) :: co2prog(pcols)      ! prognostic co2
     real(r8) :: co2diag(pcols)      ! diagnostic co2
     real(r8) :: psl(pcols)
!wangty modify
#ifdef wrf 
     real(r8) :: clflo(pcols)  ! juanxiong he
     real(r8) :: clfmi(pcols)  ! juanxiong he
     real(r8) :: clfhi(pcols)  ! juanxiong he
     real(r8) :: pblh(pcols)  ! juanxiong he
#endif
     real(r8) :: bcphiwet(pcols)     ! wet deposition of hydrophilic black carbon
     real(r8) :: bcphidry(pcols)     ! dry deposition of hydrophilic black carbon
     real(r8) :: bcphodry(pcols)     ! dry deposition of hydrophobic black carbon
     real(r8) :: ocphiwet(pcols)     ! wet deposition of hydrophilic organic carbon
     real(r8) :: ocphidry(pcols)     ! dry deposition of hydrophilic organic carbon
     real(r8) :: ocphodry(pcols)     ! dry deposition of hydrophobic organic carbon
     real(r8) :: dstwet1(pcols)      ! wet deposition of dust (bin1)
     real(r8) :: dstdry1(pcols)      ! dry deposition of dust (bin1)
     real(r8) :: dstwet2(pcols)      ! wet deposition of dust (bin2)
     real(r8) :: dstdry2(pcols)      ! dry deposition of dust (bin2)
     real(r8) :: dstwet3(pcols)      ! wet deposition of dust (bin3)
     real(r8) :: dstdry3(pcols)      ! dry deposition of dust (bin3)
     real(r8) :: dstwet4(pcols)      ! wet deposition of dust (bin4)
     real(r8) :: dstdry4(pcols)      ! dry deposition of dust (bin4)
  end type cam_out_t 

!---------------------------------------------------------------------------
! This is the merged state of sea-ice, land and ocean surface parameterizations
!---------------------------------------------------------------------------
!> \section arg_table_cam_in_t
!! \htmlinclude cam_in_t.html
!!
  type cam_in_t
     integer  :: lchnk                   ! chunk index
     integer  :: ncol                    ! number of active columns
     real(r8) :: asdir(pcols)            ! albedo: shortwave, direct
     real(r8) :: asdif(pcols)            ! albedo: shortwave, diffuse
     real(r8) :: aldir(pcols)            ! albedo: longwave, direct
     real(r8) :: aldif(pcols)            ! albedo: longwave, diffuse
     real(r8) :: lwup(pcols)             ! longwave up radiative flux
     real(r8) :: lhf(pcols)              ! latent heat flux
     real(r8) :: shf(pcols)              ! sensible heat flux
     real(r8) :: wsx(pcols)              ! surface u-stress (N)
     real(r8) :: wsy(pcols)              ! surface v-stress (N)
     real(r8) :: tref(pcols)             ! ref height surface air temp
     real(r8) :: qref(pcols)             ! ref height specific humidity 
!wangty modify
     real(r8) :: rhref(pcols)            ! ref height relative humidity
     real(r8) :: u10(pcols)              ! 10m wind speed
     real(r8) :: ts(pcols)               ! merged surface temp 
     real(r8) :: sst(pcols)              ! sea surface temp
     real(r8) :: snowhland(pcols)        ! snow depth (liquid water equivalent) over land 
     real(r8) :: snowhice(pcols)         ! snow depth over ice
     real(r8) :: fco2_lnd(pcols)         ! co2 flux from lnd
     real(r8) :: fco2_ocn(pcols)         ! co2 flux from ocn
     real(r8) :: fdms(pcols)             ! dms flux
     real(r8) :: landfrac(pcols)         ! land area fraction
     real(r8) :: icefrac(pcols)          ! sea-ice areal fraction
     real(r8) :: ocnfrac(pcols)          ! ocean areal fraction
     real(r8), pointer, dimension(:) :: ram1  !aerodynamical resistance (s/m) (pcols)
     real(r8), pointer, dimension(:) :: fv    !friction velocity (m/s) (pcols)
     real(r8) :: cflx(pcols,pcnst)      ! constituent flux (evap)
     real(r8) :: ustar(pcols)            ! atm/ocn saved version of ustar
     real(r8) :: re(pcols)               ! atm/ocn saved version of re
     real(r8) :: ssq(pcols)              ! atm/ocn saved version of ssq
     real(r8), pointer, dimension(:,:) :: depvel ! deposition velocities
!wangty modify
#ifdef wrf 
!--------------------------------------------------------------------------------------
! soil depth/height and soil temperature/moisture for wrf/cam coupling, added
! juanxiong he
!--------------------------------------------------------------------------------------

     real(r8) :: lhf_ndg_old(pcols)              ! latent heat flux
     real(r8) :: shf_ndg_old(pcols)              ! sensible heat flux
     real(r8) :: wsx_ndg_old(pcols)              ! surface u-stress (N)
     real(r8) :: wsy_ndg_old(pcols)              ! surface v-stress (N)
     real(r8) :: tref_ndg_old(pcols)             ! ref height surface air temp
     real(r8) :: qref_ndg_old(pcols)             ! ref height specific humidity
     real(r8) :: u10_ndg_old(pcols)              ! 10m wind speed
     real(r8) :: ts_ndg_old(pcols)               ! merged surface temp
     real(r8) :: ulwrf_ndg_old(pcols)            ! upward long wave radiation

     real(r8) :: lhf_ndg_new(pcols)              ! latent heat flux
     real(r8) :: shf_ndg_new(pcols)              ! sensible heat flux
     real(r8) :: wsx_ndg_new(pcols)              ! surface u-stress (N)
     real(r8) :: wsy_ndg_new(pcols)              ! surface v-stress (N)
     real(r8) :: tref_ndg_new(pcols)             ! ref height surface air temp
     real(r8) :: qref_ndg_new(pcols)             ! ref height specific humidity
     real(r8) :: u10_ndg_new(pcols)              ! 10m wind speed
     real(r8) :: ts_ndg_new(pcols)               ! merged surface temp
     real(r8) :: ulwrf_ndg_new(pcols)            ! upward long wave radiation

     real(r8) :: soildepth(pcols,4)
     real(r8) :: soilthick(pcols,4)
     real(r8) :: soilt(pcols,4)
     real(r8) :: soilm(pcols,4)

!--------------------------------------------------------------------------------------
! soil depth/height and soil temperature/moisture for wrf/cam coupling, added
! juanxiong he
!--------------------------------------------------------------------------------------
#endif
  end type cam_in_t    

  logical :: dust = .false.     ! .true. => aerosol dust package is being used

!===============================================================================
CONTAINS
!===============================================================================

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: hub2atm_alloc
!
! !DESCRIPTION:
!
!   Allocate space for the surface to atmosphere data type. And initialize
!   the values.
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine hub2atm_alloc( cam_in )
    use seq_drydep_mod, only : lnd_drydep, n_drydep
!
!!ARGUMENTS:
!
   type(cam_in_t), pointer ::  cam_in(:)     ! Merged surface state
!
!!LOCAL VARIABLES:
!
    integer :: c        ! chunk index
    integer :: ierror   ! Error code
!----------------------------------------------------------------------- 
! 
! EOP
!
    if ( .not. phys_grid_initialized() ) call endrun( "HUB2ATM_ALLOC error: phys_grid not called yet" )
    allocate (cam_in(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('HUB2ATM_ALLOC error: allocation error')
    end if

    do c = begchunk,endchunk
!wangty modify
#ifdef wrf 
!       nullify(cam_in(c)%ram1)
!       nullify(cam_in(c)%fv)
#else
       nullify(cam_in(c)%ram1)
       nullify(cam_in(c)%fv)
#endif
       nullify(cam_in(c)%depvel)
    enddo  
!wangty modify
#ifdef wrf 
!    if ( dust ) then  ! juanxiong he
#else
    if ( dust ) then
#endif
       do c = begchunk,endchunk 
          allocate (cam_in(c)%ram1(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ram1')
          allocate (cam_in(c)%fv(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error fv')
       end do
!wangty modify
#ifdef wrf 
!    endif  !dust   !juanxiong he
#else
    endif  !dust  
#endif

    if (lnd_drydep .and. n_drydep>0) then
       do c = begchunk,endchunk 
          allocate (cam_in(c)%depvel(pcols,n_drydep), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error depvel')
       end do
    endif
          
    do c = begchunk,endchunk
       cam_in(c)%lchnk = c
       cam_in(c)%ncol  = get_ncols_p(c)
       cam_in(c)%asdir    (:) = 0._r8
       cam_in(c)%asdif    (:) = 0._r8
       cam_in(c)%aldir    (:) = 0._r8
       cam_in(c)%aldif    (:) = 0._r8
       cam_in(c)%lwup     (:) = 0._r8
       cam_in(c)%lhf      (:) = 0._r8
       cam_in(c)%shf      (:) = 0._r8
       cam_in(c)%wsx      (:) = 0._r8
       cam_in(c)%wsy      (:) = 0._r8
       cam_in(c)%tref     (:) = 0._r8
       cam_in(c)%qref     (:) = 0._r8
!wangty modify
       cam_in(c)%rhref     (:) = 0._r8
       cam_in(c)%u10      (:) = 0._r8
       cam_in(c)%ts       (:) = 0._r8
       cam_in(c)%sst      (:) = 0._r8
       cam_in(c)%snowhland(:) = 0._r8
       cam_in(c)%snowhice (:) = 0._r8
       cam_in(c)%fco2_lnd (:) = 0._r8
       cam_in(c)%fco2_ocn (:) = 0._r8
       cam_in(c)%fdms     (:) = 0._r8
       cam_in(c)%landfrac (:) = inf
       cam_in(c)%icefrac  (:) = inf
       cam_in(c)%ocnfrac  (:) = inf
!wangty modify
#ifdef wrf 
!       if ( dust ) then  ! juanxiong he
          cam_in(c)%ram1  (:) = 0.1_r8
          cam_in(c)%fv    (:) = 0.1_r8
!       endif    ! juanxiong he 
#else
       if ( dust ) then  ! juanxiong he
          cam_in(c)%ram1  (:) = 0.1_r8
          cam_in(c)%fv    (:) = 0.1_r8
       endif    ! juanxiong he
#endif
       cam_in(c)%cflx   (:,:) = 0._r8
       cam_in(c)%ustar    (:) = 0._r8
       cam_in(c)%re       (:) = 0._r8
       cam_in(c)%ssq      (:) = 0._r8
       if (lnd_drydep .and. n_drydep>0) then
          cam_in(c)%depvel (:,:) = 0._r8
       endif
!wangty modify
#ifdef wrf 
!--------------------------------------------------------------------------------------
! soil depth/height and soil temperature/moisture for wrf/cam coupling, added
! juanxiong he
!--------------------------------------------------------------------------------------
       cam_in(c)%soildepth(:,:) = 0._r8
       cam_in(c)%soilthick(:,:) = 0._r8
       cam_in(c)%soilt(:,:) = 0._r8
       cam_in(c)%soilm(:,:) = 0._r8

       cam_in(c)%lhf_ndg_old(:) = 0._r8              ! latent heat flux
       cam_in(c)%shf_ndg_old(:)= 0._r8               ! sensible heat flux
       cam_in(c)%wsx_ndg_old(:)= 0._r8               ! surface u-stress (N)
       cam_in(c)%wsy_ndg_old(:)= 0._r8               ! surface v-stress (N)
       cam_in(c)%tref_ndg_old(:)= 0._r8              ! ref height surface air temp
       cam_in(c)%qref_ndg_old(:)= 0._r8              ! ref height specific humidity
       cam_in(c)%u10_ndg_old(:)= 0._r8               ! 10m wind speed
       cam_in(c)%ts_ndg_old(:)= 0._r8                ! merged surface temp
       cam_in(c)%ulwrf_ndg_old(:)= 0._r8             ! upward long wave radiation

       cam_in(c)%lhf_ndg_new(:)= 0._r8               ! latent heat flux
       cam_in(c)%shf_ndg_new(:)= 0._r8               ! sensible heat flux
       cam_in(c)%wsx_ndg_new(:)= 0._r8               ! surface u-stress (N)
       cam_in(c)%wsy_ndg_new(:)= 0._r8               ! surface v-stress (N)
       cam_in(c)%tref_ndg_new(:)= 0._r8              ! ref height surface air temp
       cam_in(c)%qref_ndg_new(:)= 0._r8              ! ref height specific humidity
       cam_in(c)%u10_ndg_new(:) = 0._r8              ! 10m wind speed
       cam_in(c)%ts_ndg_new(:) = 0._r8                ! merged surface temp
       cam_in(c)%ulwrf_ndg_new(:)= 0._r8             ! upward long wave radiation

!--------------------------------------------------------------------------------------
! soil depth/height and soil temperature/moisture for wrf/cam coupling, added
! juanxiong he
!--------------------------------------------------------------------------------------
#endif
    end do

  end subroutine hub2atm_alloc

!
!===============================================================================
!

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: atm2hub_alloc
!
! !DESCRIPTION:
!
!   Allocate space for the atmosphere to surface data type. And initialize
!   the values.
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine atm2hub_alloc( cam_out )
!
!!USES:
!
!
!!ARGUMENTS:
!
   type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
!
!!LOCAL VARIABLES:
!
    integer :: c            ! chunk index
    integer :: ierror       ! Error code
    !----------------------------------------------------------------------- 

    if ( .not. phys_grid_initialized() ) call endrun( "ATM2HUB_ALLOC error: phys_grid not called yet" )
    allocate (cam_out(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('ATM2HUB_ALLOC error: allocation error')
    end if

    do c = begchunk,endchunk
       cam_out(c)%lchnk       = c
       cam_out(c)%ncol        = get_ncols_p(c)
       cam_out(c)%tbot(:)     = 0._r8
       cam_out(c)%zbot(:)     = 0._r8
       cam_out(c)%ubot(:)     = 0._r8
       cam_out(c)%vbot(:)     = 0._r8
       cam_out(c)%qbot(:,:)   = 0._r8
       cam_out(c)%pbot(:)     = 0._r8
       cam_out(c)%rho(:)      = 0._r8
       cam_out(c)%netsw(:)    = 0._r8
       cam_out(c)%flwds(:)    = 0._r8
       cam_out(c)%precsc(:)   = 0._r8
       cam_out(c)%precsl(:)   = 0._r8
       cam_out(c)%precc(:)    = 0._r8
       cam_out(c)%precl(:)    = 0._r8
       cam_out(c)%soll(:)     = 0._r8
       cam_out(c)%sols(:)     = 0._r8
       cam_out(c)%solld(:)    = 0._r8
       cam_out(c)%solsd(:)    = 0._r8
       cam_out(c)%srfrad(:)   = 0._r8
       cam_out(c)%thbot(:)    = 0._r8
       cam_out(c)%co2prog(:)  = 0._r8
       cam_out(c)%co2diag(:)  = 0._r8
       cam_out(c)%psl(:)      = 0._r8
!wangty modify
#ifdef wrf 
       cam_out(c)%clflo(:)      = 0._r8 ! juanxiong he
       cam_out(c)%clfmi(:)      = 0._r8 ! juanxiong he
       cam_out(c)%clfhi(:)      = 0._r8 ! juanxiong he
       cam_out(c)%pblh(:)      = 0._r8 ! juanxiong he       
#endif
       cam_out(c)%bcphidry(:) = 0._r8
       cam_out(c)%bcphodry(:) = 0._r8
       cam_out(c)%bcphiwet(:) = 0._r8
       cam_out(c)%ocphidry(:) = 0._r8
       cam_out(c)%ocphodry(:) = 0._r8
       cam_out(c)%ocphiwet(:) = 0._r8
       cam_out(c)%dstdry1(:)  = 0._r8
       cam_out(c)%dstwet1(:)  = 0._r8
       cam_out(c)%dstdry2(:)  = 0._r8
       cam_out(c)%dstwet2(:)  = 0._r8
       cam_out(c)%dstdry3(:)  = 0._r8
       cam_out(c)%dstwet3(:)  = 0._r8
       cam_out(c)%dstdry4(:)  = 0._r8
       cam_out(c)%dstwet4(:)  = 0._r8
    end do

  end subroutine atm2hub_alloc

  subroutine atm2hub_deallocate(cam_out)
    type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
    if(associated(cam_out)) then
       deallocate(cam_out)
    end if
    nullify(cam_out)

  end subroutine atm2hub_deallocate
  subroutine hub2atm_deallocate(cam_in)
    type(cam_in_t), pointer :: cam_in(:)    ! Atmosphere to surface input
    integer :: c

    if(associated(cam_in)) then
       do c=begchunk,endchunk
          if(associated(cam_in(c)%ram1)) then
             deallocate(cam_in(c)%ram1)
             nullify(cam_in(c)%ram1)
          end if
          if(associated(cam_in(c)%fv)) then
             deallocate(cam_in(c)%fv)
             nullify(cam_in(c)%fv)
          end if
          if(associated(cam_in(c)%depvel)) then
             deallocate(cam_in(c)%depvel)
             nullify(cam_in(c)%depvel)
          end if
          
       enddo

       deallocate(cam_in)
    end if
    nullify(cam_in)

  end subroutine hub2atm_deallocate



!======================================================================
! 
! BOP
!
! !IROUTINE: hub2atm_setopts
!
! !DESCRIPTION:
!
!   Method for outside packages to influence what is allocated
!   (For now, just aerosol dust controls if fv & ram1 arrays are allocated
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine hub2atm_setopts( aero_dust_in )
!
!!USES:
!
!
!!ARGUMENTS:
!
    logical, intent(in),optional :: aero_dust_in


!----------------------------------------------------------------------- 
! 
! EOP
!

    if ( present (aero_dust_in ) ) then
       dust = aero_dust_in
    endif

end subroutine hub2atm_setopts

end module camsrfexch_types
