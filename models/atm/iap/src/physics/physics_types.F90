!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
#include<xjb_drag3.inc>
!#define continuous 1
module physics_types

  use shr_kind_mod, only: r8 => shr_kind_r8
!====Jinbo Xie====
  !use ppgrid,       only: pcols, pver
   use ppgrid,       only: pcols, pver,pverp,nvar_dirOA,nvar_dirOL,indexb
!====Jinbo Xie====
#ifdef CCPP
  use constituents, only: pcnst, qmin, cnst_name, cnst_get_ind
#else
  use constituents, only: pcnst, qmin, cnst_name
#endif
  use geopotential, only: geopotential_dse
  use physconst,    only: zvir, gravit, cpair, rair
  use dycore,       only: dycore_is
  use phys_grid,    only: get_ncols_p, get_rlon_all_p, get_rlat_all_p, get_gcol_all_p
  use cam_logfile,  only: iulog
  use abortutils,   only: endrun
#ifdef CCPP
  use time_manager,   only: get_step_size
  use phys_control,   only: phys_deepconv_pbl, cam_physpkg_is, phys_getopts
  use error_messages, only: alloc_err
  use ccpp_types,     only: ccpp_t
#endif

  implicit none
  private          ! Make default type private to the module

  logical, parameter :: adjust_te = .FALSE.
  real(kind=r8), parameter :: zero      = 0.0_r8
  real(kind=r8), parameter :: clear_val = zero

!> \section arg_table_physics_types
!! \htmlinclude physics_types.html
!!

! Public types:

  public physics_state
  public physics_tend
  public physics_ptend
#ifdef CCPP
  public physics_int_ephem
  public physics_int_pers
  public physics_global
#endif

! Public interfaces

  public physics_update
  public physics_ptend_reset
  public physics_ptend_init
  public physics_state_set_grid
  public physics_dme_adjust  ! adjust dry mass and energy for change in water
                             ! cannot be applied to eul or sld dycores
  public physics_state_copy  ! copy a state type
  public physics_ptend_sum   ! add 2 ptend types
  public physics_tend_init   ! initialize a tend type

  public set_state_pdry      ! calculate dry air masses in state variable
  public set_wet_to_dry
  public set_dry_to_wet
  public physics_type_alloc
!-------------------------------------------------------------------------------
!! \section arg_table_physics_state
!! \htmlinclude physics_state.html
!!
  type physics_state
     ! yhy
    integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else =  pcols
    integer                                     :: &
          lchnk,   &! chunk index
          ncol      ! number of active columns
     real(r8), dimension(pcols)                  :: &
          lat,     &! latitude (radians)
          lon,     &! longitude (radians)
          ps,      &! surface pressure
          psdry,   &! dry surface pressure
          phis,    &! surface geopotential
!wangty modify
!#ifdef wrf  !zmh
          psl,     &! seal level pressure, added by juanxiong he
!#endif
          ulat,    &! unique latitudes  (radians)
          ulon      ! unique longitudes (radians)
     real(r8), dimension(pcols,pver)             :: &
          t,       &! temperature (K)
          u,       &! zonal wind (m/s)
          v,       &! meridional wind (m/s)
!wangty modify
!#ifdef wrf  !zmh
          rh,      &! relative humidity(%), added by juanxiong he
!#endif
          s,       &! dry static energy
          omega,   &! vertical pressure velocity (Pa/s)
          pmid,    &! midpoint pressure (Pa)
          pmiddry, &! midpoint pressure dry (Pa)
          pdel,    &! layer thickness (Pa)
          pdeldry, &! layer thickness dry (Pa)
          rpdel,   &! reciprocal of layer thickness (Pa)
          rpdeldry,&! recipricol layer thickness dry (Pa)
          lnpmid,  &! ln(pmid)
          lnpmiddry,&! log midpoint pressure dry (Pa)
          exner,   &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
          uzm,     &! zonal wind for qbo (m/s)
!wangty modify
#ifdef wrf
          zm,      &! geopotential height above surface at midpoints (m)
!------------------------------------------------------------------------
! Juanxiong He
!------------------------------------------------------------------------
          taucldv3d, &
          taucldi3d, &
          u_ndg_old, &  ! u nudging value at the previous time
          v_ndg_old, &  ! v nudging value at the previous time
          t_ndg_old, &  ! t nudging value at the previous time
          q_ndg_old, &  ! q nudging value at the previous time
          u_ndg_new, &  ! u nudging value
          v_ndg_new, &  ! v nudging value
          t_ndg_new, &  ! t nudging value
          q_ndg_new     ! q nudging value
!------------------------------------------------------------------------
! Juanxiong He
!------------------------------------------------------------------------
#else
          zm        ! geopotential height above surface at midpoints (m)
#endif
     real(r8), dimension(pcols,pver,pcnst) :: &
          q         ! constituent mixing ratio (kg/kg moist or dry air depending on type)

     real(r8), dimension(pcols,pver+1)           :: &
          pint,    &! interface pressure (Pa)
          pintdry, &! interface pressure dry (Pa)
          lnpint,  &! ln(pint)
          lnpintdry,&! log interface pressure dry (Pa)
          zi        ! geopotential height above surface at interfaces (m)

     real(r8), dimension(pcols)                  :: &
          te_ini,  &! vertically integrated total (kinetic + static) energy of initial state
          te_cur,  &! vertically integrated total (kinetic + static) energy of current state
          tw_ini,  &! vertically integrated total water of initial state
          tw_cur    ! vertically integrated total water of new state
     integer :: count ! count of values with significant energy or water imbalances
     integer, dimension(pcols) :: &
          latmapback, &! map from column to unique lat for that column
          lonmapback, &! map from column to unique lon for that column
          cid        ! unique column id
     integer :: ulatcnt, &! number of unique lats in chunk
                uloncnt   ! number of unique lons in chunk
!zmh
!czy 20180119 !zmh 20180920
!zmh stored vorticity and gradients in t, q at surface, see gw_drag.F90
!#if (defined WACCM_PHYS)
     real(r8), dimension(pcols,pver)             :: &
          frontgf, &  ! frontogenesis function
          frontga     ! frontogenesis angle
!#endif
      !czy 20180119 #endif ! zmh stored vorticity and gradients in t, q at surface, see gw_drag.F90

!      ! yhy
!      real(r8), dimension(pcols,pver)        :: &
!          var1,    &! q1
!          var2,    &! q2
!          newvar1,    &! new variables for high vertical resolution
!          newvar2,    &! new variables for high vertical resolution
!
!          prevt,   &! previous t
!          prevq ! previous q

!+czy20181120==================================================
!====Jinbo Xie===========
!added for 3d GWD oro par
     real(r8), dimension(pcols)                      :: &
          var,     &!standard deviation of high-res grid height
          oc        !convexity of high-res grid height
     real(r8), dimension(pcols,nvar_dirOA)           :: &
          oadir        !orographic asymmetry in a coarse grid
#ifndef continuous
     real(r8), dimension(pcols,nvar_dirOL)           :: &
          ol        !orographic length in a coarse grid
     !real(r8), dimension(pcols,nvar_dirOL)           :: &
          !dxydir    !representative grid length in a coarse grid
#else
     real(r8), dimension(4,pcols,indexb)             :: &
          terrout       !4 data related to determining ol
#endif
!====Jinbo Xie===========
!-czy20181120==================================================

  end type physics_state

!-------------------------------------------------------------------------------
  type physics_tend
     real(r8), dimension(pcols,pver)             :: dtdt, dudt, dvdt
     real(r8), dimension(pcols     )             :: flx_net
     real(r8), dimension(pcols)                  :: &
          te_tnd,  &! cumulative boundary flux of total energy
          tw_tnd    ! cumulative boundary flux of total water
  end type physics_tend

!-------------------------------------------------------------------------------
! This is for tendencies returned from individual parameterizations
  type physics_ptend
      ! yhy
     integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else =  pcols

     character*24 :: name    ! name of parameterization which produced tendencies.

     logical ::             &
          ls,               &! true if dsdt is returned
          lu,               &! true if dudt is returned
          lv,               &! true if dvdt is returned
          lq(pcnst)          ! true if dqdt() is returned

     integer ::             &
          top_level,        &! top level index for which nonzero tendencies have been set
          bot_level          ! bottom level index for which nonzero tendencies have been set

     real(r8), dimension(pcols,pver)             :: &
          s,                &! heating rate (J/kg/s)
          u,                &! u momentum tendency (m/s/s)
          v                  ! v momentum tendency (m/s/s)
     real(r8), dimension(pcols,pver,pcnst) :: &
          q                  ! consituent tendencies (kg/kg/s)

! boundary fluxes
     real(r8), dimension(pcols) ::&
          hflux_srf,     &! net heat flux at surface (W/m2)
          hflux_top,     &! net heat flux at top of model (W/m2)
          taux_srf,      &! net zonal stress at surface (Pa)
          taux_top,      &! net zonal stress at top of model (Pa)
          tauy_srf,      &! net meridional stress at surface (Pa)
          tauy_top        ! net meridional stress at top of model (Pa)
     real(r8), dimension(pcols,pcnst) ::&
          cflx_srf,      &! constituent flux at surface (kg/m2/s)
          cflx_top        ! constituent flux top of model (kg/m2/s)
!wangty modify
#ifdef wrf
!------------------------------------------------------------------------
! Juanxiong He
!------------------------------------------------------------------------
     real(r8), dimension(pcols,pver)             ::  rundgdten, &  ! u tendency
                                                     rvndgdten, &  ! v tendency
                                                     rtndgdten, &  ! t tendency
                                                     rqndgdten  ! q tendency
#endif
  end type physics_ptend

#ifdef CCPP
!! \section arg_table_physics_int_ephem
!! \htmlinclude physics_int_ephem.html
!!
  type physics_int_ephem

    !variables previously internal to tphysbc; there need to be n_threads of this type
    real(kind=r8), pointer              :: prec(:)  => null()  !<
    real(kind=r8), pointer              :: snow(:)  => null()  !<
    type(physics_ptend),   pointer      :: ptend_deep_conv_tot => null()
    real(kind=r8), pointer              :: cmfmc(:,:) => null()       !< Convective mass flux--m sub c
    real(kind=r8), pointer              :: cmfcme(:,:)=> null()           !< cmf condensation - evaporation
    real(kind=r8), pointer              :: dlf(:,:)   => null()           !< Detraining cld H20 from shallow + deep convections
    real(kind=r8), pointer              :: pflx(:,:)  => null()           !< Conv rain flux thru out btm of lev
    real(kind=r8), pointer              :: zdu(:,:)   => null()           !< detraining mass flux from deep convection
    real(kind=r8), pointer              :: rliq(:)    => null()         !< vertical integral of liquid not yet in q(ixcldliq)

    !variables previously internal to zm_conv_intr
    real(kind=r8), pointer              :: cape(:) => null() ! convective available potential energy.
    real(kind=r8), pointer              :: pcont(:) => null()
    real(kind=r8), pointer              :: pconb(:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_qv(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_s(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_evap_qv(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_evap_s(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_momtran_s(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_momtran_u(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_momtran_v(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_convtran_q(:,:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_tot_s(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_tot_u(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_tot_v(:,:) => null()
    real(kind=r8), pointer              :: ptend_deep_conv_tot_q(:,:,:) => null()
    real(kind=r8), pointer              :: temp_state_u(:,:) => null()
    real(kind=r8), pointer              :: temp_state_v(:,:) => null()
    real(kind=r8), pointer              :: temp_state_s(:,:) => null()
    real(kind=r8), pointer              :: temp_state_q(:,:,:) => null()
    real(kind=r8), pointer              :: temp_state_t(:,:) => null()
    real(kind=r8), pointer              :: temp_state_zm(:,:) => null()
    real(kind=r8), pointer              :: temp_state_zi(:,:) => null()
    real(kind=r8), pointer              :: tend_s_snwprd  (:,:) ! Heating rate of snow production
    real(kind=r8), pointer              :: tend_s_snwevmlt(:,:) ! Heating rate of evap/melting of snow
    real(kind=r8), pointer              :: ntprprd(:,:) !net precip production in layer
    real(kind=r8), pointer              :: ntsnprd(:,:) !net snow production in layer
    real(kind=r8), pointer              :: pguall(:,:,:)
    real(kind=r8), pointer              :: pgdall(:,:,:)
    real(kind=r8), pointer              :: icwu(:,:,:)
    real(kind=r8), pointer              :: icwd(:,:,:)

    contains
      procedure :: create      => interstitial_ephemeral_create     !<   allocate array data
      procedure :: reset       => interstitial_ephemeral_reset      !<   reset array data
  end type physics_int_ephem

!! \section arg_table_physics_int_pers
!! \htmlinclude physics_int_pers.html
!!
  type physics_int_pers

   !variables from the physics buffer; memory for these should already allocated as part of pbuf, but
   !associating the memory in this DDT will make metadata easier

   !needs to have an array of nchunks of these DDTs and loop through chunks to associate pointers with the right memory


   !added from convect_deep.F90/convect_deep_register
   real(kind=r8), pointer :: jctop(:) => null()
   real(kind=r8), pointer :: jcbot(:) => null()
   real(kind=r8), pointer :: rprd(:,:) => null()
   real(kind=r8), pointer :: ql(:,:) => null()
   real(kind=r8), pointer :: slflx(:,:) => null()
   real(kind=r8), pointer :: qtflx(:,:) => null()
   real(kind=r8), pointer :: evapcdp(:,:) => null()

   !accessed from zm_conv_intr.F90
   real(kind=r8), pointer :: cld_old(:,:) => null()

   !added from aerosol_intr.F90; accessed by zm_conv_intr.F90
   real(kind=r8), pointer :: fracis(:,:,:) => null()

   !added in zm_conv_intr.F90/zm_conv_register; access in zm_conv_tend
   real(kind=r8), pointer :: flxprec(:,:) => null()
   real(kind=r8), pointer :: flxsnow(:,:) => null()
   real(kind=r8), pointer :: dp_cldliq(:,:) => null()
   real(kind=r8), pointer :: dp_cldice(:,:) => null()

   !added from zm_conv_intr.F90 module variables
   real(kind=r8), pointer :: mu(:,:) => null()
   real(kind=r8), pointer :: md(:,:) => null()
   real(kind=r8), pointer :: du(:,:) => null()
   real(kind=r8), pointer :: eu(:,:) => null()
   real(kind=r8), pointer :: ed(:,:) => null()
   real(kind=r8), pointer :: dp(:,:) => null()
   real(kind=r8), pointer :: dsubcld(:) => null()
   integer      , pointer :: jt(:) => null()
   integer      , pointer :: maxg(:) => null()
   integer      , pointer :: ideep(:) => null()
   integer                :: lengath

   contains
     procedure :: associate       => interstitial_persistent_associate
     procedure :: create          => interstitial_persistent_create
     procedure :: init            => interstitial_persistent_init
  end type physics_int_pers

!! \section arg_table_physics_global
!! \htmlinclude physics_global.html
!!
  type physics_global

    character(len=16)     :: cam_physpkg
    character(len=16)     :: cam_physpkg_cam3
    character(len=16)     :: cam_physpkg_cam4
    character(len=16)     :: cam_physpkg_cam5
    character(len=16)     :: cam_physpkg_ideal
    character(len=16)     :: cam_physpkg_adiabatic
    character(len=16)     :: microp_scheme
    logical               :: non_dilute_buoy
    logical               :: no_deep_pbl
    logical               :: fv_dycore
    logical, dimension(2) :: l_windt
    integer               :: ixcldice
    integer               :: ixcldliq
    integer               :: ixnumice
    integer               :: ixnumliq
    real(kind=r8)         :: half_ztodt

    contains
      procedure :: init      => physics_global_init
  end type physics_global
#endif

!===============================================================================
contains
!===============================================================================
#ifdef CCPP
  subroutine physics_type_alloc(phys_state, phys_tend, phys_int_ephem, phys_int_pers, begchunk, endchunk)
#else
  subroutine physics_type_alloc(phys_state, phys_tend, begchunk, endchunk)
#endif
    use infnan, only : inf
    implicit none
    type(physics_state), pointer :: phys_state(:)
    type(physics_tend), pointer :: phys_tend(:)
    integer, intent(in) :: begchunk, endchunk
#ifdef CCPP
    type(physics_int_ephem), pointer :: phys_int_ephem(:)
    type(physics_int_pers), pointer :: phys_int_pers(:)
#endif
    integer :: ierr, lchnk
    type(physics_state), pointer :: state
    type(physics_tend), pointer :: tend

    allocate(phys_state(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_types: phys_state allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_state array')
    end if

    allocate(phys_tend(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_types: phys_tend allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_tend array')
    end if

#ifdef CCPP
    allocate(phys_int_ephem(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_types: phys_int_ephem allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_int_ephem array')
    end if

    write(0,'(a,2i6)') "Calling allocate(phys_int_pers(begchunk:endchunk)) with:", begchunk, endchunk
    allocate(phys_int_pers(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_types: phys_int_pers allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_int_pers array')
    end if
#endif

   ! Set chunk id, number of columns, and coordinates
    do lchnk = begchunk,endchunk
       state => phys_state(lchnk)
       tend => phys_tend(lchnk)
       state%lat(:) = inf
       state%lon(:) = inf
       state%ulat(:) = inf
       state%ulon(:) = inf
       state%ps(:) = inf
       ! DH*
       write(0,'(a,3i6)') "XXX: begchunk, endchunk, lchnk", begchunk, endchunk, lchnk
       write(0,'(a,3i6)') "XXX: pcols, ncol, size(lat)", pcols, state%ncol, size(state%lat)
       ! *DH
!wangty modify
#ifdef wrf
       state%psl(:) = inf ! added by juanxiong he
#endif
       state%psdry(:) = inf
       state%phis(:) = inf
       state%t(:,:) = inf
       state%u(:,:) = inf
       state%v(:,:) = inf
!wangty modify
#ifdef wrf
       state%rh(:,:) = inf ! added by juanxiong he
#endif
       state%s(:,:) = inf
       state%omega(:,:) = inf
       state%pmid(:,:) = inf
       state%pmiddry(:,:) = inf
       state%pdel(:,:) = inf
       state%pdeldry(:,:) = inf
       state%rpdel(:,:) = inf
       state%rpdeldry(:,:) = inf
       state%lnpmid(:,:) = inf
       state%lnpmiddry(:,:) = inf
       state%exner(:,:) = inf
       state%uzm(:,:) = inf
       state%zm(:,:) = inf
       state%q(:,:,:) = inf

       state%pint(:,:) = inf
       state%pintdry(:,:) = inf
       state%lnpint(:,:) = inf
       state%lnpintdry(:,:) = inf
       state%zi(:,:) = inf

       state%te_ini(:) = inf
       state%te_cur(:) = inf
       state%tw_ini(:) = inf
       state%tw_cur(:) = inf
!+czy20181120==================================================
!====Jinbo Xie inf====
state%var(:)=inf
state%oc(:)=inf
state%oadir(:,:)=inf
#ifndef continuous
state%ol(:,:)=inf
!state%dxydir(:,:)=inf
#else
state%terrout=inf
#endif
!====Jinbo Xie inf====
!-czy20181120==================================================

       tend%dtdt(:,:) = inf
       tend%dudt(:,:) = inf
       tend%dvdt(:,:) = inf
       tend%flx_net(:) = inf
       tend%te_tnd(:) = inf
       tend%tw_tnd(:) = inf
!wangty modify
#ifdef wrf
!------------------------------------------------------------------------
! Juanxiong He
!------------------------------------------------------------------------
       state%taucldv3d(:,:)=0._r8
       state%taucldi3d(:,:)=0._r8
       state%u_ndg_old(:,:)=0._r8
       state%v_ndg_old(:,:)=0._r8
       state%t_ndg_old(:,:)=0._r8
       state%q_ndg_old(:,:)=0._r8
       state%u_ndg_new(:,:)=0._r8
       state%v_ndg_new(:,:)=0._r8
       state%t_ndg_new(:,:)=0._r8
       state%q_ndg_new(:,:)=0._r8
!------------------------------------------------------------------------
! Juanxiong He
!------------------------------------------------------------------------
#endif
#ifdef CCPP
       !call create for each block
       write(0,'(a,3i6)') "Calling phys_int_pers % create with:", lchnk, pcols, pver
       call phys_int_ephem(lchnk)%create(pcols, pver, pverp, pcnst)
       call phys_int_pers (lchnk)%create(pcols, pver)
#endif
    end do

  end subroutine physics_type_alloc
!===============================================================================
  subroutine physics_update(state, tend, ptend, dt)
!-----------------------------------------------------------------------
! Update the state and or tendency structure with the parameterization tendencies
!-----------------------------------------------------------------------
    use geopotential, only: geopotential_dse
    use constituents, only: cnst_get_ind
    use scamMod,      only: scm_crm_mode, single_column
    use phys_control, only: phys_getopts

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies

    type(physics_state), intent(inout)  :: state   ! Physics state variables
    type(physics_tend ), intent(inout)  :: tend    ! Physics tendencies

    real(r8), intent(in) :: dt                     ! time step
!
!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    integer :: ixnumice, ixnumliq
    integer :: ncol                                ! number of columns
    character*40 :: name    ! param and tracer name for qneg3

    character(len=16) :: microp_scheme='RK'  ! microphysics scheme
    !-----------------------------------------------------------------------

    ! The column radiation model does not update the state
    if(single_column.and.scm_crm_mode) return

    call phys_getopts(microp_scheme_out=microp_scheme)

    ncol = state%ncol

    ! Update u,v fields
    if(ptend%lu) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%u  (i,k) = state%u  (i,k) + ptend%u(i,k) * dt
             tend%dudt(i,k) = tend%dudt(i,k) + ptend%u(i,k)
          end do
       end do
    end if

    if(ptend%lv) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%v  (i,k) = state%v  (i,k) + ptend%v(i,k) * dt
             tend%dvdt(i,k) = tend%dvdt(i,k) + ptend%v(i,k)
          end do
       end do
    end if

    ! Update dry static energy
    if(ptend%ls) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%s(i,k)   = state%s(i,k)   + ptend%s(i,k) * dt
             tend%dtdt(i,k) = tend%dtdt(i,k) + ptend%s(i,k)/cpair
          end do
       end do
    end if

    ! Update constituents, all schemes use time split q: no tendency kept
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    ! Check for number concentration of cloud liquid and cloud ice (if not present
    ! the indices will be set to -1)
    call cnst_get_ind('NUMICE', ixnumice, abort=.false.)
    call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)

    do m = 1, pcnst
       if(ptend%lq(m)) then
          do k = ptend%top_level, ptend%bot_level
             do i = 1,ncol
                state%q(i,k,m) = state%q(i,k,m) + ptend%q(i,k,m) * dt
             end do
          end do

          ! now test for mixing ratios which are too small
          ! don't call qneg3 for number concentration variables
          if (m .ne. ixnumice  .and.  m .ne. ixnumliq) then
             name = trim(ptend%name) // '/' // trim(cnst_name(m))
             call qneg3(trim(name), state%lchnk, ncol, pcols, pver, m, m, qmin(m), state%q(1,1,m))
          else
             do k = ptend%top_level, ptend%bot_level
                do i = 1,ncol
                   ! checks for number concentration
                   state%q(i,k,m) = max(1.e-12_r8,state%q(i,k,m))
                   state%q(i,k,m) = min(1.e10_r8,state%q(i,k,m))
                end do
             end do
          end if

       end if
    end do

    ! special tests for cloud liquid
    if (ixcldliq > 1) then
       if(ptend%lq(ixcldliq)) then
          if (ptend%name == 'stratiform' .or. ptend%name == 'cldwat'  ) then
#ifdef PERGRO
             where (state%q(:ncol,:pver,ixcldliq) < 1.e-12_r8)
                state%q(:ncol,:pver,ixcldliq) = 0._r8
             end where

             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' .or. microp_scheme .eq. 'M2M') then
                where (state%q(:ncol,:pver,ixcldliq) < 1.e-12_r8)
                   state%q(:ncol,:pver,ixnumliq) = 0._r8
                end where
             end if
#endif
          else if (ptend%name == 'convect_deep') then
             where (state%q(:ncol,:pver,ixcldliq) < 1.e-36_r8)
                state%q(:ncol,:pver,ixcldliq) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' .or. microp_scheme .eq. 'M2M' ) then
                where (state%q(:ncol,:pver,ixcldliq) < 1.e-36_r8)
                   state%q(:ncol,:pver,ixnumliq) = 0._r8
                end where
             end if
          end if
       end if
    end if

    ! special tests for cloud ice
    if (ixcldice > 1) then
       if(ptend%lq(ixcldice)) then
          if (ptend%name == 'stratiform' .or. ptend%name == 'cldwat'  ) then
#ifdef PERGRO
             where (state%q(:ncol,:pver,ixcldice) < 1.e-12_r8)
                state%q(:ncol,:pver,ixcldice) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' .or. microp_scheme .eq. 'M2M') then
                where (state%q(:ncol,:pver,ixcldice) < 1.e-12_r8)
                   state%q(:ncol,:pver,ixnumice) = 0._r8
                end where
             end if
#endif
          else if (ptend%name == 'convect_deep') then
             where (state%q(:ncol,:pver,ixcldice) < 1.e-36_r8)
                state%q(:ncol,:pver,ixcldice) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' .or. microp_scheme .eq. 'M2M' ) then
                where (state%q(:ncol,:pver,ixcldice) < 1.e-36_r8)
                   state%q(:ncol,:pver,ixnumice) = 0._r8
                end where
             end if
          end if
       end if
    end if

    ! Derive new temperature and geopotential fields if heating or water tendency not 0.
    if (ptend%ls .or. ptend%lq(1)) then
       call geopotential_dse(                                                                    &
            state%lnpint, state%lnpmid, state%pint  , state%pmid  , state%pdel  , state%rpdel  , &
            state%s     , state%q(1,1,1),state%phis , rair        , gravit      , cpair        , &
            zvir        , state%t     , state%zi    , state%zm    , ncol         )
    end if

       ! juanxiong he
#ifdef CO2
       do i=1,ncol
       do k=1,pver
       if(isnan(state%zm(i,k))) print *,'s=',i,k,state%s(i,k)
       end do
       end do
#endif

    ! Reset all parameterization tendency flags to false
    call physics_ptend_reset(ptend)

  end subroutine physics_update

!===============================================================================

  subroutine physics_ptend_sum(ptend, ptend_sum, state)
!-----------------------------------------------------------------------
! Add ptend fields to ptend_sum for ptend logical flags = .true.
! Where ptend logical flags = .false, don't change ptend_sum
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(in)     :: ptend   ! New parameterization tendencies
    type(physics_ptend), intent(inout)  :: ptend_sum   ! Sum of incoming ptend_sum and ptend
    type(physics_state), intent(in)     :: state   ! New parameterization tendencies

!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ncol                                ! number of columns

!-----------------------------------------------------------------------
    ncol = state%ncol


! Update u,v fields
    if(ptend%lu) then
       ptend_sum%lu = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%u(i,k) = ptend_sum%u(i,k) + ptend%u(i,k)
          end do
          ptend_sum%taux_srf(i) = ptend_sum%taux_srf(i) + ptend%taux_srf(i)
          ptend_sum%taux_top(i) = ptend_sum%taux_top(i) + ptend%taux_top(i)
       end do
    end if

    if(ptend%lv) then
       ptend_sum%lv = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%v(i,k) = ptend_sum%v(i,k) + ptend%v(i,k)
          end do
          ptend_sum%tauy_srf(i) = ptend_sum%tauy_srf(i) + ptend%tauy_srf(i)
          ptend_sum%tauy_top(i) = ptend_sum%tauy_top(i) + ptend%tauy_top(i)
       end do
    end if


    if(ptend%ls) then
       ptend_sum%ls = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%s(i,k) = ptend_sum%s(i,k) + ptend%s(i,k)
          end do
          ptend_sum%hflux_srf(i) = ptend_sum%hflux_srf(i) + ptend%hflux_srf(i)
          ptend_sum%hflux_top(i) = ptend_sum%hflux_top(i) + ptend%hflux_top(i)
       end do
    end if

! Update constituents
    do m = 1, pcnst
       if(ptend%lq(m)) then
          ptend_sum%lq(m) = .true.
          do i = 1,ncol
             do k = ptend%top_level, ptend%bot_level
                ptend_sum%q(i,k,m) = ptend_sum%q(i,k,m) + ptend%q(i,k,m)
             end do
             ptend_sum%cflx_srf(i,m) = ptend_sum%cflx_srf(i,m) + ptend%cflx_srf(i,m)
             ptend_sum%cflx_top(i,m) = ptend_sum%cflx_top(i,m) + ptend%cflx_top(i,m)
          end do
       end if
    end do


  end subroutine physics_ptend_sum

!===============================================================================
  subroutine physics_ptend_reset(ptend)
!-----------------------------------------------------------------------
! Reset the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    integer :: m             ! Index for constiuent
!-----------------------------------------------------------------------

    if(ptend%ls) then
       ptend%s = 0._r8
       ptend%hflux_srf = 0._r8
       ptend%hflux_top = 0._r8
    endif
    if(ptend%lu) then
       ptend%u = 0._r8
       ptend%taux_srf = 0._r8
       ptend%taux_top = 0._r8
    endif
    if(ptend%lv) then
       ptend%v = 0._r8
       ptend%tauy_srf = 0._r8
       ptend%tauy_top = 0._r8
    endif
    do m = 1, pcnst
       if(ptend%lq(m)) then
          ptend%q(:,:,m) = 0._r8
          ptend%cflx_srf(:,m) = 0._r8
          ptend%cflx_top(:,m) = 0._r8
       endif
    end do

    ptend%name  = "none"
    ptend%lq(:) = .FALSE.
    ptend%ls    = .FALSE.
    ptend%lu    = .FALSE.
    ptend%lv    = .FALSE.

    ptend%top_level = 1
    ptend%bot_level = pver
!wangty modify
#ifdef wrf
!------------------------------------------------------------------------
! Juanxiong He
!------------------------------------------------------------------------
    ptend%rundgdten(:,:)=0
    ptend%rvndgdten(:,:)=0
    ptend%rtndgdten(:,:)=0
    ptend%rqndgdten(:,:)=0
#endif
    return
  end subroutine physics_ptend_reset

!===============================================================================
  subroutine physics_ptend_init(ptend, psetcols, name, ls, lu, lv, lq)
!-----------------------------------------------------------------------
! Initialize the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
    integer, intent(in), optional                 :: psetcols ! maximum number of columns
    character(len=*), optional                    :: name     ! optional name of parameterization which produced tendencies.
    logical, optional                   :: ls       ! if true, then fields to support dsdt are allocated
    logical, optional                   :: lu       ! if true, then fields to support dudt are allocated
    logical, optional                   :: lv       ! if true, then fields to support dvdt are allocated
    logical, dimension(pcnst),optional  :: lq       ! if true, then fields to support dqdt are allocated
!-----------------------------------------------------------------------
    if (present(psetcols)) then
        ptend%psetcols = psetcols
    else
        ptend%psetcols = 0
    end if

    if (present(name)) then
        ptend%name = name
    else
        ptend%name  = "none"
    end if

    if (present(lq)) then
        ptend%lq(:) = lq(:)
    else
        ptend%lq(:) = .true.
    end if

    if (present(ls)) then
        ptend%ls = ls
    else
        ptend%ls    = .true.
    end if

    if (present(lu)) then
        ptend%lu = lu
    else
        ptend%lu    = .true.
    end if

    if (present(lv)) then
        ptend%lv = lv
    else
        ptend%lv    = .true.
    end if

    call physics_ptend_reset(ptend)

    return
  end subroutine physics_ptend_init

!===============================================================================

  subroutine physics_state_set_grid(lchnk, phys_state)
!-----------------------------------------------------------------------
! Set the grid components of the physics_state object
!-----------------------------------------------------------------------

    integer,             intent(in)    :: lchnk
    type(physics_state), intent(inout) :: phys_state

    ! local variables
    integer  :: i, ncol
    real(r8) :: rlon(pcols)
    real(r8) :: rlat(pcols)
    !-----------------------------------------------------------------------

    ncol = get_ncols_p(lchnk)

    if(ncol<=0) then
       write(iulog,*) lchnk, ncol
       call endrun('physics_state_set_grid')
    end if

    call get_rlon_all_p(lchnk, ncol, rlon)
    call get_rlat_all_p(lchnk, ncol, rlat)
    phys_state%ncol  = ncol
    phys_state%lchnk = lchnk
    do i=1,ncol
       phys_state%lat(i) = rlat(i)
       phys_state%lon(i) = rlon(i)
    end do
    call init_geo_unique(phys_state,ncol)

  end subroutine physics_state_set_grid


  subroutine init_geo_unique(phys_state,ncol)
    integer,             intent(in)    :: ncol
    type(physics_state), intent(inout) :: phys_state
    logical :: match
    integer :: i, j, ulatcnt, uloncnt

    phys_state%ulat=-999.
    phys_state%ulon=-999.
    phys_state%latmapback=0
    phys_state%lonmapback=0
    match=.false.
    ulatcnt=0
    uloncnt=0
    match=.false.
    do i=1,ncol
       do j=1,ulatcnt
          if(phys_state%lat(i) .eq. phys_state%ulat(j)) then
             match=.true.
             phys_state%latmapback(i)=j
          end if
       end do
       if(.not. match) then
          ulatcnt=ulatcnt+1
          phys_state%ulat(ulatcnt)=phys_state%lat(i)
          phys_state%latmapback(i)=ulatcnt
       end if

       match=.false.
       do j=1,uloncnt
          if(phys_state%lon(i) .eq. phys_state%ulon(j)) then
             match=.true.
             phys_state%lonmapback(i)=j
          end if
       end do
       if(.not. match) then
          uloncnt=uloncnt+1
          phys_state%ulon(uloncnt)=phys_state%lon(i)
          phys_state%lonmapback(i)=uloncnt
       end if
       match=.false.

    end do
    phys_state%uloncnt=uloncnt
    phys_state%ulatcnt=ulatcnt

    call get_gcol_all_p(phys_state%lchnk,pcols,phys_state%cid)


  end subroutine init_geo_unique

!===============================================================================
  subroutine physics_dme_adjust(state, tend, qini, dt)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Adjust the dry mass in each layer back to the value of physics input state
    !
    ! Method: Conserve the integrated mass, momentum and total energy in each layer
    !         by scaling the specific mass of consituents, specific momentum (velocity)
    !         and specific total energy by the relative change in layer mass. Solve for
    !         the new temperature by subtracting the new kinetic energy from total energy
    !         and inverting the hydrostatic equation
    !
    !         The mass in each layer is modified, changing the relationship of the layer
    !         interfaces and midpoints to the surface pressure. The result is no longer in
    !         the original hybrid coordinate.
    !
    !         This procedure cannot be applied to the "eul" or "sld" dycores because they
    !         require the hybrid coordinate.
    !
    ! Author: Byron Boville

    ! !REVISION HISTORY:
    !   03.03.28  Boville    Created, partly from code by Lin in p_d_adjust
    !
    !-----------------------------------------------------------------------

    use constituents, only : cnst_get_type_byind

    implicit none
    !
    ! Arguments
    !
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    real(r8),            intent(in   ) :: qini(pcols,pver)    ! initial specific humidity
    real(r8),            intent(in   ) :: dt                  ! model physics timestep
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: lchnk         ! chunk identifier
    integer  :: ncol          ! number of atmospheric columns
    integer  :: i,k,m         ! Longitude, level indices
    real(r8) :: fdq(pcols)    ! mass adjustment factor
    real(r8) :: te(pcols)     ! total energy in a layer
    real(r8) :: utmp(pcols)   ! temp variable for recalculating the initial u values
    real(r8) :: vtmp(pcols)   ! temp variable for recalculating the initial v values
    !
    !-----------------------------------------------------------------------
    ! verify that the dycore is FV
    if (.not. dycore_is('LR') ) return

    lchnk = state%lchnk
    ncol  = state%ncol

    ! adjust dry mass in each layer back to input value, while conserving
    ! constituents, momentum, and total energy
    do k = 1, pver

       ! adjusment factor is just change in water vapor
       fdq(:ncol) = 1._r8 + state%q(:ncol,k,1) - qini(:ncol,k)

       ! adjust constituents to conserve mass in each layer
       do m = 1, pcnst
          state%q(:ncol,k,m) = state%q(:ncol,k,m) / fdq(:ncol)
       end do

       if (adjust_te) then
          ! compute specific total energy of unadjusted state (J/kg)
          te(:ncol) = state%s(:ncol,k) + 0.5_r8*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2)

          ! recompute initial u,v from the new values and the tendencies
          utmp(:ncol) = state%u(:ncol,k) - dt * tend%dudt(:ncol,k)
          vtmp(:ncol) = state%v(:ncol,k) - dt * tend%dvdt(:ncol,k)
          ! adjust specific total energy and specific momentum (velocity) to conserve each
          te     (:ncol)   = te     (:ncol)     / fdq(:ncol)
          state%u(:ncol,k) = state%u(:ncol,k  ) / fdq(:ncol)
          state%v(:ncol,k) = state%v(:ncol,k  ) / fdq(:ncol)
          ! compute adjusted u,v tendencies
          tend%dudt(:ncol,k) = (state%u(:ncol,k) - utmp(:ncol)) / dt
          tend%dvdt(:ncol,k) = (state%v(:ncol,k) - vtmp(:ncol)) / dt

          ! compute adjusted static energy
          state%s(:ncol,k) = te(:ncol) - 0.5_r8*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2)
       end if

! compute new total pressure variables
       state%pdel  (:ncol,k  ) = state%pdel(:ncol,k  ) * fdq(:ncol)
       state%pint  (:ncol,k+1) = state%pint(:ncol,k  ) + state%pdel(:ncol,k)
       state%lnpint(:ncol,k+1) = log(state%pint(:ncol,k+1))
       state%rpdel (:ncol,k  ) = 1._r8/ state%pdel(:ncol,k  )
    end do

! compute new T,z from new s,q,dp
    if (adjust_te) then
       call geopotential_dse(state%lnpint, state%lnpmid  , state%pint ,  &
            state%pmid  , state%pdel    , state%rpdel,  &
            state%s     , state%q(1,1,1), state%phis , rair, gravit, cpair, zvir, &
            state%t     , state%zi      , state%zm   , ncol      )
    end if

  end subroutine physics_dme_adjust
!-----------------------------------------------------------------------

!===============================================================================
  subroutine physics_state_copy(state_in, state_out)

    use ppgrid,       only: pver, pverp
    use constituents, only: pcnst

    implicit none

    !
    ! Arguments
    !
    type(physics_state), intent(in) :: state_in
    type(physics_state), intent(out) :: state_out

    !
    ! Local variables
    !
    integer i, k, m, ncol


    ncol = state_in%ncol

    state_out%lchnk = state_in%lchnk
    state_out%ncol  = state_in%ncol
    state_out%count = state_in%count

    do i = 1, ncol
       state_out%lat(i)    = state_in%lat(i)
       state_out%lon(i)    = state_in%lon(i)
       state_out%ps(i)     = state_in%ps(i)
       state_out%phis(i)   = state_in%phis(i)
       state_out%te_ini(i) = state_in%te_ini(i)
       state_out%te_cur(i) = state_in%te_cur(i)
       state_out%tw_ini(i) = state_in%tw_ini(i)
       state_out%tw_cur(i) = state_in%tw_cur(i)
    end do

    do k = 1, pver
       do i = 1, ncol
          state_out%t(i,k)         = state_in%t(i,k)
          state_out%u(i,k)         = state_in%u(i,k)
          state_out%v(i,k)         = state_in%v(i,k)
          state_out%s(i,k)         = state_in%s(i,k)
          state_out%omega(i,k)     = state_in%omega(i,k)
          state_out%pmid(i,k)      = state_in%pmid(i,k)
          state_out%pdel(i,k)      = state_in%pdel(i,k)
          state_out%rpdel(i,k)     = state_in%rpdel(i,k)
          state_out%lnpmid(i,k)    = state_in%lnpmid(i,k)
          state_out%exner(i,k)     = state_in%exner(i,k)
          state_out%zm(i,k)        = state_in%zm(i,k)
       end do
    end do

    do k = 1, pverp
       do i = 1, ncol
          state_out%pint(i,k)      = state_in%pint(i,k)
          state_out%lnpint(i,k)    = state_in%lnpint(i,k)
          state_out%zi(i,k)        = state_in% zi(i,k)
       end do
    end do


       do i = 1, ncol
          state_out%psdry(i)  = state_in%psdry(i)
       end do
       do k = 1, pver
          do i = 1, ncol
             state_out%lnpmiddry(i,k) = state_in%lnpmiddry(i,k)
             state_out%pmiddry(i,k)   = state_in%pmiddry(i,k)
             state_out%pdeldry(i,k)   = state_in%pdeldry(i,k)
             state_out%rpdeldry(i,k)  = state_in%rpdeldry(i,k)
          end do
       end do
       do k = 1, pverp
          do i = 1, ncol
             state_out%pintdry(i,k)   = state_in%pintdry(i,k)
             state_out%lnpintdry(i,k) = state_in%lnpintdry(i,k)
          end do
       end do

    do m = 1, pcnst
       do k = 1, pver
          do i = 1, ncol
             state_out%q(i,k,m) = state_in%q(i,k,m)
          end do
       end do
    end do

  end  subroutine physics_state_copy
!===============================================================================

  subroutine physics_tend_init(tend)

    implicit none

    !
    ! Arguments
    !
    type(physics_tend), intent(inout) :: tend

    !
    ! Local variables
    !

    tend%dtdt    = 0._r8
    tend%dudt    = 0._r8
    tend%dvdt    = 0._r8
    tend%flx_net = 0._r8
    tend%te_tnd  = 0._r8
    tend%tw_tnd  = 0._r8

end subroutine physics_tend_init

!===============================================================================

subroutine set_state_pdry (state,pdeld_calc)

  use ppgrid,  only: pver
  use pmgrid,  only: plev, plevp
  use hycoef,  only: hyai, hybi, ps0
  implicit none

  type(physics_state), intent(inout) :: state
  logical, optional, intent(in) :: pdeld_calc    !  .true. do calculate pdeld [default]
                                                 !  .false. don't calculate pdeld
  integer ncol
  integer i, k
  logical do_pdeld_calc

  if ( present(pdeld_calc) ) then
     do_pdeld_calc = pdeld_calc
  else
     do_pdeld_calc = .true.
  endif

  ncol = state%ncol

  state%psdry(:ncol) = ps0 * hyai(1) + state%ps(:ncol) * hybi(1)
  state%pintdry(:ncol,1) = ps0 * hyai(1) + state%ps(:ncol) * hybi(1)

  if (do_pdeld_calc)  then
     do k = 1, pver
        state%pdeldry(:ncol,k) = state%pdel(:ncol,k)*(1._r8-state%q(:ncol,k,1))
     end do
  endif
  do k = 1, pver
     state%pintdry(:ncol,k+1) = state%pintdry(:ncol,k)+state%pdeldry(:ncol,k)
     state%pmiddry(:ncol,k) = (state%pintdry(:ncol,k+1)+state%pintdry(:ncol,k))/2._r8
     state%psdry(:ncol) = state%psdry(:ncol) + state%pdeldry(:ncol,k)
  end do

  state%rpdeldry(:ncol,:) = 1._r8/state%pdeldry(:ncol,:)
  state%lnpmiddry(:ncol,:) = log(state%pmiddry(:ncol,:))
  state%lnpintdry(:ncol,:) = log(state%pintdry(:ncol,:))

end subroutine set_state_pdry

!===============================================================================

subroutine set_wet_to_dry (state)

  use constituents,  only: pcnst, cnst_type

  type(physics_state), intent(inout) :: state

  integer m, ncol

  ncol = state%ncol

  do m = 1,pcnst
     if (cnst_type(m).eq.'dry') then
        state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdel(:ncol,:)/state%pdeldry(:ncol,:)
     endif
  end do

end subroutine set_wet_to_dry

!===============================================================================

subroutine set_dry_to_wet (state)

  use constituents,  only: pcnst, cnst_type

  type(physics_state), intent(inout) :: state

  integer m, ncol

  ncol = state%ncol

  do m = 1,pcnst
     if (cnst_type(m).eq.'dry') then
        state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdeldry(:ncol,:)/state%pdel(:ncol,:)
     endif
  end do

end subroutine set_dry_to_wet

#ifdef CCPP
subroutine interstitial_ephemeral_create (int_ephem, ncol, pver, pverp, pcnst)
  implicit none

  class(physics_int_ephem)       :: int_ephem
  integer,                intent(in) :: ncol, pver, pverp, pcnst

  allocate (int_ephem%prec  (ncol))
  allocate (int_ephem%snow  (ncol))
  allocate (int_ephem%cmfmc (ncol, pverp))
  allocate (int_ephem%cmfcme(ncol, pver))
  allocate (int_ephem%cape  (ncol))
  allocate (int_ephem%dlf   (ncol, pver))
  allocate (int_ephem%pflx  (ncol, pverp))
  allocate (int_ephem%zdu   (ncol, pver))
  allocate (int_ephem%rliq  (ncol))
  allocate (int_ephem%pcont (ncol))
  allocate (int_ephem%pconb (ncol))
  allocate (int_ephem%ptend_deep_conv_qv(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_s(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_evap_qv(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_evap_s(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_momtran_s(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_momtran_u(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_momtran_v(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_convtran_q(ncol,pver,pcnst))
  allocate (int_ephem%ptend_deep_conv_tot_s(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_tot_u(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_tot_v(ncol,pver))
  allocate (int_ephem%ptend_deep_conv_tot_q(ncol,pver,pcnst))
  allocate (int_ephem%temp_state_u(ncol,pver))
  allocate (int_ephem%temp_state_v(ncol,pver))
  allocate (int_ephem%temp_state_s(ncol,pver))
  allocate (int_ephem%temp_state_q(ncol,pver,pcnst))
  allocate (int_ephem%temp_state_t(ncol,pver))
  allocate (int_ephem%temp_state_zm(ncol,pver))
  allocate (int_ephem%temp_state_zi(ncol,pverp))
  allocate (int_ephem%tend_s_snwprd(ncol,pver))
  allocate (int_ephem%tend_s_snwevmlt(ncol,pver))
  allocate (int_ephem%ntprprd(ncol,pver))
  allocate (int_ephem%ntsnprd(ncol,pver))
  allocate (int_ephem%pguall(ncol,pver,2))
  allocate (int_ephem%pgdall(ncol,pver,2))
  allocate (int_ephem%icwu(ncol,pver,2))
  allocate (int_ephem%icwd(ncol,pver,2))

end subroutine interstitial_ephemeral_create

subroutine interstitial_ephemeral_reset(int_ephem)
  implicit none

  class(physics_int_ephem)       :: int_ephem

  int_ephem%prec            = clear_val
  int_ephem%snow            = clear_val
  int_ephem%cmfmc           = clear_val
  int_ephem%cmfcme          = clear_val
  int_ephem%cape            = clear_val
  int_ephem%dlf             = clear_val
  int_ephem%pflx            = clear_val
  int_ephem%zdu             = clear_val
  int_ephem%rliq            = clear_val
  int_ephem%pcont           = clear_val
  int_ephem%pconb           = clear_val
  int_ephem%ptend_deep_conv_qv = clear_val
  int_ephem%ptend_deep_conv_s  = clear_val
  int_ephem%ptend_deep_conv_evap_qv = clear_val
  int_ephem%ptend_deep_conv_evap_s  = clear_val
  int_ephem%ptend_deep_conv_momtran_s  = clear_val
  int_ephem%ptend_deep_conv_momtran_u  = clear_val
  int_ephem%ptend_deep_conv_momtran_v  = clear_val
  int_ephem%ptend_deep_conv_convtran_q = clear_val
  int_ephem%ptend_deep_conv_tot_s  = clear_val
  int_ephem%ptend_deep_conv_tot_u  = clear_val
  int_ephem%ptend_deep_conv_tot_v  = clear_val
  int_ephem%ptend_deep_conv_tot_q = clear_val
  int_ephem%temp_state_u    = clear_val
  int_ephem%temp_state_v    = clear_val
  int_ephem%temp_state_s    = clear_val
  int_ephem%temp_state_q    = clear_val
  int_ephem%temp_state_t    = clear_val
  int_ephem%temp_state_zm   = clear_val
  int_ephem%temp_state_zi   = clear_val
  int_ephem%tend_s_snwprd   = clear_val
  int_ephem%tend_s_snwevmlt = clear_val
  int_ephem%ntprprd         = clear_val
  int_ephem%ntsnprd         = clear_val
  int_ephem%pguall          = clear_val
  int_ephem%pgdall          = clear_val
  int_ephem%icwu            = clear_val
  int_ephem%icwd            = clear_val

end subroutine interstitial_ephemeral_reset

subroutine interstitial_persistent_associate(int_pers, pcols, pver, pverp, pcnst, lchnk)
  !should be called once per block before physics init stage

  use phys_buffer,   only: pbuf, pbuf_get_fld_idx, pbuf_old_tim_idx
  use buffer,        only: pblht, tpert, qpert, tpert2, qpert2
  implicit none

  class(physics_int_pers)       :: int_pers

  integer, intent(in) :: pcols, pver, pverp, pcnst, lchnk

  integer :: var_idx, itim

  !from convect_deep.F90/convect_deep_init,convect_deep_tend
  var_idx = pbuf_get_fld_idx('CLDTOP')
  int_pers%jctop => pbuf(var_idx)%fld_ptr(1,1:pcols,1,lchnk,1)
  var_idx = pbuf_get_fld_idx('CLDBOT')
  int_pers%jcbot => pbuf(var_idx)%fld_ptr(1,1:pcols,1,lchnk,1)
  var_idx = pbuf_get_fld_idx('RPRDDP')
  int_pers%rprd => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
  var_idx = pbuf_get_fld_idx('ICWMRDP')
  int_pers%ql => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
  var_idx = pbuf_get_fld_idx('slflxdp')
  int_pers%slflx => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)
  var_idx = pbuf_get_fld_idx('qtflxdp')
  int_pers%qtflx => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)
  var_idx = pbuf_get_fld_idx('NEVAPR_DPCU')
  int_pers%evapcdp => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

  !from zm_conv_intr.F90/zm_conv_init/tend
  var_idx = pbuf_get_fld_idx('CLD')
  itim    = pbuf_old_tim_idx()
  int_pers%cld_old => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
  var_idx = pbuf_get_fld_idx('FRACIS')
  int_pers%fracis  => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1:pcnst)

  var_idx = pbuf_get_fld_idx('DP_FLXPRC')
  int_pers%flxprec => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)
  var_idx = pbuf_get_fld_idx('DP_FLXSNW')
  int_pers%flxsnow => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)
  var_idx = pbuf_get_fld_idx('DP_CLDLIQ')
  int_pers%dp_cldliq => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
  var_idx = pbuf_get_fld_idx('DP_CLDICE')
  int_pers%dp_cldice => pbuf(var_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

end subroutine interstitial_persistent_associate

subroutine interstitial_persistent_create(int_pers, pcols, pver)
  implicit none

  class(physics_int_pers)       :: int_pers

  integer, intent(in) :: pcols, pver

  integer :: istat

  allocate (int_pers%mu(pcols, pver), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'mu', pcols*pver)
  allocate (int_pers%md(pcols, pver), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'md', pcols*pver)
  allocate (int_pers%du(pcols, pver), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'du', pcols*pver)
  allocate (int_pers%eu(pcols, pver), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'eu', pcols*pver)
  allocate (int_pers%ed(pcols, pver), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'ed', pcols*pver)
  allocate (int_pers%dp(pcols, pver), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'dp', pcols*pver)
  allocate (int_pers%dsubcld(pcols), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'dsubcld', pcols)
  allocate (int_pers%jt(pcols), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'jt', pcols)
  allocate (int_pers%maxg(pcols), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'maxg', pcols)
  allocate (int_pers%ideep(pcols), stat=istat)
  call alloc_err( istat, 'interstitial_persistent_create', 'ideep', pcols)

end subroutine interstitial_persistent_create

subroutine interstitial_persistent_init(int_pers)
  implicit none

  class(physics_int_pers)       :: int_pers

  int_pers%mu = clear_val
  int_pers%md = clear_val
  int_pers%du = clear_val
  int_pers%eu = clear_val
  int_pers%ed = clear_val
  int_pers%dp = clear_val
  int_pers%dsubcld = clear_val
  int_pers%jt = clear_val
  int_pers%maxg = clear_val
  int_pers%ideep = clear_val
  int_pers%lengath = 0

end subroutine interstitial_persistent_init

subroutine physics_global_init(pglob)
  implicit none

  class(physics_global)       :: pglob

  real(kind=r8) :: ztodt
  character(len=16) :: deep_conv_scheme

  pglob%cam_physpkg_cam3 = 'cam3'
  pglob%cam_physpkg_cam4 = 'cam4'
  pglob%cam_physpkg_cam5 = 'cam5'
  pglob%cam_physpkg_ideal = 'ideal'
  pglob%cam_physpkg_adiabatic = 'adiabatic'

  if (cam_physpkg_is(pglob%cam_physpkg_cam3)) then
    pglob%cam_physpkg = pglob%cam_physpkg_cam3
  else if (cam_physpkg_is(pglob%cam_physpkg_cam4)) then
    pglob%cam_physpkg = pglob%cam_physpkg_cam4
  else if (cam_physpkg_is(pglob%cam_physpkg_cam5)) then
    pglob%cam_physpkg = pglob%cam_physpkg_cam5
  else if (cam_physpkg_is(pglob%cam_physpkg_ideal)) then
    pglob%cam_physpkg = pglob%cam_physpkg_ideal
  else if (cam_physpkg_is(pglob%cam_physpkg_adiabatic)) then
    pglob%cam_physpkg = pglob%cam_physpkg_adiabatic
  else
    pglob%cam_physpkg = 'UNSET'
  end if

  call phys_getopts(microp_scheme_out=pglob%microp_scheme)

  if (trim(pglob%cam_physpkg) == trim(pglob%cam_physpkg_cam3)) then
    pglob%non_dilute_buoy = .true.
  else
    pglob%non_dilute_buoy = .false.
  end if

  pglob%no_deep_pbl = phys_deepconv_pbl()
  pglob%fv_dycore = dycore_is ('LR')

  call cnst_get_ind('CLDICE', pglob%ixcldice, abort=.false.)
  call cnst_get_ind('CLDLIQ', pglob%ixcldliq, abort=.false.)
  ! Check for number concentration of cloud liquid and cloud ice (if not present
  ! the indices will be set to -1)
  call cnst_get_ind('NUMICE', pglob%ixnumice, abort=.false.)
  call cnst_get_ind('NUMLIQ', pglob%ixnumliq, abort=.false.)

  ztodt = get_step_size()
  pglob%half_ztodt = 0.5_r8*ztodt

  call phys_getopts(deep_scheme_out=deep_conv_scheme)
  pglob%l_windt(:) = .false.
  if ((pglob%cam_physpkg /= pglob%cam_physpkg_cam3) .and. (deep_conv_scheme == 'ZM')) then
    pglob%l_windt(:) = .true.
  end if

end subroutine physics_global_init
#endif
end module physics_types
