
module constituents

!----------------------------------------------------------------------------------------------
! 
! Purpose: Contains data and functions for manipulating advected and non-advected constituents.
!
! Revision history:
!             B.A. Boville    Original version
! June 2003   P. Rasch        Add wet/dry m.r. specifier
! 2004-08-28  B. Eaton        Add query function to allow turning off the default CAM output of
!                             constituents so that chemistry module can make the outfld calls.
!                             Allow cnst_get_ind to return without aborting when constituent not
!                             found.
! 2006-10-31  B. Eaton        Remove 'non-advected' constituent functionality.
!----------------------------------------------------------------------------------------------
!> \section arg_table_constituents
!! \htmlinclude constituents.html
!!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst,    only: r_universal
  use spmd_utils,   only: masterproc
  use abortutils,   only: endrun
  use cam_logfile,  only: iulog

  implicit none
  private
  save
!
! Public interfaces
!
  public cnst_add             ! add a constituent to the list of advected constituents
  public cnst_num_avail       ! returns the number of available slots in the constituent array
  public cnst_get_ind         ! get the index of a constituent
  public cnst_get_type_byind  ! get the type of a constituent
  public cnst_get_type_byname ! get the type of a constituent
  public cnst_read_iv         ! query whether constituent initial values are read from initial file
  public cnst_chk_dim         ! check that number of constituents added equals dimensions (pcnst)
  public cnst_cam_outfld      ! Returns true if default CAM output was specified in the cnst_add calls.

! Public data

  integer, parameter, public :: pcnst  = PCNST      ! number of advected constituents (including water vapor)
  character(len=16), public :: cnst_name(pcnst)     ! constituent names
  character(len=128),public :: cnst_longname(pcnst) ! long name of constituents

! Namelist variables
  logical, public :: readtrace = .true.             ! true => obtain initial tracer data from IC file

!
! Constants for each tracer
  real(r8),    public :: cnst_cp  (pcnst)          ! specific heat at constant pressure (J/kg/K)
  real(r8),    public :: cnst_cv  (pcnst)          ! specific heat at constant volume (J/kg/K)
  real(r8),    public :: cnst_mw  (pcnst)          ! molecular weight (kg/kmole)
  character*3, public :: cnst_type(pcnst)          ! wet or dry mixing ratio
  real(r8),    public :: cnst_rgas(pcnst)          ! gas constant ()
  real(r8),    public :: qmin     (pcnst)          ! minimum permitted constituent concentration (kg/kg)
  real(r8),    public :: qmincg   (pcnst)          ! for backward compatibility only
  logical,     public :: cnst_fixed_ubc(pcnst) = .false.  ! upper bndy condition = fixed ?

!++bee - temporary... These names should be declared in the module that makes the addfld and outfld calls.
! Lists of tracer names and diagnostics
   character(len=16), public :: apcnst    (pcnst)   ! constituents after physics  (FV core only)
   character(len=16), public :: bpcnst    (pcnst)   ! constituents before physics (FV core only)
   character(len=16), public :: hadvnam   (pcnst)   ! names of horizontal advection tendencies
   character(len=16), public :: vadvnam   (pcnst)   ! names of vertical advection tendencies
   character(len=16), public :: dcconnam  (pcnst)   ! names of convection tendencies
   character(len=16), public :: fixcnam   (pcnst)   ! names of species slt fixer tendencies
   character(len=16), public :: tendnam   (pcnst)   ! names of total tendencies of species
   character(len=16), public :: ptendnam  (pcnst)   ! names of total physics tendencies of species
   character(len=16), public :: dmetendnam(pcnst)   ! names of dme adjusted tracers (FV)
#ifdef CO2
   character(len=16), public :: sflxnam   (pcnst+4)   ! names of surface fluxes of species, juanxiong he
#else
   character(len=16), public :: sflxnam   (pcnst)   ! names of surface fluxes of species
#endif
   character(len=16), public :: tottnam   (pcnst)   ! names for horz + vert + fixer tendencies

! Private data

  integer :: padv = 0                      ! index pointer to last advected tracer
  logical :: read_init_vals(pcnst)         ! true => read initial values from initial file
  logical :: cam_outfld_(pcnst)            ! true  => default CAM output of constituents in kg/kg
                                           ! false => chemistry is responsible for making outfld
                                           !          calls for constituents

!==============================================================================================
CONTAINS
!==============================================================================================

  subroutine cnst_add (name, mwc, cpc, qminc, &
                       ind, longname, readiv, mixtype, cam_outfld, fixed_ubc)
!----------------------------------------------------------------------- 
! 
! Purpose: Register a constituent to be advected by the large scale winds and transported by
!          subgrid scale processes.
!
!---------------------------------------------------------------------------------
!
    character(len=*), intent(in) :: &
       name      ! constituent name used as variable name in history file output (8 char max)
    real(r8),intent(in)    :: mwc    ! constituent molecular weight (kg/kmol)
    real(r8),intent(in)    :: cpc    ! constituent specific heat at constant pressure (J/kg/K)
    real(r8),intent(in)    :: qminc  ! minimum value of mass mixing ratio (kg/kg)
                                     ! normally 0., except water 1.E-12, for radiation.
    integer, intent(out)   :: ind    ! global constituent index (in q array)

    character(len=*), intent(in), optional :: &
       longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)
    logical,          intent(in), optional :: &
       readiv      ! true => read initial values from initial file (default: true)
    character(len=*), intent(in), optional :: &
       mixtype     ! mixing ratio type (dry, wet)
    logical,          intent(in), optional :: &
       cam_outfld  ! true => default CAM output of constituent in kg/kg
    logical,          intent(in), optional :: &
       fixed_ubc ! true => const has a fixed upper bndy condition

!-----------------------------------------------------------------------

! set tracer index and check validity, advected tracer
    padv = padv+1
    ind  = padv
    if (padv > pcnst) then
       write(iulog,*) 'CNST_ADD: advected tracer index greater than pcnst = ', pcnst
       call endrun
    end if

! set tracer name and constants
    cnst_name(ind) = name
    if ( present(longname) )then
       cnst_longname(ind) = longname
    else
       cnst_longname(ind) = name
    end if

! set whether to read initial values from initial file
    if ( present(readiv) ) then
       read_init_vals(ind) = readiv
    else
       read_init_vals(ind) = readtrace
    end if

! set constituent mixing ratio type
    if ( present(mixtype) )then
       cnst_type(ind) = mixtype
    else
       cnst_type(ind) = 'wet'
    end if

! set outfld type 
! (false: the module declaring the constituent is responsible for outfld calls)
    if ( present(cam_outfld) ) then
       cam_outfld_(ind) = cam_outfld
    else
       cam_outfld_(ind) = .true.
    end if

! set upper boundary condition type
    if ( present(fixed_ubc) ) then
       cnst_fixed_ubc(ind) = fixed_ubc
    else
       cnst_fixed_ubc(ind) = .false.
    end if

    cnst_cp  (ind) = cpc
    cnst_mw  (ind) = mwc
    qmin     (ind) = qminc
    qmincg   (ind) = qminc
    if (ind == 1) qmincg = 0._r8  ! This crap is replicate what was there before ****

    cnst_rgas(ind) = r_universal * mwc
    cnst_cv  (ind) = cpc - cnst_rgas(ind)

    return
  end subroutine cnst_add

!==============================================================================

  function cnst_num_avail()

     ! return number of available slots in the constituent array

     integer cnst_num_avail

     cnst_num_avail = pcnst - padv

  end function cnst_num_avail

!==============================================================================

  subroutine cnst_get_ind (name, ind, abort)
!----------------------------------------------------------------------- 
! 
! Purpose: Get the index of a constituent 
! 
! Author:  B.A. Boville
! 
!-----------------------------Arguments---------------------------------
!
    character(len=*),  intent(in)  :: name  ! constituent name
    integer,           intent(out) :: ind   ! global constituent index (in q array)
    logical, optional, intent(in)  :: abort ! optional flag controlling abort

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    logical :: abort_on_error
!-----------------------------------------------------------------------

! Find tracer name in list
    do m = 1, pcnst
       if (name == cnst_name(m)) then
          ind  = m
          return
       end if
    end do

! Unrecognized name
    abort_on_error = .true.
    if ( present(abort) ) abort_on_error = abort

    if ( abort_on_error ) then
       write(iulog,*) 'CNST_GET_IND, name:', name,  ' not found in list:', cnst_name(:)
       call endrun('CNST_GET_IND: name not found')
    end if

! error return
    ind = -1

  end subroutine cnst_get_ind

!==============================================================================================

  character*3 function cnst_get_type_byind (ind)
!----------------------------------------------------------------------- 
! 
! Purpose: Get the type of a constituent 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  P. J. Rasch
! 
!-----------------------------Arguments---------------------------------
!
    integer, intent(in)   :: ind    ! global constituent index (in q array)

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index

!-----------------------------------------------------------------------

    if (ind.le.pcnst) then
       cnst_get_type_byind = cnst_type(ind)
    else
       ! Unrecognized name
       write(iulog,*) 'CNST_GET_TYPE_BYIND, ind:', ind
       call endrun
    endif


  end function cnst_get_type_byind

!==============================================================================================

  character*3 function cnst_get_type_byname (name)
!----------------------------------------------------------------------- 
! 
! Purpose: Get the type of a constituent 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  P. J. Rasch
! 
!-----------------------------Arguments---------------------------------
!
    character(len=*), intent(in) :: name ! constituent name

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index

!-----------------------------------------------------------------------

    do m = 1, pcnst
       if (name == cnst_name(m)) then
          cnst_get_type_byname = cnst_type(m)
          return
       end if
    end do

! Unrecognized name
    write(iulog,*) 'CNST_GET_TYPE_BYNAME, name:', name,  ' not found in list:', cnst_name(:)
    call endrun

  end function cnst_get_type_byname

!==============================================================================
  function cnst_read_iv(m)
!----------------------------------------------------------------------- 
! 
! Purpose: Query whether constituent initial values are read from initial file.
! 
! Author:  B. Eaton
! 
!-----------------------------Arguments---------------------------------
!
    integer, intent(in) :: m    ! constituent index

    logical :: cnst_read_iv     ! true => read initial values from inital file
!-----------------------------------------------------------------------

    cnst_read_iv = read_init_vals(m)
 end function cnst_read_iv

!==============================================================================
  subroutine cnst_chk_dim
!----------------------------------------------------------------------- 
! 
! Purpose: Check that the number of registered constituents of each type is the
!          same as the dimension
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  B.A. Boville
! 
    integer i,m
!-----------------------------------------------------------------------
!
    if (padv /= pcnst) then
       write(iulog,*)'CNST_CHK_DIM: number of advected tracer ',padv, ' not equal to pcnst = ',pcnst
       call endrun ()
    endif

    if (masterproc) then
       write(iulog,*) 'Advected constituent list:'
       do i = 1, pcnst
          write(iulog,'(i4,2x,a8,2x,a128,2x,a3)') i, cnst_name(i), cnst_longname(i), cnst_type(i)
       end do
    end if

    ! Set names of advected tracer diagnostics
    do m=1,pcnst
       apcnst    (m)  = trim(cnst_name(m))//'AP'
       bpcnst    (m)  = trim(cnst_name(m))//'BP'
       hadvnam   (m)  = 'HA'//cnst_name(m)
       vadvnam   (m)  = 'VA'//cnst_name(m)
       fixcnam   (m)  = 'DF'//cnst_name(m)
       tendnam   (m)  = 'TE'//cnst_name(m)
       ptendnam  (m)  = 'PTE'//cnst_name(m)
       dmetendnam(m)  = 'DME'//cnst_name(m)
       tottnam   (m)  = 'TA'//cnst_name(m)
       sflxnam(m)     = 'SF'//cnst_name(m)
    end do
    !-------------------------------------
    ! juanxiong he
    !-------------------------------------
#ifdef CO2
       sflxnam(pcnst+1)     = 'SFCO2_OCN'
       sflxnam(pcnst+2)     = 'SFCO2_FFF'
       sflxnam(pcnst+3)     = 'SFCO2_LND'
       sflxnam(pcnst+4)     = 'CO2_OUT'
#endif
  end subroutine cnst_chk_dim

!==============================================================================

function cnst_cam_outfld(m)
!----------------------------------------------------------------------- 
! 
! Purpose:
! Query whether default CAM outfld calls should be made.
! 
!----------------------------------------------------------------------- 
   integer, intent(in) :: m                ! constituent index
   logical             :: cnst_cam_outfld  ! true => use default CAM outfld calls
!-----------------------------------------------------------------------

   cnst_cam_outfld = cam_outfld_(m)

end function cnst_cam_outfld

!==============================================================================

end module constituents
