#include <define.h>

module colm_varctl

   use precision
   implicit none

   integer, parameter :: ndst        =   4   ! number of dust size classes (BGC only)
   integer, parameter :: nvoc        =   5   ! number of voc categories

   integer, public, parameter :: nsrStartup  = 0                          ! Startup from initial conditions
   integer, public, parameter :: nsrContinue = 1                          ! Continue from restart files
   integer, public, parameter :: nsrBranch   = 2                          ! Branch from restart files

   integer :: irad                 ! solar radiation frequency (iterations)
   logical :: wrtdia               ! true => write global average diagnostics to std out
   logical :: csm_doflxave         ! true => only communicate with flux coupler on albedo calc time steps
   logical :: doalb                ! albedo 
   integer :: nsrest = nsrStartup  ! nsrStartup, nsrContinue, nsrBranch

   logical, public :: create_glacier_mec_landunit = .false.               ! glacier_mec landunit is not created

   logical :: debug = .false.     ! control debug information output

   logical :: lcycdump  = .false. ! logical flag do enable cycdump function.
   integer :: ncycdump  = 96      ! number of steps to dump
   integer :: lucycdump = 100     ! logical unit of fcycdump
   character(len=255) :: fcycdump = "./CoLM-cycdump"  ! file name of cycdump (restart & forcing)

! history output options
   logical :: lhist_yearly  = .false. ! yearly history files
   logical :: lhist_monthly = .true.  ! monthly history files
   logical :: lhist_daily   = .false. ! daily history files
   logical :: lhist_3hourly = .false. ! 3hourly history files
   logical :: lhist_hourly  = .false. ! hourly history files
   logical :: lhist_steply  = .false. ! steply history files

   character(len=31) :: restart_freq  ! yearly, monthly, daily
   character(len=31) :: history_freq  ! yearly, monthly, daily, 6hourly, 3hourly, hourly

#ifdef RTM
   integer :: rtm_dtime = 3600*3      ! time interval(seconds) to call RTM
#endif

#ifdef SPINUP
   integer :: spinup_numyr  = 500

 !*TH
 ! integer :: forcing_begyr = 1997
 ! integer :: forcing_endyr = 2000

 !*HV
 ! integer :: forcing_begyr = 1992
 ! integer :: forcing_endyr = 1999

 !*RCN
   integer :: forcing_begyr = 2006
   integer :: forcing_endyr = 2006

   integer :: last_nstep    = 0
   integer :: curr_nstep    = 0
#endif

#ifdef COUP_CSM
   integer :: lnd_cflux_year = 0.               ! year to turn on carbon flux exchanging 
!!! Only for CPL6
   character(len=31) :: co2_option = "co2prog"  ! 'co2prog','co2diag'
!!! End for CPL6
#endif

   real(r8) :: po2  = 0.209                     ! constant atmospheric partial pressure  O2 (mol/mol)
   real(r8) :: pco2 = 355.0E-06                 ! constant atmospheric partial pressure CO2 (mol/mol)
!  real(r8) :: pco2 = 284.725E-06               ! constant atmospheric partial pressure CO2 (mol/mol)

#ifdef CPL7
!!! Only for CPL7
   real(r8) :: co2_ppmv = 367.0                 ! atmospheric CO2 molar ratio (by volume) (umol/mol)
   character(len=16) :: co2_type = 'constant'   ! values of 'prognostic','diagnostic','constant'
!!! End for CPL7
#endif

   logical :: luc_emission = .false.

end module colm_varctl
