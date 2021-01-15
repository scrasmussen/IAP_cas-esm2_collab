
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,27) = 1.20e-10_r8
      rate(:,:,28) = 2.02e-10_r8
      rate(:,:,29) = 1.204e-10_r8
      rate(:,:,30) = 1.31e-10_r8
      rate(:,:,31) = 3.50e-11_r8
      rate(:,:,32) = 9.00e-12_r8
      rate(:,:,35) = 7.20e-11_r8
      rate(:,:,36) = 6.90e-12_r8
      rate(:,:,37) = 1.60e-12_r8
      rate(:,:,41) = 1.80e-12_r8
      rate(:,:,43) = 1.80e-12_r8
      rate(:,:,58) = 1.00e-11_r8
      rate(:,:,59) = 2.20e-11_r8
      rate(:,:,60) = 3.50e-12_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,20) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,22) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      rate(:,:,23) = 3.30e-11_r8 * exp( 55._r8 * itemp(:,:) )
      rate(:,:,24) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,25) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,26) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,34) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,38) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:,39) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,40) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,44) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,89) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,45) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,47) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,49) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,65) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,68) = 7.40e-12_r8 * exp_fac(:,:)
      rate(:,:,50) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:,51) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,53) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,73) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,57) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,62) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,64) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,66) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,67) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      rate(:,:,69) = 6.00e-13_r8 * exp( 230._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,70) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,71) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,74) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,75) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,78) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,81) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,82) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,83) = 3.40e-11_r8 * exp( -1600._r8 * itemp(:,:) )
      rate(:,:,84) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,85) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,86) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,87) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,88) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,90) = 2.700E-11_r8 * exp( 390._r8 * itemp(:,:) )
      rate(:,:,91) = 5.590E-15_r8 * exp( -1814._r8 * itemp(:,:) )
      rate(:,:,96) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,98) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,33), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,42), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,48), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,52), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,54), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,56), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,61), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,72), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,76), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,80), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)


      end subroutine setrxt_hrates

      end module mo_setrxt
