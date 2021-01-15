#include <define.h>

subroutine iniDGVMVar()

   use precision, only : r8
   use paramodel, only : npftpara, numlandc
   use colm_varMod, only : numpatch, pvar
   implicit none

   real(r8) :: pftpara   (npftpara,numlandc) ! PFT parameters
   real(r8) :: vegclass           (numlandc) ! 1.tree 2.shrub 3.grass 4.crop
   real(r8) :: summergreen        (numlandc)
   real(r8) :: raingreen          (numlandc)
   real(r8) :: sla                (numlandc)
   real(r8) :: lm_sapl            (numlandc) ! sapling leafmass
   real(r8) :: sm_sapl            (numlandc) ! sapling sapwood mass
   real(r8) :: hm_sapl            (numlandc) ! sapling heartwood mass
   real(r8) :: rm_sapl            (numlandc) ! sapling rootmass
   real(r8) :: stemdiam           (numlandc) ! canopy height initialization
#ifdef DyN
   real(r8) :: cton_soil          (numlandc) ! soil C:N mass ratio
   real(r8) :: cton_pro           (numlandc) ! C:N mass ratio in production
#endif

   integer p, ivt

   CALL setDGVMpara(pftpara, vegclass, summergreen, raingreen &
#ifdef DyN
                   ,cton_pro, cton_soil &
#endif
                   ,sla, lm_sapl, sm_sapl, hm_sapl, rm_sapl, stemdiam)

   do p = 1, numpatch
      ivt = pvar%ivt(p)

      pvar%pftpara  (:,p) = pftpara  (:,ivt)
      pvar%vegclass   (p) = vegclass   (ivt)
      pvar%summergreen(p) = summergreen(ivt)
      pvar%raingreen  (p) = raingreen  (ivt)
      pvar%sla        (p) = sla        (ivt)
      pvar%lm_sapl    (p) = lm_sapl    (ivt)
      pvar%sm_sapl    (p) = sm_sapl    (ivt)
      pvar%hm_sapl    (p) = hm_sapl    (ivt)
      pvar%rm_sapl    (p) = rm_sapl    (ivt)
      pvar%stemdiam   (p) = stemdiam   (ivt)

#ifdef DyN
      pvar%cton_soil  (p) = cton_soil  (ivt)
      pvar%cton_pro   (p) = cton_pro   (ivt)
#endif  
   end do

end subroutine iniDGVMVar

subroutine setDGVMpara(pftpara,vegclass,summergreen,raingreen&
#ifdef DyN
                      ,cton_pro,cton_soil&
#endif
                      ,sla,lm_sapl,sm_sapl,hm_sapl,rm_sapl,stemdiam)

!--------------------------------------------------------------------------
!DESCRIPTION:
! If defined DGVM, initialize parameters (SLA, leafmass, rootmass,
! sapwood and heartwood mass, sterm diameter) for each PFT.
! (The vegetation parameters are the same as NCAR-DGVM)
!CALLED FROM:initialize.F90
!---------------------------------------------------------------------------
  use precision
  use nanMod, only: inf        !Set parameters for the floating point flags "inf" Infinity
  use paramodel, only : npftpara, numlandc
  implicit none

!--------------------output
  real(r8),INTENT(out)::      &!
     pftpara (npftpara,numlandc),    &! PFT parameters
     sla              (numlandc),    &! Specialized leaf area
     lm_sapl          (numlandc),    &! sapling leafmass
     sm_sapl          (numlandc),    &! sapling sapwood mass
     hm_sapl          (numlandc),    &! sapling heartwood mass
     rm_sapl          (numlandc),    &! sapling rootmass
     vegclass         (numlandc),    &! 1.tree 2.shrub 3.grass 4.crop
     summergreen      (numlandc),    &! 1 = summergreen, -1 = not summergreen
     raingreen        (numlandc),    &! 1 = raingreen, -1 = not raingreen
     stemdiam         (numlandc)      ! stemdiam initialization
#ifdef DyN
  real(r8),INTENT(out)::      &!
     cton_soil        (numlandc),    &! soil C:N mass ratio
     cton_pro         (numlandc)      ! C:N mass ratio in production
#endif
!-------------local variables:
  real(r8), parameter :: PI = 3.14159265358979323846  ! pi
  real(r8), parameter :: reinickerp = 1.6 !parameter in allometric equation
  real(r8), parameter :: wooddens = 2.0e5 !wood density (gC/m3)

! adopted from X.D.Z's shrub submodule by zhq @ 09/14/2010
!  real(r8), parameter :: latosa = 8.0e3   !ratio of leaf area to sapwood cross-sectional
                                          !area (Shinozaki et al 1964a,b)
!  real(r8), parameter :: allom1 = 100.0   !parameters in allometric
!  real(r8), parameter :: allom2 =  40.0
!  real(r8), parameter :: allom3 =   0.5

  real(r8) :: latosa(numlandc)       ! leafarea:sapwood cross-sectional area
  real(r8) :: allom1(numlandc)       ! parameters in allometric
  real(r8) :: allom2(numlandc)       ! parameters in allometric
  real(r8) :: allom3(numlandc)       ! parameters in allometric

  real(r8):: table(npftpara+2,numlandc)
  real(r8):: height_sapl(numlandc)
  integer :: nlandc, k         ! indices
  real(r8):: x

!---------------------------------------------------------------------
    ! Paramters for shrub submodule. zhq @ 09/14/2010

    latosa(1:21) = 8.0e3   ! leafarea:sapwood cross-sectional area
    allom1(1:21) = 100.0   ! parameters in allometric
    allom2(1:21) =  40.0   ! parameters in allometric
    allom3(1:21) =   0.5   ! parameters in allometric

    latosa(9:11) = 4000.0  ! 2000.0 in last test
    allom1(9:11) = 200.0
    allom2(9:11) = 10.0

!-----------------------------------------------------------------------
#ifndef IAPDGVM
    ! PFT PARAMETERS (follows LPJ & NCAR DGVM pft parameters)
    ! zhq,05/13/2010: Nitrogen PARAMETERS from DyN, Xu-Ri. et al.,2008

    !  1  fraction of roots in upper soil layer
    !  2  plants with C4 (1) or C3 (0) photosynthetic pathway
    !  3  water scalar value at which leaves shed by drought deciduous PFT
    !  4  canopy conductance component (gmin, mm/s) not associated with
    !     photosynthesis (Haxeltine & Prentice 1996, Table 4)
    !  5  maintenance respiration coefficient
    !  6  flammability threshold
    !  7  maximum foliar N content (mg/g)
    !     (Haxeltine & Prentice 1996a, Fig 4)
    !  8  fire resistance index
    !  9  leaf turnover period (years)
    ! 10  leaf longevity (years)
    ! 11  sapwood turnover period (sapwood converted to heartwood) (years)
    ! 12  root turnover period (years)
    ! 13  leaf C:N mass ratio
    ! 14  sapwood C:N mass ratio
    ! 15  root C:N mass ratio
    ! 16  leaf type: broadleaved (1), needleleaved (2) or grass (3)
    ! 17  phenology type: evergreen (1), summergreen (2), raingreen (3),
    !     any type (4)
    ! 18  leaf to root ratio under non-water stressed conditions
    ! 19  summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
    ! 20  tree maximum crown area (m2)
    ! 21  sapling (or grass on initialisation) LAI
    ! 22  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    ! 23  boreal pft (1), non-boreal pft (0)
    ! 24  low temperature limit for CO2 uptake
    ! 25  lower range of temperature optimum for photosynthesis
    ! 26  upper range of temperature optimum for photosynthesis
    ! 27  high temperature limit for CO2 unptake

    !BIOCLIMATIC LIMITS

    ! 28 minimum coldest monthly mean temperature
    ! 29 maximum coldest monthly mean temperature
    ! 30 minimum growing degree days (at or above 5 deg C)
    ! 31 upper limit of temperature of the warmest month (twmax)
    ! 32 lower limit of growth efficiency (g/m2)

    ! C:N CONSTANT PARAMETERS (from DyN, Xu-Ri. et al.,2008)

    ! 33 C:N ratios for plant production
    ! 34 C:N ratios for soil

    ! ------------------------------------------------------------------------
    !      1      2      3      4      5      6      7      8       landcover 
    ! ------------------------------------------------------------------------
    data ((table(k,nlandc),k=1,8),nlandc=1,21) /                &
         0.70,   0.0,  0.00,   0.3,  1.20,  0.15, 100.0,  0.12, &      ! P1
         0.90,   0.0,  0.00,   0.3,  1.20,  0.15, 100.0,  0.12, &      ! P2
         0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      ! P3
         0.85,   0.0,  0.00,   0.5,  0.20,  0.15, 100.0,  0.12, &      ! P4
         0.70,   0.0,  0.00,   0.5,  1.20,  0.15, 100.0,  0.50, &      ! P5
         0.70,   0.0,  0.35,   0.5,  0.20,  0.15, 100.0,  0.50, &      ! P6
         0.80,   0.0,  0.00,   0.5,  1.20,  0.15, 120.0,  0.12, &      ! P7
         0.90,   0.0,  0.00,   0.3,  1.20,  0.15, 100.0,  0.12, &      ! P8
         0.85,   0.0,  0.00,   0.5,  1.20,  0.15, 100.0,  0.12, &      ! P9
         0.80,   0.0,  0.00,   0.5,  1.20,  0.15, 120.0,  0.12, &      ! P10
         0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      ! P11
         0.90,   0.0,  0.35,   0.5,  1.00,  0.10, 100.0,  1.00, &      ! P12
         0.90,   0.0,  0.35,   0.5,  1.00,  0.10, 100.0,  1.00, &      ! P13
         0.90,   1.0,  0.35,   0.5,  1.00,  0.10, 100.0,  1.00, &      ! P14
         0.90,   1.0,  0.35,   0.5,  1.20,  0.15, 100.0,  1.00, &      ! P15
         0.90,   1.0,  0.35,   0.5,  1.20,  0.15, 100.0,  1.00, &      ! P16
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      ! 17
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      ! 18
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      ! 19
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      ! 20
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf/        ! 21


    ! ---------------------------------------------------------------------
    !     9     10     11     12     13     14     15     16     17  landcover 
    ! ---------------------------------------------------------------------
    data ((table(k,nlandc),k=9,17),nlandc=1,21) /               &
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, & ! P1
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, & ! P2
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! P3
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! P4
         1.0,  1.00,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! P5
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, & ! P6
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, & ! P7
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! P8
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! p9
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, & ! P10 (17) 2.0 -> 3.0 by X.D.Z
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! P11
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P12
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P13
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P14
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P15
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P16
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 17
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 18
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 19
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 20
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf/   ! 21

    ! ------------------------------------------------------
    !    18      19     20      21    22     23    landcover 
    ! ------------------------------------------------------
    data ((table(k,nlandc),k=18,23),nlandc=1,21) /               &
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P1
         1.0, 1000.0,  15.0,  1.500,  1.2,   1.0, &  ! P2
         1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! P3
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P4
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P5
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P6
         1.0,  200.0,  15.0,  1.500,  1.2,   0.0, &  ! P7
         1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! P8
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P9
         1.0,  200.0,   5.0,  1.500,  1.2,   0.0, &  ! P10
         1.0,  200.0,   5.0,  1.500,  1.2,   1.0, &  ! P11
        0.75,  100.0,   0.0,  0.001,  1.2,   1.0, &  ! P12
        0.75,  100.0,   0.0,  0.001,  1.2,   1.0, &  ! P13
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  ! P14
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  ! P15
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  ! P16
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 17
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 18
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 19
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 20
         inf,    inf,   inf,    inf,  inf,   inf/    ! 21

    ! -------------------------------------
    !      24     25     26      27   landcover 
    ! -------------------------------------
    data ((table(k,nlandc),k=24,27),nlandc=1,21) /               &
         -4.0,  20.0,  30.0,   42.0, & ! P1
         -4.0,  15.0,  25.0,   38.0, & ! P2
         -4.0,  15.0,  25.0,   38.0, & ! P3
          2.0,  25.0,  30.0,   55.0, & ! p4
         -4.0,  20.0,  30.0,   42.0, & ! P5
          2.0,  25.0,  30.0,   55.0, & ! P6
         -4.0,  20.0,  25.0,   38.0, & ! P7
         -4.0,  15.0,  25.0,   38.0, & ! P8
          2.0,  25.0,  30.0,   55.0, & ! P9
         -4.0,  20.0,  25.0,   38.0, & ! P10
         -4.0,  15.0,  25.0,   38.0, & ! P11
         -4.0,  10.0,  30.0,   45.0, & ! P12
         -4.0,  10.0,  30.0,   45.0, & ! P13
          6.0,  20.0,  45.0,   55.0, & ! P14
          6.0,  20.0,  45.0,   55.0, & ! P15
          6.0,  20.0,  45.0,   55.0, & ! P16
          inf,   inf,   inf,    inf, & ! 17
          inf,   inf,   inf,    inf, & ! 18
          inf,   inf,   inf,    inf, & ! 19
          inf,   inf,   inf,    inf, & ! 20
          inf,   inf,   inf,    inf/   ! 21


    ! --------------------------------------------------------
    !      28       29      30       31      32     33    34    landcover
    ! --------------------------------------------------------
    data ((table(k,nlandc),k=28,34),nlandc=1,21) /                  &
         -2.0,    22.0,  900.0,  1000.0,    0.0,  89.17,  23.86, & !  P1
        -32.5,    -2.0,  600.0,    23.0,    0.0,  52.38,  29.70, & !  P2
       9999.9,    -2.0,  350.0,    23.0,    0.0,  45.24,  18.15, & !  P3
         15.5,  1000.0,    0.0,  1000.0,    0.0,  43.75,  16.73, & !  P4
          3.0,    18.8, 1200.0,  1000.0,    0.0,  90.63,  25.78, & !  P5
         15.5,  1000.0,    0.0,  1000.0,    0.0,  32.66,   8.31, & !  P6
        -17.0,    15.5, 1200.0,  1000.0,    0.0,  65.00,  20.09, & !  P7
      -1000.0,    -2.0,  350.0,    23.0,    0.0,  45.24,  18.15, & !  P8
       9999.9,  1000.0,    0.0,  1000.0,    0.0,  29.00,  19.00, & !  P9  (was 15.5)
        -17.0,    15.5, 1200.0,  1000.0,    0.0,  29.00,  19.00, & !  P10 (was -17.0)
      -1000.0,    -2.0,  350.0,    23.0,    0.0,  29.00,  19.00, & !  P11 (was -1000.0)
      -1000.0,   -17.0,    0.0,  1000.0,    0.0,  54.29,   9.77, & !  P12
        -17.0,    15.5,    0.0,  1000.0,    0.0,  54.29,   9.77, & !  P13
         15.5,  1000.0,    0.0,  1000.0,    0.0,  69.55,  10.34, & !  P14
       9999.9,  1000.0,    0.0,  1000.0,    0.0,  29.00,  19.00, & !  P15
       9999.9,  1000.0,    0.0,  1000.0,    0.0,  29.00,  19.00, & !  P16
          inf,     inf,    inf,  1000.0,    inf,   inf ,   inf , & !  17
          inf,     inf,    inf,  1000.0,    inf,   inf ,   inf , & !  18
          inf,     inf,    inf,  1000.0,    inf,   inf ,   inf , & !  19
          inf,     inf,    inf,  1000.0,    inf,   inf ,   inf , & !  20
          inf,     inf,    inf,  1000.0,    inf,   inf ,   inf/    !  21

#else
!----------------------------------------------------------------------
!parameters for IAPDGVM
!---------------------------------------------------------------------
    ! PFT PARAMETERS (follows LPJ & NCAR DGVM pft parameters)
    ! zhq,05/13/2010: Nitrogen PARAMETERS from DyN, Xu-Ri. et al.,2008
    ! changed by zhujw, according to the parameters in IAPDGVM  20181122

    !  1  fraction of roots in upper soil layer
    !  2  plants with C4 (1) or C3 (0) photosynthetic pathway
    !  3  water scalar value at which leaves shed by drought deciduous PFT
    !  4  canopy conductance component (gmin, mm/s) not associated with
    !     photosynthesis (Haxeltine & Prentice 1996, Table 4)
    !  5  maintenance respiration coefficient
    !  6  flammability threshold
    !  7  maximum foliar N content (mg/g)
    !     (Haxeltine & Prentice 1996a, Fig 4)
    !  8  fire resistance index
    !  9  leaf turnover period (years)
    ! 10  leaf longevity (years)
    ! 11  sapwood turnover period (sapwood converted to heartwood) (years)
    ! 12  root turnover period (years)
    ! 13  leaf C:N mass ratio
    ! 14  sapwood C:N mass ratio
    ! 15  root C:N mass ratio
    ! 16  leaf type: broadleaved (1), needleleaved (2) or grass (3)
    ! 17  phenology type: evergreen (1), summergreen (2), raingreen (3),any type (4)
    ! 18  leaf to root ratio under non-water stressed conditions
    ! 19  summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
    ! 20  tree maximum crown area (m2)
    ! 21  sapling (or grass on initialisation) LAI
    ! 22  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    ! 23  boreal pft (1), non-boreal pft (0)
    ! 24  low temperature limit for CO2 uptake
    ! 25  lower range of temperature optimum for photosynthesis
    ! 26  upper range of temperature optimum for photosynthesis
    ! 27  high temperature limit for CO2 unptake

      !BIOCLIMATIC LIMITS

    ! 28 minimum coldest monthly mean temperature
    ! 29 maximum coldest monthly mean temperature
    ! 30 minimum growing degree days (at or above 5 deg C)
    ! 31 upper limit of temperature of the warmest month (twmax)
    ! 32 lower limit of growth efficiency (g/m2)

!-------------------------------------------------------------------
!add parameters(33-49) for Fire in IAPDGVM
!-------------------------------------------------------------------
! FIRE INDEX 
    ! 33 combustion completeness index for leaf
    ! 34 combustion completeness index for stem
    ! 35 combustion completeness index for root
    ! 36 combustion completeness index for above-ground litter

    ! Trace gased from burning carbon 
    ! 37 EF for CO2
    ! 38 EF for CO
    ! 39 EF for CH4
    ! 40 EF for NHMC
    ! 41 EF for H2
    ! 42 EF for NOx
    ! 43 EF for N2O

    ! Aerosol from burning carbon 
    ! 44 EF for PM25
    ! 45 EF for TPM
    ! 46 EF for TC
    ! 47 EF for OC
    ! 48 EF for BC

    ! 49 umax  m/s 
 !----------------------------------------------------------------------------
 !end of the added fire parameters for IAPDGVM
 !----------------------------------------------------------------------------
    ! C:N CONSTANT PARAMETERS (from DyN, Xu-Ri. et al.,2008)

    ! 50 C:N ratios for plant production
    ! 51 C:N ratios for soil
    ! -- ---------------------------------------------------------------------
    !      1      2      3      4      5      6      7      8       landcover 
    ! ------------------------------------------------------------------------
    data ((table(k,nlandc),k=1,8),nlandc=1,21) /               &
         0.70,   0.0,  0.00,   0.3,  0.65,  0.69, 100.0,  0.10, &      ! P1
         0.90,   0.0,  0.00,   0.3,  0.57,  0.69, 100.0,  0.12, &      ! P2
         0.90,   0.0,  0.00,   0.3,  0.60,  0.69, 100.0,  0.12, &      ! P3
         0.85,   0.0,  0.00,   0.5,  0.37,  0.69, 100.0,  0.10, &      ! P4
         0.70,   0.0,  0.00,   0.5,  0.70,  0.69, 100.0,  0.10, &      ! P5
         0.70,   0.0,  0.35,   0.5,  0.52,  0.69, 100.0,  0.07, &      ! P6
         0.80,   0.0,  0.00,   0.5,  0.92,  0.69, 120.0,  0.07, &      ! P7
         0.90,   0.0,  0.00,   0.3,  0.61,  0.69, 100.0,  0.10, &      ! P8
         0.85,   0.0,  0.00,   0.5,  1.20,  0.69, 100.0,  0.15, &      ! P9
         0.80,   0.0,  0.00,   0.5,  1.20,  0.69, 120.0,  0.15, &      ! P10
         0.90,   0.0,  0.00,   0.3,  0.38,  0.69, 100.0,  0.15, &      ! P11
         0.90,   0.0,  0.35,   0.5,  0.75,  0.69, 100.0,  0.20, &      ! P12
         0.90,   0.0,  0.35,   0.5,  0.28,  0.69, 100.0,  0.20, &      ! P13
         0.90,   1.0,  0.35,   0.5,  1.20,  0.69, 100.0,  0.20, &      ! P14
         0.90,   1.0,  0.35,   0.5,  1.20,  0.69, 100.0,  0.20, &      ! P15
         0.90,   1.0,  0.35,   0.5,  1.20,  0.69, 100.0,  0.20, &      ! P16         
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      ! 17
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      ! 18
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      ! 19
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      ! 20
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf/        ! 21
    ! ---------------------------------------------------------------------
    !     9     10     11     12     13     14     15     16     17  landcover 
    ! ---------------------------------------------------------------------
    data ((table(k,nlandc),k=9,17),nlandc=1,21) /               &
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, & ! P1
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, & ! P2
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! P3
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! P4
         1.0,  1.00,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! P5
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, & ! P6
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, & ! P7
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! P8
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! p9
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, & ! P10 (17) 2.0 ->3.0 by X.D.Z
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! P11
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P12
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P13
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P14
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P15
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & ! P16
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 17
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 18
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 19
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 20
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf/   ! 21
    ! ------------------------------------------------------
    !    18      19     20      21    22     23    landcover 
    ! ------------------------------------------------------
    data ((table(k,nlandc),k=18,23),nlandc=1,21) /               &
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P1
         1.0, 1000.0,  15.0,  1.500,  1.2,   1.0, &  ! P2
         1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! P3
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P4
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P5
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P6
         1.0,  200.0,  15.0,  1.500,  1.2,   0.0, &  ! P7
         1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! P8
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! P9
         1.0,  200.0,  15.0,  1.500,  1.2,   0.0, &  ! P10
         1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! P11
        0.75,  100.0,   0.0,  0.001,  1.2,   1.0, &  ! P12
        0.75,  100.0,   0.0,  0.001,  1.2,   1.0, &  ! P13
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  ! P14
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  ! P15
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  ! P16
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 17
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 18
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 19
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 20
         inf,    inf,   inf,    inf,  inf,   inf/    ! 21
! -------------------------------------
    !      24     25     26      27   landcover 
    ! -------------------------------------
    data ((table(k,nlandc),k=24,27),nlandc=1,21) / &
         -4.0,  20.0,  30.0,   42.0, & ! P1
         -4.0,  15.0,  25.0,   38.0, & ! P2
         -4.0,  15.0,  25.0,   38.0, & ! P3
          2.0,  25.0,  30.0,   55.0, & ! p4
         -4.0,  20.0,  30.0,   42.0, & ! P5
          2.0,  25.0,  30.0,   55.0, & ! P6
         -4.0,  20.0,  25.0,   38.0, & ! P7
         -4.0,  15.0,  25.0,   38.0, & ! P8
          2.0,  25.0,  30.0,   55.0, & ! P9
         -4.0,  20.0,  25.0,   38.0, & ! P10
         -4.0,  15.0,  25.0,   38.0, & ! P11
         -4.0,  10.0,  30.0,   45.0, & ! P12
         -4.0,  10.0,  30.0,   45.0, & ! P13
          6.0,  20.0,  45.0,   55.0, & ! P14
          6.0,  20.0,  45.0,   55.0, & ! P15
          6.0,  20.0,  45.0,   55.0, & ! P16
          inf,   inf,   inf,    inf, & ! 17
          inf,   inf,   inf,    inf, & ! 18
          inf,   inf,   inf,    inf, & ! 19
          inf,   inf,   inf,    inf, & ! 20
          inf,   inf,   inf,    inf/   ! 21
    ! --------------------------------------------------------
    !      28       29      30       31      32    landcover
    ! --------------------------------------------------------
    data ((table(k,nlandc),k=28,32),nlandc=1,21) / &
         -2.0,    22.0,  900.0,  1000.0,    0.0, & !  P1
        -32.5,    -2.0,  600.0,    23.0,    0.0, & !  P2
       9999.9,    -2.0,  350.0,    23.0,    0.0, & !  P3
         15.5,  1000.0,    0.0,  1000.0,    0.0, & !  P4
          3.0,    18.8, 1200.0,  1000.0,    0.0, & !  P5
         15.5,  1000.0,    0.0,  1000.0,    0.0, & !  P6
        -17.0,    15.5, 1200.0,  1000.0,    0.0, & !  P7
      -1000.0,    -2.0,  350.0,    23.0,    0.0, & !  P8 -1000
       9999.9,  1000.0,    0.0,  1000.0,    0.0, & !  P9  (was  15.5)
        -17.0,    15.5, 1200.0,  1000.0,    0.0, & !  P10 (was -17.0)
      -1000.0,    -2.0,  350.0,    23.0,    0.0, & !  P11 (was -1000.0)
      -1000.0,   -17.0,    0.0,  1000.0,    0.0, & !  P12 (was -1000.0)
        -17.0,    15.5,    0.0,  1000.0,    0.0, & !  P13
         15.5,  1000.0,    0.0,  1000.0,    0.0, & !  P14
       9999.9,  1000.0,    0.0,  1000.0,    0.0, & !  P15
       9999.9,  1000.0,    0.0,  1000.0,    0.0, & !  P16
          inf,     inf,    inf,  1000.0,    inf, & !  17
          inf,     inf,    inf,  1000.0,    inf, & !  18
          inf,     inf,    inf,  1000.0,    inf, & !  19
          inf,     inf,    inf,  1000.0,    inf, & !  20
          inf,     inf,    inf,  1000.0,    inf/    !  21
    ! --------------------------------------------------------
    !      33       34      35       36       
    ! --------------------------------------------------------
    data ((table(k,nlandc),k=33,36),nlandc=1,21) / &
         0.70,    0.15,    0.00,    0.50, & !  P1
         0.75,    0.20,    0.00,    0.55, & !  P2
         0.75,    0.20,    0.00,    0.55, & !  P3 
         0.70,    0.15,    0.00,    0.50, & !  P4
         0.70,    0.15,    0.00,    0.50, & !  P5
         0.70,    0.10,    0.00,    0.45, & !  P6
         0.70,    0.10,    0.00,    0.45, & !  P7
         0.70,    0.15,    0.00,    0.50, & !  P8
         0.80,    0.30,    0.00,    0.60, & !  P9 
         0.80,    0.30,    0.00,    0.60, & ! P10 
         0.80,    0.30,    0.00,    0.60, & ! P11 
         0.85,    0.00,    0.00,    0.85, & ! P12 
         0.85,    0.00,    0.00,    0.85, & ! P13 
         0.85,    0.00,    0.00,    0.85, & ! P14
         0.85,    0.00,    0.00,    0.85, & ! P15
         0.85,    0.00,    0.00,    0.85, & ! P16 
          inf,    inf,      inf,     inf, & !  17
          inf,    inf,      inf,     inf, & !  18
          inf,    inf,      inf,     inf, & !  19
          inf,    inf,      inf,     inf, & !  20
          inf,    inf,      inf,     inf/  !  21
   ! --------------------------------------------------------
   !      37    38     39    40     41    42     43   
   ! --------------------------------------------------------
    data ((table(k,nlandc),k=37,43),nlandc=1,21) / &
     1568.0,  106.0,  4.8,  5.7,  1.81,  3.00,  0.26, & !  P1
     1568.0,  106.0,  4.8,  5.7,  1.81,  3.00,  0.26, & !  P2
     1568.0,  106.0,  4.8,  5.7,  1.81,  3.00,  0.26, & !  P3 
     1580.0,  102.0,  6.8,  8.1,  3.80,  1.85,  0.20, & !  P4
     1568.0,  106.0,  4.8,  5.7,  1.81,  3.00,  0.26, & !  P5
     1568.0,  106.0,  4.8,  5.7,  1.81,  3.00,  0.26, & !  P6
     1568.0,  106.0,  4.8,  5.7,  1.81,  3.00,  0.26, & !  P7
     1568.0,  106.0,  4.8,  5.7,  1.81,  3.00,  0.26, & !  P8
     1664.0,   63.0,  2.2,  3.4,  0.99,  2.35,  0.21, & !  P9 
     1664.0,   63.0,  2.2,  3.4,  0.99,  2.35,  0.21, & ! P10 
     1664.0,   63.0,  2.2,  3.4,  0.99,  2.35,  0.21, & ! P11 
     1664.0,   63.0,  2.2,  3.4,  0.99,  2.35,  0.21, & ! P12 
     1664.0,   63.0,  2.2,  3.4,  0.99,  2.35,  0.21, & ! P13 
     1664.0,   63.0,  2.2,  3.4,  0.99,  2.35,  0.21, & ! P14
     1664.0,   63.0,  2.2,  3.4,  0.99,  2.35,  0.21, & ! P15 
     1664.0,   63.0,  2.2,  3.4,  0.99,  2.35,  0.21, & ! P16 
        inf,    inf,  inf,  inf,   inf,   inf,   inf, & ! 17
        inf,    inf,  inf,  inf,   inf,   inf,   inf, & ! 18
        inf,    inf,  inf,  inf,   inf,   inf,   inf, & ! 19
        inf,    inf,  inf,  inf,   inf,   inf,   inf, & ! 20
        inf,    inf,  inf,  inf,   inf,   inf,   inf/   ! 21
    ! --------------------------------------------------------
    !      44      45     46     47    48     49  
    ! --------------------------------------------------------
    data ((table(k,nlandc),k=44,49),nlandc=1,21) / &
         13.0,   17.6,   8.3,   9.1,   0.56,  0.13, & !  1
         13.0,   17.6,   8.3,   9.1,   0.56,  0.15, & !  2
         13.0,   17.6,   8.3,   9.1,   0.56,  0.15, & !  3 
          9.1,    8.5,   6.0,   5.2,   0.63,  0.13, & !  4
         13.0,   17.6,   8.3,   9.1,   0.56,  0.13, & !  5
         13.0,   17.6,   8.3,   9.1,   0.56,  0.11, & !  6
         13.0,   17.6,   8.3,   9.1,   0.56,  0.11, & !  7
         13.0,   17.6,   8.3,   9.1,   0.56,  0.13, & !  8
          4.9,    8.5,   3.7,   3.2,   0.46,  0.17, & !  9 
          4.9,    8.5,   3.7,   3.2,   0.46,  0.17, & ! 10 
          4.9,    8.5,   3.7,   3.2,   0.46,  0.17, & ! 11 
          4.9,    8.5,   3.7,   3.2,   0.46,  0.20, & ! 12 
          4.9,    8.5,   3.7,   3.2,   0.46,  0.20, & ! 13 
          4.9,    8.5,   3.7,   3.2,   0.46,  0.20, & ! 14
          4.9,    8.5,   3.7,   3.2,   0.46,  0.20, & ! 15
          4.9,    8.5,   3.7,   3.2,   0.46,  0.20, &  ! 16 
          inf,    inf,   inf,   inf,    inf,  0.00, & !  17
          inf,    inf,   inf,   inf,    inf,  0.00, & !  18
          inf,    inf,   inf,   inf,    inf,  0.00, & !  19
          inf,    inf,   inf,   inf,    inf,  0.00, & !  20
          inf,    inf,   inf,   inf,    inf,  0.00/   !  21
 ! --------------------------------------------------------
    !      50    51    landcover
    ! --------------------------------------------------------
    data ((table(k,nlandc),k=50,51),nlandc=1,21) /                  &
         89.17,  23.86, & !  P1
         52.38,  29.70, & !  P2
         45.24,  18.15, & !  P3
         43.75,  16.73, & !  P4
         90.63,  25.78, & !  P5
         32.66,   8.31, & !  P6
         65.00,  20.09, & !  P7
         45.24,  18.15, & !  P8
         29.00,  19.00, & !  P9 
         29.00,  19.00, & !  P10 
         29.00,  19.00, & !  P11
         54.29,   9.77, & !  P12
         54.29,   9.77, & !  P13
         69.55,  10.34, & !  P14
         29.00,  19.00, & !  P15
         29.00,  19.00, & !  P16
          inf ,   inf , & !  17
          inf ,   inf , & !  18
          inf ,   inf , & !  19
          inf ,   inf , & !  20
          inf ,   inf/    !  21
!--------------------------------------
#endif
! end of IAPDGVM
!----------------------------------------------------------------------
!other paramters of landcover
!---------------------------------------------------------------------

!1.For non-vegetated land cover:

    sla        (17:21) = inf  !set infinite value
    lm_sapl    (17:21) = inf
    rm_sapl    (17:21) = inf
    sm_sapl    (17:21) = inf
    hm_sapl    (17:21) = inf
    vegclass   (17:21) = -1
    summergreen(17:21) = -1
    raingreen  (17:21) = -1
    stemdiam   (17:21) = 0.

!2.For vegetated land cover:

    pftpara(1:npftpara,:) = table(1:npftpara,:)
#ifdef DyN
    cton_pro(:)     = table(npftpara+1,:)
    cton_soil(:)    = table(npftpara+2,:)
#endif

      do nlandc = 1,16

       ! Assign leaf and phenology logicals

       if (pftpara(16,nlandc) <= 2.0) then  !woody vegetation: trees, shrubs
          if (nlandc <= 8)then
             vegclass(nlandc) = 1.         !tree
          else
             vegclass(nlandc) = 2.         !shrub
          endif
       else                                !non woody vegetation: grasses
          if (nlandc <= 14)then
             vegclass(nlandc) = 3.         !grass
          else if (nlandc <=16)then
             vegclass(nlandc) = 4.         !crop
          endif
       endif

       if     (pftpara(17,nlandc) == 1.0) then !evergreen
          summergreen(nlandc) = -1
          raingreen(nlandc)   = -1
       else if(pftpara(17,nlandc) == 2.0) then !summergreen
          summergreen(nlandc) = 1
          raingreen(nlandc)   = -1
       else if(pftpara(17,nlandc) == 3.0) then !raingreen
          summergreen(nlandc) = -1
          raingreen(nlandc)   = 1
       else                                !any of the above
          summergreen(nlandc) = 1
          raingreen(nlandc)   = 1
       end if

       ! Calculate specific leaf area (SLA) for each PFT from leaf longevity
       ! Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
       ! Equation based on Reich et al 1997:

       ! SLA = 2e-4 * exp(6.15 - 0.46 ln (leaf_longevity * 12))

       ! SLA in m2/gC, leaf_longevity in years

       sla(nlandc) = 2.0e-4 * exp(6.15 - 0.46*log(pftpara(10,nlandc)*12.0))
       ! Define initial mass structure

       if (vegclass(nlandc) .le. 2.) then !woody PFTs

          ! Calculate leafmass for a sapling individual
          !  (1) lai = leafmass * sla / (crown area)
          !  (2) (leaf area) = latosa * (sapwood xs area)
          !         (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
          !  (3) (crown area) = allom1 * (stem diameter) ** reinickerp
          !         (Reinickes theory)
          ! From (1),
          !  (4) leafmass = lai * (crown area) / sla
          ! From (1) & (3),
          !  (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
          ! From (2),
          !  (6) leafmass = latosa * (sapwood xs area) / sla
          !  (7) (sapwood xs area) = PI * (sapwood diameter)**2 / 4
          ! From (6) and (7),
          !  (8) leafmass = latosa * PI * (sapwood diameter)**2 / 4 / sla
          ! From (8),
          !  (9) (sapwood diameter) = [ 4 * leafmass * sla / PI / latosa ]**0.5
          ! (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
          ! Define x,
          ! (11) x = [ (sapwood diameter)+(heartwood diameter) ] /
          !          (sapwood diameter)
          ! From (10) & (11),
          ! (12) (stem diameter) = x * (sapwood diameter)
          ! From (5), (9) & (12),
          ! (13) leafmass = lai * allom1 * x**reinickerp *
          !               (4*leafmass*sla/PI/latosa)**(reinickerp*0.5) / sla
          ! From (13),
          ! (14) leafmass = [ lai * allom1 * x**reinickerp *
          !      (4*sla/PI/latosa)**(reinickerp*0.5) / sla ]**(2/(2-reinickerp))

          x = pftpara(22,nlandc)

          lm_sapl(nlandc) = (pftpara(21,nlandc) * allom1(nlandc) * x**reinickerp *          &
               (4.0 * sla(nlandc) / PI / latosa(nlandc))**(reinickerp * 0.5) / sla(nlandc))** &
               (2.0/(2.0-reinickerp)) !eqn 14
          ! Calculate sapling stem diameter
          ! From (9) & (12),
          ! (15) (stem diameter) = x * [ 4 * leafmass * sla / PI / latosa ]**0.5

          stemdiam(nlandc) = x * (4.0*lm_sapl(nlandc)*sla(nlandc)/PI/latosa(nlandc))**0.5 !Eqn 15

          ! Calculate sapling height
          ! (16) height = allom2 * (stem diameter)**allom3 (source: Sitch,2003 Eqn 3)

          height_sapl(nlandc) = allom2(nlandc) * stemdiam(nlandc)**allom3(nlandc) !Eqn 16
          ! Calculate sapling sapwood mass
          ! (17) (sapwood volume) = height * (sapwood xs area)
          ! (18) (sapwood xs area) = leafmass * sla / latosa
          ! From (17) & (18),
          ! (19) (sapwood volume) = height * leafmass * sla / latosa
          ! (20) (sapwood mass) = (wood density) * (sapwood volume)
          ! From (19) & (20),
          ! (21) (sapwood mass) = (wood density) * height * leafmass * sla / latosa

          sm_sapl(nlandc)=wooddens*height_sapl(nlandc)*lm_sapl(nlandc)*sla(nlandc)/latosa(nlandc) !Eqn 21

          ! Calculate sapling heartwood mass
          ! From (11),
          ! (22) (heartwood mass) = (x-1) * (sapwood mass)

          hm_sapl(nlandc) = (x-1.0) * sm_sapl(nlandc) !Eqn 22

       else !grass PFTs

          lm_sapl(nlandc) = pftpara(21,nlandc) / sla(nlandc)
       end if

       ! Calculate sapling or initial grass rootmass
       ! (23) lmtorm = (leafmass) / (rootmass)
       ! where lmtorm=pftpara(18,nlandc)

       rm_sapl(nlandc) = lm_sapl(nlandc) / pftpara(18,nlandc) !From Eqn 23

    end do  !end land cover type loop

end subroutine setDGVMpara
