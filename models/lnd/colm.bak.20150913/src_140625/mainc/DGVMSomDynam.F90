#include <define.h>

  subroutine DGVMSomDynam( dtime, nyr, wf, tsoi25, litterag, litterbg, &
                       cton_soil, afcton_leaf, afcton_sap, afcton_root, &
                       litter_leaf, litter_wood, litter_root, litter_repr, &
                       litter_leaf_n, litter_wood_n, litter_root_n, litter_repr_n,&
                       soil_doc, soil_nh4, cpool_fast,cpool_slow, &
                       k_fast_ave, k_slow_ave, litter_decom_ave, fmicr, afmicr, &
                       cflux_litter_soil, cflux_litter_atmos)

!-----------------------------------------------------------------------
! Calculate litter and soil decomposition.
! called every tstep.
! CALLED FROM: CLMMAIN
! Adopted from DyN (Xu-Ri, 2008) by zhq. 01/22/2010
!-----------------------------------------------------------------------
    use precision
    use phycon_module, only : tfrz
    implicit none

!   implicit in arguments

    real(r8), INTENT(in) :: dtime            ! land model time step (sec)
    real(r8), INTENT(in) :: nyr              ! land model year
    real(r8), INTENT(in) :: wf               ! soil water as frac. of whc for top 0.5 m
    real(r8), INTENT(in) :: tsoi25           ! soil temperature to 0.25 m (Kelvin)
!#if(defined DyN)
    real(r8), INTENT(in) :: cton_soil        ! soil C:N mass ratio
    real(r8), INTENT(in) :: afcton_leaf      ! annual floating leaf C:N ratio
    real(r8), INTENT(in) :: afcton_root      ! annual floating root C:N ratio
    real(r8), INTENT(in) :: afcton_sap       ! annual floating sapwood C:N ratio
!
!   implicit in/out arguments
!
    real(r8), INTENT(inout) :: litter_leaf_n    ! leaf-derived N litter for PFT on modelled area basis (gN/m2)
    real(r8), INTENT(inout) :: litter_wood_n    ! heart&sapwood-derived N litter for PFT on modelled area basis(gN/m2)
    real(r8), INTENT(inout) :: litter_root_n    ! fine root-derived N litter for PFT on modelled area basis (gN/m2)
    real(r8), INTENT(inout) :: litter_repr_n    ! litter derived from allocation to reproduction N for PFT on modelled
                                                ! area basis (gN/m2)
    real(r8), INTENT(inout) :: soil_nh4         ! 
    real(r8), INTENT(inout) :: soil_doc         ! 
!#endif
    real(r8), INTENT(inout) :: litter_leaf      ! leaf-derived litter for PFT on modelled area basis (gC/m2)
    real(r8), INTENT(inout) :: litter_wood      ! heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
    real(r8), INTENT(inout) :: litter_root      ! fine root-derived litter for PFT on modelled area basis(gC/m2)
    real(r8), INTENT(inout) :: litter_repr      ! litter derived from allocation to reproduction for PFT on modelled
                                                ! area basis (gC/m2)
    real(r8), INTENT(inout) :: litterag         ! above ground litter
    real(r8), INTENT(inout) :: litterbg         ! below ground litter
    real(r8), INTENT(inout) :: cpool_fast       ! fast carbon pool
    real(r8), INTENT(inout) :: cpool_slow       ! slow carbon pool
    real(r8), INTENT(inout) :: k_fast_ave       ! decomposition rate
    real(r8), INTENT(inout) :: k_slow_ave       ! decomposition rate
    real(r8), INTENT(inout) :: litter_decom_ave ! decomposition rate
    real(r8), INTENT(inout) :: afmicr           ! annual microbial respiration (gC /m**2/year veget'd area)
!
!   implicit out arguments
!
    real(r8), INTENT(out) :: fmicr              ! microbial respiration (umol CO2 /m**2 /s veget'd area)
    real(r8), INTENT(out) :: cflux_litter_soil  ! litter decomposition flux to soil
    real(r8), INTENT(out) :: cflux_litter_atmos ! litter decomposition flux to atmosphere

! !OTHER LOCAL VARIABLES:
!
    real(r8), parameter :: soil_equil_beg = 1000.0_r8   ! 
    real(r8), parameter :: soil_equil_end = 1500.0_r8   !
    real(r8), parameter :: k_litter10 = 0.35_r8         ! litter decomp. rate at 10 deg C (/year)
    real(r8), parameter :: k_soil_fast10 = 0.03_r8      ! fast pool decomp. rate at 10 deg C (/year)
    real(r8), parameter :: k_soil_slow10 = 0.001_r8     ! slow pool decomp. rate at 10 deg C (/year)
    real(r8), parameter :: atmfrac = 0.7_r8             ! fraction of litter decomp going directly into the atmosphere
    real(r8), parameter :: soilfrac = 1.0_r8 - atmfrac  ! fraction of litter decomp. going to soil C pools
    real(r8), parameter :: fastfrac = 0.985_r8          ! fraction of litter entering fast soil decomposition pool
    real(r8), parameter :: slowfrac = 1.0_r8 - fastfrac ! fraction of litter entering slow soil decomposition pool
    real(r8), parameter :: dmcf = 24.0                  ! co2-to-biomass conversion (g biomass / mol CO2)
    real(r8) :: temp_resp          ! temperature response of decomposition
    real(r8) :: moist_resp         ! moisture response of decomposition
    real(r8) :: k_litter           ! litter decomposition rate (/tstep)
    real(r8) :: k_fast             ! fast pool decomposition rate (/tstep)
    real(r8) :: k_slow             ! slow pool decomposition rate (/tstep)
    real(r8) :: litter_decom       ! litter C decomposition
    real(r8) :: n_litter_decom     ! litter N decomposition
    real(r8) :: cton_litter        ! litter C:N ratio
    real(r8) :: nh4_inc            ! nh4 flux increment to soil 
    real(r8) :: cflux_fast_atmos   ! soil fast pool decomposition flux to atmos.
    real(r8) :: cflux_slow_atmos   ! soil slow pool decomposition flux to atmos.
    real(r8) :: soil_equil_year    ! number of years until pool sizes for soil decomposition solved analytically

    ! Determine litter and soil decomposition

       ! Temperature response function is a modified Q10 relationship
       ! (Lloyd & Taylor 1994)
       ! slevis: Original code used monthly avg soil temp (K); I use tstep value

       if (tsoi25 <= tfrz - 40.0_r8) then !avoid division by zero
          temp_resp=0.0_r8
       else                            !Lloyd & Taylor 1994
          temp_resp=exp(308.56_r8*((1.0_r8/56.02_r8)-(1.0_r8/(tsoi25-227.13_r8))))
       end if

       ! Moisture response based on soil layer 1 moisture content (Foley 1995)
       ! slevis: Orig. code used monthly soil water in upper 0.5 m (fraction of whc)
       !         I use the tstep value

       moist_resp = 0.25_r8 + 0.75_r8 * wf

       ! Original divided by 12 to get monthly decomposition rates (k, /month)
       ! as a function of temperature and moisture
       ! slevis: make rates /tstep by dividing by the number of tsteps per year

       k_litter = k_litter10    * temp_resp * moist_resp * dtime / (86400._r8 * 365._r8)
       k_fast   = k_soil_fast10 * temp_resp * moist_resp * dtime / (86400._r8 * 365._r8)
       k_slow   = k_soil_slow10 * temp_resp * moist_resp * dtime / (86400._r8 * 365._r8)

       ! Calculate monthly litter decomposition using equation
       !   (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate
       ! from (1),
       !   (2) c = c0*exp(-kt) where c0=initial pool size
       ! from (2), decomposition in any month given by
       !   (3) delta_c = c0 - c0*exp(-k)
       ! from (3)
       !   (4) delta_c = c0*(1.0-exp(-k))

       litter_decom=(litter_leaf+litter_root+litter_wood+litter_repr)*(1.0-exp(-k_litter))  !eqn 4

       litter_leaf=litter_leaf*exp(-k_litter)
       litter_root=litter_root*exp(-k_litter)
       litter_wood=litter_wood*exp(-k_litter)
       litter_repr=litter_repr*exp(-k_litter)

       ! Calculate carbon flux to atmosphere and soil

       cflux_litter_atmos = atmfrac  * litter_decom
       cflux_litter_soil  = soilfrac * litter_decom

       ! Further subdivide soil fraction between fast and slow soil pools

       cpool_fast = cpool_fast + fastfrac * cflux_litter_soil

       ! add a "fast" spin up process for cpool_slow accumulation 12/28/2009. zhq
       !   if((slowfrac*cflux_litter_soil-k_slow*cpool_slow).gt.0..and.(cpool_slow-soilfrac * slowfrac * litter_decom_ave / k_slow_ave).lt.0.)then
       !     cpool_slow = cpool_slow + 100. * slowfrac * cflux_litter_soil
       ! else
       cpool_slow = cpool_slow + slowfrac * cflux_litter_soil
       !  endif

       ! Calculate soil decomposition to the atmosphere

       cflux_fast_atmos = cpool_fast * (1.0_r8-exp(-k_fast))  !eqn 4
       cflux_slow_atmos = cpool_slow * (1.0_r8-exp(-k_slow))  !eqn 4

!#if(defined DyN)
       ! Increment doc availability
       soil_doc = soil_doc + cflux_fast_atmos + cflux_slow_atmos

       !nitrogen content in the litter of this pft

       !n_decomp_litter+=(litter_leaf/pft.cton_leaf+^M
       !        //litter_root/pft.cton_root+^M
       !        //litter_wood/pft.cton_sap+^M
       !        //litter_repr/pft.cton_leaf)*(1.0-fr_litter);^M

       n_litter_decom=(litter_leaf*afcton_leaf &
                       + litter_root*afcton_root &
                       + litter_wood*afcton_sap &
                       + litter_repr*afcton_leaf)*(1.0-exp(-k_litter))

       !/*n_decomp_litter+=(litter_leaf_n+^M
       !        litter_root_n+^M
       !        litter_wood_n+^M
       !        litter_repr_n)*(1.0-fr_litter);*/^M

       litter_leaf_n=litter_leaf_n*exp(-k_litter)
       litter_root_n=litter_root_n*exp(-k_litter)
       litter_wood_n=litter_wood_n*exp(-k_litter)
       litter_repr_n=litter_repr_n*exp(-k_litter)

       !if (date.year>=soil.soiltype.solvesom_begin)
       !soil.decomp_litter_mean+=decomp_litter

       ! Partition litter decomposition among fast and slow SOM pools
       ! and flux to atmosphere
       ! C to N ratio of litter
       if (n_litter_decom<=1.0E-6) then
           cton_litter=cton_soil
       else 
           cton_litter=litter_decom/n_litter_decom                                 
       endif

!print*,'litter',litter_decom,n_litter_decom
       ! NH4 released to soil with litter decompositon
!print*,'cflux',cflux_litter_atmos,'cton_litter',cton_litter
       nh4_inc=cflux_litter_atmos/cton_litter
!print*,'soil_nh4-2',soil_nh4,nh4_inc

       ! Increment NH4 
       nh4_inc=nh4_inc+(cflux_fast_atmos+cflux_slow_atmos)/cton_soil
!print*,'soil_nh4-3',soil_nh4,nh4_inc
                               
       ! Increment nh4 to soil inorganic N pool by decomposition
       soil_nh4 = soil_nh4 + nh4_inc
!#endif
!print*,'soil_nh4',soil_nh4,nh4_inc

       ! Update the soil pools

       cpool_fast = cpool_fast - cflux_fast_atmos
       cpool_slow = cpool_slow - cflux_slow_atmos

       ! Update the litter pools. zhq. 05/12/2010
       litterag = litter_leaf + litter_wood + litter_repr
       litterbg = litter_root

       ! Calculate heterotrophic respiration (mol CO2 /m2/s vegetated area)

       fmicr = cflux_litter_atmos + cflux_fast_atmos + cflux_slow_atmos
       afmicr = afmicr + fmicr
       fmicr = fmicr * 2.0 / dmcf / dtime

       ! Empty soil pools below a minimum threshold

       if (cpool_fast < 1.0e-5_r8) cpool_fast = 0.0_r8
       if (cpool_slow < 1.0e-5_r8) cpool_slow = 0.0_r8

       ! Update running average respiration rates and litter input
       ! slevis: had to multiply the denominator to chg units from years to tsteps

       if (nyr > soil_equil_beg .and. nyr <= soil_equil_end) then

          soil_equil_year = soil_equil_end - soil_equil_beg 

          k_fast_ave       = k_fast_ave       + k_fast / &
               (soil_equil_year * 365._r8 * 86400. / dtime)
          k_slow_ave       = k_slow_ave       + k_slow / &
               (soil_equil_year * 365._r8 * 86400. / dtime)
          litter_decom_ave = litter_decom_ave + litter_decom / &
               (soil_equil_year * 365._r8 * 86400. / dtime)

       else if (abs(nyr - soil_equil_end - 1) < 1.0e-5) then

          ! SOIL DECOMPOSITION EQUILIBRIUM CALCULATION
          ! Analytical solution of differential flux equations for fast and slow
          ! soil carbon pools.  Implemented after (soil_equil_year) simulation
          ! years, when annual litter inputs should be close to equilibrium.  Assumes
          ! average climate (temperature and soil moisture) from all years up to
          ! soil_equil_year.
          ! slevis: next could be done once

          ! Analytically calculate pool sizes this year only
          ! Rate of change of soil pool size = litter input - decomposition
          !   (5) dc/dt = litter_decom - kc
          ! At equilibrium,
          !   (6) dc/dt = 0
          ! From (5) & (6),
          !   (7) c = litter_decom / k

          cpool_fast = soilfrac * fastfrac * litter_decom_ave / k_fast_ave !eqn 7
          cpool_slow = soilfrac * slowfrac * litter_decom_ave / k_slow_ave !eqn 7

       end if

  end subroutine DGVMSomDynam
