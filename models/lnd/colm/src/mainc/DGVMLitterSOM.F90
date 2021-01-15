
  subroutine LitterSOM( dtime, nyr, z, dz, wf, tss, rootfr, litterag, litterbg, cpool_fast, &
                        cpool_slow, k_fast_ave, k_slow_ave, litter_decom_ave, fmicr, afmicr, &
                        cflux_litter_soil, cflux_litter_atmos, &
                        Isf, Iss, Ksf, Kss)

!-----------------------------------------------------------------------
! Calculate litter and soil decomposition.
! called every tstep.
! CALLED FROM: CLMMAIN
!-----------------------------------------------------------------------
    use precision
    use phycon_module, only : tfrz
    use paramodel, only : nl_soil
    use colm_varctl, only : enable_cspinup, nstep_cspinup
    use timemgr, only : get_nstep
    implicit none

!   implicit in arguments

    real(r8), INTENT(in) :: dtime            ! land model time step (sec)
    real(r8), INTENT(in) :: nyr              ! land model year
    real(r8), INTENT(in) ::  z(nl_soil)      ! node depth [m]
    real(r8), INTENT(in) :: dz(nl_soil)      ! layer thickiness [m]
    real(r8), INTENT(in) :: wf(nl_soil)      ! soil water as frac. of whc
    real(r8), INTENT(in) :: tss(nl_soil)     ! soil temperature to 0.25 m (Kelvin)
    real(r8), INTENT(in) :: rootfr(nl_soil)  ! fraction of roots in each soil layer
!
!   implicit in/out arguments
!
    real(r8), INTENT(inout) :: litterag         ! above ground litter
    real(r8), INTENT(inout) :: litterbg         ! below ground litter
    real(r8), INTENT(inout) :: cpool_fast(nl_soil)       ! fast carbon pool
    real(r8), INTENT(inout) :: cpool_slow(nl_soil)       ! slow carbon pool
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
!
!   fast spinup
!
    real(r8), INTENT(out) :: Isf(nl_soil)
    real(r8), INTENT(out) :: Iss(nl_soil)
    real(r8), INTENT(out) :: Ksf(nl_soil)
    real(r8), INTENT(out) :: Kss(nl_soil)

! !OTHER LOCAL VARIABLES:
!
    real(r8), parameter :: soil_equil_year = -800.0_r8  ! number of years until pool sizes for soil decomposition solved analytically
    real(r8), parameter :: k_litter10 = 0.35_r8         ! litter decomp. rate at 10 deg C (/year)
    real(r8), parameter :: k_soil_fast10 = 0.03_r8      ! fast pool decomp. rate at 10 deg C (/year)
    real(r8), parameter :: k_soil_slow10 = 0.001_r8     ! slow pool decomp. rate at 10 deg C (/year)
    real(r8), parameter :: atmfrac = 0.7_r8             ! fraction of litter decomp going directly into the atmosphere
    real(r8), parameter :: soilfrac = 1.0_r8 - atmfrac  ! fraction of litter decomp. going to soil C pools
    real(r8), parameter :: fastfrac = 0.985_r8          ! fraction of litter entering fast soil decomposition pool
    real(r8), parameter :: slowfrac = 1.0_r8 - fastfrac ! fraction of litter entering slow soil decomposition pool
#ifdef IAPDGVM
    real(r8), parameter :: dmcf = 28.5                  ! co2-to-biomass conversion (g biomass / mol CO2)
#else    
    real(r8), parameter :: dmcf = 24.0                  ! co2-to-biomass conversion (g biomass / mol CO2)
#endif
    real(r8) :: temp_resp(10)      ! temperature response of decomposition
    real(r8) :: moist_resp(10)     ! moisture response of decomposition
    real(r8) :: k_litter(10)       ! litter decomposition rate (/tstep)
    real(r8) :: k_fast(10)         ! fast pool decomposition rate (/tstep)
    real(r8) :: k_slow(10)         ! slow pool decomposition rate (/tstep)
    real(r8) :: litter_decom       ! litter decomposition
    real(r8) :: litter_decom_ag    ! above-ground litter decomposition
    real(r8) :: litter_decom_bg    ! below-ground litter decomposition
    real(r8) :: cflux_fast_atmos   ! soil fast pool decomposition flux to atmos.
    real(r8) :: cflux_slow_atmos   ! soil slow pool decomposition flux to atmos.

    real(r8) :: litter_decom_bg_layer(10)
    real(r8) :: cflux_fast_atmos_layer(10)
    real(r8) :: cflux_slow_atmos_layer(10)
    real(r8) :: cflux_litter_soil_layer(10)

    integer  j

    if (maxval(rootfr).gt.1.) stop 'error in rootfr 1'
    if (abs(sum(rootfr)-1).gt.1.e-6) stop 'error in rootfr 2'
    if (nl_soil.lt.10) stop 'erro in nl_soil'

    ! Determine litter and soil decomposition for top 10 soil layers

    do j = 1, 10

       ! Temperature response function is a modified Q10 relationship
       ! (Lloyd & Taylor 1994)
       ! slevis: Original code used monthly avg soil temp (K); I use tstep value

       if (tss(j) <= tfrz - 40.0_r8) then !avoid division by zero
          temp_resp(j)=0.0_r8
       else                            !Lloyd & Taylor 1994
          temp_resp(j)=exp(308.56_r8*((1.0_r8/56.02_r8)-(1.0_r8/(tss(j)-227.13_r8))))
       end if

       ! Moisture response based on soil layer 1 moisture content (Foley 1995)
       ! slevis: Orig. code used monthly soil water in upper 0.5 m (fraction of whc)
       !         I use the tstep value

       moist_resp(j) = 0.25_r8 + 0.75_r8 * wf(j)

       ! Original divided by 12 to get monthly decomposition rates (k, /month)
       ! as a function of temperature and moisture
       ! slevis: make rates /tstep by dividing by the number of tsteps per year

       k_litter(j) = k_litter10    * temp_resp(j) * moist_resp(j) * dtime / (86400._r8 * 365._r8)
       k_fast(j)   = k_soil_fast10 * temp_resp(j) * moist_resp(j) * dtime / (86400._r8 * 365._r8)
       k_slow(j)   = k_soil_slow10 * temp_resp(j) * moist_resp(j) * dtime / (86400._r8 * 365._r8)

    end do

    ! Calculate monthly litter decomposition using equation
    !   (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate
    ! from (1),
    !   (2) c = c0*exp(-kt) where c0=initial pool size
    ! from (2), decomposition in any month given by
    !   (3) delta_c = c0 - c0*exp(-k)
    ! from (3)
    !   (4) delta_c = c0*(1.0-exp(-k))

    litter_decom_ag = litterag * (1.0_r8-exp(-k_litter(1)))  !eqn 4

    do j = 1, 10
       litter_decom_bg_layer(j) = litterbg * rootfr(j) * (1.0_r8-exp(-k_litter(j)))
    end do

    litter_decom_bg = sum(litter_decom_bg_layer)

    litter_decom    = litter_decom_ag + litter_decom_bg

    ! Update the litter pools

    litterag = litterag - litter_decom_ag
    litterbg = litterbg - litter_decom_bg

    ! Empty litter pools below a minimum threshold, zhq 03/22/2010
    if (litterag < 1.0e-5_r8) litterag = 0.0_r8
    if (litterbg < 1.0e-5_r8) litterbg = 0.0_r8

    ! Calculate carbon flux to atmosphere and soil

    cflux_litter_atmos = atmfrac  * litter_decom

    cflux_litter_soil  = soilfrac * litter_decom

    do j = 1, 10
       cflux_litter_soil_layer(j) = soilfrac * litter_decom_bg_layer(j)
    end do
    cflux_litter_soil_layer(1) = cflux_litter_soil_layer(1) + soilfrac * litter_decom_ag

    if (abs(cflux_litter_soil - sum(cflux_litter_soil_layer)) > 1.e-6) stop 'error in cflux_litter_soil' 

    ! Further subdivide soil fraction between fast and slow soil pools

    do j = 1, 10
       cpool_fast(j) = cpool_fast(j) + fastfrac * cflux_litter_soil_layer(j)
       cpool_slow(j) = cpool_slow(j) + slowfrac * cflux_litter_soil_layer(j)

       Isf(j) = Isf(j) + fastfrac * cflux_litter_soil_layer(j)
       Iss(j) = Iss(j) + slowfrac * cflux_litter_soil_layer(j)
    end do

    ! add a "fast" spin up process for cpool_slow accumulation 12/28/2009. zhq
  ! if((slowfrac*cflux_litter_soil-k_slow*cpool_slow).gt.0..and.(cpool_slow-soilfrac * slowfrac * litter_decom_ave / k_slow_ave).lt.0.)then
  !    cpool_slow = cpool_slow + 100. * slowfrac * cflux_litter_soil
  ! else
  !    cpool_slow = cpool_slow + slowfrac * cflux_litter_soil_layer
  ! endif
 
    ! Calculate monthly soil decomposition to the atmosphere

    do j = 1, 10
       cflux_fast_atmos_layer(j) = cpool_fast(j) * (1.0_r8-exp(-k_fast(j)))  !eqn 4
       cflux_slow_atmos_layer(j) = cpool_slow(j) * (1.0_r8-exp(-k_slow(j)))  !eqn 4

       Ksf(j) = Ksf(j) + 1.0_r8-exp(-k_fast(j))
       Kss(j) = Kss(j) + 1.0_r8-exp(-k_slow(j))
    end do

    cflux_fast_atmos = sum(cflux_fast_atmos_layer)
    cflux_slow_atmos = sum(cflux_slow_atmos_layer)

    ! Update the soil pools

    cpool_fast(1:10) = cpool_fast(1:10) - cflux_fast_atmos_layer(1:10)
    cpool_slow(1:10) = cpool_slow(1:10) - cflux_slow_atmos_layer(1:10)

    ! Calculate heterotrophic respiration (mol CO2 /m2/s vegetated area)

    fmicr = cflux_litter_atmos + cflux_fast_atmos + cflux_slow_atmos
    afmicr = afmicr + fmicr
    fmicr = fmicr * 2.0 / dmcf / dtime

    ! Empty soil pools below a minimum threshold

    if (sum(cpool_fast) < 1.0e-5_r8) cpool_fast = 0.0_r8
    if (sum(cpool_slow) < 1.0e-5_r8) cpool_slow = 0.0_r8

    call cdiff (dtime, z, dz, tss, cpool_fast, cpool_slow)

    if (enable_cspinup .and. MOD(get_nstep(),nstep_cspinup) == 0) then
     ! print *, 'cpool_fast', sum(cpool_fast)
     ! print *, 'cpool_slow', sum(cpool_slow)
   
       where (Ksf.ne.0) cpool_fast = Isf/Ksf
       where (Kss.ne.0) cpool_slow = Iss/Kss

       Isf = 0._r8
       Iss = 0._r8
       Ksf = 0._r8
       Kss = 0._r8
   
     ! do p = 1, npatch
     !    if (isnan(cpool_fast(p)) .or. isnan(cpool_slow(p))) then
     !       write(6,*) 'pool nan', cpool_fast(p), cpool_slow(p), & 
     !                  Isf(p), Iss(p), Ksf(p), Kss(p)
     !    end if
     ! end do
   
     ! print *, 'cpool_fast', sum(cpool_fast)
     ! print *, 'cpool_slow', sum(cpool_slow)
    end if

#ifdef SOILEQUIL
    if ((nyr - soil_equil_year) <= 0.) then

       ! Update running average respiration rates and litter input
       ! slevis: had to multiply the denominator to chg units from years to tsteps
       k_fast_ave       = k_fast_ave       + k_fast / &
            (soil_equil_year * 365._r8 * 86400. / dtime)
       k_slow_ave       = k_slow_ave       + k_slow / &
            (soil_equil_year * 365._r8 * 86400. / dtime)
       litter_decom_ave = litter_decom_ave + litter_decom / &
            (soil_equil_year * 365._r8 * 86400. / dtime)

    else if (abs(nyr - soil_equil_year - 1) < 1.0e-5) then

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
#endif

  end subroutine LitterSOM
