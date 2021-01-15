  subroutine lpjave(nyr_r       ,npatch        ,npftpar       ,pftpar          ,&
                    vegclass    ,bm_inc        ,fpcgrid       ,agdd            ,&
                    t_mo_min    ,tmomin20      ,agdd20        ,afmicr          ,&
                    ifpre       ,lm_ind        ,hm_ind        ,sm_ind          ,&
                    rm_ind      ,litter_ag     ,litter_bg     ,turnover_ind    ,&
                    lai_ind     ,fpc_inc       ,crownarea     ,htop            ,&
                    nind        ,ivt           ,annpsn        ,annpsnpot       ,&
                    sla         ,agddtw        ,firelength                     ,&
#ifdef IAPDGVM
                    afirefrac1, nfireg1,&
#endif     
                    lm_sapl     ,sm_sapl       ,hm_sapl       ,rm_sapl, wt_patch,&
                    wt_column   ,prec365       ,cpool_fast    ,cpool_slow      ,&
                    vegc        ,litc          ,litc_ag       ,litc_bg         ,&
                    soic        ,soic_fast     ,soic_slow     ,&
                    litter2soil_ind, flitter2soil,&
                    litter2atmos_ind, flitter2atmos,&
                    leafc       ,woodc         ,rootc         ,dz)

!---------------------------------------------------------------------

   use precision
   use paramodel, only: maxsnl, nl_soil
   use timemgr, only: dtime
   implicit none

! INTENT IN VARIAVLES:

   integer, INTENT(in) :: npatch                   ! total patch number in this grid
   integer, INTENT(in) :: npftpar                  ! number of pft parameter
   real(r8),INTENT(in) :: nyr_r                    ! counting the model years
   real(r8),INTENT(in) :: pftpar(npftpar,npatch)   ! 32 pft parameters
   real(r8),INTENT(in) :: vegclass(npatch)         ! 1.tree 2.shrub 3.grass 4.crop -1.others
   real(r8),INTENT(in) :: agdd(npatch)             ! accumulated growing degree days above 5
   real(r8),INTENT(in) :: t_mo_min(npatch)         ! annual min of t_mo (Kelvin)
   real(r8),INTENT(in) :: annpsn(npatch)           ! annual photosynthesis (umol CO2 /m**2)
   real(r8),INTENT(in) :: annpsnpot(npatch)        ! annual potential photosynthesis (..)
   real(r8),INTENT(in) :: sla(npatch)              ! specific leaf area [m2 leaf g-1 carbon]
   real(r8),INTENT(in) :: agddtw(npatch)           ! accumulated growing degree days above twmax
   real(r8),INTENT(in) :: firelength(npatch)       ! fire season in days
#ifdef IAPDGVM
   real(r8),INTENT(in) :: afirefrac1(npatch)          ! IAPDGVM
   real(r8),INTENT(in) :: nfireg1(npatch)          ! IAPDGVM
#endif
   real(r8),INTENT(in) :: lm_sapl(npatch)          ! ecophys const - leaf mass of sapling
   real(r8),INTENT(in) :: sm_sapl(npatch)          ! ecophys const - stem mass of sapling
   real(r8),INTENT(in) :: hm_sapl(npatch)          ! ecophys const - heartwood mass of sapling
   real(r8),INTENT(in) :: rm_sapl(npatch)          ! ecophys const - root mass of saping
   real(r8),INTENT(in) :: cpool_fast(nl_soil,npatch)! fast soil C pool (gC/m2 veget'd area)
   real(r8),INTENT(in) :: cpool_slow(nl_soil,npatch)! slow soil C pool (gC/m2 veget'd area)
   real(r8),INTENT(in) :: prec365                  ! yearly running mean of precipitation [mm/s]
   real(r8),INTENT(in) :: wt_column                ! relative weight of natural vegetation to grid
   real(r8),INTENT(in) :: dz(maxsnl+1:nl_soil)     ! layer thickness [m]

! INTENT INOUT VARIABLES:
   real(r8),INTENT(in) :: wt_patch(npatch)         ! pft weight relative to grid cell
   real(r8),INTENT(in) :: tmomin20(npatch)         ! 20-yr running mean of tmomin
   real(r8),INTENT(in) :: agdd20(npatch)           ! 20-yr running mean of agdd
   real(r8),INTENT(in) :: bm_inc(npatch)           ! biomass increment
   real(r8),INTENT(in) :: fpcgrid(npatch)          ! foliar projective cover on gridcell
   real(r8),INTENT(in) :: afmicr(npatch)           ! annual microbial respiration
   real(r8),INTENT(in) :: ifpre(npatch)            ! whether pft present in this patch
   real(r8),INTENT(in) :: lm_ind(npatch)           ! individual leaf mass
   real(r8),INTENT(in) :: sm_ind(npatch)           ! individual sapwood mass
   real(r8),INTENT(in) :: hm_ind(npatch)           ! individual heartwood mass
   real(r8),INTENT(in) :: rm_ind(npatch)           ! individual root mass
   real(r8),INTENT(in) :: turnover_ind(npatch)     ! individual turnover biomass
   real(r8),INTENT(in) :: litter_ag(npatch)        ! above ground litter mass
   real(r8),INTENT(in) :: litter_bg(npatch)        ! below ground litter mass
   real(r8),INTENT(in) :: lai_ind(npatch)          ! max lai for individual
   real(r8),INTENT(in) :: fpc_inc(npatch)          ! fpc increase
   real(r8),INTENT(in) :: crownarea(npatch)        ! area each individual tree takes up (m^2)
   real(r8),INTENT(in) :: htop(npatch)             ! canopy top 
   real(r8),INTENT(in) :: nind(npatch)             ! population density over gridcell 
   integer ,INTENT(in) :: ivt(npatch)              ! land cover type 
   real(r8),INTENT(in) :: litter2soil_ind(npatch)
   real(r8),INTENT(in) :: litter2atmos_ind(npatch)
  
! INTENT OUT VARIABLES:
   real(r8),INTENT(out) :: vegc                    ! gridcell vegetaion biomass
   real(r8),INTENT(out) :: litc                    ! gridcell litter
   real(r8),INTENT(out) :: litc_ag                 ! gridcell above ground litter
   real(r8),INTENT(out) :: litc_bg                 ! gridcell below ground litter
   real(r8),INTENT(out) :: soic                    ! gridcell soil cpool
   real(r8),INTENT(out) :: soic_fast               ! gridcell fast soil cpool
   real(r8),INTENT(out) :: soic_slow               ! gridcell slow soil cpool(gC/m2 veget'd area)
   real(r8),INTENT(out) :: flitter2soil            ! 
   real(r8),INTENT(out) :: flitter2atmos           ! 
   real(r8),INTENT(out) :: leafc                   ! 
   real(r8),INTENT(out) :: woodc                   ! 
   real(r8),INTENT(out) :: rootc                   ! 

! LOCAL VARIABLES:
   integer  :: g,p,fp                              ! indices
   integer  :: nyr

    leafc = 0._r8
    woodc = 0._r8
    rootc = 0._r8
    vegc = 0._r8
    litc = 0._r8
    litc_ag = 0._r8
    litc_bg = 0._r8
    soic = 0._r8
    soic_fast = 0._r8
    soic_slow = 0._r8
    flitter2soil = 0._r8
    flitter2atmos= 0._r8

    do p = 1, npatch
       if(wt_patch(p)>0._r8 .and. wt_patch(p) <= 1._r8+1.e-6_r8) then  ! BuildFilter
          leafc = leafc + lm_ind(p)*nind(p)
          woodc = woodc + (sm_ind(p) + hm_ind(p))*nind(p)
          rootc = rootc + rm_ind(p)*nind(p)
          litc_ag = litc_ag + litter_ag(p)
          litc_bg = litc_bg + litter_bg(p)
          soic_fast = soic_fast + sum(cpool_fast(:,p))
          soic_slow = soic_slow + sum(cpool_slow(:,p))
          flitter2soil = flitter2soil + litter2soil_ind(p)
          flitter2atmos = flitter2atmos + litter2atmos_ind(p)
       end if
    end do

    vegc = leafc + woodc + rootc

    litc = litc_ag + litc_bg

    soic = soic_fast + soic_slow

    flitter2soil = flitter2soil/dtime
    flitter2atmos = flitter2atmos/dtime

  end subroutine lpjave
