#include <define.h>
 
  subroutine LPJDRIVER
 
 !=======================================================================
 !
 ! LPJ MODEL DRIVER
 !
 ! Adopted from NCAR-CLM3.0 : Ming Chen 05/10/2008
 !
 !=======================================================================
 
  use precision
  use paramodel
  use spmd
  use spmd_decomp, only: cgmap, ggmap, gmask, numgrid_proc, gxmap, gymap
  use colm_varMod, only: lon_points, lat_points, numpatch, numcolumn, numgrid, &
                         cvar, pvar, fldv_dgvm, forc, nep_residual, lnep_adjust, &
                         fldv, fLitterSoil, fLitterAtmos, fldv_dgvm_init
  use timemgr, only: idate, dtime, is_newyear
  use landuse, only: landuse_cflux
  use debug, only: c_bug
  implicit none
 
 ! ----------------------------------------------------------------
 ! I. Time invariant model variables
 ! ----------------------------------------------------------------
 
   integer , pointer :: itypwat                ! land water type
   real(r8), pointer :: dlat                   ! latitude in radians
   real(r8), pointer :: dlon                   ! longitude in radians

   real(r8), pointer :: dz  (:)                ! layer thickness [m]
!  real(r8), pointer :: wliq(:)                ! liquid water in layers [kg/m2]

   integer , pointer :: ivt(:)                 !land cover type  
   real(r8), pointer :: t10min(:)              !annual minimum of 10-day running mean (K)
   real(r8), pointer :: lai_ind(:)             !LAI per individual
   real(r8), pointer :: dphen(:)               !phenology [0 to 1]
   real(r8), pointer :: leafon(:)              !leafon days
   real(r8), pointer :: leafof(:)              !leafoff days
   real(r8), pointer :: firelength(:)          !fire season in days
#ifdef IAPDGVM
   real(r8), pointer :: afirefrac1(:)
   real(r8), pointer :: nfireg1(:)
   real(r8), pointer :: wliq6mon         ! IAPDGVM liquid water 6 mons for first3 layers
#endif
   real(r8), pointer :: litterag(:)            !above ground litter
   real(r8), pointer :: litterbg(:)            !below ground litter
   real(r8), pointer :: cpool_fast(:,:)        !fast carbon pool
   real(r8), pointer :: cpool_slow(:,:)        !slow carbon pool
   real(r8), pointer :: k_fast_ave(:)          !decomposition rate
   real(r8), pointer :: k_slow_ave(:)          !decomposition rate
   real(r8), pointer :: litter_decom_ave(:)    !decomposition rate
   real(r8), pointer :: fmicr(:)               !microbial respiration (umol CO2 /m**2 /s)
   real(r8), pointer :: nind(:)                !number of individuals (*/m**2)
   real(r8), pointer :: lm_ind(:)              !individual leaf mass
   real(r8), pointer :: sm_ind(:)              !individual sapwood mass
   real(r8), pointer :: hm_ind(:)              !individual heartwood mass
   real(r8), pointer :: rm_ind(:)              !individual root mass
   real(r8), pointer :: tmomin20(:)            !20-yr running mean of tmomin
   real(r8), pointer :: agdd0(:)               !growing dgree days above 0
   real(r8), pointer :: agdd(:)                !growing dgree days above 5
   real(r8), pointer :: agddtw(:)              !growing dgree days above twmax
   real(r8), pointer :: agdd20(:)              !20-yr running mean of agdd
   real(r8), pointer :: t_mo(:)                !30-day mean temperature of 2m (K)  
   real(r8), pointer :: t_mo_sum(:)            !30-day accumulated temperature of 2m (K)
   real(r8), pointer :: t_mo_min(:)            !annual min of t_mo (Kelvin)
   real(r8), pointer :: crownarea(:)           !area that each individual tree takes up (m2)
   real(r8), pointer :: htop(:)                !canopy top
   real(r8), pointer :: tsai(:)                !one-sided stem area index, no burying by snow
   real(r8), pointer :: fpcgrid(:)             !foliar projective cover on gridcell (fraction)
   real(r8), pointer :: bm_inc(:)              !biomass increment
   real(r8), pointer :: afmicr(:)              !annual microbial respiration
   real(r8), pointer :: annpsn(:)              !annual photosynthesis (umol CO2 /m**2)
   real(r8), pointer :: annpsnpot(:)           !annual potential photosynthesis (same units)
   real(r8), pointer :: tref10(:)              !10-day averaged temperature at 2m
   real(r8), pointer :: tref_sum(:)            !sum of temperature in current day
   real(r8), pointer :: t10(:,:)               !array ro record the 10 day temperature
   real(r8), pointer :: assimn10(:)            !10-day averaged assimilation rate
   real(r8), pointer :: assimn_sum(:)          !sum of assimn of current day
   real(r8), pointer :: an10(:,:)              !arry to record 10 day assimn
   real(r8), pointer :: anngpp(:)              !annual gpp
   real(r8), pointer :: annfrmf(:)             !annual frmf
   real(r8), pointer :: annfrms(:)             !annual frms
   real(r8), pointer :: annfrmr(:)             !annual frmr
   real(r8), pointer :: annfrg(:)              !annual frg
   real(r8), pointer :: turnover_ind(:)        !individual turnover biomass
!  real(r8), pointer :: litter_ag(:)           !above ground litter mass
!  real(r8), pointer :: litter_bg(:)           !below ground litter mass
!  real(r8), pointer :: litter_leaf(:)         !leaf-derived litter for PFT on modelled area basis (gC/m2)
!  real(r8), pointer :: litter_wood(:)         !heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
!  real(r8), pointer :: litter_root(:)         !fine root-derived litter for PFT on modelled area basis(gC/m2)
!  real(r8), pointer :: litter_repr(:)         !litter derived from allocation to reproduction for PFT on modelled
   real(r8), pointer :: fpc_inc(:)             !fpc increase
   real(r8), pointer :: pftpara(:,:)           !32 parameters of PFTs
   real(r8), pointer :: vegclass(:)            !1.tree 2.shrub 3.grass 4.crop -1.others
   real(r8), pointer :: summergreen(:)         !1. for summergreen; otherwise -1.
   real(r8), pointer :: raingreen(:)           !1. for raingreen; otherwise -1.
   real(r8), pointer :: sla(:)                 !sla
   real(r8), pointer :: lm_sapl(:)             !leafmass
   real(r8), pointer :: sm_sapl(:)             !sapwood mass
   real(r8), pointer :: hm_sapl(:)             !heartwood mass
   real(r8), pointer :: rm_sapl(:)             !rootmass
   real(r8), pointer :: ifpre(:)               !-1=no PFT present;1=PFT present in this grid

   real(r8), pointer :: nday                   !counting the model days
   real(r8), pointer :: nyr                    !counting the model years
   real(r8), pointer :: prec365                !yearly running mean of precipitation(mm/s)

   real(r8), pointer :: wt_pft(:)              !weight of each pft in a grid
   real(r8), pointer :: wt_col                 !weight of column in a grid

#if(defined DyN)
   real(r8), pointer :: cton_soil(:)        ! soil C:N mass ratio

   real(r8), pointer :: litter_leaf(:)      ! leaf-derived litter for PFT on modelled area basis (gC/m2)
   real(r8), pointer :: litter_wood(:)      ! heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
   real(r8), pointer :: litter_root(:)      ! fine root-derived litter for PFT on modelled area basis(gC/m2)
   real(r8), pointer :: litter_repr(:)      ! litter derived from allocation to reproduction for PFT on modelled

   real(r8), pointer :: litter_leaf_n(:) ! leaf-derived N litter for PFT on modelled area basis (gN/m2)
   real(r8), pointer :: litter_wood_n(:) ! heart&sapwood-derived N litter for PFT on modelled area basis(gN/m2)
   real(r8), pointer :: litter_root_n(:) ! fine root-derived N litter for PFT on modelled area basis (gN/m2)
   real(r8), pointer :: litter_repr_n(:) ! litter derived from allocation to reproduction N for PFT on modelled
                                         ! area basis (gN/m2)
   real(r8), pointer :: afcton_leaf(:)   ! annual floating leaf C:N ratio
   real(r8), pointer :: afcton_root(:)   ! annual floating root C:N ratio
   real(r8), pointer :: afcton_sap(:)    ! annual floating sapwood C:N ratio
   real(r8), pointer :: lm_ind_n(:)      ! individual leaf nitrogen mass
   real(r8), pointer :: sm_ind_n(:)      ! individual sapwood nitrogen mass
   real(r8), pointer :: hm_ind_n(:)      ! individual heartwood nitrogen mass
   real(r8), pointer :: rm_ind_n(:)      ! individual root nitrogen mass
                                         ! gN/m2 veget'd area for each pft
   real(r8), pointer :: an_up(:)         ! annual plant nitrogen uptake(gN/m2 vegt'd area)
   real(r8), pointer :: an_stress(:)     ! annual plant nitrogen stress(-)

   real(r8), pointer :: soil_no3
   real(r8), pointer :: soil_no2
   real(r8), pointer :: soil_no
   real(r8), pointer :: soil_n2o
   real(r8), pointer :: soil_n2
   real(r8), pointer :: soil_nh4
#endif

 ! -----------------------------------------------------------------
 ! V. Local declaration
 ! -----------------------------------------------------------------

   integer i,j,lb,ub,lc,uc,jm,n           ! loop/array indices
 
   integer npatch
   integer g, c, p
   integer p1, p2 

   logical newyear

   integer  :: ivt_old(17)
 
   real(r8) :: ifpre_old(17), wt_pft_old(17)

   real(r8) :: leafc, woodc, rootc, vegc
   real(r8) :: litc, litc_ag, litc_bg
   real(r8) :: soic, soic_fast, soic_slow
   real(r8) :: fveg2litter, flitter2soil, flitter2atmos
   real(r8) :: gpp, npp, nep, nbp
   real(r8) :: ra, rh, ffirec
   real(r8), pointer :: litter2soil_ind(:), litter2atmos_ind(:)

   real(r8) :: lucflux(lon_points,lat_points)

 !-------------------------------------------------------------------
 ! variables for histDGVM
 !-------------------------------------------------------------------
   real(r8) afirec, afiref, avegc, aestabc, anpp, amrh, alitc_ag, alitc_bg, asoic_fast, asoic_slow
 ! #ifdef IAPDGVM
  !  real(r8) emCO2_gcell ,emCO_gcell    ,emCH4_gcell   ,emNHMC_gcell  ,&
  !           emH2_gcell  ,emNOx_gcell   ,emN2O_gcell   ,emPM25_gcell  ,&
  !           emTPM_gcell ,emTC_gcell    ,emOC_gcell    ,emBC_gcell 
!#endif 
   real(r8) npp_ind(17)
   real(r8) gpp_ind(17)
   real(r8) frmf_ind(17)
   real(r8) frms_ind(17)
   real(r8) frmr_ind(17)
   real(r8) frg_ind(17)
#if(defined DyN) 
   real(r8) an_up_total, an_stress_total, avegn, alitn_ag, alitn_bg, asoin
#endif
 !---------------------------------------------------------------------------------------
 !LPJ is called on the gridcell level, because there are communicates between
 !patches in the anual processes.(light competetion, for example). Here we assume
 !that there are 21 patches in each gridcell. If one patch do not exist, then the
 !weight of the patch is set to zero. This part can be improved in further development.
 !---------------------------------------------------------------------------------------

    newyear = is_newyear()

    if(newyear) lnep_adjust = .true.

    call fldv_dgvm_init

    nep_residual(:) = 0.

    call landuse_cflux(idate(1),lucflux)

  ! g(C)/m2/s -> kg(C)/m2/s
    lucflux(:,:) = lucflux(:,:)*1.e-3_r8

    p1 = 1
    p2 = 0
 
    DO c = 1, numcolumn

       itypwat           => cvar%itypwat(c)
       dz                => cvar%dz(maxsnl+1:nl_soil,c) ! snowlayer + soillayer
       nday              => cvar%nday(c)
       nyr               => cvar%nyr(c)
       prec365           => cvar%prec365(c)
#ifdef IAPDGVM
       wliq6mon          => cvar%wliq6mon(c)
#endif
       wt_col            => cvar%wxy_column(c)

#ifdef DyN
       soil_no3          => cvar%soil_no3(c)
       soil_no2          => cvar%soil_no2(c)
       soil_no           => cvar%soil_no(c)
       soil_n2o          => cvar%soil_n2o(c)
       soil_n2           => cvar%soil_n2(c)
       soil_nh4          => cvar%soil_nh4(c)
#endif

       if(itypwat==0)      npatch = 16 + 1  ! natural vegetation + bare soil
       if(itypwat==1)      npatch = 1       ! urban and built-up
       if(itypwat==2)      npatch = 1       ! wetland
       if(itypwat==3)      npatch = 1       ! land ice
       if(itypwat==4)      npatch = 1       ! river or deep lake
       if(itypwat==99)then                  ! ocean
          write(6,*) 'ocean column', c
          call abort
       endif

       p2 = p2 + npatch

       wt_pft            => pvar%wxy_patch(p1:p2)
 
       pftpara           => pvar%pftpara(1:npftpara,p1:p2)
       vegclass          => pvar%vegclass(p1:p2)
       summergreen       => pvar%summergreen(p1:p2)
       raingreen         => pvar%raingreen(p1:p2)
       sla               => pvar%sla(p1:p2)
       lm_sapl           => pvar%lm_sapl(p1:p2)
       sm_sapl           => pvar%sm_sapl(p1:p2)
       hm_sapl           => pvar%hm_sapl(p1:p2)
       rm_sapl           => pvar%rm_sapl(p1:p2)

#ifdef DyN
       cton_soil         => pvar%cton_soil(p1:p2)
#endif
 
       t10min            => pvar%t10min(p1:p2)
       lai_ind           => pvar%lai_ind(p1:p2)
       dphen             => pvar%dphen(p1:p2)
       leafon            => pvar%leafon(p1:p2)
       leafof            => pvar%leafof(p1:p2)
       firelength        => pvar%firelength(p1:p2)
       litterag          => pvar%litterag(p1:p2)
       litterbg          => pvar%litterbg(p1:p2)
       cpool_fast        => pvar%cpool_fast(1:nl_soil,p1:p2)
       cpool_slow        => pvar%cpool_slow(1:nl_soil,p1:p2)
       k_fast_ave        => pvar%k_fast_ave(p1:p2)
       k_slow_ave        => pvar%k_slow_ave(p1:p2)
       litter_decom_ave  => pvar%litter_decom_ave(p1:p2)
       fmicr             => pvar%fmicr(p1:p2)
       nind              => pvar%nind(p1:p2)
       lm_ind            => pvar%lm_ind(p1:p2)
       sm_ind            => pvar%sm_ind(p1:p2)
       hm_ind            => pvar%hm_ind(p1:p2)
       rm_ind            => pvar%rm_ind(p1:p2)
       tmomin20          => pvar%tmomin20(p1:p2)
       agdd0             => pvar%agdd0(p1:p2)
       agdd              => pvar%agdd(p1:p2)
       agdd20            => pvar%agdd20(p1:p2)
       t_mo_min          => pvar%t_mo_min(p1:p2)
       crownarea         => pvar%crownarea(p1:p2)
       htop              => pvar%htop(p1:p2)
       tsai              => pvar%tsai(p1:p2)
       fpcgrid           => pvar%fpcgrid(p1:p2)
       bm_inc            => pvar%bm_inc(p1:p2)
       afmicr            => pvar%afmicr(p1:p2)
       annpsn            => pvar%annpsn(p1:p2)
       annpsnpot         => pvar%annpsnpot(p1:p2)
       tref10            => pvar%tref10(p1:p2)
       tref_sum          => pvar%tref_sum(p1:p2)
       t10               => pvar%t10(1:10,p1:p2)
       assimn10          => pvar%assimn10(p1:p2)
       assimn_sum        => pvar%assimn_sum(p1:p2)
       an10              => pvar%an10(1:10,p1:p2)
       turnover_ind      => pvar%turnover_ind(p1:p2)
       fpc_inc           => pvar%fpc_inc(p1:p2)
       ivt               => pvar%ivt(p1:p2)
       agddtw            => pvar%agddtw(p1:p2)
       ifpre             => pvar%ifpre(p1:p2)
       t_mo              => pvar%t_mo(p1:p2)
       t_mo_sum          => pvar%t_mo_sum(p1:p2)
       anngpp            => pvar%anngpp(p1:p2)
       annfrmf           => pvar%annfrmf(p1:p2)
       annfrms           => pvar%annfrms(p1:p2)
       annfrmr           => pvar%annfrmr(p1:p2)
       annfrg            => pvar%annfrg(p1:p2)
#ifdef IAPDGVM
       afirefrac1        => pvar%afirefrac1(p1:p2)
       nfireg1           => pvar%nfireg1(p1:p2)
#endif


#ifdef DyN
       litter_leaf       => pvar%litter_leaf(p1:p2)
       litter_wood       => pvar%litter_wood(p1:p2)
       litter_root       => pvar%litter_root(p1:p2)
       litter_repr       => pvar%litter_repr(p1:p2)
       litter_leaf_n     => pvar%litter_leaf_n(p1:p2)
       litter_wood_n     => pvar%litter_wood_n(p1:p2)
       litter_root_n     => pvar%litter_root_n(p1:p2)
       litter_repr_n     => pvar%litter_repr_n(p1:p2)
       afcton_leaf       => pvar%afcton_leaf(p1:p2)
       afcton_sap        => pvar%afcton_sap(p1:p2)
       afcton_root       => pvar%afcton_root(p1:p2)
       lm_ind_n          => pvar%lm_ind_n(p1:p2)
       sm_ind_n          => pvar%sm_ind_n(p1:p2)
       hm_ind_n          => pvar%hm_ind_n(p1:p2)
       rm_ind_n          => pvar%rm_ind_n(p1:p2)
       an_up             => pvar%an_up(p1:p2)
       an_stress         => pvar%an_stress(p1:p2)
#endif

       litter2soil_ind   => fLitterSoil(p1:p2)
       litter2atmos_ind  => fLitterAtmos(p1:p2)
 
       IF (itypwat==0) THEN

          g = cgmap(c)
          i = gxmap(g)
          j = gymap(g)

          leafc         = 0._r8
          woodc         = 0._r8
          rootc         = 0._r8
          vegc          = 0._r8
          litc          = 0._r8
          litc_ag       = 0._r8
          litc_bg       = 0._r8
          soic          = 0._r8
          soic_fast     = 0._r8
          soic_slow     = 0._r8
          fveg2litter   = 0._r8
          flitter2soil  = 0._r8
          flitter2atmos = 0._r8
          gpp           = 0._r8
          npp           = 0._r8
          nep           = 0._r8
          nbp           = 0._r8
          ra            = 0._r8
          rh            = 0._r8
          ffirec        = 0._r8
          afirec        = 0._r8
          afiref        = 0._r8
          avegc         = 0._r8
          aestabc       = 0._r8
          anpp          = 0._r8
          amrh          = 0._r8
          alitc_ag      = 0._r8
          alitc_bg      = 0._r8
          asoic_fast    = 0._r8
          asoic_slow    = 0._r8

          CALL lpjave(nyr         ,npatch      ,npftpara      ,pftpara          ,&
                      vegclass    ,bm_inc      ,fpcgrid       ,agdd             ,&
                      t_mo_min    ,tmomin20    ,agdd20        ,afmicr           ,&
                      ifpre       ,lm_ind      ,hm_ind        ,sm_ind           ,&
                      rm_ind      ,litterag    ,litterbg      ,turnover_ind     ,&
                      lai_ind     ,fpc_inc     ,crownarea     ,htop             ,&
                      nind        ,ivt         ,annpsn        ,annpsnpot        ,&
                      sla         ,agddtw      ,firelength                      ,&
#ifdef IAPDGVM
                      afirefrac1  ,nfireg1 ,&                 
#endif
                      lm_sapl     ,sm_sapl     ,hm_sapl     ,rm_sapl   ,wt_pft  ,&
                      wt_col      ,prec365     ,cpool_fast    ,cpool_slow       ,&
                      vegc        ,litc        ,litc_ag       ,litc_bg          ,&
                      soic        ,soic_fast   ,soic_slow     ,&
                      litter2soil_ind, flitter2soil           ,&
                      litter2atmos_ind,flitter2atmos          ,&
                      leafc       ,woodc       ,rootc         ,dz)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! the following weighting should be consistent with rules used in [fluxave] *
        ! carbon densities & fluxes should average over actual land area, not grid area.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Units of kgC/m2

          leafc     = leafc*wt_col*0.001
          woodc     = woodc*wt_col*0.001
          rootc     = rootc*wt_col*0.001
          vegc      = vegc*wt_col*0.001
          litc      = litc*wt_col*0.001
          litc_ag   = litc_ag*wt_col*0.001
          litc_bg   = litc_bg*wt_col*0.001
          soic      = soic*wt_col*0.001
          soic_fast = soic_fast*wt_col*0.001
          soic_slow = soic_slow*wt_col*0.001

        ! Units of kgC/m2/s

          flitter2soil  = flitter2soil*wt_col*0.001
          flitter2atmos = flitter2atmos*wt_col*0.001

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! GPP, NPP, NEP, NBP have already been weighted.

          gpp = fldv%assim(g)*0.012
          ra  = fldv%respc(g)*0.012
          rh  = fldv%fmicr(g)*0.012
          npp = gpp-ra
          nep = gpp-ra-rh
          nbp = gpp-ra-rh-lucflux(i,j)

          ivt_old = ivt
          ifpre_old = ifpre
          wt_pft_old = wt_pft

#ifndef VEGDATA
          IF (newyear) THEN

             gpp_ind  = anngpp
             frmf_ind = annfrmf
             frms_ind = annfrms
             frmr_ind = annfrmr
             frg_ind  = annfrg
   
#ifndef DyN
#ifndef IAPDGVM
             CALL LPJ(nyr         ,npatch      ,npftpara      ,pftpara       ,&
                      vegclass    ,bm_inc      ,fpcgrid       ,agdd          ,&
                      t_mo_min    ,tmomin20    ,agdd20        ,afmicr        ,&
                      ifpre       ,lm_ind      ,hm_ind        ,sm_ind        ,&
                      rm_ind      ,litterag    ,litterbg      ,turnover_ind  ,&
                      lai_ind     ,fpc_inc     ,crownarea     ,htop          ,&
                      nind        ,ivt         ,annpsn        ,annpsnpot     ,&
                      sla         ,agddtw      ,firelength    ,lm_sapl       ,&
                      sm_sapl     ,hm_sapl     ,rm_sapl       ,wt_pft        ,&
                      wt_col      ,prec365     ,afirec        ,afiref        ,&
                      avegc       ,anpp        ,amrh          ,alitc_ag      ,&
                      alitc_bg    ,asoic_fast  ,asoic_slow    ,cpool_fast    ,&
                      cpool_slow  ,aestabc     ,npp_ind       ,dz)
#else
        CALL IAP_DGVM(nyr         ,npatch      ,npftpara      ,pftpara       ,&
                      vegclass    ,bm_inc      ,fpcgrid       ,agdd          ,&
                      t_mo_min    ,tmomin20    ,agdd20        ,afmicr        ,&
                      ifpre       ,lm_ind      ,hm_ind        ,sm_ind        ,&
                      rm_ind      ,litterag    ,litterbg      ,turnover_ind  ,&
                      lai_ind     ,fpc_inc     ,crownarea     ,htop          ,&
                      nind        ,ivt         ,annpsn        ,annpsnpot     ,&
                      sla         ,agddtw      ,firelength    ,afirefrac1    ,&
                      nfireg1     ,lm_sapl     ,sm_sapl       ,hm_sapl       ,&
                      rm_sapl     ,wt_pft      ,wt_col        ,prec365       ,&
                      wliq6mon    ,afirec      ,afiref                       ,&
                      avegc       ,anpp        ,amrh          ,alitc_ag      ,&
                      alitc_bg    ,asoic_fast  ,asoic_slow    ,cpool_fast    ,&
                      cpool_slow  ,aestabc     ,npp_ind       ,dz)

                      !emCO2_gcell ,emCO_gcell    ,emCH4_gcell ,emNHMC_gcell  ,&
               !     emH2_gcell  ,emNOx_gcell   ,emN2O_gcell   ,emPM25_gcell  ,&
               !     emTPM_gcell ,emTC_gcell    ,emOC_gcell    ,emBC_gcell)
#endif
#else
             CALL LPJCN(nyr       ,npatch      ,npftpara      ,pftpara       ,&
                      vegclass    ,bm_inc      ,fpcgrid       ,agdd          ,&
                      t_mo_min    ,tmomin20    ,agdd20        ,afmicr        ,&
                      ifpre       ,lm_ind      ,hm_ind        ,sm_ind        ,&
                      rm_ind      ,litter_leaf ,litter_root   ,litter_wood   ,&
                      litter_repr ,litterag    ,litterbg      ,turnover_ind  ,&
                      lai_ind     ,fpc_inc     ,crownarea     ,htop          ,&
                      nind        ,ivt         ,annpsn        ,annpsnpot     ,&
                      sla         ,agddtw      ,firelength    ,lm_sapl       ,&
                      sm_sapl     ,hm_sapl     ,rm_sapl       ,wt_pft        ,&
                      wt_col      ,prec365     ,an_up         ,afcton_leaf   ,&
                      afcton_root ,afcton_sap  ,litter_leaf_n ,litter_root_n ,&
                      litter_wood_n,litter_repr_n,lm_ind_n    ,sm_ind_n      ,&
                      hm_ind_n    ,rm_ind_n    ,afirec        ,afiref        ,&
                      avegc       ,anpp        ,amrh          ,alitc_ag      ,&
                      alitc_bg    ,asoic_fast  ,asoic_slow    ,cpool_fast    ,&
                      cpool_slow  ,aestabc     ,cton_soil     ,an_up_total   ,&
                      an_stress   ,avegn       ,alitn_ag      ,alitn_bg      ,&
                      asoin       ,an_stress_total            ,npp_ind)
#endif

             CALL lpjreset(c,p1,p2,npatch,ivt,ivt_old,ifpre,ifpre_old,wt_pft,wt_pft_old,&
#if (defined IAPDGVM)
                          wliq6mon,&
#endif                  
                          prec365)
    
             nyr = nyr + 1.

           ! gC/m2/yr => molC/m2/s
             nep_residual(g) = (afirec-aestabc)*wt_col/(12.*dtime)

             nep = nep - (afirec-aestabc)*wt_col*0.001/dtime

             nbp = nbp - (afirec-aestabc)*wt_col*0.001/dtime

           ! gC/m2/yr => kgC/m2/s
             ffirec = afirec*wt_col*0.001/dtime
             fveg2litter = ((alitc_ag+alitc_bg)*wt_col*0.001-litc_ag-litc_bg)/dtime

             litc_ag = alitc_ag*wt_col*0.001
             litc_bg = alitc_bg*wt_col*0.001

          END IF
#endif

        !*Monthly variables (weighted)

          fldv_dgvm%leafc        (g) = leafc
          fldv_dgvm%woodc        (g) = woodc
          fldv_dgvm%rootc        (g) = rootc
          fldv_dgvm%vegc         (g) = vegc
          fldv_dgvm%litc_ag      (g) = litc_ag
          fldv_dgvm%litc_bg      (g) = litc_bg
          fldv_dgvm%litc         (g) = litc
          fldv_dgvm%soic_fast    (g) = soic_fast
          fldv_dgvm%soic_slow    (g) = soic_slow
          fldv_dgvm%soic         (g) = soic
          fldv_dgvm%fveg2litter  (g) = fveg2litter
          fldv_dgvm%flitter2soil (g) = flitter2soil
          fldv_dgvm%flitter2atmos(g) = flitter2atmos
          fldv_dgvm%gpp          (g) = gpp
          fldv_dgvm%npp          (g) = npp
          fldv_dgvm%nep          (g) = nep
          fldv_dgvm%nbp          (g) = nbp
          fldv_dgvm%ra           (g) = ra
          fldv_dgvm%rh           (g) = rh
          fldv_dgvm%ffirec       (g) = ffirec
          fldv_dgvm%pftFrac    (:,g) = wt_pft_old(1:numpft_nat)

#ifdef MYBUG
if(p_master) then
   if(newyear) write(6,*) 'newyear'
   if(c_bug.eq.c) write(6,*) 'wt_patc2', wt_pft_old(1:numpft_nat)
endif
#endif
   
        !*Yearly variables (unweighted)

          fldv_dgvm%bare       (g) = fpcgrid(17)
          fldv_dgvm%afirec     (g) = afirec
          fldv_dgvm%afiref     (g) = afiref
          fldv_dgvm%avegc      (g) = avegc
          fldv_dgvm%aestabc    (g) = aestabc
          fldv_dgvm%anpp       (g) = anpp
          fldv_dgvm%amrh       (g) = amrh
          fldv_dgvm%alitc_ag   (g) = alitc_ag
          fldv_dgvm%alitc_bg   (g) = alitc_bg
          fldv_dgvm%asoic_fast (g) = asoic_fast
          fldv_dgvm%asoic_slow (g) = asoic_slow
          fldv_dgvm%fpcgrid  (:,g) = fpcgrid(1:numpft_nat)
          fldv_dgvm%npp_ind  (:,g) = npp_ind(1:numpft_nat)
          fldv_dgvm%lm_ind   (:,g) = lm_ind(1:numpft_nat)
          fldv_dgvm%sm_ind   (:,g) = sm_ind(1:numpft_nat)
          fldv_dgvm%hm_ind   (:,g) = hm_ind(1:numpft_nat)
          fldv_dgvm%rm_ind   (:,g) = rm_ind(1:numpft_nat)
          fldv_dgvm%crownarea(:,g) = crownarea(1:numpft_nat)
          fldv_dgvm%htop     (:,g) = htop(1:numpft_nat)
          fldv_dgvm%nind     (:,g) = nind(1:numpft_nat)
          fldv_dgvm%lai_ind  (:,g) = lai_ind(1:numpft_nat)
          fldv_dgvm%gpp_ind  (:,g) = gpp_ind(1:numpft_nat)
          fldv_dgvm%frmf_ind (:,g) = frmf_ind(1:numpft_nat)
          fldv_dgvm%frms_ind (:,g) = frms_ind(1:numpft_nat)
          fldv_dgvm%frmr_ind (:,g) = frmr_ind(1:numpft_nat)
          fldv_dgvm%frg_ind  (:,g) = frg_ind(1:numpft_nat)

#ifdef MYBUG
          do p=1,numpft_nat
             if(fpcgrid(p)>0.)then
                print*,'bm',p,npp_ind(p)
                print*,'lai',p,lai_ind(p)    
                print*,'lm',p,lm_ind(p)
                print*,'sm',p,sm_ind(p)
                print*,'hm',p,hm_ind(p)
                print*,'rm',p,rm_ind(p)
             endif
          enddo
#endif
   
#ifdef DyN
          fldv_dgvm%an_up_total    (g) = an_up_total
          fldv_dgvm%an_stress_total(g) = an_stress_total
          fldv_dgvm%avegn          (g) = avegn
          fldv_dgvm%alitn_ag       (g) = alitn_ag
          fldv_dgvm%alitn_bg       (g) = alitn_bg
          fldv_dgvm%asoin          (g) = asoin
          fldv_dgvm%soil_no3       (g) = soil_no3
          fldv_dgvm%soil_nh4       (g) = soil_nh4
          fldv_dgvm%afcton_leaf  (:,g) = afcton_leaf(1:numpft_nat)
          fldv_dgvm%afcton_sap   (:,g) = afcton_sap(1:numpft_nat)
          fldv_dgvm%afcton_root  (:,g) = afcton_root(1:numpft_nat)
   
#ifdef MYBUG
          print*, 'npp_ind',npp_ind
          print*, 'inorganic N',soil_no3,soil_nh4
          print*, 'nitrogen',an_up_total,an_stress_total,avegn,alitn_ag,alitn_bg,asoin,soil_no3,soil_nh4
#endif

#endif
       END IF

       p1 = p2 + 1
 
    END DO  !column
 
 end subroutine LPJDRIVER
