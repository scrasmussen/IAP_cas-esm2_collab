#include <define.h>

module colm_varMod

   use precision
   use paramodel
   implicit none

   type fldv_pft_type
      real(r8), pointer :: taux(:)
      real(r8), pointer :: tauy(:)
      real(r8), pointer :: fsena(:)
      real(r8), pointer :: lfevpa(:)
      real(r8), pointer :: fevpa(:)
      real(r8), pointer :: fsenl(:)
      real(r8), pointer :: fevpl(:)
      real(r8), pointer :: etr(:)
      real(r8), pointer :: fseng(:)
      real(r8), pointer :: fevpg(:)
      real(r8), pointer :: fgrnd(:)
      real(r8), pointer :: sabvsun(:)
      real(r8), pointer :: sabvsha(:)
      real(r8), pointer :: sabg(:)
      real(r8), pointer :: olrg(:)
      real(r8), pointer :: rnet(:)
      real(r8), pointer :: zerr(:)
      real(r8), pointer :: assim(:)
      real(r8), pointer :: respc(:)
      real(r8), pointer :: fmicr(:)
      real(r8), pointer :: tlsun(:)
      real(r8), pointer :: tlsha(:)
      real(r8), pointer :: ldew(:)
      real(r8), pointer :: sigf(:)
      real(r8), pointer :: green(:)
      real(r8), pointer :: lai(:)
      real(r8), pointer :: sai(:)
      real(r8), pointer :: avsdr(:)
      real(r8), pointer :: avsdf(:)
      real(r8), pointer :: anidr(:)
      real(r8), pointer :: anidf(:)
      real(r8), pointer :: sols(:)
      real(r8), pointer :: soll(:)
      real(r8), pointer :: solsd(:)
      real(r8), pointer :: solld(:)
      real(r8), pointer :: solrs(:)
      real(r8), pointer :: solrl(:)
      real(r8), pointer :: solrsd(:)
      real(r8), pointer :: solrld(:)
      real(r8), pointer :: emis(:)
      real(r8), pointer :: z0ma(:)
      real(r8), pointer :: trad(:)
      real(r8), pointer :: ustar(:)
      real(r8), pointer :: tstar(:)
      real(r8), pointer :: qstar(:)
      real(r8), pointer :: zol(:)
      real(r8), pointer :: rib(:)
      real(r8), pointer :: fm(:)
      real(r8), pointer :: fh(:)
      real(r8), pointer :: fq(:)
      real(r8), pointer :: tref(:)
      real(r8), pointer :: qref(:)
      real(r8), pointer :: u10m(:)
      real(r8), pointer :: v10m(:)
      real(r8), pointer :: f10m(:)
      real(r8), pointer :: qsubl(:)
#ifdef DUST 
      real(r8), pointer :: dustemis_bin_1(:)
      real(r8), pointer :: dustemis_bin_2(:)
      real(r8), pointer :: dustemis_bin_3(:)
      real(r8), pointer :: dustemis_bin_4(:)
      real(r8), pointer :: dustemis_total(:)
#endif

   end type fldv_pft_type

   type fldv_col_type
      real(r8), pointer :: xerr(:)
      real(r8), pointer :: rsur(:)
      real(r8), pointer :: rnof(:)
      real(r8), pointer :: tg(:)
      real(r8), pointer :: scv(:)
      real(r8), pointer :: snowdp(:)
      real(r8), pointer :: fsno(:)
      real(r8), pointer :: us(:)
      real(r8), pointer :: vs(:)
      real(r8), pointer :: tm(:)
      real(r8), pointer :: qm(:)
      real(r8), pointer :: prc(:)
      real(r8), pointer :: prl(:)
      real(r8), pointer :: pbot(:)
      real(r8), pointer :: frl(:)
      real(r8), pointer :: solar(:)
      real(r8), pointer :: mrsos(:)
      real(r8), pointer :: mrso(:)
      real(r8), pointer :: mrfso(:)
      real(r8), pointer :: lwsnl(:)
      real(r8), pointer :: snm(:)
      real(r8), pointer :: tsn(:)
      real(r8), pointer :: nsnow(:)
      real(r8), pointer :: tss(:,:)
      real(r8), pointer :: wliq(:,:)
      real(r8), pointer :: wice(:,:)
      real(r8), pointer :: mrlsl(:,:)
      real(r8), pointer :: t_lake(:,:)
      real(r8), pointer :: lake_icefrac(:,:)
      real(r8), pointer :: savedtke1(:)
#if (defined FHNP) && (defined FTF)
 !liruichao add
      real(r8), pointer :: frostdp0(:)
      real(r8), pointer :: frostdp(:)
      real(r8), pointer :: thawdp(:)
      real(r8), pointer :: D_temperature(:)
      real(r8), pointer :: N_time(:)
      real(r8), pointer :: frost_day(:)
      real(r8), pointer :: thaw_day(:)
      !end
#endif
#if (defined FHNP) && (defined GWUSE)
      !Note that GWUSE is still in the testing phase, don't define it.
      real(r8), pointer :: qcharge_ori(:)        !aquifer recharge rate (mm/s) !wanglh 2018/11
      real(r8), pointer :: cgwintake(:)          !column groundwater intake (mm/s) !wanglh 2018/11
#endif
   end type fldv_col_type

   type fldv_dgvm_type
      real(r8), pointer :: leafc(:)
      real(r8), pointer :: woodc(:)
      real(r8), pointer :: rootc(:)
      real(r8), pointer :: vegc(:) 
      real(r8), pointer :: litc_ag(:)
      real(r8), pointer :: litc_bg(:)
      real(r8), pointer :: litc(:)
      real(r8), pointer :: soic_fast(:)
      real(r8), pointer :: soic_slow(:)
      real(r8), pointer :: soic(:)
      real(r8), pointer :: fveg2litter(:)
      real(r8), pointer :: flitter2soil(:)
      real(r8), pointer :: flitter2atmos(:)
      real(r8), pointer :: gpp(:)
      real(r8), pointer :: npp(:)
      real(r8), pointer :: nep(:)
      real(r8), pointer :: nbp(:)
      real(r8), pointer :: ra(:)
      real(r8), pointer :: rh(:)
      real(r8), pointer :: ffirec(:)

      real(r8), pointer :: bare(:)
      real(r8), pointer :: afirec(:)
      real(r8), pointer :: afiref(:)
      real(r8), pointer :: avegc(:)
      real(r8), pointer :: aestabc(:)
      real(r8), pointer :: anpp(:)
      real(r8), pointer :: amrh(:)
      real(r8), pointer :: alitc_ag(:)
      real(r8), pointer :: alitc_bg(:)
      real(r8), pointer :: asoic_fast(:)
      real(r8), pointer :: asoic_slow(:)
      real(r8), pointer :: pftFrac(:,:)
      real(r8), pointer :: fpcgrid(:,:)
      real(r8), pointer :: npp_ind(:,:)
      real(r8), pointer :: lm_ind(:,:)
      real(r8), pointer :: sm_ind(:,:)
      real(r8), pointer :: hm_ind(:,:)
      real(r8), pointer :: rm_ind(:,:)
      real(r8), pointer :: crownarea(:,:)
      real(r8), pointer :: htop(:,:)
      real(r8), pointer :: nind(:,:)
      real(r8), pointer :: lai_ind(:,:)
      real(r8), pointer :: gpp_ind(:,:)
      real(r8), pointer :: frmf_ind(:,:)
      real(r8), pointer :: frms_ind(:,:)
      real(r8), pointer :: frmr_ind(:,:)
      real(r8), pointer :: frg_ind(:,:)
#ifdef DyN
      real(r8), pointer :: afcton_leaf(:,:)
      real(r8), pointer :: afcton_sap(:,:)
      real(r8), pointer :: afcton_root(:,:)
      real(r8), pointer :: an_up_total(:)
      real(r8), pointer :: an_stress_total(:)
      real(r8), pointer :: avegn(:)
      real(r8), pointer :: alitn_ag(:)
      real(r8), pointer :: alitn_bg(:)
      real(r8), pointer :: soil_no3(:)
      real(r8), pointer :: soil_nh4(:)
#endif
   end type fldv_dgvm_type

   type fldv_type
      real(r8), pointer :: taux(:)
      real(r8), pointer :: tauy(:)
      real(r8), pointer :: fsena(:)
      real(r8), pointer :: lfevpa(:)
      real(r8), pointer :: fevpa(:)
      real(r8), pointer :: fsenl(:)
      real(r8), pointer :: fevpl(:)
      real(r8), pointer :: etr(:)
      real(r8), pointer :: fseng(:)
      real(r8), pointer :: fevpg(:)
      real(r8), pointer :: fgrnd(:)
      real(r8), pointer :: sabvsun(:)
      real(r8), pointer :: sabvsha(:)
      real(r8), pointer :: sabg(:)
      real(r8), pointer :: olrg(:)
      real(r8), pointer :: rnet(:)
      real(r8), pointer :: zerr(:)
      real(r8), pointer :: assim(:)
      real(r8), pointer :: respc(:)
      real(r8), pointer :: fmicr(:)
      real(r8), pointer :: tlsun(:)
      real(r8), pointer :: tlsha(:)
      real(r8), pointer :: ldew(:)
      real(r8), pointer :: sigf(:)
      real(r8), pointer :: green(:)
      real(r8), pointer :: lai(:)
      real(r8), pointer :: sai(:)
      real(r8), pointer :: avsdr(:)
      real(r8), pointer :: avsdf(:)
      real(r8), pointer :: anidr(:)
      real(r8), pointer :: anidf(:)
      real(r8), pointer :: sols(:)
      real(r8), pointer :: soll(:)
      real(r8), pointer :: solsd(:)
      real(r8), pointer :: solld(:)
      real(r8), pointer :: solrs(:)
      real(r8), pointer :: solrl(:)
      real(r8), pointer :: solrsd(:)
      real(r8), pointer :: solrld(:)
      real(r8), pointer :: emis(:)
      real(r8), pointer :: z0ma(:)
      real(r8), pointer :: trad(:)
      real(r8), pointer :: ustar(:)
      real(r8), pointer :: tstar(:)
      real(r8), pointer :: qstar(:)
      real(r8), pointer :: zol(:)
      real(r8), pointer :: rib(:)
      real(r8), pointer :: fm(:)
      real(r8), pointer :: fh(:)
      real(r8), pointer :: fq(:)
      real(r8), pointer :: tref(:)
      real(r8), pointer :: qref(:)
      real(r8), pointer :: u10m(:)
      real(r8), pointer :: v10m(:)
      real(r8), pointer :: f10m(:)
      real(r8), pointer :: qsubl(:)
#ifdef DUST
      real(r8), pointer :: dustemis_bin_1(:)
      real(r8), pointer :: dustemis_bin_2(:)
      real(r8), pointer :: dustemis_bin_3(:)
      real(r8), pointer :: dustemis_bin_4(:)
      real(r8), pointer :: dustemis_total(:)
#endif

      real(r8), pointer :: xerr(:)
      real(r8), pointer :: rsur(:)
      real(r8), pointer :: rnof(:)
      real(r8), pointer :: tg(:)
      real(r8), pointer :: scv(:)
      real(r8), pointer :: snowdp(:)
      real(r8), pointer :: fsno(:)
      real(r8), pointer :: us(:)
      real(r8), pointer :: vs(:)
      real(r8), pointer :: tm(:)
      real(r8), pointer :: qm(:)
      real(r8), pointer :: prc(:)
      real(r8), pointer :: prl(:)
      real(r8), pointer :: pbot(:)
      real(r8), pointer :: frl(:)
      real(r8), pointer :: solar(:)
      real(r8), pointer :: mrsos(:)
      real(r8), pointer :: mrso(:)
      real(r8), pointer :: mrfso(:)
      real(r8), pointer :: lwsnl(:)
      real(r8), pointer :: snm(:)
      real(r8), pointer :: tsn(:)
      real(r8), pointer :: nsnow(:)
      real(r8), pointer :: treeFrac(:)
      real(r8), pointer :: shrubFrac(:)
      real(r8), pointer :: grassFrac(:)
      real(r8), pointer :: baresoilFrac(:)
      real(r8), pointer :: residualFrac(:)
      real(r8), pointer :: soilFrac(:)
      real(r8), pointer :: urbanFrac(:)
      real(r8), pointer :: wetlandFrac(:)
      real(r8), pointer :: iceFrac(:)
      real(r8), pointer :: lakeFrac(:)
      real(r8), pointer :: tss(:,:)
      real(r8), pointer :: wliq(:,:)
      real(r8), pointer :: wice(:,:)
      real(r8), pointer :: mrlsl(:,:)
      real(r8), pointer :: t_lake(:,:)
      real(r8), pointer :: lake_icefrac(:,:)
      real(r8), pointer :: savedtke1(:)
#if (defined FHNP) && (defined FTF)
  !liruichao add 
      real(r8), pointer :: frostdp0(:)
      real(r8), pointer :: frostdp(:)
      real(r8), pointer :: thawdp(:)
      real(r8), pointer :: D_temperature(:)
      real(r8),  pointer :: N_time(:)
      real(r8),  pointer :: frost_day(:)
      real(r8),  pointer :: thaw_day(:)
      !end
#endif
#if (defined FHNP) && (defined GWUSE)
      !Note that GWUSE is still in the testing phase, don't define it.
      real(r8), pointer :: gw_uptake(:,:,:)!groundwater uptake from aquifer (mm/s) !wanglh 2018/11
      real(r8), pointer :: noirr_frac(:)   !industry and live use fraction of groundwater  !wanglh 2018/11
      real(r8), pointer :: ggw(:)          !for output wanglh 2018/11
#endif
   end type fldv_type

   type forc_col_type
      real(r8), pointer :: pco2m(:)
      real(r8), pointer :: po2m(:)
      real(r8), pointer :: us(:)
      real(r8), pointer :: vs(:)
      real(r8), pointer :: tm(:)
      real(r8), pointer :: qm(:)
      real(r8), pointer :: prc(:)
      real(r8), pointer :: prl(:)
      real(r8), pointer :: rain(:)
      real(r8), pointer :: snow(:)
      real(r8), pointer :: pbot(:)
      real(r8), pointer :: psrf(:)
      real(r8), pointer :: sols(:)
      real(r8), pointer :: soll(:)
      real(r8), pointer :: solsd(:)
      real(r8), pointer :: solld(:)
      real(r8), pointer :: frl(:)
      real(r8), pointer :: hu(:)
      real(r8), pointer :: ht(:)
      real(r8), pointer :: hq(:)
   end type forc_col_type

   type pvar_type
      integer , pointer :: ivt(:)                  ! land cover type
      real(r8), pointer :: wxy_patch(:)            ! weight of patch
      real(r8), pointer :: z0m(:)                  ! aerodynamic roughness length [m]
      real(r8), pointer :: displa(:)               ! displacement height [m]
      real(r8), pointer :: sqrtdi(:)               ! inverse sqrt of leaf dimension [m**-0.5]
      real(r8), pointer :: effcon(:)               ! quantum efficiency of RuBP regeneration 
      real(r8), pointer :: vmax25(:)               ! maximum carboxylation rate at 25 C at canopy top
      real(r8), pointer :: slti(:)                 ! s3: slope of low temperature inhibition function     
      real(r8), pointer :: hlti(:)                 ! s4: 1/2 point of low temperature inhibition function
      real(r8), pointer :: shti(:)                 ! s1: slope of high temperature inhibition function
      real(r8), pointer :: hhti(:)                 ! s2: 1/2 point of high temperature inhibition function 
      real(r8), pointer :: trda(:)                 ! s5: temperature coefficient in gs-a model
      real(r8), pointer :: trdm(:)                 ! s6: temperature coefficient in gs-a model
      real(r8), pointer :: trop(:)                 ! temperature coefficient in gs-a model          
      real(r8), pointer :: gradm(:)                ! conductance-photosynthesis slope parameter
      real(r8), pointer :: binter(:)               ! conductance-photosynthesis intercep
      real(r8), pointer :: extkn(:)                ! coefficient of leaf nitrogen allocation
      real(r8), pointer :: chil(:)                 ! leaf angle distribution factor
      real(r8), pointer :: refl (:,:,:)            ! leaf reflectance (iw=iband, il=life and dead)
      real(r8), pointer :: refs (:,:,:)            ! stem reflectance (iw=iband, il=life and dead)
      real(r8), pointer :: tranl(:,:,:)            ! leaf transmittance (iw=iband, il=life and dead)
      real(r8), pointer :: trans(:,:,:)            ! stem transmittance (iw=iband, il=life and dead)
      real(r8), pointer :: rootfr(:,:)             ! fraction of roots in each soil layer

      real(r8), pointer :: tlsun(:)                ! sunlit leaf temperature [K]
      real(r8), pointer :: tlsha(:)                ! shaded leaf temperature [K]
      real(r8), pointer :: ldew(:)                 ! depth of water on foliage [mm]
      real(r8), pointer :: fveg(:)                 ! fraction of vegetation cover
      real(r8), pointer :: sigf(:)                 ! fraction of veg cover, excluding snow-covered veg [-]
      real(r8), pointer :: green(:)                ! leaf greenness
      real(r8), pointer :: lai(:)                  ! leaf area index
      real(r8), pointer :: sai(:)                  ! stem area index
      real(r8), pointer :: albg(:,:,:)             ! albedo, ground [-]
      real(r8), pointer :: albv(:,:,:)             ! albedo, vegetation [-]
      real(r8), pointer :: alb (:,:,:)             ! averaged albedo [-]
      real(r8), pointer :: ssun(:,:,:)             ! sunlit canopy absorption for solar radiation (0-1)
      real(r8), pointer :: ssha(:,:,:)             ! shaded canopy absorption for solar radiation (0-1)
      real(r8), pointer :: thermk(:)               ! canopy gap fraction for tir radiation
      real(r8), pointer :: extkb(:)                ! (k, g(mu)/mu) direct solar extinction coefficient
      real(r8), pointer :: extkd(:)                ! diffuse and scattered diffuse PAR extinction coefficient

      ! Additional variables required by regional model

      real(r8), pointer :: trad(:)                 ! radiative temperature of surface [K]
      real(r8), pointer :: tref(:)                 ! 2 m height air temperature [kelvin]
      real(r8), pointer :: qref(:)                 ! 2 m height air specific humidity
      real(r8), pointer :: rst(:)                  ! canopy stomatal resistance (s/m)
      real(r8), pointer :: emis(:)                 ! averaged bulk surface emissivity
      real(r8), pointer :: z0ma(:)                 ! effective roughness [m]
      real(r8), pointer :: zol(:)                  ! dimensionless height (z/L) used in Monin-Obukhov theory
      real(r8), pointer :: rib(:)                  ! bulk Richardson number in surface layer
      real(r8), pointer :: ustar(:)                ! u* in similarity theory [m/s]
      real(r8), pointer :: qstar(:)                ! q* in similarity theory [kg/kg]
      real(r8), pointer :: tstar(:)                ! t* in similarity theory [K]
      real(r8), pointer :: fm(:)                   ! integral of profile function for momentum
      real(r8), pointer :: fh(:)                   ! integral of profile function for heat
      real(r8), pointer :: fq(:)                   ! integral of profile function for moisture
#if(defined DGVM)
      real(r8), pointer :: pftpara(:,:)            ! PFT parameters
      real(r8), pointer :: vegclass(:)             ! 1.tree 2.shrub 3.grass 4.crop
      real(r8), pointer :: summergreen(:)          ! 1 for summergreen, otherwise -1.
      real(r8), pointer :: raingreen(:)            ! 1 for raingreen, otherwise -1.
      real(r8), pointer :: sla(:)                  ! specific leaf area [m2 leaf g-1 carbon]
      real(r8), pointer :: stemdiam(:)             ! canopy height initialization
      real(r8), pointer :: lm_sapl(:)              ! sapling leafmass
      real(r8), pointer :: sm_sapl(:)              ! sapling sapwood mass
      real(r8), pointer :: hm_sapl(:)              ! sapling heartwood mass
      real(r8), pointer :: rm_sapl(:)              ! sapling rootmass
#if(defined DyN)
      real(r8), pointer :: cton_soil(:)            ! soil C:N mass ratio
      real(r8), pointer :: cton_pro(:)             ! C:N mass ratio in production
#endif
      real(r8), pointer :: t10min(:)               ! annual minimum of 10-day running mean (K)
      real(r8), pointer :: lai_ind(:)              ! LAI per individual
      real(r8), pointer :: dphen(:)                ! phenology [0 to 1]
      real(r8), pointer :: leafon(:)               ! leafon days
      real(r8), pointer :: leafof(:)               ! leafoff days
      real(r8), pointer :: firelength(:)           ! fire season in days
#if(defined IAPDGVM)
      real(r8), pointer :: afirefrac1(:)           ! fire fraction
      real(r8), pointer :: nfireg1(:)              ! fire count
#endif
      real(r8), pointer :: litterag(:)             ! above ground litter
      real(r8), pointer :: litterbg(:)             ! below ground litter
      real(r8), pointer :: cpool_fast(:,:)         ! fast carbon pool
      real(r8), pointer :: cpool_slow(:,:)         ! slow carbon pool
      real(r8), pointer :: k_fast_ave(:)           ! decomposition rate
      real(r8), pointer :: k_slow_ave(:)           ! decomposition rate
      real(r8), pointer :: litter_decom_ave(:)     ! decomposition rate
      real(r8), pointer :: fmicr(:)                ! microbial respiration (umol CO2 /m**2 /s)
      real(r8), pointer :: nind(:)                 ! number of individuals (#/m**2)
      real(r8), pointer :: lm_ind(:)               ! individual leaf mass
      real(r8), pointer :: sm_ind(:)               ! individual sapwood mass
      real(r8), pointer :: hm_ind(:)               ! individual heartwood mass
      real(r8), pointer :: rm_ind(:)               ! individual root mass
      real(r8), pointer :: tmomin20(:)             ! 20-yr running mean of tmomin
      real(r8), pointer :: agdd0(:)                ! growing dgree days above 0
      real(r8), pointer :: agdd(:)                 ! growing dgree days above 5
      real(r8), pointer :: agdd20(:)               ! 20-yr running mean of agdd
      real(r8), pointer :: agddtw(:)               ! accumulated growing degree days above twmax
      real(r8), pointer :: t_mo_min(:)             ! annual min of t_mo (Kelvin)
      real(r8), pointer :: crownarea(:)            ! area that each individual tree takes up (m^2)
      real(r8), pointer :: htop(:)                 ! canopy top
      real(r8), pointer :: tsai(:)                 ! one-sided stem area index, no burying by snow
      real(r8), pointer :: fpcgrid(:)              ! foliar projective cover over grid cell
      real(r8), pointer :: bm_inc(:)               ! biomass increment
      real(r8), pointer :: afmicr(:)               ! annual microbial respiration
      real(r8), pointer :: annpsn(:)               ! annual photosynthesis (umol CO2 /m**2)
      real(r8), pointer :: annpsnpot(:)            ! annual potential photosynthesis (same units)
      real(r8), pointer :: tref10(:)               ! 10-day averaged temperature at 2m
      real(r8), pointer :: tref_sum(:)             ! sum of tref of current day
      real(r8), pointer :: t10(:,:)                ! array to record 10 day tref
      real(r8), pointer :: assimn10(:)             ! 10-day averaged assimilation rate
      real(r8), pointer :: assimn_sum(:)           ! sum of assimn of current day
      real(r8), pointer :: an10(:,:)               ! array to record 10 day assimn
      real(r8), pointer :: turnover_ind(:)         ! individual turnover biomass
      real(r8), pointer :: fpc_inc(:)              ! fpc increase
      real(r8), pointer :: ifpre(:)                ! whether PFT present in patch (1 or -1) 
      real(r8), pointer :: t_mo(:)                 ! 30-day mean temperature of 2m
      real(r8), pointer :: t_mo_sum(:)             ! 30-day accumulated temperature of 2m
      real(r8), pointer :: anngpp(:)               ! annual gpp
      real(r8), pointer :: annfrmf(:)              ! annual frmf
      real(r8), pointer :: annfrms(:)              ! annual frms
      real(r8), pointer :: annfrmr(:)              ! annual frmr
      real(r8), pointer :: annfrg(:)               ! annual frg
#if(defined DyN)
      real(r8), pointer :: litter_leaf(:)          ! leaf-derived litter for PFT on modelled area basis (gC/m2)
      real(r8), pointer :: litter_wood(:)          ! heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
      real(r8), pointer :: litter_root(:)          ! fine root-derived litter for PFT on modelled area basis(gC/m2)
      real(r8), pointer :: litter_repr(:)          ! litter derived from allocation to reproduction for PFT on modelled
      real(r8), pointer :: litter_leaf_n(:)        ! leaf-derived N litter for PFT on modelled area basis (gN/m2)
      real(r8), pointer :: litter_wood_n(:)        ! heart&sapwood-derived N litter for PFT on modelled area basis(gN/m2)
      real(r8), pointer :: litter_root_n(:)        ! fine root-derived N litter for PFT on modelled area basis (gN/m2)
      real(r8), pointer :: litter_repr_n(:)        ! litter derived from allocation to reproduction N for PFT on modelled area basis (gN/m2)
      real(r8), pointer :: afcton_leaf(:)          ! annual floating leaf C:N ratio
      real(r8), pointer :: afcton_sap(:)           ! annual floating sapwood C:N ratio
      real(r8), pointer :: afcton_root(:)          ! annual floating root C:N ratio
      real(r8), pointer :: lm_ind_n(:)             ! individual leaf nitrogen mass
      real(r8), pointer :: sm_ind_n(:)             ! individual sapwood nitrogen mass
      real(r8), pointer :: hm_ind_n(:)             ! individual heartwood nitrogen mass
      real(r8), pointer :: rm_ind_n(:)             ! individual root nitrogen mass
                                                   ! gN/m2 veget'd area for each pft
      real(r8), pointer :: an_up(:)                ! annual plant nitrogen uptake(gN/m2 vegt'd area)
      real(r8), pointer :: an_stress(:)            ! annual plant nitrogen stress(-)
#endif
#endif
   end type pvar_type

   type cvar_type
      integer , pointer :: itypwat(:)              ! land water type
      real(r8), pointer :: dlat(:)                 ! latitude in radians
      real(r8), pointer :: dlon(:)                 ! longitude in radians
      real(r8), pointer :: wxy_column(:)           ! weight of column
      real(r8), pointer :: albsat(:,:)             ! wet soil albedo [vis/nir] for different coloured soils [-]
      real(r8), pointer :: albdry(:,:)             ! dry soil albedo [vis/nir] for different coloured soils [-]
      real(r8), pointer :: csol(:,:)               ! heat capacity of soil solids [J/(m3 K)]
      real(r8), pointer :: porsl(:,:)              ! fraction of soil that is voids [-]
      real(r8), pointer :: phi0(:,:)               ! minimum soil suction [mm]
      real(r8), pointer :: bsw(:,:)                ! clapp and hornbereger "b" parameter [-]
      real(r8), pointer :: dkmg(:,:)               ! thermal conductivity of soil minerals [W/m-K]
      real(r8), pointer :: dksatu(:,:)             ! thermal conductivity of saturated soil [W/m-K]
      real(r8), pointer :: dkdry(:,:)              ! thermal conductivity for dry soil  [W/(m-K)]
      real(r8), pointer :: hksati(:,:)             ! hydraulic conductivity at saturation [mm h2o/s]

      real(r8), pointer :: z(:,:)                  ! node depth [m]
      real(r8), pointer :: dz(:,:)                 ! layer thickness [m]
      real(r8), pointer :: tss(:,:)                ! soil temperature [K]
      real(r8), pointer :: wliq(:,:)               ! liquid water in layers [kg/m2]
      real(r8), pointer :: wice(:,:)               ! ice lens in layers [kg/m2]
      real(r8), pointer :: tg(:)                   ! ground surface temperature [K]
      real(r8), pointer :: sag(:)                  ! non dimensional snow age [-]
      real(r8), pointer :: scv(:)                  ! snow cover, water equivalent [mm]
      real(r8), pointer :: snowdp(:)               ! snow depth [meter]
      real(r8), pointer :: fsno(:)                 ! fraction of snow cover on ground
      real(r8), pointer :: coszen(:)               ! cosine of solar zenith angle
      real(r8), pointer :: lakedepth(:)            ! lake depth [m] added by Nan Wei
      real(r8), pointer :: dz_lake(:,:)            ! lake thickness [m] added by Nan Wei
      real(r8), pointer :: t_lake(:,:)             ! lake layer temperature [K] added by Nan Wei
      real(r8), pointer :: lake_icefrac(:,:)       ! lake mass fraction of lake layer that is frozen added by Nan Wei
      real(r8), pointer :: savedtke1(:)            ! added by Nan Wei

#if(defined DGVM)
      real(r8), pointer :: nday(:)                 ! counting the model days
      real(r8), pointer :: nyr(:)                  ! counting the model years
      real(r8), pointer :: prec365(:)              ! yearly running mean of precipitation(mm/s)
#if(defined IAPDGVM)
      real(r8), pointer :: wliq6mon(:)             ! IAPDGVM liquid water 6 mons for first 3 layers
#endif

#if(defined DyN)
      real(r8), pointer :: soil_no3(:)
      real(r8), pointer :: soil_no2(:)
      real(r8), pointer :: soil_no(:)
      real(r8), pointer :: soil_n2o(:)
      real(r8), pointer :: soil_n2(:)
      real(r8), pointer :: soil_nh4(:)
#endif
#endif

#if (defined FHNP) && (defined FTF)
!liruichao add
      real(r8), pointer :: frostdp0(:)             !initial frost depth       
      real(r8), pointer :: frostdp(:)              !frost depth
      real(r8), pointer :: thawdp(:)               !thaw deppth
      real(r8), pointer :: D_temperature(:)        !frost or thaw index
      real(r8), pointer :: N_time(:)               !step counter
      real(r8), pointer :: frost_day(:)            !frost days
      real(r8), pointer :: thaw_day(:)             !thaw days
      !end
#endif

#ifdef DUST 
      integer , pointer :: dustsource(:)           ! index for potential dust source
      integer , pointer :: isltyp(:)               ! dominant soil type
      real(r8), pointer :: soil_top_cat(:,:)       ! fraction for 12-categories soil
      real(r8), pointer :: mvegcov(:,:)            ! read_in monthly vegetation cover [unit: %, i.e. 0-100]
#endif
   end type cvar_type

   type(pvar_type)      :: pvar
   type(cvar_type)      :: cvar

   type(forc_col_type)  :: forc

   type(fldv_pft_type)  :: fldv_pft
   type(fldv_col_type)  :: fldv_col
   type(fldv_dgvm_type) :: fldv_dgvm
   type(fldv_type)      :: fldv

!
! basic model grid info
!
   integer :: lon_points                      ! number of longitude points on model grid
   integer :: lat_points                      ! number of latitude points on model grid

   integer :: numgrid                         ! local grids number
   integer :: numcolumn                       ! local columns number
   integer :: numpatch                        ! local patches number

   integer :: numgrid_glob                    ! global grids number
   integer :: numcolumn_glob                  ! global columns number
   integer :: numpatch_glob                   ! global patches number

! 
! land model grid location info
!
   integer , pointer :: numlon(:)             ! longitude points for each latitude strip
   real(r8), pointer :: lats(:)               ! grid cell latitude, southern edge (degrees)
   real(r8), pointer :: lonw(:,:)             ! grid cell longitude, western edge (degrees)
   real(r8), pointer :: area(:,:)             ! grid cell area (km**2)
   real(r8), pointer :: latixy(:,:)           ! latitude of grid cell (degrees)
   real(r8), pointer :: longxy(:,:)           ! longitude of grid cell (degrees)
   real(r8), pointer :: landarea              ! total land area for all gridcells (km^2)
!
! fractional land and mask
!
   real(r8), pointer :: landfrac(:,:)         ! fractional land
   integer , pointer :: landmask(:,:)         ! land mask: 1 = land. 0 = ocean
!
! patch and grid info
!
   integer,  pointer :: itypwat_glob(:)

   integer , pointer :: ixy_patch_glob(:)     ! longitude index of patch
   integer , pointer :: jxy_patch_glob(:)     ! latitude index of patch
   real(r8), pointer :: wxy_patch_glob(:)     ! weight of patch

   integer , pointer :: ixy_column_glob(:)    ! longitude index of patch
   integer , pointer :: jxy_column_glob(:)    ! latitude index of patch
   real(r8), pointer :: wxy_column_glob(:)    ! weight of column
!
! model variables
!
   real(r8), pointer :: oro(:)                ! ocean(0)/seaice(2)/ flag

#if(defined DGVM)
   integer,  pointer :: numcolumn_lat(:)      ! number of columnes of grids at lon. strip
   integer,  pointer :: numpatch_lat(:)       ! number of patches of grids at lon. strip
#if(defined IAPDGVM)
   real(r8), pointer :: iglf(:,:,:)                ! fire 
#endif
   real(r8), pointer :: nep_residual(:)       ! annual residual NEP = acfire + aestabc
   logical           :: lnep_adjust = .false.
   real(r8), pointer :: fLitterSoil(:)        ! cflux_litter_soil
   real(r8), pointer :: fLitterAtmos(:)       ! cflux_litter_atmos

   real(r8), pointer :: Isf_pft(:,:)
   real(r8), pointer :: Iss_pft(:,:)
   real(r8), pointer :: Ksf_pft(:,:)
   real(r8), pointer :: Kss_pft(:,:)
#endif

!
! temporary surface data info
!
   real(r8), allocatable :: dlat_glob    (:) ! latitude in radians
   real(r8), allocatable :: dlon_glob    (:) ! longitude in radians
   real(r8), allocatable :: rockdep_glob (:) ! depth to bedrock
   real(r8), allocatable :: sand_glob  (:,:) ! percent sand
   real(r8), allocatable :: clay_glob  (:,:) ! percent clay
   real(r8), allocatable :: soc_glob   (:,:) ! soil organic carbon density [kg/m^3]
   integer,  allocatable :: isc_glob     (:) ! color classes for soil albedos
   integer,  allocatable :: ivt_glob     (:) ! land cover type

#ifdef DUST
   integer,  allocatable :: dustsource_glob   (:)   ! index for potential dust source
   integer,  allocatable :: isltyp_glob       (:)   ! dominant soil type
   real(r8), allocatable :: soil_top_cat_glob (:,:) ! fraction for 12-categories soil
   real(r8), allocatable :: mvegcov_glob      (:,:) ! read_in monthly vegetation cover [unit: %, i.e. 0-100]
#endif
!
! tunable parameters
!
   real(r8)  ftune(nftune)              ! clm tunable constants

!
! forcing variables
!
   real(r8), pointer :: tair   (:,:)
   real(r8), pointer :: qair   (:,:)
   real(r8), pointer :: pres   (:,:)
   real(r8), pointer :: rainc  (:,:)
   real(r8), pointer :: rainl  (:,:)
   real(r8), pointer :: windu  (:,:)
   real(r8), pointer :: windv  (:,:)
   real(r8), pointer :: dswrf  (:,:)
   real(r8), pointer :: dlwrf  (:,:)
   real(r8), pointer :: tair_z (:,:)
   real(r8), pointer :: qair_z (:,:)
   real(r8), pointer :: wind_z (:,:)

   real(r8), pointer :: lai(:,:)
   real(r8), pointer :: sai(:,:)
   real(r8), pointer :: fveg(:,:)
   real(r8), pointer :: green(:,:)

#if(defined VEGDATA)
   real(r8), pointer :: mlai(:,:)            ! monthly LAI  (/12,numpatch/)
   real(r8), pointer :: msai(:,:)            ! monthly SAI  (/12,numpatch/)
   real(r8), pointer :: mhtop(:,:)           ! monthly HTOP (/12,numpatch/)
#endif

   interface colm_srfvar_alloc
      module procedure colm_srfvar_alloc
   end interface

   interface colm_srfvar_free
      module procedure colm_srfvar_free
   end interface

   interface colm_gridvar_alloc
      module procedure colm_gridvar_alloc
   end interface 

   interface colm_gridvar_free
      module procedure colm_gridvar_free
   end interface 

   interface colm_subgridvar_alloc
      module procedure colm_subgridvar_alloc
   end interface 

   interface colm_subgridvar_free
      module procedure colm_subgridvar_free
   end interface 

   interface colm_pcgvar_alloc
      module procedure colm_pcgvar_alloc
   end interface 

   interface colm_pcgvar_free
      module procedure colm_pcgvar_free
   end interface 

CONTAINS

   subroutine colm_gridvar_alloc()

      allocate (numlon                (lat_points))
      allocate (lats                (lat_points+1))
      allocate (lonw     (lon_points+1,lat_points))
      allocate (area       (lon_points,lat_points))
      allocate (latixy     (lon_points,lat_points))
      allocate (longxy     (lon_points,lat_points))
      allocate (landmask   (lon_points,lat_points))
      allocate (landfrac   (lon_points,lat_points))

   end subroutine colm_gridvar_alloc

   subroutine colm_gridvar_free()

      deallocate (numlon      )
      deallocate (lats        )
      deallocate (lonw        )
      deallocate (area        )
      deallocate (latixy      )
      deallocate (longxy      )
      deallocate (landmask    )
      deallocate (landfrac    )

   end subroutine colm_gridvar_free

   subroutine colm_subgridvar_alloc()

      allocate (itypwat_glob   (numcolumn_glob))
      allocate (ixy_column_glob(numcolumn_glob))
      allocate (jxy_column_glob(numcolumn_glob))
      allocate (wxy_column_glob(numcolumn_glob))

      allocate (ixy_patch_glob  (numpatch_glob))
      allocate (jxy_patch_glob  (numpatch_glob))
      allocate (wxy_patch_glob  (numpatch_glob))

   end subroutine colm_subgridvar_alloc

   subroutine colm_subgridvar_free()

      deallocate (   itypwat_glob)
      deallocate (ixy_column_glob)
      deallocate (jxy_column_glob)
      deallocate (wxy_column_glob)

      deallocate ( ixy_patch_glob)
      deallocate ( jxy_patch_glob)
      deallocate ( wxy_patch_glob)

   end subroutine colm_subgridvar_free

   subroutine colm_srfvar_alloc()

      allocate (dlat_glob         (numcolumn_glob))
      allocate (dlon_glob         (numcolumn_glob))
      allocate (rockdep_glob      (numcolumn_glob))
      allocate (sand_glob (nl_soil,numcolumn_glob))
      allocate (clay_glob (nl_soil,numcolumn_glob))
      allocate (soc_glob  (nl_soil,numcolumn_glob))
      allocate (isc_glob          (numcolumn_glob))
      allocate (ivt_glob           (numpatch_glob))

#ifdef DUST 
      allocate (dustsource_glob     (numcolumn_glob))
      allocate (isltyp_glob         (numcolumn_glob))
      allocate (soil_top_cat_glob(12,numcolumn_glob))
      allocate (mvegcov_glob     (12,numcolumn_glob))
#endif

   end subroutine colm_srfvar_alloc

   subroutine colm_srfvar_free()

      deallocate (dlat_glob   )
      deallocate (dlon_glob   )
      deallocate (rockdep_glob)
      deallocate (sand_glob   )
      deallocate (clay_glob   )
      deallocate (soc_glob    )
      deallocate (isc_glob    )
      deallocate (ivt_glob    )

#ifdef DUST 
      deallocate (dustsource_glob  )
      deallocate (isltyp_glob      )
      deallocate (soil_top_cat_glob)
      deallocate (mvegcov_glob     )
#endif

   end subroutine colm_srfvar_free

   subroutine colm_pcgvar_alloc

      allocate (oro                (numcolumn))

#ifdef VEGDATA
      allocate (mlai             (12,numpatch))
      allocate (msai             (12,numpatch))
      allocate (mhtop            (12,numpatch))
#endif
#ifdef IAPDGVM
allocate (iglf(lon_points,lat_points,365*8))      !IAPDGVM
#endif

#ifdef DGVM
      allocate (nep_residual         (numgrid))
      nep_residual(:) = 0.

      allocate (fLitterSoil         (numpatch))
      allocate (fLitterAtmos        (numpatch))

      allocate (Isf_pft     (nl_soil,numpatch))
      allocate (Iss_pft     (nl_soil,numpatch))
      allocate (Ksf_pft     (nl_soil,numpatch))
      allocate (Kss_pft     (nl_soil,numpatch))

      Isf_pft(:,:) = 0.
      Iss_pft(:,:) = 0.
      Ksf_pft(:,:) = 0.
      Kss_pft(:,:) = 0.
#endif

      call pvar_alloc
      call cvar_alloc
      call forc_alloc
      call fldv_alloc

   end subroutine colm_pcgvar_alloc

   subroutine colm_pcgvar_free()

      deallocate (oro         )

#ifdef VEGDATA
      deallocate (mlai        )
      deallocate (msai        )
      deallocate (mhtop       )
#endif
#ifdef IAPDGVM
      deallocate (iglf        )
#endif


#ifdef DGVM
      deallocate (nep_residual)
      deallocate (fLitterSoil )
      deallocate (fLitterAtmos)

      deallocate (Isf_pft     )
      deallocate (Iss_pft     )
      deallocate (Ksf_pft     )
      deallocate (Kss_pft     )
#endif

    ! call pvar_free
    ! call cvar_free
      call forc_free
      call fldv_free

   end subroutine colm_pcgvar_free

   subroutine pvar_alloc()
      allocate (pvar%ivt                   (numpatch))
      allocate (pvar%wxy_patch             (numpatch))
      allocate (pvar%z0m                   (numpatch))
      allocate (pvar%displa                (numpatch))
      allocate (pvar%sqrtdi                (numpatch))
      allocate (pvar%effcon                (numpatch))
      allocate (pvar%vmax25                (numpatch))
      allocate (pvar%slti                  (numpatch))
      allocate (pvar%hlti                  (numpatch))
      allocate (pvar%shti                  (numpatch))
      allocate (pvar%hhti                  (numpatch))
      allocate (pvar%trda                  (numpatch))
      allocate (pvar%trdm                  (numpatch))
      allocate (pvar%trop                  (numpatch))
      allocate (pvar%gradm                 (numpatch))
      allocate (pvar%binter                (numpatch))
      allocate (pvar%extkn                 (numpatch))
      allocate (pvar%chil                  (numpatch))
      allocate (pvar%refl              (2,2,numpatch))
      allocate (pvar%refs              (2,2,numpatch))
      allocate (pvar%tranl             (2,2,numpatch))
      allocate (pvar%trans             (2,2,numpatch))
      allocate (pvar%rootfr        (nl_soil,numpatch))
      allocate (pvar%tlsun                 (numpatch))
      allocate (pvar%tlsha                 (numpatch))
      allocate (pvar%ldew                  (numpatch))
      allocate (pvar%fveg                  (numpatch))
      allocate (pvar%sigf                  (numpatch))
      allocate (pvar%green                 (numpatch))
      allocate (pvar%lai                   (numpatch))
      allocate (pvar%sai                   (numpatch))
      allocate (pvar%albg              (2,2,numpatch))
      allocate (pvar%albv              (2,2,numpatch))
      allocate (pvar%alb               (2,2,numpatch))
      allocate (pvar%ssun              (2,2,numpatch))
      allocate (pvar%ssha              (2,2,numpatch))
      allocate (pvar%thermk                (numpatch))
      allocate (pvar%extkb                 (numpatch))
      allocate (pvar%extkd                 (numpatch))

    ! Additional variables required by reginal model
      allocate (pvar%trad                  (numpatch))
      allocate (pvar%tref                  (numpatch))
      allocate (pvar%qref                  (numpatch))
      allocate (pvar%rst                   (numpatch))
      allocate (pvar%emis                  (numpatch))
      allocate (pvar%z0ma                  (numpatch))
      allocate (pvar%zol                   (numpatch))
      allocate (pvar%rib                   (numpatch))
      allocate (pvar%ustar                 (numpatch))
      allocate (pvar%qstar                 (numpatch))
      allocate (pvar%tstar                 (numpatch))
      allocate (pvar%fm                    (numpatch))
      allocate (pvar%fh                    (numpatch))
      allocate (pvar%fq                    (numpatch))
#if(defined DGVM)
      allocate (pvar%pftpara      (npftpara,numpatch))
      allocate (pvar%vegclass              (numpatch))
      allocate (pvar%summergreen           (numpatch))
      allocate (pvar%raingreen             (numpatch))
      allocate (pvar%sla                   (numpatch))
      allocate (pvar%stemdiam              (numpatch))
      allocate (pvar%lm_sapl               (numpatch))
      allocate (pvar%sm_sapl               (numpatch))
      allocate (pvar%hm_sapl               (numpatch))
      allocate (pvar%rm_sapl               (numpatch))
#if(defined DyN)
      allocate (pvar%cton_soil             (numpatch))
      allocate (pvar%cton_pro              (numpatch))
#endif
      allocate (pvar%t10min                (numpatch))
      allocate (pvar%lai_ind               (numpatch))
      allocate (pvar%dphen                 (numpatch))
      allocate (pvar%leafon                (numpatch))
      allocate (pvar%leafof                (numpatch))
      allocate (pvar%firelength            (numpatch))
#if(defined IAPDGVM)
      allocate (pvar%afirefrac1            (numpatch))
      allocate (pvar%nfireg1               (numpatch))
#endif
      allocate (pvar%litterag              (numpatch))
      allocate (pvar%litterbg              (numpatch))
      allocate (pvar%cpool_fast    (nl_soil,numpatch))
      allocate (pvar%cpool_slow    (nl_soil,numpatch))
      allocate (pvar%k_fast_ave            (numpatch))
      allocate (pvar%k_slow_ave            (numpatch))
      allocate (pvar%litter_decom_ave      (numpatch))
      allocate (pvar%fmicr                 (numpatch))
      allocate (pvar%nind                  (numpatch))
      allocate (pvar%lm_ind                (numpatch))
      allocate (pvar%sm_ind                (numpatch))
      allocate (pvar%hm_ind                (numpatch))
      allocate (pvar%rm_ind                (numpatch))
      allocate (pvar%tmomin20              (numpatch))
      allocate (pvar%agdd0                 (numpatch))
      allocate (pvar%agdd                  (numpatch))
      allocate (pvar%agdd20                (numpatch))
      allocate (pvar%agddtw                (numpatch))
      allocate (pvar%t_mo                  (numpatch))
      allocate (pvar%t_mo_sum              (numpatch))
      allocate (pvar%t_mo_min              (numpatch))
      allocate (pvar%crownarea             (numpatch))
      allocate (pvar%htop                  (numpatch))
      allocate (pvar%tsai                  (numpatch))
      allocate (pvar%fpcgrid               (numpatch))
      allocate (pvar%bm_inc                (numpatch))
      allocate (pvar%afmicr                (numpatch))
      allocate (pvar%annpsn                (numpatch))
      allocate (pvar%annpsnpot             (numpatch))
      allocate (pvar%tref10                (numpatch))
      allocate (pvar%tref_sum              (numpatch))
      allocate (pvar%t10                (10,numpatch))
      allocate (pvar%assimn10              (numpatch))
      allocate (pvar%assimn_sum            (numpatch))
      allocate (pvar%an10               (10,numpatch))
      allocate (pvar%turnover_ind          (numpatch))
      allocate (pvar%fpc_inc               (numpatch))
      allocate (pvar%ivt                   (numpatch))
      allocate (pvar%ifpre                 (numpatch))
      allocate (pvar%anngpp                (numpatch))
      allocate (pvar%annfrmf               (numpatch))
      allocate (pvar%annfrms               (numpatch))
      allocate (pvar%annfrmr               (numpatch))
      allocate (pvar%annfrg                (numpatch))
#if(defined DyN)
      allocate (pvar%litter_leaf           (numpatch))
      allocate (pvar%litter_wood           (numpatch))
      allocate (pvar%litter_root           (numpatch))
      allocate (pvar%litter_repr           (numpatch))
      allocate (pvar%litter_leaf_n         (numpatch))
      allocate (pvar%litter_wood_n         (numpatch))
      allocate (pvar%litter_root_n         (numpatch))
      allocate (pvar%litter_repr_n         (numpatch))
      allocate (pvar%afcton_leaf           (numpatch))
      allocate (pvar%afcton_sap            (numpatch))
      allocate (pvar%afcton_root           (numpatch))
      allocate (pvar%lm_ind_n              (numpatch))
      allocate (pvar%sm_ind_n              (numpatch))
      allocate (pvar%hm_ind_n              (numpatch))
      allocate (pvar%rm_ind_n              (numpatch))
      allocate (pvar%an_up                 (numpatch))
      allocate (pvar%an_stress             (numpatch))
#endif
#endif
   end subroutine pvar_alloc

   subroutine cvar_alloc()
      allocate (cvar%dlat                 (numcolumn))
      allocate (cvar%dlon                 (numcolumn))
      allocate (cvar%itypwat              (numcolumn))
      allocate (cvar%wxy_column           (numcolumn))
      allocate (cvar%albsat             (2,numcolumn))
      allocate (cvar%albdry             (2,numcolumn))
      allocate (cvar%csol         (nl_soil,numcolumn))
      allocate (cvar%porsl        (nl_soil,numcolumn))
      allocate (cvar%phi0         (nl_soil,numcolumn))
      allocate (cvar%bsw          (nl_soil,numcolumn))
      allocate (cvar%dkmg         (nl_soil,numcolumn))
      allocate (cvar%dksatu       (nl_soil,numcolumn))
      allocate (cvar%dkdry        (nl_soil,numcolumn))
      allocate (cvar%hksati       (nl_soil,numcolumn))

      allocate (cvar%z   (maxsnl+1:nl_soil,numcolumn))
      allocate (cvar%dz  (maxsnl+1:nl_soil,numcolumn))
      allocate (cvar%tss (maxsnl+1:nl_soil,numcolumn))
      allocate (cvar%wliq(maxsnl+1:nl_soil,numcolumn))
      allocate (cvar%wice(maxsnl+1:nl_soil,numcolumn))
      allocate (cvar%tg                   (numcolumn))
      allocate (cvar%sag                  (numcolumn))
      allocate (cvar%scv                  (numcolumn))
      allocate (cvar%snowdp               (numcolumn))
      allocate (cvar%fsno                 (numcolumn))
      allocate (cvar%coszen               (numcolumn))
      allocate (cvar%lakedepth            (numcolumn))
      allocate (cvar%dz_lake      (nl_lake,numcolumn))
      allocate (cvar%t_lake       (nl_lake,numcolumn))
      allocate (cvar%lake_icefrac (nl_lake,numcolumn))
      allocate (cvar%savedtke1            (numcolumn))

#if(defined DGVM)
      allocate (cvar%nday                 (numcolumn))
      allocate (cvar%nyr                  (numcolumn))
      allocate (cvar%prec365              (numcolumn))
#if(defined IAPDGVM)
      allocate (cvar%wliq6mon             (numcolumn))
#endif

#if(defined DyN)
      allocate (cvar%soil_no3             (numcolumn))
      allocate (cvar%soil_no2             (numcolumn))
      allocate (cvar%soil_no              (numcolumn))
      allocate (cvar%soil_n2o             (numcolumn))
      allocate (cvar%soil_n2              (numcolumn))
      allocate (cvar%soil_nh4             (numcolumn))
#endif
#endif

#if (defined FHNP) && (defined FTF)
 !liruichao add
      allocate(cvar%frostdp0      (numcolumn))
      allocate(cvar%frostdp       (numcolumn))
      allocate(cvar%thawdp        (numcolumn))
      allocate(cvar%D_temperature (numcolumn))
      allocate(cvar%N_time        (numcolumn))
      allocate(cvar%frost_day     (numcolumn))
      allocate(cvar%thaw_day      (numcolumn))
      !end
#endif

#ifdef DUST 
      allocate (cvar%dustsource            (numcolumn))
      allocate (cvar%isltyp                (numcolumn))
      allocate (cvar%soil_top_cat       (12,numcolumn))
      allocate (cvar%mvegcov            (12,numcolumn))
#endif

   end subroutine cvar_alloc

   subroutine fldv_alloc()
      allocate (fldv_pft%taux   (numpatch))
      allocate (fldv_pft%tauy   (numpatch))
      allocate (fldv_pft%fsena  (numpatch))
      allocate (fldv_pft%lfevpa (numpatch))
      allocate (fldv_pft%fevpa  (numpatch))
      allocate (fldv_pft%fsenl  (numpatch))
      allocate (fldv_pft%fevpl  (numpatch))
      allocate (fldv_pft%etr    (numpatch))
      allocate (fldv_pft%fseng  (numpatch))
      allocate (fldv_pft%fevpg  (numpatch))
      allocate (fldv_pft%fgrnd  (numpatch))
      allocate (fldv_pft%sabvsun(numpatch))
      allocate (fldv_pft%sabvsha(numpatch))
      allocate (fldv_pft%sabg   (numpatch))
      allocate (fldv_pft%olrg   (numpatch))
      allocate (fldv_pft%rnet   (numpatch))
      allocate (fldv_pft%zerr   (numpatch))
      allocate (fldv_pft%assim  (numpatch))
      allocate (fldv_pft%respc  (numpatch))
      allocate (fldv_pft%fmicr  (numpatch))
      allocate (fldv_pft%tlsun  (numpatch))
      allocate (fldv_pft%tlsha  (numpatch))
      allocate (fldv_pft%ldew   (numpatch))
      allocate (fldv_pft%sigf   (numpatch))
      allocate (fldv_pft%green  (numpatch))
      allocate (fldv_pft%lai    (numpatch))
      allocate (fldv_pft%sai    (numpatch))
      allocate (fldv_pft%avsdr  (numpatch))
      allocate (fldv_pft%avsdf  (numpatch))
      allocate (fldv_pft%anidr  (numpatch))
      allocate (fldv_pft%anidf  (numpatch))
      allocate (fldv_pft%sols   (numpatch))
      allocate (fldv_pft%soll   (numpatch))
      allocate (fldv_pft%solsd  (numpatch))
      allocate (fldv_pft%solld  (numpatch))
      allocate (fldv_pft%solrs  (numpatch))
      allocate (fldv_pft%solrl  (numpatch))
      allocate (fldv_pft%solrsd (numpatch))
      allocate (fldv_pft%solrld (numpatch))
      allocate (fldv_pft%emis   (numpatch))
      allocate (fldv_pft%z0ma   (numpatch))
      allocate (fldv_pft%trad   (numpatch))
      allocate (fldv_pft%ustar  (numpatch))
      allocate (fldv_pft%tstar  (numpatch))
      allocate (fldv_pft%qstar  (numpatch))
      allocate (fldv_pft%zol    (numpatch))
      allocate (fldv_pft%rib    (numpatch))
      allocate (fldv_pft%fm     (numpatch))
      allocate (fldv_pft%fh     (numpatch))
      allocate (fldv_pft%fq     (numpatch))
      allocate (fldv_pft%tref   (numpatch))
      allocate (fldv_pft%qref   (numpatch))
      allocate (fldv_pft%u10m   (numpatch))
      allocate (fldv_pft%v10m   (numpatch))
      allocate (fldv_pft%f10m   (numpatch))
      allocate (fldv_pft%qsubl  (numpatch))
#ifdef DUST 
      allocate (fldv_pft%dustemis_bin_1  (numpatch))
      allocate (fldv_pft%dustemis_bin_2  (numpatch))
      allocate (fldv_pft%dustemis_bin_3  (numpatch))
      allocate (fldv_pft%dustemis_bin_4  (numpatch))
      allocate (fldv_pft%dustemis_total  (numpatch))
#endif

      allocate (fldv_col%xerr  (numcolumn))
      allocate (fldv_col%rsur  (numcolumn))
      allocate (fldv_col%rnof  (numcolumn))
      allocate (fldv_col%tg    (numcolumn))
      allocate (fldv_col%scv   (numcolumn))
      allocate (fldv_col%snowdp(numcolumn))
      allocate (fldv_col%fsno  (numcolumn))
      allocate (fldv_col%us    (numcolumn))
      allocate (fldv_col%vs    (numcolumn))
      allocate (fldv_col%tm    (numcolumn))
      allocate (fldv_col%qm    (numcolumn))
      allocate (fldv_col%prc   (numcolumn))
      allocate (fldv_col%prl   (numcolumn))
      allocate (fldv_col%pbot  (numcolumn))
      allocate (fldv_col%frl   (numcolumn))
      allocate (fldv_col%solar (numcolumn))
      allocate (fldv_col%mrsos (numcolumn))
      allocate (fldv_col%mrso  (numcolumn))
      allocate (fldv_col%mrfso (numcolumn))
      allocate (fldv_col%lwsnl (numcolumn))
      allocate (fldv_col%snm   (numcolumn))
      allocate (fldv_col%tsn   (numcolumn))
      allocate (fldv_col%nsnow (numcolumn))
      allocate (fldv_col%tss   (nl_soil,numcolumn))
      allocate (fldv_col%wliq  (nl_soil,numcolumn))
      allocate (fldv_col%wice  (nl_soil,numcolumn))
      allocate (fldv_col%mrlsl (nl_soil,numcolumn))
      allocate (fldv_col%t_lake      (nl_lake,numcolumn))
      allocate (fldv_col%lake_icefrac(nl_lake,numcolumn))
      allocate (fldv_col%savedtke1           (numcolumn))

#ifdef DGVM
      allocate (fldv_dgvm%leafc                 (numgrid))
      allocate (fldv_dgvm%woodc                 (numgrid))
      allocate (fldv_dgvm%rootc                 (numgrid))
      allocate (fldv_dgvm%vegc                  (numgrid))
      allocate (fldv_dgvm%litc_ag               (numgrid))
      allocate (fldv_dgvm%litc_bg               (numgrid))
      allocate (fldv_dgvm%litc                  (numgrid))
      allocate (fldv_dgvm%soic_fast             (numgrid))
      allocate (fldv_dgvm%soic_slow             (numgrid))
      allocate (fldv_dgvm%soic                  (numgrid))
      allocate (fldv_dgvm%fveg2litter           (numgrid))
      allocate (fldv_dgvm%flitter2soil          (numgrid))
      allocate (fldv_dgvm%flitter2atmos         (numgrid))
      allocate (fldv_dgvm%gpp                   (numgrid))
      allocate (fldv_dgvm%npp                   (numgrid))
      allocate (fldv_dgvm%nep                   (numgrid))
      allocate (fldv_dgvm%nbp                   (numgrid))
      allocate (fldv_dgvm%ra                    (numgrid))
      allocate (fldv_dgvm%rh                    (numgrid))
      allocate (fldv_dgvm%ffirec                (numgrid))

! annual variables
      allocate (fldv_dgvm%bare                  (numgrid))
      allocate (fldv_dgvm%afirec                (numgrid))
      allocate (fldv_dgvm%afiref                (numgrid))
      allocate (fldv_dgvm%avegc                 (numgrid))
      allocate (fldv_dgvm%aestabc               (numgrid))
      allocate (fldv_dgvm%anpp                  (numgrid))
      allocate (fldv_dgvm%amrh                  (numgrid))
      allocate (fldv_dgvm%alitc_ag              (numgrid))
      allocate (fldv_dgvm%alitc_bg              (numgrid))
      allocate (fldv_dgvm%asoic_fast            (numgrid))
      allocate (fldv_dgvm%asoic_slow            (numgrid))
      allocate (fldv_dgvm%pftFrac    (numpft_nat,numgrid))
      allocate (fldv_dgvm%fpcgrid    (numpft_nat,numgrid))
      allocate (fldv_dgvm%npp_ind    (numpft_nat,numgrid))
      allocate (fldv_dgvm%lm_ind     (numpft_nat,numgrid))
      allocate (fldv_dgvm%sm_ind     (numpft_nat,numgrid))
      allocate (fldv_dgvm%hm_ind     (numpft_nat,numgrid))
      allocate (fldv_dgvm%rm_ind     (numpft_nat,numgrid))
      allocate (fldv_dgvm%crownarea  (numpft_nat,numgrid))
      allocate (fldv_dgvm%htop       (numpft_nat,numgrid))
      allocate (fldv_dgvm%nind       (numpft_nat,numgrid))
      allocate (fldv_dgvm%lai_ind    (numpft_nat,numgrid))
      allocate (fldv_dgvm%gpp_ind    (numpft_nat,numgrid))
      allocate (fldv_dgvm%frmf_ind   (numpft_nat,numgrid))
      allocate (fldv_dgvm%frms_ind   (numpft_nat,numgrid))
      allocate (fldv_dgvm%frmr_ind   (numpft_nat,numgrid))
      allocate (fldv_dgvm%frg_ind    (numpft_nat,numgrid))

#ifdef DyN
      allocate (fldv_dgvm%afcton_leaf(numpft_nat,numgrid))
      allocate (fldv_dgvm%afcton_sap (numpft_nat,numgrid))
      allocate (fldv_dgvm%afcton_root(numpft_nat,numgrid))
      allocate (fldv_dgvm%an_up_total           (numgrid))
      allocate (fldv_dgvm%an_stress_total       (numgrid))
      allocate (fldv_dgvm%avegn                 (numgrid))
      allocate (fldv_dgvm%alitn_ag              (numgrid))
      allocate (fldv_dgvm%alitn_bg              (numgrid))
      allocate (fldv_dgvm%asoin                 (numgrid))
      allocate (fldv_dgvm%soil_no3              (numgrid))
      allocate (fldv_dgvm%soil_nh4              (numgrid))
#endif
#endif

! average from pft level
      allocate (fldv%taux      (numgrid))
      allocate (fldv%tauy      (numgrid))
      allocate (fldv%fsena     (numgrid))
      allocate (fldv%lfevpa    (numgrid))
      allocate (fldv%fevpa     (numgrid))
      allocate (fldv%fsenl     (numgrid))
      allocate (fldv%fevpl     (numgrid))
      allocate (fldv%etr       (numgrid))
      allocate (fldv%fseng     (numgrid))
      allocate (fldv%fevpg     (numgrid))
      allocate (fldv%fgrnd     (numgrid))
      allocate (fldv%sabvsun   (numgrid))
      allocate (fldv%sabvsha   (numgrid))
      allocate (fldv%sabg      (numgrid))
      allocate (fldv%olrg      (numgrid))
      allocate (fldv%rnet      (numgrid))
      allocate (fldv%zerr      (numgrid))
      allocate (fldv%assim     (numgrid))
      allocate (fldv%respc     (numgrid))
      allocate (fldv%fmicr     (numgrid))
      allocate (fldv%tlsun     (numgrid))
      allocate (fldv%tlsha     (numgrid))
      allocate (fldv%ldew      (numgrid))
      allocate (fldv%sigf      (numgrid))
      allocate (fldv%green     (numgrid))
      allocate (fldv%lai       (numgrid))
      allocate (fldv%sai       (numgrid))
      allocate (fldv%avsdr     (numgrid))
      allocate (fldv%avsdf     (numgrid))
      allocate (fldv%anidr     (numgrid))
      allocate (fldv%anidf     (numgrid))
      allocate (fldv%sols      (numgrid))
      allocate (fldv%soll      (numgrid))
      allocate (fldv%solsd     (numgrid))
      allocate (fldv%solld     (numgrid))
      allocate (fldv%solrs     (numgrid))
      allocate (fldv%solrl     (numgrid))
      allocate (fldv%solrsd    (numgrid))
      allocate (fldv%solrld    (numgrid))
      allocate (fldv%emis      (numgrid))
      allocate (fldv%z0ma      (numgrid))
      allocate (fldv%trad      (numgrid))
      allocate (fldv%ustar     (numgrid))
      allocate (fldv%tstar     (numgrid))
      allocate (fldv%qstar     (numgrid))
      allocate (fldv%zol       (numgrid))
      allocate (fldv%rib       (numgrid))
      allocate (fldv%fm        (numgrid))
      allocate (fldv%fh        (numgrid))
      allocate (fldv%fq        (numgrid))
      allocate (fldv%tref      (numgrid))
      allocate (fldv%qref      (numgrid))
      allocate (fldv%u10m      (numgrid))
      allocate (fldv%v10m      (numgrid))
      allocate (fldv%f10m      (numgrid))
      allocate (fldv%qsubl     (numgrid))
#ifdef DUST
      allocate (fldv%dustemis_bin_1      (numgrid))
      allocate (fldv%dustemis_bin_2      (numgrid))
      allocate (fldv%dustemis_bin_3      (numgrid))
      allocate (fldv%dustemis_bin_4      (numgrid))
      allocate (fldv%dustemis_total      (numgrid))
#endif

! average from column level
      allocate (fldv%xerr      (numgrid))
      allocate (fldv%rsur      (numgrid))
      allocate (fldv%rnof      (numgrid))
      allocate (fldv%tg        (numgrid))
      allocate (fldv%scv       (numgrid))
      allocate (fldv%snowdp    (numgrid))
      allocate (fldv%fsno      (numgrid))
      allocate (fldv%us        (numgrid))
      allocate (fldv%vs        (numgrid))
      allocate (fldv%tm        (numgrid))
      allocate (fldv%qm        (numgrid))
      allocate (fldv%prc       (numgrid))
      allocate (fldv%prl       (numgrid))
      allocate (fldv%pbot      (numgrid))
      allocate (fldv%frl       (numgrid))
      allocate (fldv%solar     (numgrid))
      allocate (fldv%mrsos     (numgrid))
      allocate (fldv%mrso      (numgrid))
      allocate (fldv%mrfso     (numgrid))
      allocate (fldv%lwsnl     (numgrid))
      allocate (fldv%snm       (numgrid))
      allocate (fldv%tsn       (numgrid))
      allocate (fldv%nsnow     (numgrid))

      allocate (fldv%treeFrac     (numgrid))
      allocate (fldv%shrubFrac    (numgrid))
      allocate (fldv%grassFrac    (numgrid))
      allocate (fldv%baresoilFrac (numgrid))
      allocate (fldv%residualFrac (numgrid))
      allocate (fldv%soilFrac     (numgrid))
      allocate (fldv%urbanFrac    (numgrid))
      allocate (fldv%wetlandFrac  (numgrid))
      allocate (fldv%iceFrac      (numgrid))
      allocate (fldv%lakeFrac     (numgrid))

      allocate (fldv%tss  (nl_soil,numgrid))
      allocate (fldv%wliq (nl_soil,numgrid))
      allocate (fldv%wice (nl_soil,numgrid))
      allocate (fldv%mrlsl(nl_soil,numgrid))
      allocate (fldv%t_lake      (nl_lake,numgrid))
      allocate (fldv%lake_icefrac(nl_lake,numgrid))
      allocate (fldv%savedtke1           (numgrid))

#if (defined FHNP) && (defined GWUSE)
     ! Note that GWUSE is still in the testing phase, don't define it.
      allocate (fldv_col%qcharge_ori (numcolumn))   !wanglh 2018/11
      allocate (fldv_col%cgwintake (numcolumn))     !wanglh 2018/11
      allocate (fldv%gw_uptake     (numgrid,12,46)) !wanglh 2018/11
      allocate (fldv%noirr_frac     (numgrid))      !wanglh 2018/11
      allocate (fldv%ggw     (numgrid))             !wanglh 2018/11
#endif

#if (defined FHNP) && (defined FTF)
      allocate(fldv_col%frostdp0      (numcolumn))
      allocate(fldv_col%frostdp       (numcolumn))
      allocate(fldv_col%thawdp        (numcolumn))
      allocate(fldv_col%D_temperature (numcolumn))
      allocate(fldv_col%N_time        (numcolumn))
      allocate(fldv_col%frost_day     (numcolumn))
      allocate(fldv_col%thaw_day      (numcolumn))
      allocate(fldv%frostdp0      (numgrid))
      allocate(fldv%frostdp       (numgrid))
      allocate(fldv%thawdp        (numgrid))
      allocate(fldv%D_temperature (numgrid))
      allocate(fldv%N_time        (numgrid))
      allocate(fldv%frost_day     (numgrid))
      allocate(fldv%thaw_day      (numgrid))
#endif
   end subroutine fldv_alloc

   subroutine fldv_dgvm_init()
#ifdef DGVM
      fldv_dgvm%leafc                 = 0.
      fldv_dgvm%woodc                 = 0.
      fldv_dgvm%rootc                 = 0.
      fldv_dgvm%vegc                  = 0.
      fldv_dgvm%litc_ag               = 0.
      fldv_dgvm%litc_bg               = 0.
      fldv_dgvm%litc                  = 0.
      fldv_dgvm%soic_fast             = 0.
      fldv_dgvm%soic_slow             = 0.
      fldv_dgvm%soic                  = 0.
      fldv_dgvm%fveg2litter           = 0.
      fldv_dgvm%flitter2soil          = 0.
      fldv_dgvm%flitter2atmos         = 0.
      fldv_dgvm%gpp                   = 0.
      fldv_dgvm%npp                   = 0.
      fldv_dgvm%nep                   = 0.
      fldv_dgvm%nbp                   = 0.
      fldv_dgvm%ra                    = 0.
      fldv_dgvm%rh                    = 0.
      fldv_dgvm%ffirec                = 0.

! annual variables
      fldv_dgvm%bare                  = 0.
      fldv_dgvm%afirec                = 0.
      fldv_dgvm%afiref                = 0.
      fldv_dgvm%avegc                 = 0.
      fldv_dgvm%aestabc               = 0.
      fldv_dgvm%anpp                  = 0.
      fldv_dgvm%amrh                  = 0.
      fldv_dgvm%alitc_ag              = 0.
      fldv_dgvm%alitc_bg              = 0.
      fldv_dgvm%asoic_fast            = 0.
      fldv_dgvm%asoic_slow            = 0.
      fldv_dgvm%pftFrac               = 0.
      fldv_dgvm%fpcgrid               = 0.
      fldv_dgvm%npp_ind               = 0.
      fldv_dgvm%lm_ind                = 0.
      fldv_dgvm%sm_ind                = 0.
      fldv_dgvm%hm_ind                = 0.
      fldv_dgvm%rm_ind                = 0.
      fldv_dgvm%crownarea             = 0.
      fldv_dgvm%htop                  = 0.
      fldv_dgvm%nind                  = 0.
      fldv_dgvm%lai_ind               = 0.
      fldv_dgvm%gpp_ind               = 0.
      fldv_dgvm%frmf_ind              = 0.
      fldv_dgvm%frms_ind              = 0.
      fldv_dgvm%frmr_ind              = 0.
      fldv_dgvm%frg_ind               = 0.

#ifdef DyN
      fldv_dgvm%afcton_leaf           = 0.
      fldv_dgvm%afcton_sap            = 0.
      fldv_dgvm%afcton_root           = 0.
      fldv_dgvm%an_up_total           = 0.
      fldv_dgvm%an_stress_total       = 0.
      fldv_dgvm%avegn                 = 0.
      fldv_dgvm%alitn_ag              = 0.
      fldv_dgvm%alitn_bg              = 0.
      fldv_dgvm%asoin                 = 0.
      fldv_dgvm%soil_no3              = 0.
      fldv_dgvm%soil_nh4              = 0.
#endif
#endif
   end subroutine fldv_dgvm_init

   subroutine fldv_init()
    ! Average from pft level
      fldv%taux      = 0.
      fldv%tauy      = 0.
      fldv%fsena     = 0.
      fldv%lfevpa    = 0.
      fldv%fevpa     = 0.
      fldv%fsenl     = 0.
      fldv%fevpl     = 0.
      fldv%etr       = 0.
      fldv%fseng     = 0.
      fldv%fevpg     = 0.
      fldv%fgrnd     = 0.
      fldv%sabvsun   = 0.
      fldv%sabvsha   = 0.
      fldv%sabg      = 0.
      fldv%olrg      = 0.
      fldv%rnet      = 0.
      fldv%zerr      = 0.
      fldv%assim     = 0.
      fldv%respc     = 0.
      fldv%fmicr     = 0.
      fldv%tlsun     = 0.
      fldv%tlsha     = 0.
      fldv%ldew      = 0.
      fldv%sigf      = 0.
      fldv%green     = 0.
      fldv%lai       = 0.
      fldv%sai       = 0.
      fldv%avsdr     = 0.
      fldv%avsdf     = 0.
      fldv%anidr     = 0.
      fldv%anidf     = 0.
      fldv%sols      = 0.
      fldv%soll      = 0.
      fldv%solsd     = 0.
      fldv%solld     = 0.
      fldv%solrs     = 0.
      fldv%solrl     = 0.
      fldv%solrsd    = 0.
      fldv%solrld    = 0.
      fldv%emis      = 0.
      fldv%z0ma      = 0.
      fldv%trad      = 0.
      fldv%ustar     = 0.
      fldv%tstar     = 0.
      fldv%qstar     = 0.
      fldv%zol       = 0.
      fldv%rib       = 0.
      fldv%fm        = 0.
      fldv%fh        = 0.
      fldv%fq        = 0.
      fldv%tref      = 0.
      fldv%qref      = 0.
      fldv%u10m      = 0.
      fldv%v10m      = 0.
      fldv%f10m      = 0.
      fldv%qsubl     = 0.
#ifdef DUST
      fldv%dustemis_bin_1 =0.
      fldv%dustemis_bin_2 =0.
      fldv%dustemis_bin_3 =0.
      fldv%dustemis_bin_4 =0.
      fldv%dustemis_total =0.
#endif

    ! Average from column level
      fldv%xerr      = 0.
      fldv%rsur      = 0.
      fldv%rnof      = 0.
      fldv%tg        = 0.
      fldv%scv       = 0.
      fldv%snowdp    = 0.
      fldv%fsno      = 0.
      fldv%us        = 0.
      fldv%vs        = 0.
      fldv%tm        = 0.
      fldv%qm        = 0.
      fldv%prc       = 0.
      fldv%prl       = 0.
      fldv%pbot      = 0.
      fldv%frl       = 0.
      fldv%solar     = 0.
      fldv%mrsos     = 0.
      fldv%mrso      = 0.
      fldv%mrfso     = 0.
      fldv%lwsnl     = 0.
      fldv%snm       = 0.
      fldv%tsn       = 0.
      fldv%nsnow     = 0.

      fldv%treeFrac     = 0.
      fldv%shrubFrac    = 0.
      fldv%grassFrac    = 0.
      fldv%baresoilFrac = 0.
      fldv%residualFrac = 0.
      fldv%soilFrac     = 0.
      fldv%urbanFrac    = 0.
      fldv%wetlandFrac  = 0.
      fldv%iceFrac      = 0.
      fldv%lakeFrac     = 0.

      fldv%tss          = 0.
      fldv%wliq         = 0.
      fldv%wice         = 0.
      fldv%mrlsl        = 0.
      fldv%t_lake       = 0.
      fldv%lake_icefrac = 0.
      fldv%savedtke1    = 0.
#if (defined FHNP) && (defined GWUSE)
     ! Note that GWUSE is still in the testing phase, don't define it.
      fldv_col%qcharge_ori =0.    !wanglh 2018/11
      fldv_col%cgwintake   =0.    !wanglh 2018/11
      fldv%gw_uptake       =0.    !wanglh 2018/11
      fldv%noirr_frac      =0.    !wanglh 2018/11
      fldv%ggw             =0.    !wanglh 2018/11
#endif

#if (defined FHNP) && (defined FTF)
      fldv%frostdp0      = 0.
      fldv%frostdp       = 0.
      fldv%thawdp        = 0.
      fldv%D_temperature = 0.
      fldv%N_time        = 0.
      fldv%frost_day     = 0.
      fldv%thaw_day      = 0.
#endif
   end subroutine fldv_init

   subroutine forc_alloc()
      allocate (forc%pco2m (numcolumn))
      allocate (forc%po2m  (numcolumn))
      allocate (forc%us    (numcolumn))
      allocate (forc%vs    (numcolumn))
      allocate (forc%tm    (numcolumn))
      allocate (forc%qm    (numcolumn))
      allocate (forc%prc   (numcolumn))
      allocate (forc%prl   (numcolumn))
      allocate (forc%rain  (numcolumn))
      allocate (forc%snow  (numcolumn))
      allocate (forc%pbot  (numcolumn))
      allocate (forc%psrf  (numcolumn))
      allocate (forc%sols  (numcolumn))
      allocate (forc%soll  (numcolumn))
      allocate (forc%solsd (numcolumn))
      allocate (forc%solld (numcolumn))
      allocate (forc%frl   (numcolumn))
      allocate (forc%hu    (numcolumn))
      allocate (forc%ht    (numcolumn))
      allocate (forc%hq    (numcolumn))

#ifndef COUP_CSM 
      allocate (tair   (lon_points,lat_points))
      allocate (qair   (lon_points,lat_points))
      allocate (pres   (lon_points,lat_points))
      allocate (rainc  (lon_points,lat_points))
      allocate (rainl  (lon_points,lat_points))
      allocate (windu  (lon_points,lat_points))
      allocate (windv  (lon_points,lat_points))
      allocate (dswrf  (lon_points,lat_points))
      allocate (dlwrf  (lon_points,lat_points))
      allocate (tair_z (lon_points,lat_points))
      allocate (qair_z (lon_points,lat_points))
      allocate (wind_z (lon_points,lat_points))
#endif

      allocate (lai    (lon_points,lat_points))
      allocate (sai    (lon_points,lat_points))
      allocate (fveg   (lon_points,lat_points))
      allocate (green  (lon_points,lat_points))
   end subroutine forc_alloc

   subroutine forc_free()
      deallocate (forc%pco2m)
      deallocate (forc%po2m )
      deallocate (forc%us   )
      deallocate (forc%vs   )
      deallocate (forc%tm   )
      deallocate (forc%qm   )
      deallocate (forc%prc  )
      deallocate (forc%prl  )
      deallocate (forc%rain )
      deallocate (forc%snow )
      deallocate (forc%pbot )
      deallocate (forc%psrf )
      deallocate (forc%sols )
      deallocate (forc%soll )
      deallocate (forc%solsd)
      deallocate (forc%solld)
      deallocate (forc%frl  )
      deallocate (forc%hu   )
      deallocate (forc%ht   )
      deallocate (forc%hq   )

#ifndef COUP_CSM 
      deallocate (tair      )
      deallocate (qair      )
      deallocate (pres      )
      deallocate (rainc     )
      deallocate (rainl     )
      deallocate (windu     )
      deallocate (windv     )
      deallocate (dswrf     )
      deallocate (dlwrf     )
      deallocate (tair_z    )
      deallocate (qair_z    )
      deallocate (wind_z    )
#endif

      deallocate (lai       )
      deallocate (sai       )
      deallocate (fveg      )
      deallocate (green     )
   end subroutine forc_free

   subroutine fldv_free()
      deallocate (fldv_pft%taux   )
      deallocate (fldv_pft%tauy   )
      deallocate (fldv_pft%fsena  )
      deallocate (fldv_pft%lfevpa )
      deallocate (fldv_pft%fevpa  )
      deallocate (fldv_pft%fsenl  )
      deallocate (fldv_pft%fevpl  )
      deallocate (fldv_pft%etr    )
      deallocate (fldv_pft%fseng  )
      deallocate (fldv_pft%fevpg  )
      deallocate (fldv_pft%fgrnd  )
      deallocate (fldv_pft%sabvsun)
      deallocate (fldv_pft%sabvsha)
      deallocate (fldv_pft%sabg   )
      deallocate (fldv_pft%olrg   )
      deallocate (fldv_pft%rnet   )
      deallocate (fldv_pft%zerr   )
      deallocate (fldv_pft%assim  )
      deallocate (fldv_pft%respc  )
      deallocate (fldv_pft%fmicr  )
      deallocate (fldv_pft%tlsun  )
      deallocate (fldv_pft%tlsha  )
      deallocate (fldv_pft%ldew   )
      deallocate (fldv_pft%sigf   )
      deallocate (fldv_pft%green  )
      deallocate (fldv_pft%lai    )
      deallocate (fldv_pft%sai    )
      deallocate (fldv_pft%avsdr  )
      deallocate (fldv_pft%avsdf  )
      deallocate (fldv_pft%anidr  )
      deallocate (fldv_pft%anidf  )
      deallocate (fldv_pft%sols   )
      deallocate (fldv_pft%soll   )
      deallocate (fldv_pft%solsd  )
      deallocate (fldv_pft%solld  )
      deallocate (fldv_pft%solrs  )
      deallocate (fldv_pft%solrl  )
      deallocate (fldv_pft%solrsd )
      deallocate (fldv_pft%solrld )
      deallocate (fldv_pft%emis   )
      deallocate (fldv_pft%z0ma   )
      deallocate (fldv_pft%trad   )
      deallocate (fldv_pft%ustar  )
      deallocate (fldv_pft%tstar  )
      deallocate (fldv_pft%qstar  )
      deallocate (fldv_pft%zol    )
      deallocate (fldv_pft%rib    )
      deallocate (fldv_pft%fm     )
      deallocate (fldv_pft%fh     )
      deallocate (fldv_pft%fq     )
      deallocate (fldv_pft%tref   )
      deallocate (fldv_pft%qref   )
      deallocate (fldv_pft%u10m   )
      deallocate (fldv_pft%v10m   )
      deallocate (fldv_pft%f10m   )
      deallocate (fldv_pft%qsubl  )
#ifdef DUST
      deallocate (fldv_pft%dustemis_bin_1  )
      deallocate (fldv_pft%dustemis_bin_2  )
      deallocate (fldv_pft%dustemis_bin_3  )
      deallocate (fldv_pft%dustemis_bin_4  )
      deallocate (fldv_pft%dustemis_total  )
#endif

      deallocate (fldv_col%xerr  )
      deallocate (fldv_col%rsur  )
      deallocate (fldv_col%rnof  )
      deallocate (fldv_col%tg    )
      deallocate (fldv_col%scv   )
      deallocate (fldv_col%snowdp)
      deallocate (fldv_col%fsno  )
      deallocate (fldv_col%us    )
      deallocate (fldv_col%vs    )
      deallocate (fldv_col%tm    )
      deallocate (fldv_col%qm    )
      deallocate (fldv_col%prc   )
      deallocate (fldv_col%prl   )
      deallocate (fldv_col%pbot  )
      deallocate (fldv_col%frl   )
      deallocate (fldv_col%solar )
      deallocate (fldv_col%mrsos )
      deallocate (fldv_col%mrso  )
      deallocate (fldv_col%mrfso )
      deallocate (fldv_col%lwsnl )
      deallocate (fldv_col%snm   )
      deallocate (fldv_col%tsn   )
      deallocate (fldv_col%nsnow )
      deallocate (fldv_col%tss   )
      deallocate (fldv_col%wliq  )
      deallocate (fldv_col%wice  )
      deallocate (fldv_col%mrlsl )
      deallocate (fldv_col%t_lake)
      deallocate (fldv_col%lake_icefrac)
      deallocate (fldv_col%savedtke1   )

#ifdef DGVM
      deallocate (fldv_dgvm%leafc          )
      deallocate (fldv_dgvm%woodc          )
      deallocate (fldv_dgvm%rootc          )
      deallocate (fldv_dgvm%vegc           )
      deallocate (fldv_dgvm%litc_ag        )
      deallocate (fldv_dgvm%litc_bg        )
      deallocate (fldv_dgvm%litc           )
      deallocate (fldv_dgvm%soic_fast      )
      deallocate (fldv_dgvm%soic_slow      )
      deallocate (fldv_dgvm%soic           )
      deallocate (fldv_dgvm%fveg2litter    )
      deallocate (fldv_dgvm%flitter2soil   )
      deallocate (fldv_dgvm%flitter2atmos  )
      deallocate (fldv_dgvm%gpp            )
      deallocate (fldv_dgvm%npp            )
      deallocate (fldv_dgvm%nep            )
      deallocate (fldv_dgvm%nbp            )
      deallocate (fldv_dgvm%ra             )
      deallocate (fldv_dgvm%rh             )
      deallocate (fldv_dgvm%ffirec         )

! annual variables
      deallocate (fldv_dgvm%bare           )
      deallocate (fldv_dgvm%afirec         )
      deallocate (fldv_dgvm%afiref         )
      deallocate (fldv_dgvm%avegc          )
      deallocate (fldv_dgvm%aestabc        )
      deallocate (fldv_dgvm%anpp           )
      deallocate (fldv_dgvm%amrh           )
      deallocate (fldv_dgvm%alitc_ag       )
      deallocate (fldv_dgvm%alitc_bg       )
      deallocate (fldv_dgvm%asoic_fast     )
      deallocate (fldv_dgvm%asoic_slow     )
      deallocate (fldv_dgvm%pftFrac        )
      deallocate (fldv_dgvm%fpcgrid        )
      deallocate (fldv_dgvm%npp_ind        )
      deallocate (fldv_dgvm%lm_ind         )
      deallocate (fldv_dgvm%sm_ind         )
      deallocate (fldv_dgvm%hm_ind         )
      deallocate (fldv_dgvm%rm_ind         )
      deallocate (fldv_dgvm%crownarea      )
      deallocate (fldv_dgvm%htop           )
      deallocate (fldv_dgvm%nind           )
      deallocate (fldv_dgvm%lai_ind        )
      deallocate (fldv_dgvm%gpp_ind        )
      deallocate (fldv_dgvm%frmf_ind       )
      deallocate (fldv_dgvm%frms_ind       )
      deallocate (fldv_dgvm%frmr_ind       )
      deallocate (fldv_dgvm%frg_ind        )

#ifdef DyN
      deallocate (fldv_dgvm%afcton_leaf    )
      deallocate (fldv_dgvm%afcton_sap     )
      deallocate (fldv_dgvm%afcton_root    )
      deallocate (fldv_dgvm%an_up_total    )
      deallocate (fldv_dgvm%an_stress_total)
      deallocate (fldv_dgvm%avegn          )
      deallocate (fldv_dgvm%alitn_ag       )
      deallocate (fldv_dgvm%alitn_bg       )
      deallocate (fldv_dgvm%asoin          )
      deallocate (fldv_dgvm%soil_no3       )
      deallocate (fldv_dgvm%soil_nh4       )
#endif
#endif

! average from pft level
      deallocate (fldv%taux      )
      deallocate (fldv%tauy      )
      deallocate (fldv%fsena     )
      deallocate (fldv%lfevpa    )
      deallocate (fldv%fevpa     )
      deallocate (fldv%fsenl     )
      deallocate (fldv%fevpl     )
      deallocate (fldv%etr       )
      deallocate (fldv%fseng     )
      deallocate (fldv%fevpg     )
      deallocate (fldv%fgrnd     )
      deallocate (fldv%sabvsun   )
      deallocate (fldv%sabvsha   )
      deallocate (fldv%sabg      )
      deallocate (fldv%olrg      )
      deallocate (fldv%rnet      )
      deallocate (fldv%zerr      )
      deallocate (fldv%assim     )
      deallocate (fldv%respc     )
      deallocate (fldv%fmicr     )
      deallocate (fldv%tlsun     )
      deallocate (fldv%tlsha     )
      deallocate (fldv%ldew      )
      deallocate (fldv%sigf      )
      deallocate (fldv%green     )
      deallocate (fldv%lai       )
      deallocate (fldv%sai       )
      deallocate (fldv%avsdr     )
      deallocate (fldv%avsdf     )
      deallocate (fldv%anidr     )
      deallocate (fldv%anidf     )
      deallocate (fldv%sols      )
      deallocate (fldv%soll      )
      deallocate (fldv%solsd     )
      deallocate (fldv%solld     )
      deallocate (fldv%solrs     )
      deallocate (fldv%solrl     )
      deallocate (fldv%solrsd    )
      deallocate (fldv%solrld    )
      deallocate (fldv%emis      )
      deallocate (fldv%z0ma      )
      deallocate (fldv%trad      )
      deallocate (fldv%ustar     )
      deallocate (fldv%tstar     )
      deallocate (fldv%qstar     )
      deallocate (fldv%zol       )
      deallocate (fldv%rib       )
      deallocate (fldv%fm        )
      deallocate (fldv%fh        )
      deallocate (fldv%fq        )
      deallocate (fldv%tref      )
      deallocate (fldv%qref      )
      deallocate (fldv%u10m      )
      deallocate (fldv%v10m      )
      deallocate (fldv%f10m      )
      deallocate (fldv%qsubl     )
#ifdef DUST
      deallocate (fldv%dustemis_bin_1     )
      deallocate (fldv%dustemis_bin_2     )
      deallocate (fldv%dustemis_bin_3     )
      deallocate (fldv%dustemis_bin_4     )
      deallocate (fldv%dustemis_total     )
#endif

! average from column level
      deallocate (fldv%xerr      )
      deallocate (fldv%rsur      )
      deallocate (fldv%rnof      )
      deallocate (fldv%tg        )
      deallocate (fldv%scv       )
      deallocate (fldv%snowdp    )
      deallocate (fldv%fsno      )
      deallocate (fldv%us        )
      deallocate (fldv%vs        )
      deallocate (fldv%tm        )
      deallocate (fldv%qm        )
      deallocate (fldv%prc       )
      deallocate (fldv%prl       )
      deallocate (fldv%pbot      )
      deallocate (fldv%frl       )
      deallocate (fldv%solar     )
      deallocate (fldv%mrsos     )
      deallocate (fldv%mrso      )
      deallocate (fldv%mrfso     )
      deallocate (fldv%lwsnl     )
      deallocate (fldv%snm       )
      deallocate (fldv%tsn       )
      deallocate (fldv%nsnow     )

      deallocate (fldv%treeFrac     )
      deallocate (fldv%shrubFrac    )
      deallocate (fldv%grassFrac    )
      deallocate (fldv%baresoilFrac )
      deallocate (fldv%residualFrac )
      deallocate (fldv%soilFrac     )
      deallocate (fldv%urbanFrac    )
      deallocate (fldv%wetlandFrac  )
      deallocate (fldv%iceFrac      )
      deallocate (fldv%lakeFrac     )

      deallocate (fldv%tss  )
      deallocate (fldv%wliq )
      deallocate (fldv%wice )
      deallocate (fldv%mrlsl)
      deallocate (fldv%t_lake      )
      deallocate (fldv%lake_icefrac)
      deallocate (fldv%savedtke1   )
#if (defined FHNP) && (defined GWUSE)
      !Note that GWUSE is still in the testing phase, don't define it.  
      deallocate (fldv_col%qcharge_ori)  !wanglh 2018/11
      deallocate (fldv_col%cgwintake)    !wanglh 2018/11
      deallocate (fldv%gw_uptake)        !wanglh 2018/11
      deallocate (fldv%noirr_frac)       !wanglh 2018/11
      deallocate (fldv%ggw)        !wanglh 2018/11
#endif

#if (defined FHNP) && (defined FTF)
      deallocate(fldv_col%frostdp0     )
      deallocate(fldv_col%frostdp      )
      deallocate(fldv_col%thawdp       )
      deallocate(fldv_col%D_temperature)
      deallocate(fldv_col%N_time       )
      deallocate(fldv_col%frost_day    )
      deallocate(fldv_col%thaw_day     )
      deallocate(fldv%frostdp0     )
      deallocate(fldv%frostdp      )
      deallocate(fldv%thawdp       )
      deallocate(fldv%D_temperature)
      deallocate(fldv%N_time       )
      deallocate(fldv%frost_day    )
      deallocate(fldv%thaw_day     )
#endif
   end subroutine fldv_free

end module colm_varMod
