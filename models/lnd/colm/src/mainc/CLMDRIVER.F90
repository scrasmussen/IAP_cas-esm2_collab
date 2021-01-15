#include <define.h>

 subroutine CLMDRIVER (dolai,doalb,dosst) 

!=======================================================================
!
! CLM MODEL DRIVER
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002
!
!=======================================================================

 use precision
 use phycon_module, only : tfrz, rgas, vonkar
 use paramodel
 use spmd
 use spmd_decomp, only : ppmap, cgmap, gxmap, gymap, gxmap_glob, gymap_glob,ggmap !wanglh add the last 3
 use colm_varMod, only : numpatch, numcolumn, numgrid, ftune, forc, &
                         cvar, pvar, fldv_col, fldv_pft, oro, &
                         area,&  !for IAPDGVM
                         fLitterSoil, fLitterAtmos, &
                         Isf_pft, Iss_pft, Ksf_pft, Kss_pft, latixy, longxy,fldv!wanglh add the last 3
 use timemgr, only : idate, idate_p, dtime, get_nstep, get_curr_date  ! IAPDGVM
#ifdef VEGDATA
 use vegdata, only : interp_vegdata
#endif
#ifdef IAPDGVM
 use colm_ioMod,  only : getig, getrhm
#endif
#ifdef AUTOMASK
 use metdata, only : forcmask
#endif
 use spmd_decomp, only : pgmap, cgmap, gxmap, gymap, gxmap_glob, gymap_glob, ggmap !wanglh add the last 3
 use nanMod, only : nan
 use debug
 implicit none

! ------------------- arguments  ----------------------------------

  logical, INTENT(in) :: dolai    ! true if time for time-varying vegetation paramter
  logical, INTENT(in) :: doalb    ! true if time for surface albedo calculation
  logical, INTENT(in) :: dosst    ! true if time for update sst/ice/snow

! ----------------------------------------------------------------
! I. Time invariant model variables
! ----------------------------------------------------------------

  integer , pointer :: itypwat          ! land water type
  real(r8), pointer :: dlat             ! latitude in radians
  real(r8), pointer :: dlon             ! longitude in radians
  real(r8), pointer :: area_xy          ! IAPDGVM
  real(r8), pointer :: wxy_column
  real(r8), pointer :: wxy_patch(:)

                       ! Soil physical parameters
  real(r8), pointer :: albsat(:)        ! wet soil albedo for different coloured soils [-]
  real(r8), pointer :: albdry(:)        ! dry soil albedo for different coloured soils [-]
  real(r8), pointer :: csol  (:)        ! heat capacity of soil solids [J/(m3 K)]
  real(r8), pointer :: porsl (:)        ! fraction of soil that is voids [-]
  real(r8), pointer :: phi0  (:)        ! minimum soil suction [mm]
  real(r8), pointer :: bsw   (:)        ! clapp and hornbereger "b" parameter [-]
  real(r8), pointer :: dkmg  (:)        ! thermal conductivity of soil minerals [W/m-K]
  real(r8), pointer :: dksatu(:)        ! thermal conductivity of saturated soil [W/m-K]
  real(r8), pointer :: dkdry (:)        ! thermal conductivity for dry soil  [W/(m-K)]
  real(r8), pointer :: hksati(:)        ! hydraulic conductivity at saturation [mm h2o/s]

#ifdef DUST
                       ! Column level variables
  integer,  pointer :: dustsource          ! index for potential dust source
  integer,  pointer :: isltyp              ! dominant soil type 
  real (r8),pointer :: soil_top_cat(:)     ! fraction for 12-categories soil
  real (r8),pointer :: mvegcov(:)          ! read_in monthly vegetation cover [unit: %, i.e. 0-100]            
  real (r8)         :: vegcov_r            ! vegetation cover for current timestep      
#endif  
                       ! Vegetation static parameters
  integer , pointer :: ivt(:)              ! land cover type number  add by zhq 06/02/2009
  real(r8), pointer :: z0m(:)              ! aerodynamic roughness length [m]
  real(r8), pointer :: displa(:)           ! displacement height [m]
  real(r8), pointer :: sqrtdi(:)           ! inverse sqrt of leaf dimension [m**-0.5]
  real(r8), pointer :: effcon(:)           ! quantum efficiency of RuBP regeneration (molCO2/molquanta)
  real(r8), pointer :: vmax25(:)           ! maximum carboxylation rate at 25 C at canopy top
  real(r8), pointer :: slti(:)             ! s3: slope of low temperature inhibition function     
  real(r8), pointer :: hlti(:)             ! s4: 1/2 point of low temperature inhibition function
  real(r8), pointer :: shti(:)             ! s1: slope of high temperature inhibition function  
  real(r8), pointer :: hhti(:)             ! s2: 1/2 point of high temperature inhibition function 
  real(r8), pointer :: trda(:)             ! s5: temperature coefficient in gs-a model            
  real(r8), pointer :: trdm(:)             ! s6: temperature coefficient in gs-a model           
  real(r8), pointer :: trop(:)             ! temperature coefficient in gs-a model          
  real(r8), pointer :: gradm(:)            ! conductance-photosynthesis slope parameter
  real(r8), pointer :: binter(:)           ! conductance-photosynthesis intercep
  real(r8), pointer :: extkn(:)            ! coefficient of leaf nitrogen allocation
  real(r8), pointer :: chil(:)             ! leaf angle distribution factor
  real(r8), pointer :: refl  (:,:,:)       ! leaf reflectance (iw=iband, il=life and dead)
  real(r8), pointer :: refs  (:,:,:)       ! stem reflectance (iw=iband, il=life and dead)
  real(r8), pointer :: tranl (:,:,:)       ! leaf transmittance (iw=iband, il=life and dead)
  real(r8), pointer :: trans (:,:,:)       ! stem transmittance (iw=iband, il=life and dead)
  real(r8), pointer :: rootfr(:,:)         ! fraction of roots in each soil layer
#if(defined DGVM)
  real(r8), pointer :: pftpara(:,:)        ! 32 parameters of PFTs
  real(r8), pointer :: vegclass(:)         ! 1.tree 2.shrub 3.grass 4.crop -1.others
  real(r8), pointer :: summergreen(:)      ! 1. for summergreen; otherwise -1.
  real(r8), pointer :: raingreen(:)        ! 1. for raingreen; otherwise -1.
  real(r8), pointer :: sla(:)              ! sla
  real(r8), pointer :: lm_sapl(:)          ! leafmass
  real(r8), pointer :: sm_sapl(:)          ! sapwood mass
  real(r8), pointer :: hm_sapl(:)          ! heartwood mass
  real(r8), pointer :: rm_sapl(:)          ! rootmass
#endif
#if(defined DyN)
  real(r8), pointer :: cton_soil(:)        ! soil C:N mass ratio
  real(r8), pointer :: cton_pro (:)        ! C:N mass ratio in production
#endif

              ! CLM time step and TUNABLE constants
  real(r8) :: zlnd             ! roughness length for soil [m]
  real(r8) :: zsno             ! roughness length for snow [m]
  real(r8) :: csoilc           ! drag coefficient for soil under canopy [-]
  real(r8) :: dewmx            ! maximum dew
  real(r8) :: wtfact           ! fraction of model area with high water table
  real(r8) :: capr             ! tuning factor to turn first layer T into surface T
  real(r8) :: cnfac            ! Crank Nicholson factor between 0 and 1 
  real(r8) :: ssi              ! irreducible water saturation of snow
  real(r8) :: wimp             ! water impremeable if porosity less than wimp
  real(r8) :: pondmx           ! ponding depth (mm)
  real(r8) :: smpmax           ! wilting point potential in mm
  real(r8) :: smpmin           ! restriction for min of soil poten. (mm)
  real(r8) :: trsmx0           ! max transpiration for moist soil+100% veg.  [mm/s]
  real(r8) :: tcrit            ! critical temp. to determine rain or snow

! -----------------------------------------------------------------
! II. Time-varying state variables which reaquired by restart run
! -----------------------------------------------------------------
                       ! Currrnt cal2ar 
  integer :: year             ! current year of model run
  integer :: jday             ! current julian day of model run 
  integer :: msec             ! current seconds of model run (0 - 86400)

                       ! Main land surface variables 
  real(r8), pointer :: z   (:)          ! node depth [m]
  real(r8), pointer :: dz  (:)          ! interface depth [m]
  real(r8), pointer :: tss (:)          ! soil temperature [K]
  real(r8), pointer :: wliq(:)          ! liquid water in layers [kg/m2]
  real(r8), pointer :: wice(:)          ! ice lens in layers [kg/m2]
  real(r8), pointer :: tg               ! ground surface temperature [K]
  real(r8), pointer :: tlsun(:)         ! sunlit leaf temperature [K]
  real(r8), pointer :: tlsha(:)         ! shaded leaf temperature [K]
  real(r8), pointer :: ldew(:)          ! depth of water on foliage [mm]
  real(r8), pointer :: sag              ! non dimensional snow age [-]
  real(r8), pointer :: scv              ! snow cover, water equivalent [mm]
  real(r8), pointer :: snowdp           ! snow depth [meter]
  real(r8), pointer :: lakedepth        ! lake depth [m] added by Nan Wei
  real(r8), pointer :: dz_lake(:)       ! lake thickness [m] added by Nan Wei
  real(r8), pointer :: t_lake(:)        ! lake layer temperature [K] added by Nan Wei
  real(r8), pointer :: lake_icefrac(:)  ! lake mass fraction of lake layer that is frozen added by Nan Wei
  real(r8), pointer :: savedtke1        ! added by Nan Wei

                       ! Vegetation dynamic parameters 
  real(r8), pointer :: fveg(:)             ! fraction of vegetation cover
  real(r8), pointer :: fsno                ! fraction of snow cover on ground
  real(r8), pointer :: sigf(:)             ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), pointer :: green(:)            ! leaf greenness
  real(r8), pointer :: lai(:)              ! leaf area index
  real(r8), pointer :: sai(:)              ! stem area index

#if (defined DGVM)
  real(r8), pointer :: t10min(:)              !annual minimum of 10-day running mean (K)
  real(r8), pointer :: lai_ind(:)             !LAI per individual
  real(r8), pointer :: dphen(:)               !phenology [0 to 1]
  real(r8), pointer :: leafon(:)              !leafon days
  real(r8), pointer :: leafof(:)              !leafoff days
  real(r8), pointer :: firelength(:)          !fire season in days
#if (defined IAPDGVM)  
  real(r8), pointer :: afirefrac1(:)          !IAPDGVM zjw
  real(r8), pointer :: nfireg1(:)             !IAPDGVM zjw
#endif
  real(r8), pointer :: litterag(:)            !above ground litter
  real(r8), pointer :: litterbg(:)            !below ground litter
  real(r8), pointer :: cpool_fast(:,:)        !fast carbon pool
  real(r8), pointer :: cpool_slow(:,:)        !slow carbon pool
  real(r8), pointer :: k_fast_ave(:)          !decomposition rate
  real(r8), pointer :: k_slow_ave(:)          !decomposition rate
  real(r8), pointer :: litter_decom_ave(:)    !decomposition rate
  real(r8), pointer :: fmicr(:)               !microbial respiration (mol CO2 /m**2 /s)
  real(r8), pointer :: nind(:)                !number of individuals (#/m**2)
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
  real(r8), pointer :: crownarea(:)           !area that each individual tree takes up (m^2)
  real(r8), pointer :: htop(:)                !canopy top
  real(r8), pointer :: tsai(:)                !one-sided stem area index, no burying by snow
  real(r8), pointer :: fpcgrid(:)             !foliar projective cover on gridcell (fraction)
  real(r8), pointer :: bm_inc(:)              !biomass increment
  real(r8), pointer :: afmicr(:)              !annual microbial respiration
  real(r8), pointer :: annpsn(:)              !annual photosynthesis (umol CO2 /m**2)
  real(r8), pointer :: annpsnpot(:)           !annual potential photosynthesis (same units)
  real(r8), pointer :: tref10(:)              !10-day averaged temperature at 2m
  real(r8), pointer :: tref_sum(:)            !sum of temperature in current day
  real(r8), pointer :: t10(:,:)               !array to record the 10 day temperature
  real(r8), pointer :: assimn10(:)            !10-day averaged assimilation rate
  real(r8), pointer :: assimn_sum(:)          !sum of assimn of current day
  real(r8), pointer :: an10(:,:)              !array to record 10 day assimn
  real(r8), pointer :: anngpp(:)              !annual gpp
  real(r8), pointer :: annfrmf(:)             !annual frmf
  real(r8), pointer :: annfrms(:)             !annual frms
  real(r8), pointer :: annfrmr(:)             !annual frmr
  real(r8), pointer :: annfrg(:)              !annual frg
  real(r8), pointer :: cflux_litter_soil(:)   !litter->soil
  real(r8), pointer :: cflux_litter_atmos(:)  !litter->atmos
  real(r8), pointer :: nday                   !counting the model days
  real(r8), pointer :: nyr                    !counting the model years
  real(r8), pointer :: turnover_ind(:)        !individual turnover biomass
  real(r8), pointer :: fpc_inc(:)             !fpc increase
  real(r8), pointer :: prec365                !yearly running mean of precipitation(mm/s)
#if(defined IAPDGVM)
  real(r8), pointer :: wliq6mon               ! IAPDGVM liquid water 6 mons for the first 3 layers
#endif

  real(r8), pointer :: ifpre(:)               !-1=no PFT present;1=PFT present in this grid
#if(defined DyN)
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
#endif

#if(defined VEGDATA)
  real(r8) :: lai_r  (numpft+1)           ! read-in leaf area index
  real(r8) :: sai_r  (numpft+1)           ! read-in stem+dead leaf area index
  real(r8) :: green_r(numpft+1)           ! read-in leaf greenness
  real(r8) :: fveg_r (numpft+1)           ! read-in vegetation fraction cover
  real(r8) :: htop_r (numpft+1)           ! read-in vegetation height
#endif

                       ! Radiation  related (albedoes)
  real(r8), pointer :: coszen             ! cosine of solar zenith angle
  real(r8), pointer :: albg (:,:,:)       ! albedo, ground [-]
  real(r8), pointer :: albv (:,:,:)       ! albedo, vegetation [-]
  real(r8), pointer :: alb  (:,:,:)       ! averaged albedo [-]
  real(r8), pointer :: ssun (:,:,:)       ! sunlit canopy absorption for solar radiation (0-1)
  real(r8), pointer :: ssha (:,:,:)       ! shaded canopy absorption for solar radiation (0-1)
  real(r8), pointer :: thermk(:)          ! canopy gap fraction for tir radiation
  real(r8), pointer :: extkb(:)           ! (k, g(mu)/mu) direct solar extinction coefficient
  real(r8), pointer :: extkd(:)           ! diffuse and scattered diffuse PAR extinction coefficient

                       ! Additional variables required by reginal model (WRF & RSM) 
  real(r8), pointer :: trad(:)            ! radiative temperature of surface [K]
  real(r8), pointer :: tref(:)            ! 2 m height air temperature [kelvin]
  real(r8), pointer :: qref(:)            ! 2 m height air specific humidity
  real(r8), pointer :: rst(:)             ! canopy stomatal resistance (s/m)

  real(r8), pointer :: emis(:)            ! averaged bulk surface emissivity
  real(r8), pointer :: z0ma(:)            ! effective roughness [m]
  real(r8), pointer :: zol(:)             ! dimensionless height (z/L) used in Monin-Obukhov theory
  real(r8), pointer :: rib(:)             ! bulk Richardson number in surface layer
  real(r8), pointer :: ustar(:)           ! u* in similarity theory [m/s]
  real(r8), pointer :: qstar(:)           ! q* in similarity theory [kg/kg]
  real(r8), pointer :: tstar(:)           ! t* in similarity theory [K]
  real(r8), pointer :: fm(:)              ! integral of profile function for momentum
  real(r8), pointer :: fh(:)              ! integral of profile function for heat
  real(r8), pointer :: fq(:)              ! integral of profile function for moisture

! -----------------------------------------------------------------
! III. Forcing 
! -----------------------------------------------------------------

  real(r8), pointer :: pco2m            ! CO2 concentration in atmos. (35 pa)
  real(r8), pointer :: po2m             ! O2 concentration in atmos. (20900 pa)
  real(r8), pointer :: us               ! wind in eastward direction [m/s]
  real(r8), pointer :: vs               ! wind in northward direction [m/s]
  real(r8), pointer :: tm               ! temperature at reference height [kelvin]
  real(r8), pointer :: qm               ! specific humidity at reference height [kg/kg]
  real(r8), pointer :: prc              ! convective precipitation [mm/s]
  real(r8), pointer :: prl              ! large scale precipitation [mm/s]
  real(r8), pointer :: rain             ! rainfall [mm/s]
  real(r8), pointer :: snow             ! snowfall [mm/s]
  real(r8), pointer :: pbot             ! atm bottom level pressure (or reference height) (pa)
  real(r8), pointer :: psrf             ! atmospheric pressure at the surface [pa]
  real(r8), pointer :: sols             ! atm vis direct beam solar rad onto srf [W/m2]
  real(r8), pointer :: soll             ! atm nir direct beam solar rad onto srf [W/m2]
  real(r8), pointer :: solsd            ! atm vis diffuse solar rad onto srf [W/m2]
  real(r8), pointer :: solld            ! atm nir diffuse solar rad onto srf [W/m2]
  real(r8), pointer :: frl              ! atmospheric infrared (longwave) radiation [W/m2]
  real(r8), pointer :: hu               ! observational height of wind [m]
  real(r8), pointer :: ht               ! observational height of temperature [m]
  real(r8), pointer :: hq               ! observational height of humidity [m]

  real(r8) :: rhoair                    ! air density [kg/m3]

! -----------------------------------------------------------------
! IV. Fluxes
! -----------------------------------------------------------------

  real(r8), pointer :: taux(:)             ! wind stress: E-W [kg/m/s2]
  real(r8), pointer :: tauy(:)             ! wind stress: N-S [kg/m/s2]
  real(r8), pointer :: fsena(:)            ! sensible heat from canopy height to atmosphere [W/m2]
  real(r8), pointer :: lfevpa(:)           ! latent heat flux from canopy height to atmosphere [W/m2]
  real(r8), pointer :: fevpa(:)            ! evapotranspiration from canopy height to atmosphere [mm/s]
  real(r8), pointer :: fsenl(:)            ! sensible heat from leaves [W/m2]
  real(r8), pointer :: fevpl(:)            ! evaporation+transpiration from leaves [mm/s]
  real(r8), pointer :: etr(:)              ! transpiration rate [mm/s]
  real(r8), pointer :: fseng(:)            ! sensible heat flux from ground [W/m2]
  real(r8), pointer :: fevpg(:)            ! evaporation heat flux from ground [mm/s]
  real(r8), pointer :: fgrnd(:)            ! ground heat flux [W/m2]
  real(r8), pointer :: sabvsun(:)          ! solar absorbed by sunlit vegetation [W/m2]
  real(r8), pointer :: sabvsha(:)          ! solar absorbed by shaded vegetation [W/m2]
  real(r8), pointer :: sabg(:)             ! solar absorbed by ground  [W/m2]
  real(r8), pointer :: olrg(:)             ! outgoing long-wave radiation from ground+canopy [W/m2]
! real(r8), pointer :: rnet(:)             ! net radiation by surface [W/m2]
  real(r8), pointer :: xerr                ! the error of water banace [mm/s]
  real(r8), pointer :: zerr(:)             ! the error of energy balance [W/m2]

  real(r8), pointer :: rsur                ! surface runoff (mm h2o/s)
  real(r8), pointer :: rnof                ! total runoff (mm h2o/s)
  real(r8), pointer :: assim(:)            ! canopy assimilation rate (mol m-2 s-1)
  real(r8), pointer :: respc(:)            ! canopy respiration (mol m-2 s-1)

  real(r8), pointer :: u10m(:)             ! 10m u-velocity 
  real(r8), pointer :: v10m(:)             ! 10m v-velocity 
  real(r8), pointer :: f10m(:)             ! integral of profile function for momentum at 10m

#ifdef DUST
                       ! PFT (Patch) level variables
  real(r8), pointer :: dustemis_bin_1(:)           ! dust emission flux for four size bins [kg/m2/s]
  real(r8), pointer :: dustemis_bin_2(:)           
  real(r8), pointer :: dustemis_bin_3(:)           
  real(r8), pointer :: dustemis_bin_4(:)            
  real(r8), pointer :: dustemis_total(:)         ! total dust emission flux [kg/m2/s]
#endif

! -----------------------------------------------------------------
! V. Local declaration
! -----------------------------------------------------------------

  real(r8) :: parsun(numpft+1)             ! PAR by sunlit leaves [W/m2] 
  real(r8) :: parsha(numpft+1)             ! PAR by shaded leaves [W/m2]
  real(r8) :: sabvg(numpft+1)              ! solar absorbed by ground + vegetation [W/m2]

  real(r8), pointer :: Isf(:,:)
  real(r8), pointer :: Iss(:,:)
  real(r8), pointer :: Ksf(:,:)
  real(r8), pointer :: Kss(:,:)

  real(r8) :: snm                          ! rate of snowmelt [kg/(m2 s)]
  real(r8) :: qsubl(numpft+1)              ! sublimation rate from snow pack (mm h2o /s) [+]

  integer  :: snl                          ! number of snow layers
  real(r8) :: lwsnl                        ! mass of liquid water of snow layers [kg/m2]
  real(r8) :: nsnow                        ! number of snow events [-]
  real(r8) :: tsn                          ! snow internal temperature [K]
  real(r8) :: mrsos                        ! mass of water of all phases in the upper 0.1 meters of soil [kg/m2]
  real(r8) :: mrso                         ! mass of water of all phases over all soil layers [kg/m2]
  real(r8) :: mrfso                        ! mass of frozen water over all soil layers [kg/m2]
  real(r8) :: mrlsl(1:nl_soil)             ! mass of water of all phases in each soil layer [kg/m2]
#if(defined IAPDGVM)
  real(r8) :: getig_in                                       
  real(r8) :: getrhm_in                    

  integer :: yr            ! current year
  integer :: mon           ! current month
  integer :: day           ! current day (0, 1, ...)
  integer :: sec         ! current seconds of current date (0, ..., 86400)
#endif
  integer i,j,lc,uc,lb,ub,jm               ! loop/array indices
  integer c,p,p1,p2,npatch
#if (defined FHNP) && (defined FTF)
!liruichao add
  real(r8), pointer :: frostdp          !frost depth
  real(r8), pointer :: thawdp           !thaw deppth
  real(r8), pointer :: frostdp0         !initial frost depth
  real(r8), pointer :: D_temperature    !frost or thaw index
  real(r8), pointer :: N_time           !step counter
  real(r8), pointer :: frost_day        !frost days
  real(r8), pointer :: thaw_day         !thaw days
!end
#endif

#if (defined FHNP) && (defined GWUSE)
  !Note that GWUSE is still in the testing phase, don't define it.  
  integer g_index,c_index,myrank,g1,i_index,j_index !wanglh
  real(r8), pointer :: llatixy(:,:)
  real(r8), pointer :: llongxy(:,:)
#endif

! -----------------------------------------------------------------
      parsun(:) = 0.
      parsha(:) = 0.
      sabvg(:)  = 0.
      qsubl(:)  = 0.

      fLitterSoil(:)   = 0.
      fLitterAtmos(:)  = 0.

      jm = nl_soil+abs(maxsnl) 

      ! CLM time step and TUNABLE constants

      zlnd   = ftune(1)
      zsno   = ftune(2)
      csoilc = ftune(3)
      dewmx  = ftune(4)
      wtfact = ftune(5)
      capr   = ftune(6)
      cnfac  = ftune(7)
      ssi    = ftune(8)
      wimp   = ftune(9)
      pondmx = ftune(10) 
      smpmax = ftune(11) 
      smpmin = ftune(12) 
      trsmx0 = ftune(13) 
      tcrit  = ftune(14)  

      year   = idate(1)      
      jday   = idate(2)      
      msec   = idate(3)      

#ifdef MYBUG
      xstep = xstep + 1
      if(p_master) write(6,*) 'xstep', xstep
#endif

      p1 = 1
      p2 = 0

      DO c = 1, numcolumn
#if (defined FHNP) && (defined GWUSE)
#ifdef SPMD
  call mpi_comm_rank(p_comm,myrank,p_err)
#endif
!Note that GWUSE is still in the testing phase, don't define it. 
         g_index=cgmap(c) !wanglh 2018/11
         c_index=c

         llatixy=>latixy
         llongxy=>longxy
         g1=ggmap(g_index)

         i_index = gxmap_glob(g1)
         j_index = gymap_glob(g1)
#endif
! ======================================================================
!  [1] Transfer the time invariant and time-varying variables
! ======================================================================
                            ! Time invariant model variables

         dlat               => cvar%dlat(c)
         dlon               => cvar%dlon(c)
         itypwat            => cvar%itypwat(c)
         wxy_column         => cvar%wxy_column(c)
         albsat             => cvar%albsat(:,c)
         albdry             => cvar%albdry(:,c)
         csol               => cvar%csol(1:nl_soil,c)
         porsl              => cvar%porsl(1:nl_soil,c)
         phi0               => cvar%phi0(1:nl_soil,c)
         bsw                => cvar%bsw(1:nl_soil,c)
         dkmg               => cvar%dkmg(1:nl_soil,c)
         dksatu             => cvar%dksatu(1:nl_soil,c)
         dkdry              => cvar%dkdry(1:nl_soil,c)
         hksati             => cvar%hksati(1:nl_soil,c)

#ifdef DUST
         dustsource         => cvar%dustsource(c)
         isltyp             => cvar%isltyp(c)
         soil_top_cat       => cvar%soil_top_cat(:,c)
         mvegcov            => cvar%mvegcov(:,c)
#endif

#if(defined SINGLE)
         if(itypwat==0)      npatch = 1       ! natural vegetation + bare soil
#else
         if(itypwat==0)      npatch = numpft+1! natural vegetation + bare soil
#endif
         if(itypwat==1)      npatch = 1       ! urban and built-up
         if(itypwat==2)      npatch = 1       ! wetland
         if(itypwat==3)      npatch = 1       ! land ice
         if(itypwat==4)      npatch = 1       ! river or deep lake
         if(itypwat==99)then                  ! ocean
            write(6,*) 'ocean column',c
            stop
         endif

         p2 = p2 + npatch

#ifdef AUTOMASK
         i = gxmap(cgmap(c))
         j = gymap(cgmap(c))

         if(forcmask(i,j).eq.1) then
            p1 = p2 + 1
            cycle
         endif
#endif

         ivt                => pvar%ivt(p1:p2)
         wxy_patch          => pvar%wxy_patch(p1:p2)
         z0m                => pvar%z0m(p1:p2)
         displa             => pvar%displa(p1:p2)
         sqrtdi             => pvar%sqrtdi(p1:p2)
         effcon             => pvar%effcon(p1:p2)
         vmax25             => pvar%vmax25(p1:p2)
         slti               => pvar%slti(p1:p2)
         hlti               => pvar%hlti(p1:p2)
         shti               => pvar%shti(p1:p2)
         hhti               => pvar%hhti(p1:p2)
         trda               => pvar%trda(p1:p2)
         trdm               => pvar%trdm(p1:p2)
         trop               => pvar%trop(p1:p2)
         gradm              => pvar%gradm(p1:p2)
         binter             => pvar%binter(p1:p2)
         extkn              => pvar%extkn(p1:p2)
         chil               => pvar%chil(p1:p2)
         refl               => pvar%refl(1:2,1:2,p1:p2)
         refs               => pvar%refs(1:2,1:2,p1:p2)
         tranl              => pvar%tranl(1:2,1:2,p1:p2)
         trans              => pvar%trans(1:2,1:2,p1:p2)
         rootfr             => pvar%rootfr(1:nl_soil,p1:p2)
#if(defined DGVM)
         pftpara            => pvar%pftpara(1:npftpara,p1:p2)
         vegclass           => pvar%vegclass(p1:p2)
         summergreen        => pvar%summergreen(p1:p2)
         raingreen          => pvar%raingreen(p1:p2)
         sla                => pvar%sla(p1:p2)
         lm_sapl            => pvar%lm_sapl(p1:p2)
         sm_sapl            => pvar%sm_sapl(p1:p2)
         hm_sapl            => pvar%hm_sapl(p1:p2)
         rm_sapl            => pvar%rm_sapl(p1:p2)
#endif
#if(defined DyN)
         cton_soil          => pvar%cton_soil(p1:p2)
         cton_pro           => pvar%cton_pro(p1:p2)
#endif
                            ! Time-varying variables 
         z                  => cvar%z(maxsnl+1:nl_soil,c)
         dz                 => cvar%dz(maxsnl+1:nl_soil,c)
         tss                => cvar%tss(maxsnl+1:nl_soil,c)
         wliq               => cvar%wliq(maxsnl+1:nl_soil,c)
         wice               => cvar%wice(maxsnl+1:nl_soil,c)
         tg                 => cvar%tg(c)
         sag                => cvar%sag(c)
         scv                => cvar%scv(c)
         snowdp             => cvar%snowdp(c)
         fsno               => cvar%fsno(c)
         coszen             => cvar%coszen(c)
         lakedepth          => cvar%lakedepth(c)
         dz_lake            => cvar%dz_lake     (1:nl_lake,c)
         t_lake             => cvar%t_lake      (1:nl_lake,c)
         lake_icefrac       => cvar%lake_icefrac(1:nl_lake,c)
         savedtke1          => cvar%savedtke1(c)

         if(scv<0.)then
            write(6,*),'driversnow',itypwat,tg,scv,snowdp
         endif

#if(defined DGVM)
         nday               => cvar%nday(c)
         nyr                => cvar%nyr(c)
         prec365            => cvar%prec365(c)
#if(defined IAPDGVM)
         wliq6mon           => cvar%wliq6mon(c)
#endif
#endif
#if(defined DyN)
         soil_no3           => cvar%soil_no3(c)
         soil_no2           => cvar%soil_no2(c)
         soil_no            => cvar%soil_no(c)
         soil_n2o           => cvar%soil_n2o(c)
         soil_n2            => cvar%soil_n2(c)
         soil_nh4           => cvar%soil_nh4(c)
#endif

#if (defined FHNP) && (defined FTF)
!liruichao add          
         frostdp            => cvar%frostdp(c)
         thawdp             => cvar%thawdp(c)
         frostdp0           => cvar%frostdp0(c)
         D_temperature      => cvar%D_temperature(c)
         N_time             => cvar%N_time(c)
         frost_day          => cvar%frost_day(c)
         thaw_day           => cvar%thaw_day(c)
!end
#endif

         tlsun              => pvar%tlsun(p1:p2)
         tlsha              => pvar%tlsha(p1:p2)
         ldew               => pvar%ldew(p1:p2)
         fveg               => pvar%fveg(p1:p2)
         sigf               => pvar%sigf(p1:p2)
         green              => pvar%green(p1:p2)
         lai                => pvar%lai(p1:p2)
         sai                => pvar%sai(p1:p2)
         albg               => pvar%albg(1:2,1:2,p1:p2)
         albv               => pvar%albv(1:2,1:2,p1:p2)
         alb                => pvar%alb (1:2,1:2,p1:p2)
         ssun               => pvar%ssun(1:2,1:2,p1:p2)
         ssha               => pvar%ssha(1:2,1:2,p1:p2)
         thermk             => pvar%thermk(p1:p2)
         extkb              => pvar%extkb(p1:p2)
         extkd              => pvar%extkd(p1:p2)

                            ! Additional variables required by reginal model (WRF & RSM) 
         trad               => pvar%trad(p1:p2)
         tref               => pvar%tref(p1:p2)
         qref               => pvar%qref(p1:p2)
         rst                => pvar%rst(p1:p2)
         emis               => pvar%emis(p1:p2)
         z0ma               => pvar%z0ma(p1:p2)
         zol                => pvar%zol(p1:p2)
         rib                => pvar%rib(p1:p2)
         ustar              => pvar%ustar(p1:p2)
         qstar              => pvar%qstar(p1:p2)
         tstar              => pvar%tstar(p1:p2)
         fm                 => pvar%fm(p1:p2)
         fh                 => pvar%fh(p1:p2)
         fq                 => pvar%fq(p1:p2)
#if(defined DGVM)
         t10min             => pvar%t10min(p1:p2)
         lai_ind            => pvar%lai_ind(p1:p2)
         dphen              => pvar%dphen(p1:p2)
         leafon             => pvar%leafon(p1:p2)
         leafof             => pvar%leafof(p1:p2)
         firelength         => pvar%firelength(p1:p2)
         litterag           => pvar%litterag(p1:p2)
         litterbg           => pvar%litterbg(p1:p2)
         cpool_fast         => pvar%cpool_fast(1:nl_soil,p1:p2)
         cpool_slow         => pvar%cpool_slow(1:nl_soil,p1:p2)
         k_fast_ave         => pvar%k_fast_ave(p1:p2)
         k_slow_ave         => pvar%k_slow_ave(p1:p2)
         litter_decom_ave   => pvar%litter_decom_ave(p1:p2)
         fmicr              => pvar%fmicr(p1:p2)
         nind               => pvar%nind(p1:p2)
         lm_ind             => pvar%lm_ind(p1:p2)
         sm_ind             => pvar%sm_ind(p1:p2)
         hm_ind             => pvar%hm_ind(p1:p2)
         rm_ind             => pvar%rm_ind(p1:p2)
         tmomin20           => pvar%tmomin20(p1:p2)
         agdd0              => pvar%agdd0(p1:p2)
         agdd               => pvar%agdd(p1:p2)
         agdd20             => pvar%agdd20(p1:p2)
         t_mo_min           => pvar%t_mo_min(p1:p2)
         crownarea          => pvar%crownarea(p1:p2)
         htop               => pvar%htop(p1:p2)
         tsai               => pvar%tsai(p1:p2)
         fpcgrid            => pvar%fpcgrid(p1:p2)
         bm_inc             => pvar%bm_inc(p1:p2)
         afmicr             => pvar%afmicr(p1:p2)
         annpsn             => pvar%annpsn(p1:p2)
         annpsnpot          => pvar%annpsnpot(p1:p2)
         tref10             => pvar%tref10(p1:p2)
         tref_sum           => pvar%tref_sum(p1:p2)
         t10                => pvar%t10(1:10,p1:p2)
         assimn10           => pvar%assimn10(p1:p2)
         assimn_sum         => pvar%assimn_sum(p1:p2)
         an10               => pvar%an10(1:10,p1:p2)
         turnover_ind       => pvar%turnover_ind(p1:p2)
         fpc_inc            => pvar%fpc_inc(p1:p2)
         agddtw             => pvar%agddtw(p1:p2)
         ifpre              => pvar%ifpre(p1:p2)
         t_mo               => pvar%t_mo(p1:p2)
         t_mo_sum           => pvar%t_mo_sum(p1:p2)
         anngpp             => pvar%anngpp(p1:p2)
         annfrmf            => pvar%annfrmf(p1:p2)
         annfrms            => pvar%annfrms(p1:p2)
         annfrmr            => pvar%annfrmr(p1:p2)
         annfrg             => pvar%annfrg(p1:p2)
#if(defined IAPDGVM)
         afirefrac1         => pvar%afirefrac1(p1:p2)    
         nfireg1            => pvar%nfireg1(p1:p2)       
#endif

#if(defined DyN)
         litter_leaf        => pvar%litter_leaf(p1:p2)
         litter_wood        => pvar%litter_wood(p1:p2)
         litter_root        => pvar%litter_root(p1:p2)
         litter_repr        => pvar%litter_repr(p1:p2)
         litter_leaf_n      => pvar%litter_leaf_n(p1:p2)
         litter_wood_n      => pvar%litter_wood_n(p1:p2)
         litter_root_n      => pvar%litter_root_n(p1:p2)
         litter_repr_n      => pvar%litter_repr_n(p1:p2)
         afcton_leaf        => pvar%afcton_leaf(p1:p2)
         afcton_sap         => pvar%afcton_sap(p1:p2)
         afcton_root        => pvar%afcton_root(p1:p2)
         lm_ind_n           => pvar%lm_ind_n(p1:p2)
         sm_ind_n           => pvar%sm_ind_n(p1:p2)
         hm_ind_n           => pvar%hm_ind_n(p1:p2)
         rm_ind_n           => pvar%rm_ind_n(p1:p2)
         an_up              => pvar%an_up(p1:p2)
         an_stress          => pvar%an_stress(p1:p2)
#endif
         cflux_litter_soil  => fLitterSoil (p1:p2)
         cflux_litter_atmos => fLitterAtmos(p1:p2)
!
! fast spinup
!
         Isf                => Isf_pft(:,p1:p2)
         Iss                => Iss_pft(:,p1:p2)
         Ksf                => Ksf_pft(:,p1:p2)
         Kss                => Kss_pft(:,p1:p2)
#endif

#if(defined VEGDATA)

#ifdef DUST
         !Vegetation cover added by Chenglai Wu: 04/27/2014
         call interp_vegdata(itypwat,p1,p2,lai_r,sai_r,htop_r,mvegcov,vegcov_r)
#else
         call interp_vegdata(itypwat,p1,p2,lai_r,sai_r,htop_r)
#endif

#endif

! ======================================================================
!  [2] atmospheric fields to force clm
! ======================================================================

         pco2m    => forc%pco2m(c)   !CO2 concentration in atmos. (35 pa)
         po2m     => forc%po2m (c)   !O2 concentration in atmos. (20900 pa)
         us       => forc%us   (c)   !wind in eastward direction [m/s]
         vs       => forc%vs   (c)   !wind in northward direction [m/s]
         tm       => forc%tm   (c)   !temperature at reference height [kelvin]
         qm       => forc%qm   (c)   !specific humidity at reference height [kg/kg]
         prc      => forc%prc  (c)   !convective precipitation [mm/s]
         prl      => forc%prl  (c)   !large scale precipitation [mm/s]
         rain     => forc%rain (c)   !rainfall [mm/s]
         snow     => forc%snow (c)   !snowfall [mm/s]
         pbot     => forc%pbot (c)   !atm bottom level pressure (or reference height) (pa)
         psrf     => forc%psrf (c)   !atmospheric pressure at the surface [pa]
         sols     => forc%sols (c)   !atm vis direct beam solar rad onto srf [W/m2]
         soll     => forc%soll (c)   !atm nir direct beam solar rad onto srf [W/m2]
         solsd    => forc%solsd(c)   !atm vis diffuse solar rad onto srf [W/m2]
         solld    => forc%solld(c)   !atm nir diffuse solar rad onto srf [W/m2]
         frl      => forc%frl  (c)   !atmospheric infrared (longwave) radiation [W/m2]
         hu       => forc%hu   (c)   !observational height of wind [m]
         ht       => forc%ht   (c)   !observational height of temperature [m]
         hq       => forc%hq   (c)   !observational height of humidity [m]

         rhoair   =  (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm)

! ======================================================================
!  [3] surface fluxes    add by zhq 03/06/2009
! ======================================================================

         xerr     => fldv_col%xerr(c)                !the error of water banace [mm/s]
         rsur     => fldv_col%rsur(c)                !surface runoff [mm/s]
         rnof     => fldv_col%rnof(c)                !total runoff [mm/s]

         taux     => fldv_pft%taux   (p1:p2)         !wind stress: E-W [kg/m/s2]
         tauy     => fldv_pft%tauy   (p1:p2)         !wind stress: N-S [kg/m/s2]
         fsena    => fldv_pft%fsena  (p1:p2)         !sensible heat from canopy height to atmosphere [W/m2]
         lfevpa   => fldv_pft%lfevpa (p1:p2)         !latent heat flux from canopy height to atmosphere [W/m2]
         fevpa    => fldv_pft%fevpa  (p1:p2)         !evapotranspiration from canopy height to atmosphere [mm/s]
         fsenl    => fldv_pft%fsenl  (p1:p2)         !sensible heat from leaves [W/m2]
         fevpl    => fldv_pft%fevpl  (p1:p2)         !evaporation+transpiration from leaves [mm/s]
         etr      => fldv_pft%etr    (p1:p2)         !transpiration rate [mm/s]
         fseng    => fldv_pft%fseng  (p1:p2)         !sensible heat flux from ground [W/m2]
         fevpg    => fldv_pft%fevpg  (p1:p2)         !evaporation heat flux from ground [mm/s]
         fgrnd    => fldv_pft%fgrnd  (p1:p2)         !ground heat flux [W/m2]
         sabvsun  => fldv_pft%sabvsun(p1:p2)         !solar absorbed by sunlit canopy [W/m2]
         sabvsha  => fldv_pft%sabvsha(p1:p2)         !solar absorbed by shaded [W/m2]
         sabg     => fldv_pft%sabg   (p1:p2)         !solar absorbed by ground  [W/m2]
         olrg     => fldv_pft%olrg   (p1:p2)         !outgoing long-wave radiation from ground+canopy [W/m2]
         zerr     => fldv_pft%zerr   (p1:p2)         !the error of energy balance [W/m2]
         assim    => fldv_pft%assim  (p1:p2)         !canopy assimilation rate [mol m-2 s-1]
         respc    => fldv_pft%respc  (p1:p2)         !respiration (plant+soil) [mol m-2 s-1]
         u10m     => fldv_pft%u10m   (p1:p2)         !
         v10m     => fldv_pft%v10m   (p1:p2)         !
         f10m     => fldv_pft%f10m   (p1:p2)         !

#ifdef DUST
         dustemis_bin_1 => fldv_pft%dustemis_bin_1   (p1:p2)
         dustemis_bin_2 => fldv_pft%dustemis_bin_2   (p1:p2)
         dustemis_bin_3 => fldv_pft%dustemis_bin_3   (p1:p2)
         dustemis_bin_4 => fldv_pft%dustemis_bin_4   (p1:p2)
         dustemis_total => fldv_pft%dustemis_total   (p1:p2) 
#endif

         snm      = 0._r8
         qsubl(:) = 0._r8

#ifdef MYBUG
         i = gxmap(cgmap(c))
         j = gymap(cgmap(c))

         if(i.eq.ibug .and. j.eq.jbug) then
            mybug = .true.
         else
            mybug = .false.
         endif

         if(mybug) write(6,*), 'debug', p_iam, c, p1, p2

         if(mybug) then
            write(6,*) "col", itypwat, wxy_column
            write(6,*) "wxy_patch", wxy_patch(1:npatch)
            write(6,*) "lai", lai(1:npatch)
            write(6,*) "ivt", ivt(1:npatch)
            write(6,*) "fsena", fsena(1:npatch)
            write(6,*) "lfevpa", lfevpa(1:npatch)
            write(6,*) "olrg", olrg(1:npatch)
            write(6,*) "tg", tg
         endif
#endif

#if (defined IAPDGVM) 
   i=gxmap(cgmap(c))
   j=gymap(cgmap(c))
   area_xy => area(i,j)
   getig_in=getig(c,year,jday,msec)
   getrhm_in=getrhm(tm,pbot,qm)
#endif
! ======================================================================
!  [2] Main driver for CLM
! ======================================================================
         CALL CLMMAIN (dtime, doalb, dolai, dosst, &
         nl_soil, maxsnl, dlon, dlat, itypwat, npatch, wxy_column, wxy_patch, ivt, oro(c), &
   
       ! soil information
         albsat, albdry, csol, porsl, phi0, bsw, &
         dkmg, dksatu, dkdry, hksati, &
   
       ! vegetation information
         z0m, displa, sqrtdi, &
         effcon, vmax25, slti, hlti, shti, hhti, &
         trda, trdm, trop, gradm, binter, extkn, &
         chil, refl, refs, tranl, trans, rootfr, &
#if (defined DGVM)
         pftpara, vegclass, summergreen, raingreen, &
         sla, lm_sapl, sm_sapl, hm_sapl, rm_sapl,& 
         t10min, lai_ind,&
         dphen,leafon,leafof,firelength,&
#if (defined IAPDGVM)
         afirefrac1, nfireg1,&
#endif
         litterag, &
         litterbg,cpool_fast,cpool_slow,k_fast_ave, &
         k_slow_ave,litter_decom_ave,fmicr, &
         ifpre,prec365,nind,lm_ind,sm_ind,agddtw, &
         hm_ind,rm_ind,tmomin20,agdd0,agdd,agdd20, &
         t_mo,t_mo_sum,t_mo_min,crownarea,htop,tsai, &
         fpcgrid,bm_inc,afmicr,annpsn,annpsnpot,tref10,&
         tref_sum,t10,assimn10,assimn_sum,an10,&
         anngpp,annfrmf,annfrms,annfrmr,annfrg,&
         cflux_litter_soil,cflux_litter_atmos,&
         nday,nyr,&
         Isf,Iss,Ksf,Kss,&
#if (defined IAPDGVM)
         area_xy, getig_in, getrhm_in,&
#endif

#endif   
#if (defined DyN)
         cton_soil, cton_pro, afcton_leaf,afcton_sap,afcton_root,&
         litter_leaf, litter_wood, litter_root, litter_repr,&
         litter_leaf_n, litter_wood_n, litter_root_n, litter_repr_n,&
         an_up, an_stress, soil_no3, soil_no2, soil_no, soil_n2o,&
         soil_n2, soil_nh4,&
#endif

       ! atmospheric forcing
         frl, sols, soll, solsd, solld, &
         pco2m, po2m, us, vs, tm, qm, &
         prc, prl, rain, snow, psrf, &
         rhoair, &
         hu, ht, hq, &

       ! model variables needed by restart run
         year, jday, msec, &
         z, dz, tss, wliq, wice, &
         tg, tlsun, tlsha, ldew, sag, scv, snowdp, &
         fveg, fsno, sigf, green, lai, sai, &
         coszen, albg, albv, alb, ssun, ssha, thermk, extkb, extkd, lakedepth, dz_lake, t_lake, lake_icefrac, savedtke1,&
#if (defined FHNP) && (defined FTF)
         frostdp, thawdp, frostdp0 ,D_temperature, N_time, frost_day, thaw_day, &
#endif
#if (defined FHNP) && (defined GWUSE)
       ! Note that GWUSE is still in the testing phase, don't define it.
       ! GWUSE may use the g_index wanglh 2018/11
         g_index,c_index,&

#endif
       ! fluxes
         taux, tauy, &
         fsena, fevpa, lfevpa, fsenl, fevpl, etr, &
         fseng, fevpg, olrg, fgrnd, trad, tref, qref, &
         rsur, rnof, rst, assim, respc, &
         parsun(1:npatch), parsha(1:npatch), sabvsun, sabvsha, sabg, sabvg(1:npatch), xerr, zerr, &
         snm, qsubl(1:npatch), &

       ! TUNABLE modle constants
         zlnd, zsno, csoilc, dewmx, wtfact, &
         capr, cnfac, ssi, wimp, pondmx, &
         smpmax, smpmin, trsmx0, tcrit, &

       ! time-varying vegetation from read-in file
#if (defined VEGDATA)
         lai_r(1:npatch), sai_r(1:npatch), green_r(1:npatch), fveg_r(1:npatch), htop_r(1:npatch), &
#endif

#ifdef DUST
       ! For dust emission added by Chenglai Wu: 04/27/2014
         dustsource, isltyp, soil_top_cat, vegcov_r, &
         dustemis_bin_1, dustemis_bin_2, dustemis_bin_3, dustemis_bin_4, dustemis_total, &
#endif

       ! additional variables required by coupling with regional model (WRF & RSM) 
         emis, z0ma, zol, rib, ustar, qstar, tstar, &
         u10m, v10m, f10m, fm, fh, fq )

#ifdef MYBUG
         if(mybug) then
            write(6,*) "col", itypwat, wxy_column
            write(6,*) "wxy_patch", wxy_patch(1:npatch)
            write(6,*) "lai", lai(1:npatch)
            write(6,*) "ivt", ivt(1:npatch)
            write(6,*) "fsena", fsena(1:npatch)
            write(6,*) "lfevpa", lfevpa(1:npatch)
            write(6,*) "olrg", olrg(1:npatch)
            write(6,*) "tg", tg
         endif
#endif

! ======================================================================
!  [4] Return required surface flux fields 
! ======================================================================

         fldv_col%tss (:,c) = tss (abs(maxsnl)+1:abs(maxsnl)+nl_soil)
         fldv_col%wliq(:,c) = wliq(abs(maxsnl)+1:abs(maxsnl)+nl_soil)
         fldv_col%wice(:,c) = wice(abs(maxsnl)+1:abs(maxsnl)+nl_soil)
         fldv_col%t_lake(:,c)       = t_lake
         fldv_col%lake_icefrac(:,c) = lake_icefrac

         fldv_col%tg    (c) = tg
         fldv_col%scv   (c) = scv
         fldv_col%snowdp(c) = snowdp
         fldv_col%fsno  (c) = fsno
#if (defined FHNP) && (defined FTF)
!liruichao add 
         fldv_col%frostdp     (c) = frostdp
         fldv_col%thawdp      (c) = thawdp
         fldv_col%frostdp0    (c) = frostdp0
         fldv_col%D_temperature  (c) = D_temperature
         fldv_col%N_time      (c) = N_time
         fldv_col%frost_day   (c) = frost_day
         fldv_col%thaw_day    (c) = thaw_day
!end
#endif
       ! forcing
         fldv_col%us    (c) = us
         fldv_col%vs    (c) = vs
         fldv_col%tm    (c) = tm
         fldv_col%qm    (c) = qm
         fldv_col%prc   (c) = prc
         fldv_col%prl   (c) = prl
         fldv_col%pbot  (c) = pbot
         fldv_col%frl   (c) = frl
         fldv_col%solar (c) = sols + soll + solsd + solld

         fldv_pft%rnet (p1:p2) = sabvg(1:npatch) + frl - olrg(1:npatch)
         fldv_pft%fmicr(p1:p2) = fmicr(:)
         fldv_pft%tlsun(p1:p2) = tlsun(:)
         fldv_pft%tlsha(p1:p2) = tlsha(:)
         fldv_pft%ldew (p1:p2) = ldew(:)
         fldv_pft%sigf (p1:p2) = sigf(:)
         fldv_pft%green(p1:p2) = green(:)
         fldv_pft%lai  (p1:p2) = lai(:)
         fldv_pft%sai  (p1:p2) = sai(:)
         fldv_pft%avsdr(p1:p2) = alb(1,1,:)
         fldv_pft%avsdf(p1:p2) = alb(1,2,:)
         fldv_pft%anidr(p1:p2) = alb(2,1,:)
         fldv_pft%anidf(p1:p2) = alb(2,2,:)
         fldv_pft%emis (p1:p2) = emis(:)
         fldv_pft%z0ma (p1:p2) = z0ma(:)
         fldv_pft%trad (p1:p2) = trad(:)
         fldv_pft%ustar(p1:p2) = ustar(:)
         fldv_pft%tstar(p1:p2) = tstar(:)
         fldv_pft%qstar(p1:p2) = qstar(:)
         fldv_pft%zol  (p1:p2) = zol(:)
         fldv_pft%rib  (p1:p2) = rib(:)
         fldv_pft%fm   (p1:p2) = fm(:)
         fldv_pft%fh   (p1:p2) = fh(:)
         fldv_pft%fq   (p1:p2) = fq(:)
         fldv_pft%tref (p1:p2) = tref(:)
         fldv_pft%qref (p1:p2) = qref(:)

       ! liquid water content of snow layer [lwsnl:kg/m2]
         snl = 0
         do j = 1, abs(maxsnl)
            if((wliq(j)+wice(j))>0.) snl = snl + 1
         end do

         lb = abs(maxsnl)-snl+1
         ub = abs(maxsnl)

         if (snl>0) then
            lwsnl = sum(wliq(lb:ub))
         else
            lwsnl = 0._r8
         endif

       ! snow age [agesno:day]
       ! time samples weighted by snow mass and accumulate
       ! finally divided by ave(scv)

       ! agesno = dtime*scv

       ! snow internal temperature [tsn:K]

         if (snl>0) then
            tsn = sum(tss(lb:ub)*dz(lb:ub))/sum(dz(lb:ub))
            nsnow = 1._r8
         else
            tsn = 0._r8
            nsnow = 0._r8
         endif

         lb = abs(maxsnl)+1
         ub = abs(maxsnl)+3

         mrsos = sum(wliq(lb:ub)) + sum(wice(lb:ub))

         lb = abs(maxsnl)+1
         ub = abs(maxsnl)+nl_soil

         mrso  = sum(wliq(lb:ub)) + sum(wice(lb:ub))

         mrfso = sum(wice(lb:ub))

         mrlsl = wice(lb:ub) + wliq(lb:ub)

         fldv_pft%qsubl(p1:p2) = qsubl(1:npatch)

         fldv_col%mrlsl  (:,c) = mrlsl(1:nl_soil)
         fldv_col%mrsos    (c) = mrsos
         fldv_col%mrso     (c) = mrso
         fldv_col%mrfso    (c) = mrfso
         fldv_col%lwsnl    (c) = lwsnl
         fldv_col%snm      (c) = snm
         fldv_col%tsn      (c) = tsn
         fldv_col%nsnow    (c) = nsnow
#if (defined IAPDGVM)
!====================================
!calculate wliq6mon for IAPDGVM
!===================================
call get_curr_date(yr, mon, day, sec)
     if(dlat.ge.0) then
        if(mon>=5.and.mon<=10) then
            wliq6mon = wliq6mon + wliq(6) + wliq(7) + wliq(8)
        end if
     else
        if(mon>=1.and.mon<=4.or.mon==11.or.mon==12) then
          wliq6mon = wliq6mon + wliq(6) + wliq(7) + wliq(8)
        end if
     end if
#endif
!================================


         p1 = p2 + 1

      ENDDO

 end subroutine CLMDRIVER
