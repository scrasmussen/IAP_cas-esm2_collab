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
   end type fldv_type

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
   integer :: numpatch                        ! local patches number

   integer :: numgrid_glob                    ! global grids number
   integer :: numpatch_glob                   ! global patches number

#if(defined DGVM)
   integer :: numcolumn
   integer :: numcolumn_glob
#endif

! 
! land model grid location info
!
   integer , pointer :: numlon(:)             ! longitude points for each latitude strip
   real(r8), pointer :: latixy(:,:)           ! latitude of grid cell (degrees)
   real(r8), pointer :: longxy(:,:)           ! longitude of grid cell (degrees)
   real(r8), pointer :: area(:,:)             ! grid cell area (km**2)
   real(r8), pointer :: landarea              ! total land area for all gridcells (km^2)
   real(r8), pointer :: lats(:)               ! grid cell latitude, southern edge (degrees)
   real(r8), pointer :: lonw(:,:)             ! grid cell longitude, western edge (degrees)
!
! fractional land and mask
!
   real(r8), pointer :: landfrac(:,:)         ! fractional land
   integer , pointer :: landmask(:,:)         ! land mask: 1 = land. 0 = ocean
!
! patch and grid info
!
   real(r8), pointer :: wxy_patch(:)          ! patch weight

   integer , pointer :: ixy_patch_glob(:)     ! longitude index of patch
   integer , pointer :: jxy_patch_glob(:)     ! latitude index of patch
   real(r8), pointer :: wxy_patch_glob(:)     ! weight of patch

   integer,  pointer :: itypwat_glob(:)

#if(defined DGVM)
   real(r8), pointer :: wxy_column(:)         ! patch weight

   integer , pointer :: ixy_column_glob(:)    ! longitude index of patch
   integer , pointer :: jxy_column_glob(:)    ! latitude index of patch
   real(r8), pointer :: wxy_column_glob(:)    ! weight of patch
#endif
!
! model variables
!
   real(r8), pointer :: forc(:,:)             ! forcing variables
!  real(r8), pointer :: fcon(:,:)             ! time constant variables
!  real(r8), pointer :: fvar(:,:)             ! time varying variables

   real(r8), pointer :: oro(:)                ! ocean(0)/seaice(2)/ flag
   integer , pointer :: itypwat(:)            ! land water type

#if(defined DGVM)
   integer,  pointer :: numcolumn_lat(:)      ! number of columnes of grids at lon. strip
   real(r8), pointer :: fcon_col(:,:)         ! time constant variables
   real(r8), pointer :: fvar_col(:,:)         ! time varying variables
   real(r8), pointer :: fcon_pft(:,:)         ! time constant variables
   real(r8), pointer :: fvar_pft(:,:)         ! time varying variables
   integer,  pointer :: numpatch_lat(:)       ! number of patches of grids at lon. strip
   real(r8), pointer :: nep_residual(:)       ! annual residual NEP = acfire + aestabc
   logical           :: lnep_adjust = .false.
   real(r8), pointer :: fLitterSoil(:)        ! cflux_litter_soil
   real(r8), pointer :: fLitterAtmos(:)       ! cflux_litter_atmos

   real(r8), pointer :: Isf_pft(:)
   real(r8), pointer :: Iss_pft(:)
   real(r8), pointer :: Ksf_pft(:)
   real(r8), pointer :: Kss_pft(:)
#endif

   real(r8)  ftune(nftune)                    ! clm tunable constants

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

   interface colm_var_alloc1
      module procedure colm_var_alloc1
   end interface 

   interface colm_var_alloc2
      module procedure colm_var_alloc2
   end interface 

   interface colm_var_dealloc1
      module procedure colm_var_dealloc1
   end interface 

   interface colm_var_dealloc2
      module procedure colm_var_dealloc2
   end interface 

CONTAINS

   subroutine colm_var_alloc1

      allocate (numlon                (lat_points))
      allocate (lats                (lat_points+1))
      allocate (lonw     (lon_points+1,lat_points))
      allocate (area       (lon_points,lat_points))
      allocate (latixy     (lon_points,lat_points))
      allocate (longxy     (lon_points,lat_points))
      allocate (landfrac   (lon_points,lat_points))
      allocate (landmask   (lon_points,lat_points))

      allocate (ixy_patch_glob     (numpatch_glob))
      allocate (jxy_patch_glob     (numpatch_glob))
      allocate (wxy_patch_glob     (numpatch_glob))

#if(defined DGVM)
      allocate (ixy_column_glob   (numcolumn_glob))
      allocate (jxy_column_glob   (numcolumn_glob))
      allocate (wxy_column_glob   (numcolumn_glob))
#endif

      allocate (itypwat_glob      (numcolumn_glob))

   end subroutine colm_var_alloc1

   subroutine fldv_alloc

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

   end subroutine fldv_alloc

   subroutine fldv_dgvm_init

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

   subroutine fldv_init

! average from pft level
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

! average from column level
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

   end subroutine fldv_init

   subroutine fldv_dealloc

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

   end subroutine fldv_dealloc

   subroutine colm_var_alloc2

      use paramodel, only: nflai, nforc, nfvar_col, nfvar_pft, &
                           nfcon_col, nfcon_pft, maxpatch

      allocate (oro                (numcolumn))
      allocate (itypwat            (numcolumn))
      allocate (forc         (nforc,numcolumn))

      allocate (wxy_column         (numcolumn))
      allocate (fcon_col (nfcon_col,numcolumn))
      allocate (fvar_col (nfvar_col,numcolumn))

      allocate (wxy_patch           (numpatch))
      allocate (fcon_pft  (nfcon_pft,numpatch))
      allocate (fvar_pft  (nfvar_pft,numpatch))

#if(defined VEGDATA)
      allocate (mlai             (12,numpatch))
      allocate (msai             (12,numpatch))
      allocate (mhtop            (12,numpatch))
#endif

      call fldv_alloc

#ifdef DGVM
      allocate (nep_residual         (numgrid))
      nep_residual(:) = 0.

      allocate (fLitterSoil         (numpatch))
      allocate (fLitterAtmos        (numpatch))

      allocate (Isf_pft             (numpatch))
      allocate (Iss_pft             (numpatch))
      allocate (Ksf_pft             (numpatch))
      allocate (Kss_pft             (numpatch))

      Isf_pft(:) = 0.
      Iss_pft(:) = 0.
      Ksf_pft(:) = 0.
      Kss_pft(:) = 0.
#endif

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

   end subroutine colm_var_alloc2

   subroutine colm_var_dealloc1

      deallocate (ixy_patch_glob)
      deallocate (jxy_patch_glob)
      deallocate (wxy_patch_glob)

#ifdef DGVM
      deallocate (ixy_column_glob)
      deallocate (jxy_column_glob)
      deallocate (wxy_column_glob)
#endif

      deallocate (itypwat_glob)

   end subroutine colm_var_dealloc1

   subroutine colm_var_dealloc2

      deallocate (numlon      )
      deallocate (lats        )
      deallocate (lonw        )
      deallocate (latixy      )
      deallocate (longxy      )
      deallocate (area        )
      deallocate (landmask    )
      deallocate (landfrac    )

      deallocate (oro         )
      deallocate (itypwat     )
      deallocate (wxy_patch   )
      deallocate (wxy_column  )
      deallocate (forc        )
      deallocate (fcon_col    )
      deallocate (fvar_col    )
      deallocate (fcon_pft    )
      deallocate (fvar_pft    )

#if(defined VEGDATA)
      deallocate (mlai        )
      deallocate (msai        )
      deallocate (mhtop       )
#endif

      call fldv_dealloc

#ifdef DGVM
      deallocate (nep_residual)
      deallocate (fLitterSoil )
      deallocate (fLitterAtmos)

      deallocate (Isf_pft     )
      deallocate (Iss_pft     )
      deallocate (Ksf_pft     )
      deallocate (Kss_pft     )
#endif

#ifndef COUP_CSM 
      deallocate (tair        )
      deallocate (qair        )
      deallocate (pres        )
      deallocate (rainc       )
      deallocate (rainl       )
      deallocate (windu       )
      deallocate (windv       )
      deallocate (dswrf       )
      deallocate (dlwrf       )
      deallocate (tair_z      )
      deallocate (qair_z      )
      deallocate (wind_z      )
#endif

      deallocate (lai         )
      deallocate (sai         )
      deallocate (fveg        )
      deallocate (green       )

   end subroutine colm_var_dealloc2

end module colm_varMod
