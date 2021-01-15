#include <define.h>

subroutine lpjreset(c,p1,p2,npatch,ivt_r,ivt_r_old,ifpre,ifpre_old,wt_patch,wt_patch_old,&
#if (defined IAPDGVM)                  
                    wliq6mon,&
#endif

                   prec365) 

   use paramodel, only: nl_soil, maxsnl
   use colm_varMod, only: ftune, pvar, cvar
   use precision
   implicit none

   integer , intent(in)    :: c
   integer , intent(in)    :: p1
   integer , intent(in)    :: p2
   integer , intent(in)    :: npatch
   integer , intent(in)    :: ivt_r(npatch)
   integer , intent(in)    :: ivt_r_old(npatch)
   real(r8), intent(in)    :: ifpre(npatch)
   real(r8), intent(in)    :: ifpre_old(npatch)
   real(r8), intent(in)    :: wt_patch(npatch)
   real(r8), intent(in)    :: wt_patch_old(npatch)
   real(r8), intent(inout) :: prec365
#ifdef IAPDGVM
   real(r8), intent(inout) :: wliq6mon
#endif


                       ! Main land surface variables 

   integer , pointer :: itypwat
   real(r8), pointer :: albsat(:)
   real(r8), pointer :: albdry(:)
   real(r8), pointer :: z   (:)          ! node depth [m]
   real(r8), pointer :: dz  (:)          ! interface depth [m]
   real(r8), pointer :: tss (:)          ! soil temperature [K]
   real(r8), pointer :: wliq(:)          ! liquid water in layers [kg/m2]
   real(r8), pointer :: wice(:)          ! ice lens in layers [kg/m2]
   real(r8), pointer :: tg               ! ground surface temperature [K]
   real(r8), pointer :: sag              ! non dimensional snow age [-]
   real(r8), pointer :: scv              ! snow cover, water equivalent [mm]
   real(r8), pointer :: snowdp           ! snow depth [meter]
   real(r8), pointer :: fsno             ! fraction of snow cover on ground
   real(r8), pointer :: coszen           ! cosine of solar zenith angle

   real(r8), pointer :: tlsun            ! sunlit leaf temperature [K]
   real(r8), pointer :: tlsha            ! shaded leaf temperature [K]
   real(r8), pointer :: ldew             ! depth of water on foliage [mm]

                        ! Vegetation dynamic parameters 
   real(r8), pointer :: fveg             ! fraction of vegetation cover
   real(r8), pointer :: sigf             ! fraction of veg cover, excluding snow-covered veg [-]
   real(r8), pointer :: green            ! leaf greenness
   real(r8), pointer :: lai              ! leaf area index
   real(r8), pointer :: sai              ! stem area index

   real(r8), pointer :: t10min           ! annual minimum of 10-day running mean (K)
   real(r8), pointer :: lai_ind          ! LAI per individual
   real(r8), pointer :: dphen            ! phenology [0 to 1]
   real(r8), pointer :: leafon           ! leafon days
   real(r8), pointer :: leafof           ! leafoff days
   real(r8), pointer :: firelength       ! fire season in days
#ifdef IAPDGVM
  real(r8), pointer :: afirefrac1          !IAPDGVM
  real(r8), pointer :: nfireg1          !IAPDGVM
#endif

   real(r8), pointer :: litterag         ! above ground litter
   real(r8), pointer :: litterbg         ! below ground litter
   real(r8), pointer :: cpool_fast(:)    ! fast carbon pool
   real(r8), pointer :: cpool_slow(:)    ! slow carbon pool
   real(r8), pointer :: k_fast_ave       ! decomposition rate
   real(r8), pointer :: k_slow_ave       ! decomposition rate
   real(r8), pointer :: litter_decom_ave ! decomposition rate
   real(r8), pointer :: fmicr            ! microbial respiration (umol CO2 /m**2 /s)
   real(r8), pointer :: nind             ! number of individuals (#/m**2)
   real(r8), pointer :: lm_ind           ! individual leaf mass
   real(r8), pointer :: sm_ind           ! individual sapwood mass
   real(r8), pointer :: hm_ind           ! individual heartwood mass
   real(r8), pointer :: rm_ind           ! individual root mass
   real(r8), pointer :: tmomin20         ! 20-yr running mean of tmomin
   real(r8), pointer :: agdd0            ! growing dgree days above 0
   real(r8), pointer :: agdd             ! growing dgree days above 5
   real(r8), pointer :: agddtw           ! growing dgree days above twmax
   real(r8), pointer :: agdd20           ! 20-yr running mean of agdd
   real(r8), pointer :: t_mo             ! 30-day mean temperature of 2m
   real(r8), pointer :: t_mo_sum         ! 30-day accumulated temperature of 2m
   real(r8), pointer :: t_mo_min         ! annual min of t_mo (Kelvin)
   real(r8), pointer :: crownarea        ! area that each individual tree takes up (m^2)
   real(r8), pointer :: htop             ! canopy top
   real(r8), pointer :: tsai             ! one-sided stem area index, no burying by snow
   real(r8), pointer :: fpcgrid          ! foliar projective cover on gridcell (fraction)
   real(r8), pointer :: bm_inc           ! biomass increment
   real(r8), pointer :: afmicr           ! annual microbial respiration
   real(r8), pointer :: annpsn           ! annual photosynthesis (umol CO2 /m**2)
   real(r8), pointer :: annpsnpot        ! annual potential photosynthesis (same units)
   real(r8), pointer :: tref10           ! 10-day averaged temperature at 2m
   real(r8), pointer :: tref_sum         ! sum of temperature in current day
   real(r8), pointer :: t10(:)           ! array ro record the 10 day temperature
   real(r8), pointer :: assimn10         ! 10-day averaged assimilation rate
   real(r8), pointer :: assimn_sum       ! sum of assimn of current day
   real(r8), pointer :: an10(:)          ! arry to record 10 day assimn
   real(r8), pointer :: anngpp           ! annual gpp
   real(r8), pointer :: annfrmf          ! annual frmf
   real(r8), pointer :: annfrms          ! annual frms
   real(r8), pointer :: annfrmr          ! annual frmr
   real(r8), pointer :: annfrg           ! annual frg
   real(r8), pointer :: turnover_ind     ! individual turnover biomass
!  real(r8), pointer :: litter_ag        ! above ground litter mass
!  real(r8), pointer :: litter_bg        ! below ground litter mass
   real(r8), pointer :: fpc_inc          ! fpc increase
   integer , pointer :: ivt_x            ! ivt
   real(r8), pointer :: pftpar(:)        ! 32 parameters of PFTs
   real(r8), pointer :: vegclass         ! 1.tree 2.shrub 3.grass 4.crop -1.others
   real(r8), pointer :: summergreen      ! 1. for summergreen; otherwise -1.
   real(r8), pointer :: raingreen        ! 1. for raingreen; otherwise -1.
   real(r8), pointer :: sla              ! sla
   real(r8), pointer :: lm_sapl          ! leafmass
   real(r8), pointer :: sm_sapl          ! sapwood mass
   real(r8), pointer :: hm_sapl          ! heartwood mass
   real(r8), pointer :: rm_sapl          ! rootmass

#if(defined DyN)
   real(r8), pointer :: litter_leaf      ! leaf-derived litter for PFT on modelled area basis (gC/m2)
   real(r8), pointer :: litter_wood      ! heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
   real(r8), pointer :: litter_root      ! fine root-derived litter for PFT on modelled area basis(gC/m2)
   real(r8), pointer :: litter_repr      ! litter derived from allocation to reproduction for PFT on modelled

   real(r8), pointer :: litter_leaf_n    ! leaf-derived N litter for PFT on modelled area basis (gN/m2)
   real(r8), pointer :: litter_wood_n    ! heart&sapwood-derived N litter for PFT on modelled area basis(gN/m2)
   real(r8), pointer :: litter_root_n    ! fine root-derived N litter for PFT on modelled area basis (gN/m2)
   real(r8), pointer :: litter_repr_n    ! litter derived from allocation to reproduction N for PFT on modelled
                                         ! area basis (gN/m2)
   real(r8), pointer :: afcton_leaf      ! annual floating leaf C:N ratio
   real(r8), pointer :: afcton_root      ! annual floating root C:N ratio
   real(r8), pointer :: afcton_sap       ! annual floating sapwood C:N ratio
   real(r8), pointer :: lm_ind_n         ! individual leaf nitrogen mass
   real(r8), pointer :: sm_ind_n         ! individual sapwood nitrogen mass
   real(r8), pointer :: hm_ind_n         ! individual heartwood nitrogen mass
   real(r8), pointer :: rm_ind_n         ! individual root nitrogen mass
                                         ! gN/m2 veget'd area for each pft
   real(r8), pointer :: an_up            ! annual plant nitrogen uptake(gN/m2 vegt'd area)
   real(r8), pointer :: an_stress        ! annual nitrogen stress for production 
#endif
                       ! Radiation  related (albedoes)
   real(r8), pointer :: albg (:,:)       ! albedo, ground [-]
   real(r8), pointer :: albv (:,:)       ! albedo, vegetation [-]
   real(r8), pointer :: alb  (:,:)       ! averaged albedo [-]
   real(r8), pointer :: ssun (:,:)       ! sunlit canopy absorption for solar radiation (0-1)
   real(r8), pointer :: ssha (:,:)       ! shaded canopy absorption for solar radiation (0-1)
   real(r8), pointer :: thermk           ! canopy gap fraction for tir radiation
   real(r8), pointer :: extkb            ! (k, g(mu)/mu) direct solar extinction coefficient
   real(r8), pointer :: extkd            ! diffuse and scattered diffuse PAR extinction coefficient

                        ! Additional variables required by reginal model (WRF & RSM) 
   real(r8), pointer :: trad             ! radiative temperature of surface [K]
   real(r8), pointer :: tref             ! 2 m height air temperature [kelvin]
   real(r8), pointer :: qref             ! 2 m height air specific humidity
   real(r8), pointer :: rst              ! canopy stomatal resistance (s/m)

   real(r8), pointer :: emis             ! averaged bulk surface emissivity
   real(r8), pointer :: z0ma             ! effective roughness [m]
   real(r8), pointer :: zol              ! dimensionless height (z/L) used in Monin-Obukhov theory
   real(r8), pointer :: rib              ! bulk Richardson number in surface layer
   real(r8), pointer :: ustar            ! u* in similarity theory [m/s]
   real(r8), pointer :: qstar            ! q* in similarity theory [kg/kg]
   real(r8), pointer :: tstar            ! t* in similarity theory [K]
   real(r8), pointer :: fm               ! integral of profile function for momentum
   real(r8), pointer :: fh               ! integral of profile function for heat
   real(r8), pointer :: fq               ! integral of profile function for moisture

                        ! vegetation time invariant variables
   real(r8), pointer :: z0m_r            ! aerodynamic roughness length [m]
   real(r8), pointer :: chil             ! leaf angle distribution factor
   real(r8), pointer :: refl(:,:)        ! leaf reflectance (iw=iband, il=life and dead)
   real(r8), pointer :: refs(:,:)        ! stem reflectance (iw=iband, il=life and dead)
   real(r8), pointer :: tranl(:,:)       ! leaf transmittance (iw=iband, il=life and dead)
   real(r8), pointer :: trans(:,:)       ! leaf transmittance (iw=iband, il=life and dead)

! -----------------------------------------------------------------
! Local declaration
! -----------------------------------------------------------------

   real(r8) :: ssw                       ! water volumetric content of soil surface layer [m3/m3]
   real(r8) :: wt                        ! fraction of vegetation buried (covered) by snow [-]
   real(r8) :: z0m                       ! aerodynamic roughness length [m] zhq: 07/27/2010. z0m vary with htop
   real(r8) :: zlnd
   integer lc,uc,lb,ub,jm                ! column/pft indices

   integer ivt, ivt_old                  ! land cover type  

   integer pfrom(1), p, px

   real(r8) wt_tmp(npatch), vegwt


!  vegwt = sum(wt_patch(1:npatch))
!  if(vegwt<1.0E-6) return

   wt_tmp = wt_patch_old-wt_patch

   pfrom = maxloc(wt_tmp)
   if (.not.(pfrom(1)>=1 .and. pfrom(1)<=npatch)) then
      print *, 'pfrom', pfrom
      print *, wt_patch_old
      print *, wt_patch
      print *, 'error pfrom'
      return
   end if

   pfrom = pfrom + p1 - 1

   itypwat   => cvar%itypwat(c)
   albsat    => cvar%albsat(1:2,c)
   albdry    => cvar%albdry(1:2,c)

   z         => cvar%z(maxsnl+1:nl_soil,c)
   dz        => cvar%dz(maxsnl+1:nl_soil,c)
   tss       => cvar%tss(maxsnl+1:nl_soil,c)
   wliq      => cvar%wliq(maxsnl+1:nl_soil,c)
   wice      => cvar%wice(maxsnl+1:nl_soil,c)
   tg        => cvar%tg(c)
   sag       => cvar%sag(c)
   scv       => cvar%scv(c)
   snowdp    => cvar%snowdp(c)
   fsno      => cvar%fsno(c)
   coszen    => cvar%coszen(c)

   ssw = min(1.,1.e-3*wliq(6)/dz(6))

   do p = p1, p2

      px = p-p1+1

      ivt     = ivt_r(px)
      ivt_old = ivt_r_old(px)

    ! if(ifpre(px).lt.0.) cycle       zhq: 07/27/2010. move below.   

      ivt_x              => pvar%ivt(p)
      z0m_r              => pvar%z0m(p)
      chil               => pvar%chil(p)
      refl               => pvar%refl(1:2,1:2,p)
      refs               => pvar%refs(1:2,1:2,p)
      tranl              => pvar%tranl(1:2,1:2,p)
      trans              => pvar%trans(1:2,1:2,p)

      tlsun              => pvar%tlsun(p)
      tlsha              => pvar%tlsha(p)
      ldew               => pvar%ldew(p)
      fveg               => pvar%fveg(p)
      sigf               => pvar%sigf(p)
      green              => pvar%green(p)
      lai                => pvar%lai(p)
      sai                => pvar%sai(p)
      albg               => pvar%albg(1:2,1:2,p)
      albv               => pvar%albv(1:2,1:2,p)
      alb                => pvar%alb(1:2,1:2,p)
      ssun               => pvar%ssun(1:2,1:2,p)
      ssha               => pvar%ssha(1:2,1:2,p)
      thermk             => pvar%thermk(p)
      extkb              => pvar%extkb(p)
      extkd              => pvar%extkd(p)

                 ! Additional variables required by reginal model (WRF & RSM) 
      trad               => pvar%trad(p)
      tref               => pvar%tref(p)
      qref               => pvar%qref(p)
      rst                => pvar%rst(p)
      emis               => pvar%emis(p)
      z0ma               => pvar%z0ma(p)
      zol                => pvar%zol(p)
      rib                => pvar%rib(p)
      ustar              => pvar%ustar(p)
      qstar              => pvar%qstar(p)
      tstar              => pvar%tstar(p)
      fm                 => pvar%fm(p)
      fh                 => pvar%fh(p)
      fq                 => pvar%fq(p)

      t10min             => pvar%t10min(p)
      lai_ind            => pvar%lai_ind(p)
      dphen              => pvar%dphen(p)
      leafon             => pvar%leafon(p)
      leafof             => pvar%leafof(p)
      firelength         => pvar%firelength(p)
      litterag           => pvar%litterag(p)
      litterbg           => pvar%litterbg(p)
      cpool_fast         => pvar%cpool_fast(1:nl_soil,p)
      cpool_slow         => pvar%cpool_slow(1:nl_soil,p)
      k_fast_ave         => pvar%k_fast_ave(p)
      k_slow_ave         => pvar%k_slow_ave(p)
      litter_decom_ave   => pvar%litter_decom_ave(p)
      fmicr              => pvar%fmicr(p)
      nind               => pvar%nind(p)
      lm_ind             => pvar%lm_ind(p)
      sm_ind             => pvar%sm_ind(p)
      hm_ind             => pvar%hm_ind(p)
      rm_ind             => pvar%rm_ind(p)
      tmomin20           => pvar%tmomin20(p)
      agdd0              => pvar%agdd0(p)
      agdd               => pvar%agdd(p)
      agdd20             => pvar%agdd20(p)
      t_mo_min           => pvar%t_mo_min(p)
      crownarea          => pvar%crownarea(p)
      htop               => pvar%htop(p)
      tsai               => pvar%tsai(p)
      fpcgrid            => pvar%fpcgrid(p)
      bm_inc             => pvar%bm_inc(p)
      afmicr             => pvar%afmicr(p)
      annpsn             => pvar%annpsn(p)
      annpsnpot          => pvar%annpsnpot(p)
      tref10             => pvar%tref10(p)
      tref_sum           => pvar%tref_sum(p)
      t10                => pvar%t10(1:10,p)
      assimn10           => pvar%assimn10(p)
      assimn_sum         => pvar%assimn_sum(p)
      an10               => pvar%an10(1:10,p)
      turnover_ind       => pvar%turnover_ind(p)
      fpc_inc            => pvar%fpc_inc(p)
      agddtw             => pvar%agddtw(p)
      t_mo               => pvar%t_mo(p)
      t_mo_sum           => pvar%t_mo_sum(p)
      anngpp             => pvar%anngpp(p)
      annfrmf            => pvar%annfrmf(p)
      annfrms            => pvar%annfrms(p)
      annfrmr            => pvar%annfrmr(p)
      annfrg             => pvar%annfrg(p)
#ifdef IAPDGVM
      afirefrac1         => pvar%afirefrac1(p)
      nfireg1            => pvar%nfireg1(p)
#endif

#if(defined DyN)
      litter_leaf        => pvar%litter_leaf(p)
      litter_wood        => pvar%litter_wood(p)
      litter_root        => pvar%litter_root(p)
      litter_repr        => pvar%litter_repr(p)
      litter_leaf_n      => pvar%litter_leaf_n(p)
      litter_wood_n      => pvar%litter_wood_n(p)
      litter_root_n      => pvar%litter_root_n(p)
      litter_repr_n      => pvar%litter_repr_n(p)
      afcton_leaf        => pvar%afcton_leaf(p)
      afcton_sap         => pvar%afcton_sap(p)
      afcton_root        => pvar%afcton_root(p)
      lm_ind_n           => pvar%lm_ind_n(p)
      sm_ind_n           => pvar%sm_ind_n(p)
      hm_ind_n           => pvar%hm_ind_n(p)
      rm_ind_n           => pvar%rm_ind_n(p)
      an_up              => pvar%an_up(p)
      an_stress          => pvar%an_stress(p)
#endif

!     print*,'lpjreset->',p,tlsun,tlsha 
! ======================================================================

! reset accumulated variables on pft level
      annpsn = 0.
      annpsnpot = 0.
      firelength = 0.
#ifdef IAPDGVM
      afirefrac1 = 0.
      nfireg1 = 0.
#endif

      bm_inc = 0.
      afmicr = 0.
      agdd0  = 0.
      agdd   = 0.
      agddtw = 0.
      t_mo_min = 1.0E36
      anngpp  = 0.
      annfrmf = 0.
      annfrms = 0.
      annfrmr = 0.
      annfrg  = 0.
#if(defined DyN)
      an_up=0.
      an_stress=0.
#endif

! for new pfts established, zhq. 10/09/2009
      if(ifpre(px).lt.0.) cycle       

      if(ivt.ne.ivt_x) then
         print *, ivt, ivt_x
         stop 'lpjreset, ivt check failed'
      end if

      if(ifpre(px).gt.0. .and. ifpre_old(px).lt.0.) then

         pvar%tmomin20(p) = pvar%tmomin20(pfrom(1))
         pvar%agdd20(p)   = pvar%agdd20(pfrom(1))
         pvar%t_mo_sum(p) = pvar%t_mo_sum(pfrom(1))

       ! update albedo for new pfts for the first timestep of next year, zhq. 11/11/2009

         lai = max(lai, 0.05)
         sai = lai * 0.25

         z0m = max(z0m_r*htop, 0.01)
         zlnd = ftune(1)
 
         call snowfraction (itypwat,fveg,z0m,zlnd,snowdp,scv,wt,sigf,fsno)

         call albland (itypwat,albsat,albdry,chil,refl,refs,tranl,trans,&
                       fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,tg,&
                       alb,albg,albv,ssun,ssha,thermk,extkb,extkd)
      end if

!! reset accumulated variables on pft level
!      annpsn = 0.
!      annpsnpot = 0.
!      firelength = 0.
!      bm_inc = 0.
!      afmicr = 0.
!      agdd0 = 0.
!      agdd = 0.
!      t_mo_min = 1.0E36  
!#if(defined DyN)
!      an_up=0.
!      an_stress=0.
!#endif
   end do

! reset accumulated variables on column level
      prec365 = 0.  
#ifdef IAPDGVM
      wliq6mon = 0.
#endif

end subroutine lpjreset
