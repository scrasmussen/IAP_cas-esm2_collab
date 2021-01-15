
#include <define.h>

 subroutine iniTimeVar(nl_soil,maxsnl,n_pft,itypwat,ivt&
                      ,porsl,albsat,albdry,z0m,chil,refl,refs,tranl,trans&
                      ,rootfr,z,dz,tss,wliq,wice&
                      ,tg,tlsun,tlsha,ldew,sag,scv&   
                      ,snowdp,t_lake,lake_icefrac,savedtke1,fveg,fsno,sigf,green,lai,sai,coszen,zlnd&
                      ,albg,albv,alb,ssun,ssha,thermk,extkb,extkd&
                      ,trad,tref,qref,rst,emis,z0ma,zol,rib&
                      ,ustar,qstar,tstar,fm,fh,fq&
#if(defined SOILINI)
                      ,nl_soil_ini,soil_z,soil_t,soil_w,snow_d &
#endif
#if(defined DGVM)
                      ,wxy_column,wxy_patch &
                      ,t10min,lai_ind,dphen &
                      ,leafon,leafof,firelength&
#if(defined IAPDGVM)
                      ,afirefrac1,nfireg1,wliq6mon &
#endif

                      ,litter_ag, litter_bg,cpool_fast,cpool_slow,k_fast_ave &
                      ,k_slow_ave,litter_decom_ave,fmicr &
                      ,ifpre,prec365,nind,lm_ind,sm_ind &
                      ,hm_ind,rm_ind,tmomin20,agdd0,agdd,agdd20,t_mo,t_mo_sum &
                      ,t_mo_min,crownarea,htop,tsai,fpcgrid &
                      ,bm_inc,afmicr,annpsn,annpsnpot,tref10 &
                      ,tref_sum,t10,assimn10,assimn_sum,an10,nday,nyr &
                      ,turnover_ind,fpc_inc,agddtw &
                      ,anngpp,annfrmf,annfrms,annfrmr,annfrg &
#endif
#if(defined DyN)
                      ,afcton_leaf,afcton_root,afcton_sap &
                      ,litter_leaf,litter_root,litter_wood,litter_repr &
                      ,litter_leaf_n,litter_root_n,litter_wood_n,litter_repr_n &
                      ,lm_ind_n,sm_ind_n,hm_ind_n,rm_ind_n,an_up,an_stress &
                      ,soil_no3,soil_no2,soil_no,soil_n2o,soil_n2,soil_nh4 &
#endif

#if (defined FHNP) && (defined FTF)
!liruichao add
                      ,frostdp,thawdp,frostdp0,D_temperature,N_time,frost_day,thaw_day&
!end 
#endif
                      )

!=======================================================================
! Original author : Yongjiu Dai, 08/30/2002, revised February 2004
!=======================================================================

  use precision
  use phycon_module, only : tfrz, tkwat
  use paramodel, only : nl_lake
  implicit none 

  integer, INTENT(in) ::        &! 
        nl_soil,                &! soil layer number
        maxsnl,                 &! maximum snow layer number
        n_pft,                  &! number of pfts in a single column
        itypwat,                &! index for land cover type [-]
        ivt(1:n_pft)             ! land cover type

  real(r8), INTENT(in) ::       &!
        zlnd,                   &! roughness length for soil [m]
        coszen,                 &! cosine of solar zenith angle
        albsat(2),              &! wet soil [vis/nir] albedo [-]
        albdry(2),              &! dry soil [vis/nir] albedo [-]
        z0m     (1:n_pft),      &! aerodynamic roughness length [m]
        chil    (1:n_pft),      &! leaf angle distribution factor
        refl(2,2,1:n_pft),      &! leaf reflectance (iw=iband, il=life and dead)
        refs(2,2,1:n_pft),      &! stem reflectance (iw=iband, il=life and dead)
        tranl(2,2,1:n_pft),     &! leaf transmittance (iw=iband, il=life and dead)
        trans(2,2,1:n_pft),     &! stem transmittance (iw=iband, il=life and dead)
        porsl (1:nl_soil),      &! porosity of soil
        rootfr(1:nl_soil,1:n_pft)! fraction of roots in each soil layer

  real(r8), INTENT(out)  ::     &
        lai  (1:n_pft),         &! leaf area index
        sai  (1:n_pft),         &! stem area index
        fveg (1:n_pft),         &! fraction of vegetation cover  ! changed by zhq. feb.10.09
        green(1:n_pft)           ! leaf greenness 

#if(defined SOILINI)
  integer, INTENT(in) :: nl_soil_ini
  real(r8), INTENT(in) ::       &!
        soil_z(nl_soil_ini),    &! soil layer depth for initial (m)
        soil_t(nl_soil_ini),    &! soil temperature from initial file (K)
        soil_w(nl_soil_ini),    &! soil wetness from initial file (-)
        snow_d                   ! snow depth (m)
#endif

  real(r8), INTENT(inout) ::    &!
        z (maxsnl+1:nl_soil),   &! node depth [m]
        dz(maxsnl+1:nl_soil)     ! interface depth [m]

#if(defined DGVM)
  real(r8), INTENT(inout) ::    &!
        nday,                   &! counting the model days
        nyr,                    &! counting the model years
        prec365,                &! yearly running mean of precipitation(mm/s)
        wxy_column,             &! colomn fraction
        wxy_patch (1:n_pft),    &! pft fraction
        t10min    (1:n_pft),    &! annual minimum of 10-day running mean (K)
        lai_ind   (1:n_pft),    &! LAI per individual
        dphen     (1:n_pft),    &! phenology [0 to 1]
        leafon    (1:n_pft),    &! leafon days
        leafof    (1:n_pft),    &! leafoff days
        firelength(1:n_pft),    &! fire season in days
#if(defined IAPDGVM)
        afirefrac1(1:n_pft),    &! fire  for IAPDGVM
        nfireg1   (1:n_pft),    &! fire  for IAPDGVM
        wliq6mon,               &! wliq for the first three layers
#endif
        litter_ag (1:n_pft),    &! above ground litter
        litter_bg (1:n_pft),    &! below ground litter
        cpool_fast(nl_soil,1:n_pft),    &! fast carbon pool
        cpool_slow(nl_soil,1:n_pft),    &! slow carbon pool
        k_fast_ave(1:n_pft),    &! decomposition rate
        k_slow_ave(1:n_pft),    &! decomposition rate
        litter_decom_ave(1:n_pft),&! decomposition rate
        fmicr     (1:n_pft),    &! microbial respiration (umol CO2 /m**2 /s)
        nind      (1:n_pft),    &! number of individuals (#/m**2)
        lm_ind    (1:n_pft),    &! individual leaf mass
        sm_ind    (1:n_pft),    &! individual sapwood mass
        hm_ind    (1:n_pft),    &! individual heartwood mass
        rm_ind    (1:n_pft),    &! individual root mass
        tmomin20  (1:n_pft),    &! 20-yr running mean of tmomin
        agdd0     (1:n_pft),    &! growing dgree days above 0
        agdd      (1:n_pft),    &! growing dgree days above 5
        agddtw    (1:n_pft),    &! growing dgree days above twmax
        agdd20    (1:n_pft),    &! 20-yr running mean of agdd
        t_mo      (1:n_pft),    &! 30-day mean temperature of 2m (K)
        t_mo_sum  (1:n_pft),    &! 30-day accumulated temperature of 2m (K)
        t_mo_min  (1:n_pft),    &! annual min of t_mo (Kelvin)
        crownarea (1:n_pft),    &! area that each individual tree takes up (m^2)
        htop      (1:n_pft),    &! canopy top
        tsai      (1:n_pft),    &! one-sided stem area index, no burying by snow
        fpcgrid   (1:n_pft),    &! foliar projective cover on gridcell (fraction)
        bm_inc    (1:n_pft),    &! biomass increment
        afmicr    (1:n_pft),    &! annual microbial respiration
        annpsn    (1:n_pft),    &! annual photosynthesis (umol CO2 /m**2)
        annpsnpot (1:n_pft),    &! annual potential photosynthesis (same units)
        tref10    (1:n_pft),    &! 10-day averaged temperature at 2m
        tref_sum  (1:n_pft),    &! sum of tref in current day
        t10    (10,1:n_pft),    &! arry to record the 10 day temperature
        assimn10  (1:n_pft),    &! 10-day averaged assimilation rate
        assimn_sum(1:n_pft),    &! sum of assimn of current day
        an10   (10,1:n_pft),    &! arry to record 10 day assimn
        turnover_ind(1:n_pft),  &! individual turnover biomass
        fpc_inc   (1:n_pft),    &! fpc increase
        ifpre     (1:n_pft),    &! whether PFT present in patch 1.present; -1.not present
        anngpp    (1:n_pft),    &! whether PFT present in patch
        annfrmf   (1:n_pft),    &! whether PFT present in patch
        annfrms   (1:n_pft),    &! whether PFT present in patch
        annfrmr   (1:n_pft),    &! whether PFT present in patch
        annfrg    (1:n_pft)      ! whether PFT present in patch
#endif
#if(defined DyN)
  real(r8), INTENT(out) ::   &!
        litter_leaf(1:n_pft)   ,&! leaf-derived litter for PFT on modelled area basis (gC/m2)
        litter_wood(1:n_pft)   ,&! heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
        litter_root(1:n_pft)   ,&! fine root-derived litter for PFT on modelled area basis(gC/m2)
        litter_repr(1:n_pft)   ,&! litter derived from allocation to reproduction for PFT on modelled
   
        litter_leaf_n(1:n_pft) ,&! leaf-derived N litter for PFT on modelled area basis (gN/m2)
        litter_wood_n(1:n_pft) ,&! heart&sapwood-derived N litter for PFT on modelled area basis(gN/m2)
        litter_root_n(1:n_pft) ,&! fine root-derived N litter for PFT on modelled area basis (gN/m2)
        litter_repr_n(1:n_pft) ,&! litter derived from allocation to reproduction N for PFT on modelled
                                 ! area basis (gN/m2)
        afcton_leaf(1:n_pft)   ,&! annual floating leaf C:N ratio
        afcton_root(1:n_pft)   ,&! annual floating root C:N ratio
        afcton_sap(1:n_pft)    ,&! annual floating sapwood C:N ratio
        lm_ind_n(1:n_pft)      ,&! individual leaf nitrogen mass
        sm_ind_n(1:n_pft)      ,&! individual sapwood nitrogen mass
        hm_ind_n(1:n_pft)      ,&! individual heartwood nitrogen mass
        rm_ind_n(1:n_pft)      ,&! individual root nitrogen mass
                                 ! gN/m2 veget'd area for each pft
        an_up(1:n_pft)         ,&! annual plant nitrogen uptake(gN/m2 vegt'd area)
        an_stress(1:n_pft)     ,&! annual plant nitrogen stress(-)
   
        soil_no3               ,&!
        soil_no2               ,&!
        soil_no                ,&!
        soil_n2o               ,&!
        soil_n2                ,&!
        soil_nh4                 !
#endif

  real(r8), INTENT(out) ::      &
        tss (maxsnl+1:nl_soil), &! soil temperature [K]
        wliq(maxsnl+1:nl_soil), &! liquid water in layers [kg/m2]
        wice(maxsnl+1:nl_soil), &! ice lens in layers [kg/m2]
        tg,                     &! ground surface temperature [K]
        tlsun(1:n_pft),         &! sunlit leaf temperature [K]
        tlsha(1:n_pft),         &! shaded leaf temperature [K]
        ldew(1:n_pft),          &! depth of water on foliage [mm]
        sag,                    &! non dimensional snow age [-]
        scv,                    &! snow cover, water equivalent [mm]
        snowdp,                 &! snow depth [meter]
        fsno,                   &! fraction of snow cover on ground
        sigf(1:n_pft),          &! fraction of veg cover, excluding snow-covered veg [-]
        t_lake(1:nl_lake),      &! lake layer temperature [K]
        lake_icefrac(1:nl_lake),&! lake mass fraction of lake layer that is frozen
        savedtke1,              &! weinan

        albg(2,2,1:n_pft),      &! albedo, ground [-]
        albv(2,2,1:n_pft),      &! albedo, vegetation [-]
        alb (2,2,1:n_pft),      &! averaged albedo [-]
        ssun(2,2,1:n_pft),      &! sunlit canopy absorption for solar radiation
        ssha(2,2,1:n_pft),      &! shaded canopy absorption for solar radiation
        thermk(1:n_pft),        &! canopy gap fraction for tir radiation
        extkb(1:n_pft),         &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd(1:n_pft),         &! diffuse and scattered diffuse PAR extinction coefficient

! Additional variables required by reginal model (WRF & RSM) 
        trad(1:n_pft),          &! radiative temperature of surface [K]
        tref(1:n_pft),          &! 2 m height air temperature [kelvin]
        qref(1:n_pft),          &! 2 m height air specific humidity
        rst(1:n_pft),           &! canopy stomatal resistance (s/m)
        emis(1:n_pft),          &! averaged bulk surface emissivity
        z0ma(1:n_pft),          &! effective roughness [m]
        zol(1:n_pft),           &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib(1:n_pft),           &! bulk Richardson number in surface layer
        ustar(1:n_pft),         &! u* in similarity theory [m/s]
        qstar(1:n_pft),         &! q* in similarity theory [kg/kg]
        tstar(1:n_pft),         &! t* in similarity theory [K]
        fm(1:n_pft),            &! integral of profile function for momentum
        fh(1:n_pft),            &! integral of profile function for heat
        fq(1:n_pft)              ! integral of profile function for moisture


        integer j, snl, ivtx
        integer c, p             ! indices of column and pft
        real(r8) wet(nl_soil), wt(1:n_pft), ssw, oro, rhosno_ini, a 
        real(r8),parameter :: T0 = 273.16
        real(r8), parameter :: reinickerp = 1.6 !parameter in allometric equation
        real(r8), parameter :: allom1 = 100.0   !parameters in allometric
        real(r8), parameter :: allom2 =  40.0
        real(r8), parameter :: allom3 =   0.5

!Refer to CLM-DGVM 400-year simulation results.

        real(r8), dimension(21),parameter :: &
                   nind_ini =(/  0.008,   0.014,   0.010,   0.060,   0.007,   0.065,   0.040 &
                              ,  0.093,   0.100,   0.122,   0.100,   1.000,   1.000,   1.000 &
                              ,  1.000,   1.000,   0.   ,   0.   ,   0.   ,   0.   ,   0.   /)    
        real(r8), dimension(21),parameter :: &
                  lm_ind_ini=(/ 1030.7,  1292.8,  1292.8,  1665.7,   973.6,   521.8,   626.0 &
                              ,  442.2,   440.0,   440.0,   440.0,   120.9,   135.7,   154.6 &
                              ,  200.0,   200.0,     0. ,     0. ,     0. ,     0. ,     0. /)
        real(r8), dimension(21), parameter :: &
                  rm_ind_ini=(/ 1154.0,  1546.2,  1546.2,  1478.4,  1354.2,   710.3,   759.2 &
                              ,  526.4,   526.4,   526.4,   526.4,   161.2,   181.0,   206.1 &
                              ,  150.0,   150.5,     0. ,     0. ,     0. ,     0. ,     0. /)
        real(r8), dimension(21), parameter :: &
                  sm_ind_ini=(/ 8279.2, 10892.2, 10892.2, 11274.5, 12785.9,  8575.2, 10537.1 &
                              , 6943.8,  6943.8,  6943.8,  6943.8,     0. ,     0. ,     0.  &
                              ,    0. ,     0. ,     0. ,     0. ,     0. ,     0. ,     0. /)
        real(r8), dimension(21), parameter :: &
                  hm_ind_ini=(/27853.4, 39760.8, 39760.8, 39786.6, 47792.2, 32130.6, 38717.1 &
                              ,21758.8, 21758.8, 21758.8, 21758.8,     0. ,     0. ,     0.  &
                              ,    0. ,     0. ,     0. ,     0. ,     0. ,     0. ,     0. /)  
        real(r8), dimension(21),parameter :: &  ! added by jidy@18-Mar-2014 <rough estimation based on lm_ind_ini>
                 lai_ind_ini=(/   4.12,    5.17,    5.17,    6.66,    3.89,    2.08,     2.50&
                              ,   1.77,    1.76,    1.76,    1.76,    0.48,    0.54,     0.62&
                              ,   0.8 ,    0.8 ,    0.  ,    0.  ,    0.  ,    0.  ,     0. /)
        real(r8), dimension(21), parameter :: &
               litter_ag_ini=(/   0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,     0.  &
                              ,   0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,     0.  &
                              ,   0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,     0. /)
        real(r8), dimension(21), parameter :: &
               litter_bg_ini=(/   0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,     0.  &
                              ,   0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,     0.  &
                              ,   0.  ,    0.  ,    0.  ,    0.  ,    0.  ,    0.  ,     0. /)
        real(r8), dimension(21), parameter :: &
                      height=(/  10.  ,   10.  ,   10.  ,   10.  ,   10.  ,   10.  ,    10.  &
                              ,  10.  ,    5.  ,    5.  ,    5.  ,    1.  ,    1.  ,     1.  &
                              ,   1.  ,    1.  ,    1.  ,    0.  ,    0.  ,    0.  ,     0. /)
    
#if (defined FHNP) && (defined FTF)
!liruichao add
        real(r8), INTENT(out)   :: &
        frostdp0               ,&!frost depth
        frostdp                ,&!thaw deppth
        thawdp                 ,&!initial frost depth
        D_temperature          ,&!frost or thaw index
        N_time                 ,&!step counter
        frost_day              ,&!frost days
        thaw_day                 !thaw days
!end
#endif
!-----------------------------------------------------------------------

  if(itypwat <= 5)then ! land grid

! CREAT fraction of vegetation cover, greenness, leaf area index, stem index
     lai(:)=0.0; sai(:)=0.0; green(:)=0.0; fveg(:)=0.0
#if(defined SOILINI)
     do l = 1, nl_soil
        tss(l,i) = soil_t(nl_soil_ini,i)
     enddo
#else
     tss(1:) = 283.
#endif

! Read soil data
     rhosno_ini = 250.
#if(defined SOILINI)
     do j = 1, nl_soil
        call polint(soil_z,soil_t,nl_soil_ini,z(j),tss(j))
        call polint(soil_z,soil_w,nl_soil_ini,z(j),wet(j))
        a = min(soil_t(1),soil_t(2),soil_t(3))-5.
        tss(j) = max(tss(j), a)
        a = max(soil_t(1),soil_t(2),soil_t(3))+5.
        tss(j) = min(tss(j), a)

        a = min(soil_w(1),soil_w(2),soil_w(3))
        wet(j) = max(wet(j), a, 0.1)
        a = max(soil_w(1),soil_w(2),soil_w(3))
        wet(j) = min(wet(j), a, 0.5)

        if(tss(j).ge.tfrz)then
           wliq(j) = wet(j)*dz(j)*1000.
!          wliq(j) = porsl(j)*wet(j)*dz(j)*1000.
           wice(j) = 0.
        else
           wliq(j) = 0.
           wice(j) = wet(j)*dz(j)*1000.
!          wliq(j) = porsl(j)*wet(j)*dz(j)*1000.
        endif
     enddo

     snowdp = snow_d
     sag    = 0.
     scv    = snowdp*rhosno_ini

     call snowfraction (itypwat,fveg,z0m,zlnd,snowdp,scv,wt,sigf,fsno)
     call snow_ini (itypwat,maxsnl,snowdp,snl,z,dz)

     if(snl.lt.0)then
        do j = snl+1, 0
           tss(j) = min(tfrz-1., tss(1))
           wliq(j) = 0.
           wice(j) = dz(j)*rhosno_ini         ! m * kg m-3 = kg m-2
        enddo
     endif

     if(snl>maxsnl)then
        tss (maxsnl+1:snl) = -999.
        wice(maxsnl+1:snl) = 0.
        wliq(maxsnl+1:snl) = 0.
        z   (maxsnl+1:snl) = 0.
        dz  (maxsnl+1:snl) = 0.
     endif

     ldew  = 0.
     tlsun = tss(1)
     tlsha = tss(1)
     tg    = tss(1)
#else
! soil temperature and water content
     do j = 1, nl_soil
        if(itypwat==3)then ! land ice 
           tss(j) = 253.
           wliq(j) = 0.
           wice(j) = dz(j)*1000.
        else
           tss(j) = 283.
           wliq(j) = dz(j)*porsl(j)*1000.
           wice(j) = 0.
        endif
     enddo

! snow temperature and water content
     tss(maxsnl+1:0) = -999.
     wice(maxsnl+1:0) = 0.
     wliq(maxsnl+1:0) = 0.
     z (maxsnl+1:0) = 0.
     dz(maxsnl+1:0) = 0.
     fsno   = 0.
     scv    = 0.
     sag    = 0.
     snowdp = 0.
     tg     = tss(1)

     sigf(:)   = fveg(:)     
     wt(:)     = 0.
     ldew(:)   = 0.
     tlsun(:)  = tss(1)
     tlsha(:)  = tss(1)
#endif

! lake variables
     if(itypwat==4 .or. itypwat==5)then
        t_lake      (:) = 285.
        lake_icefrac(:) = 0.
        savedtke1       = tkwat
     end if

! surface albedo
     ssw = min(1.,1.e-3*wliq(1)/dz(1))
     do p = 1, n_pft
        call albland (itypwat,albsat,albdry,chil(p),refl(:,:,p),refs(:,:,p),tranl(:,:,p),trans(:,:,p),&
                      fveg(p),green(p),lai(p),sai(p),coszen,wt(p),fsno,scv,sag,ssw,tg,&
                      alb(:,:,p),albg(:,:,p),albv(:,:,p),ssun(:,:,p),ssha(:,:,p),thermk(p),extkb(p),extkd(p))
     enddo
      
#if(defined DGVM)
     nday     = 1.0
     nyr      = 1.0
     prec365  = 0.
#if(defined IAPDGVM)
     wliq6mon  = 0.
#endif

#if(defined DyN)
     soil_no3 = 0.0
     soil_no2 = 0.0
     soil_no  = 0.0
     soil_n2o = 0.0
     soil_n2  = 0.0
     soil_nh4 = 0.0
#endif

#if (defined FHNP) && (defined FTF)
!liruichao add
        frostdp0      = 0.0
        frostdp       = 0.0
        thawdp        = 0.0
        D_temperature = 0.0
        N_time        = 0.0
        frost_day     = 0.0
        thaw_day      = 0.0
!end  
#endif
     do p = 1,n_pft   
        ivtx = ivt(p)
      
        nind(p)      = nind_ini(ivtx)
        lm_ind(p)    = lm_ind_ini(ivtx)
        sm_ind(p)    = sm_ind_ini(ivtx)
        hm_ind(p)    = hm_ind_ini(ivtx)
        rm_ind(p)    = rm_ind_ini(ivtx)
        litter_ag(p) = litter_ag_ini(ivtx)
        litter_bg(p) = litter_bg_ini(ivtx)
        htop(p)      = height(ivtx)
        if(ivtx .le. 11) then
           crownarea(p) = 3.5 !add by Ming 
        else if (ivtx .le. 16) then
           crownarea(p) = 1.0
        end if
      ! lai_ind initial value set here follows lai_empirical  
      ! lai_ind(p)   = 0.              !********************
        lai_ind(p)   = lai_ind_ini(p)  ! added by jidy@18/Mar/2014
      ! fpcgrid(p)   = 0.9 * crownarea(p) * nind(p) ! commented by jidy@12/Sep/2014
        fpcgrid(p)   = wxy_patch(p)/wxy_column      ! added by jidy@12/Sep/2014
                                                    ! Currently in CLMMAIN:
                                                    ! wxy_patch(p) = ftpgrid(p)*wxy_column(c)
                                                    ! Change this in future with removing above line from CLMMAIN
      
        t10min(p)    = 1.0e+36
      
        ! updated in Phenology
        dphen(p) = 0.0
        leafon(p) = 0.0
        leafof(p) = 0.0
        ! updated before Phenology
        agdd0(p) = 0.0
        agdd(p) = 0.0
        agddtw(p) = 0.0
        agdd20(p) = 0.0
        tref10(p) = 0.0
        tref_sum(p) = 0.0
        t10(:,p) = tref(p)
        assimn10(p) = 0.0
        assimn_sum(p) = 0.0
        an10(:,p) = 0.0

        ! accumulated in FireSeason (must reset at end of every year)
        firelength(p) = 0.0
#if(defined IAPDGVM)
        afirefrac1(p) = 0.0
        nfireg1(p)    = 0.0
#endif
      
        ! used and updated in annual portion of LPJ
        ifpre(p)     = -1.
        tmomin20(p)  = T0 - 5. !initialize this way for Phenology code
        agdd20(p)    = 0.
        t_mo_min(p)  = 1.0e+36
        t_mo(p)      = 1.0e+36
        t_mo_sum(p)  = 0.

        tsai(p) = 0.
        !
        ! accumulated in Biogeochemistry and used/reset in annual portion of LPJ
        annpsn(p) = 0.0
        annpsnpot(p) = 0.0
        bm_inc(p) = 0.0
        afmicr(p) = 0.0
      
        anngpp(p)  = 0._r8
        annfrmf(p) = 0._r8
        annfrms(p) = 0._r8
        annfrmr(p) = 0._r8
        annfrg(p)  = 0._r8
      
        cpool_fast(:,p) = 0.0
        cpool_slow(:,p) = 0.0
        k_fast_ave(p) = 0.0
        k_slow_ave(p) = 0.0
        litter_decom_ave(p) = 0.0
        fmicr(p) = 0.0   !initialize b/c use in Biogeochemistry before LitterSOM
      
        ! updated in Nitrogen module. zhq.05/13/2010
#if(defined DyN)
        afcton_leaf(p) = 0.0
        afcton_root(p) = 0.0
        afcton_sap(p) = 0.0
        litter_leaf_n(p) = 0.0
        litter_root_n(p) = 0.0
        litter_wood_n(p) = 0.0
        litter_repr_n(p) = 0.0
        litter_leaf(p) = 0.0
        litter_root(p) = 0.0
        litter_wood(p) = 0.0
        litter_repr(p) = 0.0
        lm_ind_n(p) = 0.0
        rm_ind_n(p) = 0.0
        sm_ind_n(p) = 0.0
        hm_ind_n(p) = 0.0
        an_up(p) = 0.0
        an_stress(p) = 0.0
#endif
     end do
#endif

  else                 ! ocean grid

#if(defined SOILINI)
     tss(:) = soil_t(1)
     tlsun  = soil_t(1)
     tlsha  = soil_t(1)
     tg     = soil_t(1)
#else
     tss(:)   = 283.
     tg       = 283.
     tlsun(:) = 283.
     tlsha(:) = 283.
#endif
     wice(:) = 0.
     wliq(:) = 1000.
     z (maxsnl+1:0) = 0.
     dz(maxsnl+1:0) = 0.
     sigf(:)   = 0.
     fsno   = 0.
     ldew(:)   = 0.
     scv    = 0.
     sag    = 0.
     snowdp = 0.

     oro = 0
     do p = 1, n_pft
       call albocean (oro,scv,coszen,alb(:,:,p))
       albg(:,:,p) = alb(:,:,p)
       albv(:,:,p) = 0.0
       ssun(:,:,p) = 0.0
       ssha(:,:,p) = 0.0
       thermk(p) = 0.0
       extkb(p) = 0.0
       extkd(p) = 0.0
     enddo
  endif

! Additional variables required by reginal model (WRF & RSM)
! totally arbitrarily assigned here

  trad(:)  = tg      
  tref(:)  = tg      
  qref(:)  = 0.3     
  rst(:)   = 1.e36   
  emis(:)  = 1.0     
  z0ma(:)  = 0.01    
  zol(:)   = -1.0    
  rib(:)   = -0.1    
  ustar(:) = 0.25    
  qstar(:) = 0.001   
  tstar(:) = -1.5    
  fm(:)    = alog(30.)  
  fh(:)    = alog(30.)  
  fq(:)    = alog(30.)  

 end subroutine iniTimeVar
!-----------------------------------------------------------------------
! EOP


subroutine snow_ini(itypwat,maxsnl,snowdp,snl,z,dz)

! Snow spatial discretization initially

  use precision
  implicit none

  integer, intent(in) :: maxsnl  ! maximum of snow layers
  integer, intent(in) :: itypwat ! index for land cover type [-]
  real(r8), intent(in) :: snowdp ! snow depth [m]
  real(r8), intent(out) :: z (maxsnl+1:0) ! node depth [m]
  real(r8), intent(out) :: dz(maxsnl+1:0) ! layer thickness [m]
  integer, intent(out) :: snl ! number of snow layer
  real(r8) zi
  integer i
!-----------------------------------------------------------------------

  dz(:0) = 0.
  z(:0) = 0.
  snl = 0
  if(itypwat.le.3)then ! non water bodies

     if(snowdp.lt.0.01)then
        snl = 0
     else
        if(snowdp>=0.01 .and. snowdp<=0.03)then
           snl = -1
           dz(0)  = snowdp
        else if(snowdp>0.03 .and. snowdp<=0.04)then
           snl = -2
           dz(-1) = snowdp/2.
           dz( 0) = dz(-1)
        else if(snowdp>0.04 .and. snowdp<=0.07)then
           snl = -2
           dz(-1) = 0.02
           dz( 0) = snowdp - dz(-1)
        else if(snowdp>0.07 .and. snowdp<=0.12)then
           snl = -3
           dz(-2) = 0.02
           dz(-1) = (snowdp - 0.02)/2.
           dz( 0) = dz(-1)
        else if(snowdp>0.12 .and. snowdp<=0.18)then
           snl = -3
           dz(-2) = 0.02
           dz(-1) = 0.05
           dz( 0) = snowdp - dz(-2) - dz(-1)
        else if(snowdp>0.18 .and. snowdp<=0.29)then
           snl = -4
           dz(-3) = 0.02
           dz(-2) = 0.05
           dz(-1) = (snowdp - dz(-3) - dz(-2))/2.
           dz( 0) = dz(-1)
        else if(snowdp>0.29 .and. snowdp<=0.41)then
           snl = -4
           dz(-3) = 0.02
           dz(-2) = 0.05
           dz(-1) = 0.11
           dz( 0) = snowdp - dz(-3) - dz(-2) - dz(-1)
        else if(snowdp>0.41 .and. snowdp<=0.64)then
           snl = -5
           dz(-4) = 0.02
           dz(-3) = 0.05
           dz(-2) = 0.11
           dz(-1) = (snowdp - dz(-4) - dz(-3) - dz(-2))/2.
           dz( 0) = dz(-1)
        else if(snowdp>0.64)then
           snl = -5
           dz(-4) = 0.02
           dz(-3) = 0.05
           dz(-2) = 0.11
           dz(-1) = 0.23
           dz( 0) = snowdp - dz(-4) - dz(-3) - dz(-2) - dz(-1)
        endif

        zi = 0.
        do i = 0, snl+1, -1
           z(i) = zi - dz(i)/2.
           zi = -zi-dz(i)
        enddo
     endif

  endif

end subroutine snow_ini
!-----------------------------------------------------------------------
! EOP


subroutine polint(xa,ya,n,x,y)

! Given arrays xa and ya, each of length n, and gi
! value y, and an error estimate dy. If P (x) is the p
! P (xa(i)) = ya(i), i = 1, . . . , n, then the returned value
! (from: "Numerical Recipes")

  use precision
  implicit none
  integer n,NMAX
  real(r8) dy,x,y,xa(n),ya(n)
  parameter (NMAX=10)      !Largest anticipated val
  integer i,m,ns
  real(r8) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=abs(x-xa(1))

  do i=1,n       !Here we find the index ns of the closest table entry,
     dift=abs(x-xa(i))
     if(dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)  !and initialize the tableau of c's and d's.
     d(i)=ya(i)
  enddo

  y=ya(ns)       !This is the initial approximation to y.
  ns=ns-1

  do m=1,n-1  !For each column of the tableau,
     do i=1,n-m   !we loop over the current c's and d's and update them.
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) stop 'failure in polint'  !two input xa's are identical.
           den=w/den
           d(i)=hp*den    !Here the c's and d's are updated.
           c(i)=ho*den
     enddo
     if(2*ns.lt.n-m)then  !After each column in the tableau is completed, we decide
        dy=c(ns+1)        !which correction, c or d, we want to add to our accumulating
     else                 !value of y, i.e., which path to take through
        dy=d(ns)          !the tableau-forking up or down. We do this in such a
        ns=ns-1           !way as to take the most "straight line" route through the
     endif                !tableau to its apex, updating ns accordingly to keep track
     y=y+dy               !of where we are. This route keeps the partial approximations
  enddo                   !centered (insofar as possible) on the target x. T he
                          !last dy added is thus the error indication.
end subroutine polint
!-----------------------------------------------------------------------
! EOP

