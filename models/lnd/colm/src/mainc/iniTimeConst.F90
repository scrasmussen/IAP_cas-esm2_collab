#include <define.h>

 subroutine iniTimeConst (nl_soil, n_pft , isc   , ivt   , sand  , clay  , soc   , rockdep &
                         ,itypwat, zsoi  , dzsoi , albsat, albdry, csol  , porsl   &
                         ,phi0   , bsw   , dkmg  , dksatu, dkdry , hksati, lakedepth, dz_lake &
                         ,z0m    , displa, sqrtdi, effcon, vmax25, slti    &
                         ,hlti   , shti  , hhti  , trda  , trdm  , trop    &
                         ,gradm  , binter, extkn , chil  , refl  , refs  , tranl , trans &
                         ,rootfr , zlnd  , zsno  , csoilc, dewmx , wtfact  &
                         ,capr   , cnfac , ssi   , wimp  , pondmx, smpmax  &
                         ,smpmin , trsmx0, tcrit )  
!===========================================================================
! Initialize time invariant model variables
! Original author: Yongjiu Dai, 09/15/1999; 08/30/2002
!===========================================================================

  use precision
  use paramodel, only: nsoilcolor, nl_lake
  implicit none

!--------------------------- Input 
  integer, INTENT(in) ::            &!
        nl_soil                   , &!number of model soil layers
        n_pft                     , &!number of pfts in a column
        isc                       , &!index for soil color type [-]
        ivt(1:n_pft)                 !land cover type

  real(r8), INTENT(in) ::  &!
        rockdep                      !depth to bed rock
  real(r8), INTENT(inout) ::        &!
        sand(1:nl_soil)           , &!percent of sand
        clay(1:nl_soil)           , &!percent of clay
         soc(1:nl_soil)              !density of organic carbon [kg/m^3]

!--------------------------- Output 
  real(r8), INTENT(out) :: &!--- Soil layer thickness, depths
        zsoi(1:nl_soil)          , &!soil layer depth [m]
       dzsoi(1:nl_soil)             !soil node thickness [m]

  integer, INTENT(in) ::  &!
        itypwat                     !land water type (0=soil, 1=urban, 2=wetland, 
                                    !3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(out) ::         &!--- Soil parameters
        albsat(2)                , &!wet soil albedo [vis/nir] [-]
        albdry(2)                , &!dry soil albedo [vis/nir] [-]
        csol  (1:nl_soil)        , &!heat capacity of soil solids [J/(m3 K)]
        porsl (1:nl_soil)        , &!fraction of soil that is voids [-]
        phi0  (1:nl_soil)        , &!minimum soil suction [mm]
        bsw   (1:nl_soil)        , &!clapp and hornbereger "b" parameter [-]
        dkmg  (1:nl_soil)        , &!thermal conductivity of soil minerals [W/m-K]
        dksatu(1:nl_soil)        , &!thermal conductivity of saturated soil [W/m-K]
        dkdry (1:nl_soil)        , &!thermal conductivity for dry soil  [W/(m-K)]
        hksati(1:nl_soil)        , &!hydraulic conductivity at saturation [mm h2o/s] 
        lakedepth                , &!lake depth [m] added by Nan Wei
        dz_lake(1:nl_lake)          !lake thickness [m] added by Nan Wei

  real(r8), INTENT(out) :: &!--- Vegetation static parameters
        z0m(1:n_pft)              , &!aerodynamic roughness length [m]
        displa(1:n_pft)           , &!displacement height [m]
        sqrtdi(1:n_pft)           , &!inverse sqrt of leaf dimension [m**-0.5]
        effcon(1:n_pft)           , &!quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25(1:n_pft)           , &!maximum carboxylation rate at 25 C at canopy top (mol CO2/m2s)

        shti(1:n_pft)             , &!slope of high temperature inhibition function     (s1)
        hhti(1:n_pft)             , &!1/2 point of high temperature inhibition function (s2)
        slti(1:n_pft)             , &!slope of low temperature inhibition function      (s3)
        hlti(1:n_pft)             , &!1/2 point of low temperature inhibition function  (s4)
        trda(1:n_pft)             , &!temperature coefficient in gs-a model             (s5)
        trdm(1:n_pft)             , &!temperature coefficient in gs-a model             (s6)
        trop(1:n_pft)             , &!temperature coefficient in gs-a model         (273+25)
        gradm(1:n_pft)            , &!conductance-photosynthesis slope parameter
        binter(1:n_pft)           , &!conductance-photosynthesis intercep
        extkn(1:n_pft)            , &!coefficient of leaf nitrogen allocation
        chil(1:n_pft)             , &!leaf angle distribution factor
        refl(2,2,1:n_pft)         , &!leaf reflectance (iw=iband, il=live and dead)
        refs(2,2,1:n_pft)         , &!stem reflectance (iw=iband, il=live and dead)
        tranl(2,2,1:n_pft)        , &!leaf transmittance (iw=iband, il=live and dead)
        trans(2,2,1:n_pft)        , &!stem transmittance (iw=iband, il=live and dead)
        rootfr(1:nl_soil,1:n_pft)    !fraction of roots in each soil layer 

  real(r8), INTENT(out) :: &!--- Initialize TUNABLE constants
        zlnd             , &!Roughness length for soil [m]
        zsno             , &!Roughness length for snow [m]
        csoilc           , &!Drag coefficient for soil under canopy [-]
        dewmx            , &!maximum dew
        wtfact           , &!Fraction of model area with high water table
        capr             , &!Tuning factor to turn first layer T into surface T
        cnfac            , &!Crank Nicholson factor between 0 and 1
        ssi              , &!Irreducible water saturation of snow
        wimp             , &!Water impremeable if porosity less than wimp
        pondmx           , &!Ponding depth (mm)
        smpmax           , &!Wilting point potential in mm
        smpmin           , &!Restriction for min of soil poten. (mm)
        trsmx0           , &!Max transpiration for moist soil+100% veg. [mm/s]
        tcrit               !critical temp. to determine rain or snow

!--------------------------- Local variables
  integer i, j              !indices
  integer c, p              !indices of column and pft
  integer idlak             !index (1=deep lake, 0=shallow)
  real(r8) bd(1:nl_soil)    !bulk density of dry soil material [kg/m^3]
  real(r8) dkm(1:nl_soil)   !
! real(r8) dzlak(1:nl_soil) !
! real(r8) zlak(1:nl_soil)  !
  real(r8) zsoih(0:nl_soil) !interface level below a zsoi level [m]
  real(r8) fsoc(1:nl_soil)  !fraction of soil organic carbon

!--------------------------- Data block
! Soil albedo for different colored soils (saturated soil, visible/nir beam) [-]

  real(r8) :: solour_sat(2,20) ![vis/nir,soilcolor]
  real(r8) :: solour_dry(2,20) ![vis/nir,soilcolor]

#include "vegparam.h"

     solour_sat(1,1:20) = (/0.25_r8,0.23_r8,0.21_r8,0.20_r8,0.19_r8,0.18_r8,0.17_r8,0.16_r8,&
                            0.15_r8,0.14_r8,0.13_r8,0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8,0.04_r8/)
     solour_sat(2,1:20) = (/0.50_r8,0.46_r8,0.42_r8,0.40_r8,0.38_r8,0.36_r8,0.34_r8,0.32_r8,&
                            0.30_r8,0.28_r8,0.26_r8,0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
     solour_dry(1,1:20) = (/0.36_r8,0.34_r8,0.32_r8,0.31_r8,0.30_r8,0.29_r8,0.28_r8,0.27_r8,&
                            0.26_r8,0.25_r8,0.24_r8,0.23_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
     solour_dry(2,1:20) = (/0.61_r8,0.57_r8,0.53_r8,0.51_r8,0.49_r8,0.48_r8,0.45_r8,0.43_r8,&
                            0.41_r8,0.39_r8,0.37_r8,0.35_r8,0.33_r8,0.31_r8,0.29_r8,0.27_r8,0.25_r8,0.23_r8,0.21_r8,0.16_r8/)

!-----------------------------------------------------------------------
! land water type for USGS classification
!-----------------------------------------------------------------------
!      i=nint(ivt)
!      itypwat=0                     ! soil
!#if(defined USGS)
!      if(i==1)           itypwat=1  ! urban and built-up
!      if(i==17.or.i==18) itypwat=2  ! wetland
!      if(i==24)          itypwat=3  ! land ice
!      if(i==16)          itypwat=4  ! river or deep lake
!      if(i==25)          itypwat=99 ! ocean
!#elif (defined IGBP)
!      if(i==13)          itypwat=1  ! urban and built-up
!      if(i==11)          itypwat=2  ! wetland
!      if(i==15)          itypwat=3  ! land ice
!      if(i==17)          itypwat=4  ! river or deep lake
!      if(i==18)          itypwat=99 ! ocean
!#elif (defined SIB2)
!      if(i==11)          itypwat=3  ! land ice
!      if(i==10)          itypwat=4  ! river or deep lake
!      if(i==12)          itypwat=99 ! ocean
!#elif (defined BATS)
!      if(i==13)          itypwat=2  ! wetland
!      if(i==12)          itypwat=3  ! land ice
!      if(i==14)          itypwat=4  ! river or deep lake
!      if(i==15)          itypwat=99 ! ocean
!#elif (defined OGE)
!      if(i==1)           itypwat=1  ! urban and built-up
!      if(i==13 .or. i==44 .or. i==45) itypwat=2  ! wetland
!      if(i==12)          itypwat=3  ! land ice
!      if(i==14)          itypwat=4  ! river or deep lake
!      if(i==15)          itypwat=99 ! ocean
!#elif (defined DGVM)
!      if(i==21)          itypwat=1  ! urban and built-up
!      if(i==19)          itypwat=2  ! wetland
!      if(i==20)          itypwat=3  ! land ice
!      if(i==18)          itypwat=4  ! river or deep lake
!      if(i==22)          itypwat=99 ! ocean
!#endif

!-----------------------------------------------------------------------
! soil layer thickness, depths (m)
!-----------------------------------------------------------------------

      if(itypwat<=5)then
         do j = 1, nl_soil
           zsoi(j) = 0.025*(exp(0.5*(j-0.5))-1.)  !node depths
         end do

         dzsoi(1)  = 0.5*(zsoi(1)+zsoi(2))        !=zsoih(1)
         dzsoi(nl_soil)= zsoi(nl_soil)-zsoi(nl_soil-1)
         do j = 2,nl_soil-1
           dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1))    !thickness b/n two interfaces
         end do

         zsoih(0)   = 0.
         zsoih(nl_soil) = zsoi(nl_soil) + 0.5*dzsoi(nl_soil)
         do j = 1, nl_soil-1
            zsoih(j)= 0.5*(zsoi(j)+zsoi(j+1))     !interface depths
         enddo

      ! ------ Lake ------ !
         if(itypwat==4 .or. itypwat==5)then
            idlak = 1                                !assumed all lakes are deep lake
            lakedepth = 0.0
            if(idlak == 1) then                                                 
               dz_lake = (/0.1, 1., 2., 3., 4., 5., 7., 7., 10.45, 10.45/)  ! m
               do j = 1, nl_lake
                  lakedepth = lakedepth + dz_lake(j)
               end do
!           else
!              dzlak = (/.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)   
!              zlak = (/ 0.125,  0.5,  1.125,  2.,  3.,  4.,  5.,  6.,  7.,  8./)
            end if
         end if 
!           zsoih(0) = 0.
!           do j = 1, nl_lake
!              dzsoi(j) = dzlak(j)
!               zsoi(j) =  zlak(j)
!              zsoih(j) = zsoih(j-1)+dzsoi(j)
!           enddo

      ! ------ Ocean ------ !
      else
            dzsoi(:) = -999.9
             zsoi(:) = -999.9
            zsoih(:) = -999.9
      endif

!-----------------------------------------------------------------------
! soil thermal and hydraulic properties
!-----------------------------------------------------------------------

! saturated or dry soil albedo for visible beam

      albsat = 0.
      albdry = 0.

      if(isc>=1 .and. isc<=nsoilcolor) then
         albsat(:) = solour_sat(:,isc)
         albdry(:) = solour_dry(:,isc)
      end if

! soil thermal and hydraulic properties

      fsoc = soc/130.

      do j = 1, nl_soil
         if(itypwat<2)then  ! not wetland, glacier and lake
            if(zsoi(j)<=rockdep)then  !NON ROCK
             ! for all organic case or data missing, assigned to "loam"
               if(sand(j)*clay(j).lt.0.01)then
                  sand(j) = 43.
                  clay(j) = 18.
               endif
   
                porsl(j) = 0.489 - 0.00126*sand(j)
                porsl(j) = fsoc(j)*0.9 + (1.-fsoc(j))*porsl(j)
                 phi0(j) = 10. * ( 10.**(1.88-0.0131*sand(j)) )
                 phi0(j) = fsoc(j)*10.3 + (1.-fsoc(j))*phi0(j)
                  bsw(j) = 2.91 + 0.159*clay(j)
                  bsw(j) = fsoc(j)*2.7 + (1.-fsoc(j))*bsw(j)
               hksati(j) = 0.0070556 * ( 10.**(-0.884+0.0153*sand(j)) ) ! mm/s
               hksati(j) = fsoc(j)*1.0e-4 + (1-fsoc(j))*hksati(j)
      
                   bd(j) = (1.- porsl(j))*2.7e3
                 csol(j) = (2.128*sand(j)+2.385*clay(j))/(sand(j)+clay(j))*1.e6  ! J/(m3 K)
                 csol(j) = fsoc(j)*2.5e6 + (1.-fsoc(j))*csol(j)
                 dkm(j)  = (8.80*sand(j)+2.92*clay(j)) / (sand(j)+clay(j))       ! W/(m K)
                 dkm(j)  = fsoc(j)*0.25 + (1.-fsoc(j))*dkm(j)
                 dkmg(j) = dkm(j)**(1.-porsl(j))
               dksatu(j) = dkmg(j)*0.57**porsl(j)
                dkdry(j) = (.135*bd(j) + 64.7) / (2.7e3 - 0.947*bd(j))
                dkdry(j) = fsoc(j)*0.05 + (1.-fsoc(j))*dkdry(j)
   
               if(zsoih(j)>rockdep)then
                  porsl(j) = porsl(j)*(rockdep-zsoih(j-1))/dzsoi(j)
               endif
            else                      !BEDROCK
                porsl(j) = 0.
                 phi0(j) = 1.5e-5
                  bsw(j) = 0
               hksati(j) = 0.
                   bd(j) = 2.7e3
                 csol(j) = 2700.*750. !J/(m3 K)
                 dkmg(j) = 1.0        ! not used
               dksatu(j) = 2.9        ! W/(m K)
                dkdry(j) = 2.9        ! W/(m K) [J.R. Garratt 1992, pp291]
            endif
         else                ! wetland, glacier, lake and ocean
                porsl(j) = 1.0
                 phi0(j) = 0.0
                  bsw(j) = 0.0
                 csol(j) = 4.186e06
      
                 dkmg(j) = 1.0
               dksatu(j) = 0.6
                dkdry(j) = 0.6
               hksati(j) = 0.0
         endif
      end do


!-----------------------------------------------------------------------
! vegetation static parameters
! the values for glacier and lake are assigned arbitrarily (not used)
!-----------------------------------------------------------------------

  do p = 1, n_pft

! For Lake, Glacier, Ocean
    z0m(p)       = 0.0
    displa(p)    = 0.0
    sqrtdi(p)    = 0.0
    vmax25(p)    = 0.0
    effcon(p)    = 0.0
    slti(p)      = 0.0
    hlti(p)      = 0.0
    shti(p)      = 0.0
    hhti(p)      = 0.0
    trda(p)      = 0.0
    trdm(p)      = 0.0
    trop(p)      = 0.0
    gradm(p)     = 0.0
    binter(p)    = 0.0
    extkn(p)     = 0.0
    chil(p)      = 0.0
    refl(:,:,p)  = 0.0
    refs(:,:,p)  = 0.0
    tranl(:,:,p) = 0.0
    trans(:,:,p) = 0.0
    rootfr(:,p)  = 0.0

    if(itypwat==3)then  ! glacier
       z0m(p) = 0.01
       displa(p) = 0.0
    endif

    if(itypwat<3)then   ! land grids
          i = ivt(p)
          z0m(p) = z0m_s   (i)  
       displa(p) = displa_s(i)
       sqrtdi(p) = sqrtdi_s(i)  
#ifdef IAPDGVM
       vmax25(p) = vmax0_s_iapdgvm(i)*1.e-6       
#else
       vmax25(p) = vmax0_s(i)*1.e-6       
#endif
       effcon(p) = effcon_s(i)       
       slti(p)   =   slti_s(i)       
       hlti(p)   =   hlti_s(i)       
       shti(p)   =   shti_s(i)       
       hhti(p)   =   hhti_s(i)       
       trda(p)   =   trda_s(i)       
       trdm(p)   =   trdm_s(i)       
       trop(p)   =   trop_s(i)       
       gradm(p)  =  gradm_s(i)       
       binter(p) = binter_s(i)       
       extkn(p)  =  extkn_s(i)       
  
           chil(p)  = chil_s(i)             
       refl(1,1,p)  = ref_lvis_s(i)      !live 
       refl(2,1,p)  = ref_lnir_s(i)      !live 
       refl(1,2,p)  = -9999.             !dead
       refl(2,2,p)  = -9999.             !dead

       refs(1,1,p)  = ref_svis_s(i)      !live 
       refs(2,1,p)  = ref_snir_s(i)      !live 
       refs(1,2,p)  = -9999.             !dead
       refs(2,2,p)  = -9999.             !dead

       tranl(1,1,p) = tran_lvis_s(i)     !live 
       tranl(2,1,p) = tran_lnir_s(i)     !live 
       tranl(1,2,p) = -9999.             !dead 
       tranl(2,2,p) = -9999.             !dead

       trans(1,1,p) = tran_svis_s(i)     !live 
       trans(2,1,p) = tran_snir_s(i)     !live 
       trans(1,2,p) = -9999.             !dead 
       trans(2,2,p) = -9999.             !dead

!#if(defined USGS) 
      ! The definition of global root distribution is based on
      ! Schenk and Jackson, 2002: The Global Biogeography of Roots.
      ! Ecological Monagraph 72(3): 311-328.
       rootfr(1,p)=1./(1.+(zsoih(1)*100./d50_s(i))**beta_s(i)) 
       rootfr(nl_soil,p)=1.-1./(1.+(zsoih(nl_soil-1)*100./d50_s(i))**beta_s(i)) 

       do j=2,nl_soil-1
          rootfr(j,p)=1./(1.+(zsoih(j)*100./d50_s(i))**beta_s(i)) &
                     -1./(1.+(zsoih(j-1)*100./d50_s(i))**beta_s(i))
       enddo
!#elif(defined DGVM)
      ! Initialize root fraction (computing from surface, d is depth in meter):
      ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
      ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
      ! beta & d_obs given in Zeng et al. (1998).
!
!      do j = 1, nl_soil-1
!      rootfr(j) = .5*( exp(-roota_par(i) * zsoih(j-1))  &
!                         + exp(-rootb_par(i) * zsoih(j-1))  &
!                         - exp(-roota_par(i) * zsoih(j))  &
!                         - exp(-rootb_par(i) * zsoih(j)) )
!      end do
!      rootfr(nl_soil) = .5*( exp(-roota_par(i) * zsoih(j-1))  &
!                         + exp(-rootb_par(i) * zsoih(j-1)) )
!      
!#endif
    endif

#if(defined DGVM)
    if(itypwat==0)then
         i = ivt(p)
         z0m(p) = z0m_r   (i)  
      displa(p) = displa_r(i)
    endif
#endif

  end do
!-----------------------------------------------------------------------
! Initialize TUNABLE constants
!-----------------------------------------------------------------------

      zlnd   = 0.01    !Roughness length for soil [m]
      zsno   = 0.0024  !Roughness length for snow [m]
      csoilc = 0.004   !Drag coefficient for soil under canopy [-]
      dewmx  = 0.1     !maximum dew
      wtfact = 0.3     !Fraction of model area with high water table
      capr   = 0.34    !Tuning factor to turn first layer T into surface T
      cnfac  = 0.5     !Crank Nicholson factor between 0 and 1
      ssi    = 0.033   !Irreducible water saturation of snow
      wimp   = 0.05    !Water impremeable if porosity less than wimp
      pondmx = 10.0    !Ponding depth (mm)
      smpmax = -1.5e5  !Wilting point potential in mm
      smpmin = -1.e8   !Restriction for min of soil poten. (mm)
      trsmx0 = 2.e-4   !Max transpiration for moist soil+100% veg. [mm/s]
      tcrit  = 0.      !critical temp. to determine rain or snow

 end subroutine iniTimeConst
