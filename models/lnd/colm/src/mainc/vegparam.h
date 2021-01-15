
  integer, parameter :: numlandc= 21

  real(r8), dimension(numlandc):: &!
       z0m_s, &!roughness
    displa_s, &!zero-plane-distance
    sqrtdi_s, &!inverse sqrt of leaf dimension [m**-0.5]
      chil_s, &!leaf angle distribution factor
  ref_lvis_s, &!leaf reflectance of vis.
  ref_lnir_s, &!leaf reflectance of nir.
  ref_svis_s, &!stem reflectance of vis.
  ref_snir_s, &!stem reflectance of nir.
 tran_lvis_s, &!leaf transmittance of vis.
 tran_lnir_s, &!leaf transmittance of nir.
 tran_svis_s, &!stem transmittance of vis.
 tran_snir_s, &!stem transmittance of nir.
     vmax0_s, &!maximum carboxylation rate at 25 C at canopy top
 vmax0_s_iapdgvm, &!maximum carboxylation rate at 25 C at canopy top
    effcon_s, &!quantum efficiency
     gradm_s, &!conductance-photosynthesis slope parameter
    binter_s, &!conductance-photosynthesis intercept
    respcp_s, &!respiration fraction
      shti_s, &!slope of high temperature inhibition function (s1)
      slti_s, &!slope of low temperature inhibition function (s3)
      trda_s, &!temperature coefficient in gs-a model (s5)
      trdm_s, &!temperature coefficient in gs-a model (s6)
      trop_s, &!temperature coefficient in gs-a model (273.16+25)
      hhti_s, &!1/2 point of high temperature inhibition function (s2)
      hlti_s, &!1/2 point of low temperature inhibition function (s4)
     extkn_s, &!coefficient of leaf nitrogen allocation
       d50_s, &!depth at 50% roots (cm)
       d95_s, &!depth at 95% roots (cm) 
      beta_s, &!coefficient of root profile
   roota_par, &!added by zhq.temporary use in jan20,09
   rootb_par, &!
       z0m_r, &!roughness
    displa_r   !zero-plane-distance
   
  integer  L

!--------------------------------------------------------------------------
!land cover type:
! 1 needleleaf evergreen temperate tree
! 2 needleleaf evergreen boreal tree
! 3 needleleaf deciduous boreal tree
! 4 broadleaf evergreen tropical tree
! 5 broadleaf evergreen temperate tree
! 6 broadleaf deciduous tropical tree
! 7 broadleaf deciduous temperate tree
! 8 broadleaf deciduous boreal tree
! 9 broadleaf evergreen shrub
!10 broadleaf deciduous temperate shrub
!11 broadleaf deciduous boreal shrub
!12 c3 arctic grass
!13 c3 non-arctic grass
!14 c4 grass
!15 corn
!16 wheat
!17 bare soil
!18 lake
!19 wetland
!20 glacier
!21 urban
!22 ocean
!-------------------------------------------------------------------------------------

!Plant functional type optical properties were taken from Dorman and Sellers (1989).
!*Optical properties for intercepted snow will be taken from Sellers et al.(1986)

  z0m_s    =(/0.700,  0.700,  0.700,  0.700,  3.500,  3.500,  2.000,&
              2.000,  2.000,  0.050,  0.050,  0.050,  0.100,  0.100,&
              0.100,  0.100,  0.100,  0.100,  0.100,  0.100,  0.100/)
  displa_s =(/11.333, 11.333, 11.333, 23.333, 23.333, 13.333, 13.333,&
              13.333,  0.333,  0.333,  0.333,  0.667,  0.667,  0.667,&
               0.667,  0.667,  0.667,  0.667,  0.667,  0.667,  0.667/)
  sqrtdi_s(:)=5.0
  chil_s   =(/0.010,  0.010,  0.010,  0.100,  0.100,  0.010,  0.250,&
              0.250,  0.010,  0.250,  0.250, -0.300, -0.300, -0.300,&
             -0.300, -0.300,  0.000, -0.500,  0.650,  0.650, -0.500/)
ref_lvis_s =(/0.070,  0.070,  0.070,  0.100,  0.100,  0.100,  0.100,&!rholvis
              0.100,  0.070,  0.100,  0.100,  0.110,  0.110,  0.110,&
              0.110,  0.110,  0.000,  0.110,  0.110,  0.110,  0.110/)
ref_lnir_s =(/0.350,  0.350,  0.350,  0.450,  0.450,  0.450,  0.450,&!rholnir
              0.450,  0.350,  0.450,  0.450,  0.350,  0.350,  0.350,&
              0.350,  0.350,  0.000,  0.350,  0.350,  0.350,  0.350/)
ref_svis_s =(/0.160,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160,&!rhosvis
              0.160,  0.160,  0.160,  0.160,  0.310,  0.310,  0.310,&
              0.310,  0.310,  0.000,  0.310,  0.310,  0.310,  0.310/)
ref_snir_s =(/0.390,  0.390,  0.390,  0.390,  0.390,  0.390,  0.390,&!rhosnir
              0.390,  0.390,  0.390,  0.390,  0.530,  0.530,  0.530,&
              0.530,  0.530,  0.000,  0.530,  0.530,  0.530,  0.530/)
tran_lvis_s=(/0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&!taulvis
              0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
              0.050,  0.050,  0.000,  0.050,  0.050,  0.050,  0.050/)
tran_lnir_s=(/0.100,  0.100,  0.100,  0.250,  0.250,  0.250,  0.250,&!taulnir
              0.250,  0.100,  0.250,  0.250,  0.340,  0.340,  0.340,&
              0.340,  0.340,  0.000,  0.340,  0.340,  0.340,  0.340/)
tran_svis_s=(/0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,&!tausvis
              0.001,  0.001,  0.001,  0.001,  0.120,  0.120,  0.120,&
              0.120,  0.120,  0.000,  0.120,  0.120,  0.120,  0.120/)
tran_snir_s=(/0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,&!tausnir
              0.001,  0.001,  0.001,  0.001,  0.250,  0.250,  0.250,&
              0.250,  0.250,  0.000,  0.250,  0.250,  0.250,  0.250/)

! NCAR CLM3 PFT physiology parameter, add by zhq.23/03/2009
  vmax0_s  =(/ 51.0,   43.0,   43.0,   50.0,   45.0,   40.0,   51.0,&
               51.0,   17.0,   17.0,   33.0,   43.0,   43.0,   24.0,&
               50.0,   50.0,    0.0,  100.0,  100.0,  100.0,  100.0/)
  vmax0_s_iapdgvm  =(/ 30.6,    37.41,   17.63,   60.0,  33.81,  19.6,   32.64,&
               44.37,   11.90,  11.9,   22.44,  15.05,  15.48,   33.0,&
               50.0,   50.0,    0.0,  100.0,  100.0,  100.0,  100.0/)
  z0m_r    =(/0.055,  0.055,  0.055,  0.075,  0.075,  0.055,  0.055,&
              0.055,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,&
              0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120/)
  displa_r =(/ 0.67,   0.67,   0.67,   0.67,   0.67,   0.67,   0.67,&
               0.67,   0.68,   0.68,   0.68,   0.68,   0.68,   0.68,&
               0.68,   0.68,   0.68,   0.68,   0.68,   0.68,   0.68/)

  effcon_s(1:18)=0.08;   effcon_s(19:21)=0.05
  gradm_s (1:18)=9.0;    gradm_s (19:21)=4.0
  binter_s(1:11)=0.01 
  binter_s(12:16)=0.01;  binter_s(19:21)=0.04
  respcp_s(1:18)=0.015;  respcp_s(19:21)=0.025

  shti_s(:)=0.3
  slti_s(:)=0.2
  trda_s(:)=1.3
  trdm_s(:)=328.0
  trop_s(:)=298.0

  hhti_s=(/303.0, 303.0, 303.0, 313.0, 313.0, 311.0, 311.0,&
           311.0, 313.0, 313.0, 313.0, 308.0, 308.0, 308.0,&
           308.0, 308.0, 313.0, 308.0, 308.0, 308.0, 308.0/)
  hlti_s=(/278.0, 278.0, 278.0, 288.0, 288.0, 283.0, 283.0,&
           283.0, 283.0, 283.0, 283.0, 281.0, 281.0, 281.0,&
           281.0, 281.0, 288.0, 281.0, 281.0, 281.0, 281.0/)
  extkn_s(:)=0.5

  d50_s =(/15.0,  15.0,  16.0,  15.0,  15.0,  16.0,  16.0,&
           16.0,  47.0,  47.0,  47.0,   9.3,   9.3,   9.3,&
           23.5,  23.5,   9.0,   1.0,   9.3,   1.0,   1.0/)
  d95_s =(/ 91.0,  91.0,  95.0,  91.0,  91.0,  95.0,  95.0,&
            95.0, 302.0, 302.0, 302.0,  49.0,  49.0,  49.0,&
           104.8, 104.8,  29.0,   1.0,   9.3,   1.0,   1.0/)
  beta_s=(/-1.632, -1.632, -1.681, -1.632, -1.632, -1.681,&
           -1.681, -1.681, -3.245, -3.245, -3.245, -1.359,&
           -1.359, -1.359, -1.710, -1.710, -2.621, -1.000,&
           -1.359, -1.000, -1.000/)
  ! temporary use to compare CoLM & NCAR CLM3.5's root distribution scheme - added by zhq. 30/12/2008
  roota_par =(/7.0,7.0,7.0,7.0,7.0,6.0,6.0,6.0,7.0,7.0,7.0,&
               11.0,11.0,11.0,6.0,6.0,0.0,0.0,0.0,0.0,0.0/)
  rootb_par =(/2.0,2.0,2.0,1.0,1.0,2.0,2.0,2.0,1.5,1.5,1.5,&
               2.0,2.0,2.0,3.0,3.0,0.0,0.0,0.0,0.0,0.0/)
