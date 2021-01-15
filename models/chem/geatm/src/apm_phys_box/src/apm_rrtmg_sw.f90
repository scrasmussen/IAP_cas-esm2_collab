
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: rrtmg_sw_GCAPM
!
! !DESCRIPTION: Module rrtmg\_sw\_GCAPM contains variables and routines for 
!  computing RF.
!-- Modified from AER RRTMG_SW for GEOS-Chem-APM 
!     August 2012: F. Yu, UAlbany
!\\
!\\
! !INTERFACE:
!
      MODULE rrtmg_sw_GCAPM
!
! !USES:
!
!      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!  
      PUBLIC  :: rrtmg_sw

!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw.1col.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.8 $
!     created:   $Date: 2009/05/22 22:22:21 $
!

      CONTAINS

      subroutine rrtmg_sw(II,JJ,IFCS,icld,nlayers,juldat, zenith, &
              pdp,pavel,tavel,pz,tz,tbound,coldry, wkl, &
              cldfrac, ciwp, clwp, rei, rel,&
              SALB,EXT,OMGA,G, &
              CST,FST,CSB,FSB,CST0,FST0,CSB0,FSB0, &
              TEXT, TOMGA, TG, TCST,TCSB,TFST,TFSB)

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  miacono@aer.com                             *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people: Steven J. Taubman, Patrick D. Brown,             *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.

! This routine
!    a) calls RRTMG_SW_INI to initialize data and to perform
!       the g-point interval reduction from 224 to 112
!    b) calls READPROF to read in the atmospheric profile;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPROP to set cloud optical depth based on input
!       cloud properties, or CLDPRMC to set cloud optical depth
!       for McICA
!    d) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVRT to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer
!       with or without McICA, the Monte-Carlo Independent Column
!       Approximation to represent sub-grid scale cloud variability
!    f) writes out the calculated fluxes and cooling rates
!
! Two modes of operation are possible:
!     The mode is chosen by setting flag imca below.  
!
!    1) Standard, single forward model calculation (imca = 0); this is 
!       valid only for clear sky or fully overcast clouds
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       For single column calculations, this method also requires setting flag
!       nmca below to the sample size of the Monte Carlo calculation; 
!       (nmca = 200 is recommended). This is method is valid for clear sky
!       or partial cloud conditions
!
! Two random number generators are available for use when imca = 1
!     This is chosen by setting flag irng below.
!
!    1) KISSVEC (irng = 0)
!    2) Mersenne Twister (irng = 1); the default setting
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflag, iceflag and liqflag; see text file rrtmg_sw_instructions
!     and subroutine rrtmg_sw_cldprop.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (inflgsw = 0)
!    2) Input cloud fraction and cloud physical properties (inflgsw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflgsw and liqflgsw. Ice particle size provided
!       must be appropriately defined for the ice parameterization selected. 
!
! Two methods of aerosol property input are possible:
!     Aerosol properties can be input in one of two ways (controlled by input 
!     flag iaer, see text file rrtmg_sw_instructions for further details):
!
!    1) Input aerosol optical depth, single scattering albedo and asymmetry
!       parameter directly by layer and spectral band (iaer=10)
!    2) Input aerosol optical depth and 0.55 micron directly by layer and use
!       one or more of six ECMWF aerosol types (iaer=6)
!
!
! ------- Modifications -------
!
! This version of RRTMG_SW has been modified from RRTM_SW to use a reduced
! set of g-point intervals and a two-stream model for application to GCMs. 
!
!-- Original version (derived from RRTM_SW)
!     2002: AER. Inc.
!-- Conversion to F90 formatting; addition of 2-stream radiative transfer
!     Feb 2003: J.-J. Morcrette, ECMWF
!-- Additional modifications for GCM application
!     Aug 2003: M. J. Iacono, AER Inc.
!-- Total number of g-points reduced from 224 to 112.  Original
!   set of 224 can be restored by exchanging code in module parrrsw.f90 
!   and in file rrtmg_sw_init.f90.
!     Apr 2004: M. J. Iacono, AER, Inc.
!-- Modifications to include output for direct and diffuse 
!   downward fluxes.  There are output as "true" fluxes without
!   any delta scaling applied.  Code can be commented to exclude
!   this calculation in source file rrtmg_sw_spcvrt.f90.
!     Jan 2005: E. J. Mlawer, M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability.
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modified to output direct and diffuse fluxes either with or without
!   delta scaling based on setting of idelm flag. 
!     Dec 2008: M. J. Iacono, AER, Inc.

!-- Modified for GEOS-Chem-APM
!     August 2012: F. Yu, UAlbany

! --------- Modules ---------

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrsw, only : mxlay, nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_vsn
      use mcica_subcol_gen_sw, only: mcica_subcol_sw
      use rrtmg_sw_cldprop, only: cldprop_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvrt, only: spcvrt_sw
      use rrtmg_sw_spcvmc, only: spcvmc_sw

      implicit none

! ------- Declarations

! ----- Local -----

! Control

      integer(kind=im) :: II,JJ
      integer(kind=im) :: nlayers             ! total number of layers
      integer(kind=im) :: juldat             ! total number of layers
      integer(kind=im) :: istart              ! beginning band of calculation
      integer(kind=im) :: iend                ! ending band of calculation
      integer(kind=im) :: icld                ! clear/cloud and cloud overlap flag
      integer(kind=im) :: icpr                ! cldprop/cldprmc use flag
      integer(kind=im) :: iflag               ! control flag
      integer(kind=im) :: iout                ! output option flag
      integer(kind=im) :: iaer                ! aerosol option flag
      integer(kind=im) :: idelm               ! delta-m scaling flag
                                              ! [0 = direct and diffuse fluxes are unscaled]
                                              ! [1 = direct and diffuse fluxes are scaled]
      integer(kind=im) :: isccos              ! instrumental cosine response flag
      integer(kind=im) :: i                   ! layer loop index                      ! jk
      integer(kind=im) :: ib,ib1,ib2                  ! band loop index                       ! jsw
      integer(kind=im) :: ia, ig              ! indices
      integer(kind=im) :: iplon               ! column loop index                     ! jl
      integer(kind=im) :: permuteseed        
      integer(kind=im) :: imca                ! flag for mcica [0=off, 1=on]
!      integer(kind=im) :: nmca                ! number of mcica samples (mcica mode)
      integer(kind=im) :: irng                ! flag for random number generator
                                              ! [0=kissvec, 1=mersenne twister (default)]
      integer(kind=im), parameter :: ncol = 1 ! total number of columns

      integer(kind=im) :: iout1, iout2        ! output control flags
      integer(kind=im) :: indform             ! output control flag
!Yu+
      integer(kind=im) :: IFCS, ITYP
      integer(kind=im) :: isolvar    
      integer(kind=im), parameter :: NTYP = 5 

      character page 

      real(kind=rb) :: zepsec, zepzen         ! epsilon
      real(kind=rb) :: zdpgcp                 ! flux to heating conversion ratio


! Atmosphere
      real(kind=rb) :: pavel(mxlay)           ! layer pressures (mb) 
      real(kind=rb) :: tavel(mxlay)           ! layer temperatures (K)
      real(kind=rb) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=rb) :: tz(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=rb) :: tbound                 ! surface temperature (K)
      real(kind=rb) :: pdp(mxlay)             ! layer pressure thickness (hPa, mb)
      real(kind=rb) :: coldry(mxlay)          ! 
      real(kind=rb) :: wbrodl(mxlay)          !
      real(kind=rb) :: wkl(mxmol,mxlay)       ! molecular amounts (mole/cm-2)

      real(kind=rb) :: cossza, zenith         ! cosine of solar zenith angle 
!      real(kind=rb) :: earth_sun              ! function for Earth/Sun distance factor
      real(kind=rb) :: adjflux(jpband)        ! adjustment for current Earth/Sun distance
      real(kind=rb) :: solvar(jpband)         ! solar constant scaling factor from rrtmg_sw
                                              !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=rb) :: SALB(nbndsw)           ! surface albedo  

      real(kind=rb) :: albdir(nbndsw)         ! surface albedo, direct          ! zalbp
      real(kind=rb) :: albdif(nbndsw)         ! surface albedo, diffuse         ! zalbd

      real(kind=rb) :: EXT(mxlay,nbndsw)
      real(kind=rb) :: OMGA(mxlay,nbndsw)
      real(kind=rb) :: G(mxlay,nbndsw)

      real(kind=rb) :: TEXT(mxlay,nbndsw,NTYP)
      real(kind=rb) :: TOMGA(mxlay,nbndsw,NTYP)
      real(kind=rb) :: TG(mxlay,nbndsw,NTYP)

      real(kind=rb) :: tauaer(mxlay,jpband)   ! aerosol optical depth (iaer=10 only)
                                              ! (non-delta scaled)      
      real(kind=rb) :: ssaaer(mxlay,jpband)   ! aerosol single scattering albedo (iaer=10 only)
                                              ! (non-delta scaled)      
      real(kind=rb) :: asmaer(mxlay,jpband)   ! aerosol asymmetry parameter (iaer=10 only)
                                              ! (non-delta scaled)      
                                              !   first moment of input phase function
      real(kind=rb) :: ecaer(mxlay,naerec)    ! aerosol optical thickness at 0.55 micron (iaer=6 only)
                                              ! (non-delta scaled)      

! Atmosphere - setcoef
      integer(kind=im) :: laytrop             ! tropopause layer index
      integer(kind=im) :: layswtch            ! tropopause layer index
      integer(kind=im) :: laylow              ! tropopause layer index
      integer(kind=im) :: jp(mxlay)           ! 
      integer(kind=im) :: jt(mxlay)           !
      integer(kind=im) :: jt1(mxlay)          !

      real(kind=rb) :: colh2o(mxlay)          ! column amount (h2o)
      real(kind=rb) :: colco2(mxlay)          ! column amount (co2)
      real(kind=rb) :: colo3(mxlay)           ! column amount (o3)
      real(kind=rb) :: coln2o(mxlay)          ! column amount (n2o)
      real(kind=rb) :: colch4(mxlay)          ! column amount (ch4)
      real(kind=rb) :: colo2(mxlay)           ! column amount (o2)
      real(kind=rb) :: colmol(mxlay)          ! column amount
      real(kind=rb) :: co2mult(mxlay)         ! column amount 

      integer(kind=im) :: indself(mxlay)
      integer(kind=im) :: indfor(mxlay)
      real(kind=rb) :: selffac(mxlay)
      real(kind=rb) :: selffrac(mxlay)
      real(kind=rb) :: forfac(mxlay)
      real(kind=rb) :: forfrac(mxlay)

      real(kind=rb) :: &                      !
                         fac00(mxlay), fac01(mxlay), &
                         fac10(mxlay), fac11(mxlay) 

! Atmosphere/clouds - cldprop
      integer(kind=im) :: ncbands             ! number of cloud spectral bands
      integer(kind=im) :: inflag              ! flag for cloud property method
      integer(kind=im) :: iceflag             ! flag for ice cloud properties
      integer(kind=im) :: liqflag             ! flag for liquid cloud properties

      real(kind=rb) :: cldfrac(mxlay)         ! layer cloud fraction
      real(kind=rb) :: tauc(nbndsw,mxlay)     ! in-cloud optical depth (non-delta scaled)
      real(kind=rb) :: ssac(nbndsw,mxlay)     ! in-cloud single scattering albedo (non-delta scaled)
      real(kind=rb) :: asmc(nbndsw,mxlay)     ! in-cloud asymmetry parameter (non-delta scaled)
      real(kind=rb) :: fsfc(nbndsw,mxlay)     ! in-cloud forward scattering fraction (non-delta scaled)
      real(kind=rb) :: ciwp(mxlay)            ! in-cloud ice water path
      real(kind=rb) :: clwp(mxlay)            ! in-cloud liquid water path
      real(kind=rb) :: rei(mxlay)             ! cloud ice particle effective size (microns)
                                              ! specific definition of rei depends on setting of iceflag:
                                              ! iceflag = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                              !              r_ec must be >= 10.0 microns
                                              ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                              !              r_ec range is limited to 13.0 to 130.0 microns
                                              ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                              !              r_k range is limited to 5.0 to 131.0 microns
                                              ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                              !              dge range is limited to 5.0 to 140.0 microns
                                              !              [dge = 1.0315 * r_ec]
      real(kind=rb) :: rel(mxlay)             ! cloud liquid particle effective radius (microns)

      real(kind=rb) :: taucloud(mxlay,jpband) ! in-cloud optical depth
      real(kind=rb) :: taucldorig(mxlay,jpband)! in-cloud optical depth (non-delta scaled)
      real(kind=rb) :: ssacloud(mxlay,jpband) ! in-cloud single scattering albedo
      real(kind=rb) :: asmcloud(mxlay,jpband) ! in-cloud asymmetry parameter

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=rb) :: cldfmc(ngptsw,mxlay)   ! cloud fraction [mcica]
      real(kind=rb) :: ciwpmc(ngptsw,mxlay)   ! in-cloud ice water path [mcica]
      real(kind=rb) :: clwpmc(ngptsw,mxlay)   ! in-cloud liquid water path [mcica]
      real(kind=rb) :: relqmc(mxlay)          ! liquid particle effective radius (microns)
      real(kind=rb) :: reicmc(mxlay)          ! ice particle effective radius (microns)
      real(kind=rb) :: taucmc(ngptsw,mxlay)   ! in-cloud optical depth [mcica]
      real(kind=rb) :: taormc(ngptsw,mxlay)   ! unscaled in-cloud optical depth [mcica]
      real(kind=rb) :: ssacmc(ngptsw,mxlay)   ! in-cloud single scattering albedo [mcica]
      real(kind=rb) :: asmcmc(ngptsw,mxlay)   ! in-cloud asymmetry parameter [mcica]
      real(kind=rb) :: fsfcmc(ngptsw,mxlay)   ! in-cloud forward scattering fraction [mcica]
! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=rb) :: ztauc(mxlay,nbndsw)    ! cloud optical depth
      real(kind=rb) :: ztaucorig(mxlay,nbndsw)! unscaled cloud optical depth
      real(kind=rb) :: zasyc(mxlay,nbndsw)    ! cloud asymmetry parameter 
                                              !  (first moment of phase function)
      real(kind=rb) :: zomgc(mxlay,nbndsw)    ! cloud single scattering albedo
      real(kind=rb) :: ztaua(mxlay,nbndsw)    ! total aerosol optical depth
      real(kind=rb) :: zasya(mxlay,nbndsw)    ! total aerosol asymmetry parameter 
      real(kind=rb) :: zomga(mxlay,nbndsw)    ! total aerosol single scattering albedo

      real(kind=rb) :: zcldfmc(mxlay,ngptsw)  ! cloud fraction [mcica]
      real(kind=rb) :: ztaucmc(mxlay,ngptsw)  ! cloud optical depth [mcica]
      real(kind=rb) :: ztaormc(mxlay,ngptsw)  ! unscaled cloud optical depth [mcica]
      real(kind=rb) :: zasycmc(mxlay,ngptsw)  ! cloud asymmetry parameter [mcica] 
      real(kind=rb) :: zomgcmc(mxlay,ngptsw)  ! cloud single scattering albedo [mcica]

      real(kind=rb) :: zbbfu(mxlay+1)         ! temporary upward shortwave flux (w/m2)
      real(kind=rb) :: zbbfd(mxlay+1)         ! temporary downward shortwave flux (w/m2)
      real(kind=rb) :: zbbcu(mxlay+1)         ! temporary clear sky upward shortwave flux (w/m2)
      real(kind=rb) :: zbbcd(mxlay+1)         ! temporary clear sky downward shortwave flux (w/m2)
      real(kind=rb) :: zbbfddir(mxlay+1)      ! temporary downward direct shortwave flux (w/m2)
      real(kind=rb) :: zbbcddir(mxlay+1)      ! temporary clear sky downward direct shortwave flux (w/m2)
      real(kind=rb) :: zuvfd(mxlay+1)         ! temporary UV downward shortwave flux (w/m2)
      real(kind=rb) :: zuvcd(mxlay+1)         ! temporary clear sky UV downward shortwave flux (w/m2)
      real(kind=rb) :: zuvfddir(mxlay+1)      ! temporary UV downward direct shortwave flux (w/m2)
      real(kind=rb) :: zuvcddir(mxlay+1)      ! temporary clear sky UV downward direct shortwave flux (w/m2)
      real(kind=rb) :: znifd(mxlay+1)         ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=rb) :: znicd(mxlay+1)         ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real(kind=rb) :: znifddir(mxlay+1)      ! temporary near-IR downward direct shortwave flux (w/m2)
      real(kind=rb) :: znicddir(mxlay+1)      ! temporary clear sky near-IR downward direct shortwave flux (w/m2)

! Parameters
      real(kind=rb), parameter :: cpdair = 1.004e3_rb  ! Specific heat capacity of dry air
                                                       ! at constant pressure at 273 K (J kg-1 K-1)
! Output
      real(kind=rb) :: totuflux(0:mxlay)      ! upward shortwave flux (w/m2)                  ! pfup
      real(kind=rb) :: totdflux(0:mxlay)      ! downward shortwave flux (w/m2)                ! pfdown
      real(kind=rb) :: fnet(0:mxlay)          ! net shortwave flux (w/m2)                     ! pfls
      real(kind=rb) :: htr(0:mxlay)           ! shortwave heating rate (k/day)                ! pheat
      real(kind=rb) :: totuclfl(0:mxlay)      ! clear sky upward shortwave flux (w/m2)        ! pcup 
      real(kind=rb) :: totdclfl(0:mxlay)      ! clear sky downward shortwave flux (w/m2)      ! pcdown 
      real(kind=rb) :: fnetc(0:mxlay)         ! clear sky net shortwave flux (w/m2)           ! pfcs
      real(kind=rb) :: htrc(0:mxlay)          ! clear sky shortwave heating rate (k/day)      ! pheac

      real(kind=rb) :: dirdflux(0:mxlay)      ! direct downward shortwave flux (w/m2)         ! dirdownflux
      real(kind=rb) :: difdflux(0:mxlay)      ! diffuse downward shortwave flux (w/m2)        ! difdownflux
      real(kind=rb) :: dflxuv(0:mxlay)        ! Total sky downward shortwave flux, UV/vis     ! pfdnuv
      real(kind=rb) :: dflxir(0:mxlay)        ! Total sky downward shortwave flux, near-IR    ! pfdnir 
      real(kind=rb) :: dirdnuv(0:mxlay)       ! Direct downward shortwave surface flux, UV/vis
      real(kind=rb) :: difdnuv(0:mxlay)       ! Diffuse downward shortwave surface flux, UV/vis
      real(kind=rb) :: dirdnir(0:mxlay)       ! Direct downward shortwave surface flux, near-IR
      real(kind=rb) :: difdnir(0:mxlay)       ! Diffuse downward shortwave surface flux, near-IR

! Output - inactive
!      real(kind=rb) :: zuvfu(mxlay+1)         ! temporary upward UV shortwave flux (w/m2)
!      real(kind=rb) :: zuvfd(mxlay+1)         ! temporary downward UV shortwave flux (w/m2)
!      real(kind=rb) :: zuvcu(mxlay+1)         ! temporary clear sky upward UV shortwave flux (w/m2)
!      real(kind=rb) :: zuvcd(mxlay+1)         ! temporary clear sky downward UV shortwave flux (w/m2)
!      real(kind=rb) :: zvsfu(mxlay+1)         ! temporary upward visible shortwave flux (w/m2)
!      real(kind=rb) :: zvsfd(mxlay+1)         ! temporary downward visible shortwave flux (w/m2)
!      real(kind=rb) :: zvscu(mxlay+1)         ! temporary clear sky upward visible shortwave flux (w/m2)
!      real(kind=rb) :: zvscd(mxlay+1)         ! temporary clear sky downward visible shortwave flux (w/m2)
!      real(kind=rb) :: znifu(mxlay+1)         ! temporary upward near-IR shortwave flux (w/m2)
!      real(kind=rb) :: znifd(mxlay+1)         ! temporary downward near-IR shortwave flux (w/m2)
!      real(kind=rb) :: znicu(mxlay+1)         ! temporary clear sky upward near-IR shortwave flux (w/m2)
!      real(kind=rb) :: znicd(mxlay+1)         ! temporary clear sky downward near-IR shortwave flux (w/m2)

!Yu+
      real(kind=rb) :: xpavel(ncol,mxlay)           ! layer pressures (mb) 
      real(kind=rb) :: xcldfrac(ncol,mxlay)         ! layer cloud fraction
      real(kind=rb) :: xciwp(ncol,mxlay)            ! in-cloud ice water path
      real(kind=rb) :: xclwp(ncol,mxlay)            ! in-cloud liquid water path
      real(kind=rb) :: xrei(ncol,mxlay)             ! cloud ice particle effective size (microns)
      real(kind=rb) :: xrel(ncol,mxlay)             ! cloud liquid particle effective radius (microns)
      real(kind=rb) :: xtauc(nbndsw,ncol,mxlay)     ! in-cloud optical depth (non-delta scaled)
      real(kind=rb) :: xssac(nbndsw,ncol,mxlay)     ! in-cloud single scattering albedo (non-delta scaled)
      real(kind=rb) :: xasmc(nbndsw,ncol,mxlay)     ! in-cloud asymmetry parameter (non-delta scaled)
      real(kind=rb) :: xfsfc(nbndsw,ncol,mxlay)     ! in-cloud forward scattering fraction (non-delta scaled)

      real(kind=rb) :: xcldfmc(ngptsw,ncol,mxlay)   ! cloud fraction [mcica]
      real(kind=rb) :: xciwpmc(ngptsw,ncol,mxlay)   ! in-cloud ice water path [mcica]
      real(kind=rb) :: xclwpmc(ngptsw,ncol,mxlay)   ! in-cloud liquid water path [mcica]
      real(kind=rb) :: xrelqmc(ncol,mxlay)          ! liquid particle effective radius (microns)
      real(kind=rb) :: xreicmc(ncol,mxlay)          ! ice particle effective radius (microns)
      real(kind=rb) :: xtaucmc(ngptsw,ncol,mxlay)   ! in-cloud optical depth [mcica]
      real(kind=rb) :: xssacmc(ngptsw,ncol,mxlay)   ! in-cloud single scattering albedo [mcica]
      real(kind=rb) :: xasmcmc(ngptsw,ncol,mxlay)   ! in-cloud asymmetry parameter [mcica]
      real(kind=rb) :: xfsfcmc(ngptsw,ncol,mxlay)   ! in-cloud forward scattering fraction [mcica]

      real(kind=rb) :: CST,CST0,FST,FST0
      real(kind=rb) :: CSB,CSB0,FSB,FSB0
      real(kind=rb) :: TCST(NTYP),TCSB(NTYP),TFST(NTYP),TFSB(NTYP)
      real(kind=rb) :: adjflux_jd

      LOGICAL, SAVE    :: FIRST = .TRUE.
!Yu+

!

! Initializations

      zepsec = 1.e-06_rb
      zepzen = 1.e-10_rb
      oneminus = 1.0_rb - zepsec
      pi = 2._rb * asin(1._rb)

      icpr = 0
      page = char(12)

! Set imca to select calculation type:
!  (read by subroutine readprof from input file INPUT_RRTM):  
! imca = 0, use standard forward model calculation (clear and overcast only)
! imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!           (clear, overcast or partial cloud conditions)

! Set irng to select random number generator for McICA (use when imca = 1)
! irng = 0, KISSVEC
! irng = 1, Mersenne Twister
!      irng = 0
      irng = 1

! Set icld to select of clear or cloud calculation and cloud overlap method
!  (read by subroutine readprof from input file INPUT_RRTM):  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap (McICA only)
! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
! icld = 3, with clouds using maximum cloud overlap (McICA only)

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 224 to 112 for input absorption
! coefficient data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  


      IF(FIRST) THEN
       IF(nlayers.NE.mxlay) THEN
        WRITE(6,*)"STOP: layers.NE.mxlay"
        STOP
       ENDIF

       call rrtmg_sw_ini(cpdair)
       WRITE(6,*) "run rrtmg_sw_ini"
       FIRST = .FALSE.
      ENDIF
      

! This is the main longitude/column loop within rrtmg.

!      do iplon = 1, ncol
      iplon = 1

      iout = 0
      imca = 1
!      icld = 2
      iaer = 10
      isccos = 0
      idelm = 0
      inflag = 2
      iceflag = 3
      liqflag = 1

      isolvar =  0
      solvar = 1.0
      ib1=16
      ib2=29
      if (juldat .eq. 0) then
         adjflux_jd = 1._rb
      else
         adjflux_jd = earth_sun (juldat)
      endif

      if (isolvar .eq. 0) then
         do ib = ib1,ib2
            adjflux(ib) = adjflux_jd
         enddo
      elseif (isolvar .eq. 1) then
         do ib=ib1,ib2
            adjflux(ib) = adjflux_jd * solvar(ib1)
         enddo
      elseif (isolvar .eq. 2) then
         do ib=ib1,ib2
            adjflux(ib) = adjflux_jd * solvar(ib)
         enddo
      else
         print *, 'ISOLVAR = ', isolvar, ' NOT A VALID INPUT VALUE'
         stop
      endif

      istart = jpb1                 ! jpb1 = 16
      iend = jpb2                   ! jpb2 = 29
      iflag = iout

! Return here for multiple band output
 1000 continue
      if (iflag .gt. 0 .and. iflag .le. jpb2) then
       istart = iflag
       iend = iflag
      endif


! Call sub-colum cloud generator for McICA calculations.

      tauc = 0.
      ssac = 1.
      asmc = 0.
      fsfc = 0.       

      if (imca.ne.1) then
        WRITE(6,*)"Need to check value of imca"
        STOP
      endif

      permuteseed = 1

      xpavel(1,:)=pavel(:)
      xcldfrac(1,:)=cldfrac(:)
      xciwp(1,:)=ciwp(:)
      xclwp(1,:)=clwp(:)
      xrei(1,:)=rei(:)
      xrel(1,:)=rel(:)
      xtauc(:,1,:)=tauc(:,:)
      xssac(:,1,:)=ssac(:,:)
      xasmc(:,1,:)=asmc(:,:)
      xfsfc(:,1,:)=fsfc(:,:)

      call mcica_subcol_sw(iplon,ncol, nlayers, icld, permuteseed, irng, xpavel, &
        xcldfrac, xciwp, xclwp, xrei, xrel, xtauc, xssac, xasmc, xfsfc, &
        xcldfmc, xciwpmc, xclwpmc, xreicmc, xrelqmc, xtaucmc, &
        xssacmc, xasmcmc, xfsfcmc)

      cldfmc(:,:)= xcldfmc(:,1,:)
      ciwpmc(:,:)=xciwpmc(:,1,:)
      clwpmc(:,:)= xclwpmc(:,1,:)
      reicmc(:)= xreicmc(1,:)
      relqmc(:)= xrelqmc(1,:)
      taucmc(:,:)=xtaucmc(:,1,:)
      ssacmc(:,:)=xssacmc(:,1,:)
      asmcmc(:,:)=xasmcmc(:,1,:)
      fsfcmc(:,:)=xfsfcmc(:,1,:)


!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  
!  Note: Model will be stopped if partial cloud present without McICA.

!  If McICA is requested use cloud fraction and cloud physical properties 
!  generated by sub-column cloud generator above. 

      call cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldfmc, &
                               ciwpmc, clwpmc, reicmc, relqmc, &
                               taormc, taucmc, ssacmc, asmcmc, fsfcmc)
      icpr = 1


! Calculate coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.

      call setcoef_sw(nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl, &
          laytrop, layswtch, laylow, jp, jt, jt1, &
          co2mult, colch4, colco2, colh2o, colmol, coln2o, &
          colo2, colo3, fac00, fac01, fac10, fac11, &
          selffac, selffrac, indself, forfac, forfrac, indfor)
  

! Cosine of the solar zenith angle 
!  Prevent using value of zero;

      cossza = zenith
      if (cossza.eq.0._rb) cossza = zepzen

! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 
  
! Surface albedo and aerosol optical properties
      do ib=1,nbndsw
         albdif(ib) = SALB(ib)
         albdir(ib) = SALB(ib)

         tauaer(:,ib1-1+ib)=EXT(:,ib)
         ssaaer(:,ib1-1+ib)=OMGA(:,ib)
         asmaer(:,ib1-1+ib)=G(:,ib)
      enddo

! Clouds
      if (icld.eq.0) then

         ztauc(:,:) = 0._rb
         ztaucorig(:,:) = 0._rb
         zasyc(:,:) = 0._rb
         zomgc(:,:) = 1._rb
         zcldfmc(:,:) = 0._rb
         ztaucmc(:,:) = 0._rb
         ztaormc(:,:) = 0._rb
         zasycmc(:,:) = 0._rb
         zomgcmc(:,:) = 1._rb

      elseif (icld.ge.1) then
         if (imca.eq.0) then
            do i=1,nlayers
               do ib=1,nbndsw
                  if (cldfrac(i) .ge. zepsec) then
                     ztauc(i,ib) = taucloud(i,jpb1-1+ib)
                     ztaucorig(i,ib) = taucldorig(i,jpb1-1+ib)
                     zasyc(i,ib) = asmcloud(i,jpb1-1+ib)
                     zomgc(i,ib) = ssacloud(i,jpb1-1+ib)
                  endif
               enddo
            enddo
         else
            do i=1,nlayers
               do ig=1,ngptsw
                  zcldfmc(i,ig) = cldfmc(ig,i)
                  ztaucmc(i,ig) = taucmc(ig,i)
                  ztaormc(i,ig) = taormc(ig,i)
                  zasycmc(i,ig) = asmcmc(ig,i)
                  zomgcmc(i,ig) = ssacmc(ig,i)
               enddo
            enddo
         endif   
      endif   

! Aerosol
! IAER = 0: no aerosols
      if (iaer.eq.0) then

         ztaua(:,:) = 0._rb
         zasya(:,:) = 0._rb
         zomga(:,:) = 1._rb

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer)
! or set manually here for each aerosol and layer. 
      elseif (iaer.eq.6) then
       ! deleted -- Yu

! IAER=10: Direct specification of aerosol properties from IN_AER_RRTM.
      elseif (iaer.eq.10) then

         do i = 1 ,nlayers
            do ib = 1 ,nbndsw
               ztaua(i,ib) = tauaer(i,jpb1-1+ib)
               zasya(i,ib) = asmaer(i,jpb1-1+ib)
               zomga(i,ib) = ssaaer(i,jpb1-1+ib)
            enddo
         enddo

      endif

! Call the 2-stream radiation transfer model
! total aerosol

      do i=1,nlayers+1
         zbbcu(i) = 0._rb
         zbbcd(i) = 0._rb
         zbbfu(i) = 0._rb
         zbbfd(i) = 0._rb
         zbbcddir(i) = 0._rb
         zbbfddir(i) = 0._rb
         zuvcd(i) = 0._rb
         zuvfd(i) = 0._rb
         zuvcddir(i) = 0._rb
         zuvfddir(i) = 0._rb
         znicd(i) = 0._rb
         znifd(i) = 0._rb
         znicddir(i) = 0._rb
         znifddir(i) = 0._rb
      enddo

      call spcvmc_sw &
       (nlayers, istart, iend, icpr, idelm, iout, &
        pavel, tavel, pz, tz, tbound, albdif, albdir, &
        zcldfmc, ztaucmc, zasycmc, zomgcmc, ztaormc, &
        ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
        laytrop, layswtch, laylow, jp, jt, jt1, &
        co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
        fac00, fac01, fac10, fac11, &
        selffac, selffrac, indself, forfac, forfrac, indfor, &
        zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, &
        zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir)


! Prepare output up and down, clear and total flux output
      do i = 1, nlayers+1
         totuclfl(i-1) = zbbcu(i)
         totdclfl(i-1) = zbbcd(i)
         totuflux(i-1) = zbbfu(i)
         totdflux(i-1) = zbbfd(i)
! Prepare direct/diffuse flux output
         dirdflux(i-1) = zbbfddir(i)
         difdflux(i-1) = totdflux(i-1) - dirdflux(i-1)
      enddo

! Prepare net clear and total flux output
      do i = 1, nlayers+1
         fnetc(i-1) = totdclfl(i-1) - totuclfl(i-1)
         fnet(i-1) = totdflux(i-1) - totuflux(i-1)
      enddo

! Output clear and total heating rates
      do i = 1, nlayers
         zdpgcp = heatfac / pdp(i)
         htrc(i-1) = (fnetc(i) - fnetc(i-1)) * zdpgcp
         htr(i-1) = (fnet(i) - fnet(i-1)) * zdpgcp
      enddo
      htr(nlayers) = 0._rb
      htrc(nlayers) = 0._rb

! Process output.
      CST = fnetc(nlayers)
      FST = fnet(nlayers)
      CSB = fnetc(0)
      FSB = fnet(0)

! No aeorsol RF
      IF(IFCS.GE.1) THEN
       ztaua(:,:) = 0._rb
       zasya(:,:) = 0._rb
       zomga(:,:) = 1._rb
       do i=1,nlayers+1
         zbbcu(i) = 0._rb
         zbbcd(i) = 0._rb
         zbbfu(i) = 0._rb
         zbbfd(i) = 0._rb
         zbbcddir(i) = 0._rb
         zbbfddir(i) = 0._rb
         zuvcd(i) = 0._rb
         zuvfd(i) = 0._rb
         zuvcddir(i) = 0._rb
         zuvfddir(i) = 0._rb
         znicd(i) = 0._rb
         znifd(i) = 0._rb
         znicddir(i) = 0._rb
         znifddir(i) = 0._rb
       enddo

       call spcvmc_sw &
             (nlayers, istart, iend, icpr, idelm, iout, &
              pavel, tavel, pz, tz, tbound, albdif, albdir, &
              zcldfmc, ztaucmc, zasycmc, zomgcmc, ztaormc, &
              ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
              laytrop, layswtch, laylow, jp, jt, jt1, &
              co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
              fac00, fac01, fac10, fac11, &
              selffac, selffrac, indself, forfac, forfrac, indfor, &
              zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, &
              zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir)


! Prepare output up and down, clear and total flux output
       do i = 1, nlayers+1
         totuclfl(i-1) = zbbcu(i)
         totdclfl(i-1) = zbbcd(i)
         totuflux(i-1) = zbbfu(i)
         totdflux(i-1) = zbbfd(i)
! Prepare direct/diffuse flux output
         dirdflux(i-1) = zbbfddir(i)
         difdflux(i-1) = totdflux(i-1) - dirdflux(i-1)
       enddo

! Prepare net clear and total flux output
       do i = 1, nlayers+1
         fnetc(i-1) = totdclfl(i-1) - totuclfl(i-1)
         fnet(i-1) = totdflux(i-1) - totuflux(i-1)
       enddo

! Output clear and total heating rates
       do i = 1, nlayers
         zdpgcp = heatfac / pdp(i)
         htrc(i-1) = (fnetc(i) - fnetc(i-1)) * zdpgcp
         htr(i-1) = (fnet(i) - fnet(i-1)) * zdpgcp
       enddo
       htr(nlayers) = 0._rb
       htrc(nlayers) = 0._rb

! Process output.

       CST0 = fnetc(nlayers)
       FST0 = fnet(nlayers)
       CSB0 = fnetc(0)
       FSB0 = fnet(0)

      ENDIF


      IF(IFCS.EQ.2) THEN
! RF of each type
      DO ITYP = 1,NTYP
!       do i = 1 ,nlayers
!        do ib = 1 ,nbndsw
!          ztaua(i,ib) = TEXT(i,ib,ITYP)
!          zasya(i,ib) = TOMGA(i,ib,ITYP)
!          zomga(i,ib) = TG(i,ib,ITYP)
!        enddo
!       enddo
       ztaua(:,:) = TEXT(:,:,ITYP)
       zasya(:,:) = TG(:,:,ITYP)
       zomga(:,:) = TOMGA(:,:,ITYP)

       do i=1,nlayers+1
         zbbcu(i) = 0._rb
         zbbcd(i) = 0._rb
         zbbfu(i) = 0._rb
         zbbfd(i) = 0._rb
         zbbcddir(i) = 0._rb
         zbbfddir(i) = 0._rb
         zuvcd(i) = 0._rb
         zuvfd(i) = 0._rb
         zuvcddir(i) = 0._rb
         zuvfddir(i) = 0._rb
         znicd(i) = 0._rb
         znifd(i) = 0._rb
         znicddir(i) = 0._rb
         znifddir(i) = 0._rb
       enddo

       call spcvmc_sw &
        (nlayers, istart, iend, icpr, idelm, iout, &
        pavel, tavel, pz, tz, tbound, albdif, albdir, &
        zcldfmc, ztaucmc, zasycmc, zomgcmc, ztaormc, &
        ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
        laytrop, layswtch, laylow, jp, jt, jt1, &
        co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
        fac00, fac01, fac10, fac11, &
        selffac, selffrac, indself, forfac, forfrac, indfor, &
        zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, &
        zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir)


! Prepare output up and down, clear and total flux output
       do i = 1, nlayers+1
         totuclfl(i-1) = zbbcu(i)
         totdclfl(i-1) = zbbcd(i)
         totuflux(i-1) = zbbfu(i)
         totdflux(i-1) = zbbfd(i)
! Prepare direct/diffuse flux output
         dirdflux(i-1) = zbbfddir(i)
         difdflux(i-1) = totdflux(i-1) - dirdflux(i-1)
       enddo

! Prepare net clear and total flux output
       do i = 1, nlayers+1
         fnetc(i-1) = totdclfl(i-1) - totuclfl(i-1)
         fnet(i-1) = totdflux(i-1) - totuflux(i-1)
       enddo

! Output clear and total heating rates
       do i = 1, nlayers
         zdpgcp = heatfac / pdp(i)
         htrc(i-1) = (fnetc(i) - fnetc(i-1)) * zdpgcp
         htr(i-1) = (fnet(i) - fnet(i-1)) * zdpgcp
       enddo
       htr(nlayers) = 0._rb
       htrc(nlayers) = 0._rb

! Process output.
       TCST(ITYP) = fnetc(nlayers)
       TFST(ITYP) = fnet(nlayers)
       TCSB(ITYP) = fnetc(0)
       TFSB(ITYP) = fnet(0)
      ENDDO
      ENDIF
             
!       if(II.eq.60.and.JJ.eq.33) then
!        WRITE(1003,*)'rrtmg_sw'
!        do i = 1, nlayers
!           WRITE(1003,106)i,EXT(i,6),OMGA(i,6),G(i,6), &
!              TEXT(i,6,1),TOMGA(i,6,1),TG(i,6,1)
!        enddo
!        WRITE(1003,107)'CST,CST0,CST-CST0',CST,CST0,CST-CST0
!        WRITE(1003,107)'FST,FST0,FST-FST0',FST,FST0,FST-FST0
!        WRITE(1003,107)'FST-CST,FST0-CST0',FST-CST,FST0-CST0
!        DO ITYP = 1,NTYP
!         WRITE(1003,107)'TFST,FST0,TFST-FST0',TFST(ITYP),FST0, &
!             TFST(ITYP)-FST0
!         WRITE(1003,107)'TCST,CST0,TCST-CST0',TCST(ITYP),CST0, &
!             TFST(ITYP)-FST0
!        ENDDO
! 106    FORMAT(I3,10(F10.4))
! 107    FORMAT(A17,2x,10(F10.2))
! 109    FORMAT(I4,2x,10(F10.3))
!       endif

 2000 continue

      if (iout .ge. 0 .and. iout .le. jpb2) goto 3500
      if (iflag .eq. 98) then
         iflag = jpb1
      elseif (iflag .gt. 0 .and. iflag .lt. jpb2) then
         iflag = iflag + 1
      else
         goto 3500
      endif
      goto 1000
 3500 continue

! Output module version numbers

 4000 continue

! End longitude/column loop
!      enddo


 9899 format(1x,'Wavenumbers: ',f6.0,' - ',f6.0,' cm-1, ATM ',i6)
 9900 format(1x,'LEVEL PRESSURE   UPWARD FLUX   DIFDOWN FLUX  DIRDOWN FL&
     &UX  DOWNWARD FLUX   NET FLUX    HEATING RATE')
 9901 format(1x,'         mb          W/m2          W/m2          W/m2&
     &        W/m2          W/m2       degree/day')
 9902 format(1x,i3,3x,f11.6,4x,1p,2(g12.6,2x),g13.6,3x,g16.9,0p)
 9903 format(a)
 9910 format('  Modules and versions used in this calculation:',/,/, &
              7(5x,a20,2x,a18,10x,a20,2x,a18,/))

      contains

!*************************************************************************
      real(kind=rb) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      use parkind, only : im => kind_im, rb => kind_rb
      use rrsw_con, only : pi

      implicit none

      integer(kind=im), intent(in) :: idn

      real(kind=rb) :: gamma

      gamma = 2._rb*pi*(idn-1)/365._rb

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_rb + .034221_rb * cos(gamma) + .001289_rb * sin(gamma) + &
                   .000719_rb * cos(2._rb*gamma) + .000077_rb * sin(2._rb*gamma)

      end function earth_sun


      end subroutine rrtmg_sw

!EOC
      END MODULE rrtmg_sw_GCAPM
