!  CVS: $Id: pconst_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module pconst_mod
!
#include <def-undef.h>
use precision_mod
use param_mod
!     -----------------------------------------------------------
!     Index Fields
!     -----------------------------------------------------------
!#ifdef SPMD
!YU   real,dimension(:),allocatable:: dyr_global
      integer,dimension(:),allocatable:: j_global
      integer,dimension(:),allocatable:: i_global
      integer :: ix,iy
      real(r8),dimension(:,:,:),allocatable:: vit_global
      real(r8),dimension(:,:),allocatable:: vit_global_surface
      real(r8),dimension(:,:),allocatable:: viv_global_surface
!#endif
      real(r8),dimension(:,:,:),allocatable:: vit,viv 
      real(r8),dimension(:,:,:),allocatable:: vit_1d,viv_1d
      real(r8),dimension(:),allocatable:: ahv_back
      integer,dimension(:,:),allocatable:: basin
      integer,dimension(:,:),allocatable:: itnu
!lhl1204
      integer,dimension(:,:),allocatable:: na
      real(r8) :: dfricmx,dwndmix
      ! lihuimin, 2012.7.15
      integer :: i_num   ! actual grid number in i direction in this process
      integer :: j_num   ! actual grid number in j direction in this process
      integer :: i_f_num ! formal grid number in i dircetion 
      integer :: j_f_num ! formal grid number in j direction 
!lhl1204
!
!
!     -----------------------------------------------------------
!     Grids
!     -----------------------------------------------------------
#if (defined NETCDF) || (defined ALL)
      real(r4),dimension(:),allocatable:: lon
      real(r4),dimension(:),allocatable:: lat
      real(r4),dimension(:),allocatable:: lev
      real(r4),dimension(:),allocatable:: lev1
#endif
!lhl090729
      real(r8),dimension(:),allocatable:: s_lon
      real(r8),dimension(:),allocatable:: s_lat
!lhl090729
      real(r8),dimension(:),allocatable:: zkt,dzp,odzp,odzt
      real(r8),dimension(:),allocatable:: zkp


      real(r4),dimension(:),allocatable::DYR_IN
      real(r8),dimension(:),allocatable::OUY,OTX,OUX,SOTX,SOUX, &
                    FF,CV1,CV2,SNLAT,SINT, &
                    SINU,DYT,DYR,DXDYU,DXDYT, &
                    R1A,R1B,R2A,R2B,R1C, &
                    R1D,R2C,R2D,EBEA,EBEB, &
                    EBLA,EBLB,EPEA,EPEB,EPLA,EPLB
!lhl060506
      real(r8),dimension(:),allocatable:: FF1,RRD1,RRD2
!lhl060506

!XC
#if (defined TSPAS)
      real(r8),dimension(:),allocatable:: dtdy,dtdx,RAA,RBB
#endif
!XC

#ifdef SPMD
!YU
!real,dimension(:),allocatable::OUY_global,OTX_global,OUX_global,SOTX_global,SOUX_global,
!&
!YU                 FF_global,CV1_global,CV2_global,SNLAT_global,SINT_global, &
!YU                 SINU_global,DYT_global,DYR_global,DXDYU_global,DXDYT_global,
!&
!YU                 R1A_global,R1B_global,R2A_global,R2B_global,R1C_global, &
!YU                 R1D_global,R2C_global,R2D_global,EBEA_global,EBEB_global, &
!YU
!EBLA_global,EBLB_global,EPEA_global,EPEB_global,EPLA_global,EPLB_global
      real(r4),dimension(:),allocatable::DYR_IN_global
      real(r8),dimension(:),allocatable::OUY_global,OTX_global,OUX_global,SOTX_global,SOUX_global, &
                    FF_global,CV1_global,CV2_global,SNLAT_global,SINT_global, &
                    SINU_global,DYT_global,DYR_global,DXDYU_global,DXDYT_global, &
                    R1A_global,R1B_global,R2A_global,R2B_global,R1C_global, &
                    R1D_global,R2C_global,R2D_global,EBEA_global,EBEB_global, &
                    EBLA_global,EBLB_global,EPEA_global,EPEB_global,EPLA_global,EPLB_global

!lhl060506
     real(r8),dimension(:),allocatable:: FF1_global
!lhl060506
!XC
#if (defined TSPAS)
      real(r8),dimension(:),allocatable:: dtdy_global,dtdx_global,RAA_global,RBB_global
#endif
!XC

#endif



#if ( defined SMAG)
      real(r8),dimension(:),allocatable:: CXT,CXU,CYT,CYU &
                   ,R1E,R1F,R2E,R2F &
                   ,R3E,R3F,R4E,R4F
#ifdef SPMD
      real(r8),dimension(:),allocatable:: CXT_global,CXU_global,CYT_global,CYU_global &
                   ,R1E_global,R1F_global,R2E_global,R2F_global &
                   ,R3E_global,R3F_global,R4E_global,R4F_global

#endif
#endif

!Yu

      real(r8),dimension(:,:),allocatable::ohbt,ohbu,dzph,hbx,hby
      real(r8),dimension(:),allocatable:: COSU,COST
#ifdef SPMD
       real(r8),dimension(:),allocatable::COSU_global,COST_global
       integer,dimension(:),allocatable:: i_start,j_start
#endif
      real(r8),dimension(:),allocatable::CF1,CF2,SF1,SF2
!
!
!     -----------------------------------------------------------
!     Reference T S & coefficients for calculation of d(density)
!     -----------------------------------------------------------
!YU
      REAL(r8),allocatable:: TO(:),SO(:),C(:,:),&
!lhl1204
                 PO(:)
!lhl1204

!YU
!
!
!     -----------------------------------------------------------
!     Control Parameter
!     -----------------------------------------------------------
      INTEGER:: ISOP
!
!
!     -----------------------------------------------------------
!     Phycical Parameter
!     -----------------------------------------------------------
      real(r8)::amv,ahv,ahice
!lhl1204
      real(r8),dimension(:,:,:),allocatable::akmu,akmt
      real(r8),dimension(:,:,:,:),allocatable::akt
!lhl1204
      real(r8),dimension(:),allocatable::AM,AH
      real(r8),dimension(:,:,:),allocatable::am3,ah3
!lhl
      real(r8),dimension(:,:,:),allocatable::amx,amy
!lhl
      real(r8)::gamma
#if ( defined SMAG)
      real(r8):: D0,CP,G,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
                    VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2,KARMAN,RR
#else
      real(r8)::  D0,CP,G,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
                    VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2
#endif
!
!
      CHARACTER (LEN=3):: ABMON(12)
      CHARACTER (LEN=3):: ABMON1(12)
      CHARACTER (LEN=80):: out_dir
!
      REAL(r8):: DTB,DTC,DTS,DTB2,DTC2,ONBB,ONBC,ONCC
      INTEGER:: NBB,NCC,NSS,ISB,ISC,IST,MONTH
      INTEGER:: NMONTH(12),NNMONTH(12)
      INTEGER:: number_day, number_month
!
      ! lihuimin 2012.6.18, add REFDATE
      INTEGER :: NUMBER,NSTART,IY0,IYFM,MON0,MEND,IMD,IDAY,II,JJ,IO_HIST,IO_REST,rest_freq,hist_freq,REFDATE
      integer :: klv
!
      REAL(r8):: PI,RADIUS,TORAD
!lhl20130419
    integer(kind(1))   :: curr_ymd_licom     ! Current date YYYYMMDD
    integer(kind(1))   :: curr_ymd_cpl     ! Current date YYYYMMDD
    integer(kind(1))   :: yy_licom,mm_licom,dd_licom,tod_licom     ! year, month, day
    integer(kind(1))   :: yy_cpl,mm_cpl,dd_cpl,tod_cpl     ! year, month, day
!lhl20130419
    integer:: ocn_cpl_dt

!jjb20181112  diurnal cycle
    integer:: first_step
    integer:: nstep_per_day
    integer:: totalstep   !the number of adding in a month for output
!jjb20181112 
end module pconst_mod
