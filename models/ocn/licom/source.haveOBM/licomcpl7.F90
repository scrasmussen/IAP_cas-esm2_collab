#define LOGMSG()
!write(mytid+600,'(a,i4)')"ocn",__LINE__
!#define LOGMSGCarbon()
!write(600+mytid,'(a,i4)')"ocn CARBON",__LINE__


module ocn_comp_mct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !MODULE:    ocn_COMP_MCT
!   !AUTHOR:    Huimin Li
!   !Date:      2012/5/30
!
!   !DESCRIPTION:
!               This is the main driver for the ocn coupled in CPL7
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use mct_mod
   use esmf_mod
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use shr_file_mod
   use shr_cal_mod, only : shr_cal_date2ymd
   use shr_sys_mod
   use shr_const_mod,only:SHR_CONST_SPVAL !linpf 2012Jul26
   use perf_mod
   use fluxcpl
   use POP_CplIndices
   use shr_dmodel_mod

#include <def-undef.h>
use param_mod
use pconst_mod

use shr_msg_mod
use shr_sys_mod
use control_mod
use constant_mod, only : LATVAP
use shr_cal_mod,       only: shr_cal_date2ymd


#if ( defined SPMD ) || ( defined COUP)
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
use tracer_mod
use pmix_mod
use forc_mod
#ifdef USE_OCN_CARBON
use carbon_mod
use cforce_mod
use grids_pt_mod !lyc 2014.06
#endif


  implicit none
#include <netcdf.inc>
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  SAVE
  private



  private :: ocn_export_mct
  private :: ocn_import_mct
  private :: ocn_SetGSMap_mct
  private :: ocn_domain_mct


  type(seq_infodata_type), pointer :: &
     infodata

!==========================================================================
  contains
!==========================================================================


!**************************************************************************
!   !ROUTINE:   ocn_init_mct
!   !AUTHOR:    Huimin Li
!   !Date:      2012/5/30
!
!   !INTERFACE:
  !subroutine ocn_init_mct(EClock, cdata_o, x2o_o, o2x_o,r2x_o, NLFilename)
  subroutine ocn_init_mct(EClock, cdata_o, x2o_o, o2x_o, NLFilename) !LPF 20121219
!
!   !DESCRIPTION:
!               This is the initialize routine for ocn coupled in CPL7
!   !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
    !type(mct_aVect)             , intent(inout) :: r2x_o
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!   !LOCAL VARIABLES
    type(mct_gsMap), pointer :: &
       gsMap_o

    type(mct_gGrid), pointer :: &
       dom_o

    integer (kind(1)) :: &
       nThreads

    integer (kind(1)) :: &
       OCNID,     &
       lsize

    real (r8) ::  &
       precadj

    integer (kind(1)) :: iam,ierr
    character(len=32)  :: starttype          ! infodata start type

    integer :: i_temp, j_temp, i_comp, j_comp ! used in specifying the grid number in each process

  !------------------------------------------------------
  ! ROUTINE BEGIN
  !------------------------------------------------------
  mpi_comm_ocn=0
!wangty bug
   ISB = 0
   ISC = 0
   IST = 0
!wangty
  ! lihuimin 2012.7.16
  ! OCNID, errcode
  call seq_cdata_setptrs(cdata_o, ID=OCNID, mpicom=mpi_comm_ocn, &
       gsMap=gsMap_o, dom=dom_o, infodata=infodata)

  ! five parameters initialized in msg_pass('init') in CPL6 coupled version
  ! get infodata from drv
  cdate = 0
  sec = 0
  ierr = 0
  info_time = 0
  call seq_infodata_GetData( infodata, info_debug=info_dbug)

  ! send initial state to drv
  call seq_infodata_PutData( infodata, ocn_nx=(imt_global-2), ocn_ny=jmt_global)
!  call seq_infodata_PutData( infodata, ocn_prognostic=.true.) !LPF 20120829
  ! lihuimin, 2012.7.25, not use ocnrof_p
  call seq_infodata_PutData(infodata, ocn_prognostic=.true., ocnrof_prognostic=.true., rof_present=.true.)
!open ocnrof. !LPF 20120829

  !-----------------------------
  ! namelist & logfile setup
  !-----------------------------
  ! TODO: redirect the log file


  mytid=0

! He - 2010-10-08 | NDay: Number of days since the first day when the simulation begins.
! lihuimin cancle , 2012.6.14
!     NDAY = 0



!
! mpicom_o <-> mpi_comm_ocn
! TODO
      if (mytid==0) write(6,*)"Begin mpi_comm_rank"
      call mpi_comm_rank (mpi_comm_ocn, mytid, ierr)
      if (mytid==0)  write(6,*)"End mpi_comm_rank"
      call mpi_comm_size (mpi_comm_ocn, nproc, ierr)
      write(6,*) "MYTID=",mytid,"Number of Processors is",nproc


   ! lihuimin 2012.6.14, initial indices
   call POP_CplIndicesSet()


!---------------------------------------------------------------------
!     SET THE CONSTANTS USED IN THE MODEL
!---------------------------------------------------------------------
#ifdef COUP
      call shr_sys_flush(6)
#endif
    LOGMSG()
      CALL CONST
    LOGMSG()
      if (mytid == 0) then
      write(111,*)"OK------3"
      close(111)
      end if
#ifdef COUP
      call shr_sys_flush(6)
#endif
#ifdef SHOW_TIME
      call run_time('CONST')
#endif

!***********************************************************************
!          SET SOME CONSTANTS FOR THE BOGCM
!***********************************************************************
#ifdef USE_OCN_CARBON
      CALL CTRLC
#ifdef SHOW_TIME
      call run_time('CTRLC')
#endif
#endif

!---------------------------------------------------------------------
!     SET MODEL'S RESOLUTION,TOPOGRAPHY AND THE CONSTANT
!     PARAMETERS RELATED TO LATITUDES (J)
!---------------------------------------------------------------------
    LOGMSG()
      CALL GRIDS
      if (mytid == 0) then
      write(111,*)"OK------4"
      close(111)
      end if
#ifdef SHOW_TIME
      call run_time('GRIDS')
#endif
!--------------------------------------------------------------------
!lyc 2014.06
#ifdef USE_OCN_CARBON
      CALL GRIDS_PT
#ifdef SHOW_TIME
      call run_time('GRIDS_PT')
#endif
#endif


!---------------------------------------------------------------------
!     SET SURFACE FORCING FIELDS (1: Annual mean; 0: Seasonal cycle)
!---------------------------------------------------------------------
    LOGMSG()
      CALL RDRIVER
#ifdef SHOW_TIME
      call run_time('RDRIVER')
#endif
      if (mytid == 0) then
      write(111,*)"OK------5"
      close(111)
      end if

!***********************************************************************
!      SET FORCING DATA USED IN CARBON CYCLYE
!***********************************************************************
#ifdef USE_OCN_CARBON
       CALL CFORCE
#ifdef SHOW_TIME
       call run_time('CFORCE')
#endif
#endif

!---------------------------------------------------------------------
!     INITIALIZATION
!---------------------------------------------------------------------
    LOGMSG()
      CALL INIRUN
      if (mytid == 0) then
      write(111,*)"OK------5.0"
      close(111)
      end if
#ifdef SHOW_TIME
      call run_time('INRUN')
#endif

!----------------------------------------------------------------------
!    specify the actual grid number in each process
!    i_num, j_num, i_f_num, j_f_num
!    lihuimin , 2012.7.15
!----------------------------------------------------------------------

!   ! i direciton
!   i_temp = (imt_global-2)/nx_proc
!   i_comp = (imt_global-2) - i_temp*nx_proc
!   if (i_comp == 0) then
!      i_num = i_temp
!      i_f_num = i_temp
!   else ! i_comp > 0
!      i_f_num = i_temp+1
!      if (ix == nx_proc-1) then
!         i_num = (imt_global-2) - ix*i_f_num
!      else
!         i_num = i_f_num
!      endif
!   endif
!      
!   ! j direction
!   j_temp = jmt_global/ny_proc
!   j_comp = jmt_global - j_temp*ny_proc
!   if (j_comp == 0) then
!      j_num = j_temp
!      j_f_num = j_temp
!   else
!      j_f_num = j_temp+1
!      if (iy == ny_proc-1) then
!         j_num = jmt_global - iy*j_f_num
!      else
!         j_num = j_f_num
!      endif
!   endif
   i_f_num = imt - num_overlap
   j_f_num = jmt - num_overlap
   ! i-direction
   if (ix == nx_proc-1) then
      i_num = (imt_global-2) - ix*i_f_num
   else
      i_num = i_f_num
   endif
   ! j-direction
   if (iy == ny_proc-1) then
      !j_num = jmt_global - iy*j_f_num
      j_num = jmt_global - (iy-1)*j_f_num - (j_f_num + (jst_global - 1) + 1)
   elseif (iy == 0) then
      j_num = j_f_num + (jst_global - 1) + 1
    else
      j_num = j_f_num
   endif


!----------------------------------------------------------------------
!     Inialize mct attribute vectors
!     imitate CPL7/POP
!     lihuimin 2012.7.16
!----------------------------------------------------------------------

   call ocn_SetGSMap_mct(mpi_comm_ocn, OCNID, GSMap_o)

   lsize = mct_gsMap_lsize(gsMap_o, mpi_comm_ocn)
   write(6,*) "mct_gsMap_lsize = ",lsize
   write(6,*) "ID=",mytid,"i_num=",i_num,"j_num=",j_num

   call ocn_domain_mct(lsize, gsMap_o, dom_o)

   call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
   call mct_aVect_zero(x2o_o)


   call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize)
   call mct_aVect_zero(o2x_o)
!LPF 20121219
!   call mct_aVect_init(r2x_o, rList=seq_flds_r2x_fields, lsize=lsize)
!   call mct_aVect_zero(r2x_o)
!LPF 20121219

!**********************************************************************
!      INITIALIZATION CARBON MODEL
!**********************************************************************
#ifdef USE_OCN_CARBON
      CALL INIRUN_PT
#ifdef SHOW_TIME
      call run_time('INIRUN_PT')
#endif
#endif
#ifdef CANUTO
      call turb_ini
#endif
!---------------------------------------------------------------------
!     INITIALIZATION OF ISOPYCNAL MIXING
!---------------------------------------------------------------------
    LOGMSG()
#ifdef ISO
      CALL ISOPYI
#ifdef SHOW_TIME
      call run_time('ISOPYI')
#endif
#endif

#ifdef COUP
         call shr_sys_flush(6)
#endif
      call ocn_export_mct(o2x_o)


  end subroutine ocn_init_mct


!**************************************************************************
!   !ROUTINE:   ocn_run_mct
!   !AUTHOR:    Huimin Li
!   !Date:      2012/6/2
!
!   !INTERFACE:
!  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o,r2x_o)
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
!
!   !DESCRIPTION:
!   !run ocn for a coupling interval
!
!   !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
!    type(mct_aVect)             , intent(inout) :: r2x_o

    integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD
    integer(kind(1))   :: yy,mm,dd     ! year, month, day

!----------------------------------------------------------------------
!   get Eclock data information
!   lihuimin, 2012.7.20
!----------------------------------------------------------------------
    call seq_timemgr_EClockGetData( EClock, curr_ymd=curr_ymd)
    call shr_cal_date2ymd(curr_ymd,yy,mm,dd)
!    imd = 30
    iyfm = yy  ! nwmf = iyfm
    mon0 = mm  ! month = (iyfm-1)*12 + mon0
    iday = dd  ! number_day = iday + 1
    
    if(mytid == 0) then
       write(6,*) "From CPL7-Time yy=",yy,"mm=",mm,"dd=",dd
    endif

!linpf 2012Jul27
!=====================time control for ocn output===============
!lhl20120728
       month=(IYFM-1)*12+mon0
       IMD = NMONTH (MON0)
!lhl20120728
!      month=mm
!      IY0 = (MONTH -1)/12
!      IYFM = IY0+1
!!     IYFM is the number of the current year
!
!      MON0 = MONTH - IY0*12
!      IMD = NMONTH (MON0)
!!=================================================================
!linpf 2012Jul27

!----------------------------------------------------------------------
!   call msg_pass('recv')
!----------------------------------------------------------------------
    LOGMSG()
    !call ocn_import_mct(x2o_o,r2x_o)
    call ocn_import_mct(x2o_o)

    LOGMSG()
    call post_cpl
!lyc for offline version
#if (!defined COUP)
     CALL INTFOR
!lyc 2014.6.27
     call core_daily

!lyc 2014.06
#ifdef USE_OCN_CARBON
    CALL INTFOR_PT
#endif
#endif
!---------------------------------------------------------------------
!     THERMAL CYCLE
!---------------------------------------------------------------------
    LOGMSG()
         DO II = 1,NSS


!     COMPUTE DENSITY, BAROCLINIC PRESSURE AND THE RELAVANT VARIABLES
    LOGMSG()
            CALL READYT

!----------------------------------------------------
!     BAROCLINIC & BAROTROPIC CYCLE
!---------------------------------------------------------------------
    LOGMSG()
            DO JJ = 1,NCC

!     COMPUTE MOMENTUM ADVECTION, DIFFUSION & THEIR VERTICAL INTEGRALS
    LOGMSG()
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run5-"
               CALL READYC
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run5"
!         CALL ENERGY

!     PREDICTION OF BAROTROPIC MODE
    LOGMSG()
               CALL BAROTR
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run6"
!         CALL ENERGY


!     PREDICTION OF BAROCLINIC MODE
    LOGMSG()
               CALL BCLINC

!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run7"
!         CALL ENERGY

            END DO

!*******************************************************************
!     PREDICTION OF PASSIVE TRACER
!*******************************************************************
!lyc 2014.06
!       FOR BGC
!-------------------------------------------------------------------
#ifdef USE_OCN_CARBON
      CALL PTRACER
#ifdef SHOW_TIME
      call run_time('PTRACER')
#endif
#endif

    LOGMSG()
            CALL TRACER
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run8"

    LOGMSG()
            CALL ICESNOW
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run9"

!***********************************************************************
!     PERFORM CONVECTIVE ADJUSTMENT IF UNSTABLE STRATIFICATION OCCURS
!************************************************************************
!lyc 2012.06
#ifdef USE_OCN_CARBON
       CALL CONVADJ_PT
#else
!lhl1204
!#if (!defined CANUTO)
    LOGMSG()
            CALL CONVADJ
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run10"
!#endif
#endif
!*************************************************************************
!
         END DO
    LOGMSG()
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run11-"
         CALL ENERGY
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run11"

!     COMPENSATE THE LOSS OF GROSS MASS

    LOGMSG()
         CALL ADDPS
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run12"

!     MONITOR THE MODEL INTEGRATION

!         CALL ACCUMM

!*********************************************************************
!      ACCUMULATE SOME VARIABLES IN CARBON MODEL
!*********************************************************************

#ifdef COUP
         call shr_sys_flush(6)
#endif
!
#ifdef COUP
    LOGMSG()
         call flux_cpl
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run13"
#endif

    LOGMSG()
    CALL SSAVEINS
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run14"

!     ACCUMULATE SOME VARIABLES FOR MONTHLY OUTPUT

    CALL ACCUMM
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run15"
!-----------------------------------------------------------------
!lyc 2014.06
#ifdef USE_OCN_CARBON
      CALL ACCUMM_PT
#endif

    if (iday==imd) then
      CALL SSAVEMON
!-----------------------------------------------------------------
!lyc 2014.06
#ifdef USE_OCN_CARBON
      CALL SSAVE_PT
#endif
    endif
!----------------------------------------------------------------------
!   call msg_pass('send')
!----------------------------------------------------------------------
    LOGMSG()
    call ocn_export_mct(o2x_o)
!      if (mytid == 0) write(*,*)ii,jj,NCC,"OK------run16"

    LOGMSG()

  end subroutine ocn_run_mct


  subroutine ocn_final_mct

    call cleanup
    LOGMSG()

  end subroutine ocn_final_mct


!*************************************************************************
!ROUTINE: ocn_import_mct
!         lihuimin 2012.7.16
!
  subroutine ocn_import_mct(x2o_o) !LPF 20121219
    type(mct_aVect)   , intent(inout) ::  x2o_o
    ! local
    integer :: j_begin, j_end

     n=0
      ! lihuimin, 2012.8.7, consider jst_global
    if (iy == 0) then
       j_begin = 1+1 - 1 ! no overlap to the north
       j_end = j_num - (jst_global-1) !+1 !LPF 20120817 ! j_end = jmt - 1, j_num = jmt+1, jst_global = 3 here
       do j=1,jst_global-1 
       do i=2,i_num+1 !LPF 20120818 
          n=n+1
       enddo
       enddo
    else
       j_begin = 1+1
       j_end = j_num + 1
    endif


       ! do j=1+1,j_num+1
        do j=j_begin,j_end
        do i=1+1,i_num+1 !2,i_num+1 !LPF 20120818  
           n=n+1
           !--- states ---
           ifrac(i,j) = x2o_o%rAttr(index_x2o_Si_ifrac,n) ! ice fraction
           patm (i,j) = x2o_o%rAttr(index_x2o_Sa_pslv,n)  ! sea level pressure index_x2o_Sa_pslv 
           !--- fluxes ---
           taux (i,j) = x2o_o%rAttr(index_x2o_Foxx_taux,n)  ! surface stress, zonal
           tauy (i,j) = x2o_o%rAttr(index_x2o_Foxx_tauy,n)  ! surface stress, merid
           netsw(i,j) = x2o_o%rAttr(index_x2o_Foxx_swnet,n) ! net sw rad
           sen  (i,j) = x2o_o%rAttr(index_x2o_Foxx_sen,n)   ! sensible
           lwup (i,j) = x2o_o%rAttr(index_x2o_Foxx_lwup,n)  ! long-wave up
           lwdn (i,j) = x2o_o%rAttr(index_x2o_Foxx_lwdn,n)  ! long-wave down
           melth(i,j) = x2o_o%rAttr(index_x2o_Foxx_melth,n) ! melt heat
           salt (i,j) = x2o_o%rAttr(index_x2o_Foxx_salt,n)  ! salinity flux
           prec (i,j) = x2o_o%rAttr(index_x2o_Foxx_prec,n)  !index_x2o_Foxx_prec 
           evap (i,j) = x2o_o%rAttr(index_x2o_Foxx_evap,n)  ! evaporation
           meltw(i,j) = x2o_o%rAttr(index_x2o_Foxx_meltw,n) ! melt water
           roff (i,j) = x2o_o%rAttr(index_x2o_Forr_roff,n)  ! runoff  !LPF 20121219
           duu10n(i,j) = x2o_o%rAttr(index_x2o_So_duu10n,n)  ! 10m wind speed squared
        end do
        end do

        lat1= LATVAP*evap ! latent (derive from evap)

        where(vit(:,:,1)<0.5) patm=0.0D0
        where(vit(:,:,1)<0.5) ifrac=0.0D0
        where(vit(:,:,1)<0.5) taux=0.0D0
        where(vit(:,:,1)<0.5) tauy=0.0D0
        where(vit(:,:,1)<0.5) netsw=0.0D0
        where(vit(:,:,1)<0.5) sen=0.0D0
        where(vit(:,:,1)<0.5) lwup=0.0D0
        where(vit(:,:,1)<0.5) lwdn=0.0D0
        where(vit(:,:,1)<0.5) melth=0.0D0
        where(vit(:,:,1)<0.5) salt=0.0D0
        where(vit(:,:,1)<0.5) prec=0.0D0
        where(vit(:,:,1)<0.5) evap=0.0D0
        where(vit(:,:,1)<0.5) meltw=0.0D0
        where(vit(:,:,1)<0.5) roff=0.0D0
        where(vit(:,:,1)<0.5) duu10n=0.0D0

        call exchange_2d(taux,1,1)
        call exchange_2d(tauy,1,1)

!      call chk_var2d(patm,"11",1)
!      call chk_var2d(ifrac,"22",1)
!      call chk_var2d(taux,"33",1)
!      call chk_var2d(tauy,"44",1)
!      call chk_var2d(netsw,"55",1)
!      call chk_var2d(sen,"66",1)
!      call chk_var2d(lwup,"77",1)
!      call chk_var2d(lwdn,"88",1)
!      call chk_var2d(melth,"99",1)
!      call chk_var2d(salt,"00",1)
!      call chk_var2d(prec,"12",1)
!      call chk_var2d(evap,"13",1)
!      call chk_var2d(meltw,"14",1)
!      call chk_var2d(roff,"15",1)
!      call chk_var2d(duu10n,"16",1)
  end subroutine ocn_import_mct



!*************************************************************************
!ROUTINE: ocn_export_mct
!         lihuimin 2012.7.16
!
  subroutine ocn_export_mct(o2x_o)
    type(mct_aVect)   , intent(inout) :: o2x_o
    ! local
    integer :: j_begin,j_end


       n=0
! lihuimin, 2012.8.7, consider jst_global
    if (iy == 0) then
       j_begin = 1+1 - 1 ! no overlap to the north
       j_end = j_num - (jst_global-1) ! j_end = jmt - 1, j_num = jmt+1, jst_global = 3 here
       do j=1,jst_global-1
       do i=1+1,i_num+1 
          n=n+1
          o2x_o%rAttr(index_o2x_So_t,n)    = T_cpl   (i,1) ! temperature
          o2x_o%rAttr(index_o2x_So_s,n)    = S_cpl   (i,1) ! salinity
          o2x_o%rAttr(index_o2x_So_u,n)    = U_cpl   (i,1) ! velocity, zonal
          o2x_o%rAttr(index_o2x_So_v,n)    = V_cpl   (i,1) ! velocity, meridional
          o2x_o%rAttr(index_o2x_So_dhdx,n) = dhdx(i,1) ! surface slope, zonal
          o2x_o%rAttr(index_o2x_So_dhdy,n) = dhdy(i,1) ! surface slope, meridional
          o2x_o%rAttr(index_o2x_Fioo_q,n)    = q   (i,1) ! heat of fusion xor melt pot
#ifdef USE_OCN_CARBON
!         buffs(n,cpl_fields_o2c_co2)    = co2_cpl(i,j) ! state: air-sea CO2 of ocean  ~ mol
#endif
       enddo
       enddo
    else
       j_begin = 1+1
       j_end = j_num + 1
    endif


      !do j=1+1,j_num+1
      do j=j_begin,j_end
      !do i=1+1,i_num+1
      do i=2,i_num+1 !LPF 20120819
         n=n+1
         o2x_o%rAttr(index_o2x_So_t,n)    = T_cpl   (i,j) ! temperature
         o2x_o%rAttr(index_o2x_So_s,n)    = S_cpl   (i,j) ! salinity
         o2x_o%rAttr(index_o2x_So_u,n)    = U_cpl   (i,j) ! velocity, zonal
         o2x_o%rAttr(index_o2x_So_v,n)    = V_cpl   (i,j) ! velocity, meridional
         o2x_o%rAttr(index_o2x_So_dhdx,n) = dhdx(i,j) ! surface slope, zonal
         o2x_o%rAttr(index_o2x_So_dhdy,n) = dhdy(i,j) ! surface slope, meridional
         o2x_o%rAttr(index_o2x_Fioo_q,n)    = q   (i,j) ! heat of fusion xor melt pot
#ifdef USE_OCN_CARBON
!         buffs(n,cpl_fields_o2c_co2)    = co2_cpl(i,j) ! state: air-sea CO2 of ocean  ~ mol
#endif
      end do
      end do

  end subroutine ocn_export_mct


!*************************************************************************
!ROUTINE: ocn_SetGSMap_mct
!         lihuimin 2012.7.16
!
  subroutine ocn_SetGSMap_mct(mpicom_ocn, OCNID, gsMap_ocn)

    implicit none
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: OCNID
    type(mct_gsMap), intent(inout) :: gsMap_ocn

    integer,allocatable :: &
      gindex(:)

    integer (kind(1)) ::   &
      i,j, k, n, iblock, &
      lsize, gsize,   &
      ier


!    lsize = (imt)*(jmt)
    lsize = i_num*j_num
    gsize = (imt_global-2)*jmt_global
    allocate(gindex(lsize),stat=ier)

    n = 0
!    do j=1,j_num
    do j=j_num,1,-1
       do i=1,i_num !LPF 20120819
          n=n+1
          ! TODO, use num_overlap
          ! ix,iy start from 0
          !gindex(n) = iy*(imt_global-2)*(j_f_num) + (j-1)*(imt_global-2) + ix*(i_f_num) + i
          if (iy == ny_proc-1) then
             gindex(n) = (j-1)*(imt_global-num_overlap) + ix*(i_f_num) + i
          else
             !gindex(n) = (jmt_global-(ny_proc-1)*(j_f_num))*(imt_global-2) + (ny_proc-iy-1-1)*(imt_global-2)*(j_f_num) + (j-1)*(imt_global-2) + ix*(i_f_num) + i
             !gindex(n) = (jmt_global - (iy+1)*j_f_num)*(imt_global-num_overlap) + (j-1)*(imt_global-num_overlap) + ix*(i_f_num) + i
             ! lihuimin, 2012.8.7
             ! j_num of iy == 0 is jmt+1,  j_num = j_f_num+(jst_global-1)+1=j_f_num+(3-1)+1=jmt-2+3=jmt+1
             gindex(n) = (jmt_global - iy*j_f_num - (jmt+1))*(imt_global-num_overlap) + (j-1)*(imt_global-num_overlap) + ix*(i_f_num) + i
          endif
!          if (iy == ny_proc-1) then
!             gindex(n) = (j-1)*(imt_global-2) + ix*(i_f_num) + i
!          else
!             gindex(n) = (jmt_global - (iy+1)*j_f_num)*(imt_global-2) + (j-1)*(imt_global-2) + ix*(i_f_num) + i
!          endif
       enddo
    enddo
!    write(*,*)'n=,gsmap',n
    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )

    deallocate(gindex)

  end subroutine ocn_SetGSMap_mct

!*************************************************************************
!ROUTINE: ocn_domain_mct
!         lihuimin 2012.7.16
!
  subroutine ocn_domain_mct(lsize, gsMap_o, dom_o)

    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o
    ! local
    integer :: j_begin,j_end
    integer, pointer :: &
      idata(:)

    real(r8), pointer :: &
      data(:)

    integer (kind(1)) ::   &
      i,j, k, n, iblock, &
      ier



    call mct_gGrid_init( GGrid=dom_o, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_o%data)
    allocate(data(lsize))


    call mct_gsMap_orderedPoints(gsMap_o, mytid, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)




    data(:) = shr_const_spval !-9999.0_R8 !linpf 2012Jul27
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize)
    data(:) = shr_const_spval !0.0_R8  !linpf 2012Jul27
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize)
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize)


      write(6,*)"domain_mct  special, id =  ", mytid,"   imt=",imt,"   jmt=",jmt,"   nx=",nx,"   ny=",ny,"lsize=",lsize

         ! lihuimin, 2012.8.7
    if (iy == 0) then
       j_begin = 0
       j_end = j_num - 1
    else
       j_begin = (iy-1)*j_f_num + (jmt+1)     ! jmt+1 is the j_num of iy==0
       j_end = (iy-1)*j_f_num + (jmt+1) + j_num - 1
    endif

     n = 0
!     do j=iy*j_f_num+j_num-1,iy*j_f_num,-1
      do j=j_begin,j_end
       do i=ix*i_f_num+1,i_num+ix*i_f_num
          n=n+1
          data(n) = mask(i,jmt_global-j)
          if (data(n) > 1.0_r8) data(n) = 1.0_r8
       enddo
     enddo
      write(*,*)'n=,mask',n

!      do j=iy*j_f_num,iy*j_f_num+j_num-1
!       do i=ix*i_f_num+1,i_num+ix*i_f_num
!          n=n+1
!          data(n) = mask(i,jmt_global-j)
!          if (data(n) > 1.0_r8) data(n) = 1.0_r8
!       enddo
!     enddo
     call mct_gGrid_importRattr(dom_o,"mask",data,lsize)
     call mct_gGrid_importRattr(dom_o,"frac",data,lsize)


      n = 0
     !do j=iy*j_f_num,iy*j_f_num+j_num-1
     do j=j_begin,j_end
       do i=ix*i_f_num+1,i_num+ix*i_f_num
          n=n+1
          data(n) = area(i,jmt_global-j)
       enddo
     enddo
     call mct_gGrid_importRattr(dom_o,"area",data,lsize)

     n = 0
     !do j=iy*j_f_num,iy*j_f_num+j_num-1
     do j=j_begin,j_end
       do i=ix*i_f_num+1,i_num+ix*i_f_num
          n=n+1
          data(n) = xc(i,jmt_global-j)
       enddo
     enddo
     call mct_gGrid_importRattr(dom_o,"lon",data,lsize)

     n = 0
     !do j=iy*j_f_num,iy*j_f_num+j_num-1
     do j=j_begin,j_end
       do i=ix*i_f_num+1,i_num+ix*i_f_num
          n=n+1
          data(n) = yc(i,jmt_global-j)
       enddo
     enddo
     call mct_gGrid_importRattr(dom_o,"lat",data,lsize)

    deallocate(data)
    deallocate(idata)

    ! lihuimin, 2012.7.22
    ! allocated in inirun.F90
    deallocate(xc,yc,xv,yv,mask,area)

  end subroutine ocn_domain_mct

!*****************************************************************************
!ROUTINE:   CLEANUP
!           lihuimin, 2012.7.22
!
  subroutine cleanup

    ! allocated in inirun.F90
    deallocate(t_cpl, s_cpl, u_cpl, v_cpl, dhdx, dhdy, Q)
    deallocate(taux, tauy, netsw, lat1, sen, lwup, lwdn, melth, salt, prec, evap, meltw, roff, ifrac, patm, duu10n)

    ! allocated in inirun.F90
    deallocate(h0, u, v, at, hi, itice, alead)

    ! allocated in rdriver.F90
    deallocate(su3,sv3,psa3,tsa3,qar3,uva3,swv3,cld3,sss3,sst3 ,nswv3,dqdt3,chloro3)
    deallocate(seaice3,runoff3)
    deallocate(wspd3,wspdu3,wspdv3,lwv3,rain3,snow3) 

  end subroutine cleanup



 end module ocn_comp_mct
