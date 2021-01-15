!----------------------------------------------------------------------------
!    This program bases on the original GEATM.  It refactorys the original
!    GEATM with the new modulization and stratification. The new top layer
!    program is added to embed GEATM wihtin CESM. The output method is 
!    is replaced with the new method to make sure the efficiency.
!     
!    Juanxiong He, 2013-06-30
!----------------------------------------------------------------------------

module geatm_vartype
                                     
 parameter(ndomain=5)   
  
 character*11 :: ename='naqpms_home'
 character    :: naqpms_dir*500
 integer      :: getcwd,iostatus
 integer      :: ratio
 integer      :: NEMIT,NLAY_EM
 
 integer :: igas,iopen,iprecise,ikosaline,imasskeep,igasCBM,iprocess
 integer :: NSOA,ifseacom,ifdustcom,nseacom,ndustcom
 integer :: isize,iaer

 integer :: itotalspe,isrnum,nv4poltr
 !----------------------------------------------------
 ! domain infomation, domain decomposition and data read
 !----------------------------------------------------
  integer :: ne,nzz,ntt,nest

  integer, dimension(ndomain) :: nx, ny ,nz, ntbeg,ntend,dtout,dtmet  
  integer, dimension(ndomain) :: nxlo, nylo

  integer, dimension(ndomain) :: comm2d, nbrleft, nbrright, nbrtop, nbrbottom, &
                                 sx, ex, sy, ey, stride,                       &
                                 irec, irecwet,irecdry,irecg,irecnox,irecairc, irec_as,&
                                 irec80,irecglobal,irecdust, irecsea,irec60,&
                                 irecMOZART,irecHgA,irecHgN
                                 
  integer, dimension(4,ndomain) :: bdysx, bdyex, bdysy, bdyey

  integer, dimension(4000) :: IISCPU,IIRCPU,IsLocX,IsLocY &
                              ,IsSNWE,IsNest,IrLocX,IrLocY

  integer, dimension(ndomain) :: mem_per_block

  integer :: mem2d,mem3d,mem4d,mem5d,mem2dgas,mem3daer,mem_emt2d
 
  integer, dimension(2,ndomain) :: dims, coords
  logical, dimension(2,ndomain) :: periods

  integer, allocatable, dimension(:,:) :: sxc, exc, syc, eyc ! allocated in mpi_domain_split

  integer,dimension(:,:), allocatable :: procs ! for cam, added by Juanxiong He
  integer, dimension(ndomain) :: csx, cex, csy, cey, cnx, cny ! for cam, added by Juanxiong He
  integer myid, newid, numprocs, myid_x, myid_y, local_com ! for cam, added by Juanxiong He
  integer myid2d,local_com2d

!-----------------------------------------------------------------------
  character*4 :: cdnum
  character*2 :: cdnum2

  integer :: iyear1,imonth1,idate1,ihour1,iminute1
  integer :: iyear2,imonth2,iday2, ihour2,iminute2
  integer :: iitime

!===================
! lglbrun  shun@20161206
  INTEGER MEMARK
  INTEGER,ALLOCATABLE,DIMENSION(:) :: STAMARK ! allocated in mpi_domain_split

  INTEGER,DIMENSION(ndomain,720) :: IPOLARMRK  ! for 1x1 degree is 720
  INTEGER  IPOLARNUM
  
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: ATESTR4W,ATESTS4W ! allocated in geatm_init_mct
  REAL,ALLOCATABLE,DIMENSION(:)   :: ATESTR4W0,ATESTS4W0 ! allocated in geatm_init_mct

!===================

  real,dimension(ndomain) :: time
  real :: time1hr12dt

! FOR DUST AND SEA SALT
  real, allocatable, dimension(:) :: MSIZDIS,MSIZDID ! sea salt and dust PARTICLES SIZE,allocated in geatm_init_mct
  REAL :: SFT(4)   !! SFT IS THE MASS FRACTION OF DIFFERENT DUST PARTICLES SIZE IN TOTAL 
                   !! DUST FOUR SIZE : 0.1-1.0um, 1-2.0um, 2.0-5.0um, 5.0-10.0um           

  integer :: nstep
  real :: nstep0(2000,2000)

  integer, allocatable, dimension(:) :: IPIG ! allocated in geatm_init_mct

  real, allocatable, dimension(:,:) :: atestR,atestS ! allocated in geatm_init_mct
  real, allocatable, dimension(:)   :: atestR0,atestS0 ! allocated in geatm_init_mct

  real, allocatable, dimension(:,:) :: gravel ! allocated in geatm_init_mct

  REAL, PARAMETER ::  STDATMPA = 101325.0 ! standard atmosphere  [ Pa ]
  REAL, PARAMETER ::  PA2ATM = 1.0 / STDATMPA

  logical :: lprocess
  integer,dimension(102) :: PrintTermGas  ! process by lijie
  integer,dimension(102) :: PrintGas
  integer,dimension(102) :: InitPrintGas
  real,   dimension(102) :: dryvelinit,GC_MOLWT
  character*40, dimension(102) :: GC_Unit,GC_NAME

  character*4 :: cyear
  character*2 :: cmonth,cday

  integer,dimension(ndomain) :: irecSM,ifsm

  integer :: idmSet,iSrcDefined,ismMax,iHgtLMax,ismMaxHere

  integer :: ifsmt

 integer,allocatable,dimension(:)    :: igMark  ! to mark gas   (idmSet), allocated in geatm_init_mct
 integer,allocatable,dimension(:)    :: iaMarkAer ! to mark aer type (idmSet), allocated in geatm_init_mct
 integer,allocatable,dimension(:)    :: iaMarkSiz ! to mark aer size (idmSet), allocated in geatm_init_mct

!=======================================================================
!> shun : apm variables

logical   :: lapm,laer,ignore_emtyp
character :: ctdway*4

character*2  :: ccyear*4,ccmonth,ccday,cchour,ccminute,ccsecond
character*19 :: timestr
integer      :: cur_year,cur_month,cur_day,cur_hour,cur_minute,cur_second 

logical :: lprint_time=.true.

logical :: lnaqpms_pso4
logical :: laqchem

real,parameter :: wgt_h2so4=98,wgt_nh3=17,wgt_hno3=63,wgt_hcl=36.5
real,parameter :: wgt_so4=96,wgt_nh4=18,wgt_no3=62,wgt_na=23,wgt_cl=35.5

! shun : emit sensitivity
logical   :: lnaqpms_ems
integer   :: iemsfunit,idxems
character :: emsflag*3
real      :: emsallfrc

real,allocatable,dimension(:,:)  :: emsfrc ! allocated in geatm_init_mct
integer,allocatable,dimension(:) :: ipig_ems ! allocated in geatm_init_mct

real,allocatable,dimension(:) :: ratioem1,ratioem2,ratioem3,ratioem4,ratioem5,ratioem6 ! allocated in geatm_init_mct

DATA GC_MOLWT/98.,  63.,  36.,  17.,  30.,  46.,  &  ! 1-6
              62.,  108., 47.,  79.,  48.,  16.,  &  ! 7-12
              16.,  17.,  33.,  34.,  28.,  64.,  &  ! 13-18
              16.,  30.,  47.,  28.,  30.,  46.2, &  ! 19-24
              45.,  48.,  61.,  44.,  46.,  1.,   &  ! 25-30
              72.,  121., 14.,  75.3,  72.,  28., &  ! 31-36
              27.,  27.,  92.,  106., 108., 109., &  ! 37-42
              139., 84.,  14.,  1.,   1.,   1.,   &  ! 43-48
              14.,  1.,   14.,   68.,  68.,  68., &  ! 49-54
              68.,  68.,  62.,   1.,   1.,   1.,  &  ! 55-60
              80.,  93.,  79.,  96.,  111., 125., &  ! 61-66
              132., 160., 136., 136., 168., 168., &  ! 67-72
              130., 130.,  1.,   1.,   12., 220., &  ! 73-78
              1.0,  23.0, 18.0, 35.5, 96.0, 97.0, &  ! 79-84
              62.0, 58.5, 142.0,85.0, 132.0,80.0, &  ! 85-90
              53.5, 98.0, 115.0,120.0,247.0,136., &  ! 91-96
              136., 168., 168., 130., 130., 18./     ! 97-102
              
DATA GC_Unit/'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb',                         &
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3', 'ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3', 'ugm-3'/
              
DATA GC_NAME / 'H2SO4',   'HNO3'  , 'HCL'   , 'NH3' ,   'NO'   ,  'NO2'  ,&  !1-6
               'NO3'  ,   'N2O5'  , 'HONO'  , 'HNO4',   'O3'   ,  'O1D'  ,&  !7-12
               'O3P'  ,   'OH'    , 'HO2'   , 'H2O2',   'CO'   ,  'SO2'  ,&  !13-18
               'CH4'  ,   'C2H6'  , 'CH3O2' , 'ETHP',   'HCHO' ,  'CH3OH',&  !19-24
               'ANOL' ,   'CH3OOH', 'ETHOOH', 'ALD2',   'HCOOH',  'RCOOH',&  !25-30
               'C2O3' ,   'PAN'   , 'PAR'   , 'AONE',   'MGLY' ,  'ETH'  ,&  !31-36
               'OLET' ,   'OLEI'  , 'TOL'   , 'XYL' ,   'CRES' ,  'TO2'  ,&  !37-42
               'CRO'  ,   'OPEN'  , 'ONIT'  , 'ROOH',   'RO2'  ,  'ANO2' ,&  !43-48
               'NAP'  ,   'XO2'   , 'XPAR'  , 'ISOP',   'ISOPRD', 'ISOPP',&  !49-54
               'ISOPN',   'ISOPO2', 'DMS'   , 'MSA' ,   'DMSO' ,  'DMSO2',&  !55-60
        'CH3SO2H', 'CH3SCH2OO','CH3SO2','CH3SO3','CH3SO2OO','CH3SO2CH2OO',&  !61-66
             'SULFHOX',  'TERP'   , 'SV1'   , 'SV2' ,   'SV3'  ,  'SV4'  ,&  !67-72
             'SV5'    ,  'SV6'    ,                                       &  !73-74
             'PM25'   , 'PM10',   'BC'  ,   'OC'   , 'HAQ',&                 !75-79
             'NAAQ',  'NH4AQ', 'CLAQ','SO4AQ','HSO4AQ','NO3AQ',&             !80-85
             'NACLS',  'NA2SO4','NANO3','NH42SO4','NH4NO3','NH4CL',&         !86-91
             'H2SO4S','NH4HSO4S','NAHSO4S','NH44HSO42', 'SOA1',&            !92-96
             'SOA2'  ,   'SOA3'   ,   'SOA4'   ,  'SOA5'  ,  'SOA6', 'AH2O'/ !97-102

!!!! for gas (CBM-Z) and ISORROPIA AEROSOLS--------------------------------+
!                 1        2      3      4      5      6       !
!                 H2SO4(g) HNO3   HCL    NH3    NO     NO2
!                 7        8      9      10     11     12      +
!                 NO3     N2O5   HONO   HNO4   O3     O1D
!                 13       14     15     16     17     18      +
!                 O3P     OH     HO2    H2O2   CO     SO2
!                 19       20     21     22     23     24      +
!                 CH4     C2H6    CH3O2  ETHP   HCHO  CH3OH
!                 25       26     27     28     29     30
!                 ANOL    CH3OOH  ETHOOH ALD2   HCOOH RCOOH
!                 31       32   , 33   , 34     35     36
!                 C2O3    PAN     PAR   AONE   MGLY   ETH
!                 37       38     39     40     41     42
!                 OLET  , OLEI ,  TOL   XYL    CRES   TO2    
!                 43       44     45     46     47     48
!                 CRO     OPEN    ONIT  ROOH   RO2    ANO2
!                 49       50     51    52      53     54
!                 NAP     XO2     XPAR  ISOP   ISOPRD ISOPP
!                 55       56     57    58      59     60 
!                 ISOPN   ISOPO2  DMS  MSA     DMSO   DMSO2
!                 61       62         63     64     65     66
!                 CH3SO2H CH3SCH2OO CH3SO2  CH3SO3 CH3SO2OO CH3SO2CH2OO
!                  67     68      69     70    71    72
!                 SULFH   TERP    SV1    SV2   SV3   SV4
!                  73     74 
!                  SV5    SV6
!                  75     76     77     78   79      80
!                   PM25  PM10   BC     OC   H+(AQ) NA+(AQ)
!                   81         82       83   84          85
!                 NH4+(AQ) CL-(AQ) SO4--(AQ) HSO4-(AQ)  NO3-(AQ)
!                  86      87        88         89        90
!                 NACL(S)  NA2SO4(S) NANO3(S)  NH42SO4(S) NH4NO3(S)
!                  91       92        93          94         95
!                  NH4CL(S) H2SO4(AQ) NH4HSO4(S)  NAHSO4(S) (NH4)4H(SO4)2(S)
!                  96      97        98      99    100   101, 102
!                  SOA1    SOA2     SOA3    SOA4   SOA5  SOA6 AH2O'
DATA dryvelinit /0.02 ,  0.02 , 0.01, 0.65, 0.1,  0.5 ,&
                0.2,    0.3 ,  0.05, 0.05, 0.2,  0.001,&
                0.001 , 0.01 , 0.1,  0.1,  0.1,  0.48 ,&
                0.1 ,   0.1 ,  0.02, 0.02, 0.1,  0.1  ,&
                0.1,    0.1 ,  0.1,  0.1,  0.1,  0.1  ,&
                0.02,   0.1 ,  0.1,  0.1,  0.1,  0.1  ,&
                0.1,    0.1 ,  0.1,  0.1,  0.1,  0.1  ,&
                0.1,    0.1 ,  0.1,  0.1,  0.02, 0.02 ,&
                0.1,    0.1 ,  0.1,  0.1,  0.1,  0.1  ,&
                0.1,    0.1 ,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.02,  0.03, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01,0.01 /                 

DATA PrintGas/1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1   , &
              1    , 1    , 1    , 1    , 1    , 1/

    DATA InitPrintGas/1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   /

!++++++++++++++++++++++for  process by lijie+++++++++++++++++++++++
DATA PrintTermGas/0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 1    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0/  
!++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
  data periods/10*.false./
!-----------------------------------------------------

 DATA SFT/0.03,0.10,0.75,0.12/

  contains
  
  subroutine initial_geatm_var
  
  end subroutine initial_geatm_var

  subroutine final_geatm_var
  deallocate(procs,ATESTR4W,ATESTS4W,ATESTR4W0,ATESTS4W0)
  deallocate(ratioem1,ratioem2,ratioem3,ratioem4,ratioem5,ratioem6)
  end subroutine final_geatm_var
end module geatm_vartype  
