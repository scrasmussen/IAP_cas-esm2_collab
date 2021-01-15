      PROGRAM APM_DRIV

      USE APM_PHYS_MOD, ONLY : APM_PHYS
      USE APM_INIT_MOD, ONLY : IFNUCL,IFAG,APM_INIT,APM_NTRACERS
      USE APM_INIT_MOD, ONLY : XMACID,XMLVSOG,M1ACID,M1LVSOG
      USE APM_INIT_MOD, ONLY : VDRY
      USE APM_INIT_MOD, ONLY : DENSULF
      USE APM_COAG_MOD, ONLY : READCK6DTABLE
      USE APM_NUCL_MOD, ONLY : IONRATE0

      ! Local variables
      INTEGER,PARAMETER :: IIPAR=1,JJPAR=1,LLPAR=1,NTEMPOUT1=1
      INTEGER,PARAMETER :: NGCOND=1,NSO4=40,NSEA=20,NDSTB=15
      INTEGER,PARAMETER :: NCTSO4=0,NCTBCOC=2,NCTDST=1,NCTSEA=1,NBCOCT=8
      INTEGER,PARAMETER :: NBCPIF=1,NBCPIO=1,NOCPIF=1,NOCPIO=1
      INTEGER,PARAMETER :: NBCPOF=1,NBCPOO=1,NOCPOF=1,NOCPOO=1 !luogan

      INTEGER   :: I,J,L,N,SIZENUM,IY,MDAY

      REAL*8  :: PRESS,TK,RHIN,XQ,CACID,PACID,DTAPM,DT
      REAL*8  :: MSO4,MNIT,MNH4,MMSA,SOAT
      REAL*8  :: MBCS, MOCS   ! mass of sulfate attached to BC, OC
      REAL*8  :: MDSTS, MSALTS   ! mass of sulfate attached to dust,sea salt
      REAL*8  :: MSULFT   ! total sulfate
      REAL*8  :: MBCOC8(8)
      REAL*8  :: XM1D(NSO4+NSEA), XN1D(NSO4),TEMPOUT1(NTEMPOUT1)
      REAL*8  :: XMDST(NDSTB)

      REAL*8  :: MASS1, MASS2

      REAL*8  :: VOL, CLVSOG, PLVSOG1
      REAL*8  :: MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV

      REAL*8  :: XOH, XU, XV
      INTEGER :: KYEAR,KMON,KDAY,KHOUR,KMIN,ISITE,JSITE,NSITE
      REAL*8  :: TOP, TOPP
      INTEGER :: KKOUT
      REAL*8  :: XLAT, XLON, YPSURF, YPR, YQ

      REAL*8  :: CSO2,CNH3,XN0
      REAL*8  :: CSOG1,CSOG2,CSOG3,CSOG4,CSOA1,CSOA2,CSOA3,CSOA4
      REAL*8  :: CSOG5,CSOA5
      REAL*8  :: GFTOT1,GFTOT2,DENWET1,DENWET2
      INTEGER :: IACT10, IACT20, IACT30   ! bin index for cloud act
                                          ! diameters corresponding to RDRY
      REAL*8  :: FCLOUD1(NSO4+4)
      INTEGER :: NCOAG1,NCOAG2 ! shun ??

      REAL*8  :: TEMPOUT(IIPAR,JJPAR,LLPAR,NTEMPOUT1)

      integer :: N_APMTRAC,ISURF

      LOGICAL, SAVE    :: FIRST = .TRUE.
      LOGICAL, SAVE    :: FIRST1 = .TRUE.

      character(len=255) :: DATA_DIR

      WRITE(6,*)'    - APM calculation '

      DATA_DIR='../'

      IFNUCL=1
      !initializes number of tracers associated with APM
      call APM_NTRACERS( 0, N_APMTRAC ) ! in 'apm_init_mod.f'
      !initializes APM model
      ! bins setup and emission
      call APM_INIT(DATA_DIR)           ! in 'apm_init_mod.f'

      ! Chemistry timestep [s]
      DT = 60d0
      DTAPM = DT  !s

! TEMPOUT
      IF(FIRST1)THEN
        ! Read coagulation kernel look-up tables
        CALL  READCK6DTABLE ! in 'apm_coag_mod.f'
        TEMPOUT = 0d0
        NPOUTSTEPS = 0
        FIRST1=.FALSE.
      ENDIF

      close(100)
      print*,'xminput.txt :'
      open(100,file='xminput.txt')
      read(100,*)xm1d   !XM1D(NSO4+NSEA)
      close(100)

! shun LV : low volatile organic from oxidation

        PLVSOG1 = 1.d-30
        CLVSOG  = 1.d-30
        MSULFLV = 1.d-30 ! ?
        MBCLV   = 1.d-30 ! ?
        MOCLV   = 1.d-30 ! ?
        MDSTLV  = 1.d-30 ! ?
        MSALTLV = 1.d-30 ! ?

        SOAT= 1.d-30   ! Total SV-&MV-SOA

        MSO4 = 0.d0
        DO N=1,NSO4
          write(*,*)'XM1D',N,XM1D(N)
        ENDDO

        stop 'shun_stop'

        DO N=1,NSO4
          ! luoXM1D(N) = 1.d-30  !kg/m3
          ! XM1D(N) : total mass in each bin in unit volume
          MSO4 = MSO4 + XM1D(N)   ! total bin sulfate mass
          XN1D(N)=XM1D(N)/(DENSULF*VDRY(N))*1.E-9 ! shun: number concentration
          write(*,*)'XN1D',N,XN1D(N)
        ENDDO

        DO N=1,NSEA        ! sea salt
          XM1D(NSO4+N)=1.d-30 !kg/m3
        ENDDO

        DO N=1,NDSTB      ! dust
          XMDST(N) = 1.d-30 !kg/m3
        ENDDO

! BCOC associated with fossil fuel and biomass/biofuel seperated
        DO N= 1, NBCOCT
           MBCOC8(N)=1.d-30 !kg/m3
        ENDDO

        MNIT = 1.d-30 !kg/m3
        MNH4 = 1.d-30 !kg/m3
        MMSA = 1.d-30 !kg/m3

! coated sulfate
        MBCS = 1.d-30   ! SULF on BC kg/m3
        MOCS = 1.d-30   ! SULF on OC kg/m3
        MDSTS = 1.d-30  ! SULF on DUST kg/m3
        MSALTS = 1.d-30 ! SULF on SEA-SALT kg/m3

      DO NPOUTSTEPS = 0, 1440

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

        PRESS = 1.d5   ! P at level center [Pa]
        TK    = 285.d0 ! Temperature [K]
        RHIN  = 90.d0  ! relative humidity %

        !XQ = 1.d-20    !ion-pairs/cm3s ! shun comment out
        ISURF = 1
        YPSURF = 1013.d0
        YPR = 1000.d0
        XLON = 120.d0
        XLAT = 40.d0
        !calculate ionization rate (output ZQ: ion-pairs/cm3s)
        CALL IONRATE0(ISURF, YPSURF, XLON, XLAT, YPR, YQ) ! in 'apm_nucl_mod.f'
        XQ = YQ
        !print *, 'Zifa ion-pairs/cm3 ', NPOUTSTEPS,I,J,L,XQ

        !kg/box to #/cm3
        CACID = 1.d5 ! sulfuric acid molecule number in one cubic centimeter

        !kg/(box*timestep) to #/(cm3*s)
        PACID = 1.d5 ! sulfuric acid molecule number production rate

        ! in 'apm_phys_mod.f'
        CALL APM_PHYS(I,J,L,
     &       NCOAG1,NCOAG2,IACT10,IACT20,IACT30,    ! out
     &       NTEMPOUT1,                             ! in
     &       PRESS,TK,RHIN,XQ,PLVSOG1,              ! in
     &       CACID,PACID,                           ! inout
     &       DTAPM,MMSA,MNIT,MNH4,                  ! in 
     &       MBCS, MOCS,MDSTS, MSALTS,              ! inout
     &       MBCOC8, SOAT,                          ! in
     &       CLVSOG,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV, ! in
     &       GFTOT1,GFTOT2,DENWET1,DENWET2,         ! out
     &       XM1D,XN1D,TEMPOUT1,XMDST,FCLOUD1) ! inout

      !luogan temp output
        DO N=1,NTEMPOUT1
           TEMPOUT(I,J,L,N)=TEMPOUT(I,J,L,N)+ TEMPOUT1(N)
        ENDDO

      ENDDO
      ENDDO
      ENDDO

      write(44,*)NPOUTSTEPS,TEMPOUT(1,1,1,1),TEMPOUT1(1)

      ENDDO
      write(*,*)'program is finished'

      END PROGRAM APM_DRIV
