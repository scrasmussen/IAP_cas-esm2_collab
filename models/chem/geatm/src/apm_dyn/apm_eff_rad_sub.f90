

SUBROUTINE eff_rad &
 & ( rd1d_sulf,rd1d_salt,rd1d_dust &
 &  ,RH,YGF &
 &  ,XM1D,XMDST,MBCOC8 &
 &  ,MNIT,MNH4,MMSA &
 &  ,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV &
 &  ,MBCS,MOCS,MDSTS,MSALTS &
 &  ,r1d_sulf,r1d_salt,r1d_dust ) !&

!
implicit none
include 'apm_parm.inc'
integer,parameter :: NTYP=5
integer :: ibin,N,II,JJ,LL,NI,I3nm
! in
real*8 :: rd1d_sulf(NSO4),rd1d_salt(NSEA),rd1d_dust(NDSTB)
real*8 :: RDRY(NSO4),RDST(NDSTB),RSALT(NSEA)
real*8 :: XM1D(NSO4+NSEA),XMDST(NDSTB),XMSALT(NSEA),MBCOC8(8)
real*8 :: SOAT
real*8 :: MNIT,MNH4,MMSA
real*8 :: MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV
real*8 :: MBCS,MOCS,MDSTS,MSALTS
real*8 :: YGF(99,4)

! out
real*8 :: r1d_sulf(NSO4),r1d_salt(NSEA),r1d_dust(NDSTB)
real*8 :: RSALTWET(NSEA),RSALTWETCM(NSEA),RWET(NSO4),RWETCM(NSO4)
real*8 :: REFOCWET(2),REFBCWET(2),REFBCOC(2)

! local variables
REAL*8  :: DENAER(5)
REAL*8  :: DENWATER,DENSULF,DENWET1,DENWET,DENSALTWET,DENWET2
REAL*8  :: XN(NSO4),XMA(NSO4),XVA(NSO4)
REAL*8  :: XNSALT(NSEA),XVSALT(NSEA),VSALT(NSEA)
REAL*8  :: XNDST(NDSTB),XVDST(NDSTB),VDST(NDSTB),DENDST(NDSTB),RDSTWET(NDSTB)
REAL*8  :: XN1D(NSO4+NSEA)
REAL*8  :: MSP(NTYP),MCORE(NTYP),ZK(NTYP)
REAL*8  :: MSO4

REAL*8  :: RRATIO,RRATIOWET

REAL*8  :: TOTNUMBC(2),TOTNUMOC(2),TOTAREABC(2),TOTAREAOC(2)

REAL*8  :: AREABC(2),AREAOC(2)

REAL*8  :: MSOA

REAL*8  :: MSO4B,MNITB,MNH4B,MMSAB,MSULFT,TOTMP,FSO4B

REAL*8  :: TNBC, TNOC,TNBC1,TNBC2,TNOC1,TNOC2
REAL*8  :: GF2,GF3,GF4,FRACSOA,GFWATERVOL,GFWATER

REAL*8  :: ZDACT(NTYP,3),ZCCN(NTYP,3),ZTN(NTYP),TCCN(3)

REAL*8  :: GFTOT,GFTOT1,GFGAS,GFTOT2,Y1

integer :: IRH

REAL*8  :: RH

integer :: NA

DATA (DENAER(NA),NA=1,5)/1.7,2.2,2.65,1.8,1.8/  ! density (g/cm3)
DENWATER = 1.0  ! density of water (g/cm3)


r_sulf=rd_sulf
r_salt=rd_salt
r_dust=rd_dust

return



! code from apm_phy_mod.f


!
! 1cc = 1mL = 1cm3
!******************************************************************************
! APM AEROSOL Types (N=1,NTYP)
! N=1:  Sulfate or Secondary particles (SO4 plus other species) (density: 1.7 g/cc)
! N=2:  Sea Salt (density: 2.2 g/cc)
! N=3:  DUST  (density: 2.5 g/cc for D<1 um, 2.65 g/cm3 for D>~1 um)
! N=4:  BC  (density: ? 1.0 or 1.8 g/cc)
! N=5:  OC  (density: 1.8 g/cc)
!******************************************************************************

!> just minimal values
  MCORE = 1.d-20
  MSP   = 1.d-21
  ZCCN  = 1.d-20
!!


MSO4 = 0.d0
DO N=1,NSO4
  XN(N)=XN1D(N)
  XMA(N)=XM1D(N) ! total mass in each bin in unit volume
  MSO4 = MSO4 + XM1D(N)   ! total bin sulfate mass
  IF(XMA(N).LT.1.E-40) XMA(N)=1.E-40
  ! total volume in each bin in unit volume
  XVA(N)=XMA(N)*1.d-3/DENSULF   !XVA in cm3/cm3  XMA in kg/m3
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
! Seasalt
!
DO N=1,NSEA
   XMSALT(N) = XM1D(NSO4+N) ! total seasalt mass in each bin in unit volume
   MCORE(2)=MCORE(2) + XMSALT(N)
   IF(XMSALT(N).LT.1.E-40) XMSALT(N)=1.E-40
    XVSALT(N) = XMSALT(N)*1.d-3/DENAER(2) ! XV in cm3/cm3  XM in kg/m3
    XNSALT(N) = XVSALT(N)/(1.E6*VSALT(N)) ! XN in #/cm3, VSALT in m3
ENDDO
ZTN(2) = SUM(XNSALT) ! total seasalt number concentration
!
! Dust
DO N=1,NDSTB
    MCORE(3)=MCORE(3) + XMDST(N)
    XVDST(N)=XMDST(N)*1.d-3/DENDST(N)   !XV in cm3/cm3,DENDST in g/cm3
    XNDST(N) = XVDST(N)/(1.E6*VDST(N))  ! XN in #/cm3, VDST in m3
    RDSTWET(N)=RDST(N)*100.   ! wet size in cm, use RDST for now
ENDDO
ZTN(3) = SUM(XNDST)

! BCOC
MCORE(4) = MBCOC8(1)+MBCOC8(5)+MBCOC8(2)+MBCOC8(6)  !kg/m3 BC
MCORE(5) = MBCOC8(3)+MBCOC8(7)+MBCOC8(4)+MBCOC8(8)  !kg/m3 OC

!******************************************************************************
! Sulfate particle dry size increase due to uptake of NIT, NH4, SOA via
! equilibrium/partition  (ratio same for all sizes for now)
! assume SO4, NIT, NH4, SOA have same density for now

! MSULFLV teated as a part of MSO4 for now but not involved in isoropia calculation
! MSO4 : total sulfate mass after in-cloud oxidation

MSO4B = MSO4-MSULFLV
MSO4B = MAX(MSO4B,1.d-20)
MCORE(1) = MSO4B

DO N=1,NTYP
      IF(MCORE(N).LE.0.) THEN
           WRITE(6,*)"MCORE.LE.0", N, II, JJ,LL
!           WRITE(6,100) (MBCOC8(NI),NI=1,8)
           MCORE(N)=1.d-20
      ENDIF
ENDDO

MSULFT = MSO4B+MBCS+MOCS+MDSTS+MSALTS  !Total SulfateNew  (kg/m3)
MSULFT = MAX(MSULFT,1.d-20)
FSO4B = MSO4B/MSULFT   ! fraction of SO4 in SP
MNITB = MNIT * FSO4B
MNH4B = MNH4 * FSO4B
MMSAB = MMSA * FSO4B   ! shun??


TOTMP = MCORE(5)*2.1+MSULFLV+MBCLV+MOCLV+MDSTLV+MSALTLV 
TOTMP = MAX(TOTMP,1.d-20)
MSOA  = SOAT*MSULFLV/TOTMP ! sv-mv SOA partioned into SP (kg/m3)


  MSP(1) = MSO4+MMSAB+MNITB+MNH4B+MSOA           ! total mass of SP
  MSP(2) = MSALTS*(1.+(MMSA+MNIT+MNH4)/MSULFT) & ! SP mass on sea salt 
 &        +MSALTLV*(1.+SOAT/TOTMP)               ! organic mass on sea salt
  MSP(3) = MDSTS*(1.+(MMSA+MNIT+MNH4)/MSULFT)  & ! SP mass on dust 
 &        +MDSTLV*(1.+SOAT/TOTMP)                ! organic mass on dust
  MSP(4) = MBCS*(1.+(MMSA+MNIT+MNH4)/MSULFT)   & ! SP mass on BC 
 &        +MBCLV*(1.+SOAT/TOTMP)                 ! organic mass on BC  
  MSP(5) = MOCS*(1.+(MMSA+MNIT+MNH4)/MSULFT)   & ! SP mass on OC
 &        +MOCLV*(1.+SOAT/TOTMP)               &
 &        +MCORE(5)*2.*SOAT/TOTMP                ! organic mass on OC


IF(MSO4.LE.0.) THEN
!  write(6,99)II,JJ,LL,MSO4,MSULFLV,MMSA,MNIT,MNH4,MSOA
   GFGAS = 1.
ELSE
! MSULFLV included in MSO4
   GFGAS=(MSP(1)/MSO4)**(1./3.)  ! mass increasing factor for TYP=1
ENDIF                            ! due to uptake of sv species

IRH = INT(RH+0.5)
IRH = MIN0(99,IRH)
IRH = MAX0(1,IRH)
GF2 = YGF(IRH,2)  ! growth factor of SO4+NIT+NH4 component (assuming ammonia bisulfate)
GF3 = YGF(IRH,3)  ! growth factor of SOA component 
FRACSOA = MSOA/MSP(1) ! soa mass fraction in sp(TYP=1)
GFWATERVOL = (1.-FRACSOA)*GF2**3.0 + FRACSOA*GF3**3.0
GFWATER = GFWATERVOL**(1./3.) ! radius hygrocscopic growth factor
!
! Total growth factor
GFTOT = GFGAS*GFWATER
IF(GFTOT<1.D0)THEN
    WRITE(*,*)'GFTOT < 1, set to 1',II,JJ,LL,GFTOT
    GFTOT = 1.D0
ENDIF
GFTOT1 = GFTOT

DO N=1,NSO4
   RWET(N) = RDRY(N)*GFTOT  ! Consider uptake of NIT,NH4,SOA and H20
   RWETCM(N) = RWET(N)*100. ! m -> cm
ENDDO

! Average density of wet particles
DENWET=DENSULF/GFWATERVOL+DENWATER*(1.-1./GFWATERVOL)
DENWET1 = DENWET



 GF4 = YGF(IRH,4)   
 IF(GF4.LT.1.D0)THEN
   WRITE(*,*)'GF4 < 1',II,JJ,LL,GF4
    GF4 = 1.0
 ENDIF
 GFTOT2= GF4

 DO N=1,NSEA
    RSALTWET(N) = RSALT(N)*GF4  ! different from sulfate particle
    RSALTWETCM(N) = RSALTWET(N)*100.
 ENDDO
 Y1 = 1./(GF4**3.)
 ! average density
 DENSALTWET = DENAER(2)*Y1 + DENWATER*(1. - Y1)
 DENWET2 = DENSALTWET

 ZTN(1) = SUM(XN(I3nm:NSO4))


! USE TOTNUMBC,TOTNUMOC (# per kg of BC or OC) to calculate
! TNBC, TNOC (#/cm3) from MBC, MOC (kg/m3)

 TNBC1 = (MBCOC8(1)+MBCOC8(5))*1.E-6*TOTNUMBC(1)   ! #/cm3--FF
 TNBC2 = (MBCOC8(2)+MBCOC8(6))*1.E-6*TOTNUMBC(2)   ! #/cm3--BB
 TNBC  = TNBC1 + TNBC2 
 ZTN(4) = TNBC

 TNOC1 = (MBCOC8(3)+MBCOC8(7))*1.E-6*TOTNUMOC(1)   ! #/cm3--FF
 TNOC2 = (MBCOC8(4)+MBCOC8(8))*1.E-6*TOTNUMOC(2)   ! #/cm3--BB
 TNOC  = TNOC1 + TNOC2
 ZTN(5) = TNOC


 IF(MCORE(4).GT.0.D0)THEN
    RRATIO = (1. + MSP(4)/MCORE(4))**(1./3.)
    RRATIOWET=(1.+MSP(4)*GFWATER**3./MCORE(4))**(1./3.)
 ELSE
    RRATIO=1.D0
    RRATIOWET=1.D0
 ENDIF

 AREABC(1) = 1.d-2*((MBCOC8(1)+MBCOC8(5))*TOTAREABC(1))  & !1.E-2 covert m2/m3 to cm2/cm3
  &  *RRATIOWET*RRATIOWET
 AREABC(2) = 1.d-2*((MBCOC8(2)+MBCOC8(6))*TOTAREABC(2)) &
  &   *RRATIOWET*RRATIOWET

! REFBCOC(2) : Effective radius of BCOC (m)
  REFBCWET(1) = 100.*REFBCOC(1)*RRATIOWET   ! wet size for coag of sulf with BC
  REFBCWET(2) = 100.*REFBCOC(2)*RRATIOWET   ! in cm


IF(MCORE(5).GT.0.D0)THEN
    RRATIO = (1. + MSP(5)/MCORE(5))**(1./3.)
    RRATIOWET=(1.+MOCS*(GFGAS*GFWATER)**3./MCORE(5))**(1./3.)
ELSE
   RRATIO=1.D0
   RRATIOWET=1.D0
ENDIF

 AREAOC(1) = 1.d-2*((MBCOC8(3)+MBCOC8(7))*TOTAREAOC(1)) &
   &  *RRATIOWET*RRATIOWET
 AREAOC(2) = 1.d-2*((MBCOC8(4)+MBCOC8(8))*TOTAREAOC(2)) &
   &  *RRATIOWET*RRATIOWET

 REFOCWET(1) = 100.*REFBCOC(1)*RRATIOWET   ! wet size for coag of sulf with OC
 REFOCWET(2) = 100.*REFBCOC(2)*RRATIOWET   ! in cm

END SUBROUTINE eff_rad
