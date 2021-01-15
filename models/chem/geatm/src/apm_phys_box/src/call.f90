
implicit none
include 'apm_parm.inc'
real*8  :: DT,XLON,XLAT
real*8  :: PRESS,TK,RHIN
integer :: ISURF
real*8  :: YPSURF,YPR
real*8  :: CACID,PACID
real*8  :: SOAT,MNIT,MNH4,MMSA
real*8  :: XN1D(NSO4+NSEA),XMDST(NDSTB),MBCOC8(8)
real*8  :: MBCS,MOCS,MDSTS,MSALTS
real*8  :: MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV
real*8  :: CLVSOG, PLVSOG1
real*8  :: XM1D(NSO4+NSEA)

integer :: N



DT = 60d0

XLON = 120.d0
XLAT = 40.d0

PLVSOG1 = 1.d-30
CLVSOG  = 1.d-30

MSULFLV = 1.d-30
MBCLV   = 1.d-30
MOCLV   = 1.d-30
MDSTLV  = 1.d-30
MSALTLV = 1.d-30
SOAT    = 1.d-30

MNIT = 1.d-30
MNH4 = 1.d-30
MMSA = 1.d-30

MBCS   = 1.d-30
MOCS   = 1.d-30
MDSTS  = 1.d-30
MSALTS = 1.d-30

PRESS  = 1.d5
TK     = 285.d0
RHIN   = 90.d0
ISURF  = 1
YPSURF = 1013.d0
YPR    = 1000.d0

CACID = 1.d5
PACID = 1.d5


close(100)
open(100,file='xminput.txt')
read(100,*)XM1D   !XM1D(NSO4+NSEA)
close(100)


do N=1,NSEA        ! sea salt
   XM1D(NSO4+N)=1.d-30 !kg/m3
enddo

do N=1,NDSTB      ! dust
   XMDST(N) = 1.d-30 !kg/m3
enddo

do N= 1, NBCOCT
   MBCOC8(N)=1.d-30 !kg/m3
enddo


call naqpms_apm_physics &
     &  ( DT,XLON,XLAT  & ! in
     &   ,PRESS,TK,RHIN,ISURF,YPSURF,YPR & ! in
     &   ,CACID & ! inout
     &   ,PACID & ! in
     &   ,SOAT,MNIT,MNH4,MMSA & ! in
     &   ,XM1D,XN1D & ! inout
     &   ,XMDST,MBCOC8 & ! inout
     &   ,MBCS,MOCS,MDSTS,MSALTS & ! inout
     &   ,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV ) ! inout

end
