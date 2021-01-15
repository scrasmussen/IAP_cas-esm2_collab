
subroutine naqpms_putdust &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,KDUSTTOP &
 &  ,IMONTH2 &
 &  ,SFT &
 &  ,ITT)

use naqpms_varlist
use naqpms_gridinfo
use met_fields
use work_vars, only : DUSTK,DUSTHGTF
implicit none

integer :: myid

real :: dt

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

real    :: SFT(4)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: ITT,IMONTH2,KDUSTTOP


integer :: i,j,k,is

integer :: mem3d

integer :: mem2d

integer :: I03AER,I03,I05,I02,I03_2,I05C9,I05C10
integer :: I05C1,I05C2,I05C3,I05C4,I05C5,I05C6,I05C7,I05C8

real    :: MSIZDIS(isize+1)

 DO K = 1, KDUSTTOP
    I02    = IP2MEM(NE)
    I03    = IP3MEM(K,NE)
    CALL DUSTHEIGHT ( MYID,HEIZ(I03), HGT1(I02),LAND_USE(I02),DUSTK(I03) &
                     ,TOTALDUST(I02),U(I03),V(I03) &
                     ,SX(NE), EX(NE), SY(NE), EY(NE),K,NE )
 ENDDO

 ! TO CALCULATE THE WEIGHTFACT BY DUSTK/TOTALDUST 
 DO K = 1, KDUSTTOP
   I02    = IP2MEM(NE)
   I03    = IP3MEM(K,NE)
   CALL DUSTHGTFACT ( MYID,TOTALDUST(I02),DUSTK(I03),DUSTHGTF(I03) &
                     ,SX(NE),EX(NE), SY(NE), EY(NE),NE )
 ENDDO
! ---

!!!!!!!!!!!!!!!!
! shun modified
 I03    = IP3MEM(1,NE)
 I02    = IP2MEM(NE)
 CALL GETZ0( MYID,LAND_USE(I02),Z0(I02),IMONTH2 &
            ,SX(NE), EX(NE), SY(NE), EY(NE),NE )
 CALL GETUST0( MYID,LAND_USE(I02), FSOIL(I02), FVEG(I02) &
              ,UST0(I02),SOILT(I03),SOILRH(I03), RHSFC(I02),LONGICRS(I02) &
              ,SX(NE),  EX(NE), SY(NE), EY(NE), NE)
!!!!!!!!!!!!!!!!


 DO IS = 1, ISIZE
   I03AER = IP3MEMAER (IS, 2, NE)
   I03    = IP3MEM(1,NE)
   I02    = IP2MEM(NE)
   !
   DO K = 1, KDUSTTOP
    I03_2 =  IP3MEM(K,NE)
    I05    = IP5MEM(K,IS,2,NE)
    I05C1=ip5memc(K,IS,1,NE)
    I05C2=ip5memc(K,IS,2,NE)
    I05C3=ip5memc(K,IS,3,NE)
    I05C4=ip5memc(K,IS,4,NE)
    I05C5=ip5memc(K,IS,5,NE)
    I05C6=ip5memc(K,IS,6,NE)
    I05C7=ip5memc(K,IS,7,NE)
    I05C8=ip5memc(K,IS,8,NE)
    I05C9=ip5memc(K,IS,9,NE)
    I05C10=ip5memc(K,IS,10,NE)

    CALL PUTDUST( MYID,AER(I05) &
                 ,DUSTCOMP(I05C1),DUSTCOMP(I05C2),DUSTCOMP(I05C3) &
                 ,DUSTCOMP(I05C4),DUSTCOMP(I05C5),DUSTCOMP(I05C6) &
                 ,DUSTCOMP(I05C7),DUSTCOMP(I05C8),DUSTCOMP(I05C9) &
                 ,DUSTCOMP(I05C10) &
                 ,DUSTHGTF(I03_2), FICE(I02), FSNOW(I02) &
                 ,LAND_USE(I02), FSOIL(I02) ,FVEG(I02) &
                 ,Z0(I02),UST(I02), UST0(I02), SOILT(I03) &
                 ,T2(I02),SOILRH(I03), RHSFC(I02) &
                 ,U10(I02),V10(I02) &
                 ,EMITFACT(I02),DZ(I03_2),DUSTEMISS(I02) &
                 ,SX(NE),EX(NE),SY(NE),EY(NE), DT,SFT(IS),K,NE ,IS, ITT )

   ENDDO ! K

 ENDDO ! IS





end subroutine naqpms_putdust





