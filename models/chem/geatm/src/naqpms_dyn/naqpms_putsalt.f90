
subroutine naqpms_putsalt &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,MSIZDIS )

use naqpms_varlist
use naqpms_gridinfo
use met_fields
use work_vars
implicit none

integer :: myid

real :: dt

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: i,j,k,is

integer :: mem3d

integer :: mem2d

integer :: I03AER,I03,I05,I02
integer :: I05C1,I05C2,I05C3,I05C4,I05C5,I05C6,I05C7,I05C8


real    :: MSIZDIS(isize+1)


 DO IS = 1, ISIZE

    I03AER = IP3MEMAER(IS,1,NE)
    I03=IP3MEM(1,NE)
    I05=IP5MEM(1,IS,1,NE)
    I02    = IP2MEM(NE)
    I05C1=ip5memcs(1,IS,1,NE)
    I05C2=ip5memcs(1,IS,2,NE)
    I05C3=ip5memcs(1,IS,3,NE)
    I05C4=ip5memcs(1,IS,4,NE)
    I05C5=ip5memcs(1,IS,5,NE)
    I05C6=ip5memcs(1,IS,6,NE)
    I05C7=ip5memcs(1,IS,7,NE)
    I05C8=ip5memcs(1,IS,8,NE)

   CALL PUTSALT( MYID, AER(I05),SEACOMP(I05C1),SEACOMP(I05C2),SEACOMP(I05C3) &
                ,SEACOMP(I05C4),SEACOMP(I05C5),SEACOMP(I05C6),SEACOMP(I05C7) &
                ,SEACOMP(I05C8),MSIZDIS(IS),MSIZDIS(IS+1),LAND_USE(I02) &
                ,FICE(I02),U10(I02),V10(I02),RHSFC(I02) &
                ,DZ(I03),DX(I03) &
                ,SEAEMISS(I02) &
                ,SX(NE),EX(NE), SY(NE), EY(NE), NE, DT )
 ENDDO ! IS





end subroutine naqpms_putsalt





