
subroutine naqpms_gra_dep &
 &  (myid &
 &  ,dt,numTGRV &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom  &
 &  ,gravel) ! &
! &  ,DUSTGRAV,DUSTGRAVSO4,DUSTGRAVNO3,DUSTGRAVFEII,DUSTGRAVFEIII )


use naqpms_varlist
use naqpms_gridinfo, only : dz
implicit none

integer :: myid

real    :: dt
integer :: numTGRV


integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: i,j,k,is,IT


integer :: mem3d


integer ::  mem2d

!real,dimension(mem2d) :: DUSTGRAV,DUSTGRAVSO4,DUSTGRAVNO3 &
!                        ,DUSTGRAVFEII,DUSTGRAVFEIII

integer :: ixy,i02,i03,I03_1,I03P1,iapm,iapm_1,iapm_p1,i03_kk

integer :: i05,i05_1,i05p1,ia,iduc


integer :: kk
real,dimension(mem3d) :: wk,wks


real :: rgf_tmp,rgf_tmp_p1
real :: diam,rhop

real :: gravel(isize,iaer)
real :: GVEL

integer :: ITGRV


DO ITGRV=1,numTGRV

  DO K = 1 , NZZ

    IF(K.GT.1)THEN
     I03_1=IP3MEM(K-1,NE)
    ELSE
     I03_1=IP3MEM(1,NE)
    ENDIF

     I03=IP3MEM(K,NE)

    IF(K.LT.NZZ)THEN
     I03P1=IP3MEM(K+1,NE)
    ELSE
     I03P1=IP3MEM(K,NE)
    ENDIF

  ! FOR DUST AND SEA SALT AEROSOLS
    DO IA=1,IAER
    DO IS=1,ISIZE
     I02 = IP2MEM(NE)

     IF(K.GT.1)THEN
       I05_1=IP5MEM(K-1,IS,IA,NE)
     ELSE
       I05_1=IP5MEM(1,IS,IA,NE)
     ENDIF

     I05=IP5MEM(K,IS,IA,NE)

     IF(K.LT.NZZ)THEN
       I05P1 = IP5MEM(K+1,IS,IA,NE)
     ELSE
       I05P1 = IP5MEM(K,IS,IA,NE)
     ENDIF

  ! TO ADD GRAVITY SETTLING ON JANUARY 09,2005

     GVEL=GRAVEL(IS,IA)
     CALL GRA_DEP_AER(MYID,AER(I05),AER(I05P1),GVEL,K,NZZ,&
                 DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),DT/FLOAT(numTGRV))
     IF(IA==2) THEN
       CALL DUSTGRADEP(MYID, AER(I05), AER(I05P1),DUSTGRAV(I02),GVEL,K,NZZ,&
                 DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),DT/FLOAT(numTGRV))
     ENDIF

     IF(IA==1)  THEN ! FOR SEA SALT 
      do iduc = 1, nseacom
       IF(K.GT.1)THEN
        I05_1=IP5MEMcs(K-1,IS,iduc,NE)
       ELSE
        I05_1=IP5MEMcs(1,IS,iduc,NE)
       ENDIF

        I05=IP5MEMcs(K,IS,iduc,NE)

       IF(K.LT.NZZ)THEN
        I05P1 = IP5MEMcs(K+1,IS,iduc,NE)
       ELSE
        I05P1 = IP5MEMcs(K,IS,iduc,NE)
       ENDIF

       GVEL=GRAVEL(IS,IA)
       CALL GRA_DEP_AER(MYID,SEACOMP(I05),SEACOMP(I05P1),GVEL,K,NZZ,&
             DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),DT/FLOAT(numTGRV))

      enddo ! iduc

     ELSE IF (IA==2) THEN ! DUST

      do iduc = 1, ndustcom
       IF(K.GT.1)THEN
        I05_1=IP5MEMc(K-1,IS,iduc,NE)
       ELSE
        I05_1=IP5MEMc(1,IS,iduc,NE)
       ENDIF

       I05=IP5MEMc(K,IS,iduc,NE)

       IF(K.LT.NZZ)THEN
        I05P1 = IP5MEMc(K+1,IS,iduc,NE)
       ELSE
        I05P1 = IP5MEMc(K,IS,iduc,NE)
       ENDIF

       GVEL=GRAVEL(IS,IA)

       IF(iduc==7) THEN
        CALL DUSTGRADEP(MYID, DUSTCOMP(I05), DUSTCOMP(I05P1),DUSTGRAVSO4(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),DT/FLOAT(numTGRV))
       ELSE IF (iduc==8) THEN
        CALL DUSTGRADEP(MYID, DUSTCOMP(I05), DUSTCOMP(I05P1),DUSTGRAVNO3(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),DT/FLOAT(numTGRV))
       ELSE IF(iduc==9) THEN
        CALL DUSTGRADEP(MYID, DUSTCOMP(I05),DUSTCOMP(I05P1),DUSTGRAVFEII(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),DT/FLOAT(numTGRV))
       ELSE IF(iduc==10) THEN
        CALL DUSTGRADEP(MYID, DUSTCOMP(I05), DUSTCOMP(I05P1),DUSTGRAVFEIII(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),DT/FLOAT(numTGRV))
       ENDIF

       CALL GRA_DEP_AER(MYID,DUSTCOMP(I05),DUSTCOMP(I05P1),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),DT/FLOAT(numTGRV))

       enddo ! iduc


    ENDIF


  ENDDO ! IS
  ENDDO ! IA

  ENDDO !K

  ENDDO ! ITGRV

end subroutine naqpms_gra_dep

