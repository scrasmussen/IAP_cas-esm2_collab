
subroutine apm_dust_emit &
  & ( myid &
  &  ,lapm &
  &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
  &  ,ip2mem,mem2d &
  &  ,ip3mem,mem3d &
  &  ,FICE,FSNOW,LAND_USE,FSOIL,FVEG,UST &
  &  ,SOILT,SOILRH &
  &  ,temp2m,RHSFC,U10,V10 &
  &  ,KDUSTTOP &
  &  ,DUSTHGTF &
  &  ,Z0,UST0 )

 use apm_varlist
 implicit none
 include 'apm_parm.inc'

 integer :: myid

 real    :: dt

 logical :: lapm

 integer :: ne,nest
 integer :: nx(5),ny(5),nzz
 integer :: sy(5),ey(5),sx(5),ex(5)

 integer :: mem2d
 integer :: mem3d

 integer :: ip2mem(nzz)
 integer :: ip3mem(nzz,nest)

 integer :: KDUSTTOP

 real,dimension(mem2d) :: FICE,FSNOW,LAND_USE,temp2m,RHSFC,U10,V10 &
                         ,FSOIL,FVEG,UST

 real,dimension(mem2d) :: Z0,UST0

 real,dimension(mem3d) :: SOILT,SOILRH

 real,dimension(mem3d) :: DUSTHGTF

 integer :: i03,i03_2,i02,ixy

 integer :: i,j,k

 integer :: iii

!print*,'sub-dust_emit'
!print*,'RHSFC=',RHSFC
!print*,'LAND_USE=',LAND_USE
!print*,'T2m=',temp2m
!print*,'dims=',SX(NE),EX(NE),SY(NE),EY(NE),ne


 I03 = IP3MEM(1,NE)
 I02 = IP2MEM(NE)

 DO K = 1, KDUSTTOP ! =6, in main.f90:'DUST CAN TRANSPORTED TO 1KM TOP'

    I03_2 =  IP3MEM(K,NE)

if(.false.) then
if(k.ne.1) then
do i=sx(ne),ex(ne)
do j=sy(ne),ey(ne)
 ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
 if(i.eq.50.and.j.eq.50) then
 print*,'5050'
 do iii=1,6
   print*,DUSTHGTF(IP3MEM(iii,NE)+ixy)
 enddo
 endif
enddo
enddo
endif

if(k.eq.1) then
 iii=4
! print*,'DUSTHGTF='
! print*,DUSTHGTF(IP3MEM(iii,NE):IP3MEM(iii+1,NE))
endif

endif

!cycle
    CALL apm_dust_emit_sub &
                ( MYID,lapm,DUSTHGTF(I03_2), FICE(I02), FSNOW(I02) &
                 ,LAND_USE(I02), FSOIL(I02) ,FVEG(I02) &
                 ,Z0(I02),UST(I02), UST0(I02), SOILT(I03) &
                 ,temp2m(I02),SOILRH(I03), RHSFC(I02) &
                 ,U10(I02),V10(I02) &
                 ,dust_emit(I03_2) &
                 ,SX(NE),EX(NE),SY(NE),EY(NE),K,NE )

!if(k.eq.1) then
! print*,'dust_emit=',dust_emit(IP3MEM(K,NE):IP3MEM(K+1,NE))
!endif

 ENDDO

!stop 'stop in sub dust emit'

end subroutine apm_dust_emit



SUBROUTINE apm_dust_emit_sub( MYID,lapm,DUSTHGTF,ICE,SNOW,LANDUSE,SOIL,VEG &
                         ,Z0,USTWRF,UST0, SOILT1,T2,SOILRH,RHSFC,U10,V10 &
                         ,DUSTEMISS, SX,EX,SY,EY,K,NE )
INTEGER MYID,SX,EX,SY,EY,NE
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: LANDUSE,ICE,SNOW,U10,V10,DELTZ,SOIL,VEG
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: Z0,UST0,SOILT1,T2,RHSFC,SOILRH,USTWRF
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: DUSTHGTF,DUSTEMISS ! DUSTEMISS is the emissions ug/m2/hr
REAL DT,EMITF,WIND,DTV,Gz1oz0,RB,BB,CONS1,CONS2,UST,SFT
REAL :: VEGF,SOILF ! CORRECTION FACTOR FROM VEG AND SOIL TYPE
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: EMITFACT

logical :: lapm

!print*,'dims in sub=',sx,ex,sy,ey,k,ne


DO J=SY,EY
DO I=SX,EX

  EMITF = 0.0
  DUSTEMISS(I,J)=0.0

  IF(LANDUSE(I,J).EQ.19) EMITF = 1.0
  IF(LANDUSE(I,J).EQ.7.OR.LANDUSE(I,J).EQ.10) EMITF = 0.3
  IF(LANDUSE(I,J).EQ.8.OR.LANDUSE(I,J).EQ.9)  EMITF = 0.1 


  IF(VEG(I,J)<10.) THEN
    VEGF = 1.0
  ELSE IF(VEG(I,J)<30.) THEN
    VEGF = 0.7
  ELSE IF(VEG(I,J)<50) THEN
    VEGF = 0.5
  ELSE
    VEGF = 0.0
  ENDIF

  IF(SOIL(I,J).LE.5 ) THEN
    SOILF = 1.0
  ELSE IF(SOIL(I,J).EQ.6) THEN
    SOILF = 0.7
  ELSE IF(SOIL(I,J).EQ.7) THEN
    SOILF = 0.85
  ELSE
    SOILF = 0.0
  ENDIF


  EMITF = EMITF*VEGF*SOILF

  EMITFACT(I,J)=EMITF

  IF(SOILT1(I,J)>273.15) THEN
    IF(ICE(I,J)<=0.0)THEN
      IF(SNOW(I,J)<=0.0)THEN
        IF(EMITF>0.0)THEN

          WIND=U10(I,J)*U10(I,J)+V10(I,J)*V10(I,J)
          WIND=MAX(WIND,0.1)

!WIND=10

          DTV=T2(I,J)-SOILT1(I,J)
          Gz1oz0=ALOG(10.0/Z0(I,J))
          !RB=9.8*DTV*10.0/(T2(I,J)*WIND)
          RB=9.8*DTV*(10.0-Z0(I,J))/(T2(I,J)*WIND)  !chenhs
          BB=70.0*0.4*0.4*SQRT(ABS(RB)*10.0/Z0(I,J))/(Gz1oz0*Gz1oz0)

          CONS1=0.4*SQRT(WIND)/Gz1oz0

          IF(RB>=0.0)THEN
            CONS2=1.0/(1.0+4.7*RB)
          ELSE
            CONS2=SQRT( 1.0-9.4*RB/(1.0+BB) )
          ENDIF

          UST=CONS1*CONS2

          ! USE WRF OUTPUT UST
          UST = USTWRF(I,J)



if(k.eq.1) then
! print*
! print*,'i&j',i,j
! print*,'UST',UST,UST0(I,J)
! print*,'RHSFC(I,J)=',RHSFC(I,J)
endif


          IF(UST>UST0(I,J).AND.RHSFC(I,J)<40.)THEN

!             C(I,J) = C(I,J) + &
!                      1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)* &
!                      DUSTHGTF(I,J)*(1.0+UST0(I,J)/UST)* &
!                      (1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))* &
!                      (1.0-RHSFC(I,J)/40.)*1.E09 !kg-->ug  !chenhs doctor
              
             ! unit : ug/(m2 s)
             DUSTEMISS(I,J) = 1.0e-5*1.29/9.8*EMITF*UST*UST*UST* &
                              DUSTHGTF(I,J)*(1.0+UST0(I,J)/UST)* &
                              (1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))* &
                              (1.0-RHSFC(I,J)/40.)*1.E09 !kg-->ug ! shun

if(k.eq.1) then
! print*
! print*,'i&j',i,j
! print*,'DUSTEMISS(I,J)=',DUSTEMISS(I,J)
 !print*,'EMITF=',EMITF
 !print*,'DUSTHGTF(I,J)=',DUSTHGTF(I,J)
 !print*,'UST0(I,J)=',UST0(I,J)
 !print*,'UST=',UST
 !print*,'RHSFC(I,J)=',RHSFC(I,J)
endif

  
!  THIS SUBROUTINE IS also TO ALLOCATE THE DUST COMPOSITIONS FROM THE
!  PUBLICATION Feng Yan  et al.,Global Modeling of Nitrate and Ammonium:
!  Interaction of Aerosols and Tropospheric Chemistry.
!  MASS WEIGHT CACO3: 7% MGCO3:5.5%; K2CO3 3.3%; NA2CO3 2.6% SIO2 60%; !  AL2O3
!  14.1%; Fe(III): 2.8 Fe(II) 1.2-04 ! Fe is from Zhao Phd thesis(2006)


            !IF(K==1) then
            !  DUSTEMISS(I,J) = DUSTEMISS(I,J) + &
            !                   1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT*&
            !                   (1.0+UST0(I,J)/UST)* &
            !                   (1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))* &
            !                   (1.0-RHSFC(I,J)/40.)   ! kg/m2/hr
            !ENDIF

          ENDIF

        ENDIF
      ENDIF
    ENDIF
  ENDIF ! SOILT1

  DUSTEMISS(I,J) = AMAX1( DUSTEMISS(I,J), 1.E-20 )

  !if(k.eq.1) print*,'dstem',DUSTEMISS(I,J),i,j

ENDDO ! I
ENDDO ! J

RETURN

END


