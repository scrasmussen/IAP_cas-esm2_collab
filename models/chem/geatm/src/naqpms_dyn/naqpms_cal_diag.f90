
subroutine naqpms_cal_diag &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,GC_MOLWT &
 &  ,ktop &
 &  ,kk )

use naqpms_varlist
use naqpms_gridinfo
use met_fields
implicit none

integer :: myid

real :: dt

integer :: kk

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

real    :: GC_MOLWT(igas)

integer :: i,j,k,is
!integer :: k,is

integer :: mem3d

integer :: mem2d
real,dimension(mem2d) :: ktop

integer,dimension(igas) :: PrintGas

integer :: ixy,i02,i03,i04,i0

integer :: ig,iduc,ia,i05,i05c


real,dimension(mem2d) :: trop

! local variable

real,allocatable,dimension(:) :: EXT0,EXTS,DUSTEXT0,EXTSO40,EXTNO30,EXTNH40 &
                                ,EXTBC0,EXTOC0 &
                                ,DSO2,DNO2,DO3

integer :: i04_1,i04_2,i04_3,i04_4,i04_5,i04_6,i04_7 &
          ,i04_8,i04_9,i04_10,i04_11,i04_12,i04_13,i04_14 &
          ,i04_15,i04_16,i04_17,i04_18,i04_19,i04_20,i04_21 &
          ,i04_22,i04_23,i04_24

integer :: I05_1,I05_2,I05_3,I05_4

real :: NH4,SO4,NO3,NA,BC,OC,PM10,PM25
real :: DUST01,DUST02,DUST03,DUST04,SEA01,SEA02,SEA03,SEA04
real :: AODS
real :: temp,PRESS0,RH


ALLOCATE(EXT0(NZZ),EXTS(NZZ),DUSTEXT0(NZZ))
ALLOCATE(EXTSO40(NZZ),EXTNO30(NZZ),EXTNH40(NZZ),EXTBC0(NZZ),EXTOC0(NZZ))


!print*,'kk=',kk
!stop

 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)

      ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1


      i02          =  ip2mem(ne)
      AOD(i02+ixy) =     0.0
      DUSTAOD(I02+IXY) = 0.0
      AODS         =     0.0
      PBLAOD(i02+ixy) =  0.0

      do k=1,nzz-1

        i03      =  ip3mem(k,ne)

        RH  =  rh1(i03+ixy)        ! the RH IN %

        i04_1  = ip4mem(k, 81, ne) ! NH4+(AE)
        i04_2  = ip4mem(k, 83, ne) ! SO42-(AE)
        i04_3  = ip4mem(k, 84, ne) ! HSO4-(AE)
        i04_4  = ip4mem(k, 85, ne) ! NO3-(AE)
        i04_5  = ip4mem(k, 87, ne) ! NA2SO4(S)
        i04_6  = ip4mem(k, 88, ne) ! NANO3(S)
        i04_7  = ip4mem(k, 89, ne) ! (NH4)2SO4(S)
        i04_8  = ip4mem(k, 90, ne) ! NH4NO3(S)
        i04_9  = ip4mem(k, 91, ne) ! NH4CL(S)
        i04_10 = ip4mem(k, 92, ne) ! H2SO4(AQ)
        i04_11 = ip4mem(k, 93, ne) ! NH4HSO4(S)
        i04_12 = ip4mem(k, 94, ne) ! NAHSO4(S)
        i04_13 = ip4mem(k, 95, ne) ! (NH4)3H(SO4)2(S)
        i04_14 = ip4mem(k, 96, ne) ! SOA1
        i04_15 = ip4mem(k, 97, ne) ! SOA2
        i04_16 = ip4mem(k, 98, ne) ! SOA3
        i04_17 = ip4mem(k, 99, ne) ! SOA4
        i04_18 = ip4mem(k,100, ne) ! SOA5
        i04_19 = ip4mem(k,101, ne) ! SOA6
        i04_20 = ip4mem(k, 78, ne) ! POA
        i04_21 = ip4mem(k, 77, ne) ! BC
        i04_22 = ip4mem(k, 76, ne) ! PM10
        i04_23 = ip4mem(k, 75, ne) ! PM25
        i04_24 = ip4mem(k, 80, ne) ! NA

        NH4 =  gas(i04_1+ixy) * 17.0/GC_MOLWT(81)&
             + gas(i04_7+ixy) * 17.0/GC_MOLWT(89)*2.0 &
             + gas(i04_8+ixy) * 17.0/GC_MOLWT(90)&
             + gas(i04_9+ixy) * 17.0/GC_MOLWT(91)&
             + gas(i04_11+ixy)* 17.0/GC_MOLWT(93)&
             + gas(i04_13+ixy)* 17.0/GC_MOLWT(95)*3.0

        SO4 =  gas(i04_2+ixy) * 98.0/GC_MOLWT(83)&
             + gas(i04_3+ixy) * 98.0/GC_MOLWT(84)&
             + gas(i04_5+ixy) * 98.0/GC_MOLWT(87)&
             + gas(i04_7+ixy) * 98.0/GC_MOLWT(89)&
             + gas(i04_10+ixy)* 98.0/GC_MOLWT(92)&
             + gas(i04_11+ixy)* 98.0/GC_MOLWT(93)&
             + gas(i04_12+ixy)* 98.0/GC_MOLWT(94)&
             + gas(i04_13+ixy)* 98.0/GC_MOLWT(95)*2.0

        NO3 =  gas(i04_4+ixy) * 63.0/GC_MOLWT(85)&
             + gas(i04_6+ixy) * 63.0/GC_MOLWT(88)&
             + gas(i04_8+ixy) * 63.0/GC_MOLWT(90)

        NA  =  gas(i04_23+ixy) * 23.0/GC_MOLWT(80)&
             + gas(i04_5+ixy)  * 23.0/GC_MOLWT(87)*2.0&
             + gas(i04_6+ixy)  * 23.0/GC_MOLWT(88)&
             + gas(i04_12+ixy) * 23.0/GC_MOLWT(94)

        BC  =  gas(i04_21+ixy)

        OC  =  gas(i04_14+ixy) + gas(i04_15+ixy) &
             + gas(i04_16+ixy) + gas(i04_17+ixy) &
             + gas(i04_18+ixy) + gas(i04_19+ixy) &
             + gas(i04_20+ixy)

        PM10 = gas(i04_22+ixy)

        PM25 = gas(i04_24+ixy)


        I05_1=IP5MEM(K,1,2,NE)
        I05_2=IP5MEM(K,2,2,NE)
        I05_3=IP5MEM(K,3,2,NE)
        I05_4=IP5MEM(K,4,2,NE)

        DUST01= AER(I05_1+IXY)
        DUST02= AER(I05_2+IXY)
        DUST03= AER(I05_3+IXY)
        DUST04= AER(I05_4+IXY)

        I05_1=IP5MEM(K,1,1,NE)
        I05_2=IP5MEM(K,2,1,NE)
        I05_3=IP5MEM(K,3,1,NE)
        I05_4=IP5MEM(K,4,1,NE)

        SEA01 = AER (I05_1+IXY)
        SEA02 = AER (I05_2+IXY)
        SEA03 = AER (I05_3+IXY)
        SEA04 = AER (I05_4+IXY)

        CALL GETEXT ( MYID, RH, SEA01, SEA02, SEA03, SEA04,DUST01,DUST02,DUST03,DUST04,&
                     NH4, SO4, NO3, NA, BC, OC, PM10, PM25,&
                     DUSTEXT0(K),EXTSO40(K),EXTNO30(K),EXTNH40(K),EXTBC0(K), EXTOC0(K),&
                     EXT0(K),EXTS(K),2, I, J, K)


        EXT(i03+ixy) = EXT0(K)
        EXTASO4(i03+ixy) = EXTSO40(K)
        EXTANO3(i03+ixy) = EXTNO30(K)
        EXTANH4(i03+ixy) = EXTNH40(K)
        EXTBC(i03+ixy)   = EXTBC0(K)
        EXTOC(i03+ixy)   = EXTOC0(K)

        DUSTEXT(I03+IXY) = DUSTEXT0(K)

        CALL GETVIS (MYID, RNW(i03+ixy),CLW(i03+ixy) , EXT(i03+ixy), VISIB(i03+ixy), I, J, K)

        AOD(i02+ixy) = AOD(i02+ixy) + EXT0(K)*dz(i03+ixy)/1.E03
        DUSTAOD(i02+ixy) = DUSTAOD(i02+ixy) + DUSTEXT0(K)*dz(i03+ixy)/1.E03


        AODS = AODS + EXTS(K)*dz(i03+ixy)/1.E03

        SSA(i03+ixy) = EXTS(K)/MAX(EXT0(K), 1.E-20)


        IF(K.LE.NPBL(i02+ixy)) &
        PBLAOD(i02+ixy) = PBLAOD(i02+ixy) + EXT0(K)*dz(i03+ixy)/1.E03

     enddo ! k

!     DEALLOCATE(EXT0,EXTSO40,EXTNO30,EXTNH40,EXTBC0,EXTOC0,EXTS,DUSTEXT0)

  enddo  ! i
  enddo  ! j


!return

!C*********     TO GET THE Dobson (DU) FOR SO2, NO2, O3  *************
!
  DO j = sy(ne),ey(ne)
  DO i = sx(ne),ex(ne)

      ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

      ALLOCATE(DSO2(NZZ),DNO2(NZZ),DO3(NZZ))

      i02          =  ip2mem(ne)
      DUSO2 (i02+ixy) = 0.0
      DUNO2 (i02+ixy) = 0.0
      DUO3  (i02+ixy) = 0.0

      DO k=1,int(ktop(i02+ixy))-1
        i03      =  ip3mem(k,ne)
        temp     =  t(i03+ixy)          ! the temp
        PRESS0   =  Plev(i03+ixy)       ! pressuer in hpa

        i04_1  = ip4mem(k, 6,  ne) ! NO2
        i04_2  = ip4mem(k, 11, ne) ! O3
        i04_3  = ip4mem(k, 18, ne) ! SO2

        CALL GETDU(MYID, gas(i04_1+ixy),DNO2(K), PRESS0, temp )
        CALL GETDU(MYID, gas(i04_2+ixy),DO3 (K), PRESS0, temp )
        CALL GETDU(MYID, gas(i04_3+ixy),DSO2(K), PRESS0, temp )

        DUSO2 (i02+ixy) = DUSO2 (i02+ixy) + DSO2(K)*dz(i03+ixy)/2.69E+20
        DUO3  (i02+ixy) = DUO3  (i02+ixy) + DO3 (K)*dz(i03+ixy)/2.69E+20
        DUNO2 (i02+ixy) = DUNO2 (i02+ixy) + DNO2(K)*dz(i03+ixy)/2.69E+20

      ENDDO  !K

      DEALLOCATE(DSO2,DO3,DNO2)
  ENDDO  !I
  ENDDO  !J 
!CCCCCC END  CCCCCC


end subroutine naqpms_cal_diag


