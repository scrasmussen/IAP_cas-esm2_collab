          SUBROUTINE FEEVOLUTION(MYID, FEIIIC, FEIII, FEII, SO2, HNO3,&
                         RK_HETSO2, RK_HETHNO3,FRCL,FRCM,FRCH,&
                         SWDOWN,RH1,DT,SX,EX,SY,EY,K,IS)
! ********  THIS ROUTINE IS TO CALCULATE THE FE(III) TO FE(II) IN DUST *********
! *************  PARTICLES BY LUO(2005) JGR, DOI:10.1029/2005JD006059  *********
! ****************    FAN ET AL (2006,GRL, DOI:10.1029/2005GL024852    *********
          INTEGER :: MYID, SX,EX,SY,EY
          REAL :: FAVG ! AVERAGE GLOBAL MEAN SHORTWAVE FLUX
          REAL :: TOBS ! ESTIMATED DAECAY LIFETIME
          REAL :: KSR  ! DECAY RATE BU SOLAR 
          REAL :: KFC  ! THE RATE COEFFICIENT FROM FRESH DUST TO COATED DUST BY HETEOROGENEOUS REACTIONS
          REAL :: KN, KS ! TEMPORARY VARAIABLE
          REAL :: RFE  ! RATE OF FE DISSOLUTION GRAMS OF FE IN DISSOLVED PER GRAM OF FE IN FE2O3 PER SECOND
          REAL :: RD   ! GRAMS OF FE DISSOLVED PER GRAM OF FE2O3
          REAL :: A    ! THE SPECIIFIC SURFACE AREA OF FE2O3
          REAL :: M,N,W
          REAL :: KCLD ! DECAY RATE BY CLOUD
          REAL :: FRCAVG ! GLOBAL AVERAGE CLOUD FRACTION
          REAL :: FRC   
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: FEIII,FEIIIC,FEII,SO2,HNO3 ! SO2 AND HNO3 IN PPBV
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: RK_HETSO2,RK_HETHNO3 ! HETEROGENOUS REACTIONS RATES OF SO2 AND HNO3 ON DUST
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: SWDOWN,RH1 ! in %
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: FRCL,FRCM,FRCH !LOW, MIDDLE AND HIGH CLOUD FRACTION
         

          DATA FAVG / 535.25 /  ! W/M2
          DATA TOBS / 2.59E07 / ! SECONDS          
          DATA RD   / 1.E-10 / ! MOL/M2/S
          DATA A    / 100./  ! M2/G
          DATA FRCAVG /0.05/ ! NO UNIT

          DO I = SX-1,EX+1
          DO J = SY-1,EY+1

!            IF(I==44.AND.J==48.AND.K==1.AND.IS==1) PRINT*, SWDOWN(I,J),RH1(I,J),FEIII(I,J),FEII(I,J),'before'
! --- TO ESTIMATE THE IMPACT OF SOLAR FROM LUO(JGR, 2005) 
           KSR = SWDOWN(I,J) / FAVG / TOBS
           FEII(I,J)  = FEII(I,J) + FEIII(I,J) * ( 1. -  EXP(-DT*KSR) )
           FEIII(I,J) = FEIII(I,J)* EXP( -DT * KSR )  
!            IF(I==44.AND.J==48.AND.K==1.AND.IS==1) PRINT*, SWDOWN(I,J),RH1(I,J),FEIII(I,J),FEII(I,J),'solar'

! --- TO ESTIMATE THE IMPACT OF  CLOUD  BY LUO  (2005,JGR)
!           IF ( K<=10.and.K>=4) FRC = FRCL(I,J) ! juanxiong he
           IF ( K<=10) FRC = FRCL(I,J)
           IF ( K>10.AND.K<=13.) FRC = FRCM(I,J)
           IF ( K>13) FRC = FRCH(I,J)
           KCLD =  FRC/FRCAVG/TOBS
           FEII(I,J)  = FEII(I,J) + FEIII(I,J) * ( 1. -  EXP(-DT*KCLD) )
           FEIII(I,J) = FEIII(I,J)* EXP( -DT * KCLD )
  

! --- TO ESTIMATE THE IMPACT OF HETEOROGENEOUS REACTIONS BY FAN ! (2006,GRL)
           RH = RH1(I,J)
           IF(RH < 25. ) KN = 0
           IF(RH <= 35..AND.RH >= 25.) KN = 5.E-6 * ( RH - 35.) / (35.-25.)
           IF(RH > 35. ) KN = 5.E-6
           IF(RH < 50. ) KS = 0.
           IF(RH > 60. ) KS = 3.E-6
           IF(RH>=50..AND.RH<=.60) KS = 3.E-6 * (RH -50.)/(60.-50.)

           KN = RK_HETHNO3(I,J)! 5.E-6
           KS = RK_HETSO2(I,J) ! 1.E-3

           KFC = KN * HNO3(I,J) + KS * SO2(I,J)
           FEIIIC(I,J) = FEIIIC(I,J) + FEIII(I,J) * KFC
           FEIIIC(I,J) = AMAX1(FEIIIC(I,J),1.E-20) 

           N = 2. ! MOLES FE/MOLE FE2O3
           M = 55.8 ! 55.8 G/MOLE
           W = 0.7 ! THE MASS FRACTION OF FE IN FE2O3
           RFE = RD * A * N * M / W ! g(feII)/s/g(feIII)

           FEII(I,J) = FeII(I,J) +    RFE * FEIIIC(I,J) * DT
           FEIIIC(I,J) = FEIIIC(I,J) -  RFE * FEIIIC(I,J) * DT
           FEIII(I,J) = FEIII(I,J) -  RFE * FEIIIC(I,J) * DT
 

           FEIII(I,J) = AMAX1 ( FEIII(I,J), 1.E-20)              
           FEII(I,J) = AMAX1 ( FEIII(I,J)*0.005,FEII(I,J) )

!           IF(I==44.AND.J==48.AND.K==1.AND.IS==1) PRINT*,SWDOWN(I,J),RH1(I,J),FEIII(I,J),FEII(I,J),'chemistry'
! --- TO ESTIMATE THE IMPACT OF SO4  HETEROGENEOUS REACTIONS(FE(III)+H2O+SO4-->FE(II)+SO4) BY LUO  (2005,JGR)
           



          ENDDO ! J
          ENDDO ! I

          RETURN
          ENDSUBROUTINE
