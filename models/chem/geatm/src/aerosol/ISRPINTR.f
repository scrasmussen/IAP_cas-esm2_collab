C=======================================================================
C
C ***   SUBROUTINE ISRPINTR
C ***   THIS SUBROUTINE IS A SAMPLE INTERFACE BETWEEN AN AIRSHED MODEL
C ***   AND ISORROPIA FOR INORGANIC GAS/AEROSOL PATIONING
C
C=======================================================================
C
      SUBROUTINE ISRPINTR (MYID,QI,RHI,TEMPI,CP,CGAS1,CGAS2,CGAS3,
     &                    CGAS4,AWATER,I,J,K)
C
      INTEGER  I, J, K,MYID
      DOUBLE PRECISION QI(5)
C                   ! THE INPUT PRECURSORS GAS+AEROSOL IN UG/M3
C     1: TOTAL SODIUM(GAS+AEROSOL) EXPRESSED AS EQUIVALENT  Na
C     2: TOTAL SULFATE(GAS+AEROSOL)EXPRESSED AS EQUIVALENT  H2SO4
C     3: TOTAL AMMONIUM(GAS+AEROSOL)EXPRESSED AS EQUIVALENT NH3 
C     4: TOTAL NITRATE(GAS+AEROSOL)EXPRESSED AS EQUIVALENT  HNO3 
C     5: TOTAL CHLORIDE(GAS+AEROSOL)EXPRESSED AS EQUIVALENT HCL
C      
      DOUBLE PRECISION :: RHI
C                   ! RELATIVE HUMINITY IN 0-100%
      DOUBLE PRECISION :: TEMPI
C                   ! AIR TEMPERATURE IN K
      DOUBLE PRECISION :: CP(17)
C                   ! THE OUTPUT SOLID AND LIQUID AEROSOL IN UG/M3
      DOUBLE PRECISION :: CGAS(4)
C                   ! THE OUTPUT GAS IN UG/M3
C 
C      DEFINE THE LOCAL ARRARS
C      
      DOUBLE PRECISION :: WI(5),GAS(3),AERLIQ(12),AERSLD(9),CNTRL(2),
     &                 WT(5),OTHER(6),RH
      
      DOUBLE PRECISION :: EMW(5),INTMW(17),GASMW(4)
C                   ! ENW  : THE MOLEWEIGHT FOR INPUT PRECURSORS
C                   ! INTMW: THE MOLEWEIGHT FOR OUTPUT AEROSOL              
C                   ! INTMW: THE MOLEWEIGHT FOR OUTPUT GASES
      DOUBLE PRECISION :: AWATER,CGAS1,CGAS2,CGAS3,CGAS4

      DATA EMW   /23.0D0, 98.0D0, 17.0D0, 63.0D0, 36.5D0/
      DATA INTMW /1.0D0 , 23.0D0, 18.0D0, 35.5D0, 96.0D0, 97.0D0,
     &            62.0D0, 58.5D0, 142.0D0,85.0D0, 132.0D0,80.0D0,
     &            53.5D0, 98.0D0, 115.0D0,120.0D0,247.0D0/
      DATA GASMW /17.0D0, 63.0D0, 98.0D0, 36.5D0 /

C
C**** CONVERT INPUT CONCENTRATIONS TO MOLES/M3 *************************
C
C        PRINT*,I,J,K,WI(5),'test1'

       WI(1) = QI(1)/EMW(1)*1.D-6     ! NA
       WI(2) = QI(2)/EMW(2)*1.D-6     ! H2SO4      
       WI(3) = QI(3)/EMW(3)*1.D-6     ! NH3
       WI(4) = QI(4)/EMW(4)*1.D-6     ! HNO3
       WI(5) = QI(5)/EMW(5)*1.D-6     ! HCL

       
       DO 100, N = 1, 5
         IF(WI(N).LT.1.D-10) WI(N) = 0.0D0
100    CONTINUE

C
C****  CALL ISORROPIA **************************************************

       RH       = RHI/100. ! TO CONVERT 0-100% to 0-1
       CNTRL(1) = 0D0        ! 0 = FORWARD PROBLEM ; 1 = REVERSE PROBLEM
       CNTRL(2) = 0D0        ! 0 = SOLID + LIQUID AEROSOL, 1 = METASTABLE
C         PRINT*, WI(1),WI(2),WI(3),WI(4),WI(5),TEMPI,RH
C
       CALL ISOROPIA ( WI, RH, TEMPI, CNTRL, 
     &                  WT, GAS, AERLIQ,AERSLD, SCASE, OTHER )
C
C***   CONVERT OUTPUT CPNCENTRATIONS TO UG/M3 FOR AEROSOLS
C
C AEROSOLS SPECIES
C 
C       PRINT*,I,J,K,CGAS(4),'test2'

       CP(1) = AERLIQ(1) * INTMW(1) *1.D6    ! H+(AQ)
       CP(2) = AERLIQ(2) * INTMW(2) *1.D6    ! NA+(AQ)
       CP(3) = AERLIQ(3) * INTMW(3) *1.D6    ! NH4+(AQ)
       CP(4) = AERLIQ(4) * INTMW(4) *1.D6    ! CL-(AQ)
       CP(5) = AERLIQ(5) * INTMW(5) *1.D6    ! SO4--(AQ)
       CP(6) = AERLIQ(6) * INTMW(6) *1.D6    ! HSO4-(AQ)
       CP(7) = AERLIQ(7) * INTMW(7) *1.D6    ! NO3-(AQ)
       CP(8) = AERSLD(3) * INTMW(8) *1.D6    ! NACL(S)
       CP(9) = AERSLD(5) * INTMW(9) *1.D6    ! NA2SO4(S)
       CP(10)= AERSLD(1) * INTMW(10)*1.D6    ! NANO3(S)
       CP(11)= AERSLD(6) * INTMW(11)*1.D6    ! NH42SO4(S)
       CP(12)= AERSLD(2) * INTMW(12)*1.D6    ! NH4NO3(S)
       CP(13)= AERSLD(4) * INTMW(13)*1.D6    ! NH4CL(S)
       CP(14)= 0.0       * INTMW(14)*1.D6    ! H2SO4(AQ) 
       CP(15)= AERSLD(8) * INTMW(15)*1.D6    ! NH4HSO4(S)
       CP(16)= AERSLD(7) * INTMW(16)*1.D6    ! NAHSO4
       CP(17)= AERSLD(9) * INTMW(17)*1.D6    ! (NH4)4H(SO4)2(S)
       AWATER = AERLIQ(8) * 18. * 1.E6         ! AEROSOL WATER
C        PRINT*, CP(3)
C
C GASEOUS SPECIES
C
       CGAS(1) = GAS(1) * GASMW(1) *1.D6     ! NH3(G)
       CGAS(2) = GAS(2) * GASMW(2) *1.D6     ! HNO3(G)
       CGAS(3) = 0.0D0                       ! H2SO4(G)
       CGAS(4) = GAS(3) * GASMW(4) *1.D6     ! HCL(G)
C
C       PRINT*,I,J,K,CGAS(4),'test3'
       
        IF(CP(1).LT.1.D-16) CP(1) = 0.0D0
        IF(CP(2).LT.1.D-16) CP(2) = 0.0D0
        IF(CP(3).LT.1.D-16) CP(3) = 0.0D0
        IF(CP(4).LT.1.D-16) CP(4) = 0.0D0
        IF(CP(5).LT.1.D-16) CP(5) = 0.0D0
        IF(CP(6).LT.1.D-16) CP(6) = 0.0D0  
        IF(CP(7).LT.1.D-16) CP(7) = 0.0D0
         CP(8) = 0.0D0
        IF(CP(9).LT.1.D-16)  CP(9) = 0.0D0
        IF(CP(10).LT.1.D-16) CP(10) = 0.0D0
        IF(CP(11).LT.1.D-16) CP(11) = 0.0D0
        IF(CP(12).LT.1.D-16) CP(12) = 0.0D0
        IF(CP(13).LT.1.D-16) CP(13) = 0.0D0
        IF(CP(14).LT.1.D-16) CP(14) = 0.0D0
        IF(CP(15).LT.1.D-16) CP(15) = 0.0D0
        IF(CP(16).LT.1.D-16) CP(16) = 0.0D0
        IF(CP(17).LT.1.D-16) CP(17) = 0.0D0
        IF(AWATER.LT.1.D-16) AWATER = 0.0D0

C       PRINT*,I,J,K,CGAS(4),'test4'


       DO 300, N = 1, 4
        IF(CGAS(N).LT.1.D-16) CGAS(N) = 0.0D0
        CGAS1 = CGAS(1)
        CGAS2 = CGAS(2)
        CGAS3 = CGAS(3)
        CGAS4 = CGAS(4)

  300  CONTINUE
 
C       PRINT*,I,J,K,CGAS(4),'test5'
      RETURN
C
C*** END OF SUBROUTINE ISRPINTR C****************************************
C
      END
       
       
     
      
