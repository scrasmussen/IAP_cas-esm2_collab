           SUBROUTINE HETERO_REAC(ASO4,BC,DUST01,DUST02,DUST03,DUST04,SSA01,SSA02,SSA03,SSA04,&
                               DUSTSO4,DUSTNO3,DSSASO4, DSSANO3,FSO4_DUST,FNO3_DUST,FSO4_SSA,FNO3_SSA)
            include 'chm1.inc'
            include 'gas1.inc' 
            REAL :: ASO4, BC  ! MASS CONCENTRATIONS IN UG/M3
!            REAL :: RH,TE ! IN % AND K
            REAL :: DENSI(2)   ! DENSITY OF AEROSOLS IN G/CM3 MEAN ! ASO4 AND BC
            REAL :: R(2)       ! MEAN RADIUS OF AEROSOLS IN UM
            REAL :: DENSIDUST(4), DENSISSA(4) ! DENSITY OF DUST AND SEA SALT FOR 4 BINS
            REAL :: RDUST(4), RSSA(4)         ! MEAN RADIUS OF DUST AND SEA SALT  IN UM 
            REAL :: AREA(2)    ! AEROSOL SURFACE AREA DENSITY UM2/CM3 FOR ASO4, BC
            REAL :: AREADUST(4),AREASSA(4) ! AEROSOL SURFACE AREA DENSITY UM2/CM3 FOR ASO4 AND BC
            REAL :: FRH,FRHSSA  ! FRH is the Hygroscopic growth for ASO4 and sea salt 
            REAL :: GAMMA1(28),A,B  ! 9 HETEROGENEOUS REACTIONS UPTAKE COFFICIENT 
            REAL :: DG(28)      ! gas diffusion in CM2/S  
            REAL :: C(28)       ! AVERAGE  molec speed [cm/s] 
!            REAL :: RK_HET(28)  ! 1/s HETEROGENEOUS rate 28 reactions  and have been defined in gas.inc
            REAL :: RK_HET_TMP(11:28, 4) ! 1/s HETEROGENEOUS rate FOR DUST AND SEA FOUR BINS
           !!!! 1: (NH4)2SO4 2: BC  3: DUST 4 : SEA SALT     
            REAL :: FSO4_DUST(4),FNO3_DUST(5),FSO4_SSA(4),FNO3_SSA(4) ! 4 BINS
            DATA DENSI/1.7, 1.0/
            DATA R/0.24, 0.04/
            DATA DENSIDUST / 2.5, 2.65, 2.65, 2.65 /
            DATA DENSISSA  / 2.2, 2.2 , 2.2 , 2.2  /
            DATA RDUST     / 0.15, 0.8,  1.75,  2.95  /
            DATA RSSA      / 0.15, 0.8,  1.75,  2.95  /

! 28 HETEROGENEOUS REACTIONS
!    1  N2O5 + ASO4 -> 2HNO3
!    2  NO2  + BC   -> 0.5HONO+0.5HNO3
!    3  NO3  + ASO4 -> HNO3
!    4  HO2  + ASO4 -> 0.5H2O2
!    5  HCHO + ASO4 -> PRODUCTS
!    6  OH   + ASO4 -> PRODUCTS
!    7  O3   + BC   -> PRODUCTS
!    8  NO2  + BC   -> HONO
!    9  HNO3 + BC   -> NO2
!   10  N2O5 + BC   -> 2HNO3
!   11  O3   + DUST -> PRODUCTS 
!   12  HNO3 + DUST -> ANO3 + PRODUCTS 
!   13  NO2  + DUST -> 0.5HONO + 0.5HNO3
!   14  NO3  + DUST -> HNO3
!   15  N2O5 + DUST -> 2HNO3
!   16  OH   + DUST -> PRODUCTS
!   17  HO2  + DUST -> 0.5H2O2
!   18  H2O2 + DUST -> PRODUCTS
!   19  SO2  + DUST -> ASO4 
!   20  CH3COOH + DUST -> PRODUCTS
!   21  CH3OH   + DUST -> PRODUCTS
!   22  HCHO    + DUST -> PRODUCTS
!   23  N2O5 + SSA  -> 2HNO3
!   24  NO3  + SSA  -> HNO3
!   25  HO2  + SSA  -> 0.5HONO
!   26  SO2  + SSA  -> ASO4
!   27  NO3  + SSA  -> ANO3
!   28  HNO3 + SSA  -> ANO3



            DATA GAMMA1 / 0.1    , 1.E-4,  3.E-3,  2.5E-1, 2.2E-2,&
                          0.2    , 3.3E-4, 3.3E-4, 2.1e-2, 5.E-3, &
                          2.7E-5 , 1.7E-1, 2.1E-6, 1.0E-3, 3.0E-2,&
                          0.1    , 0.2,    2.0E-3, 1.0E-4, 1.0E-3,&
                          1.0E-5 , 1.0E-5, 5.0E-3, 1.0E-3, 2.0E-1,&
                          5.0E-2 , 1.7E-2, 0.5/

            DATA DG / 0.1, 0.1, 0.1, 0.25, 0.1, 0.1, 0.1, 0.1,0.1, & 
                    ! N2O5,NO2, NO3, HO2, HCHO, OH, O3,  NO2, HNO3
                      0.1, 0.1, 0.1, 0.1,  0.1, 0.1, 0.1, 0.25, 0.1, &
                    ! N2O5, O3, HNO3,NO2, NO3, N2O5,OH, HO2,H2O2,
                      0.1, 0.1, 0.1, 0.1,  0.1, 0.1, 0.25, 0.1, 0.1,&
                    ! SO2, CH3COOH, CH3OH, HCHO, N2O5, NO3, HO2, SO2, NO3
                      0.1 /
                    ! HNO3

!      FRH IS CALCULATED FRH = 1+A*(RH/100)**B BY OBS BY PAN, X.(2009)IN  ACPD IN   BEIJING  
            FRH = 1. + 2.3 * (RH/100.)**6.27
           
        IF ( RH < 100 ) FRHSSA = 4.8
        IF ( RH < 099 ) FRHSSA = 2.9
        IF ( RH < 095 ) FRHSSA = 2.4
        IF ( RH < 090 ) FRHSSA = 2.0
        IF ( RH < 080 ) FRHSSA = 1.8
        IF ( RH < 070 ) FRHSSA = 1.6
        IF ( RH < 050 ) FRHSSA = 1.0

            AREA(1) = 3.0 * ASO4 / ( DENSI(1) * R (1)) * FRH**2
            AREA(2) = 3.0 * BC   / ( DENSI(2) * R (2)) 
            AREADUST(1) = 3.0 * DUST01 / (  DENSIDUST(1) * RDUST (1) )
            AREADUST(2) = 3.0 * DUST02 / (  DENSIDUST(2) * RDUST (2) )
            AREADUST(3) = 3.0 * DUST03 / (  DENSIDUST(3) * RDUST (3) )
            AREADUST(4) = 3.0 * DUST04 / (  DENSIDUST(4) * RDUST (4) )
            AREASSA(1) = 3.0 * SSA01 / (  DENSISSA(1) * RSSA (1) ) * FRHSSA**2
            AREASSA(2) = 3.0 * SSA02 / (  DENSISSA(2) * RSSA (2) ) * FRHSSA**2 
            AREASSA(3) = 3.0 * SSA03 / (  DENSISSA(3) * RSSA (3) ) * FRHSSA**2
            AREASSA(4) = 3.0 * SSA04 / (  DENSISSA(4) * RSSA (4) ) * FRHSSA**2

!    TO ADJUST GAMMA DEPENDING ON RH AND TEMPERATURE
            GAMMA1(7) = 1.8e-4 * EXP (-1000./TE)
            A = 2.79E-4 + 1.3E-4*RH - 3.43E-6 * RH**2.&
               + 7.52E-8 * RH**3. 
            IF(TE .GE. 282) B = 4.E-2*(TE-294)
            IF(TE .LT. 282) B = 0.48
            GAMMA1(1) = A * 10.**B

            IF (RH < 62  ) GAMMA1 (23) = 0.005
            IF (RH >= 62 ) GAMMA1 (23) = 0.03
            IF (RH < 50  ) GAMMA1 (26) = 0.005
            IF (RH >=50  ) GAMMA1 (26) = 0.05

            IF(RH <= 15) THEN
                GAMMA1(18) = 3.33E-4
            ELSE IF(RH <= 25) THEN
                GAMMA1(18) = 3.5E-4
            ELSE IF(RH <= 35) THEN
                GAMMA1(18) = 3.55E-4
            ELSE IF(RH <= 40) THEN
                GAMMA1(18) = 3.6E-4
            ELSE IF(RH <= 50) THEN
                GAMMA1(18) = 4.1E-4
            ELSE IF(RH <= 60) THEN
                GAMMA1(18) = 4.6E-4
            ELSE IF(RH <= 65) THEN
                GAMMA1(18) = 5.2E-4
            ELSE IF(RH <= 70) THEN
                GAMMA1(18) = 6.03E-4
            ELSE
                GAMMA1(18) = 6.03E-4
            ENDIF

                RH_TMP = RH / 100. - 0.15
                RH_TMP = MAX(RH_TMP,0.1)
                TMP1 = 8. * RH_TMP
                TMP2 = (1.-8.)* RH_TMP
                BET  = TMP1 / (1.-RH_TMP) / (1.-TMP2)
                GAMMA1(12) = BET * 0.033
                GAMMA1(12) = GAMMA1(12) * 0.5535 - 0.0058
! GAMMA1(12) dericed by Vlasenko (2006,ACP) and Wei (2010, Phd thesis Modelng the effect of Heteorogeneous reqction on atmospgeric chemistry and aerosol properties ) , using liear with latetr data with former datea
! ***   TEST IMPACTS of HETEOGENEOUS CHEMISTRY
!                GAMMA1(12) = 0.
                       
!    TO ESTIMATE THE avg. molec speed [cm/s]

            C(1) = 1.455e4 * sqrt(te/108.) ! N2O5
            C(2) = 1.455e4 * sqrt(te/46.)  ! NO2
            C(3) = 1.455e4 * sqrt(te/62.)  ! NO3
            C(4) = 1.455e4 * sqrt(te/33.)  ! HO2
            C(5) = 1.455e4 * sqrt(te/30.)  ! HCHO
            C(6) = 1.455e4 * sqrt(te/17.)  ! OH
            C(7) = 1.455e4 * sqrt(te/48.)  ! O3
            C(8) = 1.455e4 * sqrt(te/46.)  ! NO2
            C(9) = 1.455e4 * sqrt(te/63.)  ! HNO3
            C(10) = 1.455e4 * sqrt(te/108.)! N2O5
            C(11) = 1.455e4 * sqrt(te/48.) ! O3
            C(12) = 1.455e4 * sqrt(te/63.) ! HNO3
            C(13) = 1.455e4 * sqrt(te/46.) ! NO2
            C(14) = 1.455e4 * sqrt(te/62.) ! NO3
            C(15) = 1.455e4 * sqrt(te/108.)! N2O5
            C(16) = 1.455e4 * sqrt(te/17.) ! OH
            C(17) = 1.455e4 * sqrt(te/33.) ! HO2
            C(18) = 1.455e4 * sqrt(te/24.) ! H2O2
            C(19) = 1.455e4 * sqrt(te/80.) ! SO2
            C(20) = 1.455e4 * sqrt(te/60.) ! CH3COOH
            C(21) = 1.455e4 * sqrt(te/32.) ! CH3OH
            C(22) = 1.455e4 * sqrt(te/30.) ! HCHO
            C(23) = 1.455e4 * sqrt(te/108.)! N2O5
            C(24) = 1.455e4 * sqrt(te/62.) ! NO3
            C(25) = 1.455e4 * sqrt(te/33.) ! HO2
            C(26) = 1.455e4 * sqrt(te/80.) ! SO2
            C(27) = 1.455e4 * sqrt(te/62.) ! NO3
            C(28) = 1.455e4 * sqrt(te/63.) ! HNO3


!    TO ESTIMATE THE HETEROGENEOUS RATE CONSTANHT 
!               K = [r/Dg + 4/ (c * gamma)]-1 A 

            
            RK_HET(1) = 1./( 1.E-4*R(1)/DG(1) + 4./C(1)/GAMMA1(1)) * AREA(1) * 1.E-8
            RK_HET(2) = 1./( 1.E-4*R(2)/DG(2) + 4./C(2)/GAMMA1(2)) * AREA(2) * 1.E-8             
            RK_HET(3) = 1./( 1.E-4*R(1)/DG(3) + 4./C(3)/GAMMA1(3)) * AREA(1) * 1.E-8
            RK_HET(4) = 1./( 1.E-4*R(1)/DG(4) + 4./C(4)/GAMMA1(4)) * AREA(1) * 1.E-8
            RK_HET(5) = 1./( 1.E-4*R(1)/DG(5) + 4./C(5)/GAMMA1(5)) * AREA(1) * 1.E-8
            RK_HET(6) = 1./( 1.E-4*R(1)/DG(6) + 4./C(6)/GAMMA1(6)) * AREA(1) * 1.E-8
            RK_HET(7) = 1./( 1.E-4*R(2)/DG(7) + 4./C(7)/GAMMA1(7)) * AREA(2) * 1.E-8
            RK_HET(8) = 1./( 1.E-4*R(2)/DG(8) + 4./C(8)/GAMMA1(8)) * AREA(2) * 1.E-8
            RK_HET(9) = 1./( 1.E-4*R(2)/DG(9) + 4./C(9)/GAMMA1(9)) * AREA(2) * 1.E-8
            RK_HET(10)= 1./( 1.E-4*R(2)/DG(10) + 4./C(10)/GAMMA1(10)) * AREA(2) * 1.E-8

           DO I = 11, 22   ! FOR DUST
              RK_HET_TMP (I, 1) = 1./( 1.E-4*RDUST(1)/DG(I) + 4./C(I)/GAMMA1(I)) * AREADUST(1) * 1.E-8
              RK_HET_TMP (I, 2) = 1./( 1.E-4*RDUST(2)/DG(I) + 4./C(I)/GAMMA1(I)) * AREADUST(2) * 1.E-8
              RK_HET_TMP (I, 3) = 1./( 1.E-4*RDUST(3)/DG(I) + 4./C(I)/GAMMA1(I)) * AREADUST(3) * 1.E-8
              RK_HET_TMP (I, 4) = 1./( 1.E-4*RDUST(4)/DG(I) + 4./C(I)/GAMMA1(I)) * AREADUST(4) * 1.E-8

              RK_HET(I) = RK_HET_TMP (I, 1) + RK_HET_TMP (I, 2) + RK_HET_TMP (I, 3) + RK_HET_TMP (I, 4)
           ENDDO 

           DO I = 23, 28   ! FOR SEA SALT (SSA)
              RK_HET_TMP (I, 1) = 1./( 1.E-4*RSSA(1)/DG(I) + 4./C(I)/GAMMA1(I)) * AREASSA(1) * 1.E-8
              RK_HET_TMP (I, 2) = 1./( 1.E-4*RSSA(2)/DG(I) + 4./C(I)/GAMMA1(I)) * AREASSA(2) * 1.E-8
              RK_HET_TMP (I, 3) = 1./( 1.E-4*RSSA(3)/DG(I) + 4./C(I)/GAMMA1(I)) * AREASSA(3) * 1.E-8
              RK_HET_TMP (I, 4) = 1./( 1.E-4*RSSA(4)/DG(I) + 4./C(I)/GAMMA1(I)) * AREASSA(4) * 1.E-8

              RK_HET(I) = RK_HET_TMP (I, 1) + RK_HET_TMP (I, 2) + RK_HET_TMP (I, 3) + RK_HET_TMP (I, 4)

           ENDDO
 
           DO I = 1,28
             RK_HET(I) = AMAX1(RK_HET(I),1.E-20)
           ENDDO
              
              RK_ASO4 = RK_HET(19) + RK_HET(26)
              RK_ANO3 = RK_HET(12) + RK_HET(27) + RK_HET(28)

              FSO4_DUST (1) = RK_HET_TMP(19,1) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_DUST (2) = RK_HET_TMP(19,2) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_DUST (3) = RK_HET_TMP(19,3) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_DUST (4) = RK_HET_TMP(19,4) / AMAX1(RK_ASO4, 1.E-20)
              
              FNO3_DUST (1) = RK_HET_TMP(12,1) / AMAX1(RK_ANO3, 1.E-20)
              FNO3_DUST (2) = RK_HET_TMP(12,2) / AMAX1(RK_ANO3, 1.E-20)
              FNO3_DUST (3) = RK_HET_TMP(12,3) / AMAX1(RK_ANO3, 1.E-20)
              FNO3_DUST (4) = RK_HET_TMP(12,4) / AMAX1(RK_ANO3, 1.E-20)

              FSO4_SSA (1) = RK_HET_TMP(26,1) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_SSA (2) = RK_HET_TMP(26,2) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_SSA (3) = RK_HET_TMP(26,3) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_SSA (4) = RK_HET_TMP(26,4) / AMAX1(RK_ASO4, 1.E-20)

              FNO3_SSA (1) = ( RK_HET_TMP(27,1) + RK_HET_TMP(28,1) )/ AMAX1(RK_ANO3, 1.E-20)
              FNO3_SSA (2) = ( RK_HET_TMP(27,2) + RK_HET_TMP(28,2) )/ AMAX1(RK_ANO3, 1.E-20)
              FNO3_SSA (3) = ( RK_HET_TMP(27,3) + RK_HET_TMP(28,3) )/ AMAX1(RK_ANO3, 1.E-20)
              FNO3_SSA (4) = ( RK_HET_TMP(27,4) + RK_HET_TMP(28,4) )/ AMAX1(RK_ANO3, 1.E-20)

           IF ( DUSTNO3.GE. 0.2*(DUST01+DUST02+DUST03+DUST04) ) RK_HET(12) = 0.0
           IF ( DSSASO4.GE. (SSA01+SSA02+SSA03+SSA04) )  RK_HET(26) = 0.0
           IF ( DSSANO3.GE. (SSA01+SSA02+SSA03+SSA04) )  RK_HET(28) = 0.0

              DO IR = 1, 10 
!                RK_HET ( IR ) = 0.0           
              ENDDO 

           END
