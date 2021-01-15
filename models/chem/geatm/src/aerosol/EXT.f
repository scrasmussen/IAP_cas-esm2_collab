C=====================================================================
C
C***  SUBROUTINE GETEXT
C***  THIS SUBROUTINE IS A SAMPLE INTERFACE BETWEEN AN AIRSHED MODEL
C***  AND AEROSOL SCATTERING EXTINCTIONS MODULE USING CMAQ DESCRIPTION.
C
C=====================================================================
C
      SUBROUTINE GETEXT( MYID, RH,  AER0101,AER0102,AER0103,AER0104,
     &   AER0201,AER0202,AER0203,AER0204,NH4, 
     &   SO4, NO3,NA, BC, OC, PM10, PM25,
     &   DUSTEXT,EXTASO4, EXTANO3, EXTANH4, EXTBC, EXTOC,EXT,EXTS, 
     &   ICHOICE, I, J, K)
C     AER0101  SEASALT 0.1-1UM
C     AER0102  SEASALT 1-2.5UM
C     AER0103  SEASALT 2.5-5UM
C     AER0104  SEASALT 5-10UM
C     AER0201  DUST 0.1-1UM
C     AER0202  DUST 1-2.5UM
C     AER0203  DUST 2.5-5UM
C     AER0204  DUST 5-10UM
C     NH4     TOTAL MASS CONCENTRATION OF NH4 IN UG/M3
C     SO4     TOTAL MASS CONCENTRATION OF SO4 IN UG/M3
C     NO3     TOTAL MASS CONCENTRATION OF NO3 IN UG/M3      
C     BC      TOTAL MASS CONCENTRATION OF BLACK CARBON  IN UG/M3
C     OC      TOTAL MASS CONCENTRATION OF POA + SOA IN UG/M3      
C     PM10    TOTAL MASS CONCENTRATION OF PM10 IN UG/M3, HERE DEFINED AS
C             SOIL PARTICULARS
C     PM25    TOTAL MASS CONCENTRATION OF PM25 IN UG/M3, HERE DEFINED
C             AS SOIL PARTICULARS 
C     EXT     THE EXTINCTIONS DUE TO PARTICULAR SCATTERING IN KM-1
C     EXTS    THE EXTINCTIONS DUE TO SCTTER 
C     EXTASO4
C     EXTANO3
C     EXTANH4
C     EXTBC
C     EXTOC
C     FRH    
C     SSA     THE SINGLE SCATTER ALBEDO      
C     AOD     THE AEROSOL OPITICAL DEPTH
C     RH      RELATIVE HUMINITY IN %
C     FRH     THE hygroscopic growth factor for ANH4 ANO3 ASO4
C     FRHSSA  THE hygroscopic growth factor for SEA SALT
C     NZ      TOTAL LAYER NUMBERS   
C     HUMFAC  RELATIVE THE HUMINITY BASED AEROSOL GROWTH FACTOR ON A
C             LOOK-UP TABLE  at 1% RH values
C     ICHOICE THE SCHEDULE CHOICE
C             1: NATURAL CONDITIONS FROM OLD PROPOSAL,FIXED COFFICIENTS
C             2: CHANGED SCHEDUME , VARIED COFFICIENTS FOR AEROSOLS
C               CONCENTRATIONS
C     IREGIM  THE TYPES OF REGIONS 1: DUST, 2: CLEAN, 3: URBAN, 
C             4: MARINE 5: MIXED MIXED MEANS THE MIXED AIR OF URBAN AND MARINE

      
      INTEGER MYID,ICHOICE,IREGIM,I,J,K
      REAL    AER0201,AER0202,AER0203,AER0204,NH4,SO4,NO3,BC,OC,PM10
      REAL    PM25,EXT,RH,NA,FRH,EXTS,EXTASO4,EXTANO3,EXTANH4,EXTBC,EXTOC
      REAL    HUMFAC(99)
      REAL    MASS,AINORG,AORG
      REAL    A(5),B(5)
C
C     FRH IS CALCULATED FRH = 1+A*(RH/100)**B BY OBS BY PAN, X.(2009)IN ACPD IN
C     BEIJING      
      DATA A/0.64, 1.20, 2.30, 4.92, 7.68 /
      DATA B/5.17, 6.07, 6.27, 5.04, 7.62/
     
C***   CMAQ TABLES BUT DON'T BE RIGHT IN CHINA    
C      DATA HUMFAC/ 
C     &   1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0000E+00,
C     &   1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0000E+00,
C     &   1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0001E+00,  1.0001E+00,
C     &   1.0004E+00,  1.0006E+00,  1.0024E+00,  1.0056E+00,  1.0089E+00,
C     &   1.0097E+00,  1.0105E+00,  1.0111E+00,  1.0115E+00,  1.0118E+00,
C     &   1.0122E+00,  1.0126E+00,  1.0130E+00,  1.0135E+00,  1.0139E+00,
C     &   1.0173E+00,  1.0206E+00,  1.0254E+00,  1.0315E+00,  1.0377E+00,
C     &   1.0486E+00,  1.0596E+00,  1.0751E+00,  1.0951E+00,  1.1151E+00,
C     &   1.1247E+00,  1.1343E+00,  1.1436E+00,  1.1525E+00,  1.1615E+00,
C     &   1.1724E+00,  1.1833E+00,  1.1955E+00,  1.2090E+00,  1.2224E+00,
C     &   1.2368E+00,  1.2512E+00,  1.2671E+00,  1.2844E+00,  1.3018E+00,
C     &   1.3234E+00,  1.3450E+00,  1.3695E+00,  1.3969E+00,  1.4246E+00,
C     &   1.4628E+00,  1.5014E+00,  1.5468E+00,  1.5992E+00,  1.6516E+00,
C     &   1.6991E+00,  1.7466E+00,  1.7985E+00,  1.8549E+00,  1.9113E+00,
C     &   1.9596E+00,  2.0080E+00,  2.0596E+00,  2.1146E+00,  2.1695E+00,
C     &   2.2630E+00,  2.3565E+00,  2.4692E+00,  2.6011E+00,  2.7330E+00,
C     &   2.8461E+00,  2.9592E+00,  3.0853E+00,  3.2245E+00,  3.3637E+00,
C     &   3.5743E+00,  3.7849E+00,  4.0466E+00,  4.3594E+00,  4.6721E+00,
C     &   5.3067E+00,  5.9412E+00,  6.9627E+00,  8.3710E+00,  9.7793E+00,
C     &   1.2429E+01,  1.5078E+01,  1.8059E+01,  2.1371E+01/
      IF( (AER0201+AER0202+AER0203+AER0204).GE.200.) THEN
        IREGIM = 1
      ELSE IF((NH4+SO4+NO3+PM10+PM25+BC+OC).GE.50..AND.NA.GE.10.) THEN
        IREGIM = 5
      ELSE IF(NA.GE.10.) THEN
        IREGIM = 4
      ELSE IF((NH4+SO4+NO3+PM10+PM25+BC+OC).GE.50.) THEN 
        IREGIM = 3
      ELSE 
        IREGIM = 2      
      ENDIF  

C *** BEGIN TO RECONSTRUCTED METHOD *******************************
C
C       IRH = INT( RH  )     ! truncate relative humidity to nearest integer
C       IRH = MIN( 99, IRH)  ! set maximum value on IRH
C       IRH = MAX( 1, IRH)   ! set minimum value on IRH
        FRH = 1. + A(IREGIM)*(RH/100.)**B(IREGIM)
C       FRH = HUMFAC( IRH )  ! set humidity correction
C
! hygroscopic growth factor for sea-salt from Chin et al. (2002)

        IF ( RH < 100 ) FRHSSA = 4.8
        IF ( RH < 099 ) FRHSSA = 2.9
        IF ( RH < 095 ) FRHSSA = 2.4
        IF ( RH < 090 ) FRHSSA = 2.0
        IF ( RH < 080 ) FRHSSA = 1.8
        IF ( RH < 070 ) FRHSSA = 1.6
        IF ( RH < 050 ) FRHSSA = 1.0


      
       IF (ICHOICE == 1 ) THEN      
       EXT    =  0.003 * FRH * ( NH4 + SO4 + NO3 ) 
     &         + 0.004 * OC * 1.6
     &         + 0.001 * PM25
     &         + 0.0006* PM10
     &         + 0.01  * BC
     &         + 0.0017* FRHSSA * (AER0101+AER0102+AER0103+AER0104) ! SEA SALT
     &         + 0.001* (AER0201+AER0202) ! FINE DUST
     &         + 0.0006* (AER0203+AER0204) ! COARSE DUST
     &         + 0.01  ! 0.01 is RAYLEIGH SCATTER      
       
       EXTS   =  0.003 * FRH * ( NH4 + SO4 + NO3 )
     &         + 0.004 * OC * 1.6
     &         + 0.001 * PM25
     &         + 0.0006* PM10   
     &         + 0.0017* FRHSSA * (AER0101+AER0102+AER0103+AER0104)  ! SEA SALT  
     &         + 0.001* (AER0201+AER0202) ! FINE DUST
     &         + 0.0006* (AER0203+AER0204) ! COARSE DUST
     &         + 0.01
       DUSTEXT = 0.001 * (AER0201+AER0202) ! FINE DUST
     &         + 0.0006* (AER0203+AER0204) ! COARSE DUST
       EXTASO4 = 0.003 * FRH * SO4
       EXTANO3 = 0.003 * FRH * NO3
       EXTANH4 = 0.003 * FRH * NH4
       EXTBC   = 0.01 * BC
       EXTOC   = 0.004 * OC * 1.6
 
          
       ELSE 

       MASS   = NH4 + SO4 + NO3 + OC 
       AINORG = 0.003   * ( 0.7 + 0.008 * MASS )
       AORG   = 0.00363 * ( 0.7 + 0.008 * MASS )

       EXT    = AINORG * FRH * ( NH4 + SO4 + NO3 )
     &         + AORG  * OC * 1.6
     &         + 0.001 * PM25
     &         + 0.0006* PM10
     &         + 0.01  * BC
     &         + 0.0017* FRHSSA * (AER0101+AER0102+AER0103+AER0104) ! SEA SALT
     &         + 0.001* (AER0201+AER0202) ! FINE DUST
     &         + 0.0006* (AER0203+AER0204) ! COARSE DUST
     &         + 0.01  ! 0.01 is RAYLEIGH SCATTER

       EXTS    = AINORG * FRH * ( NH4 + SO4 + NO3 )
     &         + AORG  * OC * 1.6
     &         + 0.001 * PM25
     &         + 0.001* PM10
     &         + 0.0017* FRHSSA * (AER0101+AER0102+AER0103+AER0104) ! SEA SALT
     &         + 0.001* (AER0201+AER0202) ! FINE DUST
     &         + 0.0006* (AER0203+AER0204) ! COARSE DUST
     &         + 0.01
       
      DUSTEXT = 0.001 * (AER0201+AER0202) ! FINE DUST
     &         + 0.0006* (AER0203+AER0204) ! COARSE DUST

      EXTASO4 = AINORG * FRH * SO4
      EXTANO3 = AINORG * FRH * NO3
      EXTANH4 = AINORG * FRH * NH4
      EXTBC   = AINORG * FRH * BC
      EXTOC   = AINORG * FRH * OC 

  
C       IF(i==57.AND.J==41.AND.K==1) PRINT*,EXT,EXTASO4,SO4
       ENDIF

C       
      RETURN
      END

      
