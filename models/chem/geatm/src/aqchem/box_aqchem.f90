      SUBROUTINE AQUEOUS(T,PLEV,QVAPOR,CLW,CGAS, CAER, CPH, DT, SX,EX,SY,EY, I,J,K)
        INTEGER :: SX,EX,SY,EY,I,J,K
        REAL    :: TCELL,PCELL,TCHM,WCHM,DT
        REAL    :: CLDPH ! CLOUD PH VALUE
        REAL    :: H2  ! hydrogen concentration (ppm)
        REAL    :: O2  ! oxygen concentration (ppm)
        REAL    :: CH4 ! methane concentration (ppm)
        REAL    :: atm ! total gas concentration (M) in ppm
        REAL    :: CLIQ 
        REAL    :: convfac ! UNIT CHANGES FROM UMOL/M3 TO PPM
        REAL    :: TAMIN   ! Cloud water freezing threshold (K)
        REAL    :: CWMIN   ! Minimum cloud water threshold (g/m3)
        REAL    :: RH,PRES_PA, CW_KGM3
        REAL    :: R_GAS(11),R_AER(9)
        REAL    :: T,PLEV,CPH ! CPH IS CLOUD PH
        REAL    :: QVAPOR,WATER  ! WATER VAPOR IN KG/KG AND PPM
        REAL    :: CWC           ! cloud water content (g/m3):
        REAL    :: CLW           ! cloud water content IN KG/KG
        REAL    :: CGAS(11)
        REAL    :: CAER(9)     ! THE INPUT GAS AND AEROSOLS CONCENTRATIONS IN PPB AND UG/M3
        REAL    :: MAXBLD, MINBLD 

        DATA CWMIN / 0.05 /
        DATA TAMIN / 243./
! CGASNAME /SO2, HNO3, NXOY, CO2, NH3, HH2O2, O3, FOA, MHP, PAA, H2SO4/
! CAERNAME / ASO4, ANH4, ANO3, ACACO3, AMGCO3, ANACL, AFEIII, AMNII,POTCL/

          TCELL =    T 
          PCELL = PLEV ! HPA
          TCHM  = TCELL
          WATER = QVAPOR * 29. /18. * 1.E06 ! FROM KG/KG TO PPM
          WCHM  = WATER
                    
          O2  = 2.095e5
          CH4 = 1.75
          H2  = 0.60
          ATM = 1.E06
          convfac = 44.9 * (273./TCELL)*(PCELL/1013.)

! ****   FROM KG/KG TO G/M3
          
         CWC  = CLW  * 29. * convfac 
          
         CLIQ = CWC 

         IF ( TCELL.LT.273) THEN
          CLIQ = AMAX1 (0., CLIQ* (TCELL - TAMIN) / (273. - TAMIN) )
         ENDIF
 
           

!CCCCCCCCCCCCCC  DO RADM AQUEOUS CHEMISTRY IF CWC IS ABOVE THRESHOLD 
!CCCCCCCCCC    ALL CONC UNITS MYST BE MOL/MOL (MIXING RATIO)            
         IF ( CLIQ.GE.CWMIN .AND. TCELL.GE.TAMIN ) THEN
           PRES_PA = 100. * PCELL
           CW_KGM3 = CLIQ/1000.

           DO IG = 1 ,11
            R_GAS(IG) = CGAS(IG)*1.E-9 
           ENDDO
!            IF(I.EQ.16.AND.J.EQ.68) PRINT*, CGAS(16,68,4),R_GAS(4)
           R_AER(1) = (CAER(1)/96./CONVFAC)*1.E-6
           R_AER(2) = (CAER(2)/18./CONVFAC)*1.E-6
           R_AER(3) = (CAER(3)/62./CONVFAC)*1.E-6  
           R_AER(4) = (CAER(4)/100./CONVFAC)*1.E-6
           R_AER(5) = (CAER(5)/84./CONVFAC)*1.E-6
           R_AER(6) = (CAER(6)/58./CONVFAC)*1.E-6
           R_AER(7) = (CAER(7)/56./CONVFAC)*1.E-6
           R_AER(8) = (CAER(8)/55./CONVFAC)*1.E-6
           R_AER(9) = (CAER(9)/74./CONVFAC)*1.E-6
            
           CLDPH = CPH 

           IF(CLDPH.EQ.-1.E20) CLDPH = 5.0
              
            
           IF(R_AER(1).LT.5.0E-7.AND.R_AER(3).LT.1.5E-6) THEN

!             IF(I.EQ.46.AND.J.EQ.44.AND.K.EQ.10) PRINT*,CLDPH,'11111'
                                !R_AER(1),R_AER(2),R_AER(3),'111'
!            PRINT*, I,J,K,CLDPH,'111'
            CALL RAQCHEM(TCELL,PRES_PA,DT,CW_KGM3,R_GAS,R_AER,CLDPH,I,J,K) 
!             IF(I.EQ.46.AND.J.EQ.44.AND.K.EQ.10) PRINT*,CLDPH,R_AER(1),R_AER(2),R_AER(3),'222'     
!            PRINT*,I,J,K,CLDPH,'222'    
           ENDIF
          

           DO IG = 1, 11
            MINBLD =  CGAS(IG)*0.8
            MAXBLD =  CGAS(IG)*1.2

            IF(R_GAS (IG)* 1.E9 .GE. MAXBLD) THEN
             CGAS(IG) = MAXBLD
            ELSE  IF(R_GAS (IG)*1.E9 .LT.MINBLD) THEN
             CGAS(IG) = MINBLD
            ELSE
             CGAS(IG) = R_GAS (IG)* 1.E9
            ENDIF
             CGAS(IG) = AMAX1(CGAS(IG),1.E-20)

           ENDDO

          

           CAER(1) = AMAX1(R_AER(1)* CONVFAC*96.*1.E6, 1.E-20) ! PSO4 UG/M3
           CAER(2) = AMAX1(R_AER(2)* CONVFAC*18.*1.E6, 1.E-20) ! PNH4 UG/M3
           CAER(3) = AMAX1(R_AER(3)* CONVFAC*62.*1.E6, 1.E-20) ! PNO3 UG/M3
           CAER(6) = AMAX1(R_AER(6)* CONVFAC*58.*1.E6, 1.E-20) ! PCL  UG/M3
          
           CPH   = CLDPH

         ENDIF  ! CLIQ

         

      RETURN
      ENDSUBROUTINE
