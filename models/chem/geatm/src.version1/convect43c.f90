module convect43c
contains
!***********************************************************************
!****                       SUBROUTINE CONVECT                        **
!****                          VERSION 4.3c                           **
!****                          20 May, 2002                           **
!****                          Kerry Emanuel                          **
!***********************************************************************
!                                                                       
      SUBROUTINE CONVECT (T, Q, QS, U, V, TRA, P, PH, ND, NL, NTRA,     &
      DELT, FTRA, FTRAD, FTRAU, FTRAO, FTRAE, CBMF, IFLAG)              
!                                                                       
!-----------------------------------------------------------------------
!    *** On input:      ***                                             
!                                                                       
!     T:   Array of absolute temperature (K) of dimension ND, with first
!           index corresponding to lowest model level. Note that this ar
!           will be altered by the subroutine if dry convective adjustme
!           occurs and if IPBL is not equal to 0.                       
!                                                                       
!     Q:   Array of specific humidity (gm/gm) of dimension ND, with firs
!            index corresponding to lowest model level. Must be defined 
!            at same grid levels as T. Note that this array will be alte
!            if dry convective adjustment occurs and if IPBL is not equa
!                                                                       
!     QS:  Array of saturation specific humidity of dimension ND, with f
!            index corresponding to lowest model level. Must be defined 
!            at same grid levels as T. Note that this array will be alte
!            if dry convective adjustment occurs and if IPBL is not equa
!                                                                       
!     U:   Array of zonal wind velocity (m/s) of dimension ND, witth fir
!            index corresponding with the lowest model level. Defined at
!            same levels as T. Note that this array will be altered if  
!            dry convective adjustment occurs and if IPBL is not equal t
!                                                                       
!     V:   Same as U but for meridional velocity.                       
!                                                                       
!     TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA)
!            where NTRA is the number of different tracers. If no       
!            convective tracer transport is needed, define a dummy      
!            input array of dimension (ND,1). Tracers are defined at    
!            same vertical levels as T. Note that this array will be alt
!            if dry convective adjustment occurs and if IPBL is not equa
!                                                                       
!     P:   Array of pressure (mb) of dimension ND, with first           
!            index corresponding to lowest model level. Must be defined 
!            at same grid levels as T.                                  
!                                                                       
!     PH:  Array of pressure (mb) of dimension ND+1, with first index   
!            corresponding to lowest level. These pressures are defined 
!            levels intermediate between those of P, T, Q and QS. The fi
!            value of PH should be greater than (i.e. at a lower level t
!            the first value of the array P.                            
!                                                                       
!     ND:  The dimension of the arrays T,Q,QS,P,PH,FT and FQ            
!                                                                       
!     NL:  The maximum number of levels to which convection can         
!            penetrate, plus 1.                                         
!            NL MUST be less than or equal to ND-1.                     
!                                                                       
!     NTRA:The number of different tracers. If no tracer transport      
!            is needed, set this equal to 1. (On most compilers, setting
!            NTRA to 0 will bypass tracer calculation, saving some CPU.)
!                                                                       
!     DELT: The model time step (sec) between calls to CONVECT          
!                                                                       
!-----------------------------------------------------------------------
!    ***   On Output:         ***                                       
!                                                                       
!     IFLAG: An output integer whose value denotes the following:       
!                                                                       
!                VALUE                        INTERPRETATION            
!                -----                        --------------            
!                  0               No moist convection; atmosphere is no
!                                  unstable, or surface temperature is l
!                                  than 250 K or surface specific humidi
!                                  is non-positive.                     
!                                                                       
!                  1               Moist convection occurs.             
!                                                                       
!                  2               No moist convection: lifted condensat
!                                  level is above the 200 mb level.     
!                                                                       
!                  3               No moist convection: cloud base is hi
!                                  then the level NL-1.                 
!                                                                       
!                  4               Moist convection occurs, but a CFL co
!                                  on the subsidence warming is violated
!                                  does not cause the scheme to terminat
!                                                                       
!     FT:   Array of temperature tendency (K/s) of dimension ND, defined
!             grid levels as T, Q, QS and P.                            
!                                                                       
!     FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimensi
!             defined at same grid levels as T, Q, QS and P.            
!                                                                       
!     FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND, 
!             defined at same grid levels as T.                         
!                                                                       
!     FV:   Same as FU, but for forcing of meridional velocity.         
!                                                                       
!     FTRA: Array of forcing of tracer content, in tracer mixing ratio p
!             second, defined at same levels as T. Dimensioned (ND,NTRA)
!                                                                       
!     FTRAU: Array of forcing of tracer content from the lower layer, in
!            per second, defined at same levels as T. Dimensioned (ND,NT
!     FTRAD: same as FTRAU, but from upper layer                        
!                                                                       
!     PRECIP: Scalar convective precipitation rate (mm/day).            
!                                                                       
!     WD:    A convective downdraft velocity scale. For use in surface  
!             flux parameterizations. See convect.ps file for details.  
!                                                                       
!     TPRIME: A convective downdraft temperature perturbation scale (K).
!              For use in surface flux parameterizations. See convect.ps
!              file for details.                                        
!                                                                       
!     QPRIME: A convective downdraft specific humidity                  
!              perturbation scale (gm/gm).                              
!              For use in surface flux parameterizations. See convect.ps
!              file for details.                                        
!                                                                       
!     CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE 
!              BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT 
!              ITS NEXT CALL. That is, the value of CBMF must be "rememb
!              by the calling program between calls to CONVECT.         
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!    ***  THE PARAMETER NA SHOULD IN GENERAL BE GREATER THAN   ***      
!    ***                OR EQUAL TO  ND + 1                    ***      
!                                                                       
      PARAMETER (NA = 70) 
!                                                                       
      INTEGER NENT (NA) 
      REAL T (ND), Q (ND), QS (ND), U (ND), V (ND), TRA (ND, NTRA),     &
      P (ND), PH (ND)                                                   
      REAL FT (ND), FQ (ND), FU (ND), FV (ND), FTRA (ND, NTRA) 
      REAL FTRAU (ND, NTRA), FTRAD (ND, NTRA), FTRAO (ND, NTRA),        &
      FTRAE (ND, NTRA)                                                  
      REAL UENT (NA, NA), VENT (NA, NA), TRAENT (NA, NA, NTRA), TRATM ( &
      NA)                                                               
      REAL UP (NA), VP (NA), TRAP (NA, NTRA) 
      REAL M (NA), MP (NA), MENT (NA, NA), QENT (NA, NA), ELIJ (NA, NA) 
      REAL SIJ (NA, NA), TVP (NA), TV (NA), WATER (NA) 
      REAL QP (NA), EP (NA), TH (NA), WT (NA), EVAP (NA), CLW (NA) 
      REAL SIGP (NA), TP (NA), TOLD (NA), CPN (NA) 
      REAL LV (NA), LVCP (NA), LV0, H (NA), HP (NA), GZ (NA), HM (NA) 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!   ***                     Specify Switches                         ***
!                                                                       
!   ***   IPBL: Set to zero to bypass dry adiabatic adjustment       ***
!   ***    Any other value results in dry adiabatic adjustment       ***
!   ***     (Zero value recommended for use in models with           ***
!   ***                   boundary layer schemes)                    ***
!                                                                       
!   ***   MINORIG: Lowest level from which convection may originate  ***
!   ***     (Should be first model level at which T is defined       ***
!   ***      for models using bulk PBL schemes; otherwise, it should ***
!   ***      be the first model level at which T is defined above    ***
!   ***                      the surface layer)                      ***
!                                                                       
      IPBL = 0 
      IF (P (1) >800.) THEN 
         MINORIG = 9 
      ELSEIF (P (1) >600.) THEN 
         MINORIG = 7 
      ELSE 
         MINORIG = 6 
      ENDIF 
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!   ***                    SPECIFY PARAMETERS                        ***
!                                                                       
!   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
!   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
!   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
!   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
!   ***               BETWEEN 0 C AND TLCRIT)                        ***
!   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
!   ***                       FORMULATION                            ***
!   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
!   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
!   ***                        OF CLOUD                              ***
!   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
!   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
!   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF RAIN                             ***
!   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF SNOW                             ***
!   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
!   ***                         TRANSPORT                            ***
!   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
!   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
!   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
!   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
!   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
!   ***                   (DAMP MUST BE LESS THAN 1)                 ***
!                                                                       
      ELCRIT = .0011 
      TLCRIT = - 55.0 
      ENTP = 1.5 
      SIGD = 0.05 
      SIGS = 0.12 
      OMTRAIN = 50.0 
      OMTSNOW = 5.5 
      COEFFR = 1.0 
      COEFFS = 0.8 
      CU = 0.7 
      BETA = 10.0 
      DTMAX = 0.9 
      ALPHA = 0.2 
      DAMP = 0.1 
!                                                                       
!   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***     
!   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***     
!   ***             THESE SHOULD BE CONSISTENT WITH             ***     
!   ***              THOSE USED IN CALLING PROGRAM              ***     
!   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***     
!                                                                       
      CPD = 1005.7 
      CPV = 1870.0 
      CL = 2500.0 
      RV = 461.5 
      RD = 287.04 
      LV0 = 2.501E6 
      G = 9.8 
      ROWL = 1000.0 
!                                                                       
      CPVMCL = CL - CPV 
      EPS = RD / RV 
      EPSI = 1. / EPS 
      GINV = 1.0 / G 
      DELTI = 1.0 / DELT 
!                                                                       
!           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***           
!                                                                       
      DO 5 I = 1, ND 
         FT (I) = 0.0 
         FQ (I) = 0.0 
         FU (I) = 0.0 
         FV (I) = 0.0 
         DO 4 J = 1, NTRA 
            FTRA (I, J) = 0.0 
            FTRAU (I, J) = 0.0 
            FTRAD (I, J) = 0.0 
            FTRAO (I, J) = 0.0 
            FTRAE (I, J) = 0.0 
    4    END DO 
    5 END DO 
      DO 7 I = 1, NL + 1 
         RDCP = (RD * (1. - Q (I) ) + Q (I) * RV) / (CPD * (1. - Q (I) )&
         + Q (I) * CPV)                                                 
         TH (I) = T (I) * (1000.0 / P (I) ) **RDCP 
    7 END DO 
      PRECIP = 0.0 
      WD = 0.0 
      TPRIME = 0.0 
      QPRIME = 0.0 
      IFLAG = 0 
!                                                                       
      IF (IPBL.NE.0) THEN 
!                                                                       
!     ***            PERFORM DRY ADIABATIC ADJUSTMENT            ***    
!                                                                       
         JC = 0 
         DO 30 I = NL - 1, 1, - 1 
            JN = 0 
            SUM = TH (I) * (1. + Q (I) * EPSI - Q (I) ) 
            DO 10 J = I + 1, NL 
               SUM = SUM + TH (J) * (1. + Q (J) * EPSI - Q (J) ) 
               THBAR = SUM / FLOAT (J + 1 - I) 
               IF ( (TH (J) * (1. + Q (J) * EPSI - Q (J) ) ) .LT.THBAR) &
               JN = J                                                   
   10       END DO 
            IF (I.EQ.1) JN = MAX (JN, 2) 
            IF (JN.EQ.0) GOTO 30 
   12       CONTINUE 
            AHM = 0.0 
            RM = 0.0 
            UM = 0.0 
            VM = 0.0 
            DO K = 1, NTRA 
            TRATM (K) = 0.0 
            ENDDO 
            DO 15 J = I, JN 
               AHM = AHM + (CPD * (1. - Q (J) ) + Q (J) * CPV) * T (J)  &
               * (PH (J) - PH (J + 1) )                                 
               RM = RM + Q (J) * (PH (J) - PH (J + 1) ) 
               UM = UM + U (J) * (PH (J) - PH (J + 1) ) 
               VM = VM + V (J) * (PH (J) - PH (J + 1) ) 
               DO K = 1, NTRA 
               TRATM (K) = TRATM (K) + TRA (J, K) * (PH (J) - PH (J + 1)&
               )                                                        
               ENDDO 
   15       END DO 
            DPHINV = 1. / (PH (I) - PH (JN + 1) ) 
            RM = RM * DPHINV 
            UM = UM * DPHINV 
            VM = VM * DPHINV 
            DO K = 1, NTRA 
            TRATM (K) = TRATM (K) * DPHINV 
            ENDDO 
            A2 = 0.0 
            DO 20 J = I, JN 
               Q (J) = RM 
               U (J) = UM 
               V (J) = VM 
               DO K = 1, NTRA 
               TRA (J, K) = TRATM (K) 
               ENDDO 
               RDCP = (RD * (1. - Q (J) ) + Q (J) * RV) / (CPD *        &
               (1. - Q (J) ) + Q (J) * CPV)                             
               X = (0.001 * P (J) ) **RDCP 
               TOLD (J) = T (J) 
               T (J) = X 
               A2 = A2 + (CPD * (1. - Q (J) ) + Q (J) * CPV) * X *      &
               (PH (J) - PH (J + 1) )                                   
   20       END DO 
            DO 25 J = I, JN 
               TH (J) = AHM / A2 
               T (J) = T (J) * TH (J) 
               TC = TOLD (J) - 273.15 
               ALV = LV0 - CPVMCL * TC 
               QS (J) = QS (J) + QS (J) * (1. + QS (J) * (EPSI - 1.) )  &
               * ALV * (T (J) - TOLD (J) ) / (RV * TOLD (J) * TOLD (J) )
   25       END DO 
            IF ( (TH (JN + 1) * (1. + Q (JN + 1) * EPSI - Q (JN + 1) ) )&
            .LT. (TH (JN) * (1. + Q (JN) * EPSI - Q (JN) ) ) ) THEN     
               JN = JN + 1 
               GOTO 12 
            ENDIF 
            IF (I.EQ.1) JC = JN 
   30    END DO 
!                                                                       
!   ***   Remove any supersaturation that results from adjustment ***   
!                                                                       
         IF (JC.GT.1) THEN 
            DO 38 J = 1, JC 
               IF (QS (J) .LT.Q (J) ) THEN 
                  ALV = LV0 - CPVMCL * (T (J) - 273.15) 
                  TNEW = T (J) + ALV * (Q (J) - QS (J) ) / (CPD *       &
                  (1. - Q (J) ) + CL * Q (J) + QS (J) * (CPV - CL + ALV &
                  * ALV / (RV * T (J) * T (J) ) ) )                     
                  ALVNEW = LV0 - CPVMCL * (TNEW - 273.15) 
                  QNEW = (ALV * Q (J) - (TNEW - T (J) ) * (CPD *        &
                  (1. - Q (J) ) + CL * Q (J) ) ) / ALVNEW               
                  PRECIP = PRECIP + 24. * 3600. * 1.0E5 * (PH (J)       &
                  - PH (J + 1) ) * (Q (J) - QNEW) / (G * DELT * ROWL)   
                  T (J) = TNEW 
                  Q (J) = QNEW 
                  QS (J) = QNEW 
               ENDIF 
   38       END DO 
         ENDIF 
!                                                                       
      ENDIF 
!                                                                       
!  *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
!                                                                       
      GZ (1) = 0.0 
      CPN (1) = CPD * (1. - Q (1) ) + Q (1) * CPV 
      H (1) = T (1) * CPN (1) 
      LV (1) = LV0 - CPVMCL * (T (1) - 273.15) 
      HM (1) = LV (1) * Q (1) 
      TV (1) = T (1) * (1. + Q (1) * EPSI - Q (1) ) 
      AHMIN = 1.0E12 
      IHMIN = NL 
      DO 40 I = 2, NL + 1 
         TVX = T (I) * (1. + Q (I) * EPSI - Q (I) ) 
         TVY = T (I - 1) * (1. + Q (I - 1) * EPSI - Q (I - 1) ) 
         GZ (I) = GZ (I - 1) + 0.5 * RD * (TVX + TVY) * (P (I - 1)      &
         - P (I) ) / PH (I)                                             
         CPN (I) = CPD * (1. - Q (I) ) + CPV * Q (I) 
         H (I) = T (I) * CPN (I) + GZ (I) 
         LV (I) = LV0 - CPVMCL * (T (I) - 273.15) 
         HM (I) = (CPD * (1. - Q (I) ) + CL * Q (I) ) * (T (I) - T (1) )&
         + LV (I) * Q (I) + GZ (I)                                      
         TV (I) = T (I) * (1. + Q (I) * EPSI - Q (I) ) 
!                                                                       
!  ***  Find level of minimum moist static energy    ***                
!                                                                       
         IF (I.GE.MINORIG.AND.HM (I) .LT.AHMIN.AND.HM (I) .LT.HM (I - 1)&
         ) THEN                                                         
            AHMIN = HM (I) 
            IHMIN = I 
         ENDIF 
   40 END DO 
      IHMIN = MIN (IHMIN, NL - 1) 
!                                                                       
!  ***     Find that model level below the level of minimum moist       
!  ***  static energy that has the maximum value of moist static energy 
!                                                                       
      AHMAX = 0.0 
      NK = MINORIG ! added by juanxiong he
      DO 42 I = MINORIG, IHMIN 
         IF (HM (I) .GT.AHMAX) THEN 
            NK = I 
            AHMAX = HM (I) 
         ENDIF 
   42 END DO 
!                                                                       
!  ***  CHECK WHETHER PARCEL LEVEL TEMPERATURE AND SPECIFIC HUMIDITY   *
!  ***                          ARE REASONABLE                         *
!  ***      Skip convection if HM increases monotonically upward       *
!                                                                       
      IF (T (NK) .LT.250.0.OR.Q (NK) .LE.0.0.OR.IHMIN.EQ. (NL - 1) )    &
      THEN                                                              
         IFLAG = 0 
         CBMF = 0.0 
         RETURN 
      ENDIF 
!                                                                       
!   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT PARCEL ORIGIN LEV
!   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)   
!                                                                       
      RH = Q (NK) / QS (NK) 
      CHI = T (NK) / (1669.0 - 122.0 * RH - T (NK) ) 
      PLCL = P (NK) * (RH**CHI) 
      IF (PLCL.LT.200.0.OR.PLCL.GE.2000.0) THEN 
         IFLAG = 2 
         CBMF = 0.0 
         RETURN 
      ENDIF 
!                                                                       
!   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***                    
!                                                                       
      ICB = NL - 1 
      DO 50 I = NK + 1, NL 
         IF (P (I) .LT.PLCL) THEN 
            ICB = MIN (ICB, I) 
         ENDIF 
   50 END DO 
      IF (ICB.GE. (NL - 1) ) THEN 
         IFLAG = 3 
         CBMF = 0.0 
         RETURN 
      ENDIF 
!                                                                       
!   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY        
!                                                                       
!   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL   
!   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC          
!   ***                   LIQUID WATER CONTENT                          
!                                                                       
      CALL TLIFT (P, T, Q, QS, GZ, ICB, NK, TVP, TP, CLW, ND, NL, 1) 
      DO 54 I = NK, ICB 
         TVP (I) = TVP (I) - TP (I) * Q (NK) 
   54 END DO 
!                                                                       
!   ***  If there was no convection at last time step and parcel    *** 
!   ***       is stable at ICB then skip rest of calculation        *** 
!                                                                       
      IF (CBMF.EQ.0.0.AND.TVP (ICB) .LE. (TV (ICB) - DTMAX) ) THEN 
         IFLAG = 0 
         RETURN 
      ENDIF 
!                                                                       
!   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESS
!                                                                       
      IF (IFLAG.NE.4) IFLAG = 1 
!                                                                       
!   ***  FIND THE REST OF THE LIFTED PARCEL TEMPERATURES          ***   
!                                                                       
      CALL TLIFT (P, T, Q, QS, GZ, ICB, NK, TVP, TP, CLW, ND, NL, 2) 
!                                                                       
!   ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***   
!   ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***   
!   ***      THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)     ***   
!                                                                       
      DO 57 I = 1, NK 
         EP (I) = 0.0 
         SIGP (I) = SIGS 
   57 END DO 
      DO 60 I = NK + 1, NL 
         TCA = TP (I) - 273.15 
         IF (TCA.GE.0.0) THEN 
            ELACRIT = ELCRIT 
         ELSE 
            ELACRIT = ELCRIT * (1.0 - TCA / TLCRIT) 
         ENDIF 
         ELACRIT = MAX (ELACRIT, 0.0) 
         EPMAX = 0.999 
         EP (I) = EPMAX * (1.0 - ELACRIT / MAX (CLW (I), 1.0E-8) ) 
         EP (I) = MAX (EP (I), 0.0) 
         EP (I) = MIN (EP (I), EPMAX) 
         SIGP (I) = SIGS 
   60 END DO 
!                                                                       
!   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***   
!   ***                    VIRTUAL TEMPERATURE                    ***   
!                                                                       
      DO 64 I = ICB + 1, NL 
         TVP (I) = TVP (I) - TP (I) * Q (NK) 
   64 END DO 
      TVP (NL + 1) = TVP (NL) - (GZ (NL + 1) - GZ (NL) ) / CPD 
!                                                                       
!   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS   
!                                                                       
      DO 70 I = 1, NL + 1 
         HP (I) = H (I) 
         NENT (I) = 0 
         WATER (I) = 0.0 
         EVAP (I) = 0.0 
         WT (I) = OMTSNOW 
         MP (I) = 0.0 
         M (I) = 0.0 
         LVCP (I) = LV (I) / CPN (I) 
         DO 70 J = 1, NL + 1 
            QENT (I, J) = Q (J) 
            ELIJ (I, J) = 0.0 
            MENT (I, J) = 0.0 
            SIJ (I, J) = 0.0 
            UENT (I, J) = U (J) 
            VENT (I, J) = V (J) 
            DO 70 K = 1, NTRA 
               TRAENT (I, J, K) = TRA (J, K) 
   70 CONTINUE 
      QP (1) = Q (1) 
      UP (1) = U (1) 
      VP (1) = V (1) 
      DO 71 I = 1, NTRA 
         TRAP(1,I)=TRA(1,I)                                             
   71 CONTINUE                                                          
        DO 72 I=2,NL+1                                                  
         QP(I)=Q(I-1)                                                   
         UP(I)=U(I-1)                                                   
         VP(I)=V(I-1)                                                   
         DO 72 J=1,NTRA                                                 
          TRAP(I,J)=TRA(I-1,J)                                          
   72   CONTINUE                                                       
!                                                                       
!  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***   
!  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***   
!  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***   
!                                                                       
        CAPE=0.0                                                        
        CAPEM=0.0                                                       
        INB=ICB+1                                                       
        INB1=INB                                                        
        BYP=0.0                                                         
        DO 82 I=ICB+1,NL-1                                              
         BY=(TVP(I)-TV(I))*(PH(I)-PH(I+1))/P(I)                         
         CAPE=CAPE+BY                                                   
         IF(BY.GE.0.0)INB1=I+1                                          
         IF(CAPE.GT.0.0)THEN                                            
          INB=I+1                                                       
          BYP=(TVP(I+1)-TV(I+1))*(PH(I+1)-PH(I+2))/P(I+1)               
          CAPEM=CAPE                                                    
         END IF                                                         
   82 	CONTINUE                                                         
        INB=MAX(INB,INB1)                                               
        CAPE=CAPEM+BYP                                                  
        DEFRAC=CAPEM-CAPE                                               
        DEFRAC=MAX(DEFRAC,0.001)                                        
        FRAC=-CAPE/DEFRAC                                               
        FRAC=MIN(FRAC,1.0)                                              
        FRAC=MAX(FRAC,0.0)                                              
!                                                                       
!   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***   
!                                                                       
        DO 95 I=ICB,INB                                                 
         HP(I)=H(NK)+(LV(I)+(CPD-CPV)*T(I))*EP(I)*CLW(I)                
   95   CONTINUE                                                        
!                                                                       
!   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  *** 
!   ***                   AT EACH MODEL LEVEL                       *** 
!                                                                       
        DBOSUM=0.0                                                      
!                                                                       
!   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***   
!   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***   
!	                                                                      
        TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(P(ICB-1)-PLCL)/   &            
         (CPN(ICB-1)*P(ICB-1))                                         
        TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-P(ICB))/ &           
         (P(ICB)-P(ICB+1))                                             
        DTPBL=0.0                                                       
        DO 96 I=NK,ICB-1                                                
         DTPBL=DTPBL+(TVP(I)-TV(I))*(PH(I)-PH(I+1))                     
   96   CONTINUE                                                        
        DTPBL=DTPBL/(PH(NK)-PH(ICB))                                    
        DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL                               
        DTMA=DTMIN                                                      
!                                                                       
!   ***  ADJUST CLOUD BASE MASS FLUX   ***                              
!                                                                       
      CBMFOLD=CBMF                                                      
!   ***  by chenhs ***                                                  
      DELT0=600.0                                                       
      DAMPS=DAMP*DELT/DELT0                                             
      CBMF=(1.-DAMPS)*CBMF+0.1*ALPHA*DTMA                               
      CBMF=MAX(CBMF,0.0)                                                
!                                                                       
!   *** If cloud base mass flux is zero, skip rest of calculation  ***  
!                                                                       
      IF(CBMF.EQ.0.0.AND.CBMFOLD.EQ.0.0)THEN                            
       RETURN                                                           
      END IF                                                            
!                                                                       
!   ***   CALCULATE RATES OF MIXING,  M(I)   ***                        
!                                                                       
      M(ICB)=0.0                                                        
      DO 103 I=ICB+1,INB                                                
       K=MIN(I,INB1)                                                    
       DBO=ABS(TV(K)-TVP(K))+ &                           
       ENTP*0.02*(PH(K)-PH(K+1))                                       
       DBOSUM=DBOSUM+DBO                                                
       M(I)=CBMF*DBO                                                    
  103 CONTINUE                                                          
      DO 110 I=ICB+1,INB                                                
       M(I)=M(I)/DBOSUM                                                 
  110 CONTINUE                                                          
!                                                                       
!   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  *
!   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     *
!   ***                        FRACTION (SIJ)                          *
!                                                                       
        DO 170 I=ICB+1,INB                                              
         QTI=Q(NK)-EP(I)*CLW(I)                                         
         DO 160 J=ICB,INB                                               
          BF2=1.+LV(J)*LV(J)*QS(J)/(RV*T(J)*T(J)*CPD)                   
          ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(QTI-Q(J))                     
          DENOM=H(I)-HP(I)+(CPD-CPV)*(Q(I)-QTI)*T(J)                    
          DEI=DENOM                                                     
          IF(ABS(DEI).LT.0.01)DEI=0.01                                  
          SIJ(I,J)=ANUM/DEI                                             
          SIJ(I,I)=1.0                                                  
          ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)                   
          ALTEM=ALTEM/BF2                                               
          CWAT=CLW(J)*(1.-EP(J))                                        
          STEMP=SIJ(I,J)                                                
          IF((STEMP.LT.0.0.OR.STEMP.GT.1.0.OR.    &                      
           ALTEM.GT.CWAT).AND.J.GT.I)THEN                              
           ANUM=ANUM-LV(J)*(QTI-QS(J)-CWAT*BF2)                         
           DENOM=DENOM+LV(J)*(Q(I)-QTI)                                 
           IF(ABS(DENOM).LT.0.01)DENOM=0.01                             
           SIJ(I,J)=ANUM/DENOM                                          
           ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)                  
           ALTEM=ALTEM-(BF2-1.)*CWAT                                    
          END IF                                                        
          IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN                   
           QENT(I,J)=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI                    
           UENT(I,J)=SIJ(I,J)*U(I)+(1.-SIJ(I,J))*U(NK)                  
           VENT(I,J)=SIJ(I,J)*V(I)+(1.-SIJ(I,J))*V(NK)                  
           DO K=1,NTRA                                                  
            TRAENT(I,J,K)=SIJ(I,J)*TRA(I,K)+(1.-SIJ(I,J))* & 
            TRA(NK,K)                                                  
           END DO                                                       
           ELIJ(I,J)=ALTEM                                              
           ELIJ(I,J)=MAX(0.0,ELIJ(I,J))                                 
           MENT(I,J)=M(I)/(1.-SIJ(I,J))                                 
           NENT(I)=NENT(I)+1                                            
          END IF                                                        
          SIJ(I,J)=MAX(0.0,SIJ(I,J))                                    
          SIJ(I,J)=MIN(1.0,SIJ(I,J))                                    
  160    CONTINUE                                                       
!                                                                       
!   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS 
!   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES 
!                                                                       
         IF(NENT(I).EQ.0)THEN                                           
          MENT(I,I)=M(I)                                                
          QENT(I,I)=Q(NK)-EP(I)*CLW(I)                                  
          UENT(I,I)=U(NK)                                               
          VENT(I,I)=V(NK)                                               
          DO J=1,NTRA                                                   
           TRAENT(I,I,J)=TRA(NK,J)                                      
          END DO                                                        
          ELIJ(I,I)=CLW(I)                                              
          SIJ(I,I)=1.0                                                  
         END IF                                                         
  170   CONTINUE                                                        
        SIJ(INB,INB)=1.0                                                
!                                                                       
!   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***    
!   ***              PROBABILITIES OF MIXING                     ***    
!                                                                       
        DO 200 I=ICB+1,INB                                              
        IF(NENT(I).NE.0)THEN                                            
         QP1=Q(NK)-EP(I)*CLW(I)                                         
         ANUM=H(I)-HP(I)-LV(I)*(QP1-QS(I))                              
         DENOM=H(I)-HP(I)+LV(I)*(Q(I)-QP1)                              
         IF(ABS(DENOM).LT.0.01)DENOM=0.01                               
         SCRIT=ANUM/DENOM                                               
         ALT=QP1-QS(I)+SCRIT*(Q(I)-QP1)                                 
         IF(ALT.LT.0.0)SCRIT=1.0                                        
          SCRIT=MAX(SCRIT,0.0)                                          
         ASIJ=0.0                                                       
         SMIN=1.0                                                       
         DO 175 J=ICB,INB                                               
          IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN                   
           IF(J.GT.I)THEN                                               
            SMID=MIN(SIJ(I,J),SCRIT)                                    
            SJMAX=SMID                                                  
            SJMIN=SMID                                                  
            IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN                 
             SMIN=SMID                                                  
             SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)                       
             SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))                             
             SJMIN=MIN(SJMIN,SCRIT)                                     
            END IF                                                      
           ELSE                                                         
            SJMAX=MAX(SIJ(I,J+1),SCRIT)                                 
            SMID=MAX(SIJ(I,J),SCRIT)                                    
            SJMIN=0.0                                                   
            IF(J.GT.1)SJMIN=SIJ(I,J-1)                                  
            SJMIN=MAX(SJMIN,SCRIT)                                      
           END IF                                                       
           DELP=ABS(SJMAX-SMID)                                         
           DELM=ABS(SJMIN-SMID)                                         
           ASIJ=ASIJ+(DELP+DELM)*(PH(J)-PH(J+1))                        
           MENT(I,J)=MENT(I,J)*(DELP+DELM)*(PH(J)-PH(J+1))              
          END IF                                                        
  175    CONTINUE                                                       
         ASIJ=MAX(1.0E-21,ASIJ)                                         
         ASIJ=1.0/ASIJ                                                  
         DO 180 J=ICB,INB                                               
          MENT(I,J)=MENT(I,J)*ASIJ                                      
  180    CONTINUE                                                       
         BSUM=0.0                                                       
         DO 190 J=ICB,INB                                               
          BSUM=BSUM+MENT(I,J)                                           
  190    CONTINUE                                                       
         IF(BSUM.LT.1.0E-18)THEN                                        
          NENT(I)=0                                                     
          MENT(I,I)=M(I)                                                
          QENT(I,I)=Q(NK)-EP(I)*CLW(I)                                  
          UENT(I,I)=U(NK)                                               
          VENT(I,I)=V(NK)                                               
          DO J=1,NTRA                                                   
           TRAENT(I,I,J)=TRA(NK,J)                                      
          END DO                                                        
          ELIJ(I,I)=CLW(I)                                              
          SIJ(I,I)=1.0                                                  
         END IF                                                         
        END IF                                                          
  200   CONTINUE                                                        
!                                                                       
!   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***      
!   ***             DOWNDRAFT CALCULATION                      ***      
!                                                                       
        IF(EP(INB).LT.0.0001)GOTO 405                                   
!                                                                       
!   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***  
!   ***                AND CONDENSED WATER FLUX                    ***  
!                                                                       
        JTT=2                                                           
!                                                                       
!    ***                    BEGIN DOWNDRAFT LOOP                    *** 
!                                                                       
        DO 400 I=INB,1,-1                                               
!                                                                       
!    ***              CALCULATE DETRAINED PRECIPITATION             *** 
!                                                                       
        WDTRAIN=G*EP(I)*M(I)*CLW(I)                                     
        IF(I.GT.1)THEN                                                  
         DO 320 J=1,I-1                                                 
         AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)                               
         AWAT=MAX(0.0,AWAT)                                             
  320    WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)                               
        END IF                                                          
!                                                                       
!    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***     
!    ***              ESTIMATES OF QP(I)AND QP(I-1)             ***     
!                                                                       
!                                                                       
!  ***  Value of terminal velocity and coefficient of evaporation for sn
!                                                                       
        COEFF=COEFFS                                                    
        WT(I)=OMTSNOW                                                   
!                                                                       
!  ***  Value of terminal velocity and coefficient of evaporation for ra
!                                                                       
        IF(T(I).GT.273.0)THEN                                           
         COEFF=COEFFR                                                   
         WT(I)=OMTRAIN                                                  
        END IF                                                          
        QSM=0.5*(Q(I)+QP(I+1))                                          
        AFAC=COEFF*PH(I)*(QS(I)-QSM)/(1.0E4+2.0E3*PH(I)*QS(I))          
        AFAC=MAX(AFAC,0.0)                                              
        SIGT=SIGP(I)                                                    
        SIGT=MAX(0.0,SIGT)                                              
        SIGT=MIN(1.0,SIGT)                                              
        B6=100.*(PH(I)-PH(I+1))*SIGT*AFAC/WT(I)                         
        C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)                      
        REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))                               
        EVAP(I)=SIGT*AFAC*REVAP                                         
        WATER(I)=REVAP*REVAP                                            
!                                                                       
!    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***     
!    ***              HYDROSTATIC APPROXIMATION                 ***     
!                                                                       
        IF(I.EQ.1)GOTO 360                                              
        DHDP=(H(I)-H(I-1))/(P(I-1)-P(I))                                
        DHDP=MAX(DHDP,10.0)                                             
        MP(I)=100.*GINV*LV(I)*SIGD*EVAP(I)/DHDP                         
        MP(I)=MAX(MP(I),0.0)                                            
!                                                                       
!   ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***     
!                                                                       
        FAC=20.0/(PH(I-1)-PH(I))                                        
        MP(I)=(FAC*MP(I+1)+MP(I))/(1.+FAC)                              
!                                                                       
!    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 *** 
!    ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  *** 
!                                                                       
          IF(P(I).GT.(0.949*P(1)))THEN                                  
           JTT=MAX(JTT,I)                                               
           MP(I)=MP(JTT)*(P(1)-P(I))/(P(1)-P(JTT))                      
          END IF                                                        
  360   CONTINUE                                                        
!                                                                       
!    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***     
!                                                                       
        IF(I.EQ.INB)GOTO 400                                            
        IF(I.EQ.1)THEN                                                  
         QSTM=QS(1)                                                     
        ELSE                                                            
         QSTM=QS(I-1)                                                   
        END IF                                                          
        IF(MP(I).GT.MP(I+1))THEN                                        
          RAT=MP(I+1)/MP(I)                                             
          QP(I)=QP(I+1)*RAT+Q(I)*(1.0-RAT)+100.*GINV*  &                 
            SIGD*(PH(I)-PH(I+1))*(EVAP(I)/MP(I))                       
          UP(I)=UP(I+1)*RAT+U(I)*(1.-RAT)                               
          VP(I)=VP(I+1)*RAT+V(I)*(1.-RAT)                               
          DO J=1,NTRA                                                   
           TRAP(I,J)=TRAP(I+1,J)*RAT+TRAP(I,J)*(1.-RAT)                 
          END DO                                                        
         ELSE                                                           
          IF(MP(I+1).GT.0.0)THEN                                        
            QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+T(I+1)*(  &            
             CL-CPD))+CPD*(T(I+1)-T(I)))/(LV(I)+T(I)*(CL-CPD))         
            UP(I)=UP(I+1)                                               
            VP(I)=VP(I+1)                                               
            DO J=1,NTRA                                                 
             TRAP(I,J)=TRAP(I+1,J)                                      
            END DO                                                      
          END IF                                                        
        END IF                                                          
        QP(I)=MIN(QP(I),QSTM)                                           
        QP(I)=MAX(QP(I),0.0)                                            
  400   CONTINUE                                                        
!                                                                       
!   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***              
!                                                                       
        PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600.*24000./(ROWL*G)         
!                                                                       
  405   CONTINUE                                                        
!                                                                       
!   ***  CALCULATE DOWNDRAFT VELOCITY SCALE AND SURFACE TEMPERATURE AND 
!   ***                    WATER VAPOR FLUCTUATIONS                     
!                                                                       
      WD=BETA*ABS(MP(ICB))*0.01*RD*T(ICB)/(SIGD*P(ICB))                 
      QPRIME=0.5*(QP(1)-Q(1))                                           
      TPRIME=LV0*QPRIME/CPD                                             
!                                                                       
!   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
!   ***                      AND MIXING RATIO                        ***
!                                                                       
        DPINV=0.01/(PH(1)-PH(2))                                        
        AM=0.0                                                          
        IF(NK.EQ.1)THEN                                                 
         DO 410 K=2,INB                                                 
  410    AM=AM+M(K)                                                     
        END IF                                                          
        IF((2.*G*DPINV*AM).GE.DELTI)IFLAG=4                             
        FT(1)=FT(1)+G*DPINV*AM*(T(2)-T(1)+(GZ(2)-GZ(1))/CPN(1))         
        FT(1)=FT(1)-LVCP(1)*SIGD*EVAP(1)                                
        FT(1)=FT(1)+SIGD*WT(2)*(CL-CPD)*WATER(2)*(T(2)- &                
        T(1))*DPINV/CPN(1)                                             
        FQ(1)=FQ(1)+G*MP(2)*(QP(2)-Q(1))* &                              
         DPINV+SIGD*EVAP(1)                                            
        FQ(1)=FQ(1)+G*AM*(Q(2)-Q(1))*DPINV                              
        FU(1)=FU(1)+G*DPINV*(MP(2)*(UP(2)-U(1))+AM*(U(2)-U(1)))         
        FV(1)=FV(1)+G*DPINV*(MP(2)*(VP(2)-V(1))+AM*(V(2)-V(1)))         
        DO J=1,NTRA                                                     
         FTRA(1,J)=FTRA(1,J)+G*DPINV*(MP(2)*(TRAP(2,J)-TRA(1,J))+  &     
         AM*(TRA(2,J)-TRA(1,J)))                                       
         FTRAD(1,J)=FTRAD(1,J)+G*DPINV*(MP(2)*(TRAP(2,J))+ &             
         AM*(TRA(2,J)))                                                
         FTRAO(1,J)=FTRAO(1,J)+G*DPINV*(MP(2)*(-TRA(1,J))+ &             
         AM*(-TRA(1,J)))                                               
        END DO                                                          
        AMDE=0.0                                                        
        DO 415 J=2,INB                                                  
         FQ(1)=FQ(1)+G*DPINV*MENT(J,1)*(QENT(J,1)-Q(1))                 
         FU(1)=FU(1)+G*DPINV*MENT(J,1)*(UENT(J,1)-U(1))                 
         FV(1)=FV(1)+G*DPINV*MENT(J,1)*(VENT(J,1)-V(1))                 
         DO K=1,NTRA                                                    
          FTRA(1,K)=FTRA(1,K)+G*DPINV*MENT(J,1)*(TRAENT(J,1,K)- &        
          TRA(1,K))                                                    
          FTRAE(1,K)=FTRAE(1,K)+G*DPINV*MENT(J,1)*(TRAENT(J,1,K)- &       
          TRA(1,K))                                                    
                                                                        
         END DO                                                         
  415   CONTINUE                                                        
!                                                                       
!   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO 
!   ***               AT LEVELS ABOVE THE LOWEST LEVEL                  
!                                                                       
!   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES 
!   ***                      THROUGH EACH LEVEL                         
!                                                                       
        DO 500 I=2,INB                                                  
        DPINV=0.01/(PH(I)-PH(I+1))                                      
        CPINV=1.0/CPN(I)                                                
        AMP1=0.0                                                        
        AD=0.0                                                          
        IF(I.GE.NK)THEN                                                 
         DO 440 K=I+1,INB+1                                             
  440    AMP1=AMP1+M(K)                                                 
        END IF                                                          
        DO 450 K=1,I                                                    
        DO 450 J=I+1,INB+1                                              
         AMP1=AMP1+MENT(K,J)                                            
  450   CONTINUE                                                        
        IF((2.*G*DPINV*AMP1).GE.DELTI)IFLAG=4                           
        DO 470 K=1,I-1                                                  
        DO 470 J=I,INB                                                  
         AD=AD+MENT(J,K)                                                
  470   CONTINUE                                                        
        FT(I)=FT(I)+G*DPINV*(AMP1*(T(I+1)-T(I)+(GZ(I+1)-GZ(I))* &        
        CPINV)-AD*(T(I)-T(I-1)+(GZ(I)-GZ(I-1))*CPINV))         &        
        -SIGD*LVCP(I)*EVAP(I)                                          
        FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+  &                    
         T(I)*(CPV-CPD)*(Q(I)-QENT(I,I)))*CPINV                        
        FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)* &                   
         (T(I+1)-T(I))*DPINV*CPINV                                     
        FQ(I)=FQ(I)+G*DPINV*(AMP1*(Q(I+1)-Q(I))-     &                   
         AD*(Q(I)-Q(I-1)))                                             
        FU(I)=FU(I)+G*DPINV*(AMP1*(U(I+1)-U(I))-      &                  
         AD*(U(I)-U(I-1)))                                             
        FV(I)=FV(I)+G*DPINV*(AMP1*(V(I+1)-V(I))-       &                 
         AD*(V(I)-V(I-1)))                                             
        DO K=1,NTRA                                                     
         FTRA(I,K)=FTRA(I,K)+G*DPINV*(AMP1*(TRA(I+1,K)- &                
         TRA(I,K))-AD*(TRA(I,K)-TRA(I-1,K)))                           
         FTRAD(I,K)=FTRAD(I,K)+G*DPINV*(AMP1*(TRA(I+1,K)))              
         FTRAU(I,K)=FTRAU(I,K)+G*DPINV*(-AD*(-TRA(I-1,K)))              
         FTRAO(I,K)=FTRAO(I,K)+G*DPINV*(AMP1*(-     &
         TRA(I,K))-AD*(TRA(I,K)))                                      
        END DO                                                          
        DO 480 K=1,I-1                                                  
         AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)                               
         AWAT=MAX(AWAT,0.0)                                             
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-Q(I))            
         FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))                 
         FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))                 
         DO J=1,NTRA                                                    
          FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-   &
          TRA(I,J))                                                    
          FTRAE(I,J)=FTRAE(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-  &     
          TRA(I,J))                                                    
         END DO                                                         
  480   CONTINUE                                                        
        DO 490 K=I,INB                                                  
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-Q(I))                 
         FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))                 
         FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))                 
         DO J=1,NTRA                                                    
          FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)- &        
          TRA(I,J))                                                    
           FTRAE(I,J)=FTRAE(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)- &      
          TRA(I,J))                                                    
         END DO                                                         
  490   CONTINUE                                                        
        FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)*       &
         (QP(I+1)-Q(I))-MP(I)*(QP(I)-Q(I-1)))*DPINV                    
        FU(I)=FU(I)+G*(MP(I+1)*(UP(I+1)-U(I))-MP(I)* &                   
         (UP(I)-U(I-1)))*DPINV                                         
        FV(I)=FV(I)+G*(MP(I+1)*(VP(I+1)-V(I))-MP(I)*  &                  
         (VP(I)-V(I-1)))*DPINV                                         
        DO J=1,NTRA                                                     
         FTRA(I,J)=FTRA(I,J)+G*DPINV*(MP(I+1)*(TRAP(I+1,J)-TRA(I,J))- &  
         MP(I)*(TRAP(I,J)-TRA(I-1,J)))                                 
         FTRAD(I,J)=FTRAD(I,J)+G*DPINV*(MP(I+1)*TRAP(I+1,J))            
         FTRAU(I,J)=FTRAU(I,J)+G*DPINV*(-MP(I)*(-TRA(I-1,J)))           
         FTRAO(I,J)=FTRAO(I,J)+G*DPINV*(MP(I+1)*(-TRA(I,J))-  & 
         MP(I)*(TRAP(I,J)))                                            
        END DO                                                          
  500   CONTINUE                                                        
!                                                                       
!   *** Adjust tendencies at top of convection layer to reflect  ***    
!   ***       actual position of the level zero CAPE             ***    
!                                                                       
        FQOLD=FQ(INB)                                                   
        FQ(INB)=FQ(INB)*(1.-FRAC)                                       
        FQ(INB-1)=FQ(INB-1)+FRAC*FQOLD*((PH(INB)-PH(INB+1))/ &
        (PH(INB-1)-PH(INB)))*LV(INB)/LV(INB-1)                         
        FTOLD=FT(INB)                                                   
        FT(INB)=FT(INB)*(1.-FRAC)                                       
        FT(INB-1)=FT(INB-1)+FRAC*FTOLD*((PH(INB)-PH(INB+1))/ &           
        (PH(INB-1)-PH(INB)))*CPN(INB)/CPN(INB-1)                       
        FUOLD=FU(INB)                                                   
        FU(INB)=FU(INB)*(1.-FRAC)                                       
        FU(INB-1)=FU(INB-1)+FRAC*FUOLD*((PH(INB)-PH(INB+1))/ &           
        (PH(INB-1)-PH(INB)))                                           
        FVOLD=FV(INB)                                                   
        FV(INB)=FV(INB)*(1.-FRAC)                                       
        FV(INB-1)=FV(INB-1)+FRAC*FVOLD*((PH(INB)-PH(INB+1))/ &           
        (PH(INB-1)-PH(INB)))                                           
        DO K=1,NTRA                                                     
         FTRAOLD=FTRA(INB,K)                                            
         FTRA(INB,K)=FTRA(INB,K)*(1.-FRAC)                              
         FTRA(INB-1,K)=FTRA(INB-1,K)+FRAC*FTRAOLD*(PH(INB)-PH(INB+1))/ & 
         (PH(INB-1)-PH(INB))                                           
         FTRAD(INB,K)=FTRAD(INB,K)*(1.-FRAC)                            
         FTRAU(INB,K)=FTRAU(INB,K)*(1.-FRAC)                            
         FTRAOLDD=FTRAD(INB,K)                                          
         FTRAOLDU=FTRAU(INB,K)                                          
        FTRAD(INB-1,K)=FTRAD(INB-1,K)+FRAC*FTRAOLDD*(PH(INB)-PH(INB+1))/ &
         (PH(INB-1)-PH(INB))                                           
        FTRAU(INB-1,K)=FTRAU(INB-1,K)+FRAC*FTRAOLDU*(PH(INB)-PH(INB+1))/ &
         (PH(INB-1)-PH(INB))                                           
                                                                        
        END DO                                                          
!                                                                       
!   ***   Very slightly adjust tendencies to force exact   ***          
!   ***     enthalpy, momentum and tracer conservation     ***          
!                                                                       
        ENTS=0.0                                                        
        UAV=0.0                                                         
        VAV=0.0                                                         
        DO 680 I=1,INB                                                  
         ENTS=ENTS+(CPN(I)*FT(I)+LV(I)*FQ(I))*(PH(I)-PH(I+1))	          
         UAV=UAV+FU(I)*(PH(I)-PH(I+1))                                  
         VAV=VAV+FV(I)*(PH(I)-PH(I+1))                                  
  680 CONTINUE                                                          
        ENTS=ENTS/(PH(1)-PH(INB+1))                                     
        UAV=UAV/(PH(1)-PH(INB+1))                                       
        VAV=VAV/(PH(1)-PH(INB+1))                                       
        DO 640 I=1,INB                                                  
         FT(I)=FT(I)-ENTS/CPN(I)                                        
         FU(I)=(1.-CU)*(FU(I)-UAV)                                      
         FV(I)=(1.-CU)*(FV(I)-VAV)                                      
  640 CONTINUE                                                          
        DO 700 K=1,NTRA                                                 
         TRAAV=0.0                                                      
         DO 690 I=1,INB                                                 
          TRAAV=TRAAV+FTRA(I,K)*(PH(I)-PH(I+1))                         
  690    CONTINUE                                                       
         TRAAV=TRAAV/(PH(1)-PH(INB+1))                                  
         DO 695 I=1,INB                                                 
          FTRA(I,K)=FTRA(I,K)-TRAAV                                     
  695    CONTINUE                                                       
  700 CONTINUE                                                          
                                                                        
!   ***          Li Jie to exclude unstable   ****                      
        DO 699 K=1,NTRA                                                 
         DO 701 I=1,INB                                                 
          IF((TRA(I,K)-FTRA(I,K)*DELT).LE.0.OR.  &                       
            (TRA(I,K)-FTRA(I,K)*DELT).GE.TRA(I,K)*5.) THEN             
            FTRA  = 0.0                                                 
            FTRAD = 0.0                                                 
            FTRAU = 0.0                                                 
            FTRAO = 0.0                                                 
            FTRAE = 0.0                                                 
            GOTO 702                                                    
          ENDIF                                                         
  701    CONTINUE                                                       
  699    CONTINUE                                                       
  702    CONTINUE                                                       
!                                                                       
!   ***           RETURN           ***                                  
!                                                                       
        RETURN                                                          
!                                                                       
      END SUBROUTINE                                          
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
      SUBROUTINE TLIFT (P, T, Q, QS, GZ, ICB, NK, TVP, TPK, CLW, ND, NL,&
      KK)                                                               
      REAL GZ (ND), TPK (ND), CLW (ND), P (ND) 
      REAL T (ND), Q (ND), QS (ND), TVP (ND), LV0 
!                                                                       
!   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***              
!                                                                       
      CPD = 1005.7 
      CPV = 1870.0 
      CL = 2500.0 
      RV = 461.5 
      RD = 287.04 
      LV0 = 2.501E6 
!                                                                       
      CPVMCL = CL - CPV 
      EPS = RD / RV 
      EPSI = 1. / EPS 
!                                                                       
!   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   
!                                                                       
      AH0 = (CPD * (1. - Q (NK) ) + CL * Q (NK) ) * T (NK) + Q (NK)     &
      * (LV0 - CPVMCL * (T (NK) - 273.15) ) + GZ (NK)                   
      CPP = CPD * (1. - Q (NK) ) + Q (NK) * CPV 
      CPINV = 1. / CPP 
!                                                                       
      IF (KK.EQ.1) THEN 
!                                                                       
!   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***     
!                                                                       
         DO 50 I = 1, ICB - 1 
            CLW (I) = 0.0 
   50    END DO 
         DO 100 I = NK, ICB - 1 
            TPK (I) = T (NK) - (GZ (I) - GZ (NK) ) * CPINV 
            TVP (I) = TPK (I) * (1. + Q (NK) * EPSI) 
  100    END DO 
      ENDIF 
!                                                                       
!    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***         
!                                                                       
      NST = ICB 
      NSB = ICB 
      IF (KK.EQ.2) THEN 
         NST = NL 
         NSB = ICB + 1 
      ENDIF 
      DO 300 I = NSB, NST 
         TG = T (I) 
         QG = QS (I) 
         ALV = LV0 - CPVMCL * (T (I) - 273.15) 
         DO 200 J = 1, 2 
            S = CPD+ALV * ALV * QG / (RV * T (I) * T (I) ) 
            S = 1. / S 
            AHG = CPD * TG + (CL - CPD) * Q (NK) * T (I) + ALV * QG +   &
            GZ (I)                                                      
            TG = TG + S * (AH0 - AHG) 
            TG = MAX (TG, 35.0) 
            TC = TG - 273.15 
            DENOM = 243.5 + TC 
            IF (TC.GE.0.0) THEN 
               ES = 6.112 * EXP (17.67 * TC / DENOM) 
            ELSE 
               ES = EXP (23.33086 - 6111.72784 / TG + 0.15215 * LOG (TG)&
               )                                                        
            ENDIF 
            QG = EPS * ES / (P (I) - ES * (1. - EPS) ) 
  200    END DO 
         TPK (I) = (AH0 - (CL - CPD) * Q (NK) * T (I) - GZ (I) - ALV *  &
         QG) / CPD                                                      
         CLW (I) = Q (NK) - QG 
         CLW (I) = MAX (0.0, CLW (I) ) 
         RG = QG / (1. - Q (NK) ) 
         TVP (I) = TPK (I) * (1. + RG * EPSI) 
  300 END DO 
      RETURN 
      END SUBROUTINE TLIFT                          
end module convect43c      
