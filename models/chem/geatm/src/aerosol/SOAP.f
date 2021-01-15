C=====================================================================
C
C***  SUBROUTINE SOAP
C***  THIS SUBROUTINE IS A SAMPLE INTERFACE BETWEEN AN AIRSHED MODEL  
C***  AND A SECONDARY ORGANIC AEROSOL FORMATION MODULE.
C
C=====================================================================
C
        SUBROUTINE SOAP (NSOA,SOA,SVG,TEMP,PRESS,POA, II, JJ,KK)
C
C       NSOA   TOTAL NUMBER OF SECONDARY ORGANIC AEROSOLS                      
C       SVG    MIXING RATIOS OF SEMI-VILIATE GASES IN PPB
C          1   FROM ANTHTOPOGENIC TOL+OH AND XYL + OH
C          2   FROM ANTHTOPOGENIC TOL+OH AND XYL + OH
C          3   FROM BIOGENIC TERPB
C          4   FROM BIOGENIC TERPB
C          5   FROM BIOGENIC ISOP
C          6   FROM BIOGENIC ISOP        
C       SOA    CONCENTRATIONS OF SECONDARY ORGANIC AEROSOLS IN UG/M3
C          1   SOA FROM SVG1
C          2   SOA FROM SVG2
C          3   SOA FROM SVG3
C          4   SOA FROM SVG4
C          5   SOA FROM SVG5
C          6   SOA FROM SVG6
C       TOT    TOTAL CONCENTRATUION OF ith GAS AND AEROSOL IN UG/M3
C       CSATT  saturation concentrations of CG/SOA species(ug/m3)atTEMP        
C       TEMP   AIR TEMPERATURE IN K
C       POA    CONCENTRATION OF PRIMARY ORGANIC AEROSOL IN UG/M3
C       NWPOA  MOLECULAR WEIGHT OF POA  IN G/MOL
C       PRESS  PRESSURE IN HPA        
C
C***  LOCAL VARIABLE        
C        
C       NWSOA    molecular weights of CG/SOA species (g/mol)
C       CSAT     REFRENCE saturation concentrations of CG/SOA ug/m3)        
C       CSTEMP   REFRENCE TEMPERATURE
C       CSAT0    CSAT OF SVGI IN DIFERENT REACTIONS FOR CALCULATE CSAT                
C       DELTAH   enthalpy of vaporization of CG/SOA species (J/mol)
C       FLAGSOAP set to 1 if CG/SOA species forms solutions; 0 if not
C       ICONT    INDEXES
C       SUM      INDEXES
C       NSOL     TOTAL NUMBERS OF SOA SPECIES
C       SSOA     Concentrations of  SOA species (ug/m3)
c       SSVG     Concentrations of  GASES IN UG/M3        
C       SCSATT   saturation concentrations of solution-forming
C                SOA species (ug/m3)      
C       STOT     total concentrations of solution-forming SOA species
C               (ug/m3)   
C       SMW      MOLECULAR WEIGHT OF SOA  IN G/MOL        
C       ZNUM     counter for number of iterations
C       CONMIN   use simple solution for species below this level(ug/m3)        
C       POAMIN   no PRIMARY organics if POA < POAmin (ug/m3)        
C       XTOL     error tolerance for bi-section method
C       CONVFAC  CONVERTING FACTOR FROM PPB TO UG/M3        
        
        INTEGER NSOA,NSOL,ZNUM 
        REAL    TEMP,POA,NWPOA,SUM,CONVFAC
        REAL    SOA(NSOA),SVG(NSOA),TOT(NSOA),CSATT(NSOA)
        REAL    NWSOA(6), CSAT(6),CSTEMP(6),DELTAH(6)
        INTEGER FLAGSOAP(NSOA)
        REAL    SSOA(NSOA),SSVG(NSOA),SCSATT(NSOA),STOT(NSOA),SMW(NSOA)
        REAL    CONMIN,POAMIN,XTOL
        REAL    CSAT0(7,NSOA) ! 7 REACTION FOR PRODUCING SVG
C
        PARAMETER(CONMIN = 1.E-6)
        PARAMETER(POAMIN = 1.01E-9)
        PARAMETER(XTOL   = 5.E-5)        
C
        REAL    CPX,BB,CC,XEND,FEND,XMID,FMID,DX,IDX(NSOA)

C        DATA CSAT0 /0.053, 0.042, 0.,  0., 0., 0., 0.,          ! SVG1
C     &        0.0019,0.0014, 0.,  0., 0., 0., 0.,               ! SVG2
C     &        0., 0.,0.126386, 0.098145, 0.142002, 0.125645, 0. ! SVG3
C     &        0., 0.,0.004837, 0.003490, 0.004704, 0.004385, 0.,! SVG4
C     &        0., 0.,  0.,  0., 0., 0.,  0.00862,               ! SVG5
C     &        0., 0.,  0.,  0., 0., 0.,  1.62 /                 ! SVG6
C
        DATA CSAT  / 0.047,0.0016,0.123,0.00435,0.00862, 1.62/  ! AVERAGE 
        DATA CSTEMP/ 308., 308.,   308.,   308.,   295., 295./ 
        DATA DELTAH/  70.,  70.,    70.,    70.,    42.,  42./ 
        
        DATA NWSOA/ 136., 136., 168., 168., 130., 130. /
C        
        NWPOA = 220.
C        
C***   CONVERT GAS IN PPB TO UG/M3 ***********************************
        CONVFAC=44.9*(273./TEMP)*(PRESS/1013.)
C      
        DO 100 I =1,NSOA
         SVG (I) = SVG(I) * CONVFAC * NWSOA(I) / 1.E03
         TOT (I) = SOA(I) + SVG(I)
 100    CONTINUE        
         
        CPX = POA / NWPOA 

C***  CHANGE SATURATION CONCENTRATIONS ACCORDING TO CURRENT TEMPERATURE
        DO 101 I = 1,NSOA
         CSATT (I) = CSAT (I) * ( CSTEMP(I)/TEMP ) * 
     &          EXP( ( DELTAH(I)/8.314 ) * (1/CSTEMP(I)-1/TEMP)) 
 101    CONTINUE
C
C***  CALCULATE AEROSOL-PHASE CONCENTRATION (CAER) AND GAS-PHASE ****
C***  CONCENTRATION (CGAS) FOR NON-SOLUTION-FORMING COMPOUNDS    ****
C***  COMPOUNDS THAT HAVE A CONCENTRATION OF LESS THAN conmin ARE ***
C***  IGN0RED
         
C***  MAP COMPOUNDS THAT FORM SOLUTIONS ONTO ARRAYS ****************
         
       ICONT = 0
       DO 102 I =1, NSOA
        IF( TOT(I).LT. CONMIN ) THEN
          SOA(I) = 0.0
          SVG(I) = TOT (I)     
        ELSE
          ICONT = ICONT + 1      
          IDX(ICONT)   = I
          SMW(ICONT)   = NWSOA (I)
          SCSATT(ICONT)= CSATT(I)
          STOT(ICONT)  = TOT(I)
          SSOA(ICONT)  = SOA(I)  
        ENDIF   
 102   CONTINUE        

        NSOL = ICONT 
C
C***  CHECK FOR A TRIVIAL SOLUTION **********************************
C    
       IF(NSOL.EQ.0) GOTO 1000
C       
       IF(NSOL.EQ.1) THEN
        IF(POA.LT.POAMIN) THEN
         SSVG(1) = AMIN1(STOT(1),SCSATT(1))        
         SSOA(1) = STOT(1) - SSVG(1)
        ELSE 
         BB = SCSATT(1) - STOT(1) + CPX * SMW(1)
         CC = -STOT(1)*CPX*SMW(1) 
         SSOA(1) = AMIN1(STOT(1), 0.5*(-BB+SQRT(BB*BB-4.*CC)) )
         SSVG(1) = STOT(1) - SSOA(1)
        ENDIF 
        GOTO 900
       ENDIF 
C
C 
       
       SUM = 0.0
       DO 103 I = 1, NSOL
        SUM = SUM + STOT(I)/SCSATT(I)
 103   CONTINUE
C
C
       IF(POA.LT.POAMIN.AND.SUM.LE.1.0) THEN
         DO 104 I = 1 , NSOL       
          SSVG(I) = STOT(I)
          SSOA(I) = 0.0
 104     CONTINUE
       GOTO 900
      ENDIF
C
C *** FIND THE SOLUTION USING A BI-SECTION METHOD (APPROACH FROM MAX)       
C
      XEND = 0.0
      DO 105 I = 1, NSOL
       XEND = XEND + STOT(I)/SMW(I)    
 105  CONTINUE      
           
      XEND = XEND + CPX

      CALL  spfcn (NSOL, STOT, SCSATT, SSOA, SMW, CPX, XEND, FEND)

      IF( ABS(FEND).LE.XTOL*XEND ) GOTO 99
      IF( FEND.GT.0.0 ) THEN
       PRINT*,' ERROR in SOAP AND positive end point',II,JJ,KK
       GOTO 99
      ENDIF
C
      DX = XEND - CPX
C    
      DO 106 ZNUM = 1 , 200  
       DX = 0.5 * DX
       XMID = XEND - DX
       CALL spfcn (NSOL, STOT, SCSATT, SSOA, SMW, CPX, XMID, FMID)
       IF ( ABS(FMID).LE.XTOL*XMID.OR. DX.LE.XTOL*XMID ) GOTO 98
       IF (FMID.LT.0.) XEND = XMID
 106  CONTINUE
C
C
C
 99   XMID  = XEND 
 98   CONTINUE
C
C
      DO 107 I = 1, NSOL 
        SSOA(I) = AMIN1(STOT(I), SSOA(I) )
        SSVG(I) = STOT(I) - SSOA(I)
 107  CONTINUE
C
C
C
C***  REMAP COMPOUNDS THAT FORM SOLUTIONS BACK ONTO ORIGINAL ARRAYS
C
 900  CONTINUE
C
      DO 108 I = 1, NSOL
       SOA( IDX(I) ) = SSOA(I)
       SVG( IDX(I) ) = SSVG(I)       
 108  CONTINUE         
C
 1000 CONTINUE
C
C CONVERT TO PPBV FOR GASES ***************************************      
C
      DO 109 I = 1, NSOA
        SVG(I) = 1.E03*SVG(I)/ (CONVFAC*NWSOA(I))      
 109  CONTINUE
C    
C      IF(II==30.AND.JJ==30.AND.KK==1) 
C     &  PRINT*,(SVG(I),I=1,6),'111'
C      IF(II==30.AND.JJ==30.AND.KK==1)
C     &  PRINT*,(SOA(I),I=1,6),'222'   
      RETURN
      END


      subroutine spfcn (n,ct,cs,ca,mw,cpx,tom,fval)
      implicit none
c
c SPFCN calculates the objective function for the bi-section solver in
cSOAP
c  - revised for new objective function by bkoo (05/27/05)
c     Total Organics in Mole (TOM) = SUM_i(C_i(aer)/MW_i) + C_pre/MW_pre
c     C_i(aer) = C_i(tot) - x_i * Cstar_i
c              = C_i(tot) - (C_i(aer)/MW_i/TOM) * Cstar_i
c  => C_i(aer) = C_i(tot) * TOM / (TOM + Cstar_i/MW_i)
c  => SUM_i(C_i(tot) * TOM / (TOM*MW_i + Cstar_i)) + C_pre/MW_pre - TOM
c= 0
c
c Called by SOAP
c
      integer     n,i
      real        ct(n),cs(n),ca(n),mw(n),cpx,tom,fval
c
      fval = 0.0
      do i = 1, n
        ca(i) = ct(i) * tom / ( tom + cs(i) / mw(i) )
        fval  = fval + ca(i) / mw(i)
      enddo
      fval = fval + cpx - tom
c
      return
      end

      
