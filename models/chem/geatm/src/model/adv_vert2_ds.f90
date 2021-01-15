          SUBROUTINE ADV_VERT_DS(MYID,C,DSO4,DNO3, DFEIII,DFEIIW,U2D,V2D,DELTZ,&
                DELTX,SX,EX,SY,EY,NZZ,DT,IG)
!CCC   THIS ROUTINE IS TO CALCULATE THE ADVECTION OF DUST AND SEA SALT AEROSOLS AND AEROSOLS ccccc
!       COMPOSITIONS (ASO4, ANO3,FEII,FEIII) ASSUMING DUST (OR. SEA  SALT) INTERNAL MIXTURE cccccc

          INTEGER :: MYID,SX,EX,SY,EY,NZZ
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1,NZZ) :: C,W,U2D,V2D
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1,NZZ) :: DELTX,DELTZ
          REAL    :: DT,D0 
          REAL,DIMENSION(NZZ)                     :: Q0,QN,DXX0,DXX, &
                                                     DEN0,DEN1,U,DD0
          real, dimension(sx-1:ex+1,sy-1:ey+1,NZZ) :: dso4, dno3
          real, dimension(sx-1:ex+1,sy-1:ey+1,NZZ) ::  dfeIII, dfeII !so4, no3, FeIII and FeII on dust or sea salt particles  
          real,dimension(NZZ) :: DSSO4, DSNO3, DSFEII, DSFEIII

          DATA D0 /1.0/                                           
         
          DO J=SY,EY
           DO I=SX,EX
            DO  K=2,NZZ-1
             DXX0(K)  =  DELTX(I,J,K)
             DXX (K)  =  DELTZ(I,J,K)
             DEN0(K)  =  1.- DT/DXX0(K)*(D0*U2D(I,J,K)-D0*U2D(I-1,J,K))&
                           - DT/DXX0(K)*(D0*V2D(I,J,K)-D0*V2D(I,J-1,K))
             DEN1(K)  =  DEN0(K)-DT/DXX(K)*(D0*W(I,J,K)-D0*W(I,J,K-1)) 
             
             DD0 (K)  =  DEN0(K)
             Q0  (K)  =  C(I,J,K)
             U   (K)  =  W(I,J,K)
             DSSO4(K) = DSO4(i,j,k)
             DSNO3(K)  = DNO3(i,j,k)
             DSFEIII(K) = DFEIII (i,j,k)
             DSFEII(K) = DFEII(i,j,k)
            ENDDO  !K

             DD0 (1)  =  D0
             U   (1)  =  W(I,J,1)
             Q0  (1)  =  C(I,J,1)
             Q0  (NZZ)=  C(I,J,NZZ)

            CALL ADVEC1D1_DS(Q0,QN,,DSSO4, DSNO3, DSFEII, DSFEIII,U,DEN0,DEN1,DT,DXX,DD0,NZZ)

            DO K=2,NZZ-1
             C(I,J,K) = amax1(QN(K),0.0)
            ENDDO  !K
           ENDDO   !J
          ENDDO !I

          RETURN
          END


         SUBROUTINE ADVEC1D1_DS(Q0,QN,DSSO4, DSNO3, DSFEII, DSFEIII,U,DEN0,DEN1,DT,DXX,DD0,NZZ)
         INTEGER              :: NZZ
         REAL,DIMENSION(NZZ)  :: Q0,QN,U,DEN0,DEN1,DXX,DD0
         REAL                 :: DT,CK1,X1,FLUXI,FLUXIM1
         REAL,DIMENSION(NZZ)  :: FLUX,VCMAX,VCMIN,IMXMN
         REAL,DIMENSION(NZZ)  :: DSSO4, DSNO3, DSFEII,DSFEIII
         DATA ZR0,LSHP / 0.0,0.0/

         IMXMN = 0.0  
!!!!  INDENTIFY LOCAL MAX AND MIN,SPECIFY MIXING RATIO LIMITS AT NEW TIME
 
        DO 5 I=2,NZZ-1
         IMXMN(I) = 0
         IF(Q0(I).GE.MAX(Q0(I-1),Q0(I+1)).OR. &
            Q0(I).LE.MIN(Q0(I-1),Q0(I+1)) ) IMXMN(I)=1
         CK1 = Q0(I)
         IF(U(I).LT.ZR0)   CK1 = Q0(I+1)
         IF(U(I-1).GE.ZR0) CK1 =Q0(I-1)
         VCMAX(I) = MAX(Q0(I),CK1)
      5  VCMIN(I) = MIN(Q0(I),CK1)

!!!!    UPDATE MIXING RATIOS AND LIMITS FLUXES GOING UP WHERE U2D >0
             
        IF(U(1).GE.ZR0) FLUX(1) = Q0(1)*U(1)*DT*DD0(1)
        DO 10 I=2,NZZ-1
         IF(U(I).LT.ZR0) GOTO 10
         IF(U(I-1).LT.ZR0) THEN
            FLUX(I)=Q0(I)*U(I)*DT*DD0(I)              ! OUTFLOW-ONLY CELL    
         ELSE
            X1   = DT*U(I)/DXX(I)
            CF1  = CFCALC1_DS(X1,Q0(I),Q0(I-1),Q0(I+1),IMXMN(I-1),IMXMN(I+1),LSHP)
            FLUXI= DD0(I)*U(I)*DT*CF1/DXX(I)
            QN(I)= MAX(VCMIN(I) ,MIN(VCMAX(I), &
                  (Q0(I)*DEN0(I)-FLUXI+FLUX(I-1)/DXX(I))/DEN1(I) ) )

            FLUX(I) = DXX(I)*(Q0(I)*DEN0(I)-QN(I)*DEN1(I))+FLUX(I-1)  

        IF(FLUX(I).GT. 0. ) THEN

        DSSO4(I) = DSSO4(I) + &
                 ( FLUX(I-1)/DXX(I)/DEN1(I) * DSSO4(I-1)/Q0(I-1) -&
                   FLUX(I)/DXX(I)/DEN1(I) * DSSO4(I)/Q0(I)  )
        DSNO3(I) = DSNO3(I) + &
                 ( FLUX(I-1)/DXX(I)/DEN1(I) * DSNO3(I-1)/Q0(I-1)-&
                   FLUX(I)/DXX(I)/DEN1(I) * DSNO3(I)/Q0(I) )
        DSFEIII(I) = DSFEIII(I) + &
                 ( FLUX(I-1)/DXX(I)/DEN1(I) * DSFEIII(I-1)/Q0(I-1)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEIII(I)/Q0(I))
        DSFEII(I)  = DSFEII(I) + &
                 ( FLUX(I-1)/DXX(I)/DEN1(I) * DSFEII(I-1)/Q0(I-1)-&
                   FLUX(I)/DXX(I)/DEN1(I) * DSFEII(I)/Q0(I))
        ELSE IF( FLUX(I).LT. 0. ) THEN

        DSSO4(I) = DSSO4(I) + &
                  ( FLUX(I-1)/DXX(I)/DEN1(I) * DSSO4(I-1)/Q0(I-1) -&
                   FLUX(I)/DXX(I)/DEN1(I) * DSSO4(I+1)/Q0(I+1))
        DSNO3(I) = DSNO3(I) + &
                  ( FLUX(I-1)/DXX(I)/DEN1(I) * DSNO3(I-1)/Q0(I-1)-&
                   FLUX(I)/DXX(I)/DEN1(I) * DSNO3(I+1)/Q0(I+1))
        DSFEIII(I) = DSFEIII(I) + &
                  ( FLUX(I-1)/DXX(I)/DEN1(I) * DSFEIII(I-1)/Q0(I-1)-&
                   FLUX(I)/DXX(I)/DEN1(I) * DSFEIII(I+1)/Q0(I+1))
        DSFEII(I)  = DSFEII(I) + &
                  ( FLUX(I-1)/DXX(I)/DEN1(I) * DSFEII(I-1)/Q0(I-1)-&
                   FLUX(I)/DXX(I)/DEN1(I) * DSFEII(I+1)/Q0(I+1))
        ENDIF ! FLUX  


        DSSO4(I)   = AMAX1(DSSO4(I),   1.E-20)
        DSNO3(I)   = AMAX1(DSNO3(I),   1.E-20)
        DSFEIII(I) = AMAX1(DSFEIII(I), 1.E-20)
        DSFEII(I)  = AMAX1(DSFEII(I),  1.E-20)
 
         ENDIF        
      10 CONTINUE

!!!!    UPDATE MIXING RATIOS AND LIMIT FLUXES GOING DOWN WHERE U2D < 0
       
       IF(U(NZZ-1).LT.ZR0) FLUX(NZZ-1)= Q0(NZZ)*U(NZZ-1)*DT*DD0(NZZ-1)
       DO 20 I=NZZ-1,2,-1
        IF(U(I-1).GE.ZR0) THEN
          IF(U(I).LT.ZR0) THEN
            QN(I) =  &                   !INFLOW-ONLY CELL
           (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUX(I-1)/DXX(I))/DEN1(I)

           IF(FLUX(I).GT. 0. ) THEN
            DSSO4(I) = DSSO4(I) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DSSO4(I-1)/Q0(I-1) - &
                     FLUX(I)/DXX(I)/DEN1(I) * DSSO4(I)/Q0(I))
            DSNO3(I) = DSNO3(I) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DSNO3(I-1)/Q0(I-1)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSNO3(I)/Q0(I))
            DSFEIII(I) = DSFEIII(I) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DSFEIII(I-1)/Q0(I-1)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEIII(I)/Q0(I))
            DSFEII(I)  = DSFEII(I) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DSFEII(I-1)/Q0(I-1)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEII(I)/Q0(I))
           ELSE IF( FLUX(I).LT. 0. ) THEN

            DSSO4(I) = DSSO4(I) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DSSO4(I-1)/Q0(I-1) -&
                     FLUX(I)/DXX(I)/DEN1(I) * DSSO4(I+1)/Q0(I+1))
            DSNO3(I) = DSNO3(I) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DSNO3(I-1)/Q0(I-1)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSNO3(I+1)/Q0(I+1))
            DSFEIII(I) = DSFEIII(I) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DSFEIII(I-1)/Q0(I-1)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEIII(I+1)/Q0(I+1))
            DSFEII(I)  = DSFEII(I) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DSFEII(I-1)/Q0(I-1)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEII(I+1)/Q0(I+1))
        ENDIF ! FLUX   

            DSSO4(I)   = AMAX1(DSSO4(I),   1.E-20)
            DSNO3(I)   = AMAX1(DSNO3(I),   1.E-20)
            DSFEIII(I) = AMAX1(DSFEIII(I), 1.E-20)
            DSFEII(I)  = AMAX1(DSFEII(I),  1.E-20)

          ENDIF ! U(I).LT.ZRO

        ELSE ! U(I-1)
          X1  = DT*ABS(U(I-1))/DXX(I)
          CF1 = CFCALC1_DS(X1,Q0(I),Q0(I+1),Q0(I-1),IMXMN(I+1),IMXMN(I-1),LSHP)
          IF(U(I).GE.ZR0) CF1 = Q0(I)                  !OUTFLOW-ONLY CELL
          FLUXIM1 = DD0(I-1)*U(I-1)*DT*CF1/DXX(I)
          QN(I)= MAX( VCMIN(I), MIN( VCMAX(I), &
                 (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUXIM1)/DEN1(I) ))
          
          FLUX(I-1) = DXX(I)*(QN(I)*DEN1(I)-Q0(I)*DEN0(I)) + FLUX(I)

      IF( FLUX(I).LT. 0. ) THEN
            DSSO4(I) = DSSO4(I) -  &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DSSO4(I)/Q0(I) -&
                     FLUX(I)/DXX(I)/DEN1(I) * DSSO4(I+1)/Q0(I+1))
            DSNO3(I) = DSNO3(I) - &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DSNO3(I)/Q0(I)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSNO3(I+1)/Q0(I+1))
            DSFEIII(I) = DSFEIII(I) + &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DSFEIII(I)/Q0(I)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEIII(I+1)/Q0(I+1))
            DSFEII(I)  = DSFEII(I) + &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DSFEII(I)/Q0(I)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEII(I+1)/Q0(I+1))
       ELSE IF (  FLUX(I).GT. 0. ) THEN
            DSSO4(I) = DSSO4(I) -  &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DSSO4(I)/Q0(I) -&
                     FLUX(I)/DXX(I)/DEN1(I) * DSSO4(I)/Q0(I))
            DSNO3(I) = DSNO3(I) - &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DSNO3(I)/Q0(I)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSNO3(I)/Q0(I))
            DSFEIII(I) = DSFEIII(I) + &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DSFEIII(I)/Q0(I)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEIII(I)/Q0(I))
            DSFEII(I)  = DSFEII(I) + &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DSFEII(I)/Q0(I)-&
                     FLUX(I)/DXX(I)/DEN1(I) * DSFEII(I)/Q0(I))
       ENDIF      ! FLUX
            DSSO4(I)   = AMAX1(DSSO4(I),   1.E-20)
            DSNO3(I)   = AMAX1(DSNO3(I),   1.E-20)
            DSFEIII(I) = AMAX1(DSFEIII(I), 1.E-20)
            DSFEII(I)  = AMAX1(DSFEII(I),  1.E-20)

        ENDIF ! U(I-1)
      20 CONTINUE
      
        RETURN
       END    

      FUNCTION CFCALC1_DS(X1,VCUP1,VCUP2,VCDW1,IMXUP2,IMXDN1,LSHP)
! function to calculate mixing ratio in the fuild moving into
! neighnor
! cell durint DT; X1=Courant No. set LSHP=1 for full
! sharpening
        REAL :: CFCALC1_DS
        ALFA=1.
        IF(IMXDN1.gt.0)ALFA=1.7763-0.5*X1
        IF(IMXUP2.gt.0)ALFA=1.2578+0.58502*X1
        IF(LSHP==1)ALFA=5.
        IF(X1.lt.0.5)THEN
          ALFA=MIN(ALFA,.98*6./(1.-2.*X1))
          CF=VCUP1*(1.+ALFA/6.*(2.*X1-1.)) + &
                   VCDW1*ALFA/12.*(4-5.*X1) + VCUP2*ALFA/12.*(X1-2.)
        ELSE
          X1=MAX(X1,1.E-20)
          CF=VCUP1*(1.-ALFA*(1./X1 + 2.*X1 -3.)/6.) + &
             VCDW1*ALFA*(1./X1-X1)/12. + VCUP2*ALFA*(1./X1+5.*X1-6.)/12.
        ENDIF
! limit outflow mixing ratio to reasonable mixing ratio
        CFCALC1_DS=MIN(MAX(CF,MIN(VCUP1,VCDW1)),MAX(VCUP1,VCDW1))
        RETURN
        END        
                
