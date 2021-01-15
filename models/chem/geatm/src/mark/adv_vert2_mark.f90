          SUBROUTINE ADV_VERT_MARK(MYID,C,W,U2D,V2D,DELTZ,&
                DELTX,SX,EX,SY,EY,NZZ,DT,IG,ISMMAX,SMCONV,KKTOP)

          INTEGER :: MYID,SX,EX,SY,EY,ISMMAX,NZZ
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1,NZZ) :: C,W,U2D,V2D
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1,NZZ) :: DELTX,DELTZ
          REAL    :: DT,D0 
          REAL,DIMENSION(NZZ)                     :: Q0,QN,DXX0,DXX, &
                                                     DEN0,DEN1,U,DD0
          REAL,DIMENSION(ISMMAX,SX-1:EX+1,SY-1:EY+1,NZZ) ::  SMCONV
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1)     ::  KKTOP
          REAL,DIMENSION(ISMMAX,NZZ)    ::  TMPsmconv
          DATA D0 /1.0/                                           
          
          DO I=SX,EX
           DO J=SY,EY
            DO  K=2,NZZ-1
             DXX0(K)  =  DELTX(I,J,K)
             DXX (K)  =  DELTZ(I,J,K)
             DEN0(K)  =  1.- DT/DXX0(K)*(D0*U2D(I,J,K)-D0*U2D(I-1,J,K))&
                           - DT/DXX0(K)*(D0*V2D(I,J,K)-D0*V2D(I,J-1,K))
             DEN1(K)  =  DEN0(K)-DT/DXX(K)*(D0*W(I,J,K)-D0*W(I,J,K-1)) 
             
             DD0 (K)  =  DEN0(K)
             Q0  (K)  =  C(I,J,K)
             U   (K)  =  W(I,J,K)
            ENDDO  !K

             DD0 (1)  =  D0
             U   (1)  =  W(I,J,1)
             Q0  (1)  =  C(I,J,1)
             Q0  (NZZ)=  C(I,J,NZZ)
             KKKTOP   =  KKTOP(I,J)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!source appointent
              DO  is=1,ismMax
              DO K=1,NZZ  
               TMPsmconv(is,k)=smconv(is,i,j,k)
              ENDDO !K
              ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
            CALL ADVEC1D_MARK(Q0,QN,U,DEN0,DEN1,DT,DXX,DD0,NZZ,ISMMAX,&
                        TMPsmconv,KKKTOP,I,J,IG)

            DO K=2,NZZ-1

             C(I,J,K) = QN(K)
             DO  is=1,ismMax
              smconv(is,i,j,k)=TMPsmconv(is,k)
             ENDDO !IS 
             
            ENDDO  !K

            
           ENDDO   !J
          ENDDO !I

          RETURN
          END


         SUBROUTINE ADVEC1D_MARK(Q0,QN,U,DEN0,DEN1,DT,DXX,DD0,NZZ,&
                          ISMMAX,TMPsmconv,KKKTOP,II,JJ,IG)
         INTEGER              :: NZZ
         REAL,DIMENSION(NZZ)  :: Q0,QN,U,DEN0,DEN1,DXX,DD0
         REAL                 :: DT,CK1,X1,FLUXI,FLUXIM1,CF1,QQ0
         REAL,DIMENSION(NZZ)  :: FLUX,VCMAX,VCMIN,IMXMN
         INTEGER              :: ISMMAX,KKKTOP 
         REAL,DIMENSION(ISMMAX,NZZ) :: TMPsmconv
         REAL,DIMENSION(ISMMAX):: smthis,smother
         DATA ZR0,LSHP / 0.0,0.0/

         
         IMXMN = 0.0 
         
          
!!!!  INDENTIFY LOCAL MAX AND MIN,SPECIFY MIXING RATIO LIMITS AT NEW TIME
 
        DO 5 I=2,NZZ-1

         IMXMN(I) = 0
         IF(Q0(I).GE.MAX(Q0(I-1),Q0(I+1)).OR. &
            Q0(I).LE.MIN(Q0(I-1),Q0(I+1)) ) IMXMN(I)=1
         CK1 = Q0(I)
         IF(U(I).LT.ZR0)   CK1 = Q0(I+1)
         IF(U(I-1).GE.ZR0) CK1 = Q0(I-1)
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
            CF1  = CFCALC_MARK(X1,Q0(I),Q0(I-1),Q0(I+1),IMXMN(I-1),IMXMN(I+1),LSHP)
            FLUXI= DD0(I)*U(I)*DT*CF1/DXX(I)
            QN(I)= MAX(VCMIN(I) ,MIN(VCMAX(I), &
                  (Q0(I)*DEN0(I)-FLUXI+FLUX(I-1)/DXX(I))/DEN1(I) ) )
!!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            IF(II==53.AND.JJ==30.AND.IG==11) &
!             PRINT*,I,Q0(I)*DEN0(I)/2./DEN1(I)+Q0(I)/2.-Q0(I),FLUX(I-1)/DXX(I)/DEN1(I),(-FLUXI)/DEN1(I)
                  
            IF(FLUX(I-1)>=0.0) THEN
              DO is=1,ismMax
               smthis(is) = TMPsmconv(is,I)
               smother(is)= TMPsmconv(is,I-1) 
              ENDDO      

!              deltc1 = (Q0(I)*DEN0(I)/2.+FLUX(I-1)/DXX(I))/DEN1(I)-Q0(I)/2.
              deltc1 =  FLUX(I-1)/DXX(I)/DEN1(I)
              deltc1 = MAX(deltc1,1.E-20)
              Q0(I)  = MAX(Q0(I),1.E-20)

              IF(I>KKKTOP) THEN
               do is =1 , ismMax
               if(is==2)then
                 smthis(is)=1.0
               else
                 smthis(is) =0.0
               endif
               enddo 
              ELSE
               call SMmixing(Q0(I),smthis,deltc1,smother,ismMax)        
              ENDIF !I

              DO is=1,ismMax
               TMPsmconv(is,I)=smthis(is)
              ENDDO
       
            ENDIF     ! FLUX 
            
             QQ0 = (Q0(I)*DEN0(I)+FLUX(I-1)/DXX(I))/DEN1(I)

            IF(FLUXI<= 0.) THEN
             DO is=1,ismMax       
              smthis(is) = TMPsmconv(is,I)
              smother(is)= TMPsmconv(is,I+1)
             ENDDO 
            
            deltc1 = (-FLUXI)/DEN1(I)
            deltc1 = max(deltc1,1.E-20)
            QQ0    = max(QQ0, 1.E-20)            

            IF(I==KKKTOP) THEN
              call GetBounChange(deltc1,QQ0,smthis,2,ismMax)       
            ELSE IF(I>KKKTOP)THEN
             do is = 1 ,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) = 0.0
               endif
             enddo
            ELSE
              call SMmixing(QQ0,smthis,deltc1,smother,ismMax)
            ENDIF !I             
            
            DO is=1,ismMax
              TMPsmconv(is,I)=smthis(is)
            ENDDO
 
           ENDIF !FLUX(I)            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  
            FLUX(I) = DXX(I)*(Q0(I)*DEN0(I)-QN(I)*DEN1(I))+FLUX(I-1)   
         ENDIF        
      10 CONTINUE

!!!!    UPDATE MIXING RATIOS AND LIMIT FLUXES GOING DOWN WHERE U2D < 0
       
       IF(U(NZZ-1).LT.ZR0) FLUX(NZZ-1)= Q0(NZZ)*U(NZZ-1)*DT*DD0(NZZ-1)
       DO 20 I=NZZ-1,2,-1
        IF(U(I-1).GE.ZR0) THEN
          IF(U(I).LT.ZR0) QN(I) =  &                   !INFLOW-ONLY CELL
           (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUX(I-1)/DXX(I))/DEN1(I)

!!!!!!!!!!!!!!!!!!!!!!!!!!  SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!
!        IF(II==53.AND.JJ==30.AND.IG==11) &
!          PRINT*,I,Q0(I)*DEN0(I)/2./DEN1(I)+Q0(I)/2.-Q0(I),-FLUX(I)/DXX(I)/DEN1(I),FLUX(I-1)/DXX(I)/DEN1(I)

           IF(FLUX(I)<=0) THEN
              DO is=1,ismMax
               smthis(is) = TMPsmconv(is,I)
               smother(is)= TMPsmconv(is,I+1)
              ENDDO
              
!              deltc1=(Q0(I)*DEN0(I)/2.-FLUX(I)/DXX(I))/DEN1(I)-Q0(I)/2.
              deltc1 = -FLUX(I)/DXX(I)/DEN1(I)
              deltc1 = MAX(deltc1,1.E-20)
              Q0(I)  = max(Q0(I),1.E-20)              

              IF(I==KKKTOP) THEN
                call GetBounChange(deltc1,Q0(I),smthis,2,ismMax)    
              ELSE IF(I>KKKTOP)THEN
               do is = 1 , ismMax
                if(is==2) then
                  smthis(is)=1.0
                else
                  smthis(is) = 0.0
                endif
               enddo
              ELSE
                call SMmixing(Q0(I),smthis,deltc1,smother,ismMax)      
              ENDIF !I

              DO is=1,ismMax
                TMPsmconv(is,I)=smthis(is)
              ENDDO   
             
            ENDIF  !FLUX(I)
             
             QQ0 = Q0(I) - FLUX(I)/DXX(I)/DEN1(I)
             
            IF(FLUX(I-1)>=0) THEN
              DO is=1,ismMax            
               smthis(is) = TMPsmconv(is,I)
               smother(is)= TMPsmconv(is,I-1)
              ENDDO 

!             deltc1=(Q0(I)*DEN0(I)/2.+FLUX(I-1)/DXX(I))/DEN1(I)-Q0(I)/2.
              deltc1 = FLUX(I-1)/DXX(I)/DEN1(I)
              deltc1 = MAX(deltc1,1.E-20)
              QQ0    = MAX(QQ0,1.E-20)
             
             IF(I>KKKTOP) THEN
               do is = 1 ,ismMax
                 if(is==2) then
                   smthis(is)=1.0      
                 else 
                   smthis(is) = 0.0
                 endif
               enddo
             ELSE
               call SMmixing(QQ0,smthis,deltc1,smother,ismMax)
             ENDIF  !I

             DO is=1,ismMax
               TMPsmconv(is,I)=smthis(is)
             ENDDO  
            
            ENDIF  ! FLUX(I-1) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
        ELSE
          X1  = DT*ABS(U(I-1))/DXX(I)
          CF1 = CFCALC_MARK(X1,Q0(I),Q0(I+1),Q0(I-1),IMXMN(I+1),IMXMN(I-1),LSHP)
          IF(U(I).GE.ZR0) CF1 = Q0(I)                  !OUTFLOW-ONLY CELL
          FLUXIM1 = DD0(I-1)*U(I-1)*DT*CF1/DXX(I)
          QN(I)= MAX( VCMIN(I), MIN( VCMAX(I), &
                 (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUXIM1)/DEN1(I) ))
!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         IF(II==53.AND.JJ==30.AND.IG==11) &
!           PRINT*,I,Q0(I)*DEN0(I)/2./DEN1(I)+Q0(I)/2.-Q0(I),FLUXIM1/DEN1(I),-FLUX(I)/DXX(I)/DEN1(I)
                 
          IF(FLUXIM1>=0) THEN               
            DO is=1,ismMax      
             smthis(is) = TMPsmconv(is,I)
             smother(is)= TMPsmconv(is,I-1)
            ENDDO     

!           deltc1 = (Q0(I)*DEN0(I)/2.+FLUXIM1)/DEN1(I) - Q0(I)/2.
            deltc1 = FLUXIM1/DEN1(I)
            deltc1 = MAX(deltc1,1.E-20)
            Q0(I) = MAX(Q0(I),1.E-20)              

            IF(I>KKKTOP) THEN
              do is = 1,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) = 0.0
               endif
              enddo
            ELSE
              call SMmixing(Q0(I),smthis,deltc1,smother,ismMax)      
            ENDIF  !I  

            DO is=1,ismMax
             TMPsmconv(is,I)=smthis(is)
            ENDDO 
            
          ENDIF  ! FLUX(I-1)

           QQ0 = (Q0(I)*DEN0(I)+FLUXIM1)/DEN1(I) 

          IF(FLUX(I)<=0) THEN
            DO is=1,ismMax
             smthis(is) = TMPsmconv(is,I)
             smother(is)= TMPsmconv(is,I+1)
            ENDDO
           
           
!           deltc1 = (Q0(I)*DEN0(I)/2.-FLUX(I)/DXX(I))/DEN1(I)-Q0(I)/2.  
           deltc1 = -FLUX(I)/DXX(I)/DEN1(I)
           deltc1 = MAX(deltc1,1.E-20)
           QQ0 = MAX(QQ0,1.E-20)
             
            IF(I==KKKTOP) THEN             
             call GetBounChange(deltc1,QQ0,smthis,2,ismMax)
            ELSE IF(I>KKKTOP) THEN
             do is = 1, ismMax
              if(is==2) then
               smthis(is)=1.0
              else
               smthis(is) = 0.0
              endif
             enddo 
            ELSE 
             call SMmixing(QQ0,smthis,deltc1,smother,ismMax)
            ENDIF      !I
            
           DO is=1,ismMax          
             TMPsmconv(is,I)=smthis(is)
           ENDDO
                              
          ENDIF ! FLUX(I)        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                 
          FLUX(I-1) = DXX(I)*(QN(I)*DEN1(I)-Q0(I)*DEN0(I)) + FLUX(I)
        ENDIF
      20 CONTINUE
      
        RETURN
       END    

      FUNCTION CFCALC_MARK(X1,VCUP1,VCUP2,VCDW1,IMXUP2,IMXDN1,LSHP)
! function to calculate mixing ratio in the fuild moving into
! neighnor
! cell durint DT; X1=Courant No. set LSHP=1 for full
! sharpening
        REAL :: CFCALC_MARK   
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
        CFCALC_MARK=MIN(MAX(CF,MIN(VCUP1,VCDW1)),MAX(VCUP1,VCDW1))
        RETURN
        END        
                
