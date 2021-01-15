      SUBROUTINE GETVIS ( MYID, RNW, CLW, EXT, VIS, I, J, K)
! ccc   TO CALCULATE THE VISIBILITY IN KM

      REAL RNW, CLW, EXT, VIS, VIS1,VIS2,LWC
      INTEGER I,J,K
!      RNW       : RAINWATER IN KG/HG
!      CLW       : CLOUDWATER IN KG/KG
!      LWC       : LIQUID WATER IN G/M3
!      VIS       : VISIBILITY IN KM
!      VIS1      : VISIBILITY IN KM DUE TO FOG
!      VIS2      : VISIBILITY IN KM DUE TO AEROSOL
  
 
      LWC =  ( RNW + CLW )* 29. * 42.9    

!  IF NO FOG , VIS1 = 20.KM

      IF ( LWC .GE. 1.E-20) THEN
       VIS1 = - LOG(0.02) / ( 144.7 * LWC**0.88 ) 
      ELSE
       VIS1 = 80.
      ENDIF

! IF NO AEROSOL, VIS2 = 80 KM
      IF (EXT.GE.4.E-2) THEN
!       VIS2 = 3.5 / EXT 
        VIS2 = 10. * LOG(EXT/10.)
      ELSE
       VIS2 = 80.
      ENDIF

      VIS = MIN (VIS1, VIS2)

      RETURN
      ENDSUBROUTINE
