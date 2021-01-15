      SUBROUTINE ALLOSEA (MYID, AER, SEACOMP, IDUC,  SX,EX, SY, EY )
       INTEGER :: MYID, SX, EX, SY, EY
       REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: AER, SEACOMP
       
!  THIS SUBROUTINE IS TO ALLOCATE THE SEA SALT COMPOSITIONS FROM THE
!  PUBLICATION Athanasopoulou et al., ACP, 2008; The role of sea salt
!  emissions and heterogeneous chemistry in the air quality of polluted
!  coastal areas.
!  MASS WEIGHT CL : 55%; Na: 31%; ss-SO4 : 8%; Mg 4%; Ca:1% , K 1%
      DO J = SY-1,EY+1
       DO I = SX-1, EX+1
         IF( IDUC == 1 ) SEACOMP(I,J) = AER (I, J) * 0.55 ! CL
         IF( IDUC == 2 ) SEACOMP(I,J) = AER (I, J) * 0.31 ! NA
         IF( IDUC == 3 ) SEACOMP(I,J) = AER (I, J) * 0.08 ! SS-SO4
         IF( IDUC == 4 ) SEACOMP(I,J) = AER (I, J) * 0.04 ! Mg
         IF( IDUC == 5 ) SEACOMP(I,J) = AER (I, J) * 0.01 ! Ca
         IF( IDUC == 6 ) SEACOMP(I,J) = AER (I, J) * 0.01 ! K
          SEACOMP(I,J) = AMAX1(SEACOMP(I,J), 1.E-20)

       ENDDO
      ENDDO


      RETURN
  
      ENDSUBROUTINE


      SUBROUTINE ALLODUST (MYID, AER, DUSTCOMP, IDUC,  SX,EX, SY, EY )
       INTEGER :: MYID, SX, EX, SY, EY
       REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: AER, DUSTCOMP

!  THIS SUBROUTINE IS TO ALLOCATE THE DUST COMPOSITIONS FROM THE
!  PUBLICATION Feng Yan  et al.,Global Modeling of Nitrate and Ammonium:
!  Interaction of Aerosols and Tropospheric Chemistry.
!  MASS WEIGHT CACO3: 7% MGCO3:5.5%; K2CO3 3.3%; NA2CO3 2.6% SIO2 60%; !  AL2O3 14.1%; FE2O3 6.9
      DO J = SY-1,EY+1
       DO I = SX-1, EX+1
         IF( IDUC == 1 ) DUSTCOMP(I,J) = AER (I, J) * 0.07  ! CACO3
         IF( IDUC == 2 ) DUSTCOMP(I,J) = AER (I, J) * 0.055 ! MGCO3
         IF( IDUC == 3 ) DUSTCOMP(I,J) = AER (I, J) * 0.033 ! K2CO3
         IF( IDUC == 4 ) DUSTCOMP(I,J) = AER (I, J) * 0.026 ! NA2CO3
         IF( IDUC == 5 ) DUSTCOMP(I,J) = AER (I, J) * 0.60  ! SIO2
         IF( IDUC == 6 ) DUSTCOMP(I,J) = AER (I, J) * 0.141 ! AL2O3
        DUSTCOMP(I,J) = AMAX1(DUSTCOMP(I,J), 1.E-20)

       ENDDO
      ENDDO


      RETURN

      ENDSUBROUTINE

