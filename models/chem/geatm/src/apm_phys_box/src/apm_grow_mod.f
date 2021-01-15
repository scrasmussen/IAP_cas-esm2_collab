! $Id: apm_grow_mod.f,v 0.0 2008/08/23 11:30:00 fyu $
      MODULE APM_GROW_MOD
!
!******************************************************************************
!  Module APM_GROW_MOD contains variables and routines for computing size-
!  resolved particle growth. (fyu,8/23/08)
!
!  Module Variables:
!  ============================================================================
!   +++++
!
!  Module Routines:
!  ============================================================================
!
!  modules referenced by "apm_grow_mod.f"
!  ============================================================================
!  to be modified ++++
!
!  NOTES:
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "apm_grow_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...

      ! ... and these routines
      PUBLIC :: APM_GROW, APM_MOVEBIN

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE APM_GROW(ICOND,NSO4,TK,PRESS,CCOND,PCOND,
     &   TCSOTHER,DTNG,XN,XVA,RWETCM,TOTCONDOTHER,XMCOND)
!
!******************************************************************************
!  Subroutine APM_GROW calculates H2SO4 condensational growth of aerosol particles
!  (fyu, 9/16/08)
! Input:
! TK -- Temperature (K)
! PRESS -- Pressure (pa)
! CCOND -- Condensable vapor concentration (#/cm3)
! PCOND -- Condensable vapor production rate (#/cm3s)
! DTNG -- time detp for growth (s)
! XN -- number conc of each bin (#/cm3)
! XVA -- total acid volume of each bin (cm3/cm3)
!  
! Output:
!   updated XMA, CCOND

!  NOTES:
! XVA -- total acid volmue of each bin (cm3/cm3)
!  (1 ) 
!******************************************************************************
!
      USE APM_INIT_MOD, ONLY: DENSULF,RDRY
      USE APM_INIT_MOD, ONLY : ONEPI,BK,AVG,RGAS

      ! Local variables
      INTEGER :: NSO4,ICOND
      INTEGER :: N, NMAX

      REAL*8  :: XN(NSO4),XVA(NSO4),XMA(NSO4)
      REAL*8  :: YF(NSO4),YGR(NSO4),RWETCM(NSO4),AKELV(NSO4)

      REAL*8  :: TK,PRESS,CCOND,PCOND,DTNG,CSCOND,TCSOTHER,YFSUM1
      REAL*8  :: XMCOND,V1COND,VCOND,TEMP0,WTGAS
      REAL*8  :: WTAIR,DIAMAIR,RHOA,CDIFUS,CDIFUS2,DIFUSC,FREEPD
      REAL*8  :: YY,YKN,FCORR,TEMP1,YFSUM,YEVAP
      REAL*8  :: YEVAPV,CCOND1,CCONDA,TOTEVAP,YGRSUM,YDV,TOTCOND
      REAL*8  :: XRCM,TOTCONDOTHER
      REAL*8  :: DENCOND,AKELV0,SURFT,YTEMP1
      REAL*8  :: DVA(NSO4)
      !=================================================================
      ! APM_GROW begins here!
      !=================================================================
!
      DENCOND = DENSULF
      IF(ICOND.EQ.1) THEN !H2SO4
       CSCOND = 1.E4  ! no acid evap now
       AKELV0 = CCOND/CSCOND
       IF(AKELV0.LT.1.2) RETURN 
       AKELV =1.

      ELSE
       WRITE(6,*)"NEED to check ICOND value"
       STOP
      ENDIF

      V1COND = XMCOND/(AVG*DENCOND)  ! volume of one molecule (cm3)

      VCOND = SQRT(8.*RGAS*TK/(ONEPI*XMCOND))   ! cm/s
      TEMP0 = VCOND * XMCOND/AVG *ONEPI

! Cal. gas diffusion coef. and mean free path
!
      WTGAS = XMCOND
      WTAIR = 28.966
      DIAMAIR   = 4.5E-08
      RHOA = WTAIR*PRESS*10./(RGAS*TK)   ! AIR DENSITY (G CM-3)
      CDIFUS =3.*SQRT(0.5*RGAS*WTAIR/ONEPI)/(AVG*8.*DIAMAIR**2)
      CDIFUS2=CDIFUS * SQRT((WTGAS+WTAIR)/WTGAS)
      DIFUSC=CDIFUS2 * SQRT(TK)/RHOA
      FREEPD = 3.*DIFUSC/VCOND
!
       YF = 0.
       YFSUM = 0.
       YY = PCOND
       DO N = 1, NSO4
          XRCM = RWETCM(N)   ! RWETCM in cm
          YKN = FREEPD/XRCM
          FCORR = YKN/(0.75+YKN)
          TEMP1 = TEMP0*FCORR
          IF(AKELV(N).LT.AKELV0) THEN
           TEMP1 = TEMP1 *(1.-AKELV(N)/AKELV0)
          ELSE
           TEMP1 = 0.  ! no evapor for now
          ENDIF
          YF(N) = TEMP1 * XN(N) *XRCM*XRCM/(DENCOND*V1COND)   ! s-1
          YFSUM = YFSUM + YF(N)
!          YEVAP = CSCOND*AKELV*YF(N)    ! in #/cm3s
!          YEVAPV = YEVAP*DTNG*V1COND      ! in cm3/cm3
!          IF(YEVAPV.GT. XVA(N)) THEN
!             YEVAP = XVA(N)/(DTNG*V1COND)
!          ENDIF
!          YY = YY + YEVAP
       ENDDO
       IF(YFSUM.EQ.0.) THEN
         WRITE(6,100)AKELV0,(AKELV(N),N=1,NSO4,5)
         WRITE(6,100)FCORR,(YF(N),N=1,NSO4,5)
         STOP
       ENDIF
! Consider scavenging of H2SO4 vapor by aerosols other than sulfate
       YFSUM1 = YFSUM + TCSOTHER
       CCOND1 = YY/YFSUM1 + (CCOND - YY/YFSUM1)*exp(-YFSUM1*DTNG)

       CCONDA = 0.5 * (CCOND1 + CCOND)

       TOTEVAP = 0.
!       YGRSUM = 0.
!       DO N = 1, NSO4
!          YGR(N) = YF(N)*(CCONDA - CSCOND*AKELV)
!          IF(YGR(N).LT.0.) THEN
!            YDV = -YGR(N)*DTNG*V1COND    !cm3
!            IF(YDV.GT.XVA(N)) THEN
!               TOTEVAP = TOTEVAP + XVA(N)
!               XVA(N) = 1.E-50
!            ELSE
!               XVA(N) = XVA(N) - YDV
!               TOTEVAP = TOTEVAP + YDV
!            ENDIF
!          ELSE
!            YGRSUM = YGRSUM + YGR(N)
!          ENDIF
!       ENDDO
!
       TOTCOND=PCOND*DTNG+(CCOND-CCOND1)+TOTEVAP/V1COND   ! #/cm3
!  H2SO4 vapor condensing on particles other than sulfate
       TOTCONDOTHER = TOTCOND * (1.-YFSUM/YFSUM1)
!  H2SO4 vapor condensing on sulfate
       TOTCOND = TOTCOND * YFSUM/YFSUM1

       YGR = YF
       YGRSUM = YFSUM

       DO N = 1, NSO4
          DVA(N) = TOTCOND*V1COND*YGR(N)/YGRSUM
          XVA(N)=XVA(N)+DVA(N)
       ENDDO


 100  FORMAT(20(1PE9.2))
!
! Update gas concentration
       CCOND = CCOND1
!

      END SUBROUTINE APM_GROW

!------------------------------------------------------------------------------
! 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE APM_MOVEBIN(NSO4,XN,XVA)
!
!******************************************************************************
!  Subroutine APM_MOVEBIN moves particles across bins after growth

! XN -- number conc of each bin (#/cm3)
! XVA -- total acid mass of each bin (cm3/cm3)
!
! NC: number of core component associated with particles, set to 1 now,
! may increase later


      USE APM_INIT_MOD, ONLY: VDRY

      ! Local variables
      INTEGER :: NSO4
      INTEGER, PARAMETER   :: NC = 1
      INTEGER :: N, NMAX, IC, J, J0, J1, JS

      REAL*8  :: XN(NSO4),XVA(NSO4)
      REAL*8  :: YN(NSO4),YVA(NSO4,NC),YC(NSO4,NC)
      REAL*8  :: XVT, VDA, YFR, YFV

      NMAX = NSO4
!
!  Move the partilces across bins due to the condensation/evap
!
        DO N=1,NMAX
           YVA(N,1) = XVA(N) ! to be modified if NC>1
        ENDDO

!> shun ?? : only mathematical consideration ?
        DO n=1,NMAX
           YN(n)=1.E-20 ! shun? just a low value?
           DO IC=1,NC
              ! sulfate volume in each bin in unit volume
              YC(n,IC)=1./float(NC)*YN(n)*VDRY(n)*1.E6   ! in cm3/cm3
           ENDDO

!            IF(YN(n).GT.1.E6.or.YN(n).LE.0.) THEN
!              WRITE(6,*)"4 XN=",n,YN(n),XN(n),YC(n,1)
!            ENDIF

        ENDDO
!!!!!!!!!!!!!!


        DO n=1,NMAX
           XVT = 0.
           DO IC=1,NC
              XVT = XVT + YVA(n,IC)   ! cm3/cm3
           ENDDO
           ! particle average volume in each bin
           VDA = XVT/XN(n)*1.E-6      ! in m3
! if particles become smaller than first bin or larger than last bin, 
! scale and put in the first or last bin in a fashion that conserve mass
           IF(VDA.LT.VDRY(1)) THEN
              YN(1) = YN(1) + XN(n)*VDA/VDRY(1)  
              DO IC=1,NC               ! move compositions
                 YC(1,IC) = YC(1,IC)+YVA(n,IC)
              ENDDO

!            IF(YN(n).GT.1.E6.or.YN(n).LE.0.) THEN
!              WRITE(6,*)"41 XN=",n,YN(n),XN(n),YC(n,1),YVA(n,1),XVT,VDA
!            ENDIF

           ELSEIF(VDA.GT.VDRY(NMAX)) THEN
              YN(NMAX) = YN(NMAX) + XN(n)*VDA/VDRY(NMAX)  
              DO IC=1,NC               ! move compositions
                 YC(NMAX,IC) = YC(NMAX,IC)+YVA(n,IC)
              ENDDO
!            IF(YN(n).GT.1.E6.or.YN(n).LE.0.) THEN
!              WRITE(6,*)"42 XN=",n,YN(n),XN(n),YC(n,1),
!     &            XVA(n),YVA(n,1),XVT,VDA
!            ENDIF
           ELSE
            IF(VDA.GE.VDRY(n)) THEN
              J0=n
              J1 = NMAX-1
              JS = 1
            ELSE
              J0=n-1
              J1=1
              JS=-1
            ENDIF
            DO J=J0,J1,JS
              IF(VDA.GE.VDRY(J).AND.VDA.LT.VDRY(J+1)) THEN
                 YFR = (VDRY(J+1)-VDA)/(VDRY(J+1)-VDRY(J))
                 YFV = YFR*VDRY(J)/VDA
                 YN(J) = YN(J) + YFR*XN(n)  !move #
                 YN(J+1)=YN(J+1)+(1.-YFR)*XN(n)
                 DO IC=1,NC               ! move compositions
                    YC(J,IC) = YC(J,IC)+YFV*YVA(n,IC)
                    YC(J+1,IC)=YC(J+1,IC)
     &                                      +(1.-YFV)*YVA(n,IC)
                 ENDDO
                 GOTO 160
              ENDIF
            ENDDO
           ENDIF
 160       CONTINUE
        ENDDO

!  Update the bin values

         DO n=1,NMAX

            IF(YN(n).GT.1.E6.or.YN(n).LE.0.) THEN
!              WRITE(6,*)"5 XN=",n,YN(n),XN(n),YC(n,1)
            ENDIF

            XN(n) = YN(n)
            DO IC=1,NC
               YVA(n,IC) =YC(n,IC)
            ENDDO
            XVA(N) = YVA(n,1) ! to be modified if NC>1
         ENDDO


      END SUBROUTINE APM_MOVEBIN
!------------------------------------------------------------------------------
! 
      END MODULE APM_GROW_MOD
