! For the simulation of O2 according to OCMIP-2
! Reference:Najjar R. G., Jin X., Louanchi F. et al. Impact of circulation on export production, dissolved organic matter, 
!   and dissolved oxygen in the ocean: results from Phase II of the Ocean Carbon-cycle Model Intercomparison Project (OCMIP-2).
!   Global Biogeochemical cycles, 2007, 21, GB3007, doi: 10.1029/2006GB002857
!
! by cm
! 090330
!===================================
subroutine sgeo2

!-----------------------------------
! purpose: calculate transport velocity of o2
! kw02=b*(Sc/660)**(-1/2)*u**2
! b=0.336
!-----------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod
      USE carbon_mod
      USE tracer_mod
      USE cforce_mod
#ifdef SPMD      
!      USE msg_mod
#endif      
!
implicit none
!
real sco2(imt,jmt),ttmp!(imt,jmt)
!-------------------------------------
! ttmp:sea surface temperture (in deg C)
! sco2:Schmidt number of o2 in sea water, by Keeling et al. (1998, Global Biogeochem.
!      Cycles, 12, 141-163).
!-------------------------------------
!$OMP PARALLEL DO PRIVATE (I,J)        
DO J=1,JMT
   DO I=1,IMT
        sco2(i,j)=0.
        IF( vit(I,J,1)<0.5 ) cycle
        ttmp     =AT(I,J,1,1)
	sco2(i,j)=1638.0 - 81.83*ttmp + 1.483*ttmp**2 - 0.008004*ttmp**3
   enddo
enddo
!
!$OMP PARALLEL DO PRIVATE (I,J)        
DO J=1,JMT
   DO I=1,IMT
      sge(i,j)=0.
      IF( vit(I,J,1)<0.5 ) cycle
      sge(i,j)=0.336*(660./sco2(i,j))**0.5*w22np(i,j)**2/3.6e5
   enddo
enddo

return

end subroutine sgeo2

!===================================

subroutine flux_o2
!-----------------------------------
! purpose: calculate flux of o2
!
! o2sato: the oxygen saturation concentration at 1 atm total pressure
! in umol/kg given the temperature (t, in deg C) and the salinity (s,
! in permil). 
!
! FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
!
! Note: o2sato IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
! 0 permil <= S <= 42 permil
!-----------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod
      USE carbon_mod
      USE tracer_mod
      USE cforce_mod
      USE forc_mod
#ifdef SPMD      
!      USE msg_mod
#endif      
!
implicit none
!
real o2sato(imt,jmt),o2surf(imt,jmt)
real aO,a1,a2,a3,a4,a5,b0,b1,b2,b3,CO,cc
real ttmp(imt,jmt),stmp(imt,jmt),tt,tk,ts,ts2,ts3,ts4,ts5

!$OMP PARALLEL DO PRIVATE (I,J)        
DO J=1,JMT
   DO I=1,IMT
        ttmp(i,j)=0.
	stmp(i,j)=0.
	o2surf(i,j)=0.
        IF( vit(I,J,1)<0.5 ) cycle
!
        ttmp(i,j)=AT(I,J,1,1)
	if(ttmp(i,j)>40.) ttmp(i,j)=40.
!
        stmp(i,j)=AT(I,J,1,2)*1000.+35.
	if(stmp(i,j)<0.) stmp(i,j)=0.
	if(stmp(i,j)>42.) stmp(i,j)=42.
!
        o2surf(i,j)=ptb(I,J,1,5)
   enddo
enddo

AO=5.80818
A1=3.20684
A2=4.11890
A3=4.93845
A4=1.01567
A5=1.41575
B0=-0.00701211
B1=-0.00725958
B2=-0.00793334
B3=-0.00554491
CO=-1.32412E-7
!
!$OMP PARALLEL DO PRIVATE (I,J)        
DO J=1,JMT
   DO I=1,IMT
      o2sato(i,j) =0.
      IF( vit(I,J,1)<0.5 ) cycle
      TT  = 298.15-ttmp(i,j)
      TK  = 273.15+ttmp(i,j)
      TS  = LOG(TT/TK)
      TS2 = TS**2
      TS3 = TS**3
      TS4 = TS**4
      TS5 = TS**5
      CC  = AO + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5 &
            + stmp(i,j)*(B0 + B1*TS + B2*TS2 + B3*TS3)     &
	    + CO*(stmp(i,j)**2)
      o2sato(i,j) = exp(CC)
   enddo
enddo

call sgeo2

!$OMP PARALLEL DO PRIVATE (I,J)        
DO J=1,JMT
   DO I=1,IMT
      ssfc(i,j)=0.
      IF( vit(I,J,1)<0.5 ) cycle
      ssfc(i,j)=(1.0-iceday(i,j))*sge(I,J)*(o2sato(i,j)*pressureday(i,j)-o2surf(i,j))
   enddo
enddo

return

end subroutine flux_o2
