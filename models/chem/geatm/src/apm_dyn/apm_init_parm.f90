
! get apm parameters necessary to conduct following dynamical process :
! (1) apm vars in naqpms increasing due to emission
! (2) apm vars in naqpms variation due to horizontal and vertical advection
! (3) apm vars in naqpms variation due to horizontal and vertical diffusion
! (4) apm vars in naqpms variation due to convective cloud transport
! (5) apm vars in naqpms variation due to gravity settling
! (6) apm vars in naqpms variation due to dry deposition

subroutine apm_init_parm &
 ( CEMITSULF2, DFSALT9, CEMITBCOC2 &
 , RDRY, RSALT, RDST &
 , TOTNUMBC, TOTNUMOC )

! TOTNUMBC : *MBCOC8 is bcoc number concentration 




USE APM_INIT_MOD, ONLY : IFNUCL,IFAG,APM_INIT,APM_NTRACERS
USE APM_INIT_MOD, ONLY : XMACID,XMLVSOG,M1ACID,M1LVSOG
USE APM_INIT_MOD, ONLY : VDRY
USE APM_INIT_MOD, ONLY : DENSULF
implicit none
include 'apm_parm.inc'

! CEMITSULF2 : the particle mass conc (kg/m3) in each bin of each mode per
!              unit mass (kg/m3) of primary sulfate emitted. 
!              Multiple CEMITSULF2(NSO4,2) with total mass of primary sulfate 
!              emission to obtain sulfate mass emitted into each bin.
!
!  CEMITBCOC2(NMAX,2) is the # of BCOC particles distributed to each bin of  per
!  unit mass (kg) of primary BCOC. Multiple CEMITBCOC2(NNMAX,2)
!  with total BC or OC mass conc in the grid to obtain # conc.

call APM_NTRACERS( 0, N_APMTRAC )
call APM_INIT


end subroutine apm_init_parm
