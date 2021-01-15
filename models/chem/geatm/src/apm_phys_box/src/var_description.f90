! YGF : hygroscopic growth factor
! REAL*8, ALLOCATABLE :: YGF(:,:) 'apm_init_mod.f' 99*4

! subroutine BINSULF setup size-resolved sulfate sectional bin structure

!  MSO4  : total sulfate mass concentration 
! XMA(N) : total mass in each bin in unit volume
! XVA(N) : particle volume in bin N 

! subroutine APM_MOVEBIN(NSO4,XN,XVA) : Move particles across bins after cloud chem

!! seasalt:
! XMSALT(N) : total seasalt mass concentration in each bin
! MCORE(2)  : total seasalt mass concentration
! XNSALT(N) : number concentration in each bin
! ZTN(2)    : total seasalt number concentration
!!

!! dust:
! MCORE(3) : total seasalt mass concentration
! XNDST(N) : number concentration in each bin
! ZTN(3)   : total seasalt number concentration
!!

!! bcoc:
! MCORE(4) : 
! MCORE(5) :
!!


 
