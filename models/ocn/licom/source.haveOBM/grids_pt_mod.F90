!     ================
      MODULE GRIDS_PT_MOD
!     ================
!     TOPOGRAPHY & GRIDS
!-----------------------------------------------------------------------
!
! Purpose: Set up some constants used in carbon model.
!
! Author: LYC,2013,5
!
!-----------------------------------------------------------------------
!zb(km)
!zt(km)

#include <def-undef.h>
use param_mod
use pconst_mod
use pmix_mod
#ifdef SPMD
use msg_mod,only: mpi_comm_ocn
#endif
use carbon_mod
      IMPLICIT NONE
!#include <netcdf.inc>
    real,dimension(km)::zb,ztop,zm
    real(kind=8)::vsea2000
    real,parameter::kb0=24
    real(kind=8),dimension(km)::sdxdy,sdxdy0
    integer,dimension(imt,jmt)::kmt
!---------------------------------------------------------------------
! zb(1:km)=-zkp(2:km+1)
! dzp(1:km)=zkp(1:km)-zkp(2:km+1)
!
! 0        ====================================== zb(0)
! 1        ------------------- dzp(1)       zm(1)
! 1        ====================================== zb(1)
!          -------------------
!          ======================================
! kmmix    ------------------- dzp(kmmix)   zm(kmmix)
! kmmix    ====================================== zb(kmmix) -> ze
! kmmix+1  ------------------- dzp(kmmix+1) zm(kmmix+1)
! kmmix+1  ====================================== zb(kmmix+1)
!          -------------------
!          ======================================
! km       ------------------- dzp(km)      zm(km)
! km       ====================================== zb(km)
!---------------------------------------------------------------------

contains
!----------------------------------------------------
SUBROUTINE GRIDS_PT  
      zb(1)=dzp(1)
      do k=2,km
        zb(k)=zb(k-1)+dzp(k)
      enddo
      do k=1,km
        zm(k)=-zkt(k)
      enddo
!------------------------------------------------
    kmt=0
     do j=1,jmt
       do i=1,imt
         do k=1,km
        if(vit(i,j,k)>0.5) then
         kmt(i,j)=kmt(i,j)+1
        endif
         enddo
       enddo
      enddo
!lyc-----------------------------------------
!calculate the volume of sea water under 1500m
   if(mytid==0) then
     vsea2000=0.0
     do k=kb0,km
       do j=1,jmt_global
         do i=2,imm_global
         if(vit_global(i,j,k)<0.5) cycle
          vsea2000=vsea2000+dxdyt_global(j)*dzp(k)
         enddo
       enddo
     enddo
   endif
!---------------------------------------------------
    do k=1,km
       sdxdy(k)=0.0
       sdxdy0(k)=0.0
       do j=jsm,jem
         do i=2,imm
           if(vit(i,j,k) < 0.5) cycle
           sdxdy(k)=sdxdy(k)+dxdyt(j)
         enddo
       enddo
    enddo

#ifdef SPMD
      call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_reduce(sdxdy,sdxdy0,km,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)
#endif
!---------------------------------------------------------------------
!Fe_scav_prof=0.2(yr-1)*[1+200exp(-z/250m)]
     do k=1,km
     Fe_scav_prof(k)=0.2/86400.0/365.0*(1.0+200.0*exp(-zm(k)/250))
     enddo
     ReFe=0.0
     do k=kmmix+1,km
     ReFe(k)        =1.0-(zb(k)/zb(k-1))**(POP_k) 
     enddo
!---------------------------------------------------------------------
  RETURN
  END SUBROUTINE GRIDS_PT

!--------------------------------------------------------------------
  END MODULE GRIDS_PT_MOD


