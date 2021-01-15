!-----------------------------------------------------------------------------
!   Processing some variables from flux coupler.
!-----------------------------------------------------------------------------
!
!        By Yongqiang Yu , 16 Apr., 1999
!
!
!
      SUBROUTINE post_cpl

#include <def-undef.h>
use precision_mod
use dyn_mod
use work_mod

use param_mod
use pconst_mod
use tracer_mod
use forc_mod
use buf_mod
use control_mod
use shr_sys_mod
!use work_mod,only : wkj
use output_mod,only:spval
#ifdef SPMD
use msg_mod,only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif


!
      implicit none
      real(r8),dimension(:,:),allocatable::tmp_su,tmp_sv
!JJB20140830 
      REAL(r8)    :: ERR0,FRESH01,FRESH011,ERR1,ERR2,ERR3,ERR4,ERR5,ERR6,ERR7,sss2  
      REAL(r8), ALLOCATABLE :: FRESH_EN1(:,:)
!JJB20140830
!
        tsf=0.0D0
        ssf=0.0D0
        swv=0.0D0
        su=0.0D0
        sv=0.0D0
        mius=0.0D0
        fresh1=0.0D0
        err0=0.0D0
        err2=0.0D0
        err3=0.0D0
        err5=0.0D0
!open(5,file='out.txt3')
!open(6,file='out.txt4')

!$OMP PARALLEL DO PRIVATE (i,j)
        do j=jsm,jmm
        do i= 2,imm
    fresh1(i,j)=vit(i,j,1)*(prec(i,j)+evap(i,j)+ meltw(i,j)+roff(i,j))
    sss1(I,J)=vit(i,j,1)*(ATB (I,J,1,2)*1000.0+35.0)*1.0e-3
     end do
      end do
!JJB20140830   MAKE SURE GLOBAL AREA-MEAN FRESH FLUX IS ZERO
!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:ERROR2)  
         DO J = jsm,jmm
         DO I = 2,imm
         ERR0 = ERR0 + DYT(J)* SINT (J)* fresh1(I,J)*vit(I,J,1)
         ERR2 = ERR2 + DYT(J)* SINT (J)*sss1(I,J)*vit(I,J,1)
         END DO
         END DO
#ifdef SPMD 

      call mpi_reduce(ERR0,ERR1,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_bcast(ERR1,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_reduce(ERR2,ERR4,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_bcast(ERR4,1,MPI_PR,0,mpi_comm_ocn,ierr)
      FRESH01 = - ERR1/ ASEA
      sss2=ERR4/ASEA
#else
       FRESH01 = - ERR0/ ASEA
       sss2=ERR2/ASEA
#endif

!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JMT
         DO I = 1,IMT
            FRESH1(I,J)= (FRESH1(I,J) + FRESH01)* VIT (I,J,1)
         END DO
       END DO


!$OMP PARALLEL DO PRIVATE (i,j)
        do j=jsm,jmm
        do i= 2,imm
    fresh11(i,j)=-vit(i,j,1)*(prec(i,j)+evap(i,j)+ meltw(i,j)+roff(i,j))*((ATB (I,J,1,2)-sss2)*1000.0+35.0)*1.0e-3*OD0
    u3(I,J)=vit(i,j,1)*(SQRT(duu10n(I,J)))**3  
 !  write(5,*) fresh11(56,:)
      end do
      end do
!JJB20140830   MAKE SURE GLOBAL AREA-MEAN FRESH FLUX IS ZERO
!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:ERROR2)  
         DO J = jsm,jmm
         DO I = 2,imm
         ERR5 = ERR5 + DYT(J)* SINT (J)*fresh11(I,J)*vit(I,J,1)   
         ERR3 = ERR3 + DYT(J)* SINT (J)*u3(I,J)*vit(I,J,1)
         END DO
         END DO
#ifdef SPMD   
      call mpi_reduce(ERR5,ERR7,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_bcast(ERR7,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_reduce(ERR3,ERR6,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_bcast(ERR6,1,MPI_PR,0,mpi_comm_ocn,ierr)
       FRESH011 = ERR7/(sss2*ERR6)
  !    write(6,*) fresh011,ERR6,ERR7,ASEA,sss2,fresh01
#else
       FRESH011 =  ERR5/(ERR3*sss2)
#endif
!JJB20140830  
!jjb 20140713  miu*S+lamda*S

!open(3,file='out.txt1')
!open(4,file='out.txt2')
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = 1,JMT
               DO I = 2,IMT
! mius(I,J) = -(0.002256*(SQRT(duu10n))**3.*(ATB (I,J,1,2)*1000.0+35.0)+0.11107)*(1.0-seaice(i,j))*0.01/365./86400/1000.0*VIT(I,J,1)
 
  mius(i,j)=-FRESH011*(SQRT(duu10n(i,j)))**3.*(ATB(I,J,1,2)*1000.0+35.0)*(1.0-seaice(i,j))*VIT(I,J,1)*1.0e-3
            
!write(3,*) mius(56,:) 
           END DO
            END DO
!jjb 20140713  

!$OMP PARALLEL DO PRIVATE (i,j)
        do j=1,jmt
!        do i= 2,imt-1 ! 1,imt !LPF 20120822
        do i= 1,imt
           ! net heat flux
           TSF(i,j) = vit(i,j,1)*(lat1(i,j)+sen(i,j)+lwup(i,j)+lwdn(i,j )+netsw(i,j)+melth(i,j)) *OD0CP  
           ! net solar radiation
           SWV(i,j) = vit(i,j,1)*netsw(i,j) 
           ! none solar radiation !for BUOY
           NSWV(i,j)= vit(i,j,1)*(lat1(i,j)+sen(i,j)+lwup(i,j)+lwdn(i,j) +melth(i,j))
          !SSF(i,j) =-vit(i,j,1)*(prec(i,j)+evap(i,j)+ meltw(i,j)+roff(i,j))  *34.7*1.0e-3/DZP(1)*OD0                                     ! P+E+melting !linpf 25->DZP(1)
         
  !jjb 20140713  
        ! SSF(i,j) =-vit(i,j,1)*(prec(i,j)+evap(i,j)+ meltw(i,j)+roff(i,j))*(ATB (I,J,1,2)*1000.0+35.0)*1.0e-3/DZP(1)*OD0+mius(i,j)*vit(i,j,1)/DZP(1)*OD0    

    SSF(i,j) =-vit(i,j,1)*FRESH1(I,J)*(ATB (I,J,1,2)*1000.0+35.0)*1.0e-3/DZP(1)*OD0+vit(i,j,1)*mius(i,j)/DZP(1)    
   !  write(4,*) SSF(56,:)                      
   ! jjb 20140713  
        end do
        end do
!wangty bug
!wangty bug
        call exchange_2d(tsf,1,1)
        call exchange_2d(swv,1,1)
        call exchange_2d(ssf,1,1)
        call exchange_2d(roff,1,1)

!for ocn output
        lthf  = lat1        !latent flux
        sshf  = sen         !sensible flux
        lwv   = lwup+lwdn   !long wave flux
        fresh = ssf
        runoff= roff

        where(vit(:,:,1)<0.5) tsf=spval
        where(vit(:,:,1)<0.5) swv=spval
        where(vit(:,:,1)<0.5) nswv=spval
        where(vit(:,:,1)<0.5) ssf=spval

        where(vit(:,:,1)<0.5) lthf=spval
        where(vit(:,:,1)<0.5) sshf=spval
        where(vit(:,:,1)<0.5) lwv=spval
        where(vit(:,:,1)<0.5) fresh=spval
        where(vit(:,:,1)<0.5) runoff=spval
!
#ifdef USE_OCN_CARBON
        call exchange_2d(pco2)
#endif         
!
! surface stress
!$OMP PARALLEL DO PRIVATE (i,j)
        do j=1,jmt-1
        do i=2,imt
           SU(i,j) = (taux(i-1,j)+taux(i,j)+taux(i-1,j+1)+taux(i,j+1))*0.25
        end do
        end do
!$OMP PARALLEL DO PRIVATE (i,j)
        do j=2,jmt-1
        do i=2,imt
            SV(i,j) = (tauy(i-1,j)+tauy(i,j)+tauy(i-1,j+1)+tauy(i,j+1))*0.25 !linpf
        end do
        end do
!
        call exchange_2d(SU,1,1)
        call exchange_2d(SV,1,1)

        where(viv(:,:,1)<0.5) su=spval
        where(viv(:,:,1)<0.5) sv=spval

!calculate USTAR
      DO J = 1,jmt
         DO I = 1,imt
          USTAR(I,J)=sqrt(sqrt(taux(i,j)*taux(i,j)+tauy(i,j)*tauy(i,j))*OD0)*vit(i,j,1) 
         END DO
      END DO        
      
!$OMP PARALLEL DO PRIVATE (i,j)
       do i=1,imt
       do j=1,jmt
          licomqice (i,j)= 0.0D0
       end do
       end do

!      call chk_var2d(taux,"us",1)
!      call chk_var2d(tauy,"vs",1)
!      call chk_var2d(SU,"SU",0)
!      call chk_var2d(SV,"SV",0)
!      call chk_var2d(SWV,"sw",1)
!      call chk_var2d(NSWV,"sw",1)
!      call chk_var2d(USTAR,"US",1)

        return
        end
