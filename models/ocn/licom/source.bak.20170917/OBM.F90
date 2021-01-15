! CVS: $Id: ctrlc.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  SUBROUTINE CTRLC
!========================
! CTRLC
!---------------------------------------------------------------------
!
! purpose: set the control constants used in the CARBON model
!
! author: Zhao Liang@lapc 2004/03/02
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE carbon_mod
      USE cforce_mod
      USE coutput_mod
#ifdef SPMD
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!     
#include <netcdf.inc>
!#ifdef SPMD
!#include <mpif.h>      
!#endif
!
      NAMELIST /controlc/ NSTARTC,yearData,yearStart,yearStore,IO_out
!
#ifdef carbonBio      
! default value of r_cp is 120, r_np is 16,r_o2p is -170, sigma_b is 0.5 tau_b is 100 days, &
!  r_fec is 5.0E-6 (based on Archer&Johnson,2000) 
      NAMELIST /bioconst/ r_cp,r_np,r_fec0,sigma_b,tau_b,r_o2p,dust_data,fe_flux_data
#endif
      
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      WRITE(6,*) 'Begining------CTRLC'
#ifdef SPMD
      ENDIF
#endif
      
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      OPEN(11,FILE='controlc',FORM='formatted')
      REWIND(11)
      READ(11,nml=controlc)
      CLOSE(11)
#ifdef SPMD
      ENDIF
      call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(NSTARTC,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(yearData,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(yearStart,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(yearStore,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(IO_out,1,mpi_integer,0,mpi_comm_ocn,ierr)
#endif      
!
#ifdef carbonBio
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
!read r_cp, r_np, sigma_b,tau_b,r_o2p
      OPEN(11,FILE='bioconst',FORM='formatted')
      REWIND(11)
      READ(11,nml=bioconst)
      CLOSE(11)
#ifdef SPMD
      ENDIF
      call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(r_cp,1,mpi_real8,0,mpi_comm_ocn,ierr)
      call mpi_bcast(r_np,1,mpi_real8,0,mpi_comm_ocn,ierr)
      call mpi_bcast(r_fec0,1,mpi_real8,0,mpi_comm_ocn,ierr)
      call mpi_bcast(sigma_b,1,mpi_real8,0,mpi_comm_ocn,ierr)
      call mpi_bcast(tau_b,1,mpi_real8,0,mpi_comm_ocn,ierr)
      call mpi_bcast(r_o2p,1,mpi_real8,0,mpi_comm_ocn,ierr)
#endif
      r_fec(:,:)=r_fec0
      r_fep(:,:)=r_fec(:,:)*r_cp      
#endif
      
!
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      WRITE(6,*) 'END-----------CTRLC'
#ifdef SPMD
      ENDIF
#endif
      
  END SUBROUTINE CTRLC

! CVS: $Id: inirun_pt.F90,v 2.1 2004/06/10 07:45:17 cvsr/ot Exp $
!------------------------------------------
    SUBROUTINE INIRUN_PT
!========================
!
! purpose: initialization carbon model
!
! author: Zhao Liang@lapc 2004/03/02
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE carbon_mod
      USE cforce_mod
      USE coutput_mod
      USE pconst_mod 
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!lyc 2014.09.17
#ifdef COUP
     USE buf_mod,only:pco2,co2_cpl,pco2up
#endif
!
!---------------------------------------------------------------------
      IMPLICIT NONE
#include <netcdf.inc>      
!
      REAL*4::ddd
      integer::dd,np
      REAL*4,DIMENSION(km)::tactmp,po4tmp
!lyc 2013,07
      real*4,dimension(imt_global,jmt_global,km)::fe_in
      integer::ncid,iret
!
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      WRITE(6,*) 'Begining------INIRUN_PT'
#ifdef SPMD
      ENDIF
#endif
!lyc------------------------------------
!2008.12.26 for the calculation of normalized surface total alkalinity
!!------------------------------------
!    call mpi_bcast(basin,imt*jmt_global,mpi_integer,0,mpi_comm_ocn,ierr)  
!!-----------------------------------------------------------------
!lyc 2014.09.17
#ifdef COUP
    allocate(pco2(imt,jmt),co2_cpl(imt,jmt),pco2up(imt,jmt))
#endif
#ifdef SPMD
    IF(mytid==0) THEN
     pt_io=0.0
        
!$OMP PARALLEL DO PRIVATE (j,i)
   DO j=1,jmt_global
    DO i=1,imt_global
    totup_io(i,j)=0.0
    tpco2o_io(i,j)=0.0
    tdpco2o_io(i,j)=0.0
    ENDDO
   ENDDO
!
!$OMP PARALLEL DO PRIVATE (k,j,i)
  DO k=1,km
   DO j=1,jmt_global
    DO i=1,imt_global
      tocaco3_io(i,j,k)=0.0
      toa0_io(i,j,k)=0.0
    ENDDO
   ENDDO
  ENDDO

    ENDIF
#endif
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j=1,jmt
        DO i=1,imt
          ssfc(i,j)=0.0
          sge(i,j)=0.0
        ENDDO
      ENDDO
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
!C$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k=1,km
        DO j=1,jmt
          DO i=1,imt
            ccmon(i,j,k)=0.0
#ifdef carbonBio
            po4mon(i,j,k)=0.0
            ldocmon(i,j,k)=0.0
            tamon(i,j,k)=0.0
            prodmon(i,j,k)=0.0
            fpopmon(i,j,k)=0.0
            pldocmon(i,j,k)=0.0
            remimon(i,j,k)=0.0
            jpopmon(i,j,k)=0.0
            caco3mon(i,j,k)=0.0
	    o2mon(i,j,k)=0.
            femon(i,j,k)=0.0
#endif
          ENDDO
        ENDDO
      ENDDO
!
#ifdef carbonBio
      r_cap(1:km)=0.0
      kappa_b(1:km)=0.0
      delta_a(1:km)=0.0
!C$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k=1,km
        DO j=1,jmt
          DO i=1,imt
            a0_b(i,j,k)=0.0
            a1_b(i,j,k)=0.0
            a2_b(i,j,k)=0.0
            b0_b(i,j,k)=0.0
            b1_b(i,j,k)=0.0
            b2_b(i,j,k)=0.0
            c_b(i,j,k)=0.0
!lyc ----------------------------------
           tocaco3(i,j,k)=0.0
           toa0(i,j,k)=0.0
          ENDDO
        ENDDO
      ENDDO
#endif
!--------------------------
          monthR=1
!----------------------
      IF(NSTARTC==1) THEN
!   ------------------------------------------
!     INITIALIZE CARBON VARIABLES 
!   ------------------------------------------
#if (defined carbonBio) || (defined carbonAbio)
#ifdef SPMD
      IF (mytid==0) THEN
#ifdef preindustrial
!     initial TC data         
        OPEN(81,FILE="tcinitial.dat")
        DO k=1,km
          READ(81,*) ddd,tactmp(k)
        ENDDO
        CLOSE(81)
        DO k=1,km
          DO j=1,jmt_global
            DO i=1,imt_global 
#ifdef carbonAbio
              pt_io(i,j,k,1)=tactmp(k)-2000.
#else
              pt_io(i,j,k,1)=tactmp(k)
#endif
            ENDDO
          ENDDO
        ENDDO
#ifdef carbonBio
!-     initial PO4
      OPEN (81,file='po4initial.dat')
       do k=1,km
	   read(81,*) ddd,po4tmp(k)
       enddo
       do  k=1,km
         do j=1,jmt_global
            do i=1,imt_global
            pt_io(i,j,k,2)=po4tmp(k)
	    enddo
	 enddo
       enddo
       CLOSE(81)

!     initial TA data         
!        OPEN(81,FILE="tainitial.dat")
       OPEN(81,FILE="ta-profile.dat")
        DO k=1,km
          READ(81,*) ddd,tactmp(k)
        ENDDO
        CLOSE(81)
        DO k=1,km
          DO j=1,jmt_global
            DO i=1,imt_global 
              pt_io(i,j,k,4)=tactmp(k)
            ENDDO
          ENDDO
        ENDDO
!     initial LDOC data         
!        DO k=1,kmmix
        DO k=1,km
          DO j=1,jmt_global
            DO i=1,imt_global
              pt_io(i,j,k,3)=4.2
            ENDDO
          ENDDO
        ENDDO
!     initial o2 data
        DO k=1,km
          DO j=1,jmt_global
            DO i=1,imt_global 
              pt_io(i,j,k,5)=170.
            ENDDO
          ENDDO
        ENDDO
!     initial Fe data
      iret=nf_open('Initial_Fe.nc',nf_nowrite,ncid)
       call check_err (iret)
   
      iret=nf_get_vara_real(ncid,   4,(/1,1,1/),(/imt_global,jmt_global,km/), fe_in)
      call check_err (iret)
    
      iret=nf_close(ncid)
      call check_err (iret)

       pt_io(:,:,:,6)=fe_in
      
       
#endif
!lyc
	 write(6,*) 'tc initial:'
         write(6,*)pt_io(imt_global/2,jmt_global/2,1:km,1)
	 write(6,*) 'po4 initial:'
         write(6,*)pt_io(imt_global/2,jmt_global/2,1:km,2)
	 write(6,*)'LODC initial:'
	 write(6,*)pt_io(imt_global/2,jmt_global/2,1:km,3)
	 write(6,*) 'ta initial:'
         write(6,*)pt_io(imt_global/2,jmt_global/2,1:km,4)
	 write(6,*) 'fe initial:'
         write(6,*)pt_io(imt_global/2,jmt_global/2,1:km,6)
! write end	
#else
        OPEN(32,FILE='fort.32',FORM='unformatted')
        REWIND(32)
#ifdef carbonAbio 
        READ(32) pt_io,totup_io,tpco2o_io,tdpco2o_io,monthR
#endif     
#ifdef carbonBio 
        read(32) pt_io!,totup_io,tpco2o_io,tdpco2o_io,tocaco3_io,toa0_io,monthR
#endif
        close(32)
!        monthR=1
#endif 
  !fort.32 here is for AnthroCO2 from pre-industrial results
      ENDIF
!------------------------------------------------------------------------------------
        call mpi_bcast(monthR,1,mpi_integer,0,mpi_comm_ocn,ierr)
!        call global_to_local_4d(pt_io,pt,km,nptra)
!lyc 2014.06
        do np=1,nptra
        do k=1,km
        call global_distribute(pt_io(:,:,k,np),pt(:,:,k,np))
        enddo
        enddo
       
!       if(mytid==2)print*,'pt(1) in inirun',pt(:,:,1,1)
!       if(mytid==2)print*,'pt(2) in inirun',pt(:,:,1,6)
!-------------------------------------------------------------------------------------
#else        
!     initial TC data         
      OPEN(81,FILE="tcinital.dat")
      DO k=1,km
        READ(81,*) ddd,tactmp(k)
      ENDDO
      CLOSE(81)
      DO k=1,km
        DO j=1,jmt
          DO i=1,imt 
#ifdef carbonAbio
            pt(i,j,k,1)=tactmp(k)-2000.
#else
            pt(i,j,k,1)=tactmp(K)
#endif
          ENDDO
        ENDDO
      ENDDO
#ifdef carbonBio
!     initial PO4 data         
      OPEN(81,file="po4initial.dat")
      DO k=1,km
        DO j=1,jmt
          READ(81,*) (pt(i,j,k,2),i=1,imt)
        ENDDO
      ENDDO
!     initial TA data         
      CLOSE(81)
      OPEN(81,FILE="tainital.dat")
      DO k=1,km
        READ(81,*) ddd,tactmp(k)
      ENDDO
      CLOSE(81)
      DO k=1,km
        DO j=1,jmt
          DO i=1,imt 
            pt(i,j,k,4)=tactmp(k)
          ENDDO
        ENDDO
      ENDDO
!     initial LDOC data         
      DO k=1,km
        DO j=1,jmt
          DO i=1,imt 
            pt(i,j,k,3)=4.2
          ENDDO
        ENDDO
      ENDDO
!     initial o2 data
        DO k=1,km
          DO j=1,jmt
            DO i=1,imt 
              pt(i,j,k,5)=170.
            ENDDO
          ENDDO
        ENDDO
#endif
#endif
!--------------------------------
      DO m=1,nptra
!$OMP PARALLEL DO PRIVATE (k,j,i)
        DO k=1,km
          DO j=1,jmt
            DO i=1,imt 
             pt(i,j,k,m)=pt(i,j,k,m)*vit(i,j,k)
#ifdef Felimit
             if(pt(i,j,k,6)>1.0E+2) then
             print *,'Initial of fe is error',i_global(i),j_global(j),k,pt(i,j,k,6)
             endif
#endif
            ENDDO
          ENDDO
        ENDDO
      ENDDO
#endif
! -------------------------------   
      DO m=1,nptra
!$OMP PARALLEL DO PRIVATE (k,j,i) 
        DO k=1,km
          DO j=1,jmt
            DO i=1,imt
#if (defined carbonC)|| (defined cfc)            
              pt(i,j,k,m)=0.0
#endif              
              ptb(i,j,k,m)=pt(i,j,k,m)
            ENDDO
          ENDDO
        ENDDO
!      
!$OMP PARALLEL DO PRIVATE (j,i)
        DO j=1,jmt
          DO i=1,imt
            ptb(i,j,0,m) = 0.0
          ENDDO
        ENDDO

      ENDDO
      !
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j=1,jmt
        DO i=1,imt
          totup(i,j)=0.0
          tpco2o(i,j)=0.0
          tdpco2o(i,j)=0.0
        ENDDO
      ENDDO
#endif
!
#ifdef carbonBio
! initial variables used in biological model
#endif
!----------------------------------------
      ELSE
!     --------------------------------------------
!     READ INTERMEDIATE RESULTS (fort.32/fort.31)
!     ---------------------------------------------
#ifdef SPMD
      IF(mytid==0) THEN
        OPEN(32,FILE='fort.32',FORM='unformatted')
        REWIND(32)
#ifdef carbonC14      
        READ(32) pt_io,monthR
#endif      
!
#if(defined carbonC) || (defined carbonAbio) 
        READ(32) pt_io,totup_io,tpco2o_io,tdpco2o_io,monthR
#endif  
!xu for testing-------------------------------------------    
          DO j=1,jmt_global
            DO i=1,imt_global 
              totup_io(i,j)=0.0
          ENDDO
        ENDDO
!-----------------------------------------------
!-----------------------------------------------
#ifdef carbonBio 
! read restart data from file 32 for continuous calculation
        READ(32) pt_io,totup_io,tpco2o_io,tdpco2o_io,tocaco3_io,toa0_io,monthR
!         print*,"pt_io",pt_io(100,50,10,4)
!lyc 2014.06.08 just for one time
#endif      
        CLOSE(32)
      ENDIF

      call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(monthR,1,mpi_integer,0,mpi_comm_ocn,ierr)
!lyc 2014.06    
      do np=1,nptra
      do k=1,km
      call global_distribute(pt_io(:,:,k,np),pt(:,:,k,np))
      enddo
      enddo
!
#if(defined carbonC) || (defined carbonAbio)
      call global_distribute(totup_io,totup)
      call global_distribute(tpco2o_io,tpco2o)
      call global_distribute(tdpco2o_io,tdpco2o)
#endif      
!
#ifdef carbonBio
! special restart data of biological model
      call global_distribute(totup_io,totup)
      call global_distribute(tpco2o_io,tpco2o)
      call global_distribute(tdpco2o_io,tdpco2o)
      do k=1,km
      call global_distribute(tocaco3_io(1,1,k),tocaco3(1,1,k))
      call global_distribute(toa0_io(1,1,k),toa0(1,1,k))
      enddo
#endif      
#else
      OPEN(32,FILE='fort.32',FORM='unformatted')
      REWIND(32)
#ifdef carbonC14      
      READ(32) pt,monthR
#endif      
!
#if (defined carbonC) || (defined carbonAbio) 
      READ(32) pt,totup,tpco2o,tdpco2o,monthR
#endif      
!
#ifdef carbonBio      
      READ(32) pt,totup,tpco2o,tdpco2o,tocaco3,t,toa0,omonthR
#endif      
      CLOSE(32)
#endif      
!--------------------------------------
      DO m=1,nptra
!$OMP PARALLEL DO PRIVATE (k,j,i) 
        DO k=1,km
          DO j=1,jmt
            DO i=1,imt
              ptb(i,j,k,m)=pt(i,j,k,m)
            ENDDO
          ENDDO
        ENDDO
!      
!$OMP PARALLEL DO PRIVATE (j,i)
        DO j=1,jmt
          DO i=1,imt
            ptb(i,j,0,m) = 0.0
          ENDDO
        ENDDO
      ENDDO
!      
      ENDIF
!
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      WRITE(6,*) 'END-----------INIRUN_PT'
#ifdef SPMD
      ENDIF
#endif
!lyc20121010      deallocate(pt_r4,totup_r4,tpco2o_r4,tdpco2o_r4,tocaco3_r4,toa0_r4)
!------------------------------------------------------------------------------------
      END SUBROUTINE INIRUN_PT

! CVS: $Id: cforce.F90,v 2.2 2004/06/13 12:14:56 cvsroot Exp $
  SUBROUTINE CFORCE
!========================
! CFORCE
!---------------------------------------------------------------------
!
! purpose: set forcing data used in CARBON cycle
!
! author: Zhao Liang@lapc 2004/03/01
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE carbon_mod
      USE pconst_mod,ONLY:vit,vit_global
      USE cforce_mod
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------------
      IMPLICIT NONE
#include <netcdf.inc>      
!#ifdef SPMD      
!#include <mpif.h>
!#endif      
!      
      REAL::tmp1
      integer,dimension(4)::start4,count4
      integer::iret,ncid,t
      real*4,dimension(:,:,:,:),allocatable:: po4_r4
!wk add r8=============     
 real*4,allocatable,dimension(:)::csgn_r4,csg_r4,csgx_r4
 real*4,dimension(km,12)::taobs_r4
 real*4,dimension(imt_global,jmt_global,12)::winds_global_r4
 real*4,dimension(imt_global,jmt_global,12)::pressure_r4,fe_force
 real,dimension(imt_global,jmt_global,12)::fe_flux_io,dust_flux_io

!wk add r8=============  
      allocate(po4_r4(imt_global,jmt_global,km,12))
!lyc
#ifdef cfc
       real,dimension(ny):: cfcyear,vnerror,vserror
#endif

#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      WRITE(6,*) 'Begining------CFORCE'
#ifdef SPMD
      ENDIF
#endif
      
!cm080427--------
!#ifdef carbonC14
#if (defined carbonC14) && (!defined preindustrial)
!cm080427--------
#ifdef SPMD
      IF(mytid==0) THEN
#endif
      OPEN(11,FILE='bom14.fmt',FORM='formatted')
      READ(11,*) ny
#ifdef SPMD
      ENDIF
      CALL mpi_barrier(mpi_comm_ocn,ierr)
      CALL mpi_bcast(ny,1,mpi_integer,0,mpi_comm_ocn,ierr)
#endif
      ALLOCATE(kyear(ny),boml(ny),bomm(ny),bomh(ny))
      
#ifdef SPMD
      IF(mytid==0) THEN
#endif
      DO i=1,ny
        READ(11,*) kyear(i),boml(i),bomm(i),bomh(i)
      ENDDO
      CLOSE(11)
#ifdef SPMD
      ENDIF
      CALL mpi_barrier(mpi_comm_ocn,ierr)
      CALL mpi_bcast(kyear,ny,mpi_integer,0,mpi_comm_ocn,ierr)
      CALL mpi_bcast(boml,ny,mpi_real8,0,mpi_comm_ocn,ierr)
      CALL mpi_bcast(bomm,ny,mpi_real8,0,mpi_comm_ocn,ierr)
      CALL mpi_bcast(bomh,ny,mpi_real8,0,mpi_comm_ocn,ierr)
#endif
#endif
      
!#if (defined carbonC) ||(defined carbonAbio) || (defined carbonBio)||(defined cfc)||(defined nc14wind)
!!---------------------------------------------------------------------
!!Read wind speed
!! A.nnual mean wind spead for the calulation of exchange coefficeint sge
!! m is month number(1-12)
!!---------------------------------------------------------------------
#if (defined FRC_CORE)||(defined COUP)
       continue
#else
#ifdef SPMD
      IF(mytid==0) THEN
!!wk===================================================
       winds_global(:,:,1:12)=10.0
       open(88,file='winds.dat',form='unformatted',convert='big_endian')
       do m=1,12
       read(88) winds_global_r4(:,:,m)
       enddo
       do m=1,12
         do j=1,jmt_global
           do i=1,imt_global
            winds_global(i,j,m)=dble(winds_global_r4(i,j,m))
           enddo
         enddo
       enddo
!!wk===================================================
!      !print*,"winds=",winds_global
      ENDIF
     do m=1,12
      call global_distribute(winds_global(:,:,m),winds(:,:,m))
     enddo
#else
      winds(:,:,1:12)=10.0
      open(88,file='winds.dat',form='unformatted',convert='big_endian')
      do m=1,12
      read(88) winds(:,:,m)
      enddo
      close(88)
#endif      
!$OMP PARALLEL DO PRIVATE (j,i,m)        
        DO j=1,jmt
          DO i=1,imt
            DO m=1,12
              winds(i,j,m)=winds(i,j,m)*vit(i,j,1)
            ENDDO
          ENDDO
        ENDDO
#endif
!cm090330-------------------------------------------------------------
!read pressure
!
#if (defined FRC_CORE)||(defined COUP)
       continue
#else
#ifdef SPMD
      IF(mytid==0) THEN
!!wk====================================================
        open(88,file='pressure.dat',form='unformatted',convert='big_endian')
        do m=1,12
        read(88) pressure_r4(:,:,m)
        enddo
        do m=1,12
          do j=1,jmt_global
           do i=1,imt_global
             pressure_global(i,j,m)=dble(pressure_r4(i,j,m))
            enddo
          enddo
        enddo

!!wk====================================================
      ENDIF
      do m=1,12
      call global_distribute(pressure_global(1,1,m),pressure(1,1,m))
      enddo
#else
      open(88,file='pressure.dat',form='unformatted',convert='big_endian')
      do m=1,12
      read(88) pressure(:,:,m)
      enddo
      close(88)
#endif      
!$OMP PARALLEL DO PRIVATE (j,i,m)        
        DO j=1,jmt
          DO i=1,imt
            DO m=1,12
              pressure(i,j,m)=pressure(i,j,m)*vit(i,j,1)
            ENDDO
          ENDDO
        ENDDO
#endif
!cm090330-------------------------------------------------------------
#if (defined carbonC) || (defined carbonBio)|| (defined carbonAbio)
!---------------------------------------------------------------------
!    reading atmospheric CO2
!---------------------------------------------------------------------
#ifdef preindustrial
!        OPEN(11,FILE='prepco2atm.fmt',FORM='formatted')
#else
#ifdef SPMD      
        IF(mytid==0) THEN
#endif        
        OPEN(11,FILE='pco2atm.fmt',FORM='formatted')
        READ(11,*) nsg1
#ifdef SPMD
        ENDIF
        call mpi_barrier(mpi_comm_ocn,ierr)
        call mpi_bcast(nsg1,1,mpi_integer,0,mpi_comm_ocn,ierr)
#endif
        allocate(csgn(nsg1),csg(nsg1),csgx(nsg1))
        allocate(csgn_r4(nsg1),csg_r4(nsg1),csgx_r4(nsg1))
#ifdef SPMD
        IF(mytid==0) THEN
#endif
        DO j=1,nsg1
!wk add ============
         ! READ(11,*) csgn(j),csg(j),csgx(nsg1)
          READ(11,*) csgn_r4(j),csg_r4(j),csgx_r4(nsg1)
           csgn(j)=csgn_r4(j)
           csg(j) =csg_r4(j)
           csgx(j)=csgx_r4(nsg1)
!wk add ============
    
        ENDDO
        CLOSE(11)
        deallocate(csgn_r4,csg_r4,csgx_r4)

#ifdef SPMD        
        ENDIF
        call mpi_barrier(mpi_comm_ocn,ierr)
        call mpi_bcast(csgn,nsg1,mpi_real8,0,mpi_comm_ocn,ierr)
        call mpi_bcast(csg,nsg1,mpi_real8,0,mpi_comm_ocn,ierr)
        call mpi_bcast(csgx,nsg1,mpi_real8,0,mpi_comm_ocn,ierr)
#endif        
#endif        
!for 746
#endif
       if(mytid==0) print *, 'atm co2 is ok'


#ifdef carbonBio
#ifdef SPMD
#ifdef murnane1999
      IF(mytid==0) THEN
!----------------------------------------------------
      iret=nf_open('po4.dat',nf_nowrite,ncid)
      call check_err (iret)
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start4(1)=1 ; count4(1)=imt_global
      start4(2)=1 ; count4(2)=jmt_global
      start4(3)=1 ; count4(3)=km
      start4(4)=1 ; count4(4)=12
      iret=nf_get_vara_real(ncid,   5,start4,count4,po4_r4)
      call check_err (iret)
!
	    do m=1,12
	      do k=1,kmmix
		do j=1,jmt_global
		  do i=1,imt_global
		      po4obs_global(i,j,k,m)=dble(po4_r4(i,j,k,m))
!		      if(vit_global(i,j,k)>0.5.and.po4obs_global(i,j,k,m)==0.0) po4obs_global(i,j,k,m)=pt_io(i,j,k,2)
		  enddo
	        enddo
	      enddo
	    enddo
      ENDIF
      do m=1,12
       do k=1,kmmix
      call global_distribute(po4obs_global(:,:,k,m),po4obs(:,:,k,m))
       enddo
      enddo
!lyc
#endif
#ifdef progca
! read forcing data of TA from GEOSEC       
      IF(mytid==0) THEN
!-------cm080814------------------------------
!        OPEN(11,FILE='ta-geosec.dat',FORM='formatted')
        OPEN(11,FILE='ta-profile.dat',FORM='formatted')
!        DO m=1,12
!          READ(11,*)
          DO k=1,km
!wk=======================================
            READ(11,*) tmp1,taobs_r4(k,1)
            taobs(k,1:12)=dble(taobs_r4(k,1))
          ENDDO
!wk=======================================
!        ENDDO
        CLOSE(11)
      ENDIF
      call mpi_bcast(taobs,km*12,mpi_real8,0,mpi_comm_ocn,ierr)
#endif
!     if(mytid==0) print *,'before dust reading'
!lyc 2013,07
!-----read the iron flux-------------------------
      IF(mytid==0) THEN
      IF(Fe_FLUX_DATA) THEN
      iret=nf_open('Fe_flux.nc',nf_nowrite,ncid)
      call check_err (iret)
      
      iret=nf_get_vara_real(ncid, 4,(/1,1,1/),(/imt_global,jmt_global,12/),fe_force)
      call check_err (iret)
      print *,'reading dust flux is ok,ncid',ncid
! fe_force umol/cm2/s
      do t=1,12
       do j=1,jmt_global
         do i=1,imt_global      
!      fe_flux_io(i,j,t)=fe_bioava*fe_force(i,j,t)*10.0/1.025*vit_global(i,j,1)
!for csm1_bgc iron flux 
      print *,'i=',i,'j=',j,'t=',t
      print *,'imt_global=',imt_global,'jmt_global=',jmt_global
      print *,'fe_force(1,1,1)=',fe_force(i,j,t)
      print *,'vit_global(1,1,1)=',vit_global(i,j,1)
     fe_flux_io(i,j,t)=dble(fe_force(i,j,t))*10.0/1.025*vit_global(i,j,1)
      print *,'i=',i,'j=',j,'t=',t
      dust_flux_io(i,j,t)=(1-fe_bioava)*fe_force(i,j,t)*10.0/1.025*vit_global(i,j,1)
         if(fe_flux_io(i,j,t)>1.0D+2) fe_flux_io(i,j,t)=0.0!print*,'i,j,t,fe_flux_io',i,j,t,fe_flux_io(i,j,t)
         if(dust_flux_io(i,j,t)>1.0D+2) dust_flux_io(i,j,t)=0.0!print*,'i,j,t,fe_flux_io',i,j,t,fe_flux_io(i,j,t)
         enddo
       enddo
       enddo
      iret=nf_close(ncid)
      ENDIF

      IF(DUST_DATA) THEN
      iret=nf_open('Dust_flux.nc',nf_nowrite,ncid)
      call check_err (iret)
      
      iret=nf_get_vara_real(ncid, 4,(/1,1,1/),(/imt_global,jmt_global,12/),fe_force)
      call check_err (iret)
!for dust flux (fe_force) units ug/cm2/s
!fe_flux_io umol/kg*m/s
!dust_flux_io ug/kg*m/s
      do t=1,12
       do j=1,jmt_global
         do i=1,imt_global      
      fe_flux_io(i,j,t)=0.035/55.847*fe_bioava*fe_force(i,j,t)*10.0/1.025*vit_global(i,j,1)
      dust_flux_io(i,j,t)=(1-fe_bioava)*fe_force(i,j,t)*10.0/1.025*vit_global(i,j,1)
         enddo
       enddo
       enddo
     
      iret=nf_close(ncid)
      ENDIF
!----------------------------------------------------
      ENDIF
      call mpi_barrier(mpi_comm_ocn,ierr)
      do t=1,12
      call global_distribute(fe_flux_io(:,:,t),fe_flux(:,:,t))
      call global_distribute(dust_flux_io(:,:,t),dust_flux(:,:,t))
      enddo 
       if(mytid==2) print *,'dust flux in mytid==2',fe_flux(:,:,2)
      
#else
#ifdef murnane1999
!----------------------------------------------------
      iret=nf_open('po4.dat',nf_nowrite,ncid)
      call check_err (iret)
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start4(1)=1 ; count4(1)=imt
      start4(2)=1 ; count4(2)=jmt_global
      start4(3)=1 ; count4(3)=kmmix
      start4(4)=1 ; count4(4)=1
      iret=nf_get_vara_real(ncid,   5,start4,count4,po4_r4)
      call check_err (iret)
!
	    do m=1,12
	      do k=1,kmmix
		do j=1,jmt_global
		  do i=1,imt
		      po4obs(i,j,k,m)=po4_r4(i,j,k,1)
		  enddo
	        enddo
	      enddo
	    enddo
#endif
#ifdef progca
! read forcing data of TA from GEOSEC       
!        OPEN(11,FILE='ta-geosec.dat',FORM='formatted')
        OPEN(11,FILE='ta-profile.dat',FORM='formatted')
!        READ(11,*)
        DO m=1,12
!          READ(11,*)
          DO k=1,km
           ! READ(11,*) tmp1,taobs(k,m)
!wk========================================================
            READ(11,*) tmp1,taobs_r4(k,m)
            taobs(k,m)=dble(taobs(k,m))
          ENDDO
        ENDDO
           print*,"taobs=",taobs
!wk========================================================
        CLOSE(11)
!-------cm080814------------------------------
#endif
#ifdef murnane1999
        DO m=1,12
          DO k=1,kmmix  
            DO j=1,jmt
              DO i=1,imt
                po4obs(i,j,k,m)=po4obs(i,j,k,m)*vit(i,j,k)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
              
#endif
#endif
!for SPMD
#endif
!for carbonBio
!lyc
#ifdef cfc
#ifdef SPMD
       if(mytid==0) then
        OPEN(1,FILE='cfc-11.txt')
         DO I=1,ny
           READ(1,*) cfcyear(i),bomn(i),vnerror(i),boms(i),vserror(i)
         ENDDO
        CLOSE(1) 
	endif
      CALL mpi_barrier(mpi_comm_ocn,ierr)
      CALL mpi_bcast(bomn,ny,mpi_real8,0,mpi_comm_ocn,ierr)
      CALL mpi_bcast(boms,ny,mpi_real8,0,mpi_comm_ocn,ierr)
#endif
#endif
 !lyc
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      WRITE(6,*) 'END-----------CFORCE'
#ifdef SPMD
      ENDIF
#endif
      
  END SUBROUTINE CFORCE

! CVS: $Id: ptracer.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
!---------------------------------------------------------------
      SUBROUTINE PTRACER
!========================
!Xu 20121012
! purpose: prediction of the passive tracers of biogeochemsitry model
!
! author: Zhao Liang@lapc 2004/03/03
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod
      USE tracer_mod
      USE dyn_mod
      USE isopyc_mod
      USE pmix_mod
      USE work_mod
      USE carbon_mod
      USE cforce_mod
      USE forc_mod

#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------
      REAL    :: AIDIF,C2DTTS,AA,FAW,FIW,ALF,RNCC,ABC
      REAL    :: wt1,wt2,adv_y,adv_x,adv_z,adv_x1,adv_x2,upsh
      real*8    :: ldoctmp,ldocsum,deltata,deltapo4,deltatc
      real*8    :: tatmp,tasum,tcsum,tctmp,po4sum,po4tmp,tcsumb,tctmpb
      REAL,ALLOCATABLE,DIMENSION(:,:)  ::stf1
      REAL,ALLOCATABLE,DIMENSION(:,:,:)::wkb1,wkc1,wkd1,tf1,uwk1 !for upwell_pt
      REAL,ALLOCATABLE,DIMENSION(:,:)  ::h0f1
      REAL,ALLOCATABLE,DIMENSION(:,:,:)::utf1,vtf1,utl1,vtl1,wst
      REAL,DIMENSION(IMT,JMT,KM)::TEST
!lyc 2014.09.18
      REAL,PARAMETER::fil_latp=66.0
!------------------------------------------
#if (defined carbonC)||(defined carbonC14)
      real*8::vseac,vseap
#endif
!lyc201209
     integer::np
!------------------------------------------------------------
!     SET LOCAL CONSTANT
!----------------------------------------------------------
 
      allocate(stf1(imt,jmt),tf1(imt,jmt,km))
      allocate(wkb1(imt,jmt,km),wkc1(imt,jmt,km),wkd1(imt,jmt,km))
      allocate(h0f1(imt,jmt))
      allocate(utf1(imt,jmt,km),vtf1(imt,jmt,km),utl1(imt,jmt,km),vtl1(imt,jmt,km))
      allocate(uwk1(imt,jmt,km),wst(imt,jmt,kmp1)) ! for upwell_pt

#if (defined ISO)
      AIDIF = 1.0 
#else
      AIDIF = 0.0
#endif
 
      RNCC = 1.0/ FLOAT (NCC)
!----------------------------
!cm090302----
!      IF (ISP >= 1)THEN
      IF (mod(ISP,15)/=0) THEN
!cm090302----
         C2DTTS = DTS *2.0
         AA = 0.5
      ELSE
         C2DTTS = DTS
         AA = 0.0
      END IF

     VTL1=0.0
     TF1=0.0
 
!-------------------------------------------------------------
!     PREPARATION FOR VERTICAL ADVECTIVE TERM
!--------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            H0F1 (I,J)= H0F (I,J)* ONBC
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            STF1 (I,J)= AA * H0F1 (I,J) + (1.0- AA)* H0L (I,J)
        END DO
      END DO

    
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTF1 (I,J,K)= UTF (I,J,K)* ONCC
               VTF1 (I,J,K)= VTF (I,J,K)* ONCC
            END DO
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKD1 (I,J,K)= AA * UTF1 (I,J,K) + (1.0- AA)* UTL (I,J,K)
               uwk1 (i,j,k)= wkd1(i,j,k) !for upwell_pt 
               WKC1 (I,J,K)= AA * VTF1 (I,J,K) + (1.0- AA)* VTL (I,J,K)
               if(isnan(wkc1(i,j,k)))then
               print *,'wkc1(v) is error,i,j,k,vtf1,vtl,aa,r2b(j):',i_global(i),j_global(j),k,vtf1(i,j,k),VTF(i,j,k),vtl(i,j,k),aa,r2b(j)
               stop
	       endif
#if ( defined SMAG)
               UTL1 (I,J,K)= AA * UTF1 (I,J,K) + (1.0- AA)* UTL (I,J,K)
               VTL1 (I,J,K)= AA * VTF1 (I,J,K) + (1.0- AA)* VTL (I,J,K)
#endif
            END DO
         END DO
      END DO

!      CALL UPWELL (WKD1,WKC1,STF1)
!      wst=ws
!---------------------------------------------------------------------
!     PREPARATION FOR HORIZONAL ADVECTIVE TERM 
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 1,IMT
               UTL1 (I,J,K)= 0.25* OTX (J)* (WKD1 (I,J,K) + WKD1 (I,J -1,K))
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               WKD1 (I,J,K)= R2A (J)* (WKC1 (I,J,K) + WKC1 (I +1,J,K))
               WKB1 (I,J,K)= R2B (J)* (WKC1 (I,J -1,K) + WKC1 (I +1,J -1,K))
            END DO
         END DO
      END DO
!lyc test

!       do k=1,km
!         do j=1,jmt
!           do i=1,imt
!             if(vit(i,j,k)<0.5) cycle
!              if(isnan(pt(i,j,k,1)).or.isnan(pt(i,j,k,4))) then
!              print *, 'i,j,k,pt',i,j_global(j),k,pt(i,j,k,1),c_b(i,j,k)
!              endif
!            enddo
!           enddo
!         enddo
!       if(mytid==3) then
!        print*,'i,j,k,vit(i,j,1),pt(i,j,k,:)',101,j_global(5),1,vit(101,5,1),pt(101,5,1,:)
!       endif
!-------------------------------------------
!lyc for upwell_pt
!   if(mytid==0) print *,'call upwell_pt'
!---------------------------------------------------
        CALL UPWELL_PT(uwk1,wkc1,utl1,wkd1,wkb1,wst)
!--------------------------------------------------
!   if(mytid==0) print *,'endcall upwell_pt'
 
!-----------------------------------------------------------------------
!     PREPARATION FOR ISOPYCNAL DIFFUSION & ADVECTION
!-----------------------------------------------------------------------
#if (defined ISO)
!     Calculate K1,K2 and K3
!------------------------------
            CALL ISOPYC
!------------------------------
#endif
!--------------------------------------------
!     COMPUTING DIFFUSION COEFFICIENT
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
#if (defined ISO)

               WKC1 (I,J,K) = AHV + AHISOP * K3 (I,K,J,3)
               if(k>3.and.vit(i,j,k)>0.5) then
               if(pt(i,j,k,1)>2800.0.or.pt(i,j,k,1)<1500.0) then
               WKC1(i,j,k-1)=2*AHV+AHISOP*K3(I,K,J,3)
               WKC1(i,j,k)=2*AHV+AHISOP*K3(I,K,J,3)
              endif
              endif
!#ifdef CANUTO
!               WKC1(I,J,K)  = AKT(I,J,K,2)+AHISOP*K3(I,K,J,3)
!#endif
!               WKC1 (I,J,K) = AHV(i,j,k) + AHISOP * K3 (I,K,J,3)
#else
               WKC1 (I,J,K) = AHV
!               WKC1 (I,J,K) = AHV(i,j,k)
#endif
            END DO
         END DO
      END DO
 
 
!     AT LOW LATITUDE, DIFFUSION DEPENDs ON RICHARDSON NUMBER
 
!$OMP PARALLEL DO PRIVATE (K,J,I,ABC,ALF)
      DO K = 1,KMM1
      do j=  jsm,jem
#ifdef SPMD
         if (j_global(j)>=rtst.and.j_global(j)<=rtend) then
#else
         if (j>=rtst.and.j<=rtend) then
#endif
            DO I = 1,IMT
               IF (rit (i,j,k) < 0.0) THEN
                  ABC = diff_cbt_limit
               ELSE
                  ALF = 1.0/ (1.0+5.0* rit (i,j,k)* RNCC)
                  ABC = fricmx * ALF **3+ diff_cbt_back
               END IF
               IF (k == 1.AND.ABC < wndmix) ABC = wndmix
#if (defined ISO)
               WKC1 (I,J,K) = ABC + AHISOP * K3 (I,K,J,3)
#else
               WKC1 (I,J,K) = ABC
#endif
            END DO
         end if
      end do
      end do
!test
!       wkc1(:,:,:)=0.0
!-----------------------------------------------------------
!     prepration of the calculation of biosource
!     calculate A0_b, A1_b, A2_b, B0_b, B1_b, B2_b and C_b
!     and delta_a(km), kappa_a(km)
!--------------------------------------------------------------
#ifdef carbonBio
#ifdef printcall
#ifdef SPMD
            print*,"call readybio in ptracer, mytid=",mytid
#else
            print*,"call readybio in ptracer"
#endif
#endif
!----------------------------------           
            CALL READYBIO(C2DTTS)
!----------------------------------
#endif

!------------------------------------------------------------
!     SOLVE FOR ONE PASSIVE TRACER AT A TIME
!-------------------------------------------------------------
!     NPTRA = 1 => carbon or C14, 2 => PO4, 3 => LDOC, 4 => TA
 
      DO NP = 1,NPTRA
!
!lyc 2011.02.15
!to deal with the DIC and TA before the calculation of the advective term
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM : ZONAL COMPONENT AND MERIDIONAL COMPONENT
!11111111111111-------------------------------------------------------------
#if (defined mom_xu_pt)
!-----------------------------------
                 upsh=0.8
!-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z,wt1,wt2)
         DO K = 2,km-1
            DO J = JSM,JEM
               DO I = 2,IMM
                  adv_x1=utl1(i,j,k)*(PT (I  ,J,K,NP) + PT (I-1,J,K,NP))&
                       +upsh*abs(utl1(i,j,k))*(PT(I-1,J,K,NP)-PT(I,J,K,NP))
                  adv_x2=utl1(i+1,j,k)*(PT (I+1,J,K,NP) + PT (I,J,K,NP))&
                       +upsh*abs(utl1(i+1,j,k))*(PT(I,J,K,NP)-PT(I+1,J,K,NP))
                  adv_x = - (adv_x2-adv_x1)
                  adv_y=-(WKD1(I,J,K)*(PT(I,J+1,K,NP)+PT(I,J,K,NP))&
                      +upsh*abs(wkd1(i,j,k))*(PT(I,J,K,NP)-PT(I,J+1,K,NP))&
                      -WKB1(I,J,K)*(PT(I,J,K,NP)+PT(I,J-1,K,NP))&
                      -upsh*abs(wkb1(i,j,k))*(PT(I,J-1,K,NP)-PT(I,J,K,NP)))
                  wt1= WST(I,J,K)*(PT(I,J,K-1,NP)+PT(I,J,K,NP))&
                      +upsh*abs(wst(i,j,k))*(PT(I,J,K,NP)-PT(I,J,K-1,NP))
                  wt2= WST(I,J,K+1)*(PT(I,J,K,NP)+PT(I,J,K+1,NP))&
                     +upsh*abs(wst(i,j,k+1))*(PT(I,J,K+1,NP)-PT(I,J,K,NP))
                  adv_z=-0.5*(wt1-wt2)*ODZP(K)
                  TF1(I,J,K)= adv_x+adv_y+adv_z
                 if(isnan(tf1(i,j,k))) print *,'tf1 is error:i,j,k,np,adv_x,adv_y,adv_z,wkd1,wkb1',i_global(i),j_global(j),k,np,adv_x,adv_y,adv_z,pt(i,j,k,np),wkd1(i,j,k),wkb1(i,j,k),pt(i,j-1,k,np),vit(i,j-1,k)
               END DO
            END DO
         END DO
!--------------------------------surface level--------------
!$OMP PARALLEL DO PRIVATE (J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z)
        DO J=JSM,JEM
           DO I = 2,IMM
              adv_x1=(PT (I,J,1,NP) + PT (I-1,J,1,NP))* UTL1 (I,J,1)&
                   +upsh*abs(utl1(i,j,1))*(PT(I-1,J,1,NP)-PT(I,J,1,NP))
              adv_x2=(PT (I+1,J,1,NP) + PT (I,J,1,NP))* UTL1 (I+1,J,1)&
                   +upsh*abs(utl1(i+1,j,1))*(PT(I,J,1,NP)-PT(I+1,J,1,NP))
              adv_x = - (adv_x2-adv_x1)
              adv_y=-(WKD1(I,J,1)*(PT(I,J+1,1,NP)+PT(I,J,1,NP))&
                    +upsh*abs(wkd1(i,j,1))*(PT(I,J,1,NP)-PT(I,J+1,1,NP))&
                   -WKB1(I,J,1)*(PT(I,J,1,NP)+PT(I,J-1,1,NP))&
                   -upsh*abs(wkb1(i,j,1))*(PT(I,J-1,1,NP)-PT(I,J,1,NP)))
              adv_z= 0.5*WST(I,J,2)*(PT(I,J,1,NP)+PT(I,J,2,NP))*odzp(1)&
                   +0.5*upsh*abs(wst(i,j,2))*(PT(I,J,2,NP)-PT(I,J,1,NP))*odzp(1)
             TF1 (I,J,1)= adv_x+adv_y+adv_z
                 if(isnan(tf1(i,j,1))) print *,'tf1at k=1 is error:i,j,np,adv_x,adv_y,adv_z,pt(i,j,1,np)',i_global(i),j_global(j),np,adv_x,adv_y,adv_z,pt(i,j,1,np),wkd1(i,j,1),wkb1(i,j,1),pt(i,j-1,1,np)
           END DO
        END DO
!------------- --------bottom level---------
!$OMP PARALLEL DO PRIVATE (J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z)
        DO J=JSM,JEM
           DO I = 2,IMM
              adv_x1=(PT (I  ,J,KM,NP) + PT (I-1,J,KM,NP))* UTL1 (I,J,KM)&
                   +upsh*abs(utl1(i,j,km))*(PT(I-1,J,KM,NP)-PT(I,J,KM,NP))
              adv_x2=(PT (I+1,J,km,NP) + PT (I,J,km,NP))* UTL1 (I+1,J,KM)&
                  +upsh*abs(utl1(i+1,j,km))*(PT(i,j,KM,NP)-PT(I+1,J,KM,NP))
              adv_x = - (adv_x2-adv_x1)
              adv_y=-(WKD1(I,J,km)*(PT(I,J+1,KM,NP)+PT(I,J,KM,NP))&
                     +upsh*abs(wkd1(i,j,km))*(PT(I,J,KM,NP)-PT(I,J+1,KM,NP))&
                    -WKB1(I,J,km)*(PT(I,J,KM,NP)+PT(I,J-1,KM,NP))&
                 -upsh*abs(wkb1(i,j,km))*(PT(I,J-1,KM,NP)-PT(I,J,KM,NP)))
              adv_z= -0.5*WST(I,J,km)*(PT(I,J,KM-1,NP)+PT(I,J,KM,NP))*odzp(km)&
                   -0.5*upsh*abs(wst(i,j,km))*(PT(I,J,KM,NP)-PT(I,J,KM-1,np))*odzp(km)
              TF1 (I,J,km)= adv_x+adv_y+adv_z
           END DO
        END DO
!1111111111111111-------------------------------------------------------
! original calculation method
!
#else
!$OMP PARALLEL DO PRIVATE (K,J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z,wt1,wt2)
         DO K = 2,KM-1
            DO J = JSM,JEM
               DO I = 2,IMM
!lyc           
                  adv_x1=(PT (I  ,J,K,NP) - PT (I-1,J,K,NP))* UTL1 (I,J,K)
                  adv_x2=(PT (I+1,J,K,NP) - PT (I,J,K,NP))* UTL1 (I+1,J,K)
                  adv_x = - (adv_x2+adv_x1)
                  adv_y=-(WKD1(I,J,K)*(PT(I,J+1,K,NP)-PT(I,J,K,NP))+WKB1(I,J,K)*(PT(I,J,K,NP)-PT(I,J-1,K,NP)))
                  wt1= WST(I,J,K)*(PT(I,J,K-1,NP)-PT(I,J,K,NP))
                  wt2= WST(I,J,K+1)*(PT(I,J,K,NP)-PT(I,J,K+1,NP))
                  adv_z=-0.5*(wt1+wt2)*ODZP(K)
                  TF1 (I,J,K)= adv_x+adv_y+adv_z
               END DO
            END DO
         END DO
!$OMP PARALLEL DO PRIVATE (J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z)
        DO J=JSM,JEM
           DO I = 2,IMM
!lyc
              adv_x1=(PT (I,J,1,NP) - PT (I-1,J,1,NP))* UTL1 (I,J,1)
              adv_x2=(PT (I+1,J,1,NP) - PT (I,J,1,NP))* UTL1 (I+1,J,1)
              adv_x = - (adv_x2+adv_x1)
              adv_y=-(WKD1(I,J,1)*(PT(I,J+1,1,NP)-PT(I,J,1,NP))+WKB1(I,J,1)*(PT(I,J,1,NP)-PT(I,J-1,1,NP)))
              adv_z= -0.5*WST(I,J,2)*(PT(I,J,1,NP)-PT(I,J,2,NP))*odzp(1)
              TF1 (I,J,1)= adv_x+adv_y+adv_z
                 if(isnan(tf1(i,j,1))) then
                 print *,'tf1 before isoflux is error:i,j,k,np',i_global(i),j_global(j),np,pt(i,j,1,np),ptb(i,j,1,np),adv_x,adv_y,adv_z
                 stop
                 endif
           END DO
        END DO
!$OMP PARALLEL DO PRIVATE (J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z)
        DO J=JSM,JEM
           DO I = 2,IMM
!lyc	      
              adv_x1=(PT (I  ,J,km,NP) - PT (I-1,J,km,NP))* UTL1 (I,J,km)
              adv_x2=(PT (I+1,J,km,NP) - PT (I,J,km,NP))* UTL1 (I+1,J,km)
              adv_x = - (adv_x2+adv_x1)
              adv_y=-(WKD1(I,J,km)*(PT(I,J+1,km,NP)-PT(I,J,km,NP))+WKB1(I,J,km)*(PT(I,J,km,NP)-PT(I,J-1,km,NP)))
              adv_z= -0.5*WST(I,J,km)*(PT(I,J,km-1,NP)-PT(I,J,km,NP))*odzp(km)
              TF1 (I,J,km)= adv_x+adv_y+adv_z
           END DO
        END DO
!1111111111-----------------------------------
#endif
!-----------------
!lyc 2011.02.15
!
!-----------------------------------------------------------------------
!     COMPUTE THE ISOPYCNAL/DIPYCNAL MIXING
!-----------------------------------------------------------------------
!     XZ AND YZ ISOPYCNAL DIFFUSIVE FLUX ARE SOLVED EXPLICITLY;
!     WHILE ZZ COMPONENT WILL BE SOLVED IMPLICITLY.
 
!iso2222---------------------------------------- 
!
#if (defined ISO)
!     Calculate XZ and YZ
                   test=0.0
!---------------------------------------------------
                  CALL ISOFLUX_PT (TEST,NP)
!-----------------------------------------------------------------
!  if(mytid==0.and.n==2) print *, 'po4 isoflux term',tf1(40,34,10)
!  if(mytid==0.and.n==3) print *, 'ldoc isoflux term',tf1(40,34,10)
!   do k=1,km
!       do j=jsm,jem
!        do i=2,imm
!	    if(mytid==0.and.n==3)deltapo4(i,j,k)=tf1(i,j,k)-test(i,j,k)
!  if(mytid==0.and.n==3.and.deltapo4(i,j,k)<-0.0001.and.vit(i,j,k)>0.5) print*, &
!      'ldoc error in mixing term,delta',i,j,k,deltapo4(i,j,k),vit(i,j,k)
!        enddo
!       enddo
!   enddo
!iso2222-------------------
               do k=1,km
                  do j=2,jem
                   do i=2,imm
                 if(isnan(tf1(i,j,k))) print *,'tf1 before isoflux is error:i,j,k,np',i_global(i),j_global(j),k,np,pt(i,j,k,np),ptb(i,j,k,np)
                 if(isnan(test(i,j,k))) print *,'tf1 after isoflux is error:i,j,k,np',i_global(i),j_global(j),k,np
                  enddo
                 enddo
                enddo
!
              tf1=tf1+test
#else
 
#if ( defined SMAG)
!     Calculate AH3 
!-------------------------------------
              CALL SMAG3
!------------------------------------- 
!$OMP PARALLEL DO PRIVATE (K,J,WKI)
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0
               DO I = 2,IMT
                  WKI (I) = 0.5* (AH3 (I,J,K) + AH3 (I -1,J,K))* (PTB ( &
                           I,J,K,NP) - PTB (I -1, &
                  J,K,NP))* VIT (I,J,K)* VIT (I -1,J,K)
               END DO
               DO I = 2,IMM
                  TF1 (I,J,K) = TF1 (I,J,K) + SOTX (J)*(WKI(I+1)-WKI(I))
               END DO
            END DO
         END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
         DO J = JSM,JEM
         DO I = 2,IMM
            wt1= 0.5*(AH3(I,J,K)+AH3(I,J-1,K))*(PTB(I,J,K,NP)-PTB(I,J-1,K,NP))*VIT(I,J,K)*VIT(I,J-1,K)
            wt2= 0.5*(AH3(I,J+1,K)+AH3(I,J,K))*(PTB(I,J+1,K,NP)-PTB(I,J,K,NP))*VIT(I,J+1,K)*VIT(I,J,K)
            TF1 (I,J,K) = TF1 (I,J,K) + (R2D (J)* wt2 - R2C(J)* wt1)
         END DO
         END DO
         END DO
 
#else
 
!-----------------------------------------------------------------------
!     COMPUTE THE EDDY-DIFFUSION TERM :  ZONAL COMPONENT
!-----------------------------------------------------------------------
 
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,JMT
                  WKD1 (I,J,K) = 0.0
               END DO
            END DO
         END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,WKI)
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0
               DO I = 2,IMT
                  WKI (I) = (PTB (I,J,K,NP) - PTB (I -1,J,K,NP))* VIT (I, &
                           J,K)* VIT (I -1,J,K)
               END DO
 
               DO I = 2,IMM
                  WKD1 (I,J,K) = AH3 (I,J,K)* SOTX (J)* (WKI (I +1) - WKI (I))
               END DO
               
               IF(NX_PROC==1) THEN 
               WKD1 (1,J,K)= WKD1 (IMM,J,K)
               WKD1 (IMT,J,K)= WKD1 (2,J,K)
               ENDIF
            END DO
         END DO
#ifdef SPMD
       call exch_boundary(wkd1(1,1,1),km)
#endif
 
!$OMP PARALLEL DO PRIVATE (K,J,WKI)
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0
               DO I = 2,IMT
                  WKI (I) = (WKD1 (I,J,K) - WKD1 (I -1,J,K))* VIT (I,J,K) &
                           * VIT (I -1,J,K)
               END DO
               DO I = 2,IMM
                  TF1 (I,J,K) = TF1 (I,J,K) + SOTX (J)* (WKI (I +1) - WKI (I))
               END DO
            END DO
         END DO
 
#else
! default diffusion (ZONAL COMPONENT)
!$OMP PARALLEL DO PRIVATE (K,J,WKI)
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0
               DO I = 2,IMT
                  WKI (I) = (PTB (I,J,K,NP) - PTB (I -1,J,K,NP))* VIT (I, &
                           J,K)* VIT (I -1,J,K)
               END DO
               DO I = 2,IMM
                  TF1 (I,J,K) = TF1 (I,J,K) + AH3 (I,J,K)* SOTX (J)* (    &
                              WKI (I +1) - WKI (I))
               END DO
            END DO
         END DO
 
#endif
 
!-----------------------------------------------------------------------
!     COMPUTE THE EDDY-DIFFUSION TERM :  MERIDIONAL COMPONENT
!-----------------------------------------------------------------------
 
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,JMT
                  WKD1 (I,J,K) = 0.0
               END DO
            END DO
         END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
                  wt1= (PTB (I,J,K,NP)- PTB(I,J-1,K,NP))*VIT(I,J,K)* VIT (I,J -1,K)
                  wt2= (PTB (I,J+1,K,NP)- PTB(I,J,K,NP))*VIT(I,J,K)* VIT (I,J +1,K)
                  wkd1 (I,J,K) = ah3(i,j,k)*(R2D (J)* wt2 - R2C(J)* wt1)
               END DO
            END DO
         END DO

        IF(NX_PROC==1)THEN 
            DO J = JSM,JEM
!---xu  likely typing mistakes----------
!               WKD1 (1,J,K)= WKD1 (JMM,J,K)
!               WKD1 (JMT,J,K)= WKD1 (2,J,K)
!------------------------------------------
               WKD1 (1,J,K)= WKD1 (iMM,J,K)
               WKD1 (iMT,J,K)= WKD1 (2,J,K)
            END DO
!         END DO
         ENDIF
#ifdef SPMD
!------------------------------------
         call exch_boundary(wkd1(1,1,1),km)
!------------------------------------
#endif 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
         DO J = JSM,JEM
         DO I = 2,IMM
            wt1= (WKD1 (I,J,K) - WKD1 (I,J -1,K))* VIT (I,J,K)* VIT (I,J -1,K)
            wt2= (WKD1 (I,J+1,K) - WKD1 (I,J ,K))* VIT (I,J,K)* VIT (I,J +1,K)
            TF1 (I,J,K) = TF1 (I,J,K) + (R2D (J)* wt2 - R2C(J)* wt1)
         END DO
         END DO
         END DO

#else
!default difussion (MERIDIONAL COMPONENT)
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
         DO J = JSM,JEM
         DO I = 2,IMM
            wt1= (PTB (I,J,K,NP) - PTB (I,J -1,K,NP))* VIT (I,J,K)* VIT (I,J -1,K)
            wt2= (PTB (I,J +1,K,NP) - PTB (I,J,K,NP))* VIT (I,J,K)* VIT (I,J +1,K)
            TF1 (I,J,K) = TF1 (I,J,K) + ah3(i,j,k)* (R2D (J)* wt2 - R2C(J)* wt1)
         END DO
         END DO
         END DO

#endif
#endif
!
!iso22222----------------------
#endif
!
!-------------------------------------------------
!          VERTICAL COMPONENT
!--------------------------------------------------
 
!     EDDY-DIFFUSION
 
        wt1=0
!$OMP PARALLEL DO PRIVATE (K,J,I)
        DO K=2,KM-1
           DO J=JSM,JEM
           DO I=2,IMM
              wt1= WKC1(I,J,K-1)*(PTB(I,J,K-1,NP)-PTB(I,J,K,NP))*ODZT(K)*VIT(I,J,K)
              wt2= WKC1(I,J,K)*(PTB(I,J,K,NP)-PTB(I,J,K+1,NP))*ODZT(K+1)*VIT(I,J,K+1)
              TF1 (I,J,K)= TF1 (I,J,K)+ODZP(K)*(wt1-wt2)*(1.0-AIDIF)
                 if(isnan(tf1(i,j,k))) print *,'tf1 vertical diffusion is error:i,j,k,np',i_global(i),j_global(j),k,np
           END DO
           END DO
        END DO
!$OMP PARALLEL DO PRIVATE (J,I)
       DO J = JSM,JEM
       DO I = 2,IMM
           wt1= WKC1(I,J,1)*(PTB(I,J,1,NP)-PTB(I,J,2,NP))*ODZT(2)*VIT(I,J,2)
           wt2= WKC1(I,J,km-1)*(PTB(I,J,km-1,NP)-PTB(I,J,km,NP))*ODZT(km)*VIT(I,J,km)
           TF1(I,J,1)=TF1 (I,J,1)-ODZP(1)*wt1*(1.0-AIDIF)
           TF1(I,J,km)=TF1(I,J,km)+ODZP(km)*wt2*(1.0-AIDIF)
       END DO
       END DO
!------------------------------------------
!     SET SURFACE BOUNDARY CONDITION
!------------------------------------------
!lyc
         stf1(:,:)=0 
         ssfc(:,:)=0.
!--------------------------------------
         IF (NP == 1) THEN
!--------------------------------------
! only give surface flux of DIC ,radiocarbon, CFC
#ifdef printcall
#ifdef SPMD
            print*,"call flux_pt in ptracer, mytid=",mytid
#else
            print*,"call flux_pt in ptracer"
#endif
#endif
!------------------------------           
            CALL FLUX_PT
!--------------------------------------
#if (defined  carbonC14)||(defined cfc)
!--------------------------
            call flux_ot
!------------------------
#endif
!-------------------------------
            DO J = JSM,JEM
               DO I = 2,IMM
                IF (ITNU (I,J) > 0) THEN
#if (defined  carbonC14)||(defined cfc)
!
                 STF1 (I,J) = ssfc (I,J)*dzp(1)
#endif           
!xu-- 
!#if (defined carbonC) || (defined carbonBio)
#if (defined carbonC) || (defined carbonBio) ||(defined carbonAbio)
!
                 STF1 (I,J) = ssfc (I,J) 
#endif            
                 TF1 (I,J,1) = TF1 (I,J,1) + STF1 (I,J)*ODZP(1)* (1.0- AIDIF)
                END IF
               END DO
            END DO
        END IF
#ifdef carbonBio
    IF (NP==4) THEN
        do j=1,jmt
	   do i=1,imt
	     ssfc(i,j)=0.0
	   enddo
        enddo
!-------------------------
	  CALL FLUX_TA
!-------------------------
          do j=2,jem
            do i=2,imm
              IF (ITNU (I,J) > 0) THEN
                 stf1(i,j)= ssfc(i,j)
                 TF1 (I,J,1) = TF1 (I,J,1) + STF1 (I,J)*ODZP(1)* (1.0- AIDIF)
              ENDIF
            enddo
          enddo
      ENDIF
!cm090330--------------------------
! surface flux of o2        
         IF (NP == 5) THEN
	    ssfc(:,:)=0.
!-------------------------
            CALL FLUX_o2
!-------------------------
            DO J = JSM,JEM
               DO I = 2,IMM
                IF (ITNU (I,J) > 0) THEN
                    STF1 (I,J) = ssfc (I,J) 
                    TF1 (I,J,1) = TF1 (I,J,1) + STF1 (I,J)*ODZP(1)* (1.0- AIDIF)
                END IF
               END DO
            END DO
         END IF
!cm090330----------------------
!-------------------------------------------------------------
! for fe (lyc,2013,07)
!-------------------------------------------------------------
      IF(NP==6) THEN
       DO J=JSM, JEM
         DO I=2,IMM
          IF(ITNU(I,J)>0) THEN
           STF1(I,J)=FE_F(I,J)
           TF1(I,J,1)=TF1(I,J,1)+STF1(I,J)*ODZP(1)*(1.0-AIDIF)
          ENDIF
         ENDDO
        ENDDO
     ENDIF
       
!
#endif
!lyc end
!*************************************************************
!     Calculate biological source of DIC, PO4, LDOC and TA (biosource)
!-----------------------------------------------------------------------         
#ifdef carbonBio
#ifdef printcall
#ifdef SPMD
            print*,"call biosource in ptracer, idx=',N,'mytid=",mytid
#else
            print*,"call biosource in ptracer, idx=",N
#endif
#endif           
!---------------------------------------------
               CALL BIOSOURCE(TF1,NP)
!--------------------------------------------
#endif
!
#ifdef carbonBio
#ifdef carbonDebug
       do k=1,km
        do j=jsm,jem
         do i=2,imm
	     if(vit(i,j,k)<0.5) cycle
             if(np==3.and.tf1(i,j,k) < -0.0001) then
               print*, 'lodc error! i,j,k,lodc:',i,j,k,j_global(j),tf1(i,j,k),'mytid=',mytid   
	       print*, 'biosourc lodc!,i,j,k:', i,j,k,j_global(j),deltapo4(i,j,k),'mytid=',mytid
               print*,'lodcb!',i,j,k,j_global(j),ptb(i,j,k,3)
	       print*,'salinity',at(i,j,k,2)
	       print*,'tempreture',at(i,j,k,1)
	       print*,'u,v',utf(i,j,k),vtf(i,j,k)
	       print*,'ta',pt(i,j,k,4)
	       print*,'tc',pt(i,j,k,1)
	       print *,'error biosource!'
             endif
	 enddo
	enddo
       enddo
#endif  
#endif       
!--------------------------------------------------------------
!     SOLVE FOR "TAU+1" PASSIVE TRACER AT CENTER OF "T" CELLS
!--------------------------------------------------------------
!lyc test           if(mytid==1.and.n==6) then
!            print*,'the source of fe', tf1(91,jmt/2,:)
!           endif
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
                  VTL1 (I,J,K) = PTB (I,J,K,NP) + C2DTTS * TF1 (I,J,K)
!
!cm080418----------for decay of carbonC14------------------------------
#ifdef preindustrial 
#ifdef carbonC14
                  VTL1 (I,J,K) =VTL1 (I,J,K) - C2DTTS * rdca * PTB (I,J,K,NP)
#endif
#endif
!cm080418----------for decay of carbonC14------------------------------

               END DO
            END DO
         END DO
!
!    if(mytid==0.and.n==3) print*, 'lodc debug!',ptb(40,34,10,3),vtl1(40,34,10),tf1(40,34,10) 
!lyc end
!
#ifdef carbonBio
#ifdef carbonDebug
      do k=1,km
        do j=jsm,jem
          do i=2,imm
            if(vit(i,j,k) < 0.5) cycle
            if(vtl1(i,j,k) < 1200.0) then
              if(np==1) print*, 'TC error! i,j,k,TC:',i_global(i),j_global(j),k,vtl1(i,j,k),'mytid=',mytid   
              if(np==4) print*, 'TA error! i,j,k,TA:',i_global(i),j_global(j),k,vtl1(i,j,k),'mytid=',mytid   
            endif
	    if(vtl1(i,j,k) < -0.01) then
              if(np==2) print*, 'PO4 error! i,j,k,PO4:',i_global(i),j_global(j),k,vtl1(i,j,k),'mytid=',mytid   
              if(np==3) print*, 'LDOC error! i,j,k,LDOC:',i_global(i),j_global(j),k,vtl1(i,j,k),tf1(i,j,k),ptb(i,j,k,3),'mytid=',mytid   
            endif
          enddo
        enddo
      enddo
#endif 
#endif     
!-----------------------------------------------------------------------
!     ADD DT/DT COMPONENT DUE TO IMPLICIT VERTICAL DIFFUSION
!-----------------------------------------------------------------------
 
#if (defined ISO)
!-------------------------------------------------
         CALL INVTRI (VTL1,STF1,WKC1,AIDIF,C2DTTS)
!-------------------------------------------------
#endif

!----------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
!------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J)
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM  !lyc2014.06
                  VTL1 (I,J,K) = VTL1 (I,J,K)*VIT(I,J,K)
               END DO
            END DO
         ENDDO

!lyc 2014.06
        IF(NX_PROC==1) then
           DO K = 1,KM
            DO J = JSM,JEM
!lyc set the cyclic conditions
             VTL1(1,J,K)=VTL1(IMM,J,K)
	     VTL1(IMT,J,K)=VTL1(2,J,K)
            ENDDO
           ENDDO
        ENDIF
!   
         do k=1,km
              do i=1,imt
                 if(isnan(vtl1(i,1,k))) then
                  print *,'lyc test for the northest boundary,before exch:i,j,pt,np,tf1,vit',i_global(i),vtl1(i,1,k),j_global(1),pt(i,1,k,np),np,tf1(i,1,k),vit(i,1,k)
                endif
              enddo
          enddo
#ifdef SPMD
        call exch_boundary(vtl1(1,1,1),km)
#endif
         do k=1,km
              do i=1,imt
                 if(isnan(vtl1(i,1,k))) then
                  print *,'lyc test for the northest boundary,',vtl1(i,1,k),j_global(1),pt(i,1,k,np),np
                  stop
                endif
              enddo
          enddo
!   2003,8,2 changed by lrf for 1 X 1 model
!
!        CALL SMTS (VTL,VIT,KM,50,70,294,331)
!
!-------------------------------------
         call SMTS(VTL1,VIT,KM,fil_latp)
!---------------------------------
         do k=1,km
              do i=1,imt
                 if(isnan(vtl1(i,1,k))) then
                  print *,'lyc test for the northest boundary,after smts:i,j,k,vtl1,pt,np',i_global(i),j_global(1),k,vtl1(i,1,k),pt(i,1,k,np),np
                 stop
                endif
              enddo
          enddo

#ifdef carbonBio      
!
#ifdef carbonDebug
      do k=1,km
        do j=jsm,jem
          do i=1,imt
            if(vit(i,j,k) < 0.5) cycle
            if(np==1.or.np==4) then
            if(vtl1(i,j,k) < 1700.0) vtl1(i,j,k)=1700.0 
            if(vtl1(i,j,k) > 2500.0) vtl1(i,j,k)= 2500.0
            endif
	    if (np==2.or.np==3) then
            if(vtl1(i,j,k) < -1.0E-3)  vtl1(i,j,k)=0.001
            if(np==2.and.vtl1(i,j,k)>10.0) vtl1(i,j,k)=10.0
            if(np==3.and.vtl1(i,j,k)>110.0) vtl1(i,j,k)=110.0
            endif
          enddo
        enddo
      enddo
       do k=1,km
        do j=jsm,jem
          do i=1,imt
            if(np==2.and.vtl1(i,j,1)>5.0) then
             print *, i_global(i),j_global(j),k,vtl1(i,j,k),po4force(i,j,k),wst(i,j,k),wst(i,j,k+1),u(i,j,k),v(i,j,k)
	     stop
	    endif
	     enddo
        enddo
      enddo
#endif 
#endif
!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
#ifdef SPMD
!----------------------------------
         CALL exch_boundary(vtl1(1,1,1),km)
!----------------------------------
#endif
!
#ifdef buchang 
#if (defined carbonBio)||(defined carbonC)||(defined carbonC14)
!-------------------------------------------------------------
!lyc-------calculate the tc(n-1)+ssfc*c2dtts
!-----------------------------------------------
    if(np==1) then
         tcsumb=0.0
	 tctmpb=0.0
	    do j=2,jem
	      do i=2,imm
		if(vit(i,j,1)<0.5) cycle
		tcsumb=tcsumb+stf1(i,j)*dxdyt(j)*dts
	      enddo
	    enddo
#ifdef SPMD	    
	  call mpi_barrier(mpi_comm_ocn,ierr)
          call mpi_reduce(tcsumb,tctmpb,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    
	  if(mytid==0) print *,' the sum of stf1',tctmpb
#else
	  print *,' the sum of stf1',tcsumb
#endif
	  tctmpb=0.0
	   do k=1,km
	    do j=2,jem
	      do i=2,imm
		if(vit(i,j,k)<0.5) cycle
		 tcsumb=tcsumb+pt(i,j,k,1)*dxdyt(j)*dzp(k)
	      enddo
	    enddo
	   enddo
#ifdef SPMD	   
	  call mpi_barrier(mpi_comm_ocn,ierr)
          call mpi_reduce(tcsumb,tctmpb,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    
#else
	  tctmpb=tcsumb
#endif
	tcsum=0.0
	tctmp=0.0
!	vseasum=0.0
!	vsea_bio=0.0
	   do k=1,km
	    do j=2,jem
	      do i=2,imm
		if(vit(i,j,k)<0.5) cycle
		 tcsum=tcsum+vtl1(i,j,k)*dxdyt(j)*dzp(k)
!		 vseasum=vseasum+dxdyt(j)*dzp(k)
	      enddo
	    enddo
	   enddo
#ifdef SPMD	   
	  call mpi_barrier(mpi_comm_ocn,ierr)
          call mpi_reduce(tcsum,tctmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    
#else
	  tctmp=tcsum
#endif	 
!
#ifdef carbonBio	
#ifdef SPMD	  
	if(mytid==0) then
#endif	    
	    deltatc=(tctmpb-tctmp)/vsea
!          print *, 'the compensatory concentration of pt(1)',deltatc
#ifdef SPMD
        endif
#endif	
#else
       vseac=0.0
       do k=1,km
	   do j=2,jmm
	       do i=2,imm
	       if(vit(i,j,k)<0.5) cycle
		if(vtl1(i,j,k)<0.0001)cycle
		 vseac=vseac+dxdyt(j)*dzp(k)
		enddo
           enddo
       enddo
#ifdef SPMD       
        call mpi_barrier(mpi_comm_ocn,ierr)
        call mpi_reduce(vseac,vseap,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    
	if(mytid==0) then
#endif        
	    deltatc=(tctmpb-tctmp)/vseap
	     print *, 'the compensatory concentration of pt(1)',deltatc
#ifdef SPMD       
       endif  
#endif
#endif


#ifdef SPMD       
          call mpi_bcast(deltatc,1,mpi_real8,0,mpi_comm_ocn,ierr)
#endif
          do k=1,km
            do j=jst,jet
              do i=1,imt
	       if(vit(i,j,k)<0.5) cycle
#if (defined carbonC)||(defined carbonC14)
		if(vtl1(i,j,k)<0.0001) cycle
#endif		    
               vtl1(i,j,k)=vtl1(i,j,k)+deltatc
              enddo
            enddo
          enddo
      endif
#endif      
#endif	 
!cm090302----
!         IF (ISP >= 1)THEN
!!lyc 2011.02.15$OMP PARALLEL DO PRIVATE (K,J,I)
!            DO K = 1,KM
!               DO J = JST,JET
!                  DO I = 1,IMT
!         IF (mod(ISP,15)/=0) THEN
!                     PTB (I,J,K,N) = AFT2* PT (I,J,K,N) + AFT1* (PTB (I, &
!                                    J,K,N) + VTL1 (I,J, K))
!         ELSE
!	             PTB (I,J,K,N) =VTL1 (I,J,K)
!         END IF
!!cm090302----
!		  END DO
!               END DO
!            END DO
!!         END IF
!------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JST,JET
               DO I = 1,IMT
                  PTF(I,J,K,NP) = PT (I,J,K,NP)
                  PT (I,J,K,NP) = VTL1 (I,J,K)
               END DO
            END DO
         END DO
!!555---------------------
!
#ifdef carbonBio	 
#ifdef buchang	 
  tcsum=0.0
  tctmp=0.0
  tasum=0.0
  tatmp=0.0
  po4sum=0.0
  po4tmp=0.0
  ldocsum=0.0
  ldoctmp =0.0
     DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
		  if(vit(i,j,k)<0.5) cycle
		  if(np==1) then
		      tcsum=tcsum+pt(i,j,k,np)*dxdyt(j)*dzp(k)
		  endif
		  if(np==4) then
		      tasum=tasum+pt(i,j,k,np)*dxdyt(j)*dzp(k)
		  endif
		  if(np==2) po4sum=po4sum+pt(i,j,k,np)*dxdyt(j)*dzp(k)
               END DO
            END DO
         END DO
!
       if(np==3) then
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
		  if(vit(i,j,k)<0.5) cycle
                  if(pt(i,j,k,3)>0.0) then
		  ldocsum=ldocsum+pt(i,j,k,np)*dxdyt(j)*dzp(k)
                  endif
               END DO
            END DO
         END DO
       endif

     if(np==3) then
	 call mpi_barrier(mpi_comm_ocn,ierr)
	 call mpi_reduce(ldocsum,ldoctmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    

      	if(mytid==0) then
         delldoc=(vsea*ldocga-ldoctmp)/vsea
!        print *, 'the total storage of ldoc!', ldoctmp
!        print *, 'the compensatory concentration of ldoc',delldoc
        endif  
          call mpi_bcast(delldoc,1,mpi_real8,0,mpi_comm_ocn,ierr)
          do k=1,km
            do j=jst,jet
              do i=1,imt
	       if(vit(i,j,k)<0.5) cycle
               pt(i,j,k,np)=pt(i,j,k,np)+delldoc
              enddo
            enddo
          enddo
     endif

     if(np==1) then
	  call mpi_barrier(mpi_comm_ocn,ierr)
          call mpi_reduce(tcsum,tctmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    
	if(mytid==0) print *, 'the total storage of tc!', tctmp
      endif

     if(np==4) then
	  call mpi_barrier(mpi_comm_ocn,ierr)
          call mpi_reduce(tasum,tatmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    
	if(mytid==0) then
!         print *, 'the total storage of ta!', tatmp
         deltata=(vsea*2377-tatmp)/vsea
        endif  
          call mpi_bcast(deltata,1,mpi_real8,0,mpi_comm_ocn,ierr)
          do k=1,km
            do j=jst,jet
              do i=1,imt
	      if(vit(i,j,k)<0.5) cycle
               pt(i,j,k,np)=pt(i,j,k,np)+deltata
              enddo
            enddo
          enddo
     endif
     if(np==2) then
	  call mpi_barrier(mpi_comm_ocn,ierr)
          call mpi_reduce(po4sum,po4tmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    
      if(mytid==0) then
!	print *, 'the total storage of po4!', po4tmp
         deltapo4=(vsea*2.38-po4tmp)/vsea
      endif  
       call mpi_bcast(deltapo4,1,mpi_real8,0,mpi_comm_ocn,ierr)
       do k=1,km
            do j=jst,jet
              do i=1,imt
		if(vit(i,j,k)<0.5) cycle
               pt(i,j,k,np)=pt(i,j,k,np)+deltapo4
              enddo
            enddo
          enddo

      endif

      do k=1,km
	  do j=1,jmt
	      do i=1,imt
		if(vit(i,j,k)<0.5) cycle
		if(pt(i,j,k,np)<0.0) pt(i,j,k,np)=0.0
	      enddo
	  enddo
       enddo
!555---------------------------------
#endif      
#endif      
      END DO
!-------------------------------------
             ISP = ISP +1
!             print *,'ISP=',ISP !juanxiong he
!-------------------------------------
! 
#ifdef ISO
!zhaoliang      deallocate(K1,K2,K3,adv_vetiso,adv_vbtiso,adv_vntiso)
!lyc test
!      deallocate(K1,K2,K3,adv_vetiso,adv_vbtiso,adv_vntiso)
#endif
      deallocate(stf1,tf1)
      deallocate(wkb1,wkc1,wkd1,uwk1,wst)
      deallocate(h0f1)
      deallocate(utf1,vtf1,utl1,vtl1)
!--------------------------------------------
      RETURN
      END SUBROUTINE PTRACER
 
! CVS: $Id: biosource.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  SUBROUTINE BIOSOURCE(sourceb,idx)
!========================
! BIOSOURCE
!---------------------------------------------------------------------
!
! purpose: calculate biological source of DIC, PO4, LDOC and TA
!
! author: Zhao Liang@lapc 2004/05/21
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod
      USE tracer_mod
      USE carbon_mod
      USE cforce_mod
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!#include <netcdf.inc>      
!#ifdef SPMD      
!#include <mpif.h>
!#endif      
!      
!---------------------------------------------------------------------
!     This sub is used to calculate the FLUX of a natural carbon
!     cycle including a simple biological paramterizationi
!     1) pCO2 flux at surface
!     2) TC, TA 
!     3) PO4 
!---------------------------------------------------------------------
!
      integer::idx
      real,dimension(imt,jmt,km)::sourceb,atest
!cm---------------------
!o20:oxygen consumption is assumed to be halted below this critical oxygen level
      real::o20
!cm---------------------

#ifdef carbonBio

      real::r_cpr
!
      r_cpr=1.0/r_cp

!cm---------------------
      o20=4.0
!cm---------------------

!     calculate biological source of DIC       
       IF(idx==1) THEN
         do k=1,km
           do j=jsm,jem
             do i=2,imm
               sourceb(i,j,k)=sourceb(i,j,k) &
                            -r_cp*a0_b(i,j,k)-b0_b(i,j,k)-c_b(i,j,k)
             enddo
           enddo
         enddo
       ENDIF    
!     calculate biological source of PO4       
       IF(idx==2) THEN
         do k=1,km
           do j=jsm,jem
             do i=2,imm
               sourceb(i,j,k)=sourceb(i,j,k) &
                            -a0_b(i,j,k)-b0_b(i,j,k)*r_cpr
             enddo
           enddo
         enddo
       ENDIF    
!     calculate biological source of LDOC       
       IF(idx==3) THEN
         do k=1,km
           do j=jsm,jem
            do i=2,imm
               sourceb(i,j,k)=sourceb(i,j,k)+b0_b(i,j,k)
             enddo
           enddo
         enddo
!	if(mytid==0) print *, 'ldoc source ',b0_b(40,34,10)
       ENDIF    
!     calculate biological source of TA       
       IF(idx==4) THEN
         do k=1,km
           do j=jsm,jem
             do i=2,imm
               sourceb(i,j,k)=sourceb(i,j,k) &
                            +r_np*(a0_b(i,j,k)+b0_b(i,j,k)*r_cpr)-2.0*c_b(i,j,k)
             enddo
           enddo
         enddo
       ENDIF    
!cm---------------------
!     calculate biological source of o2   
       IF(idx==5) THEN
         do k=1,km
           do j=jsm,jem
             do i=2,imm
		 if(pt(i,j,k,5)>o20) then
                     sourceb(i,j,k)=sourceb(i,j,k) &
                            +r_o2p*(-a0_b(i,j,k)-b0_b(i,j,k)*r_cpr)
	         endif
             enddo
           enddo
         enddo
       ENDIF    
!cm---------------------
!
!lyc 2013.05------------------------------------------------------------
! calculate the biosource of Fe
      IF(idx==6) THEN
        do k=1,km
          do j=jsm,jem
            do i=2,imm
             sourceb(i,j,k)=sourceb(i,j,k)+Fe_source(i,j,k)
           enddo
          enddo
        enddo
      ENDIF       
#ifdef carbonDebug      
#ifdef SPMD
      print*, 'subroutine biosource-----------------------------------------------------------','mytid=',mytid
      print*,'idx,i,j,1,source:',idx,imt/2,jmt/2,sourceb(imt/2,jmt/2,1),'mytid=',mytid
      print*,'idx=',idx,'a0=',a0_b(imt/2,jmt/2,1),'b0=',b0_b(imt/2,jmt/2,1), &
             'c=',c_b(imt/2,jmt/2,1), &
             'a1=',a1_b(imt/2,jmt/2,1),'a2=',a2_b(imt/2,jmt/2,1), &
             'b1=',b1_b(imt/2/2,jmt/2,1),'b2=',b2_b(imt/2,jmt/2,1), &
             'mytid=',mytid,"BIOSOURCE"
#else        
      print*, 'subroutine biosource----------------------------------------------------------'
      print*,'idx,i,j,1,source:',idx,imt/2,jmt/2,sourceb(imt/2,jmt/2,1),'mytid=',mytid
      print*,'idx=',idx,'a0=',a0_b(imt/2,jmt/2,1),'b0=',b0_b(imt/2,jmt/2,1), &
             'c=',c_b(imt/2,jmt/2,1), &
             'a1=',a1_b(imt/2,jmt/2,1),'a2=',a2_b(imt/2,jmt/2,1), &
             'b1=',b1_b(imt2/2,jmt/2,1),'b2=',b2_b(imt/2,jmt/2,1), &
             "BIOSOURCE"
#endif
#endif
	     
#endif      
!      
      RETURN
      END SUBROUTINE BIOSOURCE


! CVS: $Id: flux_pt.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
!------------------------
     SUBROUTINE FLUX_PT
!========================
! FLUX_PT
!---------------------------------------------------------------------
!
! purpose: give surface flux of CCs (carbon or C14)
!          (original author: Xu Y F)
!
! author: Zhao Liang@lapc 2004/03/03
!
!---------------------------------------------------------------------
#include <def-undef.h>      
!---------------------- 
!
      USE param_mod
      USE pconst_mod
      USE tracer_mod
      USE carbon_mod
      USE coutput_mod
      USE cforce_mod
      USE forc_mod   !for sea ice
#ifdef COUP
      USE buf_mod
#endif
!----------------------------------
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!----------------------------------------------
      IMPLICIT NONE
#include <netcdf.inc>      
!#ifdef SPMD      
!#include <mpif.h>
!#endif      
!      
!for forcing data
!lyc 2015.08

#if (defined carbonC) || (defined carbonBio) ||(defined carbonAbio)
!---------------------------------------------------------------------
!     This sub is used to calculate the FLUX of a natural carbon
!     cycle including a simple biological parameterization
!     1) pCO2 flux at surface
!     2) STC, STA 
!---------------------------------------------------------------------
!
      REAL T0,ES,FACTOR,FACTX
#ifdef carbonC
      REAL C1(IMT,JMT),T1(IMT,JMT),S1(IMT,JMT)
      REAL A1,A2,A3,B1,B2,B3,Z0,Z1
#endif
#endif      
!----xu-----
#if(defined carbonBio)||(defined carbonAbio)
      REAL tt(imt,jmt),ss(imt,jmt),tatmp(imt,jmt),tctmp(imt,jmt)
#endif
!       
#ifdef carbonC14
!---------------------------------------------------------------------
!     GIVE SURFACE FLUX of C14
!---------------------------------------------------------------------
!cm 081027,use wind velosity to compute overd
#ifdef nc14wind
       real,dimension(IMT,JMT)::overdwind
       real ccs1
#else
       REAL ccs1,overd
#endif

!cm 081027
!---------------------------------------------------------------------
! one year=1.0*3.1536E7 seconds       
! a timescale of 5 year is used to relax the modelled values to
! observed ones with the surface thickness is 50m.
! if the thickness of the surface layer is 10 meter, the relax time 
! should be 1 year.
!---------------------------------------------------------------------
#ifdef nc14wind
!from Toggweiler,1989
do j=1,jmt
    do i=1,imt
	overdwind(i,j)=(3.5*(w22np(i,j)-2))/(2.0*dzp(1)*3.1536E7)
!	if(overdwind(i,j)<0.) overdwind(i,j)=0.
    enddo
enddo
!print *,((w22np(90,j),overdwind(90,j)*3.1536E7),j=1,jmt)
!stop
#else
      overd=1.0/(1.0*3.1536E7)
#endif
#endif
!----------------------------------------------------------
!      
!lyc
#ifdef cfc
INTEGER:: NN
real,dimension(2)::a1,a2,a3,a4,b1,b2,b3
real,dimension(imt,jmt)::ccs
integer::l 
#endif

#ifdef COUP
     w22np=sqrt(duu10n)
     pressureday=patm
     pressureday=patm/1.01325*1.0E-5
!     if(mytid==0) print *,'pressureday: ',pressureday(2,2),'wind: ',w22np(2,2)
          !pressureday(i,j)=psa3(i,j,1)/1.01325*1.0E-5
#endif     
!lyc 2014.09.11
#ifdef COUP
          iceday=ifrac
#else
          iceday=seaice
#endif

#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
!---------------------------------------------------------------------
!     xu, calculation of change in PCO2A with time
!     pCO2  in the atmosphere
!---------------------------------------------------------------------
!     pco2dry from INTFOR_PT.F90 from yearly mean values
!---------------------------------------------------------------------
#ifdef COUP
 !      pco2dry=pco2(1,1)
 !      pco2s=284.725*1.0E-6
       pco2dry=284.725  !for instant purpose
      !if(mytid==0) write(6,*) 'the atmospheric carbon dioxide concentration(ppm) :', pco2dry
#endif
!$OMP PARALLEL DO PRIVATE (I,J,T0,ES)
      DO J=1,JMT
        DO I=1,IMT
          IF(VIT(i,j,1) < 0.5) CYCLE
          T0 = 273.15 + AT(I,J,1,1)
          ES = EXP(20.1050 - 9.7982E-3*T0 - 6163.10/T0)
          pco2a(I,J) = pco2dry * (1.0 - ES)
#ifdef carbonC
          dpco2a(I,J)=pco2a(I,J)-pco2dry0*(1.0-ES)
#endif          
        ENDDO
      ENDDO
!---------------------------------------------------------------------
!     PCO2O: PARTIAL PRESSURE OF CO2-OCEAN (uatm)
!---------------------------------------------------------------------
!     C1(I,J) = PTB(I,J,1,1)/RHOSRF(I,J)*1.0E9
!     C1(I,J) =(PTB(I,J,1,1)+SSTC(i,j))/RHOSRF(I,J)
!     TALK(i,j)=50.26*s(i,j,1)+2317.
!---------------------------------------------------------------------
#ifdef carbonC          
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J=1,JMT
        DO I=1,IMT
          C1(I,J)=0.0
          T1(I,J)=0.0
          IF(VIT(I,J,1)==0) CYCLE
          T1(I,J)=AT(I,J,1,1)
          C1(I,J)=PTB(I,J,1,1)
        ENDDO
      ENDDO
!---------------------------------------------------------------------
!     Calculation of dpco2 from C1 using Sarmiento et al.'s equation
!     If C1 is larger than 100 umol/kg, eq is invalid (too large error).
!---------------------------------------------------------------------
      A1 =  1.7561
      A2 = -0.031618
      A3 =  0.0004444
      B1 =  0.004096
      B2 = -7.7086E-5
      B3 =  6.1000E-7
      
!$OMP PARALLEL DO PRIVATE(I,J,Z0,Z1)
      DO J=1,JMT
        DO I=1,IMT
! control       
          IF(C1(I,J) > 600.0) THEN
#ifdef SPMD
            PRINT*, "PTB is too big",I,J,C1(I,J),"MYTID=",mytid,"FLUX_PT"
#else            
            PRINT*, "PTB is too big",I,J,C1(I,J),"FLUX_PT"
#endif
          ENDIF
          Z0  = A1 + A2*T1(I,J) + A3*T1(I,J)*T1(I,J)
          Z1  = B1 + B2*T1(I,J) + B3*T1(I,J)*T1(I,J)
          pco2o(i,j) = Z0*C1(I,J) / ( 1.0-Z1*C1(I,J) )
        ENDDO
      ENDDO
#endif          
!----------------------------------------
#if(defined carbonBio) || (defined carbonAbio)          
!lyc 2014.10.28         call FLUX_ta
!------------------------
      DO J=1,JMT
        DO I=1,IMT
          tt(i,j)=at(i,j,1,1)
          ss(i,j)=at(i,j,1,2)
#ifdef carbonAbio
          tctmp(i,j)=ptb(i,j,1,1)+2000.
#else 
          tctmp(i,j)=ptb(i,j,1,1)
#endif
#ifdef carbonAbio
!-------------------------------------------
! xu, Calculation of inorganic carbon cycle
!------------------------------------------
          tatmp(i,j)=ssfc(i,j)
#else
!--------------------------
          tatmp(i,j)=ptb(i,j,1,4)
#endif
        ENDDO
      ENDDO
#ifdef carbonDebug
#ifdef SPMD
      print*, 'subroutine flux_pt-----------------------------','mytid=',mytid
      PRINT*,"i,j",imt/2,jmt/2, &
             "ta=",tatmp(imt/2,jmt/2),"tc=",tctmp(imt/2,jmt/2), &
             "tt=",tt(imt/2,jmt/2),"ss=",ss(imt/2,jmt/2), &
             "mytid=",mytid,"FLUX_PT"
!             
      do j=1,jmt
        do i=1,imt
          if (vit(i,j,1) < 0.5) cycle
          if(tt(i,j)< -1.8) print*,"abnormal tt===",i,j,tt(i,j),"mytid=",mytid,"FLUX_PT"  
          if(tatmp(i,j)< 100.0 .or. tatmp(i,j)>5000.0) &
            print*,"abnormal ta===",i,j,tatmp(i,j),"tt=",tt(i,j),"ss=",ss(i,j), &
                   "mytid=",mytid,"VIT=",vit(i,j,1),"FLUX_PT"  
          if(tctmp(i,j)< 100.0 .or. tctmp(i,j)>5000.0) &
            print*,"abnormal tc===",i,j,tctmp(i,j),"tt=",tt(i,j),"ss=",ss(i,j), &
                   "mytid=",mytid,"VIT=",vit(i,j,1),"FLUX_PT" 
        enddo
      enddo
#else
      print*, 'subroutine flux_pt--------------------------------------------------------'
      PRINT*,"i,j",imt/2,jmt/2, &
             "ta=",tatmp(imt/2,jmt/2),"tc=",tctmp(imt/2,jmt/2), &
             "tt=",tt(imt/2,jmt/2),"ss=",ss(imt/2,jmt/2), &
            "FLUX_PT"
!             
      do j=1,jmt
        do i=1,imt
          if (vit(i,j,1) < 0.5) cycle
          if(tt(i,j)< -1.8) print*,"abnormal tt===",i,j,tt(i,j),"FLUX_PT"  
          if(tatmp(i,j)< 1000.0 .or. tatmp(i,j)>3000.0) &
            print*,"abnormal ta===",i,j,tatmp(i,j),"tt=",tt(i,j),"ss=",ss(i,j), &
                   "VIT=",vit(i,j,1),"FLUX_PT"  
          if(tctmp(i,j)< 1000.0 .or. tctmp(i,j)>3000.0) &
            print*,"abnormal tc===",i,j,tctmp(i,j),"tt=",tt(i,j),"ss=",ss(i,j), &
                   "VIT=",vit(i,j,1),"FLUX_PT" 
        enddo
      enddo
#endif
#endif
#ifdef printcall
#ifdef SPMD
     print*, "call pco2 in flux_pt++++++++++++++++++++++++++++++++++", "mytid=", mytid
#else     
     print*, "call pco2 in flux_pt++++++++++++++++++++++++++++++++++"
#endif
#endif
!-----------------------------------------
      CALL OPCO2(pco2o,tatmp,tctmp,tt,ss)
!------------------------------------------
#endif
!---------------------------------------------------------------------
!     call SGEC.F90 to calculate the exchange coefficients of co2
!---------------------------------------------------------------------
        CALL SGEC
!---------------------------------------------------------------------
!     FIND CO2-FLUX FROM ATMOSPHERE TO OCEAN
!     PCODRY: PCO2-DRY ATMOSPHERE (PPM)
!     PCO2A: PARTIAL PRESSURE OF CO2-WET ATMOSPHERE (PPM)
!     ssfc in units of m/s *umol/kg (in LICOM)
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (I,J)
      DO J=1,JMT
        DO I=1,IMT
#ifdef carbonC
          dpco2o(I,J)=(dpco2a(I,J)-pco2o(I,J))*VIT(I,J,1)
#endif          
!-----xu------------------
!#ifdef carbonBio
#if(defined carbonBio)||(defined carbonAbio)
          dpco2o(I,J)=(pco2a(I,J)-pco2o(I,J))*VIT(I,J,1)
#endif          
!-----treatment of sea ice,cm080813--------
!          ssfc(I,J) = sge(I,J)*dpco2o(I,J)
          ssfc(I,J) =(1.0-iceday(i,j))*sge(I,J)*dpco2o(I,J)
! treanment of sea ice           
!          IF(AT(I,J,1,1)<=-1.8) ssfc(I,J) = 0.0
!          IF(ITICE(I,J)==1) ssfc(I,J) = 0.0
!-----cm080813-----------------------------
        ENDDO
      ENDDO
#ifdef carbonDebug      
#ifdef SPMD
      PRINT*,'i,j,pco2a:',imt/2,jmt/2,pco2a(imt/2,jmt/2),'mytid=',mytid,"FLUX_PT"
      PRINT*,'i,j,ssfc:',imt/2,jmt/2,ssfc(imt/2,jmt/2),'mytid=',mytid,"FLUX_PT"
#else        
      PRINT*,'i,j,pco2a:',imt/2,jmt/2,pco2a(imt/2,jmt/2),"FLUX_PT"
      PRINT*,'i,j,ssfc:',imt/2,jmt/2,ssfc(imt/2,jmt/2),"FLUX_PT"
#endif
#endif      
!---------------------------------------------------------------------
!     calculation of accumulative input CO2 
!       in one timestep in units of GtC
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(I,J,FACTX)
      DO J=1,JMT
        FACTX=DXDYT(J)*unitc*DTS
        DO I=1,IMT
          totup(I,J)=totup(I,J)+ssfc(i,j)*FACTX*VIT(I,J,1)
          tpco2o(I,J)=tpco2o(I,J)+pco2o(i,j)
          tdpco2o(I,J)=tdpco2o(I,J)+dpco2o(i,j)
        ENDDO
      ENDDO
#endif      
!--------------------------------------
!     GIVE SURFACE FLUX of C14
!-------------------------------------
#ifdef carbonC14
!
      DO J = 1,JMT
        DO I = 1,IMT
!lyc#ifdef SPMD          
!          IF(j_global(J)<=40) THEN
!            ccs1=catm(3)
!          ELSE IF(j_global(J)<=140) THEN
!            ccs1=catm(2)
!          ELSE
!            ccs1=catm(1)
!          ENDIF     
!#else          
!          IF(J<=40) THEN
!            ccs1=catm(3)
!          ELSE IF(J<=140) THEN
!            ccs1=catm(2)
!          ELSE
!            ccs1=catm(1)
!          ENDIF     
!#endif          
#ifdef SPMD          
!cm080415
!          IF(j_global(J)<=70) THEN
          IF(j_global(J)<=36) THEN
            ccs1=catm(3)
!          ELSE IF(j_global(J)<=111) THEN
          ELSE IF(j_global(J)<=56) THEN
            ccs1=catm(2)
          ELSE
            ccs1=catm(1)
          ENDIF     
#else          
!          IF(J<=70) THEN
          IF(J<=36) THEN
            ccs1=catm(3)
!          ELSE IF(J<=111) THEN
          ELSE IF(J<=56) THEN
            ccs1=catm(2)
          ELSE
            ccs1=catm(1)
          ENDIF     
!cm080415
#endif          
!cm080807,for the use of sea ice
!          ssfc(I,J) = overd*(ccs1-ptb(I,J,1,1))*vit(I,J,1)
#ifdef nc14wind
          ssfc(I,J) = (1.0-iceday(i,j))*overdwind(i,j)*(ccs1-ptb(I,J,1,1))*vit(I,J,1)
#else
          ssfc(I,J) = (1.0-iceday(i,j))*overd*(ccs1-ptb(I,J,1,1))*vit(I,J,1)
#endif
!cm080807
!ssfc(i,j)=0.
ENDDO
      ENDDO
#ifdef carbonDebug      
#ifdef SPMD
      PRINT*,'i,j,ssfc:',imt/2,jmt/2,ssfc(imt/2,jmt/2),'mytid=',mytid,"FLUX_PT"
#else        
      PRINT*,'i,j,ssfc:',imt/2,jmt/2,ssfc(imt/2,jmt/2),"FLUX_PT"
#endif
#endif
#endif
!     
!----xu----A--------------------------
!  checking the inorganic carbon cycle
!-------------------------------------
#ifdef carbonAbio
#ifdef SPMD
      PRINT*,'i,j,totup,ssfc,dpco2o,pco2o,dpco2a,tatmp,tctmp,tt,ss at imt/2 and jmt/2:',&
             imt/2,jmt/2,totup(imt/2,jmt/2),ssfc(imt/2,jmt/2),dpco2o(imt/2,jmt/2),pco2o(imt/2,jmt/2),&
             dpco2a(imt/2,jmt/2),tatmp(imt/2,jmt/2),tctmp(imt/2,jmt/2),tt(imt/2,jmt/2),&
             ss(imt/2,jmt/2),'mytid=',mytid,"FLUX_PT"
#else        
      PRINT*,'i,j,ssfc:',imt/2,jmt/2,ssfc(imt/2,jmt/2),"FLUX_PT"
#endif
#endif
!
! 
!lyc
!------------------------------------------------------------
!     SOLUBILITY COEFFICIENTS OF CFC11 & 12 AFTER WEISS(1974)
!                   UNIT: mol/liter/atm
!------------------------------------------------------------
#ifdef cfc
!---------------------------------------------
       DATA A1/ -229.9261   , -218.0971    /
       DATA A2/  319.6552   ,  298.9702    /
       DATA A3/  119.4471   ,  113.8049    /
       DATA A4/   -1.39165  ,   -1.39165   /
       DATA B1/   -0.142382 ,   -0.143566  /
       DATA B2/    0.091459 ,    0.091015  /
       DATA B3/   -0.0157274,   -0.0153924 /
!----------------------------------------------------------
!      COMPUTE SURFACE SATURATED CFC CONCENTRATION.
!      SOL (mol/liter/atm) : SOLUBILITY OF CFC11 & 12.
!      CCS (pico-mol/liter): SATURATED CFCs IN UPPER MOST BOXES.
!      T0(K), S0(psu)      : SURFACE TEMPERATURE & SALINITY.
!----------------------------------------------------------------
        L=1
!c$omp parallel do private (j,i,t0,s0,xsol,sol)
!$doacross local(j,i,t0,s0,xsol,sol) 
       DO  j = 1, Jmt
       DO  i = 1, Imt
        t0 = AT(i,j,1,1) + 273.15
         s0 = AT(i,j,1,2)*1000.+35.  
         if( vit(I,J,1)<0.5 ) cycle
         xsol = A1(L) + A2(L)*100.0/t0                &
            + A3(L)*LOG(t0/100.) + A4(L)*(t0/100.)**2 &
            + s0 *( B1(L) + B2(L)*t0/100.0 + B3(L)*(t0/100.)**2 )
        sol = EXP(xsol)
        if (j_global(j).le.91) then
        ccs(i,j) = sol * cfcatm(1) * vit(I,J,1)
         else 
         ccs(i,j) = sol * cfcatm(2) * vit(I,J,1)
        endif
       ENDDO
       ENDDO
!--------------------------------------
       call sgecfc
!--------------------------------------
!
! 5 year=5*3.1536E7 seconds       
      nn=1
      IF(NN.EQ.1) THEN
!14      overd=1.0/5./3.1536E7
      DO 120 J = 1,JMT
        DO 120 I = 1,IMT
!14          IF(J.LT.21) THEN
!14            ccs1=ccatm(3)
!14          ELSE IF(J.LE.69) THEN
!14            ccs1=ccatm(2)
!14          ELSE
!14            ccs1=ccatm(1)
!14          ENDIF     
!cc          ssfc(i,j) = (zk0(1)-zk0(2))/5./3.1536E7*(ccs1-ccp(i,j,1))
!cc     &                                           *vit(i,j,1)
!         if(itice(i,j).ne.1) then
         ssfc(I,J) = SGE(i,j)*ODZP(1)/360000.*(ccs(I,J)-ptb(I,J,1,1))*vit(i,j,1)
!         else
!         ssfc(i,j)=0.0
!         endif
120   CONTINUE
       write(*,*) 'ssfc=' ,ssfc(imt/2,jmt/2)
      ENDIF
#endif

      RETURN
      END SUBROUTINE FLUX_PT


SUBROUTINE FLUX_TA

!---------------------------------------------------------------------------------
!if carbonBio is defined , the sea surface TA will be modified by the relationship of TA and salinity
!the statistical relationship is from the work of Millero et al.,1998
!---------------------------------------------------------------------------------

!Author: Li Yangchun 15/05/2008


#include <def-undef.h>
#if(defined carbonBio) || (defined carbonAbio)
      USE param_mod
      USE pconst_mod
      USE tracer_mod
      USE carbon_mod
      USE forc_mod
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------------
      IMPLICIT NONE
      integer::ia,ja
      real,dimension(:,:,:),allocatable::res

!      real,dimension(imt,jmt,2)::res
!      real,dimension(imt,jmt)::nta_pt1
!      real,dimension(imt,jmt_global)::nta_io
      allocate(res(imt,jmt,2))


!--------------------------------------- 
      nta_pt(:,:)=0.0
#if (defined BOUNDARY)
    do j=1,jmt
	do i=1,imt
	    res(i,j,1)=restore(i,j,1,1)
	    res(i,j,2)=restore(i,j,1,2)*1000.0+35.0
	enddo
    enddo
#else
    do j=1,jmt
	do i=1,imt
	    res(i,j,1)=at(i,j,1,1)
	    res(i,j,2)=at(i,j,1,2)*1000.0+35.0
	enddo
    enddo
#endif
    do j=1,jmt
    do i=2,imm
#ifdef SPMD
	ja=j_global(j)
        ia=i_global(i)
#else
	ja=j
        ia=i
#endif
	if(vit(i,j,1)<0.5.or.res(i,j,2)<29.0.or.res(i,j,2)>37.0) cycle
!----------------------------------------------------------------------------------
! for the Arctic Ocean
!    if(basin(i,ja)==1) then
!	nta_pt(i,j)=2291-2.69*(res(i,j,1)-20.0)-0.046*(res(i,j,1)-20.0)**2
!    endif
!-----------------------------------------------------------------------------------
! for the Atlantic ocean and Indian ocean   
   if(basin(ia,ja)==2.and.lat(ja)>=30.0) then
	nta_pt(i,j)=2291-2.69*(res(i,j,1)-20.0)-0.046*(res(i,j,1)-20.0)**2
   endif
   if(basin(ia,ja)==2.and.lat(ja)<30.0.or.basin(ia,ja)==3) then
       nta_pt(i,j)= 2291
   endif
!-------------------------------------------------------------------------------------
!for the Pacific ocean
   if(basin(ia,ja)==4) then
     if(lat(ja)>=30.0) then
         nta_pt(i,j)=2300-7.00*(res(i,j,1)-20)-0.158*(res(i,j,1)-20)**2
     else if(lat(ja)>=-20.0.and.lat(ja)<=20.0.and.lon(ia)>=250.and.lon(ia)<=285) then
         nta_pt(i,j)=2300-2.94*(res(i,j,1)-29)-0.058*(res(i,j,1)-29)**2
     else if(lat(ja)<=10.0.and.lat(ja)>=-10.0.and.lon(ia)>=220.and.lon(ia)<=250) then
         nta_pt(i,j)=2300-2.94*(res(i,j,1)-29)-0.058*(res(i,j,1)-29)**2
     else
         nta_pt(i,j)=2300
     endif
    endif
!--------------------------------------------------------------------------------------
!for the southern ocean
     if(basin(ia,ja)==5.or.basin(ia,ja)==6) then
	  nta_pt(i,j)=2291-2.52*(res(i,j,1)-20)+0.056*(res(i,j,1)-20)**2
     endif
     if(basin(ia,ja)==7) then
	  nta_pt(i,j)=2300-2.52*(res(i,j,1)-20)+0.056*(res(i,j,1)-20)**2
     endif

 !    if(nta_pt(i,j)<1900) print *,'nta<1900', i,j_global(j),nta_pt(i,j)
     enddo
       if(nx_proc==1) then
        nta_pt(1,j)=nta_pt(imm,j)
        nta_pt(imt,j)=nta_pt(2,j)
       endif
     enddo
#ifdef SPMD
      call exch_boundary(nta_pt(1,1),1)
#endif
!    if(mytid==1) then
!	do j=1,jmt
!    write(99,*) (basin(i,ja),i=1,imt)
!        enddo
!    endif
!-----------------------------------------------------------
   do j=1,jmt
     do i=1,imt
!	 nta_pt1(i,j)=nta_pt(i,j)*vit(i,j,1)*res(i,j,2)/35.0
     IF (vit(i,j,1) > 0.5) then
     if(res(i,j,2)>29.and.res(i,j,2)<38.and.nta_pt(i,j)>0.0) THEN
!----------------------------------------------------
!xu
!----------------
#ifdef carbonAbio
!---------
!   SSFC is the observed TA for the calculation of pCO2sw
!-------------------------------------------------------
     ssfc(i,j)=nta_pt(i,j)*vit(i,j,1)*res(i,j,2)/35.0

    else
     ssfc(i,j)=2300
#else
     ssfc(i,j)=DZP(1)*3*GAMMA*(nta_pt(i,j)*vit(i,j,1)*            &
                res(i,j,2)/35.0-ptb(i,j,1,4))
#endif
    endif
     ENDIF
     enddo
    enddo

!      call local_to_global_4d(nta_pt1,nta_io,1,1)
!      if(mytid==0) then
!               write(99) ((nta_io(i,j),i=1,imt),j=1,jmt_global)
!              close(99)
!      endif

!     do j=1,jmt
!	 do i=1,imt
!	     ja=j_global(j)
!	     if(i==141.and.ja==20) print *,ssfc(i,j),nta_pt(i,j)*vit(i,j,1)*            &
!                (res(i,j,2)*1000.0+35.0)/35.0,res(i,j,1:2)
		 
!         enddo
!     enddo
!----------------------------------

     deallocate(res)
#endif  

END SUBROUTINE FLUX_TA


! CVS: $Id: pco2.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  SUBROUTINE OPCO2(pco2o,tta,ttc,tt,ss)
!========================
! PCO2
!---------------------------------------------------------------------
!
! purpose: calculate partial pressure of carbon dioxide with TA, T and S
!
! author: Zhao Liang@lapc 2004/05/26
!          (original author: Li Y CH and Xu Y F)
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod, ONLY: imt,jmt,jsm,jem,imm,mytid,nx_proc 
      USE pconst_mod, ONLY: vit,i_global,j_global 
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!#ifdef SPMD      
!#include <mpif.h>
!#endif      
!      
      REAL,DIMENSION(imt,jmt)::tta,ttc,tt,ss
      REAL::ta,tc,t,s
      REAL,DIMENSION(imt,jmt)::pco2o
!      
      REAL::tp,tsi,tb
      REAL::ion
      REAL*8::k0,pk1,pk2,k1,k2,fh,kb,ksi,kw,kp1,kp2,kp3
      REAL*8::fh1
      REAL*8::ac,ab,asi,ap,aw,ah,ah1

      INTEGER::i,j,nn

      tp=1.5
      tsi=5.0
      tb=410.6

      pco2o(:,:)=0.0
#ifdef carbonDebug
#ifdef SPMD
            print*, 'subroutine pco2 -----------------------------------------------','mytid=',mytid
#else
            print*, 'subroutine pco2 ------------------------------------------------'
#endif
#endif

      DO j=2,jem
        DO i=2,imm
          IF(vit(i,j,1) < 0.5) CYCLE
          t=tt(i,j)+273.15
          s=ss(i,j)*1000.0+35.0
          ta=tta(i,j)
          tc=ttc(i,j)

!#ifdef carbonDebug
          if(t>270.0 .and. t<313.0) then
            continue
          else
#ifdef SPMD
            print*, 'abnormal temperature i,j,t:',i,j,t,'mytid=',mytid
#else
            print*,'abnormal temperature i,j,t:',i,j,t
#endif
          endif
          if(isnan(tc))then
          print*,'DIC is nan',i_global(i),j_global(j),tc,ta
         stop
          endif
          if(ta<1600.0.or.ta>3000.0 ) then
#ifdef SPMD
            print*, 'abnormal TA i,j,ta:',i,j_global(j),ta,s
#else
            print*,'abnormal TA i,j,ta:',i,j,ta
#endif
          endif
          if(tc<1500.0 .or. tc>2500.0) then
#ifdef SPMD
            print*, 'abnormal TC i,j,tc:',i,j_global(j),tc
#else
            print*,'abnormal TC i,j,tc:',i,j,tc
#endif
          endif
!#endif
!          tb=410.6*s/35.
 
! the unit of ko is moles/(kg*atm)or umols/(kg*uatm)
          k0=exp(-60.2409+93.4517*(100./t)+23.3585*log(t/100.)+ &
                s*(0.023517-0.023656*(t/100.)+0.0047036*(t/100.)**2))

          pk1=845.0/t+3.248-0.0098*s+0.000087*s**2
          pk2=1377.3/t+4.824-0.0185*s+0.000122*s**2
          k1=10.0**(-pk1)*10**6
          k2=10.0**(-pk2)*10**6
          fh=1.29-0.00204*T+4.6E-4*S**2-1.48E-6*S**2*T
          fh1=fh*10.0**6

! the unit is molality
! dissociation of borate :
          kb=exp((-8966.90-2890.51*s**0.5-77.942*s+1.726*s**1.5- &
             0.0993*s**2)/t+(148.0248+137.194*s**0.5+1.62247*s)+ &
             (-24.4344-25.085*s**0.5-0.2474*s)*log(t)+0.053105*s**0.5*t) &
             *(1-s*0.001005)*10.0**6
! dissociation of silicic
          ion=0.7
          ksi=exp(117.40-8904.2/t-19.334*log(t)+(3.5913-458.79/t)*ion**0.5+ &
                 (-1.5998+188.74/t)*ion+ &
                 (0.07871-12.1652/t)*ion**2)*(1-s*0.001005)*10.0**6

! the unit is umol/(kg soln):
! dissociation of water
          kw=exp(148.9802-13847.26/t-23.6521*log(t)+(-5.977+118.67/t+ &
                 1.0495*log(t))*s**0.5-0.01615*s)*10.0**6

! dissociation of phosphoric
          kp1=exp(115.54-4576.752/t-18.453*log(t)+(0.69171-106.736/t)*s**0.5+ &
                 (-0.01844-0.65643/t)*s)*10.0**6
          kp2=exp(172.1033-8814.715/t-27.927*log(t)+(1.3566-160.340/t)*s**0.5+ &
                 (-0.05778+0.37335/t)*s)*10.0**6
          kp3=exp(-18.126-3070.75/t+(2.81197+17.27039/t)*s**0.5+(-0.09984- &
                  44.99486/t)*s)*10.0**6

! calculate:
          ah=0.005
	  nn=0
 
 doo:     DO
     !       ac=tc*(k1*ah+2*k1*k2)/(ah**2+k1*ah+k1*k2)
            ab=tb*kb/(ah+kb)
            asi=tsi*ksi/(ah+ksi)
            ap=tp*(kp2*ah+2*kp2*kp3)/(ah**2+kp2*ah+kp2*kp3)
            aw=kw*fh1/ah-ah/fh1
     !      ab=ta-ac-asi-ap-aw
            ac=ta-ab-asi-ap-aw
     !      ah1=tb*kb/ab-kb

            ah1=k1/(2.0*ac)*((tc-ac)+sqrt((ac-tc)**2+4.0*ac*k2*(2*tc-ac)/k1))
            IF(abs((ah1-ah)/ah) > 1.e-4) THEN
              ah=ah1

!            IF(abs(ah1-ah) > 1.e-6) THEN
!              ah=ah+(ah1-ah)/10.0
	      if(ah<0.0) then
		  print *, 'ah is error,i,j', i,j_global(j),ah,'tt',t,'ss',s,'ta',ta,'tc',tc
		stop
	      endif
              nn=nn+1
	      if(nn>50) then
		  print*,"the endless cycle",i,j_global(j),ah,"tt,", t,"ss,",s,"ta",ta,"tc",tc
		  stop
	      endif
	      CYCLE doo
            ELSE
              EXIT doo
            ENDIF
          ENDDO doo  
      
          pco2o(i,j)=tc*ah**2/(k0*(ah**2+k1*ah+k1*k2))
        
        ENDDO
       if(nx_proc==1) then
        call exchange_2d(pco2o,1,1)
         pco2o(1,j)  =pco2o(imm,j)
         pco2o(imt,j)=pco2o(2,j)
       endif
      ENDDO

#ifdef SPMD
        call exch_boundary(pco2o,1)
#endif
!      
!!!      CALL exchange_2d(pco2o)
!      
#ifdef carbonDebug
#ifdef SPMD
      print*, "i,j,pco2o:",imt/2,jmt/2,pco2o(imt/2,jmt/2), &
                  "mytid=",mytid,"PCO2"
#else
      print*, "i,j,pco2o:",imt/2,jmt/2,pco2o(imt/2,jmt/2),"PCO2"
#endif
#endif
#ifdef carbonDebug
#ifdef SPMD
            print*, 'subroutine pco2 end-----------------','mytid=',mytid
#else
            print*, 'subroutine pco2 end-----------------'
#endif
#endif
  END SUBROUTINE OPCO2

! CVS: $Id: readybio.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  SUBROUTINE READYBIO(C2DTTS)
!========================
! BIOSOURCE
!---------------------------------------------------------------------
!
! purpose: prepararion of calculating biological source
!          calculate A0, A1, A2, B0, B1, B2 and C
!
! author: Zhao Liang@lapc 2004/05/30
! the iron is added by lyc in 05/2013
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
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod
      USE tracer_mod
      USE carbon_mod
      USE cforce_mod
      USE grids_pt_mod
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
     
!
!---------------------------------------------------------------------
      IMPLICIT NONE
#include <netcdf.inc>      
!#ifdef SPMD      
!#include <mpif.h>
!#endif      
!      
#ifdef carbonBio

      real,parameter::S0=35.0
!      real,dimension(km)::zb,zm
      real::dtts,scalet,zer,r_cpr
!lyc2009.07.07      real::inta1,sr
      real::sr,ze
      real,dimension(:,:),allocatable::inta1
      real,dimension(:,:,:),allocatable::c2_b
!     kappa_a(1:km)
#ifdef anderson1995
      real(kind=8)::prodze,ldoczb,prod0,ldoc0
#endif
!     r_cap(1:km),delta_a(1:km)
      real,dimension(imt,jmt,km)::atmp
      real(kind=8),dimension(km)::intatmp,inta0
!      real(kind=8),dimension(km)::sdxdy
      real::tmpf,tmpd
      real(kind=8),dimension(km)::intatmp0,inta00
!      real(kind=8),dimension(km)::sdxdy0
!      
! add new compent by lyc
      real(kind=8):: va0tmp,va0sum,vctmp,vcsum,delcaco3,deltaa0
      real(kind=8):: vb0tmp,vb0sum,deltab0
!      integer:: kb0
!lyc 2013.05
      real::Fe_b,Fe_c,Fe_free,scb,sc
!      real,parameter::feb=0.00384 ng-1cm-1 !moore 2008
      real,parameter::feb=0.015/86400/365!kg/ug/ml
      real::nut_f  !nutrient forcing
      real::nut_b  !basic nutrient for bio-production
      real,dimension(imt,jmt,kmp1)::dust_out
      real,dimension(imt,jmt)::dust_soft,dust_hard
      real::sqfe
      real::scal

      real::c2dtts
      allocate(inta1(imt,jmt),c2_b(imt,jmt,km))
!lyc2013.05      kb0=24
      scalet=1.0/(tau_b*86400.0)
      ze=zb(kmmix)
      zer=1.0/ze
      r_cpr=1.0/r_cp

!
!      if (isp >= 1)then
      IF (mod(ISP,15)/=0) THEN
         dtts = dts * 2.0
      else
         dtts = dts
      endif
!lyc 2013.05      
!      zb(1)=dzp(1)
!      do k=2,km
!        zb(k)=zb(k-1)+dzp(k)
!      enddo
!      do k=1,km
!        zm(k)=-zkt(k)
!      enddo
!!lyc-----------------------------------------
!!calculate the volume of sea water under 1500m
!   if(mytid==0) then
!     vsea2000=0.0
!     do k=kb0,km
!       do j=1,jmt_global
!         do i=1,imm
!         if(vit_global(i,j,k)<0.5) cycle
!          vsea2000=vsea2000+dxdyt_global(j)*dzp(k)
!         enddo
!       enddo
!     enddo
!   endif
!lyc 2013.05
!---------------------------------------------------------------------
!     This subroutine is used to calculate Ax, Bx and C
!---------------------------------------------------------------------
!
!     calculate A1 - new primary production        
#ifdef murnane1999        
      do k=1,kmmix
        do j=1,jmt
          do i=2,imm
            if(vit(i,j,k) < 0.5) cycle
            if (ptb(i,j,k,2) > po4force(i,j,k).and.po4force(i,j,k)>0.1E-4) then
              a1_b(i,j,k)=scalet*(pt(i,j,k,2)-po4force(i,j,k)/1.025)
            else
              a1_b(i,j,k)=0.0
            endif
          enddo
        enddo
      enddo
!      if(mytid==0) print *, 'a1_b(85,2,1)',a1_b(85,2,1)
#ifdef carbonDebug
#ifdef SPMD
      print*, 'subroutine readybio-----------------------------------------------------------','mytid=',mytid
      print*, 'i,j,1,ptb(po4),po4force',imt/2,jmt/2,ptb(imt/2,jmt/2,1,2),po4force(imt/2,jmt/2,1), &
              'mytid=',mytid      
#else
      print*, 'subroutine readybio-----------------------------------------------------------'
      print*, 'i,j,1,ptb(po4),po4force',imt/2,jmt/2,ptb(imt/2,jmt/2,1,2),po4force(imt/2,jmt/2,1)
#endif
#endif

#else
!the new production is a funtion of I,PO4
!the new production is a funtion of I,PO4
!  EP=r0*Lf*[po4]**2/(hf+[po4])*exp(-kz)*dzp*(t+2.0)/(t+10.)
!-----------------------------------------------------------------------------------------------
      nut_f=0.0
      do k=1,kmmix
        do j=1,jmt
          do i=2,imm
            a1_b(i,j,k)=0.0
            if(vit(i,j,k) < 0.5.or.pt(i,j,k,2)<0.1E-6) cycle
#ifdef Felimit
            if (pt(i,j,k,6)>0.0) then
            if(pt(i,j,k,6)<0.2*1.0E-3) then 
            r_fec(i,j)=2.0*1.0E-6
            r_fep(i,j)=r_fec(i,j)*r_cp
            endif
            nut_f=min(pt(i,j,k,2)/(hf/1.025+pt(i,j,k,2)),pt(i,j,k,6)/(hf_Fe/1.025+pt(i,j,k,6)))
            nut_b=min(pt(i,j,k,2),pt(i,j,k,6)/r_fep(i,j))
            else
            nut_f=0.0
            nut_b=0.0
            endif
#else
            nut_f=pt(i,j,k,2)/(hf/1.025+pt(i,j,k,2))
            nut_b=pt(i,j,k,2)
#endif
            a1_b(i,j,k)=ar0*sint(j)*nut_b*nut_f*(at(i,j,k,1)+2.0)/(at(i,j,k,1)+10.0)  &
                        *exp(-zm(k)/zdc)
#ifdef Felimit
            if((pt(i,j,k,2)-a1_b(i,j,k)*c2dtts)<0.0.or.(pt(i,j,k,6)-a1_b(i,j,k)*c2dtts*r_fep(i,j))<0.0) then
            a1_b(i,j,k)=min(pt(i,j,k,2)/c2dtts,pt(i,j,k,6)/c2dtts/r_fep(i,j))
            endif
#endif
          enddo
        enddo
      enddo
#endif
#ifdef SPMD
      call exch_boundary(a1_b,km)
#endif
!     calculate A2 - flux of POP below the euphotic zone      
      do j=1,jmt
        do i=2,imm
        if(vit(i,j,1) < 0.5) cycle
          inta1(i,j)=0.0
          do k=1,kmmix
            inta1(i,j)=inta1(i,j)+a1_b(i,j,k)*dzp(k)
          enddo
          do k=kmmix,km
            a2_b(i,j,k)=(1.0-sigma_b)*inta1(i,j)*(zb(k)*zer)**(POP_k)
          enddo
        enddo
      enddo
!     calculate A0       
      do k=1,kmmix
        do j=1,jmt
          do i=2,imm
            if(vit(i,j,k) < 0.5) cycle
!lyc20140404
             if(k<ITNU(I,J))then
            a0_b(i,j,k)=(1.0-sigma_b)*a1_b(i,j,k)
             else
            a0_b(i,j,k)=-inta1(i,j)/dzp(k)
             endif
          enddo
        enddo
      enddo
      do k=kmmix+1,km
        do j=1,jmt
          do i=2,imm
            if(vit(i,j,k) < 0.5) cycle
!lyc20140404
            if(k< ITNU(i,j)) then
            a0_b(i,j,k)=-(a2_b(i,j,k-1)-a2_b(i,j,k))/dzp(k)
            else
            a0_b(i,j,k)=-a2_b(i,j,k-1)/dzp(k)
            endif
          enddo
        enddo
      enddo
#ifdef SPMD
     call exch_boundary(a0_b,km)
     call exch_boundary(a2_b,km)
#endif
!lyc--------------------------
! calculate delta a0
!lyc20140404
!      va0sum=0.0
!      do k=1,km
!        do j=2,jem
!          do  i=2,imm
!            if(vit(i,j,k)<0.5) cycle
!            va0sum=va0sum+a0_b(i,j,k)*dxdyt(j)*dzp(k)
!          enddo
!        enddo
!      enddo
!       call mpi_barrier(mpi_comm_ocn,ierr)
!       call mpi_reduce(va0sum,va0tmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)
!       if(mytid==0) then
!          deltaa0=va0tmp/vsea2000
!       endif
!       call mpi_bcast(deltaa0,1,mpi_real8,0,mpi_comm_ocn,ierr)
!       do k=kb0,km
!         do j=1,jet
!           do i=2,imm
!           if(vit(i,j,k)<0.5) cycle
!            a0_b(i,j,k)=a0_b(i,j,k)-deltaa0
!           enddo
!         enddo
!       enddo
   
!     calculate kappa_b(km)
#ifdef anderson1995
!      kappa_b(1:km)=1.0/(11.2*365.0*86400.0)
       prodze=0.0
       do k=1,kmmix
         do j=jsm,jem
           do i=2,imm
             if(vit(i,j,k) < 0.5) cycle  
             prodze=prodze+a1_b(i,j,k)*dxdyt(j)*dzp(k)
           enddo
         enddo
       enddo
       prodze=r_cp*sigma_b*prodze
       
       ldoczb=0.0
       do k=kmmix+1,km
         do j=jsm,jem
           do i=2,imm
             if(vit(i,j,k) < 0.5) cycle  
             ldoczb=ldoczb+pt(i,j,k,3)*dxdyt(j)*dzp(k)*exp((ze-zm(k))/750.0)
           enddo
         enddo
       enddo
#ifdef SPMD
       prod0=0.0
       ldoc0=0.0
       call mpi_barrier(mpi_comm_ocn,ierr)
       call mpi_reduce(prodze,prod0,1,mpi_double_precision,mpi_sum,0,mpi_comm_ocn,ierr)
       call mpi_reduce(ldoczb,ldoc0,1,mpi_double_precision,mpi_sum,0,mpi_comm_ocn,ierr)

       
       if(mytid==0) then
         kappa_b(1:km)=0.0
         do k=kmmix+1,km
           kappa_b(k)=prod0/ldoc0*exp((ze-zm(k))/750.0)
	 enddo
#ifdef carbonDebug
       print*,'prod0=',prod0,'ldoc0=',ldoc0,'mytid=',mytid
#endif
       endif
       call mpi_barrier(mpi_comm_ocn,ierr)
       call mpi_bcast(kappa_b,km,mpi_real8,0,mpi_comm_ocn,ierr)
#else       
       kappa_b(1:km)=0.0
       do k=kmmix+1,km
         kappa_b(k)=prod0/ldoc0*exp((ze-zm(k))/750.0)
       enddo
#endif
#endif
!     calculate B1 - production of LDOC in the euphotic zone       
      do k=1,kmmix
        do j=1,jmt
          do i=2,imm
            if(vit(i,j,k) < 0.5) cycle
            b1_b(i,j,k)=r_cp*sigma_b*a1_b(i,j,k)
          enddo
        enddo
      enddo
!     calculate B2 - LDOC remineralization      
      do k=1,km
        do j=1,jmt
          do i=2,imm
            if(vit(i,j,k) < 0.5) cycle
            if(pt(i,j,k,3) > 0.0) then
              b2_b(i,j,k)=kappa_b(k)*pt(i,j,k,3)
            else
              b2_b(i,j,k)=0.0  
            endif
          enddo
        enddo
      enddo
!     calculate B0       
      b0_b(:,:,:)=b1_b(:,:,:)-b2_b(:,:,:)
#ifdef SPMD
      CALL EXCH_BOUNDARY(b0_b,km)
      CALL EXCH_BOUNDARY(b1_b,km)
      CALL EXCH_BOUNDARY(b2_b,km)
#endif
!lyc--------------------------
      vb0sum=0.0
      do k=1,km
        do j=2,jem
          do  i=2,imm
            if(vit(i,j,k)<0.5) cycle
            vb0sum=vb0sum+b0_b(i,j,k)*dxdyt(j)*dzp(k)
          enddo
        enddo
      enddo
       call mpi_barrier(mpi_comm_ocn,ierr)
       call mpi_reduce(vb0sum,vb0tmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)
       if(mytid==0) then
          deltab0=vb0tmp/vsea2000
       endif
       call mpi_bcast(deltab0,1,mpi_real8,0,mpi_comm_ocn,ierr)
       do k=kb0,km
         do j=1,jet
           do i=1,imt
           if(vit(i,j,k)<0.5) cycle
            b0_b(i,j,k)=b0_b(i,j,k)-deltab0
           enddo
         enddo
       enddo
!lyc----------------------------------------
! calculate the calcite
#ifdef progca
!     calculate r_cap(1:km)
      do k=1,km
        do j=1,jmt
          do i=2,imm
            if(vit(i,j,k) < 0.5) cycle  
            sr=at(i,j,k,2)*1000.0+35.0
            atmp(i,j,k)=(pt(i,j,k,4)+r_np*dtts* &
                        (a0_b(i,j,k)+b0_b(i,j,k)*r_cpr))!*S0/sr
          enddo
        enddo
      enddo
#ifdef carbonDebug
#ifdef SPMD
      print*, '--------------------------------------------------------------------------------------------'
      print*,'i,j,1,ptb(TA),a0,b0,at(S)',imt/2,jmt/2,ptb(imt/2,jmt/2,1,4), &
              a0_b(imt/2,jmt/2,1),b0_b(imt/2,jmt/2,1),at(imt/2,jmt/2,1,2), &
              'mytid=',mytid,'READYBIO'
      print*,'i,j,10,ptb(TA),a0,b0,at(S)',imt/2,jmt/2,ptb(imt/2,jmt/2,10,4), &
              a0_b(imt/2,jmt/2,10),b0_b(imt/2,jmt/2,10),at(imt/2,jmt/2,10,2), &
              'mytid=',mytid,'READYBIO'
      print*,'i,j,1,kappa,a1,b1',imt/2,jmt/2,kappa_b(1),a1_b(imt/2,jmt/2,1), &
              b1_b(imt/2,jmt/2,1),'mytid=',mytid
      print*,'i,j,10,kappa,a2,b2',imt/2,jmt/2,kappa_b(10),a2_b(imt/2,jmt/2,10), &
              b2_b(imt/2,jmt/2,10),'mytid=',mytid
#else
      print*, '------------------------------------------------------------------------------------------'
      print*,'i,j,1,ptb(TA),a0,b0,at(S)',imt/2,jmt/2,ptb(imt/2,jmt/2,1,4), &
              a0_b(imt/2,jmt/2,1),b0_b(imt/2,jmt/2,1),at(imt/2,jmt/2,1,2), &
              'READYBIO'
      print*,'i,j,10,ptb(TA),a0,b0,at(S)',imt/2,jmt/2,ptb(imt/2,jmt/2,10,4), &
              a0_b(imt/2,jmt/2,10),b0_b(imt/2,jmt/2,10),at(imt/2,jmt/2,10,2), &
              'READYBIO'
      print*,'i,j,1,kappa,a1,b1',imt/2,jmt/2,kappa_b(1),a1_b(imt/2,jmt/2,1), &
              b1_b(imt/2,jmt/2,1)
      print*,'i,j,10,kappa,a2,b2',imt/2,jmt/2,kappa_b(10),a2_b(imt/2,jmt/2,10), &
              b2_b(imt/2,jmt/2,10)
#endif
#endif
      do k=1,km
       intatmp(k)=0.0
       inta0(k)=0.0
       sdxdy(k)=0.0
       intatmp0(k)=0.0
       inta00(k)=0.0
       sdxdy0(k)=0.0
       do j=jsm,jem
         do i=2,imm
           if(vit(i,j,k) < 0.5) cycle
           intatmp(k)=intatmp(k)+atmp(i,j,k)*dxdyt(j)
           inta0(k)=inta0(k)+a0_b(i,j,k)*dxdyt(j)
           sdxdy(k)=sdxdy(k)+dxdyt(j)
         enddo
       enddo
      enddo

#ifdef SPMD      
      call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_reduce(intatmp,intatmp0,km,mpi_double_precision,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(inta0,inta00,km,mpi_double_precision,mpi_sum,0,mpi_comm_ocn,ierr)
      if(mytid==0) then
#endif
        do k=1,km
#ifdef SPMD      
          intatmp(k)=intatmp0(k)/sdxdy0(k)  
          inta0(k)=inta00(k)/sdxdy0(k)  
#else
          intatmp(k)=intatmp(k)/sdxdy(k)  
          inta0(k)=inta0(k)/sdxdy(k)  
#endif
          if(abs(inta0(k)) > 0.0) then
            r_cap(k)=(intatmp(k)-taforce(k))/inta0(k)*0.5/dtts
          else
            r_cap(k)=0.0   
          endif
        enddo
        
!     calculate delta_a(1:km)
        delta_a(1:km)=0.0
        do k=1,kmmix
          if(intatmp(k) > taforce(k)) delta_a(k)=1.0
        enddo
        do k=kmmix+1,km
          if(intatmp(k) < taforce(k)) delta_a(k)=1.0
        enddo
        tmpf=0.0
        do k=1,kmmix
          if(delta_a(k) < 0.5) cycle
          tmpf=tmpf+(intatmp(k)-taforce(k))
        enddo
        tmpd=0.0
        do k=kmmix+1,km
          if(delta_a(k) < 0.5) cycle
          tmpd=tmpd+(taforce(k)-intatmp(k))
          if(tmpd > tmpf) then
            delta_a(k)=(tmpf-(tmpd-(taforce(k)-intatmp(k))))/(taforce(k)-intatmp(k))
#ifdef SPMD          
            if(delta_a(k) > 1.0) print*,"delta_a error!",delta_a(k),"mytid=",mytid,"READYBIO"
#else
            if(delta_a(k) > 1.0) print*,"delta_a error!",delta_a(k),"READYBIO"
#endif
            if(k<km) delta_a(k+1:km)=0.0

            exit
          endif
        enddo
#ifdef SPMD      
#ifdef carbonDebug 
#ifdef SPMD
        print*, '--------------------------------------------------------------------------------------------'
        print*, 'intatmp=',intatmp(1:km),'mytid=',mytid
        print*, 'taforce=',taforce(1:km),'mytid=',mytid
        print*, 'inta0=',inta0(1:km),'mytid=',mytid
        print*, 'r_cap=',r_cap(1:km),'mytid=',mytid
        print*, 'delta_a=',delta_a(1:km),'mytid=',mytid
#else
        print*, '--------------------------------------------------------------------------------------------'
        print*, 'intatmp=',intatmp(1:km)
        print*, 'taforce=',taforce(1:km)
        print*, 'inta0=',inta0(1:km)
        print*, 'r_cap=',r_cap(1:km)
        print*, 'delta_a=',delta_a(1:km)
#endif
#endif
      endif

      call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(r_cap,km,mpi_real8,0,mpi_comm_ocn,ierr)
      call mpi_bcast(delta_a,km,mpi_real8,0,mpi_comm_ocn,ierr)
#endif
      
!     calculate C       
      do k=1,km
        do j=1,jmt
          do i=2,imm
            if(vit(i,j,k) < 0.5) cycle
            c_b(i,j,k)=delta_a(k)*r_cap(k)*a0_b(i,j,k)
          enddo
        enddo
      enddo
#else
     do k=1,kmmix
        do j=1,jmt
           do i=2,imm
            if(vit(i,j,k)<0.5) cycle
            c_b(i,j,k)=r_rain*r_cp*a0_b(i,j,k)
           enddo
        enddo
     enddo
!----------------------------------------
! the calcite downward flux
     do k=kmmix,km
       do j=1,jmt
         do i=2,imm
            if(vit(i,j,k)<0.5) cycle
            c2_b(i,j,k)=r_rain*r_cp*(1-sigma_b)*inta1(i,j)*exp(-(zb(k)-ze)/d_cal)
         enddo
       enddo
      enddo   
     do k=kmmix+1,km
       do j=1,jmt
         do i=2,imm
           if(vit(i,j,k)<0.5) cycle
!lyc20140407
           if(k<itnu(i,j))then
           c_b(i,j,k)=-(c2_b(i,j,k-1)-c2_b(i,j,k))/dzp(k)
           else
           c_b(i,j,k)=-c2_b(i,j,k-1)/dzp(k)
           endif
         enddo
       enddo
     enddo 
#endif
#ifdef SPMD
     call exch_boundary(c_b,km)
#endif
!       
!lyc---------------------
!calculate deltacaco3
!lyc20140407
!     vcsum=0.0
!     do k=1,km
!      do j=2,jem
!        do i=2,imm
!         if(vit(i,j,k)<0.5) cycle
!          vcsum=vcsum+c_b(i,j,k)*dxdyt(j)*dzp(k)
!        enddo
!      enddo
!     enddo
!         call mpi_barrier(mpi_comm_ocn,ierr)
!	 call mpi_reduce(vcsum,vctmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)    
!      	if(mytid==0) then
!         delcaco3=vctmp/vsea2000
!        endif  
!         call mpi_bcast(delcaco3,1,mpi_real8,0,mpi_comm_ocn,ierr)
!    do k=kb0,km
!      do j=1,jet
!        do i=2,imm
!        if(vit(i,j,k)<0.5) cycle
!        c_b(i,j,k)=c_b(i,j,k)-delcaco3
!        if(isnan(c_b(i,j,k))) print *,'c_b is error',i,j_global(j),k,c_b(i,j,k)
!        enddo
!      enddo
!    enddo
!---------------------------------------------------------------------------
! Iron cycle begin
!lyc 2013.0.5
!---------------------------------------------------------------------------
!calculate the dust flux
!---------------------------------------------------------------------------
!dust_in is the dust-bioavailable

   dust_in=0.0
    do j=1,jmt
      do i=2,imm
      dust_soft(i,j)=0.03*dust_f(i,j)
      dust_hard(i,j)=0.97*dust_f(i,j)
      dust_out(i,j,1)=dust_f(i,j)
      do k=2,km+1
      dust_out(i,j,k)=dust_soft(i,j)*exp(-zb(k-1)/600.0)+dust_hard(i,j)*exp(-zb(k-1)/120000.0)
      if(dust_out(i,j,k)<0.0.or.dust_out(i,j,k)>1.0E-7) then
       write(6,*) 'dust_out is error',i,j,k,dust_out(i,j,k),'dust',dust_f(i,j)
      endif
      enddo
      do k=1,km
      dust_in(i,j,k)=(dust_out(i,j,k)-dust_out(i,j,k+1))/dzp(k)
      enddo
      enddo
    enddo

!calculate the Fe_free
!1.0E-3  is the constant ligand concentration(umol/l)
!3.33E-6 is 1/Kl, Kl=300L/nmol  (umol/l)
!Fe_free is the positive root of the quadratic
![Fe_free]**2+b*[Fe_free]+c=0
!b=L+1/Kl-[Fe],L=1.0E-3 umol/l
!c=3.33E-6*[Fe]
!umol/l=1/1.025 umol/kg
!the method is similar with NCOM1.4
    do k=1,km
     do j=jsm,jem
      do i=2,imm
        Fe_scav(i,j,k)=0.0
      enddo
     enddo
    enddo
    P_fe(:,:,:)=0.0
    
    do k=1,km
       do j=jsm,jem
         do i=2,imm
         if(vit(i,j,k)<0.5.or.pt(i,j,k,6)<0.0) cycle
#ifdef scav_moore08
!for ug/kg*m/s 
        IF(FE_FLUX_DATA) THEN
        SCAL=55.847
        ELSE
        SCAL=1
        ENDIF
        scb=Feb*(r_cp*a2_b(i,j,k)*6.0*12.01+c2_b(i,j,k)*100.0+dust_out(i,j,K+1)*SCAL)
!for umol/kg*m/s
!         scb=Feb*(a2_b(i,j,k)*6.0+c2_b(i,j,k)+dust_out(i,j,K+1))
         if(pt(i,j,k,6)>0.6*1.0E-3) then
         sc=scb+(pt(i,j,k,6)-0.6*1.0E-3)*0.00904
!         else if(pt(i,j,k,6)<0.0005) then
!         sc=scb*(pt(i,j,k,6)/0.0005)
         else
         sc=scb
         endif
         Fe_scav(i,j,k)=sc*pt(i,j,k,6)     
#else
         Fe_free=0.0
         Fe_b=(1.0E-3+3.33E-6)/1.025-pt(i,j,k,6)
         Fe_c=-3.33E-6/1.025*pt(i,j,k,6)
         sqfe=Fe_b**2-4*Fe_c
         if(sqfe<0.0) then
          print*,'i,j,k', i,j_global(j),k,'sqfe is error',sqfe,pt(i,j,k,6),ptb(i,j,k,6)
          stop
         endif 
         if(Fe_b<0.0) then 
          Fe_free=0.5*(-Fe_b+sqrt(sqfe))
         else
          Fe_free=2.0*Fe_c/(-Fe_b-sqrt(sqfe))
         endif
         if(isnan(Fe_free)) then
          print*,'i,j,k',II, i_global(i),j_global(j),k,'Fe_free is error',Fe_free,pt(i,j,k,6),pt(i,j,k,1:3),dust_in(i,j,k)
          stop
         endif 
         Fe_scav(i,j,k)=Fe_scav_prof(k)*Fe_free
#endif
         if(isnan(Fe_scav(i,j,k))) then
          print*,'i,j,k', i,j_global(j),k,'Fe_scav is error',Fe_scav,Fe_free,Fe_scav_prof(k)
          stop
         endif 
        enddo
       enddo
     enddo

!     TFe_scav=0.0
!     do j=jsm,jem
!       do i=2,imm
!       TFe_scav(i,j,1)=Fe_scav(i,j,1)*dzp(1)
!       enddo
!     enddo
!
!     do k=2,km
!        do j=jsm,jem
!          do i=2,imm
!          TFe_scav(i,j,k)=TFe_scav(i,j,k-1)+Fe_scav(i,j,k)*dzp(k)
!          enddo
!        enddo
!     enddo
!------------------------------------------------------     
     Fe_source=0.0
     do k=1,kmmix
      do j=jsm,jem
       do i=2,imm
        if(vit(i,j,k)<0.5) cycle
         Fe_source(i,j,k)=-r_fep(i,j)*a1_b(i,j,k)-Fe_scav(i,j,k)
       enddo
      enddo
     enddo
     do k=2,kmmix+1
      do j=jsm,jem
       do i=2,imm
         P_fe(i,j,k)=P_fe(i,j,k-1)+(R_FeP(i,j)*(1.0-sigma_b)*a1_b(i,j,k-1)+0.6*Fe_scav(i,j,k-1))*dzp(k-1)
       enddo
      enddo
      enddo
     do k=kmmix+1,km
       do j=jsm,jem
        do i=2,imm
        if(vit(i,j,k)<0.5) cycle
        Fe_source(i,j,k)=-r_fec(i,j)*b0_b(i,j,k)+P_fe(i,j,k)*ReFe(k)/dzp(k)-Fe_scav(i,j,k)
        P_Fe(i,j,k+1)=P_Fe(i,j,k)-P_Fe(i,j,k)*ReFe(k)+dzp(k)*0.6*Fe_scav(i,j,k)
        enddo
       enddo
     enddo     
      if(ii==nss.and.mytid==50) then
       print *,'fe_source for b0_b,flux,scav,dust_in:',imt/2,jmt/2,'vit:',vit(imt/2,jmt/2,12),Fe_source(imt/2,jmt/2,12),-r_fec(imt/2,jmt/2),b0_b(imt/2,jmt/2,12),P_fe(imt/2,jmt/2,12)*refe(12)/dzp(12),-fe_scav(imt/2,jmt/2,12),dust_in(imt/2,jmt/2,12),pt(imt/2,jmt/2,12,6)
     endif
     
     do k=kb0,km
       do j=jsm,jem
         do i=2,imm
        if(vit(i,j,k)<0.5) cycle
        Fe_source(i,j,k)=Fe_source(i,j,k)+r_feP(i,j)*deltaa0-r_feC(i,j)*deltab0
        enddo
       enddo
     enddo

!      call exchange_3d(Fe_source,km)

!-------add the dust_bio to fe-source------------
!dust remin gDust=0.035/55.847*1.0e6->mmolfe
!0.035  dust iron content is 3.5% iron by weight (Moore et al.,2004)
!55.847 g/mol fe
!1.0e6  mol->mmol
IF(DUST_DATA) THEN
    SCAL=0.035/55.847 !
ELSE
    SCAL=1
ENDIF
    dust_in =0.0  !for just csm1_bgc
     do k=1,km
        do j=1,jmt
          do i=2,imm
          Fe_source(i,j,k)=Fe_source(i,j,k)+dust_in(i,j,k)*SCAL
          enddo
        enddo
      enddo
!----source from sediments--------------------------------------------
     do j=1,jmt
       do i=2,imm
        if(kmt(i,j)>1.and.kmt(i,j)<=22)then
         k=kmt(i,j)
         Fe_source(i,j,k)=Fe_source(i,j,k)+2.0/86400.0/1.0E+3/1.025/dzp(k)
         endif
       enddo
      enddo
!lyc-----------------------------------------------------------------------------------------
! for diagnostic
!-------------------------------------------------------------------------------------------
    if(ii==nss) then
	
      va0sum=0.0
      do k=1,km
        do j=2,jem
          do  i=2,imm
            if(vit(i,j,k)<0.5) cycle
            va0sum=va0sum+a0_b(i,j,k)*dxdyt(j)*dzp(k)
          enddo
        enddo
      enddo
       call mpi_barrier(mpi_comm_ocn,ierr)
       call mpi_reduce(va0sum,va0tmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)
!lyc--------------------------
      vb0sum=0.0
      do k=1,km
        do j=2,jem
          do  i=2,imm
            if(vit(i,j,k)<0.5) cycle
            vb0sum=vb0sum+b0_b(i,j,k)*dxdyt(j)*dzp(k)
          enddo
        enddo
      enddo
       call mpi_barrier(mpi_comm_ocn,ierr)
       call mpi_reduce(vb0sum,vb0tmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)

!lyc--------------------------
      vcsum=0.0
      do k=1,km
        do j=2,jem
          do  i=2,imm
            if(vit(i,j,k)<0.5) cycle
            vcsum=vcsum+c_b(i,j,k)*dxdyt(j)*dzp(k)
          enddo
        enddo
      enddo
       call mpi_barrier(mpi_comm_ocn,ierr)
       call mpi_reduce(vcsum,vctmp,1,mpi_real8,mpi_sum,0,mpi_comm_ocn,ierr)
       if(mytid==0) print *,month,iday, 'the total va0',va0tmp/vsea, 'the total vb0',vb0tmp/vsea, 'the total c_b',vctmp/vsea
	   
     endif	   
!lyc----------------------
! cumulate caco3 and a0
      do k=1,km
        do j=1,jmt
          do i=2,imm
          tocaco3(i,j,k)=tocaco3(i,j,k)+c_b(i,j,k)*dtts*vit(i,j,k)
          toa0(i,j,k)   =toa0   (i,j,k)+a0_b(i,j,k)*dtts*vit(i,j,k)
          enddo
        enddo
      enddo

#ifdef SPMD
      call exch_boundary(tocaco3(1,1,1),km)
      call exch_boundary(toa0(1,1,1),km)
#endif
          
#ifdef carbonDebug 
#ifdef SPMD
      print*, '--------------------------------------------------------------------------------------------'
      print*,'i,j,1,c,10,c',imt/2,jmt/2,c_b(imt/2,jmt/2,1),c_b(imt/2,jmt/2,10),'mytid=',mytid,'READYBIO'
#else
      print*, '--------------------------------------------------------------------------------------------'
      print*,'i,j,1,c,10,c',imt/2,jmt/2,c_b(imt/2,jmt/2,1),c_b(imt/2,jmt/2,10),'READYBIO'
#endif
#endif
#endif
!      
      deallocate(c2_b,inta1)

      RETURN
      END SUBROUTINE READYBIO


! CVS: $Id: sgec.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  SUBROUTINE SGEC
!========================
! SGEC
!---------------------------------------------------------------------
!
! purpose: calculate transport velocity of carbon
!
! author: Zhao Liang@lapc 2004/03/03 (original author: Xu Y F)
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod
      USE carbon_mod
      USE tracer_mod
      USE cforce_mod
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------------
      IMPLICIT NONE
#include <netcdf.inc>      
!#ifdef SPMD      
!#include <mpif.h>
!#endif      
!      
      REAL tssw(IMT,JMT),sssw(IMT,JMT)
      REAL factsg,secyr,rho00

      REAL tv(IMT,JMT),sch(IMT,JMT),arphs(IMT,JMT)
      REAL TSEA,SALT,dw,unitsg
!---------------------------------------------------------------------
!     factsg    =1.0E-4/3.1536E7
!     sge       =0.05 in units of mol/m^2/yr/ppm
!     factsg    for units of umol/dm^2/s (x10, for converting dz)
!     dz        in units of cm
!     sgbrok    in units of mol/m^2/yr/ppm
!---------------------------------------------------------------------
!     sge       in units of cm/s * umol/kg/ppm after factsg
!---------------------------------------------------------------------
!     factsg    =1.0E4/secyr*10.
!---------------------------------------------------------------------
!
      rho00=1025.0
      secyr=365.0*86400.0
      factsg=1.0E8/rho00/secyr
!        
!$OMP PARALLEL DO PRIVATE (I,J)        
      DO J=1,JMT
        DO I=1,IMT
          tssw(I,J)=AT(I,J,1,1)
          IF(AT(I,J,1,1).gt.40.0.and.AT(I,J,1,1).lt.999.) then
            tssw(I,J)=40.0
#ifdef SPMD
            PRINT*, 'mytid=',mytid,'notice: temperature greater than 40.0', &
                     I,J,AT(I,J,1,1),ITICE(I,J),"SGEC"                
#else                    
            PRINT*, 'notice: temperature greater than 40.0', &
                     I,J,AT(I,J,1,1),ITICE(I,J),"SGEC"                
#endif
          ENDIF    
          sssw(I,J)=AT(I,J,1,2)*1000.0+35.0
        ENDDO
      ENDDO
!---------------------------------------------------------------------
!     if(w22np(i,j) .le. 2.0) then
!       sge(i,j)=0.0
!     else
!       sge(i,j)=sgbrok(w22np(i,j))*factsg
!     end if
!     sgea(j)=0.05
!     sgep(j)=0.05
!---------------------------------------------------------------------
!     CALCULATION OF TRANSFER VELOCITY ACCORDING TO Wanninkhof's EQUATION
!     The values of parameters are refer to the references
!---------------------------------------------------------------------
      dw=1.027e3
      unitsg=24.0*365.E-8
!
!$OMP PARALLEL DO PRIVATE (I,J,TSEA,SALT) 
      DO J=1,JMT
        DO I=1,IMT
          TSEA=tssw(I,J)+273.15
          SALT=sssw(I,J)
          sge(I,J)=0.0
!---------------------------------------------------------------------
!     Schmit number at seawater
!---------------------------------------------------------------------
          IF(SALT.LE.0.1.OR.VIT(I,J,1).EQ.0) CYCLE
          sch(I,J)=2073.1-125.62*tssw(I,J)+3.6276*tssw(I,J)**2 &
                         -0.043219*tssw(I,J)**3
#ifdef SPMD
          IF(sch(I,J).LE.0.0) THEN 
            PRINT*, 'mytid=',mytid,'notice: SCH less than 0', &
                     i,j,sch(i,j),tssw(i,j),SALT,itice(i,j),"SGEC"
          ENDIF
#else          
          IF(sch(I,J).LE.0.0) THEN 
            PRINT*, 'notice: SCH less than 0', &
                     i,j,sch(i,j),tssw(i,j),SALT,itice(i,j),"SGEC"
          ENDIF
#endif
          arphs(I,J)=EXP(-60.2409+9345.17/TSEA            &
                         +23.3585*LOG(TSEA/100.)+SALT     &
                         *(0.023517-0.023656*(TSEA/100.)  &
                         +0.0047036*(TSEA/100.)**2))
!---------------------------------------------------------------------
!     TV from Wanninkhof's equation in units of cm/h
!     Sge is an exchange coefficient for CO2 in units of mol/m^2/yr/ppm
!---------------------------------------------------------------------
!     for steady or instant wind,7.946
!     for climatological wind, 10.0
!---------------------------------------------------------------------
!          tv(I,J)=7.946*w22np(I,J)*w22np(I,J)/SQRT(sch(I,J))
          tv(I,J)=10.0*w22np(I,J)*w22np(I,J)/SQRT(sch(I,J))
          sge(I,J)=tv(I,J)*arphs(I,J)*dw*unitsg
        ENDDO
      ENDDO
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J=1,JMT
        DO I=1,IMT
!---------------------------------------------------------------------
! devided by 100.0 to transform the unit from cm/s to m/s            
!---------------------------------------------------------------------
          sge(I,J)=sge(I,J)*factsg/100.0
        ENDDO
      ENDDO
      call exch_boundary(sge(1,1),1)
!      
      RETURN
      END SUBROUTINE SGEC


! CVS: $Id: isoflux_pt.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     ==============================
      SUBROUTINE ISOFLUX_PT (TF2,MTRACE)
!     ==============================
 
!     isopycnal diffusive tracer fluxes are computed.
use param_mod
use pconst_mod
use carbon_mod
use isopyc_mod
use work_mod
#ifdef SPMD
use msg_mod,only:mpi_comm_ocn
#endif
 
      IMPLICIT NONE
!#ifdef SPMD
!#include <mpif.h>
!#endif
 
      INTEGER :: mtrace
      REAL    :: p5,p25,c1,c0,fxa,fxb,fxc,fxe
      REAL,DIMENSION(imt,jmt,km)::TF2
 
      allocate (work_1(imt,jmt,km),work_2(imt,jmt,km),work_3(imt,jmt,0:km))
      allocate (temp(imt,jmt,km))
!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
      p5 = 0.5
 
      c0 = 0.0
      c1 = 1.0
      p25 = 0.25
 
      m = mtrace
 
!-----------------------------------------------------------------------
!     first compute the vertical tracer flux "temp" at the northern
!     face of "t" cells.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 2,km -1
         DO j = 2,jmm
            DO i = 1,imt
               temp (i,j,k)= p25* dzr (k)* (ptb (i,j +1,k -1,m) - ptb ( &
                  i,j +1,k +1,m) &
                   + ptb (i,j,k -1,m) - ptb (i,j, k +1,m))    
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     now consider the top level, assuming that the surface tracer
!     values are the same as the ones at "k"=1
!-----------------------------------------------------------------------
 
      k = 1
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 1,imt
            temp (i,j,k) = 0.25* dzr (k)* (ptb (i,j +1,k,m) - ptb (i,   &
                          j +1,k +1,m)+ptb (i,j,k,m) - ptb (i,j, k +1,m))
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     finally, consider the bottom level. the extrapolative estimator
!     is used to compute the tracer values at the ocean bottom.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 1,jmt
         DO i = 1,imt
            temp (i,j,km)= 0.0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j,i,k,fxa,fxb,fxc,fxe)
      DO j = 2,jmm
         DO i = 1,imt
            k = min (ITNU (i,j),ITNU (i,j +1))
            IF (k /= 0) THEN
               fxe = dzw (k -1) + dzw (k)
               fxa = 0.5* (ptb (i,j +1,k -1,m) + ptb (i,j,k -1,m))
               fxb = 0.5* (ptb (i,j +1,k,m) + ptb (i,j,k,m))
               fxc = dzwr (k -1)* (fxb * fxe- fxa * dzw (k))
               temp (i,j,k) = dzr (k)* (0.5* (fxa + fxb) - fxc)
            END IF
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     computation of meridional tracer flux
!     first calculate the effects of purely horizontal diffusion along
!     isopycnal, the background horizontal diffusion has been computed
!     before called this subroutine. (jxz)
!     add in the effects of the along isopycnal diffusion computed
!     using "K2" component of the tensor and apply land/sea masks
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO k = 1,km
         DO i = 1,imt
            work_1 (i,1,k)= 0.0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 1,imt
               work_1 (i,j,k)= ( ahisop * dyur (j)* (ptb (i,j +1,k,m)   &
                              - ptb (i,j,k,m)) + &
               ahisop * K2 (i,k,j,3)* temp (i,j,k) )* &
               SINU (j)* vit (i,j,k)* vit (i,j +1,k) 
            END DO
         END DO
      END DO
#ifdef SPMD
      call exch_boundary(work_1,km)
#endif
 
 
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "temp" at the eastern
!     face of "t" cells.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 2,km -1
         DO j = 2,jmm
            DO i = 1,imm
               temp (i,j,k)= p25* dzr (k)* (ptb (i +1,j,k -1,m) - ptb ( &
                             i+1,j,k+1,m)+ptb(i,j,k-1,m)-ptb(i,j,k+1,m))     
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     now consider the top level, assuming that the surface tracer
!     values are the same as the ones at "k"=1
!-----------------------------------------------------------------------
 
      k = 1
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 1,imm
            temp (i,j,k)= p25* dzr (k)* (ptb (i +1,j,k,m) - ptb (i +1,j,&
               k +1,m) + ptb (i,j,k,m) - ptb (i,j,k +1,m))  
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     finally, consider the bottom level. the extrapolative estimator
!     is used to compute the tracer values at the ocean bottom.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 1,imm
            temp (i,j,km) = 0.0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j,i,k,fxa,fxb,fxc,fxe)
      DO j = 2,jmm
         DO i = 1,imm
            k = min (ITNU (i,j),ITNU (i +1,j))
            IF (k /= 0) THEN
               fxe = dzw (k -1) + dzw (k)
               fxa = p5* (ptb (i,j,k -1,m) + ptb (i +1,j,k -1,m))
               fxb = p5* (ptb (i,j,k,m) + ptb (i +1,j,k,m))
               fxc = dzwr (k -1)* (fxb * fxe- fxa * dzw (k))
               temp (i,j,k) = dzr (k)* (p5* (fxa + fxb) - fxc)
            END IF
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     computtation of zonal tracer flux
!     first calculate the effects of purely horizontal diffusion along
!     isopycnal, the background horizontal diffusion has been computed
!     before called this subroutine. (jxz)
!     add in the effects of the along isopycnal diffusion computed
!     using "K1" component of the tensor and apply land/sea masks
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 1,imm
               work_2 (i,j,k)= ( ahisop * OTX (j)* (ptb (i +1,j,k,m)    &
                              - ptb (i,j,k,m)) + &
               ahisop*K1(i,k,j,3)*temp(i,j,k) )*vit(i+1,j,k)* vit(i,j,k) 
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "work_3" containing the K31
!     and K32 components which are to be solved explicitly. The K33
!     component will be treated semi-implicitly
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 2,km
         DO j = 2,jmm
            DO i = 2,imm
               work_3 (i,j,k -1) = ahisop * p25* vit (i,j,k)* ( &
               OTX (j)* K3 (i,k -1,j,1)* &
               (vit (i -1,j,k )* (ptb (i,j,k,m) - ptb (i -1,j,k,m)) &
               + vit (i -1,j,k -1)* (ptb (i,j,k -1,m) - ptb (i -1,j,k -1,m)) &
               + vit (i +1,j,k )* (ptb (i +1,j,k,m) - ptb (i,j,k,m)) &
               + vit (i +1,j,k -1)* (ptb (i +1,j,k -1,m) - ptb (i,j,    &
                                  k -1,m))) + &
               dytr (j)* K3 (i,k -1,j,2)* &
               (vit (i,j -1,k )* (ptb (i,j,k,m) - ptb (i,j -1,k,m)) &
               + vit (i,j -1,k -1)* (ptb (i,j,k -1,m) - ptb (i,j -1,k -1,m)) &
               + vit (i,j +1,k )* (ptb (i,j +1,k,m) - ptb (i,j,k,m)) &
                                   + vit (i,j +1,k -1)* (ptb (i,j +1,   &
                                    k -1,m) - ptb (i,j,k -1,m))) )  
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     at ocean surface the flux is set to zero to reflect the no tracer
!     flux condition. Same condition is also imposed at ocean bottom.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 2,imm
            work_3 (i,j, 0)= 0.0
            work_3 (i,j,km)= 0.0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 2,imm
               tf2 (i,j,k) = tf2 (i,j,k) &
               + cstrdytr (j)* (work_1 (i,j,k) - work_1 (i,j -1,k)) &
               + OTX (j) * (work_2 (i,j,k) - work_2 (i -1,j,k)) &
                            + dzr (k) * (work_3 (i,j,k -1) - work_3 (i,j,k))
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 2,imm
               work_1 (i,j,k) = adv_vntiso (i,k,j)* (ptb (i,j +1,k,m)   &
                                + ptb (i,j,k,m))
            END DO
         END DO
      END DO
 
#ifdef SPMD 
!      if (mytid==0) then
!!$OMP PARALLEL DO PRIVATE (k,i)
 !     DO k = 1,km
 !        DO i = 2,imm
 !           work_1 (i,1,k) = 0.0
 !        END DO
 !     END DO
 !     end if

!      call exchange_3d(work_1,km)
       call exch_boundary(work_1,km)
#else
!$OMP PARALLEL DO PRIVATE (k,i)
      DO k = 1,km
         DO i = 2,imm
            work_1 (i,1,k) = 0.0
         END DO
      END DO
#endif
 
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 1,imm
               work_2 (i,j,k) = adv_vetiso (i,k,j)* (ptb (i +1,j,k,m)   &
                                + ptb (i,j,k,m))
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     compute the vertical component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 2,imm
            work_3 (i,j, 0)= 0.0
            work_3 (i,j,km)= 0.0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 2,km
         DO j = 2,jmm
            DO i = 2,imm
               work_3 (i,j,k -1)= adv_vbtiso (i,k -1,j)* (ptb (i,j,k,m) &
                                + ptb (i,j,k -1,m))
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 2,imm
               tf2(i,j,k) = tf2(i,j,k) &
               - p5* cstrdytr (j)* (work_1 (i,j,k) - work_1 (i,j -1,k)) &
               - p5* OTX (j) * (work_2 (i,j,k) - work_2 (i -1,j,k)) &
                            - p5* dzr (k) * (work_3 (i,j,k -1)          &
                              - work_3 (i,j,k)) 
            END DO
         END DO
      END DO

      deallocate (work_1,work_2,work_3,temp)
 
      RETURN
      END SUBROUTINE ISOFLUX_PT
 
 
#else
      SUBROUTINE ISOFLUX_PT ()
      RETURN
      END SUBROUTINE ISOFLUX_PT
#endif 


! CVS: $Id: intfor_pt.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  SUBROUTINE INTFOR_PT
!========================
! INTFOR_PT
!---------------------------------------------------------------------
!
! purpose: interpolate the observed yely data to certain day
!
!lyc 2013,07
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod, ONLY: nmonth,mon0,iday,vit,j_global
      USE carbon_mod
      USE cforce_mod
      USE forc_mod
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!      
      INTEGER::nn
      REAL::days

      INTEGER::ipt1,ipt2
      REAL::factor
     
      if(mytid==0) print*,' entering intfor_pt'

      days=REAL(iday)
      DO i=1,mon0-1
        days=days+REAL(nmonth(i))
      ENDDO

      IF(iday<=15) THEN
        ipt1=mon0-1
        IF(ipt1==0) ipt1=12
        ipt2=mon0
        factor=REAL(iday-15)/REAL(nmonth(ipt1))+1
      ELSE
        ipt1=mon0
        ipt2=MOD(mon0,12)+1
        factor=REAL(iday-15)/REAL(nmonth(ipt1))
      ENDIF

!lyc 2014.07.01
       do j=1,jmt
        do i=1,imt
          w22np(i,j)=sqrt(wspdu3(i,j,1)**2+wspdv3(i,j,1)**2)
          pressureday(i,j)=psa3(i,j,1)/1.01325*1.0E-5
        enddo
       enddo
!
!!$OMP PARALLEL DO PRIVATE(j,i)
!      DO j=1,jmt
!        DO i=1,imt
!          w22np(i,j)=(winds(i,j,ipt2)-winds(i,j,ipt1))*factor+winds(i,j,ipt1)
!!cm090330----------------------------------------
!          pressureday(i,j)=(pressure(i,j,ipt2)-pressure(i,j,ipt1))*factor+pressure(i,j,ipt1)
!!cm090330----------------------------------------
!        ENDDO
!      ENDDO
!#endif

#ifdef carbonBio
#ifdef murnane1999
!$OMP PARALLEL DO PRIVATE(k,j,i)
      DO k=1,kmmix
        DO j=1,jmt
          DO i=1,imt
            po4force(i,j,k)=(po4obs(i,j,k,ipt2)-po4obs(i,j,k,ipt1))*factor+po4obs(i,j,k,ipt1)
          ENDDO
        ENDDO
      ENDDO
#endif
!
#ifdef progca
!$OMP PARALLEL DO PRIVATE(k,j,i)
      DO k=1,km
        taforce(k)=(taobs(k,ipt2)-taobs(k,ipt1))*factor+taobs(k,ipt1)
!
      ENDDO
#endif     
!---------------------------------------------------------------------
! for iron flux
! Whenever FE_FLUX_DATA AND DUST_DATA are true or false, fe_f is the dissolved iron at surface water (umol/kg*m/s)
! If FE_FLUX_DATA is true, dust_f is the iron  which is released from dust in the interior ocean! (umol/kg*m/s)
! If DUST_DATA is true, dust_f is the dust which sinks in the ocean (umol/kg*m/s)
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(j,i)
        DO j=1,jmt
          DO i=1,imt
           fe_f(i,j)=((fe_flux(i,j,ipt2)-fe_flux(i,j,ipt1))*factor+fe_flux(i,j,ipt1))*vit(i,j,1)
           dust_f(i,j)=((dust_flux(i,j,ipt2)-dust_flux(i,j,ipt1))*factor+dust_flux(i,j,ipt1))*vit(i,j,1)
           if(isnan(fe_f(i,j)).or.fe_f(i,j)>1.0) then
           print*, 'fe flux is error',fe_f(i,j),ipt1,fe_flux(i,j,ipt1),fe_flux(i,j,ipt2)
           stop
           endif
          ENDDO
        ENDDO
!       if(mytid==1) then
!        print *, 'fe_f,fe_flux',fe_f(91,jmt/2),fe_flux(91,jmt/2,ipt1)
!       endif
#endif

#if (defined carbonC) || (defined carbonBio) ||(defined carbonAbio)
#ifdef preindustrial
      pco2dry=280.776
#else
!cm080415
!      nn=INT(days/365.0*2.0)+2*(yearR+yearStart-yearData)-1
      nn=INT(days/365.0*2.0)+2*(yearR+yearStart-yearData-1)-1
      if(days<=182.5) then
          pco2dry=(csg(nn+1)-csg(nn))*days/365.0*2.0+csg(nn)
      else
         pco2dry=(csg(nn+1)-csg(nn))*(days/365.0*2.0-1)+csg(nn)
      endif
!cm080415
#endif
#endif
!
!lyc
#ifdef cfc
!cm080415
!    nn=yearR+(yearStart-yearData)
    nn=yearR+(yearStart-yearData-1)
!cm080415
     cfcatm(1)=(bomn(nn+1)-bomn(nn))*days/365.0+bomn(nn)
     cfcatm(2)=(boms(nn+1)-boms(nn))*days/365.0+boms(nn)
#endif 

#ifdef carbonC14
#ifdef preindustrial
     catm(1)= 100.0
     catm(2)= 100.0
     catm(3)= 100.0
#else
!cm080415
!      nn=yearR+(yearStart-yearData)
      nn=yearR+(yearStart-yearData-1)
!cm080415
!     interpolate the carbon concentrations to a given day
      catm(1)=(boml(nn+1)-boml(nn))*days/365.0+boml(nn)
      catm(2)=(bomm(nn+1)-bomm(nn))*days/365.0+bomm(nn)
      catm(3)=(bomh(nn+1)-bomh(nn))*days/365.0+bomh(nn)
#endif
#endif      
      if(mytid==0) print*,' intfor_pt is ok'

      Return

      END SUBROUTINE INTFOR_PT

! CVS: $Id: ssave_pt.F90,v 2.1 2004/06/10 07:45:18 cvsroot Exp $
  SUBROUTINE SSAVE_PT
!========================
! SSAVE_PT
!---------------------------------------------------------------------
!
! purpose: save, smooth and reset some array for CARBON model
!
! author: Zhao Liang@lapc 2004/03/04
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod
      USE carbon_mod
      USE cforce_mod
      USE coutput_mod
      USE cdf_mod
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
#ifdef COUP
      use buf_mod,only:pco2ups
#endif
!
!---------------------------------------------------------------------
      IMPLICIT NONE
#include <netcdf.inc>      
!#ifdef SPMD      
!#include <mpif.h>
!#endif      
!!      
      CHARACTER ( LEN =   4 ) :: ftail
      CHARACTER ( LEN =  15 ) :: fname
      CHARACTER ( LEN =  24 ) :: fname32
      LOGICAL :: hist_output,rest_output
      INTEGER :: nwmf,nday1,np,t0_cdf_pt
      REAL    :: overd
      integer lon_len,lat_len,lev_len,lev1_len,time_len
#if (!defined NORMAL)
!       integer,dimension(3)::start3,count3 
!       integer,dimension(4)::start4,count4
!       character::iret,ncid
        real(kind=4),dimension(imt_global,jmt_global)::cc2
        real(kind=4),dimension(imt_global,jmt_global,km)::cc3
        real(kind=4),parameter::spvalc=1.0E+35
        integer,dimension(3)::totup_dims,tpco2o_dims,tdpco2_dims
        integer,dimension(4)::cc_dims,po4_dims,ldoc_dims,ta_dims,prod_dims,fpop_dims,pldoc_dims,remi_dims,caco3_dims,o2_dims,fe_dims
        integer::totup_id,tpco2o_id,tdpco2_id,cc_id,po4_id,ldoc_id,ta_id,fe_id,prod_id,fpop_id,pldoc_id,remi_id,caco3_id,o2_id
        character(len=10)::dd,tt
        character(len=5)::zz
        integer*4::vv(8)
#endif 
 
 
!---------------------------------------------------------------------
!     output monthly results
!---------------------------------------------------------------------
      nday1  = nmonth(mon0)
      overd = 1.0/REAL(nday1)
      lon_len=imt_global
      lat_len=jmt_global
      lev_len=km
      lev1_len=km+1
      time_len=1

      DO k=1,klv
        DO j=1,jmt
          DO i=1,imt
            IF(VIT(i,j,k)>0.0) THEN
              CCMON(i,j,k)     = CCMON(i,j,k)*overd
#ifdef carbonBio
              po4mon(i,j,k)    = po4mon(i,j,k)*overd
              ldocmon(i,j,k)   = ldocmon(i,j,k)*overd
              tamon(i,j,k)     = tamon(i,j,k)*overd
              o2mon(i,j,k)     = o2mon(i,j,k)*overd
!lyc 2013,07
              femon(i,j,k)     = femon(i,j,k)*overd
              prodmon (i,j,k)  = prodmon(i,j,k)*overd
              fpopmon (i,j,k)  = fpopmon(i,j,k)*overd
              pldocmon (i,j,k) = pldocmon(i,j,k)*overd
              remimon (i,j,k)  = remimon(i,j,k)*overd
              jpopmon (i,j,k)  = jpopmon(i,j,k)*overd
              caco3mon (i,j,k) = caco3mon(i,j,k)*overd
#endif
            ELSE
              CCMON(i,j,k)     = 0.0
#ifdef carbonBio
              po4mon(i,j,k)    = 0.0
              ldocmon(i,j,k)   = 0.0
              tamon(i,j,k)     = 0.0
!lyc 2013,07
              o2mon(i,j,k)     = 0.0
              femon(i,j,k)     = 0.0
              prodmon (i,j,k)  = 0.0
              fpopmon (i,j,k)  = 0.0
              pldocmon (i,j,k) = 0.0
              remimon (i,j,k)  = 0.0
              jpopmon (i,j,k)  = 0.0
              caco3mon (i,j,k) = 0.0
#endif
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!---------------------------------------------------------------------
!     file name
#ifdef COUP
      nwmf= iyfm
#else
!lyc20121024      nwmf = yearStart + yearR - 1
      nwmf = yearStart + yearR
#endif
      WRITE (ftail,'(i4.4)') nwmf

      fname(1:5)='CCMON'
      fname(6:9)=ftail
      fname(10:12)=abmon(mon0)
 
      IF (nwmf >= yearStore .and. mod (nwmf,IO_out)==0 ) THEN
        hist_output=.true.
      ELSE
        hist_output=.false.
      ENDIF
!
      IF (mod ((monthR-1),io_rest)==0 ) THEN
        rest_output=.true.
      ELSE
        rest_output=.false.
      ENDIF
!
#ifdef SPMD
#ifdef printcall
#ifdef SPMD
      print*, "call local to global in ssave_pt, mytid=",mytid
#else
      print*, "call local to global in ssave_pt"
#endif
#endif
!      
      CALL local_to_global_pt(hist_output,rest_output)
!
#if (defined carbonC14)||(defined cfc)
      call local_to_global_4d_double(totup,totup_io,1,1)
#endif      
      IF (hist_output.AND.mytid==0) THEN
#else
      IF (hist_output) THEN
#endif 
!************************************************************************
#if (defined NORMAL) || (defined ALL)
 
         OPEN (87,FILE = FNAME(1:12),FORM ='UNFORMATTED',STATUS ='UNKNOWN')
 
         WRITE (87) nwmf,monthR -1
#if (defined carbonC14)||(defined cfc)
#ifdef SPMD
	  write(87) ((totup_io(i,j),i=1,imt_global),j=1,jmt_global)
#else
	  write(87) ((totup(i,j),i=1,imt),j=1,jmt)
#endif
#endif
#ifdef SPMD
         DO k = 1,km
           WRITE (87) ( (ccmon_io   (i,j,k),i = 1,imt_global),j = 1,jmt_global)
#ifdef carbonBio
           WRITE (87) ( (po4mon_io  (i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (ldocmon_io (i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (tamon_io   (i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (o2mon_io   (i,j,k),i = 1,imt_global),j = 1,jmt_global)
!lyc 2013,07
           WRITE (87) ( (femon_io   (i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (prodmon_io (i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (fpopmon_io (i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (pldocmon_io(i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (remimon_io (i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (jpopmon_io (i,j,k),i = 1,imt_global),j = 1,jmt_global)
           WRITE (87) ( (caco3mon_io(i,j,k),i = 1,imt_global),j = 1,jmt_global)
#endif
         ENDDO
#else
         DO k = 1,km
           WRITE (87) ( (ccmon   (i,j,k),i = 1,imt),j = 1,jmt)
#ifdef carbonBio
           WRITE (87) ( (po4mon  (i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (ldocmon (i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (tamon   (i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (o2mon   (i,j,k),i = 1,imt),j = 1,jmt)
!lyc 2013,07
           WRITE (87) ( (femon   (i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (prodmon (i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (fpopmon (i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (pldocmon(i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (remimon (i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (jpopmon (i,j,k),i = 1,imt),j = 1,jmt)
           WRITE (87) ( (caco3mon(i,j,k),i = 1,imt),j = 1,jmt)
#endif
         ENDDO
#endif
 
         CLOSE (87)
!
#if (defined carbonC) || (defined carbonBio) 
         FNAME='Ctime'//ftail//abmon(mon0)
         OPEN (97,FILE = FNAME(1:12),FORM ='UNFORMATTED',STATUS ='UNKNOWN')
!
#ifdef carbonC 
#ifdef SPMD
         WRITE (97) totup_io,tpco2o_io,tdpco2o_io,nwmf,monthR-1
#else
         WRITE (97) totup,tpco2o,tdpco2o,nwmf,monthR-1
#endif
#endif
!
#ifdef carbonBio 
#ifdef SPMD
         WRITE (97) totup_io,tpco2o_io,tdpco2o_io,tocaco3_io,toa0_io,nwmf,monthR-1
#else
         WRITE (97) totup,tpco2o,tdpco2o,tocao3,toa0,nwmf,monthR-1
#endif
#endif
         CLOSE (97)
#endif
!**********************************************************************
#else
#ifdef carbonBio
         fname(13:15)='.nc'
! output option (netcdf)
        ! enter define mode
         iret = nf_create (fname, NF_CLOBBER, ncid)
         CALL check_err (iret)
! define dimensions
         iret = nf_def_dim (ncid, 'lat', lat_len, lat_dim)
         CALL check_err (iret)
         iret = nf_def_dim (ncid, 'lon', lon_len, lon_dim)
         CALL check_err (iret)

      IF (mon0 == 12) THEN
          iret = nf_def_dim (ncid, 'lev', lev_len, lev_dim)
          CALL check_err (iret)
      ELSE
          iret = nf_def_dim (ncid, 'lev', klv, lev_dim)
          CALL check_err (iret)
      END IF
       iret = nf_def_dim (ncid, 'time', NF_UNLIMITED, time_dim)
!      iret = nf_def_dim(ncid, 'time', time_len, time_dim)
         CALL check_err (iret)

! define variables

        lat_dims (1) = lat_dim
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, lat_dims, lat_id)
         CALL check_err (iret)
!
         lon_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, lon_dims, lon_id)
         CALL check_err (iret)
!
         lev_dims (1) = lev_dim
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, lev_dims, lev_id)
         CALL check_err (iret)
         time_dims (1) = time_dim
         iret = nf_def_var (ncid, 'time', NF_DOUBLE, time_rank, time_dims, time_id)
         CALL check_err (iret)

         totup_dims (3) = time_dim
         totup_dims (2) = lat_dim
         totup_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'totup', NF_REAL,3 , totup_dims, totup_id)
         CALL check_err (iret)
         tpco2o_dims (3) = time_dim
         tpco2o_dims (2) = lat_dim
         tpco2o_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'tpco2o', NF_REAL,3 , tpco2o_dims, tpco2o_id)
         CALL check_err (iret)
         tdpco2_dims (3) = time_dim
         tdpco2_dims (2) = lat_dim
         tdpco2_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'tdpco2', NF_REAL,3 , tdpco2_dims, tdpco2_id)
         CALL check_err (iret)
         cc_dims (4) = time_dim
         cc_dims (3) = lev_dim
         cc_dims (2) = lat_dim
         cc_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'tc', NF_REAL,4 , cc_dims, cc_id)
         CALL check_err (iret)
         po4_dims (4) = time_dim
         po4_dims (3) = lev_dim
         po4_dims (2) = lat_dim
         po4_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'po4', NF_REAL,4 , po4_dims, po4_id)
         CALL check_err (iret)
         ldoc_dims (4) = time_dim
         ldoc_dims (3) = lev_dim
         ldoc_dims (2) = lat_dim
         ldoc_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ldoc', NF_REAL,4 , ldoc_dims, ldoc_id)
         CALL check_err (iret)
         ta_dims (4) = time_dim
         ta_dims (3) = lev_dim
         ta_dims (2) = lat_dim
         ta_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ta', NF_REAL,4 , ta_dims, ta_id)
         CALL check_err (iret)
         o2_dims (4) = time_dim
         o2_dims (3) = lev_dim
         o2_dims (2) = lat_dim
         o2_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'o2', NF_REAL,4 , o2_dims, o2_id)
         CALL check_err (iret)
         fe_dims (4) = time_dim
         fe_dims (3) = lev_dim
         fe_dims (2) = lat_dim
         fe_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'Fe', NF_REAL,4 , fe_dims, fe_id)
         CALL check_err (iret)
         prod_dims (4) = time_dim
         prod_dims (3) = lev_dim
         prod_dims (2) = lat_dim
         prod_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'prod', NF_REAL,4 , prod_dims, prod_id)
         CALL check_err (iret)
         fpop_dims (4) = time_dim
         fpop_dims (3) = lev_dim
         fpop_dims (2) = lat_dim
         fpop_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'fpop', NF_REAL,4 , fpop_dims, fpop_id)
         CALL check_err (iret)
         pldoc_dims (4) = time_dim
         pldoc_dims (3) = lev_dim
         pldoc_dims (2) = lat_dim
         pldoc_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'pldoc', NF_REAL,4 , pldoc_dims, pldoc_id)
         CALL check_err (iret)
         remi_dims (4) = time_dim
         remi_dims (3) = lev_dim
         remi_dims (2) = lat_dim
         remi_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'remi', NF_REAL,4 , remi_dims, remi_id)
         CALL check_err (iret)
!         jpop_dims (4) = time_dim
!         jpop_dims (3) = lev_dim
!         jpop_dims (2) = lat_dim
!         jpop_dims (1) = lon_dim
!         iret = nf_def_var (ncid, 'jpop', NF_REAL,4 , jpop_dims, jpop_id)
!         CALL check_err (iret)
         caco3_dims (4) = time_dim
         caco3_dims (3) = lev_dim
         caco3_dims (2) = lat_dim
         caco3_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'pcaco3', NF_REAL,4 , caco3_dims, caco3_id)
         CALL check_err (iret)
! assign attributes
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 18, 'depth (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, time_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, time_id, 'units', 23, 'months since 1001-01-01')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, totup_id, 'long_name', 12, 'total uptake')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, totup_id, 'units', 3, 'GtC')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, totup_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, tpco2o_id, 'long_name', 23, 'partail pressure of CO2')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, tpco2o_id, 'units', 4, 'uatm')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, tpco2o_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, tdpco2_id, 'long_name', 37, 'deviation of pco2 between air and sea')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, tdpco2_id, 'units', 4, 'uatm')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, tdpco2_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, cc_id, 'long_name', 9, 'total DIC')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, cc_id, 'units', 7, 'umol/kg')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, cc_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, po4_id, 'long_name', 3, 'PO4')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, po4_id, 'units', 7, 'umol/kg')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, po4_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, ldoc_id, 'long_name', 3, 'DOC')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ldoc_id, 'units', 7, 'umol/kg')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ldoc_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, ta_id, 'long_name', 16, 'total alkalinity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ta_id, 'units', 7, 'umol/kg')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ta_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, o2_id, 'long_name', 16, 'Dissolved oxygen')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, o2_id, 'units', 7, 'umol/kg')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, o2_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, fe_id, 'long_name', 14, 'Dissolved iron')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, fe_id, 'units', 7, 'umol/kg')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, fe_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, prod_id, 'long_name', 18, 'the new production')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, prod_id, 'units', 9, 'umol/kg/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, prod_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, fpop_id, 'long_name', 11, 'flux of POP')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, fpop_id, 'units', 9, 'umol/kg.m')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, fpop_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, pldoc_id, 'long_name', 17, 'production of DOC')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, pldoc_id, 'units', 11, 'umol/kg/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, pldoc_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, remi_id, 'long_name', 23, 'remineralization of DOC')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, remi_id, 'units', 9, 'umol/kg/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, remi_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, caco3_id, 'long_name', 40, 'production and remineralization of caco3')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, caco3_id, 'units', 9, 'umol/kg/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, caco3_id, 'missing_value', NF_REAL, 1, spvalc)
         CALL check_err (iret)
!   define global attribute
         CALL date_and_time (dd,tt,zz,vv)

         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', 4, 'test')
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', 20, tt //'  '//dd)
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', 14, 'LAPC/IAP OBGCM')
         CALL check_err (iret)
! leave define mode
         iret = nf_enddef (ncid)
         CALL check_err (iret)

!      t0_cdf=nday1
         t0_cdf_pt = monthR -1
!      t0_cdf=month-1+12000

         iret = nf_put_var_real (ncid, lon_id, real(lon,kind=4))
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, real(lat,kind=4))
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, real(lev,kind=4))
         CALL check_err (iret)
         start1 (1)= 1
         count1 (1)= time_len
         iret = nf_put_vara_double (ncid, time_id,start1,count1,t0_cdf_pt)
         CALL check_err (iret)
! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= time_len
         cc2=real(totup_io,kind=4)
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,1)<0.5) then
               cc2(i,j)=spvalc
               endif
               enddo
             enddo
         iret = nf_put_vara_real (ncid,totup_id,start3, count3, cc2)
         CALL check_err (iret)
         cc2=real(tpco2o_io,kind=4)
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,1)<0.5) then
               cc2(i,j)=spvalc
               endif
               enddo
             enddo
         iret = nf_put_vara_real (ncid,tpco2o_id,start3, count3, cc2)
         CALL check_err (iret)
         cc2=real(tdpco2o_io,kind=4)
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,1)<0.5) then
               cc2(i,j)=spvalc
               endif
               enddo
             enddo
         iret = nf_put_vara_real (ncid,tdpco2_id,start3, count3, cc2)
         CALL check_err (iret)
         cc3=real(ccmon_io,kind=4)
         start4(1)=1
         start4(2)=1
         start4(3)=1
         start4(4)=1
         count4(1)=lon_len
         count4(2)=lat_len
         count4(3)=lev_len
         count4(4)=time_len
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,cc_id,start4, count4, cc3)
         CALL check_err (iret)
         cc3=real(po4mon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,po4_id,start4, count4, cc3)
         CALL check_err (iret)
         cc3=real(ldocmon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,ldoc_id,start4, count4, cc3)
         CALL check_err (iret)
         cc3=real(tamon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,ta_id,start4, count4, cc3)
         CALL check_err (iret)
         cc3=real(o2mon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,o2_id,start4, count4, cc3)
         CALL check_err (iret)
!lyc 2013,07
         cc3=real(femon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,fe_id,start4, count4, cc3)
         CALL check_err (iret)

         cc3=real(prodmon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,prod_id,start4, count4, cc3)
         CALL check_err (iret)
         cc3=real(fpopmon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,fpop_id,start4, count4, cc3)
         CALL check_err (iret)
         cc3=real(pldocmon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,pldoc_id,start4, count4, cc3)
         CALL check_err (iret)
         cc3=real(remimon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,remi_id,start4, count4, cc3)
         CALL check_err (iret)
         cc3=real(caco3mon_io,kind=4)
         do k=1,lev_len
             do j=1,jmt_global
               do i=1,imt_global
               if(vit_global(i,j,k)<0.5) then
               cc3(i,j,k)=spvalc
               endif
               enddo
             enddo
         enddo
         iret = nf_put_vara_real (ncid,caco3_id,start4, count4, cc3)
         CALL check_err (iret)

          iret = nf_CLOSE (ncid)
         CALL check_err (iret)
#endif
!******************************************
#endif
 
      ENDIF
!--------------------------------------------------------------------- 
! set zore at the end of the month
!$OMP PARALLEL DO PRIVATE(k,j,i)
      DO k=1,klv
        DO j=1,jmt
          DO i=1,imt
            ccmon(i,j,k)     = 0.0
#ifdef carbonBio
            po4mon(i,j,k)    = 0.0
            ldocmon(i,j,k)   = 0.0
            tamon(i,j,k)     = 0.0
            o2mon(i,j,k)     = 0.0
            femon(i,j,k)     = 0.0
            prodmon (i,j,k)  = 0.0
            fpopmon (i,j,k)  = 0.0
            pldocmon (i,j,k) = 0.0
            remimon (i,j,k)  = 0.0
            jpopmon (i,j,k)  = 0.0
            caco3mon (i,j,k) = 0.0
#endif
          ENDDO
        ENDDO
      ENDDO

!--------------------------------------------------------------------- 
! set zore at the end of the year
      IF (MOD ( (monthR -1),12) == 0) THEN
!$OMP PARALLEL DO PRIVATE(k,j,i)
      DO k=1,km
        DO j=1,jmt
          DO i=1,imt
            ccmon(i,j,k)     = 0.0
#ifdef carbonBio
            po4mon(i,j,k)    = 0.0
            ldocmon(i,j,k)   = 0.0
            tamon(i,j,k)     = 0.0
            o2mon(i,j,k)     = 0.0
            femon(i,j,k)     = 0.0
            prodmon (i,j,k)  = 0.0
            fpopmon (i,j,k)  = 0.0
            pldocmon (i,j,k) = 0.0
            remimon (i,j,k)  = 0.0
            jpopmon (i,j,k)  = 0.0
            caco3mon (i,j,k) = 0.0
#endif
          ENDDO
        ENDDO
      ENDDO
      ENDIF
!

      IF (rest_output) THEN
#ifdef SPMD
#if (defined carbonC14)||(defined cfc)
      call local_to_global_4d_double(totup,totup_io,1,1)
#endif      
      IF(mytid==0) THEN
        fname32(1:8)='fort.32.'
        fname32(9:12)=ftail
#ifdef COUP
        fname32(13:13)='-'
        if(mon0<12) then
        write(fname32(14:15),'(i2.2)')mon0+1
        else
        write(fname32(14:15),'(i2.2)')mon0-11
        write(fname32(9:12),'(i4.4)')nwmf+1
        endif
        fname32(16:24)='-01-00000'
#else
        fname32(13:15)=abmon(mon0)
#endif
        OPEN (90,FILE = fname32,FORM ='unformatted',STATUS ='unknown')
#if (defined carbonC14)||(defined cfc)        
        WRITE (90) pt_io,totup_io,monthR
#endif        
!
#ifdef carbonC        
        WRITE (90) pt_io,totup_io,tpco2o_io,tdpco2o_io,monthR
#endif        
!
#ifdef carbonBio        
#ifdef COUP
        WRITE (90) pt_io,totup_io,tpco2o_io,tdpco2o_io,tocaco3_io,toa0_io,pco2ups,monthR
#else
        WRITE (90) pt_io,totup_io,tpco2o_io,tdpco2o_io,tocaco3_io,toa0_io,monthR
#endif
#endif        
        CLOSE (90)
      ENDIF
#else
        fname32(1:8)='fort.32.'
        fname32(9:12)=ftail
        fname32(13:13)='-'
#ifdef COUP
        if(mon0<12) then
        write(fname32(14:15),'(i2.2)')mon0+1
        else
        write(fname32(14:15),'(i2.2)')mon0-11
        write(fname32(9:12),'(i4.4)')nwmf+1
        endif
        fname32(16:24)='-01-00000'
#else  
        fname32(13:15)=abmon(mon0)
#endif
        OPEN (90,FILE = fname32,FORM ='unformatted',STATUS ='unknown')
#if (defined carbonC14)||(defined cfc)        
        WRITE (90) pt,totup,monthR
#endif        
!
#ifdef carbonC        
        WRITE (90) pt,totup,tpco2o,tdpco2o,monthR
#endif        
!
#ifdef carbonBio        
#ifdef COUP
        WRITE (90) pt,totup,tpco2o,tdpco2o,tocaco3,toa0,pco2ups,monthR
#else
        WRITE (90) pt,totup,tpco2o,tdpco2o,tocaco3,toa0,monthR
#endif
      
#endif        
        CLOSE (90)
#endif
      ENDIF
 
#ifdef SPMD
#if (defined carbonC14)||(defined cfc)
      call local_to_global_4d_double(totup,totup_io,1,1)
#endif      
      IF(mytid==0) THEN
        OPEN(32,FILE='fort.32',FORM='unformatted')
!
#if (defined carbonC14)||(defined cfc)        
        WRITE (32) pt_io,totup_io,monthR
#endif        
#ifdef carbonC        
        WRITE (32) pt_io,totup_io,tpco2o_io,tdpco2o_io,monthR
#endif        
!
#ifdef carbonBio        
#ifdef COUP
        WRITE (32) pt_io,totup_io,tpco2o_io,tdpco2o_io,tocaco3_io,toa0_io,pco2ups,monthR
#else
        WRITE (32) pt_io,totup_io,tpco2o_io,tdpco2o_io,tocaco3_io,toa0_io,monthR
#endif
#endif        
        CLOSE (32)
      ENDIF
#else
      REWIND 32
#if (defined carbonC14)||(defined cfc)        
        WRITE (32) pt,totup,monthR
#endif        
!
#ifdef carbonC        
        WRITE (32) pt,totup,tpco2o,tdpco2o,monthR
#endif        
!
#ifdef carbonBio        
        WRITE (32) pt_io,totup,tpco2o,tdpco2o,tocaco3,toa0,monthR
#endif        
      CLOSE(32)
#endif

   
!---------------------------------------------------------------------
!     reset some arrays
!---------------------------------------------------------------------
      DO m = 1,nptra
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
        DO j = 1,jmt
          DO i = 1,imt
            ptb (i,j,k,m) = pt (i,j,k,m) * vit(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      ENDDO

!lyc 2014.06
      MONTHR=MONTHR+1
      
      RETURN
      END SUBROUTINE SSAVE_PT

!  CVS: $Id: upwell.F90,v 1.5 2003/08/12 09:06:39 lhl Exp $
!     ===============================
      SUBROUTINE UPWELL_PT (UWK,VWK,UTL2,WKD2,WKB2,wst)
!     ===============================
 
#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use carbon_mod
 
      IMPLICIT NONE
      REAL    :: UWK (IMT,JMT,KM),VWK (IMT,JMT,KM),H0WK (IMT,JMT)
      real    :: hb_x1,hb_x2,hb_x,hb_y
      REAL,DIMENSION(imt,jmt,km)::wkb2,wkd2,utl2
      REAL,DIMENSION(imt,jmt,kmp1)::wst
!      REAL,DIMENSION(imt,jmt_global,kmp1)::wst_io
      integer :: a
      
!      allocate(UTL2 (IMT,JMT,KM),WKD2 (IMT,JMT,KM),WKB2 (IMT,JMT,KM))
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
      allocate(uk(imt,jmt,km),vk(imt,jmt,km))
      DO J = JST,JET
         DO I = 1,IMT
            WORK (I,J) = 0.0
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K)= 0.0
            END DO
         END DO
      END DO

   wst(:,:,1:kmp1)=0.0
 
 
!!$OMP PARALLEL DO PRIVATE (J,I)
!      DO J = JSM,JEM
!         DO I = 2,IMT
!            WORK (I,J)= 0.25* (H0WK (I,J) + H0WK (I -1,J) + H0WK (I,    &
!                        J +1) + H0WK (I -1,J +1))
!         END DO
!         WORK (1,J)= WORK (IMM,J)
!      END DO

! the topography of bottom
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JSM,JEM
         DO I = 2,IMT
           if(ohbt(i,j)>0.0) then
            WORK (I,J)= 1/OHBT(I,J)
           else
            WORK (I,J)= 0.0
           endif
         END DO
         WORK (1,J)= WORK (IMM,J)
      END DO
#ifdef SPMD
      call exchange_2d(work,1,1)
#endif
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UK (I,J,K)=  UWK (I,J,K)
               VK (I,J,K)=  VWK (I,J,K)
            END DO
         END DO
      END DO
! -----------------------------
!    the partial(u)/partial(x) + partial(v)/partial(y)
!------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,K)= 0.5* OTX (J)* ( (UK (I +1,J,K) + UK (I +1,  &
                           J -1,K)) &
               - (UK (I,J,K) + UK (I,J -1,K))) &
               + R2A (J)*2.0* (VK (I,J,K) + VK (I +1,J,K)) &
               - R2B (J)*2.0* (VK (I,J -1,K) + VK (I +1,J -1,K))      
            END DO
         END DO
      END DO
    a =0
    IF (a==1) then
!calculate the ws of the last level 
         DO J=JSM,JEM
           DO I=2,IMM
!            if(vit(i,j,km)<0.5) cycle
!            HB_x1=(WORK (I  ,J) - WORK (I-1,J))* UTL2 (I,J,KM)
!            HB_x2=(WORK (I+1,J) - WORK (I  ,J))* UTL2 (I+1,J,KM)
!            HB_x = HB_x2+HB_x1
!            HB_y=WKD2(I,J,KM)*(WORK(I,J+1)-WORK(I,J))+WKB2(I,J,KM)*(WORK(I,J)-WORK(I,J-1))
!            WSt(I,J,KMP1)=(HB_X+HB_Y)*vit(i,j,km)
             WST(I,J,KMP1)=0.0
           ENDDO
         ENDDO
!!$OMP PARALLEL DO PRIVATE (J,I,HB_x,HB_x1,HB_x2,HB_y)          
         DO J=JSM,JEM
           DO I=2,IMM
           do k=km,1,-1
           if(vit(i,j,k)<0.5) cycle
! calculate ws  in the bottom           
           IF(k==ITNU(i,j)) THEN
            WST(I,J,K+1)=0.0
            HB_x1=(WORK (I  ,J) - WORK (I-1,J))* UTL2 (I,J,K)
            HB_x2=(WORK (I+1,J) - WORK (I  ,J))* UTL2 (I+1,J,K)
            HB_x = HB_x2+HB_x1
            HB_y=WKD2(I,J,K)*(WORK(I,J+1)-WORK(I,J))+WKB2(I,J,K)*(WORK(I,J)-WORK(I,J-1))
            WSt(I,J,K)=-(HB_X+HB_Y)*vit(i,j,k)
           ELSE  
 ! if vit(i,j,k) and vit(i,j,k+1) ==1
              
            WST(I,J,K)= VIT (I,J,K)* (WST (I,J,K+1 )- &
                             DZP(K)* WKA(I,J,K))
          END IF
          
          enddo
! calculate the ws of the first level          
              WST(I,J,1)=vit(i,j,1)*(WST(I,J,2)-DZP(1)*WKA(I,J,1))    
         END DO
      END DO
 
  ELSE
!      wst(:,:,1)=0.0
      DO K=2,km
        DO J=JSM,JEM
           DO I=2,IMM
           WST(I,J,K)=VIT (I,J,K)* (WST (I,J,K-1 )+ &
                             DZP(K-1)* WKA(I,J,K-1))
           END DO
        END DO
      END DO
  ENDIF      
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
     DO K = 1,KM
        DO J = JSM,JEM
           DO I = 2,IMM
              WST (I,J,K)= WST (I,J,K)* VIT (I,J,K)
!              if(wst(i,j,k)>0.001) then 
!            print *,'wst ',i,j,'k',k,wst(i,j,k),'mytid=',mytid
!              endif
           END DO
        END DO
!            print *,'wka','wka',k,wka(imt/2,jem/2,k),'mytid=',mytid
     END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      if(nx_proc==1) then
      DO K = 1,KM+1
         DO J = JSM,JEM
            WST (1,J,K) = WST (IMM,J,K)
            WST (IMT,J,K) = WST (2,J,K)
         END DO
      END DO
      endif

#ifdef SPMD
      call exch_boundary(wst,kmp1)
#endif
!      call local_to_global_4d_double(wst,wst_io,kmp1,1)
!      if(mytid==0) then
!          do k=1,km
!            do j=1,jmt_global
!              do i=1,imt
!         wstmon(i,j,k)=wstmon(i,j,k)+wst_io(i,j,k)/nss
!              enddo
!            enddo
!          enddo
!       endif
!          open (99,file='ws0.dat',form='unformatted')
!          write(99) wst_io(:,:,1:km)
!          close(99)
      deallocate(uk,vk)

      RETURN
      END SUBROUTINE UPWELL_PT


! CVS: $Id: accumm_pt.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
!========================
  SUBROUTINE ACCUMM_PT
!========================
! ACCUMM_PT
!---------------------------------------------------------------------
!
! purpose: sum up some variables used in CARBON model
!
! author: Zhao Liang@lapc 2004/03/04
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE carbon_mod
      USE coutput_mod
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!#include <netcdf.inc>      
!#ifdef SPMD      
!#include <mpif.h>
!#endif      
!      
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               CCMON (I,J,K)= CCMON (I,J,K) + PT (I,J,K,1)
#ifdef carbonC14
#endif
!
#ifdef carbonC
#endif
!
#ifdef carbonBio
               po4mon (I,J,K)   = po4mon (I,J,K) + PT (I,J,K,2)
               ldocmon (I,J,K)  = ldocmon (I,J,K) + PT (I,J,K,3)
               tamon (I,J,K)    = tamon (I,J,K) + PT (I,J,K,4)
               prodmon (I,J,K)  = prodmon (I,J,K) + a1_b (I,J,K)
               fpopmon (I,J,K)  = fpopmon (I,J,K) + a2_b (I,J,K)
               pldocmon (I,J,K) = pldocmon (I,J,K) + b1_b (I,J,K)
               remimon (I,J,K)  = remimon (I,J,K) + b2_b (I,J,K)
               jpopmon (I,J,K)  = jpopmon (I,J,K) + a0_b (I,J,K)
               caco3mon (I,J,K) = caco3mon (I,J,K) + c_b (I,J,K)
               o2mon (I,J,K)   = o2mon (I,J,K) + PT (I,J,K,5)
               femon (I,J,K)   = femon (I,J,K) + PT (I,J,K,6)
              
#endif
            END DO
         END DO
      END DO
!lyc
#if (defined cfc)||(defined carbonC14)
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = 1,JMT
           DO I= 1,IMT
           ssfcmon(i,j)=ssfcmon(i,j)+ssfc(i,j)
           ENDDO
         ENDDO
#endif         
      RETURN
      END SUBROUTINE ACCUMM_PT
 
! CVS: $Id: local_to_global_pt.F90,v 2.2 2004/06/13 12:14:56 cvsroot Exp $
!----------------------------------------------------------
  SUBROUTINE LOCAL_TO_GLOBAL_PT(hist_output,rest_output)
!========================
! LOCAL_TO_GLOBAL_PT
!---------------------------------------------------------------------
!
! purpose: collect variables from local to global in CARBON model
!
! author: Zhao Liang@lapc 2004/03/04; lyc 2012.10.10
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE pconst_mod
      USE work_mod
      USE carbon_mod
      USE cforce_mod
      USE coutput_mod
#ifdef SPMD      
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------
      IMPLICIT NONE
!#include <netcdf.inc>      
!      
!111-------------------
#ifdef SPMD
      INTEGER :: itag
      LOGICAL :: hist_output,rest_output 
!
!      INTEGER,PARAMETER :: nptram1=nptra-1 
! number of variables for output

!restart data
      IF (rest_output) THEN
! lyc------------------------------
! add the collection for c_b and a0

#ifdef carbonBio            
      call local_to_global_4d_double(tocaco3,tocaco3_io,km,1)
      call local_to_global_4d_double(toa0,toa0_io,km,1)
#endif         
! lyc ----------------------------------------------------------------
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)          
      call local_to_global_4d_double(totup,totup_io,1,1)
      call local_to_global_4d_double(tpco2o,tpco2o_io,1,1)
      call local_to_global_4d_double(tdpco2o,tdpco2o_io,1,1)
#endif
      call local_to_global_4d_double(pt,pt_io,km,nptra)
      endif


!-------------------------------------------------------------------
! output data     
!------------------------------ 
      IF (hist_output) THEN
! lyc------------------------------
! add the collection for c_b and a0

#ifdef carbonBio            
     call local_to_global_4d_double(tocaco3,tocaco3_io,km,1)
     call local_to_global_4d_double(toa0,toa0_io,km,1)
#endif         
! lyc ----------------------------
#if (defined carbonC) || (defined carbonBio)|| (defined carbonAbio) 
      call local_to_global_4d_double(totup,totup_io,1,1)
      call local_to_global_4d_double(tpco2o,tpco2o_io,1,1)
      call local_to_global_4d_double(tdpco2o,tdpco2o_io,1,1)
#endif
     call local_to_global_4d_double(ccmon,ccmon_io,km,1)

#ifdef carbonBio            
      call local_to_global_4d_double(po4mon,po4mon_io,km,1)
      call local_to_global_4d_double(ldocmon,ldocmon_io,km,1)
      call local_to_global_4d_double(tamon,tamon_io,km,1)
      call local_to_global_4d_double(o2mon,o2mon_io,km,1)
      call local_to_global_4d_double(femon,femon_io,km,1)
      call local_to_global_4d_double(prodmon,prodmon_io,km,1)
      call local_to_global_4d_double(fpopmon,fpopmon_io,km,1)
      call local_to_global_4d_double(pldocmon,pldocmon_io,km,1)
      call local_to_global_4d_double(remimon,remimon_io,km,1)
      call local_to_global_4d_double(jpopmon,jpopmon_io,km,1)
      call local_to_global_4d_double(caco3mon,caco3mon_io,km,1)
#endif         
      endif 
!-------------------------------------------------
!111***************
#endif
!-------------------
      RETURN
      END SUBROUTINE LOCAL_TO_GLOBAL_PT
 
! CVS: $Id: deallocate_pt.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  SUBROUTINE DEALLOCATE_PT
!========================
! DEALLOCATE_PT
!---------------------------------------------------------------------
!
! purpose: deallocate temporary variables used in CARBON model
!
! author: Zhao Liang@lapc 2004/03/02
!
!---------------------------------------------------------------------
#include <def-undef.h>       
!
      USE param_mod
      USE cforce_mod
#ifdef SPMD
      USE msg_mod,only:mpi_comm_ocn
#endif      
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!      
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      WRITE(6,*) 'Begining------DEALLOCATE_PT'
#ifdef SPMD
      ENDIF
#endif
      
#ifdef carbonC14
      DEALLOCATE(kyear,boml,bomm,bomh)
#endif
!      
#if (defined carbonC)|| (defined carbonBio)
      DEALLOCATE(csgn,csg,csgx)
#endif      
!  
#ifdef carbonBio

#endif
!
#ifdef SPMD
      IF(mytid==0) THEN
#endif      
      WRITE(6,*) 'END-----------DEALLOCATE_PT'
#ifdef SPMD
      ENDIF
#endif

  END SUBROUTINE DEALLOCATE_PT
