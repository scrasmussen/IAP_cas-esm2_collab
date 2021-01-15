SUBROUTINE NLITTI
!------------------------------------------------------------------------------------------------
! Purpose: NONLINEAR ITERATIVE TIME INTEGRATION & splitting method
! Original version : NLITTI.f (IAP 21L)
! Reconstructed & add splitting method in : ZhangHe
! Completed : 2006.2.20
! Update    : 2006.12.10, Zhanghe, delete CALL pstartend2 when integration for advection term  
!             2007.04.23, ZhangHe, 1) delete dumb parameter 'SETUV' & 'NCTDCB'
!                                  2) 'ptuvtend1' ==> 'tend_lin'; 'ptuvtend2' ==> 'tend_adv' 
!             2007.05.07, ZhangHe, 1) 'pstartend' ==> 'tend_pstar' ;
!                                  2) 'nliter2'   ==> 'nliter_uvtp' ;
!                                  3) 'nliter3'   ==> 'nliter_uvt'
!             2007.12.22, ZhangHe, move the calling of QPDATA from lin timestep to adv timestep
!             2008.04.10, ZhangHe, add calculation of UVW0, update calling nliter_uvtp
!                                  move statement of Istar to module flexib
!             2008.04.23, ZhangHe, delete calling QPDATA & calculation of UVW0
!             2008.06.05, WuJianping, for parallel version
!             2010.08, Juanxiong He
!             October 2012, Jiang Jinrong, for 2D parallel
! Reviewed: ZhangHe, 2011-11-19
!           ZhangHe, 2012-10-31
!           ZhangHe, 2013-03-27
! Modified: ZhangHe, 2013-03-29, use dynamic arrays 
! Modified: ZhangMinghua, 2018-08-25, use tend_lin2 and tend_adv2 for time integration
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only : NX, NY, NL
   use IAP_prog,  only : Psa, UT, VT, TT, Pstar1, T, U ,V,PT,P,Pstar2
   use flexib,    only : DTlin, DTadv, Ndt, Istar, zmh_nlitti1, zmh_nlitti2,wbd_nlitti1, wbd_nlitti2  !zmh  && wbd
   use mathconst, only : HALF
   use pmgrid,    only : beglatdyn, endlatdyn, beglev, endlev     !zhh 2012.10.31
   use spmd_utils, only: masterproc, iam ! juanxiong he, 2010.08
   use Trans_coef, only: trans_antiIAP, trans_IAP
   use smoother,   only: SMOTHP, SMOOTH
   use perf_mod, only : t_startf, t_stopf ! zhh  2013-02-06
   use tendency,  only: DU, DV, DT ,DPSA  !debug
use prognostics,      only:take_valuet !jjr               
   implicit none
!----------------------------------Local workspace-----------------------------------------------
   real(r8), allocatable :: PsaL(:,:)
   real(r8), allocatable :: PstL(:,:)
   real(r8), allocatable :: UL(:,:,:)
   real(r8), allocatable :: VL(:,:,:)
   real(r8), allocatable :: TL(:,:,:)
   integer  :: I, J, K, n, ITERNM    ! loop index

!------------------------------------------------------------------------------------------------
! allacate temporary arrays
   allocate ( PsaL(NX,beglatdyn:endlatdyn) )
   allocate ( PstL(NX,beglatdyn:endlatdyn) )
   allocate ( UL(NX,beglev:endlev,beglatdyn:endlatdyn) )
   allocate ( VL(NX,beglev:endlev,beglatdyn:endlatdyn) )
   allocate ( TL(NX,beglev:endlev,beglatdyn:endlatdyn) )

![1]  integration for linear(adaption) term  
   do n = 1, Ndt
      DO J = beglatdyn, endlatdyn
         DO I = 1, NX
            PsaL(I,J)    = Psa(I,J)
            PstL(I,J)    = Pstar1(I,J)
         END DO
      END DO
      DO J = beglatdyn, endlatdyn
         DO K = beglev,endlev 
            DO I = 1, NX
               UL(I,K,J) = UT(I,K,J)
               VL(I,K,J) = VT(I,K,J)
               TL(I,K,J) = TT(I,K,J)
            END DO
         END DO
      END DO
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, NX
!if(i==1.and.k==1.and.j==1) print*,'test jjr 1 in nlitt' &
!       ,uL(i,k,j),vL(i,k,j),TL(i,k,j),PsaL(i,j),PstL(i,j),DT(i,k,j),DPSA(i,j)
!        ,PT(I,J),TT(I,K,J),P(i,j),&
!       pstar2(i,j),pstar1(i,j),ut(i,k,j),vt(i,k,j),psa(i,j)


        !if(i==60.and.k==1.and.j==127) print*,'test jjr 1 in nlitt',&
        ! U(I,K,J+1),Ut(I,K,J+1)
!       end do
!       end do
!   end do

!
! ============================================

!zmh 2018-09-12  removed the backward integration, added embedded substepping that can replace Ndt
!wbd 2018-11-16  added the leaping format in tend_lin, tend_adv, tend_lin2, tend_adv2

 if(zmh_nlitti1) then
      CALL tend_pstar( Istar )
      if(wbd_nlitti1) then
         CALL tend_lin2_leap       ! 2007.4.23
      else
         CALL tend_lin2       ! 2007.4.23
      endif
  else
! ============================================
!-----------------------------------------------
!-----------------------------------------------
!  DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, NX
!        if(i==170.and.k==8.and.j==8) print*,'test jjr in nli0',&
!         PT(i,j),Pstar1(i,j),DPSA(i,j)
!       end do
!       end do
!   end do


      CALL tend_pstar( Istar )
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, NX
!        if(i==170.and.k==8.and.j==8) print*,'test jjr in nli1',&
!         PT(i,j),Pstar1(i,j),DPSA(I,J)
!if(i==1.and.k==1.and.j==1) print*,'test jjr 2 in nlitt' &
!,uL(i,k,j),vL(i,k,j),TL(i,k,j),PsaL(i,j),PstL(i,j),DT(i,k,j),DPSA(i,j)
!       ,uL(4,k,j),vL(4,k,j),TL(4,k,j),PsaL(4,j),PstL(4,j),DT(4,k,j),DPSA(4,j)

!       end do
!      end do
!   end do


!-----------------------------------------------
      if(wbd_nlitti1) then
         CALL tend_lin_leap       ! 2007.4.23
      else
         CALL tend_lin       ! 2007.4.23
      endif
!-----------------------------------------------
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, NX
!        if(i==1.and.k==1.and.j==1) print*,'test jjr 3 in nlitt',&
!uL(i,k,j),vL(i,k,j),TL(i,k,j),PsaL(i,j),PstL(i,j),DT(i,k,j),DPSA(i,j)
!         V(I,K,J),U(I,K,J),PT(i,j),Pstar1(i,j),DPSA(I,J)
!       end do
!       end do
!   end do
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, nx
!if(i==60.and.k==1.and.j==128) print*,'test jjr 3 in nlitt' &
     !   ,PT(I,J),TT(I,K,J),P(i,j),&
     !  pstar2(i,j),pstar1(i,j),ut(i,k,j),vt(i,k,j),psa(i,j)
!,uL(4,k,j),vL(4,k,j),TL(4,k,j),PsaL(4,j),PstL(4,j),DT(4,k,j),DPSA(4,j)

       ! if(i==60.and.k==1.and.j==127) print*,'test jjr 3 in nlitt',&
       !  U(I,K,J+1),Ut(I,K,J+1)
!       end do
!       end do
!   end do

!      call take_valuet

      CALL nliter_uvtp(UL, VL, TL, PsaL, PstL, DTlin, Istar)
! ============================================
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, NX
!        if(i==170.and.k==8.and.j==8) print*,'test jjr in nli2',&
!         V(I,K,J),U(I,K,J),PT(I,J),Pstar1(i,j)
!       end do
!       end do
!   end do
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, nx
!if(i==1.and.k==1.and.j==1) print*,'test jjr 4 in nlitt' &
!,uL(i,k,j),vL(i,k,j),TL(i,k,j),PsaL(i,j),PstL(i,j),DT(i,k,j),DPSA(i,j)
!        ,PT(I,J),TT(I,K,J),P(i,j),&
!       pstar2(i,j),pstar1(i,j),ut(i,k,j),vt(i,k,j),psa(i,j)
!
      !  if(i==60.and.k==1.and.j==127) print*,'test jjr 4 in nlitt',&
      !   U(I,K,J+1),Ut(I,K,J+1)
!       end do
!       end do
!   end do


!
      DO ITERNM = 1,1
!------------------------------------------------------------------------------------------------
         CALL tend_pstar( Istar )
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, NX
!        if(i==170.and.k==8.and.j==8) print*,'test jjr in nli3',&
!         V(I,K,J),U(I,K,J),UT(i,k,j),vt(i,k,j)
!if(i==60.and.k==1.and.j==128) print*,'test jjr 5 in nlitt' &
!        ,PT(I,J),TT(I,K,J),P(i,j),&
!       pstar2(i,j),pstar1(i,j),ut(i,k,j),vt(i,k,j),psa(i,j)

!       end do
!       end do
!   end do

! ============================================      
         if(wbd_nlitti1) then
            CALL tend_lin_leap       ! 2007.4.23
         else
            CALL tend_lin       ! 2007.4.23
         endif
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, NX
!        if(i==170.and.k==8.and.j==8) print*,'test jjr in nli4',&
!         V(I,K,J),U(I,K,J),UT(i,k,j),vt(i,k,j)
!       end do
!       end do
!   end do
!jjr test   stop
!   DO J = beglatdyn, endlatdyn
!          DO K = beglev, endlev
!            DO I = 1, nx
!if(i==60.and.k==1.and.j==128) print*,'test jjr 6 in nlitt' &
!        ,PT(I,J),TT(I,K,J),P(i,j),&
!       pstar2(i,j),pstar1(i,j),ut(i,k,j),vt(i,k,j),psa(i,j)


!       end do
!       end do
!   end do

! ======================================================================
!      call take_valuet
         CALL nliter_uvtp(UL, VL, TL, PsaL, PstL, DTlin, Istar)
! ======================================================================
!------------------------------------------------------------------------------------------------

         DO J = beglatdyn, endlatdyn
            DO I = 1, NX
               Psa(I,J)  = (Psa(I,J) + PsaL(I,J)) * HALF
               Pstar1(I,J)  = (Pstar1(I,J) + PstL(I,J)) * HALF     
            END DO
         END DO
         DO J = beglatdyn, endlatdyn
            DO K = beglev,endlev 
               DO I = 1, NX
                  UT(I,K,J) = (UT(I,K,J) + UL(I,K,J)) * HALF
                  VT(I,K,J) = (VT(I,K,J) + VL(I,K,J)) * HALF
                  TT(I,K,J) = (TT(I,K,J) + TL(I,K,J)) * HALF
               END DO
            END DO
         END DO
 

! ============================================
         CALL tend_pstar( Istar )
! ============================================           
         if(wbd_nlitti1) then
            CALL tend_lin_leap       ! 2007.4.23
         else
            CALL tend_lin       ! 2007.4.23
         endif! ====================================================================
!           call mpi_barrier(mpicom,ierr)
!             stop !jjr test
!      call take_valuet
         CALL nliter_uvtp(UL, VL, TL, PsaL, PstL, DTlin, Istar)
! ====================================================================
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
      END DO

   endif  !zmh end

   end do    ! for n = 1, Nt


!********************************************************************************
![2]  integration for advection term
   do n = 1, 1

!zmh
 if(zmh_nlitti2)then
      if(wbd_nlitti2) then
         CALL tend_adv2_leap       ! 2007.4.23
      else
         CALL tend_adv2       ! 2007.4.23
      endif
 else
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

      DO J = beglatdyn, endlatdyn
         DO K = beglev,endlev 
            DO I = 1, NX
               UL(I,K,J) = UT(I,K,J)
               VL(I,K,J) = VT(I,K,J)
               TL(I,K,J) = TT(I,K,J)
            END DO
         END DO
      END DO
!
! ============================================
      if(wbd_nlitti2) then
         CALL tend_adv_leap       ! 2007.4.23
      else
         CALL tend_adv       ! 2007.4.23
      endif
!      call take_valuet
      CALL nliter_uvt(UL, VL, TL, DTadv)
! ============================================
!------------------------------------------------------------------------------------------------
!
      DO ITERNM = 1,1
!------------------------------------------------------------------------------------------------
! ============================================
         if(wbd_nlitti2) then
            CALL tend_adv_leap       ! 2007.4.23
         else
            CALL tend_adv       ! 2007.4.23
         endif
!      call take_valuet
         CALL nliter_uvt(UL, VL, TL, DTadv)
! ============================================
!------------------------------------------------------------------------------------------------

         DO J = beglatdyn, endlatdyn
            DO K = beglev, endlev
               DO I = 1, NX
                  UT(I,K,J) = (UT(I,K,J) + UL(I,K,J)) * HALF
                  VT(I,K,J) = (VT(I,K,J) + VL(I,K,J)) * HALF
                  TT(I,K,J) = (TT(I,K,J) + TL(I,K,J)) * HALF
               END DO
            END DO
         END DO

!------------------------------------------------------------------------------------------------
! ============================================
         if(wbd_nlitti2) then
            CALL tend_adv_leap       ! 2007.4.23
         else
            CALL tend_adv       ! 2007.4.23
         endif
!      call take_valuet
         CALL nliter_uvt(UL, VL, TL, DTadv)
! ============================================
!------------------------------------------------------------------------------------------------
      END DO

endif !zmh_nlitti2
! ============================ zhh 2007.12.15 =============================
      call t_startf('smooth in nlitti')

      call trans_antiIAP
      CALL SMOTHP  ! SMOOTHING P & PT BY 2-ORDER SHAPIRO SMOOTHER(LON & LAT)
      CALL SMOOTH
      call trans_IAP

      call t_stopf('smooth in nlitti')
! ============================ zhh 2007.12.15 =============================
!zmh moved position of enddo
   end do
! deallocate

   deallocate(PsaL)
   deallocate(PstL)
   deallocate(UL)
   deallocate(VL)
   deallocate(TL)

   !end do   !original location for ITERNM loop zmh
   RETURN
END
