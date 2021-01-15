SUBROUTINE tend_adv2_leap
!------------------------------------------------------------------------------------------------
! Purpose: Compute the tendence of T, U, V ( Advection(evolution) part )
! Author : ZhangHe
! Completed : 2005.9.1
! Update    : 2006.2.22, Zhanghe, revise to advection part
!             2007.3.6, Zhanghe, 1) change subroutine's name from ptuvtend2 to tend_adv 
!                                2) delete the dumb parameter 'NCTDCB' 
!                                3) delete computation of GHI & module 'physconst' & 'stdatm'
!             2007.4.23, ZhangHe, change 'LADDSS' to 'NADDSS'
!             2008.5, WuJianping, parallel version
!             2008.6.11, ZhangHe, available for both serial & parallel
! Reviewed: ZhangHe, 2011-11-19
! Modified: Jiang Jinrong, 2012 October, for 2D parallel
! Reviewed: Zhang He, 2012-11-13
! Modified: Zhang He, 2013-01-28, new mp_send3d, dyn_state is added
!           Zhang He, 2013-04-01, use dynamic arrays UZ et. al
! Modifed   Zhang Minghua 2019-09-25: 1) implict forward-back integration together with nlitti.F90 
!                                     2) name changed to tend_adv2, added var to flexib
! Modified  Cao Hang 2018-11-16:  Added adptive leap-cofficient in high-lats
!
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only: NX, NY, NL, NZ, KE, JB, JE, IB, IE, period, NLON
   use Dyn_const 
   use IAP_prog,  only: TT, UT, VT, Ustar, Vstar, WSstar 
   use tendency,  only: DU, DV, DT, SU, SV, ST, NADDSS
   use mathconst, only: ZERO
! zmh
   !use flexib,    only: a11, a12, a13, a21, a22, a23, a31, a32, a33, kappa2, IBCFFT
   use flexib,    only: a11, a12, a13, a21, a22, a23, a31, a32, a33, kappa2, IBCFFT, &
                        zmh_nlitti2, DTadv          !zmh added flag and time step

   use pmgrid,    only: beglatdyn, endlatdyn, loc_JB, loc_JE,beglev,endlev,npr_y,&
                        myid_z,beglatdynex,endlatdynex,npr_z
   use spmd_utils, only: masterproc, iam   ! debug
   use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-01-28
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-28
   use perf_mod, only : t_startf, t_stopf              !zhh  2013-02-06
#if ( defined SPMD )
  use mod_comm, only: mp_send3d, mp_recv3d,mp_send3d_2, mp_recv3d_2
#endif


   IMPLICIT NONE
!----------------------------------Local workspace-----------------------------------------------
!zmh
   real(r8), allocatable :: UL(:,:,:)
   real(r8), allocatable :: VL(:,:,:)
   real(r8), allocatable :: TL(:,:,:)

   real(r8), allocatable :: UZ(:,:,:), VZ(:,:,:), TZ(:,:,:)   ! extended variables in vertical layer
   real(r8), allocatable :: WW(:,:)         ! wjp
   real(r8) :: TL1, TL2, TL3                ! advection terms of TT
   real(r8) :: VL1, VL2, VL3                ! advection terms of VT
   real(r8) :: UL1, UL2, UL3                ! advection terms of UT
   real(r8) :: TLP, TL3P                    ! advection terms of TT at polar
   real(r8) :: DTTP                         ! tendency of TT at polar
   real(r8) :: OYVS                         ! auxiliary variables
   real(r8) :: USVT0, USVT1, USVT2, USVT3   ! auxiliary variables
   real(r8) :: time_begin, time_end, ierr
   integer  :: KK, I, J, K,src,dest,src1,dest1,src2,dest2,src3,dest3                  ! loop index
!zmh
   integer  :: JM1,JP1,KP1,KM1, IZMH, NZMH, ITERNM !zmh                 !zmh loop index added

   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-28
   integer  :: commyz     !zhh 2013-01-28
   integer  :: I_leap1, I_leap1_1, I_leap2	, I_leap3, I_leap4, I_leap5, I_leap6	!For leap points CH 181116


!----------------------------------------------------------------------------------------------
! allocate arrays
   allocate ( UZ(NX,beglev:endlev+2,beglatdyn:endlatdyn) )
   allocate ( VZ(NX,beglev:endlev+2,beglatdyn:endlatdyn) )
   allocate ( TZ(NX,beglev:endlev+2,beglatdyn:endlatdyn) )
   allocate ( WW(NX,NY) )
!zmh
   allocate ( UL(NX,beglev:endlev,beglatdyn:endlatdyn) )
   allocate ( VL(NX,beglev:endlev,beglatdyn:endlatdyn) )
   allocate ( TL(NX,beglev:endlev,beglatdyn:endlatdyn) )

   call t_startf('tend_adv2')
!zmh
!  ====================================================
 DO ITERNM = 1, 1  ! Iteration times
!  ====================================================

   DU(IB:IE,beglev:endlev,beglatdyn:endlatdyn) = ZERO
   DV(IB:IE,beglev:endlev,beglatdyn:endlatdyn) = ZERO
   DT(IB:IE,beglev:endlev,beglatdyn:endlatdyn) = ZERO

    NZMH = 1
  if(zmh_nlitti2)then 

    NZMH = 2
      DO J = beglatdyn, endlatdyn
         DO K = beglev,endlev 
            DO I = 1, NX
               UL(I,K,J) = UT(I,K,J)
               VL(I,K,J) = VT(I,K,J)
               TL(I,K,J) = TT(I,K,J)
            END DO
         END DO
      END DO
   endif

!  ====================================================
 DO IZMH = 1, NZMH   ! forward and backward integration
!  ====================================================

   dyn_state => get_dyn_state()
   commyz = dyn_state%grid%commyz

!  *******************************************************************
!  *******************************************************************
!  *******************          COMPUTE DT        ********************
!  *******************************************************************
!  *******************************************************************
#if (defined SPMD)
!jjr   call mpi_move_left(Vstar(1,1,endlatdyn),Vstar(1,1,beglatdyn-1),NX*NL)
!jjr   call mpi_move_left(   TT(1,1,endlatdyn),   TT(1,1,beglatdyn-1),NX*NL)
!jjr   call mpi_move_right(    TT(1,1,beglatdyn),   TT(1,1,endlatdyn+1),NX*NL)
!right 
! call t_startf('tend_adv_barrcommy1')
! call MPI_BARRIER(commyz,ierr)
! call t_stopf('tend_adv_barrcommy1') 
   call t_startf('tend_adv_comm1')
   
   src = iam-1
   dest  = iam+1
   if ( mod(iam,npr_y) == 0 ) src = -1
   if ( mod(iam+1,npr_y) == 0 ) dest = -1
   call mp_send3d( commyz, dest, src, NX,NL,NY,                     &
                   1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                   1, NX,beglev,endlev, beglatdyn, beglatdyn,  tt )
   call mp_recv3d( commyz, src, NX,NL,  NY,                            &
                   1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                   1, NX,beglev,endlev, endlatdynex, endlatdynex,tt )

!left      
   src1 = iam+1
   dest1  = iam-1
   if ( mod(iam,npr_y) == 0 ) dest1 = -1
   if ( mod(iam+1,npr_y) == 0 ) src1 = -1
   call mp_send3d_2( commyz, dest1, src1, NX,NL, NY,                     &
                     1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                     1, NX,beglev,endlev, endlatdyn, endlatdyn,  TT,Vstar )
   call t_stopf('tend_adv_comm1')
#endif
!
!     COMPUTE    DTT/DT

   DO J = beglatdyn, endlatdyn

      DO K = beglev+1 ,endlev+1
         DO I = 1 ,NX
            TZ(I,K,J) = TT(I,K-1,J)
            UZ(I,K,J) = UT(I,K-1,J)
            VZ(I,K,J) = VT(I,K-1,J)
         ENDDO
      ENDDO
   ENDDO
!
   if (myid_z.eq.0) then
      do j=beglatdyn,endlatdyn
         DO K = 1 ,KE,NZ
            DO I = 1 ,NX
               UZ(I,1,J) = ZERO
               VZ(I,1,J) = ZERO
               TZ(i,1,j) = ZERO
            ENDDO
         ENDDO
      ENDDO  !jjr
   endif
   if (myid_z.eq.(npr_z-1)) then
      do j=beglatdyn,endlatdyn
         DO K = 1 ,KE,NZ
            DO I = 1 ,NX
               UZ(I,KE,J) = ZERO
               VZ(I,KE,J) = ZERO
               TZ(i,KE,j) = ZERO
            ENDDO
         ENDDO
      ENDDO  !jjr
   endif

#if (defined SPMD)
!right
!left
! call t_startf('tend_adv_barrcommy2')
! call MPI_BARRIER(commyz,ierr)
! call t_stopf('tend_adv_barrcommy2') 
   call t_startf('tend_adv_comm2')
   call mp_recv3d_2( commyz, src1, NX,NL,  NY,                             &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX,beglev,endlev, beglatdynex, beglatdynex, TT,Vstar )

!up

   src2=iam-npr_y
   dest2=iam+npr_y
   if(myid_z.eq.0) src2=-1
   if(myid_z.eq.(npr_z-1)) dest2=-1
   call mp_send3d( commyz, dest2, src2, NX,KE,NY,                     &
                   1, NX,beglev,endlev+2,beglatdyn ,endlatdyn,        &
                   1, NX,endlev+1,endlev+1, beglatdyn, endlatdyn,  TZ )
   call mp_recv3d( commyz, src2, NX,KE,NY,                     &
                   1, NX,beglev,endlev+2,beglatdyn ,endlatdyn,        &
                   1, NX,beglev,beglev, beglatdyn, endlatdyn,  TZ )
   call mp_send3d_2( commyz, dest2, src2, NX,KE,NY,                     &
                     1, NX,beglev,endlev+2,beglatdyn ,endlatdyn,        &
                     1, NX,endlev+1,endlev+1, beglatdyn,endlatdyn,UZ,VZ   )
   call mp_recv3d_2( commyz, src2, NX,KE,NY,                     &
                     1, NX,beglev,endlev+2,beglatdyn ,endlatdyn,        &
                     1, NX,beglev,beglev, beglatdyn, endlatdyn,  UZ,VZ )
 
!down
   src3=iam+npr_y
   dest3=iam-npr_y
   if(myid_z.eq.0) dest3=-1
   if(myid_z.eq.(npr_z-1)) src3=-1
   call mp_send3d( commyz, dest3, src3, NX,KE,NY,                     &
                   1, NX,beglev,endlev+2,beglatdyn ,endlatdyn,        &
                   1, NX,beglev+1,beglev+1, beglatdyn, endlatdyn,  TZ )
   call mp_recv3d( commyz, src3, NX,KE,NY,                     &
                   1, NX,beglev,endlev+2,beglatdyn ,endlatdyn,        &
                   1, NX,endlev+2,endlev+2, beglatdyn, endlatdyn,  TZ )
   call mp_send3d_2( commyz, dest3, src3, NX,KE,NY,                     &
                     1, NX,beglev,endlev+2,beglatdyn ,endlatdyn,        &
                     1, NX,beglev+1,beglev+1, beglatdyn, endlatdyn,  UZ,VZ )
   call mp_recv3d_2( commyz, src3, NX,KE,NY,                     &
                     1, NX,beglev,endlev+2,beglatdyn ,endlatdyn,        &
                     1, NX,endlev+2,endlev+2, beglatdyn, endlatdyn,UZ,VZ )
   ! call t_stopf('comm1')

!right
   call mp_send3d_2( commyz, dest, src, NX,NL,NY,                     &
                     1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                     1, NX,beglev,endlev, beglatdyn, beglatdyn,  UT,VT )
   call mp_recv3d_2( commyz, src, NX,NL,  NY,                            &
                     1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                     1, NX,beglev,endlev, endlatdynex, endlatdynex,UT,VT )

   call mp_send3d_2( commyz, dest, src, NX,NL,NY,                     &
                     1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                     1, NX,beglev,endlev, beglatdyn, beglatdyn, Ustar,Vstar  )
   call mp_recv3d_2( commyz, src, NX,NL,  NY,                            &
                     1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                     1, NX,beglev,endlev, endlatdynex, endlatdynex,Ustar,Vstar )

   call mp_send3d( commyz, dest, src, NX,NZ,NY,                     &
                   1, NX,beglev,endlev+1, beglatdynex, endlatdynex,       &
                   1, NX,beglev,endlev+1, beglatdyn, beglatdyn,  WSstar )
   call mp_recv3d( commyz, src, NX,NZ,  NY,                            &
                   1, NX,beglev,endlev+1, beglatdynex, endlatdynex,       &
                   1, NX,beglev,endlev+1, endlatdynex, endlatdynex,WSstar )

!left
   call mp_send3d_2( commyz, dest1, src1, NX,NL, NY,                     &
                     1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                     1, NX,beglev,endlev, endlatdyn, endlatdyn,  UT,VT )
   call mp_recv3d_2( commyz, src1, NX,NL,  NY,                             &
                     1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                     1, NX,beglev,endlev, beglatdynex, beglatdynex, UT,VT )
 
   call t_stopf('tend_adv_comm2')
#endif
   call t_startf('tend_adv_comput1')
   DO J = beglatdyn, endlatdyn !jjr
      IF (J.GE.JB.AND.J.LE.JE) THEN
!------------------ COMPUTE DTT/DT FOR J [JB,JE]  &  I [IB,IE]  &  K [1,NL] ---------------------
         DO K = beglev ,endlev
            DO I = IB,IE
               ! TL1  = OUXH(J)*(alpha(J)*(Ustar(I+1,K,J)*TT(I+1,K,J)-Ustar(I  ,K,J)*TT(I-1,K,J)) &   
                         ! + betaw(J)/3.0*(Ustar(I+2,K,J)*TT(I+3,K,J)-Ustar(I-1,K,J)*TT(I-3,K,J))) 
					                                                                      ! !(4.13)        
               ! TL2  = RUPH(J)*Vstar(I,K,J)*TT(I,K,J+1)-RUMH(J)*Vstar(I,K,J-1)*TT(I,K,J-1) !(4.14)
               ! TL3  = SGH(K) * (WSstar(I,K+1,J)*TZ(I,K+2,J)  &                            ! TZ(K)=TT(K-1)
                             ! -  WSstar(I,K,J)*TZ(I,K,J))                                  !(4.15)
               ! DT(I,K,J)  = - a31*TL1 - a32*TL2 - a33*TL3                         ! (3.3)
			   
			  !CH modified 181116
			   I_leap1   = mod(I+1*k_leap(J),NX) + (I+1*k_leap(J))/NX * (IB+2)			!I+1
			   I_leap1_1 = mod(I+2*k_leap(J)-1,NX) + (I+2*k_leap(J)-1)/NX * (IB+2)		!I+1  (2)
			   I_leap2   = (k_leap(J)-1+NLON)/(I+NLON) * NLON + (I-(k_leap(J)-1))		!I
			   I_leap3   = (2*k_leap(J)-1+NLON)/(I+NLON) * NLON + (I-(2*k_leap(J)-1))	!I-1

               TL1  = OUXH(J)*(Ustar(I_leap1 ,K,J)*TT(I_leap1_1 ,K,J)-Ustar(I_leap2  ,K,J)*TT(I_leap3 ,K,J)) /abs(2*k_leap(J)-1)                                                                                       
               TL2  = RUPH(J)*Vstar(I,K,J)*TT(I,K,J+1)-RUMH(J)*Vstar(I,K,J-1)*TT(I,K,J-1) !(4.14)
               TL3  = SGH(K) * (WSstar(I,K+1,J)*TZ(I,K+2,J)  &                            ! TZ(K)=TT(K-1)
                             -  WSstar(I,K,J)*TZ(I,K,J))                                  !(4.15)
               DT(I,K,J)  = - a31*TL1 - a32*TL2 - a33*TL3                         ! (3.3)
			   
            ENDDO
            call period( DT(1,K,J) )
         ENDDO
!------------------------------- FOR THE POLAR POINTS -----------------------------------
      ELSE IF(J.EQ.1) THEN         ! at north polar
!---------------------------------COMPUTE DTT/DT AT NORTH POLE------------------------------------
         DO K = beglev ,endlev
!           COMPUTE HORIZONTAL & VERTICAL ADVECTION OF TT
            TLP    = ZERO
            DO I   = IB,IE
               TLP = TLP + Vstar(I,K,1)*TT(I,K,JB)                                     !(4.14)-(2)
            ENDDO
            TLP    = RDLATI * TLP    
            TL3P   = SGH(K) * (WSstar(IB,K+1,1)*TZ(IB,K+2,1) - WSstar(IB,K,1)*TZ(IB,K,1))
            DTTP   = -a32*TLP - a33*TL3P 
            DO I = 1 ,NX
               DT(I,K,1) = DTTP
            ENDDO
         ENDDO
!
      ELSE                        ! J = NY  South polar
!-------------------------------- COMPUTE DTT/DT AT south POLES -------------------------------------
         DO K = beglev ,endlev
!           COMPUTE HORIZONTAL & VERTICAL ADVECTION OF TT
            TLP    = ZERO
            DO I   = IB,IE
               TLP = TLP + Vstar(I,K,JE)*TT(I,K,JE)                                    !(4.14)-(3)
            ENDDO
            TLP    = - RDLATI * TLP   
            TL3P   = SGH(K) * (WSstar(IB,K+1,NY)*TZ(IB,K+2,NY) - WSstar(IB,K,NY)*TZ(IB,K,NY))
            DTTP   = -a32*TLP - a33*TL3P 
            DO I = 1 ,NX
               DT(I,K,NY) = DTTP
            ENDDO
         ENDDO
      ENDIF
!=============================== zhh =====================================
!!      DO K = beglev ,endlev
!!         DO I = 1 ,NX
!!            IF (abs(DT(I,K,J)) > 1E-1 ) THEN
!!               print*, 'DT(',I,K,J,') =', DT(I,K,J)
!!               print*, 'tend_adv2--1'
!!               stop
!!            ENDIF
!!         ENDDO
!!      ENDDO
!============================ 2013.03.07 ====================================
      IF ( NADDSS == 2 ) THEN
!         ADD SOURCE & SINK
         DO K = beglev ,endlev
            DO I = 1 ,NX
               DT(I,K,J) = DT(I,K,J) + ST(I,K,J)
            ENDDO
         ENDDO
      ENDIF
   ENDDO  !(end for DO J = 1 ,NY)
   
   call t_stopf('tend_adv_comput1')
!
!     HIGH-MID LAT FILTER  FOR DP/DT & DTT/DT
!----------------------------------------------------------------
   call t_startf('tendadv_FILT2Dnew2_DUDV')
   DO KK = beglev,endlev
      WW(:,beglatdyn:endlatdyn)=DT(:,KK,beglatdyn:endlatdyn)
      ! CALL FILT2D( WW,0,1,IBCFFT )   
      CALL FILT2Dadaptive( WW,0,1,IBCFFT,3 )
      DT(:,KK,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
   ENDDO
   call t_stopf('tendadv_FILT2Dnew2_DUDV')
!----------------------------------------------------------------
!
!     *******************************************************************
!     *******************************************************************
!     ****************         COMPUTE DU & DV        *******************
!     *******************************************************************
!     *******************************************************************
!
   call t_startf('tendadv_comput2')
   DO J  = beglatdyn, endlatdyn !jjr
!
!         COMPUTE DV/DT AT J=1
!
      IF (J.EQ.1) THEN    !AT THE NORTH POLAR  
         DO K = beglev ,endlev
            DO I = IB,IE
               ! USVT0  = Ustar(I  ,K,JB) * VT(I-1,K,1)
               ! USVT1  = Ustar(I+1,K,JB) * VT(I+1,K,1)
               ! USVT2  = Ustar(I-1,K,JB) * VT(I-3,K,1)
               ! USVT3  = Ustar(I+2,K,JB) * VT(I+3,K,1)
               ! VL1 = OVXQ(1) * (alpha2(1) * ( USVT1 - USVT0 )                       & 
                   ! +      betaw2(1) / 3.0 * ( USVT3 - USVT2 ) )                      ! (4.1)
               ! VL2 = ( RVP(1)*Vstar(I,K,JB) + RDLATQ*Vstar(I,K,1) ) * VT(I,K,JB)         
               ! VL3 = SGQ(K) * ((WSstar(I,K+1,JB) + WSstar(IB,K+1,1))* VZ(I,K+2,1)   &
                   ! -           (WSstar(I,K,JB  ) + WSstar(IB,K,  1))* VZ(I,K,1) )
               ! DV(I,K,1) = - a11*VL1 - a12*VL2 - a13*VL3         ! (3.1)	
               ! DU(I,K,1) = ZERO
			   
			  !CH modified 181116
			   I_leap1   = mod(I+1*k_leap2(J),NX) + (I+1*k_leap2(J))/NX * (IB+2)			!I+1 -> I+2
			   I_leap1_1 = mod(I+2*k_leap2(J)-1,NX) + (I+2*k_leap2(J)-1)/NX * (IB+2)		!I+1 -> I+3 
			   I_leap2   = (k_leap2(J)-1+NLON)/(I+NLON) * NLON + (I-(k_leap2(J)-1))			!I 	 -> I-1
			   I_leap3   = (2*k_leap2(J)-1+NLON)/(I+NLON) * NLON + (I-(2*k_leap2(J)-1))		!I-1 -> I-3			

               USVT0  = Ustar(I_leap2  ,K,JB) * VT(I_leap3,K,1)
               USVT1  = Ustar(I_leap1,K,JB) * VT(I_leap1_1,K,1)
               USVT2  = Ustar(I-1,K,JB) * VT(I-3,K,1)
               USVT3  = Ustar(I+2,K,JB) * VT(I+3,K,1)

               VL1 = OVXQ(1)  * ( USVT1 - USVT0 )  /abs(2*k_leap2(J)-1)   
			   
               VL2 = ( RVP(1)*Vstar(I,K,JB) + RDLATQ*Vstar(I,K,1) ) * VT(I,K,JB)         
               VL3 = SGQ(K) * ((WSstar(I,K+1,JB) + WSstar(IB,K+1,1))* VZ(I,K+2,1)   &
                   -           (WSstar(I,K,JB  ) + WSstar(IB,K,  1))* VZ(I,K,1) )
               DV(I,K,1) = - a11*VL1 - a12*VL2 - a13*VL3         ! (3.1)  
               DU(I,K,1) = ZERO

            ENDDO
            call period( DU(1,K,J) )
            call period( DV(1,K,J) )
         ENDDO
!
      ELSE IF(J.LE.JE) THEN
! ----------------- COMPUTE DU/DT & DV/DT FOR J [JB,JE] & I [IB,IE] & K [1,NL] ------------------
!         
         DO K = beglev ,endlev
            DO I = IB,IE
               ! OYVS = RDLATQ * Vstar(I,K,J)
! ! for UT =================================================================================
               ! UL1  = OUXQ(J) * ( alpha(J) * ((Ustar(I  ,K,J)+Ustar(I+1,K,J))*UT(I+1,K,J)   &
                                           ! -  (Ustar(I  ,K,J)+Ustar(I-1,K,J))*UT(I-1,K,J))  &   
                            ! + betaw(J)/3.0 * ((Ustar(I+2,K,J)+Ustar(I+1,K,J))*UT(I+3,K,J)   &
	        		           ! -  (Ustar(I-1,K,J)+Ustar(I-2,K,J))*UT(I-3,K,J))  ) !(4.7)
               ! UL2  = RUPQ(J) * (Vstar(I,  K,J) + Vstar(I-1  ,K,J)) * UT(I,K,J+1)           &
                    ! - RUMQ(J) * (Vstar(I,K,J-1) + Vstar(I-1,K,J-1)) * UT(I,K,J-1)   ! (4.8)
               ! UL3  = SGQ(K) * ((WSstar(I,K+1,J)+ WSstar(I-1,K+1,J))* UZ(I,K+2,J)           &                
                    ! -           (WSstar(I,K,J)  + WSstar(I-1,K,J )) * UZ(I,K,J) )   ! (4.9)
! ! for VT =================================================================================
               ! VL1  = OVXQ(J)*( alpha2(J)*((Ustar(I+1,K,J+1)+Ustar(I+1,K,J))*VT(I+1,K,J)    &
                                        ! -  (Ustar(I  ,K,J+1)+Ustar(I  ,K,J))*VT(I-1,K,J))   &   
                    ! +      betaw2(J)/3.0 *((Ustar(I+2,K,J+1)+Ustar(I+2,K,J))*VT(I+3,K,J)    &
		        	        ! -  (Ustar(I-1,K,J+1)+Ustar(I-1,K,J))*VT(I-3,K,J))  )
! !					                                            ! (4.1)  
               ! VL2  = ( RVP(J) * Vstar(I,K,J+1) + OYVS ) * VT(I,K,J+1)                          &
                    ! - ( RVM(J) * Vstar(I,K,J-1) + OYVS ) * VT(I,K,J-1)              ! (4.2)
               ! VL3  = SGQ(K) * ( (WSstar(I,K+1,J+1) + WSstar(I,K+1,J))*VZ(I,K+2,J)          &
                    ! -            (WSstar(I,K,J+1  ) + WSstar(I,K,J  ))*VZ(I,K,J) )  ! (4.3)
               ! DU(I,K,J) = - a21*UL1 - a22*UL2 - a23*UL3      ! (3.2)    
               ! DV(I,K,J) = - a11*VL1 - a12*VL2 - a13*VL3      ! (3.1)   
			   
			!CH modified 181116
               OYVS = RDLATQ * Vstar(I,K,J)
			   I_leap1   = mod(I+1*k_leap(J),NX)   + (I+1*k_leap(J))/NX * (IB+2)			!I+1 -> I+2
			   I_leap2   = mod(I+k_leap(J)-1,NX)   + (I+k_leap(J)-1)/NX * (IB+2)		!I -> I+1  
			   I_leap3   = mod(I+2*k_leap(J)-1,NX) + (I+2*k_leap(J)-1)/NX * (IB+2)		!I+1 -> I+3  
			   
			   I_leap4   = (k_leap(J)-1+NLON)/(I+NLON)   * NLON + (I-(k_leap(J)-1))		!I -> I-1
			   I_leap5   = (k_leap(J)+NLON)/(I+NLON)     * NLON + (I-k_leap(J))		!I-1 -> I-2
			   I_leap6   = (2*k_leap(J)-1+NLON)/(I+NLON) * NLON + (I-(2*k_leap(J)-1))	!I-1 -> I-3

               UL1  = OUXQ(J)  * ((Ustar(I_leap1  ,K,J) + Ustar(I_leap2,K,J)) * UT(I_leap3,K,J)   &
                               -  (Ustar(I_leap4  ,K,J) + Ustar(I_leap5,K,J)) * UT(I_leap6,K,J))  /abs(2*k_leap(J)-1)  
										   
               UL2  = RUPQ(J) * (Vstar(I,  K,J) + Vstar(I-1  ,K,J)) * UT(I,K,J+1)           &
                    - RUMQ(J) * (Vstar(I,K,J-1) + Vstar(I-1,K,J-1)) * UT(I,K,J-1)   ! (4.8)
               UL3  = SGQ(K) * ((WSstar(I,K+1,J)+ WSstar(I-1,K+1,J))* UZ(I,K+2,J)           &                
                    -           (WSstar(I,K,J)  + WSstar(I-1,K,J )) * UZ(I,K,J) )   ! (4.9)
! for VT =================================================================================
			   !CH modified 181116
			   I_leap1   = mod(I+1*k_leap2(J),NX) + (I+1*k_leap2(J))/NX * (IB+2)			!I+1 -> I+2
			   I_leap1_1 = mod(I+2*k_leap2(J)-1,NX) + (I+2*k_leap2(J)-1)/NX * (IB+2)		!I+1 -> I+3 
			   I_leap2   = (k_leap2(J)-1+NLON)/(I+NLON) * NLON + (I-(k_leap2(J)-1))			!I 	 -> I-1
			   I_leap3   = (2*k_leap2(J)-1+NLON)/(I+NLON) * NLON + (I-(2*k_leap2(J)-1))		!I-1 -> I-3

               VL1  = OVXQ(J)*((Ustar(I_leap1  ,K,J+1)  + Ustar(I_leap1  ,K,J)) * VT(I_leap1_1 ,K,J)	&    			! equation(2.123)
                            -  (Ustar(I_leap2  ,K,J+1)  + Ustar(I_leap2  ,K,J)) * VT(I_leap3   ,K,J))  /abs(2*k_leap2(J)-1) 		   
																								
               VL2  = ( RVP(J) * Vstar(I,K,J+1) + OYVS ) * VT(I,K,J+1)                          &
                    - ( RVM(J) * Vstar(I,K,J-1) + OYVS ) * VT(I,K,J-1)              ! (4.2)
               VL3  = SGQ(K) * ( (WSstar(I,K+1,J+1) + WSstar(I,K+1,J))*VZ(I,K+2,J)          &
                    -            (WSstar(I,K,J+1  ) + WSstar(I,K,J  ))*VZ(I,K,J) )  ! (4.3)
               DU(I,K,J) = - a21*UL1 - a22*UL2 - a23*UL3      ! (3.2)    
               DV(I,K,J) = - a11*VL1 - a12*VL2 - a13*VL3      ! (3.1)  

            ENDDO
            call period( DU(1,K,J) )
            call period( DV(1,K,J) )
         ENDDO
!
      ELSE   !  J = NY  AT THE SOUTH POLAR
         DO K = beglev ,endlev
            DO I = 1 ,NX
               DU(I,K,J)= ZERO
               DV(I,K,J)= ZERO
            ENDDO
         ENDDO
      ENDIF

!=============================== zhh =====================================
      DO K = beglev ,endlev
         DO I = 1 ,NX
            IF (abs(DU(I,K,J)) > 5E-1 .OR. abs(DV(I,K,J)) > 5E-1) THEN
               print*, 'DU(',I,K,J,') =', DU(I,K,J)
               print*, 'DV(',I,K,J,') =', DV(I,K,J)
               print*, 'tend_adv2--2'
!!               stop
            ENDIF
         ENDDO
      ENDDO
!============================ 2007.8.5 ====================================
      IF ( NADDSS == 2 ) THEN
         DO K = beglev ,endlev
            DO I = 1 ,NX
               IF (DU(I,K,J).GT.1E9.OR.SU(I,K,J).GT.1E9) THEN
                  PRINT*,'SU:',SU(I,K,J),I,K,J
                  PRINT*,'DU:',DU(I,K,J),I,K,J
                  PRINT*,'DV:',DV(I,K,J),I,K,J
!!                  PAUSE 'DTPUV--1'
               ENDIF
               DU(I,K,J) = DU(I,K,J) + SU(I,K,J)
               DV(I,K,J) = DV(I,K,J) + SV(I,K,J)
!!               print*, 'SU(10,26,30) =', SU(10,26,30)
            ENDDO
         ENDDO
      ENDIF
   ENDDO        !end J = 1,NY
   call t_stopf('tendadv_comput2')
!
!     HIGH-MID LAT FILTER  FOR DU/DT & DV/DT
   call t_startf('tendadv_FILT2Dnew2_DUDV')
   DO K = beglev,endlev
      WW(:,beglatdyn:endlatdyn)=DU(:,K,beglatdyn:endlatdyn)
      ! CALL FILT2D( WW,1,1,IBCFFT )	
      CALL FILT2Dadaptive( WW,1,1,IBCFFT,3 )
      DU(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!
      WW(:,beglatdyn:endlatdyn)=DV(:,K,beglatdyn:endlatdyn)
      ! CALL FILT2D( WW,1,2,IBCFFT )	
      CALL FILT2Dadaptive( WW,1,2,IBCFFT,3 ) 
      DV(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
   ENDDO
!
   call t_stopf('tendadv_FILT2Dnew2_DUDV')
!zmh
! -----------------------------------
    IF(zmh_nlitti2)then
     CALL nliter_uvt(UL, VL, TL, DTadv)  
    ENDIF

! ===============================================
  ENDDO  ! IZMH forth and back integration

 ENDDO   ! ITERNM iteratioins
! ===============================================



   call t_stopf('tend_adv2')

! deallocate arrays
   deallocate(WW)
   deallocate(UZ)
   deallocate(VZ)
   deallocate(TZ)
!zmh
   deallocate(UL)
   deallocate(VL)
   deallocate(TL)


   RETURN
END
