SUBROUTINE FILT2Dadaptive(CH,IV,IP,IC,wf)
!------------------------------------------------------------------------------------------------
! Purpose: PERFORM 2-D FILTER
! Original version : FILT2D.f (IAP 9L)
! Reconstructed    : ZhangHe
! Completed : 2005.9.2
! Update: 2007.5.14, ZhangHe, specify which variables or subroutines to use only
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid, only: NX, NY, JB, JE, period
   use smoother, only: SOSS1D, SOSS2D
   use Filt,     only: FILTER,FILTER_9_13,FILTER_all13,FILTER_noFFT
   use pmgrid,   only: beglatdyn, endlatdyn
   use perf_mod, only : t_startf, t_stopf              !zhh 2013-02-06
   use spmd_utils, only: masterproc, iam


   implicit none
!------------------------------------Arguments---------------------------------------------------
   real(r8), intent(inout) :: CH(NX,NY)     ! input variable needed to be filtered
!!   real(r8), intent(inout) :: CH(NX,beglatdyn:endlatdyn)    ! input variable needed to be filtered
   integer , intent(in)    :: IV, IP, IC, wf    ! index
!----------------------------------Local workspace-----------------------------------------------
   integer  :: J, i             ! loop index
   integer  :: ID,JN,JS
   integer :: yr, mon, day      ! year, month, and day components of a date
   integer :: ncsec             ! current time of day [seconds]
!------------------------------------------------------------------------------------------------

!  APPLY THE SPHERICAL CYCLICITY
   DO J = beglatdyn, endlatdyn
      call period( CH(1,J) )
   END DO
   ID = IP*2 - 1
!------------------------------------------------------------------------------      

   IF (wf.eq.0) THEN        		!modified for test 180926
      call t_startf('FILTER_org')
      CALL FILTER( CH,ID )
      call t_stopf('FILTER_org')
   ELSEIF (wf.eq.1) THEN			!modified for test 180926
      call t_startf('FILTER_9+13')
      CALL FILTER_9_13( CH,ID )  !Near the poles 13 points multiple times, High latitude 9 points once，successfully runing 1.4°*1.4°
      call t_stopf('FILTER_9+13')
   ELSEIF (wf.eq.2) THEN			!modified for test 180926
      call t_startf('FILTER_all13')
      CALL FILTER_all13( CH,ID ) !Near the poles 13 points multiple times, High latitude 13 points once，successfully runing 0.5°*0.5°
      call t_stopf('FILTER_all13')
   ELSEIF (wf.eq.3) THEN			!modified for test 180926
      ! call t_startf('FILTER_all13')
      ! CALL FILTER_all13( CH,ID ) !Near the poles 13 points multiple times, High latitude 13 points once，successfully runing 0.5°*0.5°
      ! call t_stopf('FILTER_all13')

	  call t_startf('FILTER_noFFT')
      CALL FILTER_noFFT( CH,ID )
      call t_stopf('FILTER_noFFT')
   endif

!------------------------------------------------------------------------------      
   IF ( IC.GE.1 ) THEN
      IF ( IV.EQ.1 ) THEN  ! do second-order Shapiro smoother of 1D (LAT) for u & v
        JN = JB - IP + 1
         JS = JE
!------------------------------------------------------------------------------      
        CALL SOSS1D( CH,JN,JS,IC )  ! at module smoother
!------------------------------------------------------------------------------      
     ELSE  ! do second-order Shapiro smoother of 2D (LAT & LON) for PS, T & q
!------------------------------------------------------------------------------      
          CALL SOSS2D( CH,IV,IP,IC ) 
!------------------------------------------------------------------------------      
      ENDIF
   ENDIF
   RETURN
END
