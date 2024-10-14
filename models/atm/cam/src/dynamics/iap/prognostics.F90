module prognostics

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Prognostic variables held in-core for convenient access.
! q3 is specific humidity (water vapor) and other constituents.
! pcnst is advected constituents, pnats is non-advected.
! 
! Author: G. Grant
! 
! Revised by: ZhangHe, 2007.4.28
! Updated by: Juanxiong He, 2010.08
! Reviewed by: ZhangHe, 2011-11-19
! Modified by: Jiang Jinrong, October 2012, for 2D parallel
! Reviewed:  Zhang He, 2012-11-01
! Modified by: ZhangHe, 2013-01-17, added xy 2D decomposition variables
!            : Zhang He, 2013-03-21, removed n3m2
!            : Zhang He, 2013-04-15, uxy, vxy, txy --> u3xy, v3xy, t3xy, added Uxy, Vxy 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, beglat, endlat, beglev, endlev, beglonxy, endlonxy, beglatxy, endlatxy
   use infnan,       only: inf
   use constituents, only: pcnst !juanxiong he, 201008

   implicit none

   private

   public ps, u3, v3, t3, q3, qminus, omga, phis, hadv, pdeld
   public phisxy, GHSxy, psxy, omgaxy, lammpxy, phimpxy, sigmpxy
   public WSxy, u3xy, v3xy, t3xy, pdeldxy, qfcstxy, q3xy, Uxy, Vxy
   public n3, n3m1, ptimelevels
   public initialize_prognostics
   public shift_time_indices
   public take_valuet
   public take_value

   integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore
   integer :: n3   = 2                    ! current time level
   integer :: n3m1 = 1                    ! previous time level
!!   integer :: n3m2 = 1                    ! time level before previous
!-----------------------------------------------------------------------
! for CAM_EUL:  three time level  
!           n3m2     n3m1      n3
!             |________|________|     f(n3) = f(n3m2) + df * 2dt
!            n-2      n-1       n
!
! for IAP:  two time level (time n3m2 is useless) 
!           n3m2     n3m1      n3
!             |________|________|     f(n3) = f(n3m1) + df * dt
!            n-2      n-1       n
!-----------------------------------------------------------------------
   real(r8), allocatable :: ps(:,:,:)       ! surface pressure
   real(r8), allocatable :: u3(:,:,:,:)     ! u-wind component
   real(r8), allocatable :: v3(:,:,:,:)     ! v-wind component
   real(r8), allocatable :: t3(:,:,:,:)     ! temperature
   real(r8), allocatable :: pdeld(:,:,:,:)  ! layer thickness dry (Pa)
   real(r8), allocatable :: q3(:,:,:,:,:)   ! specific humidity
   real(r8), allocatable :: qminus(:,:,:,:) ! constituents
   real(r8), allocatable :: hadv  (:,:,:,:) ! horizontal advection tendency
   real(r8), allocatable :: phis(:,:)       ! surface geopotential
   real(r8), allocatable :: omga(:,:,:)     ! vertical velocity
!
   real(r8), allocatable :: phisxy(:,:)       ! surface geopotential
   real(r8), allocatable :: GHSxy(:,:)        ! surface geopotential
   real(r8), allocatable :: psxy(:,:,:)       ! surface pressure
   real(r8), allocatable :: omgaxy(:,:,:)     ! vertical velocity (pressure)
   real(r8), allocatable :: WSxy(:,:,:)       ! vertical velocity (sigma level)
   real(r8), allocatable :: Uxy(:,:,:)        ! u-wind component
   real(r8), allocatable :: Vxy(:,:,:)        ! v-wind component
   real(r8), allocatable :: lammpxy(:,:,:)    ! Lamda midpoint coordinate
   real(r8), allocatable :: phimpxy(:,:,:)    ! Phi midpoint coordinate
   real(r8), allocatable :: sigmpxy(:,:,:)    ! Sigma midpoint coordinate
   real(r8), allocatable :: u3xy(:,:,:,:)      ! u-wind component
   real(r8), allocatable :: v3xy(:,:,:,:)      ! v-wind component
   real(r8), allocatable :: t3xy(:,:,:,:)      ! temperature
   real(r8), allocatable :: pdeldxy(:,:,:,:)  ! layer thickness dry (Pa)
   real(r8), allocatable :: qfcstxy(:,:,:,:)  ! slt forecast of moisture and constituents
   real(r8), allocatable :: q3xy(:,:,:,:,:)   ! specific humidity

!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   subroutine initialize_prognostics
!-----------------------------------------------------------------------------------------------
! Purpose:  Allocate and initialize the prognostic arrays.
!-----------------------------------------------------------------------------------------------
      allocate (ps    (plon                    ,beglat:endlat    ,ptimelevels))
      allocate (u3    (plon,beglev:endlev      ,beglat:endlat,ptimelevels))
      allocate (v3    (plon,beglev:endlev      ,beglat:endlat,ptimelevels))
      allocate (t3    (plon,beglev:endlev      ,beglat:endlat,ptimelevels))
      allocate (q3    (plon,beglev:endlev,pcnst,beglat:endlat,ptimelevels)) !juanxiong he, 201008
      allocate (qminus(plon,beglev:endlev,pcnst,beglat:endlat  ))
      allocate (hadv  (plon,beglev:endlev,pcnst,beglat:endlat  ))
      allocate (phis  (plon,beglat:endlat))        
      allocate (omga  (plon,beglev:endlev,beglat:endlat))    
      allocate (pdeld (plon,beglev:endlev,beglat:endlat,ptimelevels))
!
      allocate(phisxy(beglonxy:endlonxy, beglatxy:endlatxy))
      allocate(GHSxy(beglonxy:endlonxy, beglatxy:endlatxy))
      allocate(psxy(beglonxy:endlonxy, beglatxy:endlatxy, ptimelevels))
      allocate(omgaxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
      allocate(WSxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
      allocate(Uxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
      allocate(Vxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
      allocate(lammpxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
      allocate(phimpxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
      allocate(sigmpxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
      allocate(u3xy(beglonxy:endlonxy, beglatxy:endlatxy, plev, ptimelevels))
      allocate(v3xy(beglonxy:endlonxy, beglatxy:endlatxy, plev, ptimelevels))
      allocate(t3xy(beglonxy:endlonxy, beglatxy:endlatxy, plev, ptimelevels))
      allocate(pdeldxy(beglonxy:endlonxy, beglatxy:endlatxy, plev, ptimelevels))
      allocate(qfcstxy(beglonxy:endlonxy, beglatxy:endlatxy, plev, pcnst))
      allocate(q3xy(beglonxy:endlonxy, beglatxy:endlatxy, plev, ptimelevels, pcnst))

      ps(:,:,:)       = inf
      u3(:,:,:,:)     = inf
      v3(:,:,:,:)     = inf
      t3(:,:,:,:)     = inf
      pdeld(:,:,:,:)  = inf
      q3(:,:,:,:,:)   = inf
      qminus(:,:,:,:) = inf
      hadv  (:,:,:,:) = inf
      phis  (:,:)     = inf
      omga  (:,:,:)   = inf
!
      phisxy(:,:)       = inf
      GHSxy(:,:)        = inf
      psxy(:,:,:)       = inf
      omgaxy(:,:,:)     = inf
      WSxy(:,:,:)       = inf
      Uxy(:,:,:)        = inf
      Vxy(:,:,:)        = inf
      lammpxy(:,:,:)    = inf
      phimpxy(:,:,:)    = inf
      sigmpxy(:,:,:)    = inf
      u3xy(:,:,:,:)      = inf
      v3xy(:,:,:,:)      = inf
      t3xy(:,:,:,:)      = inf
      pdeldxy(:,:,:,:)  = inf
      qfcstxy(:,:,:,:)  = inf
      q3xy(:,:,:,:,:)   = inf

      return
   end subroutine initialize_prognostics

!================================================================================================
   subroutine shift_time_indices
!-----------------------------------------------------------------------------------------------
! Purpose: 
! Shift the indices that keep track of which index stores
! the relative times (current time, previous, time before previous etc).
!-----------------------------------------------------------------------------------------------
      integer :: itmp

      itmp = n3m1

      n3m1 = n3
 
      n3   = itmp
      return       
   end subroutine shift_time_indices
  subroutine take_value
  implicit none
   integer::i,j,k,n,t
   real(kind=4):: temp
  do t=1,ptimelevels
         do j=beglat,endlat
          do k=beglev,endlev
          do i=1,plon
     temp=t3(i,k,j,t)
     t3(i,k,j,t)=dble(temp)
  end do
 end do
 end do
 end do
  do t=1,ptimelevels
         do j=beglat,endlat
          do k=beglev,endlev
          do i=1,plon
     temp=u3(i,k,j,t)
     u3(i,k,j,t)=dble(temp)
  end do
 end do
 end do
 end do
  do t=1,ptimelevels
         do j=beglat,endlat
          do k=beglev,endlev
          do i=1,plon
     temp=v3(i,k,j,t)
     v3(i,k,j,t)=dble(temp)
  end do
 end do
 end do
 end do
  do t=1,ptimelevels
         do j=beglat,endlat
            do n=1,pcnst
          do k=beglev,endlev
          do i=1,plon
     temp=q3(i,k,n,j,t)
     q3(i,k,n,j,t)=dble(temp)
  end do
 end do
 end do
 end do
end do
         do j=beglat,endlat
          do k=beglev,endlev
          do i=1,plon
     temp=omga(i,k,j)
     omga(i,k,j)=dble(temp)
  end do
 end do
 end do
         do j=beglat,endlat
          do i=1,plon
     temp=phis(i,j)
     phis(i,j)=dble(temp)
 end do
 end do
  do t=1,ptimelevels
         do j=beglat,endlat
          do i=1,plon
     temp=ps(i,j,t)
     ps(i,j,t)=dble(temp)
 end do
 end do
end do
         do j=beglat,endlat
            do n=1,pcnst
          do k=beglev,endlev
          do i=1,plon
     temp=qminus(i,k,n,j)
     qminus(i,k,n,j)=dble(temp)
  end do
 end do
 end do
 end do
         do j=beglat,endlat
            do n=1,pcnst
          do k=beglev,endlev
          do i=1,plon
     temp=hadv(i,k,n,j)
     hadv(i,k,n,j)=dble(temp)
  end do
 end do
 end do
 end do
  do t=1,ptimelevels
         do j=beglat,endlat
          do k=beglev,endlev
          do i=1,plon
     temp=pdeld(i,k,j,t)
     pdeld(i,k,j,t)=dble(temp)
  end do
 end do
 end do
 end do

  end subroutine take_value
  subroutine take_valuet
 use IAP_prog,  only : Psa, UT, VT, TT, Pstar1
use pmgrid,     only:beglatdyn,endlatdyn
  integer i,j,k
  real(kind=4) temp
  do j=beglatdyn,endlatdyn
        do k= beglev,endlev
           do i=1,plon
           temp=ut(i,k,j)
           ut(i,k,j)=dble(temp)
  end do
   end do
   end do

  do j=beglatdyn,endlatdyn
        do k= beglev,endlev
           do i=1,plon
           temp=vt(i,k,j)
           vt(i,k,j)=dble(temp)
  end do
   end do
   end do
  do j=beglatdyn,endlatdyn
        do k= beglev,endlev
           do i=1,plon
           temp=tt(i,k,j)
           tt(i,k,j)=dble(temp)
  end do
   end do
   end do
   do j=beglatdyn,endlatdyn
           do i=1,plon
           temp=psa(i,j)
           psa(i,j)=dble(temp)
  end do
   end do
   do j=beglatdyn,endlatdyn
           do i=1,plon
           temp=pstar1(i,j)
           pstar1(i,j)=dble(temp)
  end do
   end do





  end subroutine take_valuet


end module prognostics
