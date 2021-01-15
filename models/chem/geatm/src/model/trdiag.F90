      subroutine trdiag(a,b,c,d,n,m)
!
!----CAMx v4.40 061025
!
!     TRDIAG solves a tridiagnal matrix equation with multiple 
!     right-hand side arrays.
!
!     Copyright 1996-2006
!     ENVIRON International Corporation
!          
!     Modifications:
!        none
!
!     Input arguments:
!        n                  size of the coefficient matrix
!        m                  number of right-hand side arrays
!        a                  first coefficient array 
!        b                  second coefficient array 
!        c                  third coefficient array 
!        d                  right-hand side arrays 
! 
!     Output arguments:
!        d                  solution arrays 
!
!     Routines called:
!        none
!
!     Called by:
!        VRTSLV
!        VDIFFIMP
!
      dimension a(n),b(n),c(n),d(n,m)
!
!-----Entry point
!
      nm1=n-1
      do 10 i=1,nm1
        c(i)=-c(i)/b(i)
        i1=i+1
        b(i1)=c(i)*a(i1)+b(i1)
  10  continue

      do 40 k=1,m
        d(1,k)=d(1,k)/b(1)
        do 20 i=1,nm1
          i1=i+1
          d(i1,k)=(d(i1,k)-d(i,k)*a(i1))/b(i1)
  20    continue

        do 30 j=1,nm1
          i=n-j
          i1=i+1
          d(i,k)=d(i1,k)*c(i)+d(i,k)
  30    continue
  40  continue

      return
      end