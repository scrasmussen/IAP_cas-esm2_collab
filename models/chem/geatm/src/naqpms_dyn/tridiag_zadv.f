      subroutine tridiag_zadv(a,b,c,d,n,m)
c
c----CAMx v5.40 111010
c
c     TRIDIAG solves a tridiagnal matrix equation with multiple 
c     right-hand side arrays.
c
c     Copyright 1996 - 2010
c     ENVIRON International Corporation
c          
c     Modifications:
c        none
c
c     Input arguments:
c        n                  size of the coefficient matrix
c        m                  number of right-hand side arrays
c        a                  first coefficient array 
c        b                  second coefficient array 
c        c                  third coefficient array 
c        d                  right-hand side arrays 
c 
c     Output arguments:
c        d                  solution arrays 
c
c     Routines called:
c        none
c
c     Called by:
c        VRTSLV
c        VDIFFIMP
c
      implicit none
      integer n,m
      real a(n),b(n),c(n),d(n,m)
      integer nm1,i,i1,k,j
c
c-----Entry point
c
      nm1=n-1
      do 10 i=1,nm1
        c(i)=-c(i)/b(i)
        i1=i+1
        b(i1)=c(i)*a(i1)+b(i1)
  10  continue
c
      do 40 k=1,m
        d(1,k)=d(1,k)/b(1)
        do 20 i=1,nm1
          i1=i+1
          d(i1,k)=(d(i1,k)-d(i,k)*a(i1))/b(i1)
  20    continue
c
        do 30 j=1,nm1
          i=n-j
          i1=i+1
          d(i,k)=d(i1,k)*c(i)+d(i,k)
  30    continue
  40  continue
c
      return
      end
