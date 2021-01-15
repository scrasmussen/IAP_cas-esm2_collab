
!
!
  subroutine termbal( myid, c0, c ,d, k,sx, ex, sy, ey ,nx, ny, dt)
  integer myid, sx, ex, sy, ey
  real dt
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c0,c,d
  integer i, j


  do j = sy,ey
  do i = sx,ex

   d(i,j) = d(i,j) + c(i,j)-c0(i,j)
   c0(i,j)= c(i,j)
  enddo
  enddo

  return
  end

