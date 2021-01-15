
!
!
  subroutine dif_hori( myid, c, kh, deltx, delty, &
                       sx, ex, sy, ey ,nx,ny,dt,ktop,k,ig)
  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, kh, deltx, delty,ktop
  integer i, j
  integer sx1,ex1,sy1,ey1

  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c00

  sx1=sx
  ex1=ex
  sy1=sy
  ey1=ey
!  if(sx1.eq.1)sx1=2
!  if(ex1.eq.nx)ex1=nx-1
!  if(sy1.eq.1)sy1=2
!  if(ey1.eq.ny)ey1=ny-1

  dtstep=10.0*60.0

  nstep=dt/dtstep+1

  dtstep=dt/nstep

do it=1,nstep

  c00=c

  do j = sy1,ey1

    if(j.ge.(ny-4).or.j.le.5) cycle ! chenxsh@20181206

    do i = sx1,ex1
!      dt0=deltx(i,j)*deltx(i,j)/(4.*kh(i,j))
!      iii= dt/dt0+1
!      dtt=dt/iii
!      do it=1,iii
!      c(i,j)=c(i,j)+dtt*kh(i,j)*(c(i+1,j)-2*c(i,j)+c(i-1,j))/ &
      c(i,j)=c00(i,j)+dtstep*kh(i,j)*(c00(i+1,j)-2*c00(i,j)+c00(i-1,j))/ &
                    (deltx(i,j)*deltx(i,j))
!      enddo
    enddo
  enddo

  c00=c

  do j = sy1,ey1

    if(j.ge.(ny-4).or.j.le.5) cycle ! chenxsh@20181206

    do i = sx1,ex1
!      dt0=delty(i,j)*delty(i,j)/(4.*kh(i,j))
!      iii= dt/dt0+1
!      dtt=dt/iii
!      do it=1,iii
!      c(i,j)=c(i,j)+dtt*kh(i,j)*(c(i,j+1)-2*c(i,j)+c(i,j-1))/ &
      c(i,j)=c00(i,j)+dtstep*kh(i,j)*(c00(i,j+1)-2*c00(i,j)+c00(i,j-1))/ &
                    (delty(i,j)*delty(i,j))
!      enddo
    enddo
  enddo

enddo

  return
  end

