
!
!
  subroutine dif_hori( myid, c, dso4, dno3, dfeIII, dfeII,kh, deltx, delty, &
                       sx, ex, sy, ey ,nx,ny,dt,ktop,k,ig)
  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, kh, deltx, delty,ktop
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: dso4,dno3,dfeIII,dfeII
  integer i, j
  integer sx1,ex1,sy1,ey1

  sx1=sx
  ex1=ex
  sy1=sy
  ey1=ey
!  if(sx1.eq.1)sx1=2
!  if(ex1.eq.nx)ex1=nx-1
!  if(sy1.eq.1)sy1=2
!  if(ey1.eq.ny)ey1=ny-1

  do j = sy1,ey1
    do i = sx1,ex1
!      dt0=deltx(i,j)*deltx(i,j)/(4.*kh(i,j))
!      iii= dt/dt0+1
!      dtt=dt/iii
!      do it=1,iii
!      c(i,j)=c(i,j)+dtt*kh(i,j)*(c(i+1,j)-2*c(i,j)+c(i-1,j))/ &
      c(i,j)=c(i,j)+dt*kh(i,j)*(c(i+1,j)-2*c(i,j)+c(i-1,j))/ &
                    (deltx(i,j)*deltx(i,j))
      dso4(i,j) = dso4(i,j) + dt * kh(i,j)* &
                ( dso4(i+1,j)-2*dso4(i,j)+dso4(i-1),j)/(deltx(i,j)*deltx(i,j))
      dno3(i,j) = dno3(i,j) + dt * kh(i,j)* &
                ( dno3(i+1,j)-2*dno3(i,j)+dno3(i-1),j)/(deltx(i,j)*deltx(i,j))
      dfeII(i,j) = dFeII(i,j) + dt * kh(i,j)* &
                ( dFeII(i+1,j)-2*dFeII(i,j)+dfeII(i-1),j)/(deltx(i,j)*deltx(i,j))
      dfeIII(i,j) = dFeIII(i,j) + dt * kh(i,j)* &
                ( dFeIII(i+1,j)-2*dFeIII(i,j)+dFeIII(i-1),j)/(deltx(i,j)*deltx(i,j))       
!      enddo
    enddo
  enddo

  do j = sy1,ey1
    do i = sx1,ex1
!      dt0=delty(i,j)*delty(i,j)/(4.*kh(i,j))
!      iii= dt/dt0+1
!      dtt=dt/iii
!      do it=1,iii
!      c(i,j)=c(i,j)+dtt*kh(i,j)*(c(i,j+1)-2*c(i,j)+c(i,j-1))/ &
      c(i,j)=c(i,j)+dt*kh(i,j)*(c(i,j+1)-2*c(i,j)+c(i,j-1))/ &
                    (delty(i,j)*delty(i,j))
      dso4(i,j) = dso4(i,j) + dt * kh(i,j)* &
                ( dso4(i+1,j)-2*dso4(i,j)+dso4(i-1),j)/(deltx(i,j)*deltx(i,j))
      dno3(i,j) = dno3(i,j) + dt * kh(i,j)* &
                ( dno3(i+1,j)-2*dno3(i,j)+dno3(i-1),j)/(deltx(i,j)*deltx(i,j))
      dfeII(i,j) = dFeII(i,j) + dt * kh(i,j)* &
                ( dFeII(i+1,j)-2*dFeII(i,j)+dfeII(i-1),j)/(deltx(i,j)*deltx(i,j)) 
      dfeIII(i,j) = dFeIII(i,j) + dt * kh(i,j)* &
                ( dFeIII(i+1,j)-2*dFeIII(i,j)+dFeIII(i-1),j)/(deltx(i,j)*deltx(i,j))
!      enddo
    enddo
  enddo

  return
  end

