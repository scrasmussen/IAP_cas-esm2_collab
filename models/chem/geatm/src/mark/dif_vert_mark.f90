!
  subroutine dif_vert_mark( myid, c10,c0,cp10,kz,deltz,  &
                       sx, ex, sy, ey ,k, nzz,dt,ismMax, &
             sm,smb,smp)
  integer,parameter :: ismMaxSet=100
  integer myid, sx, ex, sy, ey, it
  integer i, j, k
  real smthis(ismMaxSet)
  real smother(ismMaxSet)

  real, dimension(sx-1:ex+1,sy-1:ey+1) :: kz
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c10,c0,cp10
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: deltz
  real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1):: sm,smb,smp

  real, allocatable, dimension(:,:) ::c1,c,cp1

  allocate(c1(sx-1:ex+1,sy-1:ey+1))
  allocate(c(sx-1:ex+1,sy-1:ey+1))
  allocate(cp1(sx-1:ex+1,sy-1:ey+1))


  c1=c10
  c=c0
  cp1=cp10

 ! goto 333

  do j = sy,ey
    do i = sx,ex
      deltc1=dt*kz(i,j)*(cp1(i,j)-c(i,j))/(deltz(i,j)*deltz(i,j))  !Zifa wang original
      c00=c(i,j)

      if(deltc1>0)then
        do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=smp(is,i,j)
          enddo
        if(k==nzz-1)then
          call GetBounChange(deltc1,c00,smthis,2,ismMax)
        else if(k>nzz-1) then
          smthis(2)=1.0
        else 
          call SMmixing(c00,smthis,deltc1,smother,ismMax)
        endif

      do is=1,ismMax
           sm(is,i,j)=smthis(is)
       enddo 

      endif

      deltc2=dt*kz(i,j)*(-c(i,j)+c1(i,j))/(deltz(i,j)*deltz(i,j))   !Zifa wang original
    if(deltc2>0)then
        do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=smb(is,i,j)
          enddo
        
        if(k>nzz-1) then
          smthis(2)=1.0
        else
          call SMmixing(c00,smthis,deltc2,smother,ismMax)
        endif   
       do is=1,ismMax
           sm(is,i,j)=smthis(is)
       enddo 
     endif

       c(i,j)=c(i,j)+deltc1+deltc2
    enddo
  enddo

 333 continue
  deallocate(c1,c,cp1)
  return
  end

