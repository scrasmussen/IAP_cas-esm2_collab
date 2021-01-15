! Zifa 2006/10/02 revised
! For Source Mark
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine adv_hori_mark( myid, c0, u, v, deltx, delty, &
                       sx, ex, sy, ey ,nx,ny,dt, &
                       ismMax,sm,ne,k,nzz,ktop )
  integer,parameter :: ismMaxSet=100
  integer myid, sx, ex, sy, ey, it,k,nzz
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c0, u, v, deltx, delty,ktop
  real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1) :: sm
  integer i, j
  integer sx1,ex1,sy1,ey1
  real smthis(ismMaxSet)
  real smother(ismMaxSet)
  integer, allocatable, dimension(:,:) :: IQ,JQ
  real, allocatable, dimension(:,:) :: c
  !real, allocatable, dimension(:) :: smthis,smother

!  if(ismMax.gt.ismMaxSet)stop 'ismMaxSet too small'

  allocate(IQ(sx-1:ex+1,sy-1:ey+1))
  allocate(JQ(sx-1:ex+1,sy-1:ey+1))
  allocate(c(sx-1:ex+1,sy-1:ey+1))

  c=c0

  sx1=sx
  ex1=ex
  sy1=sy
  ey1=ey
  do j = sy,ey+1
      do i = sx1,ex1+1
       if(u(i,j).gt.0.)then
          IQ(i,j)=i-1
       else
          IQ(i,j)=i
       endif
       if(v(i,j).gt.0.)then
          JQ(i,j)=j-1
       else
          JQ(i,j)=j
       endif
 
      enddo
  enddo

  dt0 = 1000.
  do j = sy,ey
  do i = sx,ex
   dt0 = min(dt0, deltx(i,j)/abs(u(i,j)),delty(i,j)/abs(v(i,j)))
  enddo
  enddo

  itt=dt/dt0+1
  dtt=dt/itt
  icontrol=1

!  print *,dt0,dtt,itt
do ii = 1, itt

  do j = sy1,ey1
       do i = sx1,ex1
       ! for two parts
         deltc= -dtt/deltx(i,j)*u(i+1,j)*c(IQ(i+1,j),j)
         c00=c(i,j)
         c(i,j)=c(i,j)+ deltc

      IF(k>ktop(i,j)) then
       if(is==2)  then    
         sm(is,i,j)=1.0
       else  
         sm(is,i,j)=0.0      
       endif  
      ENDIF
       if(icontrol==1)then
         if(deltc > 0 ) then
          do is=1,ismMax
           smthis(is) = sm(is,i,j)
           smother(is)= sm(is,i+1,j)   !!!!need check
          enddo
        
    imotherboundary =0
! to judge if it is mother domain
if(ne==1)then
  imotherboundary=0
  if(i==nx) imotherboundary=1
endif
  if(imotherboundary==1)then
!---------------------li jie modify the sm-------------     
     IF(k<=ktop(i,j)) then
       call GetBounChange(deltc,c00,smthis,1,ismMax)
     ELSE
       call GetBounChange(deltc,c00,smthis,2,ismMax)
     ENDIF  
!---------------------finish---------------------------
    else
         call SMmixing(c00,smthis,deltc,smother,ismMax)
    endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
          endif
      endif 
       !!!!!!!!!  
         deltc=  dtt/deltx(i,j)*u(i  ,j)*c(IQ(i,j),j) 
         c00=c(i,j)
         c(i,j)=c(i,j)+ deltc
if(icontrol==1)then
  if(deltc > 0 ) then
         do is=1,ismMax
           smthis(is) = sm(is,i,j)
           smother(is)= sm(is,i-1,j)
          enddo
  imotherboundary =0
! to judge if it is mother domain
if(ne==1)then
  imotherboundary=0
  if(i==1) imotherboundary=1
endif

 if(imotherboundary==1)then
!---------------------li jie modify the sm-------------     
     IF(k<=ktop(i,j) ) then
        call GetBounChange(deltc,c00,smthis,1,ismMax)   
     ELSE 
        call GetBounChange(deltc,c00,smthis,2,ismMax)
     ENDIF    
!---------------------finish---------------------------
 else
        call SMmixing(c00,smthis,deltc,smother,ismMax)
 endif

          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
 endif
endif

  enddo
  enddo

  do j = sy1,ey1
      do i = sx1,ex1
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      deltc=-dtt/delty(i,j)*v(i,j+1)*c(i,JQ(i,j+1))
      c00=c(i,j)
      c(i,j)=c(i,j)+deltc
     if(icontrol==1)then
        if(deltc > 0 ) then
         do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i,j+1)  
          enddo
          deltc1=deltc
  imotherboundary =0
! to judge if it is mother domain
if(ne==1)then
    imotherboundary=0
    if(j==ny) imotherboundary=1
endif
    if(imotherboundary==1)then
!---------------------li jie modify the sm-------------     
     IF(k<=ktop(i,j)) then
       call GetBounChange(deltc1,c00,smthis,1,ismMax)   
     ELSE
       call GetBounChange(deltc1,c00,smthis,2,ismMax)      
     ENDIF
!---------------------finish---------------------------
    else
          call SMmixing(c00,smthis,deltc1,smother,ismMax)
    endif

          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
          endif
       endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      deltc= dtt/delty(i,j)* v(i,j  )*c(i,JQ(i,j ))  
      c00=c(i,j)
      c(i,j)=c(i,j)+deltc
     if(icontrol==1)then
     if(deltc > 0 ) then
         do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i,j-1)        
          enddo
          deltc1=deltc
  imotherboundary =0
! to judge if it is mother domain
if(ne==1)then
    imotherboundary=0
    if(j==1) imotherboundary=1
endif
    if(imotherboundary==1)then
!---------------------li jie modify the sm-------------     
     IF(k<=ktop(i,j))THEN
       call GetBounChange(deltc1,c00,smthis,1,ismMax)   
     ELSE
       call GetBounChange(deltc1,c00,smthis,2,ismMax)      
     ENDIF
!---------------------finish---------------------------
    else
       call SMmixing(c00,smthis,deltc1,smother,ismMax)
    endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
          endif
       endif

      enddo
  enddo

enddo


  deallocate(IQ,JQ)
  deallocate(c)
!  deallocate(smthis,smother)
  return
  end


