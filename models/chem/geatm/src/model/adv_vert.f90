!
!
!
  subroutine adv_vert( myid, c1,c,cp1,w1,w,wp1,deltz1,deltz,deltzp1, &
                       sx, ex, sy, ey ,k, dt,ktop)
  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: w1,w,wp1,ktop
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1,c,cp1
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: deltz,deltz1,deltzp1
  real :: wtmp,deltztmp
  integer i, j, k

   
  sx1=sx
  ex1=ex
  sy1=sy
  ey1=ey
  if(ex1==nx)ex1=nx-1
  if(ey1==ny)ey1=ny-1

   dtt=dt
      do j = sy1,ey1
        do i = sx1,ex1
 IF(k<=ktop(i,j)) THEN
           wtmp=(wp1(i,j)+wp1(i,j))/2.0
!           deltztmp=(deltz(i,j)+deltz(i,j))/2.0
           if(wtmp.ge.0)then
            deltztmp=(deltz(i,j)+deltz(i,j))/2.       
           c(i,j)=c(i,j)+dtt/deltztmp*(-wtmp*c(i,j))
           else
            deltztmp=(deltz(i,j)+deltz(i,j))/2.0        
           c(i,j)=c(i,j)+dtt/deltztmp*(-wtmp*cp1(i,j))
           endif
  ENDIF
           enddo
        enddo

 if(k==1)return

      do j = sy1,ey1
        do i = sx1,ex1
 IF(k<=ktop(i,j)) THEN
           wtmp=(w(i,j)+w(i,j))/2.0
!           deltztmp=(deltz(i,j)+deltz(i,j))/2.0
           if(wtmp.ge.0)then
            deltztmp=(deltz(i,j)+deltz(i,j))/2.0       
           c(i,j)=c(i,j)+dtt/deltztmp*wtmp*c1(i,j)
           else
            deltztmp=(deltz(i,j)+deltz(i,j))/2.0       
           c(i,j)=c(i,j)+dtt/deltztmp*wtmp*c(i,j)
           endif
 ENDIF
           enddo
        enddo


  return
  end

