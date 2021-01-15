!!!!!!!!!!!!!!!!!!
!for Source Mark
!Zifa 2006/10/06
! to mention
! k+1   wp1 smp   
! k     w   sm
! k-1   w1  smb
!
  subroutine adv_vert_mark( myid, c10,c0,cp10,w1,w,wp1,deltz1,deltz,deltzp1,&
                       sx, ex, sy, ey ,k, nzz,dt,ismMax,&
                       sm,smb,smp,ktop)
  integer, parameter :: ismMaxSet=100
  real smthis(ismMaxSet)
  real smother(ismMaxSet)

  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: w1,w,wp1,ktop
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c10,c0,cp10
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: deltz,deltz1,deltzp1
  real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1) :: sm,smb,smp
  real :: wtmp,deltztmp
  real, allocatable, dimension(:,:) ::c1,c,cp1

   allocate(c1(sx-1:ex+1,sy-1:ey+1))
   allocate(c(sx-1:ex+1,sy-1:ey+1))
   allocate(cp1(sx-1:ex+1,sy-1:ey+1))

   c1=c10
   c=c0
   cp1=cp10

  sx1=sx
  ex1=ex
  sy1=sy
  ey1=ey
  if(ex1==nx)ex1=nx-1
  if(ey1==ny)ey1=ny-1

!  dt0 = 1000.

  do j = sy1,ey1
  do i = sx1,ex1
!   dt0 = min(dt0, deltz(i,j)/abs(w(i,j)))
    IF(k>ktop(i,j)) then
         if(is==2) then
             sm(is,i,j)=1.0
         else
             sm(is,i,j)=0.0
         endif
     ENDIF     
  enddo
  enddo

!  itt=dt/dt0+1
!  dtt=dt/itt
    dtt=dt
!do ii = 1, itt
      do j = sy1,ey1
        do i = sx1,ex1
           wtmp=(w(i,j)+w(i,j))/2.0
           deltztmp=(deltz(i,j)+deltz(i,j))/2.0
           if(wtmp.ge.0)then
           deltc1=dtt/deltztmp*(-wtmp*c(i,j))
           else
           deltc1=dtt/deltztmp*(-wtmp*cp1(i,j))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
       c00=c(i,j)
       do is=1,ismMax
          smthis(is)=sm(is,i,j)
          smother(is)=smp(is,i,j)
        enddo
        if(k==ktop(i,j))then
          call GetBounChange(deltc1,c00,smthis,2,ismMax)
        else if(k>ktop(i,j)) then
          smthis(2)=1.
        else if(k<ktop(i,j)) then 
          call SMmixing(c00,smthis,deltc1,smother,ismMax)
        endif
        do is=1,ismMax
           sm(is,i,j)=smthis(is)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        endif

           c(i,j)=deltc1+c(i,j)

           enddo
        enddo
!enddo

if(k>1)then

!do ii = 1, itt
      do j = sy1,ey1
        do i = sx1,ex1
            wtmp=(w1(i,j)+w1(i,j))/2.0
            deltztmp=(deltz(i,j)+deltz(i,j))/2.0
           if(wtmp.ge.0)then
           deltc2=dtt/deltztmp*wtmp*c1(i,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       c00=c(i,j)
       do is=1,ismMax           
          smthis(is)=sm(is,i,j) 
          smother(is)=smb(is,i,j) 
        enddo  
        if(k<=ktop(i,j))then
            call SMmixing(c00,smthis,deltc2,smother,ismMax)
        else
            smthis(2)=1.0
        endif
        do is=1,ismMax
           sm(is,i,j)=smthis(is)
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           else
           deltc2=dtt/deltz(i,j)*w(i,j)*c(i,j)
           endif

           c(i,j)=c(i,j)+deltc2
           enddo
        enddo
 !enddo

 endif

  deallocate(c1,c,cp1)
  return
  end

