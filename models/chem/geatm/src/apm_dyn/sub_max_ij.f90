

subroutine max_ij_location(myid,sx,ex,sy,ey,im,jm,var,vm)
implicit none
integer :: myid
integer :: sx,ex,sy,ey,im,jm,i,j
real    :: var(sx:ex,sy:ey)
real    :: vm
vm=-1.0e6
do j=sy,ey
do i=sx,ex
 if(var(i,j).ge.vm) then
   vm=var(i,j)
   im=i
   jm=j
 endif
enddo
enddo
end subroutine max_ij_location


subroutine min_ij_location(myid,sx,ex,sy,ey,im,jm,var,vm)
implicit none
integer :: myid
integer :: sx,ex,sy,ey,im,jm,i,j
real    :: var(sx:ex,sy:ey)
real    :: vm
vm=1.0e6
do j=sy,ey
do i=sx,ex
 if(var(i,j).le.vm) then
   vm=var(i,j)
   im=i
   jm=j
 endif
enddo
enddo
end subroutine min_ij_location



