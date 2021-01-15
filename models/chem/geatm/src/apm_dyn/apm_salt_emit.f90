
subroutine apm_salt_emit &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,iwb,ieb,jsb,jeb &
 &  ,land,u10,v10 &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )

 use apm_varlist
 implicit none
 include 'apm_parm.inc'

 integer :: myid

 real :: dt

 logical :: lapm

 integer :: ne,nest
 integer :: nx(5),ny(5),nzz
 integer :: sy(5),ey(5),sx(5),ex(5)
 integer :: iwb,ieb,jsb,jeb

 integer :: i,j,k,is

 integer :: mem2d
 real,dimension(mem2d)    :: u10,v10,land

 real,parameter :: pidx=3.41
 real :: u10_1d,v10_1d,ws10_1d
 
 integer :: mem3d

 integer :: ixy,i02,i03,iapm

 integer :: ip2mem(nzz)
 integer :: ip3mem(nzz,nest)


 i02=ip2mem(ne)


 loop_salt : do is=1,NSEA

   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

     if(int(land(i02+ixy)).eq.16) then ! USGS water body
       do k=1,1 ! only first layer
         iapm=ip_salt(k,is,ne)
         ws10_1d=sqrt(u10(i02+ixy)**2+v10(i02+ixy)**2)
         salt_emit(iapm+ixy)=salt_flux(is)*(ws10_1d/9.0)**pidx*kg2ug ! ug/(m2*s)
       enddo
     else
       do k=1,1
         iapm=ip_salt(k,is,ne)
         salt_emit(iapm+ixy)=0.0d0
       enddo
     endif
   enddo
   enddo

 enddo loop_salt


end subroutine apm_salt_emit



