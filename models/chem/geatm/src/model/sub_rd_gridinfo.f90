

subroutine read_gridinfo( myid,nest,sx,ex,sy,ey,nx,ny,nzz )

use naqpms_varlist, only : ip2mem,ip3mem
use naqpms_gridinfo
implicit none

real,parameter :: hztop=20000.0 ! 20 km 

integer :: myid
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

!integer :: ip2mem(nest),ip3mem(nzz,nest)

integer               :: mem2d
!real,dimension(mem2d) :: HGT1,TERRAIN,LATITCRS,LONGICRS,LAND_USE

integer               :: mem3d
!real,dimension(mem3d) :: dx,dy,dz,heiz

integer :: funit

integer :: idom,i,j,k,irec,i02,i03,ixy

character :: cdom*2


real,allocatable,dimension(:,:) :: dat2d


call get_funitnaqpms(funit)


do idom=1,nest

 write(cdom,'(i2.2)') idom
 ! juanxiong he
 open(funit,file='wrfd'//cdom//'.geatm.dat',form='unformatted' &
      ,access='direct',recl=nx(idom)*ny(idom))

 if(myid.eq.0) then
   print*,nest,'read gridinfo d'//cdom
 endif

 irec=0

 allocate(dat2d(nx(idom),ny(idom)))


 do k=1,nzz
   irec=irec+1
   read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))
   i03=ip3mem(k,idom)
   do j=sy(idom),ey(idom)
   do i=sx(idom),ex(idom)
     ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
     dx(i03+ixy)=dat2d(i,j)
   enddo
   enddo
 enddo


 do k=1,nzz
   irec=irec+1
   read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))
   i03=ip3mem(k,idom)
   do j=sy(idom),ey(idom)
   do i=sx(idom),ex(idom)
     ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
     dy(i03+ixy)=dat2d(i,j)
   enddo
   enddo
 enddo


 do k=1,nzz
   irec=irec+1
   read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))
   i03=ip3mem(k,idom)
   do j=sy(idom),ey(idom)
   do i=sx(idom),ex(idom)
     ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
     dz(i03+ixy)=dat2d(i,j)
   enddo
   enddo
 enddo

 do k=1,nzz
   irec=irec+1
   read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))
   i03=ip3mem(k,idom)
   do j=sy(idom),ey(idom)
   do i=sx(idom),ex(idom)
     ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
     heiz(i03+ixy)=dat2d(i,j)
   enddo
   enddo
 enddo


 i02=ip2mem(idom)

 irec=irec+1
 read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))

 do j=sy(idom),ey(idom)
 do i=sx(idom),ex(idom)
   ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
   TERRAIN(i02+ixy)=dat2d(i,j)
   HGT1(i02+ixy)=TERRAIN(i02+ixy)
 enddo
 enddo

 irec=irec+1
 read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))

 do j=sy(idom),ey(idom)
 do i=sx(idom),ex(idom)
   ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
   LATITCRS(i02+ixy)=dat2d(i,j)
 enddo
 enddo

 irec=irec+1
 read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))

 do j=sy(idom),ey(idom)
 do i=sx(idom),ex(idom)
   ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
   LONGICRS(i02+ixy)=dat2d(i,j)
 enddo
 enddo

 irec=irec+1
 read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))

 do j=sy(idom),ey(idom)
 do i=sx(idom),ex(idom)
   ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
   LAND_USE(i02+ixy)=dat2d(i,j)
 enddo
 enddo

 dat2d=1.0
 do j=sy(idom),ey(idom)
 do i=sx(idom),ex(idom)
   ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
   mpfac(i02+ixy)=dat2d(i,j)
   do k=1,nzz
      i03=ip3mem(k,idom)
      pzps(i03+ixy)=hztop-TERRAIN(i02+ixy)
   enddo
 enddo
 enddo


 deallocate(dat2d)

 if(myid.eq.0) print*,'irec=',irec

enddo


end



