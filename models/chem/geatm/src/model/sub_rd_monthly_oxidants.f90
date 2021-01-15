

subroutine read_monthly_oxidants( myid,nest,sx,ex,sy,ey,nx,ny,nzz )

use naqpms_varlist, only : ip2mem,ip3mem
use smpsulf_var,    only : ip4mem_oxdt,mmean_oxdt,nmoxdt,nm12
!use naqpms_gridinfo
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

integer :: idom,i,j,k,irec,i02,i03,ixy,i04
integer :: is,imon

character :: cdom*2,fname*200
logical   :: lexist

real,allocatable,dimension(:,:) :: dat2d

!return

call get_funitnaqpms(funit)


do idom=1,nest

 write(cdom,'(i2.2)') idom
 ! juanxiong he
 fname='oxdt_d'//cdom//'.dat'

!if(1==2) then
  inquire(file=fname,exist=lexist)
  if(.not.lexist) then
    print*,trim(fname)//' NOT exist'
    stop
  endif

 open(funit,file=trim(fname),form='unformatted' &
      ,access='direct',action='read',recl=nx(idom)*ny(idom))

 if(myid.eq.0) then
   print*,nest,'read oxidants d'//cdom
 endif

!endif

 irec=0

 allocate(dat2d(nx(idom),ny(idom)))

!dat2d=1.0e-2

 do imon=1,nm12

 do is=1,nmoxdt
  do k=1,nzz
   irec=irec+1
   read(funit,rec=irec) ((dat2d(i,j),i=1,nx(idom)),j=1,ny(idom))
   i04=ip4mem_oxdt(k,is,imon,idom)
   do j=sy(idom),ey(idom)
   do i=sx(idom),ex(idom)
     ixy=(ex(idom)-sx(idom)+3)*(j-sy(idom)+1)+i-sx(idom)+1
     mmean_oxdt(i04+ixy)=dat2d(i,j)*1.0e9 ! unit: vmr -> ppbv
   enddo
   enddo
  enddo
 enddo

 enddo

 deallocate(dat2d)

 if(myid.eq.0) print*,'irec=',irec

enddo


!stop

end subroutine read_monthly_oxidants



