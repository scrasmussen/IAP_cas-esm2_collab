
subroutine rd_global_conc( myid,iyear,imonth,iday,ihour,iminute &
                  ,nest,nzz,nx,ny,nz,sx,ex,sy,ey,ne )

use naqpms_varlist, only : ip2mem,ip3mem,globalno2,globalo3,globalco

implicit none 

integer :: myid
integer :: ne,nest
integer :: nx(5),ny(5),nz(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)


integer :: iyear,imonth,iday,ihour,iminute


integer :: k,i0

character :: cdnum*1
character :: date*10

integer :: irec
integer :: funit
logical :: lexist

character :: fname*100


  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)') ne
 
  call get_funitnaqpms(funit)

  fname='global/g2m/g2md0'//cdnum//'_'//date(1:4)//'-'//date(5:6)//&
       '-'//date(7:8)// '_'//date(9:10)//'.DAT'


  inquire(file=fname,exist=lexist)

  if(.not.lexist) then
    print*,trim(fname)//' NOT exist'
    stop
  endif

 open(funit,file=trim(fname),form='unformatted',access='direct',recl=nx(ne)*ny(ne))
 
 irec=1

  ! to read global concentration of no2
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,globalno2(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

  ! to read global concentration of o3
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,globalo3(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

  ! to read global concentration of co
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,globalco(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo


close(funit)

end subroutine rd_global_conc

