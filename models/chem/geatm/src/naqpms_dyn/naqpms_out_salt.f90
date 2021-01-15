
subroutine naqpms_salt_output &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,igas,iaer,isize,nseacom,ndustcom )

use naqpms_varlist !, only : ip2mem,ip2memGas,WETDEP,WETDEP2 
use met_fields, only : u,v
implicit none

integer :: myid

real :: dt


integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: iyear,imonth,iday,ihour,iminute

integer :: i,j,k,is,ia,iduc

integer :: i02,i0,i02gas

integer :: ig

character :: cdnum*1
character :: date*10

character :: cmyid*3

integer :: irec
integer :: funit

logical :: lexist

character :: fname*100

  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)') ne
  write (cmyid,'(i3.3)') myid

  call system("mkdir -p seasalt/tmp")
  call system("mkdir -p seasalt/tmp/"//date(1:8))

  fname='seasalt/tmp/'//date(1:8)//'/seasaltd'//cdnum//'.'//cmyid//'.'//date

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')

  do k=1,nzz
  i0=ip3mem(k,ne)
  call writedust_v2(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call writedust_v2(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do ia=1,1
  do is=1,isize
  do k=1,nzz
  i0=ip5mem(k,is,ia,ne)
  call writedust_v2(myid,aer(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo
  enddo

  do iduc = 1, nseacom
  do is = 1,isize
  do k =1,nzz
   i0 = ip5memcs(k,is,iduc,ne)
   call  writedust_v2(myid,seacomp(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo
  enddo

  i0 = ip2mem(ne)
  call writedust_v2(myid,SEAEMISS(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)


  close(funit)

end subroutine naqpms_salt_output





