
subroutine naqpms_dry_output &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,igas,iaer,isize,nseacom,ndustcom )

use naqpms_varlist, only : ddep_flag,ip2mem,ip2memGas,DryVelGas ! 
!use met_fields, only : RAINNON,RAINCON
use naqpms_varlist, only : iedgas

use naqpms_varlist, only : naerbin,naersp,ip2memaer,dryvelaer
use naqpms_varlist, only : ip3memaer,DryVeldust


implicit none

integer :: myid

real :: dt


integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: iyear,imonth,iday,ihour,iminute

integer :: i,j,k,is,ia

integer :: i02,i0,i02gas,i03aer,i02aer

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

  call system("mkdir -p drydep/tmp")
  call system("mkdir -p drydep/tmp/"//date(1:8))

  fname='drydep/tmp/'//date(1:8)//'/drydepd'//cdnum//'.'//cmyid//'.'//date

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')

  write(funit) trim(ddep_flag)
  do ig = 1 , iedgas
    i02gas = ip2memGas(ig,ne)
    call writedry_v2(  myid,DryVelGas(i02gas) &
                     , sx(ne), ex(ne), sy(ne), ey(ne),ig,ne,funit)
  enddo
   

  do ia=1,naersp
  do is=1,naerbin
    i02aer = ip2memaer(is,ia,ne)
    ig=(ia-1)*naerbin+is
    call writedry_v2(  myid,dryvelaer(i02aer) &
                     , sx(ne), ex(ne), sy(ne), ey(ne),ig,ne,funit)
  enddo
  enddo

  do ia=1,iaer
  do is=1,isize
    i03aer = ip3memaer(is,ia,ne)
    ig=(ia-1)*isize+is
    call writedry_v2(  myid,DryVeldust(i03aer) &
                     , sx(ne), ex(ne), sy(ne), ey(ne),ig,ne,funit)
  enddo
  enddo


   close(funit)



end subroutine naqpms_dry_output





