
subroutine naqpms_wet_output &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,igas,iaer,isize,nseacom,ndustcom )

use naqpms_varlist, only : ip2mem,ip2memGas,WETDEP,WETDEP2 
use naqpms_varlist, only : ip2memaer,aerom_wdepfld,aerom_wdepfld2
use naqpms_varlist, only : naersp,naerbin
use naqpms_varlist, only : iedgas
use met_fields, only : RAINNON,RAINCON
implicit none

integer :: myid

real :: dt


integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: iyear,imonth,iday,ihour,iminute

integer :: i,j,k,is,ia

integer :: i02,i0,i02gas,i02aer,ixy

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

  call system("mkdir -p wetdep/tmp")
  call system("mkdir -p wetdep/tmp/"//date(1:8))

  fname='wetdep/tmp/'//date(1:8)//'/wetdepd'//cdnum//'.'//cmyid//'.'//date

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')



   ! to write CUMULUS PRECIPITATION (cm/h)
   i0=ip2mem(ne)
   call writewet_v2(myid,RAINCON(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

   ! to write NON-CUMULUS PRECIPITATION (cm/h) 
   i0=ip2mem(ne)
   call writewet_v2(myid,RAINNON(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
   
   do ig = 1 , iedgas
     i02gas = ip2memGas(ig,ne)
     call writewet_v2(myid,WETDEP(i02gas),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
   enddo
   
   do ig = 1 , iedgas
     i02gas = ip2memGas(ig,ne)
     call writewet_v2(myid,WETDEP2(i02gas),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
   enddo

   do ia=1,naersp
   do is=1,naerbin
     i02aer=ip2memaer(is,ia,ne)
     call writewet_v2(myid,aerom_wdepfld(i02aer),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
   enddo
   enddo

   do ia=1,naersp
   do is=1,naerbin
     i02aer=ip2memaer(is,ia,ne)
     call writewet_v2(myid,aerom_wdepfld(i02aer),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
   enddo
   enddo

   close(funit)

! clean up

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
      ixy=(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

      do ig = 1 , igas
        i02gas = ip2memGas(ig,ne)
        WETDEP(i02gas+ixy)=0.0
        WETDEP2(i02gas+ixy)=0.0
      enddo

      do ia=1,naersp
      do is=1,naerbin
        i02aer=ip2memaer(is,ia,ne)
        aerom_wdepfld(i02aer+ixy)=0.0
        aerom_wdepfld2(i02aer+ixy)=0.0
      enddo
      enddo
   enddo
   enddo

end subroutine naqpms_wet_output





