
subroutine naqpms_term_output &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,iprocess,PrintTermGas )

use naqpms_varlist, only : IGOPos,ipGasTermBal,GasTermBal 
implicit none

integer :: myid

real :: dt


integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: iyear,imonth,iday,ihour,iminute

integer :: iprocess
integer :: PrintTermGas(igas)

integer :: i,j,k,is,ia,iduc

integer :: i05,igo,ip

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

  call system("mkdir -p term/tmp")
  call system("mkdir -p term/tmp/"//date(1:8))

  fname='term/tmp/'//date(1:8)//'/termd'//cdnum//'.'//cmyid//'.'//date

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')

  do ig=1,igas
    if(PrintTermGas(ig)==1)then
      igo = IGOPos(ig)
      do ip=1,iprocess
      do k=1,nzz
        i05=ipGasTermBal(k,ip,igo,ne)
        call writeterm_v2( myid,GasTermBal(i05) &
                          ,sx(ne),ex(ne),sy(ne),ey(ne),k,ip,ig,ne )
      enddo
      enddo
    endif
  enddo

  close(funit)

end subroutine naqpms_term_output





