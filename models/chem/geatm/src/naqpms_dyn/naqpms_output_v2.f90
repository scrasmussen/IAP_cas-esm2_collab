
subroutine naqpms_tracer_output_v2 &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d,tropp &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,PrintGas)

use naqpms_varlist
use met_fields
use naqpms_gridinfo

use smpsulf_var, only : ip4mem_ox3d,oxdt3d,nmoxdt
use smpsulf_var, only : tmp4d

implicit none

integer :: myid

real :: dt


integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: i,j,k,is
!integer :: k,is

integer :: mem3d

integer :: mem2d
real,dimension(mem2d) :: tropp

integer,dimension(igas) :: PrintGas

integer :: ixy,i02,i03,i04,i0

integer :: ig,iduc,ia,i05,i05c


integer :: iyear,imonth,iday,ihour,iminute


character :: cdnum*1
character :: date*10

character :: cmyid*3

integer :: irec
integer :: funit

logical :: lexist

character :: fname*100

integer,save :: iyearb,imonthb,idayb,ihourb,iminuteb
logical,save :: lfirst=.true.
integer,save :: nttadd

integer :: i22,nrec
integer,save :: nrec22

character :: coutfrclab*20

logical :: loutput

!return

!coutfrclab='daily'

coutfrclab=caveoutclab

if(lfirst) then
iyearb=iyear
imonthb=imonth
idayb=iday
ihourb=ihour
iminuteb=iminute
lfirst=.false.
nttadd=0
endif

loutput=.true.
if(trim(coutfrclab).eq.'monthly'.and.imonth.eq.imonthb) then
  loutput=.false.
elseif(trim(coutfrclab).eq.'daily'.and.iday.eq.idayb) then
  loutput=.false.
elseif(trim(coutfrclab).eq.'hourly'.and.ihour.eq.ihourb) then
  loutput=.false.
  return
endif


if(.not.loutput) then

 call naqpms_aveload &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d,tropp &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,PrintGas &
 &  ,nrec )

 nttadd=nttadd+1
 nrec22=nrec

else ! output

  write(date(1:4),'(i4)')iyearb
  write(date(5:6),'(i2.2)')imonthb
  write(date(7:8),'(i2.2)')idayb
  write(date(9:10),'(i2.2)')ihourb

  write(cdnum(1:1),'(i1)') ne
  write (cmyid,'(i3.3)') myid

  if(trim(coutfrclab).eq.'monthly') then
    call system("mkdir -p out/tmp/"//date(1:4))
    fname='out/tmp/'//date(1:4)//'/food'//cdnum//'.'//cmyid//'.'//date(1:6)
  elseif(trim(coutfrclab).eq.'daily') then
    call system("mkdir -p out/tmp/"//date(1:6))
    fname='out/tmp/'//date(1:6)//'/food'//cdnum//'.'//cmyid//'.'//date(1:8)
  elseif(trim(coutfrclab).eq.'hourly') then
    call system("mkdir -p out/tmp/"//date(1:8))
    fname='out/tmp/'//date(1:8)//'/food'//cdnum//'.'//cmyid//'.'//date(1:10)
  else
    stop 'coutfrclab set erro'
  endif

  if(myid.eq.0) print*,'nttadd=',nttadd,'nrec22=',nrec22

  avgvar=avgvar/nttadd

  call get_funitnaqpms(funit)
  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')
  do irec=1,nrec22
     i22=ip4aveout(irec,ne)
     call write2d_v2(myid,avgvar(i22),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  close(funit)

  avgvar=0.0
  call naqpms_aveload &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d,tropp &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,PrintGas &
 &  ,nrec )

  nttadd=1
  nrec22=nrec

endif


iyearb=iyear
imonthb=imonth
idayb=iday
ihourb=ihour
iminuteb=iminute

end subroutine naqpms_tracer_output_v2





