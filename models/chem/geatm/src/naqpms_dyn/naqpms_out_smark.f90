
subroutine naqpms_out_smark &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d,tropp &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,PrintGas &
 &  ,ifsm,idmSet,ismMax,igMark )

use naqpms_varlist
use met_fields
use naqpms_gridinfo
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

integer :: ig,iduc,ia,i05,i05c,ism,idm


integer :: iyear,imonth,iday,ihour,iminute

integer :: ifsm(5)

integer :: idmSet,ismMax

integer :: igMark(idmSet)



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


  call system("mkdir -p sm/tmp")
  call system("mkdir -p sm/tmp/"//date(1:8))


  fname='sm/tmp/'//date(1:8)//'/makd'//cdnum//'.'//cmyid//'.'//date

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')

  if(ifsm(ne)==1) then 

  do k=1,nzz
  i0=ip3mem(k,ne)
  call writesm_v2(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               k,1,1,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call writesm_v2(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               k,1,1,funit)
  enddo

  do idm=1,idmSet
    do ism=1,ismMax
    do k=1,nzz
    i0=ipSMmem(k,ism,idm,ne)
    call writesm_v2(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
               k,ism,idm,funit)
     enddo
     enddo
    do k=1,nzz
     i0=ip4mem(k,igMark(idm),ne)
     call writesm_v2(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               k,1,1,funit)
     enddo
  enddo

  close(funit)

  IF(ihour==0)THEN  !! to write down reinit source file
    !! to write source file

    fname='sm/tmp/'//date(1:8)//'/resd'//cdnum//'.'//cmyid//'.'//date

    call get_funitnaqpms(funit)
 
    open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')

    do idm=1,idmSet
      do ism=1,ismMax
      do k=1,nzz
        i0=ipSMmem(k,ism,idm,ne)
        call writesm_v2(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
               k,ism,idm,funit)
      enddo
      enddo
      do k=1,nzz
        i0=ip4mem(k,igMark(idm),ne)
        call writesm_v2(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               k,1,1,funit)
      enddo
    enddo

    close(funit)

  ENDIF

endif

end subroutine naqpms_out_smark





