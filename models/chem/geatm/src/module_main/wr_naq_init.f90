
subroutine wr_naq_init( myid,ntbeg &
                  ,iyear1,imonth1,iday1,ihour1,iminute1 &
                  ,nest,nzz,nx,ny,nz,sx,ex,sy,ey,ne &
                  ,igas,isize,iaer &
                  ,ndustcom &
                  ,ifsm,ifsmt,idmSet,ismMax,igMark )


use naqpms_varlist
use met_fields
use naqpms_gridinfo

use  work_vars, only: tropp

implicit none

integer :: myid
integer :: ne,nest
integer :: nx(5),ny(5),nz(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: ifsm(5)

integer :: ifsmt

integer :: idmSet,ismMax

integer :: igMark(idmSet)

integer :: ntbeg(5)

integer :: igas,isize,iaer

integer :: ndustcom


integer :: iyear1,imonth1,iday1,ihour1,iminute1

integer :: iyear,imonth,iday,ihour,iminute
integer :: iitime

integer :: k,i0,ig,ia,is,idug,iduc

integer :: idm,ism

character :: cdnum*1
character :: date*10

character :: cmyid*3

integer :: irec
integer :: funit

logical :: lexist

character :: fname*100





!============  Zifa to add ()write initial condition =====
!===========   2006-07-07 
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww


 
loop_nest : do ne=1,nest

  iitime = (ntbeg(ne)-1)*3600    ! in seconds
  call getnewdate(iyear1,imonth1,iday1,ihour1,iitime, & 
                iyear,imonth,iday, ihour,iminute)


  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  call system("mkdir -p out/tmp")
  call system("mkdir -p out/tmp/"//date(1:8))

  write(cdnum(1:1),'(i1)') ne
  write (cmyid,'(i3.3)') myid

  fname='out/tmp/'//date(1:8)//'/food'//cdnum//'.'//cmyid//'.'//date  

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')


  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,w(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,t(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,rh1(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,Plev(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip2mem(ne)
   call write2d_v2(myid,tropp(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  
  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,heiz(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo


  
  do ig=1,igas
!    if(PrintGas(ig)==1)then
     do k=1,nzz
     i0=ip4mem(k,ig,ne)
     call write2d_v2(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
     enddo
!    endif
  enddo

 do ia=1,naersp
 do is=1,naerbin
   do k=1,nzz
     i0=ip4mem_aer(k,is,ia,ne)
     call write2d_v2(myid,aerom(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
   enddo
 enddo
 enddo


 do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,jo1d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,jno2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,EXT(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,EXTASO4(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,EXTANO3(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,EXTANH4(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,EXTBC(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,EXTOC(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,VISIB(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,UVB(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,UVBS(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,UVA(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,VIS(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,SSA(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  i0=ip2mem(ne)
  call write2d_v2(myid,AOD(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,CLDOPD(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,DUSO2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,DUO3(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,DUNO2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)


  do k=1,nzz
   I0=ip3mem(k,ne)
   call write2d_v2(MYID,ANA(I0),SX(NE),EX(NE),SY(NE),EY(NE),k,ne,funit)
  ENDDO

  do k=1,nzz
   I0=ip3mem(k,ne)
   call write2d_v2(MYID,ASO4(I0),SX(NE),EX(NE),SY(NE),EY(NE),k,ne,funit)
  ENDDO


  do k=1,nzz
   I0=ip3mem(k,ne)
   call write2d_v2(MYID,ANH4(I0),SX(NE),EX(NE),SY(NE),EY(NE),k,ne,funit)
  ENDDO

  do k=1,nzz
   I0=ip3mem(k,ne)
   call write2d_v2(MYID,ANO3(I0),SX(NE),EX(NE),SY(NE),EY(NE),k,ne,funit)
  ENDDO

  do k=1,nzz
   I0=ip3mem(k,ne)
   call write2d_v2(MYID,ACL(I0),SX(NE),EX(NE),SY(NE),EY(NE),k,ne,funit)
  ENDDO


  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,CPH(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,OPE(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo


  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  i0=ip5mem(k,is,ia,ne)
  call write2d_v2(myid,aer(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,coefcld3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,gscav_so2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo


  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,ascav_so4(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,entrn3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,dsdt3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo



  close(funit)   !Zifa 2006/07/07


 !!!!!!!!!!!!!!!!!!!!!!!!!!
 ! for Source Mark
  if(ifsm(ne)==1)then !!!!!!!!!!!!!

  call system("mkdir -p sm")
  call system("mkdir -p sm/tmp")
  call system("mkdir -p sm/tmp/"//date(1:8))

  fname='sm/tmp/'//date(1:8)//'/makd'//cdnum//'.'//cmyid//'.'//date

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')

  !! to write source file <only save 3 levels> 1,3,5
   
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


!  IF(ihour2==0)THEN  !! to write down reinit source file
  !! to write source file
!  do idm=1,idmSet
!    do ism=1,ismMax
!    do k=1,nzz
!    i0=ipSMmem(k,ism,idm,ne)
!    call writesm_v2(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
!               k,ism,idm,1funit)
!     enddo
!     enddo
!    do k=1,nzz
!     i0=ip4mem(k,igMark(idm),ne)
!     call writesm_v2(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
!               k,1,1,1funit)
!     enddo
!  enddo
!  close(1funit)
!  ENDIF

  endif ! if(ifsm(ne)==1)


enddo loop_nest  ! write for different domain

end subroutine wr_naq_init

