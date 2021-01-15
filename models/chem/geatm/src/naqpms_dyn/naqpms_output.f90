
subroutine naqpms_tracer_output &
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


!return

  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)') ne
  write (cmyid,'(i3.3)') myid


  call system("mkdir -p out/tmp")
  call system("mkdir -p out/tmp/"//date(1:8))

!return

  fname='out/tmp/'//date(1:8)//'/food'//cdnum//'.'//cmyid//'.'//date

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')

  do k=1,nzz
  i0=ip3mem(k,ne)
  call write2d_v2(myid,dz(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
            
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


  do ig=1,iedgas
    if(PrintGas(ig)==1)then
     do k=1,nzz
     i0=ip4mem(k,ig,ne)
     call write2d_v2(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
     enddo
    endif
  enddo


! write(funit) naersp,naerbin
 do ia=1,naersp
 do is=1,naerbin
!   write(funit) ia,is
   do k=1,nzz
     i0=ip4mem_aer(k,is,ia,ne)
     call write2d_v2(myid,aerom(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
   enddo
 enddo
 enddo

! write(funit) 'edaer'

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


!if(1==2) then
  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,entrn3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,dsdt3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
!endif

if(lgaschemsmp) then
  do ig=1,nmoxdt
  do k=1,nzz
     i0=ip4mem_ox3d(k,ig,ne)
     call write2d_v2(myid,oxdt3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo
endif

   i0=ip_emit2d(57,1,ne)
   call write2d_v2(myid,emit2d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)


!  do k=1,nzz
!   i0=ip3mem(k,ne)
!   call write2d_v2(myid,tmp4d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
!  enddo



!
  close(funit)

end subroutine naqpms_tracer_output





