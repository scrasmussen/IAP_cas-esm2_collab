
subroutine naqpms_aveload &
 & ( myid &
 &  ,iyear,imonth,iday,ihour,iminute &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d,tropp &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,PrintGas &
 &  ,nrec )

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

integer :: i22,irec,nrec



  irec=0
  do k=1,nzz
    irec=irec+1
    i0=ip3mem(k,ne)
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),dz(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo
            
  do k=1,nzz
    irec=irec+1
    i0=ip3mem(k,ne)
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),u(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo
  
  do k=1,nzz
    irec=irec+1
    i0=ip3mem(k,ne)
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),v(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
    irec=irec+1
    i0=ip3mem(k,ne)
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),w(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
    irec=irec+1
    i0=ip3mem(k,ne)
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),t(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
    irec=irec+1
    i0=ip3mem(k,ne)
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),rh1(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
    irec=irec+1
    i0=ip3mem(k,ne)
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),Plev(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo


!  do k=1,nzz
    irec=irec+1
!    i0=ip3mem(k,ne) ! juanxiong he
    i0=ip2mem(ne) ! juanxiong he
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),tropp(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
!  enddo

  do k=1,nzz
    irec=irec+1
    i0=ip3mem(k,ne)
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),heiz(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo


  do ig=1,iedgas
    if(PrintGas(ig)==1)then
     do k=1,nzz
     i0=ip4mem(k,ig,ne)
     irec=irec+1
     i22=ip4aveout(irec,ne)
     call load_aveout(myid,avgvar(i22),gas(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
     enddo
    endif
  enddo


! write(funit) naersp,naerbin
 do ia=1,naersp
 do is=1,naerbin
!   write(funit) ia,is
   do k=1,nzz
     i0=ip4mem_aer(k,is,ia,ne)
     irec=irec+1
     i22=ip4aveout(irec,ne)
     call load_aveout(myid,avgvar(i22),aerom(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
   enddo
 enddo
 enddo

! write(funit) 'edaer'

 do k=1,nzz
   i0=ip3mem(k,ne)
     irec=irec+1
     i22=ip4aveout(irec,ne)
     call load_aveout(myid,avgvar(i22),jo1d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
     irec=irec+1
     i22=ip4aveout(irec,ne)
     call load_aveout(myid,avgvar(i22),jno2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
     irec=irec+1
     i22=ip4aveout(irec,ne)
     call load_aveout(myid,avgvar(i22),EXT(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
     irec=irec+1
     i22=ip4aveout(irec,ne)
     call load_aveout(myid,avgvar(i22),EXTASO4(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
    i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),EXTANO3(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),EXTANH4(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo



  do k=1,nzz
   i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),EXTBC(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),EXTOC(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),VISIB(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),UVB(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),UVBS(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),UVA(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo


  do k=1,nzz
   i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),VIS(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
    irec=irec+1
    i22=ip4aveout(irec,ne)
    call load_aveout(myid,avgvar(i22),SSA(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  i0=ip2mem(ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),AOD(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)

  i0=ip2mem(ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),CLDOPD(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)

  i0=ip2mem(ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),DUSO2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)

  i0=ip2mem(ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),DUO3(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)

  i0=ip2mem(ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),DUNO2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)


  do k=1,nzz
   I0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),ANA(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  ENDDO

  do k=1,nzz
   I0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),ASO4(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  ENDDO


  do k=1,nzz
   I0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),ANH4(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  ENDDO

  do k=1,nzz
   I0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),ANO3(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  ENDDO

  do k=1,nzz
   I0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),ACL(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  ENDDO


  do k=1,nzz
   i0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),CPH(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo


  do k=1,nzz
   i0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),OPE(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo


  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  i0=ip5mem(k,is,ia,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),aer(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo
  enddo
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),coefcld3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),gscav_so2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo


  do k=1,nzz
   i0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),ascav_so4(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo


!if(1==2) then
  do k=1,nzz
   i0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),entrn3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),dsdt3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo
!endif

if(lgaschemsmp) then
  do ig=1,nmoxdt
  do k=1,nzz
     i0=ip4mem_ox3d(k,ig,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),oxdt3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)
  enddo
  enddo
endif

   i0=ip_emit2d(57,1,ne)
  irec=irec+1
  i22=ip4aveout(irec,ne)
  call load_aveout(myid,avgvar(i22),emit2d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne)

  nrec=irec


end subroutine naqpms_aveload




subroutine load_aveout( myid, s, a, sx, ex, sy, ey, k, id)
  integer myid, sx, ex, sy, ey, k,funit
  real a(sx-1:ex+1,sy-1:ey+1),s(sx-1:ex+1,sy-1:ey+1)
  s=s+a
  return
end
