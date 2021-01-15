

subroutine cal_cldeffect_coeff &
 & ( myid &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dx,dy,dz &
 &  ,iyear2,imonth2,iday2, ihour2,iminute2 &
 &  ,CLW,RNW,temp,Plev &
 &  ,LATITCRS,LONGICRS &
 &  ,fcld
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )

implicit none

integer :: myid

real    :: dt

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: iyear2,imonth2,iday2, ihour2,iminute2

integer :: i,j,k,is,imode

integer :: mem2d,mem3d

real,dimension(mem3d) :: LATITCRS,LONGICRS

real,dimension(mem3d) :: dx,dy,dz

real,dimension(mem3d) :: CLW,RNW,temp,Plev

integer :: ixy,i02,i03,iapm,i02apm

integer :: ip3mem(nzz,nest),ip2mem(nest)

! output

real,dimension(mem3d) :: fcld


!===============================================================
! local vars
real,dimension(nzz) :: clwc,tcol,pcol,dzcol,accld
integer :: jlnd
real :: alpha

alpha=1.0

loop_j : do j = sy(ne),ey(ne)
loop_i : do i = sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

  DO k = 1, nzz
      i03 = ip3mem(k,ne)
      clwc(k) = CLW(i03+ixy)
      tcol(k) = temp(i03+ixy)
      pcol(k) = Plev(i03+ixy)
      dzcol(k)= dz(i03+ixy)
  ENDDO ! k 

  call get_cld_tr(nz,clwc,pcol,tcol,dzcol,cldflag,iztop,izbot,tr)

  if(.not.cldflag) then ! no cloud
    do k1-,nzz
      i03 = ip3mem(k,ne)
      fcld(i03+ixy)=1.0
    enddo
    cycle
  endif

  call julian(iyear2,imonth2,iday2,jlnd)

  i02=ip2mem(ne)
  XLAT=LATITCRS(i02+ixy)
  XLONG=LONGICRS(i02+ixy) 

  call cal_cosfai(GMT,jlnd,XLONG,XLAT,cosfai)

  do k=1,nzz
    i03 = ip3mem(k,ne)
    if(k.gt.iztop) then  ! above the cloud 
      accld(k)=1.0+alpha*(1.0-tr)*cosfai
    elseif(k.ge.izbot.and.k.le.iztop) then ! in the cloud
      accld(k)=1.4*cosfai 
    elseif(k.lt.izbot) then ! below the cloud
      accld(k)=1.6*tr*cosfai
    endif
    fcld(i03+ixy)=1.0+fracld(i03+ixy)*(accld(k)-1.0) 
  enddo

enddo loop_i
enddo loop_j


end subroutine cal_cldeffect_coeff



subroutine get_cld_tr(nz,clwc,pcol,tcol,dzcol,cldflag,iztop,izbot,tr)
implicit none
real,parameter :: cwmin=0.05 ! g/m3
real,parameter :: rhoh2o=1.0e6 ! g/m3
real,parameter :: radcld=1.0e-5 ! m  10um
integer :: nz
integer :: iz,iztop,izbot
real    :: dzcld,tr,taucld
logical :: cldflag
real,dimension(nz) :: clwc,pcol,tcol,dzcol

cldflag=.false.

loop1 : do iz=nz-1,1,-1
 if(clwc.ge.cwmin) then
   iztop=iz
   cldflag=.true.
   exit loop1
 endif
enddo loop1

if(.not.cldflag) then
 cldcoeff=1.0
 return
endif

if(cldflag) then
 izbot=iztop
 loop2 : do iz=iztop-1,1,-1
  if(clwc.ge.cwmin) then
   izbot=iz
  else
   exit loop2
  endif
 enddo loop2
endif

taucld=0.0
do iz=izbot,iztop
dzcld=dzcol(iz)
taucld=taucld+3.0*clwc(iz)*dzcol(iz)/(2.0*rhoh2o*radcld)
enddo

tr=( 5.0-exp(-1.0*taucld) )/( 4.0+3.0*taucld(1-0.86) )

end subroutine get_cld_tr

