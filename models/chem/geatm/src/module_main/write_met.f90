! write by winne at 2017.08.25
! weiying@mail.iap.ac.cn
! write meteorological fields 

subroutine write_met( myid,iyear,imonth,iday,ihour,iminute &
                  ,nest,nzz,nx,ny,nz,sx,ex,sy,ey,ne )

use naqpms_varlist, only: ip2mem,ip3mem
use naqpms_varlist, only: itzon
use naqpms_varlist, only: idifvert 
use naqpms_varlist, only: lrd_lai
use met_fields
use naqpms_gridinfo, only: heiz,terrain,dz,LATITCRS,LONGICRS,land_use 


implicit none 

logical,parameter :: lrain_1h=.false.

real,parameter :: deg2rad=0.01745329
real,parameter :: gamma=0.286

real,parameter :: kvmin=0.5

real,parameter :: denh2o=1.0e6 ! g/m3
real,parameter :: fspfa=0.86 ! scattering phase function asymmetry factor


integer :: myid
integer :: ne,nest
integer :: nx(5),ny(5),nz(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)


integer :: iyear,imonth,iday,ihour,iminute


integer :: i,j,k,i0,i02,i03,ixy

character :: cdnum*1
character :: date*10,fname*100
character :: cmyid*3

integer :: irec
integer :: funit
logical :: lexist

!=========

logical :: lconv
real,dimension(nzz) :: dz_1d,zm_1d,pr,tt,qv,qq,cfrac,cldtrns,qqwrf
real :: pbl,convfac,rainc,rainr


real :: zenang,zenith
real :: cldrat,coefcld,fcld
integer :: iabov

real :: cellat,cellon,tau

integer :: timec,datec
logical :: ldark

!===========
integer :: i03_z1
real :: ustar,eli,wstar,risfc

real :: tsfc,press0,wind,z0c,pblc

real :: th,qvc

real,dimension(nzz) :: tac,pac,zm,qac
real,dimension(nzz) :: zzc,thetav,uwind,vwind,ttt,qqq,rkc

real,dimension(nzz) :: tkep,elp

integer :: lu_idx

integer,parameter  :: lucats_usgs=24
real z0_USGS(lucats_usgs)
data z0_USGS /50.,15.,15.,15.,14.,20.,12.,10.,11.,15., &
              50.,50.,50.,50.,50.,0.01,20.,40.,10.,10., &
              30.,15.,10.,5./


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!===================================
! (1)  write meteorological fields

  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)') ne
  write (cmyid,'(i3.3)') myid
  call system("mkdir -p met_out/tmp")
  call system("mkdir -p met_out/tmp/"//date(1:8))

  fname='met_out/tmp/'//date(1:8)//'/metd'//cdnum//'.'//cmyid//'.'//date
call get_funitnaqpms(funit)

!  open(funit,file=trim(fname),form='unformatted',access='direct',recl=nx(ne)*ny(ne))

 open(funit,file=trim(fname),form='unformatted',status='UNKNOWN')

 !irec=1

 ! to read u
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call write2d_v2(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  
 ! to read v
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call write2d_v2(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  i0=ip2mem(ne)
 ! to read 2m temperature
  call write2d_v2(myid,T2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

   i0=ip2mem(ne)
 ! to read Syrface Pressue in Pa
  call write2d_v2(myid,PSFC(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

   i0=ip2mem(ne)
 ! to read 10m U m/s
  call write2d_v2(myid,U10(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

   i0=ip2mem(ne)
 ! to read 10 V m/s
  call write2d_v2(myid,V10(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  
 ! to read Water vapor mixing ratio (kg/kg)
  do k=1,nz(ne)
    i0=ip3mem(k,ne)
    call write2d_v2(myid,QVAPOR(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

 ! to read CLoud water mixing ratio(kg/kg)
  do k=1,nz(ne)
    i0=ip3mem(k,ne)
    call write2d_v2(myid,clw(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
 
 ! to read Rain water mixing ratio(kg/kg)
  do k=1,nz(ne)
    i0=ip3mem(k,ne)
    call write2d_v2(myid,rnw(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo  

 ! to read SOIL TEMPERATURE (K)
  do k=1,4  ! 4 layers soil
   i0=ip3mem(k,ne)
   call write2d_v2(myid,SOILT(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

 ! to read RELATIVE SOIL MOISTURE (-)
  do k=1,4  ! 4 layers soil
   i0=ip3mem(k,ne)
   call write2d_v2(myid,SOILRH(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

 ! to read  SEA ICE FLAG
   i0=ip2mem(ne)
    call write2d_v2(myid,FICE(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

 ! to read   DOMINANT SOIL CATEGORY (-)
   i0=ip2mem(ne)
    call write2d_v2(myid,FSOIL(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

 ! to read VEGETATION FRACTION (%)
    i0=ip2mem(ne)
    call write2d_v2(myid,FVEG(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

 ! to read PHYSICAL SNOW DEPTH (m)
    i0=ip2mem(ne)
    call write2d_v2(myid,FSNOW(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

    i0=ip2mem(ne)
 ! to read CUMULUS PRECIPITATION (cm/h) 
    call write2d_v2(myid,RAINCON(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  
    i0=ip2mem(ne)
 ! to read NON-CUMULUS PRECIPITATION (cm/h)
    call write2d_v2(myid,RAINNON(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

    i0=ip2mem(ne)
 ! to read sodown (w/m2) 
    call write2d_v2(myid,SWDOWN(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
              
    i0=ip2mem(ne)
 ! to read UST in m/s
    call write2d_v2(myid,UST(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

    i0=ip2mem(ne)
 ! to read 1./Monin Ob. Length
    call write2d_v2(myid,RMOL(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
     
    i0=ip2mem(ne)
 ! to read PBL HEIGHT(m) 
    call write2d_v2(myid,PBL_HGT(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
 ! to read ice cloud optical depth 

  do k=1,nz(ne)
   i0=ip3mem(k,ne)
    call write2d_v2(myid,TAUCLDI(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo 
              
 ! to read water cloud optical depth   
 
  do k=1,nz(ne)
   i0=ip3mem(k,ne)
     call write2d_v2(myid,TAUCLDC(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo                
                             
 
 ! to read model prssue hpa
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call write2d_v2(myid,Plev(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
    
 ! to read model height in km
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call write2d_v2(myid,h(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  
 ! to read Temperature in K
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call write2d_v2(myid,t(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

 ! to read Relative Humidity (%)
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call write2d_v2(myid,rh1(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

   i0=ip2mem(ne) 
 ! to read low cloud fraction 
   call write2d_v2(myid,clflo(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

   i0=ip2mem(ne)
 ! to read mid cloud fraction
   call write2d_v2(myid,clfmi(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

   i0=ip2mem(ne)
 ! to read high cloud fraction
   call write2d_v2(myid,clfhi(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

   i0=ip2mem(ne)
 ! to read relative huminity at 2m
   call write2d_v2(myid,RHSFC(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

if(1==1) then
   i0=ip2mem(ne)
 ! to read surface skin temperature
   call write2d_v2(myid,tskwrf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
else
   tskwrf=T2
endif  

   if(lrd_lai) then
     i0=ip2mem(ne)
     ! to read leaf area index
     if(.not.allocated(wrflai)) stop 'wrflai NOT allocated'
     call write2d_v2(myid,wrflai(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
   endif

   !i0=ip2mem(ne)
 ! to read sea surface temperature
   !call write2d_v2(myid,sstwrf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

   if(idifvert.eq.3) then
     ! to read turbulent kinetic energy
     do k=1,nz(ne)
       i0=ip3mem(k,ne)
!       call write2d_v2(myid,wrftke(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
!               nx(ne),ny(ne),irec,funit)
     enddo
     ! to read elp
     do k=1,nz(ne)
       i0=ip3mem(k,ne)
!       call write2d_v2(myid,wrfelp(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
!               nx(ne),ny(ne),irec,funit)
     enddo
   endif

!print*,'tsk=',tskwrf-T2
!stop

! read hourly accumulative precipitation 
if(lrain_1h) then
    i0=ip2mem(ne)
 ! to read CUMULUS PRECIPITATION (cm/h) 
    call write2d_v2(myid,RAINCON(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

    i0=ip2mem(ne)
 ! to read shallow CUMULUS PRECIPITATION (cm/h) 
    call write2d_v2(myid,RAINSH(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)


    i0=ip2mem(ne)
 ! to read NON-CUMULUS PRECIPITATION (cm/h)
    call write2d_v2(myid,RAINNON(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
endif

close(funit)



!===================================================



end subroutine write_met

