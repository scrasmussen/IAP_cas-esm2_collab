
! chenxsh@mail.iap.ac.cn 20170209
! calculate DMS emission

subroutine cal_dms_emis( myid,iyear,imonth,iday,ihour,iminute &
                        ,nest,nzz,nx,ny,nz,sx,ex,sy,ey,ne,mem2d )

use naqpms_varlist,  only: ip2mem,ip3mem,ip_emit2d,emit2d
use met_fields,      only: u10,v10,tskwrf,fice
use naqpms_gridinfo, only: LATITCRS,LONGICRS,land_use 


implicit none 


integer :: myid
integer :: ne,nest
integer :: nx(5),ny(5),nz(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: mem2d

integer :: iyear,imonth,iday,ihour,iminute


integer :: i,j,k,i0,i02,i03,ixy

character :: cdnum*1
character :: date*10,fname*100

integer :: irec
integer :: funit
logical :: lexist

!=========

real,allocatable,dimension(:) :: sdmsc

real :: cellat,cellon,tau


!===========
integer :: i0emt,ig

real :: tsfc,wspd,conc,seaice,emis

integer :: lu_idx

integer,parameter  :: lucats_usgs=24
real z0_USGS(lucats_usgs)
data z0_USGS /50.,15.,15.,15.,14.,20.,12.,10.,11.,15., &
              50.,50.,50.,50.,50.,0.01,20.,40.,10.,10., &
              30.,15.,10.,5./


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allocate(sdmsc(mem2d))


!=============================================
! (1)  read DMS concentration in the sea water

  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)') ne
 
  call get_funitnaqpms(funit)

  fname='emit/data.emit/sdmsc_d0'//cdnum//'.'//date(5:6)//'.dat'

  inquire(file=fname,exist=lexist)

  if(.not.lexist) then
    print*,trim(fname)//' NOT exist'
    stop
  endif

  open(funit,file=trim(fname),form='unformatted',access='direct',recl=nx(ne)*ny(ne))

 if(myid==0) then
    print*,'read '//trim(fname)
 endif

 irec=1

 ! to read sdmsc
 i0=ip2mem(ne)
 call read2d(myid,sdmsc(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  
 close(funit)


!=================================
! (2) calculate DMS emission flux

i02=ip2mem(ne)

do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

  conc=sdmsc(i02+ixy)

  wspd=sqrt(u10(i02+ixy)**2+v10(i0+ixy)**2)
  tsfc=tskwrf(i02+ixy)
  seaice=fice(i02+ixy)
  lu_idx=land_use(i02+ixy) ! juanxiong he
 if(conc.ge.0.and.lu_idx.eq.16) then ! chenxsh@20181206
   call dms_flux(conc,tsfc,wspd,seaice,emis,lu_idx) ! juanxiong he
   if(.not.(emis.ge.0.and.emis.le.100)) then
      print*,'dms_emis-err',i,j,emis,conc,tsfc,wspd,seaice
   endif
 else
   emis=0.0
 endif
!emis=0.0
 ig=57 ! DMS
 i0emt=ip_emit2d(ig,1,ne)
 emit2d(i0emt+ixy)=emis

enddo
enddo

!stop

!=================================

deallocate(sdmsc)

end subroutine cal_dms_emis


subroutine dms_flux(S_DMS,SST,WS10,seaice,emit_DMS,lu_idx) ! juanxiong he
! weiy@20170209
! S_DMS    : DMS conc in nmol/L
! WS10     : wind speed at 10m in m/s
! SST      : sea surface temperature in K
! emit_DMS : DMS emission flux in ug/m2/s
      implicit none
      integer,parameter  :: mth=12,nlon=360,nlat=180
      real,parameter     :: a1=0.222,a2=0.333
      real,parameter     :: b1=2674.0,b2=147.12,b3=3.726,b4=0.038
      real,parameter     :: c1=3525,c2=9.464
      real,parameter     :: d=659
      real,parameter     :: f=-0.5
      real,parameter     :: nan=-18768.20
      integer,parameter  :: mho=18,mdms=62

      character :: cmon*2 ,fname1*200 ,fname2*200 ,fname3*200
      real    :: k600,sc_dms,kw,af,ka,raf,kt
      real    :: SST,WS10,emit_DMS,S_DMS
      integer :: i,j,k,itt,imon,irec
      integer :: error, di, dj
      integer :: lu_idx ! juanxiong he
      real :: seaice

if(seaice.le.0.and.lu_idx.eq.16) then ! juanxiong he
!if(seaice.le.0) then ! juanxiong he
      if(sst.gt.321.0) sst=321 ! juanxiong he
      k600=a1*WS10**2+a2*WS10
      sc_dms=b1-b2*(SST-273.15)+b3*(SST-273.15)**2-b4*(SST-273.15)**3
      kw=k600*(sc_dms/600)**f

      af=exp(c1/SST-c2)
      ka=d*WS10*(mdms/mho*1.0)**f ! juanxiong he
      raf=1/(1+ka/(af*kw))
      kt=kw*(1-raf)
      emit_DMS=S_DMS*kt*mdms*(0.01/3600)
else
      emit_DMS=0.0
endif

end subroutine dms_flux

