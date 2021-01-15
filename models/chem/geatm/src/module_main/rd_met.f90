
! chenxsh@mail.iap.ac.cn
! read meteorological fields 
! and calculate vertical mixing coefficient, cloud parameters

subroutine rd_met( myid,iyear,imonth,iday,ihour,iminute &
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
! (1)  read meteorological fields

  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)') ne
 
  call get_funitnaqpms(funit)

  fname='input/wrfd0'//cdnum//'_'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)// '_'//date(9:10)//'.dat'

  inquire(file=fname,exist=lexist)

  if(.not.lexist) then
    print*,trim(fname)//' NOT exist'
    stop
  endif

  open(funit,file=trim(fname),form='unformatted',access='direct',recl=nx(ne)*ny(ne))

 if(myid==0) then
!    print*,'read '//trim(fname)
!    print*
 endif

 irec=1

 ! to read u
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
  
 ! to read v
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

  i0=ip2mem(ne)
 ! to read 2m temperature
  call read2d(myid,T2(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
             nx(ne),ny(ne),irec,funit)

   i0=ip2mem(ne)
 ! to read Syrface Pressue in Pa
  call read2d(myid,PSFC(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
             nx(ne),ny(ne),irec,funit)
   i0=ip2mem(ne)
 ! to read 10m U m/s
  call read2d(myid,U10(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
             nx(ne),ny(ne),irec,funit)

   i0=ip2mem(ne)
 ! to read 10 V m/s
  call read2d(myid,V10(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit) 
  
 ! to read Water vapor mixing ratio (kg/kg)
  do k=1,nz(ne)
    i0=ip3mem(k,ne)
    call read2d(myid,QVAPOR(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

 ! to read CLoud water mixing ratio(kg/kg)
  do k=1,nz(ne)
    i0=ip3mem(k,ne)
    call read2d(myid,clw(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
 
 ! to read Rain water mixing ratio(kg/kg)
  do k=1,nz(ne)
    i0=ip3mem(k,ne)
    call read2d(myid,rnw(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo  

 ! to read SOIL TEMPERATURE (K)
  do k=1,4  ! 4 layers soil
   i0=ip3mem(k,ne)
   call read2d(myid,SOILT(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
            nx(ne),ny(ne),irec,funit)
  enddo

 ! to read RELATIVE SOIL MOISTURE (-)
  do k=1,4  ! 4 layers soil
   i0=ip3mem(k,ne)
   call read2d(myid,SOILRH(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
            nx(ne),ny(ne),irec,funit)
  enddo

 ! to read  SEA ICE FLAG
   i0=ip2mem(ne)
    call read2d(myid,FICE(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)

 ! to read   DOMINANT SOIL CATEGORY (-)
   i0=ip2mem(ne)
    call read2d(myid,FSOIL(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)

 ! to read VEGETATION FRACTION (%)
    i0=ip2mem(ne)
    call read2d(myid,FVEG(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)

 ! to read PHYSICAL SNOW DEPTH (m)
    i0=ip2mem(ne)
    call read2d(myid,FSNOW(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)    

    i0=ip2mem(ne)
 ! to read CUMULUS PRECIPITATION (cm/h) 
    call read2d(myid,RAINCON(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  
    i0=ip2mem(ne)
 ! to read NON-CUMULUS PRECIPITATION (cm/h)
    call read2d(myid,RAINNON(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
    i0=ip2mem(ne)
 ! to read sodown (w/m2) 
    call read2d(myid,SWDOWN(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
              
    i0=ip2mem(ne)
 ! to read UST in m/s
    call read2d(myid,UST(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)

    i0=ip2mem(ne)
 ! to read 1./Monin Ob. Length
    call read2d(myid,RMOL(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
     
    i0=ip2mem(ne)
 ! to read PBL HEIGHT(m) 
    call read2d(myid,PBL_HGT(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)   
 ! to read ice cloud optical depth 

  do k=1,nz(ne)
   i0=ip3mem(k,ne)
    call read2d(myid,TAUCLDI(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                 nx(ne),ny(ne),irec,funit)
  enddo 
              
 ! to read water cloud optical depth   
 
  do k=1,nz(ne)
   i0=ip3mem(k,ne)
     call read2d(myid,TAUCLDC(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                nx(ne),ny(ne),irec,funit)
  enddo                
                             
 
 ! to read model prssue hpa
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,Plev(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
    
 ! to read model height in km
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,h(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
  
 ! to read Temperature in K
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,t(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

 ! to read Relative Humidity (%)
  do k=1,nz(ne)
  i0=ip3mem(k,ne)
  call read2d(myid,rh1(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

   i0=ip2mem(ne) 
 ! to read low cloud fraction 
   call read2d(myid,clflo(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)

   i0=ip2mem(ne)
 ! to read mid cloud fraction
   call read2d(myid,clfmi(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)

   i0=ip2mem(ne)
 ! to read high cloud fraction
   call read2d(myid,clfhi(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)

   i0=ip2mem(ne)
 ! to read relative huminity at 2m
   call read2d(myid,RHSFC(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)    

if(1==2) then
   i0=ip2mem(ne)
 ! to read surface skin temperature
   call read2d(myid,tskwrf(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
else
   tskwrf=T2
endif

!print*,'tskwrf=',tskwrf
!stop


   if(lrd_lai) then
     i0=ip2mem(ne)
     ! to read leaf area index
     if(.not.allocated(wrflai)) stop 'wrflai NOT allocated'
     call read2d(myid,wrflai(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                 nx(ne),ny(ne),irec,funit)

     !print*,'wrflai=',wrflai
     !stop
   endif

   if(idifvert.eq.3) then
     ! to read turbulent kinetic energy
     do k=1,nz(ne)
       i0=ip3mem(k,ne)
!       call read2d(myid,wrftke(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
!               nx(ne),ny(ne),irec,funit)
     enddo
     ! to read elp
     do k=1,nz(ne)
       i0=ip3mem(k,ne)
!       call read2d(myid,wrfelp(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
!               nx(ne),ny(ne),irec,funit)
     enddo
   endif

!print*,'tsk=',tskwrf-T2
!stop

! read hourly accumulative precipitation 
if(lrain_1h) then
    i0=ip2mem(ne)
 ! to read CUMULUS PRECIPITATION (cm/h) 
    call read2d(myid,RAINCON(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)

    i0=ip2mem(ne)
 ! to read shallow CUMULUS PRECIPITATION (cm/h) 
    call read2d(myid,RAINSH(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)


    i0=ip2mem(ne)
 ! to read NON-CUMULUS PRECIPITATION (cm/h)
    call read2d(myid,RAINNON(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
endif

close(funit)


kh=450. ! horizontal diffusion coefficient


!=================================
! (2) calculate cloud parameters

i02=ip2mem(ne)

do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

  rainc=raincon(i02+ixy)*10.0 ! cm/hr -> mm/hr
  rainr=rainnon(i02+ixy)*10.0 ! cm/hr -> mm/hr

  pbl=PBL_HGT(i02+ixy)

  do k=1,nzz
    i03=ip3mem(k,ne)
    ! Layer heights (m)
    dz_1d(k)=dz(i03+ixy)
    zm_1d(k)=heiz(i03+ixy)-terrain(i02+ixy)
    ! Layer pressure (Pa)
    pr(k)=Plev(i03+ixy)
    ! Layer temperature (K)
    tt(k)=t(i03+ixy)
    ! Layer humidity (kg/kg)
    qv(k)=qvapor(i03+ixy)
!qv(k)=0.01
    ! Layer cloud water content (g/m3)
    convfac=44.9*(273./tt(k))*(pr(k)/100/1013.)*29.0
    qq(k)=clw(i03+ixy)*convfac ! kg/kg-> g/m3 
    qqwrf(k)=clw(i03+ixy)*convfac
!qq(k)=0.05
  enddo

  lconv = .false.
  if (rainc.gt.0.2 .and. rainc.gt.rainr) lconv = .true.
  !if (scmeth .eq. 'DIAG') then
  call clddiag(nzz,lconv,pbl,rainc,dz_1d,zm_1d,pr,tt,qv,qq,cfrac)
  !endif

  do k=1,nzz
    if (cfrac(k).le.0.5) then
      cfrac(k) = 0.0
    endif
    if(qq(k).lt.0.01) then
      cfrac(k) = 0.0
    endif
  enddo

  do k=1,nzz
    i03=ip3mem(k,ne)
    cldfrc3d(i03+ixy)=cfrac(k)
    if(.not.(cfrac(k).ge.0.0.and.cfrac(k).le.1.0)) then
      print*,'cldfrc3d-erro:',i,j,k,cfrac(k)
      stop
    endif
  enddo

  cellat=LATITCRS(i02+ixy)
  cellon=LONGICRS(i02+ixy)

  timec=ihour*100+iminute

  datec=iyear*10000+imonth*100+iday
  call juldate(datec)

  call getznth(cellat,cellon,timec,datec,itzon,zenith,ldark)

  zenang = amin1(zenith,60.0)
  zenang = deg2rad*zenang

  do k=1,nzz

     tau=3.0*qqwrf(k)*dz_1d(k)/(2.0*denh2o*1.0e-5) ! 1.0e-5 is the cloud dp( in m)

     cldtrns(k)=(5.0-exp(-1.0*tau))/(4.0+3.0*tau*(1.0-fspfa))  

  enddo
     
  do k=1,nzz

     if(cldtrns(k).ne.1) then
       iabov=0
       fcld=cfrac(k)
     else
       iabov = 1
       fcld=cfrac(1)
     endif

     if (iabov.eq.1) then
        cldrat = 1. + (1. - cldtrns(k))*cos(zenang)
     else
        cldrat = 1.6*cldtrns(k)*cos(zenang)
     endif
     coefcld = 1. + fcld*(cldrat - 1.)
 
     i03=ip3mem(k,ne)
     coefcld3d(i03+ixy)=coefcld

  enddo


enddo
enddo

!stop

!===========================================
! (3) calculate Kv profiles : ysu scheme

if(idifvert.eq.1) then

  do k=1,nzz
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      CALL eddyz(myid,UST(i02),PBL_HGT(i02),RMOL(i02),&
                  heiz(i03),terrain(i02),kv(i03),sx(ne),ex(ne),sy(ne),ey(ne),k)
  enddo

elseif(idifvert.ge.2) then

i02=ip2mem(ne)

do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

  tsfc=tskwrf(i02+ixy) ! surface skin temperature
  press0=psfc(i02+ixy)/100.0 ! pa->mb

  i03_z1=ip3mem(1,ne)
  wind=sqrt(u(i03_z1+ixy)**2.0+v(i03_z1+ixy)**2.0)
  pblc=pbl_hgt(i02+ixy)

  lu_idx=int(land_USE(i02+ixy))
  z0c=z0_USGS(lu_idx)/100.0 ! Surface roughness length (m)

!  ustar=ust(i02+ixy)

  do k=1,nzz

    i03=ip3mem(k,ne)

    tac(k)=t(i03+ixy)
    pac(k)=Plev(i03+ixy) ! mb
    zm(k)=heiz(i03+ixy)-terrain(i02+ixy)

    if(k.eq.1) then
      zzc(k)=2*zm(k)
    elseif(k.gt.1.and.k.lt.nzz) then
      zzc(k)=2*(zm(k)+zm(k+1))
    elseif(k.eq.nzz) then
      zzc(k)=zzc(k-1)+5.0e3 ! not used in subroutine kv_ysu
    endif

    uwind(k)=u(i03+ixy)
    vwind(k)=v(i03+ixy)
    ttt(k)=t(i03+ixy)
    qqq(k)=qvapor(i03+ixy)

    qac(k)=qqq(k)
    th=tac(k)*(1000./pac(k))**gamma
    qvc=qac(k)*18./28.8/1.e6
    thetav(k)=th*(1. + 0.61*qvc)
  enddo

!pblc=1000  

  call wrf_micromet( tac(1),tsfc,pac(1),press0,zm(1),wind,z0c,pblc, &
               & ustar,eli,wstar,risfc )

if(1==2) then
  print*
  print*,'1111=',tac(1),pac(1),zm(1)
  print*,'surf=',tsfc,press0,z0c,wind
  print*,'pblc=',pblc
  print*,'micromet',ustar,eli,wstar,risfc
endif

  
  if(idifvert.eq.2) then ! ysu
    call kv_ysu(nzz,pblc,ustar,eli,wstar,zm,zzc &
                   ,thetav,uwind,vwind,ttt,qqq,kvmin,risfc,rkc)
  elseif(idifvert.eq.3) then ! myj < NO choosing >
    do k=1,nzz
!      tkep(k)=wrftke(i03+ixy) ! 'TKE_PBL'
!      elp(k)=wrfelp(i03+ixy)  ! 'EL_PBL'
    enddo
    call kv_tkemyj(nzz,thetav,uwind,vwind,zm,zzc &
                      ,tkep,elp,kvmin,pblc,rkc)
  elseif(idifvert.eq.4) then ! cmaq
    call kv_cmaq(nzz,pblc,ustar,eli,wstar,zm,zzc &
                    ,thetav,uwind,vwind,kvmin,rkc)
  elseif(idifvert.eq.5) then ! ob70
    call kv_ob70(nzz,kvmin,pblc,zzc,thetav,uwind,vwind,rkc)
  else
    print*,'idifvert erro in rd_met.f90'
  endif


  do k=1,nzz
    i03=ip3mem(k,ne)
    kv(i03+ixy)=rkc(k)
    if(k.eq.1) then
!      print*,i,j,rkc(k)
    endif
  enddo

enddo
enddo

endif ! idifvert.gt.2


!===================================================
!(4) calculate air density
do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

  do k=1,nzz

    i03=ip3mem(k,ne)

    ! kg/m3
    roair3d(i03+ixy)=28.973*1.0E-3*Plev(i03+ixy)*100/(8.31*t(i03+ixy))
!roair3d(i03+ixy)=1.0

  enddo

enddo
enddo
!===================================================

!===================================================
!(4) calculate air density

!===================================================


!stop

end subroutine rd_met

