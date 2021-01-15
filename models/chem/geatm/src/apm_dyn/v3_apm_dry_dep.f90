subroutine apm_dry_dep_v3 &
 & ( myid &
 &  ,lapm &
 &  ,dt &
 &  ,ddep_flag &
 &  ,lrd_lai &
 &  ,imonth2 &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dz &
 &  ,Plev,temp,u_ws,v_ws,QVAPOR,clw,rnw,hlayer,cldfrc3d,hgt1 &
 &  ,land_use,tskwrf,wrflai,LATITCRS,LONGICRS,SWDOWN,PBL_HGT &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d &
 &  ,itzon,iyear,imonth,iday,ihour,iminute )

use apm_varlist
implicit none
include 'apm_parm.inc'

character(len=*) :: ddep_flag

integer :: myid

real :: dt

integer :: imonth2

integer :: itzon,iyear,imonth,iday,ihour,iminute

logical :: lapm

logical :: lrd_lai

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is,imode

integer :: mem2d,mem3d

integer :: ixy,i02,i03,iapm,iapmode

integer :: ip3mem(nzz,nest),ip2mem(nest)

real :: dryvel=0.0 ! need to be modified


real,dimension(mem3d) :: dz,Plev,temp,QVAPOR,u_ws,v_ws,hlayer,cldfrc3d,clw,rnw
real,dimension(mem2d) :: HGT1,tskwrf,LATITCRS,LONGICRS,SWDOWN,PBL_HGT

real,dimension(mem2d) :: land_use

real,dimension(mem2d) :: wrflai

real :: cellai

real :: convfac

integer :: iwb,ieb,jsb,jeb
real   ,allocatable,dimension(:,:)   :: tsurf,xlat,xlon,SWDOWN1,pblhgt !,vdep
integer,allocatable,dimension(:,:)   :: land
real   ,allocatable,dimension(:,:)   :: gas_vdep
real   ,allocatable,dimension(:,:,:) :: ppp,ttn,uu,vv,height1,water

real   ,allocatable,dimension(:,:,:) :: pwc,cwc
real   ,allocatable,dimension(:,:,:) :: fcloud

real   ,allocatable,dimension(:,:)   :: lai2d

real   ,allocatable,dimension(:,:)   :: ddep2d


real,parameter :: rgfmin=1.0,rgfmax=100.0
real,parameter :: diamfine=1.0E-6*sqrt(0.1*2.5)
real :: rgf_tmp
real :: diam

real,dimension(mem2d) :: area0,area,ratio_area

real*8,allocatable,dimension(:) :: numcc,tmp,radius
real*8 :: totalbc,totaloc

! 1:FF 2:BB
integer,parameter :: ridx(1:8) = (/1,2,1,2,1,2,1,2/)
!integer,parameter :: bcidx(1:4) = (/1,2,5,6/)
!integer,parameter :: ocidx(1:4) = (/3,4,7,8/)

character*2,parameter :: flagbcoc(1:8) = (/'bc','bc','oc','oc' &
                                          ,'bc','bc','oc','oc'/)

real    :: vdepm
integer :: imlo,jmlo


!+++ sulferic acid vapor ++++!
logical :: lrddrt,lrdsno     !
!integer :: icddrt,icdsno     !
integer,allocatable,dimension(:,:) :: icddrt,icdsno
!++++++++++++++++++++++++++++!

integer :: iz
integer :: celland
integer :: nbin

real :: celtsf,cellat,cellon,celpblhgt
real :: vdep
real :: rhop
real :: deltz
real :: con

real,dimension(nzz) :: celheight1,celpre,windu,windv,tempk

real,allocatable,dimension(:) :: consp



!
real,parameter :: lvalue=0.0,hvalue=1.0
real :: ratiocore,ratiocoat
real :: concore1,concore2
real :: concoat1,concoat2

real :: ratiocore_bc,ratiocore_oc
real :: concore1_bc,concore2_bc
real :: concore1_oc,concore2_oc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!print*,'rd_binbc=',rd_binbc
!print*,'rd_binoc=',rd_binoc

!stop



!IF(lapm) THEN ! apm flag
!print*,'dims in sub-ddep',sx(ne),ex(ne),sy(ne),ey(ne),nzz

iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1


!allocate( ddep2d(iwb:ieb,jsb:jeb)  )
!ddep2d=-999

allocate( tsurf(iwb:ieb,jsb:jeb) )
allocate( xlat(iwb:ieb,jsb:jeb) )
allocate( xlon(iwb:ieb,jsb:jeb) )
allocate( SWDOWN1(iwb:ieb,jsb:jeb) )
allocate( pblhgt(iwb:ieb,jsb:jeb) )

allocate( land(iwb:ieb,jsb:jeb) )


allocate( ppp(iwb:ieb,jsb:jeb,nzz) )
allocate( ttn(iwb:ieb,jsb:jeb,nzz) )
allocate( uu(iwb:ieb,jsb:jeb,nzz) )
allocate( vv(iwb:ieb,jsb:jeb,nzz) )
allocate( height1(iwb:ieb,jsb:jeb,nzz) )
!allocate( QQ(iwb:ieb,jsb:jeb,nzz) )
allocate( water(iwb:ieb,jsb:jeb,nzz) )
allocate( pwc(iwb:ieb,jsb:jeb,nzz) )
allocate( cwc(iwb:ieb,jsb:jeb,nzz) )
allocate( fcloud(iwb:ieb,jsb:jeb,nzz))

allocate( gas_vdep(iwb:ieb,jsb:jeb) )

allocate( icddrt(iwb:ieb,jsb:jeb) )
allocate( icdsno(iwb:ieb,jsb:jeb) )

! get local variables

do j=sy(ne)-1,ey(ne)+1
do i=sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      ppp(i,j,k)=Plev(i03+ixy)
      ttn(i,j,k)=temp(i03+ixy)
!      QQ(i,j,k)=QVAPOR(i03+ixy)
      uu(i,j,k)=u_ws(i03+ixy)
      vv(i,j,k)=v_ws(i03+ixy)
      water(i,j,k)=QVAPOR(i03+ixy)*1.E06*29./18.!ppm


      convfac=44.9*(273./temp(i03+ixy))*(Plev(i03+ixy)/100/1013.)*29.0

      cwc(i,j,k)=clw(i03+ixy)*convfac
      pwc(i,j,k)=rnw(i03+ixy)*convfac

      fcloud(i,j,k)=cldfrc3d(i03+ixy)

      if(k==1) then
         height1(i,j,k) = 2*(hlayer(i03+ixy)-HGT1(i02+ixy))
      endif
      if(k/=1) then
         height1(i,j,k) = 2*(hlayer(i03+ixy)-HGT1(i02+ixy))-height1(i,j,k-1)
      endif
   enddo    !k     

   if(trim(ddep_flag).eq.'W89') then
     CALL getland(LAND_USE(i02+ixy),land(i,j))
   elseif(trim(ddep_flag).eq.'Z03') then
     CALL getland_usgs2zhang03 (LAND_USE(i02+ixy),land(i,j))
   else
     stop 'ddep_flag-erro in v3_apm_dry_dep.f90'
   endif


   tsurf(i,j)= tskwrf(i02+ixy)
   xlat(i,j) = LATITCRS(i02+ixy)
   xlon(i,j) = LONGICRS(i02+ixy)
   SWDOWN1(i,j)=SWDOWN(i02+ixy)
   pblhgt(i,j) = PBL_HGT(i02+ixy)

   if(lrd_lai) then
     lai2d(i,j)=wrflai(i02+ixy)
   endif

enddo
enddo

if(1==2) then
print*,'tsurf :',tsurf
print*,'xlat :',xlat
print*,'xlon :',xlon
print*,size(land)
print*,'land :',land

print*,'pblhgt:',pblhgt

print*,'height1',height1(:,:,1)
print*,'ppp',ppp(:,:,1)
print*,'uu',uu(:,:,1)
print*,'vv',vv(:,:,1)
print*,'ttn',ttn(:,:,1)
endif

!return

!stop 'kkinputs in ddep'

!========================================================================
!> sulferic acid gas
 if(lfor_h2so4) then
   lrddrt=.false.;icddrt=0  !if droughtness index
   lrdsno=.false.;icdsno=0  !if snow index

   CALL drydep_gas(myid,ddep_flag,fcloud &
               ,imonth2,tsurf,xlat,xlon,pwc,cwc,height1,ppp,uu,vv &
               ,SWDOWN1,land,water &
               ,lrd_lai,lai2d &
               ,ttn,pblhgt,lrddrt,icddrt,lrdsno,icdsno &
               ,gas_vdep &
               ,itzon,iyear,imonth,iday,ihour,iminute &
               ,sx(ne),ex(ne),sy(ne),ey(ne),nzz,1) ! ig=1

   i02=ip2mem(ne)
   do j=sy(ne)-1,ey(ne)+1
   do  i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      apmdryvel2d(i02+ixy) = gas_vdep(i,j)
if(i.eq.33.and.j.eq.33) then
!print*,'vdep-h2so4=',apmdryvel2d(i02+ixy)
endif
   enddo
   enddo

!print*,'h2so4_vdep=',gas_vdep

!h2so4_gas=10

   i03=ip3mem(1,ne) ! k=1
   call dry_dep_gas_zhu( myid, h2so4_gas(i03),apmdryvel2d(i02),dz(i03) &
                        ,sx(ne), ex(ne), sy(ne), ey(ne), dt )

   i02=ip2mem(ne)
   do j=sy(ne)-1,ey(ne)+1
   do  i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
if(i.eq.33.and.j.eq.33) then
!print*,'con-h2so4=',h2so4_gas(i03+ixy)
endif
   enddo
   enddo

 endif
!=======================================================================


!pause


!=======================================================================
!----------------------------- particles -------------------------------

k=1

i02=ip2mem(ne)
i03=ip3mem(k,ne)

loop_j : do j = sy(ne),ey(ne)
loop_i : do i = sx(ne),ex(ne)
!############################


  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

!*************************
  celtsf=tsurf(i,j)
  cellat=xlat(i,j)
  cellon=xlon(i,j)
  celland=land(i,j)
  celpblhgt=pblhgt(i,j)

  if(lrd_lai) then
    cellai=lai2d(i,j)
  endif

  do iz=1,nzz
    celheight1(iz)=height1(i,j,iz)
    celpre(iz)=ppp(i,j,iz)
    windu(iz)=uu(i,j,iz)
    windv(iz)=vv(i,j,iz)
    tempk(iz)=ttn(i,j,iz) 
  enddo
  deltz=dz(i03+ixy)
!*************************


!============================
! apm sulfate
  if(lfor_sulf) then

    nbin=NSO4

    concore1=0
    concoat1=0

    loop_so4 : do is=1,nbin
!     loop_so4 : do is=30,30     

      iapm=ip_sulf(k,is,ne)
      con=apm_sulf(iapm+ixy)

!con=10
      if(lapm_wetsize) then
          diam=rd_sulf(is)*2.0*rgf_sulf(i03+ixy)
      else
          diam=rd_sulf(is)*2.0
      endif
      rhop=densulf*1.0e6 ! g/cm3 -> g/m3

diam = diamfine

      call drydep_vel_apm_v2( myid,ddep_flag,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,lrd_lai,cellai &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)
     
      apm_sulf(iapm+ixy)=con

!      ddep2d(i,j)=vdep


      if(i.eq.33.and.j.eq.33) then
!        print*,'spbin vdep con =',is,vdep,con
      endif

    enddo loop_so4
  endif
!============================

!cycle

!============================
! apm salt
  if(lfor_salt) then

    nbin=NSEA

    concore1=0
    concoat1=0

    concore2=0
    concoat2=0

    if(lcoated_dyn) then

      allocate(consp(nbin))
      allocate(numcc(nbin))
      allocate(tmp(nbin))
      allocate(radius(nbin)) 

      do is=1,nbin
         iapm=ip_salt(k,is,ne)
         radius(is)=rd_salt(is)
!apm_salt(iapm+ixy)=10
         concore1=concore1+apm_salt(iapm+ixy)

         numcc(is)= apm_salt(iapm+ixy)/ &
                & radius(is)**3.0d0
         tmp(is)=numcc(is)*radius(is)**2
      enddo

      do is=1,nbin
        iapm=ip_salt(k,is,ne)
!msltsulf(i03+ixy)=10
        if(sum(tmp(:)).gt.0) then
           consp(is)=msltsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
           consp(is)=0
        endif
        concoat1=concoat1+consp(is)
!        consp(is)=10
      enddo

    endif

if(i.eq.33.and.j.eq.33) then
! print*,'bfdep'
! print*,'core&coat00=',concore1,concoat1
! stop
endif

    loop_salt : do is=1,nbin

      iapm=ip_salt(k,is,ne)

      con=apm_salt(iapm+ixy)

!      con=10

      if(lapm_wetsize) then
          diam=rd_salt(is)*2.0*rgf_salt(i03+ixy)
      else
          diam=rd_salt(is)*2.0
      endif
      rhop=densalt*1.0e6 ! g/cm3 -> g/m3

      call drydep_vel_apm_v2( myid,ddep_flag,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,lrd_lai,cellai &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)

      if(lcoated_dyn) then
        call dry_dep_apm_sink(myid,deltz,dt,vdep,consp(is))
      endif

      apm_salt(iapm+ixy)=con

      concore2=concore2+con
!      concoat2=concoat2+consp(is)

      if(i.eq.33.and.j.eq.33) then
!        print*,'salt_bin vdep con =',is,vdep,con,consp(is)
      endif

    enddo loop_salt

    do is=1,nbin
!        if(i.eq.100.and.j.eq.90) print*,'consp',consp(is),is
    enddo

    ratiocore=amin1(amax1(concore2/concore1,lvalue),hvalue)

if(i.eq.33.and.j.eq.33) then
!  print*,'salt_coating=',msltsulf(i03+ixy),sum(consp(1:nbin))
!print*,'afdep'
!print*,'core&coat22=',concore2,concoat2
!print*,'core&coat22=',concore2/concore1,concoat2/concoat1
!print*,'salt ratiocore=',ratiocore
!stop
endif


    if(lcoated_dyn) then
!        msltsulf(i03+ixy)=sum(consp(1:nbin))
        msltsulf(i03+ixy)=msltsulf(i03+ixy)*ratiocore
    endif

    if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

  endif
!============================

!============================
! apm dust
  if(lfor_dust) then

    nbin=NDSTB

    concore1=0
    concoat1=0

    concore2=0
    concoat2=0

    if(lcoated_dyn) then
      allocate(consp(nbin))
      allocate(numcc(nbin))
      allocate(tmp(nbin))
      allocate(radius(nbin))
      do is=1,nbin
         iapm=ip_dust(k,is,ne)
         radius(is)=rd_dust(is)
!apm_dust(iapm+ixy)=10
         concore1=concore1+apm_dust(iapm+ixy)
         numcc(is)= apm_dust(iapm+ixy)/ &
                & radius(is)**3.0d0
         tmp(is)=numcc(is)*radius(is)**2
      enddo

!      concoat1=mdstsulf(i03+ixy)

      do is=1,nbin
        iapm=ip_dust(k,is,ne)
!mdstsulf(i03+ixy)=10
        if(sum(tmp(:)).gt.0) then
           consp(is)=mdstsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
           consp(is)=0
        endif
        concoat1=concoat1+consp(is)
!        consp(is)=10
      enddo
    endif

    loop_dust : do is=1,nbin

      iapm=ip_dust(k,is,ne)

      con=apm_dust(iapm+ixy)

!      con=10

      if(lapm_wetsize) then
          diam=rd_dust(is)*2.0*rgf_dust(i03+ixy)
      else
          diam=rd_dust(is)*2.0
      endif
      rhop=dendust*1.0e6 ! g/cm3 -> g/m3

      call drydep_vel_apm_v2( myid,ddep_flag,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,lrd_lai,cellai &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)

      if(lcoated_dyn) then
        call dry_dep_apm_sink(myid,deltz,dt,vdep,consp(is))
      endif

      apm_dust(iapm+ixy)=con

      concore2=concore2+con
      concoat2=concoat2+consp(is)

      if(i.eq.33.and.j.eq.33) then
!        print*,'dust_bin vdep con =',is,vdep,con
      endif

    enddo loop_dust

    do is=1,nbin
!        if(i.eq.100.and.j.eq.90) print*,'consp',consp(is),is
    enddo

    ratiocore=amin1(amax1(concore2/concore1,lvalue),hvalue)

if(i.eq.33.and.j.eq.33) then
!  print*,'dust_coating=',mdstsulf(i03+ixy),sum(consp(1:nbin))
!print*,'core&coat00=',concore1,concoat1
!print*,'core&coat22=',concore2,concoat2
!print*,'dust core&coatch=',concore2/concore1,concoat2/concoat1
endif

    if(lcoated_dyn) then
!       mdstsulf(i03+ixy)=sum(consp(1:nbin))
      mdstsulf(i03+ixy)=mdstsulf(i03+ixy)*ratiocore
    endif

    if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

  endif
!========================================================


!============================
! apm dust
  if(lfor_bcoc) then

    nbin=NBCOCT

    concore1_bc=0.0
    concore1_oc=0.0
    concore2_bc=0.0
    concore2_oc=0.0

    if(lcoated_dyn) then

      allocate(consp(nbin))
      allocate(numcc(nbin))
      allocate(tmp(nbin))
      allocate(radius(nbin))

      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
        radius(is)=ref1d_bcoc(ridx(is))
!apm_bcoc(iapm+ixy)=10.0        
        numcc(is)= apm_bcoc(iapm+ixy)/ &
               & radius(is)**3.0d0
        tmp(is)=numcc(is)*radius(is)**2
        if(flagbcoc(is).eq.'bc') then
          concore1_bc=concore1_bc+apm_bcoc(iapm+ixy)
        elseif(flagbcoc(is).eq.'oc') then
          concore1_oc=concore1_oc+apm_bcoc(iapm+ixy)
        endif
      enddo
      totalbc=tmp(1)+tmp(5)+tmp(2)+tmp(6)
      totaloc=tmp(3)+tmp(7)+tmp(4)+tmp(8)

      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
!mbcsulf(i03+ixy)=10.0
!mocsulf(i03+ixy)=10.0
        if(flagbcoc(is).eq.'bc') then
          if(totalbc.gt.0) then
            consp(is)=mbcsulf(i03+ixy)*tmp(is)/totalbc
          else
            consp(is)=0
          endif
        elseif(flagbcoc(is).eq.'oc') then
          if(totaloc.gt.0) then
             consp(is)=mocsulf(i03+ixy)*tmp(is)/totaloc
          else
             consp(is)=0
          endif
        endif
!        consp(is)=10
      enddo

    endif

    loop_bcoc : do is=1,nbin

      iapm=ip_bcoc(k,is,ne)

      con=apm_bcoc(iapm+ixy)

!      con=10

      if(lapm_wetsize) then
         if(flagbcoc(is).eq.'bc') then
            rgf_tmp=rgf_bc(i03+ixy)
         elseif(flagbcoc(is).eq.'oc') then
            rgf_tmp=rgf_oc(i03+ixy)
         endif
      else
         rgf_tmp=1.0
      endif
      diam=ref1d_bcoc(ridx(is))*2.0*rgf_tmp

      rhop=denbcoc*1.0e6 ! g/cm3 -> g/m3

diam=diamfine

      call drydep_vel_apm_v2( myid,ddep_flag,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,lrd_lai,cellai &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)

      if(lcoated_dyn) then
        call dry_dep_apm_sink(myid,deltz,dt,vdep,consp(is))
      endif

      apm_bcoc(iapm+ixy)=con

      if(flagbcoc(is).eq.'bc') then
        concore2_bc=concore2_bc+con
      elseif(flagbcoc(is).eq.'oc') then
        concore2_oc=concore2_oc+con
      endif

      if(i.eq.33.and.j.eq.33) then
!        print*,'bcoc_bin vdep con =',is,vdep,con
      endif

    enddo loop_bcoc

    do is=1,nbin
!        if(i.eq.100.and.j.eq.90) print*,'consp',consp(is),is
    enddo

    ratiocore_bc=amin1(amax1(concore2_bc/concore1_bc,lvalue),hvalue)
    ratiocore_oc=amin1(amax1(concore2_oc/concore1_oc,lvalue),hvalue)

if(i.eq.33.and.j.eq.33) then
!  print*,'bc_coating=',mbcsulf(i03+ixy),consp(1)+consp(5)+consp(2)+consp(6)
!  print*,'oc_coating=',mocsulf(i03+ixy),consp(3)+consp(7)+consp(4)+consp(8)
!print*,'ratiocore_bc=',ratiocore_bc
!print*,'ratiocore_oc=',ratiocore_oc
endif

    if(lcoated_dyn) then
!!      mbcsulf(i03+ixy)=consp(1)+consp(5)+consp(2)+consp(6)
!!      mocsulf(i03+ixy)=consp(3)+consp(7)+consp(4)+consp(8)
!      mbcsulf(i03+ixy)=mbcsulf(i03+ixy)*ratiocore_bc
!      mocsulf(i03+ixy)=mocsulf(i03+ixy)*ratiocore_oc
    endif

    if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

  endif
!===========


    nbin=nbincb

    concore1=0
    concoat1=0

    concore2=0
    concoat2=0
      allocate(consp(nbin))
      allocate(numcc(nbin))
      allocate(tmp(nbin))
      allocate(radius(nbin))

      do is=1,nbin
         iapm=ip_cbbin(k,is,ne)
         radius(is)=rd_binbc(is)
         concore1=concore1+apm_binbc(iapm+ixy)

         numcc(is)= apm_binbc(iapm+ixy)/ &
                & radius(is)**3.0d0
         tmp(is)=numcc(is)*radius(is)**2
      enddo

      do is=1,nbin
        iapm=ip_cbbin(k,is,ne)
        if(sum(tmp(:)).gt.0) then
           consp(is)=mbcsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
           consp(is)=0
        endif
        concoat1=concoat1+consp(is)
!        consp(is)=10
      enddo

    loop_binbc : do is=1,nbincb

      iapm=ip_cbbin(k,is,ne)

      con=apm_binbc(iapm+ixy)

      diam=rd_binbc(is)*2.0

diam=diamfine

      rhop=denbcoc*1.0e6 ! g/cm3 -> g/m3

      call drydep_vel_apm_v2( myid,ddep_flag,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,lrd_lai,cellai &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)

      call dry_dep_apm_sink(myid,deltz,dt,vdep,consp(is))

      apm_binbc(iapm+ixy)=con

      concore2=concore2+con

    enddo loop_binbc

    ratiocore=amin1(amax1(concore2/concore1,lvalue),hvalue)

    mbcsulf(i03+ixy)=mbcsulf(i03+ixy)*ratiocore

    deallocate(consp,numcc,tmp,radius)
!------------------------------------------
!------------------------------------------

    nbin=nbincb

    concore1=0
    concoat1=0

    concore2=0
    concoat2=0
      allocate(consp(nbin))
      allocate(numcc(nbin))
      allocate(tmp(nbin))
      allocate(radius(nbin))

      do is=1,nbin
         iapm=ip_cbbin(k,is,ne)
         radius(is)=rd_binoc(is)
         concore1=concore1+apm_binoc(iapm+ixy)

         numcc(is)= apm_binoc(iapm+ixy)/ &
                & radius(is)**3.0d0
         tmp(is)=numcc(is)*radius(is)**2
      enddo

      do is=1,nbin
        iapm=ip_cbbin(k,is,ne)
        if(sum(tmp(:)).gt.0) then
           consp(is)=mocsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
           consp(is)=0
        endif
        concoat1=concoat1+consp(is)
!        consp(is)=10
      enddo

    loop_binoc : do is=1,nbincb

      iapm=ip_cbbin(k,is,ne)

      con=apm_binoc(iapm+ixy)

      diam=rd_binoc(is)*2.0

diam=diamfine


      rhop=densalt*1.0e6 ! g/cm3 -> g/m3

      call drydep_vel_apm_v2( myid,ddep_flag,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,lrd_lai,cellai &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)

      call dry_dep_apm_sink(myid,deltz,dt,vdep,consp(is))

      apm_binoc(iapm+ixy)=con

      concore2=concore2+con

    enddo loop_binoc

    ratiocore=amin1(amax1(concore2/concore1,lvalue),hvalue)

    mocsulf(i03+ixy)=mocsulf(i03+ixy)*ratiocore

    deallocate(consp,numcc,tmp,radius)


!#############
enddo loop_i
enddo loop_j


deallocate( ppp,ttn,land,tsurf,xlat,xlon,pwc,cwc,height1,uu,vv &
           ,water,SWDOWN1,pblhgt,fcloud,icddrt,icdsno)

!print*,ddep2d
!stop

!deallocate(ddep2d)

end subroutine apm_dry_dep_v3



