subroutine apm_dry_dep_v2 &
 & ( myid &
 &  ,lapm &
 &  ,dt &
 &  ,imonth2 &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dz &
 &  ,Plev,temp,u_ws,v_ws,QVAPOR,hlayer,hgt1 &
 &  ,land_use,temp2m,LATITCRS,LONGICRS,SWDOWN,PBL_HGT &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

real :: dt

integer :: imonth2

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is,imode

integer :: mem2d,mem3d

integer :: ixy,i02,i03,iapm,iapmode

integer :: ip3mem(nzz,nest),ip2mem(nest)

real :: dryvel=0.0 ! need to be modified


real,dimension(mem3d) :: dz,Plev,temp,QVAPOR,u_ws,v_ws,hlayer
real,dimension(mem2d) :: HGT1,temp2m,LATITCRS,LONGICRS,SWDOWN,PBL_HGT

real,dimension(mem2d) :: land_use

integer :: iwb,ieb,jsb,jeb
real   ,allocatable,dimension(:,:)   :: tsurf,xlat,xlon,SWDOWN1,pblhgt !,vdep
integer,allocatable,dimension(:,:)   :: land
real   ,allocatable,dimension(:,:)   :: gas_vdep
real   ,allocatable,dimension(:,:,:) :: ppp,ttn,uu,vv,height1,QQ,water

real,parameter :: rgfmin=1.0,rgfmax=100.0
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
integer :: icddrt,icdsno     !
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!IF(lapm) THEN ! apm flag
!print*,'dims in sub-ddep',sx(ne),ex(ne),sy(ne),ey(ne),nzz

iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1

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
allocate( QQ(iwb:ieb,jsb:jeb,nzz) )
allocate( water(iwb:ieb,jsb:jeb,nzz) )



allocate( gas_vdep(iwb:ieb,jsb:jeb) )


! get local variables

do j=sy(ne)-1,ey(ne)+1
do i=sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      ppp(i,j,k)=Plev(i03+ixy)
      ttn(i,j,k)=temp(i03+ixy)
      QQ(i,j,k)=QVAPOR(i03+ixy)
      uu(i,j,k)=u_ws(i03+ixy)
      vv(i,j,k)=v_ws(i03+ixy)
      water(i,j,k)=QVAPOR(i03+ixy)*1.E06*29./18.!ppm
      if(k==1) then
        height1(i,j,k) = 2*(hlayer(i03+ixy)*1000.-HGT1(i02+ixy))
      endif
      if(k/=1) then
        height1(i,j,k) = 2*(hlayer(i03+ixy)*1000.-HGT1(i02+ixy))- &
                       & height1(i,j,k-1)
      endif
   enddo    !k     

   CALL getland(LAND_USE(i02+ixy),land(i,j))

   tsurf(i,j)= temp2m(i02+ixy)
   xlat(i,j) = LATITCRS(i02+ixy)
   xlon(i,j) = LONGICRS(i02+ixy)
   SWDOWN1(i,j)=SWDOWN(i02+ixy)
   pblhgt(i,j) = PBL_HGT(i02+ixy)
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

   CALL drydep_gas(myid,imonth2,tsurf,xlat,xlon,QQ,QQ,height1,ppp,uu,vv &
               ,SWDOWN1,land,water,ttn,pblhgt,lrddrt,icddrt,lrdsno,icdsno &
               ,gas_vdep &
               ,sx(ne),ex(ne),sy(ne),ey(ne),nzz,1)

   i02=ip2mem(ne)
   do j=sy(ne)-1,ey(ne)+1
   do  i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      apmdryvel2d(i02+ixy) = gas_vdep(i,j)
if(i.eq.33.and.j.eq.33) then
print*,'vdep-h2so4=',apmdryvel2d(i02+ixy)
endif
   enddo
   enddo

!h2so4_gas=10

   i03=ip3mem(1,ne) ! k=1
   call dry_dep_gas_zhu( myid, h2so4_gas(i03),apmdryvel2d(i02),dz(i03) &
                        ,sx(ne), ex(ne), sy(ne), ey(ne), dt )

   i02=ip2mem(ne)
   do j=sy(ne)-1,ey(ne)+1
   do  i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
if(i.eq.33.and.j.eq.33) then
print*,'con-h2so4=',h2so4_gas(i03+ixy)
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

    loop_so4 : do is=1,nbin

      iapm=ip_sulf(k,is,ne)
      con=apm_sulf(iapm+ixy)

!con=10
      if(lapm_wetsize) then
          diam=rd_sulf(is)*2.0*rgf_sulf(i03+ixy)
      else
          diam=rd_sulf(is)*2.0
      endif
      rhop=densulf*1.0e6 ! g/cm3 -> g/m3

      call drydep_vel_apm_v2( myid,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)
     
      apm_sulf(iapm+ixy)=con

      if(i.eq.33.and.j.eq.33) then
!        print*,'spbin vdep con =',is,vdep,con
      endif

    enddo loop_so4
  endif
!============================



!============================
! apm salt
  if(lfor_salt) then

    nbin=NSEA

    if(lcoated_dyn) then

      allocate(consp(nbin))
      allocate(numcc(nbin))
      allocate(tmp(nbin))
      allocate(radius(nbin)) 

      do is=1,nbin
         iapm=ip_salt(k,is,ne)
         radius(is)=rd_salt(is)
!apm_salt(iapm+ixy)=10
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
!        consp(is)=10
      enddo

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

      call drydep_vel_apm_v2( myid,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)

      if(lcoated_dyn) then
        call dry_dep_apm_sink(myid,deltz,dt,vdep,consp(is))
      endif

      apm_salt(iapm+ixy)=con

      if(i.eq.33.and.j.eq.33) then
!        print*,'salt_bin vdep con =',is,vdep,con,consp(is)
      endif

    enddo loop_salt

    do is=1,nbin
!        if(i.eq.100.and.j.eq.90) print*,'consp',consp(is),is
    enddo

if(i.eq.33.and.j.eq.33) then
!  print*,'salt_coating=',msltsulf(i03+ixy),sum(consp(1:nbin))
endif

    if(lcoated_dyn) msltsulf(i03+ixy)=sum(consp(1:nbin))
    if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

  endif
!============================

!============================
! apm dust
  if(lfor_dust) then

    nbin=NDSTB

    if(lcoated_dyn) then
      allocate(consp(nbin))
      allocate(numcc(nbin))
      allocate(tmp(nbin))
      allocate(radius(nbin))
      do is=1,nbin
         iapm=ip_dust(k,is,ne)
         radius(is)=rd_dust(is)
!apm_dust(iapm+ixy)=10
         numcc(is)= apm_dust(iapm+ixy)/ &
                & radius(is)**3.0d0
         tmp(is)=numcc(is)*radius(is)**2
      enddo
      do is=1,nbin
        iapm=ip_dust(k,is,ne)
!mdstsulf(i03+ixy)=10
        if(sum(tmp(:)).gt.0) then
           consp(is)=mdstsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
           consp(is)=0
        endif
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

      call drydep_vel_apm_v2( myid,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)

      if(lcoated_dyn) then
        call dry_dep_apm_sink(myid,deltz,dt,vdep,consp(is))
      endif

      apm_dust(iapm+ixy)=con

      if(i.eq.33.and.j.eq.33) then
!        print*,'dust_bin vdep con =',is,vdep,con
      endif

    enddo loop_dust

    do is=1,nbin
!        if(i.eq.100.and.j.eq.90) print*,'consp',consp(is),is
    enddo

if(i.eq.33.and.j.eq.33) then
!  print*,'dust_coating=',mdstsulf(i03+ixy),sum(consp(1:nbin))
endif

    if(lcoated_dyn) mdstsulf(i03+ixy)=sum(consp(1:nbin))
    if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

  endif
!========================================================


!============================
! apm dust
  if(lfor_dust) then

    nbin=NBCOCT

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

      call drydep_vel_apm_v2( myid,imonth2 &
             & ,celtsf,cellat,cellon,celheight1,celpre,windu,windv &
             & ,celland,tempk,celpblhgt,diam,rhop,vdep &
             & ,i,j,nzz )

      call dry_dep_apm_sink(myid,deltz,dt,vdep,con)

      if(lcoated_dyn) then
        call dry_dep_apm_sink(myid,deltz,dt,vdep,consp(is))
      endif

      apm_bcoc(iapm+ixy)=con

      if(i.eq.33.and.j.eq.33) then
!        print*,'bcoc_bin vdep con =',is,vdep,con
      endif

    enddo loop_bcoc

    do is=1,nbin
        if(i.eq.100.and.j.eq.90) print*,'consp',consp(is),is
    enddo

if(i.eq.33.and.j.eq.33) then
!  print*,'bc_coating=',mbcsulf(i03+ixy),consp(1)+consp(5)+consp(2)+consp(6)
!  print*,'oc_coating=',mocsulf(i03+ixy),consp(3)+consp(7)+consp(4)+consp(8)
endif

    if(lcoated_dyn) then
      mbcsulf(i03+ixy)=consp(1)+consp(5)+consp(2)+consp(6)
      mocsulf(i03+ixy)=consp(3)+consp(7)+consp(4)+consp(8)
    endif

    if(lcoated_dyn) deallocate(consp,numcc,tmp,radius)

  endif
!===========



!#############
enddo loop_i
enddo loop_j


end subroutine apm_dry_dep_v2



