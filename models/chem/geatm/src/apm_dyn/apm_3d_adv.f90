
subroutine apm_3d_adv &
 & ( myid &
 &  ,imasskeep &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,GC_MOLWT &
 &  ,mem3d,RatioMass &
 &  ,mem2d,ktop &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,ifsm,idmSet,ismMax,igMark &
 &  ,hh )

use apm_varlist

use naqpms_varlist
use naqpms_gridinfo, only : dx,dy,dz,terrain,mpfac,pzps
use met_fields, only : u,v,w,t,Plev,roair3d,entrn3d,dilut3d,dsdt3d

use adv1d_comv

implicit none

include 'apm_parm.inc'

integer :: myid

real :: dt,dt_naqpms

logical :: lapm
integer :: imasskeep

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is

integer :: mem3d

real,dimension(igas) :: GC_MOLWT

real,dimension(mem3d) :: RatioMass,kpmass_m2

integer :: mem2d
real,dimension(mem2d) :: ktop

integer :: ifsm(5)

integer :: idmSet,ismMax

integer :: igMark(idmSet)

integer :: ixy,i03,iapm



integer :: igas,iaer,isize,nseacom,ndustcom


real,dimension(mem3d) :: wk

integer :: ig,i04,i02,ia,iduc,i05,i05c,i04aer
integer :: idm,ism,i04sm

integer :: letdoit

real,allocatable,dimension(:,:,:) :: sm

logical :: lflag 
real :: uws,vws,deltx,delty

real :: rrr


!===================
!character(len=*),parameter :: flagadv='walcek'
!character(len=*),parameter :: flagadv='camx_ppm'
!character(len=*),parameter :: fvadv='implicit' 
character(len=*),parameter :: fvadv='explicit'
integer,parameter :: iroadv=1
real :: cfl

!=============
! 
real :: hh


! 
!=====================================================
! wk arrays
 
 real    :: vmax,vmin,vbot,vtop,ulft,urgt
 real    :: tstep,ucell,dxcell
 integer :: nstep,istep
 integer :: i03_bot,i03_top

!1d -> 3d
 real,allocatable,dimension(:,:,:) :: u3d,v3d,w3d,p3d,t3d,dx3d,dy3d,dz3d,den3d,dsdtf
 real,allocatable,dimension(:,:)   :: mpfac2d,ter2d,wgtfac,ter1_2d,ter2_2d

 real,allocatable,dimension(:,:,:) :: jacobian ! (patial_z/patial_eta)

! horizontal
 real,allocatable,dimension(:) :: Q0,QN,DXX,DEN0,DEN1,u1D,DD0

 real,allocatable,dimension(:) :: p1d,t1d,rhoa,denvadv

 real,allocatable,dimension(:,:) :: denair2d

 real,allocatable,dimension(:,:,:) :: u_xstag,ro_xstag
 real,allocatable,dimension(:,:,:) :: v_ystag,ro_ystag

! vertical
 real,allocatable,dimension(:,:,:) :: ro_zstag
 real,dimension(1:nzz) :: denair1d
 real ::  DXX0
 integer :: sz,ez
!========================================================

 integer   :: izkk
 character :: label*50

! camx
 integer :: nndim
 real :: dxres,dyres
 real,allocatable,dimension(:) :: conc,vel,area,areav,conc11
 real,allocatable,dimension(:) :: depth,entrn,dilut
 real,allocatable,dimension(:) :: conc00,sl

 real :: mwgt


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(trim(flagadv).eq.'camx_ppm') then
  cfl=0.5
elseif(trim(flagadv).eq.'walcek') then
  cfl=1.0
endif


!stop



if(1==2) then
 do j = sy(ne)-1,ey(ne)+1
 do i = sx(ne)-1,ex(ne)+1

if((j.eq.sy(ne)-1.or.j.eq.ey(ne)+1).and.(i.eq.sx(ne)-1.or.i.eq.ex(ne)+1)) cycle

    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    do k=1,nzz
        i03=ip3mem(k,ne)
        uws=u(i03+ixy)
        vws=v(i03+ixy)
        deltx=dx(i03+ixy)
        delty=dy(i03+ixy)

rrr=3**(ne-1)

        lflag=(abs(uws).lt.1.0e3).and.(abs(vws).lt.1.0e3).and.&
             &(deltx.eq.27000/rrr).and.(delty.eq.27000/rrr)
        if(.not.lflag) then
          print*,'erro-hadv',sx(ne),ex(ne),sy(ne),ey(ne)
          print*,'erro-hadv11',i,j,k
          print*,'erro-hadv22',uws, vws,deltx,delty,27000/rrr 
        endif
    enddo
 enddo
 enddo
endif


!return


 allocate( u3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,v3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,w3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,p3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,t3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,dx3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,dy3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,dz3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,den3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) &
          ,dsdtf(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) )




 allocate( jacobian(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz) )

 allocate( mpfac2d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) &
          ,ter2d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) )


 allocate( ter1_2d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) &
          ,ter2_2d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) )


! allocate( wgtfac(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) )

 allocate( denair2d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) )

 allocate( u_xstag(sx(ne)-1:ex(ne),sy(ne):ey(ne),nzz) &
          ,ro_xstag(sx(ne)-1:ex(ne),sy(ne):ey(ne),nzz) &
          ,v_ystag(sx(ne):ex(ne),sy(ne)-1:ey(ne),nzz) &
          ,ro_ystag(sx(ne):ex(ne),sy(ne)-1:ey(ne),nzz) )

 allocate( ro_zstag(sx(ne):ex(ne),sy(ne):ey(ne),1:nzz-1)  )


i02=ip2mem(ne)
do j=sy(ne)-1,ey(ne)+1
do i=sx(ne)-1,ex(ne)+1

! skip corner four points
if((j.eq.sy(ne)-1.or.j.eq.ey(ne)+1).and.(i.eq.sx(ne)-1.or.i.eq.ex(ne)+1)) cycle

  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
!  wgtfac(i,j)=hh-terrain(i02+ixy)
  ter2d(i,j)=terrain(i02+ixy)
  mpfac2d(i,j)=mpfac(i02+ixy)

!ter1_2d(i,j)=etater1(i02+ixy)
!ter2_2d(i,j)=etater2(i02+ixy)

  do k=1,nzz
    i03 = ip3mem(k,ne)
    u3d(i,j,k)=u(i03+ixy)
    v3d(i,j,k)=v(i03+ixy)
    w3d(i,j,k)=w(i03+ixy)
!    w3d(i,j,k)=dsdt3d(i03+ixy)
    p3d(i,j,k)=Plev(i03+ixy)
    t3d(i,j,k)=t(i03+ixy)
    dx3d(i,j,k)=dx(i03+ixy)
    dy3d(i,j,k)=dy(i03+ixy)
    dz3d(i,j,k)=dz(i03+ixy)
    den3d(i,j,k)=roair3d(i03+ixy)
    dsdtf(i,j,k)=dsdt3d(i03+ixy) ! defined at layer interface

!    jacobian(i,j,k) = 1.0 - ter2d(i,j)*( cosh((hh-etaz(k))/scale1) &
!                                        /(scale1*sinh(hh/scale1)) )
!    jacobian(i,j,k) = 1.0 - ter1_2d(i,j)*( cosh((hh-etaz(k))/scale1) &
!                                        /(scale1*sinh(hh/scale1)) ) &
!                          - ter2_2d(i,j)*( cosh((hh-etaz(k))/scale2) &
!                                        /(scale2*sinh(hh/scale2)) ) 

    jacobian(i,j,k) = pzps(i03+ixy)
!jacobian(i,j,k) = 1.0


    if(jacobian(i,j,k).le.0) then
      print*,'jacobian erro'
      print*,'i j  k ne =',i,j,k,ne
      print*,'jacobian=',jacobian(i,j,k)
    endif
  enddo
enddo
enddo

!print*,'dsdtf=',dsdtf(:,:,1)


do k=1,nzz
   do j=sy(ne)-1,ey(ne)+1
   do i=sx(ne)-1,ex(ne)+1
!     denair2d(i,j)=den3d(i,j,k)*wgtfac(i,j) ! multiply base term
     denair2d(i,j)=den3d(i,j,k) !*jacobian(i,j,k)
   enddo
   enddo
   do j=sy(ne),ey(ne)
   do i=sx(ne)-1,ex(ne)
     u_xstag(i,j,k)=0.5*(u3d(i,j,k)+u3d(i+1,j,k))
     ro_xstag(i,j,k)=0.5*(denair2d(i,j)+denair2d(i+1,j)) ! rhoair*pzps
!u_xstag(i,j,k)=0.0
   enddo
   enddo
   do j=sy(ne)-1,ey(ne)
   do i=sx(ne),ex(ne)
     v_ystag(i,j,k)=0.5*(v3d(i,j,k)+v3d(i,j+1,k))
     ro_ystag(i,j,k)=0.5*(denair2d(i,j)+denair2d(i,j+1)) ! rhoair*pzps
!v_ystag(i,j,k)=0.0
   enddo
   enddo
enddo

do k=1,nzz-1
do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)
! multiply base term   
    ro_zstag(i,j,k)=0.5*(den3d(i,j,k)+den3d(i,j,k+1))
!   ro_zstag(i,j,k)=0.5*(den3d(i,j,k)*jacobian(i,j,k)+den3d(i,j,k+1)*jacobian(i,j,k+1)) ! rhoair*pzps
enddo
enddo
enddo


dt_naqpms=dt


if(myid.eq.0) print*,'apm hadv'

!#################################
loop_hadv : do k=1,nzz-1

   do j=sy(ne)-1,ey(ne)+1
   do i=sx(ne)-1,ex(ne)+1
! multiply base term
!      denair2d(i,j)=den3d(i,j,k)*wgtfac(i,j)
     denair2d(i,j)=den3d(i,j,k) !*jacobian(i,j,k)
   enddo
   enddo

!=================

    allocate(jcbin(sx(ne)-1:ex(ne)+1))

    allocate(Q0(sx(ne)-1:ex(ne)+1))
    allocate(QN(sx(ne):ex(ne)))

    allocate(conc(sx(ne)-1:ex(ne)+1))
    allocate(conc11(sx(ne):ex(ne)))

    allocate(p1d(sx(ne)-1:ex(ne)+1))
    allocate(t1d(sx(ne)-1:ex(ne)+1))
    allocate(rhoa(sx(ne)-1:ex(ne)+1))

    allocate( DXX(sx(ne):ex(ne)) &
             ,DEN0(sx(ne):ex(ne)) &
             ,DEN1(sx(ne):ex(ne)) )
    allocate( u1D(sx(ne)-1:ex(ne)) &
             ,DD0(sx(ne)-1:ex(ne)) )

    allocate(vel(sx(ne)-1:ex(ne)+1))
    allocate(area(sx(ne)-1:ex(ne)+1))
    allocate(areav(sx(ne)-1:ex(ne)+1))


  loop_xadv : do j=sy(ne),ey(ne)


    do i=sx(ne)-1,ex(ne)+1
      jcbin(i)=jacobian(i,j,k)
    enddo

    tstep=dt_naqpms
    do i=sx(ne),ex(ne)
      ucell=u3d(i,j,k)
      dxcell = dx3d(i,j,k)/mpfac2d(i,j)
      ulft = u_xstag(i-1,j,k)
      urgt = u_xstag(i,j,k)

      tstep = amin1 ( tstep,cfl*dxcell/(abs(ulft)+1.0e-20),cfl*dxcell/(abs(urgt)+1.0e-20) )
    enddo
    nstep=int(dt_naqpms/tstep)+1
    tstep=dt_naqpms/nstep


!print*,'tstep=',tstep,nstep

!stop

!-----------------------------------------
! for unit transfer in advection algorism
    do i=sx(ne)-1,ex(ne)+1
      p1d(i)=p3d(i,j,k)
      t1d(i)=t3d(i,j,k)
      rhoa(i)=den3d(i,j,k)
!p1d(i)=1013.25
!t1d(i)=273.15
!rhoa(i)=1.0
    enddo
!-----------------------------------------

  if(trim(flagadv).eq.'walcek') then


    do i=sx(ne)-1,ex(ne)
       u1D(i)=u_xstag(i,j,k)
!u1D(i)=dx3d(i,j,k)/tstep
       DD0(i)=ro_xstag(i,j,k)
    enddo

!print*,'u1d=',u1d
!print*,'dd0=',dd0
!stop


    do i=sx(ne),ex(ne)
       DXX(i)=dx3d(i,j,k)/mpfac2d(i,j)
       DEN0(i)=denair2d(i,j)
       DEN1(i)=DEN0(i)-tstep/DXX(i)*( ro_xstag(i,j,k)*u_xstag(i,j,k) &
                                     -ro_xstag(i-1,j,k)*u_xstag(i-1,j,k) )*iroadv
    enddo

!print*,'dxx=',DXX
!print*,'den0=',den0
!print*,'den1=',den1

!stop

  elseif(trim(flagadv).eq.'camx_ppm') then

    dxres=dx3d(sx(ne),j,k) ! resolution

    do i=sx(ne)-1,ex(ne)
        vel(i)=u_xstag(i,j,k)
!vel(i)=-dxres/tstep
    enddo
    vel(ex(ne)+1)=9999.0 ! not used

    area=9999.0 ! not used 
!    do i=sx(ne)-1,ex(ne)+1
    do i=sx(ne),ex(ne)
       area(i)=mpfac2d(i,j)*mpfac2d(i,j)/(dy3d(i,j,k)*dz3d(i,j,k))
!area(i)=1
    enddo

    do i=sx(ne)-1,ex(ne)
       areav(i)=dy3d(i,j,k)*( dz3d(i,j,k)+dz3d(i+1,j,k)) &
                             /(mpfac2d(i,j)+mpfac2d(i+1,j) )
!areav(i)=1
    enddo
    areav(ex(ne)+1)=9999.0 ! not used 

    nndim=ex(ne)-sx(ne)+3

!print*,'dxres=',dxres
!print*,'vel=',vel
!print*,'area=',area
!print*,'areav=',areav
!stop

  endif

    do istep=1,nstep

        do is=1,NSO4
           iapm=ip_sulf(k,is,ne)
           mwgt=96.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=apm_sulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek( sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_sulf(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo
        enddo

        do is=1,NSEA
           iapm=ip_salt(k,is,ne)
           mwgt=96.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=apm_salt(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_salt(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo
        enddo

        do is=1,NDSTB
           iapm=ip_dust(k,is,ne)
           mwgt=96.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=apm_dust(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_dust(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo
        enddo

        do is=1,NBCOCT
           iapm=ip_bcoc(k,is,ne)
           mwgt=96.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=apm_bcoc(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_bcoc(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo
        enddo

           iapm=ip3mem(k,ne)
           mwgt=96.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=msltsulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             msltsulf(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo

           iapm=ip3mem(k,ne)
           mwgt=96.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=mdstsulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             mdstsulf(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo


           iapm=ip3mem(k,ne)
           mwgt=96.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=mbcsulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &       
                                ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             mbcsulf(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo

           iapm=ip3mem(k,ne)
           mwgt=96.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=mocsulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &       
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             mocsulf(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo

           iapm=ip3mem(k,ne)
           mwgt=98.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=h2so4_gas(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &       
                                 ,p1d,t1d,rhoa,mwgt,'ppbv' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ppbv' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             h2so4_gas(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo

 
        do is=1,nbincb
           iapm=ip_cbbin(k,is,ne)
           mwgt=1.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=apm_binbc(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &       
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_binbc(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo
        enddo

        do is=1,nbincb
           iapm=ip_cbbin(k,is,ne)
           mwgt=1.0
           do i=sx(ne)-1,ex(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(i)=apm_binoc(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sx(ne),ex(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                 ,p1d,t1d,rhoa,mwgt,'ugom3' )
           elseif(trim(flagadv).eq.'camx_ppm') then
             call hadvppm( nndim,tstep,dxres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3' )
           endif
           do i = sx(ne),ex(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_binoc(iapm+ixy)=amax1(QN(i),1.E-20)
           enddo
        enddo


      enddo ! nstep


  enddo  loop_xadv
!=================

    deallocate(jcbin)

    deallocate(Q0,QN,DXX,DEN0,DEN1,u1D,DD0)

    deallocate(conc,conc11,vel,area,areav)

    deallocate(p1d,t1d,rhoa)

    allocate( jcbin(sy(ne)-1:ey(ne)+1) )

    allocate( Q0(sy(ne)-1:ey(ne)+1) )
    allocate( QN(sy(ne):ey(ne)) )
    allocate( conc(sy(ne)-1:ey(ne)+1) )
    allocate( conc11(sy(ne):ey(ne)) )

    allocate( p1d(sy(ne)-1:ey(ne)+1) )
    allocate( t1d(sy(ne)-1:ey(ne)+1) )
    allocate( rhoa(sy(ne)-1:ey(ne)+1) )

    allocate( DXX(sy(ne):ey(ne)) &
             ,DEN0(sy(ne):ey(ne)) &
             ,DEN1(sy(ne):ey(ne)) )
    allocate( u1D(sy(ne)-1:ey(ne)) &
             ,DD0(sy(ne)-1:ey(ne)) )

    allocate(vel(sy(ne)-1:ey(ne)+1))
    allocate(area(sy(ne)-1:ey(ne)+1))
    allocate(areav(sy(ne)-1:ey(ne)+1))




!=================
  loop_yadv : do i=sx(ne),ex(ne)

    do j=sy(ne)-1,ey(ne)+1
      jcbin(j)=jacobian(i,j,k)
    enddo

    tstep=dt_naqpms
    do j=sy(ne),ey(ne)
      ucell=v3d(i,j,k)
      dxcell = dy3d(i,j,k)/mpfac2d(i,j)
!      tstep = amin1 ( tstep, dxcell/(abs(ucell)+1.0e-20) )
      ulft = v_ystag(i,j-1,k)
      urgt = v_ystag(i,j,k)
      
      tstep = amin1 ( tstep,cfl*dxcell/(abs(ulft)+1.0e-20),cfl*dxcell/(abs(urgt)+1.0e-20) )
    enddo
    nstep=int(dt_naqpms/tstep)+1
    tstep=dt_naqpms/nstep

!-----------------------------------------
! for unit transfer in advection algorism
    do j=sy(ne)-1,ey(ne)+1
      p1d(j)=p3d(i,j,k)
      t1d(j)=t3d(i,j,k)
      rhoa(j)=den3d(i,j,k)
    enddo
!------------------------------------------

if(1==2) then
print*,'p1d=',p1d
print*,'t1d=',t1d
print*,'rrhoa=',rhoa
endif


if(trim(flagadv).eq.'walcek') then


    do j=sy(ne)-1,ey(ne)
      u1D(j)=v_ystag(i,j,k)
      DD0(j)=ro_ystag(i,j,k)
    enddo
    do j=sy(ne),ey(ne)
      DXX(j)=dy3d(i,j,k)/mpfac2d(i,j)
      DEN0(j)=denair2d(i,j)-tstep/(dx3d(i,j,k)/mpfac2d(i,j))* &
                           (ro_xstag(i,j,k)*u_xstag(i,j,k)- &
                            ro_xstag(i-1,j,k)*u_xstag(i-1,j,k))*iroadv
      DEN1(j)=DEN0(j)   -tstep/(dy3d(i,j,k)/mpfac2d(i,j))* &
                        (ro_ystag(i,j,k)*v_ystag(i,j,k)- &
                         ro_ystag(i,j-1,k)*v_ystag(i,j-1,k))*iroadv
    enddo

if(1==2) then
print*,'u1d=',u1d
print*,'DD0=',DD0
print*,'DXX=',DXX
print*,'DEN0=',DEN0
print*,'DEN1=',DEN1
endif

elseif(trim(flagadv).eq.'camx_ppm') then

    dyres=dy3d(i,sy(ne),k) ! resolution

    do j=sy(ne)-1,ey(ne)
       vel(j)=v_ystag(i,j,k)
    enddo
    vel(ey(ne)+1)=9999.0 ! not used

    area=9999.0 ! not used
!    do j=sy(ne)-1,ey(ne)+1
    do j=sy(ne),ey(ne)
       area(j)=mpfac2d(i,j)*mpfac2d(i,j)/(dx3d(i,j,k)*dz3d(i,j,k))
    enddo

    do j=sy(ne)-1,ey(ne)
       areav(j)=dx3d(i,j,k)*(dz3d(i,j,k)+dz3d(i,j+1,k)) &
                            /(mpfac2d(i,j)+mpfac2d(i,j+1))
    enddo
    areav(ey(ne)+1)=9999.0 ! not used

    nndim=ey(ne)-sy(ne)+3

endif


    do istep=1,nstep

        do is=1,NSO4
           mwgt=96.0
           iapm=ip_sulf(k,is,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=apm_sulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_sulf(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo
        enddo

        do is=1,NSEA
           mwgt=96.0
           iapm=ip_salt(k,is,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=apm_salt(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_salt(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo
        enddo

        do is=1,NDSTB
           mwgt=96.0
           iapm=ip_dust(k,is,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=apm_dust(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_dust(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo
        enddo

        do is=1,NBCOCT
           mwgt=96.0
           iapm=ip_bcoc(k,is,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=apm_bcoc(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_bcoc(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo
        enddo


           mwgt=96.0
           iapm=ip3mem(k,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=msltsulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             msltsulf(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo


           mwgt=96.0
           iapm=ip3mem(k,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=mdstsulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             mdstsulf(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo

           mwgt=96.0
           iapm=ip3mem(k,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=mbcsulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             mbcsulf(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo

           mwgt=96.0
           iapm=ip3mem(k,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=mocsulf(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep & 
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             mocsulf(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo

           mwgt=98.0
           iapm=ip3mem(k,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=h2so4_gas(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep & 
                                ,p1d,t1d,rhoa,mwgt,'ppbv')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ppbv')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             h2so4_gas(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo

        do is=1,nbincb
           mwgt=96.0
           iapm=ip_cbbin(k,is,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=apm_binbc(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_binbc(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo
        enddo

        do is=1,nbincb
           mwgt=96.0
           iapm=ip_cbbin(k,is,ne)
           do j=sy(ne)-1,ey(ne)+1
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             Q0(j)=apm_binoc(iapm+ixy)
           enddo
           if(trim(flagadv).eq.'walcek') then
             call advec1d_walcek(sy(ne),ey(ne),Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                                ,p1d,t1d,rhoa,mwgt,'ugom3')
           elseif(trim(flagadv).eq.'camx_ppm') then
              call hadvppm(nndim,tstep,dyres,Q0,vel,area,areav,QN &
                          ,p1d,t1d,rhoa,mwgt,'ugom3')
           endif
           do j = sy(ne),ey(ne)
             ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             apm_binoc(iapm+ixy)=amax1(QN(j),1.E-20)
           enddo
        enddo

      enddo ! nstep

  enddo  loop_yadv
!================

    deallocate(jcbin)

    deallocate(Q0,QN,DXX,DEN0,DEN1,u1D,DD0)

    deallocate(conc,conc11,vel,area,areav)

    deallocate(p1d,t1d,rhoa)

enddo loop_hadv
!###################################

!stop

!return

!###################################

if(myid.eq.0) print*,'apm vadv'


 sz=1
 ez=nzz-1

 nndim=nzz-1

 allocate(jcbin(sz-1:ez+1))

 allocate(Q0(sz-1:ez+1))
 allocate(QN(sz:ez))

 allocate(conc00(sz-1:ez+1))
 allocate(conc(sz:ez))

 allocate(p1d(sz-1:ez+1))
 allocate(t1d(sz-1:ez+1))
 allocate(rhoa(sz-1:ez+1))

!if(trim(flagadv).eq.'camx_ppm') then

 allocate(DXX(sz:ez))
 allocate(DEN0(sz:ez))
 allocate(DEN1(sz:ez))
 allocate(u1d(sz-1:ez))
 allocate(DD0(sz-1:ez))

 !allocate(denvadv(sz:ez))


!elseif(trim(flagadv).eq.'camx_ppm') then

 allocate(depth(sz:ez))
 allocate(entrn(sz:ez))
 allocate(dilut(sz:ez))

!endif


do j = sy(ne),ey(ne)
do i = sx(ne),ex(ne)
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

   do k=sz,ez+1
     jcbin(k)=jacobian(i,j,k)
   enddo
   jcbin(0)=jcbin(1)

   do k=1,nzz
     p1d(k)=p3d(i,j,k)
     t1d(k)=t3d(i,j,k)
     rhoa(k)=den3d(i,j,k)
   enddo
   p1d(0)=p1d(1)
   t1d(0)=t1d(1)
   rhoa(0)=rhoa(1)

if(1==2) then
print*,'p1d=',p1d
print*,'t1d=',t1d
print*,'rrhoa=',rhoa
endif


   vmin=1.0e10
   vmax=0.0

   tstep=dt_naqpms
   do k=1,nzz-1
     i03=ip3mem(k,ne)
      dxcell = dz3d(i,j,k)/jacobian(i,j,k)
      ucell  = dsdt3d(i03+ixy)

     if(k.gt.1) then
       i03_top=ip3mem(k,ne)
       i03_bot=ip3mem(k-1,ne)
       vbot=dsdt3d(i03_bot+ixy)
       vtop=dsdt3d(i03_top+ixy)
     else
       i03_top=ip3mem(k,ne)
       vtop=dsdt3d(i03_top+ixy)
       vbot=0.0
     endif

if(1==2) then
     if(k.gt.1) then
       i03_top=ip3mem(k,ne)
       i03_bot=ip3mem(k-1,ne)
       vbot=w(i03_bot+ixy)
       vtop=w(i03_top+ixy)
     else
       i03_top=ip3mem(k,ne)
       vtop=w(i03_top+ixy)
       vbot=0.0
     endif
endif




!     tstep = min ( tstep, dxcell/(abs(ucell)+1.0e-30) )
     tstep = min ( tstep, dxcell/(abs(vbot)+1.0e-30), dxcell/(abs(vtop)+1.0e-30)  )


     if(abs(ucell).ge.vmax) vmax=abs(ucell)
     if(abs(ucell).le.vmin) vmin=abs(ucell)
   enddo
   nstep=int(dt_naqpms/tstep)+1
   tstep=dt_naqpms/nstep

!   print*,'tstep=',tstep,nstep,vmin,vmax


if(trim(flagadv).eq.'walcek') then
   do k=1,nzz-1
     u1d(k)=dsdtf(i,j,k)
!u1d(k)=-100/tstep
   enddo
   u1d(0)=0.0 ! no transport between layer 0 and 1

if(1==2) then
   do k=1,nzz
     i03=ip3mem(k,ne)
     u1d(k-1)=w(i03+ixy)
   enddo
endif


   DD0(0)=ro_zstag(i,j,1)
   do k=1,nzz-1
      DD0(K)=ro_zstag(i,j,k)
   enddo

   do k=1,nzz
!      denair1d(k)=den3d(i,j,k)*wgtfac(i,j)
      denair1d(k)=den3d(i,j,k) !*jacobian(i,j,k)
   enddo

   do  k=1,nzz-1
     DXX0=dx3d(i,j,k)/mpfac2d(i,j)
     DXX(k)=dz3d(i,j,k)/jacobian(i,j,k)
!DXX(k)=100
     DEN0(K)=denair1d(k)-tstep/DXX0*( ro_xstag(i,j,k)*u_xstag(i,j,k) &
                                   -ro_xstag(i-1,j,k)*u_xstag(i-1,j,k) )*iroadv &
                      -tstep/DXX0*( ro_ystag(i,j,k)*v_ystag(i,j,k) &
                                   -ro_ystag(i,j-1,k)*v_ystag(i,j-1,k) )*iroadv
     DEN1(K)=DEN0(K)-tstep/DXX(k)*( ro_zstag(i,j,k)*u1d(k) &
                                   -ro_zstag(i,j,k-1)*u1d(k-1) )*iroadv

!     denvadv(k)=-tstep/DXX(k)*( ro_zstag(i,j,k)*u1d(k) &
!                               -ro_zstag(i,j,k-1)*u1d(k-1) )
!print*,DEN1(k),denvadv(k),amax1(abs(u1d(k)),abs(u1d(k-1)))*tstep/DXX(k)
   enddo

if(trim(fvadv).eq.'implicit') then
! Fully implicit backward Euler method
   tstep=dt_naqpms
   nstep=1

   do k=sz,ez
     i03=ip3mem(k,ne)
     depth(k)=dz3d(i,j,k)/jacobian(i,j,k)
     entrn(k)=-dsdtf(i,j,k)
     dilut(k)=0.0
   enddo
endif


!print*,'entrn=',entrn


!if(myid.eq.0) print*

if(1==2) then
print*,'u1d=',u1d
print*,'DD0=',DD0
print*,'DXX=',DXX
print*,'DEN0=',DEN0
print*,'DEN1=',DEN1
endif


elseif(trim(flagadv).eq.'camx_ppm') then

! Fully implicit backward Euler method 
   tstep=dt_naqpms
   nstep=1

   do k=sz,ez
     i03=ip3mem(k,ne)
     depth(k)=dz3d(i,j,k)
     entrn(k)=entrn3d(i03+ixy)
     dilut(k)=dilut3d(i03+ixy)
!depth(k)=100
!entrn(k)=depth(k)/tstep
   enddo

endif

   do istep=1,nstep

     do is=1,NSO4
       mwgt=96.0
       do k=1,nzz
         iapm=ip_sulf(k,is,ne)
         Q0(k)=apm_sulf(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip_sulf(k,is,ne)
         apm_sulf(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo  
     enddo

     do is=1,NSEA
       mwgt=96.0
       do k=1,nzz
         iapm=ip_salt(k,is,ne)
         Q0(k)=apm_salt(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip_salt(k,is,ne)
         apm_salt(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo
     enddo

     do is=1,NDSTB
       mwgt=96.0
       do k=1,nzz
         iapm=ip_dust(k,is,ne)
         Q0(k)=apm_dust(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip_dust(k,is,ne)
         apm_dust(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo
     enddo

     do is=1,NBCOCT
       mwgt=96.0
       do k=1,nzz
         iapm=ip_bcoc(k,is,ne)
         Q0(k)=apm_bcoc(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip_bcoc(k,is,ne)
         apm_bcoc(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo
     enddo


       mwgt=96.0
       do k=1,nzz
         iapm=ip3mem(k,ne)
         Q0(k)=msltsulf(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip3mem(k,ne)
         msltsulf(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo

       mwgt=96.0
       do k=1,nzz
         iapm=ip3mem(k,ne)
         Q0(k)=mdstsulf(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip3mem(k,ne)
         mdstsulf(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo


       mwgt=96.0
       do k=1,nzz
         iapm=ip3mem(k,ne)
         Q0(k)=mbcsulf(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip3mem(k,ne)
         mbcsulf(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo

       mwgt=96.0
       do k=1,nzz
         iapm=ip3mem(k,ne)
         Q0(k)=mocsulf(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip3mem(k,ne)
         mocsulf(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo

       mwgt=98.0
       do k=1,nzz
         iapm=ip3mem(k,ne)
         Q0(k)=h2so4_gas(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ppbv' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ppbv')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ppbv')
       endif
       do k=1,nzz-1
         iapm=ip3mem(k,ne)
         h2so4_gas(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo

     do is=1,nbincb
       mwgt=96.0
       do k=1,nzz
         iapm=ip_cbbin(k,is,ne)
         Q0(k)=apm_binbc(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip_cbbin(k,is,ne)
         apm_binbc(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo
     enddo

     do is=1,nbincb
       mwgt=96.0
       do k=1,nzz
         iapm=ip_cbbin(k,is,ne)
         Q0(k)=apm_binoc(iapm+ixy)
       enddo
       Q0(0)=Q0(1)
       if(trim(flagadv).eq.'walcek') then
        if(trim(fvadv).eq.'explicit') then
         call advec1d_walcek(sz,ez,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,tstep &
                            ,p1d,t1d,rhoa,mwgt,'ugom3' )
        elseif(trim(fvadv).eq.'implicit') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
        endif
       elseif(trim(flagadv).eq.'camx_ppm') then
         call vrtslv(nndim,i,j,tstep,entrn,dilut,depth,Q0,QN &
                    ,p1d,t1d,rhoa,mwgt,'ugom3')
       endif
       do k=1,nzz-1
         iapm=ip_cbbin(k,is,ne)
         apm_binoc(iapm+ixy)=amax1(QN(k),1.0e-20)
       enddo
     enddo


   enddo ! istep


enddo
enddo
!###################################

    deallocate(jcbin)

    deallocate(Q0,QN,DXX,DEN0,DEN1,u1D,DD0)

    deallocate(conc,conc00,depth,dilut,entrn)

    deallocate(p1d,t1d,rhoa)

!    deallocate(denair1d)

    deallocate(u3d,v3d,w3d,p3d,t3d,dx3d,dy3d,dz3d,den3d,dsdtf)
    deallocate(mpfac2d,ter2d,ter1_2d,ter2_2d)
    deallocate(jacobian,denair2d,u_xstag,ro_xstag,v_ystag,ro_ystag,ro_zstag)

end




