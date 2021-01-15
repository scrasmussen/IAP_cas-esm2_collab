module work_vars

!real,allocatable,dimension(:) :: Z0, UST0 !,FSOIL,FICE, FSNOW,FVEG
real,allocatable,dimension(:) :: DUSTK,DUSTHGTF

real, allocatable, dimension(:) :: FSO4_DUST,FNO3_DUST, FSO4_SSA,FNO3_SSA

real,allocatable, dimension(:) :: kpmass_m1,RatioMass,kpmass_m2
                                      ! kpmass_m1 to mark the copy tracer
                                      ! kpmass_m2 to remember the real tracer

real, allocatable, dimension(:) :: gasbdysouth, gasbdynorth, &
                                   gasbdyeast,  gasbdywest

  ! for 5D variables 
real, allocatable, dimension(:) :: aerbdysouth, aerbdynorth, &
                                   aerbdyeast,  aerbdywest


integer, allocatable, dimension(:,:,:) :: ib4mem,jb4mem
integer,allocatable,dimension(:,:,:,:) :: ib5mem,jb5mem


real, allocatable, dimension(:) :: RK_HETSO2_DUST, RK_HETHNO3_DUST


 !===================================================
 ! For Source Mark
! integer,dimension(5) :: irecSM,ifsm
 real,allocatable,dimension(:,:,:)   :: sm,smb,smp,sm1 !sm1 is for ACM2
 real,allocatable,dimension(:,:,:,:) :: smconv  !convective diffusion
 real,allocatable,dimension(:,:)     :: smconv1 !moist convection 
 real,allocatable,dimension(:)       :: c00,smthis,smother
! real,allocatable,dimension(:)       :: tmpMarkCon
 real,allocatable,dimension(:)       :: contem0  !!!!!!!!!!!!!!
 real,allocatable,dimension(:)       :: TmpSM    !!!!!!!!!!!!!!
! integer,allocatable,dimension(:)    :: igMark  ! to mark gas   (idmSet)
! integer,allocatable,dimension(:)    :: iaMarkAer ! to mark aer type (idmSet)
! integer,allocatable,dimension(:)    :: iaMarkSiz ! to mark aer size (idmSet)

 ! for 5D variables (for O3/SO2/SO4/et al.) <idm,ism,i,j,k>
 real :: delta,delta1,delta2,delta3,delta4
 !=====================================================


!-------------------------------cloud and convection---------------
!the tempory variables
 real,allocatable,dimension(:,:,:) :: ppp,ttn,rkv,dzz,atm,ffn,conc,concmark
 real,allocatable,dimension(:)     :: T1,Q1,QS1,U1,V1,P1,PH1
 real,allocatable,dimension(:,:)   :: TRA1,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E
 real,allocatable,dimension(:)     :: CBMF1
!****** the cbmf(cloud base mass flux,must be remerbered in convect)

 integer, allocatable, dimension(:)   :: ip2mem2dconv
 integer, allocatable, dimension(:,:) :: ip3mem3dconv
 integer :: nxx,nyy ! the maximum grids in all domains

!----------------------------------dry------------
!the tempory variables
 real,allocatable,dimension(:,:,:)  :: height1,uu,vv,ww,ddx,ddz,water,QQ
 real,allocatable,dimension(:,:)    :: tsurf,xlat,xlon,SWDOWN1,pblhgt,vdep
 integer,allocatable,dimension(:,:) :: land

!----------------to define the layer of top conditions from global model
 real :: press
 real,allocatable,dimension(:)   :: ktop !the layer
 real,allocatable,dimension(:,:) :: kktop,tp !tempory variable
 real,allocatable,dimension(:)   :: tropp  ! the pressure of tropause
!---------------------------------------


  real, allocatable, dimension(:,:,:,:) :: ratioemit,ratioemitP,ratioemitL,ratioemitB
  real, allocatable, dimension(:)   :: ratioemitt,ratioemitPt


!  REAL    ::  cppb(ngas_max),cnnzifa(ngas_max)
  LOGICAL :: LSUNLIGHT                ! Flag for daytime
  INTEGER :: ISTEPI
  INTEGER :: JDATE           ! Current date (YYYYDDD)
  INTEGER :: JTIME           ! Current time (HHMMSS)
  REAL    :: FCLD          ! THE COFFI OF CLOUD TO PHOTO
  REAL    :: CLFLO1,CLFMI1,CLFHI1,TER


real,allocatable,dimension(:) :: deltso4_naq

contains 

subroutine allo_tmpwork_var(nx,ny,nzz,nest,sx,ex,sy,ey,mem_per_block &
                          ,igas,IAER,ISIZE,ndustcom,nseacom &
                          ,ifsmt,ifsm,idmSet,ismMax,igMark &
                          ,imasskeep )
 
 implicit none


 integer :: ne,is,k,ig,ia,iduc,idm,ism
 integer :: ii,mem5d,mem4d,mem3d,mem2d,mem2dgas,mem3daer

 integer :: mem5dc

 integer :: nest
 integer :: nx(5),ny(5),nzz
 integer :: sx(5), ex(5), sy(5), ey(5)

 integer :: mem_per_block(5)

 integer :: igas

 integer :: IAER,ISIZE

 integer :: ndustcom,nseacom


 integer :: ifsmt
 integer :: ifsm(5)
 integer :: idmSet,ismMax
 integer :: igMark(idmSet)

 integer :: imasskeep

 integer :: mem2dconv,mem3dconv

 integer :: mem4bi,mem4bj

 integer :: mem5bi,mem5bj

! integer :: ip2mem2dconv(nest),ip3mem3dconv(nzz,nest)




  mem2d=0         ! number of memory for 2d variables
  do ne=1,nest
    mem2d=mem2d+mem_per_block(ne)
  enddo


  allocate( ktop(mem2d) &
           ,tropp(mem2d) &
           ,CBMF1(mem2d) )  ! the top layer of global model


  mem3d=0        ! number of memory for 3d variables
  do ne=1,nest
  do k=1,nzz
    mem3d=mem3d+mem_per_block(ne)
  enddo
  enddo

  allocate(DUSTK(mem3d),DUSTHGTF(mem3d)) ! DUST EMISSIONS

  allocate(RK_HETSO2_DUST(mem3d),RK_HETHNO3_DUST(mem3d))

  allocate(deltso4_naq(mem3d))

  if(imasskeep==1)then
    allocate(kpmass_m1(mem3d))
    allocate(kpmass_m2(mem3d))
    allocate(RatioMass(mem3d))
  endif

  mem2dconv=0
  do ne=1,nest
        mem2dconv=mem2dconv+(ex(ne)-sx(ne)+3)*(ey(ne)-sy(ne)+3)
  enddo

  allocate(ip2mem2dconv(nest))
  ii=1
  do ne=1,nest
      ip2mem2dconv(ne)=ii
      ii=ii+(ex(ne)-sx(ne)+3)*(ey(ne)-sy(ne)+3)
  enddo

  mem3dconv=0
  do ne=1,nest
     do k=1,nzz
         mem3dconv=mem3dconv+(ex(ne)-sx(ne)+3)*(ey(ne)-sy(ne)+3)
     enddo
  enddo

  allocate(ip3mem3dconv(nzz,nest))
  ii=1
  do ne=1,nest
    do k=1,nzz
        ip3mem3dconv(k,ne)=ii
        ii=ii+mem_per_block(ne)
    enddo
  enddo


  mem4bi=0        ! number of memory for north-south boundary of 4D variables 
  do ne=1,nest
  do ig=1,igas
  do k=1,nzz
  mem4bi=mem4bi+(nx(ne)+2)
  enddo
  enddo
  enddo
  allocate(gasbdynorth(mem4bi),gasbdysouth(mem4bi))

  mem4bj=0        ! number of memory for east-west boundary of 4D variables
  do ne=1,nest
  do ig=1,igas
  do k=1,nzz
  mem4bj=mem4bj+(ny(ne)+2)
  enddo
  enddo
  enddo
  allocate(gasbdywest(mem4bj),gasbdyeast(mem4bj))

  allocate(ib4mem(nzz,igas,nest))
  allocate(jb4mem(nzz,igas,nest))
  ii=1
  do ne=1,nest
  do ig=1,igas
  do k=1,nzz
  ib4mem(k,ig,ne)=ii
  ii=ii+(nx(ne)+2)
  enddo
  enddo
  enddo
  ii=1
  do ne=1,nest
  do ig=1,igas
  do k=1,nzz
  jb4mem(k,ig,ne)=ii
  ii=ii+(ny(ne)+2)
  enddo
  enddo
  enddo

  mem5bi=0        ! number of memory for north-south boundary of 5D variables 
  do ne=1,nest
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  mem5bi=mem5bi+(nx(ne)+2)
  enddo
  enddo
  enddo
  enddo
  allocate(aerbdynorth(mem5bi),aerbdysouth(mem5bi))



  mem5bj=0        ! number of memory for east-west boundary of 5D variables
  do ne=1,nest
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  mem5bj=mem5bj+(ny(ne)+2)
  enddo
  enddo
  enddo
  enddo
  allocate(aerbdywest(mem5bj),aerbdyeast(mem5bj))


  allocate(ib5mem(nzz,isize,iaer,nest))
  allocate(jb5mem(nzz,isize,iaer,nest))

  ii=1
  do ne=1,nest
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  ib5mem(k,is,ia,ne)=ii
  ii=ii+(nx(ne)+2)
  enddo
  enddo
  enddo
  enddo

  ii=1
  do ne=1,nest
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  jb5mem(k,is,ia,ne)=ii
  ii=ii+(ny(ne)+2)
  enddo
  enddo
  enddo
  enddo





end subroutine allo_tmpwork_var


end module work_vars

