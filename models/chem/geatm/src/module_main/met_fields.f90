module met_fields

real,allocatable,dimension(:) :: T2,PSFC,U10,V10,RAINCON,RAINNON,SWDOWN

real,allocatable,dimension(:) :: RHSFC

real,allocatable,dimension(:) :: TAUCLDI,TAUCLDC,clflo,clfmi,clfhi

real,allocatable,dimension(:) :: u,v,QVAPOR,clw,rnw,Plev,h,t,rh1

real,allocatable,dimension(:) :: SOILT,SOILRH,FICE,FSOIL,FVEG,FSNOW,UST0,Z0

real,allocatable,dimension(:) :: UST,RMOL,PBL_HGT

real,allocatable,dimension(:) :: tskwrf

real,allocatable,dimension(:) :: wrflai

real,allocatable,dimension(:) :: cldfrc3d,coefcld3d

real,allocatable,dimension(:) :: RAINSH

integer,allocatable,dimension(:) :: NPBL

real,allocatable,dimension(:) :: roair3d

!+++++++++

real,allocatable,dimension(:) :: w

real,allocatable,dimension(:) :: kh,kv

!real,allocatable,dimension(:) :: rkc_ysu

real,allocatable,dimension(:) :: phpt3d,dilut3d,entrn3d,dsdt3d


contains

subroutine allo_met_var(nx,ny,nzz,nest,sx,ex,sy,ey,mem_per_block,lrd_lai)
 implicit none

 integer :: ne,is,k
 integer :: ii,mem_apm,mem3d,mem2d

 integer :: nest
 integer :: nx(5),ny(5),nzz
 integer :: sx(5), ex(5), sy(5), ey(5)

 integer :: mem_per_block(5)

 logical :: lrd_lai

!===================================================================


 ! number of memory for 2d variables

 mem2d=0           
 do ne=1,nest
   mem2d=mem2d+mem_per_block(ne)
 enddo

 allocate( RAINCON(mem2d) &
          ,RAINNON(mem2d) &
          ,UST(mem2d) &
          ,U10(mem2d) &
          ,V10(mem2d) &
          ,T2(mem2d)  &
          ,SWDOWN(mem2d) &
          ,PSFC(mem2d) &
          ,PBL_HGT(mem2d) &
          ,RHSFC(mem2d) &
          ,RMOL(mem2d) &
          ,NPBL(mem2d) &
          ,clflo(mem2d) &
          ,clfmi(mem2d) &
          ,clfhi(mem2d) &
          ,RAINSH(mem2d) )

 ALLOCATE( FSOIL(MEM2D) &
          ,FICE(MEM2D) &
          ,FSNOW(MEM2D) &
          ,FVEG(MEM2D) &
          ,UST0(MEM2D) &
          ,Z0(MEM2D) )


 allocate( tskwrf(mem2d)  )

 if(lrd_lai) then
   allocate( wrflai(mem2d)  )
 endif

 mem3d=0        
 do ne=1,nest
 do k=1,nzz
    mem3d=mem3d+mem_per_block(ne)
 enddo
 enddo

 allocate( u(mem3d) &
          ,v(mem3d) &
          ,t(mem3d) &
          ,h(mem3d) &
          ,w(mem3d) &
          ,clw(mem3d) &
          ,rnw(mem3d) &
          ,rh1(mem3d) &
          ,Plev(mem3d) &
          ,QVAPOR(mem3d) &
          ,TAUCLDI(mem3d) &
          ,TAUCLDC(mem3d) &
          ,kh(mem3d) &
          ,kv(mem3d) )

! allocate(rkc_ysu(mem3d))

 allocate(cldfrc3d(mem3d))

 allocate(coefcld3d(mem3d))

 allocate(SOILT(mem3d),SOILRH(mem3d))

 allocate(roair3d(mem3d))

 allocate( phpt3d(mem3d) &
          ,dilut3d(mem3d) &
          ,entrn3d(mem3d) &
          ,dsdt3d(mem3d) )


end subroutine allo_met_var


end module met_fields
