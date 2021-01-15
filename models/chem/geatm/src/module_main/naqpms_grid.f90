module naqpms_gridinfo

real,allocatable,dimension(:) :: HGT1

real,allocatable,dimension(:) :: TERRAIN,LATITCRS,LONGICRS,LAND_USE

real,allocatable,dimension(:) :: dx,dy,dz,heiz

real,allocatable,dimension(:) :: mpfac

real,allocatable,dimension(:) :: pzps

contains

subroutine allo_grid_var(nx,ny,nzz,nest,sx,ex,sy,ey,mem_per_block)
 implicit none

 integer :: ne,is,k
 integer :: ii,mem_apm,mem3d,mem2d

 integer :: nest
 integer :: nx(5),ny(5),nzz
 integer :: sx(5), ex(5), sy(5), ey(5)

 integer :: mem_per_block(5)

!===================================================================


 ! number of memory for 2d variables

 mem2d=0
 do ne=1,nest
   mem2d=mem2d+mem_per_block(ne)
 enddo

 allocate( TERRAIN(mem2d) &
          ,LATITCRS(mem2d) &
          ,LONGICRS(mem2d) &
          ,LAND_USE(mem2d) &
          ,HGT1(mem2d) &
          ,mpfac(mem2d) )



 mem3d=0
 do ne=1,nest
 do k=1,nzz
    mem3d=mem3d+mem_per_block(ne)
 enddo
 enddo

 allocate( dx(mem3d) &
          ,dy(mem3d) &
          ,dz(mem3d) &
          ,heiz(mem3d) )

 allocate( pzps(mem3d) )


end subroutine allo_grid_var


end module naqpms_gridinfo
