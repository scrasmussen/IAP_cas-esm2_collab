
MODULE aqchem_varlist

logical,allocatable,dimension(:) :: lcloudy,lcloudy_old
logical,allocatable,dimension(:) :: lcld_app,lcld_dis,lcld_con
integer,allocatable,dimension(:) :: cld_flag
real   ,allocatable,dimension(:) :: clw_ph
real   ,allocatable,dimension(:) :: so4_aqchem

logical :: lupdt_met

contains

subroutine allo_aqchem(nx,ny,nzz,nest,sx,ex,sy,ey,mem_per_block)
 implicit none
 integer :: nx(5),ny(5),nzz,nest
 integer :: ne,k
 integer :: mem3d

 integer :: sx(5),ex(5),sy(5),ey(5)
 integer :: mem_per_block(5)
 

 mem3d=0
 do ne=1,nest
 do k=1,nzz
!   mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
   mem3d=mem3d+mem_per_block(ne)
 enddo
 enddo

 allocate( lcloudy(mem3d) )
 allocate( lcloudy_old(mem3d) )
 allocate( lcld_app(mem3d) )
 allocate( lcld_dis(mem3d) )
 allocate( cld_flag(mem3d) )
 allocate( clw_ph(mem3d) )
 allocate( lcld_con(mem3d) )
 allocate( so4_aqchem(mem3d) )

 lcloudy=.false.
 lcloudy_old=.false.
 lcld_app=.false.
 lcld_dis=.false.
 lcld_con=.false.
 cld_flag=0
 clw_ph=5.0
 so4_aqchem=0.0

end subroutine allo_aqchem

END MODULE aqchem_varlist
