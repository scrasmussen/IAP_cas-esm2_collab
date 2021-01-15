
subroutine get_cld_flag &
  & ( myid &
  &  ,dt_naqpms &
  &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
  &  ,ip3mem,mem3d &
  &  ,qvapor_3d,clw_3d,rnw_3d )

use aqchem_varlist
implicit none

integer :: myid
real    :: dt_naqpms
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer               :: mem3d
real,dimension(mem3d) :: qvapor_3d,clw_3d,rnw_3d
integer               :: ip3mem(nzz,nest)
!
! 1d vars
real    :: qvapor,clw,rnw
!logical ::  
!
integer :: ixy,i03
integer :: i,j,k

!
real,parameter :: clw_min=0.0


loop_k : do k=1,nzz-1
loop_j : do j=sy(ne),ey(ne)
loop_i : do i=sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
  i03 = ip3mem(k,ne)

  if(clw.ge.clw_min) then
    lcloudy(i03+ixy)=.true.
    clw_ph(i03+ixy)=1
  endif

  if( lcloudy(i03+ixy).and.lcloudy_old(i03+ixy) ) then
    lcld_con(i03+ixy)=.true.
    lcld_app(i03+ixy)=.not.lcld_con(i03+ixy)
    lcld_dis(i03+ixy)=.not.lcld_con(i03+ixy)
  elseif( lcloudy(i03+ixy).and.( .not.lcloudy_old(i03+ixy)) ) then
    lcld_app(i03+ixy)=.true.
    lcld_con(i03+ixy)=.not.lcld_app(i03+ixy) 
    lcld_dis(i03+ixy)=.not.lcld_app(i03+ixy)
  elseif( (.not.lcloudy(i03+ixy)).and.lcloudy_old(i03+ixy) ) then
    lcld_dis(i03+ixy)=.true.
    lcld_app(i03+ixy)=.not.lcld_dis(i03+ixy)
    lcld_con(i03+ixy)=.not.lcld_dis(i03+ixy)
  endif

  lcloudy_old(i03+ixy)=lcloudy(i03+ixy)

enddo loop_i
enddo loop_j
enddo loop_k



end subroutine get_cld_flag



