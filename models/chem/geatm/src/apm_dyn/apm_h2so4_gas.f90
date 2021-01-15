

subroutine apm_h2so4_gas & 
                         & ( myid &
                         &  ,lapm &
                         &  ,dt,dt_cbmz &
                         &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
                         &  ,ip3mem,mem3d &
                         &  ,PA2ATM,Plev,temp )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid
logical :: lapm
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

real :: dt,dt_cbmz

character :: flag*2

real    :: PA2ATM

integer :: mem3d

real,dimension(mem3d) :: Plev,temp


integer :: ip3mem(nzz,nest)

integer :: i03,i04,ixy

integer :: i,j,k,ig,idx,is


loop_k : DO k=1,nzz-1

  i03 = ip3mem(k,ne)

  loop_j : DO j = sy(ne),ey(ne)
  loop_i : DO i = sx(ne),ex(ne)

    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

    ! unit : ppb
    h2so4_gas(i03+ixy)=h2so4_gas(i03+ixy)+p_h2so4_so2_cbmz(i03+ixy)*dt 

    apm_pr_atm =  PA2ATM*Plev(i03+ixy)*100.
    apm_te = temp(i03+ixy)
    apm_cair_mlc = apm_avogad*apm_pr_atm/(82.056*apm_te)
   
    ! unit : #/(cm3 s)
!    apm_pacid(i03+ixy) = p_h2so4_so2_cbmz(i03+ixy)* &
!                       & apm_cair_mlc/ppbunit

    ! unit : #/cm3
!    apm_cacid(i03+ixy) = h2so4_gas(i03+ixy)* &
!                       & apm_cair_mlc/ppbunit

  enddo loop_i
  enddo loop_j

enddo loop_k


end subroutine apm_h2so4_gas 






