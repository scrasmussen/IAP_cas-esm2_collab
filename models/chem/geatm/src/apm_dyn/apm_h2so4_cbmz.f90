

subroutine apm_h2so4_cbmz & 
                        & ( myid &
                        &  ,lapm &
                        &  ,dt,dt_cbmz &
                        &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
                        &  ,ip3mem,mem3d &
                        &  ,ip4mem,mem4d &
                        &  ,PA2ATM,Plev,temp &
                        &  ,igas,gas,GC_MOLWT &
                        &  ,flag )

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

integer :: mem3d,mem4d

real,dimension(mem3d) :: Plev,temp

real,dimension(mem4d) :: gas ! ppb,ug/m3

real, dimension(102)  :: GC_MOLWT

integer :: ip3mem(nzz,nest),ip4mem(nzz,igas,nest)

integer :: i03,i04,ixy

integer :: i,j,k,ig,idx,is

integer :: isulf

real              :: msulf,mwght
integer,parameter :: nsulf = 8
integer,parameter :: sulf_o_gas(1:nsulf) = (/ 1, 1, 1, 1, 1, 1, 1, 1/)
integer,parameter :: sulf_index(1:nsulf) = (/83,84,87,89,92,93,94,95/)


!dt_cbmz=dt

if(flag.eq.'11') then
! print*,'h2so4_gas before CBMZ'
elseif(flag.eq.'ch') then
! print*,'h2so4_gas and its production rate after CBMZ'
endif

loop_k : DO k=1,nzz-1

  i03 = ip3mem(k,ne)

  loop_j : DO j = sy(ne),ey(ne)
  loop_i : DO i = sx(ne),ex(ne)

    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

    if(flag.eq.'11') then
      p_h2so4_so2_cbmz(i03+ixy) = gas(ip4mem(k,1,ne)+ixy)
      if(k.eq.1) then
        !print*,ch_ratio(i03+ixy),msulf,j,i
      endif
    elseif(flag.eq.'ch') then
      p_h2so4_so2_cbmz(i03+ixy) = ( gas(ip4mem(k,1,ne)+ixy) - &
                                 & p_h2so4_so2_cbmz(i03+ixy) )/dt_cbmz
    else
      stop 'flag erro in apm_h2so4_cbmz.f90 '
    endif

  enddo loop_i
  enddo loop_j

enddo loop_k

!stop

end 






