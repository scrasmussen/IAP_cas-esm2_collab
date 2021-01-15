
subroutine apm_boxphy_driver &
 & ( myid &
 &  ,lapm &
 &  ,dt_naqpms &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d &
 &  ,ip4mem,mem4d &
 &  ,igas,gas,GC_MOLWT &
 &  ,clw_3d,rh1)

use apm_varlist
use aqchem_varlist
implicit none
include 'apm_parm.inc'

real,parameter :: PA2ATM = 1.0/101325.0

integer :: myid
real    :: dt_naqpms
logical :: lapm
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,ig,idx,is,N,imode

integer               :: mem2d
real,dimension(mem2d) :: longicrs,latitcrs,land_use,PSFC

integer               :: mem3d
real,dimension(mem3d) :: Plev,t,rh1

integer               :: mem4d
real,dimension(mem4d) :: gas ! ppb(gas),ug/m3(aerosol)

real, dimension(102)  :: GC_MOLWT

integer :: ixy,i02,i03,i04,iapm,i03_00,i04_acid,iccn

integer :: ip2mem(nest),ip3mem(nzz,nest),ip4mem(nzz,igas,nest)


!>-------------     end of definitions    ------------------<!
!=============================================================




!IF(lapm) THEN ! apm flag

 i02 = ip2mem(ne)

loop_k : DO k=1,nzz-1

  i03 = ip3mem(k,ne)
  i04 = ip4mem(k,ig,ne)

  loop_j : DO j = sy(ne),ey(ne)
  loop_i : DO i = sx(ne),ex(ne)

   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

  ENDDO loop_i
  ENDDO loop_j
ENDDO loop_k

!ENDIF ! apm flag


return

end subroutine apm_boxphy_driver






