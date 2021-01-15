
subroutine apm_snstvy_so2_emit &
 & ( myid &
 &  ,lapm &
 &  ,iemittype,ig &
 &  ,ne,nx,ny,nzz,nest,sx,ex,sy,ey &
 &  ,igas &
 &  ,ip2memGas,mem2dgas &
 &  ,EmtaGas,EmtpGas,EmttGas,EmtbGas )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid
logical :: lapm

integer :: iemittype
integer :: i,j,k,ig,i02Gas

integer :: ixy

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas
integer,dimension(igas,nest) :: ip2memGas

integer :: mem2dgas
real,dimension(mem2dgas) :: EmtaGas,EmtpGas,EmttGas,EmtbGas


! 1:anthoropogenic 2:power plant 3:biomass burning 4:biogenic

if( ig.eq.18 .and. lfrac_so2_emit ) then

 !print*,'shun_kk_ig=',ig
 !print*,'lfrac_so2_emit=',lfrac_so2_emit
 !print*,'frc_so2_emit=',frc_so2_emit
 !stop 'sub_snstvy'

 i02Gas = ip2memGas(ig,ne)

 do j=sy(ne),ey(ne)
 do i=sx(ne),ex(ne)

   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

   EmtaGas(i02Gas+ixy) = EmtaGas(i02Gas+ixy)*frc_so2_emit
   EmtpGas(i02Gas+ixy) = EmtpGas(i02Gas+ixy)*frc_so2_emit
   EmttGas(i02Gas+ixy) = EmttGas(i02Gas+ixy)*frc_so2_emit
   EmtbGas(i02Gas+ixy) = EmtbGas(i02Gas+ixy)*frc_so2_emit

 enddo
 enddo

endif


end subroutine apm_snstvy_so2_emit


