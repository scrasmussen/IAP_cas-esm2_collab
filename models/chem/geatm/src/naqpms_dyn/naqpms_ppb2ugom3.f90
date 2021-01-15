
subroutine naqpms_ppb_ugom3 &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem3d &
 &  ,mem4d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,igasCBM &
 &  ,GC_MOLWT &
 &  ,flag )

use naqpms_varlist 
use met_fields, only : Plev,t

implicit none

integer :: myid

character(len=*) :: flag

real :: dt

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igasCBM
integer :: igas,iaer,isize,nseacom,ndustcom

integer :: i,j,k,ig

integer :: mem3d
!real,dimension(mem3d) :: Plev,t

integer :: ixy,i03,i04

integer :: mem4d

real :: GC_MOLWT(igas)

real :: trfac

if(flag.eq.'ppb_ugom3'.or.flag.eq.'ugom3_ppb') then

 do j=sy(ne),ey(ne)
 do i=sx(ne),ex(ne)
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    do k=1,nzz
      i03 = ip3mem(k,ne)
      do ig=1,igasCBM
        i04 = ip4mem(k,ig,ne)
        trfac=GC_MOLWT(ig)/( (0.08206*T(i03+ixy))/(Plev(i03+ixy)/1000.) )
        if(flag.eq.'ppb_ugom3') then
          gas(i04+ixy) = gas(i04+ixy)*trfac
        elseif(flag.eq.'ugom3_ppb') then
          gas(i04+ixy) = gas(i04+ixy)/trfac
        endif
      enddo
    enddo
 enddo
 enddo

else 

  print*,'in qpms_ppb2ugom3.f90'
  stop 'unit transfer flag erro '

endif

end subroutine naqpms_ppb_ugom3





