
subroutine calc_air_density &
  & ( myid &
  &  ,ne,nx,ny,nzz,nest,sx,ex,sy,ey &
  &  ,mem3d )

use naqpms_varlist, only : ip3mem
use met_fields, only : plev,rh,t

implicit none
integer :: myid
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer               :: mem3d

!
integer :: ixy,i03,idx
integer :: i,j,k


integer,parameter  :: T0 = 273.16,EPS = 0.622

logical :: lmet_erro


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! roair : kg/m3
! PV=nRT 
! Pm/ro=(m/M)RT
!print*,ne
!print*,nx(ne),ny(ne)
!print*,sy(ne),ey(ne),sx(ne),ex(ne)
!print*,nest


loop_j : do j=sy(ne),ey(ne)
loop_i : do i=sx(ne),ex(ne)

 ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

 loop_k : do k=1,nzz

    i03 = ip3mem(k,ne)

    lmet_erro=.false.

    idx=i03+ixy

    if(    (.not.(temp(idx).gt.0.and.temp(idx).le.500)) &
       .or.(.not.(pres(idx).ge.0.and.pres(idx).le.200000)) ) then
       print*,'calc_denair'
       print*,'pres=',pres(idx),i,j,k
       print*,'temp=',temp(idx),i,j,k
       lmet_erro=.true.
    endif


    if(.not.lmet_erro) then
       roair(idx)=28.973*1.0E-3*pres(idx)*100/(8.31*temp(idx))
    else
       roair(idx)=1.0
       print*,'setdenair=',roair(idx)
    endif


!roair(idx)=1.0

    if( .not.(roair(idx).gt.0.and.roair(idx).lt.5) ) then
      write(*,'(a,2x,i2,1x,i2,1x,i2,3(1x,f10.3))'),'calro',i,j,k,roair(idx),temp(idx),pres(idx)
      stop
    endif


!roair(idx)=1.0

!  tmp1     = 10.*0.6112*exp(17.67*(TK-T0)/(TK-29.65))
!  tmp2     = EPS*tmp1/(0.01 * PRES -  (1.-EPS)*tmp1)
!  tmp1     = 100.*AMAX1(AMIN1(QV/tmp2,1.0),0.0)

!   if(k.eq.1) print*,ne,i,j,roair(idx)

    if(k.eq.1.and.i.eq.33.and.j.eq.33.and..false.) then
    endif

 enddo loop_k

enddo loop_i
enddo loop_j

!stop

end subroutine calc_air_density


