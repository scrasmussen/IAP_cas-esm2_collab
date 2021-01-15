
subroutine get_apm_cacid &
 & ( myid &
 &  ,lapm &
 &  ,dt &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip3mem,mem3d &
 &  ,igas,ip4mem,mem4d &
 &  ,PA2ATM &
 &  ,gas,Plev,t)

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

real    :: dt
logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,ig

integer :: mem3d
real,dimension(mem3d) :: Plev,t

integer :: mem4d
real,dimension(mem4d) :: gas

integer :: ixy,i03,i04

integer :: igas
integer :: ip3mem(nzz,nest),ip4mem(nzz,igas,nest)

real :: PA2ATM


!print*,'cacid dims:',sx(ne),ex(ne),sy(ne),ey(ne)
!print*,'PA2ATM=',PA2ATM
!print*,'apm_avogad=',apm_avogad
!print*,'ppbunit',ppbunit

!IF(lapm) THEN ! apm flag

 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
!   do k=1,1
   do k=1,nzz
     ig=1
     i03 = ip3mem(k,ne)
     i04 = ip4mem(k,ig,ne)
     apm_pr_atm =  PA2ATM*Plev(i03+ixy)*100.
     apm_te = t(i03+ixy)
!print*
!print*,k,j,i
!print*,gas(i04+ixy)
!gas(i04+ixy)=0.1
!print*,apm_pr_atm,apm_te
     apm_cair_mlc = apm_avogad*apm_pr_atm/(82.056*apm_te)
!print*,apm_cair_mlc
     apm_cacid_old(i03+ixy) = gas(i04+ixy)* &
                            & apm_cair_mlc/ppbunit
!if(k.eq.1) print*,'acid',j,i,apm_cacid_old(i03+ixy)
!stop 'sub_acid'
     apm_cacid(i03+ixy) = apm_cacid_old(i03+ixy)
   enddo
 enddo
 enddo

!ENDIF ! apm flag

end subroutine get_apm_cacid


subroutine get_apm_pacid &
 & ( myid &
 &  ,lapm &
 &  ,dt,dt_cacid &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip3mem,mem3d &
 &  ,igas,ip4mem,mem4d &
 &  ,PA2ATM &
 &  ,gas,Plev,t)

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

real :: dt

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,ig

integer :: mem3d

real,dimension(mem3d) :: Plev,t

integer :: mem4d

real,dimension(mem4d) :: gas

integer :: ixy,i03,i04

integer :: igas
integer :: ip3mem(nzz,nest),ip4mem(nzz,igas,nest)

real :: PA2ATM

real :: dt_cacid

!IF(lapm) THEN ! apm flag

 !dt_cacid=dt

 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz
!    do k=1,1
!print*,k,j,i
     ig=1
     i03 = ip3mem(k,ne)
     i04 = ip4mem(k,ig,ne)
     apm_pr_atm =  PA2ATM*Plev(i03+ixy)*100.
!print*,'apm_pr_atm=',apm_pr_atm
     apm_te = t(i03+ixy)
!print*,'apm_te=',apm_te

     apm_cair_mlc = apm_avogad*apm_pr_atm/(82.056*apm_te)
!print*,''
     apm_cacid(i03+ixy) = gas(i04+ixy)* &
                        & apm_cair_mlc/ppbunit
!print*,'old-new',apm_cacid_old(i03+ixy),apm_cacid(i03+ixy)
     apm_pacid(i03+ixy)=(apm_cacid(i03+ixy)-apm_cacid_old(i03+ixy))/dt_cacid
     apm_pacid(i03+ixy)=amax1(apm_pacid(i03+ixy),0.0d0)
!print*,'result',apm_pacid(i03+ixy)
   enddo
 enddo
 enddo

!ENDIF ! apm flag

end subroutine get_apm_pacid

