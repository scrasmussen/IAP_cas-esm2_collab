
subroutine apm_ic(myid,lapm,ne,nx,ny,nzz,nest,sy,ey,sx,ex,ip3mem)
use apm_varlist !, only : apm_sulf,apm_salt,apm_dust,apm_bcoc
implicit none
include 'apm_parm.inc'
integer :: myid
logical :: lapm
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)
integer :: i,j,k,is
integer :: iapm,ixy
integer :: ip3mem(nzz,nest)

integer :: imode


!apm_sulf=1.0d0

!return

!print*,'kk_init',myid,ne
!print*,'kk_dim ',sx(ne),ex(ne),sy(ne),ey(ne)


!=====================================================
! initialize radius growth factor
  do j = sy(ne)-1,ey(ne)+1
  do i = sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       iapm=ip3mem(k,ne)

       rgf_sulf(iapm+ixy)=1.0
       rgf_salt(iapm+ixy)=1.0
       rgf_dust(iapm+ixy)=1.0
       rgf_bc(iapm+ixy)=1.0
       rgf_oc(iapm+ixy)=1.0


       accu_pso4_so2(iapm+ixy)=0 
     enddo
  enddo
  enddo
!=====================================================




loop_so4 : do is=1,NSO4
!loop_so4 : do is=1,1

  do j = sy(ne)-1,ey(ne)+1
  do i = sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz-1
       iapm=ip_sulf(k,is,ne)
       if(ltest_phy) then
         apm_sulf(iapm+ixy)=zero00 ! test apm ph
       elseif(ltest_dyn) then
         if(k.eq.1) apm_sulf(iapm+ixy)=15.0d0
         if(k.eq.2) apm_sulf(iapm+ixy)=10.0d0
         if(k.eq.3) apm_sulf(iapm+ixy)=5.0d0
         if(k.ge.4) apm_sulf(iapm+ixy)=1.0d0
       endif
       !rw_sulf(iapm+ixy)=rd_sulf(is)
     enddo
  enddo
  enddo

  do k=1,nzz-1
    iapm=ip_sulf(k,is,ne)
    call apm_init_bdy( myid, apm_sulf(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
  enddo

enddo loop_so4

!print*,'kkradius'
!print*,rw_sulf

!stop 'initstop'


!return


loop_salt : do is=1,NSEA
  do j = sy(ne)-1,ey(ne)+1
  do i = sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz-1
       iapm=ip_salt(k,is,ne)
       if(ltest_phy) then
         apm_salt(iapm+ixy)=zero00
       elseif(ltest_dyn) then
         if(k.eq.1) apm_salt(iapm+ixy)=15
         if(k.eq.2) apm_salt(iapm+ixy)=10
         if(k.eq.3) apm_salt(iapm+ixy)=5
         if(k.ge.4) apm_salt(iapm+ixy)=1
       endif
       !rw_salt(iapm+ixy)=rd_salt(is)
     enddo
  enddo
  enddo
  do k=1,nzz-1
    iapm=ip_salt(k,is,ne)
    call apm_init_bdy( myid, apm_salt(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
  enddo
enddo loop_salt


loop_dust : do is=1,NDSTB
  do j = sy(ne)-1,ey(ne)+1
  do i = sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz-1
       iapm=ip_dust(k,is,ne)
       if(ltest_phy) then
         apm_dust(iapm+ixy)=zero00
       elseif(ltest_dyn) then
         if(k.eq.1) apm_dust(iapm+ixy)=15
         if(k.eq.2) apm_dust(iapm+ixy)=10
         if(k.eq.3) apm_dust(iapm+ixy)=5
         if(k.ge.4) apm_dust(iapm+ixy)=1
       endif
       !rw_dust(iapm+ixy)=rd_dust(is)
     enddo
  enddo
  enddo
  do k=1,nzz-1
    iapm=ip_dust(k,is,ne)
    call apm_init_bdy( myid, apm_dust(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
  enddo
enddo loop_dust


loop_bcoc : do is=1,NBCOCT
  do j = sy(ne)-1,ey(ne)+1
  do i = sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz-1
       iapm=ip_bcoc(k,is,ne)
       if(ltest_phy) then
         apm_bcoc(iapm+ixy)=zero00
       elseif(ltest_dyn) then
         if(k.eq.1) apm_bcoc(iapm+ixy)=15
         if(k.eq.2) apm_bcoc(iapm+ixy)=10
         if(k.eq.3) apm_bcoc(iapm+ixy)=5
         if(k.ge.4) apm_bcoc(iapm+ixy)=1
       endif 
       
     enddo
  enddo
  enddo
  do k=1,nzz-1
    iapm=ip_bcoc(k,is,ne)
    call apm_init_bdy( myid, apm_bcoc(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
  enddo
enddo loop_bcoc



loop_binbc : do is=1,nbincb
  do j = sy(ne)-1,ey(ne)+1
  do i = sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz-1
       iapm=ip_cbbin(k,is,ne)
       if(ltest_phy) then
         apm_binbc(iapm+ixy)=zero00
       elseif(ltest_dyn) then
         if(k.eq.1) apm_binbc(iapm+ixy)=15
         if(k.eq.2) apm_binbc(iapm+ixy)=10
         if(k.eq.3) apm_binbc(iapm+ixy)=5
         if(k.ge.4) apm_binbc(iapm+ixy)=1
       endif

     enddo
  enddo
  enddo
  do k=1,nzz-1
    iapm=ip_cbbin(k,is,ne)
    call apm_init_bdy( myid, apm_binbc(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
  enddo
enddo loop_binbc


loop_binoc : do is=1,nbincb
  do j = sy(ne)-1,ey(ne)+1
  do i = sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz-1
       iapm=ip_cbbin(k,is,ne)
       if(ltest_phy) then
         apm_binoc(iapm+ixy)=zero00
       elseif(ltest_dyn) then
         if(k.eq.1) apm_binoc(iapm+ixy)=15
         if(k.eq.2) apm_binoc(iapm+ixy)=10
         if(k.eq.3) apm_binoc(iapm+ixy)=5
         if(k.ge.4) apm_binoc(iapm+ixy)=1
       endif

     enddo
  enddo
  enddo
  do k=1,nzz-1
    iapm=ip_cbbin(k,is,ne)
    call apm_init_bdy( myid, apm_binoc(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
  enddo
enddo loop_binoc


!do imode=1,2
!  do j = sy(ne)-1,ey(ne)+1
!  do i = sx(ne)-1,ex(ne)+1
!     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
!     do k=1,nzz-1
!       iapm=ipmode_bcoc(k,imode,ne)
!       refw_bc(iapm+ixy)=ref1d_bcoc(imode)
!       refw_oc(iapm+ixy)=ref1d_bcoc(imode)
!     enddo
!  enddo
!  enddo
!enddo





! sulfate coated on seasalt
do j = sy(ne)-1,ey(ne)+1
do i = sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz-1
     iapm=ip3mem(k,ne)
     if(ltest_phy) then
       msltsulf(iapm+ixy)=zero00
     elseif(ltest_dyn) then
       if(k.eq.1) msltsulf(iapm+ixy)=15
       if(k.eq.2) msltsulf(iapm+ixy)=10
       if(k.eq.3) msltsulf(iapm+ixy)=5
       if(k.ge.4) msltsulf(iapm+ixy)=1
     endif
   enddo
enddo
enddo
do k=1,nzz-1
  iapm=ip3mem(k,ne)
  call apm_init_bdy( myid, msltsulf(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
enddo
!!!!!


! sulfate coated on dust
do j = sy(ne)-1,ey(ne)+1
do i = sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz-1
     iapm=ip3mem(k,ne)
     if(ltest_phy) then
       mdstsulf(iapm+ixy)=zero00
     elseif(ltest_dyn) then
       if(k.eq.1) mdstsulf(iapm+ixy)=15
       if(k.eq.2) mdstsulf(iapm+ixy)=10
       if(k.eq.3) mdstsulf(iapm+ixy)=5
       if(k.ge.4) mdstsulf(iapm+ixy)=1
     endif
   enddo
enddo
enddo
do k=1,nzz-1
  iapm=ip3mem(k,ne)
  call apm_init_bdy( myid, mdstsulf(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
enddo
!!!!!


! sulfate coated on bc
do j = sy(ne)-1,ey(ne)+1
do i = sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz-1
     iapm=ip3mem(k,ne)
     if(ltest_phy) then
       mbcsulf(iapm+ixy)=zero00
     elseif(ltest_dyn) then
       if(k.eq.1) mbcsulf(iapm+ixy)=15
       if(k.eq.2) mbcsulf(iapm+ixy)=10
       if(k.eq.3) mbcsulf(iapm+ixy)=5
       if(k.ge.4) mbcsulf(iapm+ixy)=1
     endif
   enddo
enddo
enddo
do k=1,nzz-1
  iapm=ip3mem(k,ne)
  call apm_init_bdy( myid, mbcsulf(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
enddo
!!!!!


! sulfate coated on oc
do j = sy(ne)-1,ey(ne)+1
do i = sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz-1
     iapm=ip3mem(k,ne)
     if(ltest_phy) then
       mocsulf(iapm+ixy)=zero00
     elseif(ltest_dyn) then
       if(k.eq.1) mocsulf(iapm+ixy)=15
       if(k.eq.2) mocsulf(iapm+ixy)=10
       if(k.eq.3) mocsulf(iapm+ixy)=5
       if(k.ge.4) mocsulf(iapm+ixy)=1
     endif
   enddo
enddo
enddo
do k=1,nzz-1
  iapm=ip3mem(k,ne)
  call apm_init_bdy( myid, mocsulf(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
enddo
!!!!!


! sulferic vapor
do j = sy(ne)-1,ey(ne)+1
do i = sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz-1
     iapm=ip3mem(k,ne)
     if(ltest_phy) then
       h2so4_gas(iapm+ixy)=zero00
     elseif(ltest_dyn) then
       if(k.eq.1) h2so4_gas(iapm+ixy)=15
       if(k.eq.2) h2so4_gas(iapm+ixy)=10
       if(k.eq.3) h2so4_gas(iapm+ixy)=5
       if(k.ge.4) h2so4_gas(iapm+ixy)=1
     endif
   enddo
enddo
enddo
do k=1,nzz-1
  iapm=ip3mem(k,ne)
  call apm_init_bdy( myid, h2so4_gas(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
enddo
!!!!!


! sulferic vapor
do j = sy(ne)-1,ey(ne)+1
do i = sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz-1
     iapm=ip3mem(k,ne)
     if(ltest_phy) then
       sulf_inair(iapm+ixy)=zero00
       sulf_incld(iapm+ixy)=zero00
     elseif(ltest_dyn) then
       if(k.eq.1) sulf_inair(iapm+ixy)=15
       if(k.eq.2) sulf_inair(iapm+ixy)=10
       if(k.eq.3) sulf_inair(iapm+ixy)=5
       if(k.ge.4) sulf_inair(iapm+ixy)=1
       if(k.eq.1) sulf_incld(iapm+ixy)=15
       if(k.eq.2) sulf_incld(iapm+ixy)=10
       if(k.eq.3) sulf_incld(iapm+ixy)=5
       if(k.ge.4) sulf_incld(iapm+ixy)=1
     endif
   enddo
enddo
enddo
do k=1,nzz-1
  iapm=ip3mem(k,ne)
  call apm_init_bdy( myid, sulf_inair(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
  call apm_init_bdy( myid, sulf_incld(iapm) &
                    & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne) )
enddo
!!!!!





return

end subroutine apm_ic


subroutine wr_apm_init(ne)
implicit none
integer :: ne
return

end subroutine wr_apm_init



subroutine apm_init_bdy(myid, c, sx, ex, sy, ey ,nx,ny)

 implicit none
 integer :: ne
 integer myid, sx, ex, sy, ey, nx, ny
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
 integer i, j

 if(sx==1)then         ! west boundary
    do j=sy,ey
      c(sx-1,j) = c(sx,j)
    enddo
 endif

 if(ex==nx)then         ! east boundary
    do j=sy,ey
      c(ex+1,j) = c(ex,j)
    enddo
 endif

 if(sy==1)then         ! south boundary
    do i=sx,ex
      c(i,sy-1) = c(i,sy)
    enddo
 endif

 if(ey==ny)then         ! north boundary
    do i=sx,ex
      c(i,ey+1) = c(i,ey)
    enddo
 endif

 return

end subroutine apm_init_bdy



subroutine apm_init_tbdy(ne)
implicit none
integer :: ne
return

end subroutine apm_init_tbdy

