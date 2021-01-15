
subroutine naqpms_bcoc_agt &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex )

use naqpms_varlist, only: agtflag
use naqpms_varlist, only: naerbin,ip3mem,ip4mem_aer,aerom
use apm_varlist, only: tau_hb,bcagt
implicit none

integer :: myid

real :: dt

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,ia_hl,ia_hb
integer :: is

integer :: ixy,i02,i03,iapm_hb,iapm_hl

real :: org_hb,org_hl,org_hb00

real :: tau_hb_ageing

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!return
! tau_hb_ageing=tau_hb


loop_size :  do is=1,naerbin

     do k = 1,nzz-1

       i03=ip3mem(k,ne)

       do j = sy(ne),ey(ne)
       do i = sx(ne),ex(ne)

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

        if(trim(agtflag).eq.'apm3sp') then
          tau_hb_ageing=bcagt(i03+ixy)
        elseif(trim(agtflag).eq.'efolding') then
          tau_hb_ageing=129600 !tau_hb
        else
          print*,'agtflag err'
          stop
        endif
!i04aer = ip4mem_aer(k,is,ia,ne)
!192               con(k)=aerom(i04aer+ixy)
        ia_hb=13
        ia_hl=14      

        iapm_hb=ip4mem_aer(k,is,ia_hb,ne)
        iapm_hl=ip4mem_aer(k,is,ia_hl,ne)

        org_hb00=aerom(iapm_hb+ixy)
        org_hl=aerom(iapm_hl+ixy)

        org_hb=org_hb00*exp(-dt/tau_hb_ageing) ! 0.9976879
        org_hl=org_hl+org_hb00*(1.0-exp(-dt/tau_hb_ageing))

        aerom(iapm_hb+ixy)=org_hb
        aerom(iapm_hl+ixy)=org_hl

        !if(k.eq.1.and.i.eq.30.and.j.eq.30) then
        if(k.eq.1.and..false.) then
          print*,'i&j',i,j
          if(i.eq.30.and.j.eq.30) then
            print*,'i&j',i,j
            print*,'org:',org_hb00,exp(-dt/tau_hb_ageing),org_hb
            print*,'org:',org_hl
          endif
        endif

       enddo ! i
       enddo ! j
     enddo   ! k

enddo loop_size


end subroutine naqpms_bcoc_agt





