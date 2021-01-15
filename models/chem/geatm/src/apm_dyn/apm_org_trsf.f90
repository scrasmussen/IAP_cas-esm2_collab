
subroutine apm_org_trsf &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex )

use naqpms_varlist, only: agtflag, ip3mem
use apm_varlist
implicit none

integer :: myid

real :: dt

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is_hl,is_hb
!integer :: k,is

integer :: ixy,i02,i03,iapm_hb,iapm_hl

real :: org_hb,org_hl,org_hb00

real :: tau_hb_ageing

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 tau_hb_ageing=tau_hb


!IF(lapm) THEN ! apm flag


 !> shun : apm bcoc hdif
  if(lapm_trhl) then

   !print*,'run org_trhl'

   loop_bc : do is_hb=5,6 ! hydrophobic BC

     is_hl=is_hb-4

     do k = 1,nzz-1

       i03=ip3mem(k,ne)

       do j = sy(ne),ey(ne)
       do i = sx(ne),ex(ne)

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

        if(trim(agtflag).eq.'apm3sp') then
          tau_hb_ageing=bcagt(i03+ixy)
        elseif(trim(agtflag).eq.'efolding') then
          tau_hb_ageing=tau_hb
        else
          print*,'agtflag err'
          stop
        endif

        iapm_hb=ip_bcoc(k,is_hb,ne)
        iapm_hl=ip_bcoc(k,is_hl,ne)

        org_hb00=apm_bcoc(iapm_hb+ixy)
        org_hl=apm_bcoc(iapm_hl+ixy)

        org_hb=org_hb00*exp(-dt/tau_hb_ageing) ! 0.9976879
        org_hl=org_hl+org_hb00*(1.0-exp(-dt/tau_hb_ageing))

        apm_bcoc(iapm_hb+ixy)=org_hb
        apm_bcoc(iapm_hl+ixy)=org_hl

        !if(k.eq.1.and.i.eq.30.and.j.eq.30) then
        if(k.eq.1.and..false.) then
          print*,'i&j',i,j
          if(i.eq.30.and.j.eq.30) then
            print*,'i&j',i,j
            print*,'org:',org_hb00,exp(-dt/tau_hb_ageing),org_hb
            print*,'org:',org_hl
          endif
        endif

       enddo
       enddo
     enddo

   enddo loop_bc

   loop_oc : do is_hb=7,8 ! hydrophobic OC

     is_hl=is_hb-4

     do k = 1,nzz-1

       i03=ip3mem(k,ne)

       do j = sy(ne),ey(ne)
       do i = sx(ne),ex(ne)

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

      
        if(trim(agtflag).eq.'apm3sp') then
          tau_hb_ageing=ocagt(i03+ixy)
        elseif(trim(agtflag).eq.'efolding') then
          tau_hb_ageing=tau_hb
        else
          print*,'agtflag err'
          stop
        endif


        iapm_hb=ip_bcoc(k,is_hb,ne)
        iapm_hl=ip_bcoc(k,is_hl,ne)

        org_hb00=apm_bcoc(iapm_hb+ixy)
        org_hl=apm_bcoc(iapm_hl+ixy)

        org_hb=org_hb00*exp(-dt/tau_hb_ageing) ! 0.9976879
        org_hl=org_hl+org_hb00*(1.0-exp(-dt/tau_hb_ageing))

        apm_bcoc(iapm_hb+ixy)=org_hb
        apm_bcoc(iapm_hl+ixy)=org_hl

        !if(k.eq.1.and.i.eq.30.and.j.eq.30) then
        if(k.eq.1.and..false.) then
          print*,'i&j',i,j
          if(i.eq.30.and.j.eq.30) then
            print*,'i&j',i,j
            print*,'org:',org_hb00,exp(-dt/tau_hb_ageing),org_hb
            print*,'org:',org_hl
          endif
        endif

       enddo
       enddo
     enddo

   enddo loop_oc


  endif
 !< shun : end of apm bcoc hdif

!ENDIF ! apm flag


end subroutine apm_org_trsf





