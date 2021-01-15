
subroutine cal_inorganic_hgfac &
  & ( myid &
  &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
  &  ,igas,GC_MOLWT )

use naqpms_varlist, only: ip3mem,ip4mem,gas
use naqpms_varlist, only: hgfac

use naqpms_varlist, only: laerv1,laerv2
use naqpms_varlist, only: ip4mem_aer,aerom,awc3d
use naqpms_varlist, only: naerbin
use naqpms_varlist, only: idx_so4,idx_no3,idx_nh4,idx_bc,idx_oc

use met_fields, only: rh1,clw

implicit none
real,parameter :: denwater=1.0,deninsp=1.7 ! g/cm3
integer :: myid
real    :: dt_naqpms
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

real, dimension(igas) :: GC_MOLWT ! igas=102

integer :: ibin
integer :: i04_so4,i04_no3,i04_nh4,i04_bc,i04_oc

!
integer :: ixy,i03,i04,igsp,iasp,iapm
integer :: i,j,k,is

real :: kk1,kk2



! 1d variable
!====================================================================

real    :: mdry,mh2o,yspgf,mpp

integer :: i04_75,i04_76,i04_77,i04_78,i04_79 &
          ,i04_80,i04_81,i04_82,i04_83,i04_84 &
          ,i04_85,i04_86,i04_87,i04_88,i04_89 &
          ,i04_90,i04_91,i04_92,i04_93,i04_94 &
          ,i04_95,i04_96,i04_97,i04_98,i04_99 &
          ,i04_100,i04_101,i04_102

real :: celrh,celcwc



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!print*,'kk1=',kk1,'kk2=',kk2
!stop

loop_j : do j=sy(ne),ey(ne)
loop_i : do i=sx(ne),ex(ne)

 ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

 loop_k : do k=1,nzz-1

    i03 = ip3mem(k,ne)

    celrh=rh1(i03+ixy)
    celcwc=clw(i03+ixy)

if(laerv1) then
    i04_75 = ip4mem(k,75,ne)
    i04_76 = ip4mem(k,76,ne)
    i04_77 = ip4mem(k,77,ne)
    i04_78 = ip4mem(k,78,ne)

    i04_79 = ip4mem(k,79,ne) ! H+(AQ)

    i04_80 = ip4mem(k,80,ne) ! NAAQ 
    i04_81 = ip4mem(k,81,ne) ! NH4+(AQ)
    i04_82 = ip4mem(k,82,ne) ! Cl-(AQ)
    i04_83 = ip4mem(k,83,ne) ! SO4--(AQ)
    i04_84 = ip4mem(k,84,ne) ! HSO4(-AQ)
    i04_85 = ip4mem(k,85,ne) ! NO3-(AQ)
    i04_86 = ip4mem(k,86,ne) ! NaCl(S)
    i04_87 = ip4mem(k,87,ne) ! Na2SO4
    i04_88 = ip4mem(k,88,ne) ! NaNO3(S)
    i04_89 = ip4mem(k,89,ne) ! (NH4)2SO4

    i04_90 = ip4mem(k,90,ne) ! NH4NO3(S)
    i04_91 = ip4mem(k,91,ne) ! NH4CL(S)
    i04_92 = ip4mem(k,92,ne) ! H2SO4(S)
    i04_93 = ip4mem(k,93,ne) ! NH4HSO4(S)
    i04_94 = ip4mem(k,94,ne) ! NaHSO4(S)
    i04_95 = ip4mem(k,95,ne) ! (NH4)4H(SO4)2(s)

    i04_102 = ip4mem(k,102,ne) ! AH2O

! primary particles

    mpp = gas(i04_77+ixy) &
        + gas(i04_78+ixy) 

! inorganic secondary particles

    mdry = gas(i04_79+ixy) &
         + gas(i04_80+ixy) &
         + gas(i04_81+ixy) &
         + gas(i04_82+ixy) &
         + gas(i04_83+ixy) &
         + gas(i04_84+ixy) &
         + gas(i04_85+ixy) &
         + gas(i04_86+ixy) &
         + gas(i04_87+ixy) &
         + gas(i04_88+ixy) &
         + gas(i04_89+ixy) &
         + gas(i04_90+ixy) &
         + gas(i04_91+ixy) &
         + gas(i04_92+ixy) &
         + gas(i04_93+ixy) &
         + gas(i04_94+ixy) &
         + gas(i04_95+ixy) 

    mdry = mdry+mpp

    mh2o = gas(i04_102+ixy)
endif


if(laerv2) then
    mdry=0
    do ibin=1,naerbin
      i04_so4=ip4mem_aer(k,ibin,idx_so4,ne)
      i04_nh4=ip4mem_aer(k,ibin,idx_nh4,ne)
      i04_no3=ip4mem_aer(k,ibin,idx_no3,ne)
      i04_bc=ip4mem_aer(k,ibin,idx_bc,ne)
      i04_oc=ip4mem_aer(k,ibin,idx_oc,ne)
      mdry = mdry + aerom(i04_so4+ixy) &
                  + aerom(i04_nh4+ixy) &
                  + aerom(i04_no3+ixy) &
                  + aerom(i04_bc+ixy) &
                  + aerom(i04_oc+ixy)

    enddo

    mh2o = awc3d(i03+ixy)
endif

    if(mdry.gt.0) then
      yspgf = (1.0+mh2o/mdry*deninsp/denwater)**(1.0/3.0)
    else
      yspgf = 1.0
    endif

    if(k.eq.1.and.i.eq.33.and.j.eq.33.and..true.) then
!       print*,'yspgf33=',yspgf,mdry,mh2o
    endif

    if(.not.(yspgf.ge.1.0.and.yspgf.le.100) ) then
!      print*,' naq_i&j&k=',i,j,k,yspgf,mdry,mh2o,celrh,celcwc
      yspgf=1.0
    endif

    hgfac(i03+ixy) = yspgf

 enddo loop_k

enddo loop_i
enddo loop_j


end subroutine cal_inorganic_hgfac


