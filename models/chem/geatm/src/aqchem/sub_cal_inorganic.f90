
subroutine cal_inorganic_aer &
  & ( myid &
  &  ,dt_naqpms &
  &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
  &  ,ip3mem,mem3d &
  &  ,ip4mem,mem4d &
  &  ,igas,gas,GC_MOLWT &
  &  ,ANA,ASO4,ANH4,ANO3,ACL )

use naqpms_varlist, only: laerv1,laerv2
use naqpms_varlist, only: naerbin,ip4mem_aer,aerom
use naqpms_varlist, only: idx_so4,idx_no3,idx_nh4

implicit none
integer :: myid
real    :: dt_naqpms
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer               :: mem3d
real,dimension(mem3d) :: ANA,ASO4,ANH4,ANO3,ACL
integer               :: mem4d
real,dimension(mem4d) :: gas ! ppb(gas),ug/m3(aerosol)
real, dimension(igas) :: GC_MOLWT ! igas=102
integer               :: ip3mem(nzz,nest),ip4mem(nzz,igas,nest)

!
integer :: ixy,i03,i04,igsp,iasp,iapm
integer :: i,j,k,is

real :: kk1,kk2


real,parameter :: wgt_so4=96,wgt_nh4=18,wgt_no3=62,wgt_na=23,wgt_cl=35.5

! 1d variable
!====================================================================
integer :: i04_83,i04_84,i04_87,i04_89,i04_92,i04_93,i04_94
integer :: i04_81,i04_90,i04_91,i04_95
integer :: i04_85,i04_88
integer :: i04_86,i04_80
integer :: i04_03
integer :: i04_18,i04_02,i04_08,i04_04,i04_16,i04_11,i04_23,i04_01
integer :: i04_82
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!print*,'kk1=',kk1,'kk2=',kk2
!stop

loop_j : do j=sy(ne),ey(ne)
loop_i : do i=sx(ne),ex(ne)

 ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

 loop_k : do k=1,nzz-1

    i03 = ip3mem(k,ne)

if(laerv1) then
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


    ANA(i03+ixy) = gas(i04_80+ixy)*wgt_na/GC_MOLWT(80)      &
                 + gas(i04_86+ixy)*wgt_na/GC_MOLWT(86)      &
                 + gas(i04_87+ixy)*wgt_na/GC_MOLWT(87)*2.0  &
                 + gas(i04_88+ixy)*wgt_na/GC_MOLWT(88)      &
                 + gas(i04_94+ixy)*wgt_na/GC_MOLWT(94)

    ASO4(i03+ixy) = gas(i04_83+ixy)*wgt_so4/GC_MOLWT(83) &
                  + gas(i04_84+ixy)*wgt_so4/GC_MOLWT(84) &
                  + gas(i04_87+ixy)*wgt_so4/GC_MOLWT(87) &
                  + gas(i04_89+ixy)*wgt_so4/GC_MOLWT(89) &
                  + gas(i04_92+ixy)*wgt_so4/GC_MOLWT(92) &
                  + gas(i04_93+ixy)*wgt_so4/GC_MOLWT(93) &
                  + gas(i04_94+ixy)*wgt_so4/GC_MOLWT(94) &
                  + gas(i04_95+ixy)*wgt_so4/GC_MOLWT(95)*2.0 ! SO4

    ANH4(i03+ixy) = gas(i04_81+ixy)*wgt_nh4/GC_MOLWT(81)     &
                  + gas(i04_89+ixy)*wgt_nh4/GC_MOLWT(89)*2.0 &
                  + gas(i04_90+ixy)*wgt_nh4/GC_MOLWT(90)     &
                  + gas(i04_91+ixy)*wgt_nh4/GC_MOLWT(91)     &
                  + gas(i04_93+ixy)*wgt_nh4/GC_MOLWT(93)     &
                  + gas(i04_95+ixy)*wgt_nh4/GC_MOLWT(95)*3.0 ! NH4

    ANO3(i03+ixy) = gas(i04_85+ixy)*wgt_no3/GC_MOLWT(85) &
                  + gas(i04_88+ixy)*wgt_no3/GC_MOLWT(88) &
                  + gas(i04_90+ixy)*wgt_no3/GC_MOLWT(90)     ! NO3

    ACL(i03+ixy)  = gas(i04_86+ixy)*wgt_cl/GC_MOLWT(86)  &
                  + gas(i04_91+ixy)*wgt_cl/GC_MOLWT(91)  &
                  + gas(i04_82+ixy)*wgt_cl/GC_MOLWT(82)

elseif(laerv2) then
    ANA(i03+ixy) = 0.0
    ASO4(i03+ixy) = 0.0
    ANH4(i03+ixy) = 0.0
    ANO3(i03+ixy) = 0.0
    ACL(i03+ixy)  = 0.0

    do is=1,naerbin

      ANA(i03+ixy) = ANA(i03+ixy) + 0.0

      i04=ip4mem_aer(k,is,idx_so4,ne)
      ASO4(i03+ixy) = ASO4(i03+ixy) + aerom(i04+ixy)

      i04=ip4mem_aer(k,is,idx_nh4,ne)
      ANH4(i03+ixy) = ANH4(i03+ixy) + aerom(i04+ixy)

      i04=ip4mem_aer(k,is,idx_no3,ne)
      ANO3(i03+ixy) = ANO3(i03+ixy) + aerom(i04+ixy)

      ACL(i03+ixy)  = ACL(i03+ixy)  + 0.0

    enddo

endif

    if(k.eq.1.and.i.eq.33.and.j.eq.33.and..false.) then
       print*,'calino'
       print*,'ANA =',ANA(I03+IXY)
       print*,'ASO4=',ASO4(I03+IXY)
       print*,'ANH4=',ANH4(I03+IXY)
       print*,'ANO3=',ANO3(I03+IXY)
       print*,'ACL =',ACL(I03+IXY)
    endif

 enddo loop_k

enddo loop_i
enddo loop_j


end subroutine cal_inorganic_aer


