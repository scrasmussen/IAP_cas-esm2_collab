
subroutine aqchem_driver &
  & ( myid &
  &  ,dt_naqpms &
  &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
  &  ,ip3mem,mem3d &
  &  ,ip4mem,mem4d &
  &  ,igas,gas,GC_MOLWT &
  &  ,plev_3d,temp_3d,qvapor_3d,clw_3d )

use apm_varlist
use aqchem_varlist
implicit none

integer :: myid
real    :: dt_naqpms
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer               :: mem3d
real,dimension(mem3d) :: plev_3d,temp_3d 
real,dimension(mem3d) :: qvapor_3d,clw_3d ! mixing ratio kg/kg
integer               :: mem4d
real,dimension(mem4d) :: gas ! ppb(gas),ug/m3(aerosol)
real, dimension(igas) :: GC_MOLWT ! igas=102
integer               :: ip3mem(nzz,nest),ip4mem(nzz,igas,nest)

!
integer :: ixy,i03,i04,igsp,iasp
integer :: i,j,k

! 1d variable
real,parameter :: ppb2mom=1.0e-9
integer :: i04_83,i04_84,i04_87,i04_89,i04_92,i04_93,i04_94
integer :: i04_81,i04_90,i04_91,i04_95
integer :: i04_85,i04_88
integer :: i04_86
real :: o2_pp,ch4_pp,h2_pp,atm_pp
real :: cwc,cliq
real :: pcell,tcell,pres_pa
real :: clw,cw_kgom3,cldph
real :: convfac
real :: r_gas(11),r_aer(9) ! mol/mol
real :: aer_mw(9)
real :: dt_aqchem,so4_old
!real :: gas_idx(11)
real :: tamin,cwmin

! fCaCO3 = 0.065 : mass fraction of CaCO3 on dust

data    tamin /243.0/  ! K
data    cwmin /0.05/   ! g/m3
!data    gas_idx //
data    aer_mw / 96.0, 18.0, 62.0, 100.0, 84.0 &
               &,58.0, 56.0, 55.0, 74.0 /
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!print*,'aq_plev=',plev_3d

!print*,'temp_3d=',temp_3d

!print*,'qvapor_3d=',qvapor_3d

!print*,'clw_3d=',clw_3d

!print*,'dt_naqpms=',dt_naqpms

!stop


loop_k : do k=1,nzz-1
!loop_k : do k=1,1
loop_j : do j=sy(ne),ey(ne)
loop_i : do i=sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
  i03 = ip3mem(k,ne)
  !i04 = ip4mem(k,ig,ne)
  tcell = temp_3d(i03+ixy)
  pcell = plev_3d(i03+ixy) ! hPa 
  clw   = clw_3d(i03+ixy)  ! kg/kg


  o2_pp  = 2.095e5
  ch4_pp = 1.75
  h2_pp  = 0.60
  atm_pp = 1.E06
  !convfac= 44.9 * (273./tcell)*(pcell/1013.) 
  convfac = pcell*100/(8.31*tcell)  ! mol air per m3

  cwc = clw*29.0*convfac
  cliq= cwc  ! g/m3

!if(cliq.gt.0.0) then
!print*,'kk_cwc',i,j,k
!print*,tcell,pcell,cwc,cliq
!endif


  ! (1) liquid water is allowed to exist below 273 K
  ! (2) a ramp function is applied to apportion 
  !     total cloud water into liquid form between 243-273 K
  if(tcell.lt.273.0) then
    cliq = amax1( 0.0, cliq*(tcell-tamin)/(273.0-tamin) )
  endif

  !cliq=0.1 ! for test
!if(cliq.gt.0.0) then
!print*,'kk_cwc',i,j,k
!print*,tcell,pcell,cwc,cliq
!endif


!-----Do RADM aqueous chemistry if CWC is above threshold
!     All conc units must be mol/mol (mixing ratio)
  if( cliq.ge.cwmin .and. tcell.ge.tamin ) then ! enough liquid cloud water ?

    pres_pa = pcell*100  ! hPa ->Pa
    cw_kgom3 = cliq/1000 

    igsp     = ip4mem(k,18,ne) ! SO2
    r_gas(1) = gas(igsp+ixy)*ppb2mom

    igsp     = ip4mem(k,2,ne) ! HNO3
    r_gas(2) = gas(igsp+ixy)*ppb2mom

    igsp     = ip4mem(k,18,ne) ! NxOy
    r_gas(3) = gas(igsp+ixy)*ppb2mom 

    r_gas(4) = 330.0*1.0e3*ppb2mom !  CO2

    igsp     = ip4mem(k,4,ne) ! NH3
    r_gas(5) = gas(igsp+ixy)*ppb2mom

    igsp     = ip4mem(k,16,ne) ! H2O2
    r_gas(6) = gas(igsp+ixy)*ppb2mom

    igsp     = ip4mem(k,11,ne) ! O3
    r_gas(7) = gas(igsp+ixy)*ppb2mom

    igsp     = ip4mem(k,23,ne) ! NxOy
    r_gas(8) = gas(igsp+ixy)*ppb2mom

    r_gas(9) = 1.0e-3*ppb2mom ! MHP : background value

    r_gas(10) = 1.0e-3*ppb2mom ! PAA : background value
     
    igsp     = ip4mem(k,1,ne) ! H2SO4(g)
    r_gas(11) = gas(igsp+ixy)*ppb2mom

    i04_83=ip4mem(k,83,ne) ! SO4--(AQ)
    i04_84=ip4mem(k,84,ne) ! HSO4(-AQ)
    i04_87=ip4mem(k,87,ne) ! Na2SO4
    i04_89=ip4mem(k,89,ne) ! (NH4)2SO4
    i04_92=ip4mem(k,92,ne) ! H2SO4(S)
    i04_93=ip4mem(k,93,ne) ! NH4HSO4(S)
    i04_94=ip4mem(k,94,ne) ! NaHSO4(S)
    i04_95=ip4mem(k,95,ne) ! (NH4)4H(SO4)2(s)
    r_aer(1) = gas(i04_83+ixy)*96.0/GC_MOLWT(83) &
             + gas(i04_84+ixy)*96.0/GC_MOLWT(84) &
             + gas(i04_87+ixy)*96.0/GC_MOLWT(87) &
             + gas(i04_89+ixy)*96.0/GC_MOLWT(89) &
             + gas(i04_92+ixy)*96.0/GC_MOLWT(92) &
             + gas(i04_93+ixy)*96.0/GC_MOLWT(93) &
             + gas(i04_94+ixy)*96.0/GC_MOLWT(94) &
             + gas(i04_95+ixy)*96.0/GC_MOLWT(95)*2.0 ! SO4


    i04_81=ip4mem(k,81,ne)  ! NH4+(AQ)
    !i04_89=ip4mem(k,89,ne) ! (NH4)2SO4
    i04_90=ip4mem(k,90,ne)  ! NH4NO3(S)
    i04_91=ip4mem(k,91,ne)  ! NH4CL(S)
    !i04_93=ip4mem(k,93,ne) ! NH4HSO4(S)
    !i04_95=ip4mem(k,95,ne) ! (NH4)4H(SO4)2(s)
    r_aer(2) = gas(i04_81+ixy)*18.0/GC_MOLWT(81) &
             + gas(i04_89+ixy)*18.0/GC_MOLWT(89)*2.0 &
             + gas(i04_90+ixy)*18.0/GC_MOLWT(90) &
             + gas(i04_91+ixy)*18.0/GC_MOLWT(91) &
             + gas(i04_93+ixy)*18.0/GC_MOLWT(93) &
             + gas(i04_95+ixy)*18.0/GC_MOLWT(95)*3.0 ! NH4

    i04_85=ip4mem(k,85,ne)  ! NO3-(AQ) 
    i04_88=ip4mem(k,88,ne)  ! NaNO3(S)
    !i04_90=ip4mem(k,90,ne)
    r_aer(3) = gas(i04_85+ixy)*62.0/GC_MOLWT(85) &
             + gas(i04_88+ixy)*62.0/GC_MOLWT(88) &
             + gas(i04_90+ixy)*62.0/GC_MOLWT(90)     ! NO3

    r_aer(4) = r_aer(1)/1.5 ! li CaCO3 ??? 
    r_aer(4) = 0.0

    i04_86=ip4mem(k,86,ne)  ! NaCl(S)
    !i04_91=ip4mem(k,91,ne) ! NH4Cl(S) 
    !r_aer(5) = r_aer(1)/12  ! li MgCO3 ???
    r_aer(5) = 0.0 ! camx

    r_aer(6) = gas(i04_86+ixy)*35.5/GC_MOLWT(86) &
             + gas(i04_91+ixy)*35.5/GC_MOLWT(91)     ! NaCl

    r_aer(7) = 0.01  ! Fe+++

    r_aer(8) = 0.005 ! Mn++

    r_aer(9) = 0.000 ! KCl

    so4_old=r_aer(1) ! ug/m3

    do iasp=1,9
      ! ug/m3 -> mol/mol
      r_aer(iasp) =  ( r_aer(iasp)/aer_mw(iasp)/convfac )*1.0e-6
    enddo

    cldph = clw_ph(i03+ixy)
!%%======================================================%% 
!%% the key to the question : the history of cloud water %%
!%%======================================================%%
    cldph = 5.0 ! camx, default or background value ?

    !print*,'run_aq',i,j,k
    !print*,'r_gas=',r_gas
    !print*,'r_aer=',r_aer

    if( r_aer(1).lt.5.0e-7 .and. r_aer(3).lt.1.5e-6 ) then
      ! check high acid conditions
!print*,'run_aq',i,j,k
      dt_aqchem = dt_naqpms
      call raqchem(tcell,pres_pa,dt_aqchem,cw_kgom3,r_gas,r_aer,cldph,i,j,k)
    endif

    ! return the effect of aqueous chemistry 
    r_aer(1)=r_aer(1)*1.0e6*convfac*aer_mw(1) ! mol/mol -> ug/m3
    so4_aqchem(i03+ixy) = amax1( (r_aer(1)-so4_old), 0.0 )

!print*,'pso4_aq=',so4_aqchem(i03+ixy)
!print*,'ph_valu=',cldph

    !????????
    ! return the effect of absorption of gases , update gas concentration

    ! return the sulfate mass from aqueous chemistry to ??? SO4--(aq) ?
    if(lnaqpms_pso4) then
      gas(i04_83+ixy)=gas(i04_83+ixy)+so4_aqchem(i03+ixy)
    endif

    clw_ph(i03+ixy) = cldph

  endif ! enough liquid cloud water ?

enddo loop_i
enddo loop_j
enddo loop_k


!stop 'sub_aq'

end


