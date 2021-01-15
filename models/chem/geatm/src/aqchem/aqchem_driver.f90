
subroutine aqchem_driver &
  & ( myid &
  &  ,lapm &
  &  ,lnaqpms_pso4 &
  &  ,dt_cbmz &
  &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
  &  ,ip3mem,mem3d &
  &  ,ip4mem,mem4d &
  &  ,ip5mem,mem5d &
  &  ,igas,gas,GC_MOLWT &
  &  ,iaer,isize,aer &
  &  ,plev_3d,temp_3d,qvapor_3d,clw_3d &
  &  ,ANA,ASO4,ANH4,ANO3,ACL )

use naqpms_varlist, only : lgaschemsmp
use apm_varlist
use aqchem_varlist

use smpsulf_var, only : idx_smpsulf_h2so4,idx_smpsulf_h2o2 &
                       ,idx_smpsulf_so2,idx_smpsulf_dms &
                       ,idx_smpsulf_nh3
use smpsulf_var, only : idx_oxdt_oh,idx_oxdt_ho2,idx_oxdt_no3,idx_oxdt_o3
use smpsulf_var, only : ip4mem_oxdt,mmean_oxdt,nmoxdt,nm12
use smpsulf_var, only : ip4mem_ox3d,oxdt3d
use smpsulf_var, only : igassmp
use smpsulf_var, only : idx4dep_smpsulf

implicit none
include 'apm_parm.inc'
logical :: lapm
logical :: lnaqpms_pso4
integer :: myid
real    :: dt_cbmz
integer :: igas,iaer,isize
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer               :: mem3d
real,dimension(mem3d) :: plev_3d,temp_3d 
real,dimension(mem3d) :: qvapor_3d,clw_3d ! mixing ratio kg/kg
real,dimension(mem3d) :: ANA,ASO4,ANH4,ANO3,ACL
integer               :: mem4d
real,dimension(mem4d) :: gas ! ppb(gas),ug/m3(aerosol)

integer               :: mem5d
real,dimension(mem5d) :: aer

real, dimension(igas) :: GC_MOLWT ! igas=102
integer               :: ip3mem(nzz,nest),ip4mem(nzz,igas,nest) &
                        ,ip5mem(nzz,isize,iaer,nest)

!
integer :: ixy,i03,i04,igsp,iasp,iapm
integer :: i,j,k,is

integer :: ig_aq

! 1d variable
real,parameter :: ppb2mom=1.0e-9
!====================================================================
! CaCO3 and MgCO3 fraction in dust from the paper :
!  Global modeling of nitrate and ammonium : Interaction of aerosols
!  and tropospheric chemistry(Yan Feng,2005)
real,parameter :: fcaco3=0.07,fmgco3=0.055
real,parameter :: fnacl=1.0 ! NaCl mass fraction in sea salt
!====================================================================
integer :: i04_83,i04_84,i04_87,i04_89,i04_92,i04_93,i04_94
integer :: i04_81,i04_90,i04_91,i04_95
integer :: i04_85,i04_88
integer :: i04_86
integer :: i04_18,i04_02,i04_08,i04_04,i04_16,i04_11,i04_23,i04_01

real :: o2_pp,ch4_pp,h2_pp,atm_pp
real :: cwc,cliq
real :: pcell,tcell,pres_pa
real :: clw,cw_kgom3,cldph
real :: convfac
real :: r_gas(11),r_aer(9) ! mol/mol
real :: aer_mw(9)
real :: dt_aqchem,so4_old,nh4_old,no3_old
real :: so4_aqchem_1d,nh4_aqchem_1d,no3_aqchem_1d
real :: dust_tot,salt_tot
real :: caco3,mgco3,salt_nacl
real :: tamin,cwmin

! fCaCO3 = 0.065 : mass fraction of CaCO3 on dust

data    tamin /243.0/  ! K
data    cwmin /0.05/   ! g/m3
data    aer_mw / 96.0, 18.0, 62.0, 100.0, 84.0 &
               &,58.0, 56.0, 55.0, 74.0 /
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!print*,'dt_naqpms=',dt_naqpms
!stop

loop_j : do j=sy(ne),ey(ne)
loop_i : do i=sx(ne),ex(ne)

 ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

 loop_k : do k=1,nzz-1

  i03 = ip3mem(k,ne)

  i04_01 = ip4mem(k,1,ne)  ! H2SO4(g)
  i04_02 = ip4mem(k,2,ne)  ! HNO3
  i04_04 = ip4mem(k,4,ne)  ! NH3
  i04_08 = ip4mem(k,8,ne)  ! N2O5

  i04_11 = ip4mem(k,11,ne) ! O3
  i04_16 = ip4mem(k,16,ne) ! H2O2
  i04_18 = ip4mem(k,18,ne) ! SO2

  i04_23 = ip4mem(k,23,ne) ! FOA

  i04_81 = ip4mem(k,81,ne) ! NH4+(AQ)
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


  ! (1) liquid water is allowed to exist below 273 K
  ! (2) a ramp function is applied to apportion 
  !     total cloud water into liquid form between 243-273 K
  if(tcell.lt.273.0) then
    cliq = amax1( 0.0, cliq*(tcell-tamin)/(273.0-tamin) )
  endif

!-----Do RADM aqueous chemistry if CWC is above threshold
!     All conc units must be mol/mol (mixing ratio)
  if( cliq.ge.cwmin .and. tcell.ge.tamin ) then ! enough liquid cloud water ?

    pres_pa = pcell*100  ! hPa ->Pa
    cw_kgom3 = cliq/1000 

    dust_tot = 0.0
    do is=1,isize
      iapm = ip5mem(k,is,2,ne)
      dust_tot = dust_tot + aer(iapm+ixy)
    enddo
    caco3 = dust_tot*fcaco3 ! ug/m3
    mgco3 = dust_tot*fmgco3 ! ug/m3

    salt_tot = 0.0
    do is=1,isize
      iapm = ip5mem(k,is,1,ne)
      salt_tot = salt_tot + aer(iapm+ixy)
    enddo 
    salt_nacl = salt_tot*fnacl ! ug/m3


if(.not.lgaschemsmp) then
    r_gas(1)  = gas(i04_18+ixy)*ppb2mom ! SO2
    r_gas(2)  = gas(i04_02+ixy)*ppb2mom ! HNO3
    r_gas(3)  = gas(i04_08+ixy)*ppb2mom ! N2O5
    r_gas(4)  = 330.0*1.0e3*ppb2mom     ! CO2
    r_gas(5)  = gas(i04_04+ixy)*ppb2mom ! NH3
    r_gas(6)  = gas(i04_16+ixy)*ppb2mom ! H2O2
    r_gas(7)  = gas(i04_11+ixy)*ppb2mom ! O3
    r_gas(8)  = gas(i04_23+ixy)*ppb2mom ! FOA
    r_gas(9)  = 1.0e-3*ppb2mom          ! MHP : background value
    r_gas(10) = 1.0e-3*ppb2mom          ! PAA : background value
    r_gas(11) = gas(i04_01+ixy)*ppb2mom ! H2SO4(g)
else

    r_gas = 0.0

    ig_aq=idx4dep_smpsulf(idx_smpsulf_h2so4)
    i04_01 = ip4mem(k,ig_aq,ne)  ! H2SO4(g)
    r_gas(11) = gas(i04_01+ixy)*ppb2mom

!    i04_02 = ip4mem(k,ig_aq,ne)  ! HNO3

    ig_aq=idx4dep_smpsulf(idx_smpsulf_nh3)
    i04_04 = ip4mem(k,ig_aq,ne)  ! NH3
    r_gas(5)  = gas(i04_04+ixy)*ppb2mom

!    i04_08 = ip4mem(k,ig_aq,ne)  ! N2O5

    i04_11 = ip4mem_ox3d(k,idx_oxdt_o3,ne) ! O3
    r_gas(7)  = oxdt3d(i04_11+ixy)*ppb2mom

    ig_aq=idx4dep_smpsulf(idx_smpsulf_h2o2)
    i04_16 = ip4mem(k,ig_aq,ne) ! H2O2
    r_gas(6)  = gas(i04_16+ixy)*ppb2mom

    ig_aq=idx4dep_smpsulf(idx_smpsulf_so2)
    i04_18 = ip4mem(k,ig_aq,ne) ! SO2
    r_gas(1)  = gas(i04_18+ixy)*ppb2mom

!    i04_23 = ip4mem(k,ig_aq,ne) ! FOA

endif


    r_aer(1) = ASO4(i03+ixy)

    r_aer(2) = ANH4(i03+ixy)

    r_aer(3) = ANO3(i03+ixy)

    r_aer(4) = caco3 ! CaCO3 ???

    r_aer(5) = mgco3 ! MgCO3 ???

!    r_aer(6) = gas(i04_86+ixy)*35.5/GC_MOLWT(86) &
!             + gas(i04_91+ixy)*35.5/GC_MOLWT(91)     ! NaCl
    r_aer(6) = salt_nacl
    

    r_aer(7) = 0.01  ! Fe+++

    r_aer(8) = 0.005 ! Mn++

    r_aer(9) = 0.000 ! KCl

    ! reserve initial values for calculating variation quatity of aerosols
    so4_old=r_aer(1) ! ug/m3
    nh4_old=r_aer(2) ! ug/m3
    no3_old=r_aer(3) ! ug/m3 

    do iasp=1,9
      ! ug/m3 -> mol/mol
      r_aer(iasp) =  ( r_aer(iasp)/aer_mw(iasp)/convfac )*1.0e-6
    enddo

!!!!
   !:: treatment of pH value
   !%%======================================================%% 
   !%% the key to the question : the history of cloud water %%
   !%%======================================================%%

    if(lupdt_met) then ! new cloud water
      cldph = 5.0      ! camx, default or background value ?
    else               ! old cloud water
      cldph = clw_ph(i03+ixy)
    endif

    cldph = 5.0

!!!!

    if( r_aer(1).lt.5.0e-7 .and. r_aer(3).lt.1.5e-6 ) then
      ! check high acid conditions
      dt_aqchem = dt_cbmz
      call raqchem(tcell,pres_pa,dt_aqchem,cw_kgom3,r_gas,r_aer,cldph,i,j,k)
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !                                                                         !
    ! r_gas(i) is the Total concentration of species i in Air and Cloud Water !
    !                                                                         !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


    ! calculate sulfate production and transformation from gas to aerosol 
    r_aer(1) = r_aer(1)*1.0e6*convfac*aer_mw(1) ! mol/mol -> ug/m3
    so4_aqchem_1d = amax1( (r_aer(1)-so4_old), 0.0 )

    so4_aqchem(i03+ixy) = so4_aqchem_1d

    r_aer(2) = r_aer(2)*1.0e6*convfac*aer_mw(2) ! mol/mol -> ug/m3
    nh4_aqchem_1d = amax1( (r_aer(2)-nh4_old), 0.0 )

    r_aer(3)=r_aer(3)*1.0e6*convfac*aer_mw(3) ! mol/mol -> ug/m3
    no3_aqchem_1d = amax1( (r_aer(3)-no3_old), 0.0 )

    !===================================================================
    !> update tracers in NAQPMS due to aqueous chemistry
    !
    if(lnaqpms_pso4) then
      !(1) return the aqchem effect on gases , update gas concentration

if(lgaschemsmp) then
      gas(i04_18+ixy) = r_gas(1)/ppb2mom ! SO2
!      gas(i04_02+ixy) = r_gas(2)/ppb2mom ! HNO3
!      gas(i04_08+ixy) = r_gas(3)/ppb2mom ! N2O5
!      !                 r_gas(4)/ppb2mom ! CO2
      gas(i04_04+ixy) = r_gas(5)/ppb2mom ! NH3
      gas(i04_16+ixy) = r_gas(6)/ppb2mom ! H2O2
!      gas(i04_11+ixy) = r_gas(7)/ppb2mom ! O3
!      gas(i04_23+ixy) = r_gas(8)/ppb2mom ! FOA
      !                 r_gas(9)/ppb2mom ! NHP
      !                 r_gas(10)/ppb2mom !PAA
!      gas(i04_01+ixy) = r_gas(11)/ppb2mom ! H2SO4     
else
      gas(i04_18+ixy) = r_gas(1)/ppb2mom ! SO2
      gas(i04_02+ixy) = r_gas(2)/ppb2mom ! HNO3
      gas(i04_08+ixy) = r_gas(3)/ppb2mom ! N2O5
      !                 r_gas(4)/ppb2mom ! CO2
      gas(i04_04+ixy) = r_gas(5)/ppb2mom ! NH3
      gas(i04_16+ixy) = r_gas(6)/ppb2mom ! H2O2
      gas(i04_11+ixy) = r_gas(7)/ppb2mom ! O3
      gas(i04_23+ixy) = r_gas(8)/ppb2mom ! FOA
      !                 r_gas(9)/ppb2mom ! NHP
      !                 r_gas(10)/ppb2mom !PAA
      gas(i04_01+ixy) = r_gas(11)/ppb2mom ! H2SO4
endif

      !(2) return the aqchem effect on aerosols , update aerosols concentration
      ASO4(i03+ixy)=ASO4(i03+ixy)+so4_aqchem_1d
!      ANO3(i03+ixy)=ANO3(i03+ixy)+no3_aqchem_1d 
!      ANH4(i03+ixy)=ANH4(i03+ixy)+nh4_aqchem_1d 
    endif
    !===================================================================

    clw_ph(i03+ixy) = cldph

  else  ! liquid cloud water is not enough and temperature is too low

    ! (1) sulfate aqueous production is zero 
    ! (2) clean up values from last time-step
    so4_aqchem(i03+ixy)=0.0

  endif ! enough liquid cloud water ?

 enddo loop_k

enddo loop_i
enddo loop_j


!stop 'sub_aq'

end


