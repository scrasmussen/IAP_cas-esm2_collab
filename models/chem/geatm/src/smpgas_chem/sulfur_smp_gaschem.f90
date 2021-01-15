
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  chenxsh@mail.iap.ac.cn 2016DEC
! (1) HO2 + HO2 -> H2O2 + O2
! (2) H2O2 + OH -> H2O + HO2
! (3) H2O2      -> 2OH

! (4) SO2 + OH  -> H2SO4
! (5) DMS + OH  -> SO2
! (6) DMS + OH  -> 0.5*SO2 + 0.5*HO2
! (7) DMS + NO3 -> SO2 + HNO3
! (8) NH3 + OH  -> H2O

subroutine sulfur_smp_gchembox( myid,dt_cbmz,ix,iy,iz &
                               ,iyear, imonth, iday, juday, iseason &
                               ,ihour, iminute &
                               ,tbeg_dd,tbeg_hh,tbeg_mm,tbeg_ss &
                               ,xlon,xlat,zalt_m &
                               ,nsp,noxdt,conc_oxdt,conc,daso4,te,pr_atm,rh &
                               ,coxdt)

 use smpsulf_var, only : idx_smpsulf_h2so4,idx_smpsulf_h2o2 &
                        ,idx_smpsulf_so2,idx_smpsulf_dms &
                        ,idx_smpsulf_nh3

 use smpsulf_var, only : idx_oxdt_oh,idx_oxdt_ho2,idx_oxdt_no3

 implicit none

 real   ,parameter :: avogad=6.02217e+23,deg2rad=0.017453293
 real   ,parameter :: cos_sza_cut = 0.017452406
 real   ,parameter :: ppb = 1.e+9

 real   ,parameter :: gocc2ugom3=1.0e12 ! g/cc -> ug/m3

 integer,parameter :: nract =9,nm12=12

 integer,intent(in) :: myid,ix,iy,iz
 real   ,intent(in) :: dt_cbmz


 integer,intent(in) :: iyear, imonth, iday, juday, iseason, ihour, iminute
 real   ,intent(in) :: tbeg_dd,tbeg_hh,tbeg_mm,tbeg_ss
 real   ,intent(in) :: xlon,xlat,zalt_m 
 real   ,intent(in) :: te,pr_atm,rh

 integer,intent(in) :: nsp,noxdt

 real   ,intent(in) :: conc_oxdt(noxdt,nm12)

 real    :: tbeg_hr
 integer :: ihr,ihr12

 real,dimension(nsp) :: conc,cwk,cwk0

 real,dimension(noxdt) :: coxdt

 real :: avedms,daso4,aveso2

 real :: cair_mlc
 real :: rlon,rlat
 real :: tcur_sec,dt_sec,tmid_sec
 real :: cos_sza

 real :: rk(nract)

 integer :: i
 integer :: ih2so4,ih2o2,iso2,idms,inh3
 integer :: ioh,iho2,ino3

 real*8 :: aa,bb,c0,tstep
 real :: twgt(24),curwgt,avewgt,avewgt12,twgt12(24),curwgt12

 real*8,external :: ct4ode,fave4ode

 real,external :: fun_rh2h2o

 real :: h2o

 rlon=xlon*deg2rad
 rlat=xlat*deg2rad

 cair_mlc = avogad*pr_atm/(82.056*te) ! air conc [molec/cc]

 h2o=fun_rh2h2o(rh, cair_mlc, te, pr_atm)

!rlon=120*deg2rad
!rlat=36*deg2rad

if(1==2) then
 print*,xlon,xlat
 print*,rlon,rlat,zalt_m,te,pr_atm
 print*,tbeg_dd,tbeg_hh,tbeg_mm,tbeg_ss
 print*,iyear, imonth, iday, juday, iseason
!stop
endif

! tmid_sec = time in seconds from Greenich Noon March 21
 dt_sec=dt_cbmz
 call updt_smpchem_time(tbeg_dd,tbeg_hh,tbeg_mm,tbeg_ss,tcur_sec,dt_sec,tmid_sec)

 call get_cos_sza(tmid_sec,rlon,rlat,cos_sza)

!==========================================================================
! calculate diurnal variation factors of radicals (OH HO2 NO3)
 do i=1,24
 
   ihr=i-1

   if(juday>=81.) then
      if(ihr>= 12) then
        tbeg_hr = ihr- 12
      else
        tbeg_hr = ihr- 12 +24
      endif
   else
      if(ihr> 12) then
        tbeg_hr = ihr- 12-24
      else
        tbeg_hr = ihr- 12
      endif
   endif

   dt_sec=0.0
   call updt_smpchem_time(tbeg_dd,tbeg_hr,0.0,0.0,tcur_sec,dt_sec,tmid_sec)
   call get_cos_sza(tmid_sec,rlon,rlat,twgt(i))
   twgt12(i)=-twgt(i)
   if(twgt(i).lt.cos_sza_cut) twgt(i)=0.0 
   if(twgt12(i).lt.cos_sza_cut) twgt12(i)=0.0
 enddo
 avewgt=sum(twgt(:))/24.0
 curwgt=twgt(int(ihour)+1)

 avewgt12=sum(twgt12(:))/24.0
 curwgt12=twgt12(int(ihour)+1)

 if(avewgt.eq.0) then
   avewgt=9999.0
   curwgt=0.0
 endif

 if(avewgt12.eq.0) then
   avewgt12=9999.0
   curwgt12=0.0
 endif

!===========================================================================


 call get_rk(pr_atm,te,zalt_m,cos_sza,nract,rk,ix,iy,iz)

 do i=1,nract
!   print*,i,rk(i),rk(i)*1.0e10
 enddo

 call interp_oxidant( iyear,imonth,iday,juday &
                     ,tmid_sec,rlon,rlat,cos_sza &
                     ,avewgt,curwgt,avewgt12,curwgt12 &
                     ,noxdt,nm12,conc_oxdt,coxdt)
 coxdt=coxdt/ppb*cair_mlc ! ppb -> molecules/cc

! print*,'oxidant'
 do i=1,noxdt
!   print*,coxdt(i),coxdt(i)*ppb/cair_mlc
 enddo
!stop

 ih2so4=idx_smpsulf_h2so4 ! 1
 ih2o2 =idx_smpsulf_h2o2  ! 2
 iso2  =idx_smpsulf_so2   ! 3
 idms  =idx_smpsulf_dms   ! 4
 inh3  =idx_smpsulf_nh3   ! 5

 ioh =idx_oxdt_oh  ! 1
 iho2=idx_oxdt_ho2 ! 2
 ino3=idx_oxdt_no3 ! 3

! ppb -> molecules/cc
!print*,'conc00',nsp
 do i=1,nsp
   cwk0(i)=conc(i)*cair_mlc/ppb
!print*,i,cwk0(i)
 enddo
 cwk=cwk0
!stop

 do i=1,nract
   if(i.eq.4) cycle
!   rk(i)=0.0
 enddo
!rk(3)=0.0
!rk(2)=1.8e-12

!===================================!
!                                   !
! dy/dt=aa+bb*y                     !
!                                   !
! y=((aa+bb*y0)*exp(bb*dt)-aa)/bb   !
!                                   !
!===================================!
!
!

 tstep=dt_cbmz

! d[H2O2]/dt=rk(1)*[HO2]**2-(rk(2)*[OH]+rk(3))*[H2O2]
 aa=(rk(1)+rk(9)*h2o)*coxdt(iho2)**2
 bb=-(rk(2)*coxdt(ioh)+rk(3))
 c0=cwk0(ih2o2)
 cwk(ih2o2)=ct4ode(aa,bb,tstep,c0)

if(.not.(cwk(ih2o2).ge.0.0.and.cwk(ih2o2).le.1.0e30)) then
      print*,'smpchem err H2O2',ih2o2
      print*,ix,iy,iz,cwk(ih2o2)
      print*,aa,bb,dt_cbmz,cwk0(ih2o2)
      print*,rk(1),rk(2),rk(3),coxdt(ioh)
      print*,cos_sza,exp(bb*tstep)
      stop
endif


! d[DMS]/dt=-(rk(5)*[OH]+0.5*rk(6)*[OH]+rk(7)*[NO3])*[DMS]
 aa=0.0
 bb=-( rk(5)*coxdt(ioh)+0.5*rk(6)*coxdt(ioh)+rk(7)*coxdt(ino3) )
 c0=cwk0(idms)
 cwk(idms)=ct4ode(aa,bb,tstep,c0)
 avedms=fave4ode(aa,bb,tstep,c0)

! d[SO2]/dt=(rk(5)*[OH]+0.5*rk(6)*[OH]+rk(7)*[NO3])*[DMS]-rk(4)*[OH]*[SO2]
! aa=( rk(5)*coxdt(ioh)+0.5*rk(6)*coxdt(ioh)+rk(7)*coxdt(ino3) )*cwk0(idms)
!avedms=0.0
 aa=( rk(5)*coxdt(ioh)+0.5*rk(6)*coxdt(ioh)+rk(7)*coxdt(ino3) )*avedms
 bb=-rk(4)*coxdt(ioh)
 c0=cwk0(iso2)
 cwk(iso2)=ct4ode(aa,bb,tstep,c0)
 aveso2=fave4ode(aa,bb,tstep,c0)

! d[H2SO4]/dt=rk(4)*[OH]*[SO2]
! aa=rk(4)*coxdt(ioh)*cwk0(iso2)
 aa=rk(4)*coxdt(ioh)*aveso2
 bb=0.0
 c0=cwk0(ih2so4)
 cwk(ih2so4)=ct4ode(aa,bb,tstep,c0)
 daso4=cwk(ih2so4)-cwk0(ih2so4)
! print*,daso4,cwk(iso2)-cwk0(iso2),cwk0(ih2so4),cwk(ih2so4)

! d[NH3]/dt=-rk(8)*[OH]*[NH3]
 aa=0.0
 bb=-rk(8)*coxdt(ioh)
 c0=cwk0(inh3)
 cwk(inh3)=ct4ode(aa,bb,tstep,c0)

 do i=1,nsp
!   print*,i,cwk0(i),cwk(i)
 enddo
!stop

! all sufuric acid condensed or nucleated
 daso4=daso4/avogad*96.0*gocc2ugom3 ! moles/cc -> ug/m3
 cwk(ih2so4)=0.0

!daso4=0.0

 do i=1,nsp
   conc(i)=cwk(i)*ppb/cair_mlc
!   print*,i,conc(i)
   if(.not.(conc(i).ge.0.0.and.conc(i).le.1.0e6)) then
      print*,'smpchem err'
      print*,ix,iy,iz,i,conc(i)
      stop
   endif
 enddo

 coxdt=coxdt*ppb/cair_mlc ! molecules/cc -> ppb

!stop

end subroutine sulfur_smp_gchembox
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



subroutine get_cos_sza(tmid_sec,rlon,rlat,cos_sza)

! tmid_sec = time in seconds from Greenich Noon March 21

  tlocal=tmid_sec
  tdec=0.4092797*sin(1.992385E-7*tlocal)
  sidec=sin(tdec)
  codec=cos(tdec)
  tloc=7.272205E-5*tlocal                                           
  thou=cos(rlon+tloc)
  cos_sza=sin(rlat)*sidec+cos(rlat)*codec*thou
  return
end

subroutine updt_smpchem_time(tbeg_dd,tbeg_hh,tbeg_mm,tbeg_ss,tcur_sec,dt_sec,tmid_sec)
  implicit none
  real,intent(in) :: tbeg_dd,tbeg_hh,tbeg_mm,tbeg_ss

  real :: tbeg_sec,tcur_sec,tsav_sec,told_sec,tcur_min,tcur_hrs,tmid_sec
  real :: dt_sec

  tbeg_sec = ((tbeg_dd*24.+tbeg_hh)*60.+tbeg_mm)*60.+tbeg_ss
  tcur_sec = tbeg_sec

  tsav_sec = tcur_sec
  told_sec = tcur_sec
  tcur_sec = tcur_sec + dt_sec
  tcur_min = tcur_sec/60.
  tcur_hrs = tcur_min/60.
  tmid_sec = told_sec + 0.5*dt_sec
end subroutine updt_smpchem_time


subroutine get_rk(pr_atm,te,zalt_m,cos_sza,nract,rk,ix,iy,iz)

 implicit none

 real,parameter :: avogad=6.02217e+23
 real,parameter :: cos_sza_cut = 0.017452406

 real :: cair_mlc
 real :: te,pr_atm,zalt_m

 integer :: nract

 real,dimension(nract) :: rk

 real,external :: arr4smpchem,troe4smpchem

 real :: zz,alpha0,alpha8,beta0,beta8,cos_sza
 real :: rk0,rnn,rki,rmm

 integer :: ix,iy,iz

 real :: rk_2ho2

! print*,pr_atm,te,zalt_m,cos_sza,nract
!stop

 cair_mlc = avogad*pr_atm/(82.056*te) ! air conc [molec/cc]

 rk(1)=2.3e-13*exp(600./te)+1.7e-33*exp(1000./te)*cair_mlc

 rk(2)=arr4smpchem(2.9e-12,-160.0,te)

if(cos_sza .ge. cos_sza_cut) then ! daytime
 zz=amin1(1.e-4*zalt_m,1.1)
 alpha0=2.540534e-9-4.e-11*(zz-0.75)**2.
 alpha8=-3.e-5+sqrt(alpha0)
 beta0=0.539284-0.16*(zz-1.75)**2.
 beta8=-1+sqrt(beta0)
 rk(3) = alpha8*exp(beta8/cos_sza)
else ! nighttime
 rk(3) =0.0
endif

 if(.not.(rk(3).ge.0.and.rk(3).le.1.0e30)) then
   print*,'rk3 err',rk(3)
   print*,ix,iy,iz
   print*,cos_sza,beta8,alpha8
   print*,zz,zalt_m,alpha0,beta0
!   stop
 endif


 rk0 = 3.0e-31
 rnn = 3.3
 rki = 1.5e-12
 rmm = 0.0
 rk(4)=troe4smpchem(cair_mlc,te,rk0,rnn,rki,rmm)

 rk(5)=arr4smpchem(9.60e-12,-234.0,te)

 rk(6)=1.7e-42*exp(7810./te)*cair_mlc*0.21/(1+5.5E-31*exp(7460.0/te)*cair_mlc*0.21)

 rk(7)=arr4smpchem(1.40e-13,500.0,te)

 rk(8)=arr4smpchem(1.70e-12,-710.0,te)

 rk_2ho2 = 2.3e-13 * exp(600./te) + & ! ho2 + ho2 --> h2o2
         & 1.7e-33 * exp(1000./te)*cair_mlc

 rk(9) = rk_2ho2*1.4e-21*exp(2200./te)


end subroutine get_rk


subroutine interp_oxidant( iyear,imonth,iday,juday &
                           ,tmid_sec,rlon,rlat,cos_sza &
                           ,avewgt,curwgt,avewgt12,curwgt12 &
                           ,noxdt,nm12,conc_oxdt,coxdt)
 implicit none
 integer,parameter :: pi=3.141592654
 integer :: nm12,noxdt
 integer :: iyear,imonth,iday,juday,iseason
 real :: tmid_sec,rlon,rlat,cos_sza
 real :: avewgt,curwgt,avewgt12,curwgt12
 real :: conc_oxdt(noxdt,nm12)
 real :: coxdt(noxdt),conc1(noxdt),conc2(noxdt)

 integer :: dayst(nm12),dayed(nm12)

 integer :: jd1,jd2,jdm
 logical :: leap_year

 integer :: idaymd,idayaf,idaybf,jdd,imbf,imaf

 integer :: i

!iyear=2016
!imonth=2
!iday=17
!call jday(iday,imonth,iyear, juday, iseason)

 jdd=juday

 dayst(1:nm12)=(/1 , 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1/)
 dayed(1:nm12)=(/31,28,31,30, 31,30,31,31, 30,31,30,31/)
 call leapy4(iyear,leap_year)
 if(leap_year) dayed(2)=29

 idaymd=0.5*(dayst(imonth)+dayed(imonth))

 if(iday.le.idaymd) then
   imbf=imonth-1
   if(imonth.eq.1) imbf=12
   idaybf=0.5*(dayst(imbf)+dayed(imbf))
   call jday(idaybf,imbf,iyear, jd1, iseason)
   call jday(idaymd,imonth,iyear, jd2, iseason)
   conc1(:)=conc_oxdt(:,imbf)
   conc2(:)=conc_oxdt(:,imonth)
   if(imonth.eq.1) then
     jd2=jd2+365
     jdd=juday+365
   endif
!print*,'bf'
!print*,imonth,iday,juday
!print*,imbf,idaybf,jd1,jdd,jd2
 else
   imaf=imonth+1
   if(imonth.eq.12) imaf=1
   idayaf=0.5*(dayst(imaf)+dayed(imaf))
   call jday(idaymd,imonth,iyear, jd1, iseason)
   call jday(idayaf,imaf,iyear, jd2, iseason)
   conc1(:)=conc_oxdt(:,imonth)
   conc2(:)=conc_oxdt(:,imaf)
   if(imonth.eq.12) then
     jd2=jd2+365
   endif

!print*,'af'
!print*,imonth,iday,juday
!print*,imaf,idayaf,jd1,jdd,jd2

 endif

!stop

 do i=1,noxdt
   coxdt(i)=(conc2(i)-conc1(i))*(jdd-jd1)/(jd2-jd1)+conc1(i)
   if(i.le.1.or.i.eq.2) then ! OH,HO2 diurnal variation 
     coxdt(i)=coxdt(i)*curwgt/avewgt
   elseif(i.eq.3) then ! NO3 diurnal variation
     coxdt(i)=coxdt(i)*curwgt12/avewgt12
   endif
!   print*,conc1(i),conc2(i),coxdt(i)
 enddo


end subroutine interp_oxidant


subroutine leapy4(year,leap_year)
 integer :: year
 logical :: leap_year
 if((mod(year,400).eq.0).or.(mod(year,4).eq.0.and.mod(year,100).ne.0)) then
   leap_year=.true.
 else
   leap_year=.false.
 endif
end subroutine leapy4


function arr4smpchem(aa,bb,te)
  implicit none
  real :: aa,bb,te
  real :: arr4smpchem
  arr4smpchem=aa*exp(bb/te)
end

function troe4smpchem(cair_mlc,te,rk0,rnn,rki,rmm)
  implicit none
  real :: cair_mlc,te,rk0,rnn,rki,rmm
  real :: troe4smpchem
  real :: expo
  rk0 = rk0*cair_mlc*(te/300.)**(-rnn)
  rki = rki*(te/300.)**(-rmm)
  expo= 1./(1. + (ALOG10(rk0/rki))**2)
  troe4smpchem  = (rk0*rki/(rk0+rki))*.6**expo
  return
end


function ct4ode(aa,bb,dt,c0)
  implicit none
  real*8 :: aa,bb,dt,c0
  real*8 :: ct4ode
  if(bb.ne.0) then
    ct4ode=((aa+bb*c0)*exp(bb*dt)-aa)/bb
  else
    ct4ode=c0+aa*dt
  endif
end

function fave4ode(aa,bb,dt,c0)
  implicit none
  real*8 :: aa,bb,dt,c0
  real*8 :: fave4ode
  if(bb.ne.0) then
    fave4ode=( (aa+bb*c0)*(exp(bb*dt)-1)/bb-aa*dt )/(bb*dt)
  else
    fave4ode=c0+0.5*aa*dt
  endif
end


function fun_rh2h2o(rh, cair_mlc, te, pr_atm)

      t_steam = 373.15                  ! steam temperature  [K]
      pr_std   = 1.0                    ! standard pressure  [atm]

      a      = 1.0 - t_steam/te
      arg    = (((-.1299*a -.6445)*a -1.976)*a +13.3185)*a
      pr_h2o = pr_std*exp(arg)                          ! [atm]
      fun_rh2h2o = RH*(pr_h2o/pr_atm)*cair_mlc/100.     ! [molec/cc] ! juanxiong he
       
      return
end




