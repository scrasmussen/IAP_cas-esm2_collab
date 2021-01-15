
subroutine apm_cal_opt &
 & ( myid &
 &  ,lapm &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d &
 &  ,dz,rh1,clw )

use apm_varlist,  only : ZBEXT3D,vsblt2d
use apm_varlist,  only : ip_ext
use apm_varlist,  only : ifhaze,iffog
use apm_varlist,  only : ifwrffog
use apm_varlist,  only : wrfvsb

use apm_varlist,  only : YBEXT3D
use apm_varlist,  only : ipwl_type 

use apm_varlist,  only : ZW3D

use apm_varlist,  only : AOD_390nm,AAOD_390nm
use apm_varlist,  only : AOD_500nm,AAOD_500nm
use apm_varlist,  only : AOD_530nm,AAOD_530nm
use apm_varlist,  only : AOD_550nm,AAOD_550nm
use apm_varlist,  only : AOD_700nm,AAOD_700nm
use apm_varlist,  only : AOD_1010nm,AAOD_1010nm

use apm_varlist,  only : ip_2dtype
use apm_varlist,  only : TAOD

implicit none

real,parameter :: rhow=1.0 ! g/cm3
real,parameter :: clwmin=0.01 ! g/kg  <-> 0.05g/m3 <-> 1500m vis
integer,parameter :: iwl_390nm=3
integer,parameter :: iwl_500nm=4
integer,parameter :: iwl_530nm=5
integer,parameter :: iwl_550nm=6
integer,parameter :: iwl_700nm=7
integer,parameter :: iwl_1010nm=8
integer :: myid
logical :: lapm
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,iwls,is

integer :: mem2d

integer               :: mem3d
real,dimension(mem3d) :: dz,rh1,clw

integer :: ixy,i02,i03,iapm,i00_tmp

integer :: i03_01,i03_02,i03_03,i03_k

integer :: ip2mem(nest),ip3mem(nzz,nest)

real    :: ZBEXT_z1,ZBEXT_z2,ZBEXT_z3
real    :: dz_z1,dz_z2,dz_z3

real    :: dz_k,ZBEXT_k,ZW_k
real    :: YBEXT

real    :: clw_z1,clw_z2,clw_z3

real    :: beta

real    :: rh_z1

real    :: BEXTL13

real    :: yvsb

real    :: aod390_1d,aod500_1d,aod530_1d,aod550_1d,aod700_1d,aod1010_1d
real    :: aaod390_1d,aaod500_1d,aaod530_1d,aaod550_1d,aaod700_1d,aaod1010_1d
real    :: taod_1d(5)

!================================================================
!      DATA (WLS(IWL),IWL=1,MWLS)/0.23,0.30,0.39,0.50,0.53,0.55,
!     &       0.70,1.01,1.27,1.46,1.78,2.05,2.33,2.79,3.46,8.02/
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



i02 = ip2mem(ne)

i03 = ip3mem(1,ne) ! k=1

loop_j : DO j = sy(ne),ey(ne)
loop_i : DO i = sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

  i03_01=ip3mem(1,ne); dz_z1=dz(i03_01+ixy)
  i03_02=ip3mem(2,ne); dz_z2=dz(i03_02+ixy)
  i03_03=ip3mem(3,ne); dz_z3=dz(i03_03+ixy)

  rh_z1=rh1(i03_01+ixy)

! 6 -> 550nm
  iapm=ip_ext(1,6,ne); ZBEXT_z1=ZBEXT3D(iapm+ixy) ! in cm-1
  iapm=ip_ext(2,6,ne); ZBEXT_z2=ZBEXT3D(iapm+ixy)
  iapm=ip_ext(3,6,ne); ZBEXT_z3=ZBEXT3D(iapm+ixy)

  BEXTL13=(ZBEXT_z1*dz_z1+ZBEXT_z2*dz_z2+ZBEXT_z3*dz_z3)/(dz_z1+dz_z2+dz_z3)

  yvsb=3.912/(3.9E-2+BEXTL13*1.0E5)  ! in km

  vsblt2d(i02+ixy)=yvsb  ! in km

  if(rh_z1.le.95.0.and.yvsb.le.5.0) then
    ifhaze(i02+ixy)=1
    iffog(i02+ixy) =0
  elseif(rh_z1.gt.95.0.and.yvsb.le.2.0) then
    ifhaze(i02+ixy)=0
    iffog(i02+ixy) =1
  else
    ifhaze(i02+ixy)=0
    iffog(i02+ixy) =0
  endif
    
  if(i.eq.33.and.j.eq.33.and..false.) then
    print*,'kk-vis'
    print*,'ZBEXT_z1=',ZBEXT_z1,dz_z1
    print*,'ZBEXT_z2=',ZBEXT_z2,dz_z2
    print*,'ZBEXT_z3=',ZBEXT_z3,dz_z3
    print*,'BEXTL13 =',BEXTL13
    print*,'vis33&33=',vsblt2d(i02+ixy)
  endif

  clw_z1=clw(i03+ixy)*1000.0 ! kg/kg -> g/kg
  if(clw_z1.gt.0) then
    beta=144.7*(rhow*clw_z1)**0.88
    wrfvsb(i02+ixy)=3.912/beta ! km
  else
    wrfvsb(i02+ixy)=9999.0
  endif  

  if(wrfvsb(i02+ixy).le.2.0) then
    ifwrffog(i02+ixy)=1
  else
    ifwrffog(i02+ixy)=0
  endif    

  aod390_1d=1.0e-20
  aod500_1d=1.0e-20
  aod530_1d=1.0e-20
  aod550_1d=1.0e-20
  aod700_1d=1.0e-20
  aod1010_1d=1.0e-20

  aaod390_1d=1.0e-20
  aaod500_1d=1.0e-20
  aaod530_1d=1.0e-20
  aaod550_1d=1.0e-20
  aaod700_1d=1.0e-20
  aaod1010_1d=1.0e-20

  taod_1d=1.0e-20

  loop_k : do k=1,nzz-1
    i03_k   = ip3mem(k,ne)
    dz_k    = dz(i03_k+ixy)

    iapm      = ip_ext(k,iwl_390nm,ne)
    ZBEXT_k   = ZBEXT3D(iapm+ixy)
    aod390_1d = aod390_1d + ZBEXT_k*100.0*dz_k
    ZW_k      = ZW3D(iapm+ixy)
    aaod390_1d = aaod390_1d + ZBEXT_k*100.0*dz_k*(1.0-ZW_k)

    iapm      = ip_ext(k,iwl_500nm,ne)
    ZBEXT_k   = ZBEXT3D(iapm+ixy)
    aod500_1d = aod500_1d + ZBEXT_k*100.0*dz_k
    ZW_k      = ZW3D(iapm+ixy)
    aaod500_1d = aaod500_1d + ZBEXT_k*100.0*dz_k*(1.0-ZW_k)

    iapm      = ip_ext(k,iwl_530nm,ne)
    ZBEXT_k   = ZBEXT3D(iapm+ixy) ! cm-1
    aod530_1d = aod530_1d + ZBEXT_k*100.0*dz_k
    ZW_k      = ZW3D(iapm+ixy)
    aaod530_1d = aaod530_1d + ZBEXT_k*100.0*dz_k*(1.0-ZW_k)

    iapm      = ip_ext(k,iwl_550nm,ne)
    ZBEXT_k   = ZBEXT3D(iapm+ixy)
    aod550_1d = aod550_1d + ZBEXT_k*100.0*dz_k
    ZW_k      = ZW3D(iapm+ixy)
    aaod550_1d = aaod550_1d + ZBEXT_k*100.0*dz_k*(1.0-ZW_k)

    iapm      = ip_ext(k,iwl_700nm,ne)
    ZBEXT_k   = ZBEXT3D(iapm+ixy)
    aod700_1d = aod700_1d + ZBEXT_k*100.0*dz_k
    ZW_k      = ZW3D(iapm+ixy)
    aaod700_1d = aaod700_1d + ZBEXT_k*100.0*dz_k*(1.0-ZW_k)

    iapm      = ip_ext(k,iwl_1010nm,ne)
    ZBEXT_k   = ZBEXT3D(iapm+ixy)
    aod1010_1d = aod1010_1d + ZBEXT_k*100.0*dz_k
    ZW_k      = ZW3D(iapm+ixy)
    aaod1010_1d = aaod1010_1d + ZBEXT_k*100.0*dz_k*(1.0-ZW_k)

    do is=1,5
      i00_tmp = ipwl_type(k,iwl_550nm,is,ne)
      YBEXT   = YBEXT3D(i00_tmp+ixy)
      taod_1d(is)=taod_1d(is)+YBEXT*100.0*dz_k
    enddo

  enddo loop_k

  AOD_390nm(i02+ixy)  = aod390_1d
  AOD_500nm(i02+ixy)  = aod500_1d
  AOD_530nm(i02+ixy)  = aod530_1d
  AOD_550nm(i02+ixy)  = aod550_1d
  AOD_700nm(i02+ixy)  = aod700_1d
  AOD_1010nm(i02+ixy) = aod1010_1d

  AAOD_390nm(i02+ixy)  = aaod390_1d
  AAOD_500nm(i02+ixy)  = aaod500_1d
  AAOD_530nm(i02+ixy)  = aaod530_1d
  AAOD_550nm(i02+ixy)  = aaod550_1d
  AAOD_700nm(i02+ixy)  = aaod700_1d
  AAOD_1010nm(i02+ixy) = aaod1010_1d

  do is=1,5
    i00_tmp = ip_2dtype (is, ne)
    TAOD(i00_tmp+ixy) = taod_1d(is)
  enddo

ENDDO loop_i
ENDDO loop_j

return

end subroutine apm_cal_opt






