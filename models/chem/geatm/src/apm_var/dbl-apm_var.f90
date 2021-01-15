
MODULE apm_varlist

! MBCOC8
! MBCOC8(1) : hydrophilic FF BC
! MBCOC8(2) : hydrophilic BB BC
! MBCOC8(3) : hydrophilic FF OC
! MBCOC8(4) : hydrophilic BB OC
! MBCOC8(5) : hydrophobic FF BC
! MBCOC8(6) : hydrophobic BB BC
! MBCOC8(7) : hydrophobic FF OC
! MBCOC8(8) : hydrophobic BB OC

! FF: 1,3,5,7  |  BB: 2,4,6,8

! BC: 1,2,5,6  |  OC: 3,4,7,8

! hydrophilic: 1-4  |  hydrophobic: 5-8



!include 'apm_parm.inc'
!!

! constant parameters
real,parameter :: kg2ug=1.0E+9 
real,parameter :: ug2kg=1.0E-9
real,parameter :: fbc_qw=0.5,fbc_sw=1.0-fbc_qw
real,parameter :: foc_qw=0.5,foc_sw=1.0-foc_qw
real,parameter :: fs2bc=0.5,fs2oc=1.0-fs2bc
real,parameter :: apm_Pi = 3.14159265358979

real,parameter :: densulf=1.7 ! g/cm3
real,parameter :: densalt=2.2 ! g/cm3
real,parameter :: dendust=2.5 ! g/cm3
real,parameter :: denbcoc=1.8 ! g/cm3

!
integer :: apmfunit
integer :: runned_mn

logical :: ltest_box=.false.

logical,allocatable,dimension(:) :: linit_box

logical :: lapm_emit,lapm_hadv,lapm_vadv,lapm_hdif,lapm_vdif,lapm_ccld
logical :: lapm_ddep,lapm_gdep
logical :: lapm_acid,lapm_phys
logical :: lapm_wdep,lapm_init,lapm_0emt,lapm_sulf_emit,lapm_salt_emit
logical :: lapm_dust_emit,lapm_4bdy,lapm_tbdy

logical :: lapm_isor

integer :: spinup_time

logical :: ls2bcoc,ls2sulf

logical :: lapm_wr00,lapm_wrfl

integer :: dep_method

logical :: lcoated_dyn

real :: frac_sp,frac_pp

logical :: lbox_phy
logical :: ltest_dyn,ltest_phy

integer :: ncycle

logical :: lapm_chbdy,lapm_nogdt

real :: zero00

real :: upper,lower

logical :: lapm_bpho,lapm_wreach


namelist /apm_lvar/ &
 &  lapm_emit,lapm_hadv,lapm_vadv,lapm_hdif,lapm_vdif    &
 & ,lapm_ccld,lapm_ddep,lapm_gdep,lapm_acid,lapm_phys    &
 & ,ls2sulf,lapm_wdep,lapm_init,lapm_0emt,lapm_sulf_emit &
 & ,lapm_salt_emit,lapm_dust_emit,lapm_4bdy,lapm_tbdy    &
 & ,dep_method,lapm_wrfl,spinup_time,frac_sp,lcoated_dyn &
 & ,lapm_wr00,lbox_phy,ltest_dyn,ltest_phy,lapm_isor     &
 & ,ncycle,lapm_chbdy,lapm_nogdt,zero00,upper,lower      &
 & ,lapm_bpho,lapm_wreach


logical :: lfor_sulf,lfor_salt,lfor_dust,lfor_bcoc
logical :: lem_sulf,lem_salt,lem_dust,lem_bcoc
namelist /apm_ldyn/ lfor_sulf,lfor_salt,lfor_dust,lfor_bcoc &
                   ,lem_sulf,lem_salt,lem_dust,lem_bcoc



integer :: sulf2dmem,salt2dmem,dust2dmem,bc2dmem,oc2dmem,bcoc2dmem
integer,allocatable,dimension(:,:) :: ip2d_sulf,ip2d_salt,ip2d_dust &
                                   & ,ip2d_bc,ip2d_oc,ip2d_bcoc

integer :: sulfmem,saltmem,dustmem,bcmem,ocmem,bcocmem

real*8, allocatable,dimension(:) :: sulf_area,salt_area,dust_area,bcoc_area,bc_area,oc_area

! particle tracers
real*8, allocatable,dimension(:) :: apm_sulf,apm_salt,apm_dust,apm_bcoc ! ug/m3
!

real*8, allocatable,dimension(:) :: sulf_sbdy,sulf_nbdy,sulf_wbdy,sulf_ebdy
real*8, allocatable,dimension(:) :: salt_sbdy,salt_nbdy,salt_wbdy,salt_ebdy
real*8, allocatable,dimension(:) :: dust_sbdy,dust_nbdy,dust_wbdy,dust_ebdy
real*8, allocatable,dimension(:) :: bcoc_sbdy,bcoc_nbdy,bcoc_wbdy,bcoc_ebdy


real*8, allocatable,dimension(:) :: sulf_emit,bcoc_emit,bc_emit,oc_emit &
                                   ,salt_emit,dust_emit

real*8, allocatable,dimension(:) :: fdstbin


real,allocatable,dimension(:)    :: ch_ratio


! coated sulfate in bins
real*8, allocatable,dimension(:) ::  &
&  apm_mbcsulf,apm_mocsulf,apm_mdstsulf,apm_msltsulf,apm_mbcocsulf !&
!& ,apm_mbcsoa ,apm_mocsoa ,apm_mdstsoa ,apm_msltsoa  &
!& ,apm_mbcnit ,apm_mocnit ,apm_mdstnit ,apm_msltnit  &
!& ,apm_mbcnh4 ,apm_mocnh4 ,apm_mdstnh4 ,apm_msltnh4  &
!& ,apm_mbcmsa ,apm_mocmsa ,apm_mdstmsa ,apm_msltmsa  


! coated tracers
real*8, allocatable,dimension(:) ::  &
&  mbcsulf,mocsulf,mdstsulf,msltsulf !&
!& ,mbcsoa ,mocsoa ,mdstsoa ,msltsoa  &
!& ,mbcnit ,mocnit ,mdstnit ,msltnit  &
!& ,mbcnh4 ,mocnh4 ,mdstnh4 ,msltnh4  &
!& ,mbcmsa ,mocmsa ,mdstmsa ,msltmsa


real*8, allocatable,dimension(:) :: bulk_msp,tot_sulf

integer,allocatable,dimension(:,:,:) :: ip_sulf,ip_salt,ip_dust,ip_bcoc
integer,allocatable,dimension(:,:,:) :: ipmode_bcoc

integer,allocatable,dimension(:,:,:) :: ipbsn_sulf,ipbsn_salt,ipbsn_dust,ipbsn_bcoc
integer,allocatable,dimension(:,:,:) :: ipbwe_sulf,ipbwe_salt,ipbwe_dust,ipbwe_bcoc
integer,allocatable,dimension(:,:)   :: ip2_sulf,ip2_salt,ip2_dust,ip2_bcoc

integer :: mem_apm

integer,allocatable,dimension(:) :: test_data

real*8, allocatable,dimension(:) :: apm_xq3d
real*8, allocatable,dimension(:) :: apm_cacid,apm_cacid_old,apm_pacid

real,parameter :: apm_avogad=6.02217e+23

real,parameter :: ppbunit=1.0e+9

!real :: dt_cacid

real :: apm_cair_mlc

real :: apm_pr_atm,apm_te
!>------------------------

real*8, allocatable,dimension(:) :: rd_sulf,rd_salt,rd_dust
real*8, allocatable,dimension(:) :: r_sulf,r_salt,r_dust,ref_bcoc
! sea salt flux ( kg/(m2 s) ) where u10 is 9m/s
real*8, allocatable,dimension(:) :: salt_flux 

real*8, allocatable,dimension(:) :: tn_bc,tn_oc,ref1d_bcoc

real*8, allocatable,dimension(:,:) :: sulf_unit_em2,bcoc_unit_em2
!TOTNUMBC,TOTNUMOC,REFBCOC,CEMITSULF2,CEMITBCOC2

real*8, allocatable,dimension(:,:) :: grfac
!>-----------------------------------------


real, allocatable,dimension(:) :: apmdryvel2d


!============================================
!> wet deposition variables

real, allocatable,dimension(:) :: apm_wdep_sulf,apm_wdep_salt,apm_wdep_dust &
                                   ,apm_wdep_bcoc,apm_wdep2_sulf &
                                   ,apm_wdep2_salt,apm_wdep2_dust,apm_wdep2_bcoc

real, allocatable,dimension(:) :: &
 apm_wdfld_sulf,apm_wdfld_salt,apm_wdfld_dust,apm_wdfld_bcoc &
,apm_wdfld2_sulf,apm_wdfld2_salt,apm_wdfld2_dust,apm_wdfld2_bcoc


! coated sulfate deposition in bulk
real, allocatable,dimension(:) :: wdep_sltsulf,wdep_dstsulf,wdep_bcsulf &
                                   ,wdep_ocsulf &
                                   ,wdep2_sltsulf,wdep2_dstsulf,wdep2_bcsulf &
                                   ,wdep2_ocsulf

real, allocatable,dimension(:) :: wdfld_sltsulf,wdfld_dstsulf,wdfld_bcsulf &
                                   ,wdfld_ocsulf &
                                   ,wdfld2_sltsulf,wdfld2_dstsulf &
                                   ,wdfld2_bcsu,wdfld2_ocsulf

! coated sulfate deposition in bins                           
real  , allocatable,dimension(:) :: con_s

real, allocatable,dimension(:) :: apm_wdeps_salt,apm_wdeps_dust &
                                   ,apm_wdeps_bcoc  &
                                   ,apm_wdep2s_salt &
                                   ,apm_wdep2s_dust,apm_wdep2s_bcoc

real, allocatable,dimension(:) :: apm_wdflds_salt,apm_wdflds_dust &
                                   ,apm_wdflds_bcoc &
                                   ,apm_wdfld2s_salt,apm_wdfld2s_dust &
                                   ,apm_wdfld2s_bcoc

!============================================



contains

subroutine allo_apm(nx,ny,nzz,nest,ctdway)

 implicit none
 include 'apm_parm.inc'
 
 integer :: ne,is,k
 integer :: ii,mem_apm,mem3d,mem2d

 integer :: nest

 integer :: nx(5),ny(5),nzz

 character :: ctdway*4
 logical :: lbin

 allocate(test_data(2))
 test_data(1:2) = (/1,2/)

 if(ctdway.eq.'bins') then
  lbin=.true.
 elseif(ctdway.eq.'bulk') then
  lbin=.false.
 else
  stop 'ctdway value set erro in emit/input.dat'
 endif
  
 sulf2dmem=0
 do ne=1,nest
 do is=1,NSO4
  sulf2dmem=sulf2dmem+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 allocate(apm_wdep_sulf(sulf2dmem))
 allocate(apm_wdep2_sulf(sulf2dmem))
 allocate(apm_wdfld_sulf(sulf2dmem))
 allocate(apm_wdfld2_sulf(sulf2dmem))


 salt2dmem=0
 do ne=1,nest
 do is=1,NSEA
  salt2dmem=salt2dmem+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 allocate(apm_wdep_salt(salt2dmem))
 allocate(apm_wdep2_salt(salt2dmem))
 allocate(apm_wdfld_salt(salt2dmem))
 allocate(apm_wdfld2_salt(salt2dmem))


 dust2dmem=0
 do ne=1,nest
 do is=1,NDSTB
  dust2dmem=dust2dmem+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 allocate(apm_wdep_dust(dust2dmem))
 allocate(apm_wdep2_dust(dust2dmem))
 allocate(apm_wdfld_dust(dust2dmem))
 allocate(apm_wdfld2_dust(dust2dmem))


 bc2dmem=0
 do ne=1,nest
 do is=1,2 ! two modes
  bc2dmem=bc2dmem+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo

 oc2dmem=bc2dmem


 bcoc2dmem=0
 do ne=1,nest
 do is=1,NDSTB
  bcoc2dmem=bcoc2dmem+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 allocate(apm_wdep_bcoc(bcoc2dmem))
 allocate(apm_wdep2_bcoc(bcoc2dmem))
 allocate(apm_wdfld_bcoc(bcoc2dmem))
 allocate(apm_wdfld2_bcoc(bcoc2dmem))


 ! apm sulfate particle
 mem_apm=0
 do ne=1,nest
 do is=1,NSO4
 do k=1,nzz
  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 enddo
 sulfmem=mem_apm
 allocate(apm_sulf(mem_apm),r_sulf(mem_apm))

 ii = 1
 allocate(ip_sulf(nzz,NSO4,nest))
 do ne = 1,nest
 do is = 1, NSO4
 do k = 1, nzz
   ip_sulf (k, is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo
 enddo

 ii=1
 allocate(ip2_sulf(NSO4,nest))
 do ne = 1,nest
 do is = 1, NSO4
   ip2_sulf (is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo


 ! apm sea salt
 mem_apm=0
 do ne=1,nest
 do is=1,NSEA
 do k=1,nzz
  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 enddo
 saltmem=mem_apm
 allocate(apm_salt(mem_apm),r_salt(mem_apm))
 allocate(salt_emit(mem_apm)); salt_emit=0.0d0
 !if(lbin) allocate(apm_msltsulf(mem_apm))


 ii = 1
 allocate(ip_salt(nzz,NSEA,nest))
 do ne = 1,nest
 do is = 1,NSEA
 do k = 1, nzz
   ip_salt (k, is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo
 enddo

 ii=1
 allocate(ip2_salt(NSEA,nest))
 do ne = 1,nest
 do is = 1, NSEA
   ip2_salt (is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo

 ! apm dust
 mem_apm=0
 do ne=1,nest
 do is=1,NDSTB
 do k=1,nzz
  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 enddo
 dustmem=mem_apm
 allocate(apm_dust(mem_apm),r_dust(mem_apm))
 allocate(dust_emit(mem_apm)); dust_emit=0.0d0
 !if(lbin) allocate(apm_mdstsulf(mem_apm))

 ii = 1
 allocate(ip_dust(nzz,NDSTB,nest))
 do ne = 1,nest
 do is = 1,NDSTB
 do k = 1, nzz
   ip_dust (k, is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo
 enddo

 ii=1
 allocate(ip2_dust(NDSTB,nest))
 do ne = 1,nest
 do is = 1, NDSTB
   ip2_dust (is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo


 ! apm bcoc
 mem_apm=0
 do ne=1,nest
 do is=1,NBCOCT
 do k=1,nzz
  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 enddo
 bcocmem=mem_apm
 allocate(apm_bcoc(mem_apm))
 !if(lbin) allocate(apm_mbcsulf(mem_apm),apm_mocsulf(mem_apm))

 mem_apm=0
 do ne=1,nest
 do is=1,2 ! two mode
 do k=1,nzz
  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 enddo
 allocate(ref_bcoc(mem_apm))

 ii=1
 allocate(ipmode_bcoc(nzz,2,nest))
 do ne = 1,nest
 do is = 1,2 ! two mode
 do k = 1, nzz
   ipmode_bcoc (k, is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo
 enddo




 ii = 1
 allocate(ip_bcoc(nzz,NBCOCT,nest))
 do ne = 1,nest
 do is = 1,NBCOCT
 do k = 1, nzz
   ip_bcoc (k, is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo
 enddo

 ii=1
 allocate(ip2_bcoc(NBCOCT,nest))
 do ne = 1,nest
 do is = 1, NBCOCT
   ip2_bcoc (is, ne) = ii
   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
 enddo
 enddo

 !
 mem3d=0
 do ne=1,nest
 do k=1,nzz
   mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
 enddo
 enddo
 allocate(apm_xq3d(mem3d))
 allocate(apm_cacid(mem3d),apm_cacid_old(mem3d),apm_pacid(mem3d))
 allocate(sulf_emit(mem3d),bcoc_emit(mem3d),bc_emit(mem3d),oc_emit(mem3d))
 sulf_emit=0.0d0; bcoc_emit=0.0d0; bcoc_emit=0.0d0; oc_emit=0.0d0

 allocate(ch_ratio(mem3d))


 !allocate(salt_emit(mem3d))

 !if(.not.lbin) then
  allocate(msltsulf(mem3d),mdstsulf(mem3d),mbcsulf(mem3d),mocsulf(mem3d))
  allocate(bulk_msp(mem3d))
  allocate(tot_sulf(mem3d))
  allocate(linit_box(mem3d)); linit_box=.true.
 !endif

 mem2d=0
 do ne=1,nest
   mem2d=mem2d+(nx(ne)+2)*(ny(ne)+2)
 enddo
 allocate(apmdryvel2d(mem2d))





 !> boundary memory and ip

 ! sulf bdy
 mem3d=0
 do ne=1,nest
 do is=1,NSO4
 do k=1,nzz
  mem3d=mem3d+(nx(ne)+2)
 enddo
 enddo
 enddo
 allocate(sulf_sbdy(mem3d),sulf_nbdy(mem3d))

 allocate(ipbsn_sulf(nzz,NSO4,nest))
 ii=1
 do ne=1,nest
 do is=1,NSO4
 do k=1,nzz
   ipbsn_sulf(k,is,ne)=ii
   ii=ii+(nx(ne)+2)
 enddo
 enddo
 enddo

 mem3d=0
 do ne=1,nest
 do is=1,NSO4
 do k=1,nzz
  mem3d=mem3d+(ny(ne)+2)
 enddo
 enddo
 enddo
 allocate(sulf_wbdy(mem3d),sulf_ebdy(mem3d))

 allocate(ipbwe_sulf(nzz,NSO4,nest))
 ii=1
 do ne=1,nest
 do is=1,NSO4
 do k=1,nzz
   ipbwe_sulf(k,is,ne)=ii
   ii=ii+(ny(ne)+2)
 enddo
 enddo
 enddo

 ! salt bdy
 mem3d=0
 do ne=1,nest
 do is=1,NSEA
 do k=1,nzz
  mem3d=mem3d+(nx(ne)+2)
 enddo
 enddo
 enddo
 allocate(salt_sbdy(mem3d),salt_nbdy(mem3d))

 allocate(ipbsn_salt(nzz,NSEA,nest))
 ii=1
 do ne=1,nest
 do is=1,NSEA
 do k=1,nzz
   ipbsn_salt(k,is,ne)=ii
   ii=ii+(nx(ne)+2)
 enddo
 enddo
 enddo

 mem3d=0
 do ne=1,nest
 do is=1,NSEA
 do k=1,nzz
  mem3d=mem3d+(ny(ne)+2)
 enddo
 enddo
 enddo
 allocate(salt_wbdy(mem3d),salt_ebdy(mem3d))

 allocate(ipbwe_salt(nzz,NSEA,nest))
 ii=1
 do ne=1,nest
 do is=1,NSEA
 do k=1,nzz
   ipbwe_salt(k,is,ne)=ii
   ii=ii+(ny(ne)+2)
 enddo
 enddo
 enddo

 !dust bdy
 mem3d=0
 do ne=1,nest
 do is=1,NDSTB
 do k=1,nzz
  mem3d=mem3d+(nx(ne)+2)
 enddo
 enddo
 enddo
 allocate(dust_sbdy(mem3d),dust_nbdy(mem3d))

 allocate(ipbsn_dust(nzz,NSEA,nest))
 ii=1
 do ne=1,nest
 do is=1,NDSTB
 do k=1,nzz
   ipbsn_dust(k,is,ne)=ii
   ii=ii+(nx(ne)+2)
 enddo
 enddo
 enddo

 mem3d=0
 do ne=1,nest
 do is=1,NDSTB
 do k=1,nzz
  mem3d=mem3d+(ny(ne)+2)
 enddo
 enddo
 enddo
 allocate(dust_wbdy(mem3d),dust_ebdy(mem3d))

 allocate(ipbwe_dust(nzz,NSEA,nest))
 ii=1
 do ne=1,nest
 do is=1,NDSTB
 do k=1,nzz
   ipbwe_dust(k,is,ne)=ii
   ii=ii+(ny(ne)+2)
 enddo
 enddo
 enddo

 ! bcoc
 mem3d=0
 do ne=1,nest
 do is=1,NBCOCT
 do k=1,nzz
  mem3d=mem3d+(nx(ne)+2)
 enddo
 enddo
 enddo
 allocate(bcoc_sbdy(mem3d),bcoc_nbdy(mem3d))

 allocate(ipbsn_bcoc(nzz,NBCOCT,nest))
 ii=1
 do ne=1,nest
 do is=1,NBCOCT
 do k=1,nzz
   ipbsn_bcoc(k,is,ne)=ii
   ii=ii+(nx(ne)+2)
 enddo
 enddo
 enddo

 mem3d=0
 do ne=1,nest
 do is=1,NBCOCT
 do k=1,nzz
  mem3d=mem3d+(ny(ne)+2)
 enddo
 enddo
 enddo
 allocate(bcoc_wbdy(mem3d),bcoc_ebdy(mem3d))

 allocate(ipbwe_bcoc(nzz,NBCOCT,nest))
 ii=1
 do ne=1,nest
 do is=1,NBCOCT
 do k=1,nzz
   ipbwe_bcoc(k,is,ne)=ii
   ii=ii+(ny(ne)+2)
 enddo
 enddo
 enddo


end subroutine allo_apm
!======================




!===================================
subroutine apm_basic_var(naqpms_dir)

 use APM_INIT_MOD, only : APM_INIT,APM_NTRACERS
 use APM_INIT_MOD, only : RDRY,RSALT,RDST 
 use APM_INIT_MOD, only : CEMITSULF2
 use APM_INIT_MOD, only : DFMSALT9
 use APM_INIT_MOD, only : CEMITBCOC2
 use APM_INIT_MOD, only : TOTNUMBC,TOTNUMOC,REFBCOC 
 use APM_INIT_MOD, only : YGF

 use APM_INIT_MOD, only : IFNUCL

! use APM_INIT_MOD, only : IFNUCL,IFAG
! use APM_INIT_MOD, only : XMACID,XMLVSOG,M1ACID,M1LVSOG
! use APM_INIT_MOD, only : VDRY
! use APM_INIT_MOD, only : DENSULF
! use APM_COAG_MOD, only : READCK6DTABLE
! use APM_NUCL_MOD, only : IONRATE0

 implicit none
 include 'apm_parm.inc'
 integer :: ibin
 integer :: N_APMTRAC

 character :: naqpms_dir*500
 character(len=255) :: DATA_DIR_1x1a

 !character(len=255),parameter :: DATA_DIR_1x1a = &
 !& '/extra01/shun/new_naqpms/naqpms_v0.07/src/apm_phys_box'


 !allocate( r_sulf(NSO4),r_salt(NSEA),r_dust(NDSTB) )

 allocate( rd_sulf(NSO4), rd_salt(NSEA), rd_dust(NDSTB) )
 allocate( sulf_unit_em2(NSO4,2), bcoc_unit_em2(NSO4,2) )
 allocate( salt_flux(NSEA) )
 allocate( grfac(99,4) )
 allocate( tn_bc(2), tn_oc(2), ref1d_bcoc(2) )

 IFNUCL=1

 call APM_NTRACERS( 0, N_APMTRAC ) ! in 'apm_init_mod.f'

 DATA_DIR_1x1a = trim(naqpms_dir)//'/src/apm_phys_box'

 call APM_INIT(DATA_DIR_1x1a)

 if(NDSTB.eq.4) then
    allocate(fdstbin(NDSTB))
    fdstbin(1:NDSTB) = (/ 0.1, 0.3, 0.3, 0.3 /)
 else
    stop 'NDSTB set erro'
 endif

 sulf_unit_em2=CEMITSULF2
 bcoc_unit_em2=CEMITBCOC2

 tn_oc=TOTNUMOC
 tn_bc=TOTNUMBC

 ref1d_bcoc=REFBCOC

 rd_sulf=RDRY
 rd_salt=RSALT
 rd_dust=RDST

 salt_flux=DFMSALT9

 grfac=YGF

 print*,'num of dust bins=',NDSTB
 print*,'dust_radius=',rd_dust*1.0e6

 print*,'DUST',RDST*1.0e6

 print*,'salt_flux',salt_flux

 return

 DO ibin=1,99
  print*,'YGF',ibin,grfac(ibin,1)
 ENDDO


end subroutine apm_basic_var




END MODULE apm_varlist
