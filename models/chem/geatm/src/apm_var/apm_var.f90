
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
!real,parameter :: fbc_qw=0.5,fbc_sw=1.0-fbc_qw
real,parameter :: fbc_qw=0.0,fbc_sw=1.0-fbc_qw
real,parameter :: foc_qw=0.5,foc_sw=1.0-foc_qw
real,parameter :: fs2bc=0.5,fs2oc=1.0-fs2bc
real,parameter :: apm_Pi = 3.14159265358979

real,parameter :: densulf=1.7 ! g/cm3
real,parameter :: densalt=2.2 ! g/cm3
real,parameter :: dendust=2.5 ! g/cm3
real,parameter :: denbcoc=1.8 ! g/cm3

! particle diameters in two static mode assumption
real,parameter :: diam_fine=1.0E-6*sqrt(0.1*2.5)     ! meters
real,parameter :: diam_coarse=1.0E-6*sqrt(2.5*10.0)  ! meters

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real :: apmdensity

!
integer :: apmfunit
integer :: runned_mn


logical :: ltest_box=.false.

logical,allocatable,dimension(:) :: linit_box

! aqueous chemistry
!logical,allocatable,dimension(:) :: lcloudy,lcld_app,lcld_dis
!integer,allocatable,dimension(:) :: cld_flag
!real   ,allocatable,dimension(:) :: clw_ph
!

!---------------------------
! flag variables in namelist 
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

real    :: frac_sp,frac_pp

logical :: lbox_phy &
          ,lbox_nucl,lbox_cond,lbox_cond_other,lbox_coag,lbox_coag_scav

logical :: ltest_dyn,ltest_phy

integer :: ncycle

logical :: lapm_chbdy,lapm_nogdt

real    :: zero00
real    :: tau_hb

real    :: upper,lower

logical :: lapm_bpho,lapm_wreach

logical :: lapm_dust00

logical :: lapm_wetsize

logical :: lapm_trhl

logical :: lapm_pso4

logical :: lapm_restart

logical :: luse_apmaeraq
logical :: luse_apmsulf

logical :: lbincb

!----------------------

namelist /apm_lvar/ &
 &  lapm_emit,lapm_hadv,lapm_vadv,lapm_hdif,lapm_vdif    &
 & ,lapm_ccld,lapm_ddep,lapm_gdep,lapm_acid,lapm_phys    &
 & ,ls2sulf,lapm_wdep,lapm_init,lapm_0emt,lapm_sulf_emit &
 & ,lapm_salt_emit,lapm_dust_emit,lapm_4bdy,lapm_tbdy    &
 & ,dep_method,lapm_wrfl,spinup_time,frac_sp,lcoated_dyn &
 & ,lapm_wr00,lbox_phy,ltest_dyn,ltest_phy,lapm_isor     &
 & ,ncycle,lapm_chbdy,lapm_nogdt,zero00,upper,lower      &
 & ,lapm_bpho,lapm_wreach,lapm_dust00,lapm_wetsize       &
 & ,lbox_nucl,lbox_cond,lbox_cond_other,lbox_coag,lbox_coag_scav &
 & ,tau_hb,lapm_trhl,lapm_pso4,lapm_restart &
 & ,luse_apmsulf,luse_apmaeraq &
 & ,lbincb


logical :: lfor_sulf,lfor_salt,lfor_dust,lfor_bcoc,lfor_h2so4
logical :: lem_sulf,lem_salt,lem_dust,lem_bcoc
namelist /apm_ldyn/ lfor_sulf,lfor_salt,lfor_dust,lfor_bcoc &
                   ,lem_sulf,lem_salt,lem_dust,lem_bcoc,lfor_h2so4

integer :: iflag_nucleation
logical :: lfrac_so2_emit,lfrac_nh3_emit,lppsize_ch,lppmass_ch
integer :: iflag_ppsize
character :: cflag_ppsize*10
real    :: frc_so2_emit,frc_nh3_emit,frc_size_ch,frc_mass_ch
logical :: lsivgrow,lsovgrow,llovgrow

namelist /apm_sensitivity/ iflag_nucleation &
  & ,lfrac_so2_emit,lfrac_nh3_emit,lppsize_ch,lppmass_ch &
  & ,frc_so2_emit,frc_nh3_emit,frc_size_ch,frc_mass_ch &
  & ,iflag_ppsize,cflag_ppsize,lsivgrow,lsovgrow,llovgrow


integer :: sulf2dmem,salt2dmem,dust2dmem,bc2dmem,oc2dmem,bcoc2dmem
integer,allocatable,dimension(:,:) :: ip2d_sulf,ip2d_salt,ip2d_dust &
                                   & ,ip2d_bc,ip2d_oc,ip2d_bcoc

integer :: sulfmem,saltmem,dustmem,bcmem,ocmem,bcocmem

real, allocatable,dimension(:) :: sulf_area,salt_area,dust_area,bcoc_area,bc_area,oc_area

! particle tracers
real, allocatable,dimension(:) :: apm_sulf,apm_salt,apm_dust,apm_bcoc ! ug/m3
real, allocatable,dimension(:) :: apm_binbc,apm_binoc
!

real, allocatable,dimension(:) :: sulf_sbdy,sulf_nbdy,sulf_wbdy,sulf_ebdy
real, allocatable,dimension(:) :: salt_sbdy,salt_nbdy,salt_wbdy,salt_ebdy
real, allocatable,dimension(:) :: dust_sbdy,dust_nbdy,dust_wbdy,dust_ebdy
real, allocatable,dimension(:) :: bcoc_sbdy,bcoc_nbdy,bcoc_wbdy,bcoc_ebdy


real, allocatable,dimension(:) :: sulf_emit,bcoc_emit,bc_emit,oc_emit &
                                   ,salt_emit,dust_emit
real, allocatable,dimension(:) :: bc_emit_bb,bc_emit_ff,oc_emit_bb,oc_emit_ff

real, allocatable,dimension(:) :: ppmfine_emit
real, allocatable,dimension(:) :: naq_nh4,naq_no3


real, allocatable,dimension(:) :: fdstbin


real,allocatable,dimension(:)    :: ch_ratio


! coated sulfate in bins
real, allocatable,dimension(:) ::  &
&  apm_mbcsulf,apm_mocsulf,apm_mdstsulf,apm_msltsulf,apm_mbcocsulf !&
!& ,apm_mbcsoa ,apm_mocsoa ,apm_mdstsoa ,apm_msltsoa  &
!& ,apm_mbcnit ,apm_mocnit ,apm_mdstnit ,apm_msltnit  &
!& ,apm_mbcnh4 ,apm_mocnh4 ,apm_mdstnh4 ,apm_msltnh4  &
!& ,apm_mbcmsa ,apm_mocmsa ,apm_mdstmsa ,apm_msltmsa  


! coated tracers
real, allocatable,dimension(:) ::  &
&  mbcsulf,mocsulf,mdstsulf,msltsulf !&
!& ,mbcsoa ,mocsoa ,mdstsoa ,msltsoa  &
!& ,mbcnit ,mocnit ,mdstnit ,msltnit  &
!& ,mbcnh4 ,mocnh4 ,mdstnh4 ,msltnh4  &
!& ,mbcmsa ,mocmsa ,mdstmsa ,msltmsa


integer,allocatable,dimension(:) :: ncoag1_3d,ncoag2_3d,ncoag3_3d,ncoag4_3d,ncoag5_3d


real,allocatable,dimension(:) :: deltso4_apm

real, allocatable,dimension(:) :: h2so4_gas,p_h2so4_so2_cbmz
real, allocatable,dimension(:) :: oh_radical


real, allocatable,dimension(:) :: bcgetsp_rate,ocgetsp_rate,bcagt,ocagt


! tracers for aqueous chemistry 

real, allocatable,dimension(:) :: sulf_inair,sulf_incld


! opt+ diagnostic variables
real, allocatable,dimension(:) :: ZBEXT3D,YBEXT3D,ZW3D
real, allocatable,dimension(:) :: vsblt2d
real, allocatable,dimension(:) :: AOD_390nm,AOD_500nm,AOD_530nm,AOD_550nm &
                                 ,AOD_700nm,AOD_1010nm
real, allocatable,dimension(:) :: AAOD_390nm,AAOD_500nm,AAOD_530nm,AAOD_550nm &
                                 ,AAOD_700nm,AAOD_1010nm
real, allocatable,dimension(:) :: TAOD
real, allocatable,dimension(:) :: ycs3d,surf3d

real, allocatable,dimension(:) :: wrfvsb

integer,allocatable,dimension(:) :: ifwrffog
integer,allocatable,dimension(:) :: ifhaze,iffog


! CCN fraction
real, allocatable,dimension(:) :: frac_ccn,pso4_so2,accu_pso4_so2
! CCN number
real, allocatable,dimension(:) :: number_ccn1,number_ccn2,number_ccn3
real, allocatable,dimension(:) :: cn10nm,cn3nm

real, allocatable,dimension(:) :: so4_ccn2,dust_ccn2,salt_ccn2,bc_ccn2,oc_ccn2

real, allocatable,dimension(:) :: rgfdry3d
real, allocatable,dimension(:) :: ztn3d
integer,allocatable,dimension(:,:,:) :: ip_type

!opt+
integer,allocatable,dimension(:,:,:,:) :: ipwl_type
integer :: iwl
integer,allocatable,dimension(:,:) :: ip_2dtype


real, allocatable,dimension(:) :: bulk_msp,tot_sulf

integer,allocatable,dimension(:,:,:) :: ip_sulf,ip_salt,ip_dust,ip_bcoc,ip_fccn
integer,allocatable,dimension(:,:,:) :: ipmode_bcoc

integer,allocatable,dimension(:,:,:) :: ipbsn_sulf,ipbsn_salt,ipbsn_dust,ipbsn_bcoc
integer,allocatable,dimension(:,:,:) :: ipbwe_sulf,ipbwe_salt,ipbwe_dust,ipbwe_bcoc
integer,allocatable,dimension(:,:)   :: ip2_sulf,ip2_salt,ip2_dust,ip2_bcoc

integer,allocatable,dimension(:,:)   :: ip2_binbc,ip2_binoc

!opt+
integer,allocatable,dimension(:,:,:) :: ip_ext

! bincb
integer,allocatable,dimension(:,:,:) :: ip_cbbin



integer :: mem_apm

integer,allocatable,dimension(:) :: test_data

real, allocatable,dimension(:) :: apm_xq3d,npf3d
real, allocatable,dimension(:) :: apm_cacid,apm_cacid_old,apm_pacid

real, allocatable,dimension(:) :: acid_gas_1,acid_gas_2

real,parameter :: apm_avogad=6.02217e+23

real,parameter :: ppbunit=1.0e+9

!real :: dt_cacid

real :: apm_cair_mlc

real :: apm_pr_atm,apm_te
!>------------------------

real, allocatable,dimension(:) :: rd_sulf,rd_salt,rd_dust
!real, allocatable,dimension(:) :: rw_sulf,rw_salt,rw_dust,refw_bc,refw_oc

real, allocatable,dimension(:) :: rd_binbc,rd_binoc
real, allocatable,dimension(:,:) :: bc_emit2bin_mfrc,oc_emit2bin_mfrc


real, allocatable,dimension(:) :: rgf_sulf,rgf_salt,rgf_dust,rgf_bc,rgf_oc

real, allocatable,dimension(:) :: spgf_3d

! sea salt flux ( kg/(m2 s) ) where u10 is 9m/s
real, allocatable,dimension(:) :: salt_flux 

real, allocatable,dimension(:) :: tn_bc,tn_oc,ref1d_bcoc

real, allocatable,dimension(:,:) :: sulf_unit_em2,bcoc_unit_em2
!TOTNUMBC,TOTNUMOC,REFBCOC,CEMITSULF2,CEMITBCOC2

real, allocatable,dimension(:,:) :: grfac
!>-----------------------------------------


real, allocatable,dimension(:) :: apmdryvel2d


!============================================
!> wet deposition variables

real,allocatable,dimension(:,:) :: apm_depfld,apm_depfld2
real,allocatable,dimension(:,:) :: apm_depflds,apm_depfld2s

real, allocatable,dimension(:) :: apm_wdep_sulf,apm_wdep_salt,apm_wdep_dust &
                                 ,apm_wdep_bcoc &
                                 ,apm_wdep_binbc,apm_wdep_binoc &
                                 ,apm_wdep2_sulf,apm_wdep2_salt,apm_wdep2_dust &
                                 ,apm_wdep2_bcoc &
                                 ,apm_wdep2_binbc,apm_wdep2_binoc


! coated sulfate deposition in bulk
real, allocatable,dimension(:) :: wdep_sltsulf,wdep_dstsulf,wdep_bcsulf &
                                   ,wdep_ocsulf &
                                   ,wdep2_sltsulf,wdep2_dstsulf,wdep2_bcsulf &
                                   ,wdep2_ocsulf


real, allocatable,dimension(:) :: wdep_h2so4,wdep2_h2so4


! coated sulfate deposition in bins                           

real, allocatable,dimension(:) :: apm_wdeps_salt,apm_wdeps_dust &
                                   ,apm_wdeps_bcoc  &
                                   ,apm_wdep2s_salt &
                                   ,apm_wdep2s_dust,apm_wdep2s_bcoc


real, allocatable,dimension(:) :: apm_wdeps_binbc,apm_wdep2s_binbc &
                                 ,apm_wdeps_binoc,apm_wdep2s_binoc


!============================================


!============================================
! 20140604@albany
real, allocatable,dimension(:) :: naqpms_sulfate,apm_sulfate,apm_ppsulfate,apm_spsulfate
real, allocatable,dimension(:) :: ddep_naqpms_sulf,ddep_apm_sulf,ddep_apm_ppsulf,ddep_apm_spsulf


!============================================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

subroutine allo_apm(nx,ny,nzz,nest,ctdway,sx,ex,sy,ey,mem_per_block)

 implicit none
 include 'apm_parm.inc'
 
!==================================
!opt+
 INTEGER, PARAMETER :: MWLS=16

 integer :: ne,is,k
 integer :: ii,mem_apm,mem3d,mem2d

 integer :: nest

 integer :: nx(5),ny(5),nzz
 integer :: sx(5), ex(5), sy(5), ey(5)

 character :: ctdway*4
 logical :: lbin

 integer :: mem_per_block(5)


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
!  sulf2dmem=sulf2dmem+(nx(ne)+2)*(ny(ne)+2)
   sulf2dmem=sulf2dmem+mem_per_block(ne)
 enddo
 enddo
 allocate(apm_wdep_sulf(sulf2dmem))
 allocate(apm_wdep2_sulf(sulf2dmem))
 !allocate(apm_wdfld_sulf(sulf2dmem))
 !allocate(apm_wdfld2_sulf(sulf2dmem))


 salt2dmem=0
 do ne=1,nest
 do is=1,NSEA
!  salt2dmem=salt2dmem+(nx(ne)+2)*(ny(ne)+2)
   salt2dmem=salt2dmem+mem_per_block(ne)
 enddo
 enddo
 allocate(apm_wdep_salt(salt2dmem))
 allocate(apm_wdep2_salt(salt2dmem))
 !allocate(apm_wdfld_salt(salt2dmem))
 !allocate(apm_wdfld2_salt(salt2dmem))
 allocate(apm_wdeps_salt(salt2dmem))
 allocate(apm_wdep2s_salt(salt2dmem))


 dust2dmem=0
 do ne=1,nest
 do is=1,NDSTB
!  dust2dmem=dust2dmem+(nx(ne)+2)*(ny(ne)+2)
  dust2dmem=dust2dmem+mem_per_block(ne)
 enddo
 enddo
 allocate(apm_wdep_dust(dust2dmem))
 allocate(apm_wdep2_dust(dust2dmem))
 !allocate(apm_wdfld_dust(dust2dmem))
 !allocate(apm_wdfld2_dust(dust2dmem))
 allocate(apm_wdeps_dust(dust2dmem))
 allocate(apm_wdep2s_dust(dust2dmem))


 bc2dmem=0
 do ne=1,nest
 do is=1,2 ! two modes
!  bc2dmem=bc2dmem+(nx(ne)+2)*(ny(ne)+2)
   bc2dmem=bc2dmem+mem_per_block(ne)
 enddo
 enddo

 oc2dmem=bc2dmem


 bcoc2dmem=0
 do ne=1,nest
 do is=1,NBCOCT
!  bcoc2dmem=bcoc2dmem+(nx(ne)+2)*(ny(ne)+2)
  bcoc2dmem=bcoc2dmem+mem_per_block(ne)
 enddo
 enddo
 allocate(apm_wdep_bcoc(bcoc2dmem))
 allocate(apm_wdep2_bcoc(bcoc2dmem))
 !allocate(apm_wdfld_bcoc(bcoc2dmem))
 !allocate(apm_wdfld2_bcoc(bcoc2dmem))
 allocate(apm_wdeps_bcoc(bcoc2dmem))
 allocate(apm_wdep2s_bcoc(bcoc2dmem))


 ! apm sulfate particle
 mem_apm=0
 do ne=1,nest
 do is=1,NSO4
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
   mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 sulfmem=mem_apm
 allocate( apm_sulf(mem_apm) )
 !allocate( rw_sulf(mem_apm) )

 ii = 1
 allocate(ip_sulf(nzz,NSO4,nest))
 do ne = 1,nest
 do is = 1, NSO4
 do k = 1, nzz
   ip_sulf (k, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo


 ii=1
 allocate(ip2_sulf(NSO4,nest))
 do ne = 1,nest
 do is = 1, NSO4
   ip2_sulf (is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo


 ! apm sea salt
 mem_apm=0
 do ne=1,nest
 do is=1,NSEA
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 saltmem=mem_apm
 allocate( apm_salt(mem_apm) )
 !allocate( rw_salt(mem_apm) )
 allocate(salt_emit(mem_apm)); salt_emit=0.0d0
 !if(lbin) allocate(apm_msltsulf(mem_apm))


 ii = 1
 allocate(ip_salt(nzz,NSEA,nest))
 do ne = 1,nest
 do is = 1,NSEA
 do k = 1, nzz
   ip_salt (k, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo

 ii=1
 allocate(ip2_salt(NSEA,nest))
 do ne = 1,nest
 do is = 1, NSEA
   ip2_salt (is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo

 ! apm dust
 mem_apm=0
 do ne=1,nest
 do is=1,NDSTB
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
    mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 dustmem=mem_apm
 allocate( apm_dust(mem_apm) )
 !allocate( rw_dust(mem_apm) )
 allocate(dust_emit(mem_apm)); dust_emit=0.0d0
 !if(lbin) allocate(apm_mdstsulf(mem_apm))

 ii = 1
 allocate(ip_dust(nzz,NDSTB,nest))
 do ne = 1,nest
 do is = 1,NDSTB
 do k = 1, nzz
   ip_dust (k, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo

 ii=1
 allocate(ip2_dust(NDSTB,nest))
 do ne = 1,nest
 do is = 1, NDSTB
   ip2_dust (is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo


 ! apm bcoc
 mem_apm=0
 do ne=1,nest
 do is=1,NBCOCT
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
   mem_apm=mem_apm+mem_per_block(ne)
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
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 !allocate(refw_bc(mem_apm),refw_oc(mem_apm))

 ii=1
 allocate(ipmode_bcoc(nzz,2,nest))
 do ne = 1,nest
 do is = 1,2 ! two mode
 do k = 1, nzz
   ipmode_bcoc (k, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo




 ii = 1
 allocate(ip_bcoc(nzz,NBCOCT,nest))
 do ne = 1,nest
 do is = 1,NBCOCT
 do k = 1, nzz
   ip_bcoc (k, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo

 ii=1
 allocate(ip2_bcoc(NBCOCT,nest))
 do ne = 1,nest
 do is = 1, NBCOCT
   ip2_bcoc (is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo

! bincb
 mem_apm=0
 do ne=1,nest
 do is=1,nbincb
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 allocate(apm_binbc(mem_apm))
 allocate(apm_binoc(mem_apm))

 allocate(ip_cbbin(nzz,nbincb,nest))
 ii=1
 do ne = 1,nest
 do is = 1,nbincb
 do k = 1, nzz
   ip_cbbin (k, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo


 mem_apm=0
 do ne = 1,nest
 do is = 1,nbincb
!   mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)  
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
  allocate(apm_wdep_binbc(mem_apm))
  allocate(apm_wdep_binoc(mem_apm))
  allocate(apm_wdep2_binbc(mem_apm))
  allocate(apm_wdep2_binoc(mem_apm))

  allocate(apm_wdeps_binbc(mem_apm))
  allocate(apm_wdep2s_binbc(mem_apm))
  allocate(apm_wdeps_binoc(mem_apm))
  allocate(apm_wdep2s_binoc(mem_apm))


 allocate(ip2_binbc(nbincb,nest))
 allocate(ip2_binoc(nbincb,nest))
 ii=1
 do ne = 1,nest
 do is = 1,nbincb
   ip2_binbc(is,ne)=ii
   ip2_binoc(is,ne)=ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo



 mem_apm=0
 do ne=1,nest
 do is=1,NSO4+4
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 allocate(frac_ccn(mem_apm))


 ii=1
 allocate( ip_fccn(nzz,NSO4+4,nest) )
 do ne = 1,nest
 do is = 1,NSO4+4
 do k = 1,nzz
   ip_fccn (k,is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo


 ii = 1
 allocate(ip_type(nzz,5,nest))
 do ne = 1,nest
 do is = 1, 5 ! five types
 do k = 1, nzz
   ip_type (k, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo

 mem_apm=0
 do ne=1,nest
 do is=1,5  ! five types
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 allocate(number_ccn1(mem_apm))
 allocate(number_ccn2(mem_apm))
 allocate(number_ccn3(mem_apm))
 allocate(ztn3d(mem_apm))

 allocate(ycs3d(mem_apm))
 allocate(surf3d(mem_apm))

 allocate(rgfdry3d(mem_apm))


 allocate(ip_2dtype(5,nest))
 ii = 1
 do ne = 1,nest
 do is = 1, 5 ! five types
   ip_2dtype (is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo

 mem_apm=0
 do ne=1,nest
 do is=1,5  ! five types
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 allocate(TAOD(mem_apm))

 !
 mem3d=0
 do ne=1,nest
 do k=1,nzz
!   mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
  mem3d=mem3d+mem_per_block(ne)
 enddo
 enddo

 allocate(deltso4_apm(mem3d))


 allocate(apm_xq3d(mem3d))
 allocate(npf3d(mem3d))
 allocate(apm_cacid(mem3d),apm_cacid_old(mem3d),apm_pacid(mem3d))
 allocate(sulf_emit(mem3d),bcoc_emit(mem3d),bc_emit(mem3d),oc_emit(mem3d))
 allocate(bc_emit_bb(mem3d),oc_emit_bb(mem3d))
 allocate(bc_emit_ff(mem3d),oc_emit_ff(mem3d))
 sulf_emit=0.0d0; bcoc_emit=0.0d0; bcoc_emit=0.0d0; oc_emit=0.0d0
 bc_emit_bb=0.0; oc_emit_bb=0.0; bc_emit_ff=0.0; oc_emit_ff=0.0


 allocate(ppmfine_emit(mem3d)); ppmfine_emit=0.0
 
 allocate(naq_nh4(mem3d))
 allocate(naq_no3(mem3d))


 allocate(h2so4_gas(mem3d))
 allocate(p_h2so4_so2_cbmz(mem3d))
 allocate(oh_radical(mem3d))

 allocate( rgf_sulf(mem3d),rgf_salt(mem3d),rgf_dust(mem3d) &
          ,rgf_bc(mem3d),rgf_oc(mem3d) ) 

 allocate(spgf_3d(mem3d))

! allocate( number_ccn1(mem3d) )
! allocate( number_ccn2(mem3d) )
! allocate( number_ccn3(mem3d) )
 allocate( cn10nm(mem3d),cn3nm(mem3d) )

 allocate( so4_ccn2(mem3d) )
 allocate( dust_ccn2(mem3d) )
 allocate( salt_ccn2(mem3d) )
 allocate( bc_ccn2(mem3d) )
 allocate( oc_ccn2(mem3d) )


 allocate(ch_ratio(mem3d))

 allocate(pso4_so2(mem3d))
 allocate(accu_pso4_so2(mem3d))

 allocate(sulf_inair(mem3d))
 allocate(sulf_incld(mem3d))


! 20160314@bj
 allocate(bcgetsp_rate(mem3d))
 allocate(ocgetsp_rate(mem3d))

 allocate(bcagt(mem3d))
 allocate(ocagt(mem3d))
 bcagt=36.0*3600.0 ! (h) 1.5 day
 ocagt=36.0*3600.0 ! (h) 1.5 day



!20140604@albany
!naqpms_sulfate,apm_sulfate,apm_ppsulfate,apm_spsulfate
 allocate(naqpms_sulfate(mem3d))
 allocate(apm_sulfate(mem3d))
 allocate(apm_ppsulfate(mem3d))
 allocate(apm_spsulfate(mem3d))

!opt+


 ! only for test
 allocate( acid_gas_1(mem3d), acid_gas_2(mem3d) )



 !allocate(salt_emit(mem3d))

 !if(.not.lbin) then
  allocate(msltsulf(mem3d),mdstsulf(mem3d),mbcsulf(mem3d),mocsulf(mem3d))
  allocate(bulk_msp(mem3d))
  allocate(tot_sulf(mem3d))
  allocate(linit_box(mem3d)); linit_box=.true.
 !endif

  allocate( ncoag1_3d(mem3d) &
           ,ncoag2_3d(mem3d) &
           ,ncoag3_3d(mem3d) &
           ,ncoag4_3d(mem3d) &
           ,ncoag5_3d(mem3d) )
  ncoag1_3d=0
  ncoag2_3d=0
  ncoag3_3d=0
  ncoag4_3d=0
  ncoag5_3d=0



 ! aqueous chemistry flags vars
  !allocate(lcloudy(mem3d))
  !allocate(lcld_app(mem3d))
  !allocate(lcld_dis(mem3d))
  !allocate(cld_flag(mem3d))
  !allocate(clw_ph(mem3d))

 mem2d=0
 do ne=1,nest
!   mem2d=mem2d+(nx(ne)+2)*(ny(ne)+2)
   mem2d=mem2d+mem_per_block(ne)
 enddo
 allocate(apmdryvel2d(mem2d))

 allocate(wdep_sltsulf(mem2d))
 allocate(wdep2_sltsulf(mem2d))

 allocate(wdep_dstsulf(mem2d))
 allocate(wdep2_dstsulf(mem2d))

 allocate(wdep_bcsulf(mem2d))
 allocate(wdep2_bcsulf(mem2d))

 allocate(wdep_ocsulf(mem2d))
 allocate(wdep2_ocsulf(mem2d))

 allocate(wdep_h2so4(mem2d))
 allocate(wdep2_h2so4(mem2d))


!20140604@albany
! ddep_naqpms_sulf,ddep_apm_sulf,ddep_apm_ppsulf,ddep_apm_ppsulf
 allocate(ddep_naqpms_sulf(mem2d))
 allocate(ddep_apm_sulf(mem2d))
 allocate(ddep_apm_ppsulf(mem2d))
 allocate(ddep_apm_spsulf(mem2d))
 

!opt+
 allocate( vsblt2d(mem2d) )
 allocate( ifhaze(mem2d) )
 allocate( iffog(mem2d) )

 allocate( AOD_390nm(mem2d)  )
 allocate( AOD_500nm(mem2d)  )
 allocate( AOD_530nm(mem2d)  )
 allocate( AOD_550nm(mem2d)  )
 allocate( AOD_700nm(mem2d)  )
 allocate( AOD_1010nm(mem2d) )

 allocate( AAOD_390nm(mem2d)  )
 allocate( AAOD_500nm(mem2d)  )
 allocate( AAOD_530nm(mem2d)  )
 allocate( AAOD_550nm(mem2d)  )
 allocate( AAOD_700nm(mem2d)  )
 allocate( AAOD_1010nm(mem2d) )

 allocate( wrfvsb(mem2d) )
 allocate( ifwrffog(mem2d))


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


!opt+
 allocate(ip_ext(nzz,MWLS,nest)) ! 550nm
 ii=1
 do ne = 1,nest
 do is = 1,MWLS
 do k = 1, nzz
   ip_ext (k, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo

 mem_apm=0
 do ne=1,nest
 do is=1,MWLS  ! 550nm
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 allocate( ZBEXT3D(mem_apm) )
 allocate( ZW3D(mem_apm) )

!opt+



 ii = 1
 allocate(ipwl_type(nzz,MWLS,5,nest))
 do ne = 1,nest
 do is = 1, 5 ! five types
 do iwl=1,MWLS
 do k = 1, nzz
   ipwl_type (k, iwl, is, ne) = ii
!   ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
   ii = ii + mem_per_block(ne)
 enddo
 enddo
 enddo
 enddo

 mem_apm=0 
 do ne=1,nest
 do is=1,5
 do iwl=1,MWLS
 do k=1,nzz
!  mem_apm=mem_apm+(nx(ne)+2)*(ny(ne)+2)
  mem_apm=mem_apm+mem_per_block(ne)
 enddo
 enddo
 enddo
 enddo
 allocate( YBEXT3D(mem_apm) )






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

 use APM_INIT_MOD, only : radcbm,cbemt2bin_mfrc

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

 allocate(rd_binbc(nbincb))
 allocate(rd_binoc(nbincb))
 allocate( bc_emit2bin_mfrc(nbincb,2) )
 allocate( oc_emit2bin_mfrc(nbincb,2) )

 !IFNUCL=1 
 IFNUCL=iflag_nucleation

! print*,'IFNUCL=',IFNUCL

!stop

 call APM_NTRACERS( 0, N_APMTRAC ) ! in 'apm_init_mod.f'

 DATA_DIR_1x1a = trim(naqpms_dir)//'/src/apm_phys_box'

! print*,'lppsize_ch=',lppsize_ch
! print*,

 call APM_INIT( DATA_DIR_1x1a &
               ,lppsize_ch,iflag_ppsize,cflag_ppsize )

! print*,'apm_var : DATA_DIR_1x1a=',DATA_DIR_1x1a

!stop 'test apm_init'


 if(NDSTB.eq.4) then
    allocate(fdstbin(NDSTB))
    fdstbin(1:NDSTB) = (/ 0.1, 0.3, 0.3, 0.3 /)
 else
    stop 'NDSTB set erro'
 endif

 sulf_unit_em2=CEMITSULF2
 bcoc_unit_em2=CEMITBCOC2

 bc_emit2bin_mfrc=cbemt2bin_mfrc
 oc_emit2bin_mfrc=cbemt2bin_mfrc

 tn_oc=TOTNUMOC
 tn_bc=TOTNUMBC

 ref1d_bcoc=REFBCOC

 rd_sulf=RDRY
 rd_salt=RSALT
 rd_dust=RDST

 rd_binbc=radcbm
 rd_binoc=radcbm

 salt_flux=DFMSALT9

 grfac=YGF


return
print*,sum(cbemt2bin_mfrc(:,1)),sum(cbemt2bin_mfrc(:,2))
stop
return

 print*,'num of dust bins=',NDSTB
 print*,'dust_radius=',rd_dust*1.0e6

 print*,'DUST',RDST*1.0e6

 print*,'salt_flux',salt_flux

return

 open(100,file=trim(naqpms_dir)//'/apm_bins/sulfate.bins')
 do ibin=1,NSO4
  write(100,*) ibin,rd_sulf(ibin)*1.0e6
 enddo
 close(100)

 open(100,file=trim(naqpms_dir)//'/apm_bins/seasalt.bins')
 do ibin=1,NSEA
  write(100,*) ibin,rd_salt(ibin)*1.0e6
 enddo
 close(100)

 open(100,file=trim(naqpms_dir)//'/apm_bins/dust.bins')
 do ibin=1,NDSTB
  write(100,*) ibin,rd_dust(ibin)*1.0e6
 enddo
 close(100)



 return

 DO ibin=1,99
  print*,'YGF',ibin,grfac(ibin,1)
 ENDDO


end subroutine apm_basic_var




END MODULE apm_varlist
