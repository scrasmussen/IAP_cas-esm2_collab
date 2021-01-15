
subroutine apm_boxphy_driver &
 & ( myid &
 &  ,lapm &
 &  ,dt_naqpms &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d &
 &  ,ip4mem,mem4d &
 &  ,igas,gas,GC_MOLWT &
 &  ,longicrs,latitcrs,land_use &
 &  ,PSFC,Plev,t,rh1)

use apm_varlist
use aqchem_varlist
use APM_PHYS_MOD, only : APM_PHYS
use APM_COAG_MOD, only : READCK6DTABLE
use APM_NUCL_MOD, only : IONRATE0
use APM_INIT_MOD, only : VDRY,IFNUCL
use APM_INIT_MOD, only : NTYP

implicit none
include 'apm_parm.inc'

real,parameter :: PA2ATM = 1.0/101325.0

integer :: myid
real    :: dt_naqpms
logical :: lapm
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,ig,idx,is,N,imode

integer               :: mem2d
real,dimension(mem2d) :: longicrs,latitcrs,land_use,PSFC

integer               :: mem3d
real,dimension(mem3d) :: Plev,t,rh1

integer               :: mem4d
real,dimension(mem4d) :: gas ! ppb(gas),ug/m3(aerosol)

real, dimension(102)  :: GC_MOLWT

integer :: ixy,i02,i03,i04,iapm,i03_00,i04_acid,iccn

integer :: ip2mem(nest),ip3mem(nzz,nest),ip4mem(nzz,igas,nest)

real,parameter :: r_const = 8.31
real,parameter :: e_9=1.0e-9,e_3=1.0e-3

integer           :: isoa,init,inh4,imsa,isulf

integer,parameter :: nsoa = 6 ! SOA1,...,SOA6
integer,parameter :: soa_o_gas(1:nsoa) = (/ 1, 1, 1, 1,  1,  1/)
integer,parameter :: soa_index(1:nsoa) = (/96,97,98,99,100,101/) 


integer,parameter :: nnit = 3 ! NO3AQ,NANO3,NH4NO3
integer,parameter :: nit_o_gas(1:nnit) = (/ 1, 1, 1/)
integer,parameter :: nit_index(1:nnit) = (/85,88,90/)
!integer,parameter :: nnit = 2 ! NO3AQ,NANO3,NH4NO3
!integer,parameter :: nit_o_gas(1:nnit) = (/  1, 1/)
!integer,parameter :: nit_index(1:nnit) = (/ 88,90/)


integer,parameter :: nnh4 = 6 ! NH4AQ,NH42SO4,NH4NO3,NH4CL,NH4HSO4S,NH44HSO42
integer,parameter :: nh4_o_gas(1:nnh4) = (/ 1, 2, 1, 1, 1, 1/)
integer,parameter :: nh4_index(1:nnh4) = (/81,89,90,91,93,95/)

real              :: mwght

integer,parameter :: nmsa = 1
integer,parameter :: msa_o_gas(1:nmsa) = (/1/)
integer,parameter :: msa_index(1:nmsa) = (/58/)


real              :: msulf,msulf_pp,msulf_sp
integer,parameter :: nsulf = 8
integer,parameter :: sulf_o_gas(1:nsulf) = (/ 1, 1, 1, 1, 1, 1, 1, 1/)
integer,parameter :: sulf_index(1:nsulf) = (/83,84,87,89,92,93,94,95/)

! apm acid
real,parameter :: gfgas1=1.0,gfhyg1=1.0,rgf1=1.0
real,parameter :: gfgas2=50.0,gfhyg2=10.0,rgf2=100.0
real,parameter :: gasdef=gfgas1,hygdef=gfhyg1,rgfdef=rgf1
real,parameter :: mcon1=0.0,mcon2=1000*1.0e-9 ! 0-1000 ug/m3

integer :: iostate

real :: tmp

! shun ++ 2013-02-11
real*8  :: SOAT_to_APM,MSPS
! ++++++++++++++++++

! aqueous 1d
real :: pso4_so2_1d


! apm 1d vars
real*8  :: DTAPM
real*8  :: XQ
real*8  :: DT,XLON,XLAT
real*8  :: PRESS,TK,RHIN
integer :: ISURF
real*8  :: YPSURF,YPR
real*8  :: CACID,PACID
real*8  :: SOAT,MNIT,MNH4,MMSA ! kg/m3
real*8  :: XN1D(NSO4+NSEA),XMDST(NDSTB),MBCOC8(8)
real*8  :: MBCS,MOCS,MDSTS,MSALTS
real*8  :: MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV
real*8  :: CLVSOG, PLVSOG1
real*8  :: XM1D(NSO4+NSEA)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! sulfate mass in each bin before aqueous chemistry
real*8  :: XM1D_old(NSO4+NSEA) ! shun
! particle wet size 
real*8  :: rw1d_sulf(NSO4),rw1d_salt(NSEA),rw1d_dust(NDSTB) &
        & ,refw1d_bc(2),refw1d_oc(2)
! radius growth factor
real*8  :: rgf_sulf_1d,rgf_salt_1d,rgf_dust_1d,rgf_bc_1d,rgf_oc_1d
! CCN number
real*8  :: number_ccn1_1d(5),number_ccn2_1d(5),number_ccn3_1d(5)
real*8  :: ztn1d(5)
real*8  :: npf_rate
real*8  :: so4_ccn2_1d,dust_ccn2_1d,salt_ccn2_1d,bc_ccn2_1d,oc_ccn2_1d
integer :: itype,i00_type
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


REAL*8  :: FCLOUD1(NSO4+4)
INTEGER :: NCOAG1,NCOAG2 ! shun ??
REAL*8  :: TEMPOUT(IIPAR,JJPAR,LLPAR,NTEMPOUT1),TEMPOUT1(NTEMPOUT1)

REAL*8  :: GFTOT1,GFTOT2,DENWET1,DENWET2
! bin index for cloud act diameters corresponding to RDRY
INTEGER :: IACT10, IACT20, IACT30
LOGICAL,SAVE :: FIRST  = .TRUE.
LOGICAL,SAVE :: FIRST1 = .TRUE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!================================
! opt+
INTEGER, PARAMETER :: MWLS=16
INTEGER, PARAMETER :: MWLL=9
REAL*8    :: ZBEXT(MWLS),ZW(MWLS),ZG(MWLS)
REAL*8    :: ZBABS(MWLL)
REAL*8    :: YBEXT(NTYP,MWLS),YW(NTYP,MWLS),YG(NTYP,MWLS)
integer   :: iwls

real*8  :: YCS_1D(NTYP)  ! s-1
real*8  :: SURF_1D(NTYP) ! cm2/cm3


integer :: lu_1d
integer,parameter :: luidx(1:24) = (/ 1,1,1,1,1,  1,1,1,1,1 &
                                   & ,1,1,1,1,1,  0,1,1,1,1 &
                                   & ,1,1,1,0 /)
logical,save :: init_box = .true.
!real*8 :: ch_ratio

!>-------------     end of definitions    ------------------<!
!=============================================================



! may exist problem : 2012-11-27_23:00


 if(FIRST1)then
   ! Read coagulation kernel look-up tables
   CALL  READCK6DTABLE ! in 'apm_coag_mod.f'
   TEMPOUT = 0d0
   !NPOUTSTEPS = 0
   FIRST1=.FALSE.
 endif



!IF(lapm) THEN ! apm flag

 i02 = ip2mem(ne)

 loop_j : DO j = sy(ne),ey(ne)
 loop_i : DO i = sx(ne),ex(ne)

   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

   lu_1d=int(land_use(i02+ixy))

   ISURF=luidx(lu_1d) 

   loop_k : DO k=1,nzz-1

     i03 = ip3mem(k,ne)
     i04 = ip4mem(k,ig,ne)

     ! met2d to met1d
     XLON = longicrs(i02+ixy)
     XLAT = latitcrs(i02+ixy)
     !
 
     ! met3d to met1d
     PRESS = Plev(i03+ixy)*100.0  ! hPa -> Pa
     TK    = t(i03+ixy)
     RHIN  = rh1(i03+ixy)         ! RHIN : % ; rh1 : ?

     YPSURF = PSFC(i02+ixy)/100.0 ! Pa -> mb(hPa)
     YPR    = Plev(i03+ixy)       ! mb
     !

     ! species3d to species1d     
     !CACID = apm_cacid(i03+ixy) ! #/cm3

     ! calculate sulferic acid vapor production rate and concentration

     i04_acid = ip4mem(k,1,ne) ! ig=1

     apm_pr_atm = PA2ATM*Plev(i03+ixy)*100.
     apm_te = t(i03+ixy)
     apm_cair_mlc = apm_avogad*apm_pr_atm/(82.056*apm_te)
!     apm_cacid(i03+ixy) = gas(i04_acid+ixy)* &
!                        & apm_cair_mlc/ppbunit

!     PACID = apm_pacid(i03+ixy) ! #/(cm3 s)
!     CACID = apm_cacid(i03+ixy) ! #/cm3

     PACID = p_h2so4_so2_cbmz(i03+ixy)*apm_cair_mlc/ppbunit
     CACID = h2so4_gas(i03+ixy)*apm_cair_mlc/ppbunit


     SOAT   = 0.0d0
     do isoa=1,nsoa
      idx   = soa_index(isoa)
      mwght = GC_MOLWT(idx) ! SOA
      ! kg/m3
      SOAT  = SOAT + gas(ip4mem(k,idx,ne)+ixy)*soa_o_gas(isoa)*mwght* &
                     e_9/GC_MOLWT(idx)
!      SOAT  = SOAT + gas(ip4mem(k,idx,ne)+ixy)*soa_o_gas(isoa)*mwght* &
!            &        e_9*e_3*(Plev(i03+ixy)*100)/(8.31*t(i03+ixy))
     enddo


     MNIT   = 0.0d0
     do init=1,nnit
      idx   = nit_index(init)
      mwght = 62.0 ! NO3
      ! kg/m3
      MNIT  = MNIT + gas(ip4mem(k,idx,ne)+ixy)*nit_o_gas(init)*mwght* &
                     e_9/GC_MOLWT(idx)
!      MNIT  = MNIT + gas(ip4mem(k,idx,ne)+ixy)*nit_o_gas(init)*mwght* &
!            &        e_9*e_3*(Plev(i03+ixy)*100)/(8.31*t(i03+ixy))
     enddo

     MNH4   = 0.0d0
     do inh4=1,nnh4
      idx   = nh4_index(inh4)
      mwght = 18.0 ! NH4
      ! kg/m3
      MNH4  = MNH4 + gas(ip4mem(k,idx,ne)+ixy)*nh4_o_gas(inh4)*mwght* &
                     e_9/GC_MOLWT(idx)
!      MNH4  = MNH4 + gas(ip4mem(k,idx,ne)+ixy)*nh4_o_gas(inh4)*mwght* &
!            &        e_9*e_3*(Plev(i03+ixy)*100)/(8.31*t(i03+ixy))
     enddo

     MMSA   = 0.0d0
     do imsa=1,nmsa
      idx   = msa_index(imsa)
      mwght = 98.0 
      ! ppb -> kg/m3
      MMSA  = MMSA + gas(ip4mem(k,idx,ne)+ixy)*msa_o_gas(inh4)*mwght* &
            &        e_9*e_3*(Plev(i03+ixy)*100)/(8.31*t(i03+ixy))
     enddo

     msulf=0.0d0
     do isulf=1,nsulf
      idx   = sulf_index(isulf)
      mwght = 96.0 ! NH4
      ! kg/m3
      msulf = msulf + gas(ip4mem(k,idx,ne)+ixy)*sulf_o_gas(isulf)*mwght* &
                      e_9/GC_MOLWT(idx)
     enddo

     !msulf=10*ug2kg ! test

     idx=14 ! OH radical
     oh_radical(i03+ixy)=gas(ip4mem(k,idx,ne)+ixy)


     if(.not.lsivgrow) then
       MMSA=1.d-30
       MNH4=1.d-30
       MNIT=1.d-30
     endif

     if(.not.lsovgrow) then
       SOAT=1.d-30
     endif

!

     if(linit_box(i03+ixy)) then
       if(ne.eq.1.and.k.eq.1.and.j.eq.1.and.i.eq.1) print*,'shun : init box'
       ! distribute sulfate by ratios assumed
       frac_pp=1.0-frac_sp
       msulf_sp=msulf*frac_sp ! kg/m3
       msulf_pp=msulf*frac_pp ! kg/m3

       !tot_sulf(i03+ixy)=msulf*kg2ug ! kg/m3->ug/m3

       bulk_msp(i03+ixy)=msulf_sp*kg2ug ! kg/m3->ug/m3
       ! get apm_sulf
       call get_sp_distribution(myid,i,j,k,ne,ixy,msulf_sp)
       ! get coated sulfate
       call get_coated_sulfate(myid,i,j,k,ne,ixy &
                              ,msulf_pp,MSALTS,MDSTS,MBCS,MOCS) ! kg/m3

       msltsulf(i03+ixy)=MSALTS*kg2ug
       mdstsulf(i03+ixy)=MDSTS*kg2ug
       mbcsulf(i03+ixy)=MBCS*kg2ug
       mocsulf(i03+ixy)=MOCS*kg2ug

       ! init CCN fraction : assuming value is zero
       do is=1,NSO4+4
         iccn=ip_fccn(k,is,ne)
         frac_ccn(iccn+ixy)=0.0
       enddo

       linit_box(i03+ixy)=.false.
     endif

     !================================================================
     ! update sulfate mass ( SP particle and coating sulfate ) tracers
     ! in apm due to incloud process, i.e. aqueous chemistry (so2->so4)

     ! reserve sulfate particle number before aqueous chemistry 
     do n=1,NSO4
       iapm=ip_sulf(k,n,ne)
       XM1D_old(n)=apm_sulf(iapm+ixy)*ug2kg ! ug/m3 -> kg/m3
     enddo
     !!!!!

     if(lapm_pso4) then 
       pso4_so2_1d = so4_aqchem(i03+ixy)
     else
       pso4_so2_1d = 0.0
     endif


     do is=1,NSO4+4
        iccn=ip_fccn(k,is,ne)
        FCLOUD1(is)=frac_ccn(iccn+ixy)
        if(is.le.NSO4) then        ! SP
          iapm=ip_sulf(is,k,ne)
          apm_sulf(iapm+ixy) = apm_sulf(iapm+ixy) + &
                             & pso4_so2_1d*FCLOUD1(is)
        elseif(is.eq.NSO4+1) then  ! BC
          mbcsulf(i03+ixy) = mbcsulf(i03+ixy) + &
                           & pso4_so2_1d*FCLOUD1(is)
        elseif(is.eq.NSO4+2) then  ! OC
          mocsulf(i03+ixy) = mocsulf(i03+ixy) + &
                           & pso4_so2_1d*FCLOUD1(is)
        elseif(is.eq.NSO4+3) then  ! DUST
          mdstsulf(i03+ixy) = mdstsulf(i03+ixy) + & 
                            & pso4_so2_1d*FCLOUD1(is)
        elseif(is.eq.NSO4+4) then  ! SEA SALT
          msltsulf(i03+ixy) = msltsulf(i03+ixy) + &
                            & pso4_so2_1d*FCLOUD1(is)
        endif
     enddo
     !---------- end of aqueous chemistry production ---------
     !========================================================

     ! update total sulfate mass ( NAQPMS )
     tot_sulf(i03+ixy)=msulf*kg2ug


     ! apm3d to apm1d
     do n=1,NSO4
       iapm=ip_sulf(k,n,ne)
       XM1D(n)=apm_sulf(iapm+ixy)*ug2kg ! ug/m3 -> kg/m3
       XN1D(N)=XM1D_old(N)/(DENSULF*VDRY(N))*1.E-9 ! #/cm3
     enddo

     do n=1,NSEA
       iapm=ip_salt(k,n,ne)
       XM1D(NSO4+N)=apm_salt(iapm+ixy)*ug2kg
     enddo

     do n=1,NBCOCT
       iapm=ip_bcoc(k,n,ne)
       MBCOC8(N)=apm_bcoc(iapm+ixy)*ug2kg
     enddo

     do n=1,NDSTB
       iapm=ip_dust(k,n,ne)
       XMDST(n)=apm_dust(iapm+ixy)*ug2kg
     enddo

     MBCS  = mbcsulf(i03+ixy)*ug2kg  ! ug/m3 -> kg/m3
     MOCS  = mocsulf(i03+ixy)*ug2kg  ! ug/m3 -> kg/m3
     MDSTS = mdstsulf(i03+ixy)*ug2kg ! ug/m3 -> kg/m3
     MSALTS= msltsulf(i03+ixy)*ug2kg ! ug/m3 -> kg/m3


     MSPS = sum(XM1D(1:NSO4)) ! kg/m3

     !CACID = 2.0D+9
     !PACID = 1.0D+5

     !=================================================================
     !-----------------------------------------------------------------
     ! YU : these species are not important in pollutant or urban areas
     ! shun : use minimal values below
     PLVSOG1 = 1.d-30
     CLVSOG  = 1.d-30
!     MSULFLV = 1.d-30 ! ?
!     MBCLV   = 1.d-30 ! ?
!     MOCLV   = 1.d-30 ! ?
!     MDSTLV  = 1.d-30 ! ?
!     MSALTLV = 1.d-30 ! ?
     !=================================================================
     if(llovgrow.and.lsovgrow) then
       call cal_lv( SOAT,MSPS,MBCS,MOCS,MDSTS,MSALTS &
                   ,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV &
                   ,SOAT_to_APM )
     else
       MSULFLV = 1.d-30
       MBCLV   = 1.d-30
       MOCLV   = 1.d-30
       MDSTLV  = 1.d-30
       MSALTLV = 1.d-30
       SOAT_to_APM = 1.d-30
     endif
     !=================================================================


     DTAPM = dt_naqpms

     if(IFNUCL.eq.1) then ! IMN ION
       CALL IONRATE0(ISURF, YPSURF, XLON, XLAT, YPR, XQ)
       apm_xq3d(i03+ixy)=XQ
     elseif(IFNUCL.eq.2) then
       XQ=1.E-20
       apm_xq3d(i03+ixy)=XQ
     endif



if(lbox_phy) then
     ! apm box model
     CALL APM_PHYS( i,j,k &
     &      ,NCOAG1,NCOAG2,IACT10,IACT20,IACT30 &   ! out
     &      ,NTEMPOUT1                          &   ! in
     &      ,PRESS,TK,RHIN,XQ,PLVSOG1           &   ! in
     &      ,CACID,PACID                        &   ! inout
     &      ,DTAPM,MMSA,MNIT,MNH4               &   ! in 
     &      ,MBCS, MOCS,MDSTS, MSALTS           &   ! inout
     &      ,MBCOC8, SOAT_to_APM                &   ! in
     &      ,CLVSOG,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV & ! in
     &      ,GFTOT1,GFTOT2,DENWET1,DENWET2      &   ! out
     &      ,XM1D,XN1D,TEMPOUT1,XMDST,FCLOUD1   &   ! inout
     &      ,ZBEXT,ZW,ZG,ZBABS,YBEXT,YW,YG      &
     &      ,rw1d_sulf,rw1d_salt,rw1d_dust,refw1d_bc,refw1d_oc & ! out
     &      ,rgf_sulf_1d,rgf_salt_1d,rgf_dust_1d,rgf_bc_1d,rgf_oc_1d &
     &      ,number_ccn1_1d,number_ccn2_1d,number_ccn3_1d &
     &      ,ztn1d &
     &      ,lbox_nucl,lbox_cond,lbox_cond_other,lbox_coag,lbox_coag_scav &
     &      ,npf_rate &
     &      ,so4_ccn2_1d,dust_ccn2_1d,salt_ccn2_1d,bc_ccn2_1d,oc_ccn2_1d &
     &      ,YCS_1D,SURF_1D) 
     !
endif


    !=========================================
    !> apm1d to apm3d

    !-----------------------------------------
    !----------------- tracers ---------------
     do n=1,NSO4
       iapm=ip_sulf(k,n,ne)
       apm_sulf(iapm+ixy)=XM1D(n)*kg2ug
     enddo

     do n=1,NSEA
       iapm=ip_salt(k,n,ne)
       apm_salt(iapm+ixy)=XM1D(NSO4+N)*kg2ug
     enddo

     do n=1,NBCOCT
       iapm=ip_bcoc(k,n,ne)
       apm_bcoc(iapm+ixy)=MBCOC8(N)*kg2ug
     enddo

     do n=1,NDSTB
       iapm=ip_dust(k,n,ne)
       apm_dust(iapm+ixy)=XMDST(n)*kg2ug
     enddo

     mbcsulf(i03+ixy)  = MBCS*kg2ug
     mocsulf(i03+ixy)  = MOCS*kg2ug
     mdstsulf(i03+ixy) = MDSTS*kg2ug
     msltsulf(i03+ixy) = MSALTS*kg2ug
     !--------- end of tracers  ------------
     !---------------------------------------

     ! sulferic acid
     !apm_cacid(i03+ixy) = CACID
     !
     ! ccn number fraction
     do is=1,NSO4+4
        iccn=ip_fccn(k,is,ne)
        frac_ccn(iccn+ixy)=FCLOUD1(is)
     enddo
     !
     ! CCN number at three supersaturations
     do itype=1,5
       i00_type=ip_type(k,itype,ne)
       number_ccn1(i00_type+ixy)=number_ccn1_1d(itype) ! 0.8% supersaturation
       number_ccn2(i00_type+ixy)=number_ccn2_1d(itype) ! 0.4% supersaturation
       number_ccn3(i00_type+ixy)=number_ccn3_1d(itype) ! 0.2% supersaturation
       ztn3d(i00_type+ixy)=ztn1d(itype)
       ycs3d(i00_type+ixy)=YCS_1D(itype)
       surf3d(i00_type+ixy)=SURF_1D(itype)
     enddo

     so4_ccn2(i03+ixy)  = so4_ccn2_1d
     salt_ccn2(i03+ixy) = salt_ccn2_1d
     dust_ccn2(i03+ixy) = dust_ccn2_1d
     bc_ccn2(i03+ixy)   = bc_ccn2_1d
     oc_ccn2(i03+ixy)   = oc_ccn2_1d

     cn3nm(i03+ixy)=sum(ztn1d(:))
     if(NTEMPOUT1.EQ.1) then
       cn10nm(i03+ixy)=TEMPOUT1(1)
     endif
     !
     ! particle growth factor
     if(lbox_phy) then
      call force_rgfac( 'def_rgf_sulf_1d',i,j,k &
                       ,rgf1,rgf2,rgfdef,rgf_sulf_1d,iostate)
      call force_rgfac( 'def_rgf_salt_1d',i,j,k &
                       ,rgf1,rgf2,rgfdef,rgf_salt_1d,iostate)
      call force_rgfac( 'def_rgf_dust_1d',i,j,k &
                       ,rgf1,rgf2,rgfdef,rgf_dust_1d,iostate)
      call force_rgfac( 'def_rgf_bc_1d',i,j,k &
                       ,rgf1,rgf2,rgfdef,rgf_bc_1d,iostate)
      call force_rgfac( 'def_rgf_oc_1d',i,j,k &
                       ,rgf1,rgf2,rgfdef,rgf_oc_1d,iostate)

      rgf_sulf(i03+ixy) = rgf_sulf_1d
      rgf_salt(i03+ixy) = rgf_salt_1d
      rgf_dust(i03+ixy) = rgf_dust_1d 
      rgf_bc(i03+ixy)   = rgf_bc_1d
      rgf_oc(i03+ixy)   = rgf_oc_1d
     endif

     ! new particle formation
     npf3d(i03+ixy)=npf_rate

!opt+
     do iwls=1,MWLS
       iapm=ip_ext(k,iwls,ne)
       ZBEXT3D(iapm+ixy)=ZBEXT(iwls)
       ZW3D(iapm+ixy)=ZW(iwls)
       do itype=1,5
         iapm=ipwl_type(k,iwls,itype,ne)
         YBEXT3D(iapm+ixy)=YBEXT(itype,iwls)
       enddo
     enddo

     



!opt+

!     ! particle wet size
!     do is=1,NSO4
!        iapm=ip_sulf(k,is,ne)
!        rw_sulf(iapm+ixy)=rw1d_sulf(is)
!     enddo
!     do is=1,NSEA
!        iapm=ip_salt(k,is,ne)
!        rw_salt(iapm+ixy)=rw1d_salt(is)
!     enddo
!     do is=1,NDSTB
!        iapm=ip_dust(k,is,ne)
!        rw_dust(iapm+ixy)=rw1d_dust(is)
!     enddo
!     do imode=1,2
!        iapm=ipmode_bcoc(k,imode,ne)
!        refw_bc(iapm+ixy)=refw1d_bc(imode)
!        refw_oc(iapm+ixy)=refw1d_oc(imode)
!     enddo
     !
     !< end of apm1d to apm3d
     !===============================      

     !> shun : update sulferic vapor concentration 
     !         after nucleation and condensation
     apm_pr_atm = PA2ATM*Plev(i03+ixy)*100.
     apm_te     = t(i03+ixy)
     apm_cair_mlc = apm_avogad*apm_pr_atm/(82.056*apm_te)

     if(i.eq.22.and.j.eq.30.and.k.eq.1.and.ne.eq.1) then
       print*,'acid00',h2so4_gas(i03+ixy),tmp
       print*,'pacid ',PACID*ppbunit/apm_cair_mlc*dt_naqpms,PACID*dt_naqpms
     endif  

     h2so4_gas(i03+ixy)=CACID*ppbunit/apm_cair_mlc

     if(i.eq.22.and.j.eq.30.and.k.eq.1.and.ne.eq.1) then
       print*,'acid01',h2so4_gas(i03+ixy),CACID
       !stop
     endif  

   ENDDO loop_k

 ENDDO loop_i
 ENDDO loop_j

!ENDIF ! apm flag

!stop 'apm_box_ph'

!print*,'kkradius'

return

end subroutine apm_boxphy_driver


