

subroutine apm_wet_dep_v3 &
 & ( myid &
 &  ,lapm &
 &  ,itt &
 &  ,dt  &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dx,dy,dz &
 &  ,CLW,RNW,temp,Plev,CPH &
 &  ,RAINCON,RAINNON &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )

use wdepflags, only: laerincld
use apm_varlist
use aqchem_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

integer :: itt
real    :: dt,dtout

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is,imode

integer :: mem2d,mem3d

real,dimension(mem2d) :: RAINCON,RAINNON

real,dimension(mem3d) :: dx,dy,dz

real,dimension(mem3d) :: CLW,RNW,temp,Plev,CPH

integer :: ixy,i02,i03,iapm,i02apm,iapm_1,iapm_2,iapm_5,iapm_6

integer :: ip3mem(nzz,nest),ip2mem(nest)

!===============================================================
! local vars
logical :: LPREC
real    :: TMASS
integer :: iwb,ieb,jsb,jeb 
integer :: KBOTC,KTOPC
real,dimension(nzz) :: clwc,rnwc,twet,pwet,RR,VOLRAT,pwr_c,cwc_c
real,dimension(nzz) :: con,con_s,conbc_qw,conbc_sw,frcbc_qw

real,dimension(nzz) :: conqw,consw

real,parameter      :: lvalue=0.0,hvalue=2.0
real,dimension(nzz) :: concore1,concore2,ratiocore
real,dimension(nzz) :: concore1_bc,concore1_oc
real,dimension(nzz) :: concore2_bc,concore2_oc
real,dimension(nzz) :: ratiocore_bc,ratiocore_oc



real,dimension(nzz) :: diam1d,density1d

real :: cwph_1d ! 2013-03-06

real :: rgf_tmp

real,parameter :: rgfmin=1.0,rgfmax=100.0
real,parameter :: diamfine=1.0E-6*sqrt(0.1*2.5)
real :: diam
real*8,allocatable,dimension(:) :: numcc,tmp,radius
real*8 :: totalbc,totaloc

! 1:FF 2:BB
integer,parameter :: ridx(1:8) = (/1,2,1,2,1,2,1,2/)

character*2,parameter :: flagbcoc(1:8) = (/'bc','bc','oc','oc' &
                                          ,'bc','bc','oc','oc'/)
!=================================================================


integer :: nbin
real,allocatable,dimension(:,:) :: consp

real :: wcav


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!return



!=======================================================================
! allocate memory for work arrays
iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1
allocate( apm_depfld(iwb:ieb,jsb:jeb), apm_depfld2(iwb:ieb,jsb:jeb) )
allocate( apm_depflds(iwb:ieb,jsb:jeb), apm_depfld2s(iwb:ieb,jsb:jeb) )
!=======================================================================


loop_j : do j = sy(ne),ey(ne)
loop_i : do i = sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
  i02 = ip2mem(ne)

  DO k = 1, nzz
      i03 = ip3mem(k,ne)
      clwc(k) = CLW(i03+ixy)
      rnwc(k) = RNW(i03+ixy)
      TWET(k) = temp(i03+ixy)
      PWET(k) = Plev(i03+ixy)
  ENDDO ! k 

! CCC  GET THE LAYER containin precipitation bottom/top
  CALL GETCLOUDDEPTH ( MYID,TWET,PWET,clwc,rnwc,CWC_C, PWR_C,KBOTC,KTOPC,NZZ &
                      ,RR,VOLRAT,RAINCON(i02+ixy),RAINNON(i02+ixy),LPREC)

!cycle

!**************************
  IF(LPREC) THEN
   if(lfor_h2so4) then
    i02 = ip2mem(ne)
    IF(ITT.EQ.1) THEN
      apm_depfld(i,j)  = 0.0
      apm_depfld2(i,j) = 0.0
    ELSE
      apm_depfld(i,j)  = wdep_h2so4(i02+ixy)
      apm_depfld2(i,j) = wdep2_h2so4(i02+ixy)
    ENDIF

    DO K = KTOPC, KBOTC,-1
      i03 = ip3mem(k,ne)
!      cwph_1d = clw_ph(i03+ixy) 
      cwph_1d = 5.0
      con(k) = h2so4_gas(i03+ixy)
      CALL WETDEP_GAS( MYID, KBOTC, KTOPC, DT, DX(I03+IXY),DY(I03+IXY) &
                      ,DZ(I03+IXY),TWET(K), PWET(K),CWC_C(K), PWR_C(K) &
                      , 0., 0., VOLRAT(K),RR(K),cwph_1d &
                      ,con(k), TMASS, apm_depfld(i,j), apm_depfld2(i,j) &
                      ,wcav, I, J, K, 1, dt   ) ! IG=1 : H2SO4(g)
      h2so4_gas(i03+ixy)=con(k)
    ENDDO

    wdep_h2so4(i02+ixy)=apm_depfld(i,j)
    wdep2_h2so4(i02+ixy)=apm_depfld2(i,j)
   endif
  ENDIF
!***************************************

!  cycle

!***************************************
! apm sulfate

  IF(LPREC) THEN

    !> apm sulfate
    if(lfor_sulf) then
      !if(i.eq.1.and.j.eq.1) print*,'sulf '
      loop_so4 : do is=1,NSO4

        laerincld=.true.
!print*,'so4 size',is

        i02apm=ip2_sulf(is,ne)

        IF(ITT.EQ.1) THEN
          apm_depfld(i,j)  = 0.0
          apm_depfld2(i,j) = 0.0
        ELSE
          apm_depfld(i,j)  = apm_wdep_sulf(i02apm+ixy) 
          apm_depfld2(i,j) = apm_wdep2_sulf(i02apm+ixy)
        ENDIF


        DO k=KTOPC, KBOTC,-1

          i03=ip3mem(k,ne)
          iapm=ip_sulf(k,is,ne)

          con(k) = apm_sulf(iapm+ixy)

          !diam = rd_sulf(is)*2.0 ! : m  ! need to be modified
          if(lapm_wetsize) then
            diam1d(k) = rd_sulf(is)*2.0*rgf_sulf(i03+ixy)
            call check_rgf('rgf_sulf',i,j,k,rgfmin,rgfmax,rgf_sulf(i03+ixy))
          else
            diam1d(k) = rd_sulf(is)*2.0
          endif
          call check_value(i,j,k,0.0,1000.0,diam1d(k))
          density1d(k)=densulf*1.0e6 ! g/cm3 -> g/m3

          diam = diam1d(k)
diam = diamfine
          apmdensity=density1d(k)


          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k),VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)

          apm_sulf(iapm+ixy)=con(k)

        ENDDO

        apm_wdep_sulf(i02apm+ixy) =apm_depfld(i,j)
        apm_wdep2_sulf(i02apm+ixy)=apm_depfld2(i,j)

      enddo loop_so4
    endif
    !< shun : end of apm sulfate 
!*****************************************************

!cycle


!*****************************************************
    !> shun : apm seasalt

    if(lfor_salt) then

      laerincld=.true.

      nbin=NSEA

      concore1=0.0
      concore2=0.0
      ratiocore=1.0

      if(lcoated_dyn) then

        allocate(consp(nbin,nzz))
        allocate(numcc(nbin))
        allocate(tmp(nbin))
        allocate(radius(nbin))

        do k=1,nzz
          i03=ip3mem(k,ne)
          do is=1,nbin
            iapm=ip_salt(k,is,ne)
            radius(is)=rd_salt(is)
            concore1(k)=concore1(k)+apm_salt(iapm+ixy)
            numcc(is)= apm_salt(iapm+ixy)/ &
                 & radius(is)**3.0d0
!            tmp(is)=numcc(is)*radius(is)**2
tmp(is)=apm_salt(iapm+ixy)
          enddo
          do is=1,nbin
            iapm=ip_salt(k,is,ne)
            if(sum(tmp(:)).gt.0) then
              consp(is,k)=msltsulf(i03+ixy)*tmp(is)/sum(tmp(:))
            else
              consp(is,k)=0
            endif
          enddo
        enddo

      endif

      loop_salt : do is=1,NSEA

        i02apm=ip2_salt(is,ne)

        IF(ITT.EQ.1) THEN
          apm_depfld(i,j)  = 0.0
          apm_depfld2(i,j) = 0.0
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = 0.0
            apm_depfld2s(i,j) = 0.0 
          endif
        ELSE
          apm_depfld(i,j)  = apm_wdep_salt(i02apm+ixy)
          apm_depfld2(i,j) = apm_wdep2_salt(i02apm+ixy)
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = apm_wdeps_salt(i02apm+ixy)
            apm_depfld2s(i,j) = apm_wdep2s_salt(i02apm+ixy)
          endif
        ENDIF

        DO k=KTOPC, KBOTC,-1

          i03=ip3mem(k,ne)
          iapm=ip_salt(k,is,ne)

          !diam = rd_salt(is)*2.0   ! need to be modified
          if(lapm_wetsize) then
            diam1d(k)=rd_salt(is)*2.0*rgf_salt(i03+ixy)
            call check_rgf('rgf_salt',i,j,k,rgfmin,rgfmax,rgf_salt(i03+ixy))
          else
            diam1d(k)=rd_salt(is)*2.0
          endif
          call check_value(i,j,k,0.0,1000.0,diam1d(k))
          density1d(k)=densalt*1.0e6 ! g/cm3 -> g/m3

          diam = diam1d(k)
          apmdensity=density1d(k)

          con(k) = apm_salt(iapm+ixy)
          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k),VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_salt(iapm+ixy)=con(k)

!          concore2(k)=concore2(k)+con(k)

          if(lcoated_dyn.and.dep_method.eq.1) then
            con_s(k) = consp(is,k)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k) ,VOLRAT(K) &
                           ,RR(K),con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            consp(is,k) = con_s(k)
          endif

        ENDDO

        DO k=KTOPC, KBOTC,-1
          concore2(k)=concore2(k)+con(k)
        ENDDO

        apm_wdep_salt(i02apm+ixy) =apm_depfld(i,j)
        apm_wdep2_salt(i02apm+ixy)=apm_depfld2(i,j)

        apm_wdeps_salt(i02apm+ixy) =apm_depflds(i,j)
        apm_wdep2s_salt(i02apm+ixy)=apm_depfld2s(i,j)

      enddo loop_salt

       do k=KTOPC, KBOTC,-1
         ratiocore(k)=amin1(amax1(concore2(k)/concore1(k),lvalue),hvalue)
       enddo

      if(lcoated_dyn) then
        do k=1,nzz
          i03=ip3mem(k,ne)
          msltsulf(i03+ixy)=sum(consp(1:nbin,k))
        enddo
        deallocate(consp,numcc,tmp,radius)
      endif


    endif
    !< shun : end of apm seasalt
!*****************************************************

!cycle


!*****************************************************
    !> shun : apm dust

    if(lfor_dust) then

      laerincld=.true.

      nbin=NDSTB

      concore1=0.0
      concore2=0.0
      ratiocore=1.0

      if(lcoated_dyn) then

        allocate(consp(nbin,nzz))
        allocate(numcc(nbin))
        allocate(tmp(nbin))
        allocate(radius(nbin))

        do k=1,nzz
          i03=ip3mem(k,ne) 
          do is=1,nbin
            iapm=ip_dust(k,is,ne)
            radius(is)=rd_dust(is)
            concore1(k)=concore1(k)+apm_dust(iapm+ixy)
            numcc(is)= apm_dust(iapm+ixy)/ &
                 & radius(is)**3.0d0
!            tmp(is)=numcc(is)*radius(is)**2
tmp(is)=apm_dust(iapm+ixy)
          enddo
          do is=1,nbin
            iapm=ip_dust(k,is,ne)
            if(sum(tmp(:)).gt.0) then
              consp(is,k)=mdstsulf(i03+ixy)*tmp(is)/sum(tmp(:))
            else
              consp(is,k)=0
            endif
          enddo
        enddo

      endif



      loop_dust : do is=1,NDSTB

        i02apm=ip2_dust(is,ne)

        IF(ITT.EQ.1) THEN
          apm_depfld(i,j)  = 0.0
          apm_depfld2(i,j) = 0.0
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = 0.0
            apm_depfld2s(i,j) = 0.0
          endif
        ELSE
          apm_depfld(i,j)  = apm_wdep_dust(i02apm+ixy)
          apm_depfld2(i,j) = apm_wdep2_dust(i02apm+ixy)
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = apm_wdeps_dust(i02apm+ixy)
            apm_depfld2s(i,j) = apm_wdep2s_dust(i02apm+ixy)
          endif
        ENDIF

        DO k=KTOPC, KBOTC,-1

          i03=ip3mem(k,ne)
          iapm=ip_dust(k,is,ne)
         
          if(lapm_wetsize) then
            diam1d(k) = rd_dust(is)*2.0*rgf_dust(i03+ixy) ! need to be modified
            call check_rgf('rgf_dust',i,j,k,rgfmin,rgfmax,rgf_dust(i03+ixy))
          else
            diam1d(k) = rd_dust(is)*2.0
          endif
          call check_value(i,j,k,0.0,1000.0,diam1d(k))
          density1d(k)=dendust*1.0e6 ! g/cm3 -> g/m3

          diam = diam1d(k)
          apmdensity=density1d(k)

          con(k) = apm_dust(iapm+ixy)
          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k),VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_dust(iapm+ixy)=con(k)

          if(lcoated_dyn.and.dep_method.eq.1) then
            con_s(k) = consp(is,k)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K), clwc(k),rnwc(k),VOLRAT(K) &
                           ,RR(K),con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            consp(is,k)=con_s(k)
          endif

        ENDDO

DO k=KTOPC, KBOTC,-1
  concore2(k)=concore2(k)+con(k)
ENDDO

        apm_wdep_dust(i02apm+ixy) =apm_depfld(i,j)
        apm_wdep2_dust(i02apm+ixy)=apm_depfld2(i,j)

        apm_wdeps_dust(i02apm+ixy) =apm_depflds(i,j)
        apm_wdep2s_dust(i02apm+ixy)=apm_depfld2s(i,j)

      enddo loop_dust

      do k=KTOPC, KBOTC,-1
         ratiocore(k)=amin1(amax1(concore2(k)/concore1(k),lvalue),hvalue)
      enddo

      if(lcoated_dyn) then
        do k=1,nzz
          i03=ip3mem(k,ne)
          mdstsulf(i03+ixy)=sum(consp(1:nbin,k))
        enddo
        deallocate(consp,numcc,tmp,radius)
      endif

    endif
    !< shun : end of apm dust
!**************************************************************

!cycle


!**************************************************************
    !> shun : apm bcoc

    if(lfor_bcoc) then

      do k=1,nzz
          iapm_1=ip_bcoc(k,1,ne)
          iapm_2=ip_bcoc(k,2,ne)
          iapm_5=ip_bcoc(k,5,ne)
          iapm_6=ip_bcoc(k,6,ne)
          conbc_qw(k)=apm_bcoc(iapm_1+ixy)+apm_bcoc(iapm_2+ixy)
          conbc_sw(k)=apm_bcoc(iapm_5+ixy)+apm_bcoc(iapm_6+ixy)
      enddo


      nbin=NBCOCT

      concore1_bc=0.0
      concore2_bc=0.0
      ratiocore_bc=1.0

      concore1_oc=0.0
      concore2_oc=0.0
      ratiocore_oc=1.0

      if(lcoated_dyn) then

        allocate(consp(nbin,nzz))
        allocate(numcc(nbin))
        allocate(tmp(nbin))
        allocate(radius(nbin))

        do k=1,nzz
          i03=ip3mem(k,ne)
          do is=1,NBCOCT
            iapm=ip_bcoc(k,is,ne)
            radius(is)=ref1d_bcoc(ridx(is))
            numcc(is)= apm_bcoc(iapm+ixy)/ &
                 & radius(is)**3.0d0
!            tmp(is)=numcc(is)*radius(is)**2
tmp(is)=apm_bcoc(iapm+ixy)
            if(flagbcoc(is).eq.'bc') then
              concore1_bc(k)=concore1_bc(k)+apm_bcoc(iapm+ixy)
            elseif(flagbcoc(is).eq.'oc') then
              concore1_oc(k)=concore1_oc(k)+apm_bcoc(iapm+ixy)
            endif
          enddo
          totalbc=tmp(1)+tmp(5)+tmp(2)+tmp(6)
          totaloc=tmp(3)+tmp(7)+tmp(4)+tmp(8)
          do is=1,NBCOCT
            iapm=ip_bcoc(k,is,ne)
            if(flagbcoc(is).eq.'bc') then
              if(totalbc.gt.0) then
                consp(is,k)=mbcsulf(i03+ixy)*tmp(is)/totalbc
              else
                consp(is,k)=0
              endif
            elseif(flagbcoc(is).eq.'oc') then
              if(totaloc.gt.0) then
                 consp(is,k)=mocsulf(i03+ixy)*tmp(is)/totaloc
              else
                 consp(is,k)=0
              endif
            endif
          enddo
        enddo !k

      endif

      loop_bcoc : do is=1,NBCOCT

!        if(is.eq.5.or.is.eq.6) cycle ! skip hydrophobic bc 
        if(is.eq.5.or.is.eq.6) then
           laerincld=.false.
        else
           laerincld=.true.
        endif

        i02apm=ip2_bcoc(is,ne)

        IF(ITT.EQ.1) THEN
          apm_depfld(i,j)  = 0.0
          apm_depfld2(i,j) = 0.0
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = 0.0
            apm_depfld2s(i,j) = 0.0
          endif
        ELSE
          apm_depfld(i,j)  = apm_wdep_bcoc(i02apm+ixy)
          apm_depfld2(i,j) = apm_wdep2_bcoc(i02apm+ixy)
          if(lcoated_dyn.and.dep_method.eq.1) then
!            apm_depflds(i,j)  = apm_wdeps_bcoc(i02apm+ixy)
!            apm_depfld2s(i,j) = apm_wdep2s_bcoc(i02apm+ixy)
          endif
        ENDIF

        DO k=KTOPC, KBOTC,-1

          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)

          if(flagbcoc(is).eq.'bc') then
             rgf_tmp=rgf_bc(i03+ixy)
          elseif(flagbcoc(is).eq.'oc') then
             rgf_tmp=rgf_oc(i03+ixy)
          endif

          if(lapm_wetsize) then
            diam1d(k) = ref1d_bcoc(ridx(is))*2.0*rgf_tmp ! need to be modified
            call check_rgf('rgf_bcoc',i,j,k,rgfmin,rgfmax,rgf_tmp)
          else
            diam1d(k) = ref1d_bcoc(ridx(is))*2.0
          endif
          call check_value(i,j,k,0.0,1000.0,diam1d(k))
          density1d(k)=denbcoc*1.0e6        ! g/cm3 -> g/m3

          diam = diam1d(k)

diam=diamfine

          apmdensity=density1d(k)


          con(k) = apm_bcoc(iapm+ixy)
          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k) ,VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_bcoc(iapm+ixy)=con(k)

          if(lcoated_dyn.and.dep_method.eq.1) then
!print*,'iapm+ixy=',iapm+ixy
!print*,'con_s(k)=',con_s(k)
            con_s(k) = consp(is,k)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K), clwc(k),rnwc(k),VOLRAT(K) &
                           ,RR(K), con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            consp(is,k)=con_s(k)
          endif

        ENDDO

if(flagbcoc(is).eq.'bc') then
  DO k=KTOPC, KBOTC,-1
    concore2_bc(k)=concore2_bc(k)+con(k)
  ENDDO
elseif(flagbcoc(is).eq.'oc') then
  DO k=KTOPC, KBOTC,-1
    concore2_oc(k)=concore2_oc(k)+con(k)
  ENDDO
endif
        apm_wdep_bcoc(i02apm+ixy)=apm_depfld(i,j)
        apm_wdep2_bcoc(i02apm+ixy)=apm_depfld2(i,j)

!        apm_wdeps_bcoc(i02apm+ixy)=apm_depflds(i,j)
!        apm_wdep2s_bcoc(i02apm+ixy)=apm_depfld2s(i,j)

      enddo loop_bcoc

      do k=KTOPC, KBOTC,-1
          ratiocore_bc(k)=amin1(amax1(concore2_bc(k)/concore1_bc(k),lvalue),hvalue)
          ratiocore_oc(k)=amin1(amax1(concore2_oc(k)/concore1_oc(k),lvalue),hvalue)
      enddo

      if(lcoated_dyn) then
        do k=1,nzz
           i03=ip3mem(k,ne)
!           mbcsulf(i03+ixy)=consp(1,k)+consp(5,k)+consp(2,k)+consp(6,k)
!           mocsulf(i03+ixy)=consp(3,k)+consp(7,k)+consp(4,k)+consp(8,k)
        enddo
        deallocate(consp,numcc,tmp,radius)
      endif


    endif
    !< shun : end of apm bcoc
!*********************************************************************

      nbin=nbincb

      concore1=0.0
      concore2=0.0
      ratiocore=1.0

      if(lcoated_dyn) then

        allocate(consp(nbin,nzz))
        allocate(numcc(nbin))
        allocate(tmp(nbin))
        allocate(radius(nbin))

        do k=1,nzz
          i03=ip3mem(k,ne)

          if(conbc_qw(k)+conbc_sw(k).gt.0) then
             frcbc_qw(k)=conbc_qw(k)/(conbc_qw(k)+conbc_sw(k))
          else
             frcbc_qw(k)=0.0
          endif
          if(frcbc_qw(k).gt.1) frcbc_qw(k)=0.0

          do is=1,nbin
            iapm=ip_cbbin(k,is,ne)
            radius(is)=rd_binbc(is)
!radius(is)=diamfine
            concore1(k)=concore1(k)+apm_binbc(iapm+ixy)
            numcc(is)= apm_binbc(iapm+ixy)/ &
                 & radius(is)**3.0d0
!            tmp(is)=numcc(is)*radius(is)**2
tmp(is)=apm_binbc(iapm+ixy)*frcbc_qw(k)
          enddo
          do is=1,nbin
            iapm=ip_cbbin(k,is,ne)
            if(sum(tmp(:)).gt.0) then
              consp(is,k)=mbcsulf(i03+ixy)*tmp(is)/sum(tmp(:))
            else
              consp(is,k)=0
            endif
          enddo
        enddo

      endif

      loop_binbc : do is=1,nbincb

        i02apm=ip2_binbc(is,ne)

        IF(ITT.EQ.1) THEN
          apm_depfld(i,j)  = 0.0
          apm_depfld2(i,j) = 0.0
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = 0.0
            apm_depfld2s(i,j) = 0.0
          endif
        ELSE
          apm_depfld(i,j)  = apm_wdep_binbc(i02apm+ixy)
          apm_depfld2(i,j) = apm_wdep2_binbc(i02apm+ixy)
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = apm_wdeps_binbc(i02apm+ixy)
            apm_depfld2s(i,j) = apm_wdep2s_binbc(i02apm+ixy)
          endif
        ENDIF

        DO k=KTOPC, KBOTC,-1

          i03=ip3mem(k,ne)
          iapm=ip_cbbin(k,is,ne)

          if(lapm_wetsize) then
            diam1d(k) = diamfine ! need to be modified
          else
            diam1d(k) = diamfine
          endif
          density1d(k)=denbcoc*1.0e6 ! g/cm3 -> g/m3

          diam = diam1d(k)
diam=diamfine

          apmdensity=density1d(k)

if(1==2) then
          con(k) = apm_binbc(iapm+ixy)*frcbc_qw(k)


          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k) ,VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_binbc(iapm+ixy)=con(k)+apm_binbc(iapm+ixy)*(1.0-frcbc_qw(k))

          if(lcoated_dyn.and.dep_method.eq.1) then
            con_s(k) = consp(is,k)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K), clwc(k),rnwc(k),VOLRAT(K) &
                           ,RR(K),con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            consp(is,k)=con_s(k)
          endif
else

          conqw(k)=apm_binbc(iapm+ixy)*frcbc_qw(k)
          consw(k)=apm_binbc(iapm+ixy)*(1.0-frcbc_qw(k))

          laerincld=.true.
          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k),VOLRAT(K) &
                         ,RR(K),conqw(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          laerincld=.false.
          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k),VOLRAT(K) &
                         ,RR(K),consw(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_binbc(iapm+ixy)=conqw(k)+consw(k)

          if(lcoated_dyn.and.dep_method.eq.1) then
            laerincld=.true.
            con_s(k) = consp(is,k)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K), clwc(k),rnwc(k),VOLRAT(K) &
                           ,RR(K),con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            consp(is,k)=con_s(k)
          endif
endif


        ENDDO

DO k=KTOPC, KBOTC,-1
  concore2(k)=concore2(k)+con(k)
ENDDO

        apm_wdep_binbc(i02apm+ixy) =apm_depfld(i,j)
        apm_wdep2_binbc(i02apm+ixy)=apm_depfld2(i,j)

        apm_wdeps_binbc(i02apm+ixy) =apm_depflds(i,j)
        apm_wdep2s_binbc(i02apm+ixy)=apm_depfld2s(i,j)

      enddo loop_binbc

      do k=KTOPC, KBOTC,-1
         ratiocore(k)=amin1(amax1(concore2(k)/concore1(k),lvalue),hvalue)
      enddo

      if(lcoated_dyn) then
        do k=1,nzz
          i03=ip3mem(k,ne)
          mbcsulf(i03+ixy)=sum(consp(1:nbin,k))
        enddo
        deallocate(consp,numcc,tmp,radius)
      endif



      laerincld=.true.

      nbin=nbincb

      concore1=0.0
      concore2=0.0
      ratiocore=1.0

      if(lcoated_dyn) then

        allocate(consp(nbin,nzz))
        allocate(numcc(nbin))
        allocate(tmp(nbin))
        allocate(radius(nbin))

        do k=1,nzz
          i03=ip3mem(k,ne)
          do is=1,nbin
            iapm=ip_cbbin(k,is,ne)
            radius(is)=rd_binoc(is)
!radius(is)=diamfine
            concore1(k)=concore1(k)+apm_binoc(iapm+ixy)
            numcc(is)= apm_binoc(iapm+ixy)/ &
                 & radius(is)**3.0d0
!            tmp(is)=numcc(is)*radius(is)**2
tmp(is)=apm_binoc(iapm+ixy)
          enddo
          do is=1,nbin
            iapm=ip_cbbin(k,is,ne)
            if(sum(tmp(:)).gt.0) then
              consp(is,k)=mocsulf(i03+ixy)*tmp(is)/sum(tmp(:))
            else
              consp(is,k)=0
            endif
          enddo
        enddo

      endif

      loop_binoc : do is=1,nbincb

        i02apm=ip2_binoc(is,ne)

        IF(ITT.EQ.1) THEN
          apm_depfld(i,j)  = 0.0
          apm_depfld2(i,j) = 0.0
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = 0.0
            apm_depfld2s(i,j) = 0.0
          endif
        ELSE
          apm_depfld(i,j)  = apm_wdep_binoc(i02apm+ixy)
          apm_depfld2(i,j) = apm_wdep2_binoc(i02apm+ixy)
          if(lcoated_dyn.and.dep_method.eq.1) then
            apm_depflds(i,j)  = apm_wdeps_binoc(i02apm+ixy)
            apm_depfld2s(i,j) = apm_wdep2s_binbc(i02apm+ixy)
          endif
        ENDIF

        DO k=KTOPC, KBOTC,-1

          i03=ip3mem(k,ne)
          iapm=ip_cbbin(k,is,ne)

          if(lapm_wetsize) then
            diam1d(k) = diamfine ! need to be modified
          else
            diam1d(k) = diamfine
          endif
          density1d(k)=denbcoc*1.0e6 ! g/cm3 -> g/m3

          diam = diam1d(k)

diam=diamfine

          apmdensity=density1d(k)

          con(k) = apm_binoc(iapm+ixy)
          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k) ,VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_binoc(iapm+ixy)=con(k)

          if(lcoated_dyn.and.dep_method.eq.1) then
            con_s(k) = consp(is,k)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k), VOLRAT(K) &
                           ,RR(K),con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            consp(is,k)=con_s(k)
          endif

        ENDDO

DO k=KTOPC, KBOTC,-1
  concore2(k)=concore2(k)+con(k)
ENDDO

        apm_wdep_binoc(i02apm+ixy) =apm_depfld(i,j)
        apm_wdep2_binoc(i02apm+ixy)=apm_depfld2(i,j)

        apm_wdeps_binoc(i02apm+ixy) =apm_depflds(i,j)
        apm_wdep2s_binoc(i02apm+ixy)=apm_depfld2s(i,j)

      enddo loop_binoc

      do k=KTOPC, KBOTC,-1
         ratiocore(k)=amin1(amax1(concore2(k)/concore1(k),lvalue),hvalue)
      enddo

      if(lcoated_dyn) then
        do k=1,nzz
          i03=ip3mem(k,ne)
          mocsulf(i03+ixy)=sum(consp(1:nbin,k))
        enddo
        deallocate(consp,numcc,tmp,radius)
      endif





  ENDIF ! LPRC

enddo loop_i
enddo loop_j


!stop 'kk subwet sulf'


! deallocate memory for work arrays
if(allocated(apm_depfld))   deallocate(apm_depfld)
if(allocated(apm_depfld2))  deallocate(apm_depfld2)
if(allocated(apm_depflds))  deallocate(apm_depflds)
if(allocated(apm_depfld2s)) deallocate(apm_depfld2s)
!




!stop 'end of apm_wet_dep'

end subroutine apm_wet_dep_v3





