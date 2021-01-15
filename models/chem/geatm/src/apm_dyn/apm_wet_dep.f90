

subroutine apm_wet_dep &
 & ( myid &
 &  ,lapm &
 &  ,itt &
 &  ,dt &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dx,dy,dz &
 &  ,CLW,RNW,temp,Plev,CPH &
 &  ,RAINCON,RAINNON &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )


use apm_varlist
use aqchem_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

integer :: itt
real    :: dt

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is,imode

integer :: mem2d,mem3d

real,dimension(mem2d) :: RAINCON,RAINNON

real,dimension(mem3d) :: dx,dy,dz

real,dimension(mem3d) :: CLW,RNW,temp,Plev,CPH

integer :: ixy,i02,i03,iapm,i02apm

integer :: ip3mem(nzz,nest),ip2mem(nest)

!===============================================================
! local vars
logical :: LPREC
real    :: TMASS
integer :: iwb,ieb,jsb,jeb 
integer :: KBOTC,KTOPC
real,dimension(nzz) :: clwc,rnwc,twet,pwet,RR,VOLRAT,pwr_c,cwc_c
real,dimension(nzz) :: con,con_s

real,dimension(nzz) :: diam1d,density1d

real :: cwph_1d ! 2013-03-06

real :: rgf_tmp

real,parameter :: rgfmin=1.0,rgfmax=100.0

real :: diam
real*8,allocatable,dimension(:) :: numcc,tmp,radius
real*8 :: totalbc,totaloc

! 1:FF 2:BB
integer,parameter :: ridx(1:8) = (/1,2,1,2,1,2,1,2/)

character*2,parameter :: flagbcoc(1:8) = (/'bc','bc','oc','oc' &
                                          ,'bc','bc','oc','oc'/)
!=================================================================


!print*,'so4-mem',ip_sulf

!stop


! distribute bulk seasalt coated sulfate to bins
   if(.not.allocated(apm_msltsulf)) then
      allocate(apm_msltsulf(saltmem))
   else
      stop 'allocate(apm_msltsulf) erro'
   endif
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius))) then
      allocate(numcc(NSEA))
      allocate(tmp(NSEA))
      allocate(radius(NSEA))
   else
      stop 'allocate(numcc) erro'
   endif

   DO k=1,nzz
    if(k.eq.1.and.myid.eq.0) print*,'salt so4 bulk2bins'
    i03=ip3mem(k,ne)
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NSEA
        iapm=ip_salt(k,is,ne)
        radius(is)=rd_salt(is)
        numcc(is)= apm_salt(iapm+ixy)/ &
                 & radius(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius(is)**3.0d0)
        tmp(is)=numcc(is)*radius(is)**2
      enddo
      do is=1,NSEA
        iapm=ip_salt(k,is,ne)
        if(sum(tmp(:)).gt.0) then
          apm_msltsulf(iapm+ixy)=msltsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
          apm_msltsulf(iapm+ixy)=0
        endif
      enddo
    enddo
    enddo
   ENDDO
   if(allocated(numcc).and.allocated(tmp).and.allocated(radius)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius)
   else
      stop 'deallocate(numcc) erro'
   endif
   !!!!!


! distribute bulk dust coated sulfate to bins
   if(.not.allocated(apm_mdstsulf)) then
      allocate(apm_mdstsulf(dustmem))
   else
      stop 'allocate(apm_mdstsulf) erro'
   endif
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius))) then
      allocate(numcc(NDSTB))
      allocate(tmp(NDSTB))
      allocate(radius(NDSTB))
   else
      stop 'allocate(numcc) erro'
   endif

   DO k=1,nzz
    if(k.eq.1.and.myid.eq.0) print*,'dust so4 bulk2bins'
    i03=ip3mem(k,ne)
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NDSTB
        iapm=ip_dust(k,is,ne)
        radius(is)=rd_dust(is)
        numcc(is)= apm_dust(iapm+ixy)/ &
                 & radius(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius(is)**3.0d0)
        tmp(is)=numcc(is)*radius(is)**2
      enddo
      do is=1,NDSTB
        iapm=ip_dust(k,is,ne)
        if(sum(tmp(:)).gt.0) then
          apm_mdstsulf(iapm+ixy)=mdstsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
          apm_mdstsulf(iapm+ixy)=0
        endif
      enddo
    enddo
    enddo
   ENDDO
   if(allocated(numcc).and.allocated(tmp).and.allocated(radius)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius)
   else
      stop 'deallocate(numcc) erro'
   endif
   !!!!!


! distribute bulk BC coated sulfate to bins
   if(.not.allocated(apm_mbcocsulf)) then
      allocate(apm_mbcocsulf(bcocmem))
   else
      stop 'allocate(apm_mbcocsulf) erro'
   endif
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius))) then
      allocate(numcc(NBCOCT))
      allocate(tmp(NBCOCT))
      allocate(radius(NBCOCT))
   else
      stop 'allocate(numcc) erro'
   endif
   DO k=1,nzz
    if(k.eq.1.and.myid.eq.0) print*,'bcoc so4 bulk2bins'
    i03=ip3mem(k,ne)
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
        radius(is)=ref1d_bcoc(ridx(is))
        numcc(is)= apm_bcoc(iapm+ixy)/ &
                 & radius(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius(is)**3.0d0)
        tmp(is)=numcc(is)*radius(is)**2
      enddo
      totalbc=tmp(1)+tmp(5)+tmp(2)+tmp(6)
      totaloc=tmp(3)+tmp(7)+tmp(4)+tmp(8)
      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
        if(flagbcoc(is).eq.'bc') then
          if(totalbc.gt.0) then
            apm_mbcocsulf(iapm+ixy)=mbcsulf(i03+ixy)*tmp(is)/totalbc
          else
            apm_mbcocsulf(iapm+ixy)=0
          endif
        elseif(flagbcoc(is).eq.'oc') then
          if(totaloc.gt.0) then
            apm_mbcocsulf(iapm+ixy)=mocsulf(i03+ixy)*tmp(is)/totaloc
          else
            apm_mbcocsulf(iapm+ixy)=0
          endif
        else
          stop 'erro'
        endif
      enddo
    enddo
    enddo
   ENDDO
   if(allocated(numcc).and.allocated(tmp).and.allocated(radius)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius)
   else
      stop 'deallocate(numcc) erro'
   endif
   !!!!!
!< end of bulk to bins
!======================================================================



! allocate memory for work arrays
iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1
allocate( apm_depfld(iwb:ieb,jsb:jeb), apm_depfld2(iwb:ieb,jsb:jeb) )
allocate( apm_depflds(iwb:ieb,jsb:jeb), apm_depfld2s(iwb:ieb,jsb:jeb) )
!

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
      cwph_1d = clw_ph(i03+ixy) 
      con(k) = h2so4_gas(i03+ixy)
      CALL WETDEP_GAS( MYID, KBOTC, KTOPC, DT, DX(I03+IXY),DY(I03+IXY) &
                      ,DZ(I03+IXY),TWET(K), PWET(K),CWC_C(K), PWR_C(K) &
                      , 0., 0., VOLRAT(K),RR(K),cwph_1d &
                      ,con(k), TMASS, apm_depfld(i,j), apm_depfld2(i,j) &
                      , I, J, K, 1   ) ! IG=1 : H2SO4(g)
      h2so4_gas(i03+ixy)=con(k)
    ENDDO

    wdep_h2so4(i02+ixy)=apm_depfld(i,j)
    wdep2_h2so4(i02+ixy)=apm_depfld2(i,j)
   endif
  ENDIF


  IF(LPREC) THEN

    !> apm sulfate
    if(lfor_sulf) then
      !if(i.eq.1.and.j.eq.1) print*,'sulf '
      loop_so4 : do is=1,NSO4

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


!print*,'i&j&k',i,j,k
!print*,'diam=',diam1d(k),rd_sulf(is)*2.0,rgf_sulf(i03+ixy)

          !print*,'check sulf'
          !call check_value(i,j,k,0.0,1000.0,rw_sulf(iapm+ixy))
          !diam1d(k) = rw_sulf(iapm+ixy)*2.0
          !print*,'check diam1d'
          !call check_value(i,j,k,0.0,1000.0,diam1d(k))

          diam = diam1d(k)
          apmdensity=density1d(k)


          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K), VOLRAT(K) &
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

!cycle

    !> shun : apm seasalt

    if(lfor_salt) then
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
                         ,DZ(I03+IXY),TWET(K), PWET(K), VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_salt(iapm+ixy)=con(k)

          if(lcoated_dyn.and.dep_method.eq.1) then
            con_s(k) = apm_msltsulf(iapm+ixy)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K), VOLRAT(K) &
                           ,RR(K),con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            apm_msltsulf(iapm+ixy)=con_s(k)
          endif

        ENDDO

        apm_wdep_salt(i02apm+ixy) =apm_depfld(i,j)
        apm_wdep2_salt(i02apm+ixy)=apm_depfld2(i,j)

        apm_wdeps_salt(i02apm+ixy) =apm_depflds(i,j)
        apm_wdep2s_salt(i02apm+ixy)=apm_depfld2s(i,j)

      enddo loop_salt
    endif
    !< shun : end of apm seasalt

!cycle

    !> shun : apm dust

    if(lfor_dust) then
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
                         ,DZ(I03+IXY),TWET(K), PWET(K), VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_dust(iapm+ixy)=con(k)

          if(lcoated_dyn.and.dep_method.eq.1) then
            con_s(k) = apm_mdstsulf(iapm+ixy)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K), VOLRAT(K) &
                           ,RR(K),con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            apm_mdstsulf(iapm+ixy)=con_s(k)
          endif

        ENDDO

        apm_wdep_dust(i02apm+ixy) =apm_depfld(i,j)
        apm_wdep2_dust(i02apm+ixy)=apm_depfld2(i,j)

        apm_wdeps_dust(i02apm+ixy) =apm_depflds(i,j)
        apm_wdep2s_dust(i02apm+ixy)=apm_depfld2s(i,j)

      enddo loop_dust
    endif
    !< shun : end of apm dust

!cycle

    !> shun : apm bcoc

    if(lfor_bcoc) then
      loop_bcoc : do is=1,NBCOCT

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
            apm_depflds(i,j)  = apm_wdeps_bcoc(i02apm+ixy)
            apm_depfld2s(i,j) = apm_wdep2s_bcoc(i02apm+ixy)
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
          apmdensity=density1d(k)


          con(k) = apm_bcoc(iapm+ixy)
          CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                         ,DZ(I03+IXY),TWET(K), PWET(K), VOLRAT(K) &
                         ,RR(K),con(k),diam, TMASS &
                         ,apm_depfld(i,j),apm_depfld2(i,j) &
                         ,i,j,k,apmdensity)
          apm_bcoc(iapm+ixy)=con(k)

          if(lcoated_dyn.and.dep_method.eq.1) then
!print*,'iapm+ixy=',iapm+ixy
!print*,'con_s(k)=',con_s(k)
            con_s(k) = apm_mbcocsulf(iapm+ixy)
            CALL wetdep_apm(MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY) &
                           ,DZ(I03+IXY),TWET(K), PWET(K), VOLRAT(K) &
                           ,RR(K), con_s(k),diam, TMASS &
                           ,apm_depflds(i,j),apm_depfld2s(i,j) &
                           ,i,j,k,apmdensity)
            apm_mbcocsulf(iapm+ixy)=con_s(k)
          endif

        ENDDO

        apm_wdep_bcoc(i02apm+ixy)=apm_depfld(i,j)
        apm_wdep2_bcoc(i02apm+ixy)=apm_depfld2(i,j)

        apm_wdeps_bcoc(i02apm+ixy)=apm_depflds(i,j)
        apm_wdep2s_bcoc(i02apm+ixy)=apm_depfld2s(i,j)

      enddo loop_bcoc
    endif
    !< shun : end of apm bcoc

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




!==============!
! bins to bulk !
!==============!

DO k=1,nzz

 i03=ip3mem(k,ne)

 !-> salt
 if(k.eq.1.and.myid.eq.0) print*,'salt so4 bins2bulk'
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    msltsulf(i03+ixy)=0.0d0
    do is=1,NSEA
      iapm=ip_salt(k,is,ne)
      msltsulf(i03+ixy)=msltsulf(i03+ixy)+apm_msltsulf(iapm+ixy)
    enddo
 enddo
 enddo

 !-> dust
 if(k.eq.1.and.myid.eq.0) print*,'dust so4 bins2bulk'
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    mdstsulf(i03+ixy)=0.0d0
    do is=1,NDSTB
      iapm=ip_dust(k,is,ne)
      mdstsulf(i03+ixy)=mdstsulf(i03+ixy)+apm_mdstsulf(iapm+ixy)
    enddo
 enddo
 enddo

 !-> BCOC 
 if(k.eq.1.and.myid.eq.0) print*,'salt bcoc bins2bulk'
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    mbcsulf(i03+ixy)=0.0d0
    mocsulf(i03+ixy)=0.0d0
    do is=1,NBCOCT
      iapm=ip_bcoc(k,is,ne)
      if(flagbcoc(is).eq.'bc') then
        mbcsulf(i03+ixy)=mbcsulf(i03+ixy)+apm_mbcocsulf(iapm+ixy)
      elseif(flagbcoc(is).eq.'oc') then
        mocsulf(i03+ixy)=mocsulf(i03+ixy)+apm_mbcocsulf(iapm+ixy)
      else
        stop 'erro'
      endif
    enddo
 enddo
 enddo

ENDDO ! k

if(allocated(apm_msltsulf)) then
   deallocate(apm_msltsulf)
else
   stop 'deallocate(apm_msltsulf) erro'
endif

if(allocated(apm_mdstsulf)) then
   deallocate(apm_mdstsulf)
else
   stop 'deallocate(apm_mdstsulf) erro'
endif


if(allocated(apm_mbcocsulf)) then
   deallocate(apm_mbcocsulf)
else
   stop 'deallocate(apm_mbcocsulf) erro'
endif


!stop 'www'



!> shun : deposited mass from bins to bulk
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)

    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    i02 = ip2mem(ne)

    wdep_sltsulf(i02+ixy)=0
    wdep2_sltsulf(i02+ixy)=0
    wdep_dstsulf(i02+ixy)=0
    wdep2_dstsulf(i02+ixy)=0
    wdep_bcsulf(i02+ixy)=0
    wdep2_bcsulf(i02+ixy)=0
    wdep_ocsulf(i02+ixy)=0
    wdep2_ocsulf(i02+ixy)=0

    do is=1,NSEA
      iapm=ip2_salt(is,ne)
      wdep_sltsulf(i02+ixy)=wdep_sltsulf(i02+ixy)+apm_wdeps_salt(iapm+ixy)
      wdep2_sltsulf(i02+ixy)=wdep2_sltsulf(i02+ixy)+apm_wdep2s_salt(iapm+ixy)
    enddo

    do is=1,NDSTB
      iapm=ip2_dust(is,ne)
      wdep_dstsulf(i02+ixy)=wdep_dstsulf(i02+ixy)+apm_wdeps_dust(iapm+ixy)
      wdep2_dstsulf(i02+ixy)=wdep2_dstsulf(i02+ixy)+apm_wdep2s_dust(iapm+ixy)
    enddo

    do is=1,NBCOCT
      iapm=ip2_bcoc(is,ne)
      if(flagbcoc(is).eq.'bc') then
        wdep_bcsulf(i02+ixy)=wdep_bcsulf(i02+ixy)+apm_wdeps_bcoc(iapm+ixy)
        wdep2_bcsulf(i02+ixy)=wdep2_bcsulf(i02+ixy)+apm_wdep2s_bcoc(iapm+ixy)
      elseif(flagbcoc(is).eq.'oc') then
        wdep_ocsulf(i02+ixy)=wdep_ocsulf(i02+ixy)+apm_wdeps_bcoc(iapm+ixy)
        wdep2_ocsulf(i02+ixy)=wdep2_ocsulf(i02+ixy)+apm_wdep2s_bcoc(iapm+ixy)
      else
        stop 'erro'
      endif
    enddo

 enddo
 enddo
!< end of deposited mass from bins to bulk


!stop 'end of apm_wet_dep'

end subroutine apm_wet_dep





