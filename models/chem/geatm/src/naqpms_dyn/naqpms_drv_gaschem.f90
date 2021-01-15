!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccc     Gas  Chemistry with CBM-Z
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine naqpms_drv_gaschem &
 & ( myid &
 &  ,dt_cbmz &
 &  ,lprocess,iPrintTermGas &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,iyear2,imonth2,iday2,ihour2,iminute2 &
 &  ,mem2d &
 &  ,mem3d &
 &  ,mem4d &
 &  ,GC_MOLWT,GC_NAME &
 &  ,igas,igasCBM,iaer,isize,nseacom,ndustcom &
 &  ,NLAY_EM &
 &  ,ifsm,ifsmt,idmSet,ismMax,igMark )

use naqpms_varlist
use naqpms_gridinfo, only : LATITCRS,LONGICRS,TERRAIN,HGT1,heiz,dz
use met_fields, only : t,rh1,Plev,coefcld3d
use work_vars, only : RK_HETSO2_DUST,RK_HETHNO3_DUST

use apm_varlist
!implicit none

use smpsulf_var, only : idx_smpsulf_h2so4,idx_smpsulf_h2o2 &
                       ,idx_smpsulf_so2,idx_smpsulf_dms &
                       ,idx_smpsulf_nh3
use smpsulf_var, only : idx_oxdt_oh,idx_oxdt_ho2,idx_oxdt_no3,idx_oxdt_o3
use smpsulf_var, only : ip4mem_oxdt,mmean_oxdt,nmoxdt,nm12
use smpsulf_var, only : ip4mem_ox3d,oxdt3d
use smpsulf_var, only : igassmp

use smpsulf_var, only : tmp4d

include 'chm1.inc'
include 'gas1.inc'

real,parameter ::  STDATMPA = 101325.0 ! standard atmosphere  [ Pa ]
real,parameter ::  PA2ATM = 1.0 / STDATMPA

integer :: myid

real :: dt

!integer :: ichemgas

logical :: lprocess
integer :: iPrintTermGas

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: NLAY_EM

integer :: time(5)

integer :: iyear1,imonth1,idate1,ihour1,iitime
integer :: iyear2,imonth2,iday2,ihour2,iminute2 

integer :: iseason

integer :: juday

real :: dt_cbmz


integer :: ifsm(5)

integer :: ifsmt

integer :: idmSet,ismMax

integer :: igMark(idmSet)


integer :: letdoit,MapS
real :: DeltSpeMark,OrgCon

integer :: i,j,k,is

integer :: mem4d
integer :: mem3d
integer :: mem2d


real,dimension(igas) :: GC_MOLWT

character*40,dimension(igas) :: GC_NAME 


integer :: ixy,i02,i03,i04

integer :: igg,ig,iduc,ia,i05,i05c

integer :: ictg,i02emt

integer :: idm,ism,i04sm

integer :: i05c1,i05c2,i05c3

integer :: i04_bc,I05_1,I05_2,I05_3,I05_4

integer :: i0517,i0518,i0519,i0520

integer :: ibin,i04aer

real :: rlon,rlat,RH,te,pr_atm,TER

real :: DUST01,DUST02,DUST03,DUST04
real :: SEA01,SEA02,SEA03,SEA04


real :: PSO4,BC,DUST_SO4,DUST_NO3,SSA_SO4,SSA_NO3

integer :: iwrongchem
integer :: ib

real :: deltcppb

real :: delta
real :: delta1,delta2,delta3,delta4

real :: cppb(ngas_max),cnnzifa(ngas_max)

real,dimension(nzz) :: ratioemitt,ratioemitPt

real,dimension(nzz) :: ratioem1,ratioem2,ratioem3,ratioem4,ratioem5,ratioem6

real :: FSO4_DUST(isize),FNO3_DUST(isize),FSO4_SSA(isize),FNO3_SSA(isize)


real :: contem0(idmSet),TmpSM(ismMax)

integer :: IHgtLev,iHgtLMax,ISrcDefined

real    :: daso4
integer :: nsp,noxdt
real,allocatable,dimension(:) :: conc,coxdt
real,allocatable,dimension(:,:) :: conc_oxdt


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!ihour2=4


call JDAY(iday2,imonth2,iyear2, juday, iseason)

i02 = ip2mem(ne)


!if(myid.eq.0) print*,'cbmz-time11:',iyear2,imonth2,iday2,ihour2,iminute2

!------------------calculate the begin time from noon March 21 (UTC)---------------
!------------ for 2004 is 81,but for other is 80,because the eb is 29 in 2004
    IF(juday>=81.) THEN
      IF(ihour2>= 12) then
        tbeg_dd = juday - 81              !day
        tbeg_hh = ihour2- 12              !hour
        tbeg_mm = iminute2                !min   
        tbeg_ss = 0.0                     !second
      ELSE
        tbeg_dd = juday - 81 -1
        tbeg_hh = ihour2- 12 +24
        tbeg_mm = iminute2 
        tbeg_ss = 0.0
      ENDIF
    ELSE
      IF(ihour2> 12) then
        tbeg_dd = juday - 81 +1
        tbeg_hh = ihour2- 12-24
        tbeg_mm = iminute2
        tbeg_ss = 0.0
      ELSE
        tbeg_dd = juday - 81 
        tbeg_hh = ihour2- 12 
        tbeg_mm = iminute2
        tbeg_ss = 0.0
      ENDIF
    ENDIF
!------------------------------------------------------------------------
!cccccccccccccccccccccccccccc   set the run time   cccccccccccccccc
      trun_dd = 0.0       ! days
      trun_hh = 0.0       ! hours
      trun_mm = dt_cbmz/60 !20.0      ! minutes 
!      trun_mm = 5.0 ! shun_apm_each
      trun_ss = 0.0       !seconds
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       dt_min = dt_cbmz/60  !20.0      ! the dt of chemcial 
!       dt_min = 5.0 ! (5min=300s) ! shun_apm_each
       msolar = 1         !(flag) 1 = diurnally varying photolysis; 2 = fixed phot
       mphoto = 2         !(flag) 1 = Old parameterization; 2 = New parameterization
       iprint = 0         !integer = freq of output. 
                          !Every iprint*dt_min mins.as the original,
                          !but for calculate ,it's not meaning
!  allocate(ratioemitt(nzz))
!  allocate(ratioemitPt(nzz))
!  allocate(FSO4_DUST(ISIZE),FNO3_DUST(ISIZE),FSO4_SSA(ISIZE),FNO3_SSA(ISIZE))

! Zifa 2006/10/21

!if(myid.eq.0) write(*,'(a,i2.2)') 'call CBMZ in domain ',ne


Gas_Chemistry: do j = sy(ne),ey(ne)
               do i = sx(ne),ex(ne)
                  
 ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

 do k=1,nzz-1 

   i03 = ip3mem(k,ne)

if(iopt_gaschem.eq.1) then

   do ig=1,igasCBM
     i04 = ip4mem(k,ig,ne)
     i02Gas= ip2memGas(ig,ne)
     cnn(ig)= MAX(gas(i04+ixy) , 1.E-20)
     cnnzifa(ig)= gas(i04+ixy)
    !ppbv the initial concentration as the box chemical model
     species(ig)=GC_NAME(ig)

     emission(ig)=0.0
     do ictg=1,NLAY_EM
       i02emt=ip_emit2d(ig,ictg,ne) 
       emission(ig) = emission(ig) + &
                    & ( emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg) &
                    &  *3600.0/dz(i03+ixy) ) &
                    &  *(0.08206*t(i03+ixy))/(Plev(i03+ixy)/1000.)/GC_MOLWT(ig)
     enddo

                           !ug/s/m2-->ppbv/h
    !!!!!!!!!!!!!!!
    ! for Source Mark
     if(ifsm(ne)==1)then
       letdoit=0
       do idm=1,idmSet
         if(igMark(idm)==ig)letdoit=idm
       enddo
       if(letdoit>0)then
         contem0(letdoit)=gas(i04+ixy)
       endif
    endif
    !!!!!!!!!!!!!!!

   enddo  !ig

            
   !zifa for real-time forecast 
   ! to save time
   if((cnn(5)+cnn(6)).ge.0.0 )then !NO+NO2 > 0.0

      rlon   =  longicrs(i02+ixy)   !rlat the box lat
      rlat   =  latitcrs(i02+ixy)   ! as the above ,but lon 
      zalt_m =  heiz(i03+ixy)       ! the altitude(asl) of box(m)
      RH     =  rh1(i03+ixy)        ! the RH
      te     =  t(i03+ixy)          ! the temp 
      pr_atm =  PA2ATM*Plev(i03+ixy)*100.  ! the pressure but as the atm
      TER    =  HGT1(i02+ixy)     

      ! for heterogeneous chemistry 
      i04_bc = ip4mem(k,77,ne)
      I05_1=IP5MEM(K,1,2,NE)
      I05_2=IP5MEM(K,2,2,NE)
      I05_3=IP5MEM(K,3,2,NE)
      I05_4=IP5MEM(K,4,2,NE)

      DUST01= AER(I05_1+IXY)
      DUST02= AER(I05_2+IXY)
      DUST03= AER(I05_3+IXY)
      DUST04= AER(I05_4+IXY)

      I05_1=IP5MEM(K,1,1,NE)
      I05_2=IP5MEM(K,2,1,NE)
      I05_3=IP5MEM(K,3,1,NE)
      I05_4=IP5MEM(K,4,1,NE)

      SEA01 = AER (I05_1+IXY)
      SEA02 = AER (I05_2+IXY)
      SEA03 = AER (I05_3+IXY)
      SEA04 = AER (I05_4+IXY)

      PSO4 = AMAX1(ASO4(i03+ixy), 1.E-20)   

if(laerv1) then
      PBC  = AMAX1(gas(i04_bc+ixy), 1.E-20)  
elseif(laerv2) then
      PBC=0.0
      do ibin=1,naerbin
         i04aer=ip4mem_aer(k,ibin,idx_bc,ne)
         PBC=PBC+aerom(i04aer+ixy)
      enddo
endif


      DUST_SO4 = 0.
      DUST_NO3 = 0. 
      SSA_SO4  = 0.
      SSA_NO3  = 0.


      DO IS = 1, ISIZE
           i05c  =  ip5memcs (k,is,7,ne) ! SO4 on Ses Salt
           i05c1 =  ip5memcs (k,is,8,ne) ! NO3 on Sea Salt
           i05c2 =  ip5memc  (k,is,7,ne) ! SO4 on DUST
           i05c3 =  ip5memc  (k,is,8,ne) ! NO3 on DUST
           DUST_SO4 = DUST_SO4 + DUSTCOMP(i05c2+ixy)
           DUST_NO3 = DUST_NO3 + DUSTCOMP(i05c3+ixy)
           SSA_SO4  = SSA_SO4  + SEACOMP(i05c+ixy)
           SSA_NO3  = SSA_NO3  + SEACOMP(i05c1+ixy)
      ENDDO

      ! CALCULATE THE RK_HET FOR HETEROGENEOUS CHEMISTRY
      CALL HETERO_REAC( PSO4, PBC &
                       ,DUST01,DUST02,DUST03,DUST04 &
                       ,SEA01,SEA02,SEA03,SEA04 &
                       ,DUST_SO4, DUST_NO3, SSA_SO4,SSA_NO3 &
                       ,FSO4_DUST,FNO3_DUST,FSO4_SSA,FNO3_SSA )

      !!!!!!!

     FCLD = coefcld3d(i03+ixy)
!      FCLD = 1.0

      cnn(kdso4)=0.0  ! set intial value of dso4(heterogeneous) to 0.0 
      cnn(kdno3)=0.0

!C****   TO DO SENSITIVITY ANALYSIS TEST ****
!C *****
!C         RK_HET = 0.0
!C**************** SENSITIVITY ANALYSIS END

      RK_HETSO2_DUST  ( i03+ixy) = RK_HET(19) ! HETEROGENEOUS SO2 on dust
      RK_HETHNO3_DUST ( i03+ixy) = RK_HET(12) ! HETEROGENEOUS HNO3 on dust

      call cbmz(cppb,FCLD)
      jo1d(i03+ixy) = rk_photo(jphoto_o3b) 
      jno2(i03+ixy) = rk_photo(jphoto_no2)
          
         
      call chemope(myid, OPE(i03+ixy), trun_mm,i,j,k,ne)

! allocate the so4 and no3 from heterogeneous on dust and ssa 4 bins
      do is = 1, isize
        i05c  =  ip5memcs (k,is,7,ne) ! SO4 on Ses Salt
        i05c1 =  ip5memcs (k,is,8,ne) ! NO3 on Sea Salt
        i05c2 =  ip5memc  (k,is,7,ne) ! SO4 on DUST
        i05c3 =  ip5memc  (k,is,8,ne) ! NO3 on DUST

        CALL ALLD  ( MYID,cppb(kdso4),cppb(kdno3) &
                    ,SEACOMP(i05c+ixy),SEACOMP(i05c1+ixy) &
                    ,DUSTCOMP(i05c2+ixy),DUSTCOMP(i05c3+ixy) &
                    ,FSO4_SSA(is),FNO3_SSA(is),FSO4_DUST(is),FNO3_DUST(is) )
      enddo ! is  

!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(lprocess) then
      do ig=1,iPrintTermGas    ! for gas phase
        igg=IGGPOS(ig)
        if(igg==11) then  !for ozone
          i0517=ipGasTermBal(k,17,ig,ne) ! 17: production 11: ozone
          i0518=ipGasTermBal(k,18,ig,ne) ! 18: loss       11: ozone 
          i0519=ipGasTermBal(k,19,ig,ne) ! 19: radicals LOSS due to NOx 11:ozone
          i0520=ipGasTermBal(k,20,ig,ne) ! 20: radicals LOSS due to VOCs 11:ozone        

          call chemprodloss(myid,delta1,delta2,delta3,delta4,trun_mm,i,j,k,ne)
          GasTermBal(i0517+ixy)=delta1+GasTermBal(i0517+ixy)
          GasTermBal(i0518+ixy)=delta2+GasTermBal(i0518+ixy)
          GasTermBal(i0519+ixy)=delta3+GasTermBal(i0519+ixy)
          GasTermBal(i0520+ixy)=delta4+GasTermBal(i0520+ixy)

        endif !igg
      enddo
endif         
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
                         
          !  to judge the change ratio
      iwrongchem=0
      ib=11  ! o3 
      deltcppb=abs(cppb(ib)-cnnzifa(ib))
      if(deltcppb.ge.100)iwrongchem=1 ! > 50ppb 
         ib=18  ! so2 
         deltcppb=abs(cppb(ib)-cnnzifa(ib))/(cnnzifa(ib)+0.00001)
         if(deltcppb.ge.0.10)iwrongchem=1 ! >  10% decrease
         if(iwrongchem==1)goto 4190  ! not to use the output of CBZ

         do ig=1,igasCBM
            i04 = ip4mem(k,ig,ne)
            gas(i04+ixy) = cppb(ig)   !cppb the output ppbv 

            call chemprod(ig,i,j,k,delta,trun_mm,ne) !lijie to get the
!                                                   ozone  production 
              
            !!!!!!!!!!!!!!!
            ! for Source Mark
            if(ifsm(ne)==1)then
              letdoit=0
              do idm=1,idmSet
                if(igMark(idm)==ig)letdoit=idm
              enddo
              if(letdoit>0)then
                i02 = ip2mem(ne)
                MapS=int(MapSource(i02+ixy))
                IHgtLev= iHgtLMax-1  ! for chemistry reaction part Zifa/2007/03/03
                   ! this is need to be defined outside for future   
                i04 = ip4mem(k,ig,ne)
                OrgCon     = contem0(letdoit)
                if(ig==11) then
                   DeltSpeMark=delta
                else  
                   DeltSpeMark= gas(i04+ixy)-OrgCon
                endif
                if(DeltSpeMark > 0. )then
                  do ism=1,ismMax
                    i04sm=ipSMmem(k,ism,letdoit,ne)
                    TmpSM(ism)=SourceMark(i04sm+ixy)
                  enddo
                  call GetSMarkChange(DeltSpeMark,OrgCon,TmpSM,MapS, &
                  ISrcDefined,IHgtLev,ismMax)
                  do ism=1,ismMax
                    i04sm=ipSMmem(k,ism,letdoit,ne)
                    SourceMark(i04sm+ixy)=TmpSM(ism)
                  enddo
                endif
              endif
            endif
            !!!!!!!!!!!!!!!

          enddo !ig
 4190     continue
    endif  ! zifa

elseif(iopt_gaschem.eq.2) then! simplidied gas chemistry

      rlon   =  longicrs(i02+ixy)   !rlat the box lat
      rlat   =  latitcrs(i02+ixy)   ! as the above ,but lon 
      zalt_m =  heiz(i03+ixy)       ! the altitude(asl) of box(m)
      RH     =  rh1(i03+ixy)        ! the RH
      te     =  t(i03+ixy)          ! the temp 
      pr_atm =  PA2ATM*Plev(i03+ixy)*100.  ! the pressure but as the atm
      TER    =  HGT1(i02+ixy)

if(1==2) then
rlon=120
rlat=36
zalt_m=50.0
RH=60.0
te=300
pr_atm =PA2ATM*900*100.
ter=10
endif
      nsp=igassmp
      if(.not.allocated(conc)) allocate(conc(nsp))

      noxdt=nmoxdt
      if(.not.allocated(conc_oxdt)) allocate(conc_oxdt(noxdt,nm12))

      if(.not.allocated(coxdt)) allocate(coxdt(noxdt))

      do im=1,nm12
      do ig=1,noxdt
         i04=ip4mem_oxdt(k,ig,im,ne)
         conc_oxdt(ig,im)=mmean_oxdt(i04+ixy) ! ppb
      enddo
      enddo

!print*,'conc_oxdt=',conc_oxdt


      do ig=1,nsp
        i04=ip4mem(k,ig,ne)
        conc(ig)=gas(i04+ixy) ! ppb
!conc(ig)=10
      enddo

!      print*,'call sulfur_smp_gchembox',ne
      call sulfur_smp_gchembox( myid,dt_cbmz,i,j,k &
                               ,iyear2, imonth2, iday2, juday, iseason &
                               ,ihour2, iminute2 &
                               ,tbeg_dd,tbeg_hh,tbeg_mm,tbeg_ss &
                               ,rlon,rlat,zalt_m &
                               ,nsp,noxdt,conc_oxdt,conc,daso4,te,pr_atm,rh &
                               ,coxdt)


      tmp4d(i03+ixy)=coxdt(idx_oxdt_oh)


!if(i.eq.180.and.j.eq.90) then
! print*,'kk oxdt',coxdt
! stop
!endif

      do ig=1,noxdt
        i04=ip4mem_ox3d(k,ig,ne)
        oxdt3d(i04+ixy)=coxdt(ig)
      enddo

      do ig=1,nsp
        i04=ip4mem(k,ig,ne)
        gas(i04+ixy)=conc(ig)
      enddo

      i04aer=ip4mem_aer(k,1,idx_so4,ne)
      aerom(i04aer+ixy)=aerom(i04aer+ixy)+daso4

!      call simple_gaschem

else
  print*,'iopt_gaschem err set',iopt_gaschem
  stop
endif

   enddo !k
 enddo !i
 enddo Gas_Chemistry  !j

! deallocate(ratioemitt)
! deallocate(ratioemitPt)
! deallocate(FSO4_DUST,FNO3_DUST,FSO4_SSA,FNO3_SSA)



!endif ! every 20 minutes
!endif ! chem calculation flag yes/not

!> end of gas phase chemistry
!===============================================================

end subroutine naqpms_drv_gaschem
