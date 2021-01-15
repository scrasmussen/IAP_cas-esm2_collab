
subroutine naqpms_drv_isorropia &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,mem3d &
 &  ,GC_MOLWT &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,ifsm,idmSet,ismMax,igMark )

use naqpms_varlist
use met_fields, only : t,rh1,Plev
implicit none

real,parameter :: wgt_h2so4=98,wgt_nh3=17,wgt_hno3=63,wgt_hcl=36.5
real,parameter :: wgt_so4=96,wgt_nh4=18,wgt_no3=62,wgt_na=23,wgt_cl=35.5

integer :: myid

real :: dt

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: ifsm(5)

integer :: idmSet,ismMax

integer :: igMark(idmSet)


integer :: letdoit,MapS
real :: DeltSpeMark,OrgCon

integer :: i,j,k,is


integer :: mem3d
integer :: mem2d


real,dimension(igas) :: GC_MOLWT


integer :: ixy,i02,i03,i04

integer :: ig,iduc,ia,i05,i05c

integer :: i04aer,ibin

integer :: idm,ism,i04sm

integer :: i04_1,i04_2,i04_3,i04_4,i04_5,i04_6,i04_7 &
          ,i04_8,i04_9,i04_10,i04_11,i04_12,i04_13,i04_14 &
          ,i04_15,i04_16,i04_17,i04_18,i04_19,i04_20,i04_21 &
          ,i04_22

DOUBLE PRECISION                :: TEMPI,RHI
DOUBLE PRECISION,DIMENSION(5)   :: WI
DOUBLE PRECISION, DIMENSION(17) :: WO
DOUBLE PRECISION,DIMENSION(4)   :: WOG
DOUBLE PRECISION ::  AWATER

real :: contem0(idmSet)
real :: TmpSM(ismMax)

integer :: IHgtLev,iHgtLMax,ISrcDefined


integer,parameter :: nspwo=17
real,parameter :: wgt_wo(1:nspwo)=(/ 1.0D0 , 23.0D0, 18.0D0, 35.5D0, 96.0D0 &
                                   ,97.0D0, 62.0D0, 58.5D0, 142.0D0,85.0D0 &
                                   ,132.0D0,80.0D0, 53.5D0, 98.0D0, 115.0D0 &
                                   ,120.0D0,247.0D0 /)

!return

! Debug


 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)

    ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

    do k=1,nzz-1


    !!!
    !   WI(1) : TOTAL SODIUM   AS EQUIVALENT NA
    !   WI(2) : TOTAL SULFATE  AS EQUIVALENT H2SO4
    !   WI(3) : TOTAL AMMONIUM AS EQUIVALENT NH3
    !   WI(4) : TOTAL NITRATE  AS EQUIVALENT HNO3
    !   WI(5) : TOTAL CHLORIDE AS EQUIVALENT HCL
    !   WO :  17 AEROSOLS SPECIES
    !   WOG:  4  GASEOUS  SPECIES

        i03 = ip3mem(k,ne)
        TEMPI     =  DBLE(t(i03+ixy ))  ! THE TEMPERATURE IN K
        RHI       =  DBLE(rh1(i03+ixy)) ! THE RH IN 0-100%   


    !!!!!!!!!!!!!!!
    ! for SourceMark
    do ig=1, igas
      i04 = ip4mem(k,ig,ne)
      if(ifsm(ne)==1)then
        letdoit=0
        do idm=1,idmSet
         if(igMark(idm)==ig)letdoit=idm
        enddo
        if(letdoit>0)then
          contem0(letdoit)=gas(i04+ixy)
        endif
     endif! ifsm
    enddo  !ig
    !!!!!!!!!!!!!!!


!   !!! CONVERT GAS IN PPBV TO UG/M3

    do ig=1,4
       i04 = ip4mem(k,ig,ne)
       gas(i04+ixy) = gas(i04+ixy)* &
                    & GC_MOLWT(ig)/((0.08206*T(i03+ixy))/(Plev(i03+ixy)/1000.))
    enddo !ig 
!
    !!! GET WI(5)
    i04_1 = ip4mem(k,1,ne)   ! H2SO4(G)
    i04_2 = ip4mem(k,2,ne)   ! HNO3(G)
    i04_3 = ip4mem(k,3,ne)   ! HCL(G)
    i04_4 = ip4mem(k,4,ne)   ! NH3(G)

if(laerv1) then
    i04_5 = ip4mem(k,80,ne)  ! Na+(AE)
    i04_6 = ip4mem(k,81,ne)  ! NH4+(AE)
    i04_7 = ip4mem(k,82,ne)  ! CL-(AE)
    i04_8 = ip4mem(k,83,ne)  ! SO4--(AE)
    i04_9 = ip4mem(k,84,ne)  ! HSO4-(AE)
    i04_10= ip4mem(k,85,ne)  ! NO3-(AE)
    i04_11= ip4mem(k,86,ne)  ! NACL(S)
    i04_12= ip4mem(k,87,ne)  ! NA2SO4(S)
    i04_13= ip4mem(k,88,ne)  ! NANO3(S)
    i04_14= ip4mem(k,89,ne)  ! (NH4)2SO4(S)
    i04_15= ip4mem(k,90,ne)  ! NH4NO3(S)
    i04_16= ip4mem(k,91,ne)  ! NH4CL(S)
    i04_17= ip4mem(k,92,ne)  ! H2SO4(AQ)
    i04_18= ip4mem(k,93,ne)  ! NH4HSO4(S)
    i04_19= ip4mem(k,94,ne)  ! NAHSO4(S)
    i04_20= ip4mem(k,95,ne)  ! (NH4)3H(SO4)2(S)
    i04_21= ip4mem(k,79,ne)  !  H+(AE)
    i04_22= ip4mem(k,102,ne) ! AH2O 
endif

    WI    = 0.0
    WO    = 0.0
    WOG   = 0.0


        WI(1) = DBLE( ANA(i03+ixy) )

        WI(2) = DBLE( ASO4(i03+ixy)*wgt_h2so4/wgt_so4  &
                    + gas(i04_1+ixy)*wgt_h2so4/GC_MOLWT(1) )

        WI(3) = DBLE( ANH4(i03+ixy)*wgt_nh3/wgt_nh4 &
                    + gas(i04_4+ixy)*wgt_nh3/GC_MOLWT(4) )

        WI(4) = DBLE( ANO3(i03+ixy)*wgt_hno3/wgt_no3 &
                    + gas(i04_2+ixy)*wgt_hno3/GC_MOLWT(2) )

        WI(5) = DBLE( ACL(i03+ixy)*wgt_hcl/wgt_cl &
                    + gas(i04_3+ixy)*wgt_hcl/GC_MOLWT(3) )


!    WI(1)  = AMAX1(WI(1), 1.E-2)
!    WI(5)  = AMAX1(WI(1), 1.E-2)

   IF(WI(2).LE.0.1) GOTO 2880  ! shun comment out
   IF(WI(3).LE.0.1) GOTO 2880  ! shun comment out
   IF(WI(4).LE.0.1) GOTO 2880  ! shun comment out


!*****   TO CALL ISRPINTR ******
!    CALL ISRPINTR(MYID,WI,RHI,TEMPI,WO,WOG(1),WOG(2),WOG(3),WOG(4),I,J,K)
! shun@albany_20140708
    CALL ISRPINTR(MYID,WI,RHI,TEMPI,WO,WOG(1),WOG(2),WOG(3),WOG(4),AWATER,I,J,K)

    gas(i04_4+ixy)  = SNGL(WOG(1)) ! NH3(G)
    gas(i04_2+ixy)  = SNGL(WOG(2)) ! HNO3(G)
    gas(i04_1+ixy)  = SNGL(WOG(3)) ! H2SO4(G)
    gas(i04_3+ixy)  = SNGL(WOG(4)) ! HCL(G)


if(laerv1) then
    gas(i04_21+ixy) = SNGL(WO(1))  ! H+(AQ)
    gas(i04_5+ixy)  = SNGL(WO(2))  ! NA+(AQ)
    gas(i04_6+ixy)  = SNGL(WO(3))  ! NH4+(AQ)
    gas(i04_7+ixy)  = SNGL(WO(4))  ! CL-(AQ)
    gas(i04_8+ixy)  = SNGL(WO(5))  ! SO4--(AQ)
    gas(i04_9+ixy)  = SNGL(WO(6))  ! HSO4-(AQ)
    gas(i04_10+ixy) = SNGL(WO(7))  ! NO3-(AQ)
    gas(i04_11+ixy) = SNGL(WO(8))  ! NACL(S)
    gas(i04_12+ixy) = SNGL(WO(9))  ! NA2SO4(S)
    gas(i04_13+ixy) = SNGL(WO(10)) ! NANO3(S)
    gas(i04_14+ixy) = SNGL(WO(11)) ! NH42SO4(S)
    gas(i04_15+ixy) = SNGL(WO(12)) ! NH42SO4(S)
    gas(i04_16+ixy) = SNGL(WO(13)) ! NH4CL(S)
    gas(i04_17+ixy) = SNGL(WO(14)) ! H2SO4(AQ)
    gas(i04_18+ixy) = SNGL(WO(15)) ! NH4HSO4(S)
    gas(i04_19+ixy) = SNGL(WO(16)) ! NAHSO4
    gas(i04_20+ixy) = SNGL(WO(17)) ! (NH4)4H(SO4)2(S)
    gas(i04_22+IXY) = SNGL(AWATER)

elseif(laerv2) then
! update model aerosol tracers : nitrate and ammonium
    ibin=1
    i04aer=ip4mem_aer(k,ibin,idx_so4,ne)
    aerom(i04aer+ixy) = wo(5)*wgt_so4/wgt_wo(5) &
                      + wo(6)*wgt_so4/wgt_wo(6) &
                      + wo(9)*wgt_so4/wgt_wo(9) &
                      + wo(11)*wgt_so4/wgt_wo(11) &
                      + wo(14)*wgt_so4/wgt_wo(14) &
                      + wo(15)*wgt_so4/wgt_wo(15) &
                      + wo(16)*wgt_so4/wgt_wo(16) &
                      + wo(17)*wgt_so4/wgt_wo(17)*2.0 ! SO4
    ibin=1
    i04aer=ip4mem_aer(k,ibin,idx_nh4,ne)
    aerom(i04aer+ixy) = wo(3)*wgt_nh4/wgt_wo(3)        &
                      + wo(11)*wgt_nh4/wgt_wo(11)*2.0   &
                      + wo(12)*wgt_nh4/wgt_wo(12)       &
                      + wo(13)*wgt_nh4/wgt_wo(13)       &
                      + wo(15)*wgt_nh4/wgt_wo(15)       &
                      + wo(17)*wgt_nh4/wgt_wo(17)*3.0

    ibin=1
    i04aer=ip4mem_aer(k,ibin,idx_no3,ne)
    aerom(i04aer+ixy) = wo(7)*wgt_no3/wgt_wo(7)   &
                      + wo(10)*wgt_no3/wgt_wo(10)  &
                      + wo(12)*wgt_no3/wgt_wo(12)

    awc3d(i03+ixy)=SNGL(AWATER)

endif


2880  CONTINUE


if(laerv1) then
     gas(i04_17+ixy) = gas(i04_17+ixy) + gas(i04_1+ixy)
elseif(laerv2) then
    ibin=1
    i04aer=ip4mem_aer(k,ibin,idx_so4,ne)
    aerom(i04aer+ixy)=aerom(i04aer+ixy)+gas(i04_1+ixy)*wgt_so4/wgt_h2so4
endif



     gas(i04_1+ixy)=0.0


    if(k.eq.1) then
      !print*,'acid in main',i,j,WOG(3),gas(i04_1+ixy)
    endif

    !!!  CONVERT GAS IN UG/M3 TO PPBV
    do ig=1,4
       i04 = ip4mem(k,ig,ne)
       gas(i04+ixy) = gas(i04+ixy)*                          &
         (0.08206*T(i03+ixy))/(Plev(i03+ixy)/1000.)/GC_MOLWT(ig)
    enddo    !ig 


    !!!!!!!!!!!!!!!
    ! for Source Mark
    DO ig=1,igas

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
          DeltSpeMark= gas(i04+ixy)-OrgCon
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
        endif ! letdoit
      endif  ! ifsm

    ENDDO !IG
    !!!!!!


    enddo  !k

  enddo    !i
  enddo    !j


!stop 'end_iso'

end subroutine naqpms_drv_isorropia





