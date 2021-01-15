
subroutine naqpms_drv_soap &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,NSOA &
 &  ,ifsm,idmSet,ismMax,igMark )

use naqpms_varlist
use met_fields, only : t,Plev
use naqpms_gridinfo
implicit none

integer :: myid

real :: dt

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom,NSOA

integer :: i,j,k,is
!integer :: k,is

integer :: mem3d

integer :: mem2d
real,dimension(mem2d) :: ktop


integer :: ixy,i02,i03,i04

integer :: ig,iduc,ia,i05,i05c


integer :: i04_1,i04_2,i04_3,i04_4,i04_5,i04_6,i04_7 &
          ,i04_8,i04_9,i04_10,i04_11,i04_12,i04_13

real :: TEMP,PRESS0,POA

real :: SVG(NSOA),SOA(NSOA)


integer :: ifsm(5)


integer :: idmSet,ismMax

integer :: igMark(idmSet)

integer :: letdoit,idm,ism,i04sm

integer :: i04aer,ibin

real :: contem0(idmSet)
real :: TmpSM(ismMax)

integer :: IHgtLev,iHgtLMax,ISrcDefined

integer :: MapS
real :: DeltSpeMark,OrgCon

real :: frcbin(naerbin)

frcbin(1)=1.0
frcbin(2)=0.0


 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
 
   ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

   do k=1,nzz-1

     i03 = ip3mem(k,ne)


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


        TEMP     =  t(i03+ixy )    ! THE TEMPERATURE IN K
        PRESS0   =  Plev(i03+ixy ) ! THE PRESSURE IN HPA

        i04_1  = ip4mem(k,69, ne)    ! SV1
        i04_2  = ip4mem(k,70, ne)    ! SV2
        i04_3  = ip4mem(k,71, ne)    ! SV3
        i04_4  = ip4mem(k,72, ne)    ! SV4
        i04_5  = ip4mem(k,73, ne)    ! SV5
        i04_6  = ip4mem(k,74, ne)    ! SV6

if(laerv1) then
        i04_7  = ip4mem(k,96, ne)    ! SOA1
        i04_8  = ip4mem(k,97, ne)    ! SOA2
        i04_9  = ip4mem(k,98, ne)    ! SOA3
        i04_10 = ip4mem(k,99, ne)    ! SOA4
        i04_11 = ip4mem(k,100,ne)    ! SOA5
        i04_12 = ip4mem(k,101,ne)    ! SOA6
        i04_13 = ip4mem(k,78,ne)     ! POA
        POA    = gas(i04_13+ixy)     ! POA in ug/m3
endif


        SVG(1) = gas(i04_1 +ixy)     ! SV1
        SVG(2) = gas(i04_2 +ixy)     ! SV2
        SVG(3) = gas(i04_3 +ixy)     ! SV3
        SVG(4) = gas(i04_4 +ixy)     ! SV4
        SVG(5) = gas(i04_5 +ixy)     ! SV5
        SVG(6) = gas(i04_6 +ixy)     ! SV6

if(laerv1) then
        SOA(1) = gas(i04_7 +ixy)     ! SOA1
        SOA(2) = gas(i04_8 +ixy)     ! SOA2
        SOA(3) = gas(i04_9 +ixy)     ! SOA3
        SOA(4) = gas(i04_10+ixy)     ! SOA4
        SOA(5) = gas(i04_11+ixy)     ! SOA5
        SOA(6) = gas(i04_12+ixy)     ! SOA6
endif


if(laerv2) then

     POA=0.0
     do ibin=1,naerbin
        i04aer=ip4mem_aer(k,ibin,idx_oc,ne)
        POA=POA+aerom(i04aer+ixy)
     enddo

      SOA=0.0
      do ibin=1,naerbin
          i04aer=ip4mem_aer(k,ibin,idx_soa01,ne)
          SOA(1)=SOA(1)+aerom(i04aer+ixy)
          i04aer=ip4mem_aer(k,ibin,idx_soa02,ne)
          SOA(2)=SOA(2)+aerom(i04aer+ixy)
          i04aer=ip4mem_aer(k,ibin,idx_soa03,ne)
          SOA(3)=SOA(3)+aerom(i04aer+ixy)
          i04aer=ip4mem_aer(k,ibin,idx_soa04,ne)
          SOA(4)=SOA(4)+aerom(i04aer+ixy)
          i04aer=ip4mem_aer(k,ibin,idx_soa05,ne)
          SOA(5)=SOA(5)+aerom(i04aer+ixy)
          i04aer=ip4mem_aer(k,ibin,idx_soa06,ne)
          SOA(6)=SOA(6)+aerom(i04aer+ixy)
      enddo
endif



      CALL SOAP(NSOA, SOA, SVG, TEMP, PRESS0, POA, I, J, K)

       gas(i04_1 +ixy) = SVG(1)
       gas(i04_2 +ixy) = SVG(2)
       gas(i04_3 +ixy) = SVG(3)
       gas(i04_4 +ixy) = SVG(4)
       gas(i04_5 +ixy) = SVG(5)
       gas(i04_6 +ixy) = SVG(6)

if(laerv1) then
       gas(i04_7 +ixy) = SOA(1)
       gas(i04_8 +ixy) = SOA(2)
       gas(i04_9 +ixy) = SOA(3)
       gas(i04_10+ixy) = SOA(4)
       gas(i04_11+ixy) = SOA(5)
       gas(i04_12+ixy) = SOA(6)
       gas(i04_13+ixy) = POA
endif


if(laerv2) then
! POA ???
      do ibin=1,naerbin
        i04aer=ip4mem_aer(k,ibin,idx_soa01,ne)
        aerom(i04aer+ixy)=SOA(1)*frcbin(ibin)
        i04aer=ip4mem_aer(k,ibin,idx_soa02,ne)
        aerom(i04aer+ixy)=SOA(2)*frcbin(ibin)
        i04aer=ip4mem_aer(k,ibin,idx_soa03,ne)
        aerom(i04aer+ixy)=SOA(3)*frcbin(ibin)
        i04aer=ip4mem_aer(k,ibin,idx_soa04,ne)
        aerom(i04aer+ixy)=SOA(4)*frcbin(ibin)
        i04aer=ip4mem_aer(k,ibin,idx_soa05,ne)
        aerom(i04aer+ixy)=SOA(5)*frcbin(ibin)
        i04aer=ip4mem_aer(k,ibin,idx_soa06,ne)
        aerom(i04aer+ixy)=SOA(6)*frcbin(ibin)
      enddo
endif


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


    enddo ! k
 enddo    ! i
 enddo    ! j
!!!!!!!!!




end subroutine naqpms_drv_soap





